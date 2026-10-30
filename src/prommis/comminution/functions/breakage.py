#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Breakage (daughter-distribution) functions and the breakage-matrix builder.

Breakage values describe the fraction of broken parent-bin-``j`` mass finer
than a daughter size. The Austin form uses two daughter edges; its second term
uses the next-coarser edge. Forms (ascending-edge notation):

* ``luckie_austin`` (default): ``Phi*(x_k/d_j)^gamma + (1-Phi)*(x_{k+1}/d_j)^beta``
  (capped at 1 by ``breakage_matrix``)
* ``reid_stewart``: ``Phi*(x_k/d_j)^gamma + (1-Phi)*(x_k/d_j)^beta`` (same edge)
* ``vogel``: ``(x/d_j)^q * 1/2 (1 + tanh((x - d')/d'))``

The matrix builders turn cumulative values into the per-bin breakage
matrix ``b[i][j]`` (daughter ``i`` from parent ``j``). They assign the
fraction not entering any finer bin to the parent bin
(``b_jj = 1 - B_j(x_j)``) and assign material below the lowest mesh edge
to bin 0 (``b_0j = B_j(x_1)``). These assignments ensure that each column
sums to 1.

For parent bins ``j > 0``, ``vogel`` evaluates ``B_j`` at the daughter
characteristic sizes ``d_i`` instead of the edges. Its end-bin formulas are
``b_jj = 1 - B_j(d_{j-1})`` and ``b_0j = B_j(d_0)``.

References:

[1] Reid, K. J. (1965). A solution to the batch grinding equation.
    Chem. Eng. Sci., 20(11), 953-963.

[2] Austin, L. G., Klimpel, R. R., & Luckie, P. T. (1984). Process
    Engineering of Size Reduction: Ball Milling. SME-AIME, New York.

[3] Narayanan, S. S., & Whiten, W. J. (1988). Determination of comminution
    characteristics from single-particle breakage tests and its application to
    ball-mill scale-up. Trans. IMM Sect. C, 97, C115-C124.

[4] Vogel, L., & Peukert, W. (2003). Breakage behaviour of different
    materials - construction of a mastercurve for the breakage probability.
    Powder Technology, 129(1-3), 101-110.

[5] Napier-Munn, T. J., Morrell, S., Morrison, R. D., & Kojovic, T. (1996).
    Mineral Comminution Circuits: Their Operation and Optimisation. JKMRC,
    Univ. of Queensland.
"""

__author__ = "Daison Yancy Caballero"

import math
from collections.abc import Mapping
from scipy.interpolate import CubicSpline, PchipInterpolator


from prommis.comminution.core.config_utils import _is_real_number

# Allowed decrease between neighboring points of an assembled t10 curve.
_CURVE_MONOTONE_TOL = 1e-12
# The t-family markers give percent passing at fractions of the parent
# characteristic size; t10 is supplied separately.
_T10_MARKER_RATIOS = {
    "t75": 1.0 / 75,
    "t50": 1.0 / 50,
    "t25": 1.0 / 25,
    "t4": 0.25,
    "t2": 0.5,
}


def luckie_austin_cumulative(x_upper, x_coarser, d_j, phi, gamma, beta):
    """Austin/Luckie-Austin breakage value with the next-coarser edge in term two.

    ``B = Phi*(x_upper/d_j)^gamma + (1 - Phi)*(x_coarser/d_j)^beta``; all sizes
    share one length unit. See Austin et al. (1984).
    """
    return phi * (x_upper / d_j) ** gamma + (1.0 - phi) * (x_coarser / d_j) ** beta


def reid_stewart_cumulative(x, d_j, phi, gamma, beta):
    """Reid-Stewart cumulative breakage using the same edge in both terms.

    ``B = Phi*(x/d_j)^gamma + (1 - Phi)*(x/d_j)^beta``; ``x`` and ``d_j`` share
    one length unit. See Reid (1965) and Austin et al. (1984).
    """
    return phi * (x / d_j) ** gamma + (1.0 - phi) * (x / d_j) ** beta


def vogel_cumulative(x, d_j, q, dprime):
    """Vogel cumulative breakage (Vogel & Peukert, 2003).

    ``(x/d_j)^q * 1/2 (1 + tanh((x - d')/d'))``.

    All sizes share one length unit. ``dprime`` (``d'``) is an absolute length
    in the unit of ``x``, applied identically to every parent column; unlike the
    parent-relative ``(x/d_j)^q`` term it does not scale with the parent
    characteristic size ``d_j``.
    """
    return (x / d_j) ** q * 0.5 * (1.0 + math.tanh((x - dprime) / dprime))


def _finite(v, what):
    """Convert a finite, non-boolean ``numbers.Real`` value to ``float``.

    Raise ``ValueError`` for other inputs.
    """
    if not _is_real_number(v):
        raise ValueError(f"{what} must be a real number (got {type(v).__name__})")
    try:
        f = float(v)
    except OverflowError as exc:
        raise ValueError(f"{what} must be convertible to a float") from exc
    if not math.isfinite(f):
        raise ValueError(f"{what} must be a finite number (got {v!r})")
    return f


def t10_appearance_curve(appearance_table, t10):
    """Build the cumulative breakage curve ``B(ratio)`` from a t10 appearance table.

    JKMRC t10 / appearance-function method (Narayanan & Whiten, 1988; Napier-Munn
    et al., 1996).

    ``appearance_table`` maps ``T10`` rows (e.g. "10"/"20"/"30") to ``t75/t50/
    t25/t4/t2`` percent-passing values at size ratios ``1/75 ... 1/2`` of the
    parent. Interpolate each marker across the ``t10`` rows with a monotone
    cubic. Insert the supplied ``t10`` (%) directly at ratio ``1/10``, then
    build the size-ratio curve anchored at ``(0, 0)`` and ``(1, 1)``.

    ``appearance_table`` must be a dict with at least two rows. String row keys
    are converted with ``int()``; other row keys must equal their integer value.
    Keys must remain unique after conversion. Each row is a dict with exactly
    the five markers. ``t10`` and the marker values must be finite real
    numbers.

    Raises:
        ValueError: on a malformed table (non-dict table or row, invalid or
            duplicate row keys, fewer than two rows, missing or unknown
            marker), a non-real or non-finite ``t10``/marker value, a
            ``t10`` outside the table rows (no extrapolation), or a non-monotone
            assembled curve.
    """
    if not isinstance(appearance_table, dict):
        raise ValueError("appearance_table must be a dict")
    # Normalize string and integer row keys; reject keys identifying the same row.
    row_by_int = {}
    try:
        for k in appearance_table:
            r = int(k)
            if not isinstance(k, str) and r != k:
                raise ValueError("non-integral row key")
            row_by_int[r] = appearance_table[k]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("appearance_table row keys must be integer-like") from exc
    if len(row_by_int) != len(appearance_table):
        raise ValueError(
            "t10 appearance table has duplicate integer-equivalent row keys "
            "(e.g. '10' and 10)."
        )
    rows = sorted(row_by_int)
    if len(rows) < 2:
        raise ValueError("appearance_table needs at least two rows to interpolate")
    for r in rows:
        row = row_by_int[r]
        if not isinstance(row, dict):
            raise ValueError(f"appearance_table row {r} must be a dict")
        for marker in _T10_MARKER_RATIOS:
            if marker not in row:
                raise ValueError(
                    f"appearance_table row {r} is missing marker {marker!r}"
                )
        unknown = sorted((m for m in row if m not in _T10_MARKER_RATIOS), key=str)
        if unknown:
            raise ValueError(
                f"appearance_table row {r} has unknown markers {unknown!r}; "
                f"allowed markers: {list(_T10_MARKER_RATIOS)}"
            )
    t10 = _finite(t10, "t10")
    if t10 < rows[0] or t10 > rows[-1]:
        raise ValueError(
            f"t10={t10} is outside the appearance-table rows {rows} "
            "(extrapolation is not supported)."
        )
    points = {0.0: 0.0, 1.0 / 10: t10 / 100.0, 1.0: 1.0}
    for marker, ratio in _T10_MARKER_RATIOS.items():
        ys = [
            _finite(row_by_int[r][marker], f"row {r} marker {marker!r}") for r in rows
        ]
        points[ratio] = float(PchipInterpolator(rows, ys)(t10)) / 100.0
    ratios = sorted(points)
    cums = [points[r] for r in ratios]
    if any(cums[i + 1] < cums[i] - _CURVE_MONOTONE_TOL for i in range(len(cums) - 1)):
        raise ValueError("t10 appearance curve is non-monotone in size ratio.")
    pchip = PchipInterpolator(ratios, cums)

    def curve(ratio):
        return float(min(1.0, max(0.0, pchip(min(ratio, 1.0)))))

    return curve


def _cumulative(name, x_upper, x_coarser, d_j, params):
    if name == "luckie_austin":
        return luckie_austin_cumulative(
            x_upper, x_coarser, d_j, params["phi"], params["gamma"], params["beta"]
        )
    if name == "reid_stewart":
        return reid_stewart_cumulative(
            x_upper, d_j, params["phi"], params["gamma"], params["beta"]
        )
    if name == "vogel":
        return vogel_cumulative(x_upper, d_j, params["q"], params["dprime"])
    raise ValueError(f"unknown breakage function {name!r}")


def breakage_matrix(name, edges, d_char, params):
    """Build the mass-conserving per-bin breakage matrix (column ``j`` -> daughters).

    Args:
        name: ``"luckie_austin"``, ``"reid_stewart"``, ``"vogel"``, or
            ``"t10_appearance"``.
        edges: ``N+1`` ascending size edges.
        d_char: ``N`` characteristic sizes.
        params: parameter dict for the named function.

    Returns:
        ``N x N`` nested list ``b`` with ``b[i][j]`` the fraction of broken
        parent-``j`` mass landing in daughter bin ``i`` (upper-triangular,
        ``i <= j``); every column sums to 1.

    Raises:
        ValueError: for invalid appearance data.
    """
    n_int = len(d_char)
    b = [[0.0] * n_int for _ in range(n_int)]
    curve = None
    if name == "t10_appearance":
        curve = t10_appearance_curve(params["appearance_table"], params["t10"])
    for j in range(n_int):
        if j == 0:
            b[0][0] = 1.0  # All broken mass from the finest parent stays in bin 0.
            continue
        d_j = d_char[j]
        if name == "t10_appearance":
            # t10 markers use daughter upper edges over the parent characteristic
            # size. Differencing these cumulatives gives each bin's mass.
            b_cum = [curve(edges[i + 1] / d_j) for i in range(j)]
        elif name == "vogel":
            # Evaluate the curve at the finer bins' characteristic sizes.
            b_cum = [
                min(1.0, _cumulative(name, d_char[i], d_char[i], d_j, params))
                for i in range(j)
            ]
        else:
            # cumulative at daughter edges x_1..x_j (b_cum[k-1] = B_j(x_k))
            b_cum = [
                min(1.0, _cumulative(name, edges[k], edges[k + 1], d_j, params))
                for k in range(1, j + 1)
            ]
        b[0][j] = b_cum[0]  # Bin 0 also includes material below the lowest mesh edge.
        for i in range(1, j):
            b[i][j] = b_cum[i] - b_cum[i - 1]
        b[j][j] = 1.0 - b_cum[j - 1]
    return b


def ecs_from_table(ecs_table, t10, d_char_m):
    """Return specific comminution energy (kWh/t) for each size in ``d_char_m``.

    ``d_char_m`` contains positive characteristic sizes in m. ``ecs_table``
    has exactly three fields: ``sizes`` (m), ``t10_rows`` (%), and ``ecs``
    (kWh/t). Each energy row corresponds to one t10 level and has one value
    per table size. All fields must be nonempty lists or tuples. Sizes and
    t10 levels must be finite, positive, and increasing; energies must be
    finite and nonnegative, and every value a real number. Invalid tables
    raise ``ValueError``.

    For each table size, a natural cubic spline through the origin and the
    t10 rows gives the energy at ``t10``. Across table sizes, another natural
    cubic spline in log(size) gives the energy at each requested size. Both
    splines extrapolate linearly beyond their endpoints; negative results
    become zero. With one table size, its energy at ``t10`` applies to every
    requested size.

    See Napier-Munn et al. (1996), ch. 4, sec. 4.6.2, and ch. 11,
    sec. 11.2.4, for the ore-specific Ecs-t10-size relationship.
    """

    def evaluate(xs, ys, x):
        if len(xs) == 1:
            return ys[0]
        spline = CubicSpline(xs, ys, bc_type="natural")
        end = xs[0] if x < xs[0] else xs[-1] if x > xs[-1] else None
        if end is None:
            return float(spline(x))
        return float(spline(end) + spline(end, 1) * (x - end))

    keys = ["ecs", "sizes", "t10_rows"]
    if not isinstance(ecs_table, Mapping) or set(ecs_table) != set(keys):
        raise ValueError(f"ecs_table must be a mapping with exactly the keys {keys}.")
    if not all(isinstance(ecs_table[k], (list, tuple)) and ecs_table[k] for k in keys):
        raise ValueError(f"ecs_table entries {keys} must be nonempty lists or tuples.")
    if any(
        not isinstance(row, (list, tuple)) or len(row) != len(ecs_table["sizes"])
        for row in ecs_table["ecs"]
    ):
        raise ValueError(
            "each ecs_table 'ecs' row must be a list or tuple with one value per size."
        )
    if len(ecs_table["ecs"]) != len(ecs_table["t10_rows"]):
        raise ValueError(
            f"ecs_table needs one 'ecs' row per t10 row (got "
            f"{len(ecs_table['ecs'])} rows for {len(ecs_table['t10_rows'])} t10 "
            "values)."
        )
    sizes = [
        _finite(v, f"ecs_table['sizes'][{i}]") for i, v in enumerate(ecs_table["sizes"])
    ]
    levels = [
        _finite(v, f"ecs_table['t10_rows'][{i}]")
        for i, v in enumerate(ecs_table["t10_rows"])
    ]
    rows = [
        [_finite(v, f"ecs_table['ecs'][{i}][{j}]") for j, v in enumerate(row)]
        for i, row in enumerate(ecs_table["ecs"])
    ]
    for name, xs in (("sizes", sizes), ("t10_rows", levels)):
        if xs[0] <= 0.0 or any(b <= a for a, b in zip(xs, xs[1:])):
            raise ValueError(
                f"ecs_table '{name}' must be finite, positive and increasing."
            )
    if any(v < 0.0 for row in rows for v in row):
        raise ValueError("ecs_table 'ecs' entries must be finite and nonnegative.")
    at_t10 = [
        evaluate([0.0, *levels], [0.0, *(row[j] for row in rows)], t10)
        for j in range(len(sizes))
    ]
    logs = [math.log(s) for s in sizes]
    return tuple(max(0.0, evaluate(logs, at_t10, math.log(d))) for d in d_char_m)
