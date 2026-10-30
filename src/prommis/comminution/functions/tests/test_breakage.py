#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for breakage functions and matrix builders."""

import math

import pytest

from prommis.comminution.core.size_mesh import characteristic_sizes
from prommis.comminution.functions.breakage import (
    breakage_matrix,
    ecs_from_table,
    reid_stewart_cumulative,
    t10_appearance_curve,
)

EDGES4 = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]  # mm-scale mesh in meters, 4 intervals

# Synthetic appearance table: two rows so a midpoint t10 interpolates linearly
# (2-row PCHIP = secant line).  tn = percent passing 1/n of the parent, so
# markers ascend t75 < t50 < t25 < t10 < t4 < t2.
_APPEARANCE_TABLE = {
    "10": {"t75": 4, "t50": 6, "t25": 8, "t4": 22, "t2": 40},
    "30": {"t75": 12, "t50": 18, "t25": 26, "t4": 50, "t2": 70},
}


def _column(b, j):
    return [b[i][j] for i in range(len(b))]


@pytest.mark.unit
def test_luckie_austin_matrix_handcalc():
    # Phi = 0.5, gamma = 1, beta = 3, parent = top bin (j = 3)
    params = {"phi": 0.5, "gamma": 1.0, "beta": 3.0}
    d_char = characteristic_sizes(EDGES4)
    d3 = d_char[3]
    b = breakage_matrix("luckie_austin", EDGES4, d_char, params)

    def b_cum(k):  # B_3(x_k) = min(1, 0.5 (x_k/d3) + 0.5 (x_{k+1}/d3)^3)
        return min(1.0, 0.5 * (EDGES4[k] / d3) + 0.5 * (EDGES4[k + 1] / d3) ** 3)

    assert b[0][3] == pytest.approx(b_cum(1), rel=1e-12)
    assert b[1][3] == pytest.approx(b_cum(2) - b_cum(1), rel=1e-12)
    assert b[2][3] == pytest.approx(b_cum(3) - b_cum(2), rel=1e-12)
    assert b[3][3] == pytest.approx(1.0 - b_cum(3), rel=1e-12)
    for j in range(len(d_char)):
        assert sum(_column(b, j)) == pytest.approx(1.0, rel=1e-12)
    with pytest.raises(ValueError, match="unknown breakage function 'bogus'"):
        breakage_matrix("bogus", EDGES4, d_char, {})


@pytest.mark.unit
def test_breakage_variant_discrimination():
    # Austin uses the next-coarser edge in term two; Reid-Stewart uses the same edge.
    edges = [1e-3, 2e-3, 4e-3]
    params = {"phi": 0.9, "gamma": 1.0, "beta": 3.0}
    d_char = characteristic_sizes(edges)
    la = breakage_matrix("luckie_austin", edges, d_char, params)
    rs = breakage_matrix("reid_stewart", edges, d_char, params)

    # d1 = 2*sqrt(2) mm; Austin 0.9/sqrt(2) + 0.1 (sqrt(2))^3 = 1.3/sqrt(2),
    # Reid-Stewart 0.9/sqrt(2) + 0.1 (1/sqrt(2))^3 = 0.95/sqrt(2); both below 1
    la_expect = 1.3 / math.sqrt(2.0)
    rs_expect = 0.95 / math.sqrt(2.0)
    assert la[0][1] == pytest.approx(la_expect, rel=1e-12)
    assert rs[0][1] == pytest.approx(rs_expect, rel=1e-12)
    for mat in (la, rs):
        for j in range(len(d_char)):
            assert sum(_column(mat, j)) == pytest.approx(1.0, rel=1e-12)
    # Unequal weights reveal if gamma and beta are swapped; B = 0.2375.
    assert reid_stewart_cumulative(1.0, 2.0, 0.4, 1.0, 4.0) == pytest.approx(
        0.2375, rel=1e-12
    )


@pytest.mark.unit
def test_t10_appearance_curve_handcalc():
    # Two-row PCHIP is linear, so t10 = 20 gives each marker's midpoint.
    # The size-ratio interpolation passes through those markers.
    curve = t10_appearance_curve(_APPEARANCE_TABLE, 20)
    assert curve(1 / 75) == pytest.approx(0.08, rel=1e-12)  # (4 + 12) / 2 / 100
    assert curve(1 / 50) == pytest.approx(0.12, rel=1e-12)  # (6 + 18) / 2 / 100
    assert curve(1 / 25) == pytest.approx(0.17, rel=1e-12)  # (8 + 26) / 2 / 100
    assert curve(1 / 10) == pytest.approx(0.20, rel=1e-12)  # t10 anchor, 20 / 100
    assert curve(1 / 4) == pytest.approx(0.36, rel=1e-12)  # (22 + 50) / 2 / 100
    assert curve(1 / 2) == pytest.approx(0.55, rel=1e-12)  # (40 + 70) / 2 / 100
    # anchored at (0, 0) and (1, 1); ratios >= 1 clamp to 1
    assert curve(0.0) == 0.0
    assert curve(1.0) == pytest.approx(1.0, rel=1e-12)
    assert curve(2.0) == pytest.approx(1.0, rel=1e-12)
    swept = [curve(r / 100) for r in range(0, 101)]
    assert all(swept[i + 1] >= swept[i] - 1e-12 for i in range(len(swept) - 1))
    # integer-like row keys (int, float, padded string) select the same row
    for key in (10, 10.0, " 10 "):
        table = {key: _APPEARANCE_TABLE["10"], "30": _APPEARANCE_TABLE["30"]}
        assert t10_appearance_curve(table, 20)(0.25) == curve(0.25)


@pytest.mark.unit
def test_ecs_from_table_handcalc():
    cases = [
        # With one table size, the natural spline gives 0.3575 at t10=15,
        # rather than the straight-line result of 0.35.
        ([0.01], [[0.2], [0.5], [0.6]], 15.0, [0.01], (0.3575,)),
        # At t10=20, log-size interpolation gives 0.3 at 0.02 m;
        # extrapolation to 1 m is negative and is floored at zero.
        (
            [0.01, 0.04],
            [[0.2, 0.1], [0.4, 0.2], [0.6, 0.3]],
            20.0,
            [0.02, 1.0],
            (0.3, 0.0),
        ),
    ]
    for sizes, ecs, t10, d_char, expected in cases:
        table = {"sizes": sizes, "t10_rows": [10, 20, 30], "ecs": ecs}
        assert ecs_from_table(table, t10, d_char) == pytest.approx(
            expected, rel=1e-12, abs=1e-12
        )
        with pytest.raises(ValueError, match="one 'ecs' row per t10 row"):
            ecs_from_table(dict(table, ecs=ecs[:2]), t10, d_char)


NONUNIFORM = [1e-3, 2e-3, 3e-3, 9e-3, 12e-3, 30e-3]  # Size-coordinate tests.


# Standard AG/SAG appearance markers, Napier-Munn et al. (1996), ch. 4,
# Table 4.7, p. 83.
_T10_TABLE = {
    10: {"t75": 2.33, "t50": 3.06, "t25": 4.98, "t4": 23.33, "t2": 50.53},
    30: {"t75": 6.89, "t50": 9.41, "t25": 15.62, "t4": 61.58, "t2": 92.49},
    50: {"t75": 10.32, "t50": 14.71, "t25": 25.88, "t4": 82.86, "t2": 96.47},
}
PARAMS = {"appearance_table": _T10_TABLE, "t10": 30.0}


def _assemble(curve, edges, d_char, dsize, psize):
    """Assemble columns using the supplied daughter and parent size coordinates."""
    n = len(d_char)
    b = [[0.0] * n for _ in range(n)]
    for j in range(n):
        if j == 0:
            b[0][0] = 1.0
            continue
        b_cum = [min(1.0, curve(dsize[i] / psize[j])) for i in range(j)]
        b[0][j] = b_cum[0]
        for i in range(1, j):
            b[i][j] = b_cum[i] - b_cum[i - 1]
        b[j][j] = 1.0 - b_cum[j - 1]
    return b


def _max_delta(a, b):
    n = len(a)
    return max(abs(a[i][j] - b[i][j]) for i in range(n) for j in range(n))


@pytest.mark.unit
def test_t10_basis_is_edge_over_geomean_on_nonuniform_mesh():
    edges = NONUNIFORM
    d_char = list(characteristic_sizes(edges))
    curve = t10_appearance_curve(_T10_TABLE, 30.0)
    b = breakage_matrix("t10_appearance", edges, d_char, PARAMS)

    uppers = edges[1:]
    ug = _assemble(curve, edges, d_char, uppers, d_char)  # edge / geomean
    gg = _assemble(curve, edges, d_char, d_char, d_char)  # geomean / geomean
    uu = _assemble(curve, edges, d_char, uppers, uppers)  # edge / edge

    assert _max_delta(b, ug) < 1e-14
    # The other size choices differ by much more than the agreement tolerance.
    assert _max_delta(b, gg) > 1e-3
    assert _max_delta(b, uu) > 1e-3
