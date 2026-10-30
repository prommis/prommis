#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Size-dependent selection and classification functions.

``selection_value`` dispatches the numeric crusher kernels.

References:

[1] Whiten, W. J. (1972). The simulation of crushing plants with models
    developed using multiple spline regression. Proc. 10th APCOM Symp.,
    Johannesburg, SAIMM, 317-323.

[2] Austin, L. G., Klimpel, R. R., & Luckie, P. T. (1984). Process
    Engineering of Size Reduction: Ball Milling. SME-AIME, New York.

[3] Napier-Munn, T. J., Morrell, S., Morrison, R. D., & Kojovic, T. (1996).
    Mineral Comminution Circuits: Their Operation and Optimisation. JKMRC,
    Univ. of Queensland.

[4] King, R. P. (2001). Modeling and Simulation of Mineral Processing
    Systems. Butterworth-Heinemann, Oxford.

[5] Vogel, L., & Peukert, W. (2003). Breakage behaviour of different
    materials - construction of a mastercurve for the breakage probability.
    Powder Technology, 129(1-3), 101-110.
"""

__author__ = "Daison Yancy Caballero"

import math

WHITEN_K3_DEFAULT = 2.3


def whiten_selection(d, k1, k2, k3):
    """Whiten classification ``C(d; K1, K2, K3)`` (Whiten, 1972).

    ``C = 0`` for ``d <= K1``, ``C = 1`` for ``d >= K2``, and
    ``C = 1 - ((K2 - d)/(K2 - K1))^K3`` between them; ``d``, ``K1`` and
    ``K2`` share one length unit. ``K3`` shapes the ramp; ``selection_value``
    defaults it to ``WHITEN_K3_DEFAULT`` (2.3; Napier-Munn et al., 1996,
    ch. 6, eq. 6.4).
    """
    if d <= k1:
        return 0.0
    if d >= k2:
        return 1.0
    return 1.0 - ((k2 - d) / (k2 - k1)) ** k3


def austin_selection(d, s1, d1, a):
    """Austin selection ``S = min(1, S1*(d/d1)^a)`` (Austin et al., 1984).

    ``d`` and ``d1`` share a length unit. The crusher requires ``S1`` and
    ``a`` in [0, 1] and ``d1`` > 0. Austin et al. describe selection as a
    breakage rate; here the capped value is a probability per stress event.
    """
    return min(1.0, s1 * (d / d1) ** a)


def vogel_peukert_selection(d, f_mat, xw_min, v):
    """Vogel-Peukert selection from impact energy (Vogel & Peukert, 2003).

    ``S = 1 - exp(-max(d*(1/2 v^2) - xw_min, 0)*f_mat)``, with ``d`` in m,
    ``v`` in m/s, ``xw_min`` in J*m/kg and ``f_mat`` in kg/(J*m). The clip
    at zero gives ``S = 0`` at and below the threshold. This is the
    probability for one stress event; the unit applies it once per event.
    """
    w_kin = 0.5 * v**2
    excess = d * w_kin - xw_min
    return 1.0 - math.exp(-max(0.0, excess) * f_mat)


def selection_value(name, d, params):
    """Return the unsmoothed numeric selection value for a named function.

    ``"king"`` (King, 2001) is the :func:`whiten_selection` ramp between
    ``CSS*alpha1`` and ``CSS*alpha2`` with exponent ``n``; ``d`` and ``CSS``
    share one length unit.
    """
    if name == "whiten":
        return whiten_selection(
            d, params["K1"], params["K2"], params.get("K3", WHITEN_K3_DEFAULT)
        )
    if name == "king":
        css = params["CSS"]
        return whiten_selection(
            d, css * params["alpha1"], css * params["alpha2"], params["n"]
        )
    if name == "austin":
        return austin_selection(d, params["S1"], params["d1"], params["a"])
    if name == "vogel_peukert":
        return vogel_peukert_selection(
            d, params["f_mat"], params["xw_min"], params["v"]
        )
    raise ValueError(f"unknown selection function {name!r}")
