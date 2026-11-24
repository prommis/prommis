#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Smooth clamp primitives for differentiable PSD math.

These helpers wrap IDAES ``smooth_max``/``smooth_min`` so piecewise PSD
quantities can be written as smooth IPOPT-friendly expressions.  They declare no
Pyomo modeling objects and operate on either Python floats or Pyomo expressions.
"""

__author__ = "Daison Yancy Caballero"

from idaes.core.util.math import smooth_max, smooth_min

# Default smoothing parameter for clamp/bracket primitives: small enough for
# close agreement away from kinks while keeping IPOPT Jacobians well behaved.
EPS_SMOOTH = 1e-4


def smooth_clamp(value, lower, upper, eps=EPS_SMOOTH):
    """Smooth approximation of ``min(upper, max(lower, value))``.

    Args:
        value: quantity to clamp (float or Pyomo expression).
        lower: lower bound (float or Pyomo expression).
        upper: upper bound (float or Pyomo expression).
        eps: smoothing parameter (default :data:`EPS_SMOOTH`).

    Returns:
        A smooth expression for ``value`` clamped to ``[lower, upper]``.
    """
    return smooth_min(upper, smooth_max(lower, value, eps), eps)
