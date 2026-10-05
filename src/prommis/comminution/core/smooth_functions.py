#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Smooth clamp helpers for particle size distribution (PSD) models.

The helpers use IDAES ``smooth_max`` and ``smooth_min`` to smooth piecewise
expressions and accept Python floats or Pyomo expressions.
"""

__author__ = "Daison Yancy Caballero"

from idaes.core.util.math import smooth_max, smooth_min

# Default smoothing parameter: larger values round transitions more but
# increase deviation from the unsmoothed function.
EPS_SMOOTH = 1e-4


def smooth_clamp(value, lower, upper, eps=EPS_SMOOTH):
    """Smooth approximation of ``min(upper, max(lower, value))``.

    Args:
        value: quantity to clamp (float or Pyomo expression).
        lower: lower bound (float or Pyomo expression).
        upper: upper bound (float or Pyomo expression).
        eps: smoothing parameter (default :data:`EPS_SMOOTH`).

    Returns:
        A smooth approximation to clamping; finite ``eps`` does not guarantee
        exact bounds or endpoint values.
    """
    return smooth_min(upper, smooth_max(lower, value, eps), eps)


def smooth_max0(value, eps=EPS_SMOOTH):
    """Approximate ``max(0, value)`` smoothly.

    For ``eps > 0``, the result at ``value == 0`` is ``eps / 2``.
    """
    return smooth_max(0.0, value, eps)
