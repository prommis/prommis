#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Representation-independent particle-size-distribution (PSD) math.

Conventions (used everywhere in this package):

* **Ascending size order, retained-mass basis.**  A mesh of ``N+1`` strictly
  increasing edges ``x_0 < x_1 < ... < x_N`` defines ``N`` intervals; interval
  ``k`` holds particles in ``[x_k, x_{k+1})`` and ``k = 0`` is the finest bin.
* The PSD is stored as **retained mass flow per interval**; fractions and
  cumulative passing are derived.  Cumulative passing is evaluated at the
  interval upper edges: ``cum[k] = sum_{i<=k} w_i / sum_i w_i``.
* **Percentile (Pxx) interpolation is linear-in-size** on the *augmented*
  cumulative-passing curve, the point set
  ``[(x_0, 0), (x_1, cum[0]), ..., (x_N, 1)]``.  The ``(x_0, 0)`` anchor makes
  all-fine feeds and bin-0 percentiles well defined.
"""

__author__ = "Daison Yancy Caballero"

from .smooth_functions import EPS_SMOOTH, smooth_clamp

# Additive floor on cumulative-passing increments inside the smooth percentile
# form.  It avoids divide-by-zero on flat bins and bounds the local gradient for
# IPOPT when a PSD has many zero bins.
_EPS_DENOM = 1e-6


def size_at_passing_smooth(edges, cum, target, eps=EPS_SMOOTH):
    """Smooth, IPOPT-friendly linear-in-size percentile interpolation.

    Piecewise-linear interpolation on the augmented cumulative curve, written as
    a single differentiable expression via a smooth clamp of the per-interval
    fractional position.  ``edges`` and ``cum`` may be floats or Pyomo
    expressions; ``target`` and ``eps`` are floats.

    On well-formed inputs, this agrees with the exact piecewise-linear result to
    within the requested ``eps`` tolerance away from cumulative-curve nodes.
    Unlike an exact numeric form it does not reject flat cumulative intervals, so
    degenerate/near-zero states can be evaluated for diagnostics.
    """
    cum = list(cum)
    n = len(cum)
    aug = [0.0] + cum
    result = edges[0]
    for m in range(n):
        d_f = aug[m + 1] - aug[m]
        frac = (target - aug[m]) / (d_f + _EPS_DENOM)
        result = result + (edges[m + 1] - edges[m]) * smooth_clamp(frac, 0.0, 1.0, eps)
    return result
