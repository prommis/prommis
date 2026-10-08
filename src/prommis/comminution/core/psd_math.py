#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Particle-size distribution (PSD) calculations on ascending size meshes.

Conventions used by these helpers:

* **Ascending size order, retained-mass basis.**  A mesh of ``N+1`` strictly
  increasing edges ``x_0 < x_1 < ... < x_N`` defines ``N`` intervals; interval
  ``k`` holds particles in ``[x_k, x_{k+1})`` and ``k = 0`` is the finest bin.
* The PSD is stored as **retained mass flow per interval**; fractions and
  cumulative passing are derived.  Cumulative passing is evaluated at the
  interval upper edges: ``cum[k] = sum_{i<=k} w_i / sum_i w_i``.
* **Numeric percentile interpolation is linear in size** on the *augmented*
  cumulative-passing curve, the point set
  ``[(x_0, 0), (x_1, cum[0]), ..., (x_N, 1)]``.  The ``(x_0, 0)`` anchor makes
  all-fine feeds and bin-0 percentiles well defined. The smooth Pyomo form
  approximates this interpolation for model equations.

Numeric normalization assumes the sum of finite flows remains finite in float
arithmetic.
"""

__author__ = "Daison Yancy Caballero"

import math

from prommis.comminution.core.smooth_functions import EPS_SMOOTH, smooth_clamp

# Passing fraction of the P80 percentile.
P80_PASSING = 0.8

# Additive floor on cumulative-passing increments inside the smooth percentile
# form.  It avoids divide-by-zero on flat bins and bounds the local gradient for
# IPOPT when a PSD has many zero bins.
_EPS_DENOM = 1e-6


def _as_floats(values):
    return [float(v) for v in values]


def cumulative_passing(flows):
    """Cumulative passing fraction at each interval upper edge (ascending bins).

    Args:
        flows: iterable of retained mass flows per interval (finest first).

    Returns:
        list of cumulative passing fractions ``cum[k] = sum_{i<=k} w_i / total``;
        ``cum[-1]`` is exactly ``1.0`` for a finite total.

    Raises:
        ValueError: if a flow is negative or non-finite, or the total is zero
            or non-finite.
    """
    f = _as_floats(flows)
    if any(not math.isfinite(fi) or fi < 0.0 for fi in f):
        raise ValueError("cumulative_passing requires finite, non-negative flows")
    total = sum(f)
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("cumulative_passing requires a strictly positive total flow")
    cum = []
    running = 0.0
    for fi in f:
        running += fi
        cum.append(running / total)
    cum[-1] = 1.0
    return cum


def size_at_passing(edges, cum, target):
    """Interpolate a numeric passing percentile linearly in size.

    Use the first augmented-curve interval with ``F_lo < target <= F_hi``.
    If ``target == F_hi``, return the first upper edge with that passing fraction.

    Args:
        edges: ``N+1`` ascending size edges.
        cum: ``N`` cumulative passing fractions (``cum[-1]`` should equal 1).
        target: passing fraction in ``(0, 1]``.

    Returns:
        the interpolated size at the requested passing fraction.

    Raises:
        ValueError: on an inconsistent mesh length (``len(edges) != len(cum)+1``)
            or fewer than one interval, an out-of-domain ``target`` (``<= 0`` or
            ``> 1``), or if no interval brackets ``target``.
    """
    edges = _as_floats(edges)
    cum = _as_floats(cum)
    n = len(cum)
    if len(edges) != n + 1:
        raise ValueError("edges must have exactly len(cum)+1 entries")
    if n < 1:
        raise ValueError("need at least one interval")
    if target <= 0.0 or target > 1.0:
        raise ValueError(f"passing target must lie in (0, 1]; got {target}")
    aug = [0.0] + cum
    for m in range(n):
        f_lo = aug[m]
        f_hi = aug[m + 1]
        if f_lo < target <= f_hi:
            if target == f_hi:
                return edges[m + 1]
            return edges[m] + (target - f_lo) / (f_hi - f_lo) * (
                edges[m + 1] - edges[m]
            )
    raise ValueError(
        "no bracketing interval found; the cumulative curve must reach the "
        "target (check that cum is normalized to 1 at the top edge)"
    )


def finest_attainable_size(edges, target=P80_PASSING):
    """Return the minimum attainable size at passing fraction ``target``.

    For :func:`size_at_passing`, this occurs when all mass is in the finest
    interval. The default target gives the minimum P80. Requires at least
    two ascending edges and ``0 < target <= 1``; inputs are not validated.
    """
    return edges[0] + target * (edges[1] - edges[0])


def size_at_passing_smooth(edges, cum, target, eps=EPS_SMOOTH):
    """Smooth, IPOPT-friendly linear-in-size percentile interpolation.

    Approximate :func:`size_at_passing` with a smooth clamp of each interval's
    fractional position. ``edges`` and ``cum`` may be floats or Pyomo
    expressions; ``target`` and ``eps`` are floats.

    Add ``_EPS_DENOM = 1e-6`` to each cumulative increment to avoid division by
    zero on flat bins. Both this floor and clamp smoothing affect accuracy;
    ``eps`` is a smoothing parameter, not an interpolation-error bound.
    This form does not validate inputs and can evaluate flat cumulative curves
    for diagnostics. The numeric form skips flat intervals when bracketing.
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


def fractions_from_passing(passing):
    """Return retained fractions from cumulative passing at the mesh edges.

    ``passing`` holds ``Q(x_0), ..., Q(x_N)``. The finest fraction is
    ``Q(x_1)/Q(x_N)``, including material below ``x_0``. For ``k >= 1``,
    fraction ``k`` is ``(Q(x_{k+1}) - Q(x_k))/Q(x_N)``. The fractions sum to 1.
    Edge order and monotonicity are assumed, not checked.

    A Python ``int`` or ``float`` top-edge value must be finite and positive.
    For Pyomo expressions, the caller must ensure positive top-edge passing.

    Returns:
        list with one retained fraction per interval.
    """
    q = list(passing)
    n = len(q) - 1
    if n < 1:
        raise ValueError("need at least one interval")
    q_top = q[-1]
    if isinstance(q_top, (int, float)) and (not math.isfinite(q_top) or q_top <= 0.0):
        raise ValueError(
            "cumulative distribution is non-positive or non-finite at the top edge"
        )
    fracs = [q[1] / q_top]
    for k in range(1, n):
        fracs.append((q[k + 1] - q[k]) / q_top)
    return fracs


def fractions_from_cdf(cdf, edges):
    """Evaluate ``cdf`` at each edge to return one retained fraction per interval.

    Normalize by passing at the top edge. The finest fraction includes
    material below the bottom edge, as in :func:`fractions_from_passing`.

    Args:
        cdf: Callable for cumulative passing ``Q(d)`` in ``[0, 1]``; may
            return numbers or Pyomo expressions.
        edges: At least two ascending size edges. Order is assumed, not checked.
    """
    edges = _as_floats(edges)
    if len(edges) < 2:
        raise ValueError("need at least one interval")
    return fractions_from_passing(cdf(x) for x in edges)
