#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Validate PSD size meshes and calculate characteristic sizes.

A mesh has ``N+1`` strictly increasing edges defining ``N`` intervals. These
helpers use the caller's consistent length units; the PSD property package
supplies meters. If the finest edge is zero, ``bottom_size`` replaces it only
when computing characteristic sizes. Mass accounting still uses the stated
interval.
"""

__author__ = "Daison Yancy Caballero"

import math


def _finite_float(value, what):
    """Return ``value`` as a finite float, rejecting booleans and overflow."""
    if isinstance(value, bool):
        raise ValueError(f"{what} must be a number, not a boolean")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{what} must be convertible to a float") from exc
    if not math.isfinite(out):
        raise ValueError(f"{what} must be finite")
    return out


def validate_mesh(edges, bottom_size=None):
    """Validate ascending PSD size-mesh edges.

    Edges must use one consistent length unit; the PSD property package uses
    meters.

    Rules: at least 2 edges (>= 1 interval), all finite, strictly increasing,
    and strictly positive after the optional bottom-size substitution.

    ``bottom_size`` is a positive replacement for a zero finest edge
    (``x_0 == 0``) used only for characteristic-size purposes (mass accounting
    still uses the stated interval).  It is allowed only when ``x_0 == 0`` and
    must satisfy ``0 < bottom_size < x_1``; otherwise, ``ValueError`` is raised.
    Supplying it alongside ``x_0 > 0`` also raises ``ValueError``.

    Args:
        edges: iterable of size edges.
        bottom_size: optional positive replacement for a zero finest edge.

    Returns:
        the validated edges as a tuple of floats.

    Raises:
        ValueError: if any rule is violated.
    """
    try:
        raw = list(edges)
    except TypeError as exc:
        raise ValueError("size mesh edges must be a sequence of numbers") from exc
    # ``bool`` is an ``int`` subclass; reject it before float coercion.
    if any(isinstance(e, bool) for e in raw):
        raise ValueError("size mesh edges must be numbers, not booleans")
    try:
        edges = [float(e) for e in raw]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("size mesh edges must be convertible to floats") from exc
    if len(edges) < 2:
        raise ValueError("size mesh must have at least 2 edges (>= 1 interval)")
    for e in edges:
        if not math.isfinite(e):
            raise ValueError("size mesh edges must all be finite")
    for lo, hi in zip(edges, edges[1:]):
        if not lo < hi:
            raise ValueError("size mesh edges must be strictly increasing")
    if bottom_size is not None:
        bottom_size = _finite_float(bottom_size, "bottom_size")
        if edges[0] != 0.0:
            raise ValueError(
                "bottom_size may only be supplied when the finest edge x_0 == 0"
            )
        if bottom_size <= 0.0:
            raise ValueError("bottom_size must be a finite positive length")
        if not bottom_size < edges[1]:
            raise ValueError("bottom_size must be strictly less than x_1")
    elif edges[0] <= 0.0:
        raise ValueError(
            "finest edge x_0 must be > 0 (or supply bottom_size when x_0 == 0)"
        )
    return tuple(edges)


def characteristic_sizes(edges, bottom_size=None):
    """Return the geometric-mean characteristic size of each interval.

    ``d_char[k] = sqrt(x_k_eff * x_{k+1})`` with ``x_0_eff = bottom_size`` when
    the finest edge is zero.

    Returns:
        tuple of ``N`` characteristic sizes.
    """
    eff = list(validate_mesh(edges, bottom_size))
    if eff[0] == 0.0:
        eff[0] = float(bottom_size)
    return tuple(math.sqrt(eff[k] * eff[k + 1]) for k in range(len(eff) - 1))


def geometric_series(top, bottom, ratio):
    """Return ascending edges for a geometric sieve series.

    The number of intervals is the nearest integer to
    ``ln(top/bottom)/ln(ratio)``, with a minimum of one. Re-derive the ratio so
    both endpoints are hit exactly and the ratio between consecutive edges is
    constant.

    Args:
        top: largest edge (> bottom), in the same length unit as ``bottom``.
        bottom: smallest edge (> 0).
        ratio: target geometric ratio between consecutive edges (> 1).

    Returns:
        tuple of edges from ``bottom`` to ``top`` (inclusive), strictly
        increasing.

    Raises:
        ValueError: if ``top``, ``bottom``, or ``ratio`` is a boolean, non-finite,
            or an int too large to represent as a float; if ``0 < bottom < top``
            is not satisfied; or if ``ratio <= 1``.
    """
    top = _finite_float(top, "geometric_series top")
    bottom = _finite_float(bottom, "geometric_series bottom")
    ratio = _finite_float(ratio, "geometric_series ratio")
    if not (bottom > 0.0 and top > bottom):
        raise ValueError("require 0 < bottom < top")
    if not ratio > 1.0:
        raise ValueError("ratio must be > 1")
    n_intervals = max(int(round(math.log(top / bottom) / math.log(ratio))), 1)
    exact_ratio = (top / bottom) ** (1.0 / n_intervals)
    edges = [bottom * exact_ratio**i for i in range(n_intervals + 1)]
    edges[0] = bottom
    edges[-1] = top
    return tuple(edges)


def meshes_equal(edges_a, edges_b, rtol=1e-9, atol=1e-12):
    """Return True if two edge lists have equal length and match within tolerance.

    Reject non-finite or nonnumeric edges and booleans with ``ValueError``.
    This comparison does not validate minimum mesh length, edge ordering, or
    positivity. The absolute tolerance uses the same units as the edges. The
    test is symmetric: ``abs(a-b) <= atol + rtol*max(abs(a), abs(b))``.

    Args:
        edges_a: iterable of edges in one consistent length unit.
        edges_b: iterable of edges in the same length unit as ``edges_a``.
        rtol: dimensionless relative tolerance (default: 1e-9).
        atol: absolute tolerance in the edges' length unit (default: 1e-12).

    Returns:
        True if the lists have equal length and every edge pair satisfies the
        tolerance; False otherwise.

    Raises:
        ValueError: if an edge is a boolean, non-finite, or cannot be converted
            to a finite float.
    """
    a = [_finite_float(e, "mesh edge") for e in edges_a]
    b = [_finite_float(e, "mesh edge") for e in edges_b]
    if len(a) != len(b):
        return False
    return all(abs(x - y) <= atol + rtol * max(abs(x), abs(y)) for x, y in zip(a, b))
