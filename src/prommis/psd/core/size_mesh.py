#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Size-mesh validation and characteristic-size helpers.

A PSD mesh is an ascending list of ``N+1`` strictly increasing size edges in SI
meters, defining ``N`` intervals.  The optional ``bottom_size`` substitutes for a
zero finest edge only when computing characteristic sizes; mass accounting still
uses the stated interval.
"""

__author__ = "Daison Yancy Caballero"

import math

from idaes.core.util.exceptions import ConfigurationError


def _finite_float(value, what):
    """Return ``value`` as a finite float, rejecting booleans and overflow."""
    if isinstance(value, bool):
        raise ValueError(f"{what} must be a number, not a boolean")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        # Keep public mesh helpers on their documented ValueError path.
        raise ValueError(f"{what} must be a finite-representable number") from exc
    if not math.isfinite(out):
        raise ValueError(f"{what} must be finite")
    return out


def validate_mesh(edges, bottom_size=None):
    """Validate an ascending PSD size-mesh edge list (SI meters).

    Rules: at least 2 edges (>= 1 interval), all finite, strictly increasing,
    and strictly positive after the optional bottom-size substitution.

    ``bottom_size`` is a positive replacement for a zero finest edge
    (``x_0 == 0``) used only for characteristic-size purposes (mass accounting
    still uses the stated interval).  It is legal *only* when ``x_0 == 0`` and
    must satisfy ``0 < bottom_size < x_1``; supplying it alongside ``x_0 > 0``
    raises.

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
        raise ValueError(
            "size mesh edges must be finite-representable numbers"
        ) from exc
    if len(edges) < 2:
        raise ValueError("size mesh must have at least 2 edges (>= 1 interval)")
    for e in edges:
        if not math.isfinite(e):
            raise ValueError("size mesh edges must all be finite")
    for lo, hi in zip(edges, edges[1:]):
        if not lo < hi:
            raise ValueError("size mesh edges must be strictly increasing")
    if bottom_size is not None:
        if isinstance(bottom_size, bool):
            raise ValueError("bottom_size must be a number, not a boolean")
        try:
            bottom_size = float(bottom_size)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(
                "bottom_size must be a finite-representable number"
            ) from exc
        if edges[0] != 0.0:
            raise ValueError(
                "bottom_size may only be supplied when the finest edge x_0 == 0"
            )
        if not math.isfinite(bottom_size) or bottom_size <= 0.0:
            raise ValueError("bottom_size must be a finite positive length")
        if not bottom_size < edges[1]:
            raise ValueError("bottom_size must be strictly less than x_1")
    elif edges[0] <= 0.0:
        raise ValueError(
            "finest edge x_0 must be > 0 (or supply bottom_size when x_0 == 0)"
        )
    return tuple(edges)


def effective_edges(edges, bottom_size=None):
    """Return validated edges with ``x_0`` replaced for characteristic sizing."""
    eff = list(validate_mesh(edges, bottom_size))
    if eff[0] == 0.0:
        eff[0] = float(bottom_size)
    return tuple(eff)


def characteristic_sizes(edges, bottom_size=None):
    """Geometric-mean characteristic size of each interval.

    ``d_char[k] = sqrt(x_k_eff * x_{k+1})`` with ``x_0_eff = bottom_size`` when
    the finest edge is zero.

    Returns:
        tuple of ``N`` characteristic sizes.
    """
    eff = effective_edges(edges, bottom_size)
    return tuple(math.sqrt(eff[k] * eff[k + 1]) for k in range(len(eff) - 1))


def meshes_equal(edges_a, edges_b, rtol=1e-9, atol=1e-12):
    """Return True if two edge lists have equal length and match within tolerance.

    Edges are coerced through the same finite-float gate used by mesh validation,
    so malformed mesh data raises ``ValueError`` instead of being compared.
    """
    a = [_finite_float(e, "mesh edge") for e in edges_a]
    b = [_finite_float(e, "mesh edge") for e in edges_b]
    if len(a) != len(b):
        return False
    # Use a symmetric relative reference so argument order cannot affect the
    # result on small SI edges.
    return all(abs(x - y) <= atol + rtol * max(abs(x), abs(y)) for x, y in zip(a, b))


def validate_property_mesh(config, missing_edges_message):
    """Validate property-package mesh configuration."""
    if config.size_edges is None:
        raise ConfigurationError(missing_edges_message)
    raw_edges = config.size_edges
    bottom = config.bottom_size
    if bottom is not None:
        # reject bool before coercion: float(True) == 1.0 would otherwise slip
        # past validate_mesh's bool guard, which only sees the coerced float
        if isinstance(bottom, bool):
            raise ConfigurationError("bottom_size must be a number, not a boolean.")
        try:
            bottom = float(bottom)
        except (TypeError, ValueError, OverflowError) as exc:
            # no repr(bottom): a huge int trips Python's int-to-str digit limit
            raise ConfigurationError("bottom_size must be a finite number.") from exc
    # validate_mesh stays on its framework-neutral ValueError path; re-raise as
    # ConfigurationError for invalid IDAES package config.
    try:
        edges = list(validate_mesh(raw_edges, bottom))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(str(exc)) from exc
    return edges, bottom, len(edges) - 1


def assert_same_mesh(parameter_block, other, parameter_block_kind):
    """Raise unless ``other`` carries an identical size mesh."""
    if not hasattr(other, "_edge_values"):
        raise ConfigurationError(
            f"assert_same_mesh requires another {parameter_block_kind} "
            f"(got {type(other).__name__})."
        )
    if not meshes_equal(parameter_block._edge_values, other._edge_values):
        raise ConfigurationError(
            "size-mesh mismatch between parameter blocks:\n"
            f"  this : {parameter_block._edge_values}\n"
            f"  other: {other._edge_values}"
        )
