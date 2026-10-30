#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Assumed-shape cumulative distribution forms and mesh discretization.

Supported shapes are Rosin-Rammler (``Q(d) = 1 - exp(-(d/d_63)^n)``) and
Gates-Gaudin-Schuhmann (``Q(d) = (min(d, d_max)/d_max)^n``), with a smooth
approximation to the GGS cap. Per-interval retained fractions come from
``core.psd_math.fractions_from_cdf``. Material below the lowest mesh edge
is included in bin 0, and the fractions sum to 1.

``shape_size_upper_bound`` returns the largest characteristic size allowed
by the mesh. Rosin-Rammler needs at least 95% passing at the top edge; GGS
needs ``d_max`` at or below that edge.

References:

[1] Rosin, P., & Rammler, E. (1933). The laws governing the fineness of
    powdered coal. J. Inst. Fuel, 7, 29-36.

[2] Schuhmann, R., Jr. (1940). Principles of comminution, I: size
    distribution and surface calculations. AIME Technical Publication
    No. 1189.
"""

__author__ = "Daison Yancy Caballero"

import math
from collections.abc import Callable
from typing import NamedTuple

from pyomo.environ import exp, value

from prommis.comminution.core.psd_math import P80_PASSING, fractions_from_cdf
from prommis.comminution.core.smooth_functions import EPS_SMOOTH

# Minimum Rosin-Rammler passing at the top mesh edge.
RR_VALIDITY_THRESHOLD = 0.95


def rosin_rammler_cdf(d, d63, n):
    """Rosin-Rammler cumulative passing ``1 - exp(-(d/d63)^n)``.

    See Rosin and Rammler (1933).
    """
    return 1 - exp(-((d / d63) ** n))


def ggs_cdf(d, d_max, n):
    """Return Gates-Gaudin-Schuhmann passing with a smooth upper limit.

    The ideal law for ``d >= 0`` is ``Q(d) = (min(d, d_max)/d_max)^n``
    (Schuhmann, 1940). For ``d > 0``, the formula uses ``d_max/d`` so
    interval checks can bound the smoothed ratio between zero and one even
    as ``d_max`` approaches zero. For ``d <= 0``, passing is zero.
    At ``d = d_max``, passing is slightly below one.

    Args:
        d: Size where passing is evaluated, in the same units as ``d_max``.
           It must have a value at build and must not switch between
           nonpositive and positive values afterward.
        d_max: Characteristic maximum size, assumed positive.
        n: Positive shape exponent.
    """
    if value(d) <= 0:
        return 0.0
    inverse_ratio = d_max / d
    limited_ratio = 2 / (
        1
        + inverse_ratio
        + ((1 - inverse_ratio) ** 2 + (EPS_SMOOTH * inverse_ratio) ** 2) ** 0.5
    )
    return limited_ratio**n


def _rr_upper_bound(x_top, n):
    """Return the top-edge limit for Rosin-Rammler ``d_63``, adjusted for rounding."""
    if not (math.isfinite(n) and n > 0.0):
        raise ValueError("shape exponent n must be a finite positive number")
    upper = x_top / ((-math.log(1.0 - RR_VALIDITY_THRESHOLD)) ** (1.0 / n))
    while (
        math.isfinite(upper)
        and upper > 0.0
        and rosin_rammler_cdf(x_top, upper, n) < RR_VALIDITY_THRESHOLD
    ):
        upper = math.nextafter(upper, 0.0)
    return upper


def _form(shape):
    """Return the ``ShapeForm`` of ``shape``."""
    try:
        return SHAPES[shape]
    except KeyError:
        raise ValueError(f"unknown distribution shape {shape!r}") from None


def p80_ratio(shape, n):
    """Return P80 divided by the shape's characteristic size for exponent ``n``."""
    return _form(shape).p80_ratio(n)


def shape_size_from_p80(shape, p80, n):
    """Return the characteristic size whose analytic P80 is ``p80``."""
    return p80 / p80_ratio(shape, n)


def shape_size_upper_bound(shape, x_top, n):
    """Return the mesh-based upper limit for a distribution's characteristic size.

    For Rosin-Rammler, the limit must meet ``RR_VALIDITY_THRESHOLD`` (0.95)
    at the top edge. If rounding puts passing just below that threshold,
    move the limit slightly lower.
    For GGS, return the top mesh edge.

    Raises:
        ValueError: If ``shape`` is unknown or the Rosin-Rammler exponent is
            not finite and positive.
    """
    return _form(shape).upper_bound(x_top, n)


def shape_cdf(shape, d, shape_size, n):
    """Return the cumulative passing of ``shape`` at size ``d``.

    ``d`` may be a number or Pyomo expression. For GGS, it must be known when
    the expression is built because the function chooses its zero-size
    calculation then.
    """
    return _form(shape).cdf(d, shape_size, n)


def distribution_fractions(shape, edges, shape_size, n):
    """Discretize an assumed-shape cumulative form onto a mesh.

    Args:
        shape: ``"rosin_rammler"`` or ``"gates_gaudin_schuhmann"``.
        edges: ``N+1`` ascending numeric size edges, in the characteristic
            size's length unit.
        shape_size: numeric characteristic size (``d_63`` for RR, ``d_max``
            for GGS), already checked against :func:`shape_size_upper_bound`.
        n: positive finite shape exponent.

    Returns:
        list of ``N`` retained fractions that sum to 1.
    """
    form = _form(shape)
    return fractions_from_cdf(lambda x: form.cdf(x, shape_size, n), edges)


class ShapeForm(NamedTuple):
    """Calculations and characteristic-size labels for one PSD shape."""

    shape_size_name: str  # characteristic-size name: d_63 or d_max
    shape_size_doc: str  # description used for the characteristic-size Var
    cdf: Callable  # (d, characteristic_size, n) -> cumulative passing
    p80_ratio: Callable  # (n) -> P80 / characteristic size
    upper_bound: Callable  # (x_top, n) -> mesh-based characteristic-size limit


SHAPES = {
    "rosin_rammler": ShapeForm(
        shape_size_name="d_63",
        shape_size_doc="Rosin-Rammler size at 63.2% passing",
        cdf=rosin_rammler_cdf,
        p80_ratio=lambda n: (-math.log(0.2)) ** (1.0 / n),
        upper_bound=_rr_upper_bound,
    ),
    "gates_gaudin_schuhmann": ShapeForm(
        shape_size_name="d_max",
        shape_size_doc="Gates-Gaudin-Schuhmann maximum size",
        cdf=ggs_cdf,
        p80_ratio=lambda n: P80_PASSING ** (1.0 / n),
        upper_bound=lambda x_top, n: float(x_top),
    ),
}

SHAPE_SIZE_NAME = {name: form.shape_size_name for name, form in SHAPES.items()}
