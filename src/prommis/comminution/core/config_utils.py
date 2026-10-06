#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Validate configuration values for comminution models.
"""

__author__ = "Daison Yancy Caballero"

import math
import numbers

import numpy
from pyomo.core.base.units_container import UnitsError
from pyomo.environ import units, value

from idaes.core.util.exceptions import ConfigurationError


def _is_real_number(val):
    """Return whether ``val`` is a non-boolean ``numbers.Real`` scalar."""
    return isinstance(val, numbers.Real) and not isinstance(val, (bool, numpy.bool_))


def _cfg_float(val, what, lo=None, hi=None, positive=False):
    """
    Return a finite float from a real number, excluding booleans.

    Optional ``lo`` and ``hi`` bounds are inclusive; ``positive`` requires > 0.
    Invalid values raise ``ConfigurationError``.
    """
    if isinstance(val, (bool, numpy.bool_)):
        raise ConfigurationError(
            f"{what} must be a real number, not a boolean (got {val!r})."
        )
    if isinstance(val, (str, bytes)):
        raise ConfigurationError(
            f"{what} must be a number, not a string (got {val!r})."
        )
    if not _is_real_number(val):
        raise ConfigurationError(
            f"{what} must be a real int or float (got {type(val).__name__} "
            f"{val!r}); convert other numeric types before configuring."
        )
    try:
        out = float(val)
    except OverflowError as exc:
        raise ConfigurationError(f"{what} must be convertible to a float.") from exc
    if not math.isfinite(out):
        raise ConfigurationError(f"{what} must be finite (got {val!r}).")
    if positive and out <= 0.0:
        raise ConfigurationError(f"{what} must be > 0 (got {out}).")
    if lo is not None and out < lo:
        raise ConfigurationError(f"{what} must be >= {lo} (got {out}).")
    if hi is not None and out > hi:
        raise ConfigurationError(f"{what} must be <= {hi} (got {out}).")
    return out


def _cfg_quantity(val, what, to_units, lo=None, hi=None, positive=False):
    """Convert a unit-bearing config value to a finite float in ``to_units``.

    Require Pyomo units (e.g. ``25 * units.mm``) to avoid ambiguous bare numbers.
    Apply inclusive ``lo``/``hi`` bounds and strict positivity after conversion
    to ``to_units``. Invalid values raise ``ConfigurationError``.
    """
    if isinstance(val, (numbers.Real, numpy.bool_)):
        raise ConfigurationError(
            f"{what} must carry units: pass a Pyomo units expression "
            f"convertible to {to_units} (e.g. value * units.<unit>), not a "
            f"bare number (got {val!r})."
        )
    if isinstance(val, (str, bytes)):
        raise ConfigurationError(
            f"{what} must be a Pyomo units expression, not a string " f"(got {val!r})."
        )
    try:
        converted = value(units.convert(val, to_units=to_units))
    except UnitsError as exc:
        raise ConfigurationError(
            f"{what} must be convertible to {to_units} (got {val!r})."
        ) from exc
    except (TypeError, ValueError, AttributeError) as exc:
        raise ConfigurationError(
            f"{what} must be a Pyomo units expression convertible to "
            f"{to_units} (got {type(val).__name__} {val!r})."
        ) from exc
    return _cfg_float(converted, what, lo=lo, hi=hi, positive=positive)


def check_param_value(
    param,
    what,
    *,
    lo=None,
    hi=None,
    positive=False,
    finite_error="{what} must be finite (got {value!r}).",
):
    """Return a finite Pyomo Param value within the requested bounds.

    Optional ``lo`` and ``hi`` bounds are inclusive; ``positive`` requires > 0.
    ``finite_error`` accepts ``{what}`` and ``{value}`` format fields.
    """
    out = value(param)
    if not _is_real_number(out) or not math.isfinite(out):
        raise ConfigurationError(finite_error.format(what=what, value=out))
    if positive and out <= 0.0:
        raise ConfigurationError(f"{what} must be positive (got {out}).")
    if lo is not None and out < lo:
        raise ConfigurationError(f"{what} must be >= {lo} (got {out}).")
    if hi is not None and out > hi:
        raise ConfigurationError(f"{what} must be <= {hi} (got {out}).")
    return out
