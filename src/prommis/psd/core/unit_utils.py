#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Shared scaling helpers used by the bulk PSD property package.

This module declares no Pyomo modeling objects.
"""

__author__ = "Daison Yancy Caballero"

from pyomo.common.config import ConfigValue, In
from pyomo.environ import value

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum


class FactorSource(StrEnum):
    """Variable-factor sources selectable on every comminution scaler."""

    seeded = "seeded"
    baseline = "baseline"


def declare_factor_source(config):
    """Declare the ``factor_source`` option on ``config`` and return it.

    Every comminution scaler declares this on its own fresh
    ``CustomScalerBase.CONFIG()`` (never on the shared base ConfigDict), so
    ``call_submodel_scaler_method``'s ``scaler(**self.config)`` instantiation
    carries the selected source to delegated state-block scalers unchanged.
    """
    config.declare(
        "factor_source",
        ConfigValue(
            default=FactorSource.baseline,
            domain=In(FactorSource),
            description="Variable-factor source: 'baseline' (default) derives "
            "a-priori factors from declared factors, fixed feed-spec reads, "
            "and config estimates, and is safe on an uninitialized model; "
            "'seeded' reads post-initialization magnitudes (the pre-promotion "
            "default, retained for value-refined factors).",
        ),
    )
    return config


def fixed_value_or_none(var):
    """``value(var)`` if ``var`` is fixed, else ``None``.

    The single gate for baseline-mode feed reads: a fixed Var carries user
    data (safe before initialization), an unfixed Var carries a seeded or
    computed quantity a baseline path must never read.  A fixed Var without
    a value (pathological) also returns ``None`` rather than raising.

    The returned magnitude is in the variable's OWN units and must only be
    composed into factors for same-basis variables; cross-basis derivations
    use :func:`fixed_magnitude_or_none`.
    """
    return value(var, exception=False) if var.fixed else None


def required_scaling_factor(scaler, component):
    """Declared scaling factor of ``component``; raise if missing/non-positive.

    Constraint factors derive from declared variable factors, so a missing
    factor is a bug, not a fallback trigger -- run the variable scaling routine
    (or set factors manually) before scaling constraints.
    """
    sf = scaler.get_scaling_factor(component)
    if sf is not None and sf > 0:
        return sf
    raise ConfigurationError(f"Missing scaling factor for {component.name}.")
