#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Shared scaling and unit-conversion helpers for comminution models."""

__author__ = "Daison Yancy Caballero"

from pyomo.common.collections import ComponentMap
from pyomo.common.config import ConfigValue, In
from pyomo.environ import units, value

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum

# Mass- and volume-flow offsets used to avoid zero division.
MASS_FLOW_EPS = 1e-12 * units.kg / units.s
VOLUME_FLOW_EPS = 1e-15 * units.m**3 / units.s


class FactorSource(StrEnum):
    """Variable scaling sources used by comminution scalers."""

    current_values = "current_values"
    input_based = "input_based"


#: State-block Var families whose scaling factors propagate from inlet to outlet.
PROPAGATED_STATE_VAR_NAMES = (
    "flow_mass_sized_comp_size",
    "flow_mass_unsized_comp",
    "flow_mass_liquid_comp",
    "flow_mass_vapor_comp",
    "temperature",
    "pressure",
)


def declare_factor_source(config):
    """Declare ``factor_source`` on a scaler's own config and return it.

    Use a fresh ``CustomScalerBase.CONFIG()`` to avoid changing the shared base
    config.
    """
    config.declare(
        "factor_source",
        ConfigValue(
            default=FactorSource.input_based,
            domain=In(FactorSource),
            description="Variable scaling source: 'input_based' (default) "
            "derives factors from declared factors, fixed feed specifications, "
            "and configuration estimates, and never reads unfixed variable "
            "values, so it works before initialization; 'current_values' "
            "derives factors from current variable values, typically after "
            "initialization.",
        ),
    )
    return config


def fixed_value_or_none(var):
    """Return a fixed Var's value, or ``None`` if it is unfixed or unset.

    The numeric value is in the Var's declared units. Use
    :func:`fixed_magnitude_or_none` when deriving a factor in another unit basis.
    """
    return value(var, exception=False) if var.fixed else None


def unit_conversion_factor(var, target_units):
    """Multiplicative factor from ``var``'s declared units to ``target_units``.

    Pass any member of an indexed Var family; all members share the declaration.
    Dimensionally incompatible units raise ``InconsistentUnitsError``.
    """
    return value(units.convert(1.0 * units.get_units(var), to_units=target_units))


def fixed_magnitude_or_none(var, target_units):
    """Return a fixed Var's value converted to ``target_units``.

    Return ``None`` when the Var is unfixed or unset.
    """
    raw = value(var, exception=False) if var.fixed else None
    if raw is None:
        return None
    return raw * unit_conversion_factor(var, target_units)


def propagate_scaling_factors(scaler, src_block, dst_block, var_names, overwrite=False):
    """Copy factors for matching Vars from ``src_block`` to ``dst_block``.

    Skip names absent from the source. Require a positive factor for each present
    source Var. Existing destination factors are preserved unless ``overwrite``
    is true.
    """
    for name in var_names:
        src_var = getattr(src_block, name, None)
        if src_var is None:
            continue
        dst_var = getattr(dst_block, name)
        if src_var.is_indexed():
            for idx in src_var:
                scaler.set_variable_scaling_factor(
                    dst_var[idx],
                    required_scaling_factor(scaler, src_var[idx]),
                    overwrite=overwrite,
                )
        else:
            scaler.set_variable_scaling_factor(
                dst_var,
                required_scaling_factor(scaler, src_var),
                overwrite=overwrite,
            )


def delegate_state_block_scaling(
    scaler, submodels, default_scaler, submodel_scalers, method, overwrite=False
):
    """Call ``method`` on each state block's scaler.

    Use ``default_scaler`` unless ``submodel_scalers`` supplies an override.
    Instantiate scaler classes with compatible parent config options and pass
    ``overwrite`` through.
    """
    scalers = ComponentMap()
    for sub in submodels:
        scalers[sub] = default_scaler
    if submodel_scalers is not None:
        scalers.update(submodel_scalers)
    # Only pass options a scaler class declares; IDAES otherwise passes them all.
    for sub, sc in list(scalers.items()):
        if isinstance(sc, type):
            scalers[sub] = sc(
                **{k: scaler.config[k] for k in scaler.config if k in sc.CONFIG}
            )
    for sub in submodels:
        scaler.call_submodel_scaler_method(
            submodel=sub,
            submodel_scalers=scalers,
            method=method,
            overwrite=overwrite,
        )


def scale_constraints_at_time(scaler, model, name_factor_pairs, t, overwrite=False):
    """Apply per-family scaling factors to the time-``t`` members of each family.

    ``name_factor_pairs`` is an iterable of ``(constraint_name, factor)``.  A
    family absent from ``model`` is skipped; within a present family only the
    constraint data whose time index equals ``t`` are scaled (scalar- or
    tuple-indexed keys both handled).
    """
    for name, factor in name_factor_pairs:
        con = getattr(model, name, None)
        if con is None:
            continue
        for idx in con:
            if (idx[0] if isinstance(idx, tuple) else idx) != t:
                continue
            scaler.set_constraint_scaling_factor(con[idx], factor, overwrite=overwrite)


def required_scaling_factor(scaler, component):
    """Declared scaling factor of ``component``; raise if missing/non-positive.

    Constraint scaling requires declared variable factors. Run the variable
    scaling routine or set factors manually before scaling constraints.
    """
    sf = scaler.get_scaling_factor(component)
    if sf is not None and sf > 0:
        return sf
    raise ConfigurationError(
        f"Scaling factor for {component.name} is missing or non-positive."
    )


def inverse_sum_of_nominals(scaler, components, floor=1.0e-8):
    """Return ``1 / max(sum(1 / sf_i), floor)`` for component factors.

    Require a positive scaling factor for each component; do not read current
    variable values. ``floor`` is in the units of the summed nominals and keeps
    the factor finite when they are tiny.
    """
    total = sum(1.0 / required_scaling_factor(scaler, c) for c in components)
    return 1.0 / max(total, floor)


def inverse_magnitude_factor(reference_value, floor):
    """Return ``1 / max(abs(reference_value), floor)``.

    Treat ``None`` and zero as zero; callers may supply either a current value
    or an input-based estimate.
    """
    return 1.0 / max(abs(reference_value or 0.0), floor)
