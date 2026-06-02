#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Initializers and cascade helpers for the flotation bank unit model."""

import math

from pyomo.common.config import ConfigValue, In
from pyomo.environ import check_optimal_termination, value

from idaes import logger as idaeslog
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    ModularInitializerBase,
)
from idaes.core.util.exceptions import ConfigurationError, InitializationError

__author__ = "Daison Yancy Caballero"

FEED_FLOOR_KG_PER_H = 1e-12
REL_BRACKET_EPS = 1e-9
BRACKET_REFINE_MAX_ITER = 30
MARGIN_TOLERANCE = 1e-6
MARGIN_NEAR_BOUNDARY = 0.05
BRENTQ_XTOL = 1e-18
BRENTQ_RTOL = 1e-12
BOUND_PROXIMITY_FRAC = 0.10
POST_SEED_RESIDUAL_TOL = 1e-5
INLET_GUESS_SCALE_SAFETY = 1.05
INLET_GUESS_SCALE_MAX_ITER = 20
INLET_GUESS_SCALE_LIMIT = 1e12
STAGED_INLET_GUESS_SCALE_SAFETY = 2.0
STAGED_REMOVAL_GUESS_FRACTION = 0.5


class CascadeInputError(ValueError):
    """Static input or schema error in a kinetic cell-balance cascade."""

    def __init__(self, message, *, context=""):
        super().__init__(message)
        self.context = context


class CascadeInfeasibleError(ValueError):
    """Valid inputs that cannot support a positive kinetic cell-balance root."""

    def __init__(
        self,
        message,
        *,
        cell_index,
        F_zero,
        lower_margin,
        deficit,
        failure_mode,
        context="",
    ):
        details = (
            f"cell_index={cell_index}, F_zero={F_zero}, "
            f"lower_margin={lower_margin}, deficit={deficit}, "
            f"failure_mode={failure_mode}"
        )
        super().__init__(f"{message} ({details})")
        self.cell_index = cell_index
        self.F_zero = F_zero
        self.lower_margin = lower_margin
        self.deficit = deficit
        self.failure_mode = failure_mode
        self.context = context


def _is_finite_or_input_error(val, name, context):
    try:
        if not math.isfinite(val):
            raise CascadeInputError(
                f"{context}: non-finite {name}={val!r}", context=context
            )
    except TypeError as exc:
        raise CascadeInputError(
            f"{context}: non-numeric {name}={val!r} " f"(type {type(val).__name__})",
            context=context,
        ) from exc


def _adaptive_bracket(
    g,
    F_zero,
    F_total,
    cell_index,
    context,
    *,
    lower_margin=None,
    deficit=0.0,
):
    """Construct a root bracket with a strictly positive lower residual."""
    T_high = F_total
    if F_zero > FEED_FLOOR_KG_PER_H:
        lower_boundary = F_zero
        offset = REL_BRACKET_EPS * max(F_total - F_zero, F_zero)
    else:
        lower_boundary = 0.0
        offset = REL_BRACKET_EPS * F_total

    for _ in range(BRACKET_REFINE_MAX_ITER):
        T_low = lower_boundary + offset
        if T_low > 0 and g(T_low) > 0:
            return T_low, T_high
        offset *= 0.5

    raise CascadeInfeasibleError(
        f"{context}: cell {cell_index} near-singular bracket after "
        f"{BRACKET_REFINE_MAX_ITER} refinements "
        f"(failure_mode=near_singular_bracket, lower_boundary={lower_boundary}, "
        f"lower_margin={lower_margin}, deficit={deficit})",
        cell_index=cell_index,
        F_zero=F_zero,
        lower_margin=lower_margin,
        deficit=deficit,
        failure_mode="near_singular_bracket",
        context=context,
    )


def _solve_cell(
    F_in_cell_raw,
    M_total,
    k_hour,
    cell_index,
    context,
    *,
    collect_diagnostics,
):
    """Solve one cell in the fixed-inventory kinetic cascade."""
    from scipy.optimize import brentq

    F_in_cell = {
        j: (F_in_cell_raw[j] if F_in_cell_raw[j] > FEED_FLOOR_KG_PER_H else 0.0)
        for j in F_in_cell_raw
    }
    F_total = sum(F_in_cell.values())

    if F_total <= FEED_FLOOR_KG_PER_H:
        # The depleted_cell mode has no meaningful "margin deficit" (lower_margin
        # is undefined when there is no cell feed at all). We report
        # 1.0 + MARGIN_TOLERANCE as a positive sentinel so that audit summaries
        # consistently treat depleted cells as non-zero-deficit failures rather
        # than as borderline-feasible. Consumers should branch on failure_mode,
        # not on the numeric value.
        raise CascadeInfeasibleError(
            f"{context}: cell {cell_index} depleted below feed floor after "
            f"sanitization (failure_mode=depleted_cell, "
            f"F_total_sanitized={F_total:.6g}); positive M_total={M_total:.6g} "
            "cannot be supported by zero cell feed. Reduce upstream k, reduce "
            "upstream M_total, or increase upstream feed so this cell receives "
            "feed above the floor.",
            cell_index=cell_index,
            F_zero=0.0,
            lower_margin=None,
            deficit=1.0 + MARGIN_TOLERANCE,
            failure_mode="depleted_cell",
            context=context,
        )

    pos_k_with_feed = [
        j for j in F_in_cell if k_hour[j] > 0.0 and F_in_cell[j] > FEED_FLOOR_KG_PER_H
    ]
    zero_k_with_feed = [
        j for j in F_in_cell if k_hour[j] == 0.0 and F_in_cell[j] > FEED_FLOOR_KG_PER_H
    ]
    F_zero = sum(F_in_cell[j] for j in zero_k_with_feed)
    has_zero_k_carrier = bool(zero_k_with_feed)

    if not pos_k_with_feed:
        F_pulp_out = dict(F_in_cell)
        F_float = {j: 0.0 for j in F_in_cell}
        M_comp = {j: M_total * F_in_cell[j] / F_total for j in F_in_cell}
        cell_diag = None
        if collect_diagnostics:
            cell_diag = {
                "cell_index": cell_index,
                "F_total": F_total,
                "F_zero": F_zero,
                "has_zero_k_carrier": has_zero_k_carrier,
                "lower_margin": None,
                "deficit": 0.0,
                "T_low": None,
                "T_root": None,
            }
        return F_pulp_out, M_comp, F_float, cell_diag

    if F_zero <= FEED_FLOOR_KG_PER_H:
        lower_margin = sum(
            F_in_cell[j] / (k_hour[j] * M_total) for j in pos_k_with_feed
        )
        deficit = max(0.0, 1.0 + MARGIN_TOLERANCE - lower_margin)
        if lower_margin <= 1.0 + MARGIN_TOLERANCE:
            raise CascadeInfeasibleError(
                f"{context}: cell {cell_index} infeasible "
                f"(failure_mode=no_carrier_margin, no carrier, "
                f"lower_margin={lower_margin:.6g}, deficit={deficit:.6g}). "
                "Reduce k, reduce M_total, or increase F_in.",
                cell_index=cell_index,
                F_zero=F_zero,
                lower_margin=lower_margin,
                deficit=deficit,
                failure_mode="no_carrier_margin",
                context=context,
            )
    else:
        lower_margin = None
        deficit = 0.0

    def g(T):
        return sum(F_in_cell[j] / (T + k_hour[j] * M_total) for j in F_in_cell) - 1.0

    T_low, T_high = _adaptive_bracket(
        g,
        F_zero,
        F_total,
        cell_index,
        context,
        lower_margin=lower_margin,
        deficit=deficit,
    )
    T = brentq(g, T_low, T_high, xtol=BRENTQ_XTOL, rtol=BRENTQ_RTOL)

    F_pulp_out = {j: F_in_cell[j] * T / (T + k_hour[j] * M_total) for j in F_in_cell}
    M_comp = {j: M_total * F_pulp_out[j] / T for j in F_in_cell}
    F_float = {j: F_in_cell[j] - F_pulp_out[j] for j in F_in_cell}

    cell_diag = None
    if collect_diagnostics:
        cell_diag = {
            "cell_index": cell_index,
            "F_total": F_total,
            "F_zero": F_zero,
            "has_zero_k_carrier": has_zero_k_carrier,
            "lower_margin": lower_margin,
            "deficit": deficit,
            "T_low": T_low,
            "T_root": T,
        }
    return F_pulp_out, M_comp, F_float, cell_diag


def _cascade_forward_impl(
    F_in_bank,
    M_total_by_cell,
    k_hour,
    context,
    *,
    collect_diagnostics=False,
):
    """Shared cascade implementation with optional per-cell diagnostics."""
    if isinstance(M_total_by_cell, (str, int, float)):
        raise CascadeInputError(
            f"{context}: M_total_by_cell must be a sequence of length N; got "
            f"{type(M_total_by_cell).__name__}.",
            context=context,
        )
    if not hasattr(M_total_by_cell, "__len__"):
        raise CascadeInputError(
            f"{context}: M_total_by_cell must support len()", context=context
        )
    if len(M_total_by_cell) == 0:
        raise CascadeInputError(
            f"{context}: M_total_by_cell cannot be empty", context=context
        )

    if set(F_in_bank) != set(k_hour):
        raise CascadeInputError(
            f"{context}: F_in_bank and k_hour component sets differ; "
            f"missing={sorted(set(F_in_bank) - set(k_hour))}, "
            f"extra={sorted(set(k_hour) - set(F_in_bank))}",
            context=context,
        )

    for component, val in F_in_bank.items():
        _is_finite_or_input_error(val, f"F_in_bank[{component}]", context)
    for component, val in k_hour.items():
        _is_finite_or_input_error(val, f"k_hour[{component}]", context)
    for idx, M in enumerate(M_total_by_cell, start=1):
        _is_finite_or_input_error(M, f"M_total_by_cell[{idx}]", context)

    if any(M <= 0 for M in M_total_by_cell):
        raise CascadeInputError(
            f"{context}: M_total_by_cell entries must be positive; "
            f"got {list(M_total_by_cell)}",
            context=context,
        )
    if any(F < 0 for F in F_in_bank.values()):
        raise CascadeInputError(
            f"{context}: F_in_bank entries must be non-negative; "
            f"got {dict(F_in_bank)}",
            context=context,
        )
    if any(k_j < 0 for k_j in k_hour.values()):
        raise CascadeInputError(
            f"{context}: k_hour entries must be non-negative; got {dict(k_hour)}",
            context=context,
        )

    sanitized_bank_total = sum(
        F if F > FEED_FLOOR_KG_PER_H else 0.0 for F in F_in_bank.values()
    )
    if sanitized_bank_total <= FEED_FLOOR_KG_PER_H:
        raise CascadeInputError(
            f"{context}: zero total bank inlet feed after applying "
            f"FEED_FLOOR_KG_PER_H to each component "
            f"(raw_total={sum(F_in_bank.values()):.6g}, "
            f"sanitized_total={sanitized_bank_total:.6g}, "
            f"floor={FEED_FLOOR_KG_PER_H:.6g})",
            context=context,
        )

    F_pulp_out_list = []
    M_list = []
    F_float_list = []
    diagnostics_list = [] if collect_diagnostics else None
    F_in_cell = dict(F_in_bank)
    for cell_index, M_total in enumerate(M_total_by_cell, start=1):
        F_pulp_out, M_comp, F_float, cell_diag = _solve_cell(
            F_in_cell,
            M_total,
            k_hour,
            cell_index,
            context,
            collect_diagnostics=collect_diagnostics,
        )
        F_pulp_out_list.append(F_pulp_out)
        M_list.append(M_comp)
        F_float_list.append(F_float)
        if collect_diagnostics:
            diagnostics_list.append(cell_diag)
        F_in_cell = F_pulp_out

    return F_pulp_out_list, M_list, F_float_list, diagnostics_list


def cell_cascade_forward(F_in_bank, M_total_by_cell, k_hour, context=""):
    """Forward-march a fixed-inventory kinetic flotation cell cascade."""
    F_pulp_out, M_comp, F_float, _ = _cascade_forward_impl(
        F_in_bank,
        M_total_by_cell,
        k_hour,
        context,
        collect_diagnostics=False,
    )
    return F_pulp_out, M_comp, F_float


def audit_cascade_solvability(F_in_bank, M_total_by_cell, k_hour, bank_name=""):
    """Return per-cell solvability diagnostics for a kinetic cascade walk."""
    summary = {
        "bank": bank_name,
        "min_lower_margin_no_carrier": None,
        "min_margin_cell": None,
        "cells_with_zero_k_carrier": 0,
        "near_boundary_cells": 0,
        "first_infeasible_cell": None,
        "first_failure_mode": None,
        "walk_completed": False,
    }
    try:
        _, _, _, diagnostics = _cascade_forward_impl(
            F_in_bank,
            M_total_by_cell,
            k_hour,
            bank_name,
            collect_diagnostics=True,
        )
    except CascadeInfeasibleError as exc:
        summary["first_infeasible_cell"] = exc.cell_index
        summary["first_failure_mode"] = exc.failure_mode
        summary["min_margin_cell"] = {
            "cell_index": exc.cell_index,
            "lower_margin": exc.lower_margin,
            "deficit": exc.deficit,
            "failure_mode": exc.failure_mode,
        }
        if exc.failure_mode == "no_carrier_margin":
            summary["min_lower_margin_no_carrier"] = exc.lower_margin
        return summary
    except CascadeInputError:
        raise

    summary["walk_completed"] = True
    for diag in diagnostics:
        if diag["has_zero_k_carrier"]:
            summary["cells_with_zero_k_carrier"] += 1
        if diag["lower_margin"] is not None:
            cur = summary["min_lower_margin_no_carrier"]
            if cur is None or diag["lower_margin"] < cur:
                summary["min_lower_margin_no_carrier"] = diag["lower_margin"]
                summary["min_margin_cell"] = {
                    "cell_index": diag["cell_index"],
                    "lower_margin": diag["lower_margin"],
                    "deficit": diag["deficit"],
                }
            if abs(diag["lower_margin"] - 1.0) < MARGIN_NEAR_BOUNDARY:
                summary["near_boundary_cells"] += 1
    return summary


def _safe_value(obj, default):
    val = value(obj, exception=False)
    if val is None or not math.isfinite(val):
        return default
    return val


def _cell_feed_nominal(model, t, i, j):
    if i == model.cells.first():
        return _safe_value(model.properties_in[t].flow_mass_comp[j], default=1e-3)
    return _safe_value(
        model.cell_pulp_out_flow[t, model.cells.prev(i), j], default=1e-3
    )


_KINETIC_CELL_BALANCE_CONSTRAINT_NAMES = (
    "geometric_holdup_eq",
    "total_component_holdup_eq",
    "well_mixed_eq",
    "cell_material_balance_eq",
    "flotation_removal_eq",
    "concentrate_outlet_eq",
    "tails_outlet_eq",
    "recovery_eq",
)


def _kinetic_cell_balance_max_constraint_residual(model):
    """Return ``(max_residual, worst_constraint_name)`` over the kinetic
    cell-balance constraints. Used to detect drift between the helper and the
    Pyomo equations after pre-seed or staged solve."""
    max_residual = 0.0
    worst = None
    for constraint_name in _KINETIC_CELL_BALANCE_CONSTRAINT_NAMES:
        constraint_obj = getattr(model, constraint_name, None)
        if constraint_obj is None:
            continue
        for index in constraint_obj:
            con = constraint_obj[index]
            if not con.active:
                continue
            try:
                body_val = value(con.body, exception=False)
                upper_val = (
                    value(con.upper, exception=False) if con.upper is not None else 0.0
                )
            except (ValueError, TypeError):
                continue
            if body_val is None or upper_val is None:
                continue
            if not (math.isfinite(body_val) and math.isfinite(upper_val)):
                continue
            residual = abs(body_val - upper_val)
            if residual > max_residual:
                max_residual = residual
                worst = (
                    f"{constraint_name}"
                    f"{tuple(index) if isinstance(index, tuple) else (index,)}"
                )
    return max_residual, worst


def _apply_zero_feed_recovery_policy(model, time, component_list):
    """Fix zero-feed component recoveries to 0 and unfix recoveries that the
    initializer previously fixed once feed becomes positive again.

    Shared by both kinetic-cell-balance initializers. Mutates
    ``model._zero_feed_recovery_fixes``. Raises ``ConfigurationError`` if a
    component with zero feed has a user-fixed nonzero recovery.
    """
    for component in component_list:
        feed = value(model.properties_in[time].flow_mass_comp[component])
        recovery = model.recovery[time, component]
        if feed < FEED_FLOOR_KG_PER_H:
            if recovery.fixed:
                recovery_value = value(recovery)
                if abs(recovery_value) > 1e-12:
                    raise ConfigurationError(
                        f"bank={model.local_name}: zero F_in at "
                        f"(t={time}, j={component}) but recovery is fixed "
                        f"to {recovery_value}. Unfix recovery or provide "
                        "positive feed."
                    )
            else:
                recovery.fix(0.0)
                model._zero_feed_recovery_fixes.add((time, component))
        elif (time, component) in model._zero_feed_recovery_fixes:
            recovery.unfix()
            model._zero_feed_recovery_fixes.discard((time, component))


def _geometric_m_total(model):
    """Evaluate the geometric solid inventory ``M_total`` from the bank's
    operating Vars. Raises ``ConfigurationError`` if non-finite or non-positive.
    """
    rho_slurry_val = value(model.rho_slurry)
    M_total = (
        value(model.cell_volume)
        * (1 - value(model.air_holdup))
        * rho_slurry_val
        * value(model.pulp_solids_mass_fraction)
    )
    if not math.isfinite(M_total) or M_total <= 0:
        raise ConfigurationError(
            f"bank={model.local_name}: geometric M_total is {M_total}; "
            "check cell_volume, air_holdup, pulp_solids_mass_fraction, "
            "rho_solid, and rho_water are fixed and finite."
        )
    return M_total


def _cell_balance_infeasibility_message(model, exc):
    lower_margin_text = (
        "None" if exc.lower_margin is None else f"{exc.lower_margin:.6g}"
    )
    return (
        f"Kinetic-cell-balance bank={model.local_name} is algebraically "
        f"infeasible at cell {exc.cell_index} "
        f"(failure_mode={exc.failure_mode}, F_zero={exc.F_zero:.6g}, "
        f"lower_margin={lower_margin_text}, deficit={exc.deficit:.6g}). "
        "To restore feasibility, reduce k, reduce M_total (lower "
        "cell_volume, lower pulp_solids_mass_fraction, or lower "
        "rho_solid), or increase bank feed."
    )


def _cascade_guess_scale_factor(exc):
    if exc.lower_margin is not None and exc.lower_margin > 0:
        return max(
            2.0,
            (1.0 + MARGIN_TOLERANCE) / exc.lower_margin * INLET_GUESS_SCALE_SAFETY,
        )
    return 10.0


def _feasible_cell_balance_cascade_seed(
    model,
    time,
    component_list,
    M_total_by_cell,
    k_hour,
):
    """Return a feasible forward cascade, scaling caller-marked inlet guesses.

    SequentialDecomposition may call unit initializers before recycle/tear
    stream guesses have physical magnitudes. For fixed-inventory cell balances,
    those tiny guesses can be algebraically infeasible even though the connected
    flowsheet is feasible. In that case the flowsheet can mark variables that
    were originally unfixed before SequentialDecomposition temporarily fixed
    them. Only those marked inlet values are scaled. Direct user data still
    fails hard.
    """
    context = f"bank={model.local_name}"
    inlet_vars = {
        component: model.properties_in[time].flow_mass_comp[component]
        for component in component_list
    }
    F_in_bank = {
        component: value(var, exception=False) for component, var in inlet_vars.items()
    }
    try:
        F_pulp_out_cells, M_cells, F_float_cells = cell_cascade_forward(
            F_in_bank,
            M_total_by_cell,
            k_hour,
            context=context,
        )
        return F_in_bank, F_pulp_out_cells, M_cells, F_float_cells
    except CascadeInputError:
        raise
    except CascadeInfeasibleError as first_exc:
        initially_unfixed = getattr(
            model.flowsheet(),
            "_kinetic_cell_balance_initially_unfixed",
            None,
        )
        adjustable = {}
        if initially_unfixed is not None:
            adjustable = {
                component: var
                for component, var in inlet_vars.items()
                if id(var) in initially_unfixed
            }
        if not adjustable:
            raise first_exc

        trial_feed = dict(F_in_bank)
        scale_factor = _cascade_guess_scale_factor(first_exc)
        last_exc = first_exc
        for _ in range(INLET_GUESS_SCALE_MAX_ITER):
            for component, var in adjustable.items():
                base = value(var, exception=False)
                if base is None or not math.isfinite(base):
                    base = 1.0
                trial_feed[component] = (
                    max(float(base), 10.0 * FEED_FLOOR_KG_PER_H) * scale_factor
                )

            try:
                F_pulp_out_cells, M_cells, F_float_cells = cell_cascade_forward(
                    trial_feed,
                    M_total_by_cell,
                    k_hour,
                    context=context,
                )
            except CascadeInputError:
                raise
            except CascadeInfeasibleError as exc:
                last_exc = exc
                scale_factor *= _cascade_guess_scale_factor(exc)
                if scale_factor > INLET_GUESS_SCALE_LIMIT:
                    break
            else:
                for component, var in adjustable.items():
                    var.set_value(trial_feed[component])
                return trial_feed, F_pulp_out_cells, M_cells, F_float_cells

        raise last_exc


def _prepare_kinetic_cell_balance_inlet_guesses(model):
    """Scale marked inlet guesses and return the feasible cascade seed.

    ``ModularInitializerBase`` restores fixed variable values after
    ``initialization_routine``. SequentialDecomposition may temporarily fix
    inlet guesses, so those guesses must be made feasible before the IDAES
    initializer records its initial fixed-state snapshot. The returned seed is
    reused inside ``initialization_routine`` to avoid solving the analytical
    cascade twice.
    """
    time_set = list(model.flowsheet().time)
    if len(time_set) != 1:
        raise ConfigurationError(
            f"bank={model.local_name}: kinetic_cell_balance steady-state "
            f"initializer requires a singleton time set; got {len(time_set)} "
            "time points. Future dynamic support will lift this restriction."
        )
    time = time_set[0]
    component_list = list(model.config.property_package.component_list)
    M_total = _geometric_m_total(model)
    M_total_by_cell = [M_total for _ in model.cells]
    k_hour = {
        component: value(model.k_cb[time, component]) for component in component_list
    }
    try:
        (
            F_in_bank,
            F_pulp_out_cells,
            M_cells,
            F_float_cells,
        ) = _feasible_cell_balance_cascade_seed(
            model,
            time,
            component_list,
            M_total_by_cell,
            k_hour,
        )
    except CascadeInputError as exc:
        raise ConfigurationError(
            "Static input error in kinetic-cell-balance initialization for "
            f"bank={model.local_name}: {exc}"
        ) from exc
    except CascadeInfeasibleError as exc:
        raise ConfigurationError(
            _cell_balance_infeasibility_message(model, exc)
        ) from exc
    return {
        "model_id": id(model),
        "time": time,
        "component_list": component_list,
        "M_total_by_cell": M_total_by_cell,
        "k_hour": k_hour,
        "F_in_bank": F_in_bank,
        "F_pulp_out_cells": F_pulp_out_cells,
        "M_cells": M_cells,
        "F_float_cells": F_float_cells,
    }


def _get_prepared_kinetic_cell_balance_seed(initializer, model):
    seed = getattr(initializer, "_kinetic_cell_balance_prepared_seed", None)
    if seed is not None and seed.get("model_id") == id(model):
        return seed
    return _prepare_kinetic_cell_balance_inlet_guesses(model)


def _staged_inlet_guess_floor(component, M_total_by_cell, k_hour):
    return max(
        10.0 * FEED_FLOOR_KG_PER_H,
        STAGED_INLET_GUESS_SCALE_SAFETY
        * (len(M_total_by_cell) + 1)
        * max(M_total_by_cell)
        * k_hour[component],
    )


def _prepare_kinetic_cell_balance_staged_inlet_guesses(model):
    """Scale marked inlet guesses and cache static staged-initialization data.

    This intentionally does not call ``cell_cascade_forward``. The staged
    initializer should use IPOPT on the Pyomo equations to solve the cascade,
    so this pre-snapshot step only guards temporary SequentialDecomposition
    inlet guesses against being too small for fixed cell inventories.
    """
    time_set = list(model.flowsheet().time)
    if len(time_set) != 1:
        raise ConfigurationError(
            f"bank={model.local_name}: kinetic_cell_balance steady-state "
            f"initializer requires a singleton time set; got {len(time_set)} "
            "time points. Future dynamic support will lift this restriction."
        )
    time = time_set[0]
    component_list = list(model.config.property_package.component_list)
    M_total = _geometric_m_total(model)
    M_total_by_cell = [M_total for _ in model.cells]
    k_hour = {
        component: value(model.k_cb[time, component]) for component in component_list
    }

    context = f"bank={model.local_name}"
    inlet_vars = {
        component: model.properties_in[time].flow_mass_comp[component]
        for component in component_list
    }
    F_in_bank = {
        component: value(var, exception=False) for component, var in inlet_vars.items()
    }

    for component, k_value in k_hour.items():
        try:
            if not math.isfinite(k_value) or k_value < 0:
                raise ConfigurationError(
                    f"{context}: k_hour[{component}] must be finite and "
                    f"non-negative; got {k_value!r}."
                )
        except TypeError as exc:
            raise ConfigurationError(
                f"{context}: k_hour[{component}] must be numeric; " f"got {k_value!r}."
            ) from exc

    initially_unfixed = getattr(
        model.flowsheet(),
        "_kinetic_cell_balance_initially_unfixed",
        None,
    )
    if initially_unfixed is not None:
        for component, var in inlet_vars.items():
            if id(var) not in initially_unfixed:
                continue
            feed_value = F_in_bank[component]
            try:
                invalid_feed_guess = (
                    feed_value is None
                    or not math.isfinite(feed_value)
                    or feed_value < 0
                )
            except TypeError:
                invalid_feed_guess = True
            if invalid_feed_guess:
                feed_value = 0.0
            scaled_value = max(
                float(feed_value),
                _staged_inlet_guess_floor(component, M_total_by_cell, k_hour),
            )
            if scaled_value > float(feed_value):
                var.set_value(scaled_value)
                F_in_bank[component] = scaled_value

    for component, feed_value in F_in_bank.items():
        try:
            if not math.isfinite(feed_value) or feed_value < 0:
                raise ConfigurationError(
                    f"{context}: F_in_bank[{component}] must be finite and "
                    f"non-negative; got {feed_value!r}."
                )
        except TypeError as exc:
            raise ConfigurationError(
                f"{context}: F_in_bank[{component}] must be numeric; "
                f"got {feed_value!r}."
            ) from exc

    sanitized_total = sum(
        feed_value if feed_value > FEED_FLOOR_KG_PER_H else 0.0
        for feed_value in F_in_bank.values()
    )
    if sanitized_total <= FEED_FLOOR_KG_PER_H:
        raise ConfigurationError(
            f"{context}: zero total bank inlet feed after applying "
            f"FEED_FLOOR_KG_PER_H to each component "
            f"(raw_total={sum(F_in_bank.values()):.6g}, "
            f"sanitized_total={sanitized_total:.6g}, "
            f"floor={FEED_FLOOR_KG_PER_H:.6g})."
        )

    return {
        "model_id": id(model),
        "time": time,
        "component_list": component_list,
        "M_total_by_cell": M_total_by_cell,
        "k_hour": k_hour,
        "F_in_bank": F_in_bank,
    }


def _get_prepared_kinetic_cell_balance_staged_seed(initializer, model):
    seed = getattr(initializer, "_kinetic_cell_balance_prepared_seed", None)
    if seed is not None and seed.get("model_id") == id(model):
        return seed
    return _prepare_kinetic_cell_balance_staged_inlet_guesses(model)


class FlotationBankKineticClosedFormInitializer(BlockTriangularizationInitializer):
    """Initializer for closed-form kinetic flotation banks."""

    def initialization_routine(self, model):
        """Seed dependent kinetic variables before block triangularization."""
        component_list = list(model.config.property_package.component_list)
        for time in model.flowsheet().time:
            total_in = sum(
                value(model.properties_in[time].flow_mass_comp[component])
                for component in component_list
            )
            if total_in <= 0:
                raise ConfigurationError(
                    "Kinetic mode requires positive inlet dry-solids flow; "
                    "zero-flow banks cannot evaluate residence time."
                )

        number_of_cells = model.config.number_of_cells
        for time in model.flowsheet().time:
            tau_value = value(model.tau[time])
            for component in component_list:
                recovery_fixed = model.recovery[time, component].fixed
                k_fixed = model.k_cf[time, component].fixed
                if recovery_fixed or not k_fixed:
                    raise ConfigurationError(
                        "Unsupported DOF state for kinetic bank at "
                        f"(t={time}, j={component}): "
                        f"recovery.fixed={recovery_fixed}, k_cf.fixed={k_fixed}. "
                        "Closed-form kinetic initialization supports prediction "
                        "mode only; fix k_cf and leave recovery free. For "
                        "parameter estimation, estimate k_cf at the flowsheet "
                        "level instead of fixing recovery in the unit initializer."
                    )
                r_inf_value = value(model.R_inf[time, component])
                k_value = value(model.k_cf[time, component])
                recovery_value = r_inf_value * (
                    1
                    - (1 + k_value * tau_value / number_of_cells) ** (-number_of_cells)
                )
                model.recovery[time, component].set_value(recovery_value)

        return super().initialization_routine(model)


class FlotationBankKineticCellBalanceAnalyticalInitializer(ModularInitializerBase):
    """Analytical pre-seed initializer for kinetic cell-balance banks.

    Strategy
    --------
    The kinetic cell-balance equations admit an exact analytical solution per
    cell once geometry, pulp parameters, and rate constants are fixed (see
    ``cell_cascade_forward``). This initializer exploits that by:

    1. Validating inputs (singleton time, positive bank inlet, finite geometric
       ``M_total``, consistent zero-feed ``recovery`` policy).
    2. Calling ``cell_cascade_forward`` to compute, in pure Python, the
       analytical answer for every cell-level variable.
    3. Calling ``.set_value(...)`` on each Pyomo Var so the model is in a
       constraint-satisfying state to floating-point precision.
    4. Verifying that the pre-seed actually satisfies every kinetic
       cell-balance constraint within ``POST_SEED_RESIDUAL_TOL``. A failure
       here indicates the helper and Pyomo equations have drifted out of sync.

    Infeasibility is caught at step 2 (``CascadeInfeasibleError``) and surfaced
    as ``ConfigurationError`` with bank/cell context before any Pyomo solve is
    attempted. The analytical pre-seed makes a block-triangularization sweep
    unnecessary.

    This is the ``strategy="analytical"`` worker behind the default
    :class:`FlotationBankKineticCellBalanceInitializer`. It can also be used
    directly. For the IDAES staged-solve alternative, see
    :class:`FlotationBankKineticCellBalanceStagedInitializer`.
    """

    def initialize(self, model, *args, **kwargs):
        self._kinetic_cell_balance_prepared_seed = (
            _prepare_kinetic_cell_balance_inlet_guesses(model)
        )
        try:
            return super().initialize(model, *args, **kwargs)
        finally:
            self._kinetic_cell_balance_prepared_seed = None

    def initialization_routine(self, model):
        seed = _get_prepared_kinetic_cell_balance_seed(self, model)
        time = seed["time"]
        component_list = seed["component_list"]

        F_in_bank = seed["F_in_bank"]
        total_in = sum(F_in_bank.values())
        if total_in <= FEED_FLOOR_KG_PER_H:
            raise ConfigurationError(
                f"bank={model.local_name}: zero total inlet flow at t={time}."
            )

        _apply_zero_feed_recovery_policy(model, time, component_list)

        M_total_by_cell = seed["M_total_by_cell"]
        F_pulp_out_cells = seed["F_pulp_out_cells"]
        M_cells = seed["M_cells"]
        F_float_cells = seed["F_float_cells"]

        for cell_idx, cell in enumerate(model.cells):
            model.cell_total_solid_holdup[time, cell].set_value(
                M_total_by_cell[cell_idx]
            )
            for component in component_list:
                model.cell_solid_holdup[time, cell, component].set_value(
                    M_cells[cell_idx][component]
                )
                model.cell_pulp_out_flow[time, cell, component].set_value(
                    F_pulp_out_cells[cell_idx][component]
                )
                model.cell_float_flow[time, cell, component].set_value(
                    F_float_cells[cell_idx][component]
                )

        for component in component_list:
            concentrate = sum(
                F_float_cells[cell_idx][component]
                for cell_idx in range(len(F_float_cells))
            )
            feed = F_in_bank[component]
            model.properties_concentrate[time].flow_mass_comp[component].set_value(
                concentrate
            )
            model.properties_tails[time].flow_mass_comp[component].set_value(
                F_pulp_out_cells[-1][component]
            )
            if feed > FEED_FLOOR_KG_PER_H:
                model.recovery[time, component].set_value(concentrate / feed)

        # Verify the analytical pre-seed satisfies every kinetic cell-balance
        # constraint within POST_SEED_RESIDUAL_TOL. A failure here means
        # cell_cascade_forward and the Pyomo equations have drifted out of sync.
        max_residual, worst_constraint = _kinetic_cell_balance_max_constraint_residual(
            model
        )
        if max_residual > POST_SEED_RESIDUAL_TOL:
            raise ConfigurationError(
                f"bank={model.local_name}: kinetic-cell-balance pre-seed residual "
                f"{max_residual:.3g} exceeds POST_SEED_RESIDUAL_TOL="
                f"{POST_SEED_RESIDUAL_TOL:.3g} at {worst_constraint}. "
                "cell_cascade_forward and the Pyomo equations have drifted out "
                "of sync; check the helper implementation against the Section 4 "
                "constraints in the plan."
            )

        return None


class FlotationBankKineticCellBalanceStagedInitializer(ModularInitializerBase):
    """Staged-solve initializer for kinetic cell-balance banks (IDAES standard).

    Uses common IDAES fix / deactivate / solve / unfix / reactivate
    pattern with an independent, non-analytical starting point. This path does
    not call ``cell_cascade_forward`` for the cell cascade values; the staged
    IPOPT solves compute those values from the Pyomo equations themselves.

    NOTE: This initializer is an optional alternative to the default analytical
    initializer for kinetic cell-balance banks. It has not been tested
    extensively yet.

    Strategy
    --------
    1. Validate inputs (singleton time set, positive bank inlet, finite
       geometric ``M_total``) and apply the zero-feed ``recovery`` policy.
    2. Set physically motivated rough guesses for every cell-level Var from
       the current cell feed composition. If the current inlet values are
       caller-marked temporary initialization guesses and are too small for the
       fixed inventory, scale only those marked inlet guesses before solving.
    3. Deactivate the bank-level outlet equations (``concentrate_outlet_eq``,
       ``tails_outlet_eq``, ``recovery_eq``), fix the corresponding bank-level
       outlet Vars (``properties_concentrate.flow_mass_comp``,
       ``properties_tails.flow_mass_comp``, ``recovery``) to their stage-2
       guesses, and solve the resulting cell-balance sub-block with IPOPT. This
       isolates the bilinear ``well_mixed_eq`` and the per-cell material balance
       / flotation-removal equations.
    4. Unfix the bank-level outlet Vars, reactivate the outlet equations,
       seed the outlet Vars from the solved cell-level values, and solve the
       full bank with IPOPT. This stage is typically trivial because the outlet
       equations just aggregate solved cell values.
    5. Verify constraint residuals are within ``POST_SEED_RESIDUAL_TOL``.

    Any failure during stages 3-4 restores the original constraint-activation
    and Var-fix state before re-raising.

    Trade-offs vs. the default analytical initializer
    -------------------------------------------------
    - Slower: two IPOPT solves per bank versus zero.
    - More independent: every solve uses the Pyomo equations themselves after
      the initializer supplies only rough local values, not analytical cascade
      values.
    - Diagnostic clarity: a stage-3 IPOPT failure localizes the issue to the
      bilinear cell-balance block, while a stage-4 failure points at outlet
      aggregation.
    - Conforms to IDAES idioms more directly. Select it either through the
      default dispatcher, ``FlotationBankKineticCellBalanceInitializer(
      strategy="staged")``, or by using this class directly via
      ``FlotationBankKineticCellBalanceStagedInitializer().initialize(unit)``.
    - Uses the IDAES initializer solver configuration. Pass ``solver=...``,
      ``solver_options={...}``, or ``writer_config={...}`` when constructing
      the initializer to override defaults. Pass
      ``output_level=idaeslog.INFO_HIGH`` to ``initialize`` to show the major
      internal step messages.
    """

    _OUTLET_CONSTRAINT_NAMES = (
        "concentrate_outlet_eq",
        "tails_outlet_eq",
        "recovery_eq",
    )

    def initialize(self, model, *args, **kwargs):
        self._kinetic_cell_balance_prepared_seed = (
            _prepare_kinetic_cell_balance_staged_inlet_guesses(model)
        )
        try:
            return super().initialize(model, *args, **kwargs)
        finally:
            self._kinetic_cell_balance_prepared_seed = None

    def _set_initial_guesses(
        self,
        model,
        time,
        component_list,
        M_total_by_cell,
        F_in_bank,
        k_hour,
    ):
        """Stage-2 rough guesses for cell-level Vars and bank outlets.

        Sets values via ``.set_value(...)``; does not fix anything.
        """
        F_cell_in = dict(F_in_bank)
        F_float_cells = []
        F_pulp_out_cells = []

        for cell_idx, cell in enumerate(model.cells):
            M_total = M_total_by_cell[cell_idx]
            model.cell_total_solid_holdup[time, cell].set_value(M_total)
            total_cell_in = sum(
                max(F_cell_in[j], 0.0) if F_cell_in[j] > FEED_FLOOR_KG_PER_H else 0.0
                for j in component_list
            )
            if total_cell_in <= FEED_FLOOR_KG_PER_H:
                total_cell_in = sum(max(F_cell_in[j], 0.0) for j in component_list)
            if total_cell_in <= FEED_FLOOR_KG_PER_H:
                total_cell_in = len(component_list) * FEED_FLOOR_KG_PER_H

            F_float = {}
            F_pulp_out = {}
            for j in component_list:
                feed = max(F_cell_in[j], 0.0)
                holdup = M_total * feed / total_cell_in
                float_guess = min(
                    k_hour[j] * holdup,
                    STAGED_REMOVAL_GUESS_FRACTION * feed,
                )
                pulp_out = feed - float_guess

                model.cell_solid_holdup[time, cell, j].set_value(holdup)
                model.cell_pulp_out_flow[time, cell, j].set_value(pulp_out)
                model.cell_float_flow[time, cell, j].set_value(float_guess)
                F_float[j] = float_guess
                F_pulp_out[j] = pulp_out

            F_float_cells.append(F_float)
            F_pulp_out_cells.append(F_pulp_out)
            F_cell_in = F_pulp_out

        for j in component_list:
            concentrate = sum(
                F_float_cells[cell_idx][j] for cell_idx in range(len(F_float_cells))
            )
            model.properties_concentrate[time].flow_mass_comp[j].set_value(concentrate)
            model.properties_tails[time].flow_mass_comp[j].set_value(
                F_pulp_out_cells[-1][j]
            )
            if F_in_bank[j] > FEED_FLOOR_KG_PER_H:
                model.recovery[time, j].set_value(concentrate / F_in_bank[j])

    def _seed_outlet_vars_from_cells(self, model, time, component_list):
        """Stage-4 seed: aggregate solved cell-level values into bank outlets."""
        for j in component_list:
            concentrate = sum(
                value(model.cell_float_flow[time, cell, j]) for cell in model.cells
            )
            tails = value(model.cell_pulp_out_flow[time, model.cells.last(), j])
            model.properties_concentrate[time].flow_mass_comp[j].set_value(concentrate)
            model.properties_tails[time].flow_mass_comp[j].set_value(tails)
            feed = value(model.properties_in[time].flow_mass_comp[j])
            if feed > FEED_FLOOR_KG_PER_H:
                model.recovery[time, j].set_value(concentrate / feed)

    def _collect_outlet_constraints(self, model):
        """Return the list of active outlet ConstraintData objects to
        deactivate during stage 3."""
        active = []
        for name in self._OUTLET_CONSTRAINT_NAMES:
            con = getattr(model, name, None)
            if con is None:
                continue
            for idx in con:
                cdata = con[idx]
                if cdata.active:
                    active.append(cdata)
        return active

    def _collect_unfixed_outlet_vars(self, model, time, component_list):
        """Return the list of currently-unfixed bank-level outlet VarData
        objects to fix during stage 3."""
        unfixed = []
        for j in component_list:
            for var in (
                model.properties_concentrate[time].flow_mass_comp[j],
                model.properties_tails[time].flow_mass_comp[j],
                model.recovery[time, j],
            ):
                if not var.fixed:
                    unfixed.append(var)
        return unfixed

    def _run_solver(self, model, solver, solve_log, init_log, stage_name):
        init_log.info_high(f"Initialization {stage_name}: solve started.")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
        if not check_optimal_termination(results):
            raise InitializationError(
                f"bank={model.local_name}: kinetic-cell-balance staged "
                f"initializer {stage_name} did not converge: "
                f"termination_condition="
                f"{results.solver.termination_condition!r}, "
                f"status={results.solver.status!r}."
            )
        init_log.info_high(
            f"Initialization {stage_name} {idaeslog.condition(results)}."
        )
        return results

    def _verify_residuals(self, model, stage_name, init_log):
        max_residual, worst_constraint = _kinetic_cell_balance_max_constraint_residual(
            model
        )
        if max_residual > POST_SEED_RESIDUAL_TOL:
            raise InitializationError(
                f"bank={model.local_name}: staged-solve {stage_name} residual "
                f"{max_residual:.3g} exceeds POST_SEED_RESIDUAL_TOL="
                f"{POST_SEED_RESIDUAL_TOL:.3g} at {worst_constraint}."
            )
        init_log.info_high(
            f"Initialization {stage_name}: residual check complete "
            f"(max residual={max_residual:.3g})."
        )

    def initialization_routine(self, model):
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Stage 1: validate.
        seed = _get_prepared_kinetic_cell_balance_staged_seed(self, model)
        time = seed["time"]
        component_list = seed["component_list"]

        F_in_bank = seed["F_in_bank"]
        total_in = sum(F_in_bank.values())
        if total_in <= FEED_FLOOR_KG_PER_H:
            raise ConfigurationError(
                f"bank={model.local_name}: zero total inlet flow at t={time}."
            )

        _apply_zero_feed_recovery_policy(model, time, component_list)
        init_log.info_high("Initialization Step 1: validation complete.")

        M_total_by_cell = seed["M_total_by_cell"]
        k_hour = seed["k_hour"]

        # Stage 2: initial guesses.
        self._set_initial_guesses(
            model,
            time,
            component_list,
            M_total_by_cell,
            F_in_bank,
            k_hour,
        )
        init_log.info_high("Initialization Step 2: rough staged seed applied.")

        # Track state for rollback if any stage fails.
        deactivated_outlet_cons = []
        temp_fixed_outlet_vars = []
        try:
            # Stage 3: deactivate outlet equations, fix bank-level outlet Vars,
            # solve the cell-balance sub-block.
            deactivated_outlet_cons = self._collect_outlet_constraints(model)
            for cdata in deactivated_outlet_cons:
                cdata.deactivate()

            temp_fixed_outlet_vars = self._collect_unfixed_outlet_vars(
                model, time, component_list
            )
            for var in temp_fixed_outlet_vars:
                var.fix()

            solver = self._get_solver()
            self._run_solver(
                model,
                solver,
                solve_log,
                init_log,
                "Step 3 (cell-balance sub-block)",
            )
            self._verify_residuals(
                model,
                "Step 3 (cell-balance sub-block)",
                init_log,
            )

            # Stage 4: restore outlet block and solve full bank.
            for var in temp_fixed_outlet_vars:
                var.unfix()
            temp_fixed_outlet_vars = []
            for cdata in deactivated_outlet_cons:
                cdata.activate()
            deactivated_outlet_cons = []

            self._seed_outlet_vars_from_cells(model, time, component_list)

            self._run_solver(
                model,
                solver,
                solve_log,
                init_log,
                "Step 4 (full bank)",
            )

            # Stage 5: verify residuals.
            self._verify_residuals(model, "Step 5 (full bank)", init_log)
            init_log.info_high("Initialization Complete.")

            return None
        finally:
            # Roll back any state changes that did not get cleared in the
            # success path (e.g., when stage 3 raised).
            for var in temp_fixed_outlet_vars:
                var.unfix()
            for cdata in deactivated_outlet_cons:
                cdata.activate()


class FlotationBankKineticCellBalanceInitializer(ModularInitializerBase):
    """Default initializer for kinetic cell-balance banks (strategy dispatcher).

    This is the class assigned to ``default_initializer`` on a
    ``recovery_basis="kinetic_cell_balance"`` bank. It selects between the two
    concrete strategies based on the ``strategy`` configuration option and
    delegates the whole initialization to the chosen worker:

    * ``strategy="analytical"`` (default) ->
      :class:`FlotationBankKineticCellBalanceAnalyticalInitializer` -- a fast,
      solver-free analytical pre-seed (``cell_cascade_forward``).
    * ``strategy="staged"`` ->
      :class:`FlotationBankKineticCellBalanceStagedInitializer` -- the IDAES
      staged fix / deactivate / solve / unfix / reactivate pattern.

    Usage::

        # default (analytical) -- used by automated flowsheet initialization
        unit.default_initializer().initialize(unit)

        # opt into the staged strategy at the call site
        unit.default_initializer(strategy="staged").initialize(unit)

    Any constructor keyword other than ``strategy`` (e.g. ``solver``,
    ``solver_options``, ``writer_config``) is forwarded verbatim to the chosen
    worker, and ``initialize`` keywords (e.g. ``output_level``) pass straight
    through. The worker's :attr:`summary` is copied back onto this dispatcher so
    ``initializer.summary[model]`` works as usual.

    The two strategy classes can also be imported and used directly when finer
    control is wanted (for example, choosing the staged strategy for one bank in
    a flowsheet).
    """

    CONFIG = ModularInitializerBase.CONFIG()
    CONFIG.declare(
        "strategy",
        ConfigValue(
            default="analytical",
            domain=In(["analytical", "staged"]),
            doc="Which kinetic cell-balance initialization strategy to run: "
            "'analytical' (default, solver-free pre-seed) or 'staged' "
            "(IDAES staged solve).",
        ),
    )

    _STRATEGY_CLASSES = {
        "analytical": FlotationBankKineticCellBalanceAnalyticalInitializer,
        "staged": FlotationBankKineticCellBalanceStagedInitializer,
    }

    def __init__(self, **kwargs):
        # Capture the original constructor kwargs so non-strategy options can be
        # forwarded verbatim to the chosen worker initializer.
        self._init_kwargs = dict(kwargs)
        super().__init__(**kwargs)

    def initialize(self, model, *args, **kwargs):
        worker_cls = self._STRATEGY_CLASSES[self.config.strategy]
        worker_kwargs = {
            key: value for key, value in self._init_kwargs.items() if key != "strategy"
        }
        worker = worker_cls(**worker_kwargs)
        try:
            return worker.initialize(model, *args, **kwargs)
        finally:
            # Surface the worker's results on the dispatcher so callers can read
            # ``dispatcher.summary[model]`` as with any initializer.
            self.summary = worker.summary
