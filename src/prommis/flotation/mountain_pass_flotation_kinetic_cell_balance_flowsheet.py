####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Mountain Pass bastnaesite flotation flowsheet with kinetic cell-balance recoveries.

This module implements the five-bank Mountain Pass flotation circuit using the
kinetic cell-balance recovery model (``recovery_basis="kinetic_cell_balance"``).
Each bank is modeled as an ``N``-cell cascade of perfectly mixed flotation cells
with explicit per-cell solid holdup. The per-cell flotation flow is first order
in component holdup:

    F_float[i, j] = k_cb[j] * M[i, j]

where ``k_cb`` is the cell-balance first-order flotation rate constant
calibrated against a fixed-inventory CSTR cascade (geometric solid inventory
set by cell volume, air volume fraction, and pulp density), and ``M[i, j]`` is the
component solid holdup in cell ``i`` determined by the well-mixed closure and
the geometric holdup equation.

The fitted rate constants are loaded from
``mountain_pass_kinetic_cell_balance_parameters.json`` and are tied to the
``table1_product_fit`` scenario basis, bank geometry, and pulp-density
assumptions.

Reference: Pradip and Fuerstenau (2013), "Design and development of novel
flotation reagents for the beneficiation of Mountain Pass rare-earth ore,"
Minerals & Metallurgical Processing 30(1), 1-9.
"""

import functools
import json
import math
from pathlib import Path

from pyomo.core.expr.visitor import identify_variables
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TransformationFactory,
    Var,
    units,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.scaling import set_scaling_factor
from idaes.core.solvers import get_solver

from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation.bastnaesite_properties import (
    BastnaesiteParameters,
    COMPONENTS,
    SOLID_DENSITIES,
)
from prommis.flotation.flotation_bank import FlotationBank
from prommis.flotation.initializer import (
    CascadeInputError,
    CascadeInfeasibleError,
    audit_cascade_solvability,
)

__author__ = "Daison Yancy Caballero"

TABLE1_PRODUCT_FIT_SCENARIO = recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
UNIT_MODEL_DEFAULT_INITIALIZATION = recovery_flowsheet.UNIT_MODEL_DEFAULT_INITIALIZATION
DETERMINISTIC_STREAM_INITIALIZATION = (
    recovery_flowsheet.DETERMINISTIC_STREAM_INITIALIZATION
)
INITIALIZATION_METHODS = recovery_flowsheet.INITIALIZATION_METHODS
BANK_NAMES = recovery_flowsheet.BANK_NAMES
FINAL_SOLVER_OPTIONS = {
    "nlp_scaling_method": "user-scaling",
    "mu_strategy": "adaptive",
}
KINETIC_CELL_BALANCE_PARAMETERS_FILE = Path(__file__).with_name(
    "mountain_pass_kinetic_cell_balance_parameters.json"
)


@functools.lru_cache(maxsize=4)
def _load_kinetic_cell_balance_parameters_cached(path):
    with Path(path).open(encoding="utf-8") as file:
        return json.load(file)


def load_kinetic_cell_balance_parameters(path=None):
    """Load the kinetic cell-balance parameter JSON.

    Results are memoized per resolved path. The returned dict is shared across
    callers; do not mutate it.
    """
    parameter_file = (
        KINETIC_CELL_BALANCE_PARAMETERS_FILE if path is None else Path(path)
    )
    return _load_kinetic_cell_balance_parameters_cached(str(parameter_file))


def build_model(expand_arcs=True):
    """Build the Mountain Pass flowsheet with kinetic cell-balance banks."""
    parameters = load_kinetic_cell_balance_parameters()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    for component, density in parameters["property_package"][
        "rho_mass_comp_kg_per_m3"
    ].items():
        m.fs.properties.rho_mass_comp[component].set_value(density)
    m.fs.scenario = TABLE1_PRODUCT_FIT_SCENARIO
    recovery_flowsheet.add_conditioning_metadata(m.fs)

    m.fs.conditioning = recovery_flowsheet.StateJunction(
        property_package=m.fs.properties
    )
    m.fs.rougher_mixer = recovery_flowsheet.dry_solids_mixer(
        m.fs.properties, ["fresh_feed", "regrind_recycle"]
    )
    m.fs.cleaner1_mixer = recovery_flowsheet.dry_solids_mixer(
        m.fs.properties, ["rougher_concentrate", "cleaner2_tails"]
    )
    m.fs.cleaner2_mixer = recovery_flowsheet.dry_solids_mixer(
        m.fs.properties, ["cleaner1_concentrate", "cleaners34_tails"]
    )
    m.fs.final_tails_mixer = recovery_flowsheet.dry_solids_mixer(
        m.fs.properties, ["rougher_tails", "scavenger_tails"]
    )
    m.fs.regrind = recovery_flowsheet.StateJunction(property_package=m.fs.properties)

    for bank_name in BANK_NAMES:
        setattr(
            m.fs,
            bank_name,
            FlotationBank(
                property_package=m.fs.properties,
                recovery_basis="kinetic_cell_balance",
                number_of_cells=parameters["banks"][bank_name]["number_of_cells"],
            ),
        )

    m.fs.s_conditioning_to_rougher_mixer = Arc(
        source=m.fs.conditioning.outlet,
        destination=m.fs.rougher_mixer.fresh_feed,
    )
    m.fs.s_regrind_to_rougher_mixer = Arc(
        source=m.fs.regrind.outlet,
        destination=m.fs.rougher_mixer.regrind_recycle,
    )
    m.fs.s_rougher_mixer_to_rougher = Arc(
        source=m.fs.rougher_mixer.outlet,
        destination=m.fs.rougher.inlet,
    )
    m.fs.s_rougher_tails_to_final_tails = Arc(
        source=m.fs.rougher.tails,
        destination=m.fs.final_tails_mixer.rougher_tails,
    )
    m.fs.s_rougher_concentrate_to_cleaner1_mixer = Arc(
        source=m.fs.rougher.concentrate,
        destination=m.fs.cleaner1_mixer.rougher_concentrate,
    )
    m.fs.s_cleaner2_tails_to_cleaner1_mixer = Arc(
        source=m.fs.cleaner2.tails,
        destination=m.fs.cleaner1_mixer.cleaner2_tails,
    )
    m.fs.s_cleaner1_mixer_to_cleaner1 = Arc(
        source=m.fs.cleaner1_mixer.outlet,
        destination=m.fs.cleaner1.inlet,
    )
    m.fs.s_cleaner1_tails_to_scavenger = Arc(
        source=m.fs.cleaner1.tails,
        destination=m.fs.scavenger.inlet,
    )
    m.fs.s_scavenger_concentrate_to_regrind = Arc(
        source=m.fs.scavenger.concentrate,
        destination=m.fs.regrind.inlet,
    )
    m.fs.s_scavenger_tails_to_final_tails = Arc(
        source=m.fs.scavenger.tails,
        destination=m.fs.final_tails_mixer.scavenger_tails,
    )
    m.fs.s_cleaner1_concentrate_to_cleaner2_mixer = Arc(
        source=m.fs.cleaner1.concentrate,
        destination=m.fs.cleaner2_mixer.cleaner1_concentrate,
    )
    m.fs.s_cleaners34_tails_to_cleaner2_mixer = Arc(
        source=m.fs.cleaners34.tails,
        destination=m.fs.cleaner2_mixer.cleaners34_tails,
    )
    m.fs.s_cleaner2_mixer_to_cleaner2 = Arc(
        source=m.fs.cleaner2_mixer.outlet,
        destination=m.fs.cleaner2.inlet,
    )
    m.fs.s_cleaner2_concentrate_to_cleaners34 = Arc(
        source=m.fs.cleaner2.concentrate,
        destination=m.fs.cleaners34.inlet,
    )

    if expand_arcs:
        TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def _total_feed_scale(feed_scale):
    parameters = load_kinetic_cell_balance_parameters()
    return parameters["basis"]["scale_jsonbasis_to_plant_h"] * feed_scale


def _scaled_initial_streams(feed_scale=1.0):
    total_scale = _total_feed_scale(feed_scale)
    streams = {
        stream: {
            component: flow * total_scale for component, flow in component_flows.items()
        }
        for stream, component_flows in recovery_flowsheet.FEED_INPUTS.items()
    }
    streams.update(
        {
            stream: {
                component: flow * total_scale
                for component, flow in component_flows.items()
            }
            for stream, component_flows in recovery_flowsheet.SCENARIOS[
                TABLE1_PRODUCT_FIT_SCENARIO
            ]["initial_streams"].items()
        }
    )
    return streams


def _seed_deterministic_streams(model, feed_scale=1.0):
    streams = _scaled_initial_streams(feed_scale=feed_scale)
    for stream_name, ports in recovery_flowsheet.stream_port_map(model).items():
        for port in ports:
            recovery_flowsheet.set_port_values(port, 0, streams[stream_name])


def set_model_inputs(model, feed_scale=1.0):
    """Fix fresh feed, bank geometry, inventories, and rate constants."""
    parameters = load_kinetic_cell_balance_parameters()
    total_scale = parameters["basis"]["scale_jsonbasis_to_plant_h"] * feed_scale
    model.fs._feed_scale = feed_scale
    fresh_feed = {
        component: flow * total_scale
        for component, flow in recovery_flowsheet.FEED_INPUTS["fresh_feed"].items()
    }
    recovery_flowsheet.set_state_from_component_flows(
        model.fs.conditioning.properties, 0, fresh_feed
    )

    rho_water = parameters["bank_constants"]["rho_water_kg_per_m3"]
    for bank_name, bank_parameters in parameters["banks"].items():
        bank = getattr(model.fs, bank_name)
        bank.rho_water.set_value(rho_water)
        bank.rho_solid.set_value(bank_parameters["rho_solid"])
        bank.cell_volume.fix(bank_parameters["cell_volume"] * units.m**3)
        bank.volume_frac_air.fix(bank_parameters["air_holdup"])
        bank.pulp_solids_mass_fraction.fix(bank_parameters["pulp_solids_mass_fraction"])
        for time in model.fs.time:
            for component in COMPONENTS:
                bank.k_cb[time, component].fix(
                    bank_parameters["k_cb_per_hour"][component] / units.hour
                )


def initialize_model(
    model,
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Initialize the kinetic cell-balance flowsheet."""
    if initialization_method not in INITIALIZATION_METHODS:
        raise ValueError(
            f"Unknown initialization method {initialization_method!r}. "
            f"Expected one of {INITIALIZATION_METHODS}."
        )
    if (
        hasattr(model.fs, "_feed_scale")
        and abs(float(model.fs._feed_scale) - float(feed_scale)) > 1e-12
    ):
        raise ValueError(
            "initialize_model feed_scale must match the value used by "
            "set_model_inputs."
        )

    if initialization_method == DETERMINISTIC_STREAM_INITIALIZATION:
        _seed_deterministic_streams(model, feed_scale=feed_scale)

    model.fs._kinetic_cell_balance_initially_unfixed = {
        id(var)
        for var in model.component_data_objects(Var, descend_into=True)
        if not var.fixed
    }
    try:

        def initialize_unit(unit):
            if hasattr(unit, "default_initializer"):
                unit.default_initializer().initialize(unit)

        recovery_flowsheet.run_sequential_initialization(model, initialize_unit)
    finally:
        del model.fs._kinetic_cell_balance_initially_unfixed
    model.fs.initialization_method = initialization_method


def _scale_network_arc_flow_equalities(model):
    """Scale expanded network flow-equality constraints from current values."""
    # This scaling was needed to help convergence.
    scaled = 0
    for constraint in model.component_data_objects(
        Constraint,
        active=True,
        descend_into=True,
    ):
        if "_expanded.flow_mass_comp_equality" not in constraint.name:
            continue

        nominal = 1.0
        for var in identify_variables(constraint.body):
            var_value = value(var, exception=False)
            if var_value is not None and math.isfinite(var_value):
                nominal = max(nominal, abs(var_value))
        set_scaling_factor(constraint, 1 / nominal, overwrite=True)
        scaled += 1
    return scaled


def solve_model(model, solver=None, tee=False):
    """Solve the kinetic cell-balance flotation flowsheet."""
    _scale_network_arc_flow_equalities(model)
    opt = (
        get_solver("ipopt", solver_options=dict(FINAL_SOLVER_OPTIONS))
        if solver is None
        else solver
    )
    return opt.solve(model, tee=tee)


def build_and_initialize(
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Build, fix inputs, and initialize the kinetic cell-balance flowsheet."""
    model = build_model()
    set_model_inputs(model, feed_scale=feed_scale)
    initialize_model(
        model,
        feed_scale=feed_scale,
        initialization_method=initialization_method,
    )
    return model


def _safe_float(obj):
    """Evaluate ``obj`` and return a float, or ``None`` if unevaluable.

    Distinct from ``initializer._safe_value(obj, default=...)``: this variant
    deliberately returns ``None`` so that diagnostic helpers can surface
    missing values rather than substituting a default. Callers must handle
    the ``None`` case.
    """
    val = value(obj, exception=False)
    return None if val is None else float(val)


def _mixture_rho_solid(component_flows_kg_per_h, rho_mass_comp):
    """Return the inlet-composition-weighted dry-solid mixture density."""
    total_mass = sum(component_flows_kg_per_h.values())
    total_volume = sum(
        component_flows_kg_per_h[component] / rho_mass_comp[component]
        for component in COMPONENTS
    )
    return total_mass / total_volume


def cell_balance_summary(model, time=0):
    """Return per-bank cell-balance diagnostics and solvability audits."""
    summary = {}
    for bank_name in BANK_NAMES:
        bank = getattr(model.fs, bank_name)
        F_in_bank = {
            component: _safe_float(bank.inlet.flow_mass_comp[time, component])
            for component in COMPONENTS
        }
        M_total_by_cell = [
            _safe_float(bank.cell_total_solid_holdup[time, cell]) for cell in bank.cells
        ]
        k_hour = {
            component: _safe_float(bank.k_cb[time, component])
            for component in COMPONENTS
        }
        try:
            audit = audit_cascade_solvability(
                F_in_bank,
                M_total_by_cell,
                k_hour,
                bank_name=bank_name,
            )
        except (CascadeInputError, CascadeInfeasibleError) as exc:
            audit = {
                "bank": bank_name,
                "walk_completed": False,
                "first_infeasible_cell": getattr(exc, "cell_index", None),
                "first_failure_mode": getattr(exc, "failure_mode", type(exc).__name__),
            }

        current_rho_solid = _mixture_rho_solid(F_in_bank, SOLID_DENSITIES)
        reference_rho_solid = _safe_float(bank.rho_solid)
        if reference_rho_solid is None or reference_rho_solid == 0:
            rho_solid_deviation_pct = None
        else:
            rho_solid_deviation_pct = (
                100 * (current_rho_solid - reference_rho_solid) / reference_rho_solid
            )
        cells = []
        for cell in bank.cells:
            feed_total = sum(
                _safe_float(bank.cell_feed[time, cell, component])
                for component in COMPONENTS
            )
            float_total = sum(
                _safe_float(bank.cell_flotation_flow_mass_comp[time, cell, component])
                for component in COMPONENTS
            )
            pulp_out_total = sum(
                _safe_float(bank.cell_pulp_out_flow_mass_comp[time, cell, component])
                for component in COMPONENTS
            )
            cells.append(
                {
                    "cell_index": int(cell),
                    "feed_kg_per_h": feed_total,
                    "float_kg_per_h": float_total,
                    "pulp_out_kg_per_h": pulp_out_total,
                    "cell_mass_pull": (
                        None if feed_total <= 0 else float_total / feed_total
                    ),
                    "solid_holdup_kg": _safe_float(
                        bank.cell_total_solid_holdup[time, cell]
                    ),
                }
            )

        summary[bank_name] = {
            "number_of_cells": bank.config.number_of_cells,
            "rho_solid": reference_rho_solid,
            "rho_solid_current_inlet": current_rho_solid,
            "rho_solid_deviation_pct": rho_solid_deviation_pct,
            "tau_apparent_bank_h": _safe_float(bank.tau_apparent_bank[time]),
            "Q_slurry_in_m3_per_h": _safe_float(bank.Q_slurry_in[time]),
            "recovery": {
                component: _safe_float(bank.recovery[time, component])
                for component in COMPONENTS
            },
            "audit": audit,
            "cells": cells,
        }
    return summary


def report_results(model, time=0):
    """Return fixed-recovery product metrics plus cell-balance diagnostics."""
    report = recovery_flowsheet.report_results(model, time=time)
    report["cell_balance_summary"] = cell_balance_summary(model, time=time)
    return report


def print_results(model, time=0):
    """Print product metrics and per-bank cell-balance diagnostics."""
    recovery_flowsheet.print_results(model, time=time)
    print("\nKinetic cell-balance bank summary:")
    print("  Bank              N  residence (min)   audit walk")
    for bank_name, bank_summary in cell_balance_summary(model, time=time).items():
        print(
            f"  {bank_name:<14} "
            f"{bank_summary['number_of_cells']:>2} "
            f"{bank_summary['tau_apparent_bank_h'] * 60:>16.6f} "
            f"{str(bank_summary['audit']['walk_completed']):>12}"
        )


def mass_closure_residuals(model, time=0):
    """Return component dry-mass closure residuals for banks and mixers."""
    return recovery_flowsheet.mass_closure_residuals(model, time=time)


def main():
    """Build, initialize, solve, and print the kinetic cell-balance flowsheet example."""
    model = build_and_initialize(
        initialization_method=DETERMINISTIC_STREAM_INITIALIZATION
    )
    solve_model(model, tee=True)
    print_results(model)
    return model


if __name__ == "__main__":
    m = main()
