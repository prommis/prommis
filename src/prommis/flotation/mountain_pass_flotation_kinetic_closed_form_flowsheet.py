####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Mountain Pass bastnaesite flotation flowsheet with closed-form kinetic recoveries.

This module implements the five-bank Mountain Pass flotation circuit using the
closed-form tanks-in-series kinetic recovery model
(``recovery_basis="kinetic_closed_form"``). Each bank's per-component recovery
is constrained by the tanks-in-series form of the Garcia-Zuniga first-order
flotation rate law:

    R = R_inf * (1 - (1 + k_cf * tau / N)^(-N))

where ``k_cf`` is the apparent first-order flotation rate constant calibrated
against a flow-proportional cell inventory, ``tau`` is the apparent slurry
residence time inferred from bank geometry and inlet pulp density, and ``N``
is the number of perfectly mixed cells in series

The fitted rate constants are loaded from
``mountain_pass_kinetic_closed_form_parameters.json`` and are tied to the
``table1_product_fit`` scenario basis, bank-local pulp-density assumptions,
air holdup, and cell geometry.

References:

    Polat, M. and Chander, S. (2000). "First-order flotation kinetics models
    and methods for estimation of the true distribution of flotation rate
    constants." International Journal of Mineral Processing, 58(1-4), 145-166.
    (Garcia-Zuniga first-order model and its tanks-in-series extension.)

    Pradip and Fuerstenau (2013), "Design and development of novel flotation
    reagents for the beneficiation of Mountain Pass rare-earth ore," Minerals
    & Metallurgical Processing 30(1), 1-9.
"""

import json
from pathlib import Path

from pyomo.environ import ConcreteModel, TransformationFactory, units, value
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation.bastnaesite_properties import (
    BastnaesiteParameters,
    COMPONENTS,
)
from prommis.flotation.flotation_bank import FlotationBank

__author__ = "Daison Yancy Caballero"

KINETIC_PARAMETERS_FILE = Path(__file__).with_name(
    "mountain_pass_kinetic_closed_form_parameters.json"
)
TABLE1_PRODUCT_FIT_SCENARIO = recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
UNIT_MODEL_DEFAULT_INITIALIZATION = recovery_flowsheet.UNIT_MODEL_DEFAULT_INITIALIZATION
DETERMINISTIC_STREAM_INITIALIZATION = (
    recovery_flowsheet.DETERMINISTIC_STREAM_INITIALIZATION
)
INITIALIZATION_METHODS = recovery_flowsheet.INITIALIZATION_METHODS
BANK_NAMES = recovery_flowsheet.BANK_NAMES


def load_kinetic_parameters(path=None):
    """Load the preprocessed kinetic-parameter JSON."""
    parameter_file = KINETIC_PARAMETERS_FILE if path is None else Path(path)
    with parameter_file.open(encoding="utf-8") as file:
        return json.load(file)


def build_model(expand_arcs=True):
    """Build the Mountain Pass flowsheet with kinetic flotation banks."""
    kinetic_parameters = load_kinetic_parameters()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    for component, density in kinetic_parameters["property_package"][
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
                recovery_basis="kinetic_closed_form",
                number_of_cells=kinetic_parameters["banks"][bank_name][
                    "number_of_cells"
                ],
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
    kinetic_parameters = load_kinetic_parameters()
    return kinetic_parameters["basis"]["scale_jsonbasis_to_plant_h"] * feed_scale


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
    """Fix plant-basis fresh feed, bank geometry, and kinetic parameters."""
    kinetic_parameters = load_kinetic_parameters()
    total_scale = kinetic_parameters["basis"]["scale_jsonbasis_to_plant_h"] * feed_scale
    model.fs._feed_scale = feed_scale
    fresh_feed = {
        component: flow * total_scale
        for component, flow in recovery_flowsheet.FEED_INPUTS["fresh_feed"].items()
    }
    recovery_flowsheet.set_state_from_component_flows(
        model.fs.conditioning.properties, 0, fresh_feed
    )

    rho_water = kinetic_parameters["bank_constants"]["rho_water_kg_per_m3"]
    for bank_name, bank_parameters in kinetic_parameters["banks"].items():
        bank = getattr(model.fs, bank_name)
        bank.rho_water.set_value(rho_water)
        bank.cell_volume.fix(bank_parameters["cell_volume"] * units.m**3)
        bank.air_holdup.fix(bank_parameters["air_holdup"])
        bank.pulp_solids_mass_fraction.fix(bank_parameters["pulp_solids_mass_fraction"])
        for time in model.fs.time:
            for component in COMPONENTS:
                bank.k_cf[time, component].fix(
                    bank_parameters["k_cf_per_min"][component] * 60.0 / units.hour
                )


def initialize_model(
    model,
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Initialize the kinetic flowsheet using the selected initialization method."""
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

    def initialize_unit(unit):
        if hasattr(unit, "default_initializer"):
            unit.default_initializer().initialize(unit)

    recovery_flowsheet.run_sequential_initialization(model, initialize_unit)
    model.fs.initialization_method = initialization_method


def solve_model(model, solver=None, tee=False):
    """Solve the kinetic flotation flowsheet."""
    opt = get_solver("ipopt") if solver is None else solver
    return opt.solve(model, tee=tee)


def build_and_initialize(
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Build, fix inputs, and initialize the kinetic flowsheet."""
    model = build_model()
    set_model_inputs(model, feed_scale=feed_scale)
    initialize_model(
        model,
        feed_scale=feed_scale,
        initialization_method=initialization_method,
    )
    return model


def kinetic_summary(model, time=0):
    """Return per-bank kinetic operating and fitted-rate diagnostics."""
    summary = {}
    for bank_name in BANK_NAMES:
        bank = getattr(model.fs, bank_name)
        summary[bank_name] = {
            "number_of_cells": bank.config.number_of_cells,
            "tau_min": value(bank.tau[time]) * 60.0,
            "Q_slurry_m3_per_min": value(bank.flow_vol_slurry[time]) / 60.0,
            "k_cf_per_min": {
                component: value(bank.k_cf[time, component]) / 60.0
                for component in COMPONENTS
            },
            "recovery": {
                component: value(bank.recovery[time, component])
                for component in COMPONENTS
            },
        }
    return summary


def report_results(model, time=0):
    """Return fixed-recovery product metrics plus kinetic diagnostics."""
    report = recovery_flowsheet.report_results(model, time=time)
    report["kinetic_summary"] = kinetic_summary(model, time=time)
    return report


def print_results(model, time=0):
    """Print product metrics and per-bank kinetic diagnostics."""
    recovery_flowsheet.print_results(model, time=time)
    print("\nKinetic bank summary:")
    print("  Bank              N      tau (min)      k_REO (1/min)")
    for bank_name, bank_summary in kinetic_summary(model, time=time).items():
        print(
            f"  {bank_name:<14} "
            f"{bank_summary['number_of_cells']:>2} "
            f"{bank_summary['tau_min']:>14.6f} "
            f"{bank_summary['k_cf_per_min']['REO']:>16.6f}"
        )


def mass_closure_residuals(model, time=0):
    """Return component dry-mass closure residuals for banks and mixers."""
    return recovery_flowsheet.mass_closure_residuals(model, time=time)


def main():
    """Build, initialize, solve, and print the kinetic flowsheet example."""
    model = build_and_initialize()
    solve_model(model, tee=True)
    print_results(model)
    return model


if __name__ == "__main__":
    m = main()
