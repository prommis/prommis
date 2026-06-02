#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Mountain Pass bastnaesite flotation flowsheet with fixed component recoveries.

This module implements the five-bank Mountain Pass flotation circuit (rougher,
scavenger, cleaner 1, cleaner 2, cleaners 3&4) as a dry-solids recycle
flowsheet where each bank's per-component recovery is a user-fixed input.

Reference: Pradip and Fuerstenau (2013), "Design and development of novel
flotation reagents for the beneficiation of Mountain Pass rare-earth ore,"
Minerals & Metallurgical Processing 30(1), 1-9.

Scenarios
---------
Three recovery scenarios are provided in ``mountain_pass_scenario_data.json``,
each calibrated from a different subset of the published data from
Pradip and Fuerstenau (2013):

``table2_reported_with_fitted_inert``
    Published Table 2 recoveries for REO, CaO, BaO, and SrO, with
    inert_gangue recovery fitted to each bank's reported mass pull. This is
    the default scenario. It reproduces Table 2 as closely as the available
    data allow, but inert_gangue recoveries are inferred because Pradip &
    Fuerstenau did not report them.

``figure2_mass_balance``
    Recoveries back-calculated from Figure 2 stream masses and assays for
    the connected recycle circuit. This scenario uses the graphical data to
    close the circuit mass balance. It may differ slightly from Table 2
    because the mass balance data presented in Figure 2 were reconciled by
    the authors, and Figure 2 represents a different snapshot of the plant
    than Table 2.

``table1_product_fit``
    Keeps Table 2 assayed-component recoveries fixed, but recalibrates the
    unmeasured inert_gangue recoveries in Cleaner 2 and Cleaners 3&4 so the
    example flowsheet reproduces the Table 1 final-product envelope (overall
    REO recovery, concentrate REO grade, and concentrate mass pull).
"""

import json
from pathlib import Path

from pyomo.environ import ConcreteModel, Param, TransformationFactory, Var, units, value
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Mixer, StateJunction
from idaes.models.unit_models.mixer import MixingType, MomentumMixingType

from prommis.flotation.bastnaesite_properties import (
    BastnaesiteParameters,
    COMPONENTS,
)
from prommis.flotation.flotation_bank import FlotationBank

__author__ = "Daison Yancy Caballero"

SCENARIO_DATA_FILE = Path(__file__).with_name("mountain_pass_scenario_data.json")
TABLE2_REPORTED_SCENARIO = "table2_reported_with_fitted_inert"
FIGURE2_MASS_BALANCE_SCENARIO = "figure2_mass_balance"
TABLE1_PRODUCT_FIT_SCENARIO = "table1_product_fit"
UNIT_MODEL_DEFAULT_INITIALIZATION = "unit_model_default"
DETERMINISTIC_STREAM_INITIALIZATION = "deterministic_streams"
INITIALIZATION_METHODS = (
    UNIT_MODEL_DEFAULT_INITIALIZATION,
    DETERMINISTIC_STREAM_INITIALIZATION,
)


def stream_total(component_flows):
    return sum(component_flows[j] for j in COMPONENTS)


def stream_assay(component_flows, component):
    return 100 * component_flows[component] / stream_total(component_flows)


def load_scenario_data(data_file=SCENARIO_DATA_FILE):
    """Read the offline-generated Mountain Pass scenario data."""
    path = Path(data_file)
    with path.open(encoding="utf-8") as file:
        data = json.load(file)
    feed_components = tuple(data["feed"]["fresh_feed"])
    if set(feed_components) != set(COMPONENTS):
        raise ValueError(
            "Scenario data feed components do not match the property package."
        )
    for scenario in data["scenarios"].values():
        for bank_recoveries in scenario["recoveries"].values():
            if set(bank_recoveries) != set(COMPONENTS):
                raise ValueError(
                    "Scenario data recovery components do not match the property "
                    "package."
                )
    return data


SCENARIO_DATA = load_scenario_data()
FEED_INPUTS = SCENARIO_DATA["feed"]
TABLE_1_PRODUCTS = SCENARIO_DATA["table_1_products"]
SCENARIOS = SCENARIO_DATA["scenarios"]
BANK_NAMES = tuple(SCENARIOS[TABLE2_REPORTED_SCENARIO]["recoveries"])


def _scenario_record(scenario):
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario {scenario!r}.")
    return SCENARIOS[scenario]


def dry_solids_mixer(property_package, inlet_list):
    return Mixer(
        property_package=property_package,
        inlet_list=inlet_list,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_mixing_type=MomentumMixingType.none,
        energy_mixing_type=MixingType.none,
        has_holdup=False,
    )


def add_conditioning_metadata(fs):
    fs.pH_rougher_ref = Param(initialize=8.8, mutable=True, units=units.dimensionless)
    fs.pH_cleaner_ref = Param(initialize=8.8, mutable=True, units=units.dimensionless)
    fs.temperature_ref = Param(initialize=343.15, mutable=True, units=units.K)
    fs.soda_ash_dose_ref = Param(
        initialize=2.9,
        mutable=True,
        units=units.kg / units.tonne,
    )
    fs.lignin_sulfonate_dose_ref = Param(
        initialize=2.9,
        mutable=True,
        units=units.kg / units.tonne,
    )
    fs.fluosilicate_dose_ref = Param(
        initialize=0.4,
        mutable=True,
        units=units.kg / units.tonne,
    )
    fs.collector_dose_ref = Param(
        initialize=0.3,
        mutable=True,
        units=units.kg / units.tonne,
    )
    fs.frother_dose_ref = Param(
        initialize=0.0,
        mutable=True,
        units=units.kg / units.tonne,
    )


def build_model(scenario=TABLE2_REPORTED_SCENARIO, expand_arcs=True):
    """Build the Mountain Pass dry-solids flotation flowsheet."""
    _scenario_record(scenario)

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    m.fs.scenario = scenario
    add_conditioning_metadata(m.fs)

    m.fs.conditioning = StateJunction(property_package=m.fs.properties)
    m.fs.rougher_mixer = dry_solids_mixer(
        m.fs.properties, ["fresh_feed", "regrind_recycle"]
    )
    m.fs.cleaner1_mixer = dry_solids_mixer(
        m.fs.properties, ["rougher_concentrate", "cleaner2_tails"]
    )
    m.fs.cleaner2_mixer = dry_solids_mixer(
        m.fs.properties, ["cleaner1_concentrate", "cleaners34_tails"]
    )
    m.fs.final_tails_mixer = dry_solids_mixer(
        m.fs.properties, ["rougher_tails", "scavenger_tails"]
    )
    m.fs.regrind = StateJunction(property_package=m.fs.properties)

    for bank_name in BANK_NAMES:
        setattr(
            m.fs,
            bank_name,
            FlotationBank(
                property_package=m.fs.properties,
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


def set_state_from_component_flows(state_block, time, component_flows, fix=True):
    """Set or fix a state block's component mass flows."""
    state = state_block[time] if hasattr(state_block, "__getitem__") else state_block
    for component, flow in component_flows.items():
        var = state.flow_mass_comp[component]
        var.set_value(flow)
        if fix:
            var.fix(flow)


def set_port_values(port, time, component_flows):
    components = set(component_flows)
    expected = set(COMPONENTS)
    if components != expected:
        missing = sorted(expected - components)
        unexpected = sorted(components - expected)
        raise ValueError(
            "Stream initialization data must include exactly the bastnaesite "
            f"components; missing={missing}, unexpected={unexpected}."
        )
    for component, flow in component_flows.items():
        port.flow_mass_comp[time, component].set_value(flow)


def fix_bank_recoveries(model):
    """Fix bank recoveries from the selected scenario data."""
    recoveries = _scenario_record(model.fs.scenario)["recoveries"]
    for bank_name, bank_recoveries in recoveries.items():
        bank = getattr(model.fs, bank_name)
        for time in model.fs.time:
            for component in COMPONENTS:
                bank.recovery[time, component].fix(bank_recoveries[component])


def set_model_inputs(model, feed_scale=1.0):
    """Fix the fresh-feed dry component flows and scenario recoveries."""
    fresh_feed = {
        component: flow * feed_scale
        for component, flow in FEED_INPUTS["fresh_feed"].items()
    }
    set_state_from_component_flows(model.fs.conditioning.properties, 0, fresh_feed)
    fix_bank_recoveries(model)


def _scaled_initial_streams(model, feed_scale=1.0):
    streams = {
        stream: {
            component: flow * feed_scale for component, flow in component_flows.items()
        }
        for stream, component_flows in FEED_INPUTS.items()
    }
    streams.update(
        {
            stream: {
                component: flow * feed_scale
                for component, flow in component_flows.items()
            }
            for stream, component_flows in _scenario_record(model.fs.scenario)[
                "initial_streams"
            ].items()
        }
    )
    return streams


def stream_port_map(model):
    fs = model.fs
    # StateJunction inlet/outlet ports intentionally alias the same state vars.
    # Listing both ports keeps the named stream map complete for initialization
    # and reporting without adding pass-through constraints.
    return {
        "fresh_feed": [fs.conditioning.inlet, fs.conditioning.outlet],
        "rougher_feed": [fs.rougher_mixer.outlet, fs.rougher.inlet],
        "rougher_concentrate": [
            fs.rougher.concentrate,
            fs.cleaner1_mixer.rougher_concentrate,
        ],
        "rougher_tails": [fs.rougher.tails, fs.final_tails_mixer.rougher_tails],
        "cleaner1_feed": [fs.cleaner1_mixer.outlet, fs.cleaner1.inlet],
        "cleaner1_concentrate": [
            fs.cleaner1.concentrate,
            fs.cleaner2_mixer.cleaner1_concentrate,
        ],
        "cleaner1_tails": [fs.cleaner1.tails, fs.scavenger.inlet],
        "scavenger_concentrate": [
            fs.scavenger.concentrate,
            fs.regrind.inlet,
            fs.regrind.outlet,
            fs.rougher_mixer.regrind_recycle,
        ],
        "scavenger_tails": [
            fs.scavenger.tails,
            fs.final_tails_mixer.scavenger_tails,
        ],
        "cleaner2_feed": [fs.cleaner2_mixer.outlet, fs.cleaner2.inlet],
        "cleaner2_concentrate": [fs.cleaner2.concentrate, fs.cleaners34.inlet],
        "cleaner2_tails": [fs.cleaner2.tails, fs.cleaner1_mixer.cleaner2_tails],
        "cleaners34_concentrate": [fs.cleaners34.concentrate],
        "cleaners34_tails": [
            fs.cleaners34.tails,
            fs.cleaner2_mixer.cleaners34_tails,
        ],
        "final_tails": [fs.final_tails_mixer.outlet],
    }


def run_sequential_initialization(model, initialize_unit):
    """Run sequential decomposition without retaining temporary tear fixings."""
    initially_unfixed = [
        var
        for var in model.component_data_objects(Var, descend_into=True)
        if not var.fixed
    ]
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.solve_tears = True
    seq.options.tear_method = "Wegstein"
    seq.options.tol_type = "rel"
    seq.run(model, initialize_unit)

    for var in initially_unfixed:
        if var.fixed:
            var.unfix()


def _initialize_with_unit_model_defaults(model):
    def initialize_unit(unit):
        if hasattr(unit, "default_initializer"):
            unit.default_initializer().initialize(unit)

    run_sequential_initialization(model, initialize_unit)


def _initialize_with_deterministic_streams(model, feed_scale=1.0):
    """Set deterministic initial guesses from the connected dry-solids solve."""
    streams = _scaled_initial_streams(model, feed_scale=feed_scale)
    port_map = stream_port_map(model)
    for stream_name, ports in port_map.items():
        for port in ports:
            set_port_values(port, 0, streams[stream_name])


def initialize_model(
    model,
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Initialize the flowsheet using the selected initialization method."""
    if initialization_method == UNIT_MODEL_DEFAULT_INITIALIZATION:
        _initialize_with_unit_model_defaults(model)
    elif initialization_method == DETERMINISTIC_STREAM_INITIALIZATION:
        _initialize_with_deterministic_streams(model, feed_scale=feed_scale)
    else:
        raise ValueError(
            f"Unknown initialization method {initialization_method!r}. "
            f"Expected one of {INITIALIZATION_METHODS}."
        )
    model.fs.initialization_method = initialization_method


def solve_model(model, solver=None, tee=False):
    """Solve the flowsheet."""
    opt = get_solver("ipopt") if solver is None else solver
    return opt.solve(model, tee=tee)


def build_and_initialize(
    scenario=TABLE2_REPORTED_SCENARIO,
    feed_scale=1.0,
    initialization_method=UNIT_MODEL_DEFAULT_INITIALIZATION,
):
    """Build the flowsheet, fix scenario inputs, and set initial values."""
    model = build_model(scenario=scenario)
    set_model_inputs(model, feed_scale=feed_scale)
    initialize_model(
        model,
        feed_scale=feed_scale,
        initialization_method=initialization_method,
    )
    return model


def stream_component_flows(model, stream_name, time=0):
    """Read component flows from a named flowsheet stream."""
    port = stream_port_map(model)[stream_name][0]
    return {
        component: value(port.flow_mass_comp[time, component])
        for component in COMPONENTS
    }


def table1_products_within_tolerance(product_metrics):
    """Return whether product metrics match the reported Table 1 envelope."""
    return (
        abs(product_metrics["final_concentrate_yield"] - 9.6) <= 0.2
        and 64.6 <= product_metrics["final_concentrate_REO_grade"] <= 65.6
        and abs(product_metrics["final_concentrate_CaO_grade"] - 2.7) <= 0.5
        and abs(product_metrics["final_concentrate_BaO_grade"] - 0.9) <= 0.5
        and abs(product_metrics["final_concentrate_SrO_grade"] - 5.4) <= 0.5
        and abs(product_metrics["final_tails_REO_grade"] - 1.7) <= 0.2
        and abs(product_metrics["overall_REO_recovery"] - 80.1) <= 1.0
    )


def report_results(model, time=0):
    """Return model stream totals, assays, and product metrics."""
    streams = {}
    for stream_name in stream_port_map(model):
        component_flows = stream_component_flows(model, stream_name, time=time)
        streams[stream_name] = {
            "component_flows": component_flows,
            "mass": stream_total(component_flows),
            "assays": {
                component: stream_assay(component_flows, component)
                for component in COMPONENTS
            },
        }

    final_concentrate = streams["cleaners34_concentrate"]["component_flows"]
    final_tails = streams["final_tails"]["component_flows"]
    fresh_feed = streams["fresh_feed"]["component_flows"]
    product_metrics = {
        "final_concentrate_yield": 100
        * stream_total(final_concentrate)
        / stream_total(fresh_feed),
        "final_concentrate_REO_grade": stream_assay(final_concentrate, "REO"),
        "final_concentrate_CaO_grade": stream_assay(final_concentrate, "CaO"),
        "final_concentrate_BaO_grade": stream_assay(final_concentrate, "BaO"),
        "final_concentrate_SrO_grade": stream_assay(final_concentrate, "SrO"),
        "final_tails_REO_grade": stream_assay(final_tails, "REO"),
        "overall_REO_recovery": 100 * final_concentrate["REO"] / fresh_feed["REO"],
    }
    return {
        "scenario": model.fs.scenario,
        "scenario_description": _scenario_record(model.fs.scenario)["description"],
        "streams": streams,
        "product_metrics": product_metrics,
        "paper_table_1": TABLE_1_PRODUCTS,
        "table1_products_within_tolerance": table1_products_within_tolerance(
            product_metrics
        ),
    }


def print_results(model, time=0):
    """Print the flowsheet product metrics and paper comparison."""
    report = report_results(model, time=time)
    metrics = report["product_metrics"]
    paper = report["paper_table_1"]
    paper_concentrate = paper["final_concentrate"]
    paper_tails = paper["combined_tails"]

    print("\nMountain Pass bastnaesite flotation results")
    print(f"Scenario: {report['scenario']}")
    print(f"Scenario note: {report['scenario_description']}")
    print("\nFinal product metrics:")
    print(
        "  Final concentrate yield: "
        f"{metrics['final_concentrate_yield']:.3f} wt% "
        f"(paper: {paper_concentrate['weight_split']:.3f} wt%)"
    )
    print(
        "  Final concentrate REO grade: "
        f"{metrics['final_concentrate_REO_grade']:.3f} wt% "
        f"(paper: {paper_concentrate['assays']['REO']:.3f} wt%; "
        "accepted range: 64.6-65.6 wt%)"
    )
    print(
        "  Final concentrate CaO grade: "
        f"{metrics['final_concentrate_CaO_grade']:.3f} wt% "
        f"(paper: {paper_concentrate['assays']['CaO']:.3f} wt%)"
    )
    print(
        "  Final concentrate BaO grade: "
        f"{metrics['final_concentrate_BaO_grade']:.3f} wt% "
        f"(paper: {paper_concentrate['assays']['BaO']:.3f} wt%)"
    )
    print(
        "  Final concentrate SrO grade: "
        f"{metrics['final_concentrate_SrO_grade']:.3f} wt% "
        f"(paper: {paper_concentrate['assays']['SrO']:.3f} wt%)"
    )
    print(
        "  Final tails REO grade: "
        f"{metrics['final_tails_REO_grade']:.3f} wt% "
        f"(paper: {paper_tails['assays']['REO']:.3f} wt%)"
    )
    print(
        "  Overall REO recovery: "
        f"{metrics['overall_REO_recovery']:.3f}% "
        f"(paper: {paper_concentrate['reo_recovery']:.3f}%)"
    )
    print(
        "\nWithin Table 1 product tolerances: "
        f"{report['table1_products_within_tolerance']}"
    )


def mass_closure_residuals(model, time=0):
    """Return component dry-mass closure residuals for banks and mixers."""
    fs = model.fs
    residuals = {}

    for bank_name in BANK_NAMES:
        bank = getattr(fs, bank_name)
        residuals[bank_name] = {}
        for component in COMPONENTS:
            residuals[bank_name][component] = value(
                bank.properties_in[time].flow_mass_comp[component]
                - bank.properties_concentrate[time].flow_mass_comp[component]
                - bank.properties_tails[time].flow_mass_comp[component]
            )

    mixer_specs = {
        "rougher_mixer": ("fresh_feed", "regrind_recycle"),
        "cleaner1_mixer": ("rougher_concentrate", "cleaner2_tails"),
        "cleaner2_mixer": ("cleaner1_concentrate", "cleaners34_tails"),
        "final_tails_mixer": ("rougher_tails", "scavenger_tails"),
    }
    for mixer_name, inlet_names in mixer_specs.items():
        mixer = getattr(fs, mixer_name)
        residuals[mixer_name] = {}
        for component in COMPONENTS:
            residuals[mixer_name][component] = value(
                sum(
                    getattr(mixer, inlet).flow_mass_comp[time, component]
                    for inlet in inlet_names
                )
                - mixer.outlet.flow_mass_comp[time, component]
            )

    return residuals


def main():
    """Build, initialize, solve, and print the example flowsheet."""
    # Override the API default scenario (Table 2) with the Table 1 fit scenario:
    # it is the only scenario calibrated to reproduce the paper's headline product
    # envelope, so it is the one that satisfies "table1_products_within_tolerance".
    model = build_and_initialize(scenario=TABLE1_PRODUCT_FIT_SCENARIO)
    solve_model(model, tee=True)
    print_results(model)

    return model


if __name__ == "__main__":
    m = main()
