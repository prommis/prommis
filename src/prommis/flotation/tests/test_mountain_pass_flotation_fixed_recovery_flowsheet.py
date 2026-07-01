#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the Mountain Pass flotation flowsheet."""

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.flotation.bastnaesite_properties import COMPONENTS
from prommis.flotation.mountain_pass_flotation_fixed_recovery_flowsheet import (
    DETERMINISTIC_STREAM_INITIALIZATION,
    FIGURE2_MASS_BALANCE_SCENARIO,
    SCENARIOS,
    TABLE1_PRODUCT_FIT_SCENARIO,
    TABLE_1_PRODUCTS,
    TABLE2_REPORTED_SCENARIO,
    UNIT_MODEL_DEFAULT_INITIALIZATION,
    set_port_values,
    build_and_initialize,
    build_model,
    initialize_model,
    mass_closure_residuals,
    report_results,
    set_model_inputs,
    solve_model,
)

solver = get_solver("ipopt")


@pytest.fixture()
def model():
    return build_and_initialize()


@pytest.mark.component
@pytest.mark.build
def test_build_dof_units_and_diagnostics(model):
    assert model.fs.scenario == TABLE2_REPORTED_SCENARIO
    assert model.fs.initialization_method == UNIT_MODEL_DEFAULT_INITIALIZATION
    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model)
    DiagnosticsToolbox(model).assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.build
def test_recoveries_are_fixed_by_flowsheet_inputs():
    model = build_model()
    for bank_name in SCENARIOS[TABLE2_REPORTED_SCENARIO]["recoveries"]:
        assert not getattr(model.fs, bank_name).recovery[0, "REO"].fixed

    set_model_inputs(model)
    for scenario_recoveries in SCENARIOS[TABLE2_REPORTED_SCENARIO][
        "recoveries"
    ].values():
        assert set(scenario_recoveries) == set(COMPONENTS)
    for bank_name in SCENARIOS[TABLE2_REPORTED_SCENARIO]["recoveries"]:
        assert getattr(model.fs, bank_name).recovery[0, "REO"].fixed
    assert degrees_of_freedom(model) == 0


@pytest.mark.component
@pytest.mark.build
def test_deterministic_stream_initialization_option():
    model = build_and_initialize(
        initialization_method=DETERMINISTIC_STREAM_INITIALIZATION
    )
    expected = SCENARIOS[TABLE2_REPORTED_SCENARIO]["initial_streams"]["rougher_feed"][
        "REO"
    ]

    assert model.fs.initialization_method == DETERMINISTIC_STREAM_INITIALIZATION
    assert value(model.fs.rougher.inlet.flow_mass_comp[0, "REO"]) == pytest.approx(
        expected
    )
    assert degrees_of_freedom(model) == 0


@pytest.mark.component
@pytest.mark.build
def test_invalid_initialization_method():
    model = build_model()
    set_model_inputs(model)

    with pytest.raises(ValueError, match="Unknown initialization method"):
        initialize_model(model, initialization_method="unknown_method")


@pytest.mark.unit
def testset_port_values_rejects_incomplete_component_data():
    model = build_model(expand_arcs=False)
    component_flows = {component: 1.0 for component in COMPONENTS}
    component_flows.pop("inert_gangue")

    with pytest.raises(ValueError, match="missing=.*inert_gangue"):
        set_port_values(model.fs.rougher.inlet, 0, component_flows)


@pytest.mark.unit
def testset_port_values_rejects_unexpected_component_data():
    model = build_model(expand_arcs=False)
    component_flows = {component: 1.0 for component in COMPONENTS}
    component_flows["unexpected"] = 1.0

    with pytest.raises(ValueError, match="unexpected=.*unexpected"):
        set_port_values(model.fs.rougher.inlet, 0, component_flows)


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_solve_and_mass_closure(model):
    results = solve_model(model, solver=solver)
    assert_optimal_termination(results)
    DiagnosticsToolbox(model).assert_no_numerical_warnings()
    for component_residuals in mass_closure_residuals(model).values():
        for component in COMPONENTS:
            assert component_residuals[component] == pytest.approx(0, abs=1e-6)


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_table2_reported_scenario_solves_as_diagnostic_reference():
    table2_model = build_and_initialize(scenario=TABLE2_REPORTED_SCENARIO)
    results = solve_model(table2_model, solver=solver)
    assert_optimal_termination(results)

    recoveries = SCENARIOS[TABLE2_REPORTED_SCENARIO]["recoveries"]
    for bank_name, bank_recoveries in recoveries.items():
        bank = getattr(table2_model.fs, bank_name)
        for component in ("REO", "CaO", "BaO", "SrO"):
            assert value(bank.recovery[0, component]) == pytest.approx(
                bank_recoveries[component]
            )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_product_metrics_and_paper_gate(model):
    results = solve_model(model, solver=solver)
    assert_optimal_termination(results)
    report = report_results(model)
    metrics = report["product_metrics"]

    assert metrics["overall_REO_recovery"] == pytest.approx(
        TABLE_1_PRODUCTS["final_concentrate"]["reo_recovery"], abs=1.0
    )
    assert metrics["final_tails_REO_grade"] == pytest.approx(
        TABLE_1_PRODUCTS["combined_tails"]["assays"]["REO"], abs=0.2
    )

    if report["table1_products_within_tolerance"]:
        assert metrics["final_concentrate_yield"] == pytest.approx(9.6, abs=0.2)
        assert 64.6 <= metrics["final_concentrate_REO_grade"] <= 65.6
    else:
        assert report["table1_products_within_tolerance"] is False


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_table1_product_fit_example_matches_paper_envelope():
    model = build_and_initialize(scenario=TABLE1_PRODUCT_FIT_SCENARIO)
    results = solve_model(model, solver=solver)
    assert_optimal_termination(results)
    report = report_results(model)
    metrics = report["product_metrics"]

    assert report["table1_products_within_tolerance"]
    assert metrics["final_concentrate_yield"] == pytest.approx(9.6, abs=0.2)
    assert 64.6 <= metrics["final_concentrate_REO_grade"] <= 65.6
    assert metrics["final_concentrate_CaO_grade"] == pytest.approx(2.7, abs=0.5)
    assert metrics["final_concentrate_BaO_grade"] == pytest.approx(0.9, abs=0.5)
    assert metrics["final_concentrate_SrO_grade"] == pytest.approx(5.4, abs=0.5)


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_figure2_mass_balance_scenario_solves_as_diagnostic_reference():
    model = build_and_initialize(scenario=FIGURE2_MASS_BALANCE_SCENARIO)
    results = solve_model(model, solver=solver)
    assert_optimal_termination(results)
    report = report_results(model)

    assert report["scenario"] == FIGURE2_MASS_BALANCE_SCENARIO


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_feed_scale_preserves_recovery_and_grades():
    base = build_and_initialize(feed_scale=1.0)
    scaled = build_and_initialize(feed_scale=2.0)
    assert_optimal_termination(solve_model(base, solver=solver))
    assert_optimal_termination(solve_model(scaled, solver=solver))

    base_metrics = report_results(base)["product_metrics"]
    scaled_metrics = report_results(scaled)["product_metrics"]
    for metric in base_metrics:
        assert scaled_metrics[metric] == pytest.approx(base_metrics[metric])
