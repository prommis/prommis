#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from unittest.mock import MagicMock, patch

from pyomo.opt import SolverStatus, TerminationCondition

import pytest

from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice

# Import the module using the correct path for your project structure
from prommis.superstructure.report_superstructure_results import (
    report_superstructure_costing,
    report_superstructure_environmental_impacts,
    report_superstructure_results_overview,
    report_superstructure_streams,
)

# --- Fixtures for common mock objects ---


@pytest.fixture
def mock_model():
    m = MagicMock()
    m.fs = MagicMock()
    m.fs.all_opts_set = ["A", "B"]
    m.fs.option_binary_var = {"A": MagicMock(), "B": MagicMock()}
    m.fs.objective_function_choice = MagicMock()
    m.fs.costing = MagicMock()
    m.fs.environmental_impacts = MagicMock()
    m.fs.plant_life_range = [2025, 2026]
    m.fs.operational_range = [2025, 2026]
    m.fs.tracked_comps = ["comp1"]
    m.fs.f_stages = ["stage1"]

    # Complete f dictionary for all combinations
    m.fs.f = {
        ("stage1", "comp1", 2025): MagicMock(),
        ("stage1", "comp1", 2026): MagicMock(),
    }

    # Complete f_in dictionary for all combinations
    m.fs.f_in = {
        ("A", "comp1", 2025): MagicMock(),
        ("A", "comp1", 2026): MagicMock(),
        ("B", "comp1", 2025): MagicMock(),
        ("B", "comp1", 2026): MagicMock(),
    }

    # Complete f_out dictionary for all combinations
    m.fs.f_out = {
        ("A", "comp1", 2025): MagicMock(),
        ("A", "comp1", 2026): MagicMock(),
        ("B", "comp1", 2025): MagicMock(),
        ("B", "comp1", 2026): MagicMock(),
    }

    # Byproduct valorization
    m.fs.byproduct_valorization = MagicMock()
    m.fs.byproduct_valorization.byproducts_set = ["byprod1"]
    m.fs.byproduct_valorization.byproduct_produced = {
        ("byprod1", 2025): MagicMock(),
        ("byprod1", 2026): MagicMock(),
    }

    # Environmental impacts - complete dictionaries
    m.fs.environmental_impacts.total_impacts = MagicMock()
    m.fs.environmental_impacts.epsilon = MagicMock()
    m.fs.environmental_impacts.total_yearly_impacts = {
        2025: MagicMock(),
        2026: MagicMock(),
    }
    m.fs.environmental_impacts.option_yearly_impacts = {
        ("A", 2025): MagicMock(),
        ("A", 2026): MagicMock(),
        ("B", 2025): MagicMock(),
        ("B", 2026): MagicMock(),
    }

    # Costing - complete dictionaries
    m.fs.costing.net_present_value = MagicMock()
    m.fs.costing.cost_of_recovery = MagicMock()
    m.fs.costing.discount_factor = MagicMock()
    m.fs.costing.cash_flow = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.i_operating_expense_escalation = MagicMock()
    m.fs.costing.i_capital_escalation = MagicMock()
    m.fs.costing.main_product_profit = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.total_profit = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.total_byproduct_profit = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.equipment_cost = {"A": MagicMock(), "B": MagicMock()}
    m.fs.costing.total_plant_cost = MagicMock()
    m.fs.costing.lang_factor = MagicMock()
    m.fs.costing.total_overnight_cost = MagicMock()
    m.fs.costing.financing_factor = MagicMock()
    m.fs.costing.other_costs_factor = MagicMock()
    m.fs.costing.total_overnight_cost_expended = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.total_operators = MagicMock()
    m.fs.costing.aggregate_variable_operating_cost = {
        2025: MagicMock(),
        2026: MagicMock(),
    }
    m.fs.costing.aggregate_fixed_operating_cost = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.total_operating_expense = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.cost_of_labor = MagicMock()
    m.fs.costing.m_and_sm = MagicMock()
    m.fs.costing.sa_and_qa_qc = MagicMock()
    m.fs.costing.s_ip_r_and_d = {2025: MagicMock(), 2026: MagicMock()}
    m.fs.costing.a_and_sl = MagicMock()
    m.fs.costing.fb = MagicMock()
    m.fs.costing.pt_and_i = MagicMock()

    # Mock get_units for all variables
    all_vars = (
        list(m.fs.costing.cash_flow.values())
        + list(m.fs.costing.main_product_profit.values())
        + list(m.fs.costing.total_profit.values())
        + list(m.fs.costing.total_byproduct_profit.values())
        + list(m.fs.costing.equipment_cost.values())
        + list(m.fs.costing.total_overnight_cost_expended.values())
        + list(m.fs.costing.aggregate_variable_operating_cost.values())
        + list(m.fs.costing.aggregate_fixed_operating_cost.values())
        + list(m.fs.costing.total_operating_expense.values())
        + list(m.fs.f.values())
        + list(m.fs.f_in.values())
        + list(m.fs.f_out.values())
        + list(m.fs.byproduct_valorization.byproduct_produced.values())
        + list(m.fs.environmental_impacts.total_yearly_impacts.values())
        + list(m.fs.environmental_impacts.option_yearly_impacts.values())
    )

    for var in all_vars:
        var.get_units.return_value = "kg"

    return m


@pytest.fixture
def mock_results():
    results = MagicMock()
    results.solver.status = MagicMock()
    results.solver.termination_condition = MagicMock()
    return results


# --- Test all branches and outputs of each function ---


@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_overview_optimal_status(mock_pyo_value, mock_print, mock_model, mock_results):
    # Model solved optimally, NPV objective
    def value_side_effect(obj):
        if obj is mock_model.fs.option_binary_var["A"]:
            return 1
        elif obj is mock_model.fs.option_binary_var["B"]:
            return 0
        elif obj is mock_model.fs.objective_function_choice:
            return ObjectiveFunctionChoice.NET_PRESENT_VALUE.value
        elif obj is mock_model.fs.costing.net_present_value:
            return 10000.0
        elif obj is mock_model.fs.environmental_impacts.total_impacts:
            return 55.0
        return 0.0

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_results_overview(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Chosen Pathway" in arg for arg in print_args)
    assert any("NPV" in arg for arg in print_args)
    assert any("Total Impacts" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_overview_COR_objective(mock_pyo_value, mock_print, mock_model, mock_results):
    # Model solved optimally, COR objective (not NET_PRESENT_VALUE)
    # First, let's get the actual NET_PRESENT_VALUE.value and use a different value
    npv_value = ObjectiveFunctionChoice.NET_PRESENT_VALUE.value
    cor_objective_value = npv_value + 1  # Ensure it's different

    def value_side_effect(obj):
        if obj is mock_model.fs.option_binary_var["A"]:
            return 0
        elif obj is mock_model.fs.objective_function_choice:
            return (
                cor_objective_value  # Use different value than NET_PRESENT_VALUE.value
            )
        elif obj is mock_model.fs.costing.cost_of_recovery:
            return 2.0
        elif obj is mock_model.fs.environmental_impacts.total_impacts:
            return 30.0
        return 0.0

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_results_overview(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    # The actual print statement is "Cost of Recovery: {value}"
    assert any("Cost of Recovery:" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_overview_not_optimal(mock_print, mock_model, mock_results):
    # Model NOT solved optimally
    mock_results.solver.status = SolverStatus.error
    mock_results.solver.termination_condition = TerminationCondition.infeasible

    report_superstructure_results_overview(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Model was not solved optimally" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_overview_no_results(mock_print, mock_model):
    # No solver results
    report_superstructure_results_overview(mock_model, None)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("No solver results" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_costing_full(mock_pyo_value, mock_print, mock_model, mock_results):
    def value_side_effect(obj):
        value_map = {
            mock_model.fs.objective_function_choice: ObjectiveFunctionChoice.NET_PRESENT_VALUE.value,
            mock_model.fs.costing.net_present_value: 10000.0,
            mock_model.fs.costing.discount_factor: 0.09,
            mock_model.fs.costing.i_operating_expense_escalation: 0.03,
            mock_model.fs.costing.i_capital_escalation: 0.02,
            mock_model.fs.costing.main_product_profit[2025]: 42.0,
            mock_model.fs.costing.total_profit[2025]: 52.0,
            mock_model.fs.costing.total_byproduct_profit[2025]: 12.0,
            mock_model.fs.costing.equipment_cost["A"]: 100.0,
            mock_model.fs.costing.equipment_cost["B"]: 0.0,
            mock_model.fs.costing.total_plant_cost: 99.0,
            mock_model.fs.costing.lang_factor: 1.5,
            mock_model.fs.costing.total_overnight_cost: 199.0,
            mock_model.fs.costing.financing_factor: 0.05,
            mock_model.fs.costing.other_costs_factor: 0.03,
            mock_model.fs.costing.total_overnight_cost_expended[2025]: 87.0,
            mock_model.fs.costing.total_operators: 5,
            mock_model.fs.costing.aggregate_variable_operating_cost[2025]: 17.0,
            mock_model.fs.costing.aggregate_fixed_operating_cost[2025]: 6.0,
            mock_model.fs.costing.total_operating_expense[2025]: 23.0,
            mock_model.fs.costing.cost_of_labor: 999.0,
            mock_model.fs.costing.m_and_sm: 89.0,
            mock_model.fs.costing.sa_and_qa_qc: 77.0,
            mock_model.fs.costing.s_ip_r_and_d[2025]: 45.0,
            mock_model.fs.costing.a_and_sl: 61.0,
            mock_model.fs.costing.fb: 39.0,
            mock_model.fs.costing.pt_and_i: 22.0,
        }
        return value_map.get(obj, 0.0)

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_costing(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("NPV" in arg for arg in print_args)
    assert any("Cash Flows" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_costing_not_optimal(mock_print, mock_model, mock_results):
    mock_results.solver.status = SolverStatus.error
    mock_results.solver.termination_condition = TerminationCondition.infeasible

    report_superstructure_costing(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Model was not solved optimally" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_costing_no_results(mock_print, mock_model):
    report_superstructure_costing(mock_model, None)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("No solver results" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_streams_full(mock_pyo_value, mock_print, mock_model, mock_results):
    def value_side_effect(obj):
        value_map = {
            mock_model.fs.f[("stage1", "comp1", 2025)]: 10.0,
            mock_model.fs.f_in[("A", "comp1", 2025)]: 5.0,
            mock_model.fs.f_out[("A", "comp1", 2025)]: 2.0,
            mock_model.fs.byproduct_valorization.byproduct_produced[
                ("byprod1", 2025)
            ]: 4.0,
        }
        return value_map.get(obj, 0.0)

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_streams(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Superstructure Streams" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_streams_not_optimal(mock_print, mock_model, mock_results):
    mock_results.solver.status = SolverStatus.error
    mock_results.solver.termination_condition = TerminationCondition.infeasible

    report_superstructure_streams(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Model was not solved optimally" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_streams_no_results(mock_print, mock_model):
    report_superstructure_streams(mock_model, None)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("No solver results" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_envimpacts_full(mock_pyo_value, mock_print, mock_model, mock_results):
    def value_side_effect(obj):
        value_map = {
            mock_model.fs.environmental_impacts.total_impacts: 80.0,
            mock_model.fs.environmental_impacts.epsilon: 1.1,
            mock_model.fs.environmental_impacts.total_yearly_impacts[2025]: 1.7,
            mock_model.fs.environmental_impacts.option_yearly_impacts[("A", 2025)]: 2.3,
        }
        return value_map.get(obj, 0.0)

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_environmental_impacts(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Superstructure Environmental Impacts" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_envimpacts_not_optimal(mock_print, mock_model, mock_results):
    mock_results.solver.status = SolverStatus.error
    mock_results.solver.termination_condition = TerminationCondition.infeasible

    report_superstructure_environmental_impacts(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Model was not solved optimally" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_envimpacts_no_results(mock_print, mock_model):
    report_superstructure_environmental_impacts(mock_model, None)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("No solver results" in arg for arg in print_args)


@patch("prommis.superstructure.report_superstructure_results.print")
def test_envimpacts_not_tracked(mock_print, mock_model, mock_results):
    # Remove environmental_impacts to test not tracked case
    del mock_model.fs.environmental_impacts
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_environmental_impacts(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any(
        "User has not specified model to track environmental impacts" in arg
        for arg in print_args
    )


# Additional edge case tests for better coverage
@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_costing_missing_byproduct(
    mock_pyo_value, mock_print, mock_model, mock_results
):
    # Test case where total_byproduct_profit doesn't exist
    del mock_model.fs.costing.total_byproduct_profit

    def value_side_effect(obj):
        return 10.0

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_costing(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Costing" in arg for arg in print_args)


# Test for COR objective in costing function as well
@patch("prommis.superstructure.report_superstructure_results.print")
@patch("prommis.superstructure.report_superstructure_results.pyo.value")
def test_costing_COR_objective(mock_pyo_value, mock_print, mock_model, mock_results):
    # Test Cost of Recovery objective in costing function
    npv_value = ObjectiveFunctionChoice.NET_PRESENT_VALUE.value
    cor_objective_value = npv_value + 1  # Ensure it's different

    def value_side_effect(obj):
        if obj is mock_model.fs.objective_function_choice:
            return cor_objective_value  # Not NET_PRESENT_VALUE
        elif obj is mock_model.fs.costing.cost_of_recovery:
            return 3.5
        return 10.0

    mock_pyo_value.side_effect = value_side_effect
    mock_results.solver.status = SolverStatus.ok
    mock_results.solver.termination_condition = TerminationCondition.optimal

    report_superstructure_costing(mock_model, mock_results)

    assert mock_print.called
    print_args = [str(call[0][0]) for call in mock_print.call_args_list]
    assert any("Cost of Recovery:" in arg for arg in print_args)
