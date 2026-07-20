####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Tests for the kinetic cell-balance Mountain Pass flotation flowsheet."""

from pyomo.environ import Constraint, assert_optimal_termination, value
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import UnitModelBlockData
from idaes.core.scaling import get_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError

import pytest

from prommis.flotation.bastnaesite_properties import COMPONENTS
from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation import (
    mountain_pass_flotation_kinetic_cell_balance_flowsheet as cell_balance,
)

solver = get_solver("ipopt")


def _build_seeded_model(feed_scale=1.0):
    model = cell_balance.build_model()
    cell_balance.set_model_inputs(model, feed_scale=feed_scale)
    cell_balance._seed_deterministic_streams(model, feed_scale=feed_scale)
    return model


def _build_deterministically_initialized_model(feed_scale=1.0):
    return cell_balance.build_and_initialize(
        feed_scale=feed_scale,
        initialization_method=cell_balance.DETERMINISTIC_STREAM_INITIALIZATION,
    )


def _unit_model_names(fs):
    return {
        component.local_name
        for component in fs.component_objects(descend_into=False)
        if any(isinstance(data, UnitModelBlockData) for data in component.values())
    }


@pytest.mark.component
@pytest.mark.build
def test_build_dof_units_and_diagnostics():
    model = _build_seeded_model()

    assert model.fs.scenario == recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model)
    DiagnosticsToolbox(model).assert_no_structural_warnings(
        ignore_evaluation_errors=True
    )


@pytest.mark.component
def test_topology_parity_with_fixed_recovery_flowsheet():
    fixed_recovery_model = recovery_flowsheet.build_model(
        scenario=recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO,
        expand_arcs=False,
    )
    cell_balance_model = cell_balance.build_model(expand_arcs=False)
    assert recovery_flowsheet.BANK_NAMES == cell_balance.BANK_NAMES
    assert _unit_model_names(fixed_recovery_model.fs) == _unit_model_names(
        cell_balance_model.fs
    )

    fixed_recovery_arcs = {
        arc.local_name: arc for arc in fixed_recovery_model.fs.component_objects(Arc)
    }
    cell_balance_arcs = {
        arc.local_name: arc for arc in cell_balance_model.fs.component_objects(Arc)
    }
    assert set(fixed_recovery_arcs) == set(cell_balance_arcs)
    for arc_name in fixed_recovery_arcs:
        assert (
            fixed_recovery_arcs[arc_name].source.name
            == cell_balance_arcs[arc_name].source.name
        )
        assert (
            fixed_recovery_arcs[arc_name].destination.name
            == cell_balance_arcs[arc_name].destination.name
        )


@pytest.mark.component
def test_json_loading_and_audit_summary_pre_solve():
    model = _build_deterministically_initialized_model()
    parameters = cell_balance.load_kinetic_cell_balance_parameters()
    summary = cell_balance.cell_balance_summary(model)

    assert set(parameters["banks"]) == set(cell_balance.BANK_NAMES)
    for bank_name, bank_summary in summary.items():
        assert bank_summary["audit"]["walk_completed"]
        assert bank_summary["audit"]["first_failure_mode"] is None
        for component in COMPONENTS:
            assert 0 <= bank_summary["recovery"][component] <= 1


@pytest.mark.unit
def test_synthetic_infeasible_bank_initialization_reports_context():
    model = cell_balance.build_model()
    cell_balance.set_model_inputs(model)
    bank = model.fs.rougher
    for component in COMPONENTS:
        bank.k_cb[0, component].fix(1e6)
    cell_balance._seed_deterministic_streams(model)

    with pytest.raises(
        ConfigurationError, match="algebraically infeasible"
    ) as exc_info:
        bank.default_initializer().initialize(bank)

    assert "rougher" in str(exc_info.value)
    assert "failure_mode" in str(exc_info.value)


@pytest.mark.component
def test_network_arc_flow_equalities_are_nominally_scaled():
    model = _build_seeded_model()

    scaled = cell_balance._scale_network_arc_flow_equalities(model)

    assert scaled > 0
    for constraint in model.component_data_objects(
        Constraint,
        descend_into=True,
        active=True,
    ):
        if "_expanded.flow_mass_comp_equality" in constraint.name:
            assert get_scaling_factor(constraint) is not None


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_kinetic_cell_balance_reproduces_table1_product_fit():
    model = _build_deterministically_initialized_model()
    results = cell_balance.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    report = cell_balance.report_results(model)
    metrics = report["product_metrics"]

    assert report["table1_products_within_tolerance"]
    assert metrics["overall_REO_recovery"] == pytest.approx(80.065, abs=1e-3)
    assert metrics["final_concentrate_REO_grade"] == pytest.approx(65.081, abs=0.05)
    assert metrics["final_concentrate_CaO_grade"] == pytest.approx(3.033, abs=0.05)
    assert metrics["final_concentrate_BaO_grade"] == pytest.approx(0.996, abs=0.05)
    assert metrics["final_concentrate_SrO_grade"] == pytest.approx(5.293, abs=0.05)
    assert all(
        bank["audit"]["walk_completed"]
        for bank in report["cell_balance_summary"].values()
    )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_unit_model_default_initialization_solves_table1_product_fit():
    model = cell_balance.build_and_initialize(
        initialization_method=cell_balance.UNIT_MODEL_DEFAULT_INITIALIZATION
    )
    results = cell_balance.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    report = cell_balance.report_results(model)
    assert report["table1_products_within_tolerance"]
    assert all(
        bank["audit"]["walk_completed"]
        for bank in report["cell_balance_summary"].values()
    )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_feed_scale_changes_cell_balance_recoveries_deterministically():
    base = _build_deterministically_initialized_model(feed_scale=1.0)
    scaled = _build_deterministically_initialized_model(feed_scale=2.0)
    assert_optimal_termination(cell_balance.solve_model(base, solver=solver))
    assert_optimal_termination(cell_balance.solve_model(scaled, solver=solver))

    assert value(scaled.fs.rougher.recovery[0, "REO"]) < value(
        base.fs.rougher.recovery[0, "REO"]
    )
    for bank_name in cell_balance.BANK_NAMES:
        base_bank = getattr(base.fs, bank_name)
        scaled_bank = getattr(scaled.fs, bank_name)
        for component in COMPONENTS:
            assert value(scaled_bank.k_cb[0, component]) == pytest.approx(
                value(base_bank.k_cb[0, component])
            )


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_solve_and_mass_closure():
    model = _build_deterministically_initialized_model()
    results = cell_balance.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    for component_residuals in cell_balance.mass_closure_residuals(model).values():
        for component in COMPONENTS:
            assert component_residuals[component] == pytest.approx(0, abs=1e-6)


@pytest.mark.component
def test_cell_balance_summary_handles_missing_rho_solid(monkeypatch):
    # When the reference rho_solid cannot be evaluated, cell_balance_summary
    # must report rho_solid_deviation_pct as None rather than dividing by None.
    model = _build_deterministically_initialized_model()
    original_safe_float = cell_balance._safe_float
    rho_solid_ids = {
        id(getattr(model.fs, bank_name).rho_solid)
        for bank_name in cell_balance.BANK_NAMES
    }

    def patched_safe_float(obj):
        if id(obj) in rho_solid_ids:
            return None
        return original_safe_float(obj)

    monkeypatch.setattr(cell_balance, "_safe_float", patched_safe_float)
    summary = cell_balance.cell_balance_summary(model)

    for bank_summary in summary.values():
        assert bank_summary["rho_solid"] is None
        assert bank_summary["rho_solid_deviation_pct"] is None
