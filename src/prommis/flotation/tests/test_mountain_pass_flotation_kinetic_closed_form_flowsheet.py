####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Tests for the kinetic Mountain Pass flotation flowsheet."""

import logging

from pyomo.environ import assert_optimal_termination, value
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import UnitModelBlockData
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.flotation.bastnaesite_properties import COMPONENTS
from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation import (
    mountain_pass_flotation_kinetic_closed_form_flowsheet as kinetic,
)

solver = get_solver("ipopt")
_log = logging.getLogger(__name__)


def _build_seeded_model(feed_scale=1.0):
    model = kinetic.build_model()
    kinetic.set_model_inputs(model, feed_scale=feed_scale)
    kinetic._seed_deterministic_streams(model, feed_scale=feed_scale)
    return model


def _build_deterministically_initialized_model(feed_scale=1.0):
    return kinetic.build_and_initialize(
        feed_scale=feed_scale,
        initialization_method=kinetic.DETERMINISTIC_STREAM_INITIALIZATION,
    )


def _unit_model_names(fs):
    return {
        component.local_name
        for component in fs.component_objects(descend_into=False)
        if any(isinstance(data, UnitModelBlockData) for data in component.values())
    }


def _bank_outlet_water_flow(bank, outlet_name, time=0):
    inlet_solid_flow = sum(
        value(bank.inlet.flow_mass_comp[time, component]) for component in COMPONENTS
    )
    outlet_solid_flow = sum(
        value(getattr(bank, outlet_name).flow_mass_comp[time, component])
        for component in COMPONENTS
    )
    return value(bank.flow_mass_water[time]) * outlet_solid_flow / inlet_solid_flow


def _apparent_water_for_bank_inlet(bank, port, time=0):
    solid_flow = sum(
        value(port.flow_mass_comp[time, component]) for component in COMPONENTS
    )
    solids_fraction = value(bank.pulp_solids_mass_fraction)
    return solid_flow * (1 - solids_fraction) / solids_fraction


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
def test_runtime_basis_matches_calibration_basis_pre_solve():
    base = _build_seeded_model(feed_scale=1.0)
    scaled = _build_seeded_model(feed_scale=2.0)
    kinetic_parameters = kinetic.load_kinetic_parameters()

    for bank_name, bank_parameters in kinetic_parameters["banks"].items():
        base_bank = getattr(base.fs, bank_name)
        scaled_bank = getattr(scaled.fs, bank_name)
        expected_tau_h = bank_parameters["tau_required_min"] / 60

        assert value(base_bank.tau[0]) == pytest.approx(expected_tau_h, rel=1e-4)
        assert value(scaled_bank.tau[0]) == pytest.approx(expected_tau_h / 2, rel=1e-4)
        assert value(scaled_bank.k_cf[0, "REO"]) == pytest.approx(
            value(base_bank.k_cf[0, "REO"])
        )


@pytest.mark.component
def test_fixed_recovery_helpers_do_not_fix_kinetic_recoveries():
    model = kinetic.build_model()
    kinetic.set_model_inputs(model)

    for bank_name in kinetic.BANK_NAMES:
        bank = getattr(model.fs, bank_name)
        for component in COMPONENTS:
            assert not bank.recovery[0, component].fixed


@pytest.mark.component
def test_topology_parity_with_fixed_recovery_flowsheet():
    fixed_recovery_model = recovery_flowsheet.build_model(
        scenario=recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO,
        expand_arcs=False,
    )
    kinetic_model = kinetic.build_model(expand_arcs=False)
    assert recovery_flowsheet.BANK_NAMES == kinetic.BANK_NAMES
    assert _unit_model_names(fixed_recovery_model.fs) == _unit_model_names(
        kinetic_model.fs
    )

    fixed_recovery_arcs = {
        arc.local_name: arc for arc in fixed_recovery_model.fs.component_objects(Arc)
    }
    kinetic_arcs = {
        arc.local_name: arc for arc in kinetic_model.fs.component_objects(Arc)
    }
    assert set(fixed_recovery_arcs) == set(kinetic_arcs)
    for arc_name in fixed_recovery_arcs:
        assert (
            fixed_recovery_arcs[arc_name].source.name
            == kinetic_arcs[arc_name].source.name
        )
        assert (
            fixed_recovery_arcs[arc_name].destination.name
            == kinetic_arcs[arc_name].destination.name
        )

    assert set(recovery_flowsheet.stream_port_map(fixed_recovery_model)) == set(
        recovery_flowsheet.stream_port_map(kinetic_model)
    )
    assert (
        fixed_recovery_model.fs.scenario
        == recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
    )
    assert kinetic_model.fs.scenario == recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO


@pytest.mark.component
def test_initialize_model_rejects_unsupported_method():
    model = kinetic.build_model()
    kinetic.set_model_inputs(model)

    with pytest.raises(ValueError, match="Unknown initialization method"):
        kinetic.initialize_model(
            model,
            initialization_method="unknown_method",
        )


@pytest.mark.component
def test_initialize_model_rejects_feed_scale_mismatch():
    model = kinetic.build_model()
    kinetic.set_model_inputs(model, feed_scale=1.0)

    with pytest.raises(ValueError, match="feed_scale must match"):
        kinetic.initialize_model(model, feed_scale=2.0)


@pytest.mark.component
@pytest.mark.parametrize(
    "initialization_method",
    (
        kinetic.UNIT_MODEL_DEFAULT_INITIALIZATION,
        kinetic.DETERMINISTIC_STREAM_INITIALIZATION,
    ),
)
def test_flowsheet_initialization_kinetic_populates_recoveries(initialization_method):
    model = kinetic.build_model()
    kinetic.set_model_inputs(model)
    kinetic.initialize_model(model, initialization_method=initialization_method)

    assert model.fs.initialization_method == initialization_method
    scenario_recoveries = recovery_flowsheet.SCENARIOS[
        recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
    ]["recoveries"]
    for bank_name, recoveries in scenario_recoveries.items():
        bank = getattr(model.fs, bank_name)
        for component, recovery in recoveries.items():
            assert value(bank.recovery[0, component]) == pytest.approx(
                recovery, abs=1e-3
            )


@pytest.mark.component
def test_water_non_conservation_is_documented(caplog):
    model = _build_deterministically_initialized_model()

    diagnostics = {
        "rougher_mixer": (
            _apparent_water_for_bank_inlet(
                model.fs.rougher, model.fs.rougher_mixer.fresh_feed
            )
            + _bank_outlet_water_flow(model.fs.scavenger, "concentrate")
            - value(model.fs.rougher.flow_mass_water[0])
        ),
        "cleaner1_mixer": (
            _bank_outlet_water_flow(model.fs.rougher, "concentrate")
            + _bank_outlet_water_flow(model.fs.cleaner2, "tails")
            - value(model.fs.cleaner1.flow_mass_water[0])
        ),
        "cleaner2_mixer": (
            _bank_outlet_water_flow(model.fs.cleaner1, "concentrate")
            + _bank_outlet_water_flow(model.fs.cleaners34, "tails")
            - value(model.fs.cleaner2.flow_mass_water[0])
        ),
        "final_tails_mixer": (
            _bank_outlet_water_flow(model.fs.rougher, "tails")
            + _bank_outlet_water_flow(model.fs.scavenger, "tails")
        ),
    }

    with caplog.at_level(logging.INFO, logger=__name__):
        for mixer_name, residual in diagnostics.items():
            _log.info(
                "Apparent water non-conservation diagnostic for %s: %.6g kg/h",
                mixer_name,
                residual,
            )

    assert any(abs(residual) > 1e-6 for residual in diagnostics.values())
    assert "Apparent water non-conservation diagnostic" in caplog.text


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_kinetic_reproduces_table1_product_fit():
    model = _build_deterministically_initialized_model()
    results = kinetic.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    report = kinetic.report_results(model)
    metrics = report["product_metrics"]

    assert report["table1_products_within_tolerance"]
    assert metrics["overall_REO_recovery"] == pytest.approx(80.065, abs=1e-3)
    assert metrics["final_concentrate_REO_grade"] == pytest.approx(65.081, abs=0.05)
    assert metrics["final_concentrate_CaO_grade"] == pytest.approx(3.033, abs=0.05)
    assert metrics["final_concentrate_BaO_grade"] == pytest.approx(0.996, abs=0.05)
    assert metrics["final_concentrate_SrO_grade"] == pytest.approx(5.293, abs=0.05)

    scenario_recoveries = recovery_flowsheet.SCENARIOS[
        recovery_flowsheet.TABLE1_PRODUCT_FIT_SCENARIO
    ]["recoveries"]
    for bank_name, recoveries in scenario_recoveries.items():
        bank = getattr(model.fs, bank_name)
        for component, recovery in recoveries.items():
            assert value(bank.recovery[0, component]) == pytest.approx(
                recovery, abs=1e-3
            )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_derived_tau_matches_calibration():
    model = _build_deterministically_initialized_model()
    results = kinetic.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    kinetic_parameters = kinetic.load_kinetic_parameters()
    for bank_name, bank_parameters in kinetic_parameters["banks"].items():
        bank = getattr(model.fs, bank_name)
        assert value(bank.tau[0]) * 60 == pytest.approx(
            bank_parameters["tau_required_min"],
            abs=1e-3,
        )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_feed_scale_changes_tau_not_k():
    base = _build_deterministically_initialized_model(feed_scale=1.0)
    scaled = _build_deterministically_initialized_model(feed_scale=2.0)
    assert_optimal_termination(kinetic.solve_model(base, solver=solver))
    assert_optimal_termination(kinetic.solve_model(scaled, solver=solver))

    base_rougher = base.fs.rougher
    scaled_rougher = scaled.fs.rougher
    assert value(scaled_rougher.tau[0]) < value(base_rougher.tau[0])
    for component in COMPONENTS:
        assert value(scaled_rougher.recovery[0, component]) < value(
            base_rougher.recovery[0, component]
        )

    for bank_name in kinetic.BANK_NAMES:
        base_bank = getattr(base.fs, bank_name)
        scaled_bank = getattr(scaled.fs, bank_name)
        for component in COMPONENTS:
            assert value(scaled_bank.k_cf[0, component]) == pytest.approx(
                value(base_bank.k_cf[0, component])
            )


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_no_numerical_warnings_kinetic():
    model = _build_deterministically_initialized_model()
    results = kinetic.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    DiagnosticsToolbox(model).assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_solve_and_mass_closure():
    model = _build_deterministically_initialized_model()
    results = kinetic.solve_model(model, solver=solver)
    assert_optimal_termination(results)

    for component_residuals in kinetic.mass_closure_residuals(model).values():
        for component in COMPONENTS:
            assert component_residuals[component] == pytest.approx(0, abs=1e-6)
