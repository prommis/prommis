#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the generic flotation bank unit model."""

import logging

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Var,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes import logger as idaeslog
from idaes.core import FlowsheetBlock, useDefault
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.scaling import get_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError, InitializationError

import pytest

from prommis.flotation import initializer as initializer_module
from prommis.flotation.bastnaesite_properties import (
    BastnaesiteParameters,
    COMPONENTS,
)
from prommis.flotation.flotation_bank import (
    FlotationBank,
    FlotationBankScaler,
)
from prommis.flotation.initializer import (
    FEED_FLOOR_KG_PER_H,
    CascadeInputError,
    CascadeInfeasibleError,
    FlotationBankKineticCellBalanceAnalyticalInitializer,
    FlotationBankKineticCellBalanceInitializer,
    FlotationBankKineticCellBalanceStagedInitializer,
    FlotationBankKineticClosedFormInitializer,
    _adaptive_bracket,
    _apply_zero_feed_recovery_policy,
    _solve_cell,
    audit_cascade_solvability,
    cell_cascade_forward,
)

solver = get_solver("ipopt")

FEED = {
    "REO": 10.0,
    "CaO": 20.0,
    "BaO": 15.0,
    "SrO": 5.0,
    "inert_gangue": 50.0,
}
RECOVERY = {
    "REO": 0.8,
    "CaO": 0.2,
    "BaO": 0.3,
    "SrO": 0.4,
    "inert_gangue": 0.1,
}


def _fix_kinetic_inputs(unit, tau_h=2.0, k_per_h=1.5):
    unit.air_holdup.fix(0.0)
    unit.pulp_solids_mass_fraction.fix(0.5)
    unit.cell_volume.fix(
        tau_h * value(unit.flow_vol_slurry[0]) / unit.config.number_of_cells
    )
    for component in COMPONENTS:
        unit.k_cf[0, component].fix(k_per_h)


def _fix_cell_balance_inputs(unit, cell_volume=0.02, k_per_h=0.5):
    unit.air_holdup.fix(0.1)
    unit.pulp_solids_mass_fraction.fix(0.5)
    unit.rho_solid.set_value(3000.0)
    unit.cell_volume.fix(cell_volume)
    for component in COMPONENTS:
        unit.k_cb[0, component].fix(k_per_h)


def _bank_model(
    fix_recoveries=True,
    recovery_basis="fixed",
    number_of_cells=1,
    fix_kinetic=False,
):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    m.fs.unit = FlotationBank(
        property_package=m.fs.properties,
        recovery_basis=recovery_basis,
        number_of_cells=number_of_cells,
    )
    for component, flow in FEED.items():
        m.fs.unit.inlet.flow_mass_comp[0, component].fix(flow)
    if fix_recoveries:
        for component, recovery in RECOVERY.items():
            m.fs.unit.recovery[0, component].fix(recovery)
    if fix_kinetic:
        _fix_kinetic_inputs(m.fs.unit)
    return m


@pytest.fixture()
def model():
    return _bank_model()


@pytest.fixture()
def model_with_free_recoveries():
    return _bank_model(fix_recoveries=False)


@pytest.mark.unit
@pytest.mark.build
def test_build(model_with_free_recoveries):
    unit = model_with_free_recoveries.fs.unit
    assert isinstance(unit.recovery, Var)
    assert isinstance(unit.concentrate_split_eq, Constraint)
    assert isinstance(unit.tails_split_eq, Constraint)
    assert unit.default_scaler is FlotationBankScaler
    assert unit.default_initializer is BlockTriangularizationInitializer
    assert unit.config.recovery_basis == "fixed"
    assert not unit.recovery[0, "REO"].fixed
    assert "bank_name" not in unit.config
    assert "calibration_mode" not in unit.config
    for attr in (
        "cell_volume",
        "air_holdup",
        "pulp_solids_mass_fraction",
        "k_cf",
        "R_inf",
        "tau",
        "effective_volume",
        "flow_vol_slurry",
        "kinetic_recovery_eq",
    ):
        assert not hasattr(unit, attr)


@pytest.mark.unit
@pytest.mark.build
def test_build_kinetic_mode_adds_state():
    model = _bank_model(fix_recoveries=False, recovery_basis="kinetic_closed_form")
    unit = model.fs.unit

    assert unit.config.recovery_basis == "kinetic_closed_form"
    assert unit.default_initializer is FlotationBankKineticClosedFormInitializer
    assert isinstance(unit.recovery, Var)
    assert isinstance(unit.cell_volume, Var)
    assert isinstance(unit.air_holdup, Var)
    assert isinstance(unit.pulp_solids_mass_fraction, Var)
    assert isinstance(unit.k_cf, Var)
    assert isinstance(unit.R_inf, Var)
    assert isinstance(unit.effective_volume, Expression)
    assert isinstance(unit.flow_vol_slurry, Expression)
    assert isinstance(unit.tau, Expression)
    assert isinstance(unit.kinetic_recovery_eq, Constraint)


@pytest.mark.unit
@pytest.mark.build
def test_build_kinetic_cell_balance_mode_adds_cell_state():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit

    assert unit.config.recovery_basis == "kinetic_cell_balance"
    assert unit.default_initializer is FlotationBankKineticCellBalanceInitializer
    assert isinstance(unit.recovery, Var)
    assert isinstance(unit.cell_volume, Var)
    assert isinstance(unit.air_holdup, Var)
    assert isinstance(unit.pulp_solids_mass_fraction, Var)
    assert isinstance(unit.k_cb, Var)
    assert isinstance(unit.cell_total_solid_holdup, Var)
    assert isinstance(unit.cell_solid_holdup, Var)
    assert isinstance(unit.cell_pulp_out_flow, Var)
    assert isinstance(unit.cell_float_flow, Var)
    assert isinstance(unit.rho_slurry, Expression)
    assert isinstance(unit.geometric_holdup_eq, Constraint)
    assert isinstance(unit.recovery_eq, Constraint)
    assert not hasattr(unit, "concentrate_split_eq")
    assert not hasattr(unit, "tails_split_eq")
    assert not hasattr(unit, "cell_recovery")
    assert not hasattr(unit, "cell_mass_pull")
    assert unit._well_mixed_omitted_component == "inert_gangue"
    assert unit._zero_feed_recovery_fixes == set()
    assert unit.cell_volume.lb == pytest.approx(1e-6)
    assert unit.air_holdup.lb == pytest.approx(0)
    assert unit.air_holdup.ub == pytest.approx(0.5)
    assert unit.pulp_solids_mass_fraction.lb == pytest.approx(1e-3)
    assert unit.pulp_solids_mass_fraction.ub == pytest.approx(1 - 1e-3)
    assert unit.k_cb[0, "REO"].lb == pytest.approx(0)


@pytest.mark.unit
@pytest.mark.build
def test_number_of_cells_unused_in_fixed_mode_equations():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="fixed",
        number_of_cells=12,
    )
    unit = model.fs.unit

    assert unit.config.number_of_cells == 12
    for attr in (
        "cell_volume",
        "air_holdup",
        "pulp_solids_mass_fraction",
        "k_cf",
        "R_inf",
        "tau",
        "effective_volume",
        "flow_vol_slurry",
        "kinetic_recovery_eq",
    ):
        assert not hasattr(unit, attr)


@pytest.mark.unit
def test_number_of_cells_config_rejects_non_positive_int():
    for invalid_value in (0, -1, 1.5):
        with pytest.raises((ValueError, ConfigurationError)):
            _bank_model(
                fix_recoveries=False,
                recovery_basis="fixed",
                number_of_cells=invalid_value,
            )


@pytest.mark.unit
def test_kinetic_cell_balance_exception_construction():
    input_error = CascadeInputError("msg", context="bank=foo")
    assert input_error.context == "bank=foo"

    infeasible = CascadeInfeasibleError(
        "msg",
        cell_index=3,
        F_zero=0.5,
        lower_margin=0.9,
        deficit=0.1,
        failure_mode="no_carrier_margin",
        context="bank=bar",
    )
    assert infeasible.cell_index == 3
    assert infeasible.F_zero == pytest.approx(0.5)
    assert infeasible.lower_margin == pytest.approx(0.9)
    assert infeasible.deficit == pytest.approx(0.1)
    assert infeasible.failure_mode == "no_carrier_margin"
    assert infeasible.context == "bank=bar"
    assert "cell_index=3" in str(infeasible)
    assert "lower_margin=0.9" in str(infeasible)
    assert "deficit=0.1" in str(infeasible)


@pytest.mark.unit
def test_kinetic_cell_balance_rejects_dynamic_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=units.hour)
    m.fs.properties = BastnaesiteParameters()

    with pytest.raises(ConfigurationError, match="steady-state"):
        m.fs.unit = FlotationBank(
            property_package=m.fs.properties,
            recovery_basis="kinetic_cell_balance",
        )


@pytest.mark.unit
@pytest.mark.parametrize("has_holdup", [True, False, useDefault])
def test_kinetic_cell_balance_does_not_validate_has_holdup(has_holdup):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    m.fs.unit = FlotationBank(
        property_package=m.fs.properties,
        recovery_basis="kinetic_cell_balance",
        has_holdup=has_holdup,
    )

    expected = False if has_holdup is useDefault else has_holdup
    assert m.fs.unit.config.has_holdup is expected


@pytest.mark.unit
def test_dof_kinetic_mode():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model)
    DiagnosticsToolbox(model).assert_no_structural_warnings()


@pytest.mark.component
def test_dof_kinetic_cell_balance_mode_after_initialization():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit
    _fix_cell_balance_inputs(unit)

    unit.default_initializer().initialize(unit)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model)
    DiagnosticsToolbox(model).assert_no_structural_warnings()


@pytest.mark.unit
def test_kinetic_cell_balance_initializer_requires_singleton_time():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0, 1])
    m.fs.properties = BastnaesiteParameters()
    m.fs.unit = FlotationBank(
        property_package=m.fs.properties,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=1,
    )

    with pytest.raises(ConfigurationError, match="singleton time set"):
        FlotationBankKineticCellBalanceAnalyticalInitializer().initialization_routine(
            m.fs.unit
        )


def _temporary_fixed_low_feed_cell_balance_model():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit
    _fix_cell_balance_inputs(unit, cell_volume=0.02, k_per_h=0.5)
    initially_unfixed = set()
    for component in COMPONENTS:
        var = unit.inlet.flow_mass_comp[0, component]
        var.unfix()
        var.set_value(1.0)
        initially_unfixed.add(id(var))
        var.fix()
    model.fs._kinetic_cell_balance_initially_unfixed = initially_unfixed
    return model


@pytest.mark.component
def test_kinetic_cell_balance_initializer_scales_temporary_fixed_inlet_guess():
    model = _temporary_fixed_low_feed_cell_balance_model()
    unit = model.fs.unit

    initializer = FlotationBankKineticCellBalanceInitializer()
    initializer.initialize(unit)

    assert initializer.summary[unit]["status"] == InitializationStatus.Ok
    assert value(unit.inlet.flow_mass_comp[0, "REO"]) > 1.0


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_kinetic_cell_balance_staged_initializer_scales_temporary_fixed_inlet_guess(
    caplog,
    monkeypatch,
):
    model = _temporary_fixed_low_feed_cell_balance_model()
    unit = model.fs.unit

    def fail_analytical_seed(*args, **kwargs):
        raise AssertionError("staged initializer should not use analytical seed")

    monkeypatch.setattr(
        initializer_module,
        "_feasible_cell_balance_cascade_seed",
        fail_analytical_seed,
    )
    monkeypatch.setattr(
        initializer_module,
        "cell_cascade_forward",
        fail_analytical_seed,
    )

    initializer = FlotationBankKineticCellBalanceStagedInitializer()
    # idaes.init does not propagate to the root logger that caplog attaches to,
    # so re-enable propagation for the duration of this test (monkeypatch restores
    # it afterwards) to make the staged step messages visible to caplog.
    monkeypatch.setattr(logging.getLogger("idaes.init"), "propagate", True)
    with caplog.at_level(idaeslog.INFO_HIGH, logger="idaes.init"):
        initializer.initialize(unit, output_level=idaeslog.INFO_HIGH)

    assert initializer.summary[unit]["status"] == InitializationStatus.Ok
    assert value(unit.inlet.flow_mass_comp[0, "REO"]) > 1.0
    assert "Initialization Step 1: validation complete." in caplog.text
    assert "Initialization Step 2: rough staged seed applied." in caplog.text
    assert "Initialization Step 3 (cell-balance sub-block)" in caplog.text
    assert "Initialization Step 4 (full bank)" in caplog.text
    assert "Initialization Step 5 (full bank): residual check complete" in caplog.text


@pytest.mark.unit
def test_scaler_sets_kinetic_cell_balance_factors_without_init():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    bank = model.fs.unit
    _fix_cell_balance_inputs(bank)

    scaler = FlotationBankScaler()
    scaler.variable_scaling_routine(bank)
    scaler.constraint_scaling_routine(bank)

    assert get_scaling_factor(bank.cell_solid_holdup[0, 1, "REO"]) > 0
    assert get_scaling_factor(bank.cell_pulp_out_flow[0, 1, "REO"]) > 0
    assert get_scaling_factor(bank.k_cb[0, "REO"]) > 0
    assert get_scaling_factor(bank.cell_total_solid_holdup[0, 1]) > 0
    assert get_scaling_factor(bank.well_mixed_eq[0, 1, "REO"]) > 0
    assert get_scaling_factor(bank.flotation_removal_eq[0, 1, "REO"]) > 0
    assert get_scaling_factor(bank.geometric_holdup_eq[0, 1]) > 0
    assert get_scaling_factor(bank.recovery_eq[0, "REO"]) > 0


@pytest.mark.unit
def test_kinetic_recovery_one_cell():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    unit = model.fs.unit

    unit.default_initializer().initialize(unit)

    assert value(unit.tau[0]) == pytest.approx(2.0)
    expected = 1.5 * 2.0 / (1 + 1.5 * 2.0)
    assert value(unit.recovery[0, "REO"]) == pytest.approx(expected)


@pytest.mark.unit
@pytest.mark.parametrize("number_of_cells", [2, 3, 4, 5])
def test_kinetic_recovery_intermediate_cell_counts(number_of_cells):
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        number_of_cells=number_of_cells,
        fix_kinetic=True,
    )
    unit = model.fs.unit
    _fix_kinetic_inputs(unit, tau_h=1.25, k_per_h=0.8)

    unit.default_initializer().initialize(unit)

    expected = 1 - (
        1
        + value(unit.k_cf[0, "REO"]) * value(unit.tau[0]) / unit.config.number_of_cells
    ) ** (-unit.config.number_of_cells)
    assert value(unit.recovery[0, "REO"]) == pytest.approx(expected)


@pytest.mark.unit
def test_tau_responds_to_throughput():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    unit = model.fs.unit
    base_tau = value(unit.tau[0])

    for component, flow in FEED.items():
        unit.inlet.flow_mass_comp[0, component].set_value(2 * flow)

    assert value(unit.tau[0]) == pytest.approx(base_tau / 2)


@pytest.mark.unit
def test_R_inf_is_fixed_var_for_current_kinetic_model():
    model = _bank_model(fix_recoveries=False, recovery_basis="kinetic_closed_form")

    assert model.fs.unit.R_inf[0, "REO"].fixed
    assert model.fs.unit.R_inf[0, "REO"].lb == pytest.approx(0.0)
    # Upper bound is intentionally 1 + 1e-8 (not exactly 1.0) so that the
    # default fixed-at-one R_inf does not trigger a diagnostics warning.
    assert model.fs.unit.R_inf[0, "REO"].ub == pytest.approx(1 + 1e-8, abs=1e-12)
    assert value(model.fs.unit.R_inf[0, "REO"]) == pytest.approx(1.0)

    model.fs.unit.R_inf[0, "REO"].unfix()
    model.fs.unit.R_inf[0, "REO"].set_value(0.5)
    assert value(model.fs.unit.R_inf[0, "REO"]) == pytest.approx(0.5)


@pytest.mark.unit
def test_rho_water_rejects_non_positive():
    model = _bank_model(fix_recoveries=False, recovery_basis="kinetic_closed_form")

    with pytest.raises(ValueError):
        model.fs.unit.rho_water.set_value(0)
    with pytest.raises(ValueError):
        model.fs.unit.rho_water.set_value(-1)


@pytest.mark.unit
def test_zero_total_feed_raises_in_kinetic_mode():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    for component in COMPONENTS:
        model.fs.unit.inlet.flow_mass_comp[0, component].set_value(0)

    with pytest.raises(ConfigurationError, match="positive inlet dry-solids flow"):
        FlotationBankKineticClosedFormInitializer().initialization_routine(
            model.fs.unit
        )


@pytest.mark.unit
def test_unsupported_dof_state_raises():
    model = _bank_model(
        fix_recoveries=True,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )

    with pytest.raises(ConfigurationError, match="Unsupported DOF state"):
        FlotationBankKineticClosedFormInitializer().initialization_routine(
            model.fs.unit
        )


@pytest.mark.unit
def test_kinetic_initializer_rejects_fixed_recovery_free_k():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    for component in COMPONENTS:
        model.fs.unit.k_cf[0, component].unfix()
        model.fs.unit.recovery[0, component].fix(0.5)

    with pytest.raises(ConfigurationError, match="prediction mode only"):
        FlotationBankKineticClosedFormInitializer().initialization_routine(
            model.fs.unit
        )


@pytest.mark.unit
def test_kinetic_recovery_twelve_cells():
    tau_h = 5.5 / 60
    target_recovery = 0.944
    k_per_h = 12 / tau_h * ((1 - target_recovery) ** (-1 / 12) - 1)
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        number_of_cells=12,
        fix_kinetic=True,
    )
    unit = model.fs.unit
    _fix_kinetic_inputs(unit, tau_h=tau_h, k_per_h=k_per_h)

    unit.default_initializer().initialize(unit)
    expected = 1 - (1 + value(unit.k_cf[0, "REO"]) * value(unit.tau[0]) / 12) ** -12

    assert value(unit.tau[0]) == pytest.approx(tau_h)
    assert expected == pytest.approx(target_recovery)
    assert value(unit.recovery[0, "REO"]) == pytest.approx(target_recovery)


@pytest.mark.unit
def test_kinetic_plug_flow_limit():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        number_of_cells=200,
        fix_kinetic=True,
    )
    unit = model.fs.unit
    _fix_kinetic_inputs(unit, tau_h=0.25, k_per_h=2.0)

    tanks_in_series = 1 - (
        1
        + value(unit.k_cf[0, "REO"]) * value(unit.tau[0]) / unit.config.number_of_cells
    ) ** (-unit.config.number_of_cells)

    assert tanks_in_series == pytest.approx(1 - 2.718281828459045**-0.5, abs=1e-3)


@pytest.mark.unit
@pytest.mark.parametrize(
    "F_in,M_total,k_hour,message",
    [
        ({"A": 1.0}, 1.0, {"A": 0.1}, "sequence"),
        ({"A": 1.0}, [], {"A": 0.1}, "empty"),
        ({"A": 1.0}, [0.0], {"A": 0.1}, "positive"),
        ({"A": -1.0}, [1.0], {"A": 0.1}, "non-negative"),
        ({"A": 1.0}, [1.0], {"A": -0.1}, "non-negative"),
        ({"A": "1.0"}, [1.0], {"A": 0.1}, "non-numeric"),
        ({"A": float("nan")}, [1.0], {"A": 0.1}, "non-finite"),
        ({"A": 1.0}, [float("inf")], {"A": 0.1}, "non-finite"),
        (
            {"A": 0.4e-12, "B": 0.4e-12, "C": 0.4e-12},
            [1.0],
            {"A": 0, "B": 0, "C": 0},
            "sanitized_total",
        ),
        ({"A": 1.0}, [1.0], {"B": 0.1}, "missing"),
    ],
)
def test_cell_cascade_forward_static_validation(F_in, M_total, k_hour, message):
    with pytest.raises(CascadeInputError, match=message) as exc_info:
        cell_cascade_forward(F_in, M_total, k_hour, context="bank=test")

    assert "bank=test" in str(exc_info.value)


@pytest.mark.unit
def test_cell_cascade_forward_analytical_cases():
    F_pulp, M_comp, F_float = cell_cascade_forward({"A": 1.0}, [1.0], {"A": 0.5})
    assert F_pulp[0]["A"] == pytest.approx(0.5)
    assert M_comp[0]["A"] == pytest.approx(1.0)
    assert F_float[0]["A"] == pytest.approx(0.5)

    with pytest.raises(CascadeInfeasibleError, match="no_carrier_margin"):
        cell_cascade_forward({"A": 1.0}, [1.0], {"A": 1.0})

    F_pulp, _, _ = cell_cascade_forward({"A": 1.0}, [1.0], {"A": 0.999})
    assert F_pulp[0]["A"] == pytest.approx(0.001)

    F_pulp, _, _ = cell_cascade_forward(
        {"A": 1.0, "B": 0.5}, [1.0], {"A": 10.0, "B": 0.0}
    )
    assert sum(F_pulp[0].values()) == pytest.approx((-8.5 + 92.25**0.5) / 2)

    F_pulp, M_comp, F_float = cell_cascade_forward(
        {"A": 1.0, "B": 1.0}, [2.0], {"A": 0.0, "B": 0.0}
    )
    assert F_pulp[0] == {"A": 1.0, "B": 1.0}
    assert M_comp[0]["A"] == pytest.approx(1.0)
    assert F_float[0] == {"A": 0.0, "B": 0.0}

    F_pulp, _, F_float = cell_cascade_forward(
        {"A": 1.0, "B": 0.5}, [1.0], {"A": 10.0, "B": 1e-15}
    )
    assert F_float[0]["B"] == pytest.approx(1e-15 * 1.0, rel=1e-6)
    assert F_pulp[0]["B"] < 0.5

    F_pulp, _, _ = cell_cascade_forward(
        {"A": 1.0, "B": 1.0}, [1.0], {"A": 1.0, "B": 0.0}
    )
    assert sum(F_pulp[0].values()) == pytest.approx((1 + 5**0.5) / 2)

    with pytest.raises(CascadeInfeasibleError, match="no_carrier_margin"):
        cell_cascade_forward({"A": 1.0, "B": 1e-15}, [1.0], {"A": 1.0, "B": 0.0})


@pytest.mark.unit
def test_cell_cascade_forward_low_flow_root_tolerance():
    F_in = {"A": 2 * FEED_FLOOR_KG_PER_H}
    F_pulp, M_comp, F_float = cell_cascade_forward(F_in, [1e-12], {"A": 0.5})

    assert F_in["A"] == pytest.approx(F_pulp[0]["A"] + F_float[0]["A"])
    assert F_float[0]["A"] == pytest.approx(0.5 * M_comp[0]["A"])


@pytest.mark.unit
def test_adaptive_bracket_refinement_and_failure_diagnostics():
    def g_near_root(T):
        return 1e-15 / T - 1.0

    T_low, T_high = _adaptive_bracket(
        g_near_root,
        F_zero=0.0,
        F_total=1.0,
        cell_index=1,
        context="synthetic",
    )
    assert T_low > 0
    assert T_low <= 1e-15 * 100
    assert g_near_root(T_low) > 0
    assert T_high == pytest.approx(1.0)

    def g_no_root(T):
        return -1.0

    with pytest.raises(CascadeInfeasibleError) as exc_info:
        _adaptive_bracket(
            g_no_root,
            F_zero=0.0,
            F_total=1.0,
            cell_index=2,
            context="synthetic",
            lower_margin=0.5,
            deficit=0.5,
        )
    assert exc_info.value.cell_index == 2
    assert exc_info.value.lower_margin == pytest.approx(0.5)
    assert exc_info.value.deficit == pytest.approx(0.5)
    assert exc_info.value.failure_mode == "near_singular_bracket"


@pytest.mark.unit
def test_solve_cell_depleted_branch():
    with pytest.raises(CascadeInfeasibleError) as exc_info:
        _solve_cell(
            {"A": 0.4e-12, "B": 0.4e-12},
            1.0,
            {"A": 0.1, "B": 0.0},
            2,
            "bank=test",
            collect_diagnostics=False,
        )

    assert exc_info.value.cell_index == 2
    assert exc_info.value.F_zero == pytest.approx(0.0)
    assert exc_info.value.lower_margin is None
    assert exc_info.value.failure_mode == "depleted_cell"
    assert "depleted below feed floor" in str(exc_info.value)


@pytest.mark.unit
def test_audit_cascade_solvability_reports_walk_state():
    success = audit_cascade_solvability({"A": 1.0}, [1.0], {"A": 0.5}, "bank")
    assert success["walk_completed"]
    assert success["first_failure_mode"] is None

    infeasible = audit_cascade_solvability({"A": 1.0}, [1.0], {"A": 1.0}, "bank")
    assert not infeasible["walk_completed"]
    assert infeasible["first_infeasible_cell"] == 1
    assert infeasible["first_failure_mode"] == "no_carrier_margin"


@pytest.mark.component
def test_prediction_mode_pre_computes_recovery():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    unit = model.fs.unit

    unit.default_initializer().initialize(unit)

    assert value(unit.recovery[0, "REO"]) == pytest.approx(0.75)


@pytest.mark.unit
def test_kinetic_performance_contents_flatten_component_rates():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    variables = model.fs.unit._get_performance_contents(time_point=0)["vars"]

    assert "k_cf" not in variables
    for component in COMPONENTS:
        assert variables[f"k_cf[{component}]"] == pytest.approx(1.5)


@pytest.mark.unit
def test_dof_and_units(model_with_free_recoveries):
    assert degrees_of_freedom(model_with_free_recoveries) == len(COMPONENTS)
    for component, recovery in RECOVERY.items():
        model_with_free_recoveries.fs.unit.recovery[0, component].fix(recovery)
    assert degrees_of_freedom(model_with_free_recoveries) == 0
    assert_units_consistent(model_with_free_recoveries)
    DiagnosticsToolbox(model_with_free_recoveries).assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_no_numerical_warnings_kinetic_bank():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    model.fs.unit.air_holdup.fix(0.1)
    model.fs.unit.cell_volume.fix(
        2.0 * value(model.fs.unit.flow_vol_slurry[0]) / (1 - 0.1)
    )
    model.fs.unit.default_initializer().initialize(model.fs.unit)
    assert_optimal_termination(solver.solve(model))

    DiagnosticsToolbox(model).assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_no_numerical_warnings_kinetic_cell_balance_bank():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit
    _fix_cell_balance_inputs(unit)

    unit.default_initializer().initialize(unit)
    assert_optimal_termination(solver.solve(model))

    DiagnosticsToolbox(model).assert_no_numerical_warnings()
    for constraint in (
        unit.cell_material_balance_eq,
        unit.flotation_removal_eq,
        unit.concentrate_outlet_eq,
        unit.tails_outlet_eq,
        unit.recovery_eq,
    ):
        for condata in constraint.values():
            assert abs(value(condata.body - condata.lower)) <= 1e-10


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_corner_recoveries_mass_closure(model):
    unit = model.fs.unit
    for recovery_value in (0, 1):
        for component in COMPONENTS:
            unit.recovery[0, component].fix(recovery_value)
        assert_optimal_termination(solver.solve(model))
        for component in COMPONENTS:
            inlet = value(unit.inlet.flow_mass_comp[0, component])
            assert value(
                unit.concentrate.flow_mass_comp[0, component]
            ) == pytest.approx(recovery_value * inlet, abs=1e-6)
            assert value(unit.tails.flow_mass_comp[0, component]) == pytest.approx(
                (1 - recovery_value) * inlet, abs=1e-6
            )


@pytest.mark.component
def test_default_initializer(model):
    initializer = model.fs.unit.default_initializer()
    initializer.initialize(model.fs.unit)
    assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_split_regression(model):
    results = solver.solve(model)
    assert_optimal_termination(results)
    DiagnosticsToolbox(model).assert_no_numerical_warnings()

    unit = model.fs.unit
    for component in COMPONENTS:
        assert value(unit.concentrate.flow_mass_comp[0, component]) == pytest.approx(
            RECOVERY[component] * FEED[component], abs=1e-7
        )
        assert value(unit.tails.flow_mass_comp[0, component]) == pytest.approx(
            (1 - RECOVERY[component]) * FEED[component], abs=1e-7
        )
    assert value(unit.solid_mass_pull[0]) == pytest.approx(0.235)


@pytest.mark.unit
def test_cell_cascade_forward_multi_cell_marches_feed_forward():
    # Single carrier, two identical cells. For one component each cell removes
    # exactly k * M_total, so the pulp outlet steps down by k*M per cell. This
    # verifies cell i+1 is fed by cell i's pulp outlet (the forward march).
    F_pulp, M_comp, F_float = cell_cascade_forward({"A": 10.0}, [2.0, 2.0], {"A": 1.0})

    # Cell 1: T = F_in - k*M = 10 - 2 = 8.
    assert F_pulp[0]["A"] == pytest.approx(8.0)
    assert F_float[0]["A"] == pytest.approx(2.0)
    # Cell 2 is fed 8.0 (cell 1's outlet): T = 8 - 2 = 6.
    assert F_pulp[1]["A"] == pytest.approx(6.0)
    assert F_float[1]["A"] == pytest.approx(2.0)

    # F_float == k * M_comp must hold per cell.
    for cell in range(2):
        assert F_float[cell]["A"] == pytest.approx(1.0 * M_comp[cell]["A"])

    # Pulp outlet strictly decreases along the cascade, and recovery compounds:
    # two cells float more total than one.
    assert F_pulp[1]["A"] < F_pulp[0]["A"]
    _, _, one_cell_float = cell_cascade_forward({"A": 10.0}, [2.0], {"A": 1.0})
    assert F_float[0]["A"] + F_float[1]["A"] > one_cell_float[0]["A"]


@pytest.mark.unit
def test_cell_cascade_forward_propagates_downstream_infeasibility():
    # Cell 1 is feasible (margin 1.5) but strips enough feed that cell 2 can no
    # longer support a positive root (margin 0.5). The cascade walk must reach
    # cell 2 and raise with the correct 1-based index. (The depleted_cell mode
    # is covered separately by test_solve_cell_depleted_branch.)
    with pytest.raises(CascadeInfeasibleError) as exc_info:
        cell_cascade_forward({"A": 1.5}, [1.0, 1.0], {"A": 1.0})

    assert exc_info.value.cell_index == 2
    assert exc_info.value.failure_mode == "no_carrier_margin"


def _cell_balance_bank_for_policy():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(model.fs.unit)
    return model


@pytest.mark.unit
def test_zero_feed_recovery_fixed_to_zero():
    model = _cell_balance_bank_for_policy()
    unit = model.fs.unit
    unit.inlet.flow_mass_comp[0, "SrO"].fix(0.0)

    _apply_zero_feed_recovery_policy(unit, 0, list(COMPONENTS))

    assert unit.recovery[0, "SrO"].fixed
    assert value(unit.recovery[0, "SrO"]) == pytest.approx(0.0)
    assert (0, "SrO") in unit._zero_feed_recovery_fixes
    # Components with positive feed are untouched.
    assert not unit.recovery[0, "REO"].fixed


@pytest.mark.unit
def test_zero_feed_with_user_fixed_nonzero_recovery_raises():
    model = _cell_balance_bank_for_policy()
    unit = model.fs.unit
    unit.inlet.flow_mass_comp[0, "SrO"].fix(0.0)
    unit.recovery[0, "SrO"].fix(0.5)

    with pytest.raises(ConfigurationError, match="recovery is fixed"):
        _apply_zero_feed_recovery_policy(unit, 0, list(COMPONENTS))


@pytest.mark.unit
def test_zero_feed_recovery_unfixed_when_feed_returns():
    model = _cell_balance_bank_for_policy()
    unit = model.fs.unit

    unit.inlet.flow_mass_comp[0, "SrO"].fix(0.0)
    _apply_zero_feed_recovery_policy(unit, 0, list(COMPONENTS))
    assert (0, "SrO") in unit._zero_feed_recovery_fixes

    unit.inlet.flow_mass_comp[0, "SrO"].fix(5.0)
    _apply_zero_feed_recovery_policy(unit, 0, list(COMPONENTS))

    assert not unit.recovery[0, "SrO"].fixed
    assert (0, "SrO") not in unit._zero_feed_recovery_fixes


@pytest.mark.unit
def test_scaler_sets_kinetic_closed_form_factors():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    bank = model.fs.unit

    scaler = FlotationBankScaler()
    scaler.variable_scaling_routine(bank)
    scaler.constraint_scaling_routine(bank)

    assert get_scaling_factor(bank.k_cf[0, "REO"]) > 0
    assert get_scaling_factor(bank.cell_volume) > 0
    assert get_scaling_factor(bank.air_holdup) > 0
    assert get_scaling_factor(bank.pulp_solids_mass_fraction) > 0
    assert get_scaling_factor(bank.kinetic_recovery_eq[0, "REO"]) > 0


@pytest.mark.unit
def test_closed_form_scaler_is_nan_safe():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_closed_form",
        fix_kinetic=True,
    )
    bank = model.fs.unit
    # Simulate a partially-evaluated/NaN state before scaling. The closed-form
    # branch must fall back to its default nominal (via _safe_value) and still
    # produce a finite, positive scaling factor rather than a NaN factor. A NaN
    # factor would fail the "> 0" assertions below.
    bank.k_cf[0, "REO"].set_value(float("nan"), skip_validation=True)
    bank.cell_volume.set_value(float("nan"), skip_validation=True)

    scaler = FlotationBankScaler()
    scaler.variable_scaling_routine(bank)

    assert get_scaling_factor(bank.k_cf[0, "REO"]) > 0
    assert get_scaling_factor(bank.cell_volume) > 0


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_staged_and_analytical_initializers_agree():
    analytical_model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(analytical_model.fs.unit)
    FlotationBankKineticCellBalanceAnalyticalInitializer().initialize(
        analytical_model.fs.unit
    )

    staged_model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(staged_model.fs.unit)
    FlotationBankKineticCellBalanceStagedInitializer().initialize(staged_model.fs.unit)

    a_unit = analytical_model.fs.unit
    s_unit = staged_model.fs.unit
    for component in COMPONENTS:
        assert value(s_unit.recovery[0, component]) == pytest.approx(
            value(a_unit.recovery[0, component]), abs=1e-6
        )
        for cell in a_unit.cells:
            assert value(s_unit.cell_solid_holdup[0, cell, component]) == pytest.approx(
                value(a_unit.cell_solid_holdup[0, cell, component]),
                rel=1e-5,
                abs=1e-9,
            )
            assert value(
                s_unit.cell_pulp_out_flow[0, cell, component]
            ) == pytest.approx(
                value(a_unit.cell_pulp_out_flow[0, cell, component]),
                rel=1e-5,
                abs=1e-9,
            )


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_kinetic_cell_balance_staged_initializer_rolls_back_on_solver_failure(
    monkeypatch,
):
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(model.fs.unit)
    unit = model.fs.unit

    outlet_constraint_names = (
        "concentrate_outlet_eq",
        "tails_outlet_eq",
        "recovery_eq",
    )

    # Treat the staged sub-block solve as non-convergent so the initializer
    # raises and runs its rollback. The real IPOPT solve still runs; only the
    # optimal-termination verdict is overridden.
    monkeypatch.setattr(
        initializer_module, "check_optimal_termination", lambda results: False
    )

    initializer = FlotationBankKineticCellBalanceStagedInitializer()
    with pytest.raises(InitializationError, match="did not converge"):
        initializer.initialize(unit)

    # A failed staged initialization must leave the model in its original
    # structural state: outlet equations reactivated and the bank-level outlet
    # Vars unfixed.
    for name in outlet_constraint_names:
        constraint = getattr(unit, name)
        assert all(constraint[index].active for index in constraint)
    for component in COMPONENTS:
        assert not unit.properties_concentrate[0].flow_mass_comp[component].fixed
        assert not unit.properties_tails[0].flow_mass_comp[component].fixed
        assert not unit.recovery[0, component].fixed


@pytest.mark.unit
def test_audit_summary_multicell_tracking():
    # Pure positive-k 2-cell walk engineered so cell 2's lower margin (1.02) is
    # within MARGIN_NEAR_BOUNDARY of 1.0. Exercises min-margin tracking and the
    # near-boundary counter on a successful walk.
    near_boundary = audit_cascade_solvability(
        {"A": 2.02}, [1.0, 1.0], {"A": 1.0}, "bank"
    )
    assert near_boundary["walk_completed"]
    assert near_boundary["min_lower_margin_no_carrier"] == pytest.approx(1.02)
    assert near_boundary["min_margin_cell"]["cell_index"] == 2
    assert near_boundary["near_boundary_cells"] >= 1
    assert near_boundary["cells_with_zero_k_carrier"] == 0

    # A walk with a zero-k carrier present in every cell exercises the zero-k
    # carrier counter (those cells report lower_margin=None).
    zero_k = audit_cascade_solvability(
        {"A": 10.0, "B": 5.0}, [1.0, 1.0], {"A": 0.5, "B": 0.0}, "bank"
    )
    assert zero_k["walk_completed"]
    assert zero_k["cells_with_zero_k_carrier"] == 2


@pytest.mark.unit
def test_cell_balance_performance_contents_flatten_component_rates():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit
    _fix_cell_balance_inputs(unit)
    unit.default_initializer().initialize(unit)

    variables = unit._get_performance_contents(time_point=0)["vars"]

    assert "k_cb" not in variables
    for component in COMPONENTS:
        assert variables[f"k_cb[{component}]"] == pytest.approx(0.5)
    for label in (
        "Apparent Inlet Slurry Flow",
        "Apparent Cell Residence Time",
        "Apparent Bank Residence Time",
        "Cell 1 Inventory Residence Time",
    ):
        assert label in variables
        # A NaN would fail self-equality.
        assert variables[label] == variables[label]


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_kinetic_cell_balance_unit_mass_conservation():
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    unit = model.fs.unit
    _fix_cell_balance_inputs(unit)
    unit.default_initializer().initialize(unit)
    assert_optimal_termination(solver.solve(model))

    for component in COMPONENTS:
        inlet = value(unit.inlet.flow_mass_comp[0, component])
        concentrate = value(unit.concentrate.flow_mass_comp[0, component])
        tails = value(unit.tails.flow_mass_comp[0, component])
        assert inlet == pytest.approx(concentrate + tails, abs=1e-6)


@pytest.mark.unit
def test_cell_balance_dispatcher_defaults_to_analytical_strategy():
    initializer = FlotationBankKineticCellBalanceInitializer()
    assert initializer.config.strategy == "analytical"
    # The unit's default_initializer is the dispatcher, not a concrete worker.
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    assert (
        model.fs.unit.default_initializer is FlotationBankKineticCellBalanceInitializer
    )


@pytest.mark.unit
def test_cell_balance_dispatcher_rejects_unknown_strategy():
    with pytest.raises(ValueError):
        FlotationBankKineticCellBalanceInitializer(strategy="bogus")


@pytest.mark.component
def test_cell_balance_dispatcher_default_matches_analytical_worker():
    dispatched = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(dispatched.fs.unit)
    dispatcher = FlotationBankKineticCellBalanceInitializer()
    dispatcher.initialize(dispatched.fs.unit)
    # The worker's summary is surfaced on the dispatcher.
    assert dispatcher.summary[dispatched.fs.unit]["status"] == InitializationStatus.Ok

    worker_model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(worker_model.fs.unit)
    FlotationBankKineticCellBalanceAnalyticalInitializer().initialize(
        worker_model.fs.unit
    )

    d_unit = dispatched.fs.unit
    w_unit = worker_model.fs.unit
    for component in COMPONENTS:
        assert value(d_unit.recovery[0, component]) == pytest.approx(
            value(w_unit.recovery[0, component]), abs=1e-9
        )


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.skipif(not solver.available(exception_flag=False), reason="No IPOPT")
def test_cell_balance_dispatcher_staged_strategy_delegates_to_staged(monkeypatch):
    model = _bank_model(
        fix_recoveries=False,
        recovery_basis="kinetic_cell_balance",
        number_of_cells=2,
    )
    _fix_cell_balance_inputs(model.fs.unit)

    def _fail_if_analytical_used(self, *args, **kwargs):
        raise AssertionError("staged strategy must not use the analytical worker")

    monkeypatch.setattr(
        FlotationBankKineticCellBalanceAnalyticalInitializer,
        "initialize",
        _fail_if_analytical_used,
    )

    dispatcher = FlotationBankKineticCellBalanceInitializer(strategy="staged")
    dispatcher.initialize(model.fs.unit)

    assert dispatcher.summary[model.fs.unit]["status"] == InitializationStatus.Ok
    assert_optimal_termination(solver.solve(model))
