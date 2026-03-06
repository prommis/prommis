#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for the SettlerTank unit model.

Author: Douglas Allan

Portions of this file were created with the help of Google Gemini 2.5 Pro.
"""

import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ComponentMap,
    ConcreteModel,
    Constraint,
    units,
    Var,
)
from pyomo.network import Port
from idaes.core import (
    EnergyBalanceType,
    FlowsheetBlock,
)
from idaes.core.scaling.util import jacobian_cond
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.settler_tank import (
    SettlerTank,
    SettlerTankInitializer,
    SettlerTankScaler,
)
from prommis.util import assert_solution_equivalent

__author__ = "Douglas Allan"


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.dehpa_kerosene_params = REESolExOgParameters()
    m.fs.sulfate_leaching_params = SulfuricAcidLeachingParameters()

    m.fs.unit = SettlerTank(
        light_phase_alias="organic",
        heavy_phase_alias="aqueous",
        light_phase_config={
            "property_package": m.fs.dehpa_kerosene_params,
            "energy_balance_type": EnergyBalanceType.isothermal,
            "has_pressure_balance": False,
        },
        heavy_phase_config={
            "property_package": m.fs.sulfate_leaching_params,
            "energy_balance_type": EnergyBalanceType.isothermal,
            "has_pressure_balance": False,
        },
        has_holdup=True,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=4,
    )

    # Specify aqueous inlet
    m.fs.unit.aqueous_inlet.flow_vol.fix(50 * units.L / units.hr)
    m.fs.unit.aqueous_inlet.temperature.fix(300 * units.K)
    m.fs.unit.aqueous_inlet.pressure.fix(1e5 * units.Pa)

    m.fs.unit.aqueous_inlet.conc_mass_comp.fix(1e-16 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H"].fix(10.75 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(100 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e4 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Y"].fix(8.89 * units.mg / units.L)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(1752.34 * units.mg / units.L)

    # Specify organic inlet
    m.fs.unit.organic_inlet.flow_vol.fix(50 * units.L / units.hr)
    m.fs.unit.organic_inlet.temperature.fix(300 * units.K)
    m.fs.unit.organic_inlet.pressure.fix(1e5 * units.Pa)

    m.fs.unit.organic_inlet.conc_mass_comp.fix(1e-16 * units.kg / units.m**3)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(
        820e3 * units.mg / units.L
    )
    m.fs.unit.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
        0.05 * 975.8e3 * units.mg / units.L
    )
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(20 * units.mg / units.L)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Y_o"].fix(10 * units.mg / units.L)

    # Fix geometric and operational parameters
    m.fs.unit.length.fix(1 * units.m)
    m.fs.unit.settler_width.fix(0.5 * units.m)
    m.fs.unit.light_phase_weir_height.fix(0.5 * units.m)
    m.fs.unit.light_phase_weir_channel_width.fix(0.5 * units.m)

    # The heavy phase height is unspecified by the steady state model,
    # therefore fix it.
    m.fs.unit.heavy_phase_height.fix(0.25 * units.m)

    aqueous_scaler = m.fs.unit.heavy_phase.properties.default_scaler()
    aqueous_scaler.default_scaling_factors["flow_vol"] = 1 / 50
    aqueous_scaler.default_scaling_factors["conc_mass_comp[H2O]"] = 1e-6
    aqueous_scaler.default_scaling_factors["conc_mass_comp[H]"] = 1 / 50
    aqueous_scaler.default_scaling_factors["conc_mass_comp[SO4]"] = 1e-3
    aqueous_scaler.default_scaling_factors["conc_mass_comp[HSO4]"] = 3e-4
    aqueous_scaler.default_scaling_factors["conc_mass_comp[Fe]"] = 1e-3
    aqueous_scaler.default_scaling_factors["conc_mass_comp[Y]"] = 1e-1
    organic_scaler = m.fs.unit.light_phase.properties.default_scaler()
    organic_scaler.default_scaling_factors["flow_vol"] = 1 / 50

    submodel_scalers = ComponentMap()
    submodel_scalers[m.fs.unit.aqueous_phase.properties] = aqueous_scaler
    submodel_scalers[m.fs.unit.organic_phase.properties] = organic_scaler

    scaler_obj = SettlerTankScaler()
    scaler_obj.scale_model(m.fs.unit, submodel_scalers=submodel_scalers)

    return m


@pytest.mark.unit
def test_model_construction(model):
    # Number of time points, in this case 1
    nt = 1
    assert isinstance(model.fs.unit, SettlerTank)
    assert model.fs.unit.default_initializer is SettlerTankInitializer
    assert model.fs.unit.default_scaler is SettlerTankScaler
    for port in ["inlet", "outlet"]:
        for name in ["organic", "aqueous"]:
            port_obj = getattr(model.fs.unit, f"{name}_{port}")
            assert isinstance(port_obj, Port)

    assert model.fs.unit.light_phase is model.fs.unit.organic_phase
    assert model.fs.unit.heavy_phase is model.fs.unit.aqueous_phase

    assert isinstance(model.fs.unit.length, Var)
    assert len(model.fs.unit.length) == 1
    assert model.fs.unit.length is model.fs.unit.light_phase.length
    assert isinstance(model.fs.unit.settler_width, Var)
    assert len(model.fs.unit.settler_width) == 1
    assert isinstance(model.fs.unit.light_phase_weir_height, Var)
    assert len(model.fs.unit.light_phase_weir_height) == 1
    assert isinstance(model.fs.unit.heavy_phase_height, Var)
    assert len(model.fs.unit.heavy_phase_height) == nt
    assert isinstance(model.fs.unit.light_phase_height, Var)
    assert len(model.fs.unit.light_phase_height) == nt
    assert isinstance(model.fs.unit.weir_flow_eqn, Constraint)
    assert len(model.fs.unit.weir_flow_eqn) == nt

    assert number_variables(model) == 917
    assert number_total_constraints(model) == 818
    assert number_unused_variables(model) == 59


@pytest.mark.unit
def test_no_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    # smooth_max(H, 0) should never be negative, but interval
    # arithmetic does not capture that
    dt.assert_no_structural_warnings(ignore_evaluation_errors=True)


@pytest.mark.solver
@pytest.mark.component
def test_initialize(model):
    # Use the model's default initializer directly
    initializer = model.fs.unit.default_initializer()
    initializer.initialize(model.fs.unit)


@pytest.mark.solver
@pytest.mark.component
def test_solve(model):
    solver = get_solver(
        "ipopt_v2", solver_options={"bound_relax_factor": 0, "constr_viol_tol": 1e-10}
    )
    results = solver.solve(model, tee=False)
    assert_optimal_termination(results)


@pytest.mark.unit
def test_no_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    # Species with zero concentration are somehow ending up with values
    # of ~1e-21 when the have a lower bound of 1e-20.
    for vardata in model.fs.unit.heavy_phase.properties[
        0.0, 0.0
    ].conc_mol_comp.values():
        if vardata.value <= vardata.lb:
            vardata.value = vardata.lb * 1.1

    dt.assert_no_numerical_warnings()

    assert jacobian_cond(model, scaled=False) == pytest.approx(6.2849844e12)
    assert jacobian_cond(model, scaled=True) == pytest.approx(16447.48)


@pytest.mark.solver
@pytest.mark.component
def test_solution(model):
    t = 0
    _abs = 1e-8
    _rel = 1e-6
    expected_results = {
        "organic_outlet.flow_vol": {0: (50, _rel, None)},
        "organic_outlet.temperature": {0: (300, _rel, None)},
        # Presently outlet pressure is undefined
        # "organic_outlet.pressure": {
        #     0: (1e5, _rel, None)
        # },
        "aqueous_outlet.flow_vol": {
            0: (50, _rel, None),
        },
        "aqueous_outlet.temperature": {0: (300, _rel, None)},
        # Presently outlet pressure is undefined
        # "aqueous_outlet.pressure": {
        #     0: (1e5, _rel, None)
        # },
        "length": {None: (1, _rel, None)},
        "settler_width": {None: (0.5, _rel, None)},
        "light_phase_weir_height": {None: (0.5, _rel, None)},
        "light_phase_height": {0: (0.2506140, _rel, None)},
        "heavy_phase_height": {0: (0.25, _rel, None)},
    }
    # Organic outlet concentrations
    expected_results["organic_outlet.conc_mass_comp"] = {}
    for idx in model.fs.unit.organic_outlet.conc_mass_comp:
        expected_results["organic_outlet.conc_mass_comp"][idx] = (1e-16, None, _abs)
    expected_results["organic_outlet.conc_mass_comp"][(0, "Kerosene")] = (
        820e3,
        _rel,
        None,
    )
    expected_results["organic_outlet.conc_mass_comp"][(0, "DEHPA")] = (
        0.05 * 975.8e3,
        _rel,
        None,
    )
    expected_results["organic_outlet.conc_mass_comp"][(0, "Fe_o")] = (20, _rel, None)
    expected_results["organic_outlet.conc_mass_comp"][(0, "Y_o")] = (10, _rel, None)

    # Aqueous outlet concentrations
    expected_results["aqueous_outlet.conc_mass_comp"] = {}
    for idx in model.fs.unit.aqueous_outlet.conc_mass_comp:
        expected_results["aqueous_outlet.conc_mass_comp"][idx] = (1e-16, None, _abs)
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "H2O")] = (1e6, _rel, None)
    # The inlet concentrations of sulfate is not at equilibrium---therefore
    # these values are different than those at the inlet.
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "H")] = (
        3.3870004e01,
        _rel,
        None,
    )
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "SO4")] = (
        2.3195203e03,
        _rel,
        None,
    )
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "HSO4")] = (
        7.7573597e03,
        _rel,
        None,
    )
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "Y")] = (8.89, _rel, None)
    expected_results["aqueous_outlet.conc_mass_comp"][(0, "Fe")] = (1752.34, _rel, None)

    expected_results["organic_phase.area"] = {}
    for idx in model.fs.unit.organic_phase.area:
        expected_results["organic_phase.area"][idx] = (0.1253070, _rel, None)

    expected_results["aqueous_phase.area"] = {}
    for idx in model.fs.unit.aqueous_phase.area:
        expected_results["aqueous_phase.area"][idx] = (0.125, _rel, None)

    assert_solution_equivalent(model.fs.unit, expected_results=expected_results)
