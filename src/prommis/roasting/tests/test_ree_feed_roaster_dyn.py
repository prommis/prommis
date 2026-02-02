#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.common.collections import ComponentMap
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    SolverFactory,
    units,
    value,
    Var,
    TransformationFactory,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.scaling import get_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

import pytest
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)
import idaes.logger as idaeslog
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.roasting.ree_feed_properties import (
    ReeFeedParameters,
    ReeFeedPropertiesScaler,
)
from prommis.roasting.ree_roast_properties import (
    ReeRoastParameters,
    ReeRoastPropertiesScaler,
)
from prommis.roasting.ree_feed_roaster_dyn import (
    REEFeedRoaster,
    REEFeedRoasterScaler,
    REEFeedRoasterInitializer,
)


def get_model(
    dynamic=True,
    has_heat_transfer=True,
    has_pressure_change=True,
    time_set=None,
    nstep=None,
):
    m = ConcreteModel("REE_feed_roaster")
    if dynamic:
        if time_set is None:
            time_set = [0, 400]
        if nstep is None:
            nstep = 10
        m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=units.s)
    else:
        m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {"Ar", "O2", "H2O", "CO2", "N2", "SO2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )

    m.fs.ree_feed_prop = ReeFeedParameters(
        doc="ree feed property",
    )

    m.fs.ree_roast_prop = ReeRoastParameters(
        doc="ree product property",
    )

    m.fs.leach_solid_prop = CoalRefuseParameters(
        doc="leach solid property",
    )

    m.fs.roaster = REEFeedRoaster(
        dynamic=dynamic,
        gas_property_package=m.fs.prop_gas,
        solid_feed_property_package=m.fs.ree_feed_prop,
        solid_product_property_package=m.fs.ree_roast_prop,
        leach_solids_property_package=m.fs.leach_solid_prop,
        has_heat_transfer=has_heat_transfer,
        has_pressure_change=has_pressure_change,
    )

    if dynamic:
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=nstep, wrt=m.fs.time, scheme="BACKWARD")

    m.fs.roaster.volume.fix(4)
    m.fs.roaster.voidage.fix(0.4)
    if has_pressure_change:
        m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(900)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    gas_comp = {
        "Ar": 0.007,
        "O2": 0.155,
        "H2O": 0.0616,
        "CO2": 0.0235,
        "N2": 0.7528,
        "SO2": 0.0001,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[:, i].fix(v)
    # inlet flue gas mole flow rate
    m.fs.roaster.gas_inlet.flow_mol.fix(80)
    if has_heat_transfer:
        m.fs.roaster.heat_duty.fix(-3e6)
    # solid feed temperature
    m.fs.roaster.solid_inlet.temperature.fix(298.15)
    # limestone conversion
    m.fs.roaster.solid_inlet.flow_mass.fix(1030 / 3600)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "C"].fix(0.236806121)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "H"].fix(0.01639427)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "O"].fix(0.03278854)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "N"].fix(0.005464757)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "S"].fix(0.003643171)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "H2O"].fix(0.029126214)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Kaolinite"].fix(0.018784184)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Al2O3"].fix(0.182510737)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "SiO2"].fix(0.343148342)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "CaCO3"].fix(0.022581403)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "CaO"].fix(0.0029435831)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "FeS2"].fix(0.053569495)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Fe2O3"].fix(0.052239183)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Ree2X"].fix(2e-08)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Sc2X"].fix(1.57591e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Y2X"].fix(2.25028e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "La2X"].fix(3.63065e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Ce2X"].fix(7.70554e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Pr2X"].fix(1.44762e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Nd2X"].fix(3.40374e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Sm2X"].fix(7.80376e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Gd2X"].fix(6.21941e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Dy2X"].fix(3.77866e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Sc2O3"].fix(4.31618e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Y2O3"].fix(6.16319e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "La2O3"].fix(9.94383e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Ce2O3"].fix(2.11043e-05)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Pr2O3"].fix(3.96481e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Nd2O3"].fix(9.32235e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Sm2O3"].fix(2.13734e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Gd2O3"].fix(1.70341e-06)
    m.fs.roaster.solid_inlet.mass_frac_comp[:, "Dy2O3"].fix(1.03492e-06)

    if dynamic:
        m.fs.roaster.set_initial_condition()
    # Scale model
    scaler_config = {
        "zero_tolerance": 1e-18,
        "max_variable_scaling_factor": float("inf"),
        "min_variable_scaling_factor": 0,
        "max_constraint_scaling_factor": float("inf"),
        "min_constraint_scaling_factor": 0,
    }
    gas_scaler = m.fs.roaster.gas_in.default_scaler()
    gas_scaler.default_scaling_factors["flow_mol_phase"] = 1 / 80
    feed_scaler = ReeFeedPropertiesScaler(**scaler_config)
    prod_scaler = ReeRoastPropertiesScaler(**scaler_config)
    roaster_scaler = REEFeedRoasterScaler(**scaler_config)
    submodel_scalers = ComponentMap()
    submodel_scalers[m.fs.roaster.solid_in] = feed_scaler
    submodel_scalers[m.fs.roaster.solid_out] = prod_scaler
    submodel_scalers[m.fs.roaster.gas_in] = gas_scaler
    submodel_scalers[m.fs.roaster.gas_out] = gas_scaler
    roaster_scaler.scale_model(m.fs.roaster, submodel_scalers=submodel_scalers)

    if dynamic == False:
        solver = get_solver(options={"max_iter": 50})
        initializer = REEFeedRoasterInitializer()
        initializer.initialize(m.fs.roaster)
        result = solver.solve(m, tee=True)
        assert_optimal_termination(result)
    return m


@pytest.fixture(scope="module")
def model_steady_state():
    return get_model(dynamic=False, has_pressure_change=False)


@pytest.fixture(scope="module")
def model_dynamic():
    m_ss = get_model(dynamic=False)
    m_dyn = get_model(dynamic=True)
    copy_non_time_indexed_values(
        m_dyn.fs, m_ss.fs, copy_fixed=True, outlvl=idaeslog.ERROR
    )
    for t in m_dyn.fs.time:
        copy_values_at_time(
            m_dyn.fs,
            m_ss.fs,
            t,
            0.0,
            copy_fixed=True,
            outlvl=idaeslog.ERROR,
        )
    return m_dyn


@pytest.mark.unit
def test_build_steady_state(model_steady_state):
    m = model_steady_state
    assert hasattr(m.fs, "roaster")
    assert isinstance(m.fs.roaster, REEFeedRoaster)
    assert len(m.fs.roaster.config) == 12
    assert not m.fs.roaster.config.dynamic
    assert m.fs.roaster.config.has_holdup
    assert m.fs.roaster.config.has_heat_transfer
    assert not m.fs.roaster.config.has_pressure_change
    assert m.fs.roaster.config.gas_property_package is m.fs.prop_gas
    assert m.fs.roaster.config.solid_feed_property_package is m.fs.ree_feed_prop
    assert m.fs.roaster.config.solid_product_property_package is m.fs.ree_roast_prop
    assert m.fs.roaster.config.leach_solids_property_package is m.fs.leach_solid_prop
    assert len(m.fs.prop_gas.component_list) == 6
    assert isinstance(m.fs.roaster.heat_duty, Var)
    assert isinstance(m.fs.roaster.volume, Var)
    assert isinstance(m.fs.roaster.voidage, Var)
    assert isinstance(m.fs.roaster.solid_mass_balance_eqn, Constraint)
    assert isinstance(m.fs.roaster.gas_mass_balance_eqn, Constraint)
    assert isinstance(m.fs.roaster.energy_balance_eqn, Constraint)
    assert len(m.fs.roaster.solid_mass_balance_eqn) == 24
    assert len(m.fs.roaster.gas_mass_balance_eqn) == 6
    assert number_variables(m.fs.roaster) == 270
    assert number_total_constraints(m.fs.roaster) == 224
    assert number_unused_variables(m.fs.roaster) == 0
    assert_units_consistent(m.fs.roaster)


@pytest.mark.unit
def test_structural_issues_steady_state(model_steady_state):
    dt = DiagnosticsToolbox(model_steady_state)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solve_steady_state(model_steady_state):
    solver = get_solver()
    results = solver.solve(model_steady_state, tee=False)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues_steady_state(model_steady_state):
    dt = DiagnosticsToolbox(model_steady_state)
    dt.assert_no_numerical_warnings()
    dt.report_numerical_issues()


@pytest.mark.component
@pytest.mark.solver
def test_solution_steady_state(model_steady_state):
    m = model_steady_state
    assert value(m.fs.roaster.gas_out[0].flow_mol) == pytest.approx(
        81.98495217241386, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[0].temperature) == pytest.approx(
        784.0991061479093, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[0].mole_frac_comp["O2"]) == pytest.approx(
        0.0671451071666029, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[0].mole_frac_comp["H2O"]) == pytest.approx(
        0.09464375438304656, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[0].mole_frac_comp["SO2"]) == pytest.approx(
        0.00361051400612864, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].flow_mass) == pytest.approx(
        0.1847323, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["CaCO3"]) == pytest.approx(
        0.000316611070408, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["CaO"]) == pytest.approx(
        0.023976930120916822, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["Fe2O3"]) == pytest.approx(
        0.1361225, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["Ree2X"]) == pytest.approx(
        6.208683e-05, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["Sc2X"]) == pytest.approx(
        3.81932e-6, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[0].mass_frac_comp["Sc2O3"]) == pytest.approx(
        2.2049e-5, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.leach_solid_out[0].flow_mass) == pytest.approx(
        665.03639, rel=1e-5, abs=1e-6
    )
    assert value(
        m.fs.roaster.leach_solid_out[0].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2049e-5, rel=1e-5, abs=1e-6)
    assert value(
        m.fs.roaster.leach_solid_out[0].mass_frac_comp["inerts"]
    ) == pytest.approx(0.5454315, rel=1e-5, abs=1e-6)


@pytest.mark.unit
def test_build_dynamic(model_dynamic):
    m = model_dynamic
    assert hasattr(m.fs, "roaster")
    assert isinstance(m.fs.roaster, REEFeedRoaster)
    assert len(m.fs.roaster.config) == 12
    assert m.fs.roaster.config.dynamic
    assert m.fs.roaster.config.has_holdup
    assert m.fs.roaster.config.has_heat_transfer
    assert m.fs.roaster.config.has_pressure_change
    assert m.fs.roaster.config.gas_property_package is m.fs.prop_gas
    assert m.fs.roaster.config.solid_feed_property_package is m.fs.ree_feed_prop
    assert m.fs.roaster.config.solid_product_property_package is m.fs.ree_roast_prop
    assert m.fs.roaster.config.leach_solids_property_package is m.fs.leach_solid_prop
    assert len(m.fs.prop_gas.component_list) == 6
    assert isinstance(m.fs.roaster.heat_duty, Var)
    assert isinstance(m.fs.roaster.deltaP, Var)
    assert isinstance(m.fs.roaster.volume, Var)
    assert isinstance(m.fs.roaster.voidage, Var)
    assert isinstance(m.fs.roaster.solid_mass_balance_eqn, Constraint)
    assert isinstance(m.fs.roaster.gas_mass_balance_eqn, Constraint)
    assert isinstance(m.fs.roaster.energy_balance_eqn, Constraint)
    assert len(m.fs.roaster.solid_mass_balance_eqn) == 264
    assert len(m.fs.roaster.gas_mass_balance_eqn) == 66
    assert number_variables(m.fs.roaster) == 3247
    assert number_total_constraints(m.fs.roaster) == 2725
    assert number_unused_variables(m.fs.roaster) == 0


# Currently there is a structure warning related to unit of measure
# Structure assertion test is commented out
# @pytest.mark.unit
# def test_structural_issues_dynamic(model_dynamic):
#    dt = DiagnosticsToolbox(model_dynamic)
#    The unit consistency is not checked since Pyomo.Dae has not fixed the bug for time domain unit
#    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve_dynamic(model_dynamic):
    m = model_dynamic
    solver = get_solver()
    result = solver.solve(m, tee=False)
    assert_optimal_termination(result)
    # Add disturbance and solve dynamic model
    temp0 = m.fs.roaster.solid_inlet.temperature[0].value
    # solid inlet temperature ramp rate for unscaled model
    dTdt = 0.2
    for t in m.fs.time:
        if t > 30:
            m.fs.roaster.solid_inlet.temperature[t].fix(temp0 + dTdt * t)
    result = solver.solve(m, tee=False)
    assert_optimal_termination(result)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues_dynamic(model_dynamic):
    dt = DiagnosticsToolbox(model_dynamic)
    dt.assert_no_numerical_warnings()
    dt.report_numerical_issues()


@pytest.mark.component
@pytest.mark.solver
def test_solution_dynamic(model_dynamic):
    m = model_dynamic
    t = 400
    assert value(m.fs.roaster.gas_out[t].flow_mol) == pytest.approx(
        82.03953211409154, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[t].temperature) == pytest.approx(
        820.1335172303061, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[t].mole_frac_comp["O2"]) == pytest.approx(
        0.06710043631172598, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[t].mole_frac_comp["H2O"]) == pytest.approx(
        0.09458078899963354, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.gas_out[t].mole_frac_comp["SO2"]) == pytest.approx(
        0.00360817, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].flow_mass) == pytest.approx(
        0.1823290484141596, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["CaCO3"]) == pytest.approx(
        6.465804386318614e-05, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["CaO"]) == pytest.approx(
        0.024120773737993473, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["Fe2O3"]) == pytest.approx(
        0.1361374, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["Ree2X"]) == pytest.approx(
        6.233819e-05, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["Sc2X"]) == pytest.approx(
        3.74178e-6, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.solid_out[t].mass_frac_comp["Sc2O3"]) == pytest.approx(
        2.21087888e-5, rel=1e-5, abs=1e-6
    )
    assert value(m.fs.roaster.leach_solid_out[t].flow_mass) == pytest.approx(
        656.3845742909085, rel=1e-5, abs=1e-6
    )
    assert value(
        m.fs.roaster.leach_solid_out[t].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2108788e-5, rel=1e-5, abs=1e-6)
    assert value(
        m.fs.roaster.leach_solid_out[t].mass_frac_comp["inerts"]
    ) == pytest.approx(0.54524178, rel=1e-5, abs=1e-6)
