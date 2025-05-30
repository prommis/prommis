#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
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
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
import idaes.core.util.scaling as iscale
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values

import pytest
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)
import idaes.logger as idaeslog
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.roasting.ree_feed_roaster_dyn import REEFeedRoaster
from prommis.roasting.ree_feed_properties import ReeFeedParameters
from prommis.roasting.ree_roast_properties import ReeRoastParameters


def get_model(dynamic=True, time_set=None, nstep=None):
    m = ConcreteModel("REE_feed_roaster")
    if dynamic:
        if time_set is None:
            time_set = [0, 400]
        if nstep is None:
            nstep = 10
        m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=units.s)
    else:
        m.fs = FlowsheetBlock(dynamic=False)
    gas_species = {"O2", "H2O", "CO2", "N2", "SO2"}
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
    m.fs.prop_gas.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.prop_gas.set_default_scaling("pressure", 1e-5)
    m.fs.prop_gas.set_default_scaling("temperature", 1e-2)
    m.fs.prop_gas.set_default_scaling("flow_mol", 1e1)
    m.fs.prop_gas.set_default_scaling("flow_mol_phase", 1e1)
    m.fs.prop_gas.set_default_scaling("_energy_density_term", 1e-4)
    m.fs.prop_gas.set_default_scaling("phase_frac", 1)
    _mf_scale = {
        "O2": 5,
        "CO2": 10,
        "H2O": 5,
        "N2": 1,
        "SO2": 10,
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )
    _massfrac_scale_1 = {
        "Ree2X": 1e4,
        "Sc2X": 1e4,
        "Y2X": 1e4,
        "La2X": 1e4,
        "Ce2X": 1e4,
        "Pr2X": 1e4,
        "Nd2X": 1e4,
        "Sm2X": 1e4,
        "Gd2X": 1e4,
        "Dy2X": 1e4,
        "Sc2O3": 1e4,
        "Y2O3": 1e4,
        "La2O3": 1e4,
        "Ce2O3": 1e4,
        "Pr2O3": 1e4,
        "Nd2O3": 1e4,
        "Sm2O3": 1e4,
        "Gd2O3": 1e4,
        "Dy2O3": 1e4,
    }
    _massfrac_scale_2 = {
        "Sc2O3": 1e4,
        "Y2O3": 1e4,
        "La2O3": 1e4,
        "Ce2O3": 1e4,
        "Pr2O3": 1e4,
        "Nd2O3": 1e4,
        "Sm2O3": 1e4,
        "Gd2O3": 1e4,
        "Dy2O3": 1e4,
    }
    for comp, s in _massfrac_scale_1.items():
        m.fs.ree_feed_prop.set_default_scaling("mass_frac_comp", s, index=comp)
        m.fs.ree_roast_prop.set_default_scaling("mass_frac_comp", s, index=comp)
        m.fs.ree_feed_prop.set_default_scaling("flow_mol_comp", s, index=comp)
        m.fs.ree_roast_prop.set_default_scaling("flow_mol_comp", s, index=comp)
    for comp, s in _massfrac_scale_2.items():
        m.fs.leach_solid_prop.set_default_scaling("mass_frac_comp", s, index=comp)
    m.fs.ree_feed_prop.set_default_scaling("temperature", 1e-2)
    m.fs.ree_roast_prop.set_default_scaling("temperature", 1e-2)
    m.fs.ree_feed_prop.set_default_scaling("enth_mol", 1e-6)
    m.fs.ree_roast_prop.set_default_scaling("enth_mol", 1e-6)
    m.fs.ree_feed_prop.set_default_scaling("enth_mass", 1e-7)
    m.fs.ree_roast_prop.set_default_scaling("enth_mass", 1e-7)

    m.fs.roaster = REEFeedRoaster(
        dynamic=dynamic,
        gas_property_package=m.fs.prop_gas,
        solid_feed_property_package=m.fs.ree_feed_prop,
        solid_product_property_package=m.fs.ree_roast_prop,
        leach_solids_property_package=m.fs.leach_solid_prop,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    if dynamic:
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=nstep, wrt=m.fs.time, scheme="BACKWARD")

    # input data for reactor
    m.fs.roaster.volume.fix(4)
    m.fs.roaster.voidage.fix(0.4)
    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.heat_duty.fix(-3e6)
    # input data for gas inlet stream
    m.fs.roaster.gas_inlet.temperature.fix(900)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    m.fs.roaster.gas_inlet.flow_mol.fix(80)
    gas_comp = {
        "O2": 0.155,
        "H2O": 0.0616,
        "CO2": 0.0235,
        "N2": 0.7598,
        "SO2": 0.0001,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[:, i].fix(v)
    # input data for solid inlet stream
    m.fs.roaster.solid_inlet.temperature.fix(298.15)
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
    iscale.calculate_scaling_factors(m)
    if dynamic:
        m.fs.roaster.set_initial_condition()
    else:
        m.fs.roaster.initialize(outlvl=4)
    return m


@pytest.fixture(scope="module")
def model_steady_state():
    m = get_model(dynamic=False)
    return m


@pytest.fixture(scope="module")
def model_dynamic():
    m_ss = get_model(dynamic=False)
    m_dyn = get_model(dynamic=True)
    copy_non_time_indexed_values(
        m_dyn.fs, m_ss.fs, copy_fixed=True, outlvl=idaeslog.ERROR
    )
    for t in m_dyn.fs.time:
        copy_values_at_time(
            m_dyn.fs, m_ss.fs, t, 0.0, copy_fixed=True, outlvl=idaeslog.ERROR
        )
    return m_dyn


@pytest.mark.unit
def test_build_steady_state(model_steady_state):
    assert hasattr(model_steady_state.fs, "roaster")
    assert isinstance(model_steady_state.fs.roaster, REEFeedRoaster)
    assert len(model_steady_state.fs.roaster.config) == 12
    assert not model_steady_state.fs.roaster.config.dynamic
    assert model_steady_state.fs.roaster.config.has_holdup
    assert model_steady_state.fs.roaster.config.has_heat_transfer
    assert model_steady_state.fs.roaster.config.has_pressure_change
    assert (
        model_steady_state.fs.roaster.config.gas_property_package
        is model_steady_state.fs.prop_gas
    )
    assert (
        model_steady_state.fs.roaster.config.solid_feed_property_package
        is model_steady_state.fs.ree_feed_prop
    )
    assert (
        model_steady_state.fs.roaster.config.solid_product_property_package
        is model_steady_state.fs.ree_roast_prop
    )
    assert (
        model_steady_state.fs.roaster.config.leach_solids_property_package
        is model_steady_state.fs.leach_solid_prop
    )
    assert len(model_steady_state.fs.prop_gas.component_list) == 5
    assert isinstance(model_steady_state.fs.roaster.heat_duty, Var)
    assert isinstance(model_steady_state.fs.roaster.deltaP, Var)
    assert isinstance(model_steady_state.fs.roaster.volume, Var)
    assert isinstance(model_steady_state.fs.roaster.voidage, Var)
    assert isinstance(model_steady_state.fs.roaster.solid_mass_balance_eqn, Constraint)
    assert isinstance(model_steady_state.fs.roaster.gas_mass_balance_eqn, Constraint)
    assert isinstance(model_steady_state.fs.roaster.energy_balance_eqn, Constraint)
    assert len(model_steady_state.fs.roaster.solid_mass_balance_eqn) == 24
    assert len(model_steady_state.fs.roaster.gas_mass_balance_eqn) == 5
    assert number_variables(model_steady_state.fs.roaster) == 280
    assert number_total_constraints(model_steady_state.fs.roaster) == 234
    assert number_unused_variables(model_steady_state.fs.roaster) == 0
    assert_units_consistent(model_steady_state.fs.roaster)


# Currently there is a structure warning related to sqrt function in a constraint
# Structure assertion test is commented out
# @pytest.mark.unit
# def test_structural_issues_steady_state(model_steady_state):
#    dt = DiagnosticsToolbox(model_steady_state)
#    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve_steady_state(model_steady_state):
    model_steady_state.fs.roaster.initialize(outlvl=4)
    solver = SolverFactory("ipopt")
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
    assert value(model_steady_state.fs.roaster.gas_out[0].flow_mol) == pytest.approx(
        81.98495217241386, rel=1e-5, abs=1e-6
    )
    assert value(model_steady_state.fs.roaster.gas_out[0].temperature) == pytest.approx(
        784.33228, rel=1e-5, abs=1e-6
    )
    assert value(
        model_steady_state.fs.roaster.gas_out[0].mole_frac_comp["O2"]
    ) == pytest.approx(0.0671451071666029, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.gas_out[0].mole_frac_comp["H2O"]
    ) == pytest.approx(0.09464375438304656, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.gas_out[0].mole_frac_comp["SO2"]
    ) == pytest.approx(0.00361051400612864, rel=1e-5, abs=1e-6)
    assert value(model_steady_state.fs.roaster.solid_out[0].flow_mass) == pytest.approx(
        0.1847323, rel=1e-5, abs=1e-6
    )
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["CaCO3"]
    ) == pytest.approx(0.0003107265, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["CaO"]
    ) == pytest.approx(0.0239803, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["Fe2O3"]
    ) == pytest.approx(0.1361225, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["Ree2X"]
    ) == pytest.approx(6.208683e-05, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["Sc2X"]
    ) == pytest.approx(3.81932e-6, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.solid_out[0].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2049e-5, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.leach_solid_out[0].flow_mass
    ) == pytest.approx(665.03639, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.leach_solid_out[0].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2049e-5, rel=1e-5, abs=1e-6)
    assert value(
        model_steady_state.fs.roaster.leach_solid_out[0].mass_frac_comp["inerts"]
    ) == pytest.approx(0.5454315, rel=1e-5, abs=1e-6)


@pytest.mark.unit
def test_build_dynamic(model_dynamic):
    assert hasattr(model_dynamic.fs, "roaster")
    assert isinstance(model_dynamic.fs.roaster, REEFeedRoaster)
    assert len(model_dynamic.fs.roaster.config) == 12
    assert model_dynamic.fs.roaster.config.dynamic
    assert model_dynamic.fs.roaster.config.has_holdup
    assert model_dynamic.fs.roaster.config.has_heat_transfer
    assert model_dynamic.fs.roaster.config.has_pressure_change
    assert (
        model_dynamic.fs.roaster.config.gas_property_package
        is model_dynamic.fs.prop_gas
    )
    assert (
        model_dynamic.fs.roaster.config.solid_feed_property_package
        is model_dynamic.fs.ree_feed_prop
    )
    assert (
        model_dynamic.fs.roaster.config.solid_product_property_package
        is model_dynamic.fs.ree_roast_prop
    )
    assert (
        model_dynamic.fs.roaster.config.leach_solids_property_package
        is model_dynamic.fs.leach_solid_prop
    )
    assert len(model_dynamic.fs.prop_gas.component_list) == 5
    assert isinstance(model_dynamic.fs.roaster.heat_duty, Var)
    assert isinstance(model_dynamic.fs.roaster.deltaP, Var)
    assert isinstance(model_dynamic.fs.roaster.volume, Var)
    assert isinstance(model_dynamic.fs.roaster.voidage, Var)
    assert isinstance(model_dynamic.fs.roaster.solid_mass_balance_eqn, Constraint)
    assert isinstance(model_dynamic.fs.roaster.gas_mass_balance_eqn, Constraint)
    assert isinstance(model_dynamic.fs.roaster.energy_balance_eqn, Constraint)
    assert len(model_dynamic.fs.roaster.solid_mass_balance_eqn) == 264
    assert len(model_dynamic.fs.roaster.gas_mass_balance_eqn) == 55
    assert number_variables(model_dynamic.fs.roaster) == 3346
    assert number_total_constraints(model_dynamic.fs.roaster) == 2835
    assert number_unused_variables(model_dynamic.fs.roaster) == 0


# Currently there is a structure warning related to sqrt function in a constraint
# Structure assertion test is commented out
# @pytest.mark.unit
# def test_structural_issues_dynamic(model_dynamic):
#    dt = DiagnosticsToolbox(model_dynamic)
#    The unit consistency is not checked since Pyomo.Dae has not fixed the bug for time domain unit
#    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve_dynamic(model_dynamic):
    solver = SolverFactory("ipopt")
    result1 = solver.solve(model_dynamic, tee=False)
    assert_optimal_termination(result1)
    # Add disturbance and solve dynamic model
    temp0 = model_dynamic.fs.roaster.gas_inlet.temperature[0].value
    dTdt = 0.2
    for t in model_dynamic.fs.time:
        if t > 30:
            model_dynamic.fs.roaster.gas_inlet.temperature[t].fix(temp0 + dTdt * t)
    result = solver.solve(model_dynamic, tee=False)
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
    t = 400
    assert value(model_dynamic.fs.roaster.gas_out[t].flow_mol) == pytest.approx(
        81.99881, rel=1e-5, abs=1e-6
    )
    assert value(model_dynamic.fs.roaster.gas_out[t].temperature) == pytest.approx(
        791.2888, rel=1e-5, abs=1e-6
    )
    assert value(
        model_dynamic.fs.roaster.gas_out[t].mole_frac_comp["O2"]
    ) == pytest.approx(0.06713376, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.gas_out[t].mole_frac_comp["H2O"]
    ) == pytest.approx(0.09462776, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.gas_out[t].mole_frac_comp["SO2"]
    ) == pytest.approx(0.0036099, rel=1e-5, abs=1e-6)
    assert value(model_dynamic.fs.roaster.solid_out[t].flow_mass) == pytest.approx(
        0.184122, rel=1e-5, abs=1e-6
    )
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["CaCO3"]
    ) == pytest.approx(0.000261485, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["CaO"]
    ) == pytest.approx(0.0240084, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["Fe2O3"]
    ) == pytest.approx(0.136125477, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["Ree2X"]
    ) == pytest.approx(6.2129533e-05, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["Sc2X"]
    ) == pytest.approx(3.8060523e-6, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.solid_out[t].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2059358e-5, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.leach_solid_out[t].flow_mass
    ) == pytest.approx(662.8404157038169, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.leach_solid_out[t].mass_frac_comp["Sc2O3"]
    ) == pytest.approx(2.2059358e-5, rel=1e-5, abs=1e-6)
    assert value(
        model_dynamic.fs.roaster.leach_solid_out[t].mass_frac_comp["inerts"]
    ) == pytest.approx(0.545394, rel=1e-5, abs=1e-6)
