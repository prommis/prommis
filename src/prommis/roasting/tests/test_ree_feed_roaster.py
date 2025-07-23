#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    Var,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

import pytest

from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.roasting.ree_feed_roaster import REEFeedRoaster


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_solid = CoalRefuseParameters(
        doc="solid property",
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
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )

    m.fs.roaster = REEFeedRoaster(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
        ],
    )

    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(895)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    gas_comp = {
        "O2": 0.155,
        "H2O": 0.0617,
        "CO2": 0.0235,
        "N2": 0.7598,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    # inlet flue gas mole flow rate
    m.fs.roaster.gas_inlet.flow_mol.fix(80)

    # fix outlet temperature at 600 C
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature
    m.fs.roaster.temp_feed.fix(298.15)
    # limestone conversion
    m.fs.roaster.xconv_caco3.fix(0.5)

    # solid feed rate at 1030 kg/hr
    m.fs.roaster.flow_mass_feed.fix(1030 / 3600)
    # surface moisture in solid feed
    m.fs.roaster.mass_frac_moist_feed.fix(0.029126)
    # combustible organic material
    # Note that the heat of combustion will cause very large heat release
    m.fs.roaster.mass_frac_organic_feed.fix(0.294854)

    # impurity minerals modeled plus some small amount of unknown mineral that is modeled as Un2O3 here
    # note that the unknown in UK's excel sheet also include O, C, and S elements in minerals
    # assume all Al2O3 is in Kaolinite form
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Un2O3"].fix(0.00793)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "CaCO3"].fix(0.041217)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "SiO2"].fix(0.143882)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Al2O3"].fix(0.0)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Kaolinite"].fix(0.611613)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Pyrite"].fix(0.195359)

    # organic composition assumed as a typical coal
    m.fs.roaster.mass_frac_comp_organic_feed[0, "C"].fix(0.804)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "H"].fix(0.055)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "O"].fix(0.111)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "N"].fix(0.018)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "S"].fix(0.012)

    # Insoluble REE based on UK's excel file
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Sc"].fix(10.5836902)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Y"].fix(18.25358553)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "La"].fix(31.88774016)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Ce"].fix(67.77169805)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Pr"].fix(12.74306914)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Nd"].fix(30.06183493)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Sm"].fix(6.93187974)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Eu"].fix(0.902019051)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Gd"].fix(5.559717425)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Tb"].fix(0.492010392)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Dy"].fix(3.394871702)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Ho"].fix(0.984020783)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Er"].fix(3.044997646)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Tm"].fix(0.69428133)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Yb"].fix(2.20857998)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Lu"].fix(0.579478906)

    # Dissovable REE based on UK's excel file
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Sc"].fix(2.896677798)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Y"].fix(4.995871471)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "La"].fix(8.727438841)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Ce"].fix(18.54861295)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Pr"].fix(3.487683857)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Nd"].fix(8.227702072)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Sm"].fix(1.89720426)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Eu"].fix(0.246875949)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Gd"].fix(1.521653575)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Tb"].fix(0.134659608)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Dy"].fix(0.929151298)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Ho"].fix(0.269319217)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Er"].fix(0.833393354)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Tm"].fix(0.19001967)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Yb"].fix(0.60447202)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Lu"].fix(0.158599094)

    # conversion of insoluble REE to dissovable REE
    m.fs.roaster.xconv_comp_ins[0, "Sc"].fix(0.850340479)
    m.fs.roaster.xconv_comp_ins[0, "Y"].fix(0.737455969)
    m.fs.roaster.xconv_comp_ins[0, "La"].fix(0.784164924)
    m.fs.roaster.xconv_comp_ins[0, "Ce"].fix(0.706035503)
    m.fs.roaster.xconv_comp_ins[0, "Pr"].fix(0.478582677)
    m.fs.roaster.xconv_comp_ins[0, "Nd"].fix(0.674653216)
    m.fs.roaster.xconv_comp_ins[0, "Sm"].fix(0.678933333)
    m.fs.roaster.xconv_comp_ins[0, "Eu"].fix(0.1)
    m.fs.roaster.xconv_comp_ins[0, "Gd"].fix(0.023783784)
    m.fs.roaster.xconv_comp_ins[0, "Tb"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Dy"].fix(0.712053571)
    m.fs.roaster.xconv_comp_ins[0, "Ho"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Er"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Tm"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Yb"].fix(0.766633166)
    m.fs.roaster.xconv_comp_ins[0, "Lu"].fix(1)

    # recovery fraction of all REE
    m.fs.roaster.frac_comp_ree_recovery.fix(0.99)

    # recovery fraction of impurity minerals
    m.fs.roaster.frac_impurity_recovery.fix(0.99)

    return m


@pytest.mark.unit
def test_build(model):
    assert hasattr(model.fs, "roaster")
    assert isinstance(model.fs.roaster, REEFeedRoaster)
    assert len(model.fs.roaster.config) == 9
    assert not model.fs.roaster.config.dynamic
    assert not model.fs.roaster.config.has_holdup
    assert model.fs.roaster.config.has_heat_transfer
    assert model.fs.roaster.config.has_pressure_change
    assert model.fs.roaster.config.gas_property_package is model.fs.prop_gas
    assert model.fs.roaster.config.solid_property_package is model.fs.prop_solid
    assert len(model.fs.prop_gas.component_list) == 4
    assert len(model.fs.roaster.ree_list) == 16
    assert isinstance(model.fs.roaster.heat_duty, Var)
    assert isinstance(model.fs.roaster.deltaP, Var)
    assert isinstance(model.fs.roaster.flow_mol_outlet_eqn, Constraint)
    assert len(model.fs.roaster.flow_mol_outlet_eqn) == 4
    assert number_variables(model.fs.roaster) == 193
    assert number_total_constraints(model.fs.roaster) == 103
    assert number_unused_variables(model.fs.roaster) == 1
    assert_units_consistent(model.fs.roaster)


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve(model):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(model.fs.roaster)
    assert initializer.summary[model.fs.roaster]["status"] == InitializationStatus.Ok
    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=False)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    flow_mol_out_gas = value(model.fs.roaster.gas_out[0].flow_mol)
    assert flow_mol_out_gas == pytest.approx(82.01903587318019, rel=1e-5, abs=1e-6)
    mole_frac_h2o = value(model.fs.roaster.gas_out[0].mole_frac_comp["H2O"])
    assert mole_frac_h2o == pytest.approx(0.105061218702762, rel=1e-5, abs=1e-6)
    mole_frac_o2 = value(model.fs.roaster.gas_out[0].mole_frac_comp["O2"])
    assert mole_frac_o2 == pytest.approx(0.060924564, rel=1e-5, abs=1e-6)
    mole_frac_co2 = value(model.fs.roaster.gas_out[0].mole_frac_comp["CO2"])
    assert mole_frac_co2 == pytest.approx(0.0922570830375, rel=1e-5, abs=1e-6)
    heat_duty = value(model.fs.roaster.heat_duty[0])
    assert heat_duty == pytest.approx(-2541458.25, rel=1e-5, abs=1e-3)
    flow_mass_solid_out = value(model.fs.roaster.solid_out[0].flow_mass)
    assert flow_mass_solid_out == pytest.approx(579.369, rel=1e-5, abs=1e-6)
    mass_frac_comp_solid_out = {
        "inerts": 0.5438792755532115,
        "Sc2O3": 1.6043595369991272e-05,
        "Y2O3": 2.767022982988265e-05,
        "La2O3": 4.833796064532411e-05,
        "Ce2O3": 0.00010273370445831112,
        "Pr2O3": 1.931695289570546e-05,
        "Nd2O3": 4.557010896553694e-05,
        "Sm2O3": 1.050789201041692e-05,
        "Gd2O3": 8.427859759143537e-06,
        "Dy2O3": 5.146215251158445e-06,
        "Al2O3": 0.2874097993439199,
        "CaO": 0.013738348321299278,
        "Fe2O3": 0.15468882226238415,
    }
    for i in model.fs.prop_solid.component_list:
        assert value(model.fs.roaster.solid_out[0].mass_frac_comp[i]) == pytest.approx(
            mass_frac_comp_solid_out[i], rel=1e-5, abs=1e-6
        )
    ppm_insoluable_in_product = {
        "Sc": 1.8851305084719268,
        "Y": 5.703616099755223,
        "La": 8.191165847386474,
        "Ce": 23.710635917672906,
        "Pr": 7.907868653255506,
        "Nd": 11.64024057196155,
        "Sm": 2.648781847223353,
        "Eu": 0.9661802856697338,
        "Gd": 6.459496817385459,
        "Tb": 0.0,
        "Dy": 1.163415954814308,
        "Ho": 0.0,
        "Er": 0.0,
        "Tm": 0.0,
        "Yb": 0.6134119292750246,
        "Lu": 0.0,
    }
    for i in model.fs.roaster.ree_list:
        assert value(model.fs.roaster.ppm_comp_ree_ins_product[0, i]) == pytest.approx(
            ppm_insoluable_in_product[i], rel=1e-5, abs=1e-6
        )
    ppm_dissovable_in_product = {
        "Sc": 14.158464861519352,
        "Y": 21.966613730127424,
        "La": 40.146794797937645,
        "Ce": 79.02306854063823,
        "Pr": 11.409084242449953,
        "Nd": 33.92986839357539,
        "Sm": 7.859110163193565,
        "Eu": 0.4011715926573889,
        "Gd": 1.9683629417580777,
        "Tb": 0.7458282972693396,
        "Dy": 3.982799296344135,
        "Ho": 1.491656594538679,
        "Er": 4.615848461989136,
        "Tm": 1.0524465972578458,
        "Yb": 2.7345284273562336,
        "Lu": 0.8784199945616665,
    }
    for i in model.fs.roaster.ree_list:
        assert value(model.fs.roaster.ppm_comp_ree_dis_product[0, i]) == pytest.approx(
            ppm_dissovable_in_product[i], rel=1e-5, abs=1e-6
        )
