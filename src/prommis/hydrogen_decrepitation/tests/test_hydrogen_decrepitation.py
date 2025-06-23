#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    SolverFactory,
    TransformationFactory,
    Var,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.math import smooth_max
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models import Feed
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

import pytest

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)
from prommis.hydrogen_decrepitation.repm_solids_properties import REPMParameters


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {
        "H2",
    }
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_solid = REPMParameters(
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
        "H2": 1,
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )

    m.fs.hydrogen_decrepitation_furnace = REPMHydrogenDecrepitationFurnace(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Nd",
        ],
    )

    m.fs.hydrogen_decrepitation_furnace.solid_in[0].flow_mass.fix(0.0057367 * pyunits.kg/pyunits.s)
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd2Fe14B"].fix(0.99)
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd"] = 0.01

    # don't fix, already have mole frac balance so just need initial value
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.mole_frac_comp[0, "H2"].fix(1)


    # inlet flue gas mole flow rate, stoichiometric on molar basis with REPM
    @m.fs.hydrogen_decrepitation_furnace.Constraint(m.fs.time)
    def flow_mol_gas_constraint(b, t):
        return b.gas_inlet.flow_mol[t] == (
            sum(b.flow_mol_comp_impurity_feed[t, c] for c in m.fs.prop_solid.component_list)
        )

    # operating parameters
    m.fs.hydrogen_decrepitation_furnace.deltaP.fix(0)
    m.fs.hydrogen_decrepitation_furnace.operating_temperature.fix(443.15)
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.pressure.fix(101325)

    m.fs.hydrogen_decrepitation_furnace.decrepitation_duration.set_value(10800 * pyunits.s)
    m.fs.hydrogen_decrepitation_furnace.sample_density.set_value(
        m.fs.prop_solid.dens_mass
    )  # 7500 kg/m3
    m.fs.hydrogen_decrepitation_furnace.chamber_to_sample_ratio.set_value(2)

    # solid temperature, cools back to inlet temperature during shutdown
    m.fs.hydrogen_decrepitation_furnace.temp_feed.fix(298.15)
    m.fs.hydrogen_decrepitation_furnace.temp_prod.fix(298.15)

    # gas temperature, assume comes in at operating temperature
    m.fs.hydrogen_decrepitation_furnace.gas_in[0].temperature.fix(443.15)

    # no additional heat is supplied other than what's required for decrepitation
    m.fs.hydrogen_decrepitation_furnace.supplied_heat_duty.fix(0)

    return m


# TODO update these once all model equations are added
@pytest.mark.unit
def test_build(model):
     assert hasattr(model.fs, "hydrogen_decrepitation_furnace")
#    assert isinstance(
#        model.fs.hydrogen_decrepitation_furnace, REPMHydrogenDecrepitationFurnace
#    )
#    assert len(model.fs.hydrogen_decrepitation_furnace.config) == 9
#    assert not model.fs.hydrogen_decrepitation_furnace.config.dynamic
#    assert not model.fs.hydrogen_decrepitation_furnace.config.has_holdup
#    assert model.fs.hydrogen_decrepitation_furnace.config.has_heat_transfer
#    assert model.fs.hydrogen_decrepitation_furnace.config.has_pressure_change
#    assert (
#        model.fs.hydrogen_decrepitation_furnace.config.gas_property_package
#        is model.fs.prop_gas
#    )
#    assert (
#        model.fs.hydrogen_decrepitation_furnace.config.solid_property_package
#        is model.fs.prop_solid
#    )
#    assert len(model.fs.prop_gas.component_list) == 1
#    assert len(model.fs.hydrogen_decrepitation_furnace.ree_list) == 1
#    assert isinstance(model.fs.hydrogen_decrepitation_furnace.heat_duty, Var)
#    assert isinstance(model.fs.hydrogen_decrepitation_furnace.deltaP, Var)
#    assert isinstance(
#        model.fs.hydrogen_decrepitation_furnace.flow_mol_outlet_eqn, Constraint
#    )
#    assert len(model.fs.hydrogen_decrepitation_furnace.flow_mol_outlet_eqn) == 1
#    assert number_variables(model.fs.hydrogen_decrepitation_furnace) == 31
#    assert number_total_constraints(model.fs.hydrogen_decrepitation_furnace) == 23
#    assert number_unused_variables(model.fs.hydrogen_decrepitation_furnace) == 0
#    assert_units_consistent(model.fs.hydrogen_decrepitation_furnace)


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve(model):
    initializer = BlockTriangularizationInitializer()

    model.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.deactivate()  # flow mol will be fixed by initializer
    model.fs.hydrogen_decrepitation_furnace.solid_in[
        0
    ].sum_mass_frac.deactivate()  # mass frac will be fixed by initializer

    initializer.initialize(model.fs.hydrogen_decrepitation_furnace)

    # m.fs.hydrogen_decrepitation_furnace.gas_outlet.temperature.fix()
    model.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.activate()
    model.fs.hydrogen_decrepitation_furnace.solid_in[0].sum_mass_frac.activate()

    assert (
        initializer.summary[model.fs.hydrogen_decrepitation_furnace]["status"]
        == InitializationStatus.Ok
    )

    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=True)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


# TODO update these once all model equations are added, these are old results from the feed roaster test file
# @pytest.mark.component
# @pytest.mark.solver
# def test_solution(model):
#     flow_mol_out_gas = value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].flow_mol)
#     assert flow_mol_out_gas == pytest.approx(82.01903587318019, rel=1e-5, abs=1e-6)
#     mole_frac_h2o = value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].mole_frac_comp["H2O"])
#     assert mole_frac_h2o == pytest.approx(0.105061218702762, rel=1e-5, abs=1e-6)
#     mole_frac_o2 = value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].mole_frac_comp["O2"])
#     assert mole_frac_o2 == pytest.approx(0.060924564, rel=1e-5, abs=1e-6)
#     mole_frac_co2 = value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].mole_frac_comp["CO2"])
#     assert mole_frac_co2 == pytest.approx(0.0922570830375, rel=1e-5, abs=1e-6)
#     heat_duty = value(model.fs.hydrogen_decrepitation_furnace.heat_duty[0])
#     assert heat_duty == pytest.approx(-2541458.25, rel=1e-5, abs=1e-3)
#     flow_mass_solid_out = value(model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass)
#     assert flow_mass_solid_out == pytest.approx(579.369, rel=1e-5, abs=1e-6)
#     mass_frac_comp_solid_out = {
#         "inerts": 0.5438792755532115,
#         "Sc2O3": 1.6043595369991272e-05,
#         "Y2O3": 2.767022982988265e-05,
#         "La2O3": 4.833796064532411e-05,
#         "Ce2O3": 0.00010273370445831112,
#         "Pr2O3": 1.931695289570546e-05,
#         "Nd2O3": 4.557010896553694e-05,
#         "Sm2O3": 1.050789201041692e-05,
#         "Gd2O3": 8.427859759143537e-06,
#         "Dy2O3": 5.146215251158445e-06,
#         "Al2O3": 0.2874097993439199,
#         "CaO": 0.013738348321299278,
#         "Fe2O3": 0.15468882226238415,
#     }
#     for i in model.fs.prop_solid.component_list:
#         assert value(model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp[i]) == pytest.approx(
#             mass_frac_comp_solid_out[i], rel=1e-5, abs=1e-6
#         )
#     ppm_insoluable_in_product = {
#         "Sc": 1.8851305084719268,
#         "Y": 5.703616099755223,
#         "La": 8.191165847386474,
#         "Ce": 23.710635917672906,
#         "Pr": 7.907868653255506,
#         "Nd": 11.64024057196155,
#         "Sm": 2.648781847223353,
#         "Eu": 0.9661802856697338,
#         "Gd": 6.459496817385459,
#         "Tb": 0.0,
#         "Dy": 1.163415954814308,
#         "Ho": 0.0,
#         "Er": 0.0,
#         "Tm": 0.0,
#         "Yb": 0.6134119292750246,
#         "Lu": 0.0,
#     }
#     for i in model.fs.hydrogen_decrepitation_furnace.ree_list:
#         assert value(model.fs.hydrogen_decrepitation_furnace.ppm_comp_ree_ins_product[0, i]) == pytest.approx(
#             ppm_insoluable_in_product[i], rel=1e-5, abs=1e-6
#         )
#     ppm_dissovable_in_product = {
#         "Sc": 14.158464861519352,
#         "Y": 21.966613730127424,
#         "La": 40.146794797937645,
#         "Ce": 79.02306854063823,
#         "Pr": 11.409084242449953,
#         "Nd": 33.92986839357539,
#         "Sm": 7.859110163193565,
#         "Eu": 0.4011715926573889,
#         "Gd": 1.9683629417580777,
#         "Tb": 0.7458282972693396,
#         "Dy": 3.982799296344135,
#         "Ho": 1.491656594538679,
#         "Er": 4.615848461989136,
#         "Tm": 1.0524465972578458,
#         "Yb": 2.7345284273562336,
#         "Lu": 0.8784199945616665,
#     }
#     for i in model.fs.hydrogen_decrepitation_furnace.ree_list:
#         assert value(model.fs.hydrogen_decrepitation_furnace.ppm_comp_ree_dis_product[0, i]) == pytest.approx(
#             ppm_dissovable_in_product[i], rel=1e-5, abs=1e-6
#         )
