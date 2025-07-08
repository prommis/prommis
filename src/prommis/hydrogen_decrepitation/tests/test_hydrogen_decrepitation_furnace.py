#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, SolverFactory, assert_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.environ import value
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
        has_heat_transfer=False,
        has_pressure_change=False,
        ree_list=[
            "Nd",
        ],
    )

    m.fs.hydrogen_decrepitation_furnace.solid_in[0].flow_mass.fix(
        0.0057367 * pyunits.kg / pyunits.s
    )
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd2Fe14B"].fix(0.99)
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd"] = 0.01

    # don't fix, already have mole frac balance so just need initial value
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.mole_frac_comp[0, "H2"].fix(1)

    # inlet flue gas mole flow rate, stoichiometric on molar basis with REPM
    @m.fs.hydrogen_decrepitation_furnace.Constraint(m.fs.time)
    def flow_mol_gas_constraint(b, t):
        return b.gas_inlet.flow_mol[t] == (
            sum(
                b.flow_mol_comp_impurity_feed[t, c]
                for c in m.fs.prop_solid.component_list
            )
        )

    # operating parameters
    m.fs.hydrogen_decrepitation_furnace.operating_temperature.fix(443.15)
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.pressure.fix(101325)

    m.fs.hydrogen_decrepitation_furnace.decrepitation_duration.set_value(
        10800 * pyunits.s
    )
    m.fs.hydrogen_decrepitation_furnace.sample_density.set_value(
        m.fs.prop_solid.dens_mass
    )  # 7500 kg/m3
    m.fs.hydrogen_decrepitation_furnace.chamber_to_sample_ratio.set_value(2)

    # solid temperature, cools back to inlet temperature during shutdown
    m.fs.hydrogen_decrepitation_furnace.temp_feed.fix(298.15)
    m.fs.hydrogen_decrepitation_furnace.temp_prod.fix(298.15)

    # gas temperature, assume comes in at operating temperature
    m.fs.hydrogen_decrepitation_furnace.gas_in[0].temperature.fix(443.15)

    return m


# TODO update these once all model equations are added
@pytest.mark.unit
def test_build(model):
    assert hasattr(model.fs, "hydrogen_decrepitation_furnace")
    assert isinstance(
        model.fs.hydrogen_decrepitation_furnace, REPMHydrogenDecrepitationFurnace
    )
    assert len(model.fs.hydrogen_decrepitation_furnace.config) == 10
    assert not model.fs.hydrogen_decrepitation_furnace.config.dynamic
    assert not model.fs.hydrogen_decrepitation_furnace.config.has_holdup
    assert not model.fs.hydrogen_decrepitation_furnace.config.has_heat_transfer
    assert not model.fs.hydrogen_decrepitation_furnace.config.has_pressure_change
    assert (
        model.fs.hydrogen_decrepitation_furnace.config.gas_property_package
        is model.fs.prop_gas
    )
    assert (
        model.fs.hydrogen_decrepitation_furnace.config.solid_property_package
        is model.fs.prop_solid
    )
    assert len(model.fs.prop_gas.component_list) == 1
    assert len(model.fs.hydrogen_decrepitation_furnace.ree_list) == 1

    assert number_variables(model.fs.hydrogen_decrepitation_furnace) == 47
    assert number_total_constraints(model.fs.hydrogen_decrepitation_furnace) == 39
    assert number_unused_variables(model.fs.hydrogen_decrepitation_furnace) == 0
    assert_units_consistent(model.fs.hydrogen_decrepitation_furnace)


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


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    model.fs.hydrogen_decrepitation_furnace.report()
    assert value(
        model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass
    ) == pytest.approx(0.0057367, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd"]
    ) == pytest.approx(0.0100001, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd2Fe14B"]
    ) == pytest.approx(0.99, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.gas_out[0].flow_mol
    ) == pytest.approx(0.0056509, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.gas_out[0].mole_frac_comp["H2"]
    ) == pytest.approx(1.0000, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.gas_out[0].temperature
    ) == pytest.approx(443.15, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.gas_out[0].pressure
    ) == pytest.approx(1.01325e5, rel=1e-5)

    assert value(
        model.fs.hydrogen_decrepitation_furnace.flow_vol_feed[0]
    ) == pytest.approx(7.6489e-7, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.total_heat_duty[0]
    ) == pytest.approx(3771.14, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.sample_mass[0]
    ) == pytest.approx(61.956, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.sample_volume[0]
    ) == pytest.approx(0.0082608, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.furnace_chamber_volume[0]
    ) == pytest.approx(0.0165217, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.radius_chamber[0]
    ) == pytest.approx(0.095701, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.length_chamber[0]
    ) == pytest.approx(0.57421, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.volume_insulation1[0]
    ) == pytest.approx(0.091450, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.total_weight_insulation1[0]
    ) == pytest.approx(15.420, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.internal_diameter_metal1[0]
    ) == pytest.approx(0.48930, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.external_diameter_metal1[0]
    ) == pytest.approx(0.50056, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.volume_metal1[0]
    ) == pytest.approx(
        306.781, rel=1e-5
    )  # in3, report has m3
    assert value(
        model.fs.hydrogen_decrepitation_furnace.weight_metal1[0]
    ) == pytest.approx(
        82.8308, rel=1e-5
    )  # lb, report has kg
    assert value(
        model.fs.hydrogen_decrepitation_furnace.volume_insulation2[0]
    ) == pytest.approx(0.226876, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.total_weight_insulation2[0]
    ) == pytest.approx(58.708, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.internal_diameter_metal2[0]
    ) == pytest.approx(0.86812, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.external_diameter_metal2[0]
    ) == pytest.approx(0.88342, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.volume_metal2[0]
    ) == pytest.approx(
        737.2515, rel=1e-5
    )  # in3, report has m3
    assert value(
        model.fs.hydrogen_decrepitation_furnace.weight_metal2[0]
    ) == pytest.approx(
        209.379, rel=1e-5
    )  # lb, report has kg
    assert value(
        model.fs.hydrogen_decrepitation_furnace.furnace_external_surface_area[0]
    ) == pytest.approx(2.8195, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.decrepitation_duration
    ) == pytest.approx(10800, rel=1e-5)
    assert value(
        model.fs.hydrogen_decrepitation_furnace.processing_time[0]
    ) == pytest.approx(
        5.0833, rel=1e-5
    )  # h, report has s
    assert value(
        model.fs.hydrogen_decrepitation_furnace.furnace_volume[0]
    ) == pytest.approx(12.4292, rel=1e-5)  # in3, report has m3
    assert value(
        model.fs.hydrogen_decrepitation_furnace.furnace_weight[0]
    ) == pytest.approx(455.635, rel=1e-5)  # lb, report has kg
