#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    ConcreteModel,
    Param,
    SolverFactory,
    TransformationFactory,
    Var,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.math import smooth_max
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models import Feed
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)
from prommis.hydrogen_decrepitation.repm_solids_properties import REPMParameters


def main():

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

    # shredder and HDD feed to define REPM flow into hydrogen decrepitation furnace
    m.fs.shredder = Feed(property_package=m.fs.prop_solid)

    m.fs.plant_basis_year = Param(
        initialize=2023, units=pyunits.dimensionless
    )  # year of plant operation

    m.fs.flow_2_5inch_HDDs = Var(m.fs.time, initialize=600, units=pyunits.h**-1)
    m.fs.flow_2_5inch_HDDs.fix()
    m.fs.flow_3_5inch_HDDs = Var(m.fs.time, initialize=2100, units=pyunits.h**-1)
    m.fs.flow_3_5inch_HDDs.fix()

    @m.fs.shredder.Constraint(m.fs.time)
    def HDD_to_REPM_conversion_constraint(b, t):
        return b.flow_mass[t] == (
            # 2.5 g per 2.5 inch HDD
            pyunits.convert(
                2.5 * pyunits.g * m.fs.flow_2_5inch_HDDs[t],
                to_units=pyunits.kg / pyunits.s,
            )
            +
            # (17.87 - 0.35 * t) g per 3.5 inch HDD
            # t is the manufacture year of the HDDs in years since 1990, disks require less material every year
            # assume 8 year lifetime of HDDs, so manufacture year = plant basis year - 8
            # assume 3.5 inch HDDs will not get smaller than 2.5 inch HDDs, stop at 2.5 g REPM per HDD
            pyunits.convert(
                smooth_max(2.5, 17.87 - 0.35 * ((m.fs.plant_basis_year - 8) - 1990))
                * pyunits.g
                * m.fs.flow_3_5inch_HDDs[t],
                to_units=pyunits.kg / pyunits.s,
            )
        )

    m.fs.shredder.mass_frac_comp[0, "Nd2Fe14B"].fix(0.99)

    # don't fix, already have mole frac balance so just need initial value
    m.fs.shredder.mass_frac_comp[0, "Nd"] = 0.01

    m.fs.hydrogen_decrepitation_furnace = REPMHydrogenDecrepitationFurnace(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=True,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Nd",
        ],
        number_of_units=1,
    )

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
    m.fs.hydrogen_decrepitation_furnace.deltaP.fix(0)
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

    # no additional heat is supplied other than what's required for decrepitation
    m.fs.hydrogen_decrepitation_furnace.supplied_heat_duty.fix(0)

    # connect shredder and furnace
    m.fs.shredded_REPM = Arc(
        source=m.fs.shredder.outlet,
        destination=m.fs.hydrogen_decrepitation_furnace.solid_inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    return m


def initialize_and_solve(m):

    initializer = BlockTriangularizationInitializer()
    initializer.initialize(m.fs.shredder)
    propagate_state(m.fs.shredded_REPM)

    # m.fs.hydrogen_decrepitation_furnace.gas_outlet.temperature.unfix()
    m.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.deactivate()  # flow mol will be fixed by initializer
    m.fs.hydrogen_decrepitation_furnace.solid_in[
        0
    ].sum_mass_frac.deactivate()  # mass frac will be fixed by initializer

    initializer.initialize(m.fs.hydrogen_decrepitation_furnace)

    # m.fs.hydrogen_decrepitation_furnace.gas_outlet.temperature.fix()
    m.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.activate()
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].sum_mass_frac.activate()

    assert (
        initializer.summary[m.fs.hydrogen_decrepitation_furnace]["status"]
        == InitializationStatus.Ok
    )

    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    dt = DiagnosticsToolbox(m)
    dt.assert_no_numerical_warnings()

    m.fs.shredder.report()
    m.fs.hydrogen_decrepitation_furnace.report()


if __name__ == "__main__":
    m = main()

    initialize_and_solve(m)
