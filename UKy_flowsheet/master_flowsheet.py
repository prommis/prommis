#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Initial flowsheet for UKy plant

Authors: Marcus Holly
"""


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    TransformationFactory,
    units,
    Var,
    Set,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.mscontactor import (
    MSContactor,
    MSContactorInitializer,
)

from idaes.models.unit_models.separator import (
    Separator,
    SplittingType,
    SeparatorInitializer,
    MaterialBalanceType,
)

from idaes.models.unit_models.mixer import (
    Mixer,
    MaterialBalanceType,
    MixingType,
    MomentumMixingType,
    MixerInitializer,
)

from idaes.models.unit_models import Feed, Translator
from idaes.models.unit_models.product import (
    Product,
    ProductInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from workspace.UKy_flowsheet.old_leaching.leach_solution_properties import LeachSolutionParameters
from workspace.UKy_flowsheet.old_leaching.leach_solids_properties import CoalRefuseParameters
from workspace.UKy_flowsheet.old_leaching.leach_reactions import CoalRefuseLeachingReactions

from workspace.Solvent_Extraction.REESXmodel import REESX
from workspace.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from workspace.Solvent_Extraction.REEOgdistribution import REESolExOgParameters

from workspace.precipitate.precipitator import Precipitator
from workspace.precipitate.precip_prop import AqueousStateParameterBlock, PrecipitateStateParameterBlock

from workspace.roasting.ree_oxalate_roster import REEOxalateRoaster

from workspace.UKy_flowsheet.Translators.translator_leaching_SX import Translator_leaching_SX
from workspace.UKy_flowsheet.Translators.translator_SX_precipitator import Translator_SX_precipitator

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)

from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)

from idaes.core.util.initialization import propagate_state

from idaes.core.util.model_diagnostics import DiagnosticsToolbox

def main():
    m = build()

    set_operating_conditions(m)
    assert_units_consistent(m)
    # assert degrees_of_freedom(m) == 0

    print("Structural issues after setting operating conditions")
    dt = DiagnosticsToolbox(model=m)
    dt.report_structural_issues()
    dt.display_underconstrained_set()

    initialize_system(m)
    # print("Numerical issues after initialization")
    # dt.report_numerical_issues()
    # dt.display_constraints_with_large_residuals()
    # dt.display_variables_at_or_outside_bounds()

    results = solve(m)
    # print("Numerical issues after solving")
    # dt.report_numerical_issues()
    # dt.display_constraints_with_large_residuals()
    # dt.display_variables_at_or_outside_bounds()

    display_results(m)

    return m, results

def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Leaching property and unit models
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    m.fs.leach = MSContactor(
        number_of_finite_elements=1,
        streams={
            "liquid": {
                "property_package": m.fs.leach_soln,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            "solid": {
                "property_package": m.fs.coal,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        },
        heterogeneous_reactions=m.fs.leach_rxns,
    )

    m.fs.leach_filter_cake = Product(property_package=m.fs.coal)

    # ----------------------------------------------------------------------------------------------------------------
    # Solvent extraction property and unit models
    m.fs.prop_a = REESolExAqParameters()
    m.fs.prop_o = REESolExOgParameters()

    m.fs.leach_to_SX = Translator_leaching_SX(
        inlet_property_package=m.fs.leach_soln,
        outlet_property_package=m.fs.prop_a,
    )

    m.fs.solex = REESX(number_of_finite_elements=3,
                           aqueous_streams = {"Acidsoln":{"property_package":m.fs.prop_a, "flow_direction":1}},
                           organic_streams = {"Orgacid":{"property_package":m.fs.prop_o, "flow_direction":2}})

    m.fs.sx_leach_acid = Product(property_package=m.fs.prop_o)

    m.fs.sx_acid_soln = Product(property_package=m.fs.prop_a)

    # --------------------------------------------------------------------------------------------------------------
    # Precipitation property and unit models

    key_components = {
        "H^+",
        "Ce^3+",
        "Al^3+",
        "Fe^3+",
        # "Fe^2+",
        "Ca^2+",
        # "Mg^2+",
        "C2O4^2-",
        # "NO3^-",
        # "SO4^2-",
        # "Cl^-",
    }

    m.fs.properties_aq = AqueousStateParameterBlock(
        key_components=key_components,
    )
    m.fs.properties_solid = PrecipitateStateParameterBlock(
        key_components=key_components,
    )

    m.fs.SX_to_precipitator = Translator_SX_precipitator(
        inlet_property_package=m.fs.prop_a,
        outlet_property_package=m.fs.properties_aq,
    )

    m.fs.oxalate_feed = Feed(property_package=m.fs.properties_aq)

    # Mixer by itself has 3 variables but only one constraint
    m.fs.mixer = Mixer(
        property_package=m.fs.properties_aq,
        num_inlets=2,
        inlet_list=["SX_inlet", "oxalate_inlet"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        # momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.mixed_product = Product(property_package=m.fs.properties_aq)

    m.fs.precipitator = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    m.fs.precipitate_feed = Feed(property_package=m.fs.properties_solid)

    m.fs.solid_product = Product(property_package=m.fs.properties_solid)
    m.fs.liquid_product = Product(property_package=m.fs.properties_aq)

    # TODO: Replace with a filter press
    # m.fs.S01 = Separator(
    #     property_package=m.fs.leach_soln,
    #     outlet_list=["recycle", "product"],
    #     split_basis=SplittingType.phaseFlow,
    # )

    # -----------------------------------------------------------------------------------------------------------------
    # Roasting property and unit models

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )

    m.fs.prop_solid = PrecipitateStateParameterBlock(
        key_components=key_components,
    )

    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # -----------------------------------------------------------------------------------------------------------------
    # Flowsheet connections

    # Connections without recycle loops
    m.fs.s01 = Arc(source=m.fs.leach.solid_outlet, destination=m.fs.leach_filter_cake.inlet)
    m.fs.s02 = Arc(source=m.fs.leach.liquid_outlet, destination=m.fs.leach_to_SX.inlet)
    m.fs.s03 = Arc(source=m.fs.leach_to_SX.outlet, destination=m.fs.solex.Acidsoln_inlet)
    # m.fs.s03 = Arc(source=m.fs.oxalic_acid_feed.outlet, destination=m.fs.solex.Orgacid_inlet)
    m.fs.s04 = Arc(source=m.fs.solex.Orgacid_outlet, destination=m.fs.sx_leach_acid.inlet) # Should eventually convert to a recycle
    # m.fs.s05 = Arc(source=m.fs.solex.Acidsoln_outlet, destination=m.fs.sx_acid_soln.inlet)
    m.fs.s05 = Arc(source=m.fs.solex.Acidsoln_outlet, destination=m.fs.SX_to_precipitator.inlet)
    m.fs.s06 = Arc(source=m.fs.SX_to_precipitator.outlet, destination=m.fs.mixer.SX_inlet)
    m.fs.s07 = Arc(source=m.fs.oxalate_feed.outlet, destination=m.fs.mixer.oxalate_inlet)
    # m.fs.s08 = Arc(source=m.fs.mixer.outlet, destination=m.fs.mixed_product.inlet)
    m.fs.s08 = Arc(source=m.fs.mixer.outlet, destination=m.fs.precipitator.aqueous_inlet)
    m.fs.s09 = Arc(source=m.fs.precipitate_feed.outlet, destination=m.fs.precipitator.precipitate_inlet)
    m.fs.s10 = Arc(source=m.fs.precipitator.aqueous_outlet, destination=m.fs.liquid_product.inlet)  # Should eventually convert to a recycle
    # m.fs.s11 = Arc(source=m.fs.precipitator.precipitate_outlet, destination=m.fs.solid_product.inlet)
    m.fs.s11 = Arc(source=m.fs.precipitator.precipitate_outlet, destination=m.fs.roaster.solid_inlet)

    # Connections with recycle loops
    # m.fs.s01 = Arc(source=m.fs.leach.solid_outlet, destination=m.fs.leach_filter_cake.inlet)
    # m.fs.s02 = Arc(source=m.fs.leach.liquid_outlet, destination=m.fs.solex.Acidsoln_inlet_state)
    # m.fs.s03 = Arc(source=m.fs.S01.recycle, destination=m.fs.solex.Orgacid_inlet_state)
    # m.fs.s04 = Arc(source=m.fs.solex.Orgacid_outlet_state, destination=m.fs.sx_leach_acid) # Should eventually convert to a recycle
    # m.fs.s05 = Arc(source=m.fs.solex.Acidsoln_outlet_state, destination=m.fs.precipitator.aqueous_inlet)
    # m.fs.s06 = Arc(source=m.fs.oxalic_acid_feed.outlet, destination=m.fs.precipitator.precipitate_inlet)
    # m.fs.s07 = Arc(source=m.fs.precipitator.aqueous_outlet, destination=m.fs.M01.SX_acid)
    # m.fs.s08 = Arc(source=m.fs.precipitator.precipitate_outlet, destination=m.fs.M01.oxalic_acid)
    # m.fs.s09 = Arc(source=m.fs.M01.outlet, destination=m.fs.S01.inlet)
    # m.fs.s10 = Arc(source=m.fs.S01.product, destination=m.fs.roaster.inlet)
    # m.fs.s11 = Arc(source=m.fs.roaster.product, destination=m.fs.S01.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m

def set_operating_conditions(m):
    # Liquid feed to old_leaching unit
    m.fs.leach.liquid_inlet.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach.liquid_inlet.conc_mass_metals.fix(1e-10 * units.mg / units.L)    # Why fixed at lower bound?
    m.fs.leach.liquid_inlet.conc_mole_acid[0, "H"].fix(2 * 0.05 * units.mol / units.L)
    m.fs.leach.liquid_inlet.conc_mole_acid[0, "HSO4"].fix(1e-8 * units.mol / units.L)   # Why fixed at lower bound?
    m.fs.leach.liquid_inlet.conc_mole_acid[0, "SO4"].fix(0.05 * units.mol / units.L)

    # Solid feed to old_leaching unit
    m.fs.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sc2O3"].fix(2.77966E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Y2O3"].fix(3.28653E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "La2O3"].fix(6.77769E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Ce2O3"].fix(0.000156161 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Pr2O3"].fix(1.71438E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Nd2O3"].fix(6.76618E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sm2O3"].fix(1.47926E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Gd2O3"].fix(1.0405E-05 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Dy2O3"].fix(7.54827E-06 * units.kg / units.kg)

    # TODO: Replace this with a recycle loop
    # SX rougher recycle stream (treated as a product for now)
    m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01 * units.L / units.hour)

    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Al"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ca"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Fe"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sc"].fix(19.93 * units.g / units.hour)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Y"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["La"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ce"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pr"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Nd"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sm"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Gd"].fix(0)
    m.fs.solex.Orgacid_inlet_state[0].flow_mass["Dy"].fix(0)

    # Oxalic acid feed to precipitator
    #TODO: Use the appropriate value for flow_mass
    m.fs.oxalate_feed.properties[0].flow_mass.fix(1)
    m.fs.oxalate_feed.properties[0].temperature.fix(300 * units.kelvin)
    m.fs.oxalate_feed.properties[0].pressure.fix(101325 * units.Pa)
    m.fs.oxalate_feed.properties[0].log10_molality_comp["C2O4^2-"].fix(-4)
    # for key in m.fs.properties_aq.aq_dict.keys():
    #     m.fs.oxalate_feed.properties[0].log10_molality_comp[key].fix(1)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["H2C2O4"].fix(1)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["Ca^2+"].fix(0)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["Fe^3+"].fix(0)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["Al^3+"].fix(0)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["Ce^3+"].fix(0)
    # m.fs.oxalate_feed.properties[0].log10_molality_comp["H^+"].fix(0)



    # Precipitate feed to precipitator
    #TODO: What should these feed conditions be? They are assumed to be zero in the example
    m.fs.precipitate_feed.properties[0].flow_mol_comp.fix(0)
    m.fs.precipitate_feed.properties[0].temperature.fix(300)

    # Reactor volume
    m.fs.leach.volume = Var(
        m.fs.time,
        m.fs.leach.elements,
        initialize=1,
        units=units.litre,
        doc="Volume of each finite element."
    )
    m.fs.leach.volume.fix(100 * units.gallon)

    def rule_heterogeneous_reaction_extent(b, t, s, r):
        return b.heterogeneous_reaction_extent[t, s, r] == b.heterogeneous_reactions[t, s].reaction_rate[r]*b.volume[t,s]

    m.fs.leach.heterogeneous_reaction_extent_constraint = Constraint(
        m.fs.time,
        m.fs.leach.elements,
        m.fs.leach_rxns.reaction_idx,
        rule=rule_heterogeneous_reaction_extent,
    )

    # Roaster gas feed
    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # Inlet flue gas mole flow rate
    fgas = 0.00781
    # Inlet flue gas composition, typical flue gas by buring CH4 with air with stoichiometric ratio 0f 2.3
    gas_comp = {
        "O2": 0.1118,
        "H2O": 0.1005,
        "CO2": 0.0431,
        "N2": 0.7446,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

    # Fix outlet product temperature
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # Fix operating conditions
    m.fs.roaster.flow_mol_moist_feed.fix(6.75e-4)
    m.fs.roaster.frac_comp_recovery.fix(0.95)


def initialize_system(m):
    # Initialize old_leaching section
    initializer1 = MSContactorInitializer()
    initializer1.initialize(m.fs.leach)

    propagate_state(m.fs.s01)
    initializer_product = ProductInitializer()
    initializer_product.initialize(m.fs.leach_filter_cake)

    # Initialize old_leaching -> SX translator
    propagate_state(m.fs.s02)

    initializer2 = BlockTriangularizationInitializer()
    initializer2.initialize(m.fs.leach_to_SX)

    # Initialize SX section
    propagate_state(m.fs.s03)
    initializer2.initialize(m.fs.solex)

    # Initialize SX_leach_acid, which will eventually be a recycle
    propagate_state(m.fs.s04)
    initializer_product.initialize(m.fs.sx_leach_acid)

    # Initialize SX -> precipitation translator
    propagate_state(m.fs.s05)

    initializer2.initialize(m.fs.SX_to_precipitator)

    # Initialize mixer
    propagate_state(m.fs.s06)
    propagate_state(m.fs.s07)

    initializer3 = MixerInitializer()
    initializer3.initialize(m.fs.mixer)

    # Initialize precipitator
    propagate_state(m.fs.s08)
    propagate_state(m.fs.s09)

    initializer2.initialize(m.fs.precipitator)

    # Initialize roaster
    propagate_state(m.fs.s11)
    initializer2.initialize(m.fs.roaster)

def solve(m):
    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)
    return results

def display_results(m):
    # m.fs.leach.liquid_outlet.display()
    # m.fs.leach.solid_outlet.display()
    # m.fs.sx_acid_soln.display()
    # m.fs.sx_acid_soln.report()
    # print("-------SX to Precipitator--------")
    # m.fs.SX_to_precipitator.display()
    # print("-------Mixer--------")
    # m.fs.mixer.display()
    # print("-------Mixer Product--------")
    # m.fs.mixed_product.display()
    # m.fs.liquid_product.display()
    # m.fs.solid_product.display()
    # m.fs.roaster.display()
    m.fs.leach_filter_cake.display()

if __name__ == "__main__":
    m, results = main()
