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
University of Kentucky Flowsheet

Author: Marcus Holly
"""


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    Suffix,
    TransformationFactory,
    Var,
    value,
    units,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models.feed import Feed, FeedInitializer
from idaes.models.unit_models.mixer import Mixer, MixingType, MomentumMixingType
from idaes.models.unit_models.mscontactor import MSContactor, MSContactorInitializer
from idaes.models.unit_models.product import Product, ProductInitializer
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
)
from idaes.models.unit_models.solid_liquid import SLSeparator
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_aq_distribution import REESolExAqParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction


def main():
    m = build()

    set_operating_conditions(m)

    scaled_model = set_scaling(m)
    assert_units_consistent(scaled_model)
    assert degrees_of_freedom(scaled_model) == 0

    print("Structural issues after setting operating conditions")
    dt = DiagnosticsToolbox(model=scaled_model)
    dt.report_structural_issues()

    initialize_system(scaled_model)
    print("Numerical issues after initialization")
    dt.report_numerical_issues()

    results = solve(scaled_model)
    print("Numerical issues after solving")
    dt.report_numerical_issues()

    display_results(scaled_model)

    return scaled_model, results


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Leaching property and unit models
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    m.fs.leach = MSContactor(
        number_of_finite_elements=2,
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

    # Reactor volume
    m.fs.leach.volume = Var(
        m.fs.time,
        m.fs.leach.elements,
        initialize=1,
        units=units.litre,
        doc="Volume of each finite element.",
    )

    def rule_heterogeneous_reaction_extent(b, t, s, r):
        return (
            b.heterogeneous_reaction_extent[t, s, r]
            == b.heterogeneous_reactions[t, s].reaction_rate[r] * b.volume[t, s]
        )

    m.fs.leach.heterogeneous_reaction_extent_constraint = Constraint(
        m.fs.time,
        m.fs.leach.elements,
        m.fs.leach_rxns.reaction_idx,
        rule=rule_heterogeneous_reaction_extent,
    )

    m.fs.sl_sep1 = SLSeparator(
        solid_property_package=m.fs.coal,
        liquid_property_package=m.fs.leach_soln,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.leach_mixer = Mixer(
        property_package=m.fs.leach_soln,
        num_inlets=2,
        inlet_list=["recycle", "feed"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.leach_liquid_feed = Feed(property_package=m.fs.leach_soln)
    m.fs.leach_solid_feed = Feed(property_package=m.fs.coal)

    m.fs.leach_filter_cake = Product(property_package=m.fs.coal)
    m.fs.leach_filter_cake_liquid = Product(property_package=m.fs.leach_soln)
    # ----------------------------------------------------------------------------------------------------------------
    # Solvent extraction property and unit models
    m.fs.prop_a = REESolExAqParameters()
    m.fs.prop_o = REESolExOgParameters()

    m.fs.solex_rougher = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
    )

    m.fs.sep1 = Separator(
        property_package=m.fs.leach_soln,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.sx_mixer = Mixer(
        property_package=m.fs.prop_o,
        num_inlets=2,
        inlet_list=["aqueous_inlet", "organic_inlet"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.recycle1_purge = Product(property_package=m.fs.leach_soln)
    # --------------------------------------------------------------------------------------------------------------
    # Precipitation property and unit models

    key_components = {
        "H^+",
        "Ce^3+",
        "Al^3+",
        "Fe^3+",
        "Ca^2+",
        "C2O4^2-",
    }

    m.fs.properties_aq = AqueousParameter()
    m.fs.properties_solid = PrecipitateParameters()

    m.fs.solex_cleaner = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.properties_aq,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
    )

    m.fs.precipitator = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    m.fs.sl_sep2 = SLSeparator(
        solid_property_package=m.fs.properties_solid,
        liquid_property_package=m.fs.properties_aq,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.sep2 = Separator(
        property_package=m.fs.properties_aq,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.recycle2_purge = Product(property_package=m.fs.properties_aq)
    # -----------------------------------------------------------------------------------------------------------------
    # Roasting property and unit models

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )

    m.fs.prop_solid = PrecipitateParameters(
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
    # UKy flowsheet with leach recycle loop
    m.fs.sol_feed = Arc(
        source=m.fs.leach_solid_feed.outlet, destination=m.fs.leach.solid_inlet
    )
    m.fs.liq_feed = Arc(
        source=m.fs.leach_liquid_feed.outlet, destination=m.fs.leach_mixer.feed
    )
    m.fs.feed_mixture = Arc(
        source=m.fs.leach_mixer.outlet, destination=m.fs.leach.liquid_inlet
    )
    m.fs.s01 = Arc(source=m.fs.leach.solid_outlet, destination=m.fs.sl_sep1.solid_inlet)
    m.fs.s02 = Arc(
        source=m.fs.leach.liquid_outlet, destination=m.fs.sl_sep1.liquid_inlet
    )
    m.fs.sep1_solid = Arc(
        source=m.fs.sl_sep1.solid_outlet, destination=m.fs.leach_filter_cake.inlet
    )
    m.fs.sep1_retained_liquid = Arc(
        source=m.fs.sl_sep1.retained_liquid_outlet,
        destination=m.fs.leach_filter_cake_liquid.inlet,
    )
    m.fs.sep1_liquid = Arc(
        source=m.fs.sl_sep1.recovered_liquid_outlet,
        destination=m.fs.solex_rougher.mscontactor.aqueous_inlet,
    )
    m.fs.recycle1 = Arc(
        source=m.fs.solex_rougher.mscontactor.aqueous_outlet,
        destination=m.fs.sep1.inlet,
    )
    m.fs.purge1 = Arc(source=m.fs.sep1.purge, destination=m.fs.recycle1_purge.inlet)
    m.fs.recycle_feed = Arc(
        source=m.fs.sep1.recycle, destination=m.fs.leach_mixer.recycle
    )
    m.fs.s03 = Arc(
        source=m.fs.solex_rougher.mscontactor.organic_outlet,
        destination=m.fs.solex_cleaner.mscontactor.organic_inlet,
    )
    m.fs.s04 = Arc(
        source=m.fs.solex_cleaner.mscontactor.aqueous_outlet,
        destination=m.fs.sx_mixer.aqueous_inlet,
    )
    m.fs.s05 = Arc(
        source=m.fs.solex_cleaner.mscontactor.organic_outlet,
        destination=m.fs.sx_mixer.organic_inlet,
    )
    m.fs.s06 = Arc(
        source=m.fs.sx_mixer.outlet,
        destination=m.fs.precipitator.aqueous_inlet,
    )
    m.fs.s07 = Arc(
        source=m.fs.precipitator.precipitate_outlet,
        destination=m.fs.sl_sep2.solid_inlet,
    )
    m.fs.s08 = Arc(
        source=m.fs.precipitator.aqueous_outlet, destination=m.fs.sl_sep2.liquid_inlet
    )
    m.fs.sep2_solid = Arc(
        source=m.fs.sl_sep2.solid_outlet, destination=m.fs.roaster.solid_inlet
    )
    # TODO: roaster model cannot currently handle liquid inlets
    # m.fs.sep2_retained_liquid = Arc(
    #     source=m.fs.sl_sep2.retained_liquid_outlet, destination=m.fs.roaster.liquid_inlet
    # )
    m.fs.sep2_recovered_liquid = Arc(
        source=m.fs.sl_sep2.recovered_liquid_outlet, destination=m.fs.sep2.inlet
    )
    m.fs.purge2 = Arc(source=m.fs.sep2.purge, destination=m.fs.recycle2_purge.inlet)
    m.fs.recycle2 = Arc(
        source=m.fs.sep2.recycle,
        destination=m.fs.solex_cleaner.mscontactor.aqueous_inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_scaling(m):
    # Scaling
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    component_set1 = [
        "H2O",
        "H",
        "HSO4",
        "SO4",
        "Sc",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Gd",
        "Dy",
        "Al",
        "Ca",
        "Fe",
    ]

    component_set2 = [
        "Sc",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Gd",
        "Dy",
        "Al",
        "Ca",
        "Fe",
    ]

    # Leaching
    for component in component_set1:
        m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_liquid_feed.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep1.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep1.split.recovered_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep1.split.retained_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.aqueous[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[m.fs.sep1.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep1.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep1.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.recycle1_purge.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_mixer.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.leach_mixer.feed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach_mixer.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5

    for component in component_set2:
        m.scaling_factor[
            m.fs.sl_sep2.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep2.split.recovered_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep2.split.retained_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.organic[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[m.fs.sep2.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep2.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep2.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.recycle2_purge.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.aqueous[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.organic[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.sx_mixer.aqueous_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sx_mixer.organic_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.sx_mixer.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.precipitator.cv_aqueous.properties_in[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.precipitator.cv_aqueous.properties_out[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5

    m.scaling_factor[m.fs.sl_sep2.solid_state[0].temperature] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.liquid_inlet_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.split.recovered_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.split.retained_state[0].flow_vol] = 1e-2

    m.scaling_factor[m.fs.sep2.mixed_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sep2.recycle_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sep2.purge_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.recycle2_purge.properties[0].flow_vol] = 1e-2

    m.scaling_factor[m.fs.solex_cleaner.mscontactor.aqueous[0, 1].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.aqueous[0, 2].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.aqueous[0, 3].flow_vol] = 1e-2
    m.scaling_factor[
        m.fs.solex_cleaner.mscontactor.aqueous_inlet_state[0].flow_vol
    ] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 1].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 2].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 3].flow_vol] = 1e-2

    m.scaling_factor[
        m.fs.precipitator.cv_precipitate.properties_in[0].temperature
    ] = 1e2
    m.scaling_factor[
        m.fs.precipitator.cv_precipitate.properties_out[0].temperature
    ] = 1e-4

    m.scaling_factor[m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.precipitator.cv_aqueous.properties_out[0].flow_vol] = 1e-2

    m.scaling_factor[m.fs.roaster.gas_in[0].flow_mol] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_in[0].flow_mol_phase["Vap"]] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_in[0].temperature] = 1e-2
    m.scaling_factor[m.fs.roaster.gas_in[0].pressure] = 1e-5
    m.scaling_factor[m.fs.roaster.gas_out[0].flow_mol_phase["Vap"]] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_out[0].flow_mol] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_out[0].temperature] = 1e-2
    m.scaling_factor[m.fs.roaster.gas_out[0].pressure] = 1e-5
    m.scaling_factor[m.fs.roaster.solid_in[0].temperature] = 1e-2

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    return scaled_model


def set_operating_conditions(m):
    eps = 1e-7 * units.mg / units.L

    m.fs.leach_liquid_feed.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach_liquid_feed.conc_mass_comp.fix(1e-10 * units.mg / units.L)
    m.fs.leach_liquid_feed.conc_mass_comp[0, "H"].fix(
        2 * 0.05 * 1e3 * units.mg / units.L
    )
    m.fs.leach_liquid_feed.conc_mass_comp[0, "HSO4"].fix(1e-8 * units.mg / units.L)
    m.fs.leach_liquid_feed.conc_mass_comp[0, "SO4"].fix(
        0.05 * 96e3 * units.mg / units.L
    )

    m.fs.leach_solid_feed.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach_solid_feed.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Sc2O3"].fix(
        2.77966e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Y2O3"].fix(
        3.28653e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "La2O3"].fix(
        6.77769e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Ce2O3"].fix(
        0.000156161 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Pr2O3"].fix(
        1.71438e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Nd2O3"].fix(
        6.76618e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Sm2O3"].fix(
        1.47926e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Gd2O3"].fix(
        1.0405e-05 * units.kg / units.kg
    )
    m.fs.leach_solid_feed.mass_frac_comp[0, "Dy2O3"].fix(
        7.54827e-06 * units.kg / units.kg
    )

    m.fs.leach.volume.fix(100 * units.gallon)

    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].flow_vol.fix(
        62.01 * units.L / units.hour
    )
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(
        321.34 * units.mg / units.L
    )
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(eps)
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(eps)

    m.fs.sl_sep1.liquid_recovery.fix(0.7)
    m.fs.sl_sep2.liquid_recovery.fix(0.7)

    m.fs.sep1.split_fraction[:, "recycle"].fix(0.9)
    m.fs.sep2.split_fraction[:, "recycle"].fix(0.9)

    m.fs.precipitator.cv_precipitate.properties_in[0].temperature.fix(348.15 * units.K)

    # Roaster gas feed
    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # Inlet flue gas mole flow rate
    fgas = 0.00781
    # Inlet flue gas composition, typical flue gas by burning CH4 with air with stoichiometric ratio of 2.3
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
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.feed_mixture, m.fs.recycle2]

    G = seq.create_graph(m)
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    tear_guesses1 = {
        "flow_mass": {0.007},
        "conc_mass_comp": {
            ("Al"): 1493939.39,
            ("Ca"): 501864.01,
            ("Ce"): 49698.79,
            ("Dy"): 466.86,
            ("Fe"): 1685228.86,
            ("Gd"): 1624.81,
            ("La"): 32143.22,
            ("Nd"): 12552.13,
            ("Pr"): 4084.40,
            ("Sc"): 17310.17,
            ("Sm"): 931.35,
            ("Y"): 2666.95,
        },
        "flow_mol_comp": {
            ("Al"): 577.62,
            ("Ca"): 64.65,
            ("Ce"): 1.17,
            ("Dy"): 0.0086,
            ("Fe"): 231.93,
            ("Gd"): 0.054,
            ("La"): 0.63,
            ("Nd"): 0.34,
            ("Pr"): 0.099,
            ("Sc"): 303.33,
            ("Sm"): 0.021,
            ("Y"): 0.090,
        },
    }
    tear_guesses2 = {
        "flow_mol_comp": {
            ("Al2(C2O4)3(s)"): 1.76,
            ("Ce2(C2O4)3(s)"): 2.65,
            ("Dy2(C2O4)3(s)"): 0.068,
            ("Fe2(C2O4)3(s)"): 2.64,
            ("Gd2(C2O4)3(s)"): 0.27,
            ("La2(C2O4)3(s)"): 0.86,
            ("Nd2(C2O4)3(s)"): 1.35,
            ("Pr2(C2O4)3(s)"): 0.36,
            ("Sc2(C2O4)3(s)"): 0.62,
            ("Sm2(C2O4)3(s)"): 0.15,
            ("Y2(C2O4)3(s)"): 0.31,
        },
    }
    tear_guesses3 = {
        "flow_vol": {0: 747.99},
        "conc_mass_comp": {
            (0, "Al"): 180.84,
            (0, "Ca"): 28.93,
            (0, "Ce"): 5.48,
            (0, "Dy"): 4.46e-11,
            (0, "Fe"): 269.98,
            (0, "Gd"): 2.60e-7,
            (0, "H"): 20.06,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 963.06,
            (0, "La"): 0.0037,
            (0, "Nd"): 1.81e-7,
            (0, "Pr"): 3.65e-6,
            (0, "SO4"): 486.24,
            (0, "Sc"): 4.17e-11,
            (0, "Sm"): 6.30e-10,
            (0, "Y"): 7.18e-11,
        },
    }
    tear_guesses4 = {
        "flow_vol": {0: 157.14},
        "conc_mass_comp": {
            (0, "Al"): 380.67,
            (0, "Ca"): 122.97,
            (0, "Ce"): 5.67,
            (0, "Dy"): 0.15,
            (0, "Fe"): 741.24,
            (0, "Gd"): 0.56,
            (0, "H"): 2.274,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 881.19,
            (0, "La"): 2.24,
            (0, "Nd"): 2.68,
            (0, "Pr"): 0.72,
            (0, "SO4"): 3924.07,
            (0, "Sc"): 0.10,
            (0, "Sm"): 0.30,
            (0, "Y"): 0.39,
        },
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.precipitator.cv_aqueous.properties_out[0], tear_guesses1)
    seq.set_guesses_for(
        m.fs.precipitator.cv_precipitate.properties_out[0], tear_guesses2
    )
    seq.set_guesses_for(m.fs.leach.liquid_inlet, tear_guesses3)
    seq.set_guesses_for(m.fs.solex_rougher.mscontactor.aqueous_inlet, tear_guesses4)

    def function(stream):
        initializer_feed = FeedInitializer()
        initializer_product = ProductInitializer()
        initializer1 = MSContactorInitializer()
        initializer2 = BlockTriangularizationInitializer()

        propagate_state(m.fs.liq_feed)
        propagate_state(m.fs.sol_feed)

        if stream == m.fs.leach_liquid_feed:
            initializer_feed.initialize(m.fs.leach_liquid_feed)
        elif stream == m.fs.leach_solid_feed:
            initializer_feed.initialize(m.fs.leach_solid_feed)
        elif stream == m.fs.leach_filter_cake:
            print(f"Initializing {stream}")
            initializer_product.initialize(m.fs.leach_filter_cake)
        elif stream == m.fs.leach:
            print(f"Initializing {stream}")
            try:
                initializer1.initialize(m.fs.leach)
            except:
                # Fix feed states
                m.fs.leach.liquid_inlet.flow_vol.fix()
                m.fs.leach.liquid_inlet.conc_mass_comp.fix()
                m.fs.leach.solid_inlet.flow_mass.fix()
                m.fs.leach.solid_inlet.mass_frac_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.leach, tee=True)
                # Unfix feed states
                m.fs.leach_liquid_feed.flow_vol.unfix()
                m.fs.leach.liquid_inlet.conc_mass_comp.unfix()
                m.fs.leach.solid_inlet.flow_mass.unfix()
                m.fs.leach.solid_inlet.mass_frac_comp.unfix()
        elif stream == m.fs.leach_mixer:
            initializer2.initialize(m.fs.leach_mixer)
        elif stream == m.fs.solex_rougher.mscontactor:
            print(f"Initializing {stream}")
            initializer2.initialize(m.fs.solex_rougher)
        elif stream == m.fs.solex_cleaner.mscontactor:
            print(f"Initializing {stream}")
            initializer2.initialize(m.fs.solex_cleaner)
        elif stream == m.fs.precipitator:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.precipitator)
            except:
                # Fix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.fix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.fix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.precipitator, tee=True)
                # Unfix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.unfix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.unfix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.unfix()
        else:
            print(f"Initializing {stream}")
            initializer2.initialize(stream)

    seq.run(m, function)


def solve(m):
    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)

    return results


def display_results(m):
    m.fs.sl_sep1.display()
    m.fs.leach_mixer.display()

    # Total mass basis yield calculation

    # product = value(units.convert(m.fs.roaster.flow_mas_product[0], to_units=units.kg / units.hr))
    # print(f"REE product mass flow is {product} kg/hr")
    #
    # REE_metal_mass_frac = {
    #     "Al2O3": 26.98 * 2 / (26.98 * 2 + 16 * 3),
    #     "Fe2O3": 55.845 * 2 / (55.845 * 2 + 16 * 3),
    #     "CaO": 40.078 / (40.078 + 16),
    #     "Sc2O3": 44.956 * 2 / (44.956 * 2 + 16 * 3),
    #     "Y2O3": 88.906 * 2 / (88.906 * 2 + 16 * 3),
    #     "La2O3": 138.91 * 2 / (138.91 * 2 + 16 * 3),
    #     "Ce2O3": 140.12 * 2 / (140.12 * 2 + 16 * 3),
    #     "Pr2O3": 140.91 * 2 / (140.91 * 2 + 16 * 3),
    #     "Nd2O3": 144.24 * 2 / (144.24 * 2 + 16 * 3),
    #     "Sm2O3": 150.36 * 2 / (150.36 * 2 + 16 * 3),
    #     "Gd2O3": 157.25 * 2 / (157.25 * 2 + 16 * 3),
    #     "Dy2O3": 162.5 * 2 / (162.5 * 2 + 16 * 3),
    # }
    #
    # feed_REE = sum(value(m.fs.leach_solid_feed.flow_mass[0] * m.fs.leach_solid_feed.mass_frac_comp[0, molecule]) * REE_frac
    #                for molecule, REE_frac in REE_metal_mass_frac.items())
    # print(f"REE feed mass flow is {feed_REE} kg/hr")
    #
    # REE_yield = 100 * product / feed_REE
    # print(f"Yield is {REE_yield} %")

    al_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Al"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.52925 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Al2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Al"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
                * 100
        )
    )

    ca_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Ca"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.7147 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "CaO"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Ca"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    ce_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Ce"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8537 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Ce2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Ce"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    dy_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Dy"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8713 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Dy2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Dy"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    fe_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Fe"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.6994 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Fe2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Fe"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    gd_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Gd"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8676 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Gd2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Gd"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    la_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "La"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8527 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "La2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "La"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    nd_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Nd"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8574 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Nd2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Nd"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    pr_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Pr"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8545 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Pr2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Pr"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    sc_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Sc"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.6520 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Sc2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Sc"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    sm_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Sm"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.8624 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Sm2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Sm"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    yt_recovery = value(
        units.convert(m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Y"]
                      * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0], to_units=units.kg / units.hr)
        / (
                units.convert(0.7874 * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Y2O3"]
                              * m.fs.leach_solid_feed.outlet.flow_mass[0], to_units=units.kg / units.hr)
                + units.convert(m.fs.leach_mixer.outlet.conc_mass_comp[0, "Y"]
                                * m.fs.leach_mixer.outlet.flow_vol[0], to_units=units.kg / units.hr)
        )
        * 100
    )

    print(f"Aluminum recovery is {al_recovery} %")
    print(f"Calcium recovery is {ca_recovery} %")
    print(f"Cerium recovery is {ce_recovery} %")
    print(f"Dysprosium recovery is {dy_recovery} %")
    print(f"Iron recovery is {fe_recovery} %")
    print(f"Gadolinium recovery is {gd_recovery} %")
    print(f"Lanthanum recovery is {la_recovery} %")
    print(f"Neodymium recovery is {nd_recovery} %")
    print(f"Praseodymium  recovery is {pr_recovery} %")
    print(f"Scandium  recovery is {sc_recovery} %")
    print(f"Samarium  recovery is {sm_recovery} %")
    print(f"Yttrium  recovery is {yt_recovery} %")


if __name__ == "__main__":
    m, results = main()
