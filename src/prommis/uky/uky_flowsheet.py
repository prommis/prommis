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
from idaes.models.unit_models.mixer import (
    Mixer,
    MixingType,
    MomentumMixingType,
    MixerInitializer,
)
from idaes.models.unit_models.mscontactor import MSContactor, MSContactorInitializer
from idaes.models.unit_models.product import Product, ProductInitializer
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
    SeparatorInitializer,
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
        num_inlets=3,
        inlet_list=["load_recycle", "scrub_recycle", "feed"],
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
    m.fs.prop_o = REESolExOgParameters()

    m.fs.rougher_org_make_up = Feed(property_package=m.fs.prop_o)

    m.fs.solex_rougher_load = SolventExtraction(
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
        aqueous_to_organic=True,
    )

    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Al"] = (
        5.2 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Ca"] = (
        3.0 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Fe"] = (
        24.7 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Sc"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Y"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "La"] = (
        32.4 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Ce"] = (
        58.2 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Pr"] = (
        58.2 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Nd"] = (
        87.6 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Sm"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Gd"] = (
        69.8 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[1, "aqueous", "organic", "Dy"] = (
        96.6 / 100
    )

    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Al"] = (
        4.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Ca"] = (
        12.3 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Fe"] = (
        6.4 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Sc"] = (
        16.7 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Y"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "La"] = (
        23.2 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Ce"] = (
        24.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Pr"] = (
        15.1 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Nd"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Sm"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Gd"] = (
        7.6 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[2, "aqueous", "organic", "Dy"] = (
        5.0 / 100
    )

    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Al"] = (
        4.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Ca"] = (
        12.3 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Fe"] = (
        6.4 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Sc"] = (
        16.7 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Y"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "La"] = (
        23.2 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Ce"] = (
        24.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Pr"] = (
        15.1 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Nd"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Sm"] = (
        99.9 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Gd"] = (
        7.6 / 100
    )
    m.fs.solex_rougher_load.partition_coefficient[3, "aqueous", "organic", "Dy"] = (
        5.0 / 100
    )

    m.fs.acid_feed1 = Feed(property_package=m.fs.leach_soln)

    m.fs.solex_rougher_scrub = SolventExtraction(
        number_of_finite_elements=1,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Al"] = (
        100 - 0.12
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Ca"] = (
        100 - 0.55
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Fe"] = (
        100 - 0.007
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Sc"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Y"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "La"] = (
        100 - 99.8
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Ce"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Pr"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Nd"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Sm"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Gd"] = (
        100 - 99.9
    ) / 100
    m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Dy"] = (
        100 - 99.9
    ) / 100

    m.fs.acid_feed2 = Feed(property_package=m.fs.leach_soln)

    m.fs.solex_rougher_strip = SolventExtraction(
        number_of_finite_elements=2,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Al"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Ca"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Fe"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
        100 - 98.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Y"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "La"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
        100 - 0.5
    ) / 100

    m.fs.rougher_sep = Separator(
        property_package=m.fs.prop_o,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )
    m.fs.rougher_mixer = Mixer(
        property_package=m.fs.prop_o,
        num_inlets=2,
        inlet_list=["make_up", "recycle"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.load_sep = Separator(
        property_package=m.fs.leach_soln,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )
    m.fs.scrub_sep = Separator(
        property_package=m.fs.leach_soln,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.sc_circuit_purge = Product(property_package=m.fs.prop_o)

    m.fs.solex_cleaner_load = SolventExtraction(
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
        aqueous_to_organic=True,
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Al"] = (
        3.6 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Ca"] = (
        3.7 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Fe"] = (
        2.1 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
        99.9 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Y"] = (
        99.9 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "La"] = (
        75.2 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
        95.7 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
        96.5 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
        99.2 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
        99.9 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
        98.6 / 100
    )
    m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
        99.9 / 100
    )

    m.fs.solex_cleaner_strip = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Al"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Ca"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Fe"] = (
        100 - 5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
        100 - 98.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Y"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "La"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
        100 - 0.5
    ) / 100
    m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
        100 - 0.5
    ) / 100

    m.fs.cleaner_org_make_up = Feed(property_package=m.fs.prop_o)

    m.fs.cleaner_mixer = Mixer(
        property_package=m.fs.prop_o,
        num_inlets=2,
        inlet_list=["make_up", "recycle"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.cleaner_sep = Separator(
        property_package=m.fs.prop_o,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.leach_sx_mixer = Mixer(
        property_package=m.fs.leach_soln,
        num_inlets=2,
        inlet_list=["leach", "cleaner"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.acid_feed3 = Feed(property_package=m.fs.leach_soln)
    m.fs.cleaner_purge = Product(property_package=m.fs.prop_o)

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

    m.fs.precipitator = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    m.fs.sl_sep2 = SLSeparator(
        solid_property_package=m.fs.properties_solid,
        liquid_property_package=m.fs.leach_soln,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.precip_sep = Separator(
        property_package=m.fs.leach_soln,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.precip_sx_mixer = Mixer(
        property_package=m.fs.leach_soln,
        num_inlets=2,
        inlet_list=["precip", "rougher"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.precip_purge = Product(property_package=m.fs.properties_aq)
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
    # UKy flowsheet connections
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
        destination=m.fs.leach_sx_mixer.leach,
    )
    m.fs.mixed_aq_feed = Arc(
        source=m.fs.leach_sx_mixer.outlet,
        destination=m.fs.solex_rougher_load.mscontactor.aqueous_inlet,
    )
    m.fs.org_feed = Arc(
        source=m.fs.rougher_org_make_up.outlet, destination=m.fs.rougher_mixer.make_up
    )
    m.fs.mixed_org_feed = Arc(
        source=m.fs.rougher_mixer.outlet,
        destination=m.fs.solex_rougher_load.mscontactor.organic_inlet,
    )
    m.fs.s03 = Arc(
        source=m.fs.solex_rougher_load.mscontactor.aqueous_outlet,
        destination=m.fs.load_sep.inlet,
    )
    m.fs.load_recycle = Arc(
        source=m.fs.load_sep.recycle, destination=m.fs.leach_mixer.load_recycle
    )
    m.fs.s04 = Arc(
        source=m.fs.solex_rougher_load.mscontactor.organic_outlet,
        destination=m.fs.solex_rougher_scrub.mscontactor.organic_inlet,
    )
    m.fs.s05 = Arc(
        source=m.fs.acid_feed1.outlet,
        destination=m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet,
    )
    m.fs.s06 = Arc(
        source=m.fs.solex_rougher_scrub.mscontactor.aqueous_outlet,
        destination=m.fs.scrub_sep.inlet,
    )
    m.fs.scrub_recycle = Arc(
        source=m.fs.scrub_sep.recycle, destination=m.fs.leach_mixer.scrub_recycle
    )
    m.fs.s07 = Arc(
        source=m.fs.solex_rougher_scrub.mscontactor.organic_outlet,
        destination=m.fs.solex_rougher_strip.mscontactor.organic_inlet,
    )
    m.fs.s08 = Arc(
        source=m.fs.acid_feed2.outlet,
        destination=m.fs.solex_rougher_strip.mscontactor.aqueous_inlet,
    )
    m.fs.s09 = Arc(
        source=m.fs.solex_rougher_strip.mscontactor.organic_outlet,
        destination=m.fs.rougher_sep.inlet,
    )
    m.fs.s10 = Arc(
        source=m.fs.rougher_sep.purge, destination=m.fs.sc_circuit_purge.inlet
    )
    m.fs.s11 = Arc(
        source=m.fs.rougher_sep.recycle, destination=m.fs.rougher_mixer.recycle
    )
    m.fs.s12 = Arc(
        source=m.fs.solex_rougher_strip.mscontactor.aqueous_outlet,
        destination=m.fs.precip_sx_mixer.rougher,
    )
    m.fs.s13 = Arc(
        source=m.fs.precip_sx_mixer.outlet,
        destination=m.fs.solex_cleaner_load.mscontactor.aqueous_inlet,
    )
    m.fs.org_feed2 = Arc(
        source=m.fs.cleaner_org_make_up.outlet, destination=m.fs.cleaner_mixer.make_up
    )
    m.fs.s14 = Arc(
        source=m.fs.cleaner_mixer.outlet,
        destination=m.fs.solex_cleaner_load.mscontactor.organic_inlet,
    )
    m.fs.s15 = Arc(
        source=m.fs.solex_cleaner_load.mscontactor.aqueous_outlet,
        destination=m.fs.leach_sx_mixer.cleaner,
    )
    m.fs.s16 = Arc(
        source=m.fs.acid_feed3.outlet,
        destination=m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet,
    )
    m.fs.s17 = Arc(
        source=m.fs.solex_cleaner_load.mscontactor.organic_outlet,
        destination=m.fs.solex_cleaner_strip.mscontactor.organic_inlet,
    )
    m.fs.s18 = Arc(
        source=m.fs.solex_cleaner_strip.mscontactor.organic_outlet,
        destination=m.fs.cleaner_sep.inlet,
    )
    m.fs.s19 = Arc(source=m.fs.cleaner_sep.purge, destination=m.fs.cleaner_purge.inlet)
    m.fs.s20 = Arc(
        source=m.fs.cleaner_sep.recycle, destination=m.fs.cleaner_mixer.recycle
    )
    m.fs.s21 = Arc(
        source=m.fs.solex_cleaner_strip.mscontactor.aqueous_outlet,
        destination=m.fs.precipitator.aqueous_inlet,
    )
    m.fs.s22 = Arc(
        source=m.fs.precipitator.precipitate_outlet,
        destination=m.fs.sl_sep2.solid_inlet,
    )
    m.fs.s23 = Arc(
        source=m.fs.precipitator.aqueous_outlet, destination=m.fs.sl_sep2.liquid_inlet
    )
    m.fs.sep2_solid = Arc(
        source=m.fs.sl_sep2.solid_outlet, destination=m.fs.roaster.solid_inlet
    )
    # # TODO: roaster model cannot currently handle liquid inlets
    # m.fs.sep2_retained_liquid = Arc(
    #     source=m.fs.sl_sep2.retained_liquid_outlet, destination=m.fs.roaster.liquid_inlet
    # )
    m.fs.sep2_recovered_liquid = Arc(
        source=m.fs.sl_sep2.recovered_liquid_outlet, destination=m.fs.precip_sep.inlet
    )
    m.fs.s24 = Arc(source=m.fs.precip_sep.purge, destination=m.fs.precip_purge.inlet)
    m.fs.s25 = Arc(
        source=m.fs.precip_sep.recycle,
        destination=m.fs.precip_sx_mixer.precip,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_scaling(m):
    # Scaling
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    aqueous_component_set = [
        "H2O",
        "H",
        "HSO4",
        "SO4",
        "Cl",
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

    organic_component_set = [
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

    for component in aqueous_component_set:
        m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach.liquid[0, 2].conc_mol_comp[component]] = 1e5
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
        m.scaling_factor[
            m.fs.leach_mixer.load_recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_mixer.scrub_recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.leach_mixer.feed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach_mixer.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_scrub.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5

        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.aqueous[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[m.fs.acid_feed1.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_scrub.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[m.fs.acid_feed2.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.aqueous[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_sx_mixer.leach_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_sx_mixer.cleaner_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach_sx_mixer.mixed_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep2.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep2.split.retained_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sl_sep2.split.recovered_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.load_sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.load_sep.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.load_sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.scrub_sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.scrub_sep.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.scrub_sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.precip_sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.precip_sep.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.precip_sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.precip_purge.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.precip_sx_mixer.precip_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.precip_sx_mixer.rougher_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.precip_sx_mixer.mixed_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.acid_feed3.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.precip_purge.properties[0].conc_mol_comp[component]] = 1
        m.scaling_factor[
            m.fs.precipitator.cv_aqueous.properties_in[0].conc_mol_comp[component]
        ] = 1
        m.scaling_factor[
            m.fs.precipitator.cv_aqueous.properties_out[0].conc_mol_comp[component]
        ] = 1

    for component in organic_component_set:
        m.scaling_factor[
            m.fs.rougher_org_make_up.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.organic[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_load.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_scrub.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.make_up_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.mixed_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.rougher_sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.rougher_sep.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.rougher_sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.make_up_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.rougher_mixer.mixed_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sc_circuit_purge.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.cleaner_mixer.make_up_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.cleaner_mixer.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.cleaner_mixer.mixed_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.sc_circuit_purge.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.organic[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_load.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.cleaner_sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.cleaner_sep.recycle_state[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[m.fs.cleaner_sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[
            m.fs.cleaner_org_make_up.properties[0].conc_mol_comp[component]
        ] = 1e5

        m.scaling_factor[
            m.fs.cleaner_purge.properties[0].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.organic[0, 1].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.organic[0, 2].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.organic[0, 3].conc_mol_comp[component]
        ] = 1e5
        m.scaling_factor[
            m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[0].conc_mol_comp[
                component
            ]
        ] = 1e5

    m.scaling_factor[m.fs.solex_cleaner_load.mscontactor.aqueous[0, 1].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner_load.mscontactor.organic[0, 1].flow_vol] = 1e-2

    m.scaling_factor[m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 1].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 2].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner_strip.mscontactor.aqueous[0, 3].flow_vol] = 1e-2
    m.scaling_factor[
        m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[0].flow_vol
    ] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner_strip.mscontactor.organic[0, 1].flow_vol] = 1e-2

    m.scaling_factor[m.fs.sl_sep2.solid_state[0].temperature] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.liquid_inlet_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.split.recovered_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.sl_sep2.split.retained_state[0].flow_vol] = 1e-2

    m.scaling_factor[m.fs.precip_sep.mixed_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.precip_sep.recycle_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.precip_sep.purge_state[0].flow_vol] = 1e-2
    m.scaling_factor[m.fs.precip_purge.properties[0].flow_vol] = 1e-2

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

    m.fs.load_sep.split_fraction[:, "recycle"].fix(0.9)
    m.fs.scrub_sep.split_fraction[:, "recycle"].fix(0.9)

    m.fs.rougher_org_make_up.flow_vol.fix(6.201)

    m.fs.rougher_org_make_up.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "La"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.acid_feed1.flow_vol.fix(0.09)
    m.fs.acid_feed1.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed1.conc_mass_comp[0, "H"].fix(10.36)
    m.fs.acid_feed1.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Cl"].fix(359.64)
    m.fs.acid_feed1.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Dy"].fix(eps)

    # TODO: flow rate and HCl concentration are not defined in REESim
    m.fs.acid_feed2.flow_vol.fix(0.09)
    m.fs.acid_feed2.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed2.conc_mass_comp[0, "H"].fix(
        10.36 * 4
    )  # Arbitrarily choose 4x the dilute solution
    m.fs.acid_feed2.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Cl"].fix(359.64 * 4)
    m.fs.acid_feed2.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.rougher_sep.split_fraction[:, "recycle"].fix(0.9)

    # TODO: flow rate and HCl concentration are not defined in REESim
    m.fs.acid_feed3.flow_vol.fix(9)
    m.fs.acid_feed3.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed3.conc_mass_comp[0, "H"].fix(
        10.36 * 4
    )  # Arbitrarily choose 4x the dilute solution
    m.fs.acid_feed3.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Cl"].fix(359.64 * 4)
    m.fs.acid_feed3.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.cleaner_org_make_up.flow_vol.fix(6.201)

    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "La"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.cleaner_sep.split_fraction[:, "recycle"].fix(0.9)

    m.fs.sl_sep1.liquid_recovery.fix(0.7)
    m.fs.sl_sep2.liquid_recovery.fix(0.7)

    m.fs.precipitator.cv_precipitate.properties_in[0].temperature.fix(348.15 * units.K)

    m.fs.precip_sep.split_fraction[:, "recycle"].fix(0.9)

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
    seq.options.tear_set = [
        m.fs.feed_mixture,
        m.fs.mixed_aq_feed,
        m.fs.mixed_org_feed,
        m.fs.s13,
        m.fs.s14,
    ]

    G = seq.create_graph(m)
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    tear_guesses1 = {
        "flow_mass": {0.007},
        "conc_mass_comp": {
            "Al": 1493939.39,
            "Ca": 501864.01,
            "Ce": 49698.79,
            "Dy": 466.86,
            "Fe": 1685228.86,
            "Gd": 1624.81,
            "La": 32143.22,
            "Nd": 12552.13,
            "Pr": 4084.40,
            "Sc": 17310.17,
            "Sm": 931.35,
            "Y": 2666.95,
        },
        "flow_mol_comp": {
            "Al": 577.62,
            "Ca": 64.65,
            "Ce": 1.17,
            "Dy": 0.0086,
            "Fe": 231.93,
            "Gd": 0.054,
            "La": 0.63,
            "Nd": 0.34,
            "Pr": 0.099,
            "Sc": 303.33,
            "Sm": 0.021,
            "Y": 0.090,
        },
    }
    tear_guesses2 = {
        "flow_mol_comp": {
            "Al2(C2O4)3(s)": 1.76,
            "Ce2(C2O4)3(s)": 2.65,
            "Dy2(C2O4)3(s)": 0.068,
            "Fe2(C2O4)3(s)": 2.64,
            "Gd2(C2O4)3(s)": 0.27,
            "La2(C2O4)3(s)": 0.86,
            "Nd2(C2O4)3(s)": 1.35,
            "Pr2(C2O4)3(s)": 0.36,
            "Sc2(C2O4)3(s)": 0.62,
            "Sm2(C2O4)3(s)": 0.15,
            "Y2(C2O4)3(s)": 0.31,
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
            (0, "Cl"): 1e-8,
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
        "flow_vol": {0: 62.01},
        "conc_mass_comp": {
            (0, "Al"): 1e-9,
            (0, "Ca"): 1e-9,
            (0, "Ce"): 1e-4,
            (0, "Dy"): 1e-7,
            (0, "Fe"): 1e-7,
            (0, "Gd"): 1e-6,
            (0, "La"): 1e-5,
            (0, "Nd"): 1e-4,
            (0, "Pr"): 1e-6,
            (0, "Sc"): 250,
            (0, "Sm"): 1e-6,
            (0, "Y"): 1e-6,
        },
    }
    tear_guesses5 = {
        "flow_vol": {0: 520},
        "conc_mass_comp": {
            (0, "Al"): 430,
            (0, "Ca"): 99,
            (0, "Ce"): 2,
            (0, "Dy"): 0.01,
            (0, "Fe"): 660,
            (0, "Gd"): 0.1,
            (0, "H"): 2,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 900,
            (0, "Cl"): 0.1,
            (0, "La"): 1,
            (0, "Nd"): 1,
            (0, "Pr"): 0.1,
            (0, "SO4"): 4000,
            (0, "Sc"): 0.05,
            (0, "Sm"): 0.07,
            (0, "Y"): 0.1,
        },
    }
    tear_guesses6 = {
        "flow_vol": {0: 64},
        "conc_mass_comp": {
            (0, "Al"): 1e-9,
            (0, "Ca"): 1e-9,
            (0, "Ce"): 1e-5,
            (0, "Dy"): 1e-7,
            (0, "Fe"): 1e-7,
            (0, "Gd"): 1e-6,
            (0, "La"): 1e-5,
            (0, "Nd"): 1e-5,
            (0, "Pr"): 1e-6,
            (0, "Sc"): 321.34,
            (0, "Sm"): 1e-6,
            (0, "Y"): 1e-6,
        },
    }
    tear_guesses7 = {
        "flow_vol": {0: 5.7},
        "conc_mass_comp": {
            (0, "Al"): 5,
            (0, "Ca"): 16,
            (0, "Ce"): 346,
            (0, "Dy"): 6,
            (0, "Fe"): 1,
            (0, "Gd"): 22,
            (0, "H"): 14,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 1e-7,
            (0, "Cl"): 1400,
            (0, "La"): 160,
            (0, "Nd"): 121,
            (0, "Pr"): 30,
            (0, "SO4"): 1e-7,
            (0, "Sc"): 149.2,
            (0, "Sm"): 13,
            (0, "Y"): 18,
        },
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.precipitator.cv_aqueous.properties_out[0], tear_guesses1)
    seq.set_guesses_for(
        m.fs.precipitator.cv_precipitate.properties_out[0], tear_guesses2
    )
    seq.set_guesses_for(m.fs.leach.liquid_inlet, tear_guesses3)
    seq.set_guesses_for(
        m.fs.solex_rougher_load.mscontactor.organic_inlet, tear_guesses4
    )
    seq.set_guesses_for(
        m.fs.solex_rougher_load.mscontactor.aqueous_inlet, tear_guesses5
    )
    seq.set_guesses_for(
        m.fs.solex_cleaner_load.mscontactor.organic_inlet, tear_guesses6
    )
    seq.set_guesses_for(
        m.fs.solex_cleaner_load.mscontactor.aqueous_inlet, tear_guesses7
    )
    seq.set_guesses_for(m.fs.precip_sx_mixer.outlet, tear_guesses7)

    def function(stream):
        initializer_feed = FeedInitializer()
        initializer_product = ProductInitializer()
        initializer_sep = SeparatorInitializer()
        initializer_mix = MixerInitializer()

        initializer1 = MSContactorInitializer()
        initializer2 = BlockTriangularizationInitializer()

        propagate_state(m.fs.liq_feed)
        propagate_state(m.fs.sol_feed)
        propagate_state(m.fs.org_feed)
        propagate_state(m.fs.org_feed2)

        if stream == m.fs.leach_liquid_feed:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.leach_liquid_feed)
        elif stream == m.fs.leach_solid_feed:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.leach_solid_feed)
        elif stream == m.fs.rougher_org_make_up:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.rougher_org_make_up)
        elif stream == m.fs.acid_feed1:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.acid_feed1)
        elif stream == m.fs.acid_feed2:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.acid_feed2)
        elif stream == m.fs.acid_feed3:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.acid_feed3)
        elif stream == m.fs.cleaner_org_make_up:
            print(f"Initializing {stream}")
            initializer_feed.initialize(m.fs.cleaner_org_make_up)
        elif stream == m.fs.leach_filter_cake:
            print(f"Initializing {stream}")
            initializer_product.initialize(m.fs.leach_filter_cake)
        elif stream == m.fs.load_sep:
            print(f"Initializing {stream}")
            initializer_sep.initialize(m.fs.load_sep)
        elif stream == m.fs.scrub_sep:
            print(f"Initializing {stream}")
            initializer_sep.initialize(m.fs.scrub_sep)
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
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.leach, tee=True)
                # Unfix feed states
                m.fs.leach.liquid_inlet.flow_vol.unfix()
                m.fs.leach.liquid_inlet.conc_mass_comp.unfix()
                m.fs.leach.solid_inlet.flow_mass.unfix()
                m.fs.leach.solid_inlet.mass_frac_comp.unfix()
        elif stream == m.fs.leach_mixer:
            print(f"Initializing {stream}")
            initializer2.initialize(m.fs.leach_mixer)
        elif stream == m.fs.solex_rougher_load.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.solex_rougher_load)
            except:
                # Fix feed states
                m.fs.solex_rougher_load.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_load.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.fix()
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_rougher_load, tee=True)
                # Unfix feed states
                m.fs.solex_rougher_load.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_load.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.unfix()
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.unfix()
        elif stream == m.fs.solex_rougher_scrub.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.solex_rougher_scrub)
            except:
                # Fix feed states
                m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.fix()
                m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_rougher_scrub, tee=True)
                # Unfix feed states
                m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.unfix()
                m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.unfix()
        elif stream == m.fs.solex_rougher_strip.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.solex_rougher_strip)
            except:
                # Fix feed states
                m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.fix()
                m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_rougher_strip, tee=True)
                # Unfix feed states
                m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.unfix()
                m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.unfix()
        elif stream == m.fs.solex_cleaner_load.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.solex_cleaner_load)
            except:
                # Fix feed states
                m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.fix()
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_cleaner_load, tee=True)
                # Unfix feed states
                m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.unfix()
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.unfix()
        elif stream == m.fs.solex_cleaner_strip.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.solex_cleaner_strip)
            except:
                # Fix feed states
                m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.fix()
                m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.fix()
                m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_cleaner_strip, tee=True)
                # Unfix feed states
                m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[
                    0
                ].flow_vol.unfix()
                m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[
                    0
                ].conc_mass_comp.unfix()
                m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[
                    0
                ].conc_mass_comp.unfix()
        elif stream == m.fs.precipitator:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.precipitator)
            except:
                # Fix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.fix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.fix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.precipitator, tee=True)
                # Unfix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.unfix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.unfix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.unfix()
        elif stream == m.fs.sl_sep2:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.sl_sep2)
            except:
                # Fix feed states
                m.fs.sl_sep2.liquid_inlet_state[0].flow_vol.fix()
                m.fs.sl_sep2.liquid_inlet_state[0].conc_mass_comp.fix()
                m.fs.sl_sep2.solid_state[0].flow_mol_comp.fix()
                # Re-solve unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.sl_sep2, tee=True)
                # Unfix feed states
                m.fs.sl_sep2.liquid_inlet_state[0].flow_vol.unfix()
                m.fs.sl_sep2.liquid_inlet_state[0].conc_mass_comp.unfix()
                m.fs.sl_sep2.solid_state[0].flow_mol_comp.unfix()
        elif stream == m.fs.precip_sep:
            print(f"Initializing {stream}")
            initializer_sep.initialize(m.fs.precip_sep)
        elif stream == m.fs.precip_sx_mixer:
            print(f"Initializing {stream}")
            initializer_mix.initialize(m.fs.precip_sx_mixer)
        else:
            print(f"Initializing {stream}")
            initializer2.initialize(stream)

    seq.run(m, function)


def solve(m):
    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    m.fs.rougher_org_make_up.outlet.flow_vol.unfix()
    m.fs.rougher_mixer.outlet.flow_vol.fix(62.01)

    m.fs.cleaner_org_make_up.outlet.flow_vol.unfix()
    m.fs.cleaner_mixer.outlet.flow_vol.fix(62.01)

    results = solver.solve(m, tee=True)

    return results


def display_results(m):
    m.fs.solex_rougher_load.display()
    m.fs.solex_cleaner_load.display()

    metal_mass_frac = {
        "Al2O3": 26.98 * 2 / (26.98 * 2 + 16 * 3),
        "Fe2O3": 55.845 * 2 / (55.845 * 2 + 16 * 3),
        "CaO": 40.078 / (40.078 + 16),
        "Sc2O3": 44.956 * 2 / (44.956 * 2 + 16 * 3),
        "Y2O3": 88.906 * 2 / (88.906 * 2 + 16 * 3),
        "La2O3": 138.91 * 2 / (138.91 * 2 + 16 * 3),
        "Ce2O3": 140.12 * 2 / (140.12 * 2 + 16 * 3),
        "Pr2O3": 140.91 * 2 / (140.91 * 2 + 16 * 3),
        "Nd2O3": 144.24 * 2 / (144.24 * 2 + 16 * 3),
        "Sm2O3": 150.36 * 2 / (150.36 * 2 + 16 * 3),
        "Gd2O3": 157.25 * 2 / (157.25 * 2 + 16 * 3),
        "Dy2O3": 162.5 * 2 / (162.5 * 2 + 16 * 3),
    }

    molar_mass = {
        "Al2O3": (26.98 * 2 + 16 * 3) * units.g / units.mol,
        "Fe2O3": (55.845 * 2 + 16 * 3) * units.g / units.mol,
        "CaO": (40.078 + 16) * units.g / units.mol,
        "Sc2O3": (44.956 * 2 + 16 * 3) * units.g / units.mol,
        "Y2O3": (88.906 * 2 + 16 * 3) * units.g / units.mol,
        "La2O3": (138.91 * 2 + 16 * 3) * units.g / units.mol,
        "Ce2O3": (140.12 * 2 + 16 * 3) * units.g / units.mol,
        "Pr2O3": (140.91 * 2 + 16 * 3) * units.g / units.mol,
        "Nd2O3": (144.24 * 2 + 16 * 3) * units.g / units.mol,
        "Sm2O3": (150.36 * 2 + 16 * 3) * units.g / units.mol,
        "Gd2O3": (157.25 * 2 + 16 * 3) * units.g / units.mol,
        "Dy2O3": (162.5 * 2 + 16 * 3) * units.g / units.mol,
    }

    REE_mass_frac = {
        "Y2O3": 88.906 * 2 / (88.906 * 2 + 16 * 3),
        "La2O3": 138.91 * 2 / (138.91 * 2 + 16 * 3),
        "Ce2O3": 140.12 * 2 / (140.12 * 2 + 16 * 3),
        "Pr2O3": 140.91 * 2 / (140.91 * 2 + 16 * 3),
        "Nd2O3": 144.24 * 2 / (144.24 * 2 + 16 * 3),
        "Sm2O3": 150.36 * 2 / (150.36 * 2 + 16 * 3),
        "Gd2O3": 157.25 * 2 / (157.25 * 2 + 16 * 3),
        "Dy2O3": 162.5 * 2 / (162.5 * 2 + 16 * 3),
    }

    # Total mass basis yield calculation
    product = value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Y"]
            * molar_mass["Y2O3"]
            * REE_mass_frac["Y2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "La"]
            * molar_mass["La2O3"]
            * REE_mass_frac["La2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Ce"]
            * molar_mass["Ce2O3"]
            * REE_mass_frac["Ce2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Pr"]
            * molar_mass["Pr2O3"]
            * REE_mass_frac["Pr2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Nd"]
            * molar_mass["Nd2O3"]
            * REE_mass_frac["Nd2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Sm"]
            * molar_mass["Sm2O3"]
            * REE_mass_frac["Sm2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Gd"]
            * molar_mass["Gd2O3"]
            * REE_mass_frac["Gd2O3"],
            to_units=units.kg / units.hr,
        )
        + units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Dy"]
            * molar_mass["Dy2O3"]
            * REE_mass_frac["Dy2O3"],
            to_units=units.kg / units.hr,
        )
    )
    print(f"REE product mass flow is {product} kg/hr")
    feed_REE = sum(
        value(
            m.fs.leach_solid_feed.flow_mass[0]
            * m.fs.leach_solid_feed.mass_frac_comp[0, molecule]
        )
        * REE_frac
        for molecule, REE_frac in REE_mass_frac.items()
    )
    print(f"REE feed mass flow is {feed_REE} kg/hr")

    REE_recovery = 100 * product / feed_REE
    print(f"Total REE recovery is {REE_recovery} %")

    product_purity = 100 * product / value(units.convert(m.fs.roaster.flow_mas_product[0], to_units=units.kg / units.hr))
    print(f"Product purity is {product_purity} % REE")

    # Individual elemental recoveries

    total_al_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Al"]
            * molar_mass["Al2O3"]
            * metal_mass_frac["Al2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Al2O3"]
                * metal_mass_frac["Al2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_fe_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Fe"]
            * molar_mass["Fe2O3"]
            * metal_mass_frac["Fe2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Fe2O3"]
                * metal_mass_frac["Fe2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    # total_ca_recovery = 100 * value(
    #     units.convert(
    #         m.fs.roaster.flow_mol_comp_product[0, "Ca"]
    #         * molar_mass["CaO"]
    #         * metal_mass_frac["CaO"],
    #         to_units=units.kg / units.hr,
    #     )
    #     / (
    #         units.convert(
    #             m.fs.leach_solid_feed.flow_mass[0]
    #             * m.fs.leach_solid_feed.mass_frac_comp[0, "Cao"]
    #             * metal_mass_frac["CaO"],
    #             to_units=units.kg / units.hr,
    #         )
    #     )
    # )

    total_sc_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Sc"]
            * molar_mass["Sc2O3"]
            * metal_mass_frac["Sc2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Sc2O3"]
                * metal_mass_frac["Sc2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_yt_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Y"]
            * molar_mass["Y2O3"]
            * metal_mass_frac["Y2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Y2O3"]
                * metal_mass_frac["Y2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_la_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "La"]
            * molar_mass["La2O3"]
            * metal_mass_frac["La2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "La2O3"]
                * metal_mass_frac["La2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_ce_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Ce"]
            * molar_mass["Ce2O3"]
            * metal_mass_frac["Ce2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Ce2O3"]
                * metal_mass_frac["Ce2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_pr_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Pr"]
            * molar_mass["Pr2O3"]
            * metal_mass_frac["Pr2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Pr2O3"]
                * metal_mass_frac["Pr2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_nd_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Nd"]
            * molar_mass["Nd2O3"]
            * metal_mass_frac["Nd2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Nd2O3"]
                * metal_mass_frac["Nd2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_sm_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Sm"]
            * molar_mass["Sm2O3"]
            * metal_mass_frac["Sm2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Sm2O3"]
                * metal_mass_frac["Sm2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_gd_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Gd"]
            * molar_mass["Gd2O3"]
            * metal_mass_frac["Gd2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Gd2O3"]
                * metal_mass_frac["Gd2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    total_dy_recovery = 100 * value(
        units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Dy"]
            * molar_mass["Dy2O3"]
            * metal_mass_frac["Dy2O3"],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                m.fs.leach_solid_feed.flow_mass[0]
                * m.fs.leach_solid_feed.mass_frac_comp[0, "Dy2O3"]
                * metal_mass_frac["Dy2O3"],
                to_units=units.kg / units.hr,
            )
        )
    )

    al_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Al"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Al2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Al2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    # ca_recovery = value(
    #     units.convert(
    #         m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Ca"]
    #         * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
    #         to_units=units.kg / units.hr,
    #     )
    #     / (
    #         units.convert(
    #             metal_mass_frac["CaO"]
    #             * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Cao"]
    #             * m.fs.leach_solid_feed.outlet.flow_mass[0],
    #             to_units=units.kg / units.hr,
    #         )
    #     )
    #     * 100
    # )

    ce_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Ce"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Ce2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Ce2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    dy_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Dy"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Dy2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Dy2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    fe_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Fe"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Fe2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Fe2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    gd_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Gd"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Gd2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Gd2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    la_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "La"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["La2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "La2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    nd_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Nd"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Nd2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Nd2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    pr_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Pr"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Pr2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Pr2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    sc_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Sc"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Sc2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Sc2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    sm_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Sm"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Sm2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Sm2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    yt_recovery = value(
        units.convert(
            m.fs.sl_sep1.recovered_liquid_outlet.conc_mass_comp[0, "Y"]
            * m.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0],
            to_units=units.kg / units.hr,
        )
        / (
            units.convert(
                metal_mass_frac["Y2O3"]
                * m.fs.leach_solid_feed.outlet.mass_frac_comp[0, "Y2O3"]
                * m.fs.leach_solid_feed.outlet.flow_mass[0],
                to_units=units.kg / units.hr,
            )
        )
        * 100
    )

    print(f"\nLeaching Aluminum recovery is {al_recovery} %")
    print(f"Total aluminum recovery is {total_al_recovery} %")
    # print(f"\nLeaching Calcium recovery is {ca_recovery} %")
    # print(f"Total Calcium recovery is {total_ca_recovery} %")
    print(f"\nLeaching Cerium recovery is {ce_recovery} %")
    print(f"Total Cerium recovery is {total_ce_recovery} %")
    print(f"\nLeaching Dysprosium recovery is {dy_recovery} %")
    print(f"Total Dysprosium recovery is {total_dy_recovery} %")
    print(f"\nLeaching Iron recovery is {fe_recovery} %")
    print(f"Total Iron recovery is {total_fe_recovery} %")
    print(f"\nLeaching Gadolinium recovery is {gd_recovery} %")
    print(f"Total Gadolinium recovery is {total_gd_recovery} %")
    print(f"\nLeaching Lanthanum recovery is {la_recovery} %")
    print(f"Total Lanthanum recovery is {total_la_recovery} %")
    print(f"\nLeaching Neodymium recovery is {nd_recovery} %")
    print(f"Total Neodymium recovery is {total_nd_recovery} %")
    print(f"\nLeaching Praseodymium recovery is {pr_recovery} %")
    print(f"Total Praseodymium recovery is {total_pr_recovery} %")
    print(f"\nLeaching Scandium recovery is {sc_recovery} %")
    print(f"Total Scandium recovery is {total_sc_recovery} %")
    print(f"\nLeaching Samarium recovery is {sm_recovery} %")
    print(f"Total Samarium recovery is {total_sm_recovery} %")
    print(f"\nLeaching Yttrium recovery is {yt_recovery} %")
    print(f"Total Yttrium recovery is {total_yt_recovery} %")


if __name__ == "__main__":
    m, results = main()
