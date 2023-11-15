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
    Suffix,
    TransformationFactory,
    units,
    Var,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.mscontactor import (
    MSContactor,
    MSContactorInitializer,
)
from idaes.models.unit_models.separator import (
    EnergySplittingType,
)

from idaes.models.unit_models.feed import (
    Feed,
    FeedInitializer,
)
from idaes.models.unit_models.product import (
    Product,
    ProductInitializer,
)

from workspace.leaching.leach_solution_properties import (
    LeachSolutionParameters,
)
from workspace.leaching.leach_solids_properties import (
    CoalRefuseParameters,
)
from workspace.leaching.leach_reactions import (
    CoalRefuseLeachingReactions,
)

from workspace.prommis_workspace.UKy_flowsheet.autoscaling import (
    autoscale_constraints_by_jacobian_norm,
    autoscale_variables_by_magnitude,
)

from workspace.prommis_workspace.Solvent_Extraction.REESXmodel import REESX
from workspace.prommis_workspace.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from workspace.prommis_workspace.Solvent_Extraction.REEOgdistribution import REESolExOgParameters

from workspace.prommis_workspace.precipitate.precipitator import Precipitator
from workspace.prommis_workspace.precipitate.precip_prop import (
    AqueousStateParameterBlock,
    PrecipitateStateParameterBlock,
)

from workspace.prommis_workspace.roasting.ree_oxalate_roaster import REEOxalateRoaster

from workspace.prommis_workspace.UKy_flowsheet.Translators.translator_leaching_SX import (
    Translator_leaching_SX,
)
from workspace.prommis_workspace.UKy_flowsheet.Translators.translator_SX_precipitator import (
    Translator_SX_precipitator,
)
from idaes.models.unit_models.solid_liquid import SLSeparator

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)

from idaes.core.initialization import BlockTriangularizationInitializer

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_diagnostics import DiagnosticsToolbox


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

    # Reactor volume
    m.fs.leach.volume = Var(
        m.fs.time,
        m.fs.leach.elements,
        initialize=1,
        units=units.litre,
        doc="Volume of each finite element."
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

    m.fs.separator1 = SLSeparator(
        solid_property_package=m.fs.coal,
        liquid_property_package=m.fs.leach_soln,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.leach_liquid_feed = Feed(property_package=m.fs.leach_soln)
    m.fs.leach_solid_feed = Feed(property_package=m.fs.coal)
    m.fs.leach_filter_cake = Product(property_package=m.fs.coal)
    m.fs.leach_filter_cake_liquid = Product(property_package=m.fs.leach_soln)

    # ----------------------------------------------------------------------------------------------------------------
    # Solvent extraction property and unit models
    m.fs.prop_a = REESolExAqParameters()
    m.fs.prop_o = REESolExOgParameters()

    m.fs.leach_to_SX = Translator_leaching_SX(
        inlet_property_package=m.fs.leach_soln,
        outlet_property_package=m.fs.prop_a,
    )

    m.fs.solex = REESX(
        number_of_finite_elements=3,
        aqueous_streams={
            "Acidsoln": {"property_package": m.fs.prop_a, "flow_direction": 1}
        },
        organic_streams={
            "Orgacid": {"property_package": m.fs.prop_o, "flow_direction": 2}
        },
    )

    m.fs.sx_leach_acid = Product(property_package=m.fs.prop_a)

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

    m.fs.properties_aq = AqueousStateParameterBlock(
        key_components=key_components,
    )
    m.fs.properties_solid = PrecipitateStateParameterBlock(
        key_components=key_components,
    )

    m.fs.SX_to_precipitator = Translator_SX_precipitator(
        inlet_property_package=m.fs.prop_o,
        outlet_property_package=m.fs.properties_aq,
    )

    m.fs.precipitator = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    m.fs.separator2 = SLSeparator(
        solid_property_package=m.fs.properties_solid,
        liquid_property_package=m.fs.properties_aq,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
    )

    m.fs.precipitate_feed = Feed(property_package=m.fs.properties_solid)

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
    # Connections without recycle loops
    m.fs.liq_feed = Arc(
        source=m.fs.leach_liquid_feed.outlet, destination=m.fs.leach.liquid_inlet
    )
    m.fs.sol_feed = Arc(
        source=m.fs.leach_solid_feed.outlet, destination=m.fs.leach.solid_inlet
    )
    m.fs.s01 = Arc(
        source=m.fs.leach.solid_outlet, destination=m.fs.separator1.solid_inlet
    )
    m.fs.s02 = Arc(source=m.fs.leach.liquid_outlet, destination=m.fs.separator1.liquid_inlet)
    m.fs.sep1_solid = Arc(
        source=m.fs.separator1.solid_outlet, destination=m.fs.leach_filter_cake.inlet
    )
    m.fs.sep1_retained_liquid = Arc(
        source=m.fs.separator1.retained_liquid_outlet, destination=m.fs.leach_filter_cake_liquid.inlet
    )
    m.fs.sep1_liquid = Arc(
        source=m.fs.separator1.recovered_liquid_outlet, destination=m.fs.leach_to_SX.inlet
    )
    m.fs.s03 = Arc(
        source=m.fs.leach_to_SX.outlet, destination=m.fs.solex.Acidsoln_inlet
    )
    m.fs.s04 = Arc(
        source=m.fs.solex.Acidsoln_outlet, destination=m.fs.sx_leach_acid.inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.solex.Orgacid_outlet, destination=m.fs.SX_to_precipitator.inlet
    )
    m.fs.s06 = Arc(
        source=m.fs.SX_to_precipitator.outlet,
        destination=m.fs.precipitator.aqueous_inlet,
    )
    m.fs.s07 = Arc(
        source=m.fs.precipitate_feed.outlet,
        destination=m.fs.precipitator.precipitate_inlet,
    )
    m.fs.s08 = Arc(
        source=m.fs.precipitator.precipitate_outlet, destination=m.fs.separator2.solid_inlet
    )
    m.fs.s09 = Arc(
        source=m.fs.precipitator.aqueous_outlet, destination=m.fs.separator2.liquid_inlet
    )
    m.fs.sep2_solid = Arc(
        source=m.fs.separator2.solid_outlet, destination=m.fs.roaster.solid_inlet
    )
    # TODO: roaster model cannot currently handle liquid inlets
    # m.fs.sep2_retained_liquid = Arc(
    #     source=m.fs.separator2.retained_liquid_outlet, destination=m.fs.precipitate_residual_liquid.inlet
    # )
    # m.fs.sep2_liquid = Arc(
    #     source=m.fs.separator2.recovered_liquid_outlet, destination=m.fs.precipitate_free_liquid.inlet
    # )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_scaling(m):
    # Scaling
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Leaching
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp["Fe"]] = 1e5

    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp["Fe"]] = 1e5

    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.leach_liquid_feed.properties[0].conc_mol_comp["Fe"]] = 1e5

    # Separators
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.separator1.liquid_inlet_state[0].conc_mol_comp["Fe"]] = 1e5

    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.recovered_state[0].conc_mol_comp["Fe"]] = 1e5

    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.separator1.split.retained_state[0].conc_mol_comp["Fe"]] = 1e5

    m.scaling_factor[m.fs.separator2.solid_state[0].temperature] = 1e-2
    m.scaling_factor[m.fs.separator2.liquid_inlet_state[0].temperature] = 1e-2
    m.scaling_factor[m.fs.separator2.split.recovered_state[0].temperature] = 1e-2
    m.scaling_factor[m.fs.separator2.split.retained_state[0].temperature] = 1e-2

    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.leach_filter_cake_liquid.properties[0].conc_mol_comp["Fe"]] = 1e5

    # Translators
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["H2O"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["H"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["HSO4"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["SO4"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Sc"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Y"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["La"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Ce"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Pr"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Nd"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Sm"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Gd"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Dy"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Al"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Ca"]] = 1e5
    m.scaling_factor[m.fs.leach_to_SX.properties_in[0].conc_mol_comp["Fe"]] = 1e5
    m.scaling_factor[m.fs.SX_to_precipitator.properties_out[0].temperature] = 1e-2

    # Precipitator
    m.scaling_factor[m.fs.precipitator.cv_aqueous.properties_in[0].temperature] = 1e-2
    m.scaling_factor[m.fs.precipitator.cv_aqueous.properties_out[0].temperature] = 1e-2
    m.scaling_factor[m.fs.precipitator.cv_precipitate.properties_in[0].temperature] = 1e-2
    m.scaling_factor[m.fs.precipitator.cv_precipitate.properties_out[0].temperature] = 1e-2
    m.scaling_factor[m.fs.precipitate_feed.properties[0].temperature] = 1e-2

    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "Al^3+"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "H^+"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "C2O4^2-"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "Ca^2+"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "Fe^3+"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_key_comp[0, "Ce^3+"]] = 1e10

    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Ce(OH)3(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Ce2(C2O4)3(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Al(OH)3(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Fe2O3(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Ca(C2O4)*H2O(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Ca(C2O4)*3H2O(s)"]] = 1e10
    m.scaling_factor[m.fs.precipitator.molality_precipitate_comp[0, "Ca(OH)2(s)"]] = 1e10

    # Roaster
    m.scaling_factor[m.fs.roaster.gas_in[0].flow_mol] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_in[0].flow_mol_phase["Vap"]] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_in[0].temperature] = 1e-2
    m.scaling_factor[m.fs.roaster.gas_in[0].pressure] = 1e-5

    m.scaling_factor[m.fs.roaster.gas_out[0].flow_mol_phase["Vap"]] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_out[0].flow_mol] = 1e-3
    m.scaling_factor[m.fs.roaster.gas_out[0].temperature] = 1e-2
    m.scaling_factor[m.fs.roaster.gas_out[0].pressure] = 1e-5

    m.scaling_factor[m.fs.roaster.solid_in[0].temperature] = 1e-2

    scaling = TransformationFactory('core.scale_model')
    scaled_model = scaling.create_using(m, rename=False)

    return scaled_model


def set_operating_conditions(m):
    # Liquid feed to leaching unit
    m.fs.leach_liquid_feed.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach_liquid_feed.conc_mass_comp.fix(1e-10 * units.mg / units.L)

    m.fs.leach_liquid_feed.conc_mass_comp[0, "H"].fix(2 * 0.05 * 1e3 * units.mg / units.L)
    m.fs.leach_liquid_feed.conc_mass_comp[0, "HSO4"].fix(1e-8 * units.mg / units.L)
    m.fs.leach_liquid_feed.conc_mass_comp[0, "SO4"].fix(0.05 * 96e3 * units.mg / units.L)

    # Solid feed state
    m.fs.leach_solid_feed.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach_solid_feed.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Sc2O3"].fix(2.77966E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Y2O3"].fix(3.28653E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "La2O3"].fix(6.77769E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Ce2O3"].fix(0.000156161 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Pr2O3"].fix(1.71438E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Nd2O3"].fix(6.76618E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Sm2O3"].fix(1.47926E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Gd2O3"].fix(1.0405E-05 * units.kg / units.kg)
    m.fs.leach_solid_feed.mass_frac_comp[0, "Dy2O3"].fix(7.54827E-06 * units.kg / units.kg)

    m.fs.leach.volume.fix(100 * units.gallon)

    m.fs.separator1.liquid_recovery.fix(0.7)
    m.fs.separator1.split.retained_state[0].hso4_dissociation.deactivate()
    m.fs.separator1.split.recovered_state[0].hso4_dissociation.deactivate()

    # TODO: Replace this with a recycle loop
    # SX rougher recycle stream (treated as a product for now)
    m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01 * units.L / units.hour)

    eps = 1e-7 * units.mg / units.L
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Al"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ca"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Fe"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sc"].fix(321.34 * units.mg / units.L)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Y"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["La"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ce"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Pr"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Nd"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sm"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Gd"].fix(eps)
    m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Dy"].fix(eps)


    # Precipitate feed to precipitator
    # TODO: What should these feed conditions be? They are assumed to be zero in the example
    m.fs.precipitate_feed.properties[0].flow_mol_comp.fix(0)
    m.fs.precipitate_feed.properties[0].temperature.fix(300)

    m.fs.separator2.liquid_recovery.fix(0.7)

    # Roaster gas feed
    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # Inlet flue gas mole flow rate
    fgas = 0.00781
    # Inlet flue gas composition, typical flue gas by burning CH4 with air with stoichiometric ratio 0f 2.3
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
    # Initialize feeds
    initializer_feed = FeedInitializer()
    initializer_product = ProductInitializer()
    initializer1 = MSContactorInitializer()
    initializer2 = BlockTriangularizationInitializer()

    initializer_feed.initialize(m.fs.leach_liquid_feed)
    initializer_feed.initialize(m.fs.leach_solid_feed)
    initializer_feed.initialize(m.fs.precipitate_feed)

    # Initialize leaching section
    propagate_state(m.fs.liq_feed)
    propagate_state(m.fs.sol_feed)

    try:
        initializer1.initialize(
            m.fs.leach,
            initial_guesses={
                "liquid_outlet.flow_vol[0]": 224.477,
                "liquid_outlet.conc_mass_comp[0, Al]": 380.673,
                "liquid_outlet.conc_mass_comp[0, Ca]": 122.972,
                "liquid_outlet.conc_mass_comp[0, Ce]": 5.668,
                "liquid_outlet.conc_mass_comp[0, Dy]": 0.149,
                "liquid_outlet.conc_mass_comp[0, Fe]": 741.241,
                "liquid_outlet.conc_mass_comp[0, Gd]": 0.557,
                "liquid_outlet.conc_mass_comp[0, H]": 2.274,
                "liquid_outlet.conc_mass_comp[0, H2O]": 1000000,
                "liquid_outlet.conc_mass_comp[0, HSO4]": 881.193,
                "liquid_outlet.conc_mass_comp[0, La]": 2.239,
                "liquid_outlet.conc_mass_comp[0, Nd]": 2.684,
                "liquid_outlet.conc_mass_comp[0, Pr]": 0.717,
                "liquid_outlet.conc_mass_comp[0, SO4]": 3924.065,
                "liquid_outlet.conc_mass_comp[0, Sc]": 0.100,
                "liquid_outlet.conc_mass_comp[0, Sm]": 0.301,
                "liquid_outlet.conc_mass_comp[0, Y]": 0.395,
                "solid_outlet.flow_mass[0]": 22.241,
                "solid_outlet.mass_frac_comp[0, Al2O3]": 0.234,
                "solid_outlet.mass_frac_comp[0, CaO]": 0.00164,
                "solid_outlet.mass_frac_comp[0, Ce2O3]": 9.22e-5,
                "solid_outlet.mass_frac_comp[0, Dy2O3]": 5.975e-6,
                "solid_outlet.mass_frac_comp[0, Fe2O3]": 0.0548,
                "solid_outlet.mass_frac_comp[0, Gd2O3]": 4.131e-6,
                "solid_outlet.mass_frac_comp[0, La2O3]": 4.261e-5,
                "solid_outlet.mass_frac_comp[0, Nd2O3]": 3.740e-5,
                "solid_outlet.mass_frac_comp[0, Pr2O3]": 9.016e-6,
                "solid_outlet.mass_frac_comp[0, Sc2O3]": 2.679e-5,
                "solid_outlet.mass_frac_comp[0, Sm2O3]": 1.156e-5,
                "solid_outlet.mass_frac_comp[0, Y2O3]": 2.846e-5,
                "solid_outlet.mass_frac_comp[0, inerts]": 0.7089,
            }
        )
    except:
        # Fix feed states
        m.fs.leach.liquid_inlet.flow_vol.fix()
        m.fs.leach.liquid_inlet.conc_mass_comp.fix()
        m.fs.leach.solid_inlet.flow_mass.fix()
        m.fs.leach.solid_inlet.mass_frac_comp.fix()
        # Re-solve leach unit
        solver = SolverFactory("ipopt")
        solver.solve(m.fs.leach, tee=True)
        # ipopt_solve_halt_on_error(m.fs.leach)
        # Unfix feed states
        m.fs.leach_liquid_feed.flow_vol.unfix()
        m.fs.leach.liquid_inlet.conc_mass_comp.unfix()
        m.fs.leach.solid_inlet.flow_mass.unfix()
        m.fs.leach.solid_inlet.mass_frac_comp.unfix()

    # initialize first separator
    propagate_state(m.fs.s01)
    propagate_state(m.fs.s02)
    initializer2.initialize(m.fs.separator1)

    propagate_state(m.fs.sep1_solid)
    initializer_product.initialize(m.fs.leach_filter_cake)

    propagate_state(m.fs.sep1_retained_liquid)
    initializer_product.initialize(m.fs.leach_filter_cake_liquid)

    # Initialize leaching -> SX translator

    propagate_state(m.fs.sep1_liquid)
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

    # Initialize precipitator
    propagate_state(m.fs.s06)
    propagate_state(m.fs.s07)

    try:
        initializer2.initialize(
            m.fs.precipitator,
            initial_guesses={
                "cv_aqueous.properties_out[0].temperature": 348.15,
                "cv_aqueous.properties_out[0].flow_mass": 0.1147328,
                "cv_aqueous.properties_out[0].pH": 1.699646,
                "cv_aqueous.properties_out[0].log10_molality_comp[Al^3+]": -2.179,
                "cv_aqueous.properties_out[0].log10_molality_comp[C2O4^2-]": -8.344,
                "cv_aqueous.properties_out[0].log10_molality_comp[Ca^2+]": -2.964,
                "cv_aqueous.properties_out[0].log10_molality_comp[Ce^3+]": -8.757,
                "cv_aqueous.properties_out[0].log10_molality_comp[Fe^3+]": -7.764,
                "cv_aqueous.properties_out[0].log10_molality_comp[H2C2O4]": -6.645,
                "cv_aqueous.properties_out[0].log10_molality_comp[HC2O4^-]": -6.090,
                "cv_aqueous.properties_out[0].log10_molality_comp[H^+]": -1.615,
                "cv_aqueous.properties_out[0].log10_molality_comp[OH^-]": -12.212,
                "cv_precipitate.properties_out[0].temperature": 348.15,
            }
        )
    except:
        # Fix feed states
        m.fs.precipitator.cv_aqueous.properties_in[0].flow_mass.fix()
        m.fs.precipitator.cv_aqueous.properties_in[0].log10_molality_comp.fix()
        m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.fix()
        # Re-solve precipitator unit
        solver = SolverFactory("ipopt")
        solver.solve(m.fs.precipitator, tee=True)
        # Unfix feed states
        m.fs.precipitator.cv_aqueous.properties_in[0].flow_mass.unfix()
        m.fs.precipitator.cv_aqueous.properties_in[0].log10_molality_comp.unfix()
        m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.unfix()

    # initialize second separator
    propagate_state(m.fs.s08)
    propagate_state(m.fs.s09)
    initializer2.initialize(m.fs.separator2)

    # initialize roaster
    propagate_state(m.fs.sep2_solid)
    initializer2.initialize(m.fs.roaster)

def solve(m):
    scaling = TransformationFactory('core.scale_model')
    # Scale variables
    autoscale_variables_by_magnitude(m, overwrite=True)
    # Scale constraints
    autoscale_constraints_by_jacobian_norm(m, overwrite=True)
    # Scale model
    autoscaled_model = scaling.create_using(m, rename=False)
    solver = SolverFactory("ipopt")
    solver.solve(autoscaled_model, tee=True)
    results = scaling.propagate_solution(autoscaled_model, m)

    return results


def display_results(m):
    m.fs.roaster.display()

if __name__ == "__main__":
    m, results = main()
