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
r"""
University of Kentucky REE Processing Plant
===========================================

Author: Marcus Holly

The University of Kentucky (UKy) rare earth element (REE) processing plant is designed to extract salable rare earth oxides
from domestic U.S. coal and coal byproducts. While this implementation of the plant does not take into account
all of the complexities and unit processes detailed in the 2019 report, it depicts the major phenomena
by utilizing a series of conventional REE extraction techniques,
including acid tank leaching, solvent extraction, precipitation, and product roasting.

Implementation
--------------

Figure 1 shows the process flow diagram for the simplified, steady-state UKy plant where the solid and liquid feeds are
sent to a leaching tank for initial processing; then the solids (with some retained liquid) are separated out as a filter
cake while the remaining liquid is sent to the solvent extraction process which is comprised of 2 distinct circuits.
In the rougher circuit, solutes in the aqueous phase are transferred to the organic phase, and a portion of the
depleted aqueous solution is recycled back to the leaching process while the remainder is sent to the cleaner circuit.
The degree to which components are transferred from one phase to the other is dependent upon the unit's partition coefficient for that particular component.
In the cleaner circuit, solutes in the organic phase are transferred to the aqueous phase, and a portion of the loaded
aqueous solution is recycled back to the rougher circuit while the remainder is sent to the precipitator. The precipitate
(with some retained liquid) is sent to the roaster where the product rare earth oxides are generated, and the liquid from
the precipitator is recycled back to the cleaner circuit.

.. figure:: ../tutorials/uky_flowsheet.png
    :width: 800
    :align: center

    University of Kentucky flowsheet

Degrees of Freedom
------------------
The following variables must be specified by the user to run the UKy flowsheet:
    * liquid feed volumetric flow rate and component concentrations
    * solid feed mass flow and component mass fractions
    * volume of leach tank(s)
    * partition coefficients for each solvent extraction unit
    * flow rate and component concentrations for organic make-up streams and HCl feeds
    * liquid recovery fraction for solid-liquid separators
    * precipitator inlet temperature
    * roaster inlet temperature and pressure, outlet temperature, and pressure drop
    * roaster gas and moisture molar flows, vapor component mole fractions, and oxide recovery fraction
    * split fractions for each recycle loop

Default Flowsheet Specifications
--------------------------------

===================================================================== ============ ============================
Description                                                           Value        Units
===================================================================== ============ ============================
Leaching
Tank volume                                                           100          :math:`\text{m}^3`
Liquid feed volumetric flow                                           224.3        :math:`\text{L/hr}`
Liquid feed H concentration                                           100          :math:`\text{mg/L}`
Liquid feed HSO4 concentration                                        1e-8         :math:`\text{mg/L}`
Liquid feed SO4 concentration                                         4800         :math:`\text{mg/L}`
Liquid feed REE and contaminant concentrations                        1e-10        :math:`\text{mg/L}`
Solid feed mass flow                                                  22.68        :math:`\text{kg/hr}`
Solid feed inerts mass fraction                                       0.6952       :math:`\text{dimensionless}`
Solid feed Al2O3 mass fraction                                        0.237        :math:`\text{dimensionless}`
Solid feed Fe2O3 mass fraction                                        0.0642       :math:`\text{dimensionless}`
Solid feed CaO mass fraction                                          0.00331      :math:`\text{dimensionless}`
Solid feed Sc2O3 mass fraction                                        2.8e-5       :math:`\text{dimensionless}`
Solid feed Y2O3 mass fraction                                         3.3e-5       :math:`\text{dimensionless}`
Solid feed La2O3 mass fraction                                        6.8e-5       :math:`\text{dimensionless}`
Solid feed Ce2O3 mass fraction                                        1.6e-4       :math:`\text{dimensionless}`
Solid feed Pr2O3 mass fraction                                        1.7e-5       :math:`\text{dimensionless}`
Solid feed Nd2O3 mass fraction                                        6.8e-5       :math:`\text{dimensionless}`
Solid feed Sm2O3 mass fraction                                        1.5e-5       :math:`\text{dimensionless}`
Solid feed Gd2O3 mass fraction                                        1.0e-5       :math:`\text{dimensionless}`
Solid feed Dy2O3 mass fraction                                        7.5e-6       :math:`\text{dimensionless}`

Solvent Extraction Rougher
Loading section organic feed volumetric flow                          62.01        :math:`\text{L/hr}`
Organic make-up REE and contaminant concentrations                    1e-7         :math:`\text{mg/L}`
Scrubbing section acid feed volumetric flow                           0.09         :math:`\text{L/hr}`
Scrubbing section acid feed H concentration                           10.36        :math:`\text{mg/L}`
Scrubbing section acid feed Cl concentration                          359.64       :math:`\text{mg/L}`
Scrubbing section acid feed REE and contaminant concentrations        1e-7         :math:`\text{mg/L}`
Stripping section acid feed volumetric flow                           0.09         :math:`\text{L/hr}`
Stripping section acid feed H concentration                           41.44        :math:`\text{mg/L}`
Stripping section acid feed Cl concentration                          1438.56      :math:`\text{mg/L}`
Stripping section acid feed REE and contaminant concentrations        1e-7         :math:`\text{mg/L}`

Solvent Extraction Cleaner
Loading section organic feed volumetric flow                          62.01        :math:`\text{L/hr}`
Organic make-up REE and contaminant concentrations                    1e-7         :math:`\text{mg/L}`
Stripping section acid feed volumetric flow                           0.09         :math:`\text{L/hr}`
Stripping section acid feed H concentration                           41.44        :math:`\text{mg/L}`
Stripping section acid feed Cl concentration                          1438.56      :math:`\text{mg/L}`
Stripping section acid feed REE and contaminant concentrations        1e-7         :math:`\text{mg/L}`

Precipitator
Inlet temperature                                                     348.15       :math:`\text{K}`

Roaster
Pressure drop                                                         0            :math:`\text{Pa}`
Gas inlet temperature                                                 348.15       :math:`\text{K}`
Gas outlet temperature                                                873.15       :math:`\text{K}`
Gas inlet pressure                                                    101325       :math:`\text{Pa}`
Gas inlet molar flow                                                  0.00781      :math:`\text{mol/s}`
Gas inlet O2 mole fraction                                            0.1118       :math:`\text{dimensionless}`
Gas inlet H2O mole fraction                                           0.1005       :math:`\text{dimensionless}`
Gas inlet CO2 mole fraction                                           0.0431       :math:`\text{dimensionless}`
Gas inlet N2 mole fraction                                            0.7446       :math:`\text{dimensionless}`
Moisture inlet molar flow                                             6.75e-4      :math:`\text{mol/s}`
Oxide recovery fraction                                               0.95         :math:`\text{dimensionless}`

Separators
Leaching solid-liquid separator liquid recovery fraction              0.7          :math:`\text{dimensionless}`
Solvent extraction rougher load recycle split fraction                0.9          :math:`\text{dimensionless}`
Solvent extraction rougher scrub recycle split fraction               0.9          :math:`\text{dimensionless}`
Solvent extraction rougher organic recycle split fraction             0.9          :math:`\text{dimensionless}`
Solvent extraction cleaner organic recycle split fraction             0.9          :math:`\text{dimensionless}`
Precipitator solid-liquid separator liquid recovery fraction          0.7          :math:`\text{dimensionless}`
Precipitator solid-liquid separator liquid recycle split fraction     0.9          :math:`\text{dimensionless}`
===================================================================== ============ ============================

Costing
-------
Unit model costing in this flowsheet is preliminary and is based on the commercial scale unit model parameters provided in Table 4-28 :math:`^1`.
However, this flowsheet is at the pilot scale, so while some of the unit model costing parameters have been scaled down
accordingly, a more robust scale-down procedure of the costing parameters is necessary to accurately approximate the cost of this pilot scale system.


References:

[1] Steven Keim, "Production of salable rare earths products from coal and coal byproducts in the U.S.
using advanced separation processes", 2019

"""


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    SolverFactory,
    Suffix,
    TransformationFactory,
    Var,
    value,
    units,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

import numpy as np

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
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

from prommis.uky.costing.ree_plant_capcost import QGESSCosting, QGESSCostingData


def main():
    """
    Run the flowsheet by calling the appropriate functions in series.
    """
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

    costing = add_costing(scaled_model)

    display_costing(costing)

    return scaled_model, results


def build():
    """
    Build and connect the unit model blocks present in the University of Kentucky REE processing plant.
    """
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
    """
    Set the scaling factors to improve solver performance.
    """

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

    for component in component_set1:
        m.scaling_factor[m.fs.leach.liquid[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )
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
        m.scaling_factor[m.fs.recycle1_purge.properties[0].conc_mol_comp[component]] = (
            1e5
        )
        m.scaling_factor[m.fs.leach_mixer.recycle_state[0].conc_mol_comp[component]] = (
            1e5
        )
        m.scaling_factor[m.fs.leach_mixer.feed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach_mixer.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )

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
        m.scaling_factor[m.fs.recycle2_purge.properties[0].conc_mol_comp[component]] = (
            1e5
        )
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
        m.scaling_factor[m.fs.leach.liquid_inlet_state[0].conc_mol_comp[component]] = (
            1e5
        )

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
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.aqueous_inlet_state[0].flow_vol] = (
        1e-2
    )
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 1].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 2].flow_vol] = 1e-2
    m.scaling_factor[m.fs.solex_cleaner.mscontactor.organic[0, 3].flow_vol] = 1e-2

    m.scaling_factor[m.fs.precipitator.cv_precipitate[0].temperature] = 1e2

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
    """
    Set the operating conditions of the flowsheet such that the degrees of freedom are zero.
    """
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

    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Sc"] = 99.9 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "La"] = 75.2 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Ce"] = 95.7 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Pr"] = 96.5 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Nd"] = 99.2 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Gd"] = 98.6 / 100
    m.fs.solex_rougher.partition_coefficient[:, "aqueous", "organic", "Dy"] = 99.9 / 100

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

    number_of_stages = 3
    stage_number = np.arange(1, number_of_stages + 1)

    for s in stage_number:
        if s == 1:
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Al"] = (
                5.2 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
                3.0 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
                24.7 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
                99.1 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Y"] = (
                100 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "La"] = (
                32.4 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
                58.2 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
                58.2 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
                87.6 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
                100 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
                69.8 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
                96.6 / 100
            )
        else:
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Al"] = (
                4.9 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
                12.3 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
                6.4 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
                16.7 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Y"] = (
                100 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "La"] = (
                23.2 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
                24.9 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
                15.1 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
                100 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
                100 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
                7.6 / 100
            )
            m.fs.solex_cleaner.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
                5.0 / 100
            )

    m.fs.sl_sep1.liquid_recovery.fix(0.7)
    m.fs.sl_sep2.liquid_recovery.fix(0.7)

    m.fs.sep1.split_fraction[:, "recycle"].fix(0.9)
    m.fs.sep2.split_fraction[:, "recycle"].fix(0.9)

    m.fs.precipitator.cv_precipitate[0].temperature.fix(348.15 * units.K)

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

    # Touch properties that are used in the UI
    m.fs.leach.solid_inlet_state[0].flow_mass
    m.fs.leach.solid_inlet_state[0].mass_frac_comp

    m.fs.leach.liquid_inlet_state[0].flow_vol
    m.fs.leach.liquid_inlet_state[0].conc_mol_comp

    m.fs.solex_cleaner.mscontactor.organic_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp

    m.fs.solex_cleaner.mscontactor.aqueous_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp

    m.fs.precipitator.cv_aqueous.properties_out[0].flow_vol
    m.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp

    m.fs.precipitator.cv_precipitate[0].temperature
    m.fs.precipitator.cv_precipitate[0].flow_mol_comp


def initialize_system(m):
    """
    Provide initialized values for all streams in the system.
    """
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
    seq.set_guesses_for(m.fs.precipitator.cv_precipitate[0], tear_guesses2)
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
                # Re-solve precipitator unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.precipitator, tee=True)
                # Unfix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.unfix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.unfix()
        else:
            print(f"Initializing {stream}")
            initializer2.initialize(stream)

    seq.run(m, function)


def solve(m):
    """
    Solve the system with IPOPT.
    """
    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)

    return results


def display_results(m):
    """
    Print key flowsheet outputs.
    """
    m.fs.roaster.display()


def add_costing(flowsheet):
    """
    Set the costing parameters for each unit model.
    """
    # TODO: Costing is preliminary until more unit model costing metrics can be verified
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
    m.fs.costing = QGESSCosting()
    CE_index_year = "UKy_2019"

    # Leaching costs
    # 4.2 is UKy Leaching - Polyethylene Tanks
    L_pe_tanks_accounts = ["4.2"]
    m.fs.L_pe_tanks = UnitModelBlock()
    m.fs.L_pe_tanks.capacity = Var(
        initialize=value(flowsheet.fs.leach.volume[0, 2]), units=units.gal
    )
    m.fs.L_pe_tanks.capacity.fix()
    m.fs.L_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pe_tanks_accounts,
            "scaled_param": m.fs.L_pe_tanks.capacity,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.3 is UKy Leaching - Tank Mixer
    L_tank_mixer_accounts = ["4.3"]
    m.fs.L_tank_mixers = UnitModelBlock()
    m.fs.L_tank_mixers.power = Var(initialize=4.74, units=units.hp)
    m.fs.L_tank_mixers.power.fix()
    m.fs.L_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_tank_mixer_accounts,
            "scaled_param": m.fs.L_tank_mixers.power,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.4 is UKy Leaching - Process Pump
    L_pump_accounts = ["4.4"]
    m.fs.L_pump = UnitModelBlock()
    flow_4_4 = value(
        units.convert(
            flowsheet.fs.leach_liquid_feed.flow_vol[0], to_units=units.gal / units.min
        )
    )
    m.fs.L_pump.feed_rate = Var(initialize=flow_4_4, units=units.gal / units.min)
    m.fs.L_pump.feed_rate.fix()
    m.fs.L_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pump_accounts,
            "scaled_param": m.fs.L_pump.feed_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.5 is UKy Leaching - Thickener
    L_thickener_accounts = ["4.5"]
    m.fs.L_thickener = UnitModelBlock()
    m.fs.L_thickener.area = Var(initialize=225.90, units=units.ft**2)
    m.fs.L_thickener.area.fix()
    m.fs.L_thickener.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_thickener_accounts,
            "scaled_param": m.fs.L_thickener.area,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.6 is UKy Leaching - Solid Waste Filter Press
    L_filter_press_accounts = ["4.6"]
    m.fs.L_filter_press = UnitModelBlock()
    m.fs.L_filter_press.volume = Var(initialize=36.00, units=units.ft**3)
    m.fs.L_filter_press.volume.fix()
    m.fs.L_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_filter_press_accounts,
            "scaled_param": m.fs.L_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.8 is UKy Leaching - Solution Heater
    L_solution_heater_accounts = ["4.8"]
    m.fs.L_solution_heater = UnitModelBlock()
    m.fs.L_solution_heater.duty = Var(initialize=0.24, units=units.MBTU / units.hr)
    m.fs.L_solution_heater.duty.fix()
    m.fs.L_solution_heater.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_solution_heater_accounts,
            "scaled_param": m.fs.L_solution_heater.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # Solvent extraction costs
    # 5.1 is UKy Rougher Solvent Extraction - Polyethylene Tanks
    RSX_pe_tanks_accounts = ["5.1"]
    m.fs.RSX_pe_tanks = UnitModelBlock()
    m.fs.RSX_pe_tanks.capacity = Var(initialize=35.136, units=units.gal)
    m.fs.RSX_pe_tanks.capacity.fix()
    m.fs.RSX_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pe_tanks_accounts,
            "scaled_param": m.fs.RSX_pe_tanks.capacity,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.2 is UKy Rougher Solvent Extraction - Tank Mixer
    RSX_tank_mixer_accounts = ["5.2"]
    m.fs.RSX_tank_mixers = UnitModelBlock()
    m.fs.RSX_tank_mixers.power = Var(initialize=2.0, units=units.hp)
    m.fs.RSX_tank_mixers.power.fix()
    m.fs.RSX_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_tank_mixer_accounts,
            "scaled_param": m.fs.RSX_tank_mixers.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.3 is UKy Rougher Solvent Extraction - Process Pump
    RSX_pump_accounts = ["5.3"]
    m.fs.RSX_pump = UnitModelBlock()
    flow_5_3 = value(
        units.convert(
            flowsheet.fs.solex_rougher.mscontactor.aqueous_inlet.flow_vol[0],
            to_units=units.gal / units.min,
        )
    )
    m.fs.RSX_pump.feed_rate = Var(initialize=flow_5_3, units=units.gal / units.min)
    m.fs.RSX_pump.feed_rate.fix()
    m.fs.RSX_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pump_accounts,
            "scaled_param": m.fs.RSX_pump.feed_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.4 is UKy Rougher Solvent Extraction - Mixer Settler
    RSX_mixer_settler_accounts = ["5.4"]
    m.fs.RSX_mixer_settler = UnitModelBlock()
    m.fs.RSX_mixer_settler.volume = Var(initialize=61.107, units=units.gal)
    m.fs.RSX_mixer_settler.volume.fix()
    m.fs.RSX_mixer_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_mixer_settler_accounts,
            "scaled_param": m.fs.RSX_mixer_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.1 is UKy Cleaner Solvent Extraction - Polyethylene Tanks
    CSX_pe_tanks_accounts = ["6.1"]
    m.fs.CSX_pe_tanks = UnitModelBlock()
    m.fs.CSX_pe_tanks.capacity = Var(initialize=14.05, units=units.gal)
    m.fs.CSX_pe_tanks.capacity.fix()
    m.fs.CSX_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pe_tanks_accounts,
            "scaled_param": m.fs.CSX_pe_tanks.capacity,
            "source": 1,
            "n_equip": 5,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.2 is UKy Cleaner Solvent Extraction - Tank Mixer
    CSX_tank_mixer_accounts = ["6.2"]
    m.fs.CSX_tank_mixers = UnitModelBlock()
    m.fs.CSX_tank_mixers.power = Var(initialize=0.08, units=units.hp)
    m.fs.CSX_tank_mixers.power.fix()
    m.fs.CSX_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_tank_mixer_accounts,
            "scaled_param": m.fs.CSX_tank_mixers.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.3 is UKy Cleaner Solvent Extraction - Process Pump
    CSX_pump_accounts = ["6.3"]
    m.fs.CSX_pump = UnitModelBlock()
    flow_6_3 = value(
        units.convert(
            flowsheet.fs.solex_cleaner.mscontactor.aqueous_inlet.flow_vol[0],
            to_units=units.gal / units.min,
        )
    )
    m.fs.CSX_pump.feed_rate = Var(initialize=flow_6_3, units=units.gal / units.min)
    m.fs.CSX_pump.feed_rate.fix()
    m.fs.CSX_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pump_accounts,
            "scaled_param": m.fs.CSX_pump.feed_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.4 is UKy Cleaner Solvent Extraction - Mixer Settler
    CSX_mixer_settler_accounts = ["6.4"]
    m.fs.CSX_mixer_settler = UnitModelBlock()
    m.fs.CSX_mixer_settler.volume = Var(initialize=24.44, units=units.gal)
    m.fs.CSX_mixer_settler.volume.fix()
    m.fs.CSX_mixer_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_mixer_settler_accounts,
            "scaled_param": m.fs.CSX_mixer_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # Precipitation costs
    # 9.2 is UKy Rare Earth Element Precipitation - Polyethylene Tanks
    reep_pe_tanks_accounts = ["9.2"]
    m.fs.reep_pe_tanks = UnitModelBlock()
    m.fs.reep_pe_tanks.capacity = Var(initialize=15.04, units=units.gal)
    m.fs.reep_pe_tanks.capacity.fix()
    m.fs.reep_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pe_tanks_accounts,
            "scaled_param": m.fs.reep_pe_tanks.capacity,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.3 is UKy Rare Earth Element Precipitation - Tank Mixer
    reep_tank_mixer_accounts = ["9.3"]
    m.fs.reep_tank_mixers = UnitModelBlock()
    m.fs.reep_tank_mixers.power = Var(initialize=0.61, units=units.hp)
    m.fs.reep_tank_mixers.power.fix()
    m.fs.reep_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_tank_mixer_accounts,
            "scaled_param": m.fs.reep_tank_mixers.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.4 is UKy Rare Earth Element Precipitation - Process Pump
    reep_pump_accounts = ["9.4"]
    m.fs.reep_pump = UnitModelBlock()
    flow_9_4 = value(
        units.convert(
            flowsheet.fs.precipitator.aqueous_inlet.flow_vol[0],
            to_units=units.gal / units.min,
        )
    )
    m.fs.reep_pump.feed_rate = Var(initialize=flow_9_4, units=units.gal / units.min)
    m.fs.reep_pump.feed_rate.fix()
    m.fs.reep_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pump_accounts,
            "scaled_param": m.fs.reep_pump.feed_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.5 is UKy Rare Earth Element Precipitation - Filter Press
    reep_filter_press_accounts = ["9.5"]
    m.fs.reep_filter_press = UnitModelBlock()
    m.fs.reep_filter_press.volume = Var(initialize=0.405, units=units.ft**3)
    m.fs.reep_filter_press.volume.fix()
    m.fs.reep_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_filter_press_accounts,
            "scaled_param": m.fs.reep_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.8 is UKy Rare Earth Element Precipitation - Roaster
    reep_roaster_accounts = ["9.8"]
    m.fs.reep_roaster = UnitModelBlock()
    m.fs.reep_roaster.duty = Var(initialize=0.035, units=units.MBTU / units.hr)
    m.fs.reep_roaster.duty.fix()
    m.fs.reep_roaster.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_roaster_accounts,
            "scaled_param": m.fs.reep_roaster.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # Roasting costs
    # 3.1 is UKy Roasting - Storage Bins
    R_storage_bins_accounts = ["3.1"]
    m.fs.R_storage_bins = UnitModelBlock()
    m.fs.R_storage_bins.capacity = Var(initialize=10.0, units=units.ton)
    m.fs.R_storage_bins.capacity.fix()
    m.fs.R_storage_bins.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_storage_bins_accounts,
            "scaled_param": m.fs.R_storage_bins.capacity,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.2 is UKy Roasting - Conveyors
    R_conveyors_accounts = ["3.2"]
    m.fs.R_conveyors = UnitModelBlock()

    flow_3_2 = value(
        units.convert(
            flowsheet.fs.roaster.flow_mass_product[0]
            + flowsheet.fs.roaster.flow_mass_dust[0],
            to_units=units.ton / units.hr,
        )
    )
    m.fs.R_conveyors.throughput = Var(initialize=flow_3_2, units=units.ton / units.hr)
    m.fs.R_conveyors.throughput.fix()
    m.fs.R_conveyors.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_conveyors_accounts,
            "scaled_param": m.fs.R_conveyors.throughput,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.3 is UKy Roasting - Roaster
    R_roaster_accounts = ["3.3"]
    m.fs.R_roaster = UnitModelBlock()
    m.fs.R_roaster.duty = Var(initialize=73.7, units=units.MBTU / units.hr)
    m.fs.R_roaster.duty.fix()
    m.fs.R_roaster.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_roaster_accounts,
            "scaled_param": m.fs.R_roaster.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.4 is UKy Roasting - Gas Scrubber
    R_gas_scrubber_accounts = ["3.4"]
    m.fs.R_gas_scrubber = UnitModelBlock()
    m.fs.R_gas_scrubber.gas_rate = Var(initialize=11.500, units=units.ft**3 / units.min)
    m.fs.R_gas_scrubber.gas_rate.fix()
    m.fs.R_gas_scrubber.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_gas_scrubber_accounts,
            "scaled_param": m.fs.R_gas_scrubber.gas_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.5 is UKy Roasting - Spray Chamber Quencher (7000-60000 ft**3/min)
    R_spray_chamber_quencher_accounts = ["3.5"]
    m.fs.R_spray_chamber_quencher = UnitModelBlock()
    m.fs.R_spray_chamber_quencher.gas_rate = Var(
        initialize=11.500, units=units.ft**3 / units.min
    )
    m.fs.R_spray_chamber_quencher.gas_rate.fix()
    m.fs.R_spray_chamber_quencher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_spray_chamber_quencher_accounts,
            "scaled_param": m.fs.R_spray_chamber_quencher.gas_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.7 is UKy Roasting - Chiller
    R_chiller_accounts = ["3.7"]
    m.fs.R_chiller = UnitModelBlock()
    m.fs.R_chiller.duty = Var(initialize=13.1, units=units.MBTU / units.hr)
    m.fs.R_chiller.duty.fix()
    m.fs.R_chiller.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_chiller_accounts,
            "scaled_param": m.fs.R_chiller.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    feed_input = units.convert(
        flowsheet.fs.leach_solid_feed.flow_mass[0],
        to_units=units.ton / units.hr,
    )

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

    feed_REE = sum(
        flowsheet.fs.leach_solid_feed.flow_mass[0]
        * flowsheet.fs.leach_solid_feed.mass_frac_comp[0, molecule]
        * REE_frac
        for molecule, REE_frac in REE_mass_frac.items()
    )

    feed_grade = (
        units.convert(feed_REE, to_units=units.kg / units.hr)
        / flowsheet.fs.leach_solid_feed.flow_mass[0]
    )

    m.fs.feed_input = Var(initialize=feed_input, units=units.ton / units.hr)
    m.fs.feed_grade = Var(initialize=feed_grade * 1000000, units=units.ppm)

    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    # for convenience
    m.fs.annual_operating_hours = Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=True,
        units=units.hours / units.a,
    )

    recovery_rate = units.convert(
        flowsheet.fs.roaster.flow_mass_product[0], to_units=units.kg / units.hr
    )
    m.fs.recovery_rate_per_year = Var(
        initialize=recovery_rate * m.fs.annual_operating_hours,
        units=units.kg / units.yr,
    )

    # the land cost is the lease cost, or refining cost of REO produced
    m.fs.land_cost = Expression(
        expr=0.303736
        * 1e-6
        * getattr(units, "MUSD_" + CE_index_year)
        / units.ton
        * units.convert(m.fs.feed_input, to_units=units.ton / units.hr)
        * hours_per_shift
        * units.hr
        * shifts_per_day
        * units.day**-1
        * operating_days_per_year
        * units.day
    )

    solid_waste = value(
        units.convert(
            flowsheet.fs.leach_filter_cake.flow_mass[0], to_units=units.ton / units.hr
        )
    )

    m.fs.solid_waste = Var(
        m.fs.time, initialize=solid_waste, units=units.ton / units.hr
    )  # non-hazardous solid waste

    m.fs.precipitate = Var(
        m.fs.time, initialize=0, units=units.ton / units.hr
    )  # non-hazardous precipitate

    dust = value(
        units.convert(
            flowsheet.fs.roaster.flow_mass_dust[0], to_units=units.ton / units.hr
        )
    )
    m.fs.dust_and_volatiles = Var(
        m.fs.time, initialize=dust, units=units.ton / units.hr
    )  # dust and volatiles
    m.fs.power = Var(m.fs.time, initialize=14716, units=units.hp)

    resources = [
        "nonhazardous_solid_waste",
        "nonhazardous_precipitate_waste",
        "dust_and_volatiles",
        "power",
    ]

    rates = [
        m.fs.solid_waste,
        m.fs.precipitate,
        m.fs.dust_and_volatiles,
        m.fs.power,
    ]

    # define product flowrates

    REO_molar_mass = {
        "Y2O3": 88.906 * 2 + 16 * 3,
        "La2O3": 138.91 * 2 + 16 * 3,
        "Ce2O3": 140.12 * 2 + 16 * 3,
        "Pr2O3": 140.91 * 2 + 16 * 3,
        "Nd2O3": 144.24 * 2 + 16 * 3,
        "Sm2O3": 150.36 * 2 + 16 * 3,
        "Gd2O3": 157.25 * 2 + 16 * 3,
        "Dy2O3": 162.5 * 2 + 16 * 3,
        "Sc2O3": 44.96 * 2 + 16 * 3,
    }

    m.fs.Ce_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Ce"]
            * REO_molar_mass["Ce2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product cerium oxide mass flow",
    )

    m.fs.Dy_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Dy"]
            * REO_molar_mass["Dy2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product dysprosium oxide mass flow",
    )

    m.fs.Gd_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Gd"]
            * REO_molar_mass["Gd2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product gadolinium oxide mass flow",
    )

    m.fs.La_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "La"]
            * REO_molar_mass["La2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product lanthanum oxide mass flow",
    )

    m.fs.Nd_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Nd"]
            * REO_molar_mass["Nd2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product neodymium oxide mass flow",
    )

    m.fs.Pr_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Pr"]
            * REO_molar_mass["Pr2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product praseodymium oxide mass flow",
    )

    m.fs.Sc_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Sc"]
            * REO_molar_mass["Sc2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product scandium oxide mass flow",
    )

    m.fs.Sm_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Sm"]
            * REO_molar_mass["Sm2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product samarium oxide mass flow",
    )

    m.fs.Y_product = Param(
        default=units.convert(
            flowsheet.fs.roaster.flow_mol_comp_product[0, "Y"]
            * REO_molar_mass["Y2O3"]
            * units.g
            / units.mol,
            to_units=units.kg / units.hr,
        ),
        units=units.kg / units.hr,
        mutable=True,
        doc="Product yttrium oxide mass flow",
    )

    pure_product_output_rates = {}

    mixed_product_output_rates = {
        "CeO2": m.fs.Ce_product,
        "Sc2O3": m.fs.Sc_product,
        "Y2O3": m.fs.Y_product,
        "La2O3": m.fs.La_product,
        "Nd2O3": m.fs.Nd_product,
        "Pr6O11": m.fs.Pr_product,
        "Sm2O3": m.fs.Sm_product,
        "Gd2O3": m.fs.Gd_product,
        "Dy2O3": m.fs.Dy_product,
    }

    m.fs.costing.build_process_costs(
        # arguments related to installation costs
        piping_materials_and_labor_percentage=20,
        electrical_materials_and_labor_percentage=20,
        instrumentation_percentage=8,
        plants_services_percentage=10,
        process_buildings_percentage=40,
        auxiliary_buildings_percentage=15,
        site_improvements_percentage=10,
        equipment_installation_percentage=17,
        field_expenses_percentage=12,
        project_management_and_construction_percentage=30,
        process_contingency_percentage=15,
        # argument related to Fixed OM costs
        nameplate_capacity=500,  # short (US) ton/hr
        labor_types=[
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ],
        labor_rate=[24.98, 19.08, 30.39, 22.73, 21.97, 45.85],  # USD/hr
        labor_burden=25,  # % fringe benefits
        operators_per_shift=[4, 9, 2, 2, 2, 3],
        hours_per_shift=hours_per_shift,
        shifts_per_day=shifts_per_day,
        operating_days_per_year=operating_days_per_year,
        pure_product_output_rates=pure_product_output_rates,
        mixed_product_output_rates=mixed_product_output_rates,
        mixed_product_sale_price_realization_factor=0.65,  # 65% price realization for mixed products
        # arguments related to total owners costs
        land_cost=m.fs.land_cost,
        resources=resources,
        rates=rates,
        fixed_OM=True,
        variable_OM=True,
        feed_input=m.fs.feed_input,
        efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
        waste=[
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
        ],
        recovery_rate_per_year=m.fs.recovery_rate_per_year,
        CE_index_year=CE_index_year,
    )

    # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
    m.fs.costing.other_plant_costs.unfix()
    m.fs.costing.other_plant_costs_rule = Constraint(
        expr=(
            m.fs.costing.other_plant_costs
            == units.convert(
                1218.073 * units.USD_2016  # Rougher Solvent Extraction
                + 48.723 * units.USD_2016  # Cleaner Solvent Extraction
                + 182.711
                * units.USD_2016,  # Solvent Extraction Wash and Saponification
                to_units=getattr(units, "MUSD_" + CE_index_year),
            )
        )
    )

    # fix costing vars that shouldn't change
    m.fs.feed_input.fix()
    m.fs.feed_grade.fix()
    m.fs.recovery_rate_per_year.fix()
    m.fs.solid_waste.fix()
    m.fs.precipitate.fix()
    m.fs.dust_and_volatiles.fix()
    m.fs.power.fix()

    # check that the model is set up properly and has 0 degrees of freedom
    dt = DiagnosticsToolbox(model=m)
    print("Structural issues in costing")
    dt.report_structural_issues()
    assert degrees_of_freedom(m) == 0

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)

    # Solve costing
    solver = get_solver()
    solver.solve(m, tee=True)

    return m


def display_costing(m):
    """
    Print the key costing results.
    """
    QGESSCostingData.report(m.fs.costing)
    m.fs.costing.variable_operating_costs.display()  # results will be in t = 0
    QGESSCostingData.display_bare_erected_costs(m.fs.costing)
    QGESSCostingData.display_flowsheet_cost(m.fs.costing)


if __name__ == "__main__":
    m, results = main()
