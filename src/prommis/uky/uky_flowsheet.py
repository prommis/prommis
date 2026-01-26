#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
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
Unit model costing in this flowsheet is based on the commercial scale unit model parameters provided in Table 4-28 :math:`^1`.
The reference cost and capacity parameter data from [1] are at commercial scale. As this flowsheet is at the pilot scale, some
of the unit model capacity parameters have been scaled down accordingly by unit feed rate.


References:

[1] Keim, Steven Anthony and Naumann, Hans. "Production of Salable Rare Earths Products from Coal and Coal Byproducts
 in the U.S. Using Advanced Separation Processes (Final Technical Report)." , Sep. 2019. https://doi.org/10.2172/1569277

"""
import logging
from warnings import warn
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Param,
    Set,
    Suffix,
    TransformationFactory,
    Var,
    check_optimal_termination,
    units,
    value,
)
from pyomo.network import Arc, SequentialDecomposition

import idaes.logger as idaeslog
from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlock,
    UnitModelBlockData,
    UnitModelCostingBlock,
)
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.scaling import AutoScaler, CustomScalerBase, ConstraintScalingScheme
from idaes.core.scaling.util import get_scaling_factor, set_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    unscaled_variables_generator,
    unscaled_constraints_generator,
)

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    ModularPropertiesScaler,
)
from idaes.models.unit_models.feed import Feed, FeedInitializer
from idaes.models.unit_models.mixer import (
    Mixer,
    MixerInitializer,
    MixingType,
    MomentumMixingType,
)
from idaes.models.unit_models.product import Product, ProductInitializer
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SeparatorInitializer,
    SplittingType,
)
from idaes.models.unit_models.solid_liquid import SLSeparator
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
from prommis.leaching.leach_train import LeachingTrain, LeachingTrainInitializer
from prommis.properties import HClStrippingParameterBlock
from prommis.properties.hcl_stripping_properties import HClStrippingPropertiesScaler
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
)

# from prommis.solvent_extraction.translator_leach_precip import TranslatorLeachPrecip
from prommis.properties.translator_hcl_leach import TranslatorHClLeach
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)
from prommis.uky.costing.costing_dictionaries import load_REE_costing_dictionary
from prommis.uky.costing.ree_plant_capcost import QGESSCosting, QGESSCostingData

_log = idaeslog.getLogger(__name__)

# Epsilon represents near-zero component concentrations
eps = 1e-8 * units.mg / units.L


def main():
    """
    Run the flowsheet by calling the appropriate functions in series.
    """
    m = build()

    set_operating_conditions(m)

    set_scaling(m)

    if degrees_of_freedom(m) != 0:
        raise AssertionError(
            "The degrees of freedom are not equal to 0."
            "Check that the expected variables are fixed and unfixed."
            "For more guidance, run assert_no_structural_warnings from the IDAES DiagnosticToolbox "
        )

    initialize_system(m)

    solve_system(m, tee=True)

    # fixes the volumetric flow rate of the organic recycle streams and unfixes the flow of the make-up streams
    # we want to be able to adjust the total recycle flow rate, not just the make-up portion of it
    fix_organic_recycle(m)

    results = solve_system(m, tee=True)

    if not check_optimal_termination(results):
        raise RuntimeError(
            "Solver failed to terminate with an optimal solution. Please check the solver logs for more details"
        )
    add_result_expressions(m)
    display_results(m)

    add_costing(m)
    initialize_costing(m)

    # diagnostics, initialize, and solve
    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    auto = AutoScaler()
    auto.scale_variables_by_magnitude(m)
    auto.scale_constraints_by_jacobian_norm(m)

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)

    # from idaes.core.scaling.util import jacobian_cond

    # print(jacobian_cond(m))

    solve_system(m, tee=True)

    dt.assert_no_numerical_warnings()

    display_costing(m)

    return m, results


def build():
    """
    Build and connect the unit model blocks present in the University of Kentucky REE processing plant.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Leaching property and unit models
    m.fs.leach_soln = SulfuricAcidLeachingParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactionParameterBlock()
    m.fs.HCl_stripping_params = HClStrippingParameterBlock()

    m.fs.leach = LeachingTrain(
        number_of_tanks=2,
        liquid_phase={
            "property_package": m.fs.leach_soln,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        solid_phase={
            "property_package": m.fs.coal,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        reaction_package=m.fs.leach_rxns,
    )

    m.fs.sl_sep1 = SLSeparator(
        solid_property_package=m.fs.coal,
        liquid_property_package=m.fs.leach_soln,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.scrubber_HCl_leach_translator = TranslatorHClLeach(
        inlet_property_package=m.fs.HCl_stripping_params,
        outlet_property_package=m.fs.leach_soln,
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
    # Solvent extraction property, reaction and unit models
    m.fs.prop_o = REESolExOgParameters()
    m.fs.reaxn = SolventExtractionReactions()

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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=False,
        create_hydrostatic_pressure_terms=False,
    )

    m.fs.acid_feed1 = Feed(property_package=m.fs.HCl_stripping_params)

    m.fs.solex_rougher_scrub = SolventExtraction(
        number_of_finite_elements=1,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.HCl_stripping_params,
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=False,
        create_hydrostatic_pressure_terms=False,
    )

    m.fs.acid_feed2 = Feed(property_package=m.fs.HCl_stripping_params)

    m.fs.solex_rougher_strip = SolventExtraction(
        number_of_finite_elements=2,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.HCl_stripping_params,
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=False,
        create_hydrostatic_pressure_terms=False,
    )

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
        property_package=m.fs.HCl_stripping_params,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.rougher_organic_purge = Product(property_package=m.fs.prop_o)

    m.fs.solex_cleaner_load = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.HCl_stripping_params,
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=False,
        create_hydrostatic_pressure_terms=False,
    )

    m.fs.solex_cleaner_strip = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.HCl_stripping_params,
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=False,
        create_hydrostatic_pressure_terms=False,
    )

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

    m.fs.cleaner_HCl_leach_translator = TranslatorHClLeach(
        inlet_property_package=m.fs.HCl_stripping_params,
        outlet_property_package=m.fs.leach_soln,
    )

    m.fs.leach_sx_mixer = Mixer(
        property_package=m.fs.leach_soln,
        num_inlets=2,
        inlet_list=["leach", "cleaner"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.acid_feed3 = Feed(property_package=m.fs.HCl_stripping_params)
    m.fs.cleaner_organic_purge = Product(property_package=m.fs.prop_o)

    # --------------------------------------------------------------------------------------------------------------
    # Precipitation property and unit models

    m.fs.properties_solid = PrecipitateParameters()

    m.fs.precipitator = Precipitator(
        property_package_aqueous=m.fs.HCl_stripping_params,
        property_package_precipitate=m.fs.properties_solid,
        make_volume_balance_constraint=False,
    )

    m.fs.sl_sep2 = SLSeparator(
        solid_property_package=m.fs.properties_solid,
        liquid_property_package=m.fs.HCl_stripping_params,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.precip_sep = Separator(
        property_package=m.fs.HCl_stripping_params,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.precip_sx_mixer = Mixer(
        property_package=m.fs.HCl_stripping_params,
        num_inlets=2,
        inlet_list=["precip", "rougher"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.precip_purge = Product(property_package=m.fs.HCl_stripping_params)
    # -----------------------------------------------------------------------------------------------------------------
    # Roasting property and unit models

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )

    m.fs.prop_solid = PrecipitateParameters()

    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate_solid=m.fs.prop_solid,
        property_package_precipitate_liquid=m.fs.HCl_stripping_params,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # -----------------------------------------------------------------------------------------------------------------
    # Translator blocks

    # m.fs.translator_leaching_to_precipitate = TranslatorLeachPrecip(
    #     inlet_property_package=m.fs.leach_soln,
    #     outlet_property_package=m.fs.HCl_stripping_params,
    # )
    # m.fs.translator_precipitate_to_leaching = TranslatorLeachPrecip(
    #     inlet_property_package=m.fs.HCl_stripping_params,
    #     outlet_property_package=m.fs.leach_soln,
    # )
    # m.fs.translator_sep_to_roast = TranslatorLeachPrecip(
    #     inlet_property_package=m.fs.leach_soln,
    #     outlet_property_package=m.fs.HCl_stripping_params,
    # )
    # m.fs.translator_precip_sep_to_purge = TranslatorLeachPrecip(
    #     inlet_property_package=m.fs.leach_soln,
    #     outlet_property_package=m.fs.HCl_stripping_params,
    # )

    # -----------------------------------------------------------------------------------------------------------------

    # UKy flowsheet connections
    m.fs.leaching_sol_feed = Arc(
        source=m.fs.leach_solid_feed.outlet, destination=m.fs.leach.solid_inlet
    )
    m.fs.leaching_liq_feed = Arc(
        source=m.fs.leach_liquid_feed.outlet, destination=m.fs.leach_mixer.feed
    )
    m.fs.leaching_feed_mixture = Arc(
        source=m.fs.leach_mixer.outlet, destination=m.fs.leach.liquid_inlet
    )
    m.fs.leaching_solid_outlet = Arc(
        source=m.fs.leach.solid_outlet, destination=m.fs.sl_sep1.solid_inlet
    )
    m.fs.leaching_liquid_outlet = Arc(
        source=m.fs.leach.liquid_outlet, destination=m.fs.sl_sep1.liquid_inlet
    )
    m.fs.sl_sep1_solid_outlet = Arc(
        source=m.fs.sl_sep1.solid_outlet, destination=m.fs.leach_filter_cake.inlet
    )
    m.fs.sl_sep1_retained_liquid_outlet = Arc(
        source=m.fs.sl_sep1.retained_liquid_outlet,
        destination=m.fs.leach_filter_cake_liquid.inlet,
    )
    m.fs.sl_sep1_liquid_outlet = Arc(
        source=m.fs.sl_sep1.recovered_liquid_outlet,
        destination=m.fs.leach_sx_mixer.leach,
    )
    m.fs.sx_rougher_load_aq_feed = Arc(
        source=m.fs.leach_sx_mixer.outlet,
        destination=m.fs.solex_rougher_load.aqueous_inlet,
    )
    m.fs.sx_rougher_org_feed = Arc(
        source=m.fs.rougher_org_make_up.outlet, destination=m.fs.rougher_mixer.make_up
    )
    m.fs.sx_rougher_mixed_org_recycle = Arc(
        source=m.fs.rougher_mixer.outlet,
        destination=m.fs.solex_rougher_load.organic_inlet,
    )
    m.fs.sx_rougher_load_aq_outlet = Arc(
        source=m.fs.solex_rougher_load.aqueous_outlet,
        destination=m.fs.load_sep.inlet,
    )
    m.fs.sx_rougher_load_aq_recycle = Arc(
        source=m.fs.load_sep.recycle, destination=m.fs.leach_mixer.load_recycle
    )
    m.fs.sx_rougher_load_org_outlet = Arc(
        source=m.fs.solex_rougher_load.organic_outlet,
        destination=m.fs.solex_rougher_scrub.organic_inlet,
    )
    m.fs.sx_rougher_scrub_acid_feed = Arc(
        source=m.fs.acid_feed1.outlet,
        destination=m.fs.solex_rougher_scrub.aqueous_inlet,
    )
    m.fs.sx_rougher_scrub_aq_outlet = Arc(
        source=m.fs.solex_rougher_scrub.aqueous_outlet,
        destination=m.fs.scrub_sep.inlet,
    )
    m.fs.sx_rougher_scrub_aq_translator = Arc(
        source=m.fs.scrub_sep.recycle,
        destination=m.fs.scrubber_HCl_leach_translator.inlet,
    )
    m.fs.translator_scrub_recycle = Arc(
        source=m.fs.scrubber_HCl_leach_translator.outlet,
        destination=m.fs.leach_mixer.scrub_recycle,
    )
    m.fs.sx_rougher_scrub_org_outlet = Arc(
        source=m.fs.solex_rougher_scrub.organic_outlet,
        destination=m.fs.solex_rougher_strip.organic_inlet,
    )
    m.fs.sx_rougher_strip_acid_feed = Arc(
        source=m.fs.acid_feed2.outlet,
        destination=m.fs.solex_rougher_strip.aqueous_inlet,
    )
    m.fs.sx_rougher_strip_org_outlet = Arc(
        source=m.fs.solex_rougher_strip.organic_outlet,
        destination=m.fs.rougher_sep.inlet,
    )
    m.fs.sx_rougher_strip_org_purge = Arc(
        source=m.fs.rougher_sep.purge, destination=m.fs.rougher_organic_purge.inlet
    )
    m.fs.sx_rougher_strip_org_recycle = Arc(
        source=m.fs.rougher_sep.recycle, destination=m.fs.rougher_mixer.recycle
    )
    m.fs.sx_rougher_strip_aq_outlet = Arc(
        source=m.fs.solex_rougher_strip.aqueous_outlet,
        destination=m.fs.precip_sx_mixer.rougher,
    )
    m.fs.sx_cleaner_load_aq_feed = Arc(
        source=m.fs.precip_sx_mixer.outlet,
        destination=m.fs.solex_cleaner_load.aqueous_inlet,
    )
    m.fs.sx_cleaner_org_feed = Arc(
        source=m.fs.cleaner_org_make_up.outlet, destination=m.fs.cleaner_mixer.make_up
    )
    m.fs.sx_cleaner_mixed_org_recycle = Arc(
        source=m.fs.cleaner_mixer.outlet,
        destination=m.fs.solex_cleaner_load.organic_inlet,
    )
    m.fs.sx_cleaner_load_aq_outlet_translator = Arc(
        source=m.fs.solex_cleaner_load.aqueous_outlet,
        destination=m.fs.cleaner_HCl_leach_translator.inlet,
    )
    m.fs.sx_cleaner_load_translator_leach_sx_mixer = Arc(
        source=m.fs.cleaner_HCl_leach_translator.outlet,
        destination=m.fs.leach_sx_mixer.cleaner,
    )
    m.fs.sx_cleaner_strip_acid_feed = Arc(
        source=m.fs.acid_feed3.outlet,
        destination=m.fs.solex_cleaner_strip.aqueous_inlet,
    )
    m.fs.sx_cleaner_load_org_outlet = Arc(
        source=m.fs.solex_cleaner_load.organic_outlet,
        destination=m.fs.solex_cleaner_strip.organic_inlet,
    )
    m.fs.sx_cleaner_strip_org_outlet = Arc(
        source=m.fs.solex_cleaner_strip.organic_outlet,
        destination=m.fs.cleaner_sep.inlet,
    )
    m.fs.sx_cleaner_strip_org_purge = Arc(
        source=m.fs.cleaner_sep.purge, destination=m.fs.cleaner_organic_purge.inlet
    )
    m.fs.sx_cleaner_strip_org_recycle = Arc(
        source=m.fs.cleaner_sep.recycle, destination=m.fs.cleaner_mixer.recycle
    )
    m.fs.sx_cleaner_strip_aq_precip = Arc(
        source=m.fs.solex_cleaner_strip.aqueous_outlet,
        destination=m.fs.precipitator.aqueous_inlet,
    )
    m.fs.precip_solid_outlet = Arc(
        source=m.fs.precipitator.precipitate_outlet,
        destination=m.fs.sl_sep2.solid_inlet,
    )
    m.fs.precip_aq_sl_sep2 = Arc(
        source=m.fs.precipitator.aqueous_outlet,
        destination=m.fs.sl_sep2.liquid_inlet,
    )
    m.fs.sl_sep2_solid_outlet = Arc(
        source=m.fs.sl_sep2.solid_outlet, destination=m.fs.roaster.solid_inlet
    )
    m.fs.sl_sep2_liquid_outlet = Arc(
        source=m.fs.sl_sep2.recovered_liquid_outlet, destination=m.fs.precip_sep.inlet
    )
    m.fs.sl_sep2_retained_liquid_roaster = Arc(
        source=m.fs.sl_sep2.retained_liquid_outlet,
        destination=m.fs.roaster.liquid_inlet,
    )
    m.fs.sl_sep2_aq_purge = Arc(
        source=m.fs.precip_sep.purge,
        destination=m.fs.precip_purge.inlet,
    )
    m.fs.sl_sep2_aq_recycle = Arc(
        source=m.fs.precip_sep.recycle,
        destination=m.fs.precip_sx_mixer.precip,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_scaling(m):
    """
    Set the scaling factors to improve solver performance.

    Args:
        m: pyomo model
    """

    # Changing the default scaling factor in the class dictionary
    # applies this scaling factor globally. This sort of global
    # mutation is potentially dangerous, but we'll use it here until
    # there is a better way to set global default scaling factors
    HClStrippingPropertiesScaler.DEFAULT_SCALING_FACTORS["flow_vol"] = 1
    ModularPropertiesScaler.DEFAULT_SCALING_FACTORS["flow_mol_phase"] = 1 / 0.00781

    # Also use global mutation to change the max and min scaling factors
    # allowed from objects derived from CustomScalerBase, i.e., all the
    # model scaler objects
    CustomScalerBase.CONFIG["max_variable_scaling_factor"] = 1e20
    CustomScalerBase.CONFIG["max_constraint_scaling_factor"] = 1e20
    CustomScalerBase.CONFIG["max_expression_scaling_hint"] = 1e20

    CustomScalerBase.CONFIG["min_variable_scaling_factor"] = 1e-20
    CustomScalerBase.CONFIG["min_constraint_scaling_factor"] = 1e-20
    CustomScalerBase.CONFIG["min_expression_scaling_hint"] = 1e-20

    csb = CustomScalerBase()

    for blk in m.fs.block_data_objects(descend_into=True):
        if blk.parent_block() is m.fs:
            if isinstance(blk, UnitModelBlockData):
                if hasattr(blk, "default_scaler") and blk.default_scaler is not None:
                    print(f"Scaling {blk.name}")
                    scaler = blk.default_scaler()
                    scaler.scale_model(blk)
                else:
                    print(f"No default scaler for unit model {blk.name}")
            elif "_expanded" in blk.name:
                print(f"Scaling {blk.name}")
                # Expanded arc block
                for con in blk.component_data_objects(Constraint):
                    csb.scale_constraint_by_nominal_value(
                        con, scheme=ConstraintScalingScheme.inverseMaximum
                    )


def set_operating_conditions(m):
    """
    Set the operating conditions of the flowsheet such that the degrees of freedom are zero.

    Args:
        m: pyomo model
    """
    # Constants
    # Assume a 5% volume-by-volume ratio
    dosage = 5 / 100
    dehpa_conc = 975.8e3 * dosage * units.mg / units.L
    kerosene_conc = 8.2e5 * units.mg / units.L
    Temp_room = 303 * units.K
    P_atm = 101235 * units.Pa

    m.fs.leach_liquid_feed.properties[0.0].pressure.fix(P_atm)
    m.fs.leach_liquid_feed.properties[0.0].temperature.fix(Temp_room)
    m.fs.leach_liquid_feed.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach_liquid_feed.conc_mass_comp.fix(1e-10 * units.mg / units.L)
    m.fs.leach_liquid_feed.conc_mass_comp[0, "H2O"].fix(1e6 * units.mg / units.L)
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

    # Fix all temperatures and pressure
    m.fs.scrubber_HCl_leach_translator.outlet.temperature.fix(Temp_room)
    m.fs.scrubber_HCl_leach_translator.outlet.pressure.fix(P_atm)

    m.fs.cleaner_HCl_leach_translator.outlet.temperature.fix(Temp_room)
    m.fs.cleaner_HCl_leach_translator.outlet.pressure.fix(P_atm)

    m.fs.solex_rougher_load.mscontactor.aqueous[:, :].temperature.fix(Temp_room)
    m.fs.solex_rougher_load.mscontactor.aqueous[:, :].pressure.fix(P_atm)
    m.fs.solex_rougher_load.mscontactor.organic[:, :].temperature.fix(Temp_room)
    m.fs.solex_rougher_load.mscontactor.organic[:, :].pressure.fix(P_atm)

    m.fs.solex_rougher_scrub.mscontactor.organic[:, :].temperature.fix(Temp_room)
    m.fs.solex_rougher_scrub.mscontactor.organic[:, :].pressure.fix(P_atm)

    m.fs.solex_rougher_strip.mscontactor.organic[:, :].temperature.fix(Temp_room)
    m.fs.solex_rougher_strip.mscontactor.organic[:, :].pressure.fix(P_atm)

    m.fs.solex_cleaner_load.mscontactor.organic[:, :].temperature.fix(Temp_room)
    m.fs.solex_cleaner_load.mscontactor.organic[:, :].pressure.fix(P_atm)

    m.fs.solex_cleaner_strip.mscontactor.organic[:, :].temperature.fix(Temp_room)
    m.fs.solex_cleaner_strip.mscontactor.organic[:, :].pressure.fix(P_atm)

    m.fs.cleaner_sep.recycle_state[0.0].temperature.fix(Temp_room)
    m.fs.cleaner_sep.recycle_state[0.0].pressure.fix(P_atm)

    m.fs.cleaner_sep.purge_state[0.0].temperature.fix(Temp_room)
    m.fs.cleaner_sep.purge_state[0.0].pressure.fix(P_atm)

    m.fs.leach.mscontactor.liquid[0.0, 2].temperature.fix(Temp_room)
    m.fs.leach.mscontactor.liquid[0.0, 2].pressure.fix(P_atm)
    m.fs.leach_mixer.mixed_state[0.0].pressure.fix(P_atm)
    m.fs.leach_mixer.mixed_state[0.0].temperature.fix(Temp_room)
    m.fs.load_sep.recycle_state[0.0].pressure.fix(P_atm)
    m.fs.load_sep.recycle_state[0.0].temperature.fix(Temp_room)
    m.fs.leach_sx_mixer.mixed_state[0.0].pressure.fix(P_atm)
    m.fs.leach_sx_mixer.mixed_state[0.0].temperature.fix(Temp_room)

    m.fs.load_sep.split_fraction[:, "recycle"].fix(0.9)
    m.fs.scrub_sep.split_fraction[:, "recycle"].fix(0.9)

    m.fs.rougher_org_make_up.flow_vol.fix(6.201)

    m.fs.rougher_org_make_up.properties[0.0].pressure.fix(P_atm)
    m.fs.rougher_org_make_up.properties[0.0].temperature.fix(Temp_room)
    m.fs.rougher_mixer.mixed_state[0.0].pressure.fix(P_atm)
    m.fs.rougher_mixer.mixed_state[0.0].temperature.fix(Temp_room)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Al_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Ca_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Fe_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Sc_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Y_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "La_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Ce_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Pr_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Nd_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Sm_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Gd_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Dy_o"].fix(eps)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "DEHPA"].fix(dehpa_conc)
    m.fs.rougher_org_make_up.conc_mass_comp[0, "Kerosene"].fix(kerosene_conc)

    m.fs.acid_feed1.flow_vol.fix(90)  # Increased by 1000x from REESim
    m.fs.acid_feed1.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed1.conc_mass_comp[0, "H"].fix(10.36)
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

    m.fs.acid_feed2.flow_vol.fix(9)
    m.fs.acid_feed2.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed2.conc_mass_comp[0, "H"].fix(
        10.36 * 4
    )  # Arbitrarily choose 4x the dilute solution
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
    m.fs.rougher_sep.purge_state[0.0].pressure.fix(P_atm)
    m.fs.rougher_sep.purge_state[0.0].temperature.fix(Temp_room)
    m.fs.rougher_sep.recycle_state[0.0].pressure.fix(P_atm)
    m.fs.rougher_sep.recycle_state[0.0].temperature.fix(Temp_room)

    m.fs.acid_feed3.flow_vol.fix(9)
    m.fs.acid_feed3.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed3.conc_mass_comp[0, "H"].fix(
        10.36 * 4
    )  # Arbitrarily choose 4x the dilute solution
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

    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Al_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ca_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Fe_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sc_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Y_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "La_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ce_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Pr_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Nd_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sm_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Gd_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Dy_o"].fix(eps)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "DEHPA"].fix(dehpa_conc)
    m.fs.cleaner_org_make_up.conc_mass_comp[0, "Kerosene"].fix(kerosene_conc)
    m.fs.cleaner_org_make_up.properties[0.0].pressure.fix(P_atm)
    m.fs.cleaner_org_make_up.properties[0.0].temperature.fix(Temp_room)
    m.fs.cleaner_mixer.mixed_state[0.0].pressure.fix(P_atm)
    m.fs.cleaner_mixer.mixed_state[0.0].temperature.fix(Temp_room)

    m.fs.cleaner_sep.split_fraction[:, "recycle"].fix(0.9)

    m.fs.sl_sep1.liquid_recovery.fix(0.7)
    m.fs.sl_sep1.split.recovered_state[0.0].pressure.fix(P_atm)
    m.fs.sl_sep1.split.recovered_state[0.0].temperature.fix(Temp_room)
    m.fs.sl_sep1.split.retained_state[0.0].pressure.fix(P_atm)
    m.fs.sl_sep1.split.retained_state[0.0].temperature.fix(Temp_room)

    m.fs.sl_sep2.liquid_recovery.fix(0.9)

    m.fs.precip_sep.split_fraction[:, "recycle"].fix(0.9)

    # Fix preciptator outlet temperature
    m.fs.precipitator.precipitate_state_block[0].temperature.fix(348.15 * units.K)

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
    m.fs.roaster.frac_comp_recovery.fix(0.95)

    # Touch properties that are used in the UI
    m.fs.leach.mscontactor.solid_inlet_state[0].flow_mass
    m.fs.leach.mscontactor.solid_inlet_state[0].mass_frac_comp

    m.fs.leach.mscontactor.liquid_inlet_state[0].flow_vol
    m.fs.leach.mscontactor.liquid_inlet_state[0].conc_mol_comp

    m.fs.solex_cleaner_load.mscontactor.organic_inlet_state[0].conc_mass_comp
    m.fs.solex_cleaner_strip.mscontactor.organic_inlet_state[0].conc_mass_comp

    m.fs.solex_rougher_load.mscontactor.organic_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher_strip.mscontactor.organic_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher_scrub.mscontactor.organic_inlet_state[0].conc_mass_comp

    m.fs.solex_cleaner_load.mscontactor.aqueous_inlet_state[0].conc_mass_comp
    m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet_state[0].conc_mass_comp

    m.fs.solex_rougher_load.mscontactor.aqueous_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher_strip.mscontactor.aqueous_inlet_state[0].conc_mass_comp
    m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet_state[0].conc_mass_comp

    m.fs.precipitator.cv_aqueous.properties_out[0].flow_vol
    m.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp

    # TODO is temperature actually used in the AI?
    # This state block got switched to the HCl properties
    # which does not have temperature as a state variable
    # m.fs.precipitator.precipitate_state_block[0].temperature
    m.fs.precipitator.precipitate_state_block[0].flow_mol_comp


def initialize_system(m):
    """
    Provide initialized values for all streams in the system.

    Args:
        m: pyomo model
    """
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [
        m.fs.leaching_feed_mixture,
        m.fs.sx_rougher_load_aq_feed,
        m.fs.sx_rougher_mixed_org_recycle,
        m.fs.sx_cleaner_load_aq_feed,
        m.fs.sx_cleaner_mixed_org_recycle,
    ]

    G = seq.create_graph(m)
    order = seq.calculation_order(G)
    if _log.isEnabledFor(logging.INFO):
        _init_ord = ", ".join([o[0].name for o in order])
        _log.info("Initialization Order: {_init_ord}")

    tear_guesses1 = {
        "flow_vol": {0: 866.06},
        "conc_mass_comp": {
            (0, "Al"): 207.46,
            (0, "Ca"): 40.23,
            (0, "Ce"): 2.11,
            (0, "Cl"): 158.36,
            (0, "Dy"): 1.13e-2,
            (0, "Fe"): 292.56,
            (0, "Gd"): 0.24,
            (0, "H"): 13.66,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 1940.93,
            (0, "La"): 0.76,
            (0, "Nd"): 1.06,
            (0, "Pr"): 0.26,
            (0, "SO4"): 1438.92,
            (0, "Sc"): 2.07e-3,
            (0, "Sm"): 0.10,
            (0, "Y"): 2.02e-2,
        },
    }
    tear_guesses2 = {
        "flow_vol": {0: 62.01},
        "conc_mass_comp": {
            (0, "Al_o"): 0.048,
            (0, "Ca_o"): 1.98e-2,
            (0, "Ce_o"): 5.71e-3,
            (0, "Dy_o"): 1.077,
            (0, "Fe_o"): 1.954,
            (0, "Gd_o"): 0.14,
            (0, "La_o"): 4.03e-3,
            (0, "Nd_o"): 3.37e-3,
            (0, "Pr_o"): 1.04e-3,
            (0, "Sc_o"): 1.74,
            (0, "Sm_o"): 4.91e-3,
            (0, "Y_o"): 4.17,
            (0, "DEHPA"): 9.8e5 * 0.05,
            (0, "Kerosene"): 8.2e5,
        },
    }
    tear_guesses3 = {
        "flow_vol": {0: 623.07},
        "conc_mass_comp": {
            (0, "Al"): 320.46,
            (0, "Ca"): 62.14,
            (0, "Ce"): 3.26,
            (0, "Cl"): 192.63,
            (0, "Dy"): 4.6e-2,
            (0, "Fe"): 452.28,
            (0, "Gd"): 0.40,
            (0, "H"): 2.92,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 732.71,
            (0, "La"): 1.18,
            (0, "Nd"): 1.63,
            (0, "Pr"): 0.41,
            (0, "SO4"): 2543.95,
            (0, "Sc"): 2.25e-2,
            (0, "Sm"): 0.16,
            (0, "Y"): 0.11,
        },
    }
    tear_guesses4 = {
        "flow_vol": {0: 62},
        "conc_mass_comp": {
            (0, "Al_o"): 3.64e-3,
            (0, "Ca_o"): 2.13e-3,
            (0, "Ce_o"): 5.93e-4,
            (0, "Dy_o"): 0.33,
            (0, "Fe_o"): 0.75,
            (0, "Gd_o"): 4.00e-2,
            (0, "La_o"): 4.08e-4,
            (0, "Nd_o"): 3.76e-4,
            (0, "Pr_o"): 1.47e-4,
            (0, "Sc_o"): 3.97e-3,
            (0, "Sm_o"): 7.87e-4,
            (0, "Y_o"): 1.03,
            (0, "DEHPA"): 9.8e5 * 0.05,
            (0, "Kerosene"): 8.2e5,
        },
    }
    tear_guesses5 = {
        "flow_vol": {0: 16.70},
        "conc_mass_comp": {
            (0, "Al"): 2.42,
            (0, "Ca"): 0.68,
            (0, "Ce"): 0.16,
            (0, "Cl"): 1438.56,
            (0, "Dy"): 0.64,
            (0, "Fe"): 22.67,
            (0, "Gd"): 1.01,
            (0, "H"): 39.81,
            (0, "H2O"): 1000000,
            # (0, "HSO4"): 2.88e-6,
            (0, "La"): 0.13,
            (0, "Nd"): 8.52e-2,
            (0, "Pr"): 2.10e-2,
            # (0, "SO4"): 2.54e-6,
            (0, "Sc"): 1.65e-3,
            (0, "Sm"): 7.88e-2,
            (0, "Y"): 1.17,
        },
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.leach.liquid_inlet, tear_guesses1)
    seq.set_guesses_for(m.fs.solex_rougher_load.organic_inlet, tear_guesses2)
    seq.set_guesses_for(m.fs.solex_rougher_load.aqueous_inlet, tear_guesses3)
    seq.set_guesses_for(m.fs.solex_cleaner_load.organic_inlet, tear_guesses4)
    seq.set_guesses_for(m.fs.solex_cleaner_load.aqueous_inlet, tear_guesses5)

    initializer_feed = FeedInitializer()
    feed_units = [
        m.fs.leach_liquid_feed,
        m.fs.leach_solid_feed,
        m.fs.rougher_org_make_up,
        m.fs.acid_feed1,
        m.fs.acid_feed2,
        m.fs.acid_feed3,
        m.fs.cleaner_org_make_up,
    ]

    initializer_product = ProductInitializer()
    product_units = [
        m.fs.leach_filter_cake,
        m.fs.leach_filter_cake_liquid,
        m.fs.cleaner_organic_purge,
        m.fs.rougher_organic_purge,
        m.fs.precip_purge,
    ]

    initializer_sep = SeparatorInitializer()
    sep_units = [
        m.fs.scrub_sep,
        m.fs.precip_sep,
    ]

    initializer_mix = MixerInitializer()
    mix_units = [
        m.fs.precip_sx_mixer,
    ]

    initializer_leach = LeachingTrainInitializer()
    leach_units = [
        m.fs.leach,
    ]

    initializer_sx = SolventExtractionInitializer()
    sx_units = [
        m.fs.solex_rougher_load,
        m.fs.solex_rougher_scrub,
        m.fs.solex_rougher_strip,
        m.fs.solex_cleaner_load,
        m.fs.solex_cleaner_strip,
    ]

    initializer_bt = BlockTriangularizationInitializer()

    def function(unit):
        if unit in feed_units:
            _log.info(f"Initializing {unit}")
            initializer_feed.initialize(unit)
        elif unit in product_units:
            _log.info(f"Initializing {unit}")
            initializer_product.initialize(unit)
        elif unit in sep_units:
            _log.info(f"Initializing {unit}")
            initializer_sep.initialize(unit)

        elif unit in mix_units:
            _log.info(f"Initializing {unit}")
            initializer_mix.initialize(unit)
        elif unit in leach_units:
            _log.info(f"Initializing {unit}")
            initializer_leach.initialize(unit)
        elif unit in sx_units:
            _log.info(f"Initializing {unit}")
            initializer_sx.initialize(unit)
        else:
            _log.info(f"Initializing {unit}")
            initializer_bt.initialize(unit)

    seq.run(m, function)


def solve_system(m, solver=None, tee=False):
    """
    Solve the model.

    Args:
        m: pyomo model
        solver: optimization solver
        tee: boolean indicator to stream IPOPT solution
    """
    if hasattr(solver, "solve"):
        solver = solver
    else:
        # Why isn't it getting ipopt_v2 automatically?
        solver = get_solver("ipopt_v2")
    solver.options.constr_viol_tol = 1e-8

    results = solver.solve(m, tee=tee)

    return results


def fix_organic_recycle(m):
    """
    Fix the volumetric flow rate of the organic recycle streams and unfix the flow of make-up streams.

    Args:
        m: pyomo model
    """

    m.fs.rougher_org_make_up.outlet.flow_vol.unfix()
    m.fs.rougher_mixer.outlet.flow_vol.fix(62.01)

    m.fs.cleaner_org_make_up.outlet.flow_vol.unfix()
    m.fs.cleaner_mixer.outlet.flow_vol.fix(62.01)


def add_result_expressions(m):
    fs = m.fs
    ree_list = ["Sc", "Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]
    gangue_list = ["Al", "Fe", "Ca"]
    fs.ree_set = Set(initialize=ree_list)
    fs.gangue_set = Set(initialize=gangue_list)
    fs.metal_set = Set(initialize=ree_list + gangue_list)

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

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals into leaching. "
        "Includes metals from both the coal refuse solids "
        "and the solvent extraction recycle stream",
    )
    def leaching_metal_inlet_flow(b, t, j):
        if j == "Ca":
            j_oxide = "CaO"
        else:
            j_oxide = f"{j}2O3"

        return units.convert(
            b.leach_solid_feed.flow_mass[t]
            * b.leach_solid_feed.mass_frac_comp[t, j_oxide]
            * metal_mass_frac[j_oxide],
            to_units=units.kg / units.hr,
        ) + units.convert(
            b.leach.mscontactor.liquid_inlet_state[t].conc_mass_comp[j]
            * b.leach.mscontactor.liquid_inlet_state[t].flow_vol,
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals fed into the overall process",
    )
    def metal_feed_flow(b, t, j):
        if j == "Ca":
            j_oxide = "CaO"
        else:
            j_oxide = f"{j}2O3"

        return units.convert(
            b.leach_solid_feed.flow_mass[t]
            * b.leach_solid_feed.mass_frac_comp[t, j_oxide]
            * metal_mass_frac[j_oxide],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving sl_sep1 to go " "to solvent extraction.",
    )
    def leaching_metal_recovery_flow(b, t, j):
        return units.convert(
            b.sl_sep1.recovered_liquid_outlet.conc_mass_comp[t, j]
            * b.sl_sep1.recovered_liquid_outlet.flow_vol[t],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the flowsheet product stream.",
    )
    def metal_product_flow(b, t, j):
        if j == "Ca":
            j_oxide = "CaO"
        else:
            j_oxide = f"{j}2O3"

        return units.convert(
            b.roaster.flow_mol_comp_product[t, j]
            * molar_mass[j_oxide]
            * metal_mass_frac[j_oxide],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the leach filter cake solids",
    )
    def metal_leach_filter_cake_solids_flow(b, t, j):
        if j == "Ca":
            j_oxide = "CaO"
        else:
            j_oxide = f"{j}2O3"
        return units.convert(
            +b.leach_filter_cake.flow_mass[t]
            * b.leach_filter_cake.mass_frac_comp[t, j_oxide]
            * metal_mass_frac[j_oxide],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "liquid entrained in the filter cake",
    )
    def metal_leach_filter_cake_liquid_flow(b, t, j):
        return units.convert(
            b.leach_filter_cake_liquid.flow_vol[t]
            * b.leach_filter_cake_liquid.conc_mass_comp[t, j],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the rougher load aqueous purge stream",
    )
    def metal_rougher_load_aqueous_purge_flow(b, t, j):
        return units.convert(
            b.load_sep.purge.flow_vol[t] * b.load_sep.purge.conc_mass_comp[t, j],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the rougher scrub aqueous purge stream",
    )
    def metal_rougher_scrub_aqueous_purge_flow(b, t, j):
        return units.convert(
            b.scrub_sep.purge.flow_vol[t] * b.scrub_sep.purge.conc_mass_comp[t, j],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the rougher organic purge stream",
    )
    def metal_rougher_organic_purge_flow(b, t, j):
        return units.convert(
            b.rougher_sep.purge.flow_vol[t]
            * b.rougher_sep.purge.conc_mass_comp[t, f"{j}_o"],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the cleaner aqueous purge stream",
    )
    def metal_cleaner_aqueous_purge_flow(b, t, j):
        return units.convert(
            b.precip_sep.purge.flow_vol[t] * b.precip_sep.purge.conc_mass_comp[t, j],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Mass flow rate of metals leaving process in "
        "the cleaner organic purge stream",
    )
    def metal_cleaner_organic_purge_flow(b, t, j):
        return units.convert(
            b.cleaner_sep.purge.flow_vol[t]
            * b.cleaner_sep.purge.conc_mass_comp[t, f"{j}_o"],
            to_units=units.kg / units.hr,
        )

    @fs.Expression(
        fs.time, fs.metal_set, doc="Single-pass recovery for each metal in leaching"
    )
    def leaching_metal_recovery_percentage(b, t, j):
        return (
            100
            * b.leaching_metal_recovery_flow[t, j]
            / b.leaching_metal_inlet_flow[t, j]
        )

    @fs.Expression(
        fs.time, fs.metal_set, doc="Overall recovery of each metal for the flowsheet"
    )
    def overall_metal_recovery_percentage(b, t, j):
        return 100 * b.metal_product_flow[t, j] / b.metal_feed_flow[t, j]

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected in the " "leach filter cake solids.",
    )
    def leach_filter_cake_solids_percent_rejected(b, t, j):
        return (
            100 * b.metal_leach_filter_cake_solids_flow[t, j] / b.metal_feed_flow[t, j]
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected in the "
        "liquid entrained in leach filter cake.",
    )
    def leach_filter_cake_liquid_percent_rejected(b, t, j):
        return (
            100 * b.metal_leach_filter_cake_liquid_flow[t, j] / b.metal_feed_flow[t, j]
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected through the rougher "
        "load aqueous purge streams",
    )
    def metal_rougher_load_aqueous_percent_rejected(b, t, j):
        return (
            100
            * b.metal_rougher_load_aqueous_purge_flow[t, j]
            / b.metal_feed_flow[t, j]
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected through the rougher "
        "scrub aqueous purge streams",
    )
    def metal_rougher_scrub_aqueous_percent_rejected(b, t, j):
        return (
            100
            * b.metal_rougher_scrub_aqueous_purge_flow[t, j]
            / b.metal_feed_flow[t, j]
        )

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected through the rougher " "organic purge stream",
    )
    def metal_rougher_organic_percent_rejected(b, t, j):
        return 100 * b.metal_rougher_organic_purge_flow[t, j] / b.metal_feed_flow[t, j]

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected through the cleaner " "aqueous purge stream",
    )
    def metal_cleaner_aqueous_percent_rejected(b, t, j):
        return 100 * b.metal_cleaner_aqueous_purge_flow[t, j] / b.metal_feed_flow[t, j]

    @fs.Expression(
        fs.time,
        fs.metal_set,
        doc="Percent of metal rejected through the cleaner " "organic purge stream",
    )
    def metal_cleaner_organic_percent_rejected(b, t, j):
        return 100 * b.metal_cleaner_organic_purge_flow[t, j] / b.metal_feed_flow[t, j]

    @fs.Expression(fs.time, doc="Total mass flow rate of all REEs into process.")
    def ree_feed_flow(b, t):
        return sum(b.metal_feed_flow[t, i] for i in b.ree_set)

    @fs.Expression(
        fs.time, doc="Total mass flow rate of all REEs in flowsheet product stream."
    )
    def ree_product_flow(b, t):
        return sum(b.metal_product_flow[t, i] for i in b.ree_set)

    @fs.Expression(fs.time)
    def overall_ree_recovery_percentage(b, t):
        return 100 * b.ree_product_flow[t] / b.ree_feed_flow[t]

    @fs.Expression(fs.time)
    def ree_product_purity_percentage(b, t):
        return (
            100
            * b.ree_product_flow[t]
            / units.convert(
                b.roaster.flow_mass_product[0], to_units=units.kg / units.hr
            )
        )


def display_results(m):
    """
    Print key flowsheet outputs.

    Args:
        m: pyomo model
    """
    metal_name_dict = {
        "Sc": "scandium",
        "Y": "yttrium",
        "La": "lanthanum",
        "Ce": "cerium",
        "Pr": "praseodymium",
        "Nd": "neodynium",
        "Sm": "samarium",
        "Gd": "gadolinium",
        "Dy": "dysprosium",
        "Al": "aluminum",
        "Fe": "iron",
        "Ca": "calcium",
    }

    def print_element_report(j):
        name = metal_name_dict[j]
        leach_recovery_percent = value(m.fs.leaching_metal_recovery_percentage[0, j])
        overall_recovery_percent = value(m.fs.overall_metal_recovery_percentage[0, j])
        lfcs_reject = value(m.fs.leach_filter_cake_solids_percent_rejected[0, j])
        lfcl_reject = value(m.fs.leach_filter_cake_liquid_percent_rejected[0, j])
        ral_reject = value(m.fs.metal_rougher_load_aqueous_percent_rejected[0, j])
        ras_reject = value(m.fs.metal_rougher_scrub_aqueous_percent_rejected[0, j])
        ro_reject = value(m.fs.metal_rougher_organic_percent_rejected[0, j])
        ca_reject = value(m.fs.metal_cleaner_aqueous_percent_rejected[0, j])
        co_reject = value(m.fs.metal_cleaner_organic_percent_rejected[0, j])
        total = (
            overall_recovery_percent
            + lfcs_reject
            + lfcl_reject
            + ral_reject
            + ras_reject
            + ro_reject
            + ca_reject
            + co_reject
        )
        print(f"\nLeaching {name} recovery is {leach_recovery_percent} %")
        print(f"Overall {name} recovery is {overall_recovery_percent} %")
        print(
            f"{name.capitalize()} rejected in leach filter cake solids is {lfcs_reject} %"
        )
        print(
            f"{name.capitalize()} rejected in leach filter cake entrained liquid is {lfcl_reject} %"
        )
        print(
            f"{name.capitalize()} rejected in rougher load aqueous purge is {ral_reject} %"
        )
        print(
            f"{name.capitalize()} rejected in rougher scrub aqueous purge is {ras_reject} %"
        )
        print(f"{name.capitalize()} rejected in rougher organic purge is {ro_reject} %")
        print(f"{name.capitalize()} rejected in cleaner aqueous purge is {ca_reject} %")
        print(f"{name.capitalize()} rejected in cleaner organic purge is {co_reject} %")
        print(f"{name.capitalize()} accounted for is {total} %")
        print("")

    m.fs.roaster.report()
    print(f"REE product mass flow is {value(m.fs.ree_product_flow[0])} kg/hr")
    print(f"REE feed mass flow is {value(m.fs.ree_feed_flow[0])} kg/hr")
    print(f"Total REE recovery is {value(m.fs.overall_ree_recovery_percentage[0])} %")
    print(f"Product purity is {value(m.fs.ree_product_purity_percentage[0])} % REE")
    print("")
    print("")
    print("REE element recovery:")
    for j in m.fs.ree_set:
        print_element_report(j)
    print("")
    print("Gangue recovery:")
    for j in m.fs.gangue_set:
        print_element_report(j)
    print("")


def calculate_results(fs):
    """
    Calculate key flowsheet output for use in the flowsheet UI.

    Args:
        fs: Flowsheet
    """

    data = {}  # put all results here

    # Total mass basis yield calculation
    data["REE-product"] = value(fs.ree_product_flow[0])
    data["REE-feed"] = value(fs.ree_feed_flow[0])
    data["REE-recovery"] = value(fs.overall_ree_recovery_percentage[0])
    data["product-purity"] = value(fs.ree_product_purity_percentage[0])

    # Individual elemental recoveries
    for metal in fs.metal_set:
        data[f"total-{metal}-recovery"] = value(
            fs.overall_metal_recovery_percentage[0, metal]
        )
    for metal in fs.metal_set:
        data[f"{metal}-recovery"] = value(
            fs.leaching_metal_recovery_percentage[0, metal]
        )

    return data


def add_costing(m):
    """
    Set the costing parameters for each unit model.

    Args:
        m: pyomo model
    """

    m.fs.scaling_constraints = ComponentMap()

    m.fs.costing = QGESSCosting()
    CE_index_year = "UKy_2019"

    # define reference values for empirical scaling to estimate balance of
    # plant unit operation process parameters
    # scaled_parameter = reference_parameter * (scaled_basis_flow/reference_basis_flow)

    # reference values from UKy study - Table 4-7 p. 351
    REE_costing_params = load_REE_costing_dictionary()
    reference_basis_flow = {
        "leach_sol_flow_mass": 495 * units.ton / units.hr,  # p. 273, 351
        "rougher_solex_aqueous_flow_vol": 23131 * units.L / units.min,
        "cleaner_solex_aqueous_flow_vol": 925 * units.L / units.min,
        "precipitator_solex_aqueous_flow_vol": 231 * units.L / units.min,
    }

    # Leaching costs
    # 4.2 is UKy Leaching - Polyethylene Tanks
    L_pe_tanks_accounts = ["4.2"]
    m.fs.leach.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pe_tanks_accounts,
            "scaled_param": m.fs.leach.volume[0, 1],
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.3 is UKy Leaching - Tank Mixer
    L_tank_mixer_accounts = ["4.3"]
    m.fs.leach_mixer.power = Var(initialize=4.74, units=units.hp, bounds=(0, None))

    @m.fs.leach_mixer.Constraint(L_tank_mixer_accounts)
    def power_scaling_constraint(c, k):
        return m.fs.leach_mixer.power == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.hp
            * (
                m.fs.leach_solid_feed.flow_mass[0]
                / reference_basis_flow["leach_sol_flow_mass"]
            ),
            to_units=units.hp,
        )

    m.fs.scaling_constraints[m.fs.leach_mixer.power] = (
        m.fs.leach_mixer.power_scaling_constraint
    )

    m.fs.leach_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_tank_mixer_accounts,
            "scaled_param": m.fs.leach_mixer.power,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.4 is UKy Leaching - Process Pump
    L_pump_accounts = ["4.4"]
    m.fs.leach_pump = UnitModelBlock()
    m.fs.leach_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pump_accounts,
            "scaled_param": m.fs.leach_liquid_feed.flow_vol[0],
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.5 is UKy Leaching - Thickener
    L_thickener_accounts = ["4.5"]
    m.fs.leach_sx_mixer.area = Var(
        initialize=225.90, units=units.ft**2, bounds=(0, None)
    )

    @m.fs.leach_sx_mixer.Constraint(L_thickener_accounts)
    def area_scaling_constraint(c, k):
        return m.fs.leach_sx_mixer.area == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.ft**2
            * (
                m.fs.leach_solid_feed.flow_mass[0]
                / reference_basis_flow["leach_sol_flow_mass"]
            ),
            to_units=units.ft**2,
        )

    m.fs.scaling_constraints[m.fs.leach_sx_mixer.area] = (
        m.fs.leach_sx_mixer.area_scaling_constraint
    )

    m.fs.leach_sx_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_thickener_accounts,
            "scaled_param": m.fs.leach_sx_mixer.area,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.6 is UKy Leaching - Solid Waste Filter Press
    L_filter_press_accounts = ["4.6"]
    m.fs.sl_sep1.volume = Var(initialize=36.00, units=units.ft**3, bounds=(0, None))

    @m.fs.sl_sep1.Constraint(L_filter_press_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.sl_sep1.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.ft**3
            * (
                m.fs.leach_solid_feed.flow_mass[0]
                / reference_basis_flow["leach_sol_flow_mass"]
            ),
            to_units=units.ft**3,
        )

    m.fs.scaling_constraints[m.fs.sl_sep1.volume] = (
        m.fs.sl_sep1.volume_scaling_constraint
    )

    m.fs.sl_sep1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_filter_press_accounts,
            "scaled_param": m.fs.sl_sep1.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.8 is UKy Leaching - Solution Heater
    L_solution_heater_accounts = ["4.8"]
    m.fs.leach_solution_heater = UnitModelBlock()
    m.fs.leach_solution_heater.duty = Var(
        initialize=0.24, units=units.MBTU / units.hr, bounds=(0, None)
    )

    @m.fs.leach_solution_heater.Constraint(L_solution_heater_accounts)
    def duty_scaling_constraint(c, k):
        return m.fs.leach_solution_heater.duty == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.MBTU
            / units.hr
            * (
                m.fs.leach_solid_feed.flow_mass[0]
                / reference_basis_flow["leach_sol_flow_mass"]
            ),
            to_units=units.MBTU / units.hr,
        )

    m.fs.scaling_constraints[m.fs.leach_solution_heater.duty] = (
        m.fs.leach_solution_heater.duty_scaling_constraint
    )

    m.fs.leach_solution_heater.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_solution_heater_accounts,
            "scaled_param": m.fs.leach_solution_heater.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # Solvent extraction costs
    # 5.1 is UKy Rougher Solvent Extraction - Polyethylene Tanks
    RSX_pe_tanks_accounts = ["5.1"]
    m.fs.rougher_solex_tank = UnitModelBlock()
    m.fs.rougher_solex_tank.volume = Var(
        initialize=35.136, units=units.gal, bounds=(0, None)
    )

    @m.fs.rougher_solex_tank.Constraint(RSX_pe_tanks_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.rougher_solex_tank.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.gal
            * (
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["rougher_solex_aqueous_flow_vol"]
            ),
            to_units=units.gal,
        )

    m.fs.scaling_constraints[m.fs.rougher_solex_tank.volume] = (
        m.fs.rougher_solex_tank.volume_scaling_constraint
    )

    m.fs.rougher_solex_tank.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pe_tanks_accounts,
            "scaled_param": m.fs.rougher_solex_tank.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.2 is UKy Rougher Solvent Extraction - Tank Mixer
    RSX_tank_mixer_accounts = ["5.2"]
    m.fs.rougher_mixer.power = Var(initialize=2.0, units=units.hp, bounds=(0, None))

    @m.fs.rougher_mixer.Constraint(RSX_tank_mixer_accounts)
    def power_scaling_constraint(c, k):
        return m.fs.rougher_mixer.power == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.hp
            * (
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["rougher_solex_aqueous_flow_vol"]
            ),
            to_units=units.hp,
        )

    m.fs.scaling_constraints[m.fs.rougher_mixer.power] = (
        m.fs.rougher_mixer.power_scaling_constraint
    )

    m.fs.rougher_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_tank_mixer_accounts,
            "scaled_param": m.fs.rougher_mixer.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.3 is UKy Rougher Solvent Extraction - Process Pump
    RSX_pump_accounts = ["5.3"]
    m.fs.rougher_pump = UnitModelBlock()
    m.fs.rougher_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pump_accounts,
            "scaled_param": m.fs.solex_rougher_load.mscontactor.aqueous_inlet.flow_vol[
                0
            ],
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.4 is UKy Rougher Solvent Extraction - Mixer Settler
    RSX_mixer_settler_accounts = ["5.4"]
    m.fs.rougher_solex_settler = UnitModelBlock()
    m.fs.rougher_solex_settler.volume = Var(
        initialize=61.107, units=units.gal, bounds=(0, None)
    )

    @m.fs.rougher_solex_settler.Constraint(RSX_mixer_settler_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.rougher_solex_settler.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.gal
            * (
                m.fs.solex_rougher_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["rougher_solex_aqueous_flow_vol"]
            ),
            to_units=units.gal,
        )

    m.fs.scaling_constraints[m.fs.rougher_solex_settler.volume] = (
        m.fs.rougher_solex_settler.volume_scaling_constraint
    )

    m.fs.rougher_solex_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_mixer_settler_accounts,
            "scaled_param": m.fs.rougher_solex_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.1 is UKy Cleaner Solvent Extraction - Polyethylene Tanks
    CSX_pe_tanks_accounts = ["6.1"]
    m.fs.cleaner_solex_tank = UnitModelBlock()
    m.fs.cleaner_solex_tank.volume = Var(
        initialize=14.05, units=units.gal, bounds=(0, None)
    )

    @m.fs.cleaner_solex_tank.Constraint(CSX_pe_tanks_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.cleaner_solex_tank.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.gal
            * (
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["cleaner_solex_aqueous_flow_vol"]
            ),
            to_units=units.gal,
        )

    m.fs.scaling_constraints[m.fs.cleaner_solex_tank.volume] = (
        m.fs.cleaner_solex_tank.volume_scaling_constraint
    )

    m.fs.cleaner_solex_tank.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pe_tanks_accounts,
            "scaled_param": m.fs.cleaner_solex_tank.volume,
            "source": 1,
            "n_equip": 5,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.2 is UKy Cleaner Solvent Extraction - Tank Mixer
    CSX_tank_mixer_accounts = ["6.2"]
    m.fs.cleaner_mixer.power = Var(initialize=0.08, units=units.hp, bounds=(0, None))

    @m.fs.cleaner_mixer.Constraint(CSX_tank_mixer_accounts)
    def power_scaling_constraint(c, k):
        return m.fs.cleaner_mixer.power == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.hp
            * (
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["cleaner_solex_aqueous_flow_vol"]
            ),
            to_units=units.hp,
        )

    m.fs.scaling_constraints[m.fs.cleaner_mixer.power] = (
        m.fs.cleaner_mixer.power_scaling_constraint
    )

    m.fs.cleaner_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_tank_mixer_accounts,
            "scaled_param": m.fs.cleaner_mixer.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.3 is UKy Cleaner Solvent Extraction - Process Pump
    CSX_pump_accounts = ["6.3"]
    m.fs.cleaner_pump = UnitModelBlock()
    m.fs.cleaner_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pump_accounts,
            "scaled_param": m.fs.solex_cleaner_load.mscontactor.aqueous_inlet.flow_vol[
                0
            ],
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.4 is UKy Cleaner Solvent Extraction - Mixer Settler
    CSX_mixer_settler_accounts = ["6.4"]
    m.fs.cleaner_solex_settler = UnitModelBlock()
    m.fs.cleaner_solex_settler.volume = Var(
        initialize=24.44, units=units.gal, bounds=(0, None)
    )

    @m.fs.cleaner_solex_settler.Constraint(CSX_mixer_settler_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.cleaner_solex_settler.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.gal
            * (
                m.fs.solex_cleaner_load.mscontactor.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["cleaner_solex_aqueous_flow_vol"]
            ),
            to_units=units.gal,
        )

    m.fs.scaling_constraints[m.fs.cleaner_solex_settler.volume] = (
        m.fs.cleaner_solex_settler.volume_scaling_constraint
    )

    m.fs.cleaner_solex_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_mixer_settler_accounts,
            "scaled_param": m.fs.cleaner_solex_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # Precipitation costs
    # 10.1 is UKy Oxalate Precipitation - Polyethylene Tanks
    reep_pe_tanks_accounts = ["10.1"]
    m.fs.precipitator.volume = Var(initialize=15.04, units=units.gal, bounds=(0, None))

    @m.fs.precipitator.Constraint(reep_pe_tanks_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.precipitator.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.gal
            * (
                m.fs.precipitator.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["precipitator_solex_aqueous_flow_vol"]
            ),
            to_units=units.gal,
        )

    m.fs.scaling_constraints[m.fs.precipitator.volume] = (
        m.fs.precipitator.volume_scaling_constraint
    )

    m.fs.precipitator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pe_tanks_accounts,
            "scaled_param": m.fs.precipitator.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 10.2 is UKy Oxalate Precipitation - Tank Mixer
    reep_tank_mixer_accounts = ["10.2"]
    m.fs.precipitator_mixer = UnitModelBlock()
    m.fs.precipitator_mixer.power = Var(
        initialize=0.61, units=units.hp, bounds=(0, None)
    )

    @m.fs.precipitator_mixer.Constraint(reep_tank_mixer_accounts)
    def power_scaling_constraint(c, k):
        return m.fs.precipitator_mixer.power == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.hp
            * (
                m.fs.precipitator.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["precipitator_solex_aqueous_flow_vol"]
            ),
            to_units=units.hp,
        )

    m.fs.scaling_constraints[m.fs.precipitator_mixer.power] = (
        m.fs.precipitator_mixer.power_scaling_constraint
    )

    m.fs.precipitator_mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_tank_mixer_accounts,
            "scaled_param": m.fs.precipitator_mixer.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 10.3 is UKy Oxalate Precipitation - Process Pump
    reep_pump_accounts = ["10.3"]
    m.fs.precipitator_pump = UnitModelBlock()
    m.fs.precipitator_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pump_accounts,
            "scaled_param": m.fs.precipitator.aqueous_inlet.flow_vol[0],
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 10.4 is UKy Oxalate Precipitation - Filter Press
    reep_filter_press_accounts = ["10.4"]
    m.fs.sl_sep2.volume = Var(initialize=0.405, units=units.ft**3, bounds=(0, None))

    @m.fs.sl_sep2.Constraint(reep_filter_press_accounts)
    def volume_scaling_constraint(c, k):
        return m.fs.sl_sep2.volume == units.convert(
            REE_costing_params["1"][k]["RP Value"]
            * units.ft**3
            * (
                m.fs.precipitator.aqueous_inlet.flow_vol[0]
                / reference_basis_flow["precipitator_solex_aqueous_flow_vol"]
            ),
            to_units=units.ft**3,
        )

    m.fs.scaling_constraints[m.fs.sl_sep2.volume] = (
        m.fs.sl_sep2.volume_scaling_constraint
    )

    m.fs.sl_sep2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_filter_press_accounts,
            "scaled_param": m.fs.sl_sep2.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 10.5 is UKy Oxalate Precipitation - Roaster
    reep_roaster_accounts = ["10.5"]
    m.fs.roaster.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_roaster_accounts,
            "scaled_param": abs(m.fs.roaster.heat_duty[0]),
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
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
        m.fs.leach_solid_feed.flow_mass[0]
        * m.fs.leach_solid_feed.mass_frac_comp[0, molecule]
        * REE_frac
        for molecule, REE_frac in REE_mass_frac.items()
    )

    m.fs.feed_input = Var(
        initialize=0.025, units=units.ton / units.hr, bounds=(0, None)
    )
    m.fs.feed_input_constraint = Constraint(
        expr=m.fs.feed_input
        == units.convert(
            m.fs.leach_solid_feed.flow_mass[0], to_units=units.ton / units.hr
        )
    )

    m.fs.scaling_constraints[m.fs.feed_input] = m.fs.feed_input_constraint

    m.fs.feed_grade = Var(initialize=318.015, units=units.ppm, bounds=(0, None))
    m.fs.feed_grade_constraint = Constraint(
        expr=m.fs.feed_grade
        == units.convert(
            feed_REE / m.fs.leach_solid_feed.flow_mass[0], to_units=units.ppm
        )
    )

    m.fs.scaling_constraints[m.fs.feed_grade] = m.fs.feed_grade_constraint

    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    # for convenience
    m.fs.annual_operating_hours = Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=True,
        units=units.hours / units.a,
    )

    m.fs.recovery_rate_per_year = Var(
        initialize=13.306, units=units.kg / units.yr, bounds=(0, None)
    )
    m.fs.recovery_rate_per_year_constraint = Constraint(
        expr=m.fs.recovery_rate_per_year
        == units.convert(
            m.fs.roaster.flow_mass_product[0] * m.fs.annual_operating_hours,
            to_units=units.kg / units.yr,
        )
    )

    m.fs.scaling_constraints[m.fs.recovery_rate_per_year] = (
        m.fs.recovery_rate_per_year_constraint
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

    m.fs.solid_waste = Var(
        m.fs.time, initialize=0.0245, units=units.ton / units.hr, bounds=(0, None)
    )
    m.fs.solid_waste_constraint = Constraint(
        expr=m.fs.solid_waste[0]
        == units.convert(
            m.fs.leach_filter_cake.flow_mass[0], to_units=units.ton / units.hr
        )
    )

    m.fs.scaling_constraints[m.fs.solid_waste] = m.fs.solid_waste_constraint

    # TODO where is the corresponding constraint for this Var?
    m.fs.precipitate = Var(
        m.fs.time, initialize=1e-8, units=units.ton / units.hr, bounds=(0, None)
    )  # non-hazardous precipitate

    m.fs.dust_and_volatiles = Var(
        m.fs.time, initialize=9.5e-8, units=units.ton / units.hr, bounds=(0, None)
    )
    m.fs.dust_and_volatiles_constraint = Constraint(
        expr=m.fs.dust_and_volatiles[0]
        == units.convert(m.fs.roaster.flow_mass_dust[0], to_units=units.ton / units.hr)
    )
    m.fs.scaling_constraints[m.fs.dust_and_volatiles] = (
        m.fs.dust_and_volatiles_constraint
    )

    m.fs.power = Var(m.fs.time, initialize=7, units=units.hp, bounds=(0, None))
    m.fs.power_constraint = Constraint(
        expr=m.fs.power[0]
        == units.convert(
            m.fs.precipitator_mixer.power
            + m.fs.cleaner_mixer.power
            + m.fs.rougher_mixer.power
            + m.fs.leach_mixer.power,
            to_units=units.hp,
        )
    )

    m.fs.scaling_constraints[m.fs.power] = m.fs.power_constraint

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

    # TODO Why are these Params and not Expressions?
    m.fs.Ce_product = Param(
        default=units.convert(
            m.fs.roaster.flow_mol_comp_product[0, "Ce"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Dy"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Gd"]
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
            m.fs.roaster.flow_mol_comp_product[0, "La"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Nd"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Pr"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Sc"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Sm"]
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
            m.fs.roaster.flow_mol_comp_product[0, "Y"]
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
    m.fs.costing.other_plant_costs_eq = Constraint(
        expr=(
            m.fs.costing.other_plant_costs
            == units.convert(
                1218.073
                * units.USD_2016  # Rougher Solvent Extraction
                * (
                    m.fs.solex_rougher_load.mscontactor.aqueous_inlet.flow_vol[0]
                    / reference_basis_flow["rougher_solex_aqueous_flow_vol"]
                )
                ** 0.7
                + 48.723
                * units.USD_2016  # Cleaner Solvent Extraction
                * (
                    m.fs.solex_cleaner_load.mscontactor.aqueous_inlet.flow_vol[0]
                    / reference_basis_flow["cleaner_solex_aqueous_flow_vol"]
                )
                ** 0.7
                + 182.711
                * units.USD_2016  # Solvent Extraction Wash and Saponification
                * (
                    m.fs.precipitator.aqueous_inlet.flow_vol[0]
                    / reference_basis_flow["precipitator_solex_aqueous_flow_vol"]
                )
                ** 0.7,
                to_units=getattr(units, "MUSD_" + CE_index_year),
            )
        )
    )

    # fix costing vars that shouldn't change
    m.fs.precipitate.fix()

    return m


def initialize_costing(m):
    """
    Initializes costing by calling block triangularization on the entire flowsheet.

    Args:
        m: Model containing flowsheet with already-initialized unit models.
    """
    init = BlockTriangularizationInitializer()
    init.initialize(m)


def display_costing(m):
    """
    Print the key costing results.

    Args:
        m: pyomo model
    """
    QGESSCostingData.report(m.fs.costing)
    m.fs.costing.variable_operating_costs.display()
    QGESSCostingData.display_bare_erected_costs(m.fs.costing)
    QGESSCostingData.display_flowsheet_cost(m.fs.costing)


def optimize_model(m):
    m.obj = Objective(
        expr=(
            0.01
            * (
                m.fs.leach_liquid_feed.flow_vol[0]
                + m.fs.acid_feed1.flow_vol[0]
                + m.fs.acid_feed2.flow_vol[0]
                + m.fs.acid_feed3.flow_vol[0]
            )
            + 1e5
            * (
                -m.fs.ree_product_flow[0]
                + m.fs.metal_product_flow[0, "Al"]
                + m.fs.metal_product_flow[0, "Ca"]
                + m.fs.metal_product_flow[0, "Fe"]
            )
        )
    )
    # Unfix the H2SO4 feed rate and feed concentration
    m.fs.leach_liquid_feed.flow_vol.unfix()

    m.fs.leach_liquid_feed.conc_mass_comp[0, "H"].unfix()
    m.fs.leach_liquid_feed.conc_mass_comp[0, "HSO4"].unfix()
    m.fs.leach_liquid_feed.conc_mass_comp[0, "SO4"].unfix()

    @m.fs.leach_liquid_feed.Constraint(m.fs.time)
    def H2SO4_stoich_eqn(b, t):
        return (
            b.properties[t].conc_mol_comp["H"]
            == 2 * b.properties[t].conc_mol_comp["SO4"]
            + b.properties[t].conc_mol_comp["HSO4"]
        )

    for condata in m.fs.leach_liquid_feed.H2SO4_stoich_eqn.values():
        set_scaling_factor(condata, 10)

    # Because we have defined_state=True for the feed
    # block, we need to create the dissociation
    # equilibrium manually
    # TODO maybe we should convert it into a FeedFlash?
    @m.fs.leach_liquid_feed.Constraint(m.fs.time)
    def HSO4_dissociation(b, t):
        return (
            b.properties[t].params.Ka2 * b.properties[t].conc_mol_comp["HSO4"]
            == b.properties[t].conc_mol_comp["SO4"] * b.properties[t].conc_mol_comp["H"]
        )

    sf = get_scaling_factor(m.fs.leach.mscontactor.liquid[0, 1].hso4_dissociation)
    for condata in m.fs.leach_liquid_feed.HSO4_dissociation.values():
        set_scaling_factor(condata, sf)

    # According to Wikipedia, pH values of less than 0 are
    # impossible in an aqueous solution
    m.fs.leach_liquid_feed.properties[0].pH_phase["liquid"].setlb(0)

    # IPOPT wanted to use a liquid feed of 13.7 L/hr, which would be
    # far too low to let the leach solids move as a slurry.
    # The original flow rate is 224.3 L/hr, so use 100 as a rough
    # estimate for a lower bound.
    m.fs.leach_liquid_feed.properties[0].flow_vol.setlb(100)

    # The true figure of merit is probably the volumetric flow rate
    # out of the mixer, which is 68.3 L/hr in the case without a lower
    # bound. In the case before optimization, the mixer outlet flow
    # rate is 865 L/hr. In the case with a lower bound of 100, the
    # mixer outlet flow rate is 291 L/hr.
    #
    # We probably also need a performance equation for the filter
    # press  (sl_sep1) to determine the fraction of fluid entrained
    # as a function of the liquid to solid ratio.

    # Unfix HCl feed flow rates and concentrations
    for feed in [m.fs.acid_feed1, m.fs.acid_feed2, m.fs.acid_feed3]:
        feed.flow_vol.unfix()
        feed.conc_mass_comp[0, "H"].unfix()
        feed.conc_mass_comp[0, "Cl"].unfix()
        # According to Wikipedia, pH values of less than 0 are
        # impossible in an aqueous solution
        feed.properties[0].pH_phase["liquid"].setlb(0)

        @feed.Constraint(m.fs.time)
        def HCl_stoich_eqn(b, t):
            return (
                b.properties[t].conc_mol_comp["H"]
                == b.properties[t].conc_mol_comp["Cl"]
            )

        for condata in feed.HCl_stoich_eqn.values():
            set_scaling_factor(condata, 10)

    # Make extractant dosage a decision variable
    m.fs.rougher_org_make_up.conc_mass_comp[0, "DEHPA"].unfix()
    m.fs.rougher_org_make_up.properties[0].extractant_dosage.bounds = (3, 10)
    # m.fs.rougher_org_make_up.properties[0].extractant_dosage.fix(5)

    m.fs.cleaner_org_make_up.conc_mass_comp[0, "DEHPA"].unfix()
    m.fs.cleaner_org_make_up.properties[0].extractant_dosage.bounds = (3, 10)
    # m.fs.cleaner_org_make_up.properties[0].extractant_dosage.fix(5)

    # If the pH in the rougher scrub goes above 5 or 6, the equations get
    # extremely ill conditioned.
    m.fs.solex_rougher_scrub.mscontactor.aqueous[0.0, 1].pH_phase["liquid"].setub(6)

    solver = get_solver("ipopt_v2")
    solver.options.constr_viol_tol = 1e-8
    solver.options.max_iter = 300

    results = solver.solve(m, tee=True)

    if check_optimal_termination(results):
        display_results(m)
    else:
        print("Flowsheet optimization did not converge.")


def data_reconcilliation(m):
    m.fs.acid_feed1.flow_vol.unfix()
    m.obj = Objective(expr=m.fs.acid_feed1.flow_vol[0])
    m.fs.solex_rougher_scrub.mscontactor.aqueous[0.0, 1].pH_phase["liquid"].setub(3)

    solver = get_solver("ipopt_v2")
    solver.options.constr_viol_tol = 1e-8
    solver.options.max_iter = 300

    results = solver.solve(m, tee=True)

    if check_optimal_termination(results):
        display_results(m)
    else:
        print("Data reconcilliation optimization did not converge.")


if __name__ == "__main__":
    m, results = main()
    # warn(
    #     "Recent changes to this UKy flowsheet have made the underlying process more realistic, but the REE recovery values have fallen as a result."
    # )
    # warn(
    #     "Efforts are ongoing to increase the REE recovery while keeping the system as realistic as possible. https://github.com/prommis/prommis/issues/152 in the PrOMMiS repository is tracking the status of this issue."
    # )
    # optimize_model(m)
    data_reconcilliation(m)
