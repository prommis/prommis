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
    TransformationFactory,
    Var,
    check_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.environ import (
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core.util.scaling import constraint_scaling_transform, set_scaling_factor
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Separator,
    StoichiometricReactor,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

import prommis.examples.cmi_process_flowsheet.cmi_process_adjustment_rxn_prop_pack as adjustment_reaction_props
import prommis.examples.cmi_process_flowsheet.cmi_process_dissolution_rxn_prop_pack as dissolution_reaction_props
import prommis.examples.cmi_process_flowsheet.cmi_process_precipitation_rxn_prop_pack as precipitation_reaction_props
import prommis.examples.cmi_process_flowsheet.cmi_process_prop_pack as thermo_props
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster


def main():
    """
    Run the flowsheet by calling the appropriate functions in series.
    """
    m = build()

    set_operation_conditions(m)

    set_scaling(m)

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    if degrees_of_freedom(scaled_model) != 0:
        raise AssertionError(
            "The degrees of freedom are not equal to 0."
            "Check that the expected variables are fixed and unfixed."
            "For more guidance, run assert_no_structural_warnings from the IDAES DiagnosticToolbox "
        )

    # structural diagnostics check
    dt = DiagnosticsToolbox(scaled_model)
    dt.assert_no_structural_warnings()

    initialize_system(scaled_model)

    scaled_results = solve_system(scaled_model, False)

    if not check_optimal_termination(scaled_results):
        raise RuntimeError(
            "Solver failed to terminate with an optimal solution. Please check the solver logs for more details"
        )

    # numerical diagnostics test
    dt.assert_no_numerical_warnings()

    res = scaling.propagate_solution(scaled_model, m)

    display_results(m)

    return m, res


def build():
    """
    Build and connect the unit model blocks present in the CMI Process.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.thermo_params = GenericParameterBlock(**thermo_props.thermo_config)
    m.fs.dissolution_reaction_params = GenericReactionParameterBlock(
        property_package=m.fs.thermo_params, **dissolution_reaction_props.config_dict
    )
    m.fs.adjustment_reaction_params = GenericReactionParameterBlock(
        property_package=m.fs.thermo_params, **adjustment_reaction_props.config_dict
    )
    m.fs.precipitation_reaction_params = GenericReactionParameterBlock(
        property_package=m.fs.thermo_params, **precipitation_reaction_props.config_dict
    )

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_solid = PrecipitateParameters()
    m.fs.prop_liquid = AqueousParameter()

    ### Feed Stream (entering stage 1)
    m.fs.FEED = Feed(property_package=m.fs.thermo_params)

    ### Copper(II) Nitrate Dissolution Stage (stage 1)
    m.fs.Dissolution = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.dissolution_reaction_params,
        has_heat_of_reaction=False,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    ### S/L Filter following dissolution (following stage 1)
    m.fs.S101 = Separator(
        property_package=m.fs.thermo_params,
        split_basis=SplittingType.phaseFlow,
        outlet_list=["liq_outlet", "sol_outlet"],
        ideal_separation=False,
        has_phase_equilibrium=False,
    )

    ### pH Adjustment Feed (entering stage 2)
    m.fs.AdjFeed = Feed(property_package=m.fs.thermo_params)

    ### pH Adjustment Mixer (stage 2)
    m.fs.AdjMixer = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["reactant_feed", "separator_stream"],
    )

    ### pH Adjustment Stage Reactor (stage 2)
    m.fs.Adjustment = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.adjustment_reaction_params,
        has_heat_of_reaction=False,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    ### Precipitation Reactant Feed Stream (stage 3)
    m.fs.PrecipFeed = Feed(property_package=m.fs.thermo_params)

    ### Precipitation Mixer (stage 3)
    m.fs.PrecipMixer = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["reactant_feed", "adjustment_stream"],
    )

    ### Precipitation reactor (stage 3)
    m.fs.Precipitation = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.precipitation_reaction_params,
        has_heat_of_reaction=False,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    ### S/L Filter (following stage 3)
    m.fs.S102 = Separator(
        property_package=m.fs.thermo_params,
        split_basis=SplittingType.phaseFlow,
        outlet_list=["liq_outlet", "sol_outlet"],
        ideal_separation=False,
        has_phase_equilibrium=False,
    )

    ### Calcination unit
    m.fs.Calcination = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate_solid=m.fs.prop_solid,
        property_package_precipitate_liquid=m.fs.prop_liquid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        metal_list=["Nd"],
    )

    # Connect arc from Feed to Dissolution Stage
    m.fs.FEED_Diss = Arc(source=m.fs.FEED.outlet, destination=m.fs.Dissolution.inlet)

    # Connect arc from Dissolution Stage to S101
    m.fs.Diss_S101 = Arc(source=m.fs.Dissolution.outlet, destination=m.fs.S101.inlet)

    # Connect arc from S101 to pH adjustment stage
    m.fs.S101_AdjMixer = Arc(
        source=m.fs.S101.liq_outlet, destination=m.fs.AdjMixer.separator_stream
    )

    # Connect arc from pH adjustment feed to mixer
    m.fs.AdjFeed_AdjMixer = Arc(
        source=m.fs.AdjFeed.outlet, destination=m.fs.AdjMixer.reactant_feed
    )

    # Connect arc from pH adjustment mixer to reactor
    m.fs.AdjMixer_Adjustment = Arc(
        source=m.fs.AdjMixer.outlet, destination=m.fs.Adjustment.inlet
    )

    # Connect arc from pH adjustment reactor to precipitation stage mixer
    m.fs.Adjustment_PrecipMixer = Arc(
        source=m.fs.Adjustment.outlet, destination=m.fs.PrecipMixer.adjustment_stream
    )

    # Connect arc from precipitation feed to mixer
    m.fs.PrecipFeed_PrecipMixer = Arc(
        source=m.fs.PrecipFeed.outlet, destination=m.fs.PrecipMixer.reactant_feed
    )

    # Connect arc from precipitation mixer to reactor
    m.fs.PrecipMixer_Precipitation = Arc(
        source=m.fs.PrecipMixer.outlet, destination=m.fs.Precipitation.inlet
    )

    # Connect arc from precipitation reactor to filter
    m.fs.Precipitation_S102 = Arc(
        source=m.fs.Precipitation.outlet, destination=m.fs.S102.inlet
    )

    # Connect filter to Calcinator
    m.Calc_feed_con = Constraint(
        expr=m.fs.Calcination.flow_mol_comp_feed[0, "Nd"]
        == m.fs.S102.sol_outlet.flow_mol_phase_comp[0, "Sol", "Nd2(C2O4)3 * 10H2O"]
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_operation_conditions(m):
    """
    Set the operating conditions of the flowsheet such that the degrees of freedom are zero.

    Args:
        m: pyomo model
    """

    #########  FEED Specification (entering Stage 1)
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"].fix(1 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "Cu_2+"].fix(34 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "NO3_-"].fix(68 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Vap", "O2"].fix(50 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Liq", "H2O"].fix(300 * pyunits.mol / pyunits.s)

    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "Nd_3+"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "Fe_2+"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "Fe_3+"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "NH4_+"].fix(1e-5 * pyunits.mol / pyunits.s)

    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "OH_-"].fix(1e-5 * pyunits.mol / pyunits.s)

    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Cu3(BO3)2"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Cu2O"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Cu"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "H2C2O4"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.FEED.flow_mol_phase_comp[0, "Aq", "C2O4_2-"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Nd2(C2O4)3 * 10H2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )

    # set pressure of feed
    m.fs.FEED.pressure.fix(1 * pyunits.atm)

    # set temperature of feed
    m.fs.FEED.temperature.fix(298.15 * pyunits.K)
    ######### Copper Nitrate Dissolution Stoichiometric Reactor Specifications (Stage 1)
    m.fs.Dissolution.Nd_magnet_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Dissolution.Nd_magnet_conv_constraint = Constraint(
        expr=m.fs.Dissolution.Nd_magnet_conversion
        * m.fs.Dissolution.inlet.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"]
        == (
            m.fs.Dissolution.inlet.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"]
            - m.fs.Dissolution.outlet.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"]
        )
    )
    m.fs.Dissolution.Nd_magnet_conversion.fix(0.999999)

    m.fs.Dissolution.Iron2_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Dissolution.Iron2_conv_constraint = Constraint(
        expr=m.fs.Dissolution.Iron2_conversion
        * m.fs.Dissolution.inlet.flow_mol_phase_comp[0, "Aq", "Fe_2+"]
        == (
            m.fs.Dissolution.inlet.flow_mol_phase_comp[0, "Aq", "Fe_2+"]
            - m.fs.Dissolution.outlet.flow_mol_phase_comp[0, "Aq", "Fe_2+"]
        )
    )
    m.fs.Dissolution.Iron2_conversion.fix(0.999999)

    # set temperature
    m.fs.Dissolution.outlet.temperature.fix(343.15 * pyunits.K)
    # set upper bound
    m.fs.Dissolution.control_volume.properties_in[0.0].temperature.setub(600)

    ######### S101 Filter Specifications (S/L filtration following Stage 1)
    m.fs.S101.split_fraction[0, "liq_outlet", "Liq"].fix(0.9999999)
    m.fs.S101.split_fraction[0, "liq_outlet", "Aq"].fix(0.9999999)
    m.fs.S101.split_fraction[0, "sol_outlet", "Sol"].fix(0.99999)
    m.fs.S101.split_fraction[0, "sol_outlet", "Vap"].fix(0.999999)

    ######### pH Adjustment stage FEED Specifications (entering Stage 2)
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "Cu_2+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "NO3_-"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Vap", "O2"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )

    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "Nd_3+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "Fe_2+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "Fe_3+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "NH4_+"].fix(35 * pyunits.mol / pyunits.s)

    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "OH_-"].fix(35 * pyunits.mol / pyunits.s)

    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Cu3(BO3)2"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Cu2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Cu"].fix(1e-5 * pyunits.mol / pyunits.s)
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "H2C2O4"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Aq", "C2O4_2-"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.flow_mol_phase_comp[0, "Sol", "Nd2(C2O4)3 * 10H2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.AdjFeed.temperature.fix(298.15 * pyunits.K)
    m.fs.AdjFeed.pressure.fix(1 * pyunits.atm)

    ######## pH Adjustment stage mixer specifications (Stage 2)
    # set upper bounds on temperature
    m.fs.AdjMixer.reactant_feed_state[0].temperature.setub(600)
    m.fs.AdjMixer.separator_stream_state[0].temperature.setub(600)
    m.fs.AdjMixer.mixed_state[0].temperature.setub(600)

    ######### pH Adjustment stage reactor specifications (Stage 2)
    m.fs.Adjustment.Fe_nitrate_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Adjustment.Fe_nitrate_conv_constraint = Constraint(
        expr=m.fs.Adjustment.Fe_nitrate_conversion
        * m.fs.Adjustment.inlet.flow_mol_phase_comp[0, "Aq", "Fe_3+"]
        == (
            m.fs.Adjustment.inlet.flow_mol_phase_comp[0, "Aq", "Fe_3+"]
            - m.fs.Adjustment.outlet.flow_mol_phase_comp[0, "Aq", "Fe_3+"]
        )
    )
    m.fs.Adjustment.Fe_nitrate_conversion.fix(0.999999)

    m.fs.Adjustment.Nd_nitrate_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Adjustment.Nd_nitrate_conv_constraint = Constraint(
        expr=m.fs.Adjustment.Nd_nitrate_conversion
        * m.fs.Adjustment.inlet.flow_mol_phase_comp[0, "Aq", "Nd_3+"]
        == (
            m.fs.Adjustment.inlet.flow_mol_phase_comp[0, "Aq", "Nd_3+"]
            - m.fs.Adjustment.outlet.flow_mol_phase_comp[0, "Aq", "Nd_3+"]
        )
    )
    m.fs.Adjustment.Nd_nitrate_conversion.fix(0.999999)

    m.fs.Adjustment.outlet.temperature.fix(333.15 * pyunits.K)
    # set upper bound
    m.fs.Adjustment.control_volume.properties_in[0].temperature.setub(600)

    ######### Precipitator Feed specifications (Stage 3)
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "Cu_2+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "NO3_-"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Vap", "O2"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "Nd_3+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "Fe_2+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "Fe_3+"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "NH4_+"].fix(
        35 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "OH_-"].fix(
        35 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Cu3(BO3)2"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Cu2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Cu"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "H2C2O4"].fix(
        35 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Aq", "C2O4_2-"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )
    m.fs.PrecipFeed.flow_mol_phase_comp[0, "Sol", "Nd2(C2O4)3 * 10H2O"].fix(
        1e-5 * pyunits.mol / pyunits.s
    )

    m.fs.PrecipFeed.temperature.fix(298.15 * pyunits.K)
    m.fs.PrecipFeed.pressure.fix(1 * pyunits.atm)

    ######## Precipitator mixer (stage 3)
    m.fs.PrecipMixer.reactant_feed_state[0].temperature.setub(600)
    m.fs.PrecipMixer.adjustment_stream_state[0].temperature.setub(600)
    m.fs.PrecipMixer.mixed_state[0].temperature.setub(600)

    ######## Precipitator reactor (stage 3)
    m.fs.Precipitation.Fe_hydroxide_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Precipitation.Fe_hydroxide_conv_constraint = Constraint(
        expr=m.fs.Precipitation.Fe_hydroxide_conversion
        * m.fs.Precipitation.inlet.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"]
        == (
            m.fs.Precipitation.inlet.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"]
            - m.fs.Precipitation.outlet.flow_mol_phase_comp[0, "Sol", "Fe(OH)3"]
        )
    )
    m.fs.Precipitation.Fe_hydroxide_conversion.fix(0.999999)

    m.fs.Precipitation.Nd_hydroxide_conversion = Var(
        initialize=1, bounds=(0, 1), units=pyunits.dimensionless
    )

    m.fs.Precipitation.Nd_hydroxide_conv_constraint = Constraint(
        expr=m.fs.Precipitation.Nd_hydroxide_conversion
        * m.fs.Precipitation.inlet.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"]
        == (
            m.fs.Precipitation.inlet.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"]
            - m.fs.Precipitation.outlet.flow_mol_phase_comp[0, "Sol", "Nd(OH)3"]
        )
    )
    m.fs.Precipitation.Nd_hydroxide_conversion.fix(0.999999)

    m.fs.Precipitation.outlet.temperature.fix(333.15 * pyunits.K)
    # set upper bound
    m.fs.Precipitation.control_volume.properties_in[0].temperature.setub(600)

    ######### S102 Filter Specifications (S/L Filtration following Stage 3)
    m.fs.S102.split_fraction[0, "liq_outlet", "Liq"].fix(0.999999)
    m.fs.S102.split_fraction[0, "liq_outlet", "Aq"].fix(0.999999)
    m.fs.S102.split_fraction[0, "sol_outlet", "Sol"].fix(0.99999)
    m.fs.S102.split_fraction[0, "sol_outlet", "Vap"].fix(0.9999)

    ######### Calcinator
    m.fs.Calcination.deltaP.fix(0)
    m.fs.Calcination.gas_inlet.temperature.fix(873.15)
    m.fs.Calcination.gas_inlet.pressure.fix(101325)

    # Gas mole flow rate
    fgas = 100
    # Pure O2 Gas
    gas_comp = {
        "O2": 0.2095,
        "H2O": 0.01,
        "CO2": 0.0004,
        "N2": 0.7801,
    }
    for i, v in gas_comp.items():
        m.fs.Calcination.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.Calcination.gas_inlet.flow_mol.fix(fgas)

    # fix outlet product temperature
    m.fs.Calcination.gas_outlet.temperature.fix(873.15)

    # defined by the solid precipitate property package
    m.fs.Calcination.solid_in[0].temperature.fix(333.15)

    # no liquid entering calcinator. Set to lower bound.
    m.fs.Calcination.liquid_in[0].flow_vol.fix(1e-5)  # in L/hr
    m.fs.Calcination.liquid_in[0].conc_mass_comp.fix(1e-5)
    m.fs.Calcination.liquid_in[0].conc_mass_comp["H2O"].fix(1e-5)  # mg/L
    m.fs.Calcination.frac_comp_recovery.fix(1)


def set_scaling(m):
    """
    Set the scaling factors to improve solver performance.

    Args:
        m: pyomo model
    """
    set_scaling_factor(
        m.fs.PrecipMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "H2C2O4"], 1e4
    )
    set_scaling_factor(
        m.fs.PrecipMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "C2O4_2-"], 1e4
    )
    set_scaling_factor(
        m.fs.PrecipMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "Nd_3+"], 1e4
    )

    constraint_scaling_transform(
        m.fs.Dissolution.control_volume.enthalpy_balances[0.0], 1e-6
    )
    constraint_scaling_transform(
        m.fs.Adjustment.control_volume.enthalpy_balances[0.0], 1e-6
    )
    constraint_scaling_transform(m.fs.AdjMixer.enthalpy_mixing_equations[0.0], 1e-6)
    constraint_scaling_transform(m.fs.PrecipMixer.enthalpy_mixing_equations[0.0], 1e-6)
    constraint_scaling_transform(
        m.fs.Precipitation.control_volume.enthalpy_balances[0.0], 1e-6
    )

    set_scaling_factor(
        m.fs.AdjMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "H2C2O4"], 1e4
    )
    set_scaling_factor(
        m.fs.AdjMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "C2O4_2-"], 1e4
    )
    set_scaling_factor(
        m.fs.AdjMixer.mixed_state[0.0].mole_frac_phase_comp["Aq", "Nd_3+"], 1e4
    )

    set_scaling_factor(
        m.fs.Adjustment.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "H2C2O4"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Adjustment.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "C2O4_2-"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Adjustment.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "Nd_3+"
        ],
        1e4,
    )

    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "H2C2O4"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "C2O4_2-"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Aq", "H2C2O4"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Aq", "C2O4_2-"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_in[0.0].mole_frac_phase_comp[
            "Aq", "Nd_3+"
        ],
        1e4,
    )
    set_scaling_factor(
        m.fs.Precipitation.control_volume.properties_out[0.0].mole_frac_phase_comp[
            "Aq", "Nd_3+"
        ],
        1e4,
    )

    set_scaling_factor(
        m.fs.S102.liq_outlet_state[0.0].mole_frac_phase_comp_apparent["Vap", "O2"], 1e-4
    )


def initialize_system(m):
    """
    Initialize system.

    Args:
        m: pyomo model
    """

    ### Initialize Feed and propagate state to Dissolution Stage
    m.fs.FEED.initialize()
    propagate_state(arc=m.fs.FEED_Diss)

    ### Initialize Dissolution Stage and propagate state to S101
    m.fs.Dissolution.initialize()
    propagate_state(arc=m.fs.Diss_S101)

    ### Initialize S101 and propagate state to pH adjustment stage mixer
    m.fs.S101.initialize()
    propagate_state(arc=m.fs.S101_AdjMixer)

    # # Initialize pH Adjustment feed and propagate state to mixer
    m.fs.AdjFeed.initialize()
    propagate_state(arc=m.fs.AdjFeed_AdjMixer)

    # # Initialize pH Adjustment mixer and propagate state to reactor
    m.fs.AdjMixer.initialize()
    propagate_state(arc=m.fs.AdjMixer_Adjustment)

    # # Initialize pH Adjustment reactor and propagate state to precipitation stage mixer
    m.fs.Adjustment.initialize()
    propagate_state(arc=m.fs.Adjustment_PrecipMixer)

    # Initialize precipitation stage feed and propagate state to precipitation stage mixer
    m.fs.PrecipFeed.initialize()
    propagate_state(arc=m.fs.PrecipFeed_PrecipMixer)

    # Initialize precipitation stage feed and propagate state to precipitation stage reactor
    m.fs.PrecipMixer.initialize()
    propagate_state(arc=m.fs.PrecipMixer_Precipitation)

    # Initialize precipitation stage reactor
    m.fs.Precipitation.initialize()
    propagate_state(arc=m.fs.Precipitation_S102)

    ### Initialize S102
    m.fs.S102.initialize()

    # Initialize precipitation stage calcinator
    m.fs.Calcination.initialize()


def solve_system(m, tee=False):
    """
    Args:
        m: pyomo model
        tee: boolean indicator to stream IPOPT solution
    """
    # Solve flowsheet
    solver_obj = SolverFactory(
        "ipopt",
        options={
            "nlp_scaling_method": "user-scaling",
            "tol": 1e-6,
            "max_iter": 1000,
        },
    )

    results = solver_obj.solve(m, tee=tee)

    return results


def display_results(m):
    """
    Print key flowsheet outputs.

    Args:
        m: pyomo model
    """

    print(
        "Molar flowrate of REPM entering process: {:0.2f} mol/s".format(
            value(m.fs.FEED.flow_mol_phase_comp[0, "Sol", "Nd2Fe14B"])
        )
    )

    m.fs.FEED.report()
    m.fs.Dissolution.report()
    m.fs.S101.report()
    m.fs.AdjFeed.report()
    m.fs.AdjMixer.report()
    m.fs.Adjustment.report()
    m.fs.PrecipFeed.report()
    m.fs.PrecipMixer.report()
    m.fs.Precipitation.report()
    m.fs.S102.report()

    print(
        "\nMolar Flowrate of product Nd2O3 recovered: {:0.4f} mol/s".format(
            value(m.fs.Calcination.flow_mol_comp_product[0, "Nd"])
        )
    )


if __name__ == "__main__":
    main()
