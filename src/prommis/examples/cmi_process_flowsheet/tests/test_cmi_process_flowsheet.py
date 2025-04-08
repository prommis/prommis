#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    TransformationFactory,
    assert_optimal_termination,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
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

import pytest

from prommis.examples.cmi_process_flowsheet.cmi_process_flowsheet import (
    build,
    display_results,
    initialize_system,
    set_operation_conditions,
    set_scaling,
    solve_system,
)
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster

# global var for storing scaled model
_scaled_model = None


@pytest.fixture(scope="module")
def system_frame():
    m = build()

    return m


@pytest.mark.unit
def test_build_flowsheet(system_frame):
    m = system_frame
    assert isinstance(m.fs, FlowsheetBlock)

    # Check property packages
    assert isinstance(m.fs.thermo_params, GenericParameterBlock)
    assert isinstance(m.fs.dissolution_reaction_params, GenericReactionParameterBlock)
    assert isinstance(m.fs.adjustment_reaction_params, GenericReactionParameterBlock)
    assert isinstance(m.fs.precipitation_reaction_params, GenericReactionParameterBlock)

    assert isinstance(m.fs.prop_gas, GenericParameterBlock)
    assert isinstance(m.fs.prop_solid, PrecipitateParameters)
    assert isinstance(m.fs.prop_liquid, AqueousParameter)

    # Feed Stream (entering stage 1)
    assert isinstance(m.fs.FEED, Feed)

    # Copper(II) Nitrate Dissolution Stage (stage 1)
    assert isinstance(m.fs.Dissolution, StoichiometricReactor)

    # S/L Filter following dissolution (following stage 1)
    assert isinstance(m.fs.S101, Separator)

    # pH Adjustment Feed (entering stage 2)
    assert isinstance(m.fs.AdjFeed, Feed)

    # pH Adjustment Mixer (stage 2)
    assert isinstance(m.fs.AdjMixer, Mixer)

    # pH Adjustment Stage Reactor (stage 2)
    assert isinstance(m.fs.Adjustment, StoichiometricReactor)

    # Precipitation Reactant Feed Stream (stage 3)
    assert isinstance(m.fs.PrecipFeed, Feed)

    # Precipitation Mixer (stage 3)
    assert isinstance(m.fs.PrecipMixer, Mixer)

    # Precipitation Reactor (stage 3)
    assert isinstance(m.fs.Precipitation, StoichiometricReactor)

    # S/L Filter (following stage 3)
    assert isinstance(m.fs.S102, Separator)

    # Calcination Unit
    assert isinstance(m.fs.Calcination, REEOxalateRoaster)

    # Flowsheet connections
    assert isinstance(m.fs.FEED_Diss, Arc)
    assert isinstance(m.fs.Diss_S101, Arc)
    assert isinstance(m.fs.S101_AdjMixer, Arc)
    assert isinstance(m.fs.AdjFeed_AdjMixer, Arc)
    assert isinstance(m.fs.AdjMixer_Adjustment, Arc)
    assert isinstance(m.fs.Adjustment_PrecipMixer, Arc)
    assert isinstance(m.fs.PrecipFeed_PrecipMixer, Arc)
    assert isinstance(m.fs.PrecipMixer_Precipitation, Arc)
    assert isinstance(m.fs.Precipitation_S102, Arc)


@pytest.mark.component
def test_structural_issues(system_frame):
    global _scaled_model
    model = system_frame
    set_operation_conditions(model)

    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.fixture(scope="module")
def scaled_model(system_frame):
    model = system_frame
    set_operation_conditions(model)
    set_scaling(model)

    # Apply scaling transformation
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model, rename=False)

    # Initialize the scaled model
    initialize_system(scaled_model)

    return scaled_model


@pytest.mark.component
def test_initialized_system(scaled_model):
    initialize_system(scaled_model)


@pytest.mark.component
@pytest.mark.solver
def test_solve(scaled_model, system_frame):
    model = system_frame

    results = solve_system(scaled_model)

    scaling = TransformationFactory("core.scale_model")
    scaling.propagate_solution(scaled_model, model)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(scaled_model):
    # Use the already scaled model from test_solve
    dt = DiagnosticsToolbox(scaled_model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(system_frame):
    m = system_frame

    #  FEED (entering Stage 1)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Liq", "H2O"].value == pytest.approx(
        300, rel=1e-4
    )
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Sol", "Cu"].value == pytest.approx(
        1e-5, rel=1e-4
    )
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Aq", "OH_-"].value == pytest.approx(
        1e-5, rel=1e-4
    )
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(34, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Vap", "O2"].value == pytest.approx(
        50, rel=1e-4
    )

    # Copper Nitrate Dissolution Stoichiometric Reactor Specifications (Stage 1)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(293, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-6, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(0.50001, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(7.5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(0.50001, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.6667, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(2, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1.0001e-11, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3333, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(41.25, rel=1e-4)

    # S101 Filter Specifications (S/L filtration following Stage 1)
    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(293, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(2.93e-5, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1.0001e-11, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(9.9999e-7, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(5.0001e-6, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(0.5, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(7.5e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(7.4999, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(5.0001e-6, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(0.5000, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.6667e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.6666, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-10, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(9.9999e-6, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-10, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(9.9999e-6, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-12, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(1.0015e-12, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(1.7e-6, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(6.8e-6, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(2, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(2e-7, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1.0001e-11, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1.4995e-15, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3333, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3333e-7, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(1.0015e-12, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1.0015e-12, rel=1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(4.125e-5, rel=1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(41.25, rel=1e-4)

    # pH Adjustment stage FEED Specifications (entering Stage 2)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(1e-5, rel=1e-4)

    # pH Adjustment stage mixer specifications (Stage 2)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(293, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1.5e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.5e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1.5e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(5.6667e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(2, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3333, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(5.1250e-5, rel=1e-4)

    # pH Adjustment stage reactor specifications (Stage 2)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(293, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1.5e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.5e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1.5e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(9.3334, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(2, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(0.99994, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(2e-6, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3333e-6, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(5.1250e-5, rel=1e-4)

    # Precipitator Feed specifications (Stage 3)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(1e-5, rel=1e-4)

    # Precipitator mixer (stage 3)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(293, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.5e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(9.5e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.5e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(9.3334, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(2, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(35, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(36, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1.2e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(1.9333e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(70, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(3e-5, rel=1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(6.1250e-5, rel=1e-4)

    # Precipitator reactor (stage 3)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(345, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.5e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(9.5e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.5e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(9.3334e-6, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(2e-6, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(3.9998, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(26.667, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1.2e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3334, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(70, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(28, rel=1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(6.1250e-5, rel=1e-4)

    # S102 Filter Specifications (S/L Filtration following Stage 3)
    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(345, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(0.000345, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-10, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.5e-10, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.5e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(9.5e-10, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(9.4999e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.5e-10, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.5e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(9.3336e-11, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(9.3333e-6, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(2.0002e-11, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(2e-6, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(3.9998, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(3.9998e-6, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(26.667, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "OH_-"
    ].value == pytest.approx(2.6667e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(17, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu_2+"
    ].value == pytest.approx(1.7e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(68, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NO3_-"
    ].value == pytest.approx(6.8e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1.2e-5, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd_3+"
    ].value == pytest.approx(1.2002e-11, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(2e-5, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_2+"
    ].value == pytest.approx(2.0002e-11, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3334, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe_3+"
    ].value == pytest.approx(9.3334e-6, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(70, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4_+"
    ].value == pytest.approx(7e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(28, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "C2O4_2-"
    ].value == pytest.approx(2.8e-5, rel=1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(6.125e-9, rel=1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(6.1244e-5, rel=1e-4)

    # Calcinator
    assert m.fs.Calcination.flow_mol_comp_product[0, "Nd"].value == pytest.approx(
        1, 1e-4
    )


@pytest.mark.unit
def test_display(system_frame):
    model = system_frame
    display_results(model)
