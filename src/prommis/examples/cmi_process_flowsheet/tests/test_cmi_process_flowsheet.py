#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    SolverFactory,
    TransformationFactory,
    assert_optimal_termination,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Product,
    Separator,
    StoichiometricReactor,
)

import pytest

from prommis.examples.cmi_process_flowsheet.cmi_process_flowsheet import (
    build,
    initialize_system,
    main,
    set_operation_conditions,
    solve_system,
)
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster


@pytest.fixture(scope="module")
def system_frame():
    m = build()
    set_operation_conditions(m)
    initialize_system(m)

    return m


@pytest.mark.component
def test_structural_issues(system_frame):
    model = system_frame
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


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
@pytest.mark.solver
def test_degrees_of_freedom(system_frame):
    m = system_frame
    if degrees_of_freedom(m) != 0:
        raise AssertionError("The degrees of freedom are not equal to 0.")


@pytest.mark.component
@pytest.mark.solver
def test_solve(system_frame):
    m = system_frame

    results = solve_system(m)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(system_frame):
    m = system_frame

    dt = DiagnosticsToolbox(m)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(system_frame):
    m = system_frame

    #  FEED (entering Stage 1)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Liq", "H2O"].value == pytest.approx(
        100, 1e-4
    )
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(0.92497, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Sol", "Cu"].value == pytest.approx(
        1e-5, 1e-4
    )
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(100, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.FEED.outlet.flow_mol_phase_comp[0, "Vap", "O2"].value == pytest.approx(
        100, 1e-4
    )

    # Copper Nitrate Dissolution Stoichiometric Reactor Specifications (Stage 1)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(93.525, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(9.2497e-7, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(0.46249, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(6.9373, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(0.46249, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.3165, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.8499, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-11, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(8.6330, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Dissolution.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(91.907, 1e-4)

    # S101 Filter Specifications (S/L filtration following Stage 1)
    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(93.525, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(9.3525e-06, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(9.2497e-12, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(9.2496e-07, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(4.6249e-06, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(0.46249, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(6.9373e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(6.9372, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(4.6249e-06, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(0.46249, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.3165e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(4.3165, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.0000e-10, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(9.9999e-06, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1.0000e-10, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(9.9999e-06, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(8.4276e-06, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.8499, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.8499e-07, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1.0000e-11, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1.0000e-18, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(8.6330, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(8.6330e-07, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(1.0000e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(1.0000e-12, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1.0000e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1.0000e-12, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1.0000e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1.0000e-12, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1.0000e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1.0000e-12, 1e-4)

    assert m.fs.S101.liq_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(9.1907e-05, 1e-4)
    assert m.fs.S101.sol_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(91.906, 1e-4)

    # pH Adjustment stage FEED Specifications (entering Stage 2)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(100, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjFeed.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(1e-5, 1e-4)

    # pH Adjustment stage mixer specifications (Stage 2)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(93.525, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1.4625e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(7.9373e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1.4625e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(5.3165e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.8500, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(8.6330, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(100, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.AdjMixer.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(0.00010191, 1e-4)

    # pH Adjustment stage reactor specifications (Stage 2)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(93.525, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1.4625e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(7.9373e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1.4625e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(8.6331, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.8500, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.8500e-06, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(8.6330e-06, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(68.551, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(31.449, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.Adjustment.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(0.00010191, 1e-4)

    # Precipitator Feed specifications (Stage 3)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(1e-05, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(100, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(100.00, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(1e-5, 1e-4)
    assert m.fs.PrecipFeed.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(1e-5, 1e-4)

    # Precipitator mixer (stage 3)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(93.525, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.4625e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.9373e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.4625e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(8.6331, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.8500, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(2e-5, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.1850e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(2e-5, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1.8633e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(168.55, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(31.449, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(100.00, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(3.0000e-05, 1e-4)
    assert m.fs.PrecipMixer.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(0.00011191, 1e-4)

    # Precipitator reactor (stage 3)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(141.62, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.4625e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.9373e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.4625e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(8.6331e-06, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.8500e-06, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(0.92500, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.1850e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(2e-5, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1.8633e-05, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(142.65, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(31.449, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(71.326, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(8.6331, 1e-4)
    assert m.fs.Precipitation.outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(0.00011191, 1e-4)

    # S102 Filter Specifications (S/L Filtration following Stage 3)
    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(141.62, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Liq", "H2O"
    ].value == pytest.approx(1.4162e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2.0000e-10, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2Fe14B"
    ].value == pytest.approx(2.0000e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.4625e-10, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu3(BO3)2"
    ].value == pytest.approx(2.4625e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.9373e-10, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu2O"
    ].value == pytest.approx(8.9372e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.4625e-10, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Cu"
    ].value == pytest.approx(2.4625e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(8.6331e-11, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Fe(OH)3"
    ].value == pytest.approx(8.6330e-06, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.8500e-11, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd(OH)3"
    ].value == pytest.approx(1.8500e-06, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(9.2500e-06, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Sol", "Nd2(C2O4)3 * 10H2O"
    ].value == pytest.approx(0.92499, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(84.276, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Cu(NO3)2"
    ].value == pytest.approx(8.4276e-06, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.1850e-05, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Nd(NO3)3"
    ].value == pytest.approx(1.1850e-12, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(2.0000e-05, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)2"
    ].value == pytest.approx(2.0000e-12, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1.8633e-05, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "Fe(NO3)3"
    ].value == pytest.approx(1.8633e-12, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(142.65, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4OH"
    ].value == pytest.approx(1.4265e-05, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(31.449, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "NH4NO3"
    ].value == pytest.approx(3.1449e-06, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(71.326, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "H2C2O4"
    ].value == pytest.approx(7.1326e-06, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(8.6331, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Aq", "(NH4)3[Fe(C2O4)3]"
    ].value == pytest.approx(8.6331e-07, 1e-4)

    assert m.fs.S102.liq_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(1.1191e-08, 1e-4)
    assert m.fs.S102.sol_outlet.flow_mol_phase_comp[
        0, "Vap", "O2"
    ].value == pytest.approx(0.00011190, 1e-4)

    # Calcinator
    assert m.fs.Calcination.flow_mol_comp_product[0, "Nd"].value == pytest.approx(
        0.92499, 1e-4
    )
