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
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.unit_models import (
    Feed, 
    StoichiometricReactor, 
    Separator,
    Mixer,
    Product
)


import pytest

from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster

from prommis.examples.CMI_Process_Flowsheet.CMI_Process_Flowsheet import (
    main,
    build,
    set_operation_conditions,
    initialize_system,
    solve_system,
)


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



