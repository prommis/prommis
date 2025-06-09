#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import os
from pyomo.environ import (
    SolverFactory,
    TransformationFactory,
    assert_optimal_termination,
    value,
    ConcreteModel,
    units,
)

from idaes.core.util import DiagnosticsToolbox, from_json
from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_train import LeachingTrain
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.leaching.leach_flowsheet_dynamic import (
    build_model,
    discretization,
    copy_first_steady_state,
    set_inputs,
)
from prommis.leaching.leach_flowsheet import (
    set_inputs as set_inputs_steady_state,
    set_scaling,
)


@pytest.fixture(scope="module")
def model():
    m = build_model(time_duration=24, number_of_tanks=1)
    discretization(m)
    current_directory = os.path.dirname(__file__)
    parent_directory = os.path.dirname(current_directory)
    json_file_path = os.path.join(parent_directory, "leaching.json")
    from_json(m, fname=json_file_path)
    copy_first_steady_state(m)
    set_inputs(m, perturb_time=12)

    return m


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_solve(model):

    # Solve scaled model
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=False)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    time_duration = 24
    conversion = {
        "Al2O3": 0.03003765546930494,
        "CaO": 0.5145052459223426,
        "Ce2O3": 0.4207435931354908,
        "Dy2O3": 0.22365403062817607,
        "Fe2O3": 0.16338507097541013,
        "Gd2O3": 0.6106358012228825,
        "La2O3": 0.38342494088341084,
        "Nd2O3": 0.45795409973702533,
        "Pr2O3": 0.4842509098265882,
        "Sc2O3": 0.05486783341762388,
        "Sm2O3": 0.23377407833771635,
        "Y2O3": 0.15092538416884166,
        "inerts": 0.0,
    }
    recovery = {
        "inerts": 2.2679023103933105e-06,
        "Sc2O3": 5.848925760639934,
        "Y2O3": 15.06074131348737,
        "La2O3": 37.21499032227469,
        "Ce2O3": 41.940618340188195,
        "Pr2O3": 47.719654567828115,
        "Nd2O3": 46.549362626788906,
        "Sm2O3": 23.64501254412712,
        "Gd2O3": 65.48567880986579,
        "Dy2O3": 21.76629835153545,
        "Al2O3": 4.374322459626034,
        "CaO": 45.117626139151476,
        "Fe2O3": 17.958786657478853,
    }

    for k, v in model.fs.leach.mscontactor.solid[0, 1].conversion.items():
        f_in = model.fs.leach.solid_inlet.flow_mass[time_duration]
        f_out = model.fs.leach.solid_outlet.flow_mass[time_duration]
        x_in = model.fs.leach.solid_inlet.mass_frac_comp[time_duration, k]
        x_out = model.fs.leach.solid_outlet.mass_frac_comp[time_duration, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        if k[0] == time_duration:
            assert value(v) == pytest.approx(conversion[k], rel=1e-5, abs=1e-6)

        if k[0] == time_duration:
            assert r == pytest.approx(recovery[k[1]], rel=1e-5, abs=1e-6)


@pytest.fixture(scope="module")
def model2():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    m.fs.leach = LeachingTrain(
        number_of_tanks=1,
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
        has_holdup=True,
    )

    set_inputs_steady_state(m)

    m.fs.leach.mscontactor.volume.fix(100 * units.gallon)
    m.fs.leach.mscontactor.volume_frac_stream[:, :, "liquid"].fix(0.5)

    return m


@pytest.mark.unit
def test_structural_issues2(model2):
    dt = DiagnosticsToolbox(model2)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solve2(model2):

    set_scaling(model2)
    # Create a scaled version of the model to solve
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model2, rename=False)

    initializer = model2.fs.leach.default_initializer()
    try:
        initializer.initialize(scaled_model.fs.leach)
    except:
        pass

    # Solve scaled model
    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=False)

    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, model2)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues2(model2):
    dt = DiagnosticsToolbox(model2)
    dt.assert_no_numerical_warnings()
