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
from idaes.core.solvers import get_solver

import pytest

from prommis.leaching.leach_train import LeachingTrain
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.leaching.leach_flowsheet_dynamic import (
    build_model,
    discretization,
    copy_first_steady_state,
    set_inputs,
    set_scaling as set_scaling_dynamic,
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

    set_scaling_dynamic(m)

    return m


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_solve(model):

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model, rename=False)

    # Solve scaled model
    solver = get_solver("ipopt_v2")
    results = solver.solve(scaled_model, tee=False)
    scaling.propagate_solution(scaled_model, model)

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
    conversion_comp = {
        "Al2O3": 0.0355181940608198,
        "CaO": 0.4995144924025545,
        "Ce2O3": 0.43107111613774735,
        "Dy2O3": 0.22437005073278507,
        "Fe2O3": 0.17529336473370566,
        "Gd2O3": 0.6565721684072603,
        "La2O3": 0.3851408193127902,
        "Nd2O3": 0.4749501283359911,
        "Pr2O3": 0.48959926690601413,
        "Sc2O3": 0.05777096543627356,
        "Sm2O3": 0.23932037012182222,
        "Y2O3": 0.15465402409264634,
        "inerts": 1e-6,
    }
    recovery = {
        "inerts": -8.35084156178796,
        "Sc2O3": -2.0914109302357486,
        "Y2O3": 8.405960505624089,
        "La2O3": 33.37942371004168,
        "Ce2O3": 38.35601500079665,
        "Pr2O3": 44.697595733104656,
        "Nd2O3": 43.110347653703165,
        "Sm2O3": 17.579639523490776,
        "Gd2O3": 62.789268219398494,
        "Dy2O3": 15.959758216112995,
        "Al2O3": -4.502519846955044,
        "CaO": 45.7719198343315,
        "Fe2O3": 10.642252669700325,
    }

    for k, v in model.fs.leach.mscontactor.solid[time_duration, 1].conversion_comp.items():
        f_in = model.fs.leach.solid_inlet.flow_mass[time_duration]
        f_out = model.fs.leach.solid_outlet.flow_mass[time_duration]
        x_in = model.fs.leach.solid_inlet.mass_frac_comp[time_duration, k]
        x_out = model.fs.leach.solid_outlet.mass_frac_comp[time_duration, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        assert value(v) == pytest.approx(conversion_comp[k], rel=1e-5, abs=1e-6)
        assert r == pytest.approx(recovery[k], rel=1e-5, abs=1e-6)


@pytest.fixture(scope="module")
def model2():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactionParameterBlock()

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
    initializer.initialize(scaled_model.fs.leach)

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


@pytest.mark.component
@pytest.mark.solver
def test_solution2(model2):
    conversion_comp = {
        "Al2O3": 0.03003765546937117,
        "CaO": 0.5145052459214582,
        "Ce2O3": 0.4207435931312086,
        "Dy2O3": 0.2236540306278243,
        "Fe2O3": 0.16338507097514096,
        "Gd2O3": 0.6106358012135612,
        "La2O3": 0.3834249408816641,
        "Nd2O3": 0.45795409972930473,
        "Pr2O3": 0.48425090967291895,
        "Sc2O3": 0.05486783341758024,
        "Sm2O3": 0.23377407833714964,
        "Y2O3": 0.15092538416861428,
        "inerts": 0.0,
    }
    recovery = {
        "inerts": 0.0,
        "Sc2O3": 5.486783341758061,
        "Y2O3": 15.092538416861444,
        "La2O3": 38.342494088166426,
        "Ce2O3": 42.074359313120866,
        "Pr2O3": 48.42509096729191,
        "Nd2O3": 45.79540997293049,
        "Sm2O3": 23.37740783371498,
        "Gd2O3": 61.06358012135613,
        "Dy2O3": 22.365403062782462,
        "Al2O3": 3.0037655469370983,
        "CaO": 51.45052459214581,
        "Fe2O3": 16.33850709751409,
    }

    for k, v in model2.fs.leach.mscontactor.solid[0, 1].conversion_comp.items():
        f_in = model2.fs.leach.solid_inlet.flow_mass[0]
        f_out = model2.fs.leach.solid_outlet.flow_mass[0]
        x_in = model2.fs.leach.solid_inlet.mass_frac_comp[0, k]
        x_out = model2.fs.leach.solid_outlet.mass_frac_comp[0, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        assert value(v) == pytest.approx(conversion_comp[k], rel=1e-5, abs=1e-6)
        assert r == pytest.approx(recovery[k], rel=1e-5, abs=1e-6)


@pytest.fixture(scope="module")
def model3():

    m = build_model(time_duration=24, number_of_tanks=2)
    discretization(m)
    current_directory = os.path.dirname(__file__)
    parent_directory = os.path.dirname(current_directory)
    json_file_path = os.path.join(parent_directory, "leaching2.json")
    from_json(m, fname=json_file_path)
    copy_first_steady_state(m)
    set_inputs(m, perturb_time=12)

    # Fixing the volume of the leach reactor
    m.fs.leach.volume.fix(50 * units.gallon)
    m.fs.leach.mscontactor.volume.fix(50 * units.gallon)

    return m


@pytest.mark.unit
def test_structural_issues3(model3):
    dt = DiagnosticsToolbox(model3)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_solve3(model3):

    set_scaling_dynamic(model3)

    # Solve scaled model
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model3, rename=False)
    solver = SolverFactory("ipopt_v2")
    results = solver.solve(scaled_model, tee=False)
    scaling.propagate_solution(scaled_model, model3)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues3(model3):
    dt = DiagnosticsToolbox(model3)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution3(model3):
    time_duration = 24
    conversion_comp = {
        "Al2O3": 0.040760147592520515,
        "CaO": 0.5239553356985602,
        "Ce2O3": 0.4381391585076152,
        "Dy2O3": 0.21831467846123956,
        "Fe2O3": 0.1754365137602976,
        "Gd2O3": 0.7013644011136486,
        "La2O3": 0.3876511363515543,
        "Nd2O3": 0.4879977553351574,
        "Pr2O3": 0.503579628144801,
        "Sc2O3": 0.05582841670114907,
        "Sm2O3": 0.23427259102997697,
        "Y2O3": 0.14909005595098204,
        "inerts": 0.0,
    }
    recovery = {
        "inerts": -8.285461265002091,
        "Sc2O3": -2.2400554108932935,
        "Y2O3": 7.858824213957128,
        "La2O3": 33.69152084600201,
        "Ce2O3": 39.15863961302317,
        "Pr2O3": 46.24489105368862,
        "Nd2O3": 44.557600768029815,
        "Sm2O3": 17.082854316882045,
        "Gd2O3": 67.66210642262544,
        "Dy2O3": 15.35484439386654,
        "Al2O3": -3.871729882517183,
        "CaO": 48.45128393345267,
        "Fe2O3": 10.711762549599014,
    }

    for k, v in model3.fs.leach.mscontactor.solid[time_duration, 2].conversion_comp.items():
        f_in = model3.fs.leach.solid_inlet.flow_mass[time_duration]
        f_out = model3.fs.leach.solid_outlet.flow_mass[time_duration]
        x_in = model3.fs.leach.solid_inlet.mass_frac_comp[time_duration, k]
        x_out = model3.fs.leach.solid_outlet.mass_frac_comp[time_duration, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        assert value(v) == pytest.approx(conversion_comp[k], rel=1e-5, abs=1e-6)
        assert r == pytest.approx(recovery[k], rel=1e-5, abs=1e-6)
