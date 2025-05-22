#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    SolverFactory,
    TransformationFactory,
    assert_optimal_termination,
    value,
)

from idaes.core.util import DiagnosticsToolbox

import pytest

from prommis.leaching.leach_flowsheet import build_model, set_inputs, set_scaling


@pytest.fixture(scope="module")
def model():
    m = build_model()
    set_inputs(m)

    return m


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solve(model):
    set_scaling(model)
    # Create a scaled version of the model to solve
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model, rename=False)

    initializer = model.fs.leach.default_initializer()
    try:
        initializer.initialize(scaled_model.fs.leach)
    except:
        pass

    # Solve scaled model
    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=False)

    # Propagate results back to unscaled model
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
        "inerts": 3.3306690738754696e-14,
        "Sc2O3": 5.486783341761903,
        "Y2O3": 15.092538416882373,
        "La2O3": 38.342494088335066,
        "Ce2O3": 42.074359313539645,
        "Pr2O3": 48.42509098259142,
        "Nd2O3": 45.79540997369025,
        "Sm2O3": 23.377407833768427,
        "Gd2O3": 61.06358012227954,
        "Dy2O3": 22.365403062815247,
        "Al2O3": 3.0037655469309033,
        "CaO": 51.45052459223044,
        "Fe2O3": 16.33850709753919,
    }

    for k, v in model.fs.leach.mscontactor.solid[0, 1].conversion.items():
        f_in = model.fs.leach.solid_inlet.flow_mass[0]
        f_out = model.fs.leach.solid_outlet.flow_mass[0]
        x_in = model.fs.leach.solid_inlet.mass_frac_comp[0, k]
        x_out = model.fs.leach.solid_outlet.mass_frac_comp[0, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        assert value(v) == pytest.approx(conversion[k], rel=1e-5, abs=1e-6)
        assert r == pytest.approx(recovery[k], rel=1e-5, abs=1e-6)
