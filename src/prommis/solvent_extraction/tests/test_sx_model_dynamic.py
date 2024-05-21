from pyomo.environ import (
    SolverFactory,
    assert_optimal_termination,
    value,
)

from idaes.core.util import DiagnosticsToolbox

from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)

import pytest

from prommis.solvent_extraction.sx_model_flowsheet_dynamic import (
    build_model,
    set_inputs,
)


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
    initializer = BlockTriangularizationInitializer(constraint_tolerance=1e-4)
    initializer.initialize(model.fs.solex)

    try:
        initializer.initialize(model.fs.solex)
    except:
        pass

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
    number_of_stages = 3
    aqueous_outlet = {
        "H2O": 1000000,
        "H": 1.85700,
        "SO4": 3758.32573,
        "HSO4": 689.13997,
        "Al": 343.1403,
        "Ca": 77.43828,
        "Ce": 0.50865,
        "Dy": 0.0013665,
        "Fe": 430.23639,
        "Gd": 0.06313,
        "La": 0.37252,
        "Nd": 1.1133e-07,
        "Pr": 0.08650,
        "Sc": 0.000189363,
        "Sm": 2.70164e-10,
        "Y": 2.95748e-10,
    }

    organic_outlet = {
        "Al": 57.76735,
        "Ca": 26.28296,
        "Ce": 1.71612,
        "Dy": 0.04555,
        "Fe": 230.08492,
        "Gd": 0.19107,
        "La": 0.57677,
        "Nd": 0.94441,
        "Pr": 0.20955,
        "Sc": 1.75925,
        "Sm": 0.09701,
        "Y": 0.12401,
    }

    for k, v in model.fs.solex.mscontactor.organic[
        time_duration, 1
    ].conc_mass_comp.items():
        assert value(v) == pytest.approx(organic_outlet[k], rel=1e-4)

    for k, v in model.fs.solex.mscontactor.aqueous[
        time_duration, number_of_stages
    ].conc_mass_comp.items():
        assert value(v) == pytest.approx(aqueous_outlet[k], rel=1e-4)
