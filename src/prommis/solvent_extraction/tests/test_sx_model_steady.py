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

from prommis.solvent_extraction.sx_model_flowsheet_steady_state import (
main,
)


@pytest.fixture(scope="module")
def model():
    m = main()

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
    number_of_stages = 3
    aqueous_outlet = {
        "H2O": 1000000,
        "H": 1.75563,
        "SO4": 3999.885,
        "HSO4": 693.3903,
        "Al": 362.132,
        "Ca": 81.724,
        "Ce": 0.5368,
        "Dy": 0.0014421,
        "Fe": 454.049,
        "Gd": 0.066625,
        "La": 0.39313,
        "Nd": 1.173e-07,
        "Pr": 0.091292,
        "Sc": 0.00019984,
        "Sm": 9.6992e-11,
        "Y": 1.239e-10,
    }

    organic_outlet = {
        "Al": 60.242,
        "Ca": 27.817,
        "Ce": 1.7405,
        "Dy": 0.045565,
        "Fe": 234.216,
        "Gd": 0.1918,
        "La": 0.59296,
        "Nd": 0.9461,
        "Pr": 0.211744,
        "Sc": 1.7658,
        "Sm": 0.097017,
        "Y": 0.12402,
    }

    for k, v in model.fs.solex.mscontactor.organic[0, 1].conc_mass_comp.items():
        assert value(v) == pytest.approx(organic_outlet[k], rel=1e-4)

    for k, v in model.fs.solex.mscontactor.aqueous[
        0, number_of_stages
    ].conc_mass_comp.items():
        assert value(v) == pytest.approx(aqueous_outlet[k], rel=1e-4)
