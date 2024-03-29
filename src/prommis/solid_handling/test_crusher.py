import pytest
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

# Assuming these imports are adjusted to your project's structure
from prommis.solid_handling.crusher import CrushAndBreakageUnit
from prommis.leaching.leach_solids_properties import CoalRefuseParameters

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_solid = CoalRefuseParameters(
        doc="solid property",
    )

    m.fs.unit = CrushAndBreakageUnit(property_package=m.fs.properties_solid)

    # Assert specific config options as per your model's requirements
    # Example:
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.property_package is m.fs.properties_solid


# -----------------------------------------------------------------------------
class TestSolidHandling(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties_solid = CoalRefuseParameters(
            doc="solid property",
        )

        m.fs.unit = CrushAndBreakageUnit(property_package=m.fs.properties_solid)
        # Set up your model initialization here
        m.fs.unit.inlet.mass_frac_comp[0, :].fix(0.1)  # set mass frac value. 13
        m.fs.unit.control_volume.properties_in[0].flow_mass.fix(2000)  # kg/hr
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        model.display()
        # Assertions related to the build status of the model

        # More assertions as needed for your model

        assert number_variables(model.fs.unit) == 57
        assert number_total_constraints(model.fs.unit) == 43
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, model):
        # assert_units_consistent(model.fs.unit)

        dt = DiagnosticsToolbox(model=model)
        dt.report_structural_issues()
        dt.display_underconstrained_set()
        dt.display_components_with_inconsistent_units()
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)
        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

    def test_solve(self, model):
        results = solver.solve(model)
        assert_optimal_termination(results)

    def test_solution(self, model):
        assert pytest.approx(0.8, abs=1e-0) == value(model.fs.unit.probfeed80[0])
        assert pytest.approx(0.8, abs=1e-0) == value(model.fs.unit.probfeed80[0])
        assert pytest.approx(3.95, abs=1e-0) == value(model.fs.unit.crushpower[0])
