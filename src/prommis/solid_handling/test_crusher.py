import pytest
from pyomo.environ import ConcreteModel, assert_optimal_termination
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
        # Example: m.fs.unit.some_inlet_variable.fix(some_value)
        m.fs.unit.inlet.fow_mass.fix(2)  # tonne/hr
        m.fs.unit.inlet.feed50size.fix(60)  # eed particle size, micrometer
        m.fs.unit.inlet.Prod50size.fix(
            30
        )  # Product Particle Size that allows 80% passing, micrometer

        # m.fs.unit.soliddistribution.fix(0.48)
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        # Assertions related to the build status of the model

        # More assertions as needed for your model

        assert number_variables(model.fs.unit) == 3
        assert number_total_constraints(model.fs.unit) == 3
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model.fs.unit)

        dt = DiagnosticsToolbox(model=model)
        dt.report_structural_issues()
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)
        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)
        assert_optimal_termination(results)
