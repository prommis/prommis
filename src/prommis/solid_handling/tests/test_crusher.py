#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Var,
    assert_optimal_termination,
    value,
)

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

import pytest

from prommis.solid_handling.crusher import Crusher
from prommis.solid_handling.crusher_solids_properties import CoalRefuseParameters

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

    m.fs.unit = Crusher(property_package=m.fs.properties_solid)

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

        m.fs.unit = Crusher(property_package=m.fs.properties_solid)
        # Set up your model initialization here
        m.fs.unit.inlet.mass_frac_comp[0, :].fix(
            0.1
        )  # set mass frac value. There are 13 component in property package
        m.fs.unit.properties_in[0].flow_mass.fix(2000)  # kg/hr
        m.fs.unit.properties_in[0].particle_size_median.fix(80)  # micrometer
        m.fs.unit.properties_in[0].particle_size_width.fix(1.5)  # dimensionless
        m.fs.unit.properties_out[0].particle_size_median.fix(58)  # micrometer
        m.fs.unit.properties_out[0].particle_size_width.fix(1.5)  # dimensionless
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        # More assertions as needed for your model
        unit = model.fs.unit
        assert isinstance(
            unit.work, Var
        ), "Variable 'work' is not properly initialized as a Pyomo Var."
        assert isinstance(
            unit.feed_p80, Expression
        ), "Expression 'feed_p80' is not properly initialized as a Pyomo Expression."
        assert isinstance(
            unit.prod_p80, Expression
        ), "Expression 'prod_p80' is not properly initialized as a Pyomo Expression."
        assert isinstance(
            unit.crush_work_eq, Constraint
        ), "Constraint 'crush_work_eq' is not properly initialized as a Pyomo Constraint."

        assert number_variables(model.fs.unit) == 59
        assert number_total_constraints(model.fs.unit) == 41
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)
        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        assert pytest.approx(
            value(model.fs.unit.properties_in[0].flow_mass), rel=1e-5
        ) == value(model.fs.unit.properties_out[0].flow_mass)
        for j in model.fs.properties_solid.component_list:
            assert pytest.approx(
                value(model.fs.unit.properties_in[0].mass_frac_comp[j]), rel=1e-5
            ) == value(model.fs.unit.properties_out[0].mass_frac_comp[j])

        assert pytest.approx(114.31301, rel=1e-5) == value(
            model.fs.unit.feed_p80[0]
        )  # Test feed size expressions.
        assert pytest.approx(82.87693, rel=1e-5) == value(
            model.fs.unit.prod_p80[0]
        )  # Test prod size expressions.

        assert pytest.approx(3915.710575, rel=1e-5) == value(model.fs.unit.work[0])
