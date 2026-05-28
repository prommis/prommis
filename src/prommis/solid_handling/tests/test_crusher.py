#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import math

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Var,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.util.check_units import assert_units_consistent


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
    degrees_of_freedom,
)

import pytest

from prommis.solid_handling.crusher import Crusher
from prommis.solid_handling.crusher_solids_properties import CoalRefuseParameters

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def _build_model(**crusher_kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_solid = CoalRefuseParameters(doc="solid property")
    m.fs.unit = Crusher(property_package=m.fs.properties_solid, **crusher_kwargs)
    return m


def _set_base_state(unit, feed_median_um=8e4, prod_median_um=5.8e4):
    unit.inlet.mass_frac_comp[0, :].fix(0.1)
    unit.properties_in[0].flow_mass.fix(2000)
    unit.properties_in[0].particle_size_median.fix(feed_median_um)
    unit.properties_in[0].particle_size_width.fix(1.5)
    unit.properties_out[0].particle_size_median.fix(prod_median_um)
    unit.properties_out[0].particle_size_width.fix(1.5)


def _fix_product_p80(unit, product_p80_target):
    width = value(unit.properties_out[0].particle_size_width)
    p80_factor = (-math.log(1 - 0.8)) ** (width / 2)
    product_p80_um = value(
        pyunits.convert(product_p80_target, to_units=pyunits.um) / pyunits.um
    )
    unit.properties_out[0].particle_size_median.fix(product_p80_um / p80_factor)


def _cm_value(native_size, unit):
    return value(
        pyunits.convert(
            native_size * unit.properties_in[0].particle_size_median.get_units(),
            to_units=pyunits.cm,
        )
        / pyunits.cm
    )


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = _build_model()

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.property_package is m.fs.properties_solid
    assert m.fs.unit.config.crusher_stage is None
    assert m.fs.unit.config.crusher_equipment is None
    assert not m.fs.unit.config.enforce_stage_size_limits
    assert m.fs.unit.crusher_stage == "secondary"
    assert m.fs.unit.stage_selection_basis == "default"


@pytest.mark.unit
@pytest.mark.parametrize(
    "stage_name, expected_lower_cm, expected_upper_cm",
    [
        ("primary", 10.0, None),
        ("secondary", 2.0, 10.0),
        ("tertiary", 0.5, 2.0),
    ],
)
def test_stage_selection(stage_name, expected_lower_cm, expected_upper_cm):
    m = _build_model(crusher_stage=stage_name)

    assert m.fs.unit.crusher_stage == stage_name
    assert m.fs.unit.crusher_equipment is None
    assert m.fs.unit.stage_selection_basis == "stage"
    assert pytest.approx(expected_lower_cm, rel=1e-8) == _cm_value(
        m.fs.unit._applicable_product_p80_lower, m.fs.unit
    )
    if expected_upper_cm is None:
        assert m.fs.unit._applicable_product_p80_upper is None
    else:
        assert pytest.approx(expected_upper_cm, rel=1e-8) == _cm_value(
            m.fs.unit._applicable_product_p80_upper, m.fs.unit
        )


@pytest.mark.unit
@pytest.mark.parametrize(
    "equipment_name, expected_stage",
    [
        ("jaw", "primary"),
        ("gyratory", "primary"),
        ("cone", "secondary"),
        ("roll1", "secondary"),
        ("roll2", "tertiary"),
        ("short_head_cone", "tertiary"),
        ("hammer_mill", "tertiary"),
    ],
)
def test_equipment_selection_without_stage(equipment_name, expected_stage):
    m = _build_model(crusher_equipment=equipment_name)

    assert m.fs.unit.config.crusher_stage is None
    assert m.fs.unit.crusher_equipment == equipment_name
    assert m.fs.unit.crusher_stage == expected_stage
    assert m.fs.unit.stage_selection_basis == "equipment"


@pytest.mark.unit
def test_consistent_stage_and_equipment_allowed():
    m = _build_model(crusher_stage="primary", crusher_equipment="jaw")

    assert m.fs.unit.crusher_stage == "primary"
    assert m.fs.unit.crusher_equipment == "jaw"
    assert m.fs.unit.stage_selection_basis == "stage_and_equipment"


@pytest.mark.unit
def test_inconsistent_stage_and_equipment_raises():
    with pytest.raises(ValueError, match="Inconsistent crusher specification"):
        _build_model(crusher_stage="tertiary", crusher_equipment="jaw")


@pytest.mark.unit
@pytest.mark.parametrize(
    "product_p80_target, expected_stage, in_modeled_range",
    [
        (12 * pyunits.cm, "primary", True),
        (10 * pyunits.cm, "primary", True),
        (5 * pyunits.cm, "secondary", True),
        (2 * pyunits.cm, "secondary", True),
        (1 * pyunits.cm, "tertiary", True),
        (0.3 * pyunits.cm, "tertiary", False),
    ],
)
def test_recommend_stage_for_product_p80(
    product_p80_target, expected_stage, in_modeled_range
):
    m = _build_model(enforce_stage_size_limits=False)
    _set_base_state(m.fs.unit, feed_median_um=200000, prod_median_um=108000)
    _fix_product_p80(m.fs.unit, product_p80_target)

    recommended_stage, modeled_range = m.fs.unit.recommend_stage_for_product_p80()

    assert recommended_stage == expected_stage
    assert modeled_range is in_modeled_range


# -----------------------------------------------------------------------------
class TestSolidHandling(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = _build_model(crusher_stage="secondary", enforce_stage_size_limits=False)
        _set_base_state(m.fs.unit)

        # No unfixed variable exist for crushing_direction_constraint,
        # deactivate it to avoid assert_no_numerical_warning failure
        m.fs.unit.crushing_direction_constraint[0].deactivate()
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        unit = model.fs.unit
        assert isinstance(unit.work, Var)
        assert isinstance(unit.feed_p80, Expression)
        assert isinstance(unit.prod_p80, Expression)
        assert isinstance(unit.crush_work_eq, Constraint)

        assert number_variables(model.fs.unit) == 59
        assert number_total_constraints(model.fs.unit) == 42
        assert number_unused_variables(model.fs.unit) == 0

    
    @pytest.mark.component
    def test_structural_issues(self, model):
        assert_units_consistent(model.fs.unit)
        assert degrees_of_freedom(model.fs.unit) == 0
        
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
        # The assert_no_numerical_warnings will fail if there is an active constraint with no Jacobian
        # entries and that all variables in that constraint (e.g. crushing_direction_constraint) are fixed
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

        assert pytest.approx(114313.009, rel=1e-5) == value(
            model.fs.unit.feed_p80[0]
        )  # Test feed size expressions.
        assert pytest.approx(82876.93, rel=1e-5) == value(
            model.fs.unit.prod_p80[0]
        )  # Test prod size expressions.

        assert pytest.approx(123.826, rel=1e-5) == value(model.fs.unit.work[0])


@pytest.mark.unit
def test_check_applicability_no_warning_for_in_range_secondary(caplog):
    m = _build_model(crusher_stage="secondary", enforce_stage_size_limits=True)
    _set_base_state(m.fs.unit, feed_median_um=150000, prod_median_um=58)
    _fix_product_p80(m.fs.unit, 5 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert warning_messages == []
    assert caplog.text == ""


@pytest.mark.unit
def test_check_applicability_warns_when_secondary_target_belongs_to_primary(caplog):
    m = _build_model(crusher_stage="secondary", enforce_stage_size_limits=False)
    _set_base_state(m.fs.unit, feed_median_um=200000, prod_median_um=108000)
    _fix_product_p80(m.fs.unit, 12 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert len(warning_messages) == 1
    assert "falls in the 'primary' crushing range" in warning_messages[0]
    assert "[2, 10) cm" in warning_messages[0]
    assert "falls in the 'primary' crushing range" in caplog.text


@pytest.mark.unit
def test_check_applicability_warns_when_primary_target_belongs_to_secondary(caplog):
    m = _build_model(crusher_stage="primary", enforce_stage_size_limits=False)
    _set_base_state(m.fs.unit, feed_median_um=150000, prod_median_um=108000)
    _fix_product_p80(m.fs.unit, 9.999 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert len(warning_messages) == 1
    assert "falls in the 'secondary' crushing range" in warning_messages[0]
    assert ">= 10 cm" in warning_messages[0]
    assert "falls in the 'secondary' crushing range" in caplog.text


@pytest.mark.unit
def test_check_applicability_warning_mentions_equipment(caplog):
    m = _build_model(crusher_equipment="hammer_mill")
    _set_base_state(m.fs.unit, feed_median_um=50000)
    _fix_product_p80(m.fs.unit, 1 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert warning_messages == []
    assert "equipment='hammer_mill'" not in caplog.text


@pytest.mark.unit
def test_check_applicability_warning_mentions_equipment_when_out_of_range(caplog):
    m = _build_model(crusher_equipment="hammer_mill", enforce_stage_size_limits=False)
    _set_base_state(m.fs.unit, feed_median_um=150000, prod_median_um=58000)
    _fix_product_p80(m.fs.unit, 5 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert len(warning_messages) == 1
    assert "equipment='hammer_mill'" in warning_messages[0]
    assert "falls in the 'secondary' crushing range" in warning_messages[0]
    assert "equipment='hammer_mill'" in caplog.text


@pytest.mark.unit
def test_check_applicability_warns_when_target_is_finer_than_modeled_range(caplog):
    m = _build_model(crusher_stage="tertiary", enforce_stage_size_limits=False)
    _set_base_state(m.fs.unit, feed_median_um=50000, prod_median_um=48000)
    _fix_product_p80(m.fs.unit, 0.3 * pyunits.cm)

    with caplog.at_level("WARNING"):
        warning_messages = m.fs.unit.check_applicability()

    assert len(warning_messages) == 1
    assert (
        "finer than the minimum modeled crushing range (0.5 cm)" in warning_messages[0]
    )
    assert "finer than the minimum modeled crushing range (0.5 cm)" in caplog.text

@pytest.mark.component
@pytest.mark.solver
def test_solve_with_stage_bounds_and_crushing_direction_active():
    m = _build_model(crusher_stage="secondary", enforce_stage_size_limits=True)

    _set_base_state(m.fs.unit, feed_median_um=150000, prod_median_um=58000)

    # Let the solver determine the outlet median particle size.
    m.fs.unit.properties_out[0].particle_size_median.unfix()

    # Close the model by specifying the desired product P80.
    m.fs.unit.product_p80_spec = Constraint(
        expr=m.fs.unit.prod_p80[0]
        == value(pyunits.convert(5 * pyunits.cm, to_units=pyunits.um) / pyunits.um)
    )

    # Following constraints should remain active.
    assert m.fs.unit.crushing_direction_constraint[0].active
    assert m.fs.unit.product_p80_lower_bound[0].active
    assert m.fs.unit.product_p80_upper_bound[0].active
    
    assert_units_consistent(m.fs.unit)
    assert degrees_of_freedom(m.fs.unit) == 0
        
    results = solver.solve(m)
    assert_optimal_termination(results)

    assert pytest.approx(5.0, rel=1e-6) == _cm_value(
        value(m.fs.unit.prod_p80[0]), m.fs.unit
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()
    dt.assert_no_numerical_warnings()