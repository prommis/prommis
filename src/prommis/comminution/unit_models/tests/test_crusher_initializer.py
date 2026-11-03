#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Test crusher initialization, state restoration, and parameter estimation."""

import logging
import math
from unittest.mock import patch

import pytest

from pyomo.common.collections import ComponentMap
from pyomo.contrib.solver.common.util import NoSolutionError
from pyomo.core.base.var import VarData
from pyomo.environ import (
    Block,
    Constraint,
    Var,
    value,
)
from pyomo.opt import SolverResults, SolverStatus, TerminationCondition

from idaes.core.initialization import InitializationStatus
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom

from prommis.comminution.functions import distributions
from prommis.comminution.functions.power_laws import numeric_specific_energy
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models._crusher import initializer as crusher_initializer
from prommis.comminution.unit_models.tests.crusher_test_support import (
    DENSITY_TWO,
    FEED_REF,
    PATH_CONFIGS,
    RR,
    SINGLE_STRESS_REF,
    crusher,
    feed_percentile,
    fix_feed,
    flowsheet,
    one_mineral_flowsheet,
    two_time_flowsheet,
)

# Use an 8 mm target whose characteristic size lies inside the mesh.
INTERIOR_P80_REF = 8.0e-3


class _NonOptimalSolver:
    """A solver stub whose result reports an infeasible termination."""

    def __init__(self, *args, **kwargs):
        self.options = {}

    def solve(self, model, **kwargs):
        res = SolverResults()
        res.solver.status = SolverStatus.warning
        res.solver.termination_condition = TerminationCondition.infeasible
        return res


class _RaisingSolver(_NonOptimalSolver):
    """A solver stub whose solve call itself raises."""

    def solve(self, model, **kwargs):
        raise RuntimeError("planted")


class _NoSolutionSolver(_NonOptimalSolver):
    """A solver stub whose solve call returns no solution."""

    def solve(self, model, **kwargs):
        raise NoSolutionError("planted")


def _solver_setting_value(target, value_to_set):
    """Return a solver stub class that sets ``target`` to ``value_to_set``.

    The stub reports optimal termination without solving the model.
    """

    class _ValueSettingSolver(_NonOptimalSolver):
        def solve(self, model, **kwargs):
            target.set_value(value_to_set, skip_validation=True)
            res = SolverResults()
            res.solver.status = SolverStatus.ok
            res.solver.termination_condition = TerminationCondition.optimal
            return res

    return _ValueSettingSolver


def _flags(block):
    """Return Var fixed flags and Constraint active flags."""
    fixed = {v.name: v.fixed for v in block.component_data_objects(Var)}
    active = {
        c.name: c.active for c in block.component_data_objects(Constraint, active=None)
    }
    return fixed, active


@pytest.mark.unit
def test_legacy_initialize_is_refused_before_any_side_effect():
    m, cr = _ok_unit()
    # Register a block to check that the rejected legacy call leaves it active.
    cr.registered = Block()
    cr._initialization_order.append(cr.registered)
    with pytest.raises(NotImplementedError, match="CrusherSolidPSDInitializer"):
        cr.initialize()
    with pytest.raises(NotImplementedError, match="CrusherSolidPSDInitializer"):
        cr.initialize(None, outlvl=logging.DEBUG)
    assert cr.registered.active


# Required summary fields for a successful initialization.
OK_SUMMARY = {"status": InitializationStatus.Ok, "DoF": 0, "solver_status": True}


def _ok_unit():
    """Return a single-stress crusher with the reference feed fixed."""
    m = flowsheet()
    cr = crusher(m, bond_work_index=14.0, **SINGLE_STRESS_REF)
    fix_feed(cr.properties_in[0], FEED_REF)
    return m, cr


@pytest.mark.unit
def test_initializer_raises_on_nonoptimal_solve(monkeypatch):
    m = one_mineral_flowsheet()
    cr = crusher(m, psd_method="tabular", tabular_psd=[0.25, 0.25, 0.25, 0.25])
    fix_feed(cr.properties_in[0], {"OreA": [0.0, 0.0, 0.0, 10.0]})
    # A nonoptimal solve must restore this inlet variable's unfixed status.
    cr.properties_in[0].temperature.unfix()
    entry_flags = _flags(cr)
    entry_fixed = {v.name: v.value for v in cr.component_data_objects(Var) if v.fixed}
    init = native.CrusherSolidPSDInitializer()
    assert {
        "tol": 1e-8,
        "bound_push": 1e-8,
    }.items() <= init.config.solver_options.value().items()
    for stub, match in (
        (_NonOptimalSolver, "did not return optimal termination"),
        (_NoSolutionSolver, "returned no solution"),
    ):
        monkeypatch.setattr(crusher_initializer, "get_solver", stub)
        with pytest.raises(InitializationError, match=match):
            init.initialize(cr)
        assert {
            "status": InitializationStatus.Failed,
            "DoF": 0,
            "solver_status": False,
        }.items() <= init.summary[cr].items()
        assert {
            v.name: v.value for v in cr.component_data_objects(Var) if v.fixed
        } == entry_fixed
        assert _flags(cr) == entry_flags
    # A reused initializer must not keep the failed run's solver status.
    monkeypatch.setattr(crusher_initializer, "get_solver", _RaisingSolver)
    with pytest.raises(RuntimeError, match="planted"):
        init.initialize(cr)
    assert {"status": InitializationStatus.Error, "DoF": 0}.items() <= init.summary[
        cr
    ].items()
    assert "solver_status" not in init.summary[cr]
    # Optimal termination still requires finite variables and acceptable residuals.
    out = cr.properties_out[0].flow_mass_sized_comp_size
    for unit, target, value_to_set, match in [
        (cr, out["OreA", 1], math.nan, r"\(1\): \S*comp_size\[OreA,1\]\.$"),
        (cr, cr.power[0], math.nan, r"\(1\): fs\.cr\.power\[0\.0\]\.$"),
        (cr, out["OreA", 1], None, "could not be evaluated"),
        # This small flow change must fail the residual check (constr_viol_tol=1e-8).
        (
            cr,
            out["OreA", 1],
            2.5 + 1e-6,
            r"constr_viol_tol = 1e-08 \(\d+, largest .*product_psd_eqn",
        ),
    ]:
        monkeypatch.setattr(
            crusher_initializer,
            "get_solver",
            _solver_setting_value(target, value_to_set),
        )
        with pytest.raises(InitializationError, match=match):
            init.initialize(unit)
        assert init.summary[unit]["status"] == InitializationStatus.Error


@pytest.mark.component
@pytest.mark.solver
def test_no_entry_fixed_var_changes_value_during_initialization():
    m = flowsheet(
        unsized=["Inert"], density=dict(DENSITY_TWO, Inert=2000.0), vapors=["Air"]
    )
    cr = crusher(m, bond_work_index=14.0, **PATH_CONFIGS["classification_recycle"])
    feed = cr.properties_in[0]
    fix_feed(feed, FEED_REF, liquid=2.0, unsized={"Inert": 0.5})
    feed.flow_mass_vapor_comp["Air"].fix(0.25)
    # Fix measured power at the value a solved twin with Wi = 20 computes;
    # initialization must leave it unchanged.
    m_twin = flowsheet(
        unsized=["Inert"], density=dict(DENSITY_TWO, Inert=2000.0), vapors=["Air"]
    )
    twin = crusher(
        m_twin, bond_work_index=20.0, **PATH_CONFIGS["classification_recycle"]
    )
    fix_feed(twin.properties_in[0], FEED_REF, liquid=2.0, unsized={"Inert": 0.5})
    twin.properties_in[0].flow_mass_vapor_comp["Air"].fix(0.25)
    native.CrusherSolidPSDInitializer().initialize(twin)
    datum = value(twin.power[0])
    cr.bond_work_index.unfix()
    cr.bond_work_index.set_value(20.0)
    cr.power[0].fix(datum)
    entry_fixed = ComponentMap(
        (v, v.value) for v in cr.component_data_objects(Var) if v.fixed
    )
    assert cr.power[0] in entry_fixed
    real_set_value = VarData.set_value

    def _guarded(var, val, skip_validation=False):
        if var in entry_fixed and val != entry_fixed[var]:
            raise AssertionError(
                f"attempted to set entry-fixed Var {var.name} "
                f"from {entry_fixed[var]!r} to {val!r}"
            )
        return real_set_value(var, val, skip_validation=skip_validation)

    init = native.CrusherSolidPSDInitializer()
    with patch.object(VarData, "set_value", _guarded):
        init.initialize(cr)
    assert init.summary[cr]["status"] == InitializationStatus.Ok
    for var, entry_value in entry_fixed.items():
        assert var.fixed and var.value == entry_value
    assert cr.power[0].fixed and value(cr.power[0]) == datum
    assert value(cr.bond_work_index) == pytest.approx(20.0, rel=1e-8)
    out = cr.properties_out[0]
    assert value(out.flow_mass_unsized_comp["Inert"]) == pytest.approx(0.5, rel=1e-8)
    assert value(out.flow_mass_liquid_comp["H2O"]) == pytest.approx(2.0, rel=1e-8)
    assert value(out.flow_mass_vapor_comp["Air"]) == pytest.approx(0.25, rel=1e-8)
    # Writing the calculated power value would evade the guard. The second
    # unit fixes outlet temperature away from the inlet to expose a write,
    # then deactivates its pass-through equation to keep zero DoF.
    m2 = flowsheet(
        unsized=["Inert"], density=dict(DENSITY_TWO, Inert=2000.0), vapors=["Air"]
    )
    cr2 = crusher(m2, bond_work_index=14.0, **PATH_CONFIGS["classification_recycle"])
    feed2 = cr2.properties_in[0]
    fix_feed(feed2, FEED_REF, liquid=2.0, unsized={"Inert": 0.5})
    feed2.flow_mass_vapor_comp["Air"].fix(0.25)
    outlet_temperature = cr2.properties_out[0].temperature
    outlet_temperature.fix(310.0)
    cr2.temperature_eqn[0].deactivate()
    entry_fixed.update((v, v.value) for v in cr2.component_data_objects(Var) if v.fixed)
    assert outlet_temperature in entry_fixed
    init2 = native.CrusherSolidPSDInitializer()
    with patch.object(VarData, "set_value", _guarded):
        init2.initialize(cr2)
    assert init2.summary[cr2]["status"] == InitializationStatus.Ok
    assert outlet_temperature.fixed and value(outlet_temperature) == 310.0


# Test recovery of a work index near its lower bound from two initial guesses.
_WI_LO, _WI_HI = native.WORK_INDEX_BOUNDS
IMPLIED_WORK_INDEX = _WI_LO + 0.001


def _estimation_twin(wi_entry):
    """Return a model, crusher, and fixed size measurement for two time points.

    ``IMPLIED_WORK_INDEX`` sets measured power; ``wi_entry`` sets the initial work
    index.
    """
    m = two_time_flowsheet(species=["OreA", "OreB"], density=DENSITY_TWO)
    cfg = dict(PATH_CONFIGS["distribution_power_spec"], distribution_shape=RR)
    cr = crusher(m, bond_work_index=14.0, **cfg)
    for t in (0.0, 1.0):
        fix_feed(cr.properties_in[t], FEED_REF)
    # Compute the measured duty directly; power has not been fixed yet.
    measurement = cr.d_63[0.0]
    measurement.fix(distributions.shape_size_from_p80(RR, INTERIOR_P80_REF, 1.5))
    duty = (
        3.6
        * 10.0
        * numeric_specific_energy(
            "bond", IMPLIED_WORK_INDEX, INTERIOR_P80_REF, feed_percentile(cr)
        )
    )
    for t in (0.0, 1.0):
        cr.power[t].fix(duty)
    cr.bond_work_index.unfix()
    cr.bond_work_index.set_value(wi_entry)
    return m, cr, measurement


@pytest.mark.component
@pytest.mark.solver
def test_two_time_studies_recover_the_planted_work_index():
    twins = {}
    for wi_entry in (14.0, 60.0):
        m, cr, measurement = _estimation_twin(wi_entry)
        entry_fixed = {
            v.name: v.value for v in cr.component_data_objects(Var) if v.fixed
        }
        measured_value = value(measurement)
        init = native.CrusherSolidPSDInitializer()
        init.initialize(cr)
        twins[wi_entry] = {
            "m": m,
            "cr": cr,
            "measurement": measurement,
            "measured_value": measured_value,
            "entry_fixed": entry_fixed,
            "init": init,
        }
    for twin in twins.values():
        cr, init = twin["cr"], twin["init"]
        assert {
            v.name: v.value for v in cr.component_data_objects(Var) if v.fixed
        } == twin["entry_fixed"]
        assert OK_SUMMARY.items() <= init.summary[cr].items()
        assert value(cr.bond_work_index) == pytest.approx(IMPLIED_WORK_INDEX, rel=1e-6)
    cr_a, cr_b = twins[14.0]["cr"], twins[60.0]["cr"]
    got = value(cr_b.d_63[1.0])
    exp = value(cr_a.d_63[1.0])
    assert got == pytest.approx(exp, rel=1e-8)
    # Each twin's solved size at t=1 matches the planted measurement at t=0.
    for twin in twins.values():
        assert value(twin["cr"].d_63[1.0]) == pytest.approx(
            twin["measured_value"], rel=1e-8
        )
    pp = cr_a.config.property_package
    for s in pp.sized_solid_list:
        for k in pp.size_interval_set:
            row_a = value(cr_a.properties_out[1.0].flow_mass_sized_comp_size[s, k])
            row_b = value(cr_b.properties_out[1.0].flow_mass_sized_comp_size[s, k])
            assert row_b == pytest.approx(row_a, rel=1e-8, abs=1e-12)
    # Reinitialize after doubling the feed and power, then check the outlet bin total.
    out = cr_a.properties_out[1.0]
    before = value(out.flow_mass_size[1])
    for component in cr_a.properties_in[1.0].define_state_vars().values():
        if component.local_name.startswith("flow_mass"):
            for var in component.values():
                var.fix(2.0 * value(var))
    cr_a.power[1.0].fix(2.0 * value(cr_a.power[1.0]))
    assert degrees_of_freedom(cr_a) == 0
    twins[14.0]["init"].initialize(cr_a)
    assert value(out.flow_mass_size[1]) == pytest.approx(2.0 * before, rel=1e-8)


@pytest.mark.component
@pytest.mark.solver
def test_target_product_p80_changed_after_build_is_used_and_rechecked():
    m = flowsheet()
    cr = crusher(m, **PATH_CONFIGS["distribution_size_spec"])
    fix_feed(cr.properties_in[0], FEED_REF)
    cr.target_product_p80 = 4e-3
    native.CrusherSolidPSDInitializer().initialize(cr)
    assert value(cr.product_p80[0]) == pytest.approx(4e-3, rel=1e-8)
    cr.target_product_p80 = 20e-3
    init = native.CrusherSolidPSDInitializer()
    with pytest.raises(ConfigurationError, match="attainable mesh band"):
        init.initialize(cr)
    assert init.summary[cr]["status"] == InitializationStatus.PrecheckFailed
