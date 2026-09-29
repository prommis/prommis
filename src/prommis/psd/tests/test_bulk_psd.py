#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the bulk PSD property package."""

import math

import pytest

from pyomo.environ import ConcreteModel, units, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.scaling import get_scaling_factor, set_scaling_factor

from prommis.psd.properties.bulk_psd import (
    BulkPSDInitializer,
    BulkPSDParameterBlock,
    BulkPSDScaler,
)

# Hand-calc mesh [1, 2, 4, 8, 16] mm in meters, 4 intervals.
HANDCALC_EDGES_M = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]


def _handcalc_params():
    return BulkPSDParameterBlock(
        size_edges=HANDCALC_EDGES_M, component_list=["Ore1", "Ore2", "Ore3"]
    )


def _fix_handcalc_feed(b):
    """Fix the hand-calc state: size flows [1,2,3,4] kg/s, comps summing to 10."""
    for k, val in zip([0, 1, 2, 3], [1.0, 2.0, 3.0, 4.0]):
        b.flow_mass_size[k].fix(val)
    b.flow_mass_comp["Ore1"].fix(5.0)
    b.flow_mass_comp["Ore2"].fix(3.0)
    b.flow_mass_comp["Ore3"].fix(2.0)


# -----------------------------------------------------------------------------
# Build / structure
@pytest.mark.build
@pytest.mark.unit
def test_build_structure():
    m = ConcreteModel()
    m.params = _handcalc_params()
    assert list(m.params.size_interval_set) == [0, 1, 2, 3]
    assert list(m.params.size_edge_index) == [0, 1, 2, 3, 4]
    assert list(m.params.component_list) == ["Ore1", "Ore2", "Ore3"]

    m.state = m.params.build_state_block([0], defined_state=True)
    b = m.state[0]
    assert len(b.flow_mass_comp) == 3
    assert len(b.flow_mass_size) == 4
    assert_units_consistent(b.flow_mass_comp["Ore1"])
    assert str(units.get_units(b.flow_mass_comp["Ore1"])) == str(units.kg / units.s)
    assert set(b.define_state_vars().keys()) == {"flow_mass_comp", "flow_mass_size"}
    # defined_state=True has no consistency constraint
    assert not hasattr(b, "flow_consistency_eqn")


@pytest.mark.unit
def test_param_block_rejects_bad_config():
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(size_edges=HANDCALC_EDGES_M, component_list=[])
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(
            size_edges=HANDCALC_EDGES_M, component_list=["Ore1", "Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(
            size_edges=HANDCALC_EDGES_M, component_list=["Ore1(s)"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(component_list=["Ore1"])  # missing edges
    # invalid size meshes surface as ConfigurationError (not raw core ValueError)
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # not strictly increasing
            size_edges=[1e-3, 1e-3, 2e-3], component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # < 2 edges
            size_edges=[1e-3], component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # non-numeric bottom_size
            size_edges=[0.0, 1e-3, 2e-3], bottom_size="x", component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # non-numeric mesh edge (float(None))
            size_edges=[1e-3, None], component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # boolean mesh edge (would be 1.0)
            size_edges=[True, 2e-3], component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError, match="boolean"):
        m = ConcreteModel()
        # coarse mesh where True -> 1.0 would pass the bottom_size < x_1 range rule,
        # so this exercises the bool type guard rather than the range check
        m.params = BulkPSDParameterBlock(  # boolean bottom_size (would be 1.0)
            size_edges=[0.0, 2.0, 4.0], bottom_size=True, component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # bare string splits into char components
            size_edges=HANDCALC_EDGES_M, component_list="ABC"
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # int too large to convert to float
            size_edges=[1, 10**10000], component_list=["Ore1"]
        )
    with pytest.raises(ConfigurationError):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(  # non-sequence size_edges (list(1.0) fails)
            size_edges=1.0, component_list=["Ore1"]
        )


@pytest.mark.unit
@pytest.mark.parametrize("reserved", ["solid", "size_edges", "d_char", "flow_eps"])
def test_param_block_rejects_reserved_component_names(reserved):
    # A valid identifier that collides with the solid phase name or an internal
    # parameter-block attribute must fail fast with a clear ConfigurationError
    # rather than silently replacing the internal object during build.
    with pytest.raises(ConfigurationError, match="reserved"):
        m = ConcreteModel()
        m.params = BulkPSDParameterBlock(
            size_edges=HANDCALC_EDGES_M, component_list=["Ore1", reserved]
        )


@pytest.mark.unit
def test_assert_same_mesh_wrong_type_raises_config_error():
    m = ConcreteModel()
    m.params = _handcalc_params()
    # a wrong-type argument yields a clear ConfigurationError, not AttributeError
    with pytest.raises(ConfigurationError, match="bulk-PSD parameter block"):
        m.params.assert_same_mesh("not a parameter block")


@pytest.mark.unit
def test_defined_state_false_has_consistency_constraint():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.state = m.params.build_state_block([0], defined_state=False)
    assert hasattr(m.state[0], "flow_consistency_eqn")
    assert m.state[0].flow_consistency_eqn.active


# -----------------------------------------------------------------------------
# Hand-calc state
@pytest.mark.unit
def test_handcalc_state():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.state = m.params.build_state_block([0], defined_state=True)
    b = m.state[0]
    _fix_handcalc_feed(b)

    assert value(b.flow_mass) == pytest.approx(10.0, rel=1e-12)
    assert value(b.flow_mass_sized) == pytest.approx(10.0, rel=1e-12)
    cum = [value(b.cum_passing[k]) for k in m.params.size_interval_set]
    assert cum == pytest.approx([0.1, 0.3, 0.6, 1.0], rel=1e-9)
    # P80 = 12.0 mm (the linear-in-size value; smooth form at default eps)
    assert value(units.convert(b.P80, units.mm)) == pytest.approx(12.0, rel=1e-4)
    assert value(units.convert(b.P50, units.mm)) == pytest.approx(20.0 / 3.0, rel=1e-4)
    # characteristic size of top bin (not a Pxx value)
    assert value(units.convert(m.params.d_char[3], units.mm)) == pytest.approx(
        math.sqrt(128.0), rel=1e-9
    )
    assert b.validate_feed()


# -----------------------------------------------------------------------------
# Sum-invariant and units
@pytest.mark.unit
def test_sum_invariant_and_units():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.state = m.params.build_state_block([0], defined_state=True)
    b = m.state[0]
    _fix_handcalc_feed(b)

    sum_fc = sum(value(b.mass_frac_comp[j]) for j in m.params.component_list)
    sum_fs = sum(value(b.mass_frac_size[k]) for k in m.params.size_interval_set)
    cum_top = value(b.cum_passing[m.params.size_interval_set.last()])
    assert abs(sum_fc - 1.0) <= 1e-12
    assert abs(sum_fs - 1.0) <= 1e-12
    assert abs(cum_top - 1.0) <= 1e-12
    # exercises every floored Expression for unit consistency
    assert_units_consistent(b)


# -----------------------------------------------------------------------------
# Expression safety on degenerate states
@pytest.mark.unit
def test_expression_safety_on_degenerate_states():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.state = m.params.build_state_block([0], defined_state=True)
    b = m.state[0]

    def _eval_all():
        vals = []
        vals.append(value(b.flow_mass))
        vals.append(value(b.flow_mass_sized))
        vals += [value(b.mass_frac_comp[j]) for j in m.params.component_list]
        vals += [value(b.mass_frac_size[k]) for k in m.params.size_interval_set]
        vals += [value(b.cum_passing[k]) for k in m.params.size_interval_set]
        vals.append(value(b.P80))
        vals.append(value(b.P50))
        return vals

    # default-initialized (small positive), zero, and near-zero states
    for state in ("default", "zero", "near_zero"):
        if state == "zero":
            for v in list(b.flow_mass_comp.values()) + list(b.flow_mass_size.values()):
                v.set_value(0.0)
        elif state == "near_zero":
            for v in list(b.flow_mass_comp.values()) + list(b.flow_mass_size.values()):
                v.set_value(1e-15)
        for got in _eval_all():
            assert math.isfinite(got)


# -----------------------------------------------------------------------------
# Feed validation edge cases
@pytest.mark.unit
def test_validate_feed_rejects_invalid_classes():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.state = m.params.build_state_block([0], defined_state=True)
    b = m.state[0]

    # baseline consistent feed
    _fix_handcalc_feed(b)
    assert b.validate_feed()

    # 1. negative flow
    b.flow_mass_size[1].fix(-1.0)
    with pytest.raises(ConfigurationError):
        b.validate_feed()
    b.flow_mass_size[1].fix(2.0)

    # 2. inconsistent totals (size 1% high)
    b.flow_mass_size[3].fix(4.1)  # size sum 10.1 vs comp 10
    with pytest.raises(ConfigurationError):
        b.validate_feed()
    b.flow_mass_size[3].fix(4.0)

    # 3. fractions summing to 0.98 scaled to flows -> consistency path
    for k, frac in zip([0, 1, 2, 3], [0.1, 0.2, 0.3, 0.38]):  # size sum 9.8
        b.flow_mass_size[k].fix(10.0 * frac)
    with pytest.raises(ConfigurationError):
        b.validate_feed()

    # 4. zero total
    for v in list(b.flow_mass_comp.values()) + list(b.flow_mass_size.values()):
        v.fix(0.0)
    with pytest.raises(ConfigurationError):
        b.validate_feed()

    # 5. non-finite
    _fix_handcalc_feed(b)
    b.flow_mass_comp["Ore1"].fix(float("inf"))
    with pytest.raises(ConfigurationError):
        b.validate_feed()


# -----------------------------------------------------------------------------
# Initialization lifecycle
@pytest.mark.component
def test_initialization_lifecycle_restores_constraint():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    b = m.out[0]
    # consistent values (as would be propagated from upstream); not fixed, so the
    # Initializer fixes them itself and postcheck sees the restored constraint
    # satisfied
    for j, val in zip(["Ore1", "Ore2", "Ore3"], [5.0, 3.0, 2.0]):
        b.flow_mass_comp[j].set_value(val)
    for k, val in zip([0, 1, 2, 3], [1.0, 2.0, 3.0, 4.0]):
        b.flow_mass_size[k].set_value(val)
    assert b.flow_consistency_eqn.active

    BulkPSDInitializer().initialize(m.out)
    # restored active and vars unfixed after init
    assert b.flow_consistency_eqn.active
    assert not b.flow_mass_comp["Ore1"].fixed


@pytest.mark.component
def test_fix_initialization_states_deactivates_consistency():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    b = m.out[0]
    m.out.fix_initialization_states()
    assert not b.flow_consistency_eqn.active
    assert b.flow_mass_comp["Ore1"].fixed


@pytest.mark.component
def test_initialization_does_not_reactivate_caller_deactivated():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    b = m.out[0]
    b.flow_consistency_eqn.deactivate()  # caller deactivated before init
    BulkPSDInitializer().initialize(m.out)
    assert not b.flow_consistency_eqn.active  # stays deactivated


@pytest.mark.component
def test_initialization_raise_path_restores_state():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.feed = m.params.build_state_block([0], defined_state=True)
    b = m.feed[0]
    # deliberately inconsistent fixed feed: comp 10, size 11
    b.flow_mass_comp["Ore1"].fix(5.0)
    b.flow_mass_comp["Ore2"].fix(3.0)
    b.flow_mass_comp["Ore3"].fix(2.0)
    for k, val in zip([0, 1, 2, 3], [1.0, 2.0, 3.0, 5.0]):
        b.flow_mass_size[k].fix(val)

    with pytest.raises(ConfigurationError):
        BulkPSDInitializer().initialize(m.feed)
    # fixedness restored (they were fixed before init)
    assert b.flow_mass_comp["Ore1"].fixed
    assert b.flow_mass_size[0].fixed


# -----------------------------------------------------------------------------
# Scaler zero-bin floor
@pytest.mark.unit
def test_scaler_zero_bin_floor():
    # tertiary product vector: 6 of 25 bins nonzero
    q = [0.0] * 8 + [0.42, 0.13, 0.21, 0.04, 0.14, 0.06] + [0.0] * 11
    assert len(q) == 25
    edges = [1e-6] + [(i + 1) * 1e-3 for i in range(25)]  # 26 edges
    m = ConcreteModel()
    m.params = BulkPSDParameterBlock(size_edges=edges, component_list=["Ore1"])
    m.s = m.params.build_state_block([0], defined_state=True)
    b = m.s[0]
    total = 100.0
    for k in m.params.size_interval_set:
        b.flow_mass_size[k].set_value(total * q[k])
    b.flow_mass_comp["Ore1"].set_value(total)

    # the zero-bin floor is a SEEDED-source contract (value reads); the
    # promoted baseline default would give unfixed set_value bins the coarse
    # 1.0 default instead
    scaler = BulkPSDScaler(factor_source="seeded")
    scaler.variable_scaling_routine(b)
    floor = scaler.EPS_REL * total

    factors = []
    for k in m.params.size_interval_set:
        v = b.flow_mass_size[k]
        expected = 1.0 / max(abs(value(v)), floor)
        got = scaler.get_scaling_factor(v)
        assert got == pytest.approx(expected, rel=1e-10)
        factors.append(got)
    # no 1/eps_abs blow-ups: the factor spread stays within the eps_rel band
    assert max(factors) / min(factors) <= 1.0 / scaler.EPS_REL + 1.0


@pytest.mark.unit
def test_scaler_consistency_row_requires_declared_factors():
    # PR-264 contract: the consistency-row factor derives from the declared
    # flow-Var factors; scaling constraints without them raises.
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    b = m.out[0]
    scaler = BulkPSDScaler()
    with pytest.raises(ConfigurationError, match="Missing scaling factor"):
        scaler.constraint_scaling_routine(b)
    # after the variable routine the row factor is assigned normally
    scaler.variable_scaling_routine(b)
    scaler.constraint_scaling_routine(b)
    assert scaler.get_scaling_factor(b.flow_consistency_eqn) > 0


# -----------------------------------------------------------------------------
# Diagnostics
@pytest.mark.component
def test_diagnostics_defined_state_block():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.feed = m.params.build_state_block([0], defined_state=True)
    _fix_handcalc_feed(m.feed[0])
    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()


@pytest.mark.component
def test_diagnostics_defined_state_false_after_fix():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    m.out.fix_initialization_states()  # deactivates consistency, fixes vars
    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()


# -----------------------------------------------------------------------------
# Baseline factor source (factor_source="baseline")
@pytest.mark.unit
def test_baseline_matches_seeded_on_fixed_feed():
    # a fully fixed feed block is data, so the gated baseline reads reproduce
    # the seeded factors exactly
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.a = m.params.build_state_block([0], defined_state=True)
    m.b = m.params.build_state_block([0], defined_state=True)
    _fix_handcalc_feed(m.a[0])
    _fix_handcalc_feed(m.b[0])
    BulkPSDScaler().variable_scaling_routine(m.a[0])
    BulkPSDScaler(factor_source="baseline").variable_scaling_routine(m.b[0])
    for name in ("flow_mass_comp", "flow_mass_size"):
        va, vb = getattr(m.a[0], name), getattr(m.b[0], name)
        for idx in va:
            assert get_scaling_factor(vb[idx]) == pytest.approx(
                get_scaling_factor(va[idx]), rel=1e-12
            )


@pytest.mark.unit
def test_baseline_unfixed_block_coarse_default_and_consistency_row():
    # an unfixed (outlet-role) block standalone gets the documented coarse
    # default 1.0 per flow Var (no value reads), and the consistency row then
    # derives from those declared factors
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.out = m.params.build_state_block([0], defined_state=False)
    b = m.out[0]
    scaler = BulkPSDScaler(factor_source="baseline")
    scaler.variable_scaling_routine(b)
    for var in list(b.flow_mass_comp.values()) + list(b.flow_mass_size.values()):
        assert get_scaling_factor(var) == pytest.approx(1.0)
    scaler.constraint_scaling_routine(b)
    # nominal = max(3 comps, 4 bins) x 1.0
    assert get_scaling_factor(b.flow_consistency_eqn) == pytest.approx(0.25)


@pytest.mark.unit
def test_baseline_overwrite_false_never_clobbers():
    m = ConcreteModel()
    m.params = _handcalc_params()
    m.feed = m.params.build_state_block([0], defined_state=True)
    b = m.feed[0]
    _fix_handcalc_feed(b)
    set_scaling_factor(b.flow_mass_comp["Ore1"], 123.0)
    BulkPSDScaler(factor_source="baseline").variable_scaling_routine(b)
    assert get_scaling_factor(b.flow_mass_comp["Ore1"]) == pytest.approx(123.0)
    # every other Var still received a baseline factor
    assert get_scaling_factor(b.flow_mass_comp["Ore2"]) == pytest.approx(1.0 / 3.0)
