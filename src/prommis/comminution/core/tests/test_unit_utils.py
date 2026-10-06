#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for shared comminution unit-model helpers."""

import pytest

from pyomo.common.collections import ComponentMap
from pyomo.environ import Block, ConcreteModel, Constraint, Set, Var

from idaes.core.scaling import CustomScalerBase
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import get_scaling_factor, set_scaling_factor

from prommis.comminution.core.unit_utils import (
    FactorSource,
    declare_factor_source,
    delegate_state_block_scaling,
    inverse_magnitude_factor,
    inverse_sum_of_nominals,
    propagate_scaling_factors,
    required_scaling_factor,
    scale_constraints_at_time,
)


class _NoOpScaler(CustomScalerBase):
    """Scaler stub using the inherited factor getters and setters."""

    def variable_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        pass

    def constraint_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        pass


def _two_time_model():
    m = ConcreteModel()
    m.time = Set(initialize=[0.0, 1.0])
    m.comp = Set(initialize=["a", "b"])
    m.x = Var(m.time, m.comp, initialize=1.0)
    m.y = Var(m.time, initialize=1.0)
    m.tuple_con = Constraint(m.time, m.comp, rule=lambda b, t, j: b.x[t, j] == 1.0)
    m.scalar_con = Constraint(m.time, rule=lambda b, t: b.y[t] == 2.0)
    return m


@pytest.mark.unit
def test_scale_constraints_at_time_scales_only_target_time():
    # Two times with different factors expose factors reused across time points.
    m = _two_time_model()
    scaler = _NoOpScaler()
    scale_constraints_at_time(
        scaler,
        m,
        [("tuple_con", 0.5), ("scalar_con", 0.25), ("absent_con", 9.9)],
        0.0,
    )
    for j in m.comp:
        assert get_scaling_factor(m.tuple_con[0.0, j]) == pytest.approx(0.5)
        assert get_scaling_factor(m.tuple_con[1.0, j]) is None
    assert get_scaling_factor(m.scalar_con[0.0]) == pytest.approx(0.25)
    assert get_scaling_factor(m.scalar_con[1.0]) is None
    scale_constraints_at_time(
        scaler, m, [("tuple_con", 0.125), ("scalar_con", 0.0625)], 1.0
    )
    for j in m.comp:
        assert get_scaling_factor(m.tuple_con[0.0, j]) == pytest.approx(0.5)
        assert get_scaling_factor(m.tuple_con[1.0, j]) == pytest.approx(0.125)
    assert get_scaling_factor(m.scalar_con[0.0]) == pytest.approx(0.25)
    assert get_scaling_factor(m.scalar_con[1.0]) == pytest.approx(0.0625)
    scale_constraints_at_time(scaler, m, [("scalar_con", 0.5)], 0.0)
    assert get_scaling_factor(m.scalar_con[0.0]) == pytest.approx(0.25)
    scale_constraints_at_time(scaler, m, [("scalar_con", 0.5)], 0.0, overwrite=True)
    assert get_scaling_factor(m.scalar_con[0.0]) == pytest.approx(0.5)


def _two_block_model():
    m = ConcreteModel()
    m.comp = Set(initialize=["a", "b"])
    m.src = Block()
    m.dst = Block()
    for blk in (m.src, m.dst):
        blk.x = Var(m.comp, initialize=1.0)
        blk.y = Var(initialize=1.0)
    return m


@pytest.mark.unit
def test_inverse_sum_of_nominals_sums_and_floors():
    m = _two_time_model()
    scaler = _NoOpScaler()
    set_scaling_factor(m.x[0.0, "a"], 0.5)  # nominal 2.0
    set_scaling_factor(m.x[0.0, "b"], 0.25)  # nominal 4.0
    assert inverse_sum_of_nominals(
        scaler, (m.x[0.0, "a"], m.x[0.0, "b"])
    ) == pytest.approx(1.0 / 6.0, rel=1e-12)
    # The 1e-8 nominal floor caps the inverse at 1e8.
    set_scaling_factor(m.y[0.0], 1.0e12)  # nominal 1e-12
    assert inverse_sum_of_nominals(scaler, (m.y[0.0],)) == pytest.approx(
        1.0e8, rel=1e-12
    )
    # Missing factors raise rather than using current variable values.
    with pytest.raises(ConfigurationError, match="is missing or non-positive"):
        inverse_sum_of_nominals(scaler, (m.x[0.0, "a"], m.x[1.0, "a"]))
    with pytest.raises(ConfigurationError, match="is missing or non-positive"):
        required_scaling_factor(scaler, m.y[1.0])
    set_scaling_factor(m.y[1.0], 0.5)
    assert required_scaling_factor(scaler, m.y[1.0]) == 0.5
    set_scaling_factor(m.y[1.0], 0.0)
    with pytest.raises(ConfigurationError, match="is missing or non-positive"):
        required_scaling_factor(scaler, m.y[1.0])
    assert inverse_magnitude_factor(None, 0.5) == 2.0
    assert inverse_magnitude_factor(0.0, 0.5) == 2.0
    assert inverse_magnitude_factor(-4.0, 1.0) == 0.25
    assert inverse_magnitude_factor(1e-9, 1e-6) == pytest.approx(1e6, rel=1e-12)
    # Every variable in a source group must have a positive scaling factor.
    m = _two_block_model()
    set_scaling_factor(m.src.x["a"], 0.5)
    with pytest.raises(ConfigurationError, match="is missing or non-positive"):
        propagate_scaling_factors(scaler, m.src, m.dst, ("x",))
    set_scaling_factor(m.src.x["b"], 0.25)
    set_scaling_factor(m.dst.x["a"], 9.0)
    propagate_scaling_factors(scaler, m.src, m.dst, ("x", "absent"))
    assert get_scaling_factor(m.dst.x["a"]) == pytest.approx(9.0)
    assert get_scaling_factor(m.dst.x["b"]) == pytest.approx(0.25)
    propagate_scaling_factors(scaler, m.src, m.dst, ("x",), overwrite=True)
    assert get_scaling_factor(m.dst.x["a"]) == pytest.approx(0.5)


class _ConfiguredScaler(CustomScalerBase):
    """Scaler stub declaring factor_source on its own CONFIG."""

    CONFIG = declare_factor_source(CustomScalerBase.CONFIG())

    def variable_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        pass

    def constraint_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        pass


@pytest.mark.unit
def test_declare_factor_source_default_and_construction():
    # Declaring an option on a subclass must leave the shared base CONFIG unchanged.
    assert _ConfiguredScaler().config.factor_source == FactorSource.input_based
    assert (
        _ConfiguredScaler(factor_source="current_values").config.factor_source
        == FactorSource.current_values
    )
    assert (
        _ConfiguredScaler(
            factor_source=FactorSource.current_values
        ).config.factor_source
        == FactorSource.current_values
    )
    assert "factor_source" not in CustomScalerBase.CONFIG()
    with pytest.raises(ValueError):
        _ConfiguredScaler(factor_source="bogus")

    # Pass the parent's mode only to scalers that declare factor_source.
    used = ComponentMap()

    class _RecordingNoOp(_NoOpScaler):
        def variable_scaling_routine(
            self, model, overwrite=False, submodel_scalers=None
        ):
            used[model] = self

    class _RecordingConfigured(_ConfiguredScaler):
        def variable_scaling_routine(
            self, model, overwrite=False, submodel_scalers=None
        ):
            used[model] = self

    m = ConcreteModel()
    m.blk_a = Block()
    m.blk_b = Block()
    parent = _ConfiguredScaler(factor_source="current_values")
    overrides = ComponentMap()
    overrides[m.blk_b] = _RecordingConfigured
    delegate_state_block_scaling(
        parent,
        [m.blk_a, m.blk_b],
        _RecordingNoOp,
        overrides,
        "variable_scaling_routine",
    )
    assert isinstance(used[m.blk_a], _RecordingNoOp)
    assert isinstance(used[m.blk_b], _RecordingConfigured)
    assert used[m.blk_b].config.factor_source == FactorSource.current_values
    inst = _RecordingNoOp()
    overrides[m.blk_a] = inst
    delegate_state_block_scaling(
        parent, [m.blk_a], _RecordingNoOp, overrides, "variable_scaling_routine"
    )
    assert used[m.blk_a] is inst
