#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the mineral-by-size PSD property package."""

import math

from pyomo.environ import ConcreteModel, units, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.scaling import get_scaling_factor, set_scaling_factor

import pytest

from prommis.psd.properties.mineral_size_psd import (
    MineralSizePSDInitializer,
    MineralSizePSDParameterBlock,
    MineralSizePSDScaler,
)

HANDCALC_EDGES_M = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]
HANDCALC_COMPONENTS = ["Ore1", "Ore2", "Ore3"]
HANDCALC_JOINT_FLOWS = (
    (0.5, 0.3, 0.2),
    (1.0, 0.6, 0.4),
    (1.5, 0.9, 0.6),
    (2.0, 1.2, 0.8),
)


def _params(**kwargs):
    return MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3], component_list=["Ore"], **kwargs
    )


def _handcalc_params(**kwargs):
    return MineralSizePSDParameterBlock(
        size_edges=HANDCALC_EDGES_M,
        component_list=HANDCALC_COMPONENTS,
        **kwargs,
    )


def _set_handcalc_feed(state):
    for size, row in enumerate(HANDCALC_JOINT_FLOWS):
        for mineral, flow in zip(HANDCALC_COMPONENTS, row):
            state.flow_mass_size_comp[size, mineral].set_value(flow)


def _fix_complete_handcalc_feed(state):
    _set_handcalc_feed(state)
    for flow in state.flow_mass_size_comp.values():
        flow.fix()
    for flow in state.flow_mass_liquid.values():
        flow.fix()
    state.temperature.fix()
    state.pressure.fix()


def test_liquid_components_default_and_vapor_is_optional():
    model = ConcreteModel()
    model.params = _params()
    model.state = model.params.build_state_block([0], defined_state=True)

    state = model.state[0]
    assert list(model.params.liquid_component_set) == ["H2O"]
    assert list(model.params.vapor_component_set) == []
    assert list(model.params.phase_component_set) == [("Sol", "Ore"), ("Liq", "H2O")]
    assert list(state.flow_mass_liquid) == ["H2O"]
    assert not hasattr(state, "flow_mass_vapor")
    assert state.get_material_flow_terms("Liq", "H2O") is state.flow_mass_liquid["H2O"]


def test_component_indexed_liquid_and_optional_vapor_components():
    model = ConcreteModel()
    model.params = _params(
        liquid_component_list=["H2O", "NaCl"],
        vapor_component_list=["H2O", "CO2"],
    )
    model.state = model.params.build_state_block([0], defined_state=True)

    state = model.state[0]
    assert list(model.params.vapor_component_set) == ["H2O", "CO2"]
    assert list(model.params.phase_component_set) == [
        ("Sol", "Ore"),
        ("Liq", "H2O"),
        ("Liq", "NaCl"),
        ("Vap", "H2O"),
        ("Vap", "CO2"),
    ]
    assert list(state.flow_mass_liquid) == ["H2O", "NaCl"]
    assert list(state.flow_mass_vapor) == ["H2O", "CO2"]
    assert (
        state.get_material_flow_terms("Liq", "NaCl") is state.flow_mass_liquid["NaCl"]
    )
    assert state.get_material_flow_terms("Vap", "H2O") is state.flow_mass_vapor["H2O"]


def test_mineral_specific_cumulative_passing_and_percentile():
    model = ConcreteModel()
    model.params = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 4e-3], component_list=["Fine", "Coarse"]
    )
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]

    state.flow_mass_size_comp[0, "Fine"].set_value(3.0)
    state.flow_mass_size_comp[1, "Fine"].set_value(1.0)
    state.flow_mass_size_comp[0, "Coarse"].set_value(1.0)
    state.flow_mass_size_comp[1, "Coarse"].set_value(3.0)

    assert value(state.cum_passing_mineral[0, "Fine"]) == pytest.approx(0.75)
    assert value(state.cum_passing_mineral[0, "Coarse"]) == pytest.approx(0.25)
    assert value(state.P80_mineral["Fine"]) == pytest.approx(
        value(state.mineral_Pxx("Fine", 0.8))
    )
    assert value(state.P80_mineral["Fine"]) < value(state.P80_mineral["Coarse"])


@pytest.mark.build
@pytest.mark.unit
def test_build_structure():
    model = ConcreteModel()
    model.params = _handcalc_params()
    assert list(model.params.size_interval_set) == [0, 1, 2, 3]
    assert list(model.params.size_edge_index) == [0, 1, 2, 3, 4]
    assert list(model.params.solid_component_set) == HANDCALC_COMPONENTS

    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]
    assert len(state.flow_mass_size_comp) == 12
    assert len(state.flow_mass_liquid) == 1
    assert not hasattr(state, "flow_mass_vapor")
    assert set(state.define_state_vars()) == {
        "flow_mass_size_comp",
        "flow_mass_liquid",
        "temperature",
        "pressure",
    }
    assert_units_consistent(state)
    assert str(units.get_units(state.flow_mass_size_comp[0, "Ore1"])) == str(
        units.kg / units.s
    )


@pytest.mark.unit
def test_param_block_rejects_bad_config():
    invalid_configs = (
        {"size_edges": HANDCALC_EDGES_M, "component_list": []},
        {"size_edges": HANDCALC_EDGES_M, "component_list": ["Ore1", "Ore1"]},
        {"size_edges": HANDCALC_EDGES_M, "component_list": ["Ore1(s)"]},
        {"component_list": ["Ore1"]},
        {"size_edges": [1e-3, 1e-3, 2e-3], "component_list": ["Ore1"]},
        {"size_edges": [1e-3], "component_list": ["Ore1"]},
        {
            "size_edges": [0.0, 1e-3, 2e-3],
            "bottom_size": "x",
            "component_list": ["Ore1"],
        },
        {"size_edges": [1e-3, None], "component_list": ["Ore1"]},
        {"size_edges": [True, 2e-3], "component_list": ["Ore1"]},
        {
            "size_edges": [0.0, 2.0, 4.0],
            "bottom_size": True,
            "component_list": ["Ore1"],
        },
        {"size_edges": HANDCALC_EDGES_M, "component_list": "ABC"},
        {"size_edges": [1, 10**10000], "component_list": ["Ore1"]},
        {"size_edges": 1.0, "component_list": ["Ore1"]},
    )
    for config in invalid_configs:
        with pytest.raises(ConfigurationError):
            model = ConcreteModel()
            model.params = MineralSizePSDParameterBlock(**config)


@pytest.mark.unit
@pytest.mark.parametrize(
    "reserved",
    ["Sol", "size_edges", "d_char", "flow_eps", "flow_mass_size_comp"],
)
def test_param_block_rejects_reserved_component_names(reserved):
    with pytest.raises(ConfigurationError, match="reserved"):
        model = ConcreteModel()
        model.params = MineralSizePSDParameterBlock(
            size_edges=HANDCALC_EDGES_M,
            component_list=["Ore1", reserved],
        )


@pytest.mark.unit
def test_assert_same_mesh_wrong_type_raises_config_error():
    model = ConcreteModel()
    model.params = _handcalc_params()
    with pytest.raises(ConfigurationError, match="mineral-by-size PSD parameter block"):
        model.params.assert_same_mesh("not a parameter block")


@pytest.mark.unit
def test_handcalc_state():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]
    _set_handcalc_feed(state)

    assert value(state.flow_mass) == pytest.approx(10.0, rel=1e-12)
    assert [
        value(state.flow_mass_comp[j]) for j in HANDCALC_COMPONENTS
    ] == pytest.approx([5.0, 3.0, 2.0], rel=1e-12)
    assert [value(state.flow_mass_size[k]) for k in range(4)] == pytest.approx(
        [1.0, 2.0, 3.0, 4.0], rel=1e-12
    )
    cum = [value(state.cum_passing[k]) for k in model.params.size_interval_set]
    assert cum == pytest.approx([0.1, 0.3, 0.6, 1.0], rel=1e-9)
    assert value(units.convert(state.P80, units.mm)) == pytest.approx(12.0, rel=1e-4)
    assert value(units.convert(state.P50, units.mm)) == pytest.approx(
        20.0 / 3.0, rel=1e-4
    )
    assert value(units.convert(model.params.d_char[3], units.mm)) == pytest.approx(
        math.sqrt(128.0), rel=1e-9
    )
    assert state.validate_feed()


@pytest.mark.unit
def test_sum_invariant_and_units():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]
    _set_handcalc_feed(state)

    sum_component_fractions = sum(
        value(state.mass_frac_comp[j]) for j in model.params.solid_component_set
    )
    sum_size_fractions = sum(
        value(state.mass_frac_size[k]) for k in model.params.size_interval_set
    )
    for size in model.params.size_interval_set:
        grade_sum = sum(
            value(state.grade_by_size[size, mineral])
            for mineral in model.params.solid_component_set
        )
        assert grade_sum == pytest.approx(1.0, abs=2e-12)

    assert sum_component_fractions == pytest.approx(1.0, abs=1e-12)
    assert sum_size_fractions == pytest.approx(1.0, abs=1e-12)
    assert value(
        state.cum_passing[model.params.size_interval_set.last()]
    ) == pytest.approx(1.0, abs=1e-12)
    assert_units_consistent(state)


@pytest.mark.unit
def test_expression_safety_on_degenerate_states():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]

    def evaluate_all():
        values = [value(state.flow_mass)]
        values += [value(state.flow_mass_comp[j]) for j in HANDCALC_COMPONENTS]
        values += [
            value(state.flow_mass_size[k]) for k in model.params.size_interval_set
        ]
        values += [value(state.mass_frac_comp[j]) for j in HANDCALC_COMPONENTS]
        values += [
            value(state.mass_frac_size[k]) for k in model.params.size_interval_set
        ]
        values += [
            value(state.grade_by_size[k, j]) for k, j in state.flow_mass_size_comp
        ]
        values += [value(state.cum_passing[k]) for k in model.params.size_interval_set]
        values += [
            value(state.cum_passing_mineral[k, j]) for k, j in state.flow_mass_size_comp
        ]
        values += [value(state.P80), value(state.P50)]
        values += [value(state.P80_mineral[j]) for j in HANDCALC_COMPONENTS]
        values += [value(state.P50_mineral[j]) for j in HANDCALC_COMPONENTS]
        return values

    for case in ("default", "zero", "near_zero"):
        if case == "zero":
            for flow in state.flow_mass_size_comp.values():
                flow.set_value(0.0)
        elif case == "near_zero":
            for flow in state.flow_mass_size_comp.values():
                flow.set_value(1e-15)
        assert all(math.isfinite(result) for result in evaluate_all())


@pytest.mark.unit
def test_validate_feed_rejects_invalid_values():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]
    _set_handcalc_feed(state)
    assert state.validate_feed()

    state.flow_mass_size_comp[1, "Ore1"].fix(-1.0)
    with pytest.raises(ConfigurationError, match="negative"):
        state.validate_feed()

    for flow in state.flow_mass_size_comp.values():
        flow.fix(0.0)
    with pytest.raises(ConfigurationError, match="zero"):
        state.validate_feed()

    _set_handcalc_feed(state)
    state.flow_mass_size_comp[0, "Ore1"].fix(float("inf"))
    with pytest.raises(ConfigurationError, match="non-finite"):
        state.validate_feed()


@pytest.mark.component
def test_initialization_lifecycle_validates_defined_feed():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.feed = model.params.build_state_block([0], defined_state=True)
    state = model.feed[0]
    _fix_complete_handcalc_feed(state)

    MineralSizePSDInitializer().initialize(model.feed)

    assert state.validate_feed()
    assert all(flow.fixed for flow in state.flow_mass_size_comp.values())


@pytest.mark.component
def test_fix_initialization_states_fixes_joint_state():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.out = model.params.build_state_block([0], defined_state=False)

    model.out.fix_initialization_states()

    state = model.out[0]
    assert all(flow.fixed for flow in state.flow_mass_size_comp.values())
    assert all(flow.fixed for flow in state.flow_mass_liquid.values())
    assert state.temperature.fixed
    assert state.pressure.fixed


@pytest.mark.component
def test_initialization_failure_preserves_fixed_feed():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.feed = model.params.build_state_block([0], defined_state=True)
    state = model.feed[0]
    _fix_complete_handcalc_feed(state)
    for flow in state.flow_mass_size_comp.values():
        flow.fix(0.0)

    with pytest.raises(ConfigurationError, match="zero"):
        MineralSizePSDInitializer().initialize(model.feed)

    assert all(flow.fixed for flow in state.flow_mass_size_comp.values())


@pytest.mark.unit
def test_scaler_zero_cell_floor():
    model = ConcreteModel()
    model.params = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 3e-3, 4e-3, 5e-3],
        component_list=["Ore1", "Ore2"],
    )
    model.state = model.params.build_state_block([0], defined_state=True)
    state = model.state[0]
    fractions = ((0.0, 0.0), (0.42, 0.13), (0.21, 0.04), (0.14, 0.06))
    for size, row in enumerate(fractions):
        for mineral, fraction in zip(("Ore1", "Ore2"), row):
            state.flow_mass_size_comp[size, mineral].set_value(100.0 * fraction)

    scaler = MineralSizePSDScaler(factor_source="seeded")
    scaler.variable_scaling_routine(state)
    floor = scaler.EPS_REL * 100.0

    factors = []
    for flow in state.flow_mass_size_comp.values():
        expected = 1.0 / max(abs(value(flow)), floor)
        actual = scaler.get_scaling_factor(flow)
        assert actual == pytest.approx(expected, rel=1e-10)
        factors.append(actual)
    assert max(factors) / min(factors) <= 1.0 / scaler.EPS_REL + 1.0
    scaler.constraint_scaling_routine(state)


@pytest.mark.unit
def test_baseline_matches_seeded_on_fixed_joint_feed():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.seeded = model.params.build_state_block([0], defined_state=True)
    model.baseline = model.params.build_state_block([0], defined_state=True)
    _fix_complete_handcalc_feed(model.seeded[0])
    _fix_complete_handcalc_feed(model.baseline[0])

    MineralSizePSDScaler(factor_source="seeded").variable_scaling_routine(
        model.seeded[0]
    )
    MineralSizePSDScaler(factor_source="baseline").variable_scaling_routine(
        model.baseline[0]
    )

    for name in (
        "flow_mass_size_comp",
        "flow_mass_liquid",
        "temperature",
        "pressure",
    ):
        seeded_vars = getattr(model.seeded[0], name)
        baseline_vars = getattr(model.baseline[0], name)
        if hasattr(seeded_vars, "is_indexed") and not seeded_vars.is_indexed():
            assert get_scaling_factor(baseline_vars) == pytest.approx(
                get_scaling_factor(seeded_vars), rel=1e-12
            )
        else:
            for index in seeded_vars:
                assert get_scaling_factor(baseline_vars[index]) == pytest.approx(
                    get_scaling_factor(seeded_vars[index]), rel=1e-12
                )


@pytest.mark.unit
def test_baseline_unfixed_block_uses_coarse_flow_default():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.out = model.params.build_state_block([0], defined_state=False)
    state = model.out[0]
    scaler = MineralSizePSDScaler(factor_source="baseline")

    scaler.variable_scaling_routine(state)

    for flow in state.flow_mass_size_comp.values():
        assert scaler.get_scaling_factor(flow) == pytest.approx(1.0)
    for flow in state.flow_mass_liquid.values():
        assert scaler.get_scaling_factor(flow) == pytest.approx(1.0)


@pytest.mark.unit
def test_baseline_overwrite_false_preserves_declared_factor():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.feed = model.params.build_state_block([0], defined_state=True)
    state = model.feed[0]
    _fix_complete_handcalc_feed(state)
    set_scaling_factor(state.flow_mass_size_comp[0, "Ore1"], 123.0)

    MineralSizePSDScaler(factor_source="baseline").variable_scaling_routine(state)

    assert get_scaling_factor(state.flow_mass_size_comp[0, "Ore1"]) == pytest.approx(
        123.0
    )
    assert get_scaling_factor(state.flow_mass_size_comp[0, "Ore2"]) == pytest.approx(
        1.0 / 0.3
    )


@pytest.mark.component
def test_diagnostics_defined_state_block():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.feed = model.params.build_state_block([0], defined_state=True)
    _fix_complete_handcalc_feed(model.feed[0])

    DiagnosticsToolbox(model).assert_no_structural_warnings()


@pytest.mark.component
def test_diagnostics_after_fixing_initialization_states():
    model = ConcreteModel()
    model.params = _handcalc_params()
    model.out = model.params.build_state_block([0], defined_state=False)
    model.out.fix_initialization_states()

    DiagnosticsToolbox(model).assert_no_structural_warnings()
