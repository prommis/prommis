#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Public tests for the per-component slurry PSD property package."""

import math

import pytest

from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TransformationFactory,
    units,
    value,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import ControlVolume0DBlock, FlowsheetBlock, MaterialBalanceType
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import get_scaling_factor

from prommis.comminution.properties.solid_psd_properties import (
    SolidPSDInitializer,
    SolidPSDParameterBlock,
    SolidPSDScaler,
)

# Ascending mesh edges in meters; densities in kg/m^3.
EDGES = [1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3, 16.0e-3]
SIZED = ["OreA", "OreB"]
UNSIZED = ["Inert"]
DENSITY = {"OreA": 3000.0, "OreB": 5000.0, "Inert": 2000.0}


def _params(**overrides):
    kwargs = dict(
        size_edges=EDGES,
        sized_solid_component_list=list(SIZED),
        unsized_solid_component_list=list(UNSIZED),
        solid_density=dict(DENSITY),
    )
    kwargs.update(overrides)
    m = ConcreteModel()
    m.params = SolidPSDParameterBlock(**kwargs)
    return m


CONFIG_REJECTION_CASES = [
    ("size_edges_missing", dict(size_edges=None), "requires a 'size_edges' list"),
    ("size_edges_too_short", dict(size_edges=[1.0e-3]), "at least 2 edges"),
    (
        "size_edges_decreasing",
        dict(size_edges=[2.0e-3, 1.0e-3]),
        "strictly increasing",
    ),
    (
        "zero_edge_without_bottom_size",
        dict(size_edges=[0.0, 1.0e-3, 2.0e-3]),
        "supply bottom_size when x_0 == 0",
    ),
    (
        "solid_list_missing",
        dict(sized_solid_component_list=None),
        "sized_solid_component_list is required",
    ),
    (
        "solid_list_empty",
        dict(sized_solid_component_list=[]),
        "sized_solid_component_list must be non-empty",
    ),
    (
        "solid_list_not_sequence",
        dict(sized_solid_component_list="OreA"),
        "must be a list or tuple",
    ),
    (
        "solid_list_duplicate",
        dict(sized_solid_component_list=["OreA", "OreA"]),
        "duplicate component name",
    ),
    (
        "solid_name_not_identifier",
        dict(sized_solid_component_list=["Ore1(s)"]),
        "valid Python identifier",
    ),
    (
        "sized_and_unsized_overlap",
        dict(unsized_solid_component_list=["OreA"]),
        "both sized and unsized solid lists",
    ),
    (
        "unsized_list_duplicate",
        dict(unsized_solid_component_list=["Inert", "Inert"]),
        "duplicate component name",
    ),
    (
        "liquid_list_duplicate",
        dict(liquid_component_list=["H2O", "H2O"]),
        "duplicate component name",
    ),
    (
        "vapor_list_duplicate",
        dict(vapor_component_list=["H2O", "H2O"]),
        "duplicate component name",
    ),
    (
        "solid_density_missing",
        dict(solid_density=None),
        "solid_density is required as a dict",
    ),
    (
        "solid_density_key_missing",
        dict(solid_density={"OreA": 3000.0, "OreB": 5000.0}),
        r"missing: \['Inert'\]",
    ),
    (
        "solid_density_key_unexpected",
        dict(solid_density={**DENSITY, "Extra": 1000.0}),
        r"unexpected: \['Extra'\]",
    ),
    (
        "solid_density_negative",
        dict(solid_density={**DENSITY, "OreA": -3000.0}),
        "positive finite density",
    ),
    (
        "solid_density_nonfinite",
        dict(solid_density={**DENSITY, "OreA": float("nan")}),
        "positive finite density",
    ),
    (
        "liquid_density_missing",
        dict(liquid_component_list=["Brine"]),
        "liquid_density is required",
    ),
    (
        "reserved_name_liquid",
        dict(liquid_component_list=["Liq"], liquid_density={"Liq": 1000.0}),
        "'Liq' .* is reserved",
    ),
    (
        "reserved_name_solid",
        dict(
            sized_solid_component_list=["Sol", "OreB"],
            solid_density={"Sol": 3000.0, "OreB": 5000.0, "Inert": 2000.0},
        ),
        "'Sol' .* is reserved",
    ),
    (
        "reserved_name_vapor",
        dict(vapor_component_list=["Vap"]),
        "'Vap' .* is reserved",
    ),
    (
        "reserved_name_component_list",
        dict(
            sized_solid_component_list=["component_list", "OreB"],
            solid_density={"component_list": 3000.0, "OreB": 5000.0, "Inert": 2000.0},
        ),
        "'component_list' .* is reserved",
    ),
    (
        "block_attribute_collision",
        dict(
            sized_solid_component_list=["config", "OreB"],
            solid_density={"config": 3000.0, "OreB": 5000.0, "Inert": 2000.0},
        ),
        "collides with an existing block attribute",
    ),
    (
        "liquid_density_bool",
        dict(liquid_component_list=["Brine"], liquid_density={"Brine": True}),
        r"liquid_density\['Brine'\] must be an int or float \(not bool\)",
    ),
    (
        "percentile_targets_empty",
        dict(percentile_targets=[]),
        "non-empty list or tuple",
    ),
    ("percentile_targets_zero", dict(percentile_targets=[0.0]), "must be > 0"),
    ("percentile_targets_above_one", dict(percentile_targets=[1.1]), "must be <= 1.0"),
    (
        "percentile_targets_duplicate",
        dict(percentile_targets=[0.8, 0.8]),
        "must not contain duplicates",
    ),
    ("percentile_targets_boolean", dict(percentile_targets=[True]), "not a boolean"),
    ("percentile_targets_nan", dict(percentile_targets=[float("nan")]), "finite"),
]


@pytest.mark.unit
def test_invalid_configurations_rejected():
    failures = []
    for case_id, overrides, match in CONFIG_REJECTION_CASES:
        try:
            with pytest.raises(ConfigurationError, match=match):
                _params(**overrides)
        except (Exception, pytest.fail.Exception) as exc:
            failures.append(f"{case_id}: {type(exc).__name__}: {exc}")
    try:
        m = _params(
            liquid_component_list=["Brine"],
            liquid_density={"Brine": 1100.0},
            vapor_component_list=["Air"],
        )
        assert m.params.liquid_list == ("Brine",)
        assert m.params.vapor_list == ("Air",)
        assert value(m.params.dens_mass_liquid_comp["Brine"]) == pytest.approx(1100.0)
    except (Exception, pytest.fail.Exception) as exc:
        failures.append(
            f"custom_liquid_and_vapor_accepted: {type(exc).__name__}: {exc}"
        )
    assert not failures, "\n".join(failures)


@pytest.mark.unit
def test_bottom_size_and_assert_same_mesh():
    m1 = _params(size_edges=[0.0] + EDGES[1:], bottom_size=0.5e-3)
    # Use bottom_size as the lower edge of the first geometric mean.
    assert value(m1.params.size_char[0]) == pytest.approx(math.sqrt(0.5e-3 * 2.0e-3))
    m2 = _params()
    with pytest.raises(ConfigurationError):
        m1.params.assert_same_mesh(m2.params)
    m3 = _params()
    m2.params.assert_same_mesh(m3.params)
    with pytest.raises(ConfigurationError):
        m2.params.assert_same_mesh(object())
    # Shared edges do not make meshes of different lengths equal.
    with pytest.raises(ConfigurationError, match="size-mesh mismatch"):
        m2.params.assert_same_mesh(_params(size_edges=EDGES[:-1]).params)
    # A 1e-13 m edge shift remains within the mesh tolerance.
    m2.params.assert_same_mesh(
        _params(size_edges=EDGES[:-1] + [EDGES[-1] + 1e-13]).params
    )


# Feed flows in kg/s: sized solids 2.0, Inert 0.5, H2O 2.5.
FLOWS_A = [0.1, 0.2, 0.3, 0.4]
FLOWS_B = [0.4, 0.3, 0.2, 0.1]


def _state(m=None, **overrides):
    m = m if m is not None else _params(**overrides)
    m.state = m.params.build_state_block([0], defined_state=True)
    return m, m.state[0]


def _fix_hand_calc(b):
    for k, (va, vb) in enumerate(zip(FLOWS_A, FLOWS_B)):
        b.flow_mass_sized_comp_size["OreA", k].fix(va)
        b.flow_mass_sized_comp_size["OreB", k].fix(vb)
    b.flow_mass_unsized_comp["Inert"].fix(0.5)
    b.flow_mass_liquid_comp["H2O"].fix(2.5)
    b.temperature.fix(298.15)
    b.pressure.fix(101325.0)


DEFERRED_PROPERTY_CASES = [
    ("percentile_size_comp", ("OreA", 0.8)),
    ("cum_passing_comp_size", ("OreA", 2)),
    ("mass_frac_size_comp", ("OreA", 3)),
    ("mass_frac_comp_size", ("OreA", 3)),
    ("mass_frac_size", 3),
    ("mass_frac_solid_comp", "OreA"),
    ("dens_mass_solid_mass_weighted", None),
]


@pytest.mark.component
@pytest.mark.parametrize("name,index", DEFERRED_PROPERTY_CASES)
def test_optional_property_constructed_on_access(name, index):
    m = _params(vapor_component_list=["Air"])
    m.state = m.params.build_state_block([0, 1], defined_state=True)
    deferred = {name for name, _ in DEFERRED_PROPERTY_CASES}
    for b in m.state.values():
        _fix_hand_calc(b)
        b.flow_mass_vapor_comp["Air"].fix(0.1)
        assert not any(b.is_property_constructed(prop) for prop in deferred)

    SolidPSDInitializer().initialize(m.state)
    for b in m.state.values():
        for factor_source in ("current_values", "input_based"):
            SolidPSDScaler(factor_source=factor_source).variable_scaling_routine(b)
        assert not any(b.is_property_constructed(prop) for prop in deferred)

    b = m.state[0]
    component = getattr(b, name)
    assert b.is_property_constructed(name)
    assert not any(m.state[1].is_property_constructed(prop) for prop in deferred)

    output = component if index is None else component[index]
    before = value(output)
    b.flow_mass_sized_comp_size["OreA", 3].fix(0.8)
    SolidPSDInitializer().initialize(m.state)
    assert value(output) != pytest.approx(before)


@pytest.mark.unit
def test_hand_calc_aggregates():
    m, b = _state(percentile_targets=(0.5, 0.8))
    _fix_hand_calc(b)
    assert value(b.flow_mass_phase["Sol"]) == pytest.approx(2.5, rel=1e-9)
    assert value(b.flow_mass_phase["Liq"]) == pytest.approx(2.5, rel=1e-9)
    assert set(b.flow_mass_phase) == {"Sol", "Liq"}
    assert value(b.flow_mass_solid_comp["OreA"]) == pytest.approx(1.0, rel=1e-9)
    assert value(b.flow_mass_solid_comp["Inert"]) == pytest.approx(0.5, rel=1e-9)
    for k in range(4):
        assert value(b.flow_mass_size[k]) == pytest.approx(0.5, rel=1e-9)
        assert value(b.mass_frac_size[k]) == pytest.approx(0.25, rel=1e-6)
    assert value(b.mass_frac_solid_comp["OreA"]) == pytest.approx(0.4, rel=1e-6)
    # The solid phase includes unsized solids, so it exceeds the sized total.
    assert value(b.flow_mass_phase["Sol"]) > sum(
        value(b.flow_mass_size[k]) for k in range(4)
    )
    # Bulk PSDs exclude unsized solids.
    assert [value(b.cum_passing_size[k]) for k in range(4)] == pytest.approx(
        [0.25, 0.5, 0.75, 1.0], rel=1e-6
    )
    # Linear interpolation gives bulk P80 = 8 + 0.2 * (16 - 8) = 9.6 mm.
    assert value(units.convert(b.percentile_size[0.8], units.mm)) == pytest.approx(
        9.6, rel=1e-3
    )
    assert value(units.convert(b.percentile_size[0.5], units.mm)) == pytest.approx(
        4.0, rel=1e-3
    )
    assert list(m.params.percentile_targets) == [0.5, 0.8]
    assert list(b.percentile_size) == [0.5, 0.8]
    # Size fractions use the component total; composition uses the bin total.
    assert value(b.mass_frac_size_comp["OreA", 3]) == pytest.approx(0.4, rel=1e-6)
    assert [value(b.mass_frac_comp_size["OreA", k]) for k in range(4)] == pytest.approx(
        [0.2, 0.4, 0.6, 0.8], rel=1e-6
    )
    assert [value(b.mass_frac_comp_size["OreB", k]) for k in range(4)] == pytest.approx(
        [0.8, 0.6, 0.4, 0.2], rel=1e-6
    )
    assert [
        value(b.cum_passing_comp_size["OreA", k]) for k in range(4)
    ] == pytest.approx([0.1, 0.3, 0.6, 1.0], rel=1e-6)
    # OreA P80 = 8 + 0.5 * (16 - 8) = 12 mm.
    assert value(
        units.convert(b.percentile_size_comp["OreA", 0.8], units.mm)
    ) == pytest.approx(12.0, rel=1e-3)
    # OreB P80 = 4 + 0.5 * (8 - 4) = 6 mm.
    assert value(
        units.convert(b.percentile_size_comp["OreB", 0.8], units.mm)
    ) == pytest.approx(6.0, rel=1e-3)
    # P50: OreA = 4 + (2/3) * (8 - 4); OreB = 2 + (1/3) * (4 - 2), in mm.
    assert value(
        units.convert(b.percentile_size_comp["OreA", 0.5], units.mm)
    ) == pytest.approx(20.0 / 3.0, rel=1e-3)
    assert value(
        units.convert(b.percentile_size_comp["OreB", 0.5], units.mm)
    ) == pytest.approx(8.0 / 3.0, rel=1e-3)
    for s in SIZED:
        b.flow_mass_sized_comp_size[s, 1].fix(0.0)
    assert [value(b.mass_frac_comp_size[s, 1]) for s in SIZED] == pytest.approx(
        [0.0, 0.0], abs=1e-12
    )


@pytest.mark.unit
def test_configured_percentile_sizes():
    m, b = _state(percentile_targets=(0.5, 0.9))
    _fix_hand_calc(b)

    assert list(m.params.percentile_targets) == [0.5, 0.9, 0.8]
    assert list(b.percentile_size) == [0.5, 0.9, 0.8]
    assert value(units.convert(b.percentile_size[0.5], units.mm)) == pytest.approx(
        4.0, rel=1e-3
    )
    assert value(units.convert(b.percentile_size[0.9], units.mm)) == pytest.approx(
        12.8, rel=1e-3
    )
    assert value(
        units.convert(b.percentile_size_comp["OreA", 0.5], units.mm)
    ) == pytest.approx(20.0 / 3.0, rel=1e-3)
    assert value(
        units.convert(b.percentile_size_comp["OreA", 0.9], units.mm)
    ) == pytest.approx(14.0, rel=1e-3)
    assert value(
        units.convert(b.percentile_size_comp["OreB", 0.9], units.mm)
    ) == pytest.approx(8.0, rel=1e-3)
    assert value(units.convert(b.percentile_size[0.8], units.mm)) == pytest.approx(
        9.6, rel=1e-3
    )


def _percentile_report_state():
    m, b = _state(
        liquid_component_list=["H2O", "Brine"],
        liquid_density={"H2O": 997.048, "Brine": 1100.0},
        vapor_component_list=["Air"],
    )
    _fix_hand_calc(b)
    for comp, row in {"OreA": [1, 1, 0, 0], "OreB": [0, 0, 0, 3]}.items():
        for k, flow in enumerate(row):
            b.flow_mass_sized_comp_size[comp, k].fix(flow)
    b.flow_mass_liquid_comp["H2O"].fix(2.0)
    b.flow_mass_liquid_comp["Brine"].fix(1.0)
    b.flow_mass_vapor_comp["Air"].fix(0.2)
    return m, b


@pytest.mark.unit
@pytest.mark.parametrize(
    "comp,target,reference_mm",
    [
        (None, 0.5, 28.0 / 3.0),
        (None, 0.8, 40.0 / 3.0),
        (None, 1.0, 16.0),
        ("OreA", 0.5, 2.0),
        ("OreA", 0.8, 3.2),
        ("OreA", 1.0, 4.0),
    ],
)
def test_percentile_report_reference_and_state(comp, target, reference_mm):
    _, b = _percentile_report_state()
    b.flow_mass_sized_comp_size["OreA", 0].unfix()
    b.temperature.unfix()
    before = {
        var.name: (var.value, var.fixed)
        for component in b.define_state_vars().values()
        for var in component.values()
    }
    report = b.percentile_size_report(target, comp=comp, rel_tol=0.01)
    assert report["status"] == "reliable"
    assert report["reference_size_m"] * 1e3 == pytest.approx(reference_mm, rel=1e-12)
    assert report["size_m"] * 1e3 == pytest.approx(reference_mm, rel=0.01)
    assert report["smooth_size_m"] == report["size_m"]
    assert report["relative_error"] == pytest.approx(
        abs(report["smooth_size_m"] / report["reference_size_m"] - 1.0), abs=1e-14
    )
    assert {
        var.name: (var.value, var.fixed)
        for component in b.define_state_vars().values()
        for var in component.values()
    } == before


@pytest.mark.unit
@pytest.mark.parametrize("comp", [None, "OreB"])
@pytest.mark.parametrize(
    "flow", [0.0, 1e-12, 1.0], ids=["empty", "near-floor", "healthy"]
)
def test_percentile_report_reliability(comp, flow):
    _, b = _percentile_report_state()
    row = [0.4 * flow, 0.6 * flow, 0.0, 0.0]
    rows = (
        {"OreA": row, "OreB": [0.0] * 4}
        if comp is None
        else {"OreA": [0.4, 0.6, 0.0, 0.0], "OreB": row}
    )
    for component, flows in rows.items():
        for k, bin_flow in enumerate(flows):
            b.flow_mass_sized_comp_size[component, k].fix(bin_flow)
    for var in b.flow_mass_liquid_comp.values():
        var.fix(2.0)
    report = b.percentile_size_report(0.8, comp=comp, rel_tol=0.01)
    assert math.isfinite(report["smooth_size_m"])
    if flow == 0.0:
        assert report["status"] == "no_particles"
        assert report["size_m"] is None
        assert report["reference_size_m"] is None
        assert report["relative_error"] is None
    else:
        # Normalization without a flow floor makes this reference scale invariant.
        assert report["reference_size_m"] * 1e3 == pytest.approx(10.0 / 3.0)
        reference_m = 1.0 / 300.0
        expected_error = abs(report["smooth_size_m"] - reference_m) / reference_m
        assert report["relative_error"] == pytest.approx(expected_error, abs=1e-14)
        reliable = expected_error <= 0.01
        assert report["status"] == ("reliable" if reliable else "unreliable")
        if reliable:
            assert report["size_m"] * 1e3 == pytest.approx(10.0 / 3.0, rel=0.01)
        else:
            assert report["size_m"] is None
    if comp is not None:
        assert b.percentile_size_report(0.8, rel_tol=0.01)["status"] == "reliable"


@pytest.mark.unit
def test_size_at_passing_from_flows_matches_state_p80_at_small_flow():
    m, sb = _state()
    pp = m.params
    order = list(pp.size_interval_set)
    size_bin_flows = {}
    for i, s in enumerate(pp.sized_solid_list):
        size_bin_flows[s] = [1e-11 * (i + 1) * (k + 1) for k in order]
        for k in order:
            sb.flow_mass_sized_comp_size[s, k].fix(size_bin_flows[s][k])
    assert pp.size_at_passing_from_flows(size_bin_flows) == pytest.approx(
        value(sb.percentile_size[0.8]), rel=1e-12
    )


@pytest.mark.unit
def test_slurry_volume_and_density_expressions():
    m, b = _state()
    _fix_hand_calc(b)
    vol_solid = 1.0 / 3000.0 + 1.0 / 5000.0 + 0.5 / 2000.0
    vol_liquid = 2.5 / 997.048
    assert value(b.flow_vol_solid) == pytest.approx(vol_solid, rel=1e-9)
    assert value(b.flow_vol_liquid) == pytest.approx(vol_liquid, rel=1e-9)
    assert value(b.flow_vol_slurry) == pytest.approx(vol_solid + vol_liquid, rel=1e-9)
    assert value(b.vol_frac_solid_slurry) == pytest.approx(
        vol_solid / (vol_solid + vol_liquid), rel=1e-6
    )
    assert value(b.dens_mass_solid) == pytest.approx(2.5 / vol_solid, rel=1e-6)
    mw = (1.0 * 3000.0 + 1.0 * 5000.0 + 0.5 * 2000.0) / 2.5
    assert value(b.dens_mass_solid_mass_weighted) == pytest.approx(mw, rel=1e-6)
    assert value(b.dens_mass_slurry) == pytest.approx(
        5.0 / (vol_solid + vol_liquid), rel=1e-6
    )


@pytest.mark.unit
def test_solid_density_remains_accurate_at_low_flow():
    _, b = _state(solid_density={**DENSITY, "OreA": 7500.0})
    for var in (
        b.flow_mass_sized_comp_size,
        b.flow_mass_unsized_comp,
        b.flow_mass_liquid_comp,
    ):
        for item in var.values():
            item.fix(0.0)
    b.flow_mass_sized_comp_size["OreA", 0].fix(1e-6)
    assert b.validate_feed()
    assert value(b.dens_mass_solid) == pytest.approx(7500.0, rel=1e-5)


@pytest.mark.unit
def test_water_only_feed_admitted_and_zero_rejected():
    m, b = _state()
    for k in range(4):
        b.flow_mass_sized_comp_size["OreA", k].fix(0.0)
        b.flow_mass_sized_comp_size["OreB", k].fix(0.0)
    b.flow_mass_unsized_comp["Inert"].fix(0.0)
    b.flow_mass_liquid_comp["H2O"].fix(1.0)
    b.temperature.fix(298.15)
    b.pressure.fix(101325.0)
    assert b.validate_feed()
    b.flow_mass_liquid_comp["H2O"].fix(0.0)
    with pytest.raises(ConfigurationError):
        b.validate_feed()
    b.flow_mass_liquid_comp["H2O"].fix(1e-13)
    assert b.validate_feed()
    for spoil, match in (
        (
            lambda b: b.flow_mass_sized_comp_size["OreA", 0].fix(-1.0),
            "negative flow value",
        ),
        (
            lambda b: b.flow_mass_liquid_comp["H2O"].fix(float("nan")),
            "non-finite flow value",
        ),
        (
            lambda b: b.flow_mass_unsized_comp["Inert"].fix(float("inf")),
            "non-finite flow value",
        ),
        (lambda b: b.temperature.fix(0.0), "temperature must be a positive finite"),
        (lambda b: b.pressure.fix(-101325.0), "pressure must be a positive finite"),
    ):
        m, b = _state()
        _fix_hand_calc(b)
        spoil(b)
        with pytest.raises(ConfigurationError, match=match):
            b.validate_feed()


@pytest.mark.unit
def test_expression_safety_on_zero_state():
    m, b = _state(vapor_component_list=["Air"], percentile_targets=(0.5,))
    for var in b.define_state_vars().values():
        for idx in var:
            var[idx].fix(0.0)
    names = [
        "flow_mass_phase",
        "flow_vol_solid",
        "flow_vol_liquid",
        "flow_vol_slurry",
        "vol_frac_solid_slurry",
        "dens_mass_solid",
        "dens_mass_solid_mass_weighted",
        "dens_mass_slurry",
        "percentile_size",
    ]
    for name in names:
        for expr in getattr(b, name).values():
            assert math.isfinite(value(expr))
    # Additive denominator floors keep zero-flow fractions and densities at zero.
    for name in (
        "flow_mass_phase",
        "flow_vol_solid",
        "flow_vol_liquid",
        "flow_vol_slurry",
        "vol_frac_solid_slurry",
        "dens_mass_solid",
        "dens_mass_solid_mass_weighted",
        "dens_mass_slurry",
    ):
        for expr in getattr(b, name).values():
            assert value(expr) == pytest.approx(0.0, abs=1e-12)
    for k in range(4):
        assert math.isfinite(value(b.mass_frac_size[k]))
        assert value(b.mass_frac_size[k]) == pytest.approx(0.0, abs=1e-12)
        assert value(b.cum_passing_size[k]) == pytest.approx(0.0, abs=1e-12)
    for s in SIZED:
        assert math.isfinite(value(b.percentile_size_comp[s, 0.5]))
        assert math.isfinite(value(b.percentile_size_comp[s, 0.8]))
        for k in range(4):
            assert value(b.mass_frac_size_comp[s, k]) == pytest.approx(0.0, abs=1e-12)
            assert value(b.mass_frac_comp_size[s, k]) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.unit
def test_material_flow_terms_and_units():
    m, b = _state(vapor_component_list=["Air"])
    _fix_hand_calc(b)
    b.flow_mass_vapor_comp["Air"].fix(0.1)
    totals = {
        ("Sol", "OreA"): 1.0,
        ("Sol", "OreB"): 1.0,
        ("Sol", "Inert"): 0.5,
        ("Liq", "H2O"): 2.5,
        ("Vap", "Air"): 0.1,
    }
    for phase, expected in (
        ("Sol", (True, False, False)),
        ("Liq", (False, True, False)),
        ("Vap", (False, False, True)),
    ):
        obj = m.params.get_phase(phase)
        assert (
            obj.is_solid_phase(),
            obj.is_liquid_phase(),
            obj.is_vapor_phase(),
        ) == expected
    # Only declared phase-component pairs belong in material balances.
    assert set(m.params.get_phase_component_set()) == set(totals)
    assert b.get_material_flow_basis().name == "mass"
    for pair, total in totals.items():
        term = b.get_material_flow_terms(*pair)
        assert_units_equivalent(term, units.kg / units.s)
        assert value(term) == pytest.approx(total, rel=1e-12)
    assert set(b.flow_mass_phase) == {"Sol", "Liq", "Vap"}
    for phase, total in {"Sol": 2.5, "Liq": 2.5, "Vap": 0.1}.items():
        assert_units_equivalent(b.flow_mass_phase[phase], units.kg / units.s)
        assert value(b.flow_mass_phase[phase]) == pytest.approx(total, rel=1e-12)
    assert_units_consistent(b)
    _, b2 = _state()
    with pytest.raises(KeyError):
        b2.get_material_flow_terms("Vap", "Air")
    assert "flow_mass_vapor_comp" not in b2.define_state_vars()
    m3 = _params()
    m3.fs = FlowsheetBlock(dynamic=False)
    m3.fs.cv = ControlVolume0DBlock(property_package=m3.params)
    m3.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m3.fs.cv.add_material_balances(balance_type=MaterialBalanceType.useDefault)
    assert {j for _, j in m3.fs.cv.material_balances} == {*SIZED, *UNSIZED, "H2O"}


@pytest.mark.component
def test_shared_components_have_independent_phase_flows_and_balances():
    m = _params(
        sized_solid_component_list=["NaCl"],
        unsized_solid_component_list=[],
        solid_density={"NaCl": 2160.0},
        liquid_component_list=["H2O", "NaCl"],
        liquid_density={"H2O": 997.048, "NaCl": 1200.0},
        vapor_component_list=["H2O"],
    )
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.cv = ControlVolume0DBlock(property_package=m.params)
    cv = m.fs.cv
    cv.add_state_blocks(has_phase_equilibrium=False)
    cv.add_material_balances(balance_type=MaterialBalanceType.useDefault)
    for state in (cv.properties_in[0], cv.properties_out[0]):
        for k, flow in enumerate(FLOWS_A):
            state.flow_mass_sized_comp_size["NaCl", k].fix(flow)
        state.flow_mass_liquid_comp["H2O"].fix(2.5)
        state.flow_mass_liquid_comp["NaCl"].fix(0.6)
        state.flow_mass_vapor_comp["H2O"].fix(0.1)
    assert set(m.params.component_list) == {"NaCl", "H2O"}
    totals = {
        ("Sol", "NaCl"): 1.0,
        ("Liq", "H2O"): 2.5,
        ("Liq", "NaCl"): 0.6,
        ("Vap", "H2O"): 0.1,
    }
    assert set(m.params.get_phase_component_set()) == set(totals)
    for pair, total in totals.items():
        assert value(
            cv.properties_in[0].get_material_flow_terms(*pair)
        ) == pytest.approx(total)
    for state in (cv.properties_in[0], cv.properties_out[0]):
        for phase, total in {"Sol": 1.0, "Liq": 3.1, "Vap": 0.1}.items():
            assert value(state.flow_mass_phase[phase]) == pytest.approx(total)
    cv.properties_out[0].flow_mass_vapor_comp["H2O"].fix(0.2)
    cv.properties_out[0].flow_mass_liquid_comp["H2O"].fix(2.4)
    assert value(cv.properties_out[0].flow_mass_phase["Liq"]) == pytest.approx(3.0)
    assert value(cv.properties_out[0].flow_mass_phase["Vap"]) == pytest.approx(0.2)
    assert value(cv.material_balances[0, "H2O"].body) == pytest.approx(0.0, abs=1e-12)
    cv.properties_out[0].flow_mass_liquid_comp["NaCl"].fix(0.4)
    assert abs(value(cv.material_balances[0, "NaCl"].body)) == pytest.approx(0.2)
    cv.properties_out[0].flow_mass_sized_comp_size["NaCl", 0].fix(0.3)
    assert value(cv.properties_out[0].flow_mass_phase["Sol"]) == pytest.approx(1.2)
    assert value(cv.material_balances[0, "NaCl"].body) == pytest.approx(0.0, abs=1e-12)


flow_mass = units.kg / units.s
flow_vol = units.m**3 / units.s
density = units.kg / units.m**3

METADATA_UNITS = {
    "temperature": units.K,
    "pressure": units.Pa,
    "flow_mass_sized_comp_size": flow_mass,
    "flow_mass_unsized_comp": flow_mass,
    "flow_mass_liquid_comp": flow_mass,
    "flow_mass_vapor_comp": flow_mass,
    "flow_mass_size": flow_mass,
    "flow_mass_sized": flow_mass,
    "flow_mass_phase": flow_mass,
    "flow_mass_sized_comp": flow_mass,
    "flow_mass_solid_comp": flow_mass,
    "mass_frac_solid_comp": units.dimensionless,
    "mass_frac_size": units.dimensionless,
    "mass_frac_size_comp": units.dimensionless,
    "mass_frac_comp_size": units.dimensionless,
    "cum_passing_size": units.dimensionless,
    "cum_passing_comp_size": units.dimensionless,
    "percentile_size": units.m,
    "percentile_size_comp": units.m,
    "flow_vol_solid": flow_vol,
    "flow_vol_liquid": flow_vol,
    "flow_vol_slurry": flow_vol,
    "vol_frac_solid_slurry": units.dimensionless,
    "dens_mass_solid": density,
    "dens_mass_solid_mass_weighted": density,
    "dens_mass_slurry": density,
}


@pytest.mark.build
@pytest.mark.unit
def test_metadata_units_supported_flags_and_indices():
    m, b = _state(vapor_component_list=["Air"])
    p = m.params
    assert b.config.parameters.get_metadata() is p.get_metadata()
    props = p.get_metadata().properties
    for name, expected_units in METADATA_UNITS.items():
        meta = props[name]
        assert meta.name == name
        assert meta.supported
        assert_units_equivalent(meta.units, expected_units)
        comp = getattr(b, name)
        data = next(iter(comp.values())) if comp.is_indexed() else comp
        assert_units_equivalent(data, expected_units)
    # Custom solids properties leave the standard IDAES flow and fraction
    # entries unsupported.
    assert not props.flow_mass["comp"].supported
    assert not props.flow_mass["none"].supported
    assert not props.mass_frac["comp"].supported
    assert p.sized_solid_list == ("OreA", "OreB")
    assert p.unsized_solid_list == ("Inert",)
    assert p.all_solid_list == ("OreA", "OreB", "Inert")
    assert p.liquid_list == ("H2O",)
    assert p.vapor_list == ("Air",)
    assert value(p.dens_mass_liquid_comp["H2O"]) == pytest.approx(997.048)
    assert set(b.define_state_vars()) == {
        "flow_mass_sized_comp_size",
        "flow_mass_unsized_comp",
        "flow_mass_liquid_comp",
        "flow_mass_vapor_comp",
        "temperature",
        "pressure",
    }
    solids = {"OreA", "OreB", "Inert"}
    assert set(b.flow_mass_solid_comp) == solids
    assert set(b.mass_frac_solid_comp) == solids
    assert set(b.flow_mass_sized_comp_size) == {(s, k) for s in SIZED for k in range(4)}
    assert set(b.mass_frac_comp_size) == {(s, k) for s in SIZED for k in range(4)}
    assert set(b.percentile_size_comp) == {(s, 0.8) for s in SIZED}
    assert set(b.flow_mass_unsized_comp) == set(UNSIZED)
    assert set(b.flow_mass_liquid_comp) == {"H2O"}
    assert set(b.flow_mass_vapor_comp) == {"Air"}
    assert set(b.flow_mass_size) == set(range(4))
    with pytest.raises(KeyError):
        b.flow_mass_solid_comp["H2O"]
    assert m.state.default_initializer is SolidPSDInitializer
    assert m.state.default_scaler is SolidPSDScaler


@pytest.mark.component
def test_port_arc_expansion_carries_full_state():
    m = _params(vapor_component_list=["Air"])
    m.sb = m.params.build_state_block([0, 1], defined_state=True)
    m.port_out, _ = m.sb.build_port(index=0, slice_index=0)
    m.port_in, _ = m.sb.build_port(index=1, slice_index=1)
    for t, port in ((0, m.port_out), (1, m.port_in)):
        members = m.sb[t].define_state_vars()
        assert set(port.vars) == set(members)
        for name, component in members.items():
            assert set(port.vars[name]) == set(component)
            for index, var in component.items():
                assert port.vars[name][index] is var
    m.arc = Arc(source=m.port_out, destination=m.port_in)
    TransformationFactory("network.expand_arcs").apply_to(m)
    linked = ComponentSet()
    for con in m.arc_expanded.component_data_objects(Constraint):
        linked.update(identify_variables(con.body))
    for t in (0, 1):
        for component in m.sb[t].define_state_vars().values():
            for var in component.values():
                assert var in linked, var.name


@pytest.mark.component
def test_initializer_validates_and_restores():
    m, b = _state()
    _fix_hand_calc(b)
    init = SolidPSDInitializer()
    init.initialize(m.state)  # Defined feed states are validated without a solver.
    b.flow_mass_sized_comp_size["OreA", 0].fix(-1.0)
    b.temperature.unfix()
    with pytest.raises(ConfigurationError, match="negative flow value"):
        init.initialize(m.state)
    # InitializerBase restores the fixedness present before the failed call.
    assert b.flow_mass_sized_comp_size["OreA", 0].fixed
    assert not b.temperature.fixed


@pytest.mark.unit
def test_scaler_stream_relative_floor_and_tp():
    m, b = _state()
    _fix_hand_calc(b)
    # A zero-flow bin uses the stream-relative floor.
    b.flow_mass_sized_comp_size["OreA", 0].fix(0.0)
    SolidPSDScaler().variable_scaling_routine(b)
    total = 1.9 + 0.5 + 2.5  # sized solids (2.0 - 0.1) + unsized + liquid
    floor = 1e-8 * total
    assert get_scaling_factor(b.flow_mass_sized_comp_size["OreA", 0]) == pytest.approx(
        1.0 / floor
    )
    assert get_scaling_factor(b.flow_mass_sized_comp_size["OreB", 0]) == pytest.approx(
        1.0 / 0.4
    )
    assert get_scaling_factor(b.flow_mass_unsized_comp["Inert"]) == pytest.approx(
        1.0 / 0.5
    )
    assert get_scaling_factor(b.flow_mass_liquid_comp["H2O"]) == pytest.approx(
        1.0 / 2.5
    )
    assert get_scaling_factor(b.temperature) == pytest.approx(1.0 / 300.0)
    assert get_scaling_factor(b.pressure) == pytest.approx(1.0e-5)
    for var in b.define_state_vars().values():
        for idx in var:
            assert get_scaling_factor(var[idx]) is not None


@pytest.mark.unit
def test_factor_sources_fixed_block_agree_unfixed_block_differ():
    # Fixed flows give both factor sources the same reference values.
    ma, ba = _state()
    mb, bb = _state()
    _fix_hand_calc(ba)
    _fix_hand_calc(bb)
    SolidPSDScaler(factor_source="current_values").variable_scaling_routine(ba)
    SolidPSDScaler(factor_source="input_based").variable_scaling_routine(bb)
    for name, va in ba.define_state_vars().items():
        vb = getattr(bb, name)
        for idx in va:
            assert get_scaling_factor(vb[idx]) == pytest.approx(
                get_scaling_factor(va[idx]), rel=1e-12
            )
    # Unfixed flows use the input-based default or their current magnitude;
    # temperature and pressure use the same fixed nominals in both modes.
    for factor_source, expected in (("input_based", 1.0), ("current_values", 4.0)):
        m = _params()
        m.state = m.params.build_state_block([0], defined_state=False)
        b = m.state[0]
        for var in (
            b.flow_mass_sized_comp_size,
            b.flow_mass_unsized_comp,
            b.flow_mass_liquid_comp,
        ):
            for idx in var:
                var[idx].set_value(0.25)
        SolidPSDScaler(factor_source=factor_source).variable_scaling_routine(b)
        for var in (
            b.flow_mass_sized_comp_size,
            b.flow_mass_unsized_comp,
            b.flow_mass_liquid_comp,
        ):
            for idx in var:
                assert get_scaling_factor(var[idx]) == pytest.approx(
                    expected
                ), f"{factor_source}: {var[idx].name}"
        assert get_scaling_factor(b.temperature) == pytest.approx(1.0 / 300.0)
        assert get_scaling_factor(b.pressure) == pytest.approx(1.0e-5)
