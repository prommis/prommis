#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Crusher selection and breakage tests, including per-component product rows."""

import math

import numpy
import pytest

from pyomo.environ import value

from idaes.core.initialization import InitializationStatus
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom

from prommis.comminution.functions.breakage import breakage_matrix
from prommis.comminution.functions.selection import selection_value
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    B_A,
    B_B,
    EDGES_X,
    FEED_X,
    HC_EDGES,
    LUCKIE_AUSTIN_REF,
    S_A,
    S_B,
    ZERO_BOTTOM,
    ZERO_EDGES,
    ZERO_VOGEL_BY_COMP,
    ZERO_WHITEN,
    ZERO_WHITEN_BY_COMP,
    crusher,
    fix_feed,
    flowsheet,
    one_mineral_flowsheet,
    x_individual,
)

# Hand-calculated single-pass selection and breakage data for HC_EDGES.
HC_B = [[1.0, 1.0, 0.4], [0.0, 0.0, 0.6], [0.0, 0.0, 0.0]]


# Build-time parameter cases on HC_EDGES: (case ID, crusher options, error match).
_VOGEL_B = {"q": 0.05, "dprime": 0.005}
_RULE_WHITEN = {"K1": 2e-3, "K2": 8e-3}
_RULE_KING = {"CSS": 0.025, "alpha1": 0.5, "alpha2": 1.2, "n": 1.0}
_RULE_AUSTIN = {"S1": 0.5, "d1": 4e-3, "a": 0.5}
_RULE_VP = {"f_mat": 1.0, "xw_min": 0.0, "v": 10.0}
_RULE_APPEARANCE = {
    "10": {"t75": 1.0, "t50": 1.5, "t25": 3.0, "t4": 20.0, "t2": 40.0},
    "20": {"t75": 2.0, "t50": 3.0, "t25": 6.0, "t4": 35.0, "t2": 60.0},
}
_RULE_ECS = {"sizes": [2e-3, 8e-3], "t10_rows": [10, 20], "ecs": [[0.4, 0.3]] * 2}


def _single_pass(sel_fn, sel, brk_fn="vogel", brk=None):
    """Return single-pass crusher options, with Vogel breakage by default."""
    return dict(
        psd_method="selection_breakage",
        recycle="single_pass",
        selection_function=sel_fn,
        selection_params=sel,
        breakage_function=brk_fn,
        breakage_params=_VOGEL_B if brk is None else brk,
    )


def _vogel_peukert(params):
    """Return options with Vogel-Peukert selection and Luckie-Austin breakage."""
    return _single_pass("vogel_peukert", params, "luckie_austin", LUCKIE_AUSTIN_REF)


def _t10_recycle(table, t10=15.0, ecs_table=None):
    """Return Whiten-recycle options with t10 appearance breakage.

    Use the pendulum power law when ecs_table is provided.
    """
    kwargs = dict(
        psd_method="selection_breakage",
        recycle="classification_recycle",
        selection_function="whiten",
        selection_params={"K1": 2e-3, "K2": 8e-3, "K3": 2.3},
        breakage_function="t10_appearance",
        breakage_params={"appearance_table": table, "t10": t10},
    )
    if ecs_table is not None:
        kwargs["breakage_params"]["ecs_table"] = ecs_table
        kwargs["power_law"] = "pendulum"
        kwargs["pendulum_power_params"] = {"power_factor": 1.3, "no_load_power": 20.0}
    return kwargs


def _missing(table, row, marker):
    """Return a copy of ``table`` without ``marker`` in ``row``."""
    edited = {key: dict(markers) for key, markers in table.items()}
    del edited[row][marker]
    return edited


def _assert_build_rules(rejections, controls):
    """Build each case on a fresh model; collect every rule that fails."""
    failures = []
    for case_id, kwargs, match in rejections:
        m = one_mineral_flowsheet(edges=HC_EDGES)
        try:
            with pytest.raises(ConfigurationError, match=match):
                crusher(m, **kwargs)
        except (Exception, pytest.fail.Exception) as exc:
            failures.append(f"{case_id}: {exc}")
    for case_id, kwargs in controls:
        m = one_mineral_flowsheet(edges=HC_EDGES)
        try:
            crusher(m, **kwargs)
        except (Exception, pytest.fail.Exception) as exc:
            failures.append(f"{case_id}: {exc}")
    assert not failures, "\n".join(failures)


SELECTION_REJECTIONS = [
    (
        "whiten_K1_zero",
        _single_pass("whiten", dict(_RULE_WHITEN, K1=0.0)),
        r"'K[12]' in selection_params must be > 0",
    ),
    (
        "whiten_K2_negative",
        _single_pass("whiten", dict(_RULE_WHITEN, K2=-1e-3)),
        r"'K[12]' in selection_params must be > 0",
    ),
    (
        "whiten_K3_zero",
        _single_pass("whiten", dict(_RULE_WHITEN, K3=0.0)),
        "'K3' in selection_params must be > 0",
    ),
    (
        "whiten_K2_below_1.01_K1",
        _single_pass("whiten", {"K1": 0.2, "K2": 0.2}),
        "K2 >= 1.01",
    ),
    (
        "whiten_K2_below_1.01_K1_recycle",
        dict(
            psd_method="selection_breakage",
            recycle="classification_recycle",
            selection_function="whiten",
            selection_params={"K1": 4e-3, "K2": 4e-3, "K3": 2.3},
            breakage_function="user",
            breakage_params=[[1.0, 0.0, 0.4], [0.0, 1.0, 0.6], [0.0, 0.0, 0.0]],
        ),
        "K2 >= 1.01",
    ),
    (
        "whiten_string",
        _single_pass("whiten", {"K1": "0.002", "K2": "0.008"}),
        "not a string",
    ),
    (
        "whiten_boolean",
        _single_pass("whiten", {"K1": 2e-3, "K2": True}),
        "boolean",
    ),
    (
        "whiten_stray_k3",
        _single_pass("whiten", {"K1": 2e-3, "K2": 8e-3, "k3": 2.0}),
        "unknown selection_params",
    ),
    (
        "whiten_params_none",
        _single_pass("whiten", None),
        r"selection_params with keys \['K1', 'K2'\]; missing \['K1', 'K2'\]",
    ),
    (
        "king_CSS_zero",
        _single_pass("king", dict(_RULE_KING, CSS=0.0)),
        "'CSS' in selection_params must be > 0",
    ),
    (
        "king_alpha_order",
        _single_pass("king", dict(_RULE_KING, alpha1=1.2, alpha2=0.5)),
        "alpha1 < alpha2",
    ),
    (
        "king_n_zero",
        _single_pass("king", dict(_RULE_KING, n=0.0)),
        "'n' in selection_params must be > 0",
    ),
    (
        "king_missing_alpha",
        _single_pass("king", {"CSS": 0.025, "n": 1.0}),
        r"selection_params with keys \['CSS', 'alpha1', 'alpha2', 'n'\]; "
        r"missing \['alpha1', 'alpha2'\]",
    ),
    (
        "king_stray_m",
        _single_pass("king", dict(_RULE_KING, m=2.0)),
        "unknown selection_params",
    ),
    (
        "austin_S1_negative",
        _single_pass("austin", dict(_RULE_AUSTIN, S1=-0.1)),
        "'S1' in selection_params must be [<>]=",
    ),
    (
        "austin_S1_above_one",
        _single_pass("austin", dict(_RULE_AUSTIN, S1=1.1)),
        "'S1' in selection_params must be [<>]=",
    ),
    (
        "austin_a_negative",
        _single_pass("austin", dict(_RULE_AUSTIN, a=-0.1)),
        "'a' in selection_params must be [<>]=",
    ),
    (
        "austin_a_above_one",
        _single_pass("austin", dict(_RULE_AUSTIN, a=1.5)),
        "'a' in selection_params must be [<>]=",
    ),
    (
        "austin_d1_zero",
        _single_pass("austin", {"S1": 0.5, "d1": 0.0, "a": 0.5}),
        "'d1' in selection_params must be > 0",
    ),
    *(
        (
            f"vogel_peukert_missing_{key}",
            _vogel_peukert({k: v for k, v in _RULE_VP.items() if k != key}),
            "f_mat.*xw_min.*v",
        )
        for key in _RULE_VP
    ),
    (
        "vogel_peukert_stray_key",
        _vogel_peukert(dict(_RULE_VP, k=1.0)),
        "allowed keys",
    ),
    ("vogel_peukert_f_mat_zero", _vogel_peukert(dict(_RULE_VP, f_mat=0.0)), "f_mat"),
    (
        # Any negative xw_min is rejected.
        "vogel_peukert_xw_min_negative",
        _vogel_peukert(dict(_RULE_VP, xw_min=-1e-12)),
        "xw_min",
    ),
    (
        "vogel_peukert_v_zero",
        _single_pass("vogel_peukert", dict(_RULE_VP, xw_min=0.1, v=0.0)),
        "'v' in selection_params must be > 0",
    ),
    (
        "whiten_by_species_missing_K1",
        dict(
            psd_method="selection_breakage",
            recycle="single_pass",
            selection_function="whiten",
            selection_params_by_comp={"OreA": {"K2": 8e-3}},
            breakage_function="vogel",
            breakage_params=_VOGEL_B,
        ),
        r"selection_params_by_comp\['OreA'\] with keys \['K1', 'K2'\]; "
        r"missing \['K1'\]",
    ),
]

SELECTION_CONTROLS = [
    ("whiten", _single_pass("whiten", dict(_RULE_WHITEN, K3=2.0))),
    ("king", _single_pass("king", _RULE_KING)),
    ("austin", _single_pass("austin", _RULE_AUSTIN)),
    # Any positive f_mat and v is accepted, and xw_min may be zero.
    ("vogel_peukert_f_mat_tiny", _vogel_peukert(dict(_RULE_VP, f_mat=1e-12))),
    ("vogel_peukert_xw_min_zero", _vogel_peukert(dict(_RULE_VP, xw_min=0.0))),
    ("vogel_peukert_v_tiny", _vogel_peukert(dict(_RULE_VP, v=1e-12))),
]


@pytest.mark.unit
def test_selection_params_rejected_at_build():
    _assert_build_rules(SELECTION_REJECTIONS, SELECTION_CONTROLS)


BREAKAGE_REJECTIONS = [
    (
        "luckie_austin_phi_negative",
        _single_pass(
            "whiten", _RULE_WHITEN, "luckie_austin", dict(LUCKIE_AUSTIN_REF, phi=-0.1)
        ),
        "'phi' in breakage_params must be [<>]=",
    ),
    (
        "luckie_austin_phi_above_one",
        _single_pass(
            "whiten", _RULE_WHITEN, "luckie_austin", dict(LUCKIE_AUSTIN_REF, phi=1.1)
        ),
        "'phi' in breakage_params must be [<>]=",
    ),
    (
        "luckie_austin_gamma_negative",
        _single_pass(
            "user",
            [0.0, 0.5, 1.0],
            "luckie_austin",
            {"phi": 0.5, "gamma": -1.0, "beta": 3.0},
        ),
        "gamma",
    ),
    (
        "reid_stewart_phi_negative",
        _single_pass(
            "whiten", _RULE_WHITEN, "reid_stewart", dict(LUCKIE_AUSTIN_REF, phi=-0.1)
        ),
        "'phi' in breakage_params must be [<>]=",
    ),
    (
        "reid_stewart_phi_above_one",
        _single_pass(
            "whiten", _RULE_WHITEN, "reid_stewart", dict(LUCKIE_AUSTIN_REF, phi=1.1)
        ),
        "'phi' in breakage_params must be [<>]=",
    ),
    (
        "vogel_dprime_zero",
        _single_pass("user", [0.0, 0.5, 1.0], "vogel", {"q": 0.05, "dprime": 0.0}),
        "dprime",
    ),
    (
        "vogel_dprime_string",
        _single_pass("user", [0.0, 0.5, 1.0], "vogel", {"q": 0.05, "dprime": "0.005"}),
        "not a string",
    ),
    (
        "vogel_stray_Q",
        _single_pass("user", [0.0, 0.5, 1.0], "vogel", dict(_VOGEL_B, Q=0.1)),
        "unknown breakage_params",
    ),
    (
        "vogel_missing_q",
        _single_pass("user", [0.0, 0.5, 1.0], "vogel", {"dprime": 0.005}),
        r"breakage_params with keys \['q', 'dprime'\]; missing \['q'\]",
    ),
    (
        "t10_appearance_single_pass",
        _single_pass(
            "whiten",
            _RULE_WHITEN,
            "t10_appearance",
            {"appearance_table": _RULE_APPEARANCE, "t10": 15.0},
        ),
        "valid only in the Whiten",
    ),
    (
        "appearance_table_list",
        _t10_recycle(list(_RULE_APPEARANCE.values())),
        "appearance_table must be a dict",
    ),
    (
        "t10_below_range",
        _t10_recycle(_RULE_APPEARANCE, 5.0),
        "outside the appearance-table rows",
    ),
    ("t10_nan", _t10_recycle(_RULE_APPEARANCE, float("nan")), "t10 must be a finite"),
    (
        "missing_t2",
        _t10_recycle(_missing(_RULE_APPEARANCE, "10", "t2")),
        "missing marker",
    ),
    (
        "row_none",
        _t10_recycle({"10": None, "20": _RULE_APPEARANCE["20"]}),
        "must be a dict",
    ),
    (
        "duplicate_key",
        _t10_recycle(
            {
                **_RULE_APPEARANCE,
                10: {"t75": 9.0, "t50": 9.0, "t25": 9.0, "t4": 9.0, "t2": 9.0},
            }
        ),
        "duplicate",
    ),
    (
        "non_integer_key",
        _t10_recycle({10.9: _RULE_APPEARANCE["10"], "20": _RULE_APPEARANCE["20"]}),
        "integer-like",
    ),
    (
        "unknown_marker",
        _t10_recycle(
            {
                "10": dict(_RULE_APPEARANCE["10"], t20=10.0),
                "20": _RULE_APPEARANCE["20"],
            }
        ),
        r"unknown markers \['t20'\]",
    ),
    (
        "ecs_empty",
        _t10_recycle(_RULE_APPEARANCE, ecs_table=dict(_RULE_ECS, ecs=[])),
        "must be nonempty lists or tuples",
    ),
    (
        "ecs_row_width",
        _t10_recycle(_RULE_APPEARANCE, ecs_table=dict(_RULE_ECS, ecs=[[0.4], [0.8]])),
        "one value per size",
    ),
    (
        "ecs_row_count",
        _t10_recycle(_RULE_APPEARANCE, ecs_table=dict(_RULE_ECS, ecs=[[0.4, 0.3]])),
        "one 'ecs' row per t10 row",
    ),
    (
        "marker_string",
        _t10_recycle(
            {
                "10": dict(_RULE_APPEARANCE["10"], t2="40.0"),
                "20": _RULE_APPEARANCE["20"],
            }
        ),
        "must be a real number",
    ),
    (
        "ecs_cell_numpy_bool",
        _t10_recycle(
            _RULE_APPEARANCE,
            ecs_table=dict(_RULE_ECS, ecs=[[numpy.bool_(True), 0.3], [0.4, 0.3]]),
        ),
        "must be a real number",
    ),
]

BREAKAGE_CONTROLS = [
    ("vogel", _single_pass("user", [0.0, 0.5, 1.0], "vogel", _VOGEL_B)),
    (
        "luckie_austin",
        _single_pass("whiten", _RULE_WHITEN, "luckie_austin", LUCKIE_AUSTIN_REF),
    ),
    ("t10_appearance_recycle", _t10_recycle(_RULE_APPEARANCE)),
]


@pytest.mark.unit
def test_breakage_params_rejected_at_build():
    _assert_build_rules(BREAKAGE_REJECTIONS, BREAKAGE_CONTROLS)


# Compare generated operators with the named kernels at characteristic sizes.
_ANALYTIC_SELECTION = {
    "king": {"CSS": 4.0e-3, "alpha1": 0.5, "alpha2": 1.2, "n": 2.0},
    "austin": {"S1": 0.5, "d1": 4.0e-3, "a": 0.5},
    "vogel_peukert": {"f_mat": 1.0, "xw_min": 0.1, "v": 10.0},
}

# Hand values of each selection row at size_char = 2**(k + 0.5) mm on HC_EDGES.
_ANALYTIC_SELECTION_ROWS = {
    "king": [0.0, 1.0 - ((4.8 - 2.0 * math.sqrt(2.0)) / 2.8) ** 2, 1.0],
    "austin": [2.0 ** ((k - 3.5) / 2.0) for k in range(3)],
    "vogel_peukert": [
        0.0,
        1.0 - math.exp(-(50.0 * 2.0**1.5 * 1e-3 - 0.1)),
        1.0 - math.exp(-(50.0 * 2.0**2.5 * 1e-3 - 0.1)),
    ],
}


def _stored_breakage(cr, s, n_int):
    return [[value(cr.breakage[s, i, j]) for j in range(n_int)] for i in range(n_int)]


def _assert_breakage_matrix_matches_expected(stored, expected):
    """Compare breakage matrices and check that stored columns sum to one."""
    for row, expected_row in zip(stored, expected, strict=True):
        assert row == pytest.approx(expected_row, rel=1e-12, abs=0)
    for column in zip(*stored):
        assert math.fsum(column) == pytest.approx(1.0, rel=0, abs=1e-15)


@pytest.mark.unit
def test_generated_operators_follow_named_kernels_at_d_char():
    models = []  # Keep each parent model alive while comparing its unit.

    def build(m, **cfg):
        models.append(m)
        return crusher(m, **cfg)

    def zero_edge_flowsheet(species=("OreA",)):
        return flowsheet(
            edges=ZERO_EDGES,
            species=list(species),
            density={s: {"OreA": 2800.0, "OreB": 4200.0}[s] for s in species},
            bottom_size=ZERO_BOTTOM,
        )

    # With a zero finest edge, selection still uses each interval's
    # characteristic size rather than zero.
    vogel = {"q": 0.5, "dprime": 2.0e-3}
    base = dict(
        psd_method="selection_breakage",
        selection_function="whiten",
        selection_params=ZERO_WHITEN,
        breakage_function="vogel",
        breakage_params=vogel,
    )
    configs = {
        "single_stress": dict(base, recycle="single_pass", n_stress_events=1),
        "repeated_stress": dict(base, recycle="single_pass", n_stress_events=3),
        "classification_recycle": dict(base, recycle="classification_recycle"),
    }
    built = {path: build(zero_edge_flowsheet(), **cfg) for path, cfg in configs.items()}
    size_char = built["single_stress"].config.property_package.size_char_m
    expected_sel = [selection_value("whiten", d, ZERO_WHITEN) for d in size_char]
    # Whiten C at the bin-0 size sqrt(0.5 x 1) mm: 1 - ((4 - 0.7071)/3.5)**2.3.
    assert expected_sel[0] == pytest.approx(0.13090, abs=1e-5)
    for path, cr in built.items():
        row = cr.classification if path == "classification_recycle" else None
        for k in range(3):
            stored = row[k] if row is not None else cr.selection["OreA", k]
            assert value(stored) == expected_sel[k]
        # Vogel B(x) = (x/d_2)**0.5 (1 + tanh((x - 2 mm)/2 mm))/2 at characteristic
        # sizes, parent d_2 = 2.828 mm: B(0.707 mm) and B(1.414 mm) - B(0.707 mm).
        assert value(cr.breakage["OreA", 0, 2]) == pytest.approx(0.10768, abs=1e-5)
        assert value(cr.breakage["OreA", 1, 2]) == pytest.approx(0.14518, abs=1e-5)

    # Compare each selection kernel with the hand-calculated values above.
    for name, params in _ANALYTIC_SELECTION.items():
        m = one_mineral_flowsheet(edges=HC_EDGES)
        cr = build(m, **_single_pass(name, params, "user", [list(r) for r in HC_B]))
        assert [value(cr.selection["OreA", k]) for k in range(3)] == pytest.approx(
            _ANALYTIC_SELECTION_ROWS[name], rel=1e-12
        )

    # Stored breakage matrices must match their kernels, each column normalized
    # to sum to one, and differ by kernel.
    kernel_params = {"phi": 0.5, "gamma": 1.5, "beta": 4.0}
    stored = {}
    for name in ("luckie_austin", "reid_stewart"):
        m = one_mineral_flowsheet(edges=HC_EDGES)
        cr = build(m, **_single_pass("whiten", _RULE_WHITEN, name, kernel_params))
        expected = breakage_matrix(name, HC_EDGES, m.fs.pp.size_char_m, kernel_params)
        stored[name] = _stored_breakage(cr, "OreA", 3)
        _assert_breakage_matrix_matches_expected(stored[name], expected)
    assert stored["luckie_austin"] != stored["reid_stewart"]

    # Each component uses its own parameters. Compare selection entries exactly.
    sel, brk = ZERO_WHITEN_BY_COMP, ZERO_VOGEL_BY_COMP
    m = zero_edge_flowsheet(("OreA", "OreB"))
    cr = build(
        m,
        **dict(
            base,
            recycle="single_pass",
            selection_params=None,
            selection_params_by_comp=sel,
        ),
    )
    rows = {s: [value(cr.selection[s, k]) for k in range(3)] for s in ("OreA", "OreB")}
    for s in rows:
        assert rows[s] == [selection_value("whiten", d, sel[s]) for d in size_char]
    assert rows["OreA"] != rows["OreB"]
    m = zero_edge_flowsheet(("OreA", "OreB"))
    cr = build(
        m,
        **dict(
            base,
            recycle="classification_recycle",
            breakage_params=None,
            breakage_params_by_comp=brk,
        ),
    )
    matrices = {s: _stored_breakage(cr, s, 3) for s in ("OreA", "OreB")}
    for s, matrix in matrices.items():
        expected = breakage_matrix("vogel", ZERO_EDGES, size_char, brk[s])
        _assert_breakage_matrix_matches_expected(matrix, expected)
    assert matrices["OreA"] != matrices["OreB"]


@pytest.mark.component
@pytest.mark.solver
def test_independent_individuality_four_combinations():
    combos = {
        ("common", "common"): {"OreA": [2, 4, 2], "OreB": [2, 4, 2]},
        ("individual", "common"): {"OreA": [2, 4, 2], "OreB": [1, 2, 5]},
        ("common", "individual"): {"OreA": [2, 4, 2], "OreB": [1, 1, 6]},
        ("individual", "individual"): {"OreA": [2, 4, 2], "OreB": [0.5, 0.5, 7]},
    }
    for (sel_form, brk_form), oracle in combos.items():
        cfg = dict(
            psd_method="selection_breakage",
            recycle="single_pass",
            n_stress_events=1,
            selection_function="user",
            breakage_function="user",
        )
        if sel_form == "common":
            cfg["selection_params"] = list(S_A)
        else:
            cfg["selection_params_by_comp"] = {"OreA": list(S_A), "OreB": list(S_B)}
        if brk_form == "common":
            cfg["breakage_params"] = [list(r) for r in B_A]
        else:
            cfg["breakage_params_by_comp"] = {
                "OreA": [list(r) for r in B_A],
                "OreB": [list(r) for r in B_B],
            }
        m = flowsheet(edges=EDGES_X)
        cr = crusher(m, bond_work_index=14.0, **cfg)
        fix_feed(cr.properties_in[0], FEED_X)
        assert degrees_of_freedom(m) == 0
        assert (
            native.CrusherSolidPSDInitializer().initialize(cr)
            == InitializationStatus.Ok
        )
        product = {
            idx: value(var)
            for idx, var in cr.properties_out[0].flow_mass_sized_comp_size.items()
        }
        for s, rows in oracle.items():
            for k, exp in enumerate(rows):
                assert product[(s, k)] == pytest.approx(exp, rel=1e-8)


@pytest.mark.component
@pytest.mark.solver
def test_per_mineral_conservation_and_separation():
    m = flowsheet(edges=EDGES_X)
    cr = crusher(m, bond_work_index=14.0, **x_individual("repeated_stress"))
    fix_feed(cr.properties_in[0], FEED_X)
    native.CrusherSolidPSDInitializer().initialize(cr)
    # These products follow two stress events; stage 1 still holds the
    # single-event result.
    for s, product, stage1 in (
        ("OreA", [7 / 2, 4, 1 / 2], [2, 4, 2]),
        ("OreB", [31 / 32, 29 / 32, 49 / 8], [1 / 2, 1 / 2, 7]),
    ):
        for k in range(3):
            assert value(
                cr.properties_out[0].flow_mass_sized_comp_size[s, k]
            ) == pytest.approx(product[k], rel=1e-8)
            assert value(cr.psd_stage[0, 1, s, k]) == pytest.approx(stage1[k], rel=1e-8)
