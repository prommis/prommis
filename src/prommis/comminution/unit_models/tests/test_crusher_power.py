#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Test crusher power laws, zero-power behavior, and pendulum duty."""

from io import StringIO

import pytest

from pyomo.environ import check_optimal_termination, value

from idaes.core.initialization import InitializationStatus
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom

from prommis.comminution.functions import distributions
from prommis.comminution.functions.distributions import ggs_cdf
from prommis.comminution.functions.power_laws import numeric_specific_energy
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    FEED_REF,
    GGS,
    PENDULUM_APPEARANCE,
    PENDULUM_ECS_A,
    PENDULUM_REF,
    RR,
    SINGLE_STRESS_REF,
    crusher,
    fix_feed,
    flowsheet,
    one_mineral_flowsheet,
)

LAW_MESHES = {
    "rittinger": ([1e-6, 2e-6, 4e-6, 8e-6, 16e-6, 32e-6], 10e-6),
    "kick": ([0.05, 0.1, 0.2, 0.4, 0.8], 0.15),
}


def _law_unit(law, power_kw=None):
    """Build a size- or power-specified crusher and return its model and unit."""

    edges, d80 = LAW_MESHES[law]
    n_int = len(edges) - 1
    m = one_mineral_flowsheet(edges=edges)
    cfg = dict(
        psd_method="distribution_function",
        shape_exponent=1.5,
        power_law=law,
        bond_work_index=12.0,
    )
    if power_kw is None:
        cfg["target_product_p80"] = d80
    cr = crusher(m, **cfg)
    rows = [0.0] * (n_int - 1) + [10.0]
    fix_feed(cr.properties_in[0], {"OreA": rows})
    if power_kw is not None:
        cr.power[0].fix(power_kw)
    # Keep the parent model alive while using the crusher.
    return m, cr


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.parametrize("law", sorted(LAW_MESHES))
def test_power_law_power_spec_round_trip(law):
    _, d80 = LAW_MESHES[law]
    m_size, size_spec = _law_unit(law)
    native.CrusherSolidPSDInitializer().initialize(size_spec)
    # For the 10 kg/s sized feed, 3.6 converts kg/s * kWh/t to kW.
    assert value(size_spec.power[0]) == pytest.approx(
        10.0
        * numeric_specific_energy(
            law,
            12.0,
            value(size_spec.product_p80[0]),
            value(size_spec.properties_in[0].percentile_size[0.8]),
        )
        * 3.6,
        rel=1e-8,
    )
    power_target = value(size_spec.power[0])
    m_power, power_spec = _law_unit(law, power_kw=power_target)
    native.CrusherSolidPSDInitializer().initialize(power_spec)
    assert value(power_spec.product_p80[0]) == pytest.approx(d80, rel=1e-8)


@pytest.mark.component
@pytest.mark.solver
def test_power_excludes_unsized_solids():
    m_control = flowsheet()
    control = crusher(m_control, bond_work_index=14.0, **SINGLE_STRESS_REF)
    fix_feed(control.properties_in[0], FEED_REF)
    native.CrusherSolidPSDInitializer().initialize(control)
    m = flowsheet(
        unsized=["Inert"], density={"OreA": 2800.0, "OreB": 4200.0, "Inert": 2000.0}
    )
    cr = crusher(m, bond_work_index=14.0, **SINGLE_STRESS_REF)
    fix_feed(cr.properties_in[0], FEED_REF, unsized={"Inert": 10.0})
    native.CrusherSolidPSDInitializer().initialize(cr)
    assert value(cr.power[0]) > 0
    assert value(cr.power[0]) == pytest.approx(value(control.power[0]), rel=1e-8)
    # Only the 10 kg/s sized feed contributes to power; total solid flow is 20 kg/s.
    f80 = value(cr.properties_in[0].percentile_size[0.8])
    p80 = value(cr.product_p80[0])
    e_spec = numeric_specific_energy(cr.config.power_law, 14.0, p80, f80)
    assert value(cr.power[0]) == pytest.approx(10.0 * e_spec * 3.6, rel=1e-8)


# Feed and breakage data for the Bond zero-power cases.
SIGN_EDGES = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]
SIGN_FEED = {"OreA": [0.0, 79.98, 20.02, 0.0]}


@pytest.mark.component
@pytest.mark.solver
def test_power_is_zero_for_nonpositive_bond_response():
    # Cumulative dominance can admit a product whose smoothed P80 is larger
    # than the feed's. Both that case and an unchanged PSD require zero duty.
    for table, feed, raw_sign in (
        ([0.69983, 0.09997, 0.20020, 0.0], SIGN_FEED, -1),
        ([1.0, 0.0, 0.0, 0.0], {"OreA": [10.0, 0.0, 0.0, 0.0]}, 0),
    ):
        m = one_mineral_flowsheet(edges=SIGN_EDGES)
        cr = crusher(
            m,
            psd_method="tabular",
            tabular_psd=table,
            power_law="bond",
            bond_work_index=14.0,
        )
        fix_feed(cr.properties_in[0], feed)
        m_s = sum(feed["OreA"])
        assert degrees_of_freedom(m) == 0
        init = native.CrusherSolidPSDInitializer()
        assert init.initialize(cr) == InitializationStatus.Ok
        for k in range(4):
            assert value(
                cr.properties_out[0].flow_mass_sized_comp_size["OreA", k]
            ) == pytest.approx(m_s * table[k], rel=1e-8, abs=1e-12)
        f80 = value(cr.properties_in[0].percentile_size[0.8])
        p80 = value(cr.product_p80[0])
        if raw_sign < 0:
            assert p80 > f80
        else:
            assert p80 == pytest.approx(f80, rel=1e-8)
        assert value(cr.power[0]) == pytest.approx(0.0, abs=1e-8)


@pytest.mark.component
@pytest.mark.solver
def test_distribution_size_spec_admitted_by_either_route():
    # A distribution target passes when either cumulative passing improves
    # or its analytic P80 is below the feed F80.
    # The GGS target P80 is 12.8 mm, above the feed F80. Moving
    # 1e-4 kg/s from the finest to the coarsest feed interval makes the
    # product finer by cumulative passing.
    m = one_mineral_flowsheet(edges=SIGN_EDGES)
    cr = crusher(
        m,
        psd_method="distribution_function",
        distribution_shape=GGS,
        shape_exponent=1.0,
        target_product_p80=0.0128,
        power_law="bond",
        bond_work_index=14.0,
    )
    shape_size = distributions.shape_size_from_p80(GGS, 0.0128, 1.0)
    assert shape_size == pytest.approx(0.016, rel=1e-12, abs=0)
    cdf = [ggs_cdf(e, shape_size, 1.0) for e in SIGN_EDGES]
    q_top = cdf[-1]
    fracs = [cdf[1] / q_top] + [(cdf[k + 1] - cdf[k]) / q_top for k in range(1, 4)]
    rows = [0.0] + [100.0 * f for f in fracs[1:]]
    rows[0] = 100.0 - sum(rows[1:])
    rows[0] -= 0.0001
    rows[3] += 0.0001
    fix_feed(cr.properties_in[0], {"OreA": rows})
    f80_live = value(cr.properties_in[0].percentile_size[0.8])
    assert f80_live < 0.0128
    init = native.CrusherSolidPSDInitializer()
    assert degrees_of_freedom(m) == 0
    assert init.initialize(cr) == InitializationStatus.Ok
    assert abs(value(cr.power[0])) <= 1e-8
    # For n=1 and d_max=16 mm, ideal GGS product fractions are
    # [1/8, 1/8, 1/4, 1/2].
    total = value(cr.properties_out[0].flow_mass_sized)
    for k, fraction in enumerate([0.125, 0.125, 0.25, 0.5]):
        assert value(
            cr.properties_out[0].flow_mass_sized_comp_size["OreA", k]
        ) / total == pytest.approx(fraction, abs=1e-4)
    # A finest-only feed makes the public cumulative product coarser.
    for k, v in enumerate([100.0, 0.0, 0.0, 0.0]):
        cr.properties_in[0].flow_mass_sized_comp_size["OreA", k].fix(v)
    feed = cr.properties_in[0].flow_mass_sized_comp_size
    product = cr.properties_out[0].flow_mass_sized_comp_size
    assert (
        min(
            sum(
                value(product["OreA", j]) - value(feed["OreA", j]) for j in range(k + 1)
            )
            for k in range(3)
        )
        < 0
    )
    with pytest.raises(
        ConfigurationError,
        match=(
            "no-coarsening.*cumulative product fraction"
            ".*below cumulative feed fraction"
        ),
    ):
        init.initialize(cr)
    # This feed has 3 kg/s more in the finest bin than the product, so
    # the no-coarsening constraint prevents a direct solve.
    for k, v in enumerate([rows[0] + 3.0001, rows[1], rows[2], rows[3] - 3.0001]):
        cr.properties_in[0].flow_mass_sized_comp_size["OreA", k].fix(v)
    assert not check_optimal_termination(get_solver("ipopt_v2").solve(m))
    # The Rosin-Rammler product fails the cumulative check, but its
    # 5 mm analytic P80 is below the feed F80, so it passes with positive duty.
    m_p80 = one_mineral_flowsheet(edges=SIGN_EDGES)
    cr_p80 = crusher(
        m_p80,
        psd_method="distribution_function",
        distribution_shape=RR,
        shape_exponent=1.5,
        target_product_p80=5.0e-3,
        bond_work_index=14.0,
    )
    fix_feed(cr_p80.properties_in[0], {"OreA": [50.0, 0.0, 0.0, 50.0]})
    f80_p80 = value(cr_p80.properties_in[0].percentile_size[0.8])
    assert f80_p80 > 5.0e-3
    init_p80 = native.CrusherSolidPSDInitializer()
    assert degrees_of_freedom(m_p80) == 0
    assert init_p80.initialize(cr_p80) == InitializationStatus.Ok
    duty = 3.6 * 14.0 * numeric_specific_energy("bond", 1.0, 5.0e-3, f80_p80) * 100.0
    assert value(cr_p80.power[0]) == pytest.approx(duty, rel=1e-8)
    row0 = value(cr_p80.properties_out[0].flow_mass_sized_comp_size["OreA", 0])
    # 100 F(2 mm)/F(16 mm), with F(d) = 1 - exp(-(d/d_63)**1.5) for the 5 mm target.
    assert row0 == pytest.approx(33.44951, abs=1e-5)


# Pendulum power configuration and solved duty
_ECS_B = {
    "sizes": [2e-3, 8e-3],
    "t10_rows": [10, 20, 30],
    "ecs": [[0.6, 0.5], [1.2, 1.0], [1.6, 1.4]],
}

# Pendulum power P = A*P_p + P_0 requires Whiten recycle, t10 breakage,
# and a nonnegative Ecs table (Napier-Munn et al., 1996, ch. 11).
# Each case is (ID, crusher options, expected error).
PENDULUM_CONFIG_REJECTIONS = [
    (
        "single_pass",
        dict(PENDULUM_REF, recycle="single_pass"),
        "requires psd_method.*recycle='classification_recycle'",
    ),
    (
        "no_power_params",
        dict(PENDULUM_REF, pendulum_power_params=None),
        "requires pendulum_power_params",
    ),
    (
        "work_index_set",
        dict(PENDULUM_REF, bond_work_index=14.0),
        "does not use bond_work_index",
    ),
    (
        "no_ecs_table",
        dict(
            PENDULUM_REF,
            breakage_params={"appearance_table": PENDULUM_APPEARANCE, "t10": 15.0},
        ),
        r"missing \['ecs_table'\]",
    ),
    (
        "ecs_table_not_a_mapping",
        dict(
            PENDULUM_REF,
            breakage_params=dict(PENDULUM_REF["breakage_params"], ecs_table=[1]),
        ),
        "breakage_params: ecs_table must be a mapping",
    ),
    (
        "ecs_sizes_string",
        dict(
            PENDULUM_REF,
            breakage_params=dict(
                PENDULUM_REF["breakage_params"],
                ecs_table=dict(PENDULUM_ECS_A, sizes=["0.002", "0.008"]),
            ),
        ),
        r"ecs_table\['sizes'\]\[0\] must be a real number \(got str\)",
    ),
    (
        "negative_ecs",
        dict(
            PENDULUM_REF,
            breakage_params=dict(
                PENDULUM_REF["breakage_params"],
                ecs_table=dict(
                    PENDULUM_ECS_A, ecs=[[0.4, -0.3], [0.8, 0.6], [1.1, 0.9]]
                ),
            ),
        ),
        "breakage_params: ecs_table 'ecs' entries must be finite and nonnegative",
    ),
    (
        "pendulum_params_under_bond",
        dict(PENDULUM_REF, power_law="bond"),
        "consumed only by power_law='pendulum'",
    ),
    (
        "breakage_from_jaw_preset",
        dict(
            {k: v for k, v in PENDULUM_REF.items() if k != "breakage_function"},
            crusher_equipment="jaw",
        ),
        "breakage_function='vogel' came from crusher_equipment='jaw'",
    ),
]


@pytest.mark.unit
def test_pendulum_config_rejected():
    failures = []
    for case_id, kwargs, match in PENDULUM_CONFIG_REJECTIONS:
        try:
            with pytest.raises(ConfigurationError, match=match):
                crusher(flowsheet(), **kwargs)
        except (Exception, pytest.fail.Exception) as exc:
            failures.append(f"{case_id}: {exc}")
    assert not failures, "\n".join(failures)


@pytest.mark.component
@pytest.mark.solver
def test_pendulum_power_closes_on_solved_loads():
    m = flowsheet()
    by_species = {
        "OreA": PENDULUM_REF["breakage_params"],
        "OreB": dict(PENDULUM_REF["breakage_params"], ecs_table=_ECS_B),
    }
    cr = crusher(
        m,
        **dict(PENDULUM_REF, breakage_params=None, breakage_params_by_comp=by_species),
    )
    fix_feed(cr.properties_in[0], FEED_REF)
    sset = list(m.fs.pp.size_interval_set)
    native.CrusherSolidPSDInitializer().initialize(cr)
    broken = 3.6 * sum(  # 3.6 converts kg/s * kWh/t to kW.
        value(cr.ecs[s, k])
        * value(cr.classification[k])
        * value(cr.recycle_load[0, s, k])
        for s in ("OreA", "OreB")
        for k in sset
    )
    assert value(cr.ecs["OreA", 1]) != value(cr.ecs["OreB", 1])
    # At t10=15, OreA's Ecs values at 2 and 8 mm are 0.6075 and
    # 0.45 kWh/t. Interpolate in log(size) at size_char=2**(k+0.5) mm.
    assert [value(cr.ecs["OreA", k]) for k in sset] == pytest.approx(
        [0.64688, 0.56813, 0.48938, 0.41063], abs=1e-5
    )
    assert value(cr.power[0]) == pytest.approx(20.0 + 1.3 * broken, rel=1e-8)
    # Fix power and solve for the released pendulum power factor.
    cr.power[0].fix()
    cr.pendulum_power_factor.unfix()
    cr.pendulum_power_factor.set_value(2.0)
    assert degrees_of_freedom(cr) == 0
    assert check_optimal_termination(get_solver("ipopt_v2").solve(cr))
    assert value(cr.pendulum_power_factor) == pytest.approx(1.3, rel=1e-6)
    report = StringIO()
    cr.report(ostream=report)
    assert "Pendulum power" in report.getvalue()
