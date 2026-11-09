#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Test crusher operator arithmetic and component flows."""

import pytest

from pyomo.environ import units, value

from idaes.core.initialization import InitializationStatus

from prommis.comminution.functions.selection import selection_value
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    COMMON_B,
    COMMON_S,
    WHITEN_REF,
    crusher,
    fix_feed,
    flowsheet,
    recycle_reference,
    stress_reference,
)

# Different feed shapes expose an incorrect aggregate product split by
# mass share, even when the total product mass is correct.
COMMON_RECYCLE_B = [
    [1.0, 0.6, 0.6, 0.6],
    [0.0, 0.4, 0.0, 0.0],
    [0.0, 0.0, 0.4, 0.0],
    [0.0, 0.0, 0.0, 0.4],
]
COMMON_FEED = {"OreA": [0.0, 0.5, 1.5, 4.0], "OreB": [0.2, 0.0, 0.8, 3.0]}
COMMON_EVENTS = {"single_stress": 1, "repeated_stress": 3}


def common_operator_unit(path):
    """Build a crusher for ``path`` and fix the shared feed."""
    m = flowsheet()
    if path == "classification_recycle":
        cr = crusher(
            m,
            psd_method="selection_breakage",
            recycle="classification_recycle",
            selection_function="whiten",
            selection_params=WHITEN_REF,
            breakage_function="user",
            breakage_params=COMMON_RECYCLE_B,
        )
    else:
        cr = crusher(
            m,
            psd_method="selection_breakage",
            recycle="single_pass",
            selection_function="user",
            selection_params=COMMON_S,
            breakage_function="user",
            breakage_params=COMMON_B,
            n_stress_events=COMMON_EVENTS[path],
        )
    fix_feed(cr.properties_in[0], COMMON_FEED)
    return m, cr


def common_whiten_classification(cr):
    """Return Whiten classification fractions from the shared selection inputs."""
    pp = cr.config.property_package
    return [
        selection_value(
            "whiten",
            value(units.convert(pp.size_char[k], to_units=units.m)),
            WHITEN_REF,
        )
        for k in pp.size_interval_set
    ]


def common_operator_reference(cr, path, feed):
    """Return reference working rows and product flows for one feed.

    Working rows are indexed by variable name and stage index, if any.
    """
    working = {}
    if path == "classification_recycle":
        load, product = recycle_reference(
            COMMON_RECYCLE_B, common_whiten_classification(cr), feed
        )
        working["recycle_load"] = {(): load}
        return working, product
    events = COMMON_EVENTS[path]
    if events > 1:
        working["psd_stage"] = {
            (r,): stress_reference(COMMON_S, COMMON_B, feed, r)
            for r in range(1, events)
        }
    return working, stress_reference(COMMON_S, COMMON_B, feed, events)


@pytest.mark.component
@pytest.mark.solver
@pytest.mark.parametrize(
    "path", ["single_stress", "repeated_stress", "classification_recycle"]
)
def test_common_operator_keeps_minerals_separate(path):
    m, cr = common_operator_unit(path)
    init = native.CrusherSolidPSDInitializer()
    assert init.initialize(cr) == InitializationStatus.Ok
    order = list(cr.config.property_package.size_interval_set)
    references = {}
    for s, feed in COMMON_FEED.items():
        working, product = common_operator_reference(cr, path, feed)
        references[s] = product
        targets = [
            (cr.properties_out[0].flow_mass_sized_comp_size[s, k], product[k])
            for k in order
        ]
        for name, rows in working.items():
            component = getattr(cr, name)
            targets += [
                (component[(0, *index, s, k)], row[k])
                for index, row in rows.items()
                for k in order
            ]
        for target, expected in targets:
            solved = value(units.convert(target, to_units=units.kg / units.s))
            assert solved == pytest.approx(expected, rel=1e-8, abs=1e-12)
    aggregate = [sum(feed[k] for feed in COMMON_FEED.values()) for k in order]
    _, blind = common_operator_reference(cr, path, aggregate)
    for s, feed in COMMON_FEED.items():
        share = sum(feed) / sum(aggregate)
        assert max(abs(share * blind[k] - references[s][k]) for k in order) > 0.1
