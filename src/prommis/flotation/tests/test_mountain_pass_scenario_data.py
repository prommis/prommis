####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Tests for Mountain Pass flotation scenario data."""

import pytest

from prommis.flotation.bastnaesite_properties import COMPONENTS
from prommis.flotation.mountain_pass_flotation_fixed_recovery_flowsheet import (
    TABLE1_PRODUCT_FIT_SCENARIO,
    TABLE2_REPORTED_SCENARIO,
    SCENARIO_DATA,
    SCENARIOS,
    stream_total,
)


@pytest.mark.unit
def test_generated_data_has_expected_scenarios():
    assert set(SCENARIO_DATA) == {"feed", "scenarios", "table_1_products"}
    assert set(SCENARIO_DATA["feed"]) == {"fresh_feed"}
    assert set(SCENARIO_DATA["feed"]["fresh_feed"]) == set(COMPONENTS)
    assert SCENARIO_DATA["table_1_products"]["final_concentrate"][
        "reo_recovery"
    ] == pytest.approx(80.1)
    assert SCENARIO_DATA["table_1_products"]["combined_tails"]["assays"][
        "REO"
    ] == pytest.approx(1.7)
    assert set(SCENARIOS) == {
        "table2_reported_with_fitted_inert",
        "figure2_mass_balance",
        "table1_product_fit",
    }
    for scenario in SCENARIOS.values():
        assert set(scenario) == {
            "description",
            "initial_streams",
            "recoveries",
        }
        assert "fresh_feed" not in scenario["initial_streams"]


@pytest.mark.unit
def test_generated_data_closes_connected_stream_mass_balance():
    streams = {
        "fresh_feed": SCENARIO_DATA["feed"]["fresh_feed"],
        **SCENARIOS[TABLE2_REPORTED_SCENARIO]["initial_streams"],
    }
    for component in COMPONENTS:
        assert streams["rougher_feed"][component] == pytest.approx(
            streams["fresh_feed"][component]
            + streams["scavenger_concentrate"][component]
        )
        assert streams["cleaner1_feed"][component] == pytest.approx(
            streams["rougher_concentrate"][component]
            + streams["cleaner2_tails"][component]
        )
        assert streams["cleaner2_feed"][component] == pytest.approx(
            streams["cleaner1_concentrate"][component]
            + streams["cleaners34_tails"][component]
        )
        assert streams["final_tails"][component] == pytest.approx(
            streams["rougher_tails"][component] + streams["scavenger_tails"][component]
        )


@pytest.mark.unit
def test_table1_product_fit_scenario_hits_product_envelope():
    scenario = SCENARIOS[TABLE1_PRODUCT_FIT_SCENARIO]
    feed = SCENARIO_DATA["feed"]["fresh_feed"]
    concentrate = scenario["initial_streams"]["cleaners34_concentrate"]
    assert 100 * stream_total(concentrate) / stream_total(feed) == pytest.approx(
        9.6, abs=0.2
    )
