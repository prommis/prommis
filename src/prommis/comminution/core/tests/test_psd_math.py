#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the PSD math core (cumulative passing, percentiles, discretization)."""

import math

import pytest

from prommis.comminution.core.psd_math import (
    EPS_SMOOTH,
    cumulative_passing,
    finest_attainable_size,
    fractions_from_cdf,
    size_at_passing,
    size_at_passing_smooth,
)

HANDCALC_EDGES = [1.0, 2.0, 4.0, 8.0, 16.0]
HANDCALC_FLOWS = [0.1, 0.2, 0.3, 0.4]
HANDCALC_CUM = [0.1, 0.3, 0.6, 1.0]


@pytest.mark.unit
def test_size_at_passing_hand_values():
    cum = cumulative_passing(HANDCALC_FLOWS)
    assert cum == pytest.approx(HANDCALC_CUM, rel=1e-12)
    assert cum[-1] == 1.0
    # P80 = 8 + (0.8-0.6)/(1.0-0.6)*(16-8) = 12.0 ; bracket [8, 16]
    assert size_at_passing(HANDCALC_EDGES, HANDCALC_CUM, 0.8) == pytest.approx(
        12.0, rel=1e-9
    )
    # P50 = 4 + (0.5-0.3)/(0.6-0.3)*(8-4) = 20/3 ; bracket [4, 8]
    assert size_at_passing(HANDCALC_EDGES, HANDCALC_CUM, 0.5) == pytest.approx(
        20.0 / 3.0, rel=1e-9
    )
    assert size_at_passing(HANDCALC_EDGES, HANDCALC_CUM, 1.0) == 16.0
    # Treat the lower edge as 0% passing; P5 lies in the first bin.
    assert size_at_passing(HANDCALC_EDGES, HANDCALC_CUM, 0.05) == pytest.approx(
        1.5, rel=1e-12
    )
    # At an exact plateau hit, use the first edge with that passing fraction.
    cum = cumulative_passing([1.0, 2.0, 0.0, 1.0])
    assert size_at_passing(HANDCALC_EDGES, cum, 0.75) == 4.0
    edges = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    cum = cumulative_passing([1.0, 0.0, 0.0, 1.0, 0.0])
    # After the flat section, P80 lies between 8 and 16: 8 + 0.6 * 8 = 12.8.
    assert size_at_passing(edges, cum, 0.8) == pytest.approx(12.8, rel=1e-12)
    assert size_at_passing([2.0, 10.0], [1.0], 0.5) == pytest.approx(6.0, rel=1e-12)
    # All mass is in the finest bin; minimum P80 = 1 + 0.8 * (2 - 1) = 1.8.
    cum = cumulative_passing([10.0, 0.0, 0.0, 0.0])
    assert size_at_passing(HANDCALC_EDGES, cum, 0.8) == pytest.approx(1.8, rel=1e-12)
    assert finest_attainable_size(HANDCALC_EDGES) == pytest.approx(1.8, rel=1e-12)


@pytest.mark.unit
def test_size_at_passing_smooth_is_finite_on_step_curve():
    # Empty bins create flat sections that can cause division by zero.
    cum = cumulative_passing([0.0, 10.0, 0.0, 0.0])
    val = size_at_passing_smooth(HANDCALC_EDGES, cum, 0.8, eps=1e-5)
    assert math.isfinite(val)
    # P80 lies between 2 and 4: 2 + 0.8 * (4 - 2) = 3.6.
    assert val == pytest.approx(3.6, rel=2e-3)
    # Away from bin boundaries, smoothing should closely match exact percentiles.
    assert size_at_passing_smooth(HANDCALC_EDGES, HANDCALC_CUM, 0.8) == pytest.approx(
        12.0, rel=1e-4
    )
    assert size_at_passing_smooth(HANDCALC_EDGES, HANDCALC_CUM, 0.5) == pytest.approx(
        20.0 / 3.0, rel=1e-4
    )
    # A plateau at 0.8 should return the first matching edge (4.0), within
    # 0.001 for smoothing; the later matching edge is 8.0.
    assert size_at_passing_smooth(
        HANDCALC_EDGES, [0.2, 0.8, 0.8, 1.0], 0.8, EPS_SMOOTH
    ) == pytest.approx(4.0, abs=1e-3)


@pytest.mark.unit
def test_fractions_from_cdf_sums_to_one_and_bottom_closure():
    # Use a Rosin-Rammler distribution.
    d63, n = 4.0, 1.5
    cdf = lambda d: 1.0 - math.exp(-((d / d63) ** n))
    edges = [0.5, 1.0, 2.0, 4.0, 8.0, 40.0]
    fracs = fractions_from_cdf(cdf, edges)
    assert sum(fracs) == pytest.approx(1.0, rel=1e-12)
    # bottom closure folds the sub-grid tail into bin 0: frac[0] = Q(x_1)/Q(x_N)
    q_top = cdf(edges[-1])
    assert fracs[0] == pytest.approx(cdf(edges[1]) / q_top, rel=1e-12)
    assert fracs[2] == pytest.approx((cdf(edges[3]) - cdf(edges[2])) / q_top, rel=1e-12)
    # a zero or NaN top CDF cannot normalize the bins
    with pytest.raises(ValueError):
        fractions_from_cdf(lambda d: 0.0, [1.0, 2.0, 4.0])
    with pytest.raises(ValueError):
        fractions_from_cdf(lambda d: float("nan"), [1.0, 2.0, 4.0])
