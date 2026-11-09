#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the smooth clamp primitives."""

import pytest

from prommis.comminution.core.smooth_functions import smooth_clamp, smooth_max0


@pytest.mark.unit
def test_smooth_max0():
    assert smooth_max0(5.0, eps=1e-6) == pytest.approx(5.0, abs=1e-5)
    assert smooth_max0(-5.0, eps=1e-6) == pytest.approx(0.0, abs=1e-5)


@pytest.mark.unit
@pytest.mark.parametrize(
    "input_value,expected",
    [(-2.0, 1.0), (1.0, 1.0), (2.5, 2.5), (4.0, 4.0), (7.0, 4.0)],
)
def test_smooth_clamp(input_value, expected):
    assert smooth_clamp(input_value, 1.0, 4.0, eps=1e-6) == pytest.approx(
        expected, rel=0.0, abs=1e-5
    )
