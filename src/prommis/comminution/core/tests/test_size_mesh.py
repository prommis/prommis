#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the size-mesh core helpers."""

import math

import pytest

from prommis.comminution.core.size_mesh import (
    characteristic_sizes,
    validate_mesh,
)

# Sizes are in mm; the helpers require consistent units but do not attach units.
HANDCALC_EDGES = [1.0, 2.0, 4.0, 8.0, 16.0]


@pytest.mark.unit
def test_validate_mesh_rejections():
    assert validate_mesh(HANDCALC_EDGES) == tuple(HANDCALC_EDGES)
    assert validate_mesh([1.0, 2.0]) == (1.0, 2.0)
    # A zero lower edge requires 0 < bottom_size < the next edge.
    assert validate_mesh([0.0, 1.3, 2.6], bottom_size=0.1) == (0.0, 1.3, 2.6)

    with pytest.raises(ValueError, match="at least 2 edges"):
        validate_mesh([1.0])
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_mesh([1.0, 2.0, 2.0, 4.0])
    with pytest.raises(ValueError, match="strictly increasing"):
        validate_mesh([1.0, 4.0, 2.0])
    with pytest.raises(ValueError, match="must all be finite"):
        validate_mesh([1.0, float("inf")])
    # Unconvertible inputs must raise ValueError, not TypeError.
    with pytest.raises(ValueError, match="convertible to floats"):
        validate_mesh([1.0, None])
    with pytest.raises(ValueError, match="must be a sequence of numbers"):
        validate_mesh(1.0)
    with pytest.raises(ValueError, match="finest edge x_0 must be > 0"):
        validate_mesh([0.0, 1.3, 2.6])

    with pytest.raises(
        ValueError, match="^bottom_size must be a number, not a boolean$"
    ):
        validate_mesh([0.0, 1.3], True)
    with pytest.raises(
        ValueError, match="^bottom_size must be convertible to a float$"
    ):
        validate_mesh([0.0, 1.3], "x")
    with pytest.raises(ValueError, match="^bottom_size must be finite$"):
        validate_mesh([0.0, 1.3], float("nan"))
    with pytest.raises(
        ValueError,
        match="^bottom_size may only be supplied when the finest edge x_0 == 0$",
    ):
        validate_mesh([1.0, 2.0], 0.5)
    with pytest.raises(
        ValueError, match="^bottom_size must be a finite positive length$"
    ):
        validate_mesh([0.0, 1.3], 0.0)
    with pytest.raises(
        ValueError, match="^bottom_size must be strictly less than x_1$"
    ):
        validate_mesh([0.0, 1.3], 1.3)
    with pytest.raises(
        ValueError, match="^bottom_size must be strictly less than x_1$"
    ):
        validate_mesh([0.0, 1.3, 2.6], 2.0)


@pytest.mark.unit
def test_characteristic_sizes_handcalc():
    d_char = characteristic_sizes(HANDCALC_EDGES)
    expected = [math.sqrt(2.0), math.sqrt(8.0), math.sqrt(32.0), math.sqrt(128.0)]
    assert len(d_char) == 4
    for got, exp in zip(d_char, expected):
        assert got == pytest.approx(exp, rel=1e-12)
    # bottom_size replaces zero only in the first bin's geometric mean.
    d_char = characteristic_sizes([0.0, 1.3, 2.6], bottom_size=0.1)
    assert d_char[0] == pytest.approx(math.sqrt(0.1 * 1.3), rel=1e-12)
    assert d_char[1] == pytest.approx(math.sqrt(1.3 * 2.6), rel=1e-12)
    with pytest.raises(ValueError, match="finest edge x_0 must be > 0"):
        characteristic_sizes([0.0, 1.3, 2.6])
