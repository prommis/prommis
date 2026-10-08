#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the config-coercion helpers ``_cfg_float`` and ``_cfg_quantity``."""

import math
from decimal import Decimal

import numpy
import pytest

from pyomo.environ import units

from idaes.core.util.exceptions import ConfigurationError

from prommis.comminution.core.config_utils import _cfg_float, _cfg_quantity


@pytest.mark.unit
def test_cfg_float_validation():
    def check(val):
        return _cfg_float(val, "k", lo=0.5, hi=2.0)

    for val in (True, numpy.bool_(True)):
        with pytest.raises(ConfigurationError, match="not a boolean"):
            check(val)
    with pytest.raises(ConfigurationError, match="not a string"):
        check("1")
    for bad in (Decimal("1"), 1 + 0j, 0.4 * units.dimensionless):
        with pytest.raises(ConfigurationError, match="must be a real int or float"):
            check(bad)
    with pytest.raises(ConfigurationError, match="must be convertible to a float"):
        check(10**400)
    with pytest.raises(ConfigurationError, match="must be finite"):
        check(float("nan"))
    with pytest.raises(ConfigurationError, match="must be finite"):
        check(float("inf"))
    with pytest.raises(ConfigurationError, match=r"must be >= 0\.5"):
        check(math.nextafter(0.5, -math.inf))
    with pytest.raises(ConfigurationError, match=r"must be <= 2\.0"):
        check(math.nextafter(2.0, math.inf))

    assert check(0.5) == 0.5
    assert check(2.0) == 2.0
    for val in (1, numpy.float64(1), numpy.float32(1), numpy.int64(1)):
        out = check(val)
        assert type(out) is float and out == 1.0


@pytest.mark.unit
def test_quantity_same_unit_is_the_bare_number():
    assert _cfg_quantity(25.0 * units.mm, "k", units.mm) == 25.0


@pytest.mark.unit
def test_quantity_converts_across_units():
    assert _cfg_quantity(0.025 * units.m, "k", units.mm) == pytest.approx(25.0)
    assert _cfg_quantity(2.5 * units.cm, "k", units.mm) == pytest.approx(25.0)


@pytest.mark.unit
def test_quantity_converts_compound_units():
    rho = _cfg_quantity(
        5400.0 * units.kg / units.m**3, "k", units.metric_ton / units.m**3
    )
    assert rho == pytest.approx(5.4)
    wi = _cfg_quantity(
        14.0 * units.kWh / units.metric_ton, "k", units.kWh / units.metric_ton
    )
    assert wi == pytest.approx(14.0)


@pytest.mark.unit
@pytest.mark.parametrize("bad", [25, 25.0, 0, -3.5])
def test_bare_number_is_rejected(bad):
    with pytest.raises(ConfigurationError, match="must carry units"):
        _cfg_quantity(bad, "k", units.mm)


@pytest.mark.unit
def test_bool_is_rejected():
    with pytest.raises(ConfigurationError, match="must carry units"):
        _cfg_quantity(True, "k", units.mm)


@pytest.mark.unit
def test_string_is_rejected():
    with pytest.raises(ConfigurationError, match="not a string"):
        _cfg_quantity("25", "k", units.mm)


@pytest.mark.unit
@pytest.mark.parametrize(
    "bad",
    [
        25.0 * units.kg,
        25.0 * units.dimensionless,
        25.0 * units.s,
    ],
)
def test_wrong_dimension_is_rejected(bad):
    with pytest.raises(ConfigurationError, match="must be convertible to"):
        _cfg_quantity(bad, "k", units.mm)


@pytest.mark.unit
def test_domain_checks_apply_after_conversion():
    # Bounds use the target unit after conversion.
    with pytest.raises(ConfigurationError, match="must be > 0"):
        _cfg_quantity(-1.0 * units.mm, "k", units.mm, positive=True)
    with pytest.raises(ConfigurationError, match="must be >= 2"):
        _cfg_quantity(0.001 * units.m, "k", units.mm, lo=2.0)
    with pytest.raises(ConfigurationError, match="must be <= 5"):
        _cfg_quantity(0.01 * units.m, "k", units.mm, hi=5.0)
    assert _cfg_quantity(1.0 * units.m, "k", units.mm, lo=2.0) == pytest.approx(1000.0)


@pytest.mark.unit
def test_non_finite_quantity_is_rejected():
    with pytest.raises(ConfigurationError, match="must be finite"):
        _cfg_quantity(float("nan") * units.mm, "k", units.mm)
    with pytest.raises(ConfigurationError, match="must be finite"):
        _cfg_quantity(float("inf") * units.mm, "k", units.mm)


@pytest.mark.unit
def test_error_message_names_the_target_unit_and_the_key():
    with pytest.raises(ConfigurationError) as exc:
        _cfg_quantity(3.0, "geometry[diameter]", units.m)
    assert "geometry[diameter]" in str(exc.value)
    assert "m" in str(exc.value)
