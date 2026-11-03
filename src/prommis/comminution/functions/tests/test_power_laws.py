#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the comminution power laws."""

import math

import pytest

from pyomo.environ import units, value

from prommis.comminution.functions.power_laws import (
    invert_specific_energy,
    numeric_specific_energy,
    specific_energy,
)


def _kwh_per_tonne(expr):
    # W*h/kg equivalent to kWh/tonne
    return value(units.convert(expr, units.W * units.hour / units.kg))


@pytest.mark.unit
def test_specific_energy_values_and_coefficients():
    """Check Bond, Rittinger, and Kick energy values and coefficients.

    The Rittinger and Kick coefficients follow Zogg (1993, sec. 2.1.2.2).
    Expected values are calculated independently of implementation constants.
    """
    # Bond at clean round inputs:
    # 10*15*(1/sqrt(100) - 1/sqrt(10000)) = 150*(0.1 - 0.01) = 13.5 kWh/t.
    wi = 15 * units.W * units.hour / units.kg
    f80 = 10000e-6 * units.m  # 10000 um
    p80 = 100e-6 * units.m  # 100 um
    assert _kwh_per_tonne(specific_energy("bond", wi, p80, f80)) == pytest.approx(
        13.5, rel=1e-12
    )

    wi = 12 * units.W * units.hour / units.kg
    # Rittinger: C_R = 0.5*(10*E_i)*sqrt(d_Bl) ~ 424.26407 (Zogg 2.15).
    # W = C_R*(1/20 - 1/40) ~ 10.60660 kWh/t.
    got = _kwh_per_tonne(
        specific_energy("rittinger", wi, 20e-6 * units.m, 40e-6 * units.m)
    )
    assert got == pytest.approx(10.60660, abs=1e-5)
    # Kick: C_K = 1.151*(10*E_i)/sqrt(d_Bu) ~ 0.61769 (Zogg 2.14).
    # W = C_K*log10(200/100) ~ 0.18594 kWh/t.
    got = _kwh_per_tonne(
        specific_energy("kick", wi, 100e-3 * units.m, 200e-3 * units.m)
    )
    assert got == pytest.approx(0.18594, abs=1e-5)

    # Zogg's boundary relations, evaluated independently
    e_i = 12.0  # kWh/tonne == W*h/kg
    d_lower_um, d_upper_um = 50.0, 50000.0
    c_bond = 10.0 * e_i
    c_rittinger = 0.5 * c_bond * math.sqrt(d_lower_um)
    c_kick = 1.151 * c_bond / math.sqrt(d_upper_um)
    ritt = _kwh_per_tonne(
        specific_energy("rittinger", wi, 10e-6 * units.m, 40e-6 * units.m)
    )
    kick = _kwh_per_tonne(
        specific_energy("kick", wi, 60000e-6 * units.m, 200000e-6 * units.m)
    )
    assert ritt == pytest.approx(c_rittinger * (1 / 10.0 - 1 / 40.0), rel=1e-12)
    assert kick == pytest.approx(c_kick * math.log10(200000.0 / 60000.0), rel=1e-12)

    with pytest.raises(NotImplementedError):
        specific_energy("nonsense", wi, 80e-6 * units.m, 100e-6 * units.m)


@pytest.mark.unit
@pytest.mark.parametrize("law", ["bond", "rittinger", "kick"])
def test_invert_specific_energy_round_trip(law):
    # the inverse followed by the forward law returns the input energy, across
    # the crusher's work-index bounds and a nominal value
    f80 = 0.01
    for wi in (1.0, 14.0, 150.0):
        for e_spec in (0.0, 0.05, 12.0):
            p80 = invert_specific_energy(law, e_spec, wi, f80)
            assert 0.0 < p80 <= f80
            if e_spec > 0.0:
                assert p80 < f80
            else:
                assert abs(p80 - f80) <= 2 * math.ulp(f80)
            assert numeric_specific_energy(law, wi, p80, f80) == pytest.approx(
                e_spec, rel=1e-12, abs=1e-12
            )
    for e_spec in (-1e-12, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="e_spec must be finite and >= 0"):
            invert_specific_energy(law, e_spec, 14.0, f80)
