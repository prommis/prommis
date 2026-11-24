#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for crusher distribution sizes, power, and invalid specifications."""

import math

import pytest

from pyomo.environ import value


from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    EDGES_REF,
    RR,
    crusher,
    fix_feed,
    one_mineral_flowsheet,
)


@pytest.mark.component
@pytest.mark.solver
def test_distribution_function_zero_power_preserves_p80_not_full_psd():
    m = one_mineral_flowsheet()
    cr = crusher(
        m,
        psd_method="distribution_function",
        distribution_shape=RR,
        shape_exponent=1.5,
        bond_work_index=12.0,
    )
    feed = [1.0, 1.0, 7.0, 1.0]  # F80 about 7.4 mm; the RR shape is mesh-valid
    fix_feed(cr.properties_in[0], {"OreA": feed})
    cr.power[0].fix(0.0)
    native.CrusherSolidPSDInitializer().initialize(cr)
    assert value(cr.product_p80[0]) == pytest.approx(
        value(cr.properties_in[0].percentile_size[0.8]), rel=1e-8
    )
    total = value(cr.properties_out[0].flow_mass_sized)
    prod = [
        value(cr.properties_out[0].flow_mass_sized_comp_size["OreA", k]) / total
        for k in range(4)
    ]
    feed_frac = [v / sum(feed) for v in feed]
    assert any(abs(p - f) > 0.02 for p, f in zip(prod, feed_frac))


@pytest.mark.component
@pytest.mark.solver
def test_distribution_function_rr_size_spec_bins_match_analytic():
    n, d80 = 1.5, 5e-3
    m = one_mineral_flowsheet()
    cr = crusher(
        m,
        psd_method="distribution_function",
        distribution_shape=RR,
        shape_exponent=n,
        target_product_p80=d80,
        bond_work_index=12.0,
    )
    fix_feed(cr.properties_in[0], {"OreA": [0.0, 0.0, 0.0, 10.0]})
    native.CrusherSolidPSDInitializer().initialize(cr)
    d63 = d80 / ((-math.log(0.2)) ** (1.0 / n))
    edges = list(EDGES_REF)
    cdf = [1.0 - math.exp(-((e / d63) ** n)) for e in edges]
    q_top = cdf[-1]
    expected = [cdf[1] / q_top] + [(cdf[k + 1] - cdf[k]) / q_top for k in range(1, 4)]
    total = value(cr.properties_out[0].flow_mass_sized)
    got = [
        value(cr.properties_out[0].flow_mass_sized_comp_size["OreA", k]) / total
        for k in range(4)
    ]
    for g, exp in zip(got, expected):
        assert g == pytest.approx(exp, rel=1e-8)
