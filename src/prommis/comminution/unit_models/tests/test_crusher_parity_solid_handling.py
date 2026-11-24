#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Compare the PSD crusher with the solid-handling crusher reference case.

At 2000 kg/h, the reference feed P80 is 114313.01 um, product P80 is
82876.93 um, and power is 123.83 W (``TestSolidHandling.test_solution``).
The solid-handling ``particle_size_median`` is the Rosin-Rammler
characteristic size ``d_63``. Its width of 1.5 maps to ``n = 4/3``.
At fixed P80, Bond power does not depend on ``n``, so the P50/P80 check
tests product shape separately.
"""

import math

import pytest

from pyomo.environ import value

from prommis.comminution.core.size_mesh import geometric_series
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    crusher,
    fix_feed,
    one_mineral_flowsheet,
)

SOLID_HANDLING_F80_UM = 114313.01
SOLID_HANDLING_P80_UM = 82876.93
SOLID_HANDLING_WORK_W = 123.83
SOLID_HANDLING_WI = 12.0
MDOT_KGPS = 2000.0 / 3600.0
SOLID_HANDLING_WIDTH = 1.5
SOLID_HANDLING_N = 2.0 / SOLID_HANDLING_WIDTH  # = 4/3


def _bond_work_w(wi_whkg, mdot_kgps, p80_m, f80_m):
    e_spec = 10 * wi_whkg * (1 / math.sqrt(p80_m * 1e6) - 1 / math.sqrt(f80_m * 1e6))
    return mdot_kgps * e_spec * 3600.0


def _two_bin_feed(edges_m, target_f80_m, total_kgps):
    """Set adjacent feed bins to give the reference F80 by linear interpolation."""
    ki = next(
        k
        for k in range(len(edges_m) - 1)
        if edges_m[k] <= target_f80_m < edges_m[k + 1]
    )
    xlo, xhi = edges_m[ki], edges_m[ki + 1]
    g = (target_f80_m - xlo) / (xhi - xlo)
    f = (0.8 - g) / (1.0 - g)
    size = [0.0] * (len(edges_m) - 1)
    size[ki - 1] = f * total_kgps
    size[ki] = (1.0 - f) * total_kgps
    return size


def _parity_unit(target_product_p80=None, power_kw=None):
    edges = list(geometric_series(1.0, 1e-3, 1000 ** (1.0 / 30)))  # 30 intervals
    m = one_mineral_flowsheet(edges=edges)
    cfg = dict(
        psd_method="distribution_function",
        distribution_shape="rosin_rammler",
        shape_exponent=SOLID_HANDLING_N,
        bond_work_index=SOLID_HANDLING_WI,
    )
    if target_product_p80 is not None:
        cfg["target_product_p80"] = target_product_p80
    cr = crusher(m, **cfg)
    rows = _two_bin_feed(edges, SOLID_HANDLING_F80_UM * 1e-6, MDOT_KGPS)
    fix_feed(cr.properties_in[0], {"OreA": rows})
    if power_kw is not None:
        cr.power[0].fix(power_kw)
    return m


@pytest.mark.component
@pytest.mark.solver
def test_value_parity():
    m = _parity_unit(target_product_p80=SOLID_HANDLING_P80_UM * 1e-6)
    cr = m.fs.cr
    # Initialization already solves the model.
    native.CrusherSolidPSDInitializer().initialize(cr)
    f80_m = value(cr.properties_in[0].percentile_size[0.8])
    p80_m = value(cr.product_p80[0])
    assert f80_m == pytest.approx(SOLID_HANDLING_F80_UM * 1e-6, rel=1e-4)
    assert p80_m == pytest.approx(SOLID_HANDLING_P80_UM * 1e-6, rel=1e-4)
    # Convert crusher power from kW to W for the reference comparison.
    power_w = value(cr.power[0]) * 1000.0
    recompute = _bond_work_w(SOLID_HANDLING_WI, MDOT_KGPS, p80_m, f80_m)
    assert power_w == pytest.approx(recompute, rel=1e-8)
    assert power_w == pytest.approx(SOLID_HANDLING_WORK_W, rel=1e-3)
    total_out = sum(
        value(cr.properties_out[0].flow_mass_sized_comp_size["OreA", k])
        for k in m.fs.pp.size_interval_set
    )
    assert total_out == pytest.approx(MDOT_KGPS, rel=1e-8)
    # The P50/P80 ratio checks the exponent n = 2/width.
    out = cr.properties_out[0]
    assert value(out.size_at_passing(0.5)) / value(
        out.percentile_size[0.8]
    ) == pytest.approx((math.log(2.0) / math.log(5.0)) ** 0.75, rel=1e-2)


@pytest.mark.component
@pytest.mark.solver
def test_power_spec_round_trip():
    m = _parity_unit(power_kw=SOLID_HANDLING_WORK_W / 1000.0)
    native.CrusherSolidPSDInitializer().initialize(m.fs.cr)
    assert value(m.fs.cr.product_p80[0]) == pytest.approx(
        SOLID_HANDLING_P80_UM * 1e-6, rel=1e-4
    )
