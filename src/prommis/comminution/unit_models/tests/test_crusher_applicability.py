#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for crusher equipment presets, stage selection, and diagnostics."""

import pytest


from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    crusher,
    fix_feed,
    one_mineral_flowsheet,
)


def _tabular_with_stage(q, **stage):
    edges = [0.005, 0.02, 0.05, 0.1, 0.15, 0.3]  # m, 0.5 to 30 cm
    m = one_mineral_flowsheet(edges=edges)
    cr = crusher(
        m,
        psd_method="tabular",
        tabular_psd=q,
        bond_work_index=12.0,
        **stage,
    )
    fix_feed(cr.properties_in[0], {"OreA": [0.0, 0.0, 0.0, 0.0, 10.0]})
    native.CrusherSolidPSDInitializer().initialize(cr)
    # Return the parent model alongside the unit to keep it alive.
    return m, cr


@pytest.mark.component
@pytest.mark.solver
def test_check_applicability_in_range_silent_and_out_warns(caplog):
    m_in, in_range = _tabular_with_stage(
        [0.0, 1.0, 0.0, 0.0, 0.0], crusher_stage="secondary"
    )
    assert in_range.check_applicability(emit_warning=False) == []
    m_out, out = _tabular_with_stage(
        [0.0, 0.0, 0.0, 0.0, 1.0], crusher_equipment="cone"
    )
    assert out.config.crusher_stage == "secondary"
    caplog.set_level("WARNING", logger="idaes.prommis.comminution.unit_models.crusher")
    caplog.clear()
    messages = out.check_applicability()
    assert len(messages) == 1
    assert [r.getMessage() for r in caplog.records] == [
        f"CrusherSolidPSD applicability: {messages[0]}"
    ]
