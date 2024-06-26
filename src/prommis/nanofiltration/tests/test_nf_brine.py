#################################################################################
# adapted from test_nf.py
#################################################################################

"""
Tests the nanofiltration flowsheet
"""

from pyomo.environ import value

import pytest

pytest.importorskip("watertap", reason="WaterTAP dependency not available")
# pylint: disable-next=wrong-import-position
from prommis.nanofiltration.nf_brine import main


@pytest.mark.component
def test_main():
    """
    Tests the execution of the main function in nf_brine.py
    """
    m = main()
    test_dict = {
        "pressure": [m.fs.pump.outlet.pressure[0] * 1e-5, 11.362344380368109],
        "area": [m.fs.unit.area, 999.9990553265868],
        "nf_recovery": [
            m.fs.unit.recovery_vol_phase[0.0, "Liq"] * 100,
            94.9999989334488,
        ],
        "li_rejection": [
            m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Li_+"].value * 100,
            1.8379272407186817,
        ],
        "feed_ion_ratio": [
            (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
            / (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
            0.5077446970643547,
        ],
        "perm_ion_ratio": [
            (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
            / (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
            0.4972950260825848,
        ],
    }
    for model_result, testval in test_dict.values():
        assert pytest.approx(testval, rel=1e-5) == value(model_result)
