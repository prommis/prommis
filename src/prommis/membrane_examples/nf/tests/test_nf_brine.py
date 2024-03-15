#################################################################################
# adapted from test_nf.py
#################################################################################

"""
Tests the nanofiltration flowsheet
"""

from pyomo.environ import value

import pytest

from prommis.membrane_examples.nf.nf_brine import main


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_main():
    """
    Tests the execution of the main function in nf_brine.py
    """
    m = main()
    test_dict = {
        "pressure": [m.fs.pump.outlet.pressure[0] * 1e-5, 2.0],
        "area": [m.fs.unit.area, 985.283490],
        "nf_recovery": [
            m.fs.unit.recovery_vol_phase[0.0, "Liq"] * 100,
            9.983259,
        ],
        "li_rejection": [
            m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Li_+"].value * 100,
            3.694042,
        ],
        "feed_ion_ratio": [
            (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
            / (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
            0.507745,
        ],
        "perm_ion_ratio": [
            (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
            / (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
            0.493056,
        ],
    }
    for model_result, testval in test_dict.values():
        assert pytest.approx(testval, rel=1e-3) == value(model_result)
