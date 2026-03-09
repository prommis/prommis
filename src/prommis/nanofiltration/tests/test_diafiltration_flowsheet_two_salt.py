#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests the two-salt diafiltration flowsheet
"""

from pyomo.environ import value

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.unit_models import Feed, Product

import matplotlib.pyplot as plt

import pytest

from prommis.nanofiltration.diafiltration_flowsheet_two_salt import main


@pytest.mark.component
def test_main():
    """
    Tests the execution of the main function in diafiltration.py
    """
    m, overall_results_plot, membrane_results_plot = main()
    dt = DiagnosticsToolbox(m)
    dt.assert_no_numerical_warnings()

    # verify flowsheet structure
    assert isinstance(m.fs.feed_block, Feed)
    assert isinstance(m.fs.diafiltrate_block, Feed)
    assert isinstance(m.fs.retentate_block, Product)
    assert isinstance(m.fs.permeate_block, Product)
    assert hasattr(m.fs, "feed_stream")
    assert hasattr(m.fs, "diafiltrate_stream")
    assert hasattr(m.fs, "retentate_stream")
    assert hasattr(m.fs, "permeate_stream")

    # verify necessary variables are fixed
    assert m.fs.membrane.total_module_length.fixed
    assert m.fs.membrane.total_membrane_length.fixed
    for t in m.fs.membrane.time:
        assert m.fs.membrane.applied_pressure[t].fixed
        assert m.fs.membrane.feed_flow_volume[t].fixed
        assert m.fs.membrane.diafiltrate_flow_volume[t].fixed
        for j in m.fs.membrane.solutes:
            assert m.fs.membrane.feed_conc_mol_comp[t, j].fixed
            assert m.fs.membrane.diafiltrate_conc_mol_comp[t, j].fixed

    # verify plots exist
    assert isinstance(overall_results_plot, plt.Figure)
    assert isinstance(membrane_results_plot, plt.Figure)

    test_dict = {
        "retentate_final": [
            value(m.fs.membrane.retentate_flow_volume[0, 1]),
            6.0465,
        ],
        "Li_retentate_final": [
            value(m.fs.membrane.retentate_conc_mol_comp[0, 1, "Li"]),
            188.88,
        ],
        "Co_retentate_final": [
            value(m.fs.membrane.retentate_conc_mol_comp[0, 1, "Co"]),
            246.67,
        ],
        "Cl_retentate_final": [
            value(m.fs.membrane.retentate_conc_mol_comp[0, 1, "Cl"]),
            682.22,
        ],
        "permeate_final": [
            value(m.fs.membrane.permeate_flow_volume[0, 1]),
            10.033,
        ],
        "Li_permeate_final": [
            value(m.fs.membrane.permeate_conc_mol_comp[0, 1, "Li"]),
            191.55,
        ],
        "Co_permeate_final": [
            value(m.fs.membrane.permeate_conc_mol_comp[0, 1, "Co"]),
            222.76,
        ],
        "Cl_permeate_final": [
            value(m.fs.membrane.permeate_conc_mol_comp[0, 1, "Cl"]),
            637.06,
        ],
    }

    for model_result, test_val in test_dict.values():
        assert pytest.approx(test_val, rel=1e-4) == value(model_result)
