#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests the two-salt diafiltration flowsheet
"""

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.unit_models import Feed, Product

import matplotlib.pyplot as plt

import pytest

from prommis.nanofiltration.diafiltration_flowsheet_two_salt import main
from prommis.util import assert_solution_equivalent


@pytest.mark.component
def test_main():
    """
    Tests the execution of the main function in diafiltration.py
    """
    (
        m,
        overall_results_plot,
        boundary_layer_results_plot,
        membrane_results_plot,
        rejection_plot,
    ) = main()
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
    assert isinstance(boundary_layer_results_plot, plt.Figure)
    assert isinstance(membrane_results_plot, plt.Figure)
    assert isinstance(rejection_plot, plt.Figure)

    test_dict = {
        "retentate_flow_volume": {(0, 1): (6.0854, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (190.89, 1e-4, None),
            (0, 1, "Co"): (239.83, 1e-4, None),
            (0, 1, "Cl"): (670.55, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (10.035, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Li"): (191.70, 1e-4, None),
            (0, 1, "Co"): (222.48, 1e-4, None),
            (0, 1, "Cl"): (636.67, 1e-4, None),
        },
    }

    assert_solution_equivalent(m.fs.membrane, test_dict)
