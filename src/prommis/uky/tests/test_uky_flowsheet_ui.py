#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import pytest

ui = pytest.importorskip(
    "prommis.uky.uky_flowsheet_ui", reason="Flowsheet UI components not available"
)


@pytest.mark.unit
def test_uky_ui_export():
    interface = ui.export_to_ui()


@pytest.mark.component
def test_uky_ui_build():
    ui.build_flowsheet()


@pytest.mark.component
def test_uky_ui_solve():
    interface = ui.export_to_ui()
    interface.build()
    interface.solve()
