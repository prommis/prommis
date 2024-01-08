import pytest
from prommis.uky import uky_flowsheet_ui as ui


@pytest.mark.unit
def test_export():
    interface = ui.export_to_ui()


@pytest.mark.component
def test_build():
    ui.build_flowsheet()


@pytest.mark.component
def test_solve():
    flowsheet = ui.build_flowsheet()
    ui.solve_flowsheet(flowsheet)

