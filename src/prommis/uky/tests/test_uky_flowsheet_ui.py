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
    flowsheet = ui.build_flowsheet()
    ui.solve_flowsheet(flowsheet)
