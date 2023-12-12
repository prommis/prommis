import pytest
from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from pyomo.environ import ConcreteModel, check_optimal_termination, value

from prommis.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from prommis.Solvent_Extraction.REEOgdistribution import REESolExOgParameters
from prommis.Solvent_Extraction.SolventExtraction import SolventExtraction

solver = get_solver()


class TestSXmodel:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.prop_a = REESolExAqParameters()
        m.fs.prop_o = REESolExOgParameters()

        m.fs.solex = SolventExtraction(
            number_of_finite_elements=3,
            dynamic=False,
            aqueous_stream={
                "property_package": m.fs.prop_a,
                "flow_direction": FlowDirection.forward,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            organic_stream={
                "property_package": m.fs.prop_o,
                "flow_direction": FlowDirection.backward,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        )

        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(820)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(270)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)

        m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol.fix(4.4)

        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(7.54e-10)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(4.955e-9)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(1.491e-7)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(321.34)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(5.67e-6)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(1.78e-05)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(4.019e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(6.73e-6)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(1.82e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(3.285e-6)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(1.55e-6)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(9e-7)

        m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

        return m

    @pytest.mark.known_issue(6)
    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.report_structural_issues()
        dt.display_underconstrained_set()
        dt.assert_no_structural_warnings()

    @pytest.mark.known_issue(6)
    @pytest.mark.component
    def test_block_triangularization(self, SolEx_frame):
        model = SolEx_frame
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.solex)

        assert initializer.summary[model.fs.solex]["status"] == InitializationStatus.Ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, SolEx_frame):
        m = SolEx_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.known_issue(6)
    @pytest.mark.component
    def test_solution(self, SolEx_frame):
        m = SolEx_frame
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Al"]
        ) == pytest.approx(730, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Ca"]
        ) == pytest.approx(4680, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Fe"]
        ) == pytest.approx(250, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Sc"]
        ) == pytest.approx(9.99e-03, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Y"]
        ) == pytest.approx(9.99e-03, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["La"]
        ) == pytest.approx(30.84, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Ce"]
        ) == pytest.approx(0.36, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Pr"]
        ) == pytest.approx(0.0312, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Nd"]
        ) == pytest.approx(9.99e-03, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Sm"]
        ) == pytest.approx(9.99e-03, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Gd"]
        ) == pytest.approx(9.99e-03, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Dy"]
        ) == pytest.approx(9.99e-03, rel=1e-2)

        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Al"]
        ) == pytest.approx(6.03, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Ca"]
        ) == pytest.approx(39.64, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Fe"]
        ) == pytest.approx(1.19, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Sc"]
        ) == pytest.approx(336.25, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Y"]
        ) == pytest.approx(45.41, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["La"]
        ) == pytest.approx(142.56, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Ce"]
        ) == pytest.approx(321.59, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Pr"]
        ) == pytest.approx(53.88, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Nd"]
        ) == pytest.approx(145.83, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Sm"]
        ) == pytest.approx(26.28, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Gd"]
        ) == pytest.approx(12.42, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Dy"]
        ) == pytest.approx(7.20, rel=1e-2)
