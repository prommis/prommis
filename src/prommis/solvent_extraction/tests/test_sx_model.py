from pyomo.environ import ConcreteModel, check_optimal_termination, value

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox

import pytest

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.ree_aq_distribution import REESolExAqParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

solver = get_solver()


class TestSXmodel:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.leach_soln = LeachSolutionParameters()
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
            aqueous_to_organic=True,
        )

        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sc"] = 99.9 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "La"] = 75.2 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ce"] = 95.7 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Pr"] = 96.5 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Nd"] = 99.2 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sm"] = 99.9 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Gd"] = 98.6 / 100
        m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Dy"] = 99.9 / 100

        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e-9)
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

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.report_structural_issues()
        dt.display_underconstrained_set()
        dt.assert_no_structural_warnings()

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

    @pytest.mark.component
    def test_solution(self, SolEx_frame):
        m = SolEx_frame
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["H2O"]
        ) == pytest.approx(1e-9, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["H"]
        ) == pytest.approx(1e-9, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["SO4"]
        ) == pytest.approx(1e-9, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["HSO4"]
        ) == pytest.approx(1e-9, rel=1e-2)
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
        ) == pytest.approx(2.093e-7, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Y"]
        ) == pytest.approx(6.377e-07, rel=1e-2)
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
        ) == pytest.approx(1.084e-3, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Sm"]
        ) == pytest.approx(3.69e-07, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Gd"]
        ) == pytest.approx(4.784e-04, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Dy"]
        ) == pytest.approx(1.011e-07, rel=1e-2)

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
        assert value(m.fs.solex.clean_sx_pe_tank_cap) == pytest.approx(1405, rel=1e-2)
        assert value(m.fs.solex.num_sx_pe_tank) == pytest.approx(5, rel=1e-2)
        assert value(m.fs.solex.clean_sx_tank_mixer_pow) == pytest.approx(1, rel=1e-2)
        assert value(m.fs.solex.num_sx_tank_mixer) == pytest.approx(2, rel=1e-2)
        assert value(m.fs.solex.clean_sx_process_pump_feed) == pytest.approx(
            281, rel=1e-2
        )
        assert value(m.fs.solex.num_sx_process_pump) == pytest.approx(3, rel=1e-2)
        assert value(m.fs.solex.clean_sx_mix_set_cap) == pytest.approx(2444, rel=1e-2)
        assert value(m.fs.solex.num_sx_mix_set) == pytest.approx(6, rel=1e-2)
