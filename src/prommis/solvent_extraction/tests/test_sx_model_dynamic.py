from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
    TransformationFactory,
    units,
)

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

        m = ConcreteModel()

        time_duration = 60

        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

        m.fs.prop_o = REESolExOgParameters()
        m.fs.prop_a = REESolExAqParameters()
        m.fs.leach_soln = LeachSolutionParameters()

        number_of_stages = 3

        m.fs.solex = SolventExtraction(
            number_of_finite_elements=number_of_stages,
            aqueous_stream={
                "property_package": m.fs.leach_soln,
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

        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=5, wrt=m.fs.time, scheme="BACKWARD")

        m.fs.solex.mscontactor.volume[:].fix(0.4)

        m.fs.solex.mscontactor.volume_frac_stream[:, :, "aqueous"].fix(0.5)

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

        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(820)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(5230)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(270)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(209.31)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(637.74)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(2032.77)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(4516.13)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(756.64)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(2047.85)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(369.1)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(174.38)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(101.12)

        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(1e-4)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(1e-4)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(1e-9)
        m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(1e-9)

        m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0,:,'Ka2'].fix(0)

        m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(4.4)

        m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(4.4)

        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(7.54e-10)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(4.955e-9)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(1.491e-7)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(321.34)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(5.67e-6)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(1.78e-05)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(4.019e-5)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(6.73e-6)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(1.82e-5)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(3.285e-6)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(1.55e-6)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(9e-7)

        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Y"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ce"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Nd"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sm"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Gd"].fix(1e-9)
        m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Dy"].fix(1e-9)

        m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

        m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)

        return m

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
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

    @pytest.mark.integration
    def test_solution(self, SolEx_frame):
        m = SolEx_frame
        time_duration = 60
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["HSO4"]
        ) == pytest.approx(7.651e-5, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["H"]
        ) == pytest.approx(8.4708, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["SO4"]
        ) == pytest.approx(9.0714e-5, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Al"]
        ) == pytest.approx(120.183, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Ca"]
        ) == pytest.approx(764.155, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Fe"]
        ) == pytest.approx(41.448, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Sc"]
        ) == pytest.approx(3.455e-8, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Y"]
        ) == pytest.approx(1.046e-07, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["La"]
        ) == pytest.approx(5.072, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Ce"]
        ) == pytest.approx(0.0587, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Pr"]
        ) == pytest.approx(0.0053, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Nd"]
        ) == pytest.approx(0.000171, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Sm"]
        ) == pytest.approx(7.0718e-08, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Gd"]
        ) == pytest.approx(7.828e-05, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[time_duration, 3].conc_mass_comp["Dy"]
        ) == pytest.approx(2.7748e-08, rel=1e-1)

        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Al"]
        ) == pytest.approx(3.969, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Ca"]
        ) == pytest.approx(26.001, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Fe"]
        ) == pytest.approx(0.769, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Sc"]
        ) == pytest.approx(334.225, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Y"]
        ) == pytest.approx(45.216, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["La"]
        ) == pytest.approx(127.656, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Ce"]
        ) == pytest.approx(315.161, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Pr"]
        ) == pytest.approx(52.968, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Nd"]
        ) == pytest.approx(144.825, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Sm"]
        ) == pytest.approx(26.169, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Gd"]
        ) == pytest.approx(12.305, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp["Dy"]
        ) == pytest.approx(7.169, rel=1e-1)
