#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, check_optimal_termination, value
import numpy as np

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
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

solver = get_solver()


class TestSXmodel:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.leach_soln = LeachSolutionParameters()
        m.fs.prop_o = REESolExOgParameters()

        number_of_stages = 3

        m.fs.solex = SolventExtraction(
            number_of_finite_elements=number_of_stages,
            dynamic=False,
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

        stage_number = np.arange(1, number_of_stages + 1)

        for s in stage_number:
            if s == 1:
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = (
                    5.2 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
                    3 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
                    24.7 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
                    99.1 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = (
                    99.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = (
                    32.4 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
                    58.2 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
                    58.2 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
                    87.6 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
                    99.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
                    69.8 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
                    96.6 / 100
                )
            else:
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = (
                    4.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
                    12.3 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
                    6.4 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
                    16.7 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = (
                    99.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = (
                    23.2 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
                    24.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
                    15.1 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
                    99.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
                    99.9 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
                    7.6 / 100
                )
                m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
                    5 / 100
                )

        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e6)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1.755)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(
            3999.818
        )
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(
            693.459
        )
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(422.375)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(109.542)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(688.266)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(0.032)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(0.124)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(0.986)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(2.277)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(0.303)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(0.946)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(0.097)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(0.2584)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(0.047)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Cl"].fix(1e-8)

        m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol.fix(62.01)

        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(1.267e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(2.684e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(2.873e-6)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(1.734)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(2.179e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(0.000105)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(0.00031)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(3.711e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(0.000165)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(1.701e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(3.357e-5)
        m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(8.008e-6)

        m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

        return m

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_block_triangularization(self, SolEx_frame):
        model = SolEx_frame
        initializer = BlockTriangularizationInitializer(constraint_tolerance=1e-4)
        initializer.initialize(model.fs.solex)

        assert initializer.summary[model.fs.solex]["status"] == InitializationStatus.Ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, SolEx_frame):
        m = SolEx_frame
        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, SolEx_frame):
        m = SolEx_frame
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["H2O"]
        ) == pytest.approx(1e6, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["H"]
        ) == pytest.approx(1.755, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["SO4"]
        ) == pytest.approx(3999.885, abs=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["HSO4"]
        ) == pytest.approx(693.39, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Al"]
        ) == pytest.approx(362.132, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Ca"]
        ) == pytest.approx(81.724, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Fe"]
        ) == pytest.approx(454.049, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Sc"]
        ) == pytest.approx(0.000199, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Y"]
        ) == pytest.approx(1.239e-10, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["La"]
        ) == pytest.approx(0.393, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Ce"]
        ) == pytest.approx(0.536, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Pr"]
        ) == pytest.approx(0.091, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Nd"]
        ) == pytest.approx(1.17e-7, rel=1e-1)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Sm"]
        ) == pytest.approx(9.69e-11, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Gd"]
        ) == pytest.approx(0.066, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp["Dy"]
        ) == pytest.approx(0.00144, rel=1e-2)

        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Al"]
        ) == pytest.approx(60.242, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Ca"]
        ) == pytest.approx(27.817, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Fe"]
        ) == pytest.approx(234.216, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Sc"]
        ) == pytest.approx(1.765, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Y"]
        ) == pytest.approx(0.124, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["La"]
        ) == pytest.approx(0.592, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Ce"]
        ) == pytest.approx(1.7405, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Pr"]
        ) == pytest.approx(0.211, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Nd"]
        ) == pytest.approx(0.946, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Sm"]
        ) == pytest.approx(0.097, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Gd"]
        ) == pytest.approx(0.1918, rel=1e-2)
        assert value(
            m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp["Dy"]
        ) == pytest.approx(0.0455, rel=1e-2)

    @pytest.mark.component
    def test_numerical_issues(self, SolEx_frame):
        m = SolEx_frame
        dt = DiagnosticsToolbox(m)
        dt.assert_no_numerical_warnings()
