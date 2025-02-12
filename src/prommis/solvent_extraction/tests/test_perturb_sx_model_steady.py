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
from idaes.core.util import DiagnosticsToolbox
from idaes.core.solvers import get_solver

import pytest

from prommis.solvent_extraction.solvent_extraction import SolventExtractionInitializer

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

        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(1.755)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(3999.818)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(693.459)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-8)

        m.fs.solex.aqueous_inlet.flow_vol.fix(62.01)

        m.fs.solex.organic_inlet.conc_mass_comp[0, "Al"].fix(1.267e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Ca"].fix(2.684e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Fe"].fix(2.873e-6)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Sc"].fix(1.734)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Y"].fix(2.179e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "La"].fix(0.000105)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Ce"].fix(0.00031)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Pr"].fix(3.711e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Nd"].fix(0.000165)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Sm"].fix(1.701e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Gd"].fix(3.357e-5)
        m.fs.solex.organic_inlet.conc_mass_comp[0, "Dy"].fix(8.008e-6)

        m.fs.solex.organic_inlet.flow_vol.fix(62.01)

        return m

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_initialization(self, SolEx_frame):
        model = SolEx_frame
        initializer = SolventExtractionInitializer()
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

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_ferroglobe_plant(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(7.943)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.391)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(1.297)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.546)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.497)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(1.278)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.173)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.191)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.098)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_kentucky_13_amd(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(220.761)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.802)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(48.190)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(10.65)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(71.197)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(45.272)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(78.684)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(60.564)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_kentucky_13_amd_5_dehpa_10_tbp(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(100)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(220.086)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.329)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(49.882)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(11.512)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(71.526)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(46.930)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(78.413)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(60.704)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_kentucky_13_amd_5_dehpa(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(100)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(218.61)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(3.219)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(49.157)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(11.235)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(70.224)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(44.943)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(78.651)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(61.095)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_kentucky_13_amd_5_cyanex(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(150)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(237.795)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.808)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(52.668)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(11.938)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(73.033)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(47.752)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(82.162)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(63.202)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_kentucky_13_distribution_isotherm(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(100)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(203.468)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.293)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(44.151)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(9.174)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(61.926)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(41.284)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(71.1)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(54.472)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def final_report_various_conditions(self, SolEx_frame):
        m = SolEx_frame
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(3.162)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.278)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.344)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(2.075)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(5.013)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.714)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(2.088)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.229)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.566)
        m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.080)

        results = solver.solve(m, tee=True)
        m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mol_comp["H2O"].pprint()

        # Check for optimal solution
        assert check_optimal_termination(results)
