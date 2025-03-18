#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, check_optimal_termination, value, units
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
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)

solver = get_solver()


class TestSXmodel:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.prop_o = REESolExOgParameters()
        m.fs.leach_soln = LeachSolutionParameters()
        m.fs.reaxn = SolventExtractionReactions()

        number_of_stages = 3
        stage_number = np.arange(1, number_of_stages + 1)
        dosage = 5

        m.fs.reaxn.extractant_dosage = dosage

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
            reaction_package=m.fs.reaxn,
            has_holdup=True,
        )

        m.fs.solex.mscontactor.volume[:].fix(0.4 * units.m**3)
        m.fs.solex.mscontactor.volume_frac_stream[:, :, "aqueous"].fix(0.5)
        m.fs.solex.cross_sec_area[:] = 1
        m.fs.solex.elevation[:] = 0

        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(10.75)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(100)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(1e4)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(422.375)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(109.542)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Cl"].fix(1e-7)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(688.266)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(0.032)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(0.124)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(0.986)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(2.277)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(0.303)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(0.946)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(0.097)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(0.2584)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

        m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Kerosene"].fix(
            820e3
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["DEHPA"].fix(
            975.8e3 * dosage / 100
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al_o"].fix(
            1.267e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca_o"].fix(
            2.684e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe_o"].fix(
            2.873e-6
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc_o"].fix(1.734)
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y_o"].fix(
            2.179e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La_o"].fix(
            0.000105
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce_o"].fix(
            0.00031
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr_o"].fix(
            3.711e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd_o"].fix(
            0.000165
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm_o"].fix(
            1.701e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd_o"].fix(
            3.357e-5
        )
        m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy_o"].fix(
            8.008e-6
        )

        m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

        m.fs.solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
        m.fs.solex.mscontactor.aqueous_inlet_state[:].temperature.fix(305.15 * units.K)
        m.fs.solex.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)
        m.fs.solex.mscontactor.organic_inlet_state[:].temperature.fix(305.15 * units.K)

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

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, SolEx_frame):

        model = SolEx_frame
        aqueous_outlet = {
            "H2O": 1000000,
            "H": 42.3638,
            "SO4": 1945.036,
            "HSO4": 8135.841,
            "Al": 399.2865,
            "Ca": 105.7483,
            "Cl": 9.9999e-8,
            "Ce": 2.1173,
            "Dy": 0.001774,
            "Fe": 509.07864,
            "Gd": 0.19972,
            "La": 0.91903,
            "Nd": 0.8842,
            "Pr": 0.27728,
            "Sc": 0.026568,
            "Sm": 0.08757,
            "Y": 8.3131e-06,
        }

        organic_outlet = {
            "Al_o": 23.0673,
            "Ca_o": 3.7901,
            "Ce_o": 0.15997,
            "DEHPA": 44793.2243,
            "Dy_o": 0.045233,
            "Fe_o": 179.1681,
            "Gd_o": 0.058708,
            "Kerosene": 820000,
            "La_o": 0.067067,
            "Nd_o": 0.061947,
            "Pr_o": 0.02575,
            "Sc_o": 1.73943,
            "Sm_o": 0.009440,
            "Y_o": 0.12401,
        }

        for k, v in model.fs.solex.organic_outlet.conc_mass_comp.items():
            assert value(v) == pytest.approx(organic_outlet[k[1]], rel=1e-4)

        for k, v in model.fs.solex.aqueous_outlet.conc_mass_comp.items():
            assert value(v) == pytest.approx(aqueous_outlet[k[1]], rel=1e-4)
