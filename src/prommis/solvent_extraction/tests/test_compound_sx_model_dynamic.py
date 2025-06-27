#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import os
from pyomo.environ import check_optimal_termination, value

from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox, from_json

import pytest

from prommis.solvent_extraction.compound_sx_flowsheet_dynamic import (
    build_model,
    discretization_scheme,
    set_initial_guess,
    set_inputs,
    copy_first_steady_state,
)

solver = get_solver()


class TestSXdynamicmodel:
    @pytest.fixture(scope="class")
    def Comp_SolEx_frame(self):
        dosage = 5
        number_of_stages = 3
        time_duration = 12
        perturb_time = 4

        m = build_model(dosage, number_of_stages, time_duration)
        discretization_scheme(m)
        current_directory = os.path.dirname(__file__)
        parent_directory = os.path.dirname(current_directory)
        json_file_path = os.path.join(
            parent_directory, "compound_solvent_extraction.json"
        )
        from_json(m, fname=json_file_path)
        copy_first_steady_state(m)
        set_inputs(m, dosage, perturb_time)
        set_initial_guess(m)
        return m

    @pytest.mark.component
    def test_structural_issues(self, Comp_SolEx_frame):
        model = Comp_SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings(ignore_unit_consistency=True)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, Comp_SolEx_frame):
        m = Comp_SolEx_frame
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, Comp_SolEx_frame):
        model = Comp_SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, Comp_SolEx_frame):

        model = Comp_SolEx_frame
        time_duration = 12
        aqueous_outlet = {
            "H2O": 1000000,
            "H": 39.47196,
            "SO4": 2058.1002,
            "HSO4": 8021.5031,
            "Al": 400.13261,
            "Ca": 102.39299,
            "Cl": 6.3421e-7,
            "Ce": 2.10673,
            "Dy": 0.0010188,
            "Fe": 586.26075,
            "Gd": 0.19143,
            "La": 0.91468,
            "Nd": 0.88054,
            "Pr": 0.27649,
            "Sc": 0.0027415,
            "Sm": 0.086752,
            "Y": 4.67305e-06,
        }

        organic_outlet = {
            "Al_o": 22.42494,
            "Ca_o": 7.20595,
            "Ce_o": 0.17169,
            "DEHPA": 46086.679,
            "Dy_o": 0.046208,
            "Fe_o": 102.67241,
            "Gd_o": 0.06724,
            "Kerosene": 820000,
            "La_o": 0.071891,
            "Nd_o": 0.066063,
            "Pr_o": 0.026727,
            "Sc_o": 1.76353,
            "Sm_o": 0.010324,
            "Y_o": 0.12511,
        }

        for k, v in model.fs.compound_solex.organic_outlet.conc_mass_comp.items():
            if k[0] == time_duration:
                assert value(v) == pytest.approx(organic_outlet[k[1]], rel=1e-4)

        for k, v in model.fs.compound_solex.aqueous_outlet.conc_mass_comp.items():
            if k[0] == time_duration:
                assert value(v) == pytest.approx(aqueous_outlet[k[1]], rel=1e-4)
