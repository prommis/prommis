#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import os
from pyomo.environ import check_optimal_termination, value

from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox

import pytest

from prommis.solvent_extraction.mixer_settler_ex_flowsheet_dynamic import (
    build_model_and_discretize,
    initialize_set_input_and_initial_conditions,
    import_steady_value,
    main,
)

solver = get_solver()


class Test_Mixer_Settler_Ex_dynamic_model:
    @pytest.fixture(scope="class")
    def Mix_Settle_Ex_frame(self):
        dosage = 5
        number_of_stages = 3
        time_duration = 12
        perturb_time = 4

        m = build_model_and_discretize(dosage, number_of_stages, time_duration)
        current_directory = os.path.dirname(__file__)
        parent_directory = os.path.dirname(current_directory)
        json_file_path = os.path.join(parent_directory, "mixer_settler_extraction.json")
        import_steady_value(m, json_file_path)
        initialize_set_input_and_initial_conditions(m, dosage, perturb_time)

        return m

    @pytest.mark.component
    def test_structural_issues(self, Mix_Settle_Ex_frame):
        model = Mix_Settle_Ex_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings(ignore_unit_consistency=True)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, Mix_Settle_Ex_frame):
        m = Mix_Settle_Ex_frame
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, Mix_Settle_Ex_frame):
        model = Mix_Settle_Ex_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, Mix_Settle_Ex_frame):

        model = Mix_Settle_Ex_frame
        time_duration = 12
        aqueous_outlet = {
            "H2O": 1000000,
            "H": 39.3626,
            "SO4": 2062.6383,
            "HSO4": 8016.9174,
            "Al": 400.6008,
            "Ca": 102.5401,
            "Cl": 7.23224e-07,
            "Ce": 2.10975,
            "Dy": 1.03533e-03,
            "Fe": 588.068,
            "Gd": 0.19219,
            "La": 0.91596,
            "Nd": 0.88173,
            "Pr": 0.27698,
            "Sc": 2.7415e-03,
            "Sm": 8.6918e-02,
            "Y": 9.91505e-06,
        }

        organic_outlet = {
            "Al_o": 22.42494,
            "Ca_o": 7.20595,
            "Ce_o": 0.17169,
            "DEHPA": 46086.679,
            "Dy_o": 4.68083e-02,
            "Fe_o": 102.67241,
            "Gd_o": 6.72757e-02,
            "Kerosene": 820000,
            "La_o": 7.1891e-02,
            "Nd_o": 6.6063e-02,
            "Pr_o": 2.6727e-02,
            "Sc_o": 1.76417,
            "Sm_o": 1.03253e-02,
            "Y_o": 0.127708,
        }

        for k, v in model.fs.mixer_settler_ex.organic_outlet.conc_mass_comp.items():
            if k[0] == time_duration:
                assert value(v) == pytest.approx(organic_outlet[k[1]], rel=1e-4)

        for k, v in model.fs.mixer_settler_ex.aqueous_outlet.conc_mass_comp.items():
            if k[0] == time_duration:
                assert value(v) == pytest.approx(aqueous_outlet[k[1]], rel=1e-4)

    @pytest.fixture(scope="class")
    def Mix_Settle_Ex_total_flowsheet(self):
        dosage = 5
        number_of_stages = 3
        time_duration = 12
        perturb_time = 4
        current_directory = os.path.dirname(__file__)
        parent_directory = os.path.dirname(current_directory)
        json_file_path = os.path.join(parent_directory, "mixer_settler_extraction.json")

        m, results = main(
            dosage, number_of_stages, time_duration, perturb_time, json_file_path
        )

        return m, results

    @pytest.mark.component
    def test_solve_total(self, Mix_Settle_Ex_total_flowsheet):
        m, results = Mix_Settle_Ex_total_flowsheet
        assert check_optimal_termination(results)
