#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import os
from pyomo.environ import check_optimal_termination

from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox

import pytest

from prommis.solvent_extraction.mixer_settler_ex_flowsheet_dynamic import (
    build_model_and_discretize,
    initialize_set_input_and_initial_conditions,
    import_steady_value,
    main,
)
from prommis.util import assert_solution_equivalent

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
        expected_results = {
            "organic_outlet.conc_mass_comp": {
                (12, "Kerosene"): (8.2000e05, 1e-4, None),
                (12, "DEHPA"): (4.6087e04, 1e-4, None),
                (12, "Al_o"): (2.2425e01, 1e-4, None),
                (12, "Ca_o"): (7.2059e00, 1e-4, None),
                (12, "Fe_o"): (1.0267e02, 1e-4, None),
                (12, "Sc_o"): (1.76417e00, 1e-4, None),
                (12, "Y_o"): (1.27672e-01, 1e-4, None),
                (12, "La_o"): (7.1891e-02, 1e-4, None),
                (12, "Ce_o"): (1.6778e-01, 1e-4, None),
                (12, "Pr_o"): (2.6727e-02, 1e-4, None),
                (12, "Nd_o"): (6.5178e-02, 1e-4, None),
                (12, "Sm_o"): (1.04399e-02, 1e-4, None),
                (12, "Gd_o"): (6.60131e-02, 1e-4, None),
                (12, "Dy_o"): (4.66250e-02, 1e-4, None),
            },
            "aqueous_outlet.conc_mass_comp": {
                (12, "H2O"): (1.000e06, 1e-4, None),
                (12, "H"): (3.93625e01, 1e-4, None),
                (12, "HSO4"): (8.01691e03, 1e-4, None),
                (12, "SO4"): (2.06264e03, 1e-4, None),
                (12, "Cl"): (1.0000e-07, 1e-4, 1e-6),
                (12, "Sc"): (2.7415e-03, 1e-4, None),
                (12, "Y"): (8.81556e-06, 1e-4, None),
                (12, "La"): (9.15961e-01, 1e-4, None),
                (12, "Ce"): (2.11362e00, 1e-4, None),
                (12, "Pr"): (2.76983e-01, 1e-4, None),
                (12, "Nd"): (8.82605e-01, 1e-4, None),
                (12, "Sm"): (8.68056e-02, 1e-4, None),
                (12, "Gd"): (1.93454e-01, 1e-4, None),
                (12, "Dy"): (1.19306e-03, 1e-4, None),
                (12, "Al"): (4.00601e02, 1e-4, None),
                (12, "Ca"): (1.02540e02, 1e-4, None),
                (12, "Fe"): (5.88068e02, 1e-4, None),
            },
        }
        assert_solution_equivalent(model.fs.mixer_settler_ex, expected_results)

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
