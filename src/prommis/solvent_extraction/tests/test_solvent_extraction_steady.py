#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import check_optimal_termination

from idaes.core.initialization import InitializationStatus
from idaes.core.scaling.util import jacobian_cond
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox

import pytest

from prommis.solvent_extraction.solvent_extraction import SolventExtractionInitializer
from prommis.solvent_extraction.solvent_extraction_steady import (
    main,
    model_buildup_and_set_inputs,
)
from prommis.util import assert_solution_equivalent

solver = get_solver()


class Test_Solvent_Extraction_steady_model:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        dosage = 5
        number_of_stages = 3
        m = model_buildup_and_set_inputs(dosage, number_of_stages, has_holdup=False)

        return m

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_initialization(self, SolEx_frame):
        model = SolEx_frame
        initializer = model.fs.solex.default_initializer()
        assert model.fs.solex.default_initializer is SolventExtractionInitializer
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

        assert jacobian_cond(model, scaled=False) == pytest.approx(2.46261e14, rel=1e-3)
        assert jacobian_cond(model, scaled=True) == pytest.approx(1.1119e7, rel=1e-3)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, SolEx_frame):

        model = SolEx_frame
        expected_results = {
            "organic_outlet.conc_mass_comp": {
                (0.0, "Kerosene"): (8.2000e05, 1e-4, None),
                (0.0, "DEHPA"): (4.6087e04, 1e-4, None),
                (0.0, "Al_o"): (2.2425e01, 1e-4, None),
                (0.0, "Ca_o"): (7.2059e00, 1e-4, None),
                (0.0, "Fe_o"): (1.0267e02, 1e-4, None),
                (0.0, "Sc_o"): (1.7633e00, 1e-4, None),
                (0.0, "Y_o"): (1.2402e-01, 1e-4, None),
                (0.0, "La_o"): (7.1891e-02, 1e-4, None),
                (0.0, "Ce_o"): (1.6778e-01, 1e-4, None),
                (0.0, "Pr_o"): (2.6727e-02, 1e-4, None),
                (0.0, "Nd_o"): (6.5178e-02, 1e-4, None),
                (0.0, "Sm_o"): (1.0438e-02, 1e-4, None),
                (0.0, "Gd_o"): (6.5974e-02, 1e-4, None),
                (0.0, "Dy_o"): (4.5838e-02, 1e-4, None),
            },
            "aqueous_outlet.conc_mass_comp": {
                (0.0, "H2O"): (1.000e06, 1e-4, None),
                (0.0, "H"): (3.9513e01, 1e-4, None),
                (0.0, "HSO4"): (8.0232e03, 1e-4, None),
                (0.0, "SO4"): (2.0564e03, 1e-4, None),
                (0.0, "Cl"): (1.0000e-07, 1e-4, None),
                (0.0, "Sc"): (2.7415e-03, 1e-4, None),
                (0.0, "Y"): (6.2927e-06, 1e-4, None),
                (0.0, "La"): (9.1421e-01, 1e-4, None),
                (0.0, "Ce"): (2.1095e00, 1e-4, None),
                (0.0, "Pr"): (2.7631e-01, 1e-4, None),
                (0.0, "Nd"): (8.8099e-01, 1e-4, None),
                (0.0, "Sm"): (8.6579e-02, 1e-4, None),
                (0.0, "Gd"): (1.9246e-01, 1e-4, None),
                (0.0, "Dy"): (1.1699e-03, 1e-4, None),
                (0.0, "Al"): (3.9995e02, 1e-4, None),
                (0.0, "Ca"): (1.0234e02, 1e-4, None),
                (0.0, "Fe"): (5.8559e02, 1e-4, None),
            },
        }
        assert_solution_equivalent(model.fs.solex, expected_results)

    @pytest.fixture(scope="class")
    def SolEx_total_flowsheet(self):
        dosage = 5
        number_of_stages = 3
        model, results = main(dosage, number_of_stages, has_holdup=False)

        return model, results

    @pytest.mark.component
    def test_solve_total(self, SolEx_total_flowsheet):
        # TODO this test is redundant with the above solve
        m, results = SolEx_total_flowsheet
        assert check_optimal_termination(results)

        dt = DiagnosticsToolbox(m)
        dt.assert_no_numerical_warnings()
        assert jacobian_cond(m, scaled=False) == pytest.approx(2.46261e14, rel=1e-3)
        assert jacobian_cond(m, scaled=True) == pytest.approx(1.1119e7, rel=1e-3)


class Test_Solvent_Extraction_steady_model_hydrostatic_pressure:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):
        dosage = 5
        number_of_stages = 3
        m = model_buildup_and_set_inputs(dosage, number_of_stages, has_holdup=True)

        return m

    @pytest.mark.component
    def test_structural_issues(self, SolEx_frame):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_initialization(self, SolEx_frame):
        model = SolEx_frame
        initializer = model.fs.solex.default_initializer()
        assert model.fs.solex.default_initializer is SolventExtractionInitializer
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

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            8.415018e12, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(1.2842e7, rel=1e-3)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, SolEx_frame):

        model = SolEx_frame
        expected_results = {
            "organic_outlet.conc_mass_comp": {
                (0.0, "Kerosene"): (8.2000e05, 1e-4, None),
                (0.0, "DEHPA"): (4.6087e04, 1e-4, None),
                (0.0, "Al_o"): (2.2425e01, 1e-4, None),
                (0.0, "Ca_o"): (7.2059e00, 1e-4, None),
                (0.0, "Fe_o"): (1.0267e02, 1e-4, None),
                (0.0, "Sc_o"): (1.7633e00, 1e-4, None),
                (0.0, "Y_o"): (1.2402e-01, 1e-4, None),
                (0.0, "La_o"): (7.1891e-02, 1e-4, None),
                (0.0, "Ce_o"): (1.6778e-01, 1e-4, None),
                (0.0, "Pr_o"): (2.6727e-02, 1e-4, None),
                (0.0, "Nd_o"): (6.5178e-02, 1e-4, None),
                (0.0, "Sm_o"): (1.0438e-02, 1e-4, None),
                (0.0, "Gd_o"): (6.5974e-02, 1e-4, None),
                (0.0, "Dy_o"): (4.5838e-02, 1e-4, None),
            },
            "aqueous_outlet.conc_mass_comp": {
                (0.0, "H2O"): (1.000e06, 1e-4, None),
                (0.0, "H"): (3.9513e01, 1e-4, None),
                (0.0, "HSO4"): (8.0232e03, 1e-4, None),
                (0.0, "SO4"): (2.0564e03, 1e-4, None),
                (0.0, "Cl"): (1.0000e-07, 1e-4, None),
                (0.0, "Sc"): (2.7415e-03, 1e-4, None),
                (0.0, "Y"): (6.2927e-06, 1e-4, None),
                (0.0, "La"): (9.1421e-01, 1e-4, None),
                (0.0, "Ce"): (2.1095e00, 1e-4, None),
                (0.0, "Pr"): (2.7631e-01, 1e-4, None),
                (0.0, "Nd"): (8.8099e-01, 1e-4, None),
                (0.0, "Sm"): (8.6579e-02, 1e-4, None),
                (0.0, "Gd"): (1.9246e-01, 1e-4, None),
                (0.0, "Dy"): (1.1699e-03, 1e-4, None),
                (0.0, "Al"): (3.9995e02, 1e-4, None),
                (0.0, "Ca"): (1.0234e02, 1e-4, None),
                (0.0, "Fe"): (5.8559e02, 1e-4, None),
            },
        }
        assert_solution_equivalent(model.fs.solex, expected_results)

    @pytest.fixture(scope="class")
    def SolEx_total_flowsheet(self):
        dosage = 5
        number_of_stages = 3
        model, results = main(dosage, number_of_stages, has_holdup=True)

        return model, results

    @pytest.mark.component
    def test_solve_total(self, SolEx_total_flowsheet):
        # TODO this test is redundant with the above solve
        m, results = SolEx_total_flowsheet
        assert check_optimal_termination(results)

        dt = DiagnosticsToolbox(m)
        dt.assert_no_numerical_warnings()
        assert jacobian_cond(m, scaled=False) == pytest.approx(8.415018e12, rel=1e-3)
        assert jacobian_cond(m, scaled=True) == pytest.approx(1.2842e7, rel=1e-3)
