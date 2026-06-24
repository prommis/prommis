#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import pandas as pd

from pyomo.environ import check_optimal_termination, units as pyunits, value

from idaes.core.initialization import InitializationStatus
from idaes.core.scaling.util import jacobian_cond
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.testing import assert_solution_equivalent

import pytest

from prommis.solvent_extraction.solvent_extraction import (
    ExtractionDirection,
    SolventExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_steady import (
    model_buildup_and_set_inputs,
)

solver = get_solver()


def _validate_distribution_coefficients(unit, m_list, expr_dict):
    for j in m_list:
        D_mean = 1
        for e in unit.mscontactor.elements:
            D_mean *= value(
                unit.mscontactor.heterogeneous_reactions[0, e].distribution_coefficient[
                    j
                ]
            )
        D_mean = D_mean ** (1 / len(unit.mscontactor.elements))
        assert value(
            expr_dict[f"Geometric mean distribution coefficient {j}"]
        ) == pytest.approx(D_mean, rel=1e-12, abs=1e-12)


metal_list = ["La", "Y", "Pr", "Ce", "Nd", "Sm", "Gd", "Dy", "Al", "Ca", "Fe", "Sc"]


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

    def test_get_stream_table_contents(self, SolEx_frame):
        # TODO this is essentially the same data as in the previous test.
        # There should be a single data structure that can cover both tests.
        nan = float("NaN")
        expected = {
            "Units": {
                "flow_vol": getattr(pyunits.pint_registry, "meter ** 3 / second"),
                "conc_mass_comp H2O": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp H": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp HSO4": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp SO4": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Cl": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sc": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Y": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp La": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ce": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Pr": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Nd": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sm": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Gd": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Dy": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Al": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ca": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Fe": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "temperature": getattr(pyunits.pint_registry, "kelvin"),
                "pressure": getattr(pyunits.pint_registry, "pascal"),
                "conc_mass_comp Kerosene": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp DEHPA": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Al_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ca_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Fe_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sc_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Y_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp La_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ce_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Pr_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Nd_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sm_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Gd_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Dy_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
            },
            "aqueous Inlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": 999.9999999999998,
                "conc_mass_comp H": 0.010749999999999997,
                "conc_mass_comp HSO4": 9.999999999999998,
                "conc_mass_comp SO4": 0.09999999999999998,
                "conc_mass_comp Cl": 9.999999999999998e-11,
                "conc_mass_comp Sc": 3.199999999999999e-05,
                "conc_mass_comp Y": 0.00012399999999999998,
                "conc_mass_comp La": 0.000986,
                "conc_mass_comp Ce": 0.0022769999999999995,
                "conc_mass_comp Pr": 0.00030299999999999994,
                "conc_mass_comp Nd": 0.0009459999999999998,
                "conc_mass_comp Sm": 9.699999999999999e-05,
                "conc_mass_comp Gd": 0.0002584,
                "conc_mass_comp Dy": 4.699999999999999e-05,
                "conc_mass_comp Al": 0.42237499999999994,
                "conc_mass_comp Ca": 0.10954199999999999,
                "conc_mass_comp Fe": 0.6882659999999998,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": nan,
                "conc_mass_comp DEHPA": nan,
                "conc_mass_comp Al_o": nan,
                "conc_mass_comp Ca_o": nan,
                "conc_mass_comp Fe_o": nan,
                "conc_mass_comp Sc_o": nan,
                "conc_mass_comp Y_o": nan,
                "conc_mass_comp La_o": nan,
                "conc_mass_comp Ce_o": nan,
                "conc_mass_comp Pr_o": nan,
                "conc_mass_comp Nd_o": nan,
                "conc_mass_comp Sm_o": nan,
                "conc_mass_comp Gd_o": nan,
                "conc_mass_comp Dy_o": nan,
            },
            "aqueous Outlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": 999.9999999999998,
                "conc_mass_comp H": 0.039513087794080776,
                "conc_mass_comp HSO4": 8.023222048810739,
                "conc_mass_comp SO4": 2.0563987970532893,
                "conc_mass_comp Cl": 9.999999999999987e-11,
                "conc_mass_comp Sc": 2.7415124943096402e-06,
                "conc_mass_comp Y": 6.292748432495276e-09,
                "conc_mass_comp La": 0.0009142143743736765,
                "conc_mass_comp Ce": 0.0021095278998357163,
                "conc_mass_comp Pr": 0.0002763100455024723,
                "conc_mass_comp Nd": 0.0008809872254832235,
                "conc_mass_comp Sm": 8.657862645327003e-05,
                "conc_mass_comp Gd": 0.00019245933405617006,
                "conc_mass_comp Dy": 1.1699115491985285e-06,
                "conc_mass_comp Al": 0.39995007985144854,
                "conc_mass_comp Ca": 0.1023360816022842,
                "conc_mass_comp Fe": 0.5855947170037661,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": nan,
                "conc_mass_comp DEHPA": nan,
                "conc_mass_comp Al_o": nan,
                "conc_mass_comp Ca_o": nan,
                "conc_mass_comp Fe_o": nan,
                "conc_mass_comp Sc_o": nan,
                "conc_mass_comp Y_o": nan,
                "conc_mass_comp La_o": nan,
                "conc_mass_comp Ce_o": nan,
                "conc_mass_comp Pr_o": nan,
                "conc_mass_comp Nd_o": nan,
                "conc_mass_comp Sm_o": nan,
                "conc_mass_comp Gd_o": nan,
                "conc_mass_comp Dy_o": nan,
            },
            "organic Inlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": nan,
                "conc_mass_comp H": nan,
                "conc_mass_comp HSO4": nan,
                "conc_mass_comp SO4": nan,
                "conc_mass_comp Cl": nan,
                "conc_mass_comp Sc": nan,
                "conc_mass_comp Y": nan,
                "conc_mass_comp La": nan,
                "conc_mass_comp Ce": nan,
                "conc_mass_comp Pr": nan,
                "conc_mass_comp Nd": nan,
                "conc_mass_comp Sm": nan,
                "conc_mass_comp Gd": nan,
                "conc_mass_comp Dy": nan,
                "conc_mass_comp Al": nan,
                "conc_mass_comp Ca": nan,
                "conc_mass_comp Fe": nan,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": 819.9999999999999,
                "conc_mass_comp DEHPA": 48.78999999999999,
                "conc_mass_comp Al_o": 1.2669999999999999e-08,
                "conc_mass_comp Ca_o": 2.6839999999999994e-08,
                "conc_mass_comp Fe_o": 2.8729999999999995e-09,
                "conc_mass_comp Sc_o": 0.0017339999999999996,
                "conc_mass_comp Y_o": 2.1789999999999995e-08,
                "conc_mass_comp La_o": 1.0499999999999999e-07,
                "conc_mass_comp Ce_o": 3.0999999999999994e-07,
                "conc_mass_comp Pr_o": 3.710999999999999e-08,
                "conc_mass_comp Nd_o": 1.6499999999999996e-07,
                "conc_mass_comp Sm_o": 1.7009999999999998e-08,
                "conc_mass_comp Gd_o": 3.356999999999999e-08,
                "conc_mass_comp Dy_o": 8.007999999999998e-09,
            },
            "organic Outlet": {
                "flow_vol": 1.7225000000000002e-05,
                "conc_mass_comp H2O": nan,
                "conc_mass_comp H": nan,
                "conc_mass_comp HSO4": nan,
                "conc_mass_comp SO4": nan,
                "conc_mass_comp Cl": nan,
                "conc_mass_comp Sc": nan,
                "conc_mass_comp Y": nan,
                "conc_mass_comp La": nan,
                "conc_mass_comp Ce": nan,
                "conc_mass_comp Pr": nan,
                "conc_mass_comp Nd": nan,
                "conc_mass_comp Sm": nan,
                "conc_mass_comp Gd": nan,
                "conc_mass_comp Dy": nan,
                "conc_mass_comp Al": nan,
                "conc_mass_comp Ca": nan,
                "conc_mass_comp Fe": nan,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": 819.9999999999999,
                "conc_mass_comp DEHPA": 46.08675988668226,
                "conc_mass_comp Al_o": 0.022424932818551315,
                "conc_mass_comp Ca_o": 0.007205945237715807,
                "conc_mass_comp Fe_o": 0.10267128586923359,
                "conc_mass_comp Sc_o": 0.0017632584875056905,
                "conc_mass_comp Y_o": 0.0001240154972515675,
                "conc_mass_comp La_o": 7.189062562632331e-05,
                "conc_mass_comp Ce_o": 0.0001677821001642839,
                "conc_mass_comp Pr_o": 2.67270644975277e-05,
                "conc_mass_comp Nd_o": 6.517777451677624e-05,
                "conc_mass_comp Sm_o": 1.043838354672995e-05,
                "conc_mass_comp Gd_o": 6.59742359438299e-05,
                "conc_mass_comp Dy_o": 4.5838096450801465e-05,
            },
        }

        out = SolEx_frame.fs.solex._get_stream_table_contents()

        pd.testing.assert_frame_equal(
            pd.DataFrame(expected), out, rtol=1e-4, atol=1e-12
        )

    def test_get_performance_contents(self, SolEx_frame):
        unit = SolEx_frame.fs.solex

        out = unit._get_performance_contents()
        assert len(out) == 3
        assert len(out["vars"]) == 0
        assert len(out["params"]) == 0

        out = out["exprs"]
        assert len(out) == len(metal_list)

        _validate_distribution_coefficients(unit, metal_list, out)


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

    def test_get_stream_table_contents(self, SolEx_frame):
        # TODO this is essentially the same data as in the previous test.
        # There should be a single data structure that can cover both tests.
        nan = float("NaN")
        expected = {
            "Units": {
                "flow_vol": getattr(pyunits.pint_registry, "meter ** 3 / second"),
                "conc_mass_comp H2O": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp H": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp HSO4": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp SO4": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Cl": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sc": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Y": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp La": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ce": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Pr": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Nd": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sm": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Gd": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Dy": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Al": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ca": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Fe": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "temperature": getattr(pyunits.pint_registry, "kelvin"),
                "pressure": getattr(pyunits.pint_registry, "pascal"),
                "conc_mass_comp Kerosene": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp DEHPA": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Al_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ca_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Fe_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sc_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Y_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp La_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Ce_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Pr_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Nd_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Sm_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Gd_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
                "conc_mass_comp Dy_o": getattr(
                    pyunits.pint_registry, "kilogram / meter ** 3"
                ),
            },
            "aqueous Inlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": 999.9999999999998,
                "conc_mass_comp H": 0.010749999999999997,
                "conc_mass_comp HSO4": 9.999999999999998,
                "conc_mass_comp SO4": 0.09999999999999998,
                "conc_mass_comp Cl": 9.999999999999998e-11,
                "conc_mass_comp Sc": 3.199999999999999e-05,
                "conc_mass_comp Y": 0.00012399999999999998,
                "conc_mass_comp La": 0.000986,
                "conc_mass_comp Ce": 0.0022769999999999995,
                "conc_mass_comp Pr": 0.00030299999999999994,
                "conc_mass_comp Nd": 0.0009459999999999998,
                "conc_mass_comp Sm": 9.699999999999999e-05,
                "conc_mass_comp Gd": 0.0002584,
                "conc_mass_comp Dy": 4.699999999999999e-05,
                "conc_mass_comp Al": 0.42237499999999994,
                "conc_mass_comp Ca": 0.10954199999999999,
                "conc_mass_comp Fe": 0.6882659999999998,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": nan,
                "conc_mass_comp DEHPA": nan,
                "conc_mass_comp Al_o": nan,
                "conc_mass_comp Ca_o": nan,
                "conc_mass_comp Fe_o": nan,
                "conc_mass_comp Sc_o": nan,
                "conc_mass_comp Y_o": nan,
                "conc_mass_comp La_o": nan,
                "conc_mass_comp Ce_o": nan,
                "conc_mass_comp Pr_o": nan,
                "conc_mass_comp Nd_o": nan,
                "conc_mass_comp Sm_o": nan,
                "conc_mass_comp Gd_o": nan,
                "conc_mass_comp Dy_o": nan,
            },
            "aqueous Outlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": 999.9999999999998,
                "conc_mass_comp H": 0.039513087794080776,
                "conc_mass_comp HSO4": 8.023222048810739,
                "conc_mass_comp SO4": 2.0563987970532893,
                "conc_mass_comp Cl": 9.999999999999987e-11,
                "conc_mass_comp Sc": 2.7415124943096402e-06,
                "conc_mass_comp Y": 6.292748432495276e-09,
                "conc_mass_comp La": 0.0009142143743736765,
                "conc_mass_comp Ce": 0.0021095278998357163,
                "conc_mass_comp Pr": 0.0002763100455024723,
                "conc_mass_comp Nd": 0.0008809872254832235,
                "conc_mass_comp Sm": 8.657862645327003e-05,
                "conc_mass_comp Gd": 0.00019245933405617006,
                "conc_mass_comp Dy": 1.1699115491985285e-06,
                "conc_mass_comp Al": 0.39995007985144854,
                "conc_mass_comp Ca": 0.1023360816022842,
                "conc_mass_comp Fe": 0.5855947170037661,
                "temperature": 305.15,
                "pressure": 104894.6206,
                "conc_mass_comp Kerosene": nan,
                "conc_mass_comp DEHPA": nan,
                "conc_mass_comp Al_o": nan,
                "conc_mass_comp Ca_o": nan,
                "conc_mass_comp Fe_o": nan,
                "conc_mass_comp Sc_o": nan,
                "conc_mass_comp Y_o": nan,
                "conc_mass_comp La_o": nan,
                "conc_mass_comp Ce_o": nan,
                "conc_mass_comp Pr_o": nan,
                "conc_mass_comp Nd_o": nan,
                "conc_mass_comp Sm_o": nan,
                "conc_mass_comp Gd_o": nan,
                "conc_mass_comp Dy_o": nan,
            },
            "organic Inlet": {
                "flow_vol": 1.7225000000000005e-05,
                "conc_mass_comp H2O": nan,
                "conc_mass_comp H": nan,
                "conc_mass_comp HSO4": nan,
                "conc_mass_comp SO4": nan,
                "conc_mass_comp Cl": nan,
                "conc_mass_comp Sc": nan,
                "conc_mass_comp Y": nan,
                "conc_mass_comp La": nan,
                "conc_mass_comp Ce": nan,
                "conc_mass_comp Pr": nan,
                "conc_mass_comp Nd": nan,
                "conc_mass_comp Sm": nan,
                "conc_mass_comp Gd": nan,
                "conc_mass_comp Dy": nan,
                "conc_mass_comp Al": nan,
                "conc_mass_comp Ca": nan,
                "conc_mass_comp Fe": nan,
                "temperature": 305.15,
                "pressure": 101300.0,
                "conc_mass_comp Kerosene": 819.9999999999999,
                "conc_mass_comp DEHPA": 48.78999999999999,
                "conc_mass_comp Al_o": 1.2669999999999999e-08,
                "conc_mass_comp Ca_o": 2.6839999999999994e-08,
                "conc_mass_comp Fe_o": 2.8729999999999995e-09,
                "conc_mass_comp Sc_o": 0.0017339999999999996,
                "conc_mass_comp Y_o": 2.1789999999999995e-08,
                "conc_mass_comp La_o": 1.0499999999999999e-07,
                "conc_mass_comp Ce_o": 3.0999999999999994e-07,
                "conc_mass_comp Pr_o": 3.710999999999999e-08,
                "conc_mass_comp Nd_o": 1.6499999999999996e-07,
                "conc_mass_comp Sm_o": 1.7009999999999998e-08,
                "conc_mass_comp Gd_o": 3.356999999999999e-08,
                "conc_mass_comp Dy_o": 8.007999999999998e-09,
            },
            "organic Outlet": {
                "flow_vol": 1.7225000000000002e-05,
                "conc_mass_comp H2O": nan,
                "conc_mass_comp H": nan,
                "conc_mass_comp HSO4": nan,
                "conc_mass_comp SO4": nan,
                "conc_mass_comp Cl": nan,
                "conc_mass_comp Sc": nan,
                "conc_mass_comp Y": nan,
                "conc_mass_comp La": nan,
                "conc_mass_comp Ce": nan,
                "conc_mass_comp Pr": nan,
                "conc_mass_comp Nd": nan,
                "conc_mass_comp Sm": nan,
                "conc_mass_comp Gd": nan,
                "conc_mass_comp Dy": nan,
                "conc_mass_comp Al": nan,
                "conc_mass_comp Ca": nan,
                "conc_mass_comp Fe": nan,
                "temperature": 305.15,
                "pressure": 102933.2906,
                "conc_mass_comp Kerosene": 819.9999999999999,
                "conc_mass_comp DEHPA": 46.08675988668226,
                "conc_mass_comp Al_o": 0.022424932818551315,
                "conc_mass_comp Ca_o": 0.007205945237715807,
                "conc_mass_comp Fe_o": 0.10267128586923359,
                "conc_mass_comp Sc_o": 0.0017632584875056905,
                "conc_mass_comp Y_o": 0.0001240154972515675,
                "conc_mass_comp La_o": 7.189062562632331e-05,
                "conc_mass_comp Ce_o": 0.0001677821001642839,
                "conc_mass_comp Pr_o": 2.67270644975277e-05,
                "conc_mass_comp Nd_o": 6.517777451677624e-05,
                "conc_mass_comp Sm_o": 1.043838354672995e-05,
                "conc_mass_comp Gd_o": 6.59742359438299e-05,
                "conc_mass_comp Dy_o": 4.5838096450801465e-05,
            },
        }

        out = SolEx_frame.fs.solex._get_stream_table_contents()

        pd.testing.assert_frame_equal(
            pd.DataFrame(expected), out, rtol=1e-4, atol=1e-12
        )

    def test_get_performance_contents(self, SolEx_frame):
        unit = SolEx_frame.fs.solex

        out = unit._get_performance_contents()
        assert len(out) == 3
        assert len(out["vars"]) == 3
        for e in unit.mscontactor.elements:
            assert (
                out["vars"][f"Aqueous phase frac stage {e}"]
                is unit.mscontactor.volume_frac_stream[0, e, "aqueous"]
            )

        assert len(out["params"]) == 2
        assert out["params"]["Stage base area"] is unit.area_cross_stage
        assert out["params"]["Elevation"] is unit.elevation

        out = out["exprs"]
        assert len(out) == len(metal_list)

        _validate_distribution_coefficients(unit, metal_list, out)
