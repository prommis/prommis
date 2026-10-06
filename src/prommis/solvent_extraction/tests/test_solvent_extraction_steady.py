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

from prommis.properties.mixed_acid_properties import get_aliases

from prommis.solvent_extraction.solvent_extraction import (
    SolventExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_steady import (
    model_buildup_and_set_inputs,
)

solver = get_solver()

metal_list = ["La", "Y", "Pr", "Ce", "Nd", "Sm", "Gd", "Dy", "Al", "Ca", "Fe", "Sc"]


@pytest.mark.parametrize(
    "has_holdup",
    [False, True],
    scope="class",
)
@pytest.mark.parametrize(
    "use_mixed_acid",
    [False, True],
    scope="class",
)
class Test_Solvent_Extraction_steady_model:

    @pytest.fixture(scope="class")
    @classmethod
    def SolEx_frame(cls, has_holdup, use_mixed_acid):
        dosage = 5
        number_of_stages = 3
        m = model_buildup_and_set_inputs(
            dosage,
            number_of_stages,
            has_holdup=has_holdup,
            use_mixed_acid=use_mixed_acid,
        )
        if use_mixed_acid:
            # The sulfuric acid leaching properties use molecular weights
            # with fewer significant figures than the mixed acid properties,
            # which then leads to a discrepancy in the model solution.
            # Use versions with fewer significant figures here to make sure
            # that the models are structurally equivalent.
            m.fs.leach_soln.mw["H2O"] = 18e-3
            m.fs.leach_soln.mw["H_+"] = 1e-3
            m.fs.leach_soln.mw["HSO4_-"] = 97e-3
            m.fs.leach_soln.mw["SO4_2-"] = 96e-3

        return m

    @pytest.fixture(scope="class")
    @classmethod
    def expected_results(cls, has_holdup, use_mixed_acid):
        out = {
            "aqueous_inlet.flow_vol": {(0,): (62.01, 1e-4, None)},
            "aqueous_inlet.temperature": {(0,): (305.15, 1e-4, None)},
            "aqueous_inlet.pressure": {(0,): (101300.0, 1e-4, None)},
            "aqueous_inlet.conc_mass_comp": {
                (0.0, "H2O"): (1e6, 1e-4, None),
                (0.0, "H"): (10.75, 1e-4, None),
                (0.0, "HSO4"): (1e4, 1e-4, None),
                (0.0, "SO4"): (1e2, 1e-4, None),
                (0.0, "Cl"): (1e-7, 1e-4, None),
                (0.0, "Sc"): (3.2e-2, 1e-4, None),
                (0.0, "Y"): (0.124, 1e-4, None),
                (0.0, "La"): (0.986, 1e-4, None),
                (0.0, "Ce"): (2.277, 1e-4, None),
                (0.0, "Pr"): (0.303, 1e-4, None),
                (0.0, "Nd"): (0.946, 1e-4, None),
                (0.0, "Sm"): (9.7e-02, 1e-4, None),
                (0.0, "Gd"): (0.2584, 1e-4, None),
                (0.0, "Dy"): (4.7e-02, 1e-4, None),
                (0.0, "Al"): (4.22375e2, 1e-4, None),
                (0.0, "Ca"): (1.09542e2, 1e-4, None),
                (0.0, "Fe"): (6.88266e2, 1e-4, None),
            },
            "organic_inlet.flow_vol": {(0,): (62.01, 1e-4, None)},
            "organic_inlet.temperature": {(0,): (305.15, 1e-4, None)},
            "organic_inlet.pressure": {(0,): (101300.0, 1e-4, None)},
            "organic_inlet.conc_mass_comp": {
                (0.0, "Kerosene"): (8.2e5, 1e-4, None),
                (0.0, "DEHPA"): (4.879e4, 1e-4, None),
                (0.0, "Al_o"): (1.267e-5, 1e-4, None),
                (0.0, "Ca_o"): (2.684e-5, 1e-4, None),
                (0.0, "Fe_o"): (2.873e-6, 1e-4, None),
                (0.0, "Sc_o"): (1.734, 1e-4, None),
                (0.0, "Y_o"): (2.179e-5, 1e-4, None),
                (0.0, "La_o"): (1.05e-4, 1e-4, None),
                (0.0, "Ce_o"): (3.1e-4, 1e-4, None),
                (0.0, "Pr_o"): (3.711e-5, 1e-4, None),
                (0.0, "Nd_o"): (1.65e-4, 1e-4, None),
                (0.0, "Sm_o"): (1.701e-5, 1e-4, None),
                (0.0, "Gd_o"): (3.357e-5, 1e-4, None),
                (0.0, "Dy_o"): (8.008e-6, 1e-4, None),
            },
            "aqueous_outlet.flow_vol": {(0,): (62.01, 1e-4, None)},
            "aqueous_outlet.temperature": {(0,): (305.15, 1e-4, None)},
            "aqueous_outlet.pressure": {(0,): (101300.0, 1e-4, None)},
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
            "organic_outlet.flow_vol": {(0,): (62.01, 1e-4, None)},
            "organic_outlet.temperature": {(0,): (305.15, 1e-4, None)},
            "organic_outlet.pressure": {(0,): (101300.0, 1e-4, None)},
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
        }
        if has_holdup:
            # Uses hydrostatic pressure, not the inlet pressure
            out["aqueous_outlet.pressure"][(0,)] = (1.04895e05, 1e-4, None)
            out["organic_outlet.pressure"][(0,)] = (1.02933e05, 1e-4, None)

        if use_mixed_acid:
            for port in ["aqueous_inlet", "aqueous_outlet"]:
                old_dict = out[port + ".conc_mass_comp"]
                new_dict = {}
                for j1, j2 in get_aliases(include_sulfates=True).items():
                    new_dict[(0.0, j2)] = old_dict[(0.0, j1)]

                out[port + ".conc_mass_comp"] = new_dict

        return out

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
    def test_numerical_issues(self, SolEx_frame, has_holdup, use_mixed_acid):
        model = SolEx_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # Why does adding holdup *reduce* the unscaled condition number?
        if has_holdup:
            if use_mixed_acid:
                assert jacobian_cond(model, scaled=False) == pytest.approx(
                    8.4210e12, rel=1e-3
                )
                assert jacobian_cond(model, scaled=True) == pytest.approx(
                    1.3238e7, rel=1e-3
                )
            else:
                assert jacobian_cond(model, scaled=False) == pytest.approx(
                    8.415018e12, rel=1e-3
                )
                assert jacobian_cond(model, scaled=True) == pytest.approx(
                    1.2842e7, rel=1e-3
                )
        else:
            if use_mixed_acid:
                assert jacobian_cond(model, scaled=False) == pytest.approx(
                    2.4643e14, rel=1e-3
                )
                assert jacobian_cond(model, scaled=True) == pytest.approx(
                    1.1392e7, rel=1e-3
                )
            else:
                assert jacobian_cond(model, scaled=False) == pytest.approx(
                    2.46261e14, rel=1e-3
                )
                assert jacobian_cond(model, scaled=True) == pytest.approx(
                    1.1119e7, rel=1e-3
                )

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, SolEx_frame, expected_results):

        model = SolEx_frame
        assert_solution_equivalent(model.fs.solex, expected_results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_get_stream_table_contents(
        self, SolEx_frame, expected_results, use_mixed_acid
    ):
        nan = float("NaN")
        aq_components = [
            "H2O",
            "H",
            "HSO4",
            "SO4",
            "Cl",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Gd",
            "Dy",
            "Al",
            "Ca",
            "Fe",
        ]
        if use_mixed_acid:
            aliases = get_aliases(include_sulfates=True)
            for k, j in enumerate(aq_components):
                aq_components[k] = aliases[j]
        org_components = [
            "Kerosene",
            "DEHPA",
            "Al_o",
            "Ca_o",
            "Fe_o",
            "Sc_o",
            "Y_o",
            "La_o",
            "Ce_o",
            "Pr_o",
            "Nd_o",
            "Sm_o",
            "Gd_o",
            "Dy_o",
        ]
        expected = {
            "Units": {
                "flow_vol": getattr(pyunits.pint_registry, "meter ** 3 / second"),
                "temperature": getattr(pyunits.pint_registry, "kelvin"),
                "pressure": getattr(pyunits.pint_registry, "pascal"),
            }
        }
        for j in aq_components + org_components:
            expected["Units"][f"conc_mass_comp {j}"] = getattr(
                pyunits.pint_registry, "kilogram / meter ** 3"
            )

        for phase in ["aqueous", "organic"]:
            for direction in ["inlet", "outlet"]:
                # e.g. port_name = "aqueous_inlet."
                port_name = phase + "_" + direction + "."
                stream_info = {}
                # flow_vol is reported in m**3/s
                stream_info["flow_vol"] = (
                    expected_results[port_name + "flow_vol"][(0,)][0] / 3600 / 1e3
                )
                stream_info["temperature"] = expected_results[
                    port_name + "temperature"
                ][(0,)][0]
                stream_info["pressure"] = expected_results[port_name + "pressure"][
                    (0,)
                ][0]
                conc_mass_comp = expected_results[port_name + "conc_mass_comp"]
                if phase == "aqueous":
                    for j in aq_components:
                        # Concentration is reported in kg/m**3
                        stream_info[f"conc_mass_comp {j}"] = (
                            conc_mass_comp[(0, j)][0] * 1e-3
                        )
                    for j in org_components:
                        stream_info[f"conc_mass_comp {j}"] = nan
                else:
                    for j in org_components:
                        # Concentration is reported in kg/m**3
                        stream_info[f"conc_mass_comp {j}"] = (
                            conc_mass_comp[(0, j)][0] * 1e-3
                        )
                    for j in aq_components:
                        stream_info[f"conc_mass_comp {j}"] = nan
                expected[phase + " " + direction.capitalize()] = stream_info

        out = SolEx_frame.fs.solex._get_stream_table_contents()

        pd.testing.assert_frame_equal(
            pd.DataFrame(expected),
            out,
            rtol=1e-4,
            atol=1e-12,
            # Ignore ordering since we're constructing the dataframe
            # from a dictionary
            check_like=True,
        )

    @pytest.mark.component
    @pytest.mark.solver
    def test_get_performance_contents(self, SolEx_frame, has_holdup):
        unit = SolEx_frame.fs.solex
        out = unit._get_performance_contents()
        assert len(out) == 3

        if has_holdup:
            assert len(out["vars"]) == 3
            for e in unit.mscontactor.elements:
                assert (
                    out["vars"][f"Aqueous phase frac stage {e}"]
                    is unit.mscontactor.volume_frac_stream[0, e, "aqueous"]
                )

            assert len(out["params"]) == 2
            assert out["params"]["Stage base area"] is unit.area_cross_stage
            assert out["params"]["Elevation"] is unit.elevation
        else:
            assert len(out["vars"]) == 0
            assert len(out["params"]) == 0

        out = out["exprs"]
        assert len(out) == len(metal_list)

        # Validate distribution coefficients
        for j in metal_list:
            D_mean = 1
            for e in unit.mscontactor.elements:
                D_mean *= value(
                    unit.mscontactor.heterogeneous_reactions[
                        0, e
                    ].distribution_coefficient[j]
                )
            D_mean = D_mean ** (1 / len(unit.mscontactor.elements))
            assert value(
                out[f"Geometric mean distribution coefficient {j}"]
            ) == pytest.approx(D_mean, rel=1e-12, abs=1e-12)
