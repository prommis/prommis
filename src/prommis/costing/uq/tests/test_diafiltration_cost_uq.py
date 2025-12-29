#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

r"""
Author: Lingyan Deng, Brandon Paul

Tests for diafiltration_cost_UQ.

This test suite checks that:

1. The seiving coefficient are reset value for two different value.

2. The uncertain costing parameters and their distributions are as expected,
   and are consistent with the model (e.g., lognormal specs reproduce the
   nominal mean from the model parameters).

3. The diafiltration flowsheet and QGESS costing structure/constraints are
   present (adapted from `prommis.nanofiltration.diafiltration` and QGESS
   costing).  
"""

import math

import pyomo.environ as pyo

from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import numpy as np
import pytest

import prommis.costing.uq.diafiltration_cost_uq as uq
from prommis.costing.uq.diafiltration_cost_uq import (
    analyze_sensitivity,
    build_diafiltration_model,
    build_uncertainty_specs,
    decision_variables_bounds,
    estimate_lognormal_params_from_data,
    identify_uncertain_params,
    lhs_unit,
    main,
    plot_distributions_by_technology,
    triangular_icdf,
)
from prommis.uky.costing.ree_plant_capcost import QGESSCosting


class TestDiafiltrationCostUQStructure:
    @pytest.fixture(scope="class")
    def model(self):
        """Build a fresh diafiltration UQ model for all tests."""
        m = build_diafiltration_model(sieving_coeffs=(1.3, 0.5), technology_name=None)
        return m

    # 0. Sieving coefficient options (technology comparison)
    @pytest.mark.unit
    def test_two_sieving_coeff_options_are_applied(self):
        """Build two models with different (Li, Co) sieving coefficients and verify they differ."""
        m_a = build_diafiltration_model(
            sieving_coeffs=(1.3, 0.5), technology_name="tech_a"
        )
        m_b = build_diafiltration_model(
            sieving_coeffs=(1.5, 0.8), technology_name="tech_b"
        )

        # Both should be fixed (build_model fixes defaults; UQ script overrides and fixes)
        assert m_a.fs.sieving_coefficient["Li"].fixed
        assert m_a.fs.sieving_coefficient["Co"].fixed
        assert m_b.fs.sieving_coefficient["Li"].fixed
        assert m_b.fs.sieving_coefficient["Co"].fixed

        # Values should match requested options
        assert pyo.value(m_a.fs.sieving_coefficient["Li"]) == pytest.approx(1.3)
        assert pyo.value(m_a.fs.sieving_coefficient["Co"]) == pytest.approx(0.5)
        assert pyo.value(m_b.fs.sieving_coefficient["Li"]) == pytest.approx(1.5)
        assert pyo.value(m_b.fs.sieving_coefficient["Co"]) == pytest.approx(0.8)

        # The two technologies should not be identical
        assert pyo.value(m_a.fs.sieving_coefficient["Li"]) != pyo.value(
            m_b.fs.sieving_coefficient["Li"]
        )
        assert pyo.value(m_a.fs.sieving_coefficient["Co"]) != pyo.value(
            m_b.fs.sieving_coefficient["Co"]
        )

    # Decision variable bounds
    @pytest.mark.unit
    def test_decision_variables_bounds(self):
        m = build_diafiltration_model(
            sieving_coeffs=(1.3, 0.5), technology_name="tech_test"
        )
        decision_variables_bounds(m)
        dv = [
            m.fs.stage1.length,
            m.fs.stage2.length,
            m.fs.stage3.length,
        ]

        for v in dv:
            assert isinstance(v, pyo.Var)
            assert v.lb is not None
            assert v.ub is not None
            assert pyo.value(v.lb) <= pyo.value(v.ub)
            assert pyo.value(v.lb) == pytest.approx(0.1)
            assert pyo.value(v.ub) == pytest.approx(10000)

    # 1. Uncertain parameters and distribution specs
    @pytest.mark.unit
    def test_estimate_lognormal_params_from_data_basic(self):
        data = np.array([1.0, np.e, np.e**2], dtype=float)  # logs = [0,1,2]
        mu, sigma = estimate_lognormal_params_from_data(data)

        assert mu == pytest.approx(1.0158, rel=1e-3)
        assert sigma == pytest.approx(0.7657, rel=1e-3)

    # Test LHS sampling method
    @pytest.mark.unit
    def test_lhs_unit_shape_bounds(self):
        n, d = 10, 3
        rng = np.random.default_rng(123)

        U = lhs_unit(n_samples=n, n_dim=d, rng=rng)
        X2 = lhs_unit(n_samples=n, n_dim=d, rng=rng)

        assert isinstance(U, np.ndarray)
        assert U.shape == (n, d)

        # In [0,1]
        assert np.all(U >= 0.0)
        assert np.all(U <= 1.0)

    @pytest.mark.unit
    def test_lhs_unit_reproducible_with_seed(self):
        n, d, seed = 10, 3, 123

        X1 = lhs_unit(n_samples=n, n_dim=d, rng=np.random.default_rng(seed))
        X2 = lhs_unit(n_samples=n, n_dim=d, rng=np.random.default_rng(seed))
        assert np.allclose(X1, X2)

    @pytest.mark.unit
    def test_triangular_icdf_endpoints_and_mode(self):
        low, mode, high = 2.0, 5.0, 10.0
        Fc = (mode - low) / (high - low)

        # endpoints
        assert triangular_icdf(0.0, low, mode, high) == pytest.approx(low)
        assert triangular_icdf(1.0, low, mode, high) == pytest.approx(high)

        # exactly at CDF at mode -> returns mode
        assert triangular_icdf(Fc, low, mode, high) == pytest.approx(mode)

    @pytest.mark.unit
    def test_triangular_icdf_monotone_in_u(self):
        low, mode, high = 0.0, 1.0, 3.0
        us = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])

        xs = np.array([triangular_icdf(float(u), low, mode, high) for u in us])

        # nondecreasing and within bounds
        assert np.all(xs >= low - 1e-12)
        assert np.all(xs <= high + 1e-12)
        assert np.all(np.diff(xs) >= -1e-12)

    def test_uncertain_params_and_distributions(self, model):
        m = model
        cp = m.fs.costing

        uncertain_params = identify_uncertain_params(m)

        # Basic checks: non-empty, all are mutable Pyomo Params
        assert len(uncertain_params) > 0
        for p in uncertain_params:
            assert isinstance(p, pyo.Param)
            assert p.mutable

        def _contains_param(param_list, target):
            return any(p is target for p in param_list)

        # Check that all key costing parameters are in the uncertain set
        assert _contains_param(uncertain_params, cp.membrane_cost)
        assert _contains_param(uncertain_params, cp.factor_membrane_replacement)
        assert _contains_param(uncertain_params, cp.electricity_cost)
        assert _contains_param(uncertain_params, cp.pump_efficiency)
        assert _contains_param(uncertain_params, cp.operating_days_per_year)
        assert _contains_param(uncertain_params, cp.Lang_factor)
        assert _contains_param(uncertain_params, cp.Li_price)
        assert _contains_param(uncertain_params, cp.Co_price)
        assert _contains_param(uncertain_params, cp.income_tax_percentage)
        assert _contains_param(
            uncertain_params, cp.royalty_charge_percentage_of_revenue
        )

        # Build distribution specs using the fallback lognormal construction and a small dummy empirical sample
        specs = build_uncertainty_specs(
            m, lognormal_params=None, income_tax_samples=[21.0, 24.0, 30.5]
        )

        # Keys should match the uncertain params' names
        names_from_params = {p.getname() for p in uncertain_params}
        names_from_specs = set(specs.keys())
        assert names_from_specs == names_from_params

        # --- Triangular distributions ---
        mc_spec = specs[cp.membrane_cost.getname()]
        assert mc_spec["type"] == "triangular"
        assert mc_spec["low"] == pytest.approx(36.0)
        assert mc_spec["mode"] == pytest.approx(50.0)
        assert mc_spec["high"] == pytest.approx(450.0)

        eff_spec = specs[cp.pump_efficiency.getname()]
        assert eff_spec["type"] == "triangular"
        assert eff_spec["low"] == pytest.approx(0.1)
        assert eff_spec["mode"] == pytest.approx(0.7)
        assert eff_spec["high"] == pytest.approx(1.0)

        op_days_spec = specs[cp.operating_days_per_year.getname()]
        assert op_days_spec["type"] == "triangular"
        assert op_days_spec["low"] == pytest.approx(300.0)
        assert op_days_spec["mode"] == pytest.approx(336.0)
        assert op_days_spec["high"] == pytest.approx(365.0)

        # --- Uniform distributions ---
        repl_spec = specs[cp.factor_membrane_replacement.getname()]
        assert repl_spec["type"] == "uniform"
        assert repl_spec["low"] == pytest.approx(0.1)
        assert repl_spec["high"] == pytest.approx(0.2)

        lang_spec = specs[cp.Lang_factor.getname()]
        assert lang_spec["type"] == "uniform"
        assert lang_spec["low"] == pytest.approx(2.0)
        assert lang_spec["high"] == pytest.approx(5.93)

        tax_spec = specs[cp.income_tax_percentage.getname()]
        assert tax_spec["type"] == "discrete"
        assert "values" in tax_spec and len(tax_spec["values"]) == 3

        royalty_spec = specs[cp.royalty_charge_percentage_of_revenue.getname()]
        assert royalty_spec["type"] == "uniform"
        assert royalty_spec["low"] == pytest.approx(1.0)
        assert royalty_spec["high"] == pytest.approx(7.0)

        # --- Lognormal distributions ---
        # The fallback construction chooses (mu, sigma) so that
        #   E[X] = exp(mu + 0.5 * sigma^2) = nominal value
        def _check_lognormal_mean_matches_nominal(param):
            spec = specs[param.getname()]
            assert spec["type"] == "lognormal"
            assert "mu" in spec and "sigma" in spec
            mu = float(spec["mu"])
            sigma = float(spec["sigma"])
            assert sigma > 0.0

            mean_from_spec = math.exp(mu + 0.5 * sigma * sigma)
            nominal = float(pyo.value(param))

            # Very tight relative tolerance: they should match by construction
            assert mean_from_spec == pytest.approx(nominal, rel=1e-10)

        _check_lognormal_mean_matches_nominal(cp.electricity_cost)
        _check_lognormal_mean_matches_nominal(cp.Li_price)
        _check_lognormal_mean_matches_nominal(cp.Co_price)

    # Customized lognormal params
    @pytest.mark.unit
    def test_build_undertainty_specs_uses_custom_lognormal_params(self, model):
        m = model
        cp = m.fs.costing

        custom = {
            cp.electricity_cost.getname(): {"mu": -3.0, "sigma": 0.2},
            cp.Li_price.getname(): {"mu": 1.5, "sigma": 0.1},
            cp.Co_price.getname(): {"mu": 2.0, "sigma": 0.3},
        }

        specs = build_uncertainty_specs(
            m,
            lognormal_params=custom,
            income_tax_samples=[21.0, 24.0, 30.5],
        )

        for pname, ms in custom.items():
            assert specs[pname]["type"] == "lognormal"
            assert float(specs[pname]["mu"]) == pytest.approx(ms["mu"])
            assert float(specs[pname]["sigma"]) == pytest.approx(ms["sigma"])

    # Income-tax fallback behavior (None / empty samples)
    @pytest.mark.unit
    @pytest.mark.parametrize("samples", [None, []])
    def test_income_tax_samples_behavior(self, model, samples):
        m = model
        cp = m.fs.costing

        specs = build_uncertainty_specs(
            m,
            lognormal_params=None,
            income_tax_samples=samples,
        )

        tax_name = cp.income_tax_percentage.getname()
        tax_spec = specs[tax_name]

        if samples is None:
            # Fallback uniform path (explicit in implementation)
            assert tax_spec["type"] == "uniform"
            assert tax_spec["low"] == pytest.approx(21.0)
            assert tax_spec["high"] == pytest.approx(37.5)
            assert tax_spec["low"] < tax_spec["high"]
        else:
            # Discrete path, even if empty
            assert tax_spec["type"] == "discrete"
            assert "values" in tax_spec
            assert isinstance(tax_spec["values"], np.ndarray)
            assert tax_spec["values"].size == 0

    # Identify_uncertain_params is stable / repeatable
    @pytest.mark.unit
    def test_identify_uncertain_params_repeatable(self, model):
        p1 = identify_uncertain_params(model)
        p2 = identify_uncertain_params(model)

        # Same objects, same order (or at least same set)
        assert [id(p) for p in p1] == [id(p) for p in p2]
        assert {p.getname() for p in p1} == {p.getname() for p in p2}

    # 2. Flowsheet & costing structure (constraints from diafiltration + QGESS)
    def test_flowsheet_and_costing_structure(self, model):
        m = model

        # Diagnostics: ensure no obvious structural issues are raised
        dt = DiagnosticsToolbox(m)
        dt.report_structural_issues()
        dt.display_potential_evaluation_errors()

        # QGESS flowsheet costing block exists and is the right type
        assert hasattr(m.fs, "costing")
        assert isinstance(m.fs.costing, QGESSCosting)

        # UnitModelCostingBlocks for all key units (membranes, cascade, pumps)
        for blk_name in [
            "stage1",
            "stage2",
            "stage3",
            "cascade",
            "feed_pump",
            "diafiltrate_pump",
        ]:
            blk = getattr(m.fs, blk_name)
            assert hasattr(blk, "costing")
            assert isinstance(blk.costing, UnitModelCostingBlock)

        # QGESS costing should have cost_of_recovery and other high-level quantities
        assert hasattr(m.fs.costing, "cost_of_recovery")
        assert hasattr(m.fs.costing, "total_plant_cost")
        assert hasattr(m.fs.costing, "total_variable_OM_cost")
        assert hasattr(m.fs.costing, "total_fixed_OM_cost")

        # Diafiltration recovery-related expressions from the flowsheet
        assert hasattr(m.fs, "Li_product")
        assert hasattr(m.fs, "Co_product")
        assert hasattr(m.fs, "Li_feed")
        assert hasattr(m.fs, "Co_feed")
        assert hasattr(m.fs, "Li_recovery")
        assert hasattr(m.fs, "Co_recovery")
        assert hasattr(m.fs, "recovery_rate_per_year")
        assert hasattr(m.fs, "annual_operating_hours")

        assert isinstance(m.fs.Li_product, pyo.Expression)
        assert isinstance(m.fs.Co_product, pyo.Expression)
        assert isinstance(m.fs.Li_feed, pyo.Expression)
        assert isinstance(m.fs.Co_feed, pyo.Expression)
        assert isinstance(m.fs.Li_recovery, pyo.Expression)
        assert isinstance(m.fs.Co_recovery, pyo.Expression)
        assert isinstance(m.fs.recovery_rate_per_year, pyo.Expression)
        assert isinstance(m.fs.annual_operating_hours, pyo.Param)

        # Objective: by construction the script uses cost_of_recovery
        objs = list(m.component_objects(pyo.Objective, active=True))
        assert len(objs) == 1
        obj = objs[0]
        assert obj is m.obj
        assert obj.expr is m.fs.costing.cost_of_recovery

    @pytest.mark.unit
    def test_plot_distributions_by_technology_raises_on_empty(self):
        with pytest.raises(ValueError, match="results_by_technology is empty"):
            plot_distributions_by_technology({})

    @pytest.mark.unit
    def test_plot_distributions_by_technology_smoke(
        self, monkeypatch, tmp_path, capsys
    ):
        # Avoid writing files / GUI side effects
        monkeypatch.setattr(uq, "get_script_dir", lambda: str(tmp_path))
        monkeypatch.setattr(uq.plt, "savefig", lambda *a, **k: None)
        monkeypatch.setattr(uq.plt, "close", lambda *a, **k: None)

        results_by_technology = {
            "tech_A": {
                "samples_first_param": np.array([1.0, 2.0, 3.0]),
                "recovery_cost_samples": np.array(
                    [10.0, 11.0, 12.0, np.nan]
                ),  # valid.size>0
            },
            "tech_B": {
                "samples_first_param": np.array(
                    [9.0, 9.5, 10.0]
                ),  # not used (only first case used)
                "recovery_cost_samples": np.array([20.0, 21.0, 22.0]),  # valid.size>0
            },
        }

        uq.plot_distributions_by_technology(results_by_technology)

        out = capsys.readouterr().out
        assert "Saved plot to:" in out

    @pytest.mark.unit
    def test_analyze_sensitivity_covers_main_path(self, capsys):
        X = np.array(
            [
                [1.0, 0.0, 10.0],
                [1.0, 1.0, 11.0],
                [1.0, 2.0, 12.0],
                [1.0, 3.0, 13.0],
                [1.0, 4.0, 14.0],
                [1.0, 5.0, 15.0],
            ],
            dtype=float,
        )

        y = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)

        names = ["const_col", "x1", "x2"]

        analyze_sensitivity(X, y, names, technology_name="tech_ok")

        # Captures everything printed to the terminal
        out = capsys.readouterr().out
        assert "Sensitivity ranking (standardized linear coefficients)" in out
        assert "Sensitivity ranking (|Pearson correlation|)" in out

    @pytest.mark.unit
    def test_plot_stage_length_histograms_by_technology_raises_on_empty(self):
        with pytest.raises(ValueError, match="results_by_technology is empty"):
            uq.plot_stage_length_histograms_by_technology({})

    @pytest.mark.unit
    def test_plot_stage_length_histograms_by_technology_smoke(
        self, monkeypatch, tmp_path, capsys
    ):
        monkeypatch.setattr(uq, "get_script_dir", lambda: str(tmp_path))
        monkeypatch.setattr(uq.plt, "savefig", lambda *a, **k: None)
        monkeypatch.setattr(uq.plt, "close", lambda *a, **k: None)

        results_by_technology = {
            "tech_A": {
                "recovery_cost_samples": np.array([1.0, np.nan, 2.0, 3.0]),
                "stage1_len": np.array([10.0, 11.0, np.nan, 13.0]),
                "stage2_len": np.array([20.0, np.nan, 22.0, 23.0]),
                "stage3_len": np.array([30.0, 31.0, 32.0, np.nan]),
            },
            "tech_B": {
                "recovery_cost_samples": np.array([5.0, 6.0, np.nan]),
                "stage1_len": np.array([1.0, np.nan, 3.0]),
                "stage2_len": np.array([4.0, 5.0, 6.0]),
                "stage3_len": np.array([7.0, 8.0, np.nan]),
            },
        }

        uq.plot_stage_length_histograms_by_technology(results_by_technology)

        out = capsys.readouterr().out
        assert "Saved stage-length histogram plot to:" in out

    @pytest.mark.unit
    def test_main_monte_carlo(self, tmp_path):
        main(
            n_samples=3,
            use_lhs=False,
            run_plots=True,
            run_stage1_cost=True,
            save_plots=False,
            output_dir=str(tmp_path),
        )

    @pytest.mark.component
    def test_main_smoke_lhs(self, tmp_path):
        main(
            n_samples=3,
            use_lhs=True,
            run_plots=True,
            run_stage1_cost=True,
            save_plots=False,
            output_dir=str(tmp_path),
        )
