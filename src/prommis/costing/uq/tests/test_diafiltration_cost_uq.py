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

from prommis.costing.uq.diafiltration_cost_uq import (
    build_diafiltration_model,
    build_uncertainty_specs,
    decision_variables_bounds,
    estimate_lognormal_params_from_data,
    identify_uncertain_params,
    main,
)
from prommis.uky.costing.ree_plant_capcost import QGESSCosting


class TestDiafiltrationCostUQStructure:
    @pytest.fixture(scope="class")
    def model(self):
        """Build a fresh diafiltration UQ model for all tests."""
        m = build_diafiltration_model()
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
        m = build_diafiltration_model(sieving_coeffs=(1.3, 0.5), technology_name="tech_test")
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
        data = np.array([1.0, np.e, np.e**2], dtype=float)   # logs = [0,1,2]
        mu, sigma = estimate_lognormal_params_from_data(data)

        assert mu == pytest.approx(1.0158, rel=1e-3)
        assert sigma == pytest.approx(0.7657, rel=1e-3)
    
    # Test load income tax from csv
    
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
    def test_main(self):
        main(
            n_samples=3,
            use_lhs=False,
            run_plots=False,
            run_stage1_cost=False,
        )