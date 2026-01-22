#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Uncertainty quantification of the cost of recovery for a diafiltration process.

This script evaluates uncertainty in the cost_of_recovery for a diafiltration flowsheet under two
alternative technologies, defined by distinct combinations of lithium and cobalt sieving coefficients.

Overview of the workflow:
- Solve a deterministic diafiltration + QGESS cost flowsheet to get an optimal baseline design.
- Identify a subset of QGESS cost parameters and model them as uncertain.
- Samples uncertain parameters from specified probability distributions (triangular, uniform, discrete, 
   or lognormal) using Monte Carlo or Latin Hypercube sampling.
- Re-solve the flowsheet for each sampled parameter set and record the resulting cost_of_recovery.
- Produces summary outputs, including:
    (1) Stage1 membrane bare-erected cost (BEC) versus membrane length, with confidence intervals;
    (2) The distribution of cost_of_recovery;
    (3) The distribution of optimized membrane lengths;
    (4) Summary statistics and selected sensitivity analysis results.
"""

__author__ = "L. Deng, B. Paul, A. Fritz, A. Garciadiego, A. Ostace, and M. Zamarripa"
__version__ = "1.0.0"

import math
import os

import pyomo.environ as pyo
from pyomo.core.base.param import ScalarParam

from idaes.core import UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox, degrees_of_freedom

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCostingData,
)
from prommis.nanofiltration.diafiltration import (
    add_product_constraints,
    build_model,
    initialize_model,
    unfix_opt_variables,
)
from prommis.uky.costing.costing_dictionaries import load_default_sale_prices
from prommis.uky.costing.ree_plant_capcost import QGESSCosting

# Make sure ScalarParam has attribute units
if not hasattr(ScalarParam, "units"):
    ScalarParam.units = property(lambda self: pyo.units.get_units(self))


def get_script_dir():
    """Return the directory where this script lives."""
    return os.path.dirname(os.path.abspath(__file__))


# 0 Bounds for decision variables
def decision_variables_bounds(m):
    m.fs.stage1.length.setlb(0.1)
    m.fs.stage1.length.setub(10000)

    m.fs.stage2.length.setlb(0.1)
    m.fs.stage2.length.setub(10000)

    m.fs.stage3.length.setlb(0.1)
    m.fs.stage3.length.setub(10000)


# Sieving coefficient pairs for technology comparison:
# Case A = (1.3, 0.5), Case B = (1.5, 0.8)
def set_sieving_coefficients(m, li_sc, co_sc):
    """
    Override the Li/Co sieving coefficients on an already constructed diafiltration model.

    Note:
    The function prommis.nanofiltration.diafiltration.build_model() assigns default
    sieving coefficient values within diafiltration.py. These defaults are intentionally
    overridden here, to avoid modifying the underlying flowsheet implementation.
    """
    sc = m.fs.sieving_coefficient
    sc["Li"].fix(li_sc)
    sc["Co"].fix(co_sc)


# 1. Build the flowsheet + costing
def build_diafiltration_model(sieving_coeffs=(1.3, 0.5), technology_name=None):
    # Build the base diafiltration model
    m = build_model()

    dt = DiagnosticsToolbox(m)
    assert degrees_of_freedom(m) == 0

    # Override sieving coefficients (Li, Co) for this run
    li_sc, co_sc = sieving_coeffs
    set_sieving_coefficients(m, li_sc, co_sc)
    # Initialize the flowsheet
    initialize_model(m)
    # Add recovery constraint and unfix optimization variables.
    unfix_opt_variables(m)
    add_product_constraints(m, Li_recovery_bound=0.945, Co_recovery_bound=0.635)

    # Create dummy variables to store the UnitModelCostingBlocks
    m.fs.cascade = UnitModelBlock()  # to cost the pressure drop
    m.fs.feed_pump = UnitModelBlock()  # to cost feed pump
    m.fs.diafiltrate_pump = UnitModelBlock()  # to cost diafiltrate pump
    m.fs.costing = QGESSCosting()
    cp = m.fs.costing

    sale_prices_default = load_default_sale_prices()
    Li_default = sale_prices_default["Li"]
    Co_default = sale_prices_default["Co"]

    cp.Li_price = pyo.Param(
        mutable=True,
        initialize=pyo.value(Li_default),
        units=pyo.units.USD_2021 / pyo.units.kg,
        doc="Lithium sale price (2021_$/kg)",
    )

    cp.Co_price = pyo.Param(
        mutable=True,
        initialize=pyo.value(Co_default),
        units=pyo.units.USD_2021 / pyo.units.kg,
        doc="Cobalt sale price (2021_$/kg)",
    )

    # Same defaults / units as in diafiltration_cost_model.py
    cp.membrane_cost = pyo.Param(
        mutable=True,
        initialize=50,
        doc="Membrane unit price",
        units=pyo.units.USD_2021 / (pyo.units.m**2),
    )

    cp.factor_membrane_replacement = pyo.Param(
        mutable=True,
        initialize=0.2,
        doc="Annual membrane replacement fraction",
        units=pyo.units.year**-1,
    )

    cp.electricity_cost = pyo.Param(
        mutable=True,
        initialize=0.141,
        doc="Electricity price",
        units=pyo.units.USD_2021 / pyo.units.kWh,
    )

    cp.pump_correlation_factor = pyo.Param(
        mutable=False,
        initialize=622.59,
        doc="Pump reference cost",
        units=pyo.units.USD_1996
        / (pyo.units.kPa * (pyo.units.m**3 / pyo.units.hr)) ** 0.39,
    )

    cp.pump_exponential_factor = pyo.Param(
        mutable=False,
        initialize=0.39,
        doc="Cost-scaling exponent for pump capital cost",
        units=pyo.units.dimensionless,
    )

    cp.pump_efficiency = pyo.Param(
        mutable=True,
        initialize=0.7,
        doc="Pump efficiency",
        units=pyo.units.dimensionless,
    )

    # Costing for stages 1-3 (membrane costs)
    m.fs.stage1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage1.length,
            "membrane_width": m.w,
        },
    )
    m.fs.stage2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage2.length,
            "membrane_width": m.w,
        },
    )
    m.fs.stage3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage3.length,
            "membrane_width": m.w,
        },
    )
    # Costing of pressure drop for cascade
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop_utility,
        costing_method_arguments={
            "water_flux": m.Jw,
            "vol_flow_feed": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # cascade feed
            "vol_flow_perm": m.fs.stage3.permeate_outlet.flow_vol[
                0
            ],  # cascade permeate
        },
    )
    m.fs.feed_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,  # 14.7 psia
            "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
            * pyo.units.psi,  # this should make m.fs.feed_pump.costing.variable_operating_cost ~0
            "inlet_vol_flow": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # feed
        },
    )
    m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,  # 14.7 psia
            "outlet_pressure": m.operating_pressure,
            "inlet_vol_flow": m.fs.stage3.retentate_inlet.flow_vol[0],  # diafiltrate
        },
    )

    # Product/feed expressions and recoveries
    m.fs.Li_product = pyo.Expression(
        expr=pyo.units.convert(
            m.fs.stage3.permeate_outlet.flow_vol[0]
            * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"],
            to_units=pyo.units.kg / pyo.units.h,
        )
    )

    m.fs.Co_product = pyo.Expression(
        expr=pyo.units.convert(
            m.fs.stage1.retentate_outlet.flow_vol[0]
            * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"],
            to_units=pyo.units.kg / pyo.units.h,
        )
    )
    # Total feed = Diafiltrate + Feed
    m.fs.Li_feed = pyo.Expression(
        expr=pyo.units.convert(
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Li"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"],
            to_units=pyo.units.kg / pyo.units.h,
        )
    )

    m.fs.Co_feed = pyo.Expression(
        expr=pyo.units.convert(
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Co"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"],
            to_units=pyo.units.kg / pyo.units.h,
        )
    )

    m.fs.Li_recovery = pyo.Expression(expr=m.fs.Li_product / m.fs.Li_feed)
    m.fs.Co_recovery = pyo.Expression(expr=m.fs.Co_product / m.fs.Co_feed)

    # Operation parameters
    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    m.fs.annual_operating_hours = pyo.Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=True,
        units=pyo.units.hours / pyo.units.a,
    )

    # Define the recovery rate
    m.fs.recovery_rate_per_year = pyo.Expression(
        expr=pyo.units.convert(
            (m.fs.Li_product + m.fs.Co_product) * m.fs.annual_operating_hours,
            to_units=pyo.units.kg / pyo.units.year,
        )
    )

    # Build QGESS process costs
    m.fs.costing.build_process_costs(
        Lang_factor=2,  # includes installation, material, construction
        labor_types=[],  # labor costs already included in maintenance, admin
        fixed_OM=True,
        variable_OM=True,
        resources=[],  # List of strings of resources.
        rates=[],  # Resource consumption rate.
        prices={},  # Resource prices is not considered.
        # Assume the product is pure enough and sale price is pure product
        pure_product_output_rates={
            "Li": m.fs.Li_product,
            "Co": m.fs.Co_product,
        },
        # Override default sale prices with mutable Params.
        sale_prices={
            "Li": m.fs.costing.Li_price,
            "Co": m.fs.costing.Co_price,
        },
        recovery_rate_per_year=m.fs.recovery_rate_per_year,
        CE_index_year="2021",
        consider_taxes=True,
    )

    # Objective: minimize cost_of_recovery (for deterministic solve)
    if not any(
        isinstance(c, pyo.Objective) for c in m.component_objects(pyo.Objective)
    ):
        m.obj = pyo.Objective(expr=m.fs.costing.cost_of_recovery, sense=pyo.minimize)

    return m


# 2. Identify uncertain parameters
def identify_uncertain_params(m):
    """
    Return the set of mutable costing parameters included in the uncertainty analysis.
    """
    cp = m.fs.costing
    # membrane cost 2021_$/m2; factor membrane replacement /yr; electricity cost 2021_$/kWh;
    # pump correlation factor 1996_$/(kPa*m3/hr)**pump exponential factor; Li/Co price 2021_$/kg
    uncertain_params = [
        cp.membrane_cost,
        cp.factor_membrane_replacement,
        cp.electricity_cost,
        cp.pump_efficiency,
        cp.operating_days_per_year,
        cp.Lang_factor,
        cp.Li_price,
        cp.Co_price,
        cp.income_tax_percentage,
        cp.royalty_charge_percentage_of_revenue,
    ]

    return uncertain_params


def estimate_lognormal_params_from_data(data):
    """
    Given a 1D numpy array of strictly positive historical prices,
    estimate the underlying normal parameters (mu, sigma) such that:

        X ~ LogNormal(mu, sigma^2)

    where numpy.random.lognormal(mu, sigma) uses the same convention.

    Returns
    -------
    mu : float
    sigma : float
    """
    data = np.asarray(data, dtype=float)
    data = data[data > 0]  # drop non-positive, just in case

    if data.size == 0:
        raise ValueError("No positive data available for lognormal fit.")

    m = np.mean(data)
    v = np.var(data, ddof=1)

    # lognormal moment relationships:
    # sigma^2 = ln(1 + v / m^2)
    # mu = ln(m) - 0.5 * sigma^2
    sigma2 = math.log(1.0 + v / (m * m))
    sigma = math.sqrt(sigma2)
    mu = math.log(m) - 0.5 * sigma2

    return mu, sigma


def load_income_tax_samples_from_csv(csv_path, column_name="2021_tax"):
    """
    Load income tax percentage samples from a CSV file using a named column.
    Returns a 1D numpy array of floats with NaNs removed.
    """
    df = pd.read_csv(csv_path)

    # Be robust to extra spaces in the header (e.g., "2021_tax ")
    cols_stripped = {c: c.strip() if isinstance(c, str) else c for c in df.columns}
    df = df.rename(columns=cols_stripped)

    if column_name not in df.columns:
        raise KeyError(
            f"Column '{column_name}' not found in {csv_path}. "
            f"Available columns: {list(df.columns)}"
        )

    s = pd.to_numeric(df[column_name], errors="coerce")
    vals = s.dropna().to_numpy(dtype=float)

    if vals.size == 0:
        raise ValueError(
            f"No valid numeric values found in column '{column_name}' of: {csv_path}"
        )

    return vals


# 3. Build uncertainty specs for uncertain parameters
def build_uncertainty_specs(m, lognormal_params=None, income_tax_samples=None):
    """
    Build distribution specifications for each uncertain parameter.

    Parameters
    ----------
    lognormal_params : dict, optional
        Optional mapping from parameter names to lognormal distribution
        parameters of the form {"mu": ..., "sigma": ...}. This option is
        intended for parameters whose distributions are estimated from
        historical data. For example::

            lognormal_params = {
                "fs.costing.electricity_cost": {"mu": ..., "sigma": ...},
                "fs.costing.Li_price": {"mu": ..., "sigma": ...},
                "fs.costing.Co_price": {"mu": ..., "sigma": ...},
            }

    Notes
    -----
    If `lognormal_params` is None or does not contain an entry for a given
    parameter, a lognormal distribution is approximated using the nominal
    parameter value and an assumed coefficient of variation.
    """
    cp = m.fs.costing
    specs = {}

    # --- 1. Triangular distributions ---

    # membrane_cost: (low, mode, high) = (36, 50, 450)
    specs[cp.membrane_cost.getname()] = {
        "type": "triangular",
        "low": 36.0,
        "mode": 50.0,
        "high": 450.0,
    }

    # pump_efficiency: (0.1, 0.7, 1.0)
    specs[cp.pump_efficiency.getname()] = {
        "type": "triangular",
        "low": 0.1,
        "mode": 0.7,
        "high": 1.0,
    }

    # operating_days_per_year: (300, 336, 365)
    specs[cp.operating_days_per_year.getname()] = {
        "type": "triangular",
        "low": 300.0,
        "mode": 336.0,
        "high": 365.0,
    }

    # --- 2. Uniform distributions ---

    # factor_membrane_replacement: [0.1, 0.2]
    specs[cp.factor_membrane_replacement.getname()] = {
        "type": "uniform",
        "low": 0.1,
        "high": 0.2,
    }

    # Lang_factor: [2, 5.93]
    specs[cp.Lang_factor.getname()] = {
        "type": "uniform",
        "low": 2.0,
        "high": 5.93,
    }

    # income_tax_percentage for Pennsylvania
    if income_tax_samples is not None:
        specs[cp.income_tax_percentage.getname()] = {
            "type": "discrete",
            "values": np.asarray(income_tax_samples, dtype=float),
        }
    else:
        # Fallback to uniform assumption if no CSV samples were provided
        specs[cp.income_tax_percentage.getname()] = {
            "type": "uniform",
            "low": 21.0,
            "high": 37.5,
        }

    # royalty_charge_percentage_of_revenue: [1, 7]
    specs[cp.royalty_charge_percentage_of_revenue.getname()] = {
        "type": "uniform",
        "low": 1.0,
        "high": 7.0,
    }

    # --- 3. Lognormal distributions ---

    # Helper to build a lognormal spec either from provided (mu, sigma)
    # or from a nominal value with assumed coefficient of variation (CV).
    def _lognormal_spec_for_param(param):
        name = param.getname()
        nominal = float(pyo.value(param))

        if lognormal_params is not None and name in lognormal_params:
            mu = float(lognormal_params[name]["mu"])
            sigma = float(lognormal_params[name]["sigma"])
        else:
            # Fallback: construct a lognormal with nominal mean and an assumed CV
            # Say CV = 0.2 (20%), can adjust this if needed.
            cv = 0.2
            m = nominal
            v = (cv * m) ** 2
            sigma2 = math.log(1.0 + v / (m * m))
            sigma = math.sqrt(sigma2)
            mu = math.log(m) - 0.5 * sigma2

        return {"type": "lognormal", "mu": mu, "sigma": sigma}

    specs[cp.electricity_cost.getname()] = _lognormal_spec_for_param(
        cp.electricity_cost
    )

    specs[cp.Li_price.getname()] = _lognormal_spec_for_param(cp.Li_price)

    specs[cp.Co_price.getname()] = _lognormal_spec_for_param(cp.Co_price)

    return specs


# 4. Sampling and propagation
# 4.1 LHS sampling. equally stratify the uncertain parameters into bins
# and take one sample from each bin, and randomize them.
def lhs_unit(n_samples, n_dim, rng):
    """
    Generate a Latin Hypercube sample in [0,1]^n_dim.

    Returns an array U of shape (n_samples, n_dim) where each column
    is stratified into n_samples bins with one sample per bin.
    """
    U = np.zeros((n_samples, n_dim))
    for j in range(n_dim):
        # Stratified cut points
        cut = np.linspace(0.0, 1.0, n_samples + 1)
        # Sample one point uniformly in each bin
        u = rng.uniform(low=cut[:-1], high=cut[1:])
        # Randomly permute within the column
        rng.shuffle(u)
        U[:, j] = u
    return U


# Inverse CDF (quantile function) for a triangular distribution
def triangular_icdf(u, low, mode, high):
    """
    Inverse CDF for a triangular(low, mode, high) distribution at scalar u in [0,1].
    """
    Fc = (mode - low) / (high - low)  # CDF at the mode

    if u < Fc:  # return a single sampled value of the triangle
        return low + math.sqrt(u * (high - low) * (mode - low))
    else:
        return high - math.sqrt((1.0 - u) * (high - low) * (high - mode))


def run_LHS(
    m,
    uncertain_params,
    uncertainty_specs,
    n_samples=200,
    solver=None,
    random_seed=1,
):
    """
    Latin Hypercube Sampling propagation of uncertainty through the model.

    For each sample:
      - Draw one LHS point in [0,1]^dim for all uncertain params.
      - Map each dimension through the appropriate distribution:
          * triangular: inverse CDF
          * uniform: linear map
          * lognormal: fall back to standard lognormal MC
            (using mu, sigma in uncertainty_specs)
      - Set the Param to that value.
      - Re-solve the model (objective = cost_of_recovery).
      - Record the resulting cost_of_recovery, and optimized membrane lengths.

    Returns
    -------
    samples_first_param : (n_samples,) array of the first uncertain param.
    recovery_cost_samples : (n_samples,) array of cost_of_recovery (NaN on failure).
    param_samples : (n_samples, n_params) matrix of all parameter samples.
    membrane lengths: array of optimized membrane lengths for the three stages.
    """
    if solver is None:
        solver = pyo.SolverFactory("ipopt")

    rng = np.random.default_rng(seed=random_seed)

    n = n_samples
    n_params = len(uncertain_params)

    samples_first_param = np.zeros(n)
    recovery_cost_samples = np.full(n, np.nan)
    stage1_len = np.full(n_samples, np.nan, dtype=float)
    stage2_len = np.full(n_samples, np.nan, dtype=float)
    stage3_len = np.full(n_samples, np.nan, dtype=float)
    param_samples = np.zeros((n, n_params))

    first_param = uncertain_params[0]
    first_param_name = first_param.getname()

    # Generate LHS points in [0,1]^n_params
    U = lhs_unit(n, n_params, rng)

    # ------------------------------------------------------------------
    # Loop over samples, assign Params, solve model, record cost
    # ------------------------------------------------------------------
    for i in range(n):

        for j, p in enumerate(uncertain_params):
            name = p.getname()
            spec = uncertainty_specs[name]
            dist_type = spec["type"]
            u = float(U[i, j])

            if dist_type == "triangular":
                val = float(
                    triangular_icdf(
                        u,
                        low=spec["low"],
                        mode=spec["mode"],
                        high=spec["high"],
                    )
                )
            elif dist_type == "uniform":
                val = float(spec["low"] + u * (spec["high"] - spec["low"]))
            elif dist_type == "lognormal":
                # To avoid pulling in SciPy for inverse normal,
                # we keep lognormal as standard MC using (mu, sigma):
                val = float(
                    rng.lognormal(
                        mean=spec["mu"],
                        sigma=spec["sigma"],
                    )
                )
            elif dist_type == "discrete":
                # Use the LHS u in [0,1] as an empirical quantile index
                vals = np.asarray(spec["values"], dtype=float)
                vals = vals[~np.isnan(vals)]
                if vals.size == 0:
                    raise ValueError(
                        f"Discrete distribution for {name} has no valid values."
                    )

                vals_sorted = np.sort(vals)
                idx = int(np.floor(u * vals_sorted.size))
                if idx >= vals_sorted.size:
                    idx = vals_sorted.size - 1
                val = float(vals_sorted[idx])
            else:
                raise ValueError(f"Unknown distribution type for {name}: {dist_type}")

            p.set_value(val)
            param_samples[i, j] = val

            if name == first_param_name:
                samples_first_param[i] = val

        # Solve the model for this realization
        try:
            res = solver.solve(m, tee=False, load_solutions=False)
        except Exception as e:
            print(f"WARNING: LHS solve failed for sample {i} with exception: {e}")
            continue

        status = res.solver.status
        term = res.solver.termination_condition

        if (status == pyo.SolverStatus.ok) and (
            term
            in (
                pyo.TerminationCondition.optimal,
                pyo.TerminationCondition.locallyOptimal,
            )
        ):
            m.solutions.load_from(res)
            recovery_cost_samples[i] = pyo.value(m.fs.costing.cost_of_recovery)
            stage1_len[i] = pyo.value(m.fs.stage1.length)
            stage2_len[i] = pyo.value(m.fs.stage2.length)
            stage3_len[i] = pyo.value(m.fs.stage3.length)

        else:
            print(
                f"WARNING: non-optimal solve in LHS at sample {i}: "
                f"status={status}, term={term}"
            )

    return (
        samples_first_param,
        recovery_cost_samples,
        param_samples,
        stage1_len,
        stage2_len,
        stage3_len,
    )


def run_monte_carlo(
    m,
    uncertain_params,
    uncertainty_specs,
    n_samples=200,
    solver=None,
    random_seed=1,
):
    """
    Monte Carlo propagation of uncertainty through the model.

    For each sample:
      - Draw each uncertain param from its distribution.
      - Set the Param to that value.
      - Re-solve the model (keeping the objective as cost_of_recovery).
      - Record the resulting cost_of_recovery, and optimized membrane lengths.

    Returns:
      - samples_first_param: (n_samples,) array of the first uncertain param.
      - recovery_cost_samples: (n_samples,) array of cost_of_recovery.
      - membrane lengths: array of optimized membrane lengths for the three stages.
    """
    if solver is None:
        solver = pyo.SolverFactory("ipopt")

    rng = np.random.default_rng(seed=random_seed)

    n = n_samples
    n_params = len(uncertain_params)

    samples_first_param = np.zeros(n)
    recovery_cost_samples = np.full(n, np.nan)  # start everything as NaN
    stage1_len = np.full(n_samples, np.nan, dtype=float)
    stage2_len = np.full(n_samples, np.nan, dtype=float)
    stage3_len = np.full(n_samples, np.nan, dtype=float)

    param_samples = np.zeros((n, n_params))  # store all draws

    first_param = uncertain_params[0]
    first_param_name = first_param.getname()

    for i in range(n):
        # --- sample all uncertain parameters ---
        for j, p in enumerate(uncertain_params):
            name = p.getname()
            spec = uncertainty_specs[name]
            dist_type = spec["type"]

            if dist_type == "triangular":
                val = float(
                    rng.triangular(
                        left=spec["low"],
                        mode=spec["mode"],
                        right=spec["high"],
                    )
                )
            elif dist_type == "uniform":
                val = float(
                    rng.uniform(
                        low=spec["low"],
                        high=spec["high"],
                    )
                )
            elif dist_type == "lognormal":
                val = float(
                    rng.lognormal(
                        mean=spec["mu"],
                        sigma=spec["sigma"],
                    )
                )
            elif dist_type == "discrete":
                vals = np.asarray(spec["values"], dtype=float)
                vals = vals[~np.isnan(vals)]
                if vals.size == 0:
                    raise ValueError(
                        f"Discrete distribution for {name} has no valid values."
                    )
                val = float(rng.choice(vals))
            else:
                raise ValueError(f"Unknown distribution type for {name}: {dist_type}")

            p.set_value(val)

            # store full parameter sample
            param_samples[i, j] = val

            # track the first uncertain parameter
            if name == first_param_name:
                samples_first_param[i] = val

        # --- solve the model for this realization ---
        try:
            # IMPORTANT: don't auto-load a bad solution into the model
            res = solver.solve(m, tee=False, load_solutions=False)
        except Exception as e:
            # Something went really wrong at the solver level; leave this sample as NaN
            print(f"WARNING: solve failed for sample {i} with exception: {e}")
            continue

        status = res.solver.status
        term = res.solver.termination_condition

        if (status == pyo.SolverStatus.ok) and (
            term
            in (
                pyo.TerminationCondition.optimal,
                pyo.TerminationCondition.locallyOptimal,
            )
        ):
            # Now it is safe to load the solution
            m.solutions.load_from(res)
            recovery_cost_samples[i] = pyo.value(m.fs.costing.cost_of_recovery)
            stage1_len[i] = pyo.value(m.fs.stage1.length)
            stage2_len[i] = pyo.value(m.fs.stage2.length)
            stage3_len[i] = pyo.value(m.fs.stage3.length)

        else:
            # leave as NaN
            print(
                f"WARNING: non-optimal solve at sample {i}: "
                f"status={status}, term={term}"
            )

    return (
        samples_first_param,
        recovery_cost_samples,
        param_samples,
        stage1_len,
        stage2_len,
        stage3_len,
    )


# 5. Sensitivity analysis
def analyze_sensitivity(
    param_samples, recovery_cost_samples, param_names, technology_name=None
):
    """
    Simple global sensitivity analysis using:
      - standardized linear regression coefficients
      - Pearson correlations

    param_samples: (n_samples, n_params)
    recovery_cost_samples: (n_samples,)
    param_names: list of strings (length n_params)
    """
    if technology_name is None:
        technology_name = "technology"

    print(
        f"\n\n#############################\nSensitivity analysis: {technology_name}\n#############################"
    )
    # Keep only successful solves
    valid = ~np.isnan(recovery_cost_samples)
    X = param_samples[valid, :]
    y = recovery_cost_samples[valid]

    if X.shape[0] < 5:
        print("Not enough valid samples for sensitivity analysis.")
        return

    # Standardize columns of X
    means = X.mean(axis=0)
    stds = X.std(axis=0, ddof=1)
    # avoid division by zero
    stds[stds == 0.0] = 1.0
    X_std = (X - means) / stds

    # Center y
    y_centered = y - y.mean()

    # Solve least squares: X_std * beta = y_centered
    beta, *_ = np.linalg.lstsq(X_std, y_centered, rcond=None)
    abs_beta = np.abs(beta)

    # Pearson correlations
    corrs = []
    for j in range(X.shape[1]):
        xj = X[:, j]
        if xj.std(ddof=1) == 0:
            corrs.append(0.0)
        else:
            rho = np.corrcoef(xj, y)[0, 1]
            corrs.append(rho)
    corrs = np.array(corrs)
    abs_corrs = np.abs(corrs)

    # Ranking by absolute standardized effect
    idx_sorted_beta = np.argsort(-abs_beta)
    idx_sorted_corr = np.argsort(-abs_corrs)

    print("\n=== Sensitivity ranking (standardized linear coefficients) ===")
    for idx in idx_sorted_beta:
        print(f"{param_names[idx]:40s}  |beta| = {abs_beta[idx]:.3g}")

    print("\n=== Sensitivity ranking (|Pearson correlation|) ===")
    for idx in idx_sorted_corr:
        print(f"{param_names[idx]:40s}  |rho|  = {abs_corrs[idx]:.3g}")


# 6. Plotting inputs and outputs
def plot_distributions_by_technology(
    results_by_technology, save_plot=True, output_dir=None
):
    """
    - Histogram of the fist uncertain parameter
    - Histogram of cost_of_recovery (ignoring NaNs)
    - Compare different sieving coefficient technology
    """
    if not results_by_technology:
        raise ValueError("results_by_technology is empty")

    technology_names = list(results_by_technology.keys())

    # Use the first case to show the first-uncertain-parameter distribution
    first_case = results_by_technology[technology_names[0]]
    membrane_cost_samples = np.asarray(first_case["samples_first_param"], dtype=float)

    n_technology = len(technology_names)
    fig, axes = plt.subplots(
        1, 1 + n_technology, figsize=(6 * (1 + n_technology), 4), squeeze=False
    )
    axes = axes[0]

    # (0) membrane_cost samples
    axes[0].hist(membrane_cost_samples, bins=30, density=True, alpha=0.5)
    axes[0].set_title("First uncertain parameter: membrane_cost")
    axes[0].set_xlabel("Value")
    axes[0].set_ylabel("Density")

    # (1) cost_of_recovery histogram per case
    for k, technology_name in enumerate(technology_names, start=1):
        rec = np.asarray(
            results_by_technology[technology_name]["recovery_cost_samples"], dtype=float
        )
        valid = rec[~np.isnan(rec)]

        axes[k].hist(valid, bins=30, density=True, alpha=0.7)
        axes[k].set_title(f"cost_of_recovery\n{technology_name}")
        axes[k].set_xlabel("cost_of_recovery (USD/kg)")
        axes[k].set_ylabel("Density")

        # Mark key quantiles (valid only)
        if valid.size > 0:
            q5, q50, q95 = np.percentile(valid, [5, 50, 95])
            axes[k].axvline(q5, linestyle="--", linewidth=2, color="tab:blue")
            axes[k].axvline(q50, linestyle="-.", linewidth=2, color="tab:orange")
            axes[k].axvline(q95, linestyle="--", linewidth=2, color="tab:red")
            axes[k].legend(
                [f"5%: {q5:.4g}", f"50%: {q50:.4g}", f"95%: {q95:.4g}"],
                fontsize=8,
            )

    fig.tight_layout()

    script_dir = output_dir if output_dir is not None else get_script_dir()
    suffix = f"_{technology_name}" if technology_name else ""

    out_file = os.path.join(
        script_dir, f"uncertainty_distributions_by_sieving_case ({suffix}).png"
    )

    if save_plot:
        plt.savefig(out_file, dpi=200)
        print(f"Saved plot to: {out_file}", flush=True)
    plt.close(fig)


# 7. Single membrane (stage1) with Lang factor considered confidence interval analysis
def analyze_stage1_membrane_cost(
    m,
    technology_name=None,
    n_lengths=50,  # number of length points between 0.1 and 10000
    n_samples_per_length=50,  # number of uncertainty samples per length
    random_seed=1,
    save_plot=True,
    output_dir=None,
):
    """Explore stage 1 membrane capital cost vs membrane length.

    For each membrane length L in [0.1, 10000] m:
      - Compute membrane area = L * m.w.
      - Compute the *nominal* capital cost for that membrane:
            cost_nominal = area * membrane_cost_nominal * Lang_factor_nominal
        where membrane_cost_nominal = 50 USD/m^2
              Lang_factor_nominal   = 4 (user-specified nominal value)
      - Sample membrane_cost and Lang_factor from their uncertainty ranges:
            membrane_cost ~ Triangular(36, 50, 450)
            Lang_factor   ~ Uniform(2, 5.93)
        and compute sampled capital costs:
            cost_sample = area * membrane_cost * Lang_factor
      - Store the nominal and sampled costs.

    A scatter plot is generated with:
        x-axis: membrane length (m)
        y-axis: capital cost = membrane_cost * Lang_factor * area (USD)
        blue dots: sampled (uncertain) costs
        black line: nominal costs

    The figure is saved in the same directory as this script as
        "stage1_membrane_cost_vs_length.png".
    """
    rng = np.random.default_rng(random_seed)

    # Membrane length grid
    lengths = np.linspace(0.1, 10000.0, n_lengths)

    # Width of the membrane from the model
    width = float(pyo.value(m.w))

    # Nominal values specified by the user
    membrane_cost_nominal = 50.0  # USD_2021 / m^2
    Lang_factor_nominal = 4.0  # dimensionless

    # Uncertainty ranges
    membrane_cost_low, membrane_cost_mode, membrane_cost_high = 36.0, 50.0, 450.0
    Lang_factor_low, Lang_factor_high = 2.0, 5.93

    nominal_costs = []
    scatter_lengths = []
    scatter_costs = []

    # For CI curves
    ci_lower = []
    ci_upper = []

    for L in lengths:
        area = L * width

        # Nominal capital cost for this length
        cost_nominal = area * membrane_cost_nominal * Lang_factor_nominal
        nominal_costs.append(cost_nominal)

        # Sample uncertainty for this membrane length
        memmbrane_samples = rng.triangular(
            membrane_cost_low,
            membrane_cost_mode,
            membrane_cost_high,
            size=n_samples_per_length,
        )
        lang_factor_samples = rng.uniform(
            Lang_factor_low, Lang_factor_high, size=n_samples_per_length
        )

        cost_samples = area * memmbrane_samples * lang_factor_samples

        scatter_lengths.extend([L] * n_samples_per_length)
        scatter_costs.extend(cost_samples)

        # 95% confidence interval (2.5%-97.5%)
        ci_l = np.percentile(cost_samples, 2.5)
        ci_u = np.percentile(cost_samples, 97.5)
        ci_lower.append(ci_l)
        ci_upper.append(ci_u)

    scatter_lengths = np.array(scatter_lengths)
    scatter_costs = np.array(scatter_costs)
    nominal_costs = np.array(nominal_costs)
    ci_lower = np.array(ci_lower)
    ci_upper = np.array(ci_upper)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        scatter_lengths,
        scatter_costs,
        s=10,
        alpha=0.3,
        color="tab:blue",
        label="Uncertain membrane BEC",
    )

    ax.plot(
        lengths,
        nominal_costs,
        "k-",
        linewidth=2,
        label="Nominal membrane BEC",
    )

    # 95% CI bounds
    ax.plot(lengths, ci_lower, "r--", linewidth=1.5, label="95% CI (lower)")
    ax.plot(lengths, ci_upper, "r--", linewidth=1.5, label="95% CI (upper)")

    ax.set_xlabel("Stage 1 membrane length (m)")
    ax.set_ylabel("Membrane BEC (USD)")
    ax.set_title(
        "Stage 1 membrane BEC vs membrane length\nwith 95% confidence interval"
    )

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend()

    script_dir = output_dir if output_dir is not None else get_script_dir()
    suffix = f"_{technology_name}" if technology_name else ""

    out_file = os.path.join(script_dir, f"stage1_membrane_BEC_vs_length({suffix}).png")
    if save_plot:
        plt.savefig(out_file, dpi=200)
        print(f"Saved stage1 membrane BEC plot to: {out_file}", flush=True)
    plt.close(fig)


# Plot membrane length histogram
def plot_stage_length_histograms_by_technology(
    results_by_technology, save_plot=True, output_dir=None
):
    if not results_by_technology:
        raise ValueError("results_by_technology is empty")

    techs = list(results_by_technology.keys())
    ntech = len(techs)
    fig, axes = plt.subplots(ntech, 3, figsize=(6 * 3, 4 * ntech), squeeze=False)

    for r, tech in enumerate(techs):
        d = results_by_technology[tech]

        # Use only the valid samples
        valid = ~np.isnan(np.asarray(d["recovery_cost_samples"], dtype=float))

        L1 = np.asarray(d["stage1_len"], dtype=float)[valid]
        L2 = np.asarray(d["stage2_len"], dtype=float)[valid]
        L3 = np.asarray(d["stage3_len"], dtype=float)[valid]

        axes[r, 0].hist(L1[~np.isnan(L1)], bins=30, density=True, alpha=0.7)
        axes[r, 0].set_title(f"Stage 1 length\n{tech}")
        axes[r, 0].set_xlabel("Length (m)")
        axes[r, 0].set_ylabel("Density")

        axes[r, 1].hist(L2[~np.isnan(L2)], bins=30, density=True, alpha=0.7)
        axes[r, 1].set_title(f"Stage 2 length\n{tech}")
        axes[r, 1].set_xlabel("Length (m)")
        axes[r, 1].set_ylabel("Density")

        axes[r, 2].hist(L3[~np.isnan(L3)], bins=30, density=True, alpha=0.7)
        axes[r, 2].set_title(f"Stage 3 length\n{tech}")
        axes[r, 2].set_xlabel("Length (m)")
        axes[r, 2].set_ylabel("Density")

    fig.tight_layout()

    script_dir = output_dir if output_dir is not None else get_script_dir()
    out_file = os.path.join(script_dir, "stage_lengths_histogram_by_sieving_case.png")
    if save_plot:
        plt.savefig(out_file, dpi=200)
        print(f"Saved stage-length histogram plot to: {out_file}", flush=True)
    plt.close(fig)


# 8. Main driver
def main(
    n_samples=200,
    use_lhs=False,
    random_seed=1,
    run_plots=True,
    run_stage1_cost=True,
    solver_name="ipopt",
    max_iter=5000,
    tol=1e-6,
    acceptable_tol=1e-5,
    save_plots=True,
    output_dir=None,
):
    technologies = {
        "Li_sc=1.3, Co_sc=0.5": (1.3, 0.5),
        "Li_sc=1.5, Co_sc=0.8": (1.5, 0.8),
    }

    script_dir = output_dir if output_dir is not None else get_script_dir()

    results_by_technology = {}

    N_SAMPLES = n_samples  # Sample size
    USE_LHS = use_lhs  # Sampling method

    # Set up solver
    local_solver = pyo.SolverFactory(solver_name)
    local_solver.options["max_iter"] = max_iter
    local_solver.options["tol"] = tol
    local_solver.options["acceptable_tol"] = acceptable_tol

    for technology_name, sieving_coeffs in technologies.items():
        print(f"\n===Running technology {technology_name} ===")

        # Build model for each technology
        m = build_diafiltration_model(
            sieving_coeffs=sieving_coeffs, technology_name=technology_name
        )

        # Bound decision variables
        decision_variables_bounds(m)

        # Deterministic baseline optimization
        det_results = local_solver.solve(m)
        pyo.assert_optimal_termination(det_results)
        det_obj = pyo.value(m.fs.costing.cost_of_recovery)
        print(
            f"Deterministic optimal cost_of_recovery of membrane with sieving at ({technology_name}) is: {det_obj:.6f}"
        )

        # Identify uncertain parameters
        uncertain_params = identify_uncertain_params(m)

        # Lognormal params from historical data
        # Lithium price USD/kg
        script_dir = get_script_dir()
        li_file = os.path.join(
            script_dir, "lithium_price_USA_dailymetalprice_yr2021.csv"
        )
        li_df = pd.read_csv(li_file)
        li_data = li_df["2021_$/kg"].values
        mu_li, sigma_li = estimate_lognormal_params_from_data(li_data)

        # Cobalt price USD/kg
        co_file = os.path.join(
            script_dir, "cobalt_price_USA_tradingeconomics_yr2021.csv"
        )
        co_df = pd.read_csv(co_file)
        co_data = co_df["2021_$/kg"].values
        mu_co, sigma_co = estimate_lognormal_params_from_data(co_data)

        # Electricity price of Pennsylvania
        elec_file = os.path.join(
            script_dir, "iea_1990_2025_industry_elec_monthly_price_PA.csv"
        )
        elec_df = pd.read_csv(elec_file)
        elec_data = elec_df["PA_2021_$/kWh"].values
        mu_elec, sigma_elec = estimate_lognormal_params_from_data(elec_data)

        # Map to the actual Param names used in the model
        cp = m.fs.costing
        lognormal_params = {
            cp.electricity_cost.getname(): {"mu": mu_elec, "sigma": sigma_elec},
            cp.Li_price.getname(): {"mu": mu_li, "sigma": sigma_li},
            cp.Co_price.getname(): {"mu": mu_co, "sigma": sigma_co},
        }

        print("Lognormal parameter estimates:")
        print(f"  electricity_cost: mu={mu_elec:.4g}, sigma={sigma_elec:.4g}")
        print(f"  Li_price:         mu={mu_li:.4g},   sigma={sigma_li:.4g}")
        print(f"  Co_price:         mu={mu_co:.4g},   sigma={sigma_co:.4g}")

        # Income tax empirical samples
        income_tax_file = os.path.join(script_dir, "PA_income_tax.csv")
        income_tax_samples = load_income_tax_samples_from_csv(
            income_tax_file, column_name="2021_tax"
        )
        print(
            f"Loaded {income_tax_samples.size} income tax samples from: {income_tax_file}"
        )

        # Build mixed distribution specs
        uncertainty_specs = build_uncertainty_specs(
            m,
            lognormal_params=lognormal_params,
            income_tax_samples=income_tax_samples,
        )

        # Run sampling
        if use_lhs:
            print(f"Running LHS with {N_SAMPLES} samples...")
            (
                samples_first_param,
                recovery_cost_samples,
                param_samples,
                stage1_len,
                stage2_len,
                stage3_len,
            ) = run_LHS(
                m,
                uncertain_params,
                uncertainty_specs,
                n_samples=N_SAMPLES,
                solver=local_solver,
                random_seed=random_seed,
            )
        else:
            print(f"Running Monte Carlo with {N_SAMPLES} samples...")
            (
                samples_first_param,
                recovery_cost_samples,
                param_samples,
                stage1_len,
                stage2_len,
                stage3_len,
            ) = run_monte_carlo(
                m,
                uncertain_params,
                uncertainty_specs,
                n_samples=N_SAMPLES,
                solver=local_solver,
                random_seed=random_seed,
            )

        # Build parameter names
        param_names = [p.getname() for p in uncertain_params]
        # Sensitivity analysis
        analyze_sensitivity(
            param_samples,
            recovery_cost_samples,
            param_names,
            technology_name=technology_name,
        )

        # Summary statistics
        valid_cost = recovery_cost_samples[~np.isnan(recovery_cost_samples)]
        print(f"\n=== Summary statistics: {technology_name} ===")
        print(f"Valid samples: {valid_cost.size}/{N_SAMPLES}")
        if valid_cost.size > 0:
            mean_cost = np.mean(valid_cost)
            std_cost = np.std(valid_cost, ddof=1)
            q5, q50, q95 = np.percentile(valid_cost, [5, 50, 95])
            print(f"Mean cost_of_recovery: {mean_cost:.6f}")
            print(f"Std  cost_of_recovery: {std_cost:.6f}")
            print(f"5% / 50% / 95% quantiles: {q5:.6f}, {q50:.6f}, {q95:.6f}")

        results_by_technology[technology_name] = {
            "samples_first_param": samples_first_param,
            "recovery_cost_samples": recovery_cost_samples,
            "param_samples": param_samples,
            "stage1_len": stage1_len,
            "stage2_len": stage2_len,
            "stage3_len": stage3_len,
            "m": m,
        }

        # Stage 1 membrane BEC:
        if run_stage1_cost:
            analyze_stage1_membrane_cost(
                m,
                technology_name=technology_name.replace("=", "")
                .replace(",", "_")
                .replace(" ", ""),
                n_lengths=50,
                n_samples_per_length=50,
                random_seed=random_seed,
                save_plot=save_plots,
                output_dir=script_dir,
            )

    # Cross-case comparison plot
    if run_plots:
        plot_distributions_by_technology(
            results_by_technology,
            save_plot=save_plots,
            output_dir=script_dir,
        )
        plot_stage_length_histograms_by_technology(
            results_by_technology,
            save_plot=save_plots,
            output_dir=script_dir,
        )


if __name__ == "__main__":
    main()
