#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This test is a modified version from
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/tests/test_ion_exchange_clark.py

# WaterTAP License Agreement

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import os
import pytest
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.scaling.scaling_base import ScalerBase
from idaes.core.util.testing import initialization_tester
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from prommis.ion_exchange.ion_exchange_multicomponent import (
    IonExchangeMultiComp,
)

"""
modified by: Soraya Rawlings
"""


@pytest.fixture(scope="module")
def m(flow_factor=1):

    path = os.path.dirname(os.path.realpath(__file__))
    resin_file = os.path.join(path, "..", "data", "resin_data.json")
    resin = "S950"
    target_component = "La"
    regenerant = "single_use"
    list_solvent = ["H2O"]
    list_reactive_ions = ["La", "Dy", "Ho", "Er", "Yb", "Sm"]

    # Add sets for solvent and ion species
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.set_solvent = pyo.Set(initialize=list_solvent)
    m.fs.set_reactive_ions = pyo.Set(initialize=list_reactive_ions)
    m.fs.set_all = pyo.Set(initialize=list_solvent + list_reactive_ions)

    ion_props = {
        "solute_list": [],
        "diffusivity_data": {},
        "molar_volume_data": {},
        "mw_data": {},  # in kg/mol
        "charge": {},
    }

    ion_props["solute_list"] = list_reactive_ions
    ion_props["diffusivity_data"] = {
        ("Liq", "La"): 6.19e-10,
        ("Liq", "Dy"): 5.82e-10,
        ("Liq", "Sm"): 6.08e-10,
    }
    ion_props["molar_volume_data"] = {
        ("Liq", "Ho"): 1.8753e-05,
        ("Liq", "Er"): 1.8449e-05,
        ("Liq", "Yb"): 2.484e-05,
    }
    ion_props["mw_data"] = {
        "H2O": 0.018,
        "La": 0.138905,
        "Dy": 0.1625,
        "Ho": 0.16493,
        "Er": 0.16726,
        "Yb": 0.17305,
        "Sm": 0.15036,
    }
    ion_props["charge"] = {"La": +3, "Dy": +3, "Ho": +3, "Er": +3, "Yb": +3, "Sm": +3}
    ion_props["diffus_calculation"] = "HaydukLaudie"
    m.fs.properties = MCASParameterBlock(**ion_props)

    ix_config = {
        "property_package": m.fs.properties,
        "regenerant": regenerant,
        "target_component": target_component,
        "reactive_ions": list_reactive_ions,
        "number_traps": 25,
        "c_trap_min": 1e-3,  # initial concentration in kg/m3
        "resin_data_path": resin_file,
        "resin": resin,
    }

    m.fs.unit = ix = IonExchangeMultiComp(**ix_config)

    @m.fs.unit.Expression(
        m.fs.set_reactive_ions, doc="Percentage of recovery for all the REEs"
    )
    def recovery_comp(b, s):
        conc_in = b.process_flow.properties_in[0].conc_mass_phase_comp["Liq", s]
        conc_out = b.process_flow.properties_out[0].conc_mass_phase_comp["Liq", s]
        return ((conc_in - conc_out) / conc_in) * 100

    # Add bounds
    ix.bed_diameter.setlb(0.01)
    ix.bed_diameter.setub(100)
    ix.bed_depth.setlb(0.01)
    ix.bed_depth.setub(100)
    ix.service_flow_rate.setlb(1e-10)
    ix.process_flow.properties_in[0.0].visc_k_phase["Liq"].setlb(1e-16)
    ix.process_flow.properties_in[0].flow_vol_phase["Liq"].setlb(1e-16)
    ix.N_Sc.setlb(1e-16)
    for c in m.fs.set_reactive_ions:
        ix.process_flow.properties_in[0.0].diffus_phase_comp["Liq", c].setlb(1e-16)
        ix.bv_50[c].setlb(1e-3)
        ix.loading_rate.setlb(1e-16)
        ix.process_flow.properties_in[0].flow_mass_phase_comp["Liq", c].setlb(1e-16)
        ix.process_flow.properties_out[0].flow_mass_phase_comp["Liq", c].setlb(1e-16)
        for i in range(1, ix.num_traps + 1):
            ix.c_traps[c, i].setlb(1e-16)

    # Set operating conditions
    flow_mol = {  # in mol/s
        "H2O": 0.0009277777777777778,
        "La": 1.2022605377776176e-10,
        "Dy": 3.288615384615384e-11,
        "Ho": 5.2652640514157516e-12,
        "Er": 1.1981346406791819e-11,
        "Yb": 7.720312048540882e-12,
        "Sm": 4.553737696195796e-11,
    }

    ix.process_flow.properties_in[0].pressure.fix(101325)
    ix.process_flow.properties_in[0].temperature.fix(298.15)
    for i in m.fs.set_all:
        ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", i].fix(
            flow_mol[i] * flow_factor
        )

    ix.resin_diam.fix(0.00075)
    ix.resin_density.fix(1126.61)
    ix.bed_porosity.fix(0.80)
    ix.bed_depth.fix(1)
    ix.bed_diameter.fix(0.015)
    ix.number_columns.fix(1)
    ix.number_columns_redundant.fix(1)

    # Equilibrium parameters
    parmest_data = {
        "mass_transfer_coeff": {
            "La": 83.22416280992222,
            "Dy": 59.484016710115064,
            "Er": 1.7286360467604134,
            "Ho": 75.60499429888984,
            "Sm": 77.67086243578649,
            "Yb": 41.26897106783112,
        },
        "freundlich_n": {
            "La": 1.0000026107532496,
            "Dy": 1.000003681926298,
            "Er": 1.0001257362841591,
            "Ho": 1.0000030377701503,
            "Sm": 1.0000027974184122,
            "Yb": 1.0000052673068462,
        },
        "bv_50": {
            "La": 56.32137778092139,
            "Dy": 64.1931100064633,
            "Er": 79.5742013052232,
            "Ho": 54.08209908870865,
            "Sm": 67.48274394895874,
            "Yb": 53.16220908663995,
        },
    }
    for c in m.fs.set_reactive_ions:
        ix.bv[c].setub(400)
        ix.bv[c].set_value(10)

        ix.freundlich_n[c].fix(parmest_data["freundlich_n"][c])
        ix.mass_transfer_coeff[c].fix(parmest_data["mass_transfer_coeff"][c])
        ix.bv_50[c].fix(parmest_data["bv_50"][c])

        if c == target_component:
            ix.c_norm[c].fix(0.99)
        else:
            ix.c_norm[c].fix(0.99)

    # Set initial values
    ix.loading_rate.set_value(9.43138198040347e-05)
    ix.ebct.set_value(10581.711782300948)  # in seconds

    # Add scaling factors
    scale_from_value = [
        "bed_volume_total",
        "bed_diameter",
        "column_height",
        "service_flow_rate",
        "ebct",
        "N_Re",
        "N_Sc",
        "N_Sh",
        "N_Pe_bed",
        "N_Pe_particle",
        "bv",
    ]

    def get_comp_list(blk, comp=pyo.Var):
        cs = []
        split_name = blk.name + "."
        skip_list = ["ref", "process_flow", "regeneration"]
        for c in blk.component_objects(comp):
            if any(s in c.name for s in skip_list):
                continue
            cs.append(c.name.split(split_name)[1])
        return cs

    # Scaling
    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

    sb = ScalerBase()

    # Apply scaling to variables
    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "temperature" in var.name:
            sb.set_variable_scaling_factor(var, 1e-1, overwrite=True)
        if "pressure" in var.name:
            sb.set_variable_scaling_factor(var, 1e-5)
        if "flow_vol" in var.name:
            sb.set_variable_scaling_factor(var, 1e8)
        if "flow_mol_phase_comp" in var.name:
            if "H2O" in var.name:
                sb.set_variable_scaling_factor(var, 1e4)
            else:
                sb.set_variable_scaling_factor(var, 1e8)

    for v in get_comp_list(ix):
        ixv = getattr(ix, v)
        if ixv.is_indexed():
            idx = [*ixv.index_set()]
            for i in idx:
                if ixv[i].is_fixed():
                    if ixv[i]() == 0:
                        continue
                    sb.set_variable_scaling_factor(ixv[i], 1 / ixv[i]())
        else:
            if ixv.is_fixed():
                if ixv() == 0:
                    continue
                sb.set_variable_scaling_factor(ixv, 1 / ixv())

    for v in scale_from_value:
        ixv = getattr(ix, v)
        if ixv.is_indexed():
            idx = [*ixv.index_set()]
            for i in idx:
                if ixv[i]() == 0:
                    continue
                sb.set_variable_scaling_factor(ixv[i], 1 / ixv[i]())
        else:
            if ixv() == 0:
                continue
            sb.set_variable_scaling_factor(ixv, 1 / ixv())

    return m


def build_clark_with_costing(m, flow_factor=1):

    ix = m.fs.unit
    flow_out = ix.process_flow.properties_out[0].flow_vol_phase["Liq"]

    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyo.units.USD_2021
    m.fs.costing.base_period = pyo.units.year
    m.fs.costing.utilization_factor.fix(1)

    # Add costs related to IX unit
    ix.costing = ix_cost = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Touch properties
    m.fs.unit.costing.capital_cost
    m.fs.unit.costing.capital_cost_resin

    # Calculate costs of entireprocess
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(flow_out)
    m.fs.costing.add_annual_water_production(flow_out)
    m.fs.costing.add_specific_energy_consumption(flow_out)

    # Add costs for REEs. References are: [a]
    # https://www.metal.com/Rare-Earth-Metals/ and [b]
    # https://www.metal.com/price/Rare%20Earth/Rare-Earth-Oxides.
    market_prices = {
        # From ref[a]
        "La": 2.642,
        "Dy": 247.07,
        "Ho": 63.610,
        # From ref[b]
        "Er": 39.641,
        "Yb": 12.291,
        "Sm": 8.973,
    }
    m.fs.market_price = pyo.Param(
        m.fs.set_reactive_ions,
        initialize=market_prices,
        units=pyo.units.USD_2021 / pyo.units.kg,
        doc="Market price for REEs",
    )

    m.fs.operational_daily_hours = pyo.Param(
        initialize=8, units=pyo.units.hours / pyo.units.day, doc="IX operational hours"
    )
    m.fs.operational_yearly_days = pyo.Param(
        initialize=365, units=pyo.units.day / pyo.units.year, doc="IX operational hours"
    )
    m.fs.ix_lifetime = pyo.Param(
        initialize=15, units=pyo.units.years, doc="IX lifetime"
    )

    @m.fs.unit.Expression()
    def expected_annual_profit_ree(b):
        return (
            sum(
                m.fs.market_price[s]
                * pyo.units.convert(
                    (
                        b.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", s]
                        - b.process_flow.properties_out[0.0].flow_mass_phase_comp[
                            "Liq", s
                        ]
                    ),
                    to_units=pyo.units.kg / pyo.units.hour,
                )
                for s in m.fs.set_reactive_ions
            )
            * m.fs.operational_daily_hours
            * m.fs.operational_yearly_days
        )

    return m


@pytest.mark.unit
def test_structural_issues(m):
    dt = DiagnosticsToolbox(m)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings()


@pytest.mark.component
def test_resin_specific_data(m):

    ix = m.fs.unit
    resin = "S950"

    # Test that data coming from .json file is correct
    test_json_data = {
        "S950": {
            "bed_expansion": {"param_a": -0.0004, "param_b": 0.0587, "param_c": 0},
            "pressure_drop": {"param_a": -0.2332, "param_b": 0.1714, "param_c": 0},
        }
    }

    assert test_json_data[resin]["bed_expansion"]["param_a"] == pyo.value(
        ix.bed_expansion_frac_param_A
    )
    assert test_json_data[resin]["bed_expansion"]["param_b"] == pyo.value(
        ix.bed_expansion_frac_param_B
    )
    assert test_json_data[resin]["bed_expansion"]["param_c"] == pyo.value(
        ix.bed_expansion_frac_param_C
    )

    assert test_json_data[resin]["pressure_drop"]["param_a"] == pyo.value(
        ix.pressure_drop_param_A
    )
    assert test_json_data[resin]["pressure_drop"]["param_b"] == pyo.value(
        ix.pressure_drop_param_B
    )
    assert test_json_data[resin]["pressure_drop"]["param_c"] == pyo.value(
        ix.pressure_drop_param_C
    )

    # Scale model
    init_scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = init_scaling.create_using(m, rename=False)
    initialization_tester(scaled_model)

    # Solve scaled model with zero degrees of freedom
    solver = get_solver()
    scaled_results = solver.solve(scaled_model)
    pyo.assert_optimal_termination(scaled_results)

    # Propagate the solution back to the original model
    results = init_scaling.propagate_solution(scaled_model, m)

    # Test expected results for ix variables. NOTE: bed_expansion_frac
    # and pressure_drop are the ones calculated with resin-specific
    # data from .json file
    ix_vars_results = {
        "loading_rate": 9.450285573574311e-05,
        "bed_depth": 1,
        "bed_diameter": 0.015,
        "bed_volume_total": 0.00017671458676442585,
        "target_breakthrough_time": 1693392.3126900392,
        "service_flow_rate": 0.3402102806486752,
        "column_height": 2.2931,
    }

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit, v)
        print(f"mv={mv}")
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)


@pytest.mark.component
def test_ion_exchange_clark_optimization(m):

    target_component = "La"

    flow_factor = 1
    m = build_clark_with_costing(m, flow_factor=flow_factor)

    check_dof(m, fail_flag=True)

    # Scale model
    init_scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = init_scaling.create_using(m, rename=False)
    initialization_tester(scaled_model)

    # Solve scaled model
    solver = get_solver()
    results = solver.solve(scaled_model)
    pyo.assert_optimal_termination(results)

    # Propagate the solution back to the original model
    init_scaling.propagate_solution(scaled_model, m)

    # Unfix variables to solve optimization model
    ix = m.fs.unit
    ix.bed_depth.unfix()
    ix.bed_diameter.unfix()
    for c in m.fs.set_reactive_ions:
        ix.c_norm[c].unfix()

    # For this example, we are optimizing the model to have an
    # effluent with a very small concentration of the multiple REEs
    @m.fs.unit.Constraint(m.fs.set_reactive_ions)
    def components_specifications(b, c):
        if c == target_component:
            return b.c_norm[c] == 0.999
        else:
            return b.c_norm[c] >= 0.9

    @m.fs.unit.costing.Expression()
    def annualized_capital_cost(b):
        return b.capital_cost / m.fs.ix_lifetime

    m.fs.obj = pyo.Objective(
        expr=(
            ix.costing.annualized_capital_cost
            + ix.costing.fixed_operating_cost
            - ix.expected_annual_profit_ree
        )
    )

    solver = get_solver()
    results_opt = solver.solve(m)
    pyo.assert_optimal_termination(results_opt)

    # Test expected results for ix variables. NOTE: bed_expansion_frac
    # and pressure_drop are the ones calculated with resin-specific
    # data from .json file
    ix_vars_results = {
        "resin_diam": 0.00075,
        "resin_density": 1126.61,
        "bed_volume": 0.0001767144932687618,
        "bed_volume_total": 0.000176714586764426,
        "bed_depth": 2.249999764674673,
        "column_height": 3.9094746956963786,
        "bed_porosity": 0.8,
        "bed_diameter": 0.010000000000242665,
        "ebct": 10581.690809389767,
        "loading_rate": 0.00021263142539510235,
        "service_flow_rate": 0.3402102806486752,
        "target_breakthrough_time": 2291425.7789889984,
    }

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)

    # Expected results for system costs
    sys_cost_results = {
        "wacc": 0.09307339771758533,
        "aggregate_capital_cost": 2396.1687043450247,
        "aggregate_fixed_operating_cost": 15.615588474344836,
        "aggregate_variable_operating_cost": 0,
        "aggregate_flow_electricity": 0,
        "total_operating_cost": 87.50064960469557,
    }

    for v, r in sys_cost_results.items():
        mv = getattr(m.fs.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)

    # Expected results for ion exchange costs considering "single_use"
    # for regenerant use
    ix_cost_results = {
        "capital_cost": 2396.1687043450247,
        "fixed_operating_cost": 15.615588474344836,
        "capital_cost_resin": 1.1338614721473923,
        "capital_cost_regen_tank": 0,
        "flow_vol_resin": 0.0024337195258175713,
        "single_use_resin_replacement_cost": 15.615588474344836,
    }

    for v, r in ix_cost_results.items():
        mv = getattr(m.fs.unit.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)
