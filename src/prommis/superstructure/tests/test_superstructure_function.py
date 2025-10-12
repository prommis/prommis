#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Param,
    RangeSet,
    Set,
    SolverFactory,
    Var,
    assert_optimal_termination,
    value,
)

from idaes.core.scaling import get_scaling_factor
from idaes.core.solvers import get_solver

import pytest

from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice
from prommis.superstructure.report_superstructure_results import (
    report_superstructure_costing,
    report_superstructure_environmental_impacts,
    report_superstructure_results_overview,
    report_superstructure_streams,
)
from prommis.superstructure.superstructure_function import (
    SuperstructureScaler,
    build_model,
    define_custom_units,
)

solver_available = SolverFactory("gurobi").available()
if solver_available:
    solver = get_solver(solver="gurobi")
else:
    solver = None

### Define common parameters
# define custom units
define_custom_units()
obj_func = ObjectiveFunctionChoice.NET_PRESENT_VALUE

plant_start = 2024
plant_lifetime = 5

available_feed = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
}
collection_rate = 0.1
tracked_comps = ["Nd", "Fe"]
prod_comp_mass = {
    "Nd": 0.206 * 3,
    "Fe": 0.691 * 3,
}

num_stages = 3
options_in_stage = {1: 1, 2: 2, 3: 3}
option_outlets = {
    (1, 1): [1, 2],
    (2, 1): [1],
    (2, 2): [2, 3],
}

option_efficiencies = {
    (1, 1): {"Nd": 1, "Fe": 1},
    (2, 1): {"Nd": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Fe": 1},
    (3, 1): {"Nd": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Fe": 0},
    (3, 3): {"Nd": 1, "Fe": 0},
}

profit = {
    (3, 1): {"Nd": 45.4272, "Fe": 0},
    (3, 2): {"Nd": 69.888, "Fe": 0},
    (3, 3): {"Nd": 45.4272, "Fe": 0},
}
opt_var_oc_params = {
    (2, 1): {"a": 0.0053, "b": 7929.7},
    (2, 2): {"a": 0.0015, "b": 2233.16},
    (3, 1): {"a": 15.594, "b": 4e6},
    (3, 2): {"a": 35.58463, "b": 4e6},
    (3, 3): {"a": 1.58, "b": 0},
}
operators_per_discrete_unit = {(1, 1): 1}
yearly_cost_per_unit = {(1, 1): 0}
capital_cost_per_unit = {(1, 1): 0}
processing_rate = {(1, 1): 7868}
num_operators = {
    (2, 1): 0.65,
    (2, 2): 0.65,
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
}
labor_rate = 8000 * 38.20

discretized_purchased_equipment_cost = {
    (2, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            10130.08515,
            31353.21173,
            48788.84678,
            60305.81927,
            77063.4884,
            117214.7546,
            151018.0699,
            195698.5419,
        ],
    },
    (2, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            11702.08515,
            39023.21173,
            62134.84678,
            77539.81927,
            100113.4884,
            154792.7546,
            201326.0699,
            263374.5419,
        ],
    },
    (3, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            343228.652,
            482425.4684,
            618182.0594,
            743750.2902,
            844443.0443,
            978479.5225,
            1183834.522,
            1440660.587,
        ],
    },
    (3, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            643228.652,
            782425.4684,
            918182.0594,
            1043750.2902,
            1144443.0443,
            1278479.5225,
            1483834.522,
            1740660.587,
        ],
    },
    (3, 3): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            4274.7216,
            30479.121,
            53459.01,
            69261.68,
            92510.61,
            143803.33,
            197644.75,
            261513.79,
        ],
    },
}

options_environmental_impacts = {
    (1, 1): 0,
    (2, 1): 1000,
    (2, 2): 0,
    (3, 1): 600,
    (3, 2): 600,
    (3, 3): 0,
}
epsilon = 1

byproduct_values = {
    "Jarosite": -0.17,
    "Iron oxide": 10,
    "Residue": -0.17,
}
byproduct_opt_conversions = {
    (3, 1): {"Jarosite": 0.75},
    (3, 2): {"Iron oxide": 1},
    (3, 3): {"Residue": 0.25},
}


class TestNPV(object):
    @pytest.fixture(scope="class")
    def NPV_model(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        return m

    def test_build(self, NPV_model):
        assert isinstance(NPV_model, ConcreteModel)

        # add_objective_function_choice_param
        assert isinstance(NPV_model.fs, Block)
        assert isinstance(NPV_model.fs.objective_function_choice, Param)

        # add_plant_lifetime_params
        assert isinstance(NPV_model.fs.plant_start, Param)
        assert isinstance(NPV_model.fs.plant_end, Param)
        assert isinstance(NPV_model.fs.plant_life_range, RangeSet)
        assert isinstance(NPV_model.fs.operational_range, RangeSet)

        # add_feed_params
        assert isinstance(NPV_model.fs.tracked_comps, Set)
        assert isinstance(NPV_model.fs.prod_comp_mass, Param)
        assert isinstance(NPV_model.fs.feed_entering, Param)
        assert isinstance(NPV_model.fs.max_feed_entering, Param)

        # add_supe_formulation_params
        assert isinstance(NPV_model.fs.max_options, RangeSet)
        assert isinstance(NPV_model.fs.num_stages, Param)
        assert isinstance(NPV_model.fs.stages_set, RangeSet)
        assert isinstance(NPV_model.fs.options_in_stage, Param)
        assert isinstance(NPV_model.fs.all_opts_set, Set)
        assert isinstance(NPV_model.fs.discrete_opts_set, Set)
        assert isinstance(NPV_model.fs.continuous_opts_set, Set)
        assert isinstance(NPV_model.fs.option_outlet_pairs, Set)
        assert isinstance(NPV_model.fs.option_efficiencies, Param)
        assert isinstance(NPV_model.fs.final_opts_set, Set)

        # add_operating_params
        assert isinstance(NPV_model.fs.costing, Block)
        assert isinstance(NPV_model.fs.costing.profit, Param)
        assert isinstance(NPV_model.fs.costing.opt_var_oc_param_A, Param)
        assert isinstance(NPV_model.fs.costing.opt_var_oc_param_B, Param)
        assert isinstance(NPV_model.fs.costing.operators_per_discrete_unit, Param)
        assert isinstance(NPV_model.fs.costing.yearly_cost_per_unit, Param)
        assert isinstance(NPV_model.fs.costing.capital_cost_per_unit, Param)
        assert isinstance(NPV_model.fs.costing.num_operators, Param)
        assert isinstance(NPV_model.fs.costing.labor_rate, Param)
        assert isinstance(NPV_model.fs.costing.discrete_units_per_option, Param)

        # add_discretized_costing_params
        assert isinstance(NPV_model.fs.costing.discrete_units_per_option, Param)
        assert isinstance(
            NPV_model.fs.continuous_opts_discretized_costing_data_points, Set
        )
        assert isinstance(NPV_model.fs.costing.flowrates_data, Param)
        assert isinstance(NPV_model.fs.costing.costs_data, Param)

        # add_mass_balance_params
        assert isinstance(NPV_model.fs.f_stages, RangeSet)
        assert isinstance(NPV_model.fs.big_m_val, Param)
        assert isinstance(NPV_model.fs.max_flow_upper_bound, Param)

        # add_mass_balance_vars
        assert isinstance(NPV_model.fs.f, Var)
        assert isinstance(NPV_model.fs.f_in, Var)
        assert isinstance(NPV_model.fs.f_out, Var)
        assert isinstance(NPV_model.fs.option_binary_var, Var)
        assert isinstance(NPV_model.fs.piecewise_flow_entering, Var)

        # add_mass_balance_cons
        assert isinstance(NPV_model.fs.inlet_flow_cons, Constraint)
        assert isinstance(NPV_model.fs.init_flow_cons, Constraint)
        assert isinstance(NPV_model.fs.intermediate_flow_cons, Constraint)
        assert isinstance(NPV_model.fs.outlet_flow_cons, Constraint)
        assert isinstance(NPV_model.fs.stage_binary_cons, Constraint)
        assert isinstance(NPV_model.fs.connection_binary_cons, Constraint)
        assert isinstance(NPV_model.fs.f_in_big_m_cons, Constraint)
        assert isinstance(NPV_model.fs.f_out_big_m_cons, Constraint)
        assert isinstance(NPV_model.fs.max_flow_entering_cons, Constraint)

        # add_costing_params
        assert isinstance(NPV_model.fs.costing.lang_factor, Param)
        assert isinstance(NPV_model.fs.costing.i_operating_expense_escalation, Param)
        assert isinstance(NPV_model.fs.costing.i_capital_escalation, Param)
        assert isinstance(NPV_model.fs.costing.discount_factor, Param)
        assert isinstance(NPV_model.fs.costing.financing_factor, Param)
        assert isinstance(NPV_model.fs.costing.other_costs_factor, Param)
        assert isinstance(NPV_model.fs.costing.m_and_sm_costing_factor, Param)
        assert isinstance(NPV_model.fs.costing.sa_and_qa_qc_costing_factor, Param)
        assert isinstance(NPV_model.fs.costing.s_ip_r_and_d_costing_factor, Param)
        assert isinstance(NPV_model.fs.costing.a_and_sl_costing_factor, Param)
        assert isinstance(NPV_model.fs.costing.fb_costing_factor, Param)
        assert isinstance(NPV_model.fs.costing.pt_and_i_costing_factor, Param)
        assert isinstance(
            NPV_model.fs.costing.total_overnight_capital_fraction_expended, Param
        )
        assert isinstance(NPV_model.fs.costing.plant_overhead_factor, Param)

        # add_costing_vars
        assert isinstance(NPV_model.fs.costing.opt_profit, Var)
        assert isinstance(NPV_model.fs.costing.net_present_value, Var)
        assert isinstance(NPV_model.fs.costing.main_product_profit, Var)
        assert isinstance(NPV_model.fs.costing.total_profit, Var)
        assert isinstance(NPV_model.fs.costing.piecewise_equipment_cost, Var)
        assert isinstance(NPV_model.fs.costing.equipment_cost, Var)
        assert isinstance(NPV_model.fs.costing.total_plant_cost, Var)
        assert isinstance(NPV_model.fs.costing.financing, Var)
        assert isinstance(NPV_model.fs.costing.other_costs, Var)
        assert isinstance(NPV_model.fs.costing.total_overnight_cost, Var)
        assert isinstance(NPV_model.fs.costing.opt_variable_operating_cost, Var)
        assert isinstance(NPV_model.fs.costing.aggregate_variable_operating_cost, Var)
        assert isinstance(NPV_model.fs.costing.operators_per_option, Var)
        assert isinstance(NPV_model.fs.costing.total_operators, Var)
        assert isinstance(NPV_model.fs.costing.cost_of_labor, Var)
        assert isinstance(NPV_model.fs.costing.m_and_sm, Var)
        assert isinstance(NPV_model.fs.costing.sa_and_qa_qc, Var)
        assert isinstance(NPV_model.fs.costing.s_ip_r_and_d, Var)
        assert isinstance(NPV_model.fs.costing.a_and_sl, Var)
        assert isinstance(NPV_model.fs.costing.fb, Var)
        assert isinstance(NPV_model.fs.costing.pt_and_i, Var)
        assert isinstance(NPV_model.fs.costing.aggregate_fixed_operating_cost, Var)
        assert isinstance(NPV_model.fs.costing.total_overnight_cost_expended, Var)
        assert isinstance(NPV_model.fs.costing.plant_overhead, Var)
        assert isinstance(NPV_model.fs.costing.total_operating_expense, Var)
        assert isinstance(NPV_model.fs.costing.cash_flow, Var)

        # add_profit_cons
        assert isinstance(
            NPV_model.fs.costing.calculate_final_opts_profit_cons, Constraint
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_main_product_profit_cons, Constraint
        )
        assert isinstance(NPV_model.fs.costing.calculate_total_profit_cons, Constraint)

        # add_capital_cost_cons
        assert isinstance(
            NPV_model.fs.costing.discrete_opts_equipment_cost_cons, Constraint
        )
        assert isinstance(NPV_model.fs.costing.piecewise_cons, Block)
        assert isinstance(NPV_model.fs.costing.add_units_to_piecewise_costs, Constraint)
        assert isinstance(
            NPV_model.fs.costing.calculate_total_plant_cost_con, Constraint
        )
        assert isinstance(NPV_model.fs.costing.calculate_financing_cost_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_other_costs_con, Constraint)
        assert isinstance(
            NPV_model.fs.costing.calculate_total_overnight_cost_con, Constraint
        )

        # add_operating_cost_cons
        assert isinstance(
            NPV_model.fs.costing.calculate_opt_yearly_variable_expense_cons, Constraint
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_total_yearly_variable_operating_costs_cons,
            Constraint,
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_operators_per_option_cons, Constraint
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_total_operators_cons, Constraint
        )
        assert isinstance(NPV_model.fs.costing.calculate_cost_of_labor_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_m_and_sm_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_sa_and_qa_qc_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_s_ip_r_and_d_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_a_and_sl_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_fb_con, Constraint)
        assert isinstance(NPV_model.fs.costing.calculate_pt_and_i_con, Constraint)
        assert isinstance(
            NPV_model.fs.costing.calculate_total_yearly_fixed_operating_costs_cons,
            Constraint,
        )

        # add_cash_flow_cons
        assert isinstance(
            NPV_model.fs.costing.calculate_total_overnight_cost_expended_cons,
            Constraint,
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_plant_overhead_cons, Constraint
        )
        assert isinstance(
            NPV_model.fs.costing.calculate_total_operating_expense_cons, Constraint
        )
        assert isinstance(NPV_model.fs.costing.calculate_cash_flows, Constraint)

        # add_costing_objective_functions
        assert isinstance(
            NPV_model.fs.costing.calculate_net_present_value_con, Constraint
        )

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, NPV_model):
        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(NPV_model)

        results = solver.solve(NPV_model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, NPV_model):
        ## Store correct values for all variables
        # flow entering each stage (except the last stage).
        f_dict = {
            (1, "Fe", 2025): 60173.5929,
            (1, "Fe", 2026): 56934.5304,
            (1, "Fe", 2027): 59393.9376,
            (1, "Fe", 2028): 101124.8787,
            (1, "Nd", 2025): 17938.8714,
            (1, "Nd", 2026): 16973.2464,
            (1, "Nd", 2027): 17706.441600000002,
            (1, "Nd", 2028): 30147.214200000002,
            (2, "Fe", 2025): 60173.5929,
            (2, "Fe", 2026): 56934.5304,
            (2, "Fe", 2027): 59393.9376,
            (2, "Fe", 2028): 101124.8787,
            (2, "Nd", 2025): 17938.8714,
            (2, "Nd", 2026): 16973.2464,
            (2, "Nd", 2027): 17706.441600000002,
            (2, "Nd", 2028): 30147.214200000002,
        }
        # flow entering each option.
        f_in_dict = {
            (1, 1, "Fe", 2025): 60173.5929,
            (1, 1, "Fe", 2026): 56934.5304,
            (1, 1, "Fe", 2027): 59393.9376,
            (1, 1, "Fe", 2028): 101124.8787,
            (1, 1, "Nd", 2025): 17938.8714,
            (1, 1, "Nd", 2026): 16973.2464,
            (1, 1, "Nd", 2027): 17706.441600000002,
            (1, 1, "Nd", 2028): 30147.214200000002,
            (2, 1, "Fe", 2025): 0.0,
            (2, 1, "Fe", 2026): 0.0,
            (2, 1, "Fe", 2027): 0.0,
            (2, 1, "Fe", 2028): 0.0,
            (2, 1, "Nd", 2025): 0.0,
            (2, 1, "Nd", 2026): 0.0,
            (2, 1, "Nd", 2027): 0.0,
            (2, 1, "Nd", 2028): 0.0,
            (2, 2, "Fe", 2025): 60173.592899999996,
            (2, 2, "Fe", 2026): 56934.530399999996,
            (2, 2, "Fe", 2027): 59393.937600000005,
            (2, 2, "Fe", 2028): 101124.8787,
            (2, 2, "Nd", 2025): 17938.871400000004,
            (2, 2, "Nd", 2026): 16973.246400000004,
            (2, 2, "Nd", 2027): 17706.441600000006,
            (2, 2, "Nd", 2028): 30147.214200000002,
            (3, 1, "Fe", 2025): 0.0,
            (3, 1, "Fe", 2026): 0.0,
            (3, 1, "Fe", 2027): 0.0,
            (3, 1, "Fe", 2028): 0.0,
            (3, 1, "Nd", 2025): 0.0,
            (3, 1, "Nd", 2026): 0.0,
            (3, 1, "Nd", 2027): 0.0,
            (3, 1, "Nd", 2028): 0.0,
            (3, 2, "Fe", 2025): 0.0,
            (3, 2, "Fe", 2026): 0.0,
            (3, 2, "Fe", 2027): 0.0,
            (3, 2, "Fe", 2028): 0.0,
            (3, 2, "Nd", 2025): 0.0,
            (3, 2, "Nd", 2026): 0.0,
            (3, 2, "Nd", 2027): 0.0,
            (3, 2, "Nd", 2028): 0.0,
            (3, 3, "Fe", 2025): 60173.592899999974,
            (3, 3, "Fe", 2026): 56934.530399999974,
            (3, 3, "Fe", 2027): 59393.937599999976,
            (3, 3, "Fe", 2028): 101124.87869999999,
            (3, 3, "Nd", 2025): 17938.871399999993,
            (3, 3, "Nd", 2026): 16973.246399999996,
            (3, 3, "Nd", 2027): 17706.441599999995,
            (3, 3, "Nd", 2028): 30147.214200000002,
        }
        # flow exiting each option.
        f_out_dict = {
            (1, 1, "Fe", 2025): 60173.5929,
            (1, 1, "Fe", 2026): 56934.5304,
            (1, 1, "Fe", 2027): 59393.9376,
            (1, 1, "Fe", 2028): 101124.8787,
            (1, 1, "Nd", 2025): 17938.8714,
            (1, 1, "Nd", 2026): 16973.2464,
            (1, 1, "Nd", 2027): 17706.441600000002,
            (1, 1, "Nd", 2028): 30147.214200000002,
            (2, 1, "Fe", 2025): 0.0,
            (2, 1, "Fe", 2026): 0.0,
            (2, 1, "Fe", 2027): 0.0,
            (2, 1, "Fe", 2028): 0.0,
            (2, 1, "Nd", 2025): 0.0,
            (2, 1, "Nd", 2026): 0.0,
            (2, 1, "Nd", 2027): 0.0,
            (2, 1, "Nd", 2028): 0.0,
            (2, 2, "Fe", 2025): 60173.592899999996,
            (2, 2, "Fe", 2026): 56934.530399999996,
            (2, 2, "Fe", 2027): 59393.937600000005,
            (2, 2, "Fe", 2028): 101124.8787,
            (2, 2, "Nd", 2025): 17938.871400000004,
            (2, 2, "Nd", 2026): 16973.246400000004,
            (2, 2, "Nd", 2027): 17706.441600000006,
            (2, 2, "Nd", 2028): 30147.214200000002,
            (3, 1, "Fe", 2025): 0.0,
            (3, 1, "Fe", 2026): 0.0,
            (3, 1, "Fe", 2027): 0.0,
            (3, 1, "Fe", 2028): 0.0,
            (3, 1, "Nd", 2025): 0.0,
            (3, 1, "Nd", 2026): 0.0,
            (3, 1, "Nd", 2027): 0.0,
            (3, 1, "Nd", 2028): 0.0,
            (3, 2, "Fe", 2025): 0.0,
            (3, 2, "Fe", 2026): 0.0,
            (3, 2, "Fe", 2027): 0.0,
            (3, 2, "Fe", 2028): 0.0,
            (3, 2, "Nd", 2025): 0.0,
            (3, 2, "Nd", 2026): 0.0,
            (3, 2, "Nd", 2027): 0.0,
            (3, 2, "Nd", 2028): 0.0,
            (3, 3, "Fe", 2025): 0.0,
            (3, 3, "Fe", 2026): 0.0,
            (3, 3, "Fe", 2027): 0.0,
            (3, 3, "Fe", 2028): 0.0,
            (3, 3, "Nd", 2025): 17938.871399999993,
            (3, 3, "Nd", 2026): 16973.246399999996,
            (3, 3, "Nd", 2027): 17706.441599999995,
            (3, 3, "Nd", 2028): 30147.214200000002,
        }
        # Binary variables to indicate whether or not an option has been selected.
        option_binary_var_dict = {
            (1, 1): 1.0,
            (2, 1): 0.0,
            (2, 2): 1.0,
            (3, 1): 0.0,
            (3, 2): 0.0,
            (3, 3): 1.0,
        }
        # The max total flow that enters each continuous option over the lifetime of the plant.
        piecewise_flow_entering_dict = {
            (2, 1): 0.0,
            (2, 2): 131272.0929,
            (3, 1): 0.0,
            (3, 2): 0.0,
            (3, 3): 131272.09289999993,
        }
        # Profit generated by each option in final processing stage
        opt_profit_dict = {
            (3, 1, 2025): 0.0,
            (3, 1, 2026): 0.0,
            (3, 1, 2027): 0.0,
            (3, 1, 2028): 0.0,
            (3, 2, 2025): 0.0,
            (3, 2, 2026): 0.0,
            (3, 2, 2027): 0.0,
            (3, 2, 2028): 0.0,
            (3, 3, 2025): 814912.6988620796,
            (3, 3, 2026): 771047.0588620799,
            (3, 3, 2027): 804354.0638515197,
            (3, 3, 2028): 1369503.5289062401,
        }
        # Total profit from main product each year
        main_product_profit_dict = {
            2025: 814912.6988620803,
            2026: 771047.0588620801,
            2027: 804354.0638515205,
            2028: 1369503.5289062406,
        }
        # Total profit generated by plant each year
        total_profit_dict = {
            2025: 814912.6988620803,
            2026: 771047.0588620801,
            2027: 804354.0638515205,
            2028: 1369503.5289062406,
        }
        # Cost of purchased equipment (dimensionless for piecewise constraints)
        piecewise_equipment_cost_dict = {
            (2, 1): 0.0,
            (2, 2): 16034.63796200485,
            (3, 1): 0.0,
            (3, 2): 0.0,
            (3, 3): 8430.18510092596,
        }
        # Cost of purchased equipment with units
        equipment_cost_dict = {
            (1, 1): 0.0,
            (2, 1): 0.0,
            (2, 2): 16034.63796200485,
            (3, 1): 0.0,
            (3, 2): 0.0,
            (3, 3): 8430.18510092596,
        }
        # Yearly variable operating expense for each option
        opt_variable_operating_cost_dict = {
            (1, 1, 2025): 0.0,
            (1, 1, 2026): 0.0,
            (1, 1, 2027): 0.0,
            (1, 1, 2028): 0.0,
            (2, 1, 2025): 0.0,
            (2, 1, 2026): 0.0,
            (2, 1, 2027): 0.0,
            (2, 1, 2028): 0.0,
            (2, 2, 2025): 2350.32869645,
            (2, 2, 2026): 2344.0216652,
            (2, 2, 2027): 2348.8105688,
            (2, 2, 2028): 2430.06813935,
            (3, 1, 2025): 0.0,
            (3, 1, 2026): 0.0,
            (3, 1, 2027): 0.0,
            (3, 1, 2028): 0.0,
            (3, 2, 2025): 0.0,
            (3, 2, 2026): 0.0,
            (3, 2, 2027): 0.0,
            (3, 2, 2028): 0.0,
            (3, 3, 2025): 123417.69359399995,
            (3, 3, 2026): 116774.28734399995,
            (3, 3, 2027): 121818.59913599996,
            (3, 3, 2028): 207409.90678199998,
        }
        # Total yearly variable operating expense
        aggregate_variable_operating_cost_dict = {
            2025: 125768.02229044959,
            2026: 119118.30900919996,
            2027: 124167.4097048007,
            2028: 209839.9749213513,
        }
        # Total yearly fixed operating cost
        aggregate_fixed_operating_cost_dict = {
            2025: 4273448.942723528,
            2026: 4273010.286323528,
            2027: 4273343.356373422,
            2028: 4278994.851023969,
        }
        # Number of operators needed for each option
        operators_per_option_dict = {
            (1, 1): 7.0,
            (2, 1): 0.0,
            (2, 2): 0.65,
            (3, 1): 0.0,
            (3, 2): 0.0,
            (3, 3): 1.3,
        }
        # Sales, IP, R&D expenses
        s_ip_r_and_d_dict = {
            2025: 8149.126988620804,
            2026: 7710.470588620801,
            2027: 8043.5406385152055,
            2028: 13695.035289062405,
        }
        # Total overnight cost expended each year
        total_overnight_cost_expended_dict = {
            2024: 8552.14373328566,
            2025: 51312.86239971396,
            2026: 25656.43119985698,
            2027: 0.0,
            2028: 0.0,
        }
        # Yearly plant overhead
        plant_overhead_dict = {
            2025: 879843.3930027956,
            2026: 878425.7190665458,
            2027: 879502.1532156446,
            2028: 897766.9651890642,
        }
        # Total operating expense each year
        total_operating_expense_dict = {
            2025: 5279060.358016773,
            2026: 5270554.314399274,
            2027: 5277012.919293867,
            2028: 5386601.791134384,
        }
        # Yearly cash flow (negative = investment/cost, positive = profit)
        cash_flow_dict = {
            2024: -8552.14373328566,
            2025: -4651232.214375437,
            2026: -4801064.19238049,
            2027: -4887395.093130949,
            2028: -4521279.484773466,
        }
        # Net present value
        net_present_value = -16440490.52395555
        # Total plant cost as defined by QGESS
        total_plant_cost = 72660.52449690452
        # Total financing cost of the plant
        financing = 1961.834161416422
        # Other costs associated with the plant
        other_costs = 10899.078674535676
        # Total overnight cost of the plant
        total_overnight_cost = 85521.4373328566
        # Total number of operators needed for the process
        total_operators = 9.0
        # Cost of labor for the process
        cost_of_labor = 2750400.0
        # Maintenance & Supply Materials (M&SM)
        m_and_sm = 1453.2104899380904
        # Sample Analysis & Quality Assurance/Quality Control (SA&QA/QC)
        sa_and_qa_qc = 275040.0
        # Administrative & Supporting Labor (A&SL)
        a_and_sl = 550080.0
        # Fringe Benefits (FB)
        fb = 687600.0
        # Property Taxes & Insurance (PT&I)
        pt_and_i = 726.6052449690452

        ## Test Variables
        # test flow entering each stage (except the last stage)
        for key, val in f_dict.items():
            assert value(NPV_model.fs.f[key]) == pytest.approx(val, rel=1e-8)

        # test flow entering each option
        for key, val in f_in_dict.items():
            assert value(NPV_model.fs.f_in[key]) == pytest.approx(val, rel=1e-8)

        # test flow exiting each option
        for key, val in f_out_dict.items():
            assert value(NPV_model.fs.f_out[key]) == pytest.approx(val, rel=1e-8)

        # test binary variables to indicate whether or not an option has been selected
        for key, val in option_binary_var_dict.items():
            assert value(NPV_model.fs.option_binary_var[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test max total flow that enters each continuous option over the lifetime of the plant.
        for key, val in piecewise_flow_entering_dict.items():
            assert value(NPV_model.fs.piecewise_flow_entering[key]) == pytest.approx(
                val, rel=1e-8
            )

            # test profit generated by each option in final processing stage
        for key, val in opt_profit_dict.items():
            assert value(NPV_model.fs.costing.opt_profit[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test total profit from main product each year
        for key, val in main_product_profit_dict.items():
            assert value(
                NPV_model.fs.costing.main_product_profit[key]
            ) == pytest.approx(val, rel=1e-8)

        # test total profit generated by plant each year
        for key, val in total_profit_dict.items():
            assert value(NPV_model.fs.costing.total_profit[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test cost of purchased equipment (dimensionless for piecewise constraints)
        for key, val in piecewise_equipment_cost_dict.items():
            assert value(
                NPV_model.fs.costing.piecewise_equipment_cost[key]
            ) == pytest.approx(val, rel=1e-8)

        # test cost of purchased equipment with units
        for key, val in equipment_cost_dict.items():
            assert value(NPV_model.fs.costing.equipment_cost[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test yearly variable operating expense for each option
        for key, val in opt_variable_operating_cost_dict.items():
            assert value(
                NPV_model.fs.costing.opt_variable_operating_cost[key]
            ) == pytest.approx(val, rel=1e-8)

        # test total yearly variable operating expense
        for key, val in aggregate_variable_operating_cost_dict.items():
            assert value(
                NPV_model.fs.costing.aggregate_variable_operating_cost[key]
            ) == pytest.approx(val, rel=1e-8)

        # test total yearly fixed operating cost
        for key, val in aggregate_fixed_operating_cost_dict.items():
            assert value(
                NPV_model.fs.costing.aggregate_fixed_operating_cost[key]
            ) == pytest.approx(val, rel=1e-8)

        # test number of operators needed for each option
        for key, val in operators_per_option_dict.items():
            assert value(
                NPV_model.fs.costing.operators_per_option[key]
            ) == pytest.approx(val, rel=1e-8)

        # test sales, IP, R&D expenses
        for key, val in s_ip_r_and_d_dict.items():
            assert value(NPV_model.fs.costing.s_ip_r_and_d[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test total overnight cost expended each year
        for key, val in total_overnight_cost_expended_dict.items():
            assert value(
                NPV_model.fs.costing.total_overnight_cost_expended[key]
            ) == pytest.approx(val, rel=1e-8)

        # test yearly plant overhead
        for key, val in plant_overhead_dict.items():
            assert value(NPV_model.fs.costing.plant_overhead[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test total operating expense each year
        for key, val in total_operating_expense_dict.items():
            assert value(
                NPV_model.fs.costing.total_operating_expense[key]
            ) == pytest.approx(val, rel=1e-8)

        # test yearly cash flow
        for key, val in cash_flow_dict.items():
            assert value(NPV_model.fs.costing.cash_flow[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test net present value
        assert value(NPV_model.fs.costing.net_present_value) == pytest.approx(
            net_present_value, rel=1e-8
        )

        # test total plant cost
        assert value(NPV_model.fs.costing.total_plant_cost) == pytest.approx(
            total_plant_cost, rel=1e-8
        )

        # test financing cost
        assert value(NPV_model.fs.costing.financing) == pytest.approx(
            financing, rel=1e-8
        )

        # test other costs
        assert value(NPV_model.fs.costing.other_costs) == pytest.approx(
            other_costs, rel=1e-8
        )

        # test total overnight cost
        assert value(NPV_model.fs.costing.total_overnight_cost) == pytest.approx(
            total_overnight_cost, rel=1e-8
        )

        # test total number of operators
        assert value(NPV_model.fs.costing.total_operators) == pytest.approx(
            total_operators, rel=1e-8
        )

        # test cost of labor
        assert value(NPV_model.fs.costing.cost_of_labor) == pytest.approx(
            cost_of_labor, rel=1e-8
        )

        # test maintenance & supply materials
        assert value(NPV_model.fs.costing.m_and_sm) == pytest.approx(m_and_sm, rel=1e-8)

        # test sample analysis & quality assurance/quality control
        assert value(NPV_model.fs.costing.sa_and_qa_qc) == pytest.approx(
            sa_and_qa_qc, rel=1e-8
        )

        # test administrative & supporting labor
        assert value(NPV_model.fs.costing.a_and_sl) == pytest.approx(a_and_sl, rel=1e-8)

        # test fringe benefits
        assert value(NPV_model.fs.costing.fb) == pytest.approx(fb, rel=1e-8)

        # test property taxes & insurance
        assert value(NPV_model.fs.costing.pt_and_i) == pytest.approx(pt_and_i, rel=1e-8)

        # test objective function
        assert value(NPV_model.fs.costing.obj) == pytest.approx(
            net_present_value, rel=1e-8
        )


class TestByproductValorization(object):
    @pytest.fixture(scope="class")
    def BV_model(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        return m

    def test_build(self, BV_model):
        # add_byproduct_valorization_params
        assert isinstance(BV_model.fs.byproduct_valorization, Block)
        assert isinstance(BV_model.fs.byproduct_valorization.byproducts_set, Set)
        assert isinstance(BV_model.fs.byproduct_valorization.byproduct_opts_set, Set)
        assert isinstance(BV_model.fs.byproduct_valorization.opt_byproduct_set, Set)
        assert isinstance(BV_model.fs.costing.byproduct_values, Param)

        # add_byproduct_valorization_vars
        assert isinstance(BV_model.fs.byproduct_valorization.byproduct_produced, Var)
        assert isinstance(BV_model.fs.costing.byproduct_profit, Var)
        assert isinstance(BV_model.fs.costing.total_byproduct_profit, Var)

        # add_byproduct_valorization_cons
        assert isinstance(
            BV_model.fs.byproduct_valorization.calculate_byproduct_produced_cons,
            Constraint,
        )

        # add_profit_cons
        assert isinstance(
            BV_model.fs.costing.calculate_byproduct_profit_cons, Constraint
        )
        assert isinstance(BV_model.fs.costing.calculate_opt_byprod_val_cons, Constraint)
        assert isinstance(BV_model.fs.costing.calculate_total_profit_cons, Constraint)
        assert isinstance(
            BV_model.fs.byproduct_valorization.calculate_byproduct_produced_cons,
            Constraint,
        )

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, BV_model):
        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(BV_model)

        results = solver.solve(BV_model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, BV_model):
        ## Store correct values for all variables
        # Amount of each byproduct produced each year (kg/a)
        byproduct_produced_dict = {
            ("Iron oxide", 2025): 0.0,
            ("Iron oxide", 2026): 0.0,
            ("Iron oxide", 2027): 0.0,
            ("Iron oxide", 2028): 0.0,
            ("Jarosite", 2025): 0.0,
            ("Jarosite", 2026): 0.0,
            ("Jarosite", 2027): 0.0,
            ("Jarosite", 2028): 0.0,
            ("Residue", 2025): 19528.11607499999,
            ("Residue", 2026): 18476.94419999999,
            ("Residue", 2027): 19275.09479999999,
            ("Residue", 2028): 32818.023225,
        }
        # Amount of profit generated from each byproduct each year (USD/a)
        byproduct_profit_dict = {
            ("Iron oxide", 2025): 0.0,
            ("Iron oxide", 2026): 0.0,
            ("Iron oxide", 2027): 0.0,
            ("Iron oxide", 2028): 0.0,
            ("Jarosite", 2025): 0.0,
            ("Jarosite", 2026): 0.0,
            ("Jarosite", 2027): 0.0,
            ("Jarosite", 2028): 0.0,
            ("Residue", 2025): -3319.779732749999,
            ("Residue", 2026): -3141.080513999999,
            ("Residue", 2027): -3276.766115999999,
            ("Residue", 2028): -5579.06394825,
        }
        # Total profit generated by the plant each year from byproduct valorization (USD/a)
        total_byproduct_profit_dict = {
            2025: -3319.7797327500302,
            2026: -3141.0805140000302,
            2027: -3276.7661160000134,
            2028: -5579.06394824991,
        }
        # objective function value
        objective_value = -16454574.208341438

        ## Test Variables
        # test amount of each byproduct produced each year
        for key, val in byproduct_produced_dict.items():
            assert value(
                BV_model.fs.byproduct_valorization.byproduct_produced[key]
            ) == pytest.approx(val, rel=1e-8)

        # test amount of profit generated from each byproduct each year
        for key, val in byproduct_profit_dict.items():
            assert value(BV_model.fs.costing.byproduct_profit[key]) == pytest.approx(
                val, rel=1e-8
            )

        # test objective function
        assert value(BV_model.fs.costing.obj) == pytest.approx(
            objective_value, rel=1e-8
        )


class TestEnvironmentalImpacts(object):
    @pytest.fixture(scope="class")
    def EI_model(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        return m

    def test_build(self, EI_model):
        # add_environmental_impact_params
        assert isinstance(EI_model.fs.environmental_impacts, Block)
        assert isinstance(
            EI_model.fs.environmental_impacts.options_environmental_impacts, Param
        )
        assert isinstance(EI_model.fs.environmental_impacts.epsilon, Param)

        # add_environmental_impact_vars
        assert isinstance(EI_model.fs.environmental_impacts.option_yearly_impacts, Var)
        assert isinstance(EI_model.fs.environmental_impacts.total_yearly_impacts, Var)
        assert isinstance(EI_model.fs.environmental_impacts.total_impacts, Var)

        # add_environmental_impact_cons
        assert isinstance(
            EI_model.fs.environmental_impacts.calculate_opt_yearly_impacts_con,
            Constraint,
        )
        assert isinstance(
            EI_model.fs.environmental_impacts.calculate_yearly_impacts_con, Constraint
        )
        assert isinstance(
            EI_model.fs.environmental_impacts.calculate_total_impacts_con, Constraint
        )
        assert isinstance(EI_model.fs.environmental_impacts.epsilon_con, Constraint)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, EI_model):
        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(EI_model)

        results = solver.solve(EI_model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, EI_model):
        ## Store correct values for all variables
        # Yearly environmental impacts generated by each option (1/a)
        option_yearly_impacts_dict = {
            (1, 1, 2025): 0.0,
            (1, 1, 2026): 0.0,
            (1, 1, 2027): 0.0,
            (1, 1, 2028): 0.0,
            (2, 1, 2025): 0.0,
            (2, 1, 2026): 0.0,
            (2, 1, 2027): 0.0,
            (2, 1, 2028): 0.0,
            (2, 2, 2025): 0.0,
            (2, 2, 2026): 0.0,
            (2, 2, 2027): 0.0,
            (2, 2, 2028): 0.0,
            (3, 1, 2025): 0.0,
            (3, 1, 2026): 0.0,
            (3, 1, 2027): 0.0,
            (3, 1, 2028): 0.0,
            (3, 2, 2025): 0.0,
            (3, 2, 2026): 0.0,
            (3, 2, 2027): 0.0,
            (3, 2, 2028): 0.0,
            (3, 3, 2025): 0.0,
            (3, 3, 2026): 0.0,
            (3, 3, 2027): 0.0,
            (3, 3, 2028): 0.0,
        }
        # Total yearly environmental impacts generated by the process (1/a)
        total_yearly_impacts_dict = {2025: 0.0, 2026: 0.0, 2027: 0.0, 2028: 0.0}
        # The total environmental impacts generated by the process over its entire lifetime.
        total_impacts = 0
        # objective function value
        objective_value = -16454574.208341438

        ## Test Variables
        # test yearly environmental impacts generated by each option
        for key, val in option_yearly_impacts_dict.items():
            assert value(
                EI_model.fs.environmental_impacts.option_yearly_impacts[key]
            ) == pytest.approx(val, rel=1e-8)

        # test total yearly environmental impacts generated by the process
        for key, val in total_yearly_impacts_dict.items():
            assert value(
                EI_model.fs.environmental_impacts.total_yearly_impacts[key]
            ) == pytest.approx(val, rel=1e-8)

        # test total environmental impacts generated by the process over its entire lifetime
        assert value(EI_model.fs.environmental_impacts.total_impacts) == pytest.approx(
            total_impacts, rel=1e-8
        )

        # test objective function
        assert value(EI_model.fs.costing.obj) == pytest.approx(
            objective_value, rel=1e-8
        )

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_report_superstructure_results(self, EI_model):
        results = solver.solve(EI_model, tee=True)

        report_superstructure_results_overview(EI_model, results)
        report_superstructure_costing(EI_model, results)
        report_superstructure_streams(EI_model, results)
        report_superstructure_environmental_impacts(EI_model, results)


class TestCOR(object):
    @pytest.fixture(scope="class")
    def COR_model(self):
        m = build_model(
            ### Choice of objective function
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        return m

    def test_build(self, COR_model):
        # add_costing_vars
        assert isinstance(COR_model.fs.costing.cost_of_recovery, Var)

        # add_profit_cons
        assert isinstance(
            COR_model.fs.costing.calculate_main_product_profit_cons, Constraint
        )
        assert isinstance(COR_model.fs.costing.calculate_total_profit_cons, Constraint)

        # add_costing_objective_functions
        assert isinstance(
            COR_model.fs.costing.set_net_present_value_to_zero_con, Constraint
        )

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, COR_model):
        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(COR_model)

        results = solver.solve(COR_model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, COR_model):
        ## Store correct values for all variables
        # Cost of recovery
        cost_of_recovery = 261.6415514638381
        # objective function value
        objective_value = 261.6415514638381

        ## Test Variables
        assert value(COR_model.fs.costing.cost_of_recovery) == pytest.approx(
            cost_of_recovery, rel=1e-8
        )
        assert value(COR_model.fs.costing.obj) == pytest.approx(
            objective_value, rel=1e-8
        )

    @pytest.fixture(scope="class")
    def COR_model2(self):
        m = build_model(
            ### Choice of objective function
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        return m

    def test_build2(self, COR_model2):
        # add_profit_cons
        assert isinstance(COR_model2.fs.costing.calculate_total_profit_cons, Constraint)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve2(self, COR_model2):
        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(COR_model2)

        results2 = solver.solve(COR_model2, tee=True)
        assert_optimal_termination(results2)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution2(self, COR_model2):
        ## Store correct values for all variables
        # Cost of recovery
        cost_of_recovery = 261.45649078422645
        # objective function value
        objective_value = 261.45649078422645

        ## Test Variables
        assert value(COR_model2.fs.costing.cost_of_recovery) == pytest.approx(
            cost_of_recovery, rel=1e-8
        )
        assert value(COR_model2.fs.costing.obj) == pytest.approx(
            objective_value, rel=1e-8
        )


class TestSuperstructureScaler:
    """Test class for SuperstructureScaler functionality."""

    @pytest.fixture
    def sample_model(self):
        """Create a sample model for testing scaling."""
        # Use your existing model building function with minimal parameters
        model = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        return model

    def test_scaler_initialization(self):
        """Test that SuperstructureScaler can be initialized."""
        scaler = SuperstructureScaler()
        assert isinstance(scaler, SuperstructureScaler)

    def test_scale_model_runs_without_error(self, sample_model):
        """Test that scale_model runs without throwing errors."""
        scaler = SuperstructureScaler()

        # This should not raise any exceptions
        scaler.scale_model(sample_model)

    def test_flow_variable_scaling(self, sample_model):
        """Test that flow variables are scaled correctly."""
        scaler = SuperstructureScaler()
        scaler.scale_model(sample_model)

        # Check that flow variables have the expected scaling factors
        for var in sample_model.fs.f.values():
            scaling_factor = get_scaling_factor(var)
            assert (
                scaling_factor == 1e-4
            ), f"Flow variable scaling should be 1e-4, got {scaling_factor}"

        for var in sample_model.fs.f_in.values():
            scaling_factor = get_scaling_factor(var)
            assert (
                scaling_factor == 1e-4
            ), f"Flow input variable scaling should be 1e-4, got {scaling_factor}"

        for var in sample_model.fs.f_out.values():
            scaling_factor = get_scaling_factor(var)
            assert (
                scaling_factor == 1e-4
            ), f"Flow output variable scaling should be 1e-4, got {scaling_factor}"

    def test_piecewise_flow_scaling(self, sample_model):
        """Test that piecewise flow variables are scaled correctly."""
        scaler = SuperstructureScaler()
        scaler.scale_model(sample_model)

        for var in sample_model.fs.piecewise_flow_entering.values():
            scaling_factor = get_scaling_factor(var)
            assert (
                scaling_factor == 1e-5
            ), f"Piecewise flow scaling should be 1e-5, got {scaling_factor}"

    def test_costing_variable_scaling(self, sample_model):
        """Test that costing variables are scaled correctly."""
        scaler = SuperstructureScaler()
        scaler.scale_model(sample_model)

        # Test specific costing variable scaling factors
        expected_scalings = {
            "net_present_value": 1e-6,
            "main_product_profit": 1e-6,
            "total_profit": 1e-6,
            "piecewise_equipment_cost": 1e-5,
            "equipment_cost": 1e-4,
            "total_plant_cost": 1e-6,
            "financing": 1e-4,
            "other_costs": 1e-5,
            "total_overnight_cost": 1e-5,
            "opt_variable_operating_cost": 1e-3,
            "aggregate_variable_operating_cost": 1e-5,
            "operators_per_option": 1e1,
            "cost_of_labor": 1e-5,
        }

        for var_name, expected_scaling in expected_scalings.items():
            if hasattr(sample_model.fs.costing, var_name):
                var = getattr(sample_model.fs.costing, var_name)
                if hasattr(var, "values"):  # Indexed variable
                    for v in var.values():
                        scaling_factor = get_scaling_factor(v)
                        assert (
                            scaling_factor == expected_scaling
                        ), f"{var_name} scaling should be {expected_scaling}, got {scaling_factor}"
                else:  # Scalar variable
                    scaling_factor = get_scaling_factor(var)
                    assert (
                        scaling_factor == expected_scaling
                    ), f"{var_name} scaling should be {expected_scaling}, got {scaling_factor}"

    def test_constraint_scaling(self, sample_model):
        """Test that constraints are scaled using the nominal value method."""
        scaler = SuperstructureScaler()
        scaler.scale_model(sample_model)

        # Check that constraints have scaling factors applied
        # Test a few key constraints
        constraint_sets = [
            "inlet_flow_cons",
            "init_flow_cons",
            "intermediate_flow_cons",
            "outlet_flow_cons",
        ]

        for constraint_name in constraint_sets:
            if hasattr(sample_model.fs, constraint_name):
                constraint_set = getattr(sample_model.fs, constraint_name)
                for constraint in constraint_set.values():
                    scaling_factor = get_scaling_factor(constraint)
                    assert (
                        scaling_factor is not None
                    ), f"Constraint {constraint_name} should have a scaling factor"

    def test_environmental_impacts_scaling(self):
        """Test scaling when environmental impacts are considered."""
        model_with_env = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_with_env)

        # Test environmental impact variable scaling
        if hasattr(model_with_env.fs, "environmental_impacts"):
            for (
                var
            ) in model_with_env.fs.environmental_impacts.option_yearly_impacts.values():
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1e-7

    def test_byproduct_valorization_scaling(self):
        """Test scaling when byproduct valorization is considered."""
        model_with_byproducts = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_with_byproducts)

        # Test byproduct variable scaling
        if hasattr(model_with_byproducts.fs, "byproduct_valorization"):
            for (
                var
            ) in (
                model_with_byproducts.fs.byproduct_valorization.byproduct_produced.values()
            ):
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1e-4

    def test_cost_of_recovery_variable_scaling(self):
        """Test scaling when Cost of Recovery is chosen as objective function."""
        model_with_cor = build_model(
            ### Choice of objective function - Use COR instead of NPV
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,  # This triggers the COR-specific scaling
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_with_cor)

        # Test that cost_of_recovery variable is scaled with factor 1
        if hasattr(model_with_cor.fs.costing, "cost_of_recovery"):
            for var in model_with_cor.fs.costing.cost_of_recovery.values():
                scaling_factor = get_scaling_factor(var)
                assert (
                    scaling_factor == 1
                ), f"Cost of recovery scaling should be 1, got {scaling_factor}"

    def test_npv_zero_constraint_scaling(self):
        """Test scaling of NPV zero constraint when COR is objective function."""
        model_with_cor = build_model(
            ### Choice of objective function - Use COR instead of NPV
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,  # This triggers the NPV zero constraint
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_with_cor)

        # Test that set_net_present_value_to_zero_con constraint is scaled
        if hasattr(model_with_cor.fs.costing, "set_net_present_value_to_zero_con"):
            for (
                constraint
            ) in model_with_cor.fs.costing.set_net_present_value_to_zero_con.values():
                scaling_factor = get_scaling_factor(constraint)
                assert (
                    scaling_factor is not None
                ), "NPV zero constraint should have a scaling factor"

    def test_combined_cor_with_environmental_impacts(self):
        """Test COR objective function combined with environmental impacts."""
        model_cor_env = build_model(
            ### Choice of objective function
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,  # COR objective
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_cor_env)

        # Test both COR-specific and environmental impact scaling
        if hasattr(model_cor_env.fs.costing, "cost_of_recovery"):
            for var in model_cor_env.fs.costing.cost_of_recovery.values():
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1

        if hasattr(model_cor_env.fs, "environmental_impacts"):
            for (
                var
            ) in model_cor_env.fs.environmental_impacts.option_yearly_impacts.values():
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1e-7

    def test_combined_cor_with_byproducts(self):
        """Test COR objective function combined with byproduct valorization."""
        model_cor_byproduct = build_model(
            ### Choice of objective function
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,  # COR objective
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        scaler = SuperstructureScaler()
        scaler.scale_model(model_cor_byproduct)

        # Test both COR-specific and byproduct scaling
        if hasattr(model_cor_byproduct.fs.costing, "cost_of_recovery"):
            for var in model_cor_byproduct.fs.costing.cost_of_recovery.values():
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1

        if hasattr(model_cor_byproduct.fs, "byproduct_valorization"):
            for (
                var
            ) in (
                model_cor_byproduct.fs.byproduct_valorization.byproduct_produced.values()
            ):
                scaling_factor = get_scaling_factor(var)
                assert scaling_factor == 1e-4
