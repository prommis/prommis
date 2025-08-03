# import pyomo.environ as pyo

# from prommis.superstructure.add_superstructure_blocks import (
#     add_byproduct_valorization_cons,
#     add_byproduct_valorization_params,
#     add_byproduct_valorization_vars,
#     add_capital_cost_cons,
#     add_cash_flow_cons,
#     add_costing_objective_functions,
#     add_costing_params,
#     add_costing_vars,
#     add_discretized_costing_params,
#     add_environmental_impact_cons,
#     add_environmental_impact_params,
#     add_environmental_impact_vars,
#     add_feed_params,
#     add_mass_balance_cons,
#     add_mass_balance_params,
#     add_mass_balance_vars,
#     add_objective_function_choice_param,
#     add_operating_cost_cons,
#     add_operating_params,
#     add_plant_lifetime_params,
#     add_profit_cons,
#     add_supe_formulation_params,
# )
# from prommis.superstructure.check_superstructure_inputs import (
#     check_byproduct_valorization_params,
#     check_discretized_costing_params,
#     check_environmental_impact_params,
#     check_feed_params,
#     check_objective_function_choice,
#     check_operating_params,
#     check_plant_lifetime_params,
#     check_supe_formulation_params,
# )
# from prommis.superstructure.superstructure_function import ObjectiveFunction

npv_func = ObjectiveFunction.net_present_value

check_objective_function_choice(npv_func)

m = pyo.ConcreteModel()

add_objective_function_choice_param(m, npv_func)

m.fs.objective_function_choice.display()
