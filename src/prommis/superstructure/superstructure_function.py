#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Superstructure Code
===================

Author: Chris Laliwala
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from prommis.superstructure.add_superstructure_blocks import (
    add_byproduct_valorization_cons,
    add_byproduct_valorization_params,
    add_byproduct_valorization_vars,
    add_capital_cost_cons,
    add_cash_flow_cons,
    add_costing_objective_functions,
    add_costing_params,
    add_costing_vars,
    add_discretized_costing_params,
    add_environmental_impact_cons,
    add_environmental_impact_params,
    add_environmental_impact_vars,
    add_feed_params_block,
    add_mass_balance_cons,
    add_mass_balance_params,
    add_mass_balance_vars,
    add_operating_cost_cons,
    add_operating_params,
    add_plant_lifetime_params_block,
    add_profit_cons,
    add_supe_formulation_params,
)
from prommis.superstructure.check_superstructure_inputs import (
    check_byproduct_valorization_params,
    check_discretized_costing_params,
    check_environmental_impact_params,
    check_feed_params,
    check_objective_function_choice,
    check_operating_params,
    check_plant_lifetime_params,
    check_supe_formulation_params,
)

# from idaes.core.solvers import get_solver


def define_custom_units():
    """
    This function defines custom units that are needed throughout the model.
    """

    # Define custom units.
    pyunits.load_definitions_from_strings(["year = [time]"])
    pyunits.load_definitions_from_strings(["EOL_Product = [item]"])
    pyunits.load_definitions_from_strings(["Operator = [item]"])
    pyunits.load_definitions_from_strings(["Disassembly_Unit = [item]"])
    pyunits.load_definitions_from_strings(["USD = [currency]"])
    pyunits.load_definitions_from_strings(["KUSD = 1000 * USD"])
    pyunits.load_definitions_from_strings(["MUSD = 1000 * KUSD"])


def build_model(
    ### Choice of objective function
    obj_func,
    ### Plant lifetime parameters
    plant_start,
    plant_lifetime,
    ### Feed parameters
    available_feed,
    collection_rate,
    tracked_comps,
    prod_comp_mass,
    ### Superstructure formulation parameters
    num_stages,
    options_in_stage,
    option_outlets,
    option_efficiencies,
    ### Operating parameters
    profit,
    opt_var_oc_params,
    operators_per_discrete_unit,
    yearly_cost_per_unit,
    capital_cost_per_unit,
    processing_rate,
    num_operators,
    labor_rate,
    ### Discretized costing parameters
    discretized_purchased_equipment_cost,
    ### Environmental impacts parameters
    consider_environmental_impacts,
    options_environmental_impacts,
    epsilon,
    ### Byproduct valorization parameters
    consider_byproduct_valorization,
    byproduct_values,
    byproduct_opt_conversions,
):
    """
    This function builds a superstructure model based on specifications from the user.

    Args:
        obj_func: (str) Choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.

        plant_start: (int) The year that plant construction begins.
        plant_lifetime: (int) The total lifetime of the plant, including plant construction. Must be at least three years.

        available_feed: (dict) Total feedstock available (number of EOL products) for recycling each year.
        collection_rate: (float) Collection rate for how much available feed is processed by the plant each year.
        tracked_comps: (list) List of tracked components.
        prod_comp_mass: (dict) Mass of tracked components per EOL product (kg / EOL product).

        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        option_efficiencies: (dict) Tracked component retention efficiency for each option.

        profit: (dict) Profit per unit of product in terms of tracked components ($/kg tracked component).
        opt_var_oc_params: (dict) Holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        operators_per_discrete_unit: (dict) Number of operators needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit ($/year) for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit ($) for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of end-of-life products disassembled per year per unit (number of EOL products / year).
        num_operators: (dict) Number of operators needed for each continuous option.
        labor_rate: (float) Yearly wage per operator ($ / year).

        discretized_equipment_cost: (dict) Discretized cost by flows entering for each continuous option ($/kg).



        consider_byproduct_valorization: (bool) Decide whether or not to consider the valorization of byproducts.
    """

    #################################################################################################
    ### Define custom units
    define_custom_units()

    ### Build model
    m = pyo.ConcreteModel()

    ### Objective function
    check_objective_function_choice(obj_func)

    ### Plant lifetime parameters
    # Check that plant lifetime parameters are feasible.
    check_plant_lifetime_params(plant_start, plant_lifetime)
    # Create separate block to hold plant lifetime parameters.
    add_plant_lifetime_params_block(m, plant_start, plant_lifetime)

    ### Feed parameters
    # Check that feed parameters are feasible.
    check_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # Create separate block to hold feed parameters.
    add_feed_params_block(
        m, available_feed, collection_rate, tracked_comps, prod_comp_mass
    )

    ### Superstructure formulation parameters
    # Check that superstructure formulation parameters are feasible.
    check_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )
    # Create separate block to hold superstructure formulation parameters.
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )

    ### Operating parameters
    # Check that operating parameters are feasible.
    check_operating_params(
        m,
        profit,
        opt_var_oc_params,
        operators_per_discrete_unit,
        yearly_cost_per_unit,
        capital_cost_per_unit,
        processing_rate,
        num_operators,
        labor_rate,
    )
    # Create a costing block, and add operating parameters to it.
    add_operating_params(
        m,
        profit,
        opt_var_oc_params,
        operators_per_discrete_unit,
        yearly_cost_per_unit,
        capital_cost_per_unit,
        processing_rate,
        num_operators,
        labor_rate,
    )

    ### Costing parameters
    # Check that costing parameters are feasible.
    check_discretized_costing_params(m, discretized_purchased_equipment_cost)
    # Create add parameters to costing block
    add_discretized_costing_params(m, discretized_purchased_equipment_cost)

    ### Mass balances
    # Generate mass balance parameters.
    add_mass_balance_params(m)
    # Generate mass balance variables.
    add_mass_balance_vars(m)
    # Generate mass balance constraints.
    add_mass_balance_cons(m)

    ### Costing
    # Generate costing parameters.
    add_costing_params(m)
    # Generate costing variables.
    add_costing_vars(m, obj_func)

    ### Environmental impacts
    check_environmental_impact_params(
        m, consider_environmental_impacts, options_environmental_impacts, epsilon
    )
    # only add params, vars, and cons if environmental impacts are considered.
    if consider_environmental_impacts:
        add_environmental_impact_params(
            m, consider_environmental_impacts, options_environmental_impacts, epsilon
        )
        add_environmental_impact_vars(m)
        add_environmental_impact_cons(m)

    ### Byproduct valorization
    check_byproduct_valorization_params(
        m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
    )
    # only add vars and cons if byproduct valorization is considered
    if consider_byproduct_valorization:
        add_byproduct_valorization_params(
            m, byproduct_values, byproduct_opt_conversions
        )
        add_byproduct_valorization_vars(m)
        add_byproduct_valorization_cons(m)

    # Generate costing constraints.
    add_profit_cons(m, obj_func, consider_byproduct_valorization)
    add_capital_cost_cons(m)
    add_operating_cost_cons(m)
    add_cash_flow_cons(m)
    add_costing_objective_functions(m, obj_func)

    ### Return model
    return m
