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

from idaes.core.solvers import get_solver

from prommis.superstructure.add_superstructure_blocks import (
    add_byproduct_valorization_cons,
    add_byproduct_valorization_params,
    add_byproduct_valorization_vars,
    add_costing_cons,
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
    add_operating_params,
    add_plant_lifetime_params_block,
    add_supe_formulation_params,
)
from prommis.superstructure.check_superstructure_inputs import (
    check_byproduct_valorization_params,
    check_discretized_costing_params,
    check_environmental_impact_params,
    check_feed_params,
    check_operating_params,
    check_plant_lifetime_params,
    check_supe_formulation_params,
)


def configure_model(m, obj_func):
    """
    The configures the model based on the specifications of the user by activating and deactivating different
    blocks to ensure the correct constraints are considered.

    Args:
        m: pyomo model.
        obj_func: (str) Choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive.
    """
    ### Check types and structure.
    ## Check that obj_fun is of type str.
    if not isinstance(obj_func, str):
        raise TypeError("obj_func is not of type str.")

    ### Run tests
    if (obj_func != "NPV") and (obj_func != "COR"):
        raise ValueError(
            "Invalid choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive."
        )

    # Check the objective function.
    if obj_func == "NPV":
        # deactivate the cost of recovery constraints and objective function if the NPV objective function is chosen.
        m.fs.costing.cost_of_recovery.deactivate()
    else:
        # deactivate the net present value constraints and objective function if the cost of recovery objective function is chosen.
        m.fs.costing.net_present_value.deactivate()

    # Check if environmental impacts are considered.
    if pyo.value(m.fs.environmental_impacts.consider_environmental_impacts) == False:
        # deactivate environmental impact constraints if they are not considered.
        m.fs.environmental_impacts.deactivate()

    # Check if byproduct valorization is considered.  Different constraints are considered depending on if they are considered or not.
    if (
        pyo.value(m.fs.byproduct_valorization.consider_byproduct_valorization)
        == True
    ):
        # deactivate the constraints that're associated with no byproduct valorization if it is considered.
        m.fs.no_byproduct_valorization.deactivate()
    else:
        # deactivate the constraints that're associated with byproduct valorization if it is not considered.
        m.fs.byproduct_valorization.deactivate()


def build_model(
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
    ### Choice of objective function
    obj_func,
):

    #################################################################################################
    ### Build model
    m = pyo.ConcreteModel()

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
    add_costing_vars(m)
    # Generate costing constraints.
    add_costing_cons(m)

    ### Environmental impacts
    check_environmental_impact_params(
        m, consider_environmental_impacts, options_environmental_impacts, epsilon
    )
    add_environmental_impact_params(
        m, consider_environmental_impacts, options_environmental_impacts, epsilon
    )
    add_environmental_impact_vars(m)
    add_environmental_impact_cons(m)

    ### Byproduct valorization
    check_byproduct_valorization_params(
        m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
    )
    add_byproduct_valorization_params(
        m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
    )
    add_byproduct_valorization_vars(m)
    add_byproduct_valorization_cons(m)

    ### Configure model
    configure_model(m, obj_func)

    ### Return model
    return m
