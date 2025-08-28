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

from idaes.core.scaling import CustomScalerBase

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
    add_feed_params,
    add_mass_balance_cons,
    add_mass_balance_params,
    add_mass_balance_vars,
    add_objective_function_choice_param,
    add_operating_cost_cons,
    add_operating_params,
    add_plant_lifetime_params,
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


class SuperstructureScaler(CustomScalerBase):

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale flowrate variables
        for c, v in model.fs.f.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.f_in.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.f_out.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.piecewise_flow_entering.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)

        # scale costing variables
        for c, v in model.fs.costing.net_present_value.items():
            self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)
        for c, v in model.fs.costing.main_product_profit.items():
            self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)
        for c, v in model.fs.costing.total_profit.items():
            self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)
        for c, v in model.fs.costing.piecewise_equipment_cost.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.equipment_cost.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.costing.total_plant_cost.items():
            self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)
        for c, v in model.fs.costing.financing.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.costing.other_costs.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.total_overnight_cost.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.opt_variable_operating_cost.items():
            self.set_variable_scaling_factor(v, 1e-3, overwrite=overwrite)
        for c, v in model.fs.costing.aggregate_variable_operating_cost.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.operators_per_option.items():
            self.set_variable_scaling_factor(v, 1e1, overwrite=overwrite)
        for c, v in model.fs.costing.cost_of_labor.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.m_and_sm.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.costing.sa_and_qa_qc.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.costing.s_ip_r_and_d.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for c, v in model.fs.costing.a_and_sl.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.fb.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.pt_and_i.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.aggregate_fixed_operating_cost.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.total_overnight_cost_expended.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.plant_overhead.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for c, v in model.fs.costing.total_operating_expense.items():
            self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)
        for c, v in model.fs.costing.cash_flow.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)

        # scale vars only if environmental impacts are considered
        if hasattr(model.fs, "environmental_impacts"):
            for c, v in model.fs.environmental_impacts.option_yearly_impacts.items():
                self.set_variable_scaling_factor(v, 1e-7, overwrite=overwrite)
            for c, v in model.fs.environmental_impacts.total_yearly_impacts.items():
                self.set_variable_scaling_factor(v, 1e-8, overwrite=overwrite)
            for c, v in model.fs.environmental_impacts.total_impacts.items():
                self.set_variable_scaling_factor(v, 1e-9, overwrite=overwrite)

        # scale vars only if byproduct valorization is being considered
        if hasattr(model.fs, "byproduct_valorization"):
            for c, v in model.fs.byproduct_valorization.byproduct_produced.items():
                self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
            for c, v in model.fs.byproduct_valorization.byproduct_produced.items():
                self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        if hasattr(model.fs.costing, "byproduct_profit"):
            for c, v in model.fs.costing.byproduct_profit.items():
                self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        if hasattr(model.fs.costing, "total_byproduct_profit"):
            for c, v in model.fs.costing.total_byproduct_profit.items():
                self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)

        # scale var only if NPV is chosen as the objective function
        if hasattr(model.fs.costing, "opt_profit"):
            for c, v in model.fs.costing.opt_profit.items():
                self.set_variable_scaling_factor(v, 1e-6, overwrite=overwrite)

        # scale var only if COR is chosen as the objective function
        if hasattr(model.fs.costing, "cost_of_recovery"):
            for c, v in model.fs.costing.cost_of_recovery.items():
                self.set_variable_scaling_factor(v, 1, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale flowrate constraints
        for j, c in model.fs.inlet_flow_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.init_flow_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.intermediate_flow_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.outlet_flow_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.init_flow_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.f_in_big_m_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.f_out_big_m_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.max_flow_entering_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )

        # Scale costing constraints
        for j, c in model.fs.costing.calculate_total_profit_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.discrete_opts_equipment_cost_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.add_units_to_piecewise_costs.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_total_plant_cost_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_financing_cost_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_other_costs_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_total_overnight_cost_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_opt_yearly_variable_expense_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for (
            j,
            c,
        ) in (
            model.fs.costing.calculate_total_yearly_variable_operating_costs_cons.items()
        ):
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_cost_of_labor_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_m_and_sm_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_sa_and_qa_qc_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_s_ip_r_and_d_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_a_and_sl_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_fb_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_pt_and_i_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for (
            j,
            c,
        ) in model.fs.costing.calculate_total_yearly_fixed_operating_costs_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for (
            j,
            c,
        ) in model.fs.costing.calculate_total_overnight_cost_expended_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_plant_overhead_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_total_operating_expense_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_cash_flows.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_net_present_value_con.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )
        for j, c in model.fs.costing.calculate_main_product_profit_cons.items():
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )

        # Scale constraint only if npv is objective function
        if hasattr(model.fs.costing, "calculate_final_opts_profit_cons"):
            for j, c in model.fs.costing.calculate_final_opts_profit_cons.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )

        # Scale constraint only if cor is objective function
        if hasattr(model.fs.costing, "set_net_present_value_to_zero_con"):
            for j, c in model.fs.costing.set_net_present_value_to_zero_con.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )

        # Scale constraints only if byproduct valorization is considered
        if hasattr(model.fs, "byproduct_valorization"):
            for (
                j,
                c,
            ) in (
                model.fs.byproduct_valorization.calculate_byproduct_produced_cons.items()
            ):
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
        if hasattr(model.fs.costing, "calculate_byproduct_profit_cons"):
            for j, c in model.fs.costing.calculate_byproduct_profit_cons.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
        if hasattr(model.fs.costing, "calculate_opt_byprod_val_cons"):
            for j, c in model.fs.costing.calculate_opt_byprod_val_cons.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )

        # Scale constraints only if environmental impacts are considered
        if hasattr(model.fs, "environmental_impacts"):
            for (
                j,
                c,
            ) in (
                model.fs.environmental_impacts.calculate_opt_yearly_impacts_con.items()
            ):
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
            for (
                j,
                c,
            ) in model.fs.environmental_impacts.calculate_yearly_impacts_con.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
            for (
                j,
                c,
            ) in model.fs.environmental_impacts.calculate_total_impacts_con.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )
            for j, c in model.fs.environmental_impacts.epsilon_con.items():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )


def define_custom_units():
    """
    This function defines custom units that are needed throughout the model.
    """
    # Define custom units.
    pyunits.load_definitions_from_strings(["EOL_Product = [item]"])
    pyunits.load_definitions_from_strings(["Operator = [item]"])
    pyunits.load_definitions_from_strings(["Disassembly_Unit = [item]"])
    pyunits.load_definitions_from_strings(["USD = [currency]"])


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

        # Plant lifetime parameters
        plant_start: (int) The year that plant construction begins.
        plant_lifetime: (int) The total lifetime of the plant, including plant construction. Must be at least three years.

        # Feed parameters
        available_feed: (dict) Total feedstock available (number of EOL products) for recycling each year.
        collection_rate: (float) Collection rate for how much available feed is processed by the plant each year.
        tracked_comps: (list) List of tracked components.
        prod_comp_mass: (dict) Mass of tracked components per EOL product (kg / EOL product).

        # Superstructure formulation parameters
        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        option_efficiencies: (dict) Tracked component retention efficiency for each option.

        # Operating parameters
        profit: (dict) Profit per unit of product in terms of tracked components ($/kg tracked component).
        opt_var_oc_params: (dict) Holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        operators_per_discrete_unit: (dict) Number of operators needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit ($/year) for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit ($) for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of end-of-life products disassembled per year per unit (number of EOL products / year).
        num_operators: (dict) Number of operators needed for each continuous option.
        labor_rate: (float) Yearly wage per operator ($ / year).

        # Discretized costing parameters
        discretized_equipment_cost: (dict) Discretized cost by flows entering for each continuous option ($/kg).

        # Environmental impacts parameters
        consider_environmental_impacts: (bool) Choice of whether or not to consider environmental impacts.
        options_environmental_impacts: (dict) Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.
        epsilon: (float) Epsilon factor for generating the Pareto front.

        # Byproduct valorization parameters
        consider_byproduct_valorization: (bool) Decide whether or not to consider the valorization of byproducts.
        byproduct_values: (dict) Byproducts considered, and their value ($/kg).
        byproduct_opt_conversions: (dict) Conversion factors for each byproduct for each option.
    """

    #################################################################################################
    ### Define custom units
    define_custom_units()

    ### Build model
    m = pyo.ConcreteModel()

    ### Objective function
    # Check that choice of objective function is feasible
    check_objective_function_choice(obj_func)
    # Add choice of objective function as parameter to model
    add_objective_function_choice_param(m, obj_func)

    ### Plant lifetime parameters
    # Check that plant lifetime parameters are feasible.
    check_plant_lifetime_params(plant_start, plant_lifetime)
    # Create separate block to hold plant lifetime parameters.
    add_plant_lifetime_params(m, plant_start, plant_lifetime)

    ### Feed parameters
    # Check that feed parameters are feasible.
    check_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # Create separate block to hold feed parameters.
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)

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

    ### Environmental impacts
    check_environmental_impact_params(
        m, consider_environmental_impacts, options_environmental_impacts, epsilon
    )
    # only add params, vars, and cons if environmental impacts are considered.
    if consider_environmental_impacts:
        add_environmental_impact_params(m, options_environmental_impacts, epsilon)
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
    add_profit_cons(m, consider_byproduct_valorization)
    add_capital_cost_cons(m)
    add_operating_cost_cons(m)
    add_cash_flow_cons(m)
    add_costing_objective_functions(m, obj_func)

    ### Return model
    return m
