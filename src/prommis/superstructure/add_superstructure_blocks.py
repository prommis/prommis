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

import copy
import math

import pyomo.environ as pyo


def add_plant_lifetime_params_block(m, plant_start, plant_lifetime):
    """
    This function builds the rest of the plant lifetime parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model.
        plant_start: (int) The year that plant construction begins.
        plant_lifetime: (int) The total lifetime of the plant, including plant construction. Must be at least three years.
    """
    ### Define parameters from user input.
    # Define the start of plant production. Assume it starts one year after plant construction.
    prod_start = plant_start + 1
    # Define the final year of the plant's lifetime.
    plant_end = plant_start + plant_lifetime - 1

    ### Define necessary pyomo parameters.
    ## Create block
    m.fs = pyo.Block(doc="Main flowsheet.")
    ## Pyomo parameters
    m.fs.plant_start = pyo.Param(
        initialize=plant_start, doc="The year that plant construction begins."
    )
    m.fs.plant_lifetime = pyo.Param(
        initialize=plant_lifetime,
        doc="The total lifetime of the plant, including plant construction. Must be at least three years.",
    )
    m.fs.prod_start = pyo.Param(
        initialize=prod_start, doc="The first year of plant production."
    )
    m.fs.plant_end = pyo.Param(
        initialize=plant_end, doc="The final year of plant production."
    )
    m.fs.plant_life_range = pyo.RangeSet(
        plant_start, plant_end, doc="Lifetime of the plant."
    )
    m.fs.operational_range = pyo.RangeSet(
        prod_start, plant_end, doc="Operational lifetime of the plant."
    )


def add_feed_params_block(
    m, available_feed, collection_rate, tracked_comps, prod_comp_mass
):
    """
    This function builds the rest of the feed parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model.
        available_feed: (dict) Total feedstock available for recycling each year.
        collection_rate: (float) Collection rate for how much available feed is processed by the plant each year.
        tracked_comps: (list) List of tracked components.
        prod_comp_mass: (dict) Mass of tracked components per EOL product.
    """
    ### Define parameters from user input.
    # Define feed entering each year of plant operation in terms of available feed and collection rate.
    feed_entering = copy.deepcopy(available_feed)
    for key in feed_entering:
        feed_entering[key] = available_feed[key] * collection_rate
    # Define max feed entering plant over production period.
    max_feed_entering = max(feed_entering.values())
    # Define year in which max feed enters plant.
    max_feed_entering_year = max(feed_entering, key=feed_entering.get)

    ### Define necessary pyomo parameters.
    ## Pyomo parameters
    m.fs.available_feed = pyo.Param(
        m.fs.operational_range,
        initialize=available_feed,
        doc="The total feedstock available for recycling each year.",
    )
    m.fs.collection_rate = pyo.Param(
        initialize=collection_rate,
        doc="The fraction of available feed that is processed by the plant each year.",
    )
    m.fs.tracked_comps = pyo.Set(
        initialize=tracked_comps, doc="Tracked components."
    )
    m.fs.prod_comp_mass = pyo.Param(
        m.fs.tracked_comps,
        initialize=prod_comp_mass,
        doc="The mass of each tracked component per EOL product.",
    )
    m.fs.feed_entering = pyo.Param(
        m.fs.operational_range,
        initialize=feed_entering,
        doc="The amount of feed entering the plant each year.",
    )
    m.fs.max_feed_entering = pyo.Param(
        initialize=max_feed_entering,
        doc="The max yearly feed that enters the plant over the production period.",
    )
    m.fs.max_feed_entering_year = pyo.Param(
        initialize=max_feed_entering_year,
        doc="The year that the max feed enters the plant.",
    )


def add_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, option_efficiencies
):
    """
    This function builds the rest of the superstructure formulation parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model.
        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        option_efficiencies: (dict) Tracked component retention efficiency for each option.
    """
    ### Define parameters from user input.
    # Define a parameter for the max number of options in any of the stages.
    max_options = max(options_in_stage.values())
    # Define a set of all the discrete options. It is assumed that all options in the first stage are discrete.
    discrete_opts_set = set((1, k) for k in range(1, options_in_stage[1] + 1))
    # Define a set containing all the options in the superstructure.
    all_opts_set = set(
        (j, k)
        for j in range(1, num_stages + 1)
        for k in range(1, options_in_stage[j] + 1)
    )
    # Define a set containing all the continuous options in the superstructure. All options after the first stage are continuous.
    continuous_opts_set = all_opts_set - discrete_opts_set
    # Define a set containing all the  options in the final stage.
    final_opts_list = [
        (num_stages, k) for k in range(1, options_in_stage[num_stages] + 1)
    ]

    ### Define functions needed to initialize pyomo parameters.
    # Define a function for initializing option efficiency pyomo parameter.
    def option_efficiencies_initialize(m, j, k, c):
        return option_efficiencies[(j, k)][c]

    ### Define necessary pyomo parameters.

    ## Pyomo parameters
    m.fs.num_stages = pyo.Param(
        initialize=num_stages, doc="The total number of stages in the superstructure."
    )
    m.fs.stages_set = pyo.RangeSet(
        1, num_stages, doc="A set of all the stages in the superstructure."
    )
    m.fs.options_in_stage = pyo.Param(
        m.fs.stages_set,
        initialize=options_in_stage,
        doc="The number of options in each stage.",
    )
    m.fs.all_opts_set = pyo.Set(
        initialize=(
            (j, k)
            for j in m.fs.stages_set
            for k in range(1, options_in_stage[j] + 1)
        ),
        doc="Set containing all options in the superstructure.",
    )
    m.fs.discrete_opts_set = pyo.Set(
        initialize=((opt) for opt in discrete_opts_set),
        doc="Set containing all the options which utilize discrete units (discrete options). These are the options in the first stage, "
        "which utilize discrete units (or operators) to disassemble the incoming end-of-life products.",
    )
    m.fs.continuous_opts_set = pyo.Set(
        initialize=((opt) for opt in continuous_opts_set),
        doc="Set containing all the options in the stage which don't utilize discrete units (continuous options). "
        "All options after the first stage are assumed to be continuous. The sizing of these options is calculated "
        "based on the sum of the tracked components entering the option.",
    )
    m.fs.option_outlets = pyo.Param(
        m.fs.all_opts_set,
        initialize=option_outlets,
        doc="Defines the set of options k' in stage j+1 connected to option k in stage j.",
    )
    m.fs.option_efficiencies = pyo.Param(
        m.fs.all_opts_set,
        m.fs.tracked_comps,
        initialize=option_efficiencies_initialize,
        doc="Defines the tracked component efficiencies for each option in the superstructure. "
        "Efficiency defined for each tracked component 'c' for each option 'k' in each stage 'j'",
    )
    m.fs.final_opts_set = pyo.Set(
        initialize=final_opts_list,
        doc="Set containing all of the options in the final stage.",
    )


def add_operating_params(
    m,
    profit,
    opt_var_oc_params,
    operators_per_discrete_unit,
    yearly_cost_per_unit,
    capital_cost_per_unit,
    processing_rate,
    num_operators,
    labor_rate,
):
    """
    This function builds the rest of the operating parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model.
        profit: (dict) Profit per unit of product in terms of tracked components.
        opt_var_oc_params: (dict) Holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        operators_per_discrete_unit: (dict) Number of operators needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of end-of-life products disassembled per year per unit.
        num_operators: (dict) Number of operators needed for each continuous option.
        labor_rate: (float) Yearly wage per operator.
    """
    ### Define parameters from user input.
    ## Calculate the number of discrete units needed for each option. There must be enough discrete units to handle the maximum flow of end-of-life products
    # that enters the plant over its operational lifetime.
    discrete_units_per_option = copy.deepcopy(processing_rate)
    for key in discrete_units_per_option.keys():
        discrete_units_per_option[key] = math.ceil(
            m.fs.max_feed_entering / processing_rate[key]
        )

    ### Define functions needed to initialize pyomo parameters.
    # Define a function for initializing profit pyomo parameter.
    def profit_initialization(m, j, k, c):
        return profit[(j, k)][c]

    # Define a function for initializing opt_var_oc_params pyomo parameter.
    def opt_var_oc_params_initialization(m, j, k, var_oc_param):
        return opt_var_oc_params[(j, k)][var_oc_param]

    ### Define necessary pyomo parameters.
    ## Create block
    m.fs.costing = pyo.Block(
        doc="Block to hold costing parameters, variables, and constraints."
    )
    ## Pyomo parameters
    m.fs.costing.profit = pyo.Param(
        m.fs.final_opts_set,
        m.fs.tracked_comps,
        initialize=profit_initialization,
        doc="Holds the profit for all options in the final stage in terms of the tracked components.",
    )
    m.fs.costing.var_oc_params_set = pyo.Set(
        initialize=["a", "b"],
        doc="Set containing the necessary parameters for calculating the variable operating costs for all continuous options.",
    )
    m.fs.costing.opt_var_oc_params = pyo.Param(
        m.fs.continuous_opts_set,
        m.fs.costing.var_oc_params_set,
        initialize=opt_var_oc_params_initialization,
        doc="Holds all the variable operating costs parameter values for all continuous options.",
    )
    m.fs.costing.operators_per_discrete_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=operators_per_discrete_unit,
        doc="The number of operators needed per discrete unit for discrete options.",
    )
    m.fs.costing.yearly_cost_per_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=yearly_cost_per_unit,
        doc="The operating costs per discrete unit for discrete options.",
    )
    m.fs.costing.capital_cost_per_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=capital_cost_per_unit,
        doc="The capital cost per discrete unit for discrete options.",
    )
    m.fs.costing.processing_rate = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=processing_rate,
        doc="The processing rate per discrete unit for discrete options. In terms of number of end-of-life products "
        "disassembled per discrete unit per year.",
    )
    m.fs.costing.num_operators = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=num_operators,
        doc="The number of operators per continuous options.",
    )
    m.fs.costing.labor_rate = pyo.Param(
        initialize=labor_rate, doc="The yearly wage per operator."
    )
    m.fs.costing.discrete_units_per_option = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=discrete_units_per_option,
        doc="The number of discrete units per option needed to disassemble all the incoming end-of-life products over the operational "
        "lifetime of the plant.",
    )


def add_discretized_costing_params(m, discretized_purchased_equipment_cost):
    """
    This function adds all the discretized costing parameters to a block.

    Args:
        m: pyomo model.
        discretized_purchased_equipment_cost: (dict) Discretized cost by flows entering for each continuous option
    """
    ### Define parameters from user input.
    # Define a dict to hold discretized flowrate data for each option.
    flowrates_data = {}
    for opt in m.fs.continuous_opts_set:
        flowrates_data[opt] = discretized_purchased_equipment_cost[opt]["Flowrates"]
    # Define a dict to hold discretized cost data for each option.
    costs_data = {}
    for opt in m.fs.continuous_opts_set:
        costs_data[opt] = discretized_purchased_equipment_cost[opt]["Costs"]

    ### Define necessary pyomo parameters.
    ## Pyomo parameters
    m.fs.costing.flowrates_data = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=flowrates_data,
        doc="Discretized flowrate data for all continuous options.",
    )
    m.fs.costing.costs_data = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=costs_data,
        doc="Discretized costing data for all continuous options.",
    )


def add_mass_balance_params(m):
    """
    This function builds the mass balance parameters.

    Args:
        m: pyomo model.
    """
    ## Define parameters
    # Calculate Big-M values for each tracked component.
    m_val = {}
    for c in m.fs.tracked_comps:
        # Max value assumes 100% efficiency (no losses) over the plant's lifetime. Rounded up to ensure validity
        # in constraints.
        m_val[c] = math.ceil(
            m.fs.max_feed_entering * m.fs.prod_comp_mass[c]
        )

    # Create an upper bound for the total flow entering a continuous option.
    max_flow_upper_bound = sum(
        m_val[c] for c in m.fs.tracked_comps
    )

    ## Pyomo parameters
    # Flow entering each stage (except the last stage).
    m.fs.f_stages = pyo.RangeSet(
        1,
        m.fs.num_stages - 1,
        doc="Set of all stages except the last. Used to define the flow variable: 'f'.",
    )
    m.fs.flow_set = pyo.Set(
        initialize=m.fs.all_opts_set
        * m.fs.tracked_comps
        * m.fs.operational_range,
        doc="Set of all options, tracked components, and operational years of the plant. Used to define the flow variables: 'f_in' and 'f_out'.",
    )
    m.fs.big_m_val = pyo.Param(
        m.fs.tracked_comps,
        initialize=m_val,
        doc="Big-M parameters used in Equations (7) and (8) from the documentation.",
    )
    m.fs.max_flow_upper_bound = pyo.Param(
        initialize=max_flow_upper_bound,
        doc="Upper bound for the total flow entering a continuous option.",
    )


def add_mass_balance_vars(m):
    """
    This function builds the mass balance variables.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.fs.f = pyo.Var(
        m.fs.f_stages
        * m.fs.tracked_comps
        * m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each stage (except the last stage). See documentation for more details.",
    )
    m.fs.f_in = pyo.Var(
        m.fs.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )
    m.fs.f_out = pyo.Var(
        m.fs.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )
    m.fs.option_binary_var = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.Binary,
        doc="Binary variables to indicate whether or not an option has been selected.",
    )
    m.fs.flow_entering = pyo.Var(
        m.fs.continuous_opts_set,
        bounds=(0, m.fs.max_flow_upper_bound),
        doc="The max total flow that enters each continuous option over the lifetime"
        "of the plant.",
    )


def add_mass_balance_cons(m):
    """
    This function builds the mass balance constraints.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    @m.fs.Constraint(
        m.fs.stages_set,
        m.fs.tracked_comps,
        m.fs.operational_range,
        doc="Equation (1) from the documentation. This constraint ensures that the sum of the inlet streams is equal to the upstream flow from the previous stage.",
    )
    def inlet_flow_cons(b, j, c, t):
        if j == 1:
            return pyo.Constraint.Skip
        else:
            # Extract all the options available in stage 'j'
            num_options = range(1, b.options_in_stage[j] + 1)
            return b.f[j - 1, c, t] == sum(
                b.f_in[j, k, c, t] for k in num_options
            )

    @m.fs.Constraint(
        m.fs.stages_set,
        m.fs.tracked_comps,
        m.fs.operational_range,
        doc="Equation (2) from the documentation. The flow of each component entering the plant each year is equal to the number of EOL products the plant processes multiplied by the amount of the component per EOL product.",
    )
    def init_flow_cons(b, j, c, t):
        if j == 1:
            # Extract all the options available in stage 'j'
            num_options = range(1, b.options_in_stage[j] + 1)
            return b.feed_entering[t] * b.prod_comp_mass[
                c
            ] == sum(b.f_in[j, k, c, t] for k in num_options)
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.flow_set, doc="Equation (3) from the documentation. The steam exiting an option is related to the inlet stream by an efficiency parameter."
    )
    def intermediate_flow_cons(b, j, k, c, t):
        alpha = b.option_efficiencies[j, k, c]
        return (
            b.f_in[j, k, c, t] * alpha
            == b.f_out[j, k, c, t]
        )

    @m.fs.Constraint(
        m.fs.stages_set,
        m.fs.tracked_comps,
        m.fs.operational_range,
        doc="Equation (4) from the documentation. The flow exiting a stage is equal to the sume of the outlet streams.",
    )
    def outlet_flow_cons(b, j, c, t):
        if j != b.num_stages:
            # Extract all the options available in stage 'j'
            num_options = range(1, b.options_in_stage[j] + 1)
            return (
                sum(b.f_out[j, k, c, t] for k in num_options)
                == b.f[j, c, t]
            )
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.stages_set, doc="Equation (5) from the documentation. Ensures that only one option per stage is chosen."
    )
    def stage_binary_cons(b, j):
        # Extract all the options available in stage 'j'
        num_options = range(1, b.options_in_stage[j] + 1)
        return sum(b.option_binary_var[j, k] for k in num_options) == 1

    @m.fs.Constraint(
        m.fs.all_opts_set, doc="Equation (6) from the documentation. Esures that if an option in stage $j + 1$ can only be chosen if it is connected to an option that was chosen in stage $j$."
    )
    def connection_binary_cons(b, j, k):
        if j != m.fs.num_stages:
            # Extract the set of options k' in stage j+1 connected to option k in stage j.
            opt_connections = b.option_outlets[j, k]
            return b.option_binary_var[j, k] <= sum(
                b.option_binary_var[j + 1, kp] for kp in opt_connections
            )
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.flow_set, doc="Equation (7) from the documentation. Big-M constraint for inlet streams. Inlet streams must be zero if the option is not chosen."
    )
    def f_in_big_m_cons(b, j, k, c, t):
        return (
            b.f_in[j, k, c, t]
            <= b.option_binary_var[j, k] * b.big_m_val[c]
        )

    @m.fs.Constraint(
        m.fs.flow_set, doc="Equation (8) from the documentation. Big-M constraint for outlet streams. Outlet streams must be zero if the option is not chosen."
    )
    def f_out_big_m_cons(b, j, k, c, t):
        return (
            b.f_out[j, k, c, t]
            <= b.option_binary_var[j, k] * b.big_m_val[c]
        )
    
    @m.fs.costing.Constraint(
        m.fs.continuous_opts_set,
        m.fs.operational_range,
        doc="Constraint to determine the max flow entering each continuous option over the lifetime of the plant.",
    )
    def max_flow_entering_cons(b, j, k, t):
        return m.fs.flow_entering[j, k] >= sum(
            m.fs.f_in[j, k, c, t] for c in m.fs.tracked_comps
        )


def add_costing_params(m):
    """
    This function builds the costing parameters. 
    
    To calculate the capital cost parameters, the methodology presented in NETL's 
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate 
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed.

    **References:**
    
    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.
    
    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    Args:
        m: pyomo model.
    """
    ## Define parameters
    ## Create blocks
    m.fs.costing.net_present_value = pyo.Block(doc="Block for holding net present value.")
    m.fs.costing.cost_of_recovery = pyo.Block(doc="Block for holding cost of recovery.")

    ## Pyomo parameters
    m.fs.costing.lang_factor = pyo.Param(initialize=2.97)
    m.fs.costing.i_operating_expense_escalation = pyo.Param(
        initialize=0.03, doc="Operating expenses escalation rate."
    )
    m.fs.costing.i_capital_escalation = pyo.Param(
        initialize=0.036, doc="Capital expenses escalation rate."
    )
    m.fs.costing.discount_factor = pyo.Param(
        initialize=0.0577, doc="Discount factor for calculate the net present value."
    )


def add_costing_vars(m):
    """
    This function builds the costing variables.

    To calculate the capital cost parameters, the methodology presented in NETL's 
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate 
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed.

    **References:**
    
    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.
    
    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.fs.costing.opt_profit = pyo.Var(
        m.fs.final_opts_set,
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The profit generated by each option in the final processing stage each year.",
    )
    m.fs.costing.total_profit = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total profit generated by the plant each year.",
    )
    m.fs.costing.total_byproduct_profit = pyo.Var(
        m.fs.operational_range,
        domain=pyo.Reals,
        doc="The total profit generated by the plant each year from the valorization of byproducts.",
    )
    m.fs.costing.npv = pyo.Var(domain=pyo.Reals, doc="The net present value.")
    m.fs.costing.cor = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The cost of recovery."
    )
    m.fs.costing.purchased_equipment_cost = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The cost of purchased equipment.",
    )
    m.fs.costing.total_plant_cost = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The total plant cost."
    )
    m.fs.costing.financing = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The total financing cost of the plant."
    )
    m.fs.costing.other_costs = pyo.Var(
        domain=pyo.NonNegativeReals, 
        doc="'Other costs' associated with the plant as defined by QGESS[1]."
    )
    m.fs.costing.total_overnight_cost = pyo.Var(
        domain=pyo.NonNegativeReals, 
        doc="""The total overnight cost of the plant."""
    )
    m.fs.costing.opt_var_operating_expense = pyo.Var(
        m.fs.all_opts_set * m.fs.operational_range,
        domain=pyo.Reals,
        doc="Yearly variable operating expense for each option.",
    )
    m.fs.costing.total_var_operating_expense = pyo.Var(
        m.fs.operational_range,
        domain=pyo.Reals,
        doc="Total yearly variable operating expense.",
    )
    m.fs.costing.operators_per_option = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The number of operators needed for each option.",
    )
    m.fs.costing.total_operators = pyo.Var(
        domain=pyo.NonNegativeIntegers,
        doc="The total number of operators needed for the process. Must be an integer value.",
    )
    m.fs.costing.cost_of_labor = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The cost of labor for the process."
    )
    m.fs.costing.m_and_sm = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Maintenance & Supply Materials (M&SM)."
    )
    m.fs.costing.sa_and_qa_qc = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Sample Analysis & Quality Assurance/Quality Control (SA&QA/QC).",
    )
    m.fs.costing.s_ip_r_and_d = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Sales, Intellectual Property, and Research & Development (S,IP,R&D).",
    )
    m.fs.costing.a_and_sl = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Administrative & Supporting Labor (A&SL)."
    )
    m.fs.costing.fb = pyo.Var(domain=pyo.NonNegativeReals, doc="Fringe Benefits (FB).")
    m.fs.costing.pt_and_i = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Property Taxes & Insurance (PT&I)."
    )
    m.fs.costing.total_fixed_operating_expense = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Total yearly fixed operating expense.",
    )
    m.fs.costing.total_overnight_cost_expended = pyo.Var(
        m.fs.plant_life_range,
        domain=pyo.NonNegativeReals,
        doc="The total overnight cost expended each year.",
    )
    m.fs.costing.plant_overhead = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The yearly plant overhead.",
    )
    m.fs.costing.total_operating_expense = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total operating expense each year.",
    )
    m.fs.costing.cash_flow = pyo.Var(
        m.fs.plant_life_range,
        domain=pyo.Reals,
        doc="The yearly cash flow.",
    )


def add_costing_cons(m):
    """
    This function builds the costing constraints.

    To calculate the capital cost parameters, the methodology presented in NETL's 
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate 
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed.

    **References:**
    
    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.
    
    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    @m.fs.costing.net_present_value.Constraint(
        m.fs.final_opts_set,
        m.fs.operational_range,
        doc="Calculates the profit from each option each year in the final stage when the net present value is selected as "
        "the objective function.",
    )
    def calculate_final_opts_profit_cons(b, j, k, t):
        return m.fs.costing.opt_profit[j, k, t] == sum(
            m.fs.f_out[j, k, c, t] * m.fs.costing.profit[j, k, c]
            for c in m.fs.tracked_comps
        )

    @m.fs.costing.net_present_value.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly profit when the net present value is selected as the objective function.",
    )
    def calculate_total_profit_cons(b, t):
        return (
            m.fs.costing.total_profit[t]
            == sum(
                m.fs.costing.opt_profit[opt, t]
                for opt in m.fs.final_opts_set
            )
            + m.fs.costing.total_byproduct_profit[t]
        )

    @m.fs.costing.cost_of_recovery.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly profit when the cost of recovery is selected as the objective function.",
    )
    def calculate_total_profit_cons(b, t):
        return (
            m.fs.costing.total_profit[t]
            == m.fs.costing.cor
            * sum(
                sum(
                    m.fs.f_out[opt, c, t]
                    for c in m.fs.tracked_comps
                )
                for opt in m.fs.final_opts_set
            )
            + m.fs.costing.total_byproduct_profit[t]
        )

    @m.fs.costing.Constraint(
        m.fs.discrete_opts_set,
        doc="Calculates the purchased cost of equipment for all discrete options. "
        "Done by multiplying the number of discrete units by the capital cost per unit.",
    )
    def discrete_opts_purchased_equipment_cost_cons(b, j, k):
        return (
            m.fs.costing.purchased_equipment_cost[j, k]
            == m.fs.costing.discrete_units_per_option[j, k]
            * m.fs.costing.capital_cost_per_unit[j, k]
        )

    def piecewise_rule(b, j, k):
        flow_data = m.fs.costing.flowrates_data[j, k]
        purchased_equipment_cost_data = m.fs.costing.costs_data[j, k]

        # use m.add_component to generate all piecewise functions
        # piecewise = Piecewise(yval, xval, *kwargs)
        piecewise = pyo.Piecewise(
            m.fs.costing.purchased_equipment_cost[j, k],
            m.fs.flow_entering[j, k],
            pw_pts=flow_data,
            pw_constr_type="EQ",
            f_rule=purchased_equipment_cost_data,
            pw_repn="SOS2",
        )
        b.add_component("Option_" + str((j, k)) + "_Piecewise_Constraint", piecewise)

    m.fs.costing.piecewise_cons = pyo.Block(
        m.fs.continuous_opts_set,
        rule=piecewise_rule,
        doc="This block holds all the piecewise constraints for calculating the "
        "purchased equipment costs for all the continuous option.",
    )

    @m.fs.costing.Constraint(
        doc="Calculates the total plant cost of the plant. See Equation (x) in the documentation."
    )
    def calculate_total_plant_cost_con(b):
        return m.fs.costing.total_plant_cost == sum(
            m.fs.costing.purchased_equipment_cost[opt]
            for opt in m.fs.discrete_opts_set
        ) + m.fs.costing.lang_factor * sum(
            m.fs.costing.purchased_equipment_cost[opt]
            for opt in m.fs.continuous_opts_set
        )

    @m.fs.costing.Constraint(
        doc="Calculates the total financing cost of the plant. Assumed to be 2.7% of the total plant cost."
    )
    def calculate_financing_cost_con(b):
        return m.fs.costing.financing == 0.027 * m.fs.costing.total_plant_cost

    @m.fs.costing.Constraint(
        doc="Calculates the 'other costs' of the plant. Assumed to be 15% of the total plant cost."
    )
    def calculate_other_costs_con(b):
        return m.fs.costing.other_costs == 0.15 * m.fs.costing.total_plant_cost

    @m.fs.costing.Constraint(
        doc="Calculates the total overnight cost of the plant. "
        "Equal to the total plant cost + financing + 'other costs'."
    )
    def calculate_total_overnight_cost_con(b):
        return (
            m.fs.costing.total_overnight_cost
            == m.fs.costing.total_plant_cost + m.fs.costing.financing + m.fs.costing.other_costs
        )

    @m.fs.costing.Constraint(
        m.fs.all_opts_set,
        m.fs.operational_range,
        doc="Calculates yearly operating expense for each option.",
    )
    def calculate_opt_yearly_variable_expense_cons(b, j, k, t):
        if (j, k) in m.fs.discrete_opts_set:
            return (
                m.fs.costing.opt_var_operating_expense[j, k, t]
                == m.fs.costing.discrete_units_per_option[j, k]
                * m.fs.costing.yearly_cost_per_unit[j, k]
                * m.fs.option_binary_var[j, k]
            )
        else:
            return (
                m.fs.costing.opt_var_operating_expense[j, k, t]
                == m.fs.costing.opt_var_oc_params[j, k, "a"]
                * sum(
                    m.fs.f_in[j, k, c, t]
                    for c in m.fs.tracked_comps
                )
                + m.fs.costing.opt_var_oc_params[j, k, "b"]
                * m.fs.option_binary_var[j, k]
            )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly variable operating cost.",
    )
    def calculate_total_yearly_variable_operating_costs_cons(b, t):
        return m.fs.costing.total_var_operating_expense[t] == sum(
            m.fs.costing.opt_var_operating_expense[opt, t]
            for opt in m.fs.all_opts_set
        )

    @m.fs.costing.Constraint(
        m.fs.all_opts_set,
        doc="Calculate the number of operators needed for each option.",
    )
    def calculate_operators_per_option_cons(b, j, k):
        if (j, k) in m.fs.discrete_opts_set:
            return (
                m.fs.costing.operators_per_option[j, k]
                == m.fs.costing.discrete_units_per_option[j, k]
                * m.fs.costing.operators_per_discrete_unit[j, k]
                * m.fs.option_binary_var[j, k]
            )
        else:
            return (
                m.fs.costing.operators_per_option[j, k]
                == m.fs.costing.num_operators[j, k]
                * m.fs.option_binary_var[j, k]
            )

    @m.fs.costing.Constraint(
        doc="Calculate the number of operators needed for the entire process."
    )
    def calculate_total_operators_cons(b):
        return m.fs.costing.total_operators >= sum(
            m.fs.costing.operators_per_option[opt] for opt in m.fs.all_opts_set
        )

    @m.fs.costing.Constraint(doc="Calculates the cost of labor for the process.")
    def calculate_cost_of_labor_con(b):
        return (
            m.fs.costing.cost_of_labor == m.fs.costing.total_operators * m.fs.costing.labor_rate
        )

    @m.fs.costing.Constraint(
        doc="Calculate M&SM. Assumed to be 2% of the total plant cost."
    )
    def calculate_m_and_sm_con(b):
        return m.fs.costing.m_and_sm == 0.02 * m.fs.costing.total_plant_cost

    @m.fs.costing.Constraint(
        doc="Calculate SA&QA/QC. Assumed to be 10% of the cost of labor."
    )
    def calculate_sa_and_qa_qc_con(b):
        return m.fs.costing.sa_and_qa_qc == 0.1 * m.fs.costing.cost_of_labor

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculate S,IP,R&D. Assumed to be 1% of the revenue each year.",
    )
    def calculate_s_ip_r_and_d_con(b, t):
        return m.fs.costing.s_ip_r_and_d[t] == 0.01 * m.fs.costing.total_profit[t]

    @m.fs.costing.Constraint(doc="Calculate A&SL. Assumed to be 20% of the cost of labor.")
    def calculate_a_and_sl_con(b):
        return m.fs.costing.a_and_sl == 0.2 * m.fs.costing.cost_of_labor

    @m.fs.costing.Constraint(doc="Calculate FB. Assumed to be 25% of the cost of labor.")
    def calculate_fb_con(b):
        return m.fs.costing.fb == 0.25 * m.fs.costing.cost_of_labor

    @m.fs.costing.Constraint(
        doc="Calculate PT&I. Assumed to be 1% of the total plant cost."
    )
    def calculate_pt_and_i_con(b):
        return m.fs.costing.pt_and_i == 0.01 * m.fs.costing.total_plant_cost

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly fixed operating cost.",
    )
    def calculate_total_yearly_fixed_operating_costs_cons(b, t):
        return (
            m.fs.costing.total_fixed_operating_expense[t]
            == m.fs.costing.cost_of_labor
            + m.fs.costing.m_and_sm
            + m.fs.costing.sa_and_qa_qc
            + m.fs.costing.s_ip_r_and_d[t]
            + m.fs.costing.a_and_sl
            + m.fs.costing.fb
            + m.fs.costing.pt_and_i
        )

    @m.fs.costing.Constraint(
        m.fs.plant_life_range,
        doc="Calculates the total overnight cost expended in a given year. It is assumed that 10% is expended "
        "in the first year, 60% in the second year, and 30% in the third year.",
    )
    def calculate_total_overnight_cost_expended_cons(b, t):
        if t == m.fs.plant_start:
            return (
                m.fs.costing.total_overnight_cost_expended[t]
                == 0.1 * m.fs.costing.total_overnight_cost
            )
        elif t == (m.fs.plant_start + 1):
            return (
                m.fs.costing.total_overnight_cost_expended[t]
                == 0.6 * m.fs.costing.total_overnight_cost
            )
        elif t == (m.fs.plant_start + 2):
            return (
                m.fs.costing.total_overnight_cost_expended[t]
                == 0.3 * m.fs.costing.total_overnight_cost
            )
        else:
            return m.fs.costing.total_overnight_cost_expended[t] == 0

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the plant overhead each year. Assumed to be 20% of the fixed and variable "
        "operating costs.",
    )
    def calculate_plant_overhead_cons(b, t):
        return m.fs.costing.plant_overhead[t] == 0.2 * (
            m.fs.costing.total_var_operating_expense[t]
            + m.fs.costing.total_fixed_operating_expense[t]
        )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total operating expenses each year. Includes variable, fixed, and plant overhead.",
    )
    def calculate_total_operating_expense_cons(b, t):
        return (
            m.fs.costing.total_operating_expense[t]
            == m.fs.costing.total_var_operating_expense[t]
            + m.fs.costing.total_fixed_operating_expense[t]
            + m.fs.costing.plant_overhead[t]
        )

    @m.fs.costing.Constraint(
        m.fs.plant_life_range,
        doc="Calculate the cash flow for each year of the plant's lifetime.",
    )
    def calculate_yearly_costing(b, t):
        if t == m.fs.plant_start:
            return m.fs.costing.cash_flow[t] == -m.fs.costing.total_overnight_cost_expended[t]
        else:
            return m.fs.costing.cash_flow[t] == (
                m.fs.costing.total_profit[t] - m.fs.costing.total_operating_expense[t]
            ) * (1 + m.fs.costing.i_operating_expense_escalation) ** (
                t - m.fs.plant_start
            ) - m.fs.costing.total_overnight_cost_expended[
                t
            ] * (
                1 + m.fs.costing.i_capital_escalation
            ) ** (
                t - m.fs.plant_start
            )

    @m.fs.costing.Constraint(doc="Calculates the net present value of the plant.")
    def calculate_net_present_value_con(b):
        return m.fs.costing.npv == sum(
            m.fs.costing.cash_flow[t]
            / (
                (1 + m.fs.costing.discount_factor)
                ** (t - m.fs.plant_start)
            )
            for t in m.fs.plant_life_range
        )

    @m.fs.costing.cost_of_recovery.Constraint(
        doc="Sets the net present value to zero when the cost of recovery is chosen as an "
        "objective function."
    )
    def set_net_present_value_to_zero_con(b):
        return m.fs.costing.npv == 0

    ## Set objective function
    m.fs.costing.net_present_value.obj = pyo.Objective(
        expr=m.fs.costing.npv,
        sense=pyo.maximize,
        doc="Objective function is maximizing the net present value when the NPV is chosen.",
    )
    m.fs.costing.cost_of_recovery.obj = pyo.Objective(
        expr=m.fs.costing.cor,
        sense=pyo.minimize,
        doc="Objective function is minimizing the cost of recovery when the COR is chosen.",
    )


def add_environmental_impact_params(
    m, consider_environmental_impacts, options_environmental_impacts, epsilon
):
    """
    This function builds the environmental impact parameters.

    Args:
        m: pyomo model.
        consider_environmental_impacts: (bool) Choice of whether or not to consider environmental impacts.
        options_environmental_impacts: (dict) Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.
        epsilon: (float) Epsilon factor for generating the Pareto front.
    """
    ## Create block
    m.fs.environmental_impacts = pyo.Block(doc="Block to hold environmental impacts.")

    ## Pyomo parameters
    m.fs.environmental_impacts.consider_environmental_impacts = pyo.Param(
        initialize=consider_environmental_impacts,
        within=pyo.Boolean,
        doc="Choice of whether or not to consider environmental impacts.",
    )

    # Only add parameters if user wants to consider environmental impacts. Otherwise, add defaults.
    if consider_environmental_impacts:
        m.fs.environmental_impacts.options_environmental_impacts = pyo.Param(
            m.fs.all_opts_set,
            initialize=options_environmental_impacts,
            doc="Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.",
        )
        m.fs.environmental_impacts.epsilon = pyo.Param(
            initialize=epsilon,
            domain=pyo.NonNegativeReals,
            doc="Epsilon factor for generating the Pareto front.",
        )
    else:
        m.fs.environmental_impacts.options_environmental_impacts = pyo.Param(
            m.fs.all_opts_set,
            initialize={},
            doc="Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.",
        )
        m.fs.environmental_impacts.epsilon = pyo.Param(
            initialize=1,
            domain=pyo.NonNegativeReals,
            doc="Epsilon factor for generating the Pareto front.",
        )


def add_environmental_impact_vars(m):
    """
    This function builds the environmental impact parameters.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.fs.environmental_impacts.option_yearly_impacts = pyo.Var(
        m.fs.all_opts_set,
        m.fs.operational_range,
        doc="Yearly environmental impacts generated by each option.",
    )
    m.fs.environmental_impacts.total_yearly_impacts = pyo.Var(
        m.fs.operational_range,
        doc="Total yearly environmental impacts generated by the process.",
    )
    m.fs.environmental_impacts.total_impacts = pyo.Var(
        doc="The total environmental impacts generated by the process over its entire lifetime."
    )


def add_environmental_impact_cons(m):
    """
    This function builds the environmental impact constraints.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    # only add if user wants to consider environmental impacts. Otherwise use default values to not throw errors. Constraints will be deactivated anyways.
    if m.fs.environmental_impacts.consider_environmental_impacts:

        @m.fs.environmental_impacts.Constraint(
            m.fs.all_opts_set,
            m.fs.operational_range,
            doc="Calculate the yearly environmental impacts generated by each option.",
        )
        def calculate_opt_yearly_impacts_con(b, j, k, t):
            return (
                m.fs.environmental_impacts.option_yearly_impacts[j, k, t]
                == sum(
                    m.fs.f_in[j, k, c, t]
                    for c in m.fs.tracked_comps
                )
                * m.fs.environmental_impacts.options_environmental_impacts[j, k]
            )

        @m.fs.environmental_impacts.Constraint(
            m.fs.operational_range,
            doc="Calculate the total yearly environmental impacts generated by the process.",
        )
        def calculate_yearly_impacts_con(b, t):
            return m.fs.environmental_impacts.total_yearly_impacts[t] == sum(
                m.fs.environmental_impacts.option_yearly_impacts[opt, t]
                for opt in m.fs.all_opts_set
            )

        @m.fs.environmental_impacts.Constraint(
            doc="Calculate the total environmental impacts generated by the process throughout it's lifetime."
        )
        def calculate_total_impacts_con(b):
            return m.fs.environmental_impacts.total_impacts == sum(
                m.fs.environmental_impacts.total_yearly_impacts[t]
                for t in m.fs.operational_range
            )

        @m.fs.environmental_impacts.Constraint(
            doc="Epsilon constraint. Total process impacts must be less than epsilon."
        )
        def epsilon_con(b):
            return (
                m.fs.environmental_impacts.total_impacts <= m.fs.environmental_impacts.epsilon
            )

    else:

        @m.fs.environmental_impacts.Constraint(
            m.fs.all_opts_set,
            m.fs.operational_range,
            doc="Calculate the yearly environmental impacts generated by each option.",
        )
        def calculate_opt_yearly_impacts_con(b, j, k, t):
            return m.fs.environmental_impacts.option_yearly_impacts[j, k, t] == 0

        @m.fs.environmental_impacts.Constraint(
            m.fs.operational_range,
            doc="Calculate the total yearly environmental impacts generated by the process.",
        )
        def calculate_yearly_impacts_con(b, t):
            return m.fs.environmental_impacts.total_yearly_impacts[t] == 0

        @m.fs.environmental_impacts.Constraint(
            doc="Calculate the total environmental impacts generated by the process throughout it's lifetime."
        )
        def calculate_total_impacts_con(b):
            return m.fs.environmental_impacts.total_impacts == 0

        @m.fs.environmental_impacts.Constraint(
            doc="Epsilon constraint. Total process impacts must be less than epsilon."
        )
        def epsilon_con(b):
            return m.fs.environmental_impacts.total_impacts <= 1


def add_byproduct_valorization_params(
    m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
):
    """
    This function builds the byproduct valorization parameters.

    Args:
        m: pyomo model.
        consider_byproduct_valorization: (bool) Decide whether or not to consider the valorization of byproducts.
        byproduct_values: (dict) Byproducts considered, and their value ($/kg).
        byproduct_opt_conversions: (dict) Conversion factors for each byproduct for each option.
    """
    ## Create blocks
    m.fs.byproduct_valorization = pyo.Block()

    ### Define parameters from user input.
    m.fs.byproduct_valorization.consider_byproduct_valorization = pyo.Param(
        initialize=consider_byproduct_valorization,
        within=pyo.Boolean,
        doc="Choice of whether or not to consider byproduct valorization.",
    )
    ## Only add parameters if user wants to consider byproduct valorization.
    if consider_byproduct_valorization:
        # Create a set of the byproducts considered.
        byproducts = set(byproduct_values.keys())
        # Define a dictionary that maps all byproducts to the options that produce it.
        byproduct_producing_opts = {}
        # First, iterate through the byproduct_opt_conversion dict. Extract the options, and the innter dict
        # which contains the byproducts it produces as keys, and the corresponding conversion factor.
        for outer_key, inner_dict in byproduct_opt_conversions.items():
            # Extract the byproducts from the inner dict.
            for inner_key in inner_dict:
                # Append opt to the list associated with the byproduct. Create a new list if the byproduct entry doesn't exist.
                byproduct_producing_opts.setdefault(inner_key, []).append(outer_key)

        ### Define functions needed to initialize pyomo parameters.
        # Define a function for initializing byproduct option conversion pyomo parameter
        def byproduct_opt_conversion_initialize(b, j, k, byproduct):
            return byproduct_opt_conversions[(j, k)][byproduct]

        ### Define necessary pyomo parameters.
        ## Pyomo parameters
        m.fs.byproduct_valorization.byproducts_set = pyo.Set(
            initialize=byproducts, doc="Set of byproducts considered."
        )
        m.fs.byproduct_valorization.byproduct_values = pyo.Param(
            m.fs.byproduct_valorization.byproducts_set,
            initialize=byproduct_values,
            doc="Defines the value of each byproduct ($/kg).",
        )
        m.fs.byproduct_valorization.byproduct_opts_set = pyo.Set(
            initialize=byproduct_opt_conversions.keys(),
            doc="Set of options that produce byproducts.",
        )
        m.fs.byproduct_valorization.opt_byproduct_set = pyo.Set(
            initialize=(
                (opt, byproduct)
                for opt in byproduct_opt_conversions.keys()
                for byproduct in byproduct_opt_conversions[opt].keys()
            ),
            doc="Set of options that produce byproducts, and the byproducts that they produce.",
        )

        m.fs.byproduct_valorization.byproduct_opt_conversion = pyo.Param(
            m.fs.byproduct_valorization.opt_byproduct_set,
            initialize=byproduct_opt_conversion_initialize,
            doc="Defines the conversion factors for all byproducts for all options that produce them.",
        )
        m.fs.byproduct_valorization.byproduct_producing_opts = pyo.Param(
            m.fs.byproduct_valorization.byproducts_set,
            initialize=byproduct_producing_opts,
            doc="Pyomo parameter which defines what options produce each byproduct.",
        )


def add_byproduct_valorization_vars(m):
    """
    This function builds the byproduct valorization variables.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables

    # Only build variables if user wants to consider byproduct valorization. Otherwise, will build non-indexed vars to not throw errors. Cons will be deactivated anyways.
    if m.fs.byproduct_valorization.consider_byproduct_valorization:
        m.fs.byproduct_valorization.byproduct_produced = pyo.Var(
            m.fs.byproduct_valorization.byproducts_set
            * m.fs.operational_range,
            domain=pyo.NonNegativeReals,
            doc="The amount of each byproduct produced each year.",
        )
        m.fs.byproduct_valorization.byproduct_profit = pyo.Var(
            m.fs.byproduct_valorization.byproducts_set
            * m.fs.operational_range,
            domain=pyo.Reals,
            doc="The amount of profit generated from each byproduct each year.",
        )
    else:
        m.fs.byproduct_valorization.byproduct_produced = pyo.Var(
            domain=pyo.NonNegativeReals,
            doc="The amount of each byproduct produced each year.",
        )
        m.fs.byproduct_valorization.byproduct_profit = pyo.Var(
            domain=pyo.Reals,
            doc="The amount of profit generated from each byproduct each year.",
        )


def add_byproduct_valorization_cons(m):
    """
    This function builds the byproduct valorization constraints.

    Args:
        m: pyomo model.
    """
    ## Create blocks
    m.fs.no_byproduct_valorization = pyo.Block()

    ## Pyomo constraints
    @m.fs.no_byproduct_valorization.Constraint(
        m.fs.operational_range,
        doc="Sets the profit from the valorization of byproducts to zero.",
    )
    def calculate_opt_byprod_val_cons(b, t):
        return m.fs.costing.total_byproduct_profit[t] == 0

    # only build constraints if user wants to consider byproduct valorization.
    if m.fs.byproduct_valorization.consider_byproduct_valorization:

        @m.fs.byproduct_valorization.Constraint(
            m.fs.byproduct_valorization.byproducts_set,
            m.fs.operational_range,
            doc="Calculates the amount of each byproduct produced each year.",
        )
        def calculate_byproduct_produced_cons(b, byprod, t):
            return m.fs.byproduct_valorization.byproduct_produced[byprod, t] == sum(
                sum(
                    m.fs.f_in[opt, c, t] for c in m.fs.tracked_comps
                )
                * m.fs.byproduct_valorization.byproduct_opt_conversion[opt, byprod]
                for opt in m.fs.byproduct_valorization.byproduct_producing_opts[
                    byprod
                ]
            )

        @m.fs.byproduct_valorization.Constraint(
            m.fs.byproduct_valorization.byproducts_set,
            m.fs.operational_range,
            doc="Calculate the yearly profit from each byproduct.",
        )
        def calculate_byproduct_profit_cons(b, byprod, t):
            return (
                m.fs.byproduct_valorization.byproduct_profit[byprod, t]
                == m.fs.byproduct_valorization.byproduct_produced[byprod, t]
                * m.fs.byproduct_valorization.byproduct_values[byprod]
            )

        @m.fs.byproduct_valorization.Constraint(
            m.fs.operational_range,
            doc="Calculates the total yearly profit from the valorization of byproducts.",
        )
        def calculate_opt_byprod_val_cons(b, t):
            return m.fs.costing.total_byproduct_profit[t] == sum(
                m.fs.byproduct_valorization.byproduct_profit[byprod, t]
                for byprod in m.fs.byproduct_valorization.byproducts_set
            )
