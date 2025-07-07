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
from pyomo.environ import units as pyunits


def add_plant_lifetime_params_block(m, plant_start, plant_lifetime):
    """
    This function builds the rest of the plant lifetime parameters from the ones provided by the user.

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
        initialize=plant_start,
        doc="The year that plant construction begins.",
        units=pyunits.year,
    )
    m.fs.plant_lifetime = pyo.Param(
        initialize=plant_lifetime,
        doc="The total lifetime of the plant, including plant construction. Must be at least three years.",
        units=pyunits.year,
    )
    m.fs.prod_start = pyo.Param(
        initialize=prod_start,
        doc="The first year of plant production.",
        units=pyunits.year,
    )
    m.fs.plant_end = pyo.Param(
        initialize=plant_end,
        doc="The final year of plant production.",
        units=pyunits.year,
    )
    m.fs.plant_life_range = pyo.RangeSet(
        pyo.value(plant_start),
        pyo.value(plant_end),
        doc="Lifetime of the plant.",
    )
    m.fs.operational_range = pyo.RangeSet(
        pyo.value(prod_start),
        pyo.value(plant_end),
        doc="Operational lifetime of the plant.",
    )


def add_feed_params_block(
    m, available_feed, collection_rate, tracked_comps, prod_comp_mass
):
    """
    This function builds the rest of the feed parameters from the ones provided by the user.

    Args:
        m: pyomo model.
        available_feed: (dict) Total feedstock available (nunber of EOL products) for recycling each year.
        collection_rate: (float) Collection rate for how much available feed is processed by the plant each year.
        tracked_comps: (list) List of tracked components.
        prod_comp_mass: (dict) Mass of tracked components per EOL product (kg per EOL product).
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
        units=pyunits.EOL_Product / pyunits.year,
    )
    m.fs.collection_rate = pyo.Param(
        initialize=collection_rate,
        doc="The fraction of available feed that is processed by the plant each year.",
    )
    m.fs.tracked_comps = pyo.Set(initialize=tracked_comps, doc="Tracked components.")
    m.fs.prod_comp_mass = pyo.Param(
        m.fs.tracked_comps,
        initialize=prod_comp_mass,
        doc="The mass of each tracked component per EOL product.",
        units=pyunits.kg / pyunits.EOL_Product,
    )
    m.fs.feed_entering = pyo.Param(
        m.fs.operational_range,
        initialize=feed_entering,
        doc="The amount of EOL Products entering the plant each year.",
        units=pyunits.EOL_Product / pyunits.year,
    )
    m.fs.max_feed_entering = pyo.Param(
        initialize=max_feed_entering,
        doc="The max yearly feed that enters the plant over the production period.",
        units=pyunits.EOL_Product / pyunits.year,
    )
    m.fs.max_feed_entering_year = pyo.Param(
        initialize=max_feed_entering_year,
        doc="The year that the max feed enters the plant.",
        units=pyunits.year,
    )


def add_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, option_efficiencies
):
    """
    This function builds the rest of the superstructure formulation parameters from the ones provided by the user.

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
    m.fs.max_options = pyo.RangeSet(
        1,
        max_options,
        doc="Set of the max number of options that exist in a single stage in the superstructure.",
    )
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
            (j, k) for j in m.fs.stages_set for k in range(1, options_in_stage[j] + 1)
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
    # m.fs.option_outlets = pyo.Param(
    #     m.fs.all_opts_set,
    #     initialize=option_outlets,
    #     doc="Defines the set of options k' in stage j+1 connected to option k in stage j.",
    # )
    m.fs.option_outlet_pairs = pyo.Set(
        initialize=[(j, k, kp) for (j, k), kps in option_outlets.items() for kp in kps],
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
    This function builds the rest of the operating parameters from the ones provided by the user.

    Args:
        m: pyomo model.
        profit: (dict) Profit per unit of product in terms of tracked components ($/kg tracked component).
        opt_var_oc_params: (dict) Holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        operators_per_discrete_unit: (dict) Number of operators needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit ($/year) for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit ($) for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of end-of-life products disassembled per year per unit (number of EOL products / year).
        num_operators: (dict) Number of operators needed for each continuous option.
        labor_rate: (float) Yearly wage per operator ($ / year).
    """
    ### Define parameters from user input.
    ## Calculate the number of discrete units needed for each option. There must be enough discrete units to handle the maximum flow of end-of-life products
    # that enters the plant over its operational lifetime.
    discrete_units_per_option = copy.deepcopy(processing_rate)
    for key in discrete_units_per_option.keys():
        discrete_units_per_option[key] = math.ceil(
            pyo.value(m.fs.max_feed_entering) / processing_rate[key]
        )

    ### Define functions needed to initialize pyomo parameters.
    # Define a function for initializing profit pyomo parameter.
    def profit_initialization(m, j, k, c):
        return profit[(j, k)][c]

    # Define a function for initializing opt_var_oc_params pyomo A parameter.
    def opt_var_oc_A_param_initialization(m, j, k):
        return opt_var_oc_params[(j, k)]["a"]

    # Define a function for initializing opt_var_oc_params pyomo B parameter.
    def opt_var_oc_B_param_initialization(m, j, k):
        return opt_var_oc_params[(j, k)]["b"]

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
        units=pyunits.USD / pyunits.kg,
    )
    m.fs.costing.opt_var_oc_param_A = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=opt_var_oc_A_param_initialization,
        doc="Holds all the 'A' variable operating costs parameter values for all continuous options.",
        units=pyunits.USD / pyunits.kg,
    )
    m.fs.costing.opt_var_oc_param_B = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=opt_var_oc_B_param_initialization,
        doc="Holds all the 'B' variable operating costs parameter values for all continuous options.",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.operators_per_discrete_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=operators_per_discrete_unit,
        doc="The number of operators needed per discrete unit for discrete options.",
        units=pyunits.Operator / pyunits.Disassembly_Unit,
    )
    m.fs.costing.yearly_cost_per_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=yearly_cost_per_unit,
        doc="The operating costs per discrete unit for discrete options.",
        units=pyunits.USD / pyunits.Disassembly_Unit / pyunits.year,
    )
    m.fs.costing.capital_cost_per_unit = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=capital_cost_per_unit,
        doc="The capital cost per discrete unit for discrete options.",
        units=pyunits.USD / pyunits.Disassembly_Unit,
    )
    m.fs.costing.processing_rate = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=processing_rate,
        doc="The processing rate per discrete unit for discrete options. In terms of number of end-of-life products "
        "disassembled per discrete unit per year.",
        units=pyunits.EOL_Product / pyunits.Disassembly_Unit / pyunits.year,
    )
    m.fs.costing.num_operators = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=num_operators,
        doc="The number of operators per continuous options.",
        units=pyunits.Operator,
    )
    m.fs.costing.labor_rate = pyo.Param(
        initialize=labor_rate,
        doc="The yearly wage per operator.",
        units=pyunits.USD / pyunits.Operator / pyunits.year,
    )
    m.fs.costing.discrete_units_per_option = pyo.Param(
        m.fs.discrete_opts_set,
        initialize=discrete_units_per_option,
        doc="The number of discrete units per option needed to disassemble all the incoming end-of-life products over the operational "
        "lifetime of the plant.",
        units=pyunits.Disassembly_Unit,
    )


def add_discretized_costing_params(m, discretized_equipment_cost):
    """
    This function adds all the discretized costing parameters.

    Args:
        m: pyomo model.
        discretized_equipment_cost: (dict) Discretized cost by flows entering for each continuous option ($/kg).
    """

    ### Define functions needed to initialize pyomo parameters.
    # Define a function that returns the number of discretized data points that exist for a given option.
    def opt_discretized_costing_datapoints_initialization(m, j, k):
        return len(discretized_equipment_cost[(j, k)]["Flowrates"])

    # Define a function that returns the corresponding flowrate for a given option and datapoint number.
    def opt_discretized_costing_flowrates_initialization(m, j, k, dp):
        return discretized_equipment_cost[(j, k)]["Flowrates"][dp]

    # Define a function that returns the corresponding purchased equipment cost for a given option and datapoint number.
    def opt_discretized_costing_costs_initialization(m, j, k, dp):
        return discretized_equipment_cost[(j, k)]["Costs"][dp]

    ### Define necessary pyomo parameters.
    ## Pyomo parameters
    m.fs.costing.opt_discretized_costing_datapoints = pyo.Param(
        m.fs.continuous_opts_set,
        initialize=opt_discretized_costing_datapoints_initialization,
        doc="Holds the number of discretized costing datapoints that exist for all continuous options.",
    )
    m.fs.continuous_opts_discretized_costing_data_points = pyo.Set(
        initialize=(
            opt + (dp,)
            for opt in m.fs.continuous_opts_set
            for dp in range(len(discretized_equipment_cost[opt]["Flowrates"]))
        ),
        doc="Stores the discretized costing datapoints that exist for each continuous option.",
    )
    m.fs.costing.flowrates_data = pyo.Param(
        m.fs.continuous_opts_discretized_costing_data_points,
        initialize=opt_discretized_costing_flowrates_initialization,
        doc="Discretized flowrate data for all continuous options.",
        units=pyunits.kg / pyunits.year,
    )
    m.fs.costing.costs_data = pyo.Param(
        m.fs.continuous_opts_discretized_costing_data_points,
        initialize=opt_discretized_costing_costs_initialization,
        doc="Discretized costing data for all continuous options.",
        units=pyunits.USD,
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
            pyo.value(m.fs.max_feed_entering) * pyo.value(m.fs.prod_comp_mass[c])
        )

    # Create an upper bound for the total flow entering a continuous option.
    max_flow_upper_bound = sum(m_val[c] for c in m.fs.tracked_comps)

    ## Pyomo parameters
    # Flow entering each stage (except the last stage).
    m.fs.f_stages = pyo.RangeSet(
        1,
        m.fs.num_stages - 1,
        doc="Set of all stages except the last. Used to define the flow variable: 'f'.",
    )
    m.fs.flow_set = pyo.Set(
        initialize=m.fs.all_opts_set * m.fs.tracked_comps * m.fs.operational_range,
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
        units=pyunits.kg / pyunits.year,
    )


def add_mass_balance_vars(m):
    """
    This function builds the mass balance variables.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.fs.f = pyo.Var(
        m.fs.f_stages * m.fs.tracked_comps * m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each stage (except the last stage). See documentation for more details.",
        units=pyunits.kg / pyunits.year,
    )
    m.fs.f_in = pyo.Var(
        m.fs.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
        units=pyunits.kg / pyunits.year,
    )
    m.fs.f_out = pyo.Var(
        m.fs.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
        units=pyunits.kg / pyunits.year,
    )
    m.fs.option_binary_var = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.Binary,
        doc="Binary variables to indicate whether or not an option has been selected.",
    )
    m.fs.piecewise_flow_entering = pyo.Var(
        m.fs.continuous_opts_set,
        bounds=(0, m.fs.max_flow_upper_bound),
        doc="The max total flow that enters each continuous option over the lifetime"
        "of the plant. Used in piecewise constraints and thus must be unitless.",
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
            return b.f[j - 1, c, t] == sum(b.f_in[j, k, c, t] for k in num_options)

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
            return b.feed_entering[t] * b.prod_comp_mass[c] == sum(
                b.f_in[j, k, c, t] for k in num_options
            )
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.flow_set,
        doc="Equation (3) from the documentation. The steam exiting an option is related to the inlet stream by an efficiency parameter.",
    )
    def intermediate_flow_cons(b, j, k, c, t):
        alpha = b.option_efficiencies[j, k, c]
        return b.f_in[j, k, c, t] * alpha == b.f_out[j, k, c, t]

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
            return sum(b.f_out[j, k, c, t] for k in num_options) == b.f[j, c, t]
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.stages_set,
        doc="Equation (5) from the documentation. Ensures that only one option per stage is chosen.",
    )
    def stage_binary_cons(b, j):
        # Extract all the options available in stage 'j'
        num_options = range(1, b.options_in_stage[j] + 1)
        return sum(b.option_binary_var[j, k] for k in num_options) == 1

    # @m.fs.Constraint(
    #     m.fs.all_opts_set,
    #     doc="Equation (6) from the documentation. Ensures that if an option in stage $j + 1$ can only be chosen if it is connected to an option that was chosen in stage $j$.",
    # )
    # def connection_binary_cons(b, j, k):
    #     if j != m.fs.num_stages:
    #         # Extract the set of options k' in stage j+1 connected to option k in stage j.
    #         opt_connections = b.option_outlets[j, k]
    #         return b.option_binary_var[j, k] <= sum(
    #             b.option_binary_var[j + 1, kp] for kp in opt_connections
    #         )
    #     else:
    #         return pyo.Constraint.Skip
    @m.fs.Constraint(
        m.fs.all_opts_set,
    )
    def connetion_binary_cons(b, j, k):
        if j != m.fs.num_stages:
            return b.option_binary_var[j, k] <= sum(
                b.option_binary_var[j + 1, kp]
                for kp in m.fs.max_options
                if (j, k, kp) in m.fs.option_outlet_pairs
            )
        else:
            return pyo.Constraint.Skip

    @m.fs.Constraint(
        m.fs.flow_set,
        doc="Equation (7) from the documentation. Big-M constraint for inlet streams. Inlet streams must be zero if the option is not chosen.",
    )
    def f_in_big_m_cons(b, j, k, c, t):
        return b.f_in[j, k, c, t] <= b.option_binary_var[j, k] * b.big_m_val[c]

    @m.fs.Constraint(
        m.fs.flow_set,
        doc="Equation (8) from the documentation. Big-M constraint for outlet streams. Outlet streams must be zero if the option is not chosen.",
    )
    def f_out_big_m_cons(b, j, k, c, t):
        return b.f_out[j, k, c, t] <= b.option_binary_var[j, k] * b.big_m_val[c]

    @m.fs.Constraint(
        m.fs.continuous_opts_set,
        m.fs.operational_range,
        doc="Constraint to determine the max flow entering each continuous option over the lifetime of the plant.",
    )
    def max_flow_entering_cons(b, j, k, t):
        return b.piecewise_flow_entering[j, k] >= sum(
            b.f_in[j, k, c, t] for c in b.tracked_comps
        )


def add_costing_params(m):
    """
    This function builds the costing parameters.

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed. For the objective functions,
    the user has a choice of either maximizing the net present value, or minimizing the cost of recovery.

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    Args:
        m: pyomo model.
    """
    ## Define parameters
    # Define a Pyomo Param for fraction of total overnight cost expended
    total_overnight_capital_fraction_expended = {
        pyo.value(m.fs.plant_start): 0.1,
        pyo.value(m.fs.plant_start) + 1: 0.6,
        pyo.value(m.fs.plant_start) + 2: 0.3,
    }

    ## Pyomo parameters
    m.fs.costing.lang_factor = pyo.Param(
        initialize=2.97,
        mutable=True,
        doc="Lang factor for calculating total plant cost from equipment costs. Assumed to be 2.97[2].",
    )
    m.fs.costing.i_operating_expense_escalation = pyo.Param(
        initialize=0.03,
        mutable=True,
        doc="Operating expenses escalation rate. Assumed to be 3%[2].",
    )
    m.fs.costing.i_capital_escalation = pyo.Param(
        initialize=0.036,
        mutable=True,
        doc="Capital expenses escalation rate. Assumed to be 3.6%[2].",
    )
    m.fs.costing.discount_factor = pyo.Param(
        initialize=0.0577,
        mutable=True,
        doc="Discount factor for calculate the net present value. The after-tax weighted average cost of capital was used (ATWACC)[1].",
    )
    m.fs.costing.financing_factor = pyo.Param(
        initialize=0.027,
        mutable=True,
        doc="Factor for calculating financing costs from total plant cost. Assumed to be 2.7%[1].",
    )
    m.fs.costing.other_costs_factor = pyo.Param(
        initialize=0.15,
        mutable=True,
        doc="Factor for calculating the 'other costs', as defined by QGESS[1] from the total plant cost. Assumed to be 15%[1].",
    )
    m.fs.costing.m_and_sm_costing_factor = pyo.Param(
        initialize=0.02,
        mutable=True,
        doc="Factor for calculating Maintenance & Supply Materials. Assumed to be 2% of the total plant cost[2].",
    )
    m.fs.costing.sa_and_qa_qc_costing_factor = pyo.Param(
        initialize=0.1,
        mutable=True,
        doc="Factor for calculating Sample Analysis & Quality Assurance/Quality Control. Assumed to be 10% of the cost of labor[2].",
    )
    m.fs.costing.s_ip_r_and_d_costing_factor = pyo.Param(
        initialize=0.01,
        mutable=True,
        doc="Factor for calculating Sales, Intellectual Property, and Research & Development. Assumed to be 1% of the total profit[2].",
    )
    m.fs.costing.a_and_sl_costing_factor = pyo.Param(
        initialize=0.2,
        mutable=True,
        doc="Factor for calculating Administrative & Supporting Labor. Assumed to be 20% of the cost of labor[2].",
    )
    m.fs.costing.fb_costing_factor = pyo.Param(
        initialize=0.25,
        mutable=True,
        doc="Factor for calculating Fringe Benefits. Assumed to be 25% of the cost of labor[2].",
    )
    m.fs.costing.pt_and_i_costing_factor = pyo.Param(
        initialize=0.01,
        mutable=True,
        doc="Factor for calculating Property Taxes & Insurance. Assumed to be 1% of the total plant cost[2].",
    )
    m.fs.costing.total_overnight_capital_fraction_expended = pyo.Param(
        m.fs.plant_life_range,
        initialize=lambda m, t: total_overnight_capital_fraction_expended.get(t, 0),
        mutable=True,
        doc="Fraction of overnight cost expended in each year. It is assumed that 10% is expended "
        "in the first year, 60% in the second year, and 30% in the third year.",
        units=1 / pyunits.year,
    )
    m.fs.costing.plant_overhead_factor = pyo.Param(
        initialize=0.2,
        mutable=True,
        doc="Factor for calculating the plant overhead. Assumed to be 20% of the total operating costs each year[2].",
    )


def add_costing_vars(m, obj_func: str):
    """
    This function builds the costing variables.

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
    """
    ## Pyomo variables
    # if NPV is objective function, build relevant vars
    if obj_func == "NPV":
        m.fs.costing.opt_profit = pyo.Var(
            m.fs.final_opts_set,
            m.fs.operational_range,
            domain=pyo.NonNegativeReals,
            doc="The profit generated by each option in the final processing stage each year.",
            units=pyunits.USD / pyunits.year,
        )

    # if cor is objective function, build relevant vars
    else:
        m.fs.costing.cost_of_recovery = pyo.Var(
            domain=pyo.NonNegativeReals,
            doc="The cost of recovery.",
            units=pyunits.USD / pyunits.kg,
        )

    # build vars that're needed regardless of the objective function.
    m.fs.costing.net_present_value = pyo.Var(
        domain=pyo.Reals, doc="The net present value.", units=pyunits.USD
    )
    m.fs.costing.total_profit = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total profit generated by the plant each year.",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.piecewise_equipment_cost = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The cost of purchased equipment. Calculated with piecewise constraints and thus must be dimensionless.",
    )
    m.fs.costing.equipment_cost = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The cost of purchased equipment with units.",
        units=pyunits.USD,
    )
    m.fs.costing.total_plant_cost = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="The total plant cost, as defined by QGESS[1].",
        units=pyunits.USD,
    )
    m.fs.costing.financing = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="The total financing cost of the plant, as defined by QGESS[1].",
        units=pyunits.USD,
    )
    m.fs.costing.other_costs = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="'Other costs' associated with the plant as defined by QGESS[1].",
        units=pyunits.USD,
    )
    m.fs.costing.total_overnight_cost = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="The total overnight cost of the plant, as defined by QGESS[1].",
        units=pyunits.USD,
    )
    m.fs.costing.opt_variable_operating_cost = pyo.Var(
        m.fs.all_opts_set * m.fs.operational_range,
        domain=pyo.Reals,
        doc="Yearly variable operating expense for each option[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.aggregate_variable_operating_cost = pyo.Var(
        m.fs.operational_range,
        domain=pyo.Reals,
        doc="Total yearly variable operating expense[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.operators_per_option = pyo.Var(
        m.fs.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The number of operators needed for each option[3].",
        units=pyunits.Operator,
    )
    m.fs.costing.total_operators = pyo.Var(
        domain=pyo.NonNegativeIntegers,
        doc="The total number of operators needed for the process[3]. Must be an integer value.",
        units=pyunits.Operator,
    )
    m.fs.costing.cost_of_labor = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="The cost of labor for the process[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.m_and_sm = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Maintenance & Supply Materials (M&SM)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.sa_and_qa_qc = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Sample Analysis & Quality Assurance/Quality Control (SA&QA/QC)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.s_ip_r_and_d = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Sales, Intellectual Property, and Research & Development (S,IP,R&D)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.a_and_sl = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Administrative & Supporting Labor (A&SL)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.fb = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Fringe Benefits (FB)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.pt_and_i = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Property Taxes & Insurance (PT&I)[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.aggregate_fixed_operating_cost = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Total yearly fixed operating cost[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.total_overnight_cost_expended = pyo.Var(
        m.fs.plant_life_range,
        domain=pyo.NonNegativeReals,
        doc="The total overnight cost expended each year[1].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.plant_overhead = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The yearly plant overhead[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.total_operating_expense = pyo.Var(
        m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total operating expense each year[2].",
        units=pyunits.USD / pyunits.year,
    )
    m.fs.costing.cash_flow = pyo.Var(
        m.fs.plant_life_range,
        domain=pyo.Reals,
        doc="The yearly cash flow[4].",
        units=pyunits.USD / pyunits.year,
    )


def add_byproduct_valorization_params(m, byproduct_values, byproduct_opt_conversions):
    """
    This function builds the byproduct valorization parameters.

    In byproduct valorization, the profit generated from byproducts is considered in addition to the profit generated from the main product.

    Args:
        m: pyomo model.
        byproduct_values: (dict) Byproducts considered, and their value ($/kg).
        byproduct_opt_conversions: (dict) Conversion factors for each byproduct for each option.
    """
    ## Create blocks
    m.fs.byproduct_valorization = pyo.Block()

    ### Define parameters from user input.
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
    ## Add params related to calculating amount of byproduct produced to byproduct_valorization block.
    m.fs.byproduct_valorization.byproducts_set = pyo.Set(
        initialize=byproducts, doc="Set of byproducts considered."
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

    ## Add params related to profit generated from byproducts to costing block.
    m.fs.costing.byproduct_values = pyo.Param(
        m.fs.byproduct_valorization.byproducts_set,
        initialize=byproduct_values,
        doc="Defines the value of each byproduct ($/kg).",
    )


def add_byproduct_valorization_vars(m):
    """
    This function builds the byproduct valorization variables.

    In byproduct valorization, the profit generated from byproducts is considered in addition to the profit generated from the main product.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    ## Add vars related to calculating amount of byproduct produced to byproduct_valorization block.
    m.fs.byproduct_valorization.byproduct_produced = pyo.Var(
        m.fs.byproduct_valorization.byproducts_set * m.fs.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The amount of each byproduct produced each year.",
    )

    ## Add vars related to profit generated from byproducts to costing block.
    m.fs.costing.byproduct_profit = pyo.Var(
        m.fs.byproduct_valorization.byproducts_set * m.fs.operational_range,
        domain=pyo.Reals,
        doc="The amount of profit generated from each byproduct each year.",
    )
    m.fs.costing.total_byproduct_profit = pyo.Var(
        m.fs.operational_range,
        domain=pyo.Reals,
        doc="The total profit generated by the plant each year from the valorization of byproducts.",
    )


def add_byproduct_valorization_cons(m):
    """
    This function builds the byproduct valorization constraints.

    In byproduct valorization, the profit generated from byproducts is considered in addition to the profit generated from the main product.

    Args:
        m: pyomo model.
    """

    ## Add cons related to calculating amount of byproduct produced to byproduct_valorization block.
    @m.fs.byproduct_valorization.Constraint(
        m.fs.byproduct_valorization.byproducts_set,
        m.fs.operational_range,
        doc="Calculates the amount of each byproduct produced each year.",
    )
    def calculate_byproduct_produced_cons(b, byprod, t):
        return m.fs.byproduct_valorization.byproduct_produced[byprod, t] == sum(
            sum(m.fs.f_in[opt, c, t] for c in m.fs.tracked_comps)
            * m.fs.byproduct_valorization.byproduct_opt_conversion[opt, byprod]
            for opt in m.fs.byproduct_valorization.byproduct_producing_opts[byprod]
        )


def add_profit_cons(m, obj_func: str, consider_byproduct_valorization: bool):
    """
    This function builds the profit constraints.

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
        consider_byproduct_valorization: (bool) Decide whether or not to consider the valorization of byproducts.
    """
    ## Pyomo constraints
    # Check if NPV is objective function
    if obj_func == "NPV":

        @m.fs.costing.Constraint(
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

    # Check if profit from byproducts is considered
    if consider_byproduct_valorization:
        # add cons for calculating profit from byproducts if so
        @m.fs.costing.Constraint(
            m.fs.byproduct_valorization.byproducts_set,
            m.fs.operational_range,
            doc="Calculate the yearly profit from each byproduct.",
        )
        def calculate_byproduct_profit_cons(b, byprod, t):
            return (
                m.fs.costing.byproduct_profit[byprod, t]
                == m.fs.byproduct_valorization.byproduct_produced[byprod, t]
                * m.fs.costing.byproduct_values[byprod]
            )

        @m.fs.costing.Constraint(
            m.fs.operational_range,
            doc="Calculates the total yearly profit from the valorization of byproducts.",
        )
        def calculate_opt_byprod_val_cons(b, t):
            return m.fs.costing.total_byproduct_profit[t] == sum(
                m.fs.costing.byproduct_profit[byprod, t]
                for byprod in m.fs.byproduct_valorization.byproducts_set
            )

        # Check if NPV is objective function
        if obj_func == "NPV":
            # Add profit constraint for when byproduct valorization is considered if so.
            @m.fs.costing.Constraint(
                m.fs.operational_range,
                doc="Calculates the total yearly profit when the net present value is selected as the objective function, and byproduct valorization is considered.",
            )
            def calculate_total_profit_cons(b, t):
                return (
                    m.fs.costing.total_profit[t]
                    == sum(
                        m.fs.costing.opt_profit[opt, t] for opt in m.fs.final_opts_set
                    )
                    + m.fs.costing.total_byproduct_profit[t]
                )

        # COR is objective function
        else:

            @m.fs.costing.Constraint(
                m.fs.operational_range,
                doc="Calculates the total yearly profit when the cost of recovery is selected as the objective function, and byproduct valorization is considered.",
            )
            def calculate_total_profit_cons(b, t):
                return (
                    m.fs.costing.total_profit[t]
                    == m.fs.costing.cost_of_recovery
                    * sum(
                        sum(m.fs.f_out[opt, c, t] for c in m.fs.tracked_comps)
                        for opt in m.fs.final_opts_set
                    )
                    + m.fs.costing.byproduct_valorization.total_byproduct_profit[t]
                )

    # don't include profit from byproducts
    else:
        # Check if NPV is objective function
        if obj_func == "NPV":

            @m.fs.costing.Constraint(
                m.fs.operational_range,
                doc="Calculates the total yearly profit when the net present value is selected as the objective function, and byproduct valorization is not considered.",
            )
            def calculate_total_profit_cons(b, t):
                return m.fs.costing.total_profit[t] == sum(
                    m.fs.costing.opt_profit[opt, t] for opt in m.fs.final_opts_set
                )

        # COR is objective function
        else:

            @m.fs.costing.Constraint(
                m.fs.operational_range,
                doc="Calculates the total yearly profit when the cost of recovery is selected as the objective function, and byproduct valorization is not considered.",
            )
            def calculate_total_profit_cons(b, t):
                return m.fs.costing.total_profit[
                    t
                ] == m.fs.costing.cost_of_recovery * sum(
                    sum(m.fs.f_out[opt, c, t] for c in m.fs.tracked_comps)
                    for opt in m.fs.final_opts_set
                )


def add_capital_cost_cons(m):
    """
    This function builds the capital cost constraints.

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
    """

    @m.fs.costing.Constraint(
        m.fs.discrete_opts_set,
        doc="Calculates the purchased cost of equipment for all discrete options. "
        "Done by multiplying the number of discrete units by the capital cost per unit.",
    )
    def discrete_opts_equipment_cost_cons(b, j, k):
        return (
            m.fs.costing.equipment_cost[j, k]
            == m.fs.costing.discrete_units_per_option[j, k]
            * m.fs.costing.capital_cost_per_unit[j, k]
        )

    def piecewise_rule(b, j, k):
        # Create a range encompassing the number of discretized costing data points that exist for a particular option.
        dp_range = range(m.fs.costing.opt_discretized_costing_datapoints[j, k])
        # Create a list of the flow data for the option. Necessary to use pyo.Piecewise.
        flow_data = [
            pyo.value(m.fs.costing.flowrates_data[j, k, dp]) for dp in dp_range
        ]
        # Create a list of the cost data for the option. Necessary to use pyo.Piecewise.
        equipment_cost_data = [
            pyo.value(m.fs.costing.costs_data[j, k, dp]) for dp in dp_range
        ]

        # use m.add_component to generate all piecewise functions
        # piecewise = Piecewise(yval, xval, *kwargs)
        piecewise = pyo.Piecewise(
            m.fs.costing.piecewise_equipment_cost[j, k],
            m.fs.piecewise_flow_entering[j, k],
            pw_pts=flow_data,
            pw_constr_type="EQ",
            f_rule=equipment_cost_data,
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
        m.fs.continuous_opts_set,
        doc="Adds units to equipment costs calculated from piecewise constraints.",
    )
    def add_units_to_piecewise_costs(b, j, k):
        return (
            m.fs.costing.equipment_cost[j, k]
            == m.fs.costing.piecewise_equipment_cost[j, k] * pyunits.USD
        )

    @m.fs.costing.Constraint(doc="Calculates the total plant cost[2].")
    def calculate_total_plant_cost_con(b):
        return m.fs.costing.total_plant_cost == sum(
            m.fs.costing.equipment_cost[opt] * pyunits.USD
            for opt in m.fs.discrete_opts_set
        ) + m.fs.costing.lang_factor * sum(
            m.fs.costing.equipment_cost[opt] * pyunits.USD
            for opt in m.fs.continuous_opts_set
        )

    @m.fs.costing.Constraint(doc="Calculates the total financing cost of the plant[1].")
    def calculate_financing_cost_con(b):
        return (
            m.fs.costing.financing
            == m.fs.costing.financing_factor * m.fs.costing.total_plant_cost
        )

    @m.fs.costing.Constraint(doc="Calculates the 'other costs' of the plant[1].")
    def calculate_other_costs_con(b):
        return (
            m.fs.costing.other_costs
            == m.fs.costing.other_costs_factor * m.fs.costing.total_plant_cost
        )

    @m.fs.costing.Constraint(
        doc="Calculates the total overnight cost of the plant. "
        "Equal to the total plant cost + financing + 'other costs'[1]."
    )
    def calculate_total_overnight_cost_con(b):
        return (
            m.fs.costing.total_overnight_cost
            == m.fs.costing.total_plant_cost
            + m.fs.costing.financing
            + m.fs.costing.other_costs
        )


def add_operating_cost_cons(m):
    """
    This function builds the variable and fixed operating cost constraints.

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
    """

    @m.fs.costing.Constraint(
        m.fs.all_opts_set,
        m.fs.operational_range,
        doc="Calculates yearly operating expense for each option[2].",
    )
    def calculate_opt_yearly_variable_expense_cons(b, j, k, t):
        if (j, k) in m.fs.discrete_opts_set:
            return (
                m.fs.costing.opt_variable_operating_cost[j, k, t]
                == m.fs.costing.discrete_units_per_option[j, k]
                * m.fs.costing.yearly_cost_per_unit[j, k]
                * m.fs.option_binary_var[j, k]
            )
        else:
            return (
                m.fs.costing.opt_variable_operating_cost[j, k, t]
                == m.fs.costing.opt_var_oc_param_A[j, k]
                * sum(m.fs.f_in[j, k, c, t] for c in m.fs.tracked_comps)
                + m.fs.costing.opt_var_oc_param_B[j, k] * m.fs.option_binary_var[j, k]
            )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly variable operating cost[2].",
    )
    def calculate_total_yearly_variable_operating_costs_cons(b, t):
        return m.fs.costing.aggregate_variable_operating_cost[t] == sum(
            m.fs.costing.opt_variable_operating_cost[opt, t]
            for opt in m.fs.all_opts_set
        )

    @m.fs.costing.Constraint(
        m.fs.all_opts_set,
        doc="Calculate the number of operators needed for each option[3].",
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
                == m.fs.costing.num_operators[j, k] * m.fs.option_binary_var[j, k]
            )

    @m.fs.costing.Constraint(
        doc="Calculate the number of operators needed for the entire process[3]."
    )
    def calculate_total_operators_cons(b):
        return m.fs.costing.total_operators >= sum(
            m.fs.costing.operators_per_option[opt] for opt in m.fs.all_opts_set
        )

    @m.fs.costing.Constraint(doc="Calculates the cost of labor for the process.")
    def calculate_cost_of_labor_con(b):
        return (
            m.fs.costing.cost_of_labor
            == m.fs.costing.total_operators * m.fs.costing.labor_rate
        )

    @m.fs.costing.Constraint(doc="Calculate Maintenance & Supply Materials[2].")
    def calculate_m_and_sm_con(b):
        return (
            m.fs.costing.m_and_sm
            == m.fs.costing.m_and_sm_costing_factor * m.fs.costing.total_plant_cost
        )

    @m.fs.costing.Constraint(
        doc="Calculate Sample Analysis & Quality Assurance/Quality Control[2]."
    )
    def calculate_sa_and_qa_qc_con(b):
        return (
            m.fs.costing.sa_and_qa_qc
            == m.fs.costing.sa_and_qa_qc_costing_factor * m.fs.costing.cost_of_labor
        )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculate Sales, Intellectual Property, and Research & Development[2].",
    )
    def calculate_s_ip_r_and_d_con(b, t):
        return (
            m.fs.costing.s_ip_r_and_d[t]
            == m.fs.costing.s_ip_r_and_d_costing_factor * m.fs.costing.total_profit[t]
        )

    @m.fs.costing.Constraint(doc="Calculate Administrative & Supporting Labor[2].")
    def calculate_a_and_sl_con(b):
        return (
            m.fs.costing.a_and_sl
            == m.fs.costing.a_and_sl_costing_factor * m.fs.costing.cost_of_labor
        )

    @m.fs.costing.Constraint(doc="Calculate Fringe Benefits[2].")
    def calculate_fb_con(b):
        return (
            m.fs.costing.fb
            == m.fs.costing.fb_costing_factor * m.fs.costing.cost_of_labor
        )

    @m.fs.costing.Constraint(doc="Calculate Property Taxes & Insurance[2].")
    def calculate_pt_and_i_con(b):
        return (
            m.fs.costing.pt_and_i
            == m.fs.costing.pt_and_i_costing_factor * m.fs.costing.total_plant_cost
        )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total yearly fixed operating cost[2].",
    )
    def calculate_total_yearly_fixed_operating_costs_cons(b, t):
        return (
            m.fs.costing.aggregate_fixed_operating_cost[t]
            == m.fs.costing.cost_of_labor
            + m.fs.costing.m_and_sm
            + m.fs.costing.sa_and_qa_qc
            + m.fs.costing.s_ip_r_and_d[t]
            + m.fs.costing.a_and_sl
            + m.fs.costing.fb
            + m.fs.costing.pt_and_i
        )


def add_cash_flow_cons(m):
    """
    This function builds the cash flow constraints

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
    """

    @m.fs.costing.Constraint(
        m.fs.plant_life_range,
        doc="Calculates the total overnight cost expended in a given year[1].",
    )
    def calculate_total_overnight_cost_expended_cons(b, t):
        return (
            m.fs.costing.total_overnight_cost_expended[t]
            == m.fs.costing.total_overnight_capital_fraction_expended[t]
            * m.fs.costing.total_overnight_cost
        )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the plant overhead each year[2].",
    )
    def calculate_plant_overhead_cons(b, t):
        return m.fs.costing.plant_overhead[t] == m.fs.costing.plant_overhead_factor * (
            m.fs.costing.aggregate_variable_operating_cost[t]
            + m.fs.costing.aggregate_fixed_operating_cost[t]
        )

    @m.fs.costing.Constraint(
        m.fs.operational_range,
        doc="Calculates the total operating expenses each year. Includes variable, fixed, and plant overhead[2].",
    )
    def calculate_total_operating_expense_cons(b, t):
        return (
            m.fs.costing.total_operating_expense[t]
            == m.fs.costing.aggregate_variable_operating_cost[t]
            + m.fs.costing.aggregate_fixed_operating_cost[t]
            + m.fs.costing.plant_overhead[t]
        )

    @m.fs.costing.Constraint(
        m.fs.plant_life_range,
        doc="Calculate the cash flow for each year of the plant's lifetime[4].",
    )
    def calculate_cash_flows(b, t):
        start_year = pyo.value(m.fs.plant_start)
        # no profit or operational expenses in first year.
        if t == start_year:
            return (
                m.fs.costing.cash_flow[t]
                == -m.fs.costing.total_overnight_cost_expended[t]
            )
        else:
            return m.fs.costing.cash_flow[t] == (
                m.fs.costing.total_profit[t] - m.fs.costing.total_operating_expense[t]
            ) * (1 + m.fs.costing.i_operating_expense_escalation) ** (
                t - start_year
            ) - m.fs.costing.total_overnight_cost_expended[
                t
            ] * (
                1 + m.fs.costing.i_capital_escalation
            ) ** (
                t - start_year
            )


def add_costing_objective_functions(m, obj_func: str):
    """
    This function builds the costing objective function (either maximizing the NPV, or minimizing the COR).

    To calculate the capital cost parameters, the methodology presented in NETL's
    Quality Guidelines for Energy System Studies (QGESS)[1] was followed. However, since the methodology was created for costing
    electric power plants, modifications were made to match the economic assumptions made in Keim et al. (2019)[2]. To calculate
    the fixed and variable operating costs, the methodology followed by Kiem et al. (2019) was followed[2]. To calculate the number
    of operators needed, the methodology described by Ulrich et al. (2004) was followed[3]. The cash flows were calculated using
    the methodology described in Seider et al. (2017)[4].

    **References:**

    [1] Theis, J. *Cost Estimation Methodology for NETL Assessments of Power Plant Performance*; 2021.

    [2] Keim, S. A. *Production of Salable Rare Earths Products from Coal and Coal Byproducts in the U.S. Using Advanced
    Processes: Phase 1*; 2019.

    [3] Ulrich, G. D.; Vasudevan, P. T. *Chemical Engineering Process Design and Economics A Practical Guide*;
    Process Publishing, 2004.

    [4] Seider, W. D.; Lewin, D. R,; Seader, J.; Gani, S. W. R.; Ng, K. M. *Product and Process Design Principles*;
    Wiley, 2017.

    Args:
        m: pyomo model.
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
    """

    @m.fs.costing.Constraint(doc="Calculates the net present value of the plant[4].")
    def calculate_net_present_value_con(b):
        start_year = pyo.value(m.fs.plant_start)
        return m.fs.costing.net_present_value == sum(
            (m.fs.costing.cash_flow[t] * pyunits.year)
            / ((1 + m.fs.costing.discount_factor) ** (t - start_year))
            for t in m.fs.plant_life_range
        )

    # Check if NPV is objective function
    if obj_func == "NPV":
        # set objective function
        m.fs.costing.obj = pyo.Objective(
            expr=m.fs.costing.net_present_value,
            sense=pyo.maximize,
            doc="Objective function is maximizing the net present value when the NPV is chosen.",
        )

    # COR is objective function
    else:

        @m.fs.costing.Constraint(
            doc="Sets the net present value to zero when the cost of recovery is chosen as an "
            "objective function."
        )
        def set_net_present_value_to_zero_con(b):
            return m.fs.costing.net_present_value == 0

        # set objective function
        m.fs.costing.obj = pyo.Objective(
            expr=m.fs.costing.cost_of_recovery,
            sense=pyo.minimize,
            doc="Objective function is minimizing the cost of recovery when the COR is chosen.",
        )


def add_environmental_impact_params(
    m, consider_environmental_impacts, options_environmental_impacts, epsilon
):
    """
    This function builds the environmental impact parameters.

    When environmental impacts are considered, a multi-objective optimization problem is formulated in
    which some environmental metric, specified by the user, is minimized. A pareto front can be generated
    by varying the epsilon parameter value[1]. Although the global warming potential is the most common
    environmental metric investigated, many others also exist. Refer to the EPA's Tool for
    Reduction and Assessment of Chemicals and Other Environmental Impacts (TRACI)[2] for more information.

    **References:**

    [1] Mavrotas, G. Effective implementation of the epsilon-constraint method in multi-objective programming problems.
    *Applied Mathematics and Computation*, 213:455-465, 2009.

    [2] https://www.epa.gov/chemical-research/tool-reduction-and-assessment-chemicals-and-other-environmental-impacts-traci

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

    m.fs.environmental_impacts.options_environmental_impacts = pyo.Param(
        m.fs.all_opts_set,
        initialize=options_environmental_impacts,
        doc="Environmental impacts matrix. Units of chosen indicator per unit of incoming flowrate.",
    )
    m.fs.environmental_impacts.epsilon = pyo.Param(
        initialize=epsilon,
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
    @m.fs.environmental_impacts.Constraint(
        m.fs.all_opts_set,
        m.fs.operational_range,
        doc="Calculate the yearly environmental impacts generated by each option.",
    )
    def calculate_opt_yearly_impacts_con(b, j, k, t):
        return (
            m.fs.environmental_impacts.option_yearly_impacts[j, k, t]
            == sum(m.fs.f_in[j, k, c, t] for c in m.fs.tracked_comps)
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
            m.fs.environmental_impacts.total_impacts
            <= m.fs.environmental_impacts.epsilon
        )
