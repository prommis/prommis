#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Superstructure Code version 2
=============================

Author: Chris Laliwala
"""

import copy
import math
import sys
import warnings

import pyomo.environ as pyo


###################################################################################################
### Plant Lifetime Parameters
def check_plant_lifetime_params(plant_lifetime):
    """
    This function checks that the lifetime parameters are feasible.

    Args:
        plant_lifetime: (int) The total lifetime of the plant, including plant construction. Must be at least three years.
    """
    ## Check that plant lifetime is at least three years.
    if plant_lifetime < 3:
        raise ValueError("Plant lifetime must be a minimum of three years.")


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
    m.plant_lifetime_params = pyo.Block(doc="Block to hold plant lifetime parameters.")
    ## Pyomo parameters
    m.plant_lifetime_params.plant_start = pyo.Param(
        initialize=plant_start, doc="The year that plant construction begins."
    )
    m.plant_lifetime_params.plant_lifetime = pyo.Param(
        initialize=plant_lifetime,
        doc="The total lifetime of the plant, including plant construction. Must be at least three years.",
    )
    m.plant_lifetime_params.prod_start = pyo.Param(
        initialize=prod_start, doc="The first year of plant production."
    )
    m.plant_lifetime_params.plant_end = pyo.Param(
        initialize=plant_end, doc="The final year of plant production."
    )
    m.plant_lifetime_params.plant_life_range = pyo.RangeSet(
        plant_start, plant_end, doc="Lifetime of the plant."
    )
    m.plant_lifetime_params.operational_range = pyo.RangeSet(
        prod_start, plant_end, doc="Operational lifetime of the plant."
    )


###################################################################################################
### Feed Parameters
def check_feed_params(
    m, available_feed, collection_rate, tracked_comps, prod_comp_mass
):
    """
    The function checks that the feed parameter inputs are feasible.

    Args:
        m: pyomo model.
        available_feed: (dict) Total feedstock available for recycling each year.
        collection_rate: (float) How much available feed is processed by the plant each year.
        tracked_comps: (list) List of tracked components.
        prod_comp_mass: (dict) Mass of tracked components per EOL product.
    """
    ### Define parameters necessary for tests.
    # Define a set for the years in which the amount of available feed is defined.
    feed_years = set(available_feed.keys())
    # Define a set for the years in which the plant is in operation.
    operational_years = set(m.plant_lifetime_params.operational_range.data())
    # Define a set for the keys in prod_comp_mass.
    prod_comp_mass_keys = set(prod_comp_mass.keys())
    # Define a set for the tracked components.
    tracked_comps_set = set(tracked_comps)

    ### Define necessary data structures for tests.
    # Define a list to a list to keep track of the years in which the amount of available feed
    # passed by the user is negative.
    negative_years = []
    # Define a list to track the tracked components for which the amount of them within an EOL
    # product is defined to be negative.
    negative_prods = []
    # Define a list to track the tracked components for which the amount of them within an EOL
    # product is defined to be zero.
    zero_prods = []

    ### Run tests
    ## Check that available feed is provided for each year of plant operation.
    if feed_years != operational_years:
        raise ValueError(
            "Years of available_feed do not match the plant's operational period. "
            f"Expected years: {sorted(operational_years)}, "
            f"but got: {feed_years}"
        )
    ## Check that none of the available feeds passed are negative.
    if any(v < 0 for v in available_feed.values()):
        negative_years = [year for year, value in available_feed.items() if value < 0]
        raise ValueError(
            f"available_feed contains negative values for years: {negative_years}. "
            "Feedstock availability cannot be negative."
        )
    ## Check that available feed is not all zero.
    if not any(value != 0 for value in available_feed.values()):
        raise ValueError(
            "All values in available_feed are zero. At least one year must have non-zero feedstock available."
        )
    ## Check that collection rate is positive value.
    if collection_rate <= 0:
        raise ValueError("Collection rate must be a positive value.")
    ## Check that at least one tracked component is specified.
    if not tracked_comps:
        raise ValueError(
            "tracked_comps list is empty. At least one component must be tracked."
        )
    ## Check an amount contained within an EOL product is specified for each tracked component.
    if prod_comp_mass_keys != tracked_comps_set:
        raise ValueError(
            f"prod_comp_mass keys don't match up with the set of tracked components."
        )
    ## Check that amounts contained within EOL product for each tracked component is non-negative.
    if any(v < 0 for v in prod_comp_mass.values()):
        negative_prods = [prod for prod, amount in prod_comp_mass.items() if amount < 0]
        raise ValueError(
            f"prod_comp_mass contains negative values for the tracked components: {negative_prods}. "
            "Amounts cannot be negative"
        )
    ## Raise warning if amounts contained within EOL product for a tracked component is zero.
    if any(value == 0 for value in prod_comp_mass.values()):
        zero_prods = [prod for prod, amount in prod_comp_mass.items() if amount == 0]
        warnings.warn(
            f"prod_comp_mass contains zero values for the tracked components: {zero_prods}. "
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
    ## Create block
    m.feed_params = pyo.Block(doc="Block to hold feed parameters.")
    ## Pyomo parameters
    m.feed_params.available_feed = pyo.Param(
        m.plant_lifetime_params.operational_range,
        initialize=available_feed,
        doc="The total feedstock available for recycling each year.",
    )
    m.feed_params.collection_rate = pyo.Param(
        initialize=collection_rate,
        doc="The fraction of available feed that is processed by the plant each year.",
    )
    m.feed_params.tracked_comps = pyo.Set(
        initialize=tracked_comps, doc="Tracked components."
    )
    m.feed_params.prod_comp_mass = pyo.Param(
        m.feed_params.tracked_comps,
        initialize=prod_comp_mass,
        doc="The mass of each tracked component per EOL product.",
    )
    m.feed_params.feed_entering = pyo.Param(
        m.plant_lifetime_params.operational_range,
        initialize=feed_entering,
        doc="The amount of feed entering the plant each year.",
    )
    m.feed_params.max_feed_entering = pyo.Param(
        initialize=max_feed_entering,
        doc="The max yearly feed that enters the plant over the production period.",
    )
    m.feed_params.max_feed_entering_year = pyo.Param(
        initialize=max_feed_entering_year,
        doc="The year that the max feed enters the plant.",
    )


###################################################################################################
### Superstructure Formulation Parameters
def check_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
):
    """
    The function checks that the superstructure formulation parameter inputs are feasible.

    Args:
        m: pyomo model.
        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        discrete_opts: (list) List of options that utilize discrete units.
        option_eff: (dict) Tracked component retention efficiency for each option.
    """
    ### Define parameters necessary for tests.
    # Define a set of the number of stages in the superstructure.
    num_stages_set = set(pyo.RangeSet(num_stages).data())
    # Define a set for the keys of the options_in_stage dict.
    options_in_stage_keys_set = set(options_in_stage.keys())
    # Define a set for all the options in the superstructure.
    all_opts_set = {
        (j, k)
        for j in range(1, num_stages + 1)
        for k in range(1, options_in_stage[j] + 1)
    }
    discr_opts_set = set(discrete_opts)
    # Define a set for all the tracked components.
    tracked_comps = set(m.feed_params.tracked_comps.data())
    # Define a set containing the keys of the opt_var_oc_params dict.

    ### Define necessary data structures for tests.
    # Define a list to track the options in stage j which are not connected to an option in
    # stage j+1.
    missing_values = []
    # Define a list to track the options which don't have a connection from the previous
    # stage.
    disconnected_options = []
    # Define a list to track the missing efficiencies in the option_eff dict.
    missing_effs = []
    # Define a list to track the efficiencies defined as negative in the option_eff dict.
    negative_effs = []

    ### Run tests
    ## Check that there is at least 2 stages.
    if num_stages < 2:
        raise ValueError("There must be at least 2 processing stages.")
    ## Check that the stages must be numbered starting at 1 and counting up.
    if num_stages_set != options_in_stage_keys_set:
        raise ValueError(
            "Stages must start at 1 and count up. Each stage must contain at least 1 option. options_in_stage does not follow this convention. "
            f"Expected keys: {num_stages_set}, "
            f"but got: {options_in_stage_keys_set}"
        )
    ## Check that connections between options in superstructure are feasible.
    # Check that each option in stage j is connected to an option in stage j+1.
    missing_values = [key for key, value in option_outlets.items() if value is None]
    if missing_values:
        raise ValueError(
            f"Options {missing_values} are missing connections in the next stage, as defined by option_outlets."
        )
    # Check that each option in stage j is connected to an option in the preceding stage, stage j-1.
    # Iterate over stages starting from 2 (since stage 1 has no predecessors).
    for current_stage in range(2, num_stages + 1):
        # Track the previous stage.
        previous_stage = current_stage - 1
        # Get number of options in current stage.
        num_current_options = options_in_stage[current_stage]
        # Get number of options in previous.
        num_previous_options = options_in_stage[previous_stage]
        # Check each option in the current stage.
        for current_option in range(1, num_current_options + 1):
            # Define a boolean to track if current_option has any connections to options from previous stage.
            connected = False
            # Check if option in current stage is connected to any options in the previous stage.
            for previous_option in range(1, num_previous_options + 1):
                # Get outlets for option: (previous_stage, previous_option).
                outlets = option_outlets.get((previous_stage, previous_option), [])
                # If the current option is connected to an option from the previous stage, break from loop.
                if current_option in outlets:
                    connected = True
                    break  # No need to check further
            # If the option is not connected to any of the options in the previous option, track it. This is an error.
            if not connected:
                disconnected_options.append((current_stage, current_option))
    # Raise error if there are any options that are missing connections from the previous stage.
    if disconnected_options:
        error_msg = "The following options are not connected from the previous stage:\n"
        error_msg += "\n".join(
            f"  - Stage {stage}, Option {option}"
            for stage, option in disconnected_options
        )
        raise ValueError(error_msg)
    ## Check that all discrete options listed are feasible.
    # This means they must be contained within the superstructure.
    if not (discr_opts_set <= all_opts_set):
        raise ValueError("Discrete options listed are not feasible.")
    ## Check that an option efficiency is defined for each option and is nonnegative.
    for j in range(1, num_stages + 1):
        for k in range(1, options_in_stage[j] + 1):
            # Check for missing efficiencies.
            option_eff_tracked_comp_keys = set(option_eff[(j, k)].keys())
            if option_eff_tracked_comp_keys != tracked_comps:
                # Keep track of the components for which an efficiency is not defined.
                missing = tracked_comps - option_eff_tracked_comp_keys
                missing_effs.append((j, k, missing))
            # Check for negative efficiencies.
            # Keep track of the components for which a negative efficiency is defined.
            negative_comps = [c for c, eff in option_eff[j, k].items() if eff < 0]
            if negative_comps:
                negative_effs.append((j, k, negative_comps))
    # Raise an error if there are missing efficiencies.
    if missing_effs:
        msg = "Efficiencies not specified for all tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing components: {missing}"
            for j, k, missing in missing_effs
        )
        raise ValueError(msg)
    # Raise an error if there are negative efficiencies.
    if negative_effs:
        msg = "Negative efficiencies specified for some tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) negative efficiencies specified for: {negative}"
            for j, k, negative in negative_effs
        )
        raise ValueError(msg)


def add_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
):
    """
    This function builds the rest of the superstructure formulation parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model.
        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        discrete_opts: (list) List of options that utilize discrete units.
        option_eff: (dict) Tracked component retention efficiency for each option.
    """
    ### Define parameters from user input.
    # Define a parameter for the max number of options in any of the stages.
    max_options = max(options_in_stage.values())
    # Define a set of all the discrete options.
    discrete_opts_set = set(discrete_opts)
    # Define a set containing all the options in the superstructure.
    all_opts_set = set(
        (j, k)
        for j in range(1, num_stages + 1)
        for k in range(1, options_in_stage[j] + 1)
    )
    # Define a set containing all the continuous options in the superstructure.
    continuous_opts_set = all_opts_set - discrete_opts_set
    # Define a set containing all the  options in the final stage.
    final_opts_list = [
        (num_stages, k) for k in range(1, options_in_stage[num_stages] + 1)
    ]

    ### Define functions needed to initialize pyomo parameters.
    # Define a function for initializing option efficiency pyomo parameter.
    def option_eff_initialize(m, j, k, c):
        return option_eff[(j, k)][c]

    ### Define necessary pyomo parameters.
    ## Create block
    m.supe_form_params = pyo.Block(
        doc="Block to hold superstructure formulation parameters."
    )
    ## Pyomo parameters
    m.supe_form_params.num_stages = pyo.Param(
        initialize=num_stages, doc="The total number of stages in the superstructure."
    )
    m.supe_form_params.stages_set = pyo.RangeSet(
        1, num_stages, doc="A set of all the stages in the superstructure."
    )
    m.supe_form_params.options_in_stage = pyo.Param(
        m.supe_form_params.stages_set,
        initialize=options_in_stage,
        doc="The number of options in each stage.",
    )
    m.supe_form_params.max_options = pyo.Param(
        initialize=max_options, doc="The max number of options in any of the stages."
    )
    m.supe_form_params.max_options_set = pyo.RangeSet(
        1, max_options, doc="Set containing max number of options in any of the stages."
    )
    m.supe_form_params.all_opts_set = pyo.Set(
        initialize=(
            (j, k)
            for j in m.supe_form_params.stages_set
            for k in range(1, options_in_stage[j] + 1)
        ),
        doc="Set containing all options in the superstructure.",
    )
    m.supe_form_params.discrete_opts_set = pyo.Set(
        initialize=((opt) for opt in discrete_opts_set),
        doc="Set containing all the options which utilize discrete units (discrete options).",
    )
    m.supe_form_params.continuous_opts_set = pyo.Set(
        initialize=((opt) for opt in continuous_opts_set),
        doc="Set containing all the options in the stage which don't utilize discrete units (continuous options).",
    )
    m.supe_form_params.option_outlets = pyo.Param(
        m.supe_form_params.all_opts_set,
        initialize=option_outlets,
        doc="Defines the set of options k' in stage j+1 connected to option k in stage j.",
    )
    m.supe_form_params.option_eff = pyo.Param(
        m.supe_form_params.all_opts_set,
        m.feed_params.tracked_comps,
        initialize=option_eff_initialize,
        doc="Defines the tracked component efficiencies for each option in the superstructure. "
        "Efficiency defined for each tracked component 'c' for each option 'k' in each stage 'j'",
    )
    m.supe_form_params.final_opts_set = pyo.Set(
        initialize=final_opts_list,
        doc="Set containing all of the options in the final stage.",
    )


###################################################################################################
### Operating Parameters
def check_operating_params(
    m,
    profit,
    opt_var_oc_params,
    workers_per_discr_unit,
    yearly_cost_per_unit,
    capital_cost_per_unit,
    processing_rate,
    num_operators,
    labor_rate,
):
    """
    This function checks that all the operating parameters are feasible.

    Args:
        m: pyomo model.
        profit: (dict) Profit per unit of product in terms of tracked components.
        opt_var_oc_params: (dict) Holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        workers_per_discr_unit: (dict) Number of workers needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of units of incoming feed processed per year per unit.
        num_operators: (dict) Number of operators needed for each option.
        labor_rate: (float) Yearly wage per operator.
    """
    ### Define parameters necessary for tests.
    # Define a set of all the keys in profit dict.
    profit_opt_keys = set(profit.keys())
    # Define a set containing all the options in the final stage.
    final_opts_set = set(m.supe_form_params.final_opts_set)
    # Define a set of all the tracked components.
    tracked_comps = set(m.feed_params.tracked_comps.data())
    # Define a set containing the keys in the opt_var_oc_params dict.
    opt_var_oc_params_keys = set(opt_var_oc_params.keys())
    # Define a set containing all the discrete options.
    discr_opts_set = set(m.supe_form_params.discrete_opts_set.data())
    # Define a set containing all the continuous options.
    continuous_opts_set = set(m.supe_form_params.continuous_opts_set.data())
    # Define a set containing the necessary variable operating parameters for all continuous options.
    var_oc_params = set(["a", "b"])
    # Define a set containing the keys of the workers_per_discr_unit dict.
    workers_per_discr_unit_keys = set(workers_per_discr_unit.keys())
    # Define a set containing the keys from the yearly_cost_per_unit dict.
    yearly_cost_per_unit_keys = set(yearly_cost_per_unit.keys())
    # Define a set containing all the keys in the capital_cost_per_unit dict.
    capital_cost_per_unit_set = set(capital_cost_per_unit.keys())
    # Define a set containing all the keys from the processing_rate dict.
    processing_rate_keys = set(processing_rate.keys())
    # Define a set containing all the keys from the num_operators dict.
    num_operators_keys = set(num_operators.keys())
    ### Define necessary data structures for tests.
    # Define a list to track the final options which don't define profits for all tracked components.
    missing_profit_comps = []
    # Define a list to track the final options which define negative profits for all tracked components.
    negative_profit_comps = []
    # Define a list to track the options for which all the necessary variable operating cost parameters
    # are not defined.
    missing_var_oc_params = []
    # Define a list for tracking the options which define negative workers per discrete unit.
    negative_workers_per_discr_units = []
    # Define a list for tracking the options which define a negative yearly cost per unit
    negative_yearly_cost_per_unit = []
    # Define a list for tracking the discrete options which define a negative capital cost per unit.
    negative_capital_cost_per_unit = []
    # Define a list for tracking the discrete options which define a non-positive processing rate for a unit.
    nonpositive_processing_rate = []
    # Define a list for tracking the continuous options which define a negative number of operators.
    negative_num_operators = []

    ### Run tests
    ## Check that profit per product is defined for all options in the final stage
    if profit_opt_keys != final_opts_set:
        raise ValueError(
            "Must include profit per unit of product for all options in the final stage."
        )
    ## Check that profit per product is in terms of the tracked components for all options in the final stage and are all nonnegative.
    for opt in m.supe_form_params.final_opts_set:
        # Check for missing profits.
        profit_tracked_comps = set(profit[opt].keys())
        if profit_tracked_comps != tracked_comps:
            # Keep track of the tracked components for which no profit per component is defined.
            missing = tracked_comps - profit_tracked_comps
            missing_profit_comps.append(opt + (missing,))
        # Check for negative profits.
        # Keep track of the tracked components for which a negative profit per component is defined.
        negative_comps = [c for c, profit in profit[opt].items() if profit < 0]
        if negative_comps:
            negative_profit_comps.append(opt + (negative_comps,))
    # Raise error if there are options missing profits per tracked component.
    if missing_profit_comps:
        msg = "Profits not specified for all tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing components: {missing}"
            for j, k, missing in missing_profit_comps
        )
        raise ValueError(msg)
    # Raise error if there are options with negative profits per tracked component defined.
    if negative_profit_comps:
        msg = "Profits some tracked components are listed as negative in the following options\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) negative components: {negative}"
            for j, k, negative in negative_profit_comps
        )
        raise ValueError(msg)
    ## Check that variable operating cost params are defined for all continuous options.
    if opt_var_oc_params_keys != continuous_opts_set:
        raise ValueError(
            "Variable operating cost params not defined for all continuous options."
        )
    ## Check that both necessary variable operating cost parameters ('a' and 'b') are defined for all continuous options.
    for opt in continuous_opts_set:
        # Check for missing variable operating cost parameters.
        params = set(opt_var_oc_params[opt].keys())
        if params != var_oc_params:
            # Keep tracking of the missing variable operating cost parameters for the option.
            missing = var_oc_params - params
            missing_var_oc_params.append(opt + (missing,))
    # Raise an error if there are continuous options missing variable operating cost parameters.
    if missing_var_oc_params:
        msg = "not all variable operating cost parameters defined in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing parameters: {missing}"
            for j, k, missing in missing_var_oc_params
        )
        raise ValueError(msg)
    ## Check that workers per discrete unit is defined for all options that utilize discrete units.
    if workers_per_discr_unit_keys != discr_opts_set:
        raise ValueError("workers_per_discr_unit not defined for all discrete options.")
    ## Check that workers per discrete unit are all defined to be non-negative.
    negative_workers_per_discr_units = [
        opt for opt, w in workers_per_discr_unit.items() if w < 0
    ]
    # If there are negative workers defined per discrete unit for any discrete options, raise an error.
    if negative_workers_per_discr_units:
        raise ValueError("Workers per discrete unit must all be non-negative.")
    ## Check that yearly cost per unit is defined for all discrete options.
    if yearly_cost_per_unit_keys != discr_opts_set:
        raise ValueError(
            "yearly_cost_per_unit must be defined for all discrete options."
        )
    ## Check that the yearly cost per unit values are all non-negative.
    negative_yearly_cost_per_unit = [
        opt for opt, cost in yearly_cost_per_unit.items() if cost < 0
    ]
    # If there any discrete options that define a negative yearly cost per unit, raise an error.
    if negative_yearly_cost_per_unit:
        raise ValueError("yearly_cost_per_unit values must all be non-negative.")
    ## Check that cost per unit is defined for all discrete options.
    if capital_cost_per_unit_set != discr_opts_set:
        raise ValueError(
            "capital_cost_per_unit must be defined for all discrete options."
        )
    ## Check that the cost per unit values are all non-negative.
    negative_capital_cost_per_unit = [
        opt for opt, cost in capital_cost_per_unit.items() if cost < 0
    ]
    # If there are any discrete options which define a negative capital cost per unit, raise an error.
    if negative_capital_cost_per_unit:
        raise ValueError("capital_cost_per_unit values must all be non-negative.")
    ## Check that prosessing rate is defined for all discrete options.
    if processing_rate_keys != discr_opts_set:
        raise ValueError("processing_rate must be defined for all discrete options.")
    ## Check that all processing rates are positive.
    nonpositive_processing_rate = [
        opt for opt, rate in processing_rate.items() if rate <= 0
    ]
    # If there are any discrete options which define a non-positive processing rate for a unit, raise an error.
    if nonpositive_processing_rate:
        raise ValueError("processing rates must all be positive.")
    ## Check that num_operators is defined for all continuous options
    if num_operators_keys != continuous_opts_set:
        raise ValueError(
            "num_operators must be defined for all continuous options, and must not be defined for discrete options."
        )
    ## Check that num_operators values are all non-negative
    negative_num_operators = [
        opt for opt, workers in num_operators.items() if workers < 0
    ]
    if negative_num_operators:
        raise ValueError("number of workers must all be non-negative.")
    ## Check that labor rate is non-negative
    if labor_rate < 0:
        raise ValueError("labor rate must be non-negative.")


def add_operating_params(
    m,
    profit,
    opt_var_oc_params,
    workers_per_discr_unit,
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
        workers_per_discr_unit: (dict) Number of workers needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of units of incoming feed processed per year per unit.
        num_operators: (dict) Number of operators needed for each option.
        labor_rate: (float) Yearly wage per operator.
    """

    ### Define functions needed to initialize pyomo parameters.
    # Define a function for initializing profit pyomo parameter.
    def profit_initialization(m, j, k, c):
        return profit[(j, k)][c]

    # Define a function for initializing opt_var_oc_params pyomo parameter.
    def opt_var_oc_params_initialization(m, j, k, var_oc_param):
        return opt_var_oc_params[(j, k)][var_oc_param]

    ### Define necessary pyomo parameters.
    ## Create block
    m.operating_params = pyo.Block(doc="Block to hold operating parameters.")
    ## Pyomo parameters
    m.operating_params.profit = pyo.Param(
        m.supe_form_params.final_opts_set,
        m.feed_params.tracked_comps,
        initialize=profit_initialization,
        doc="Holds the profit for all options in the final stage in terms of the tracked components.",
    )
    m.operating_params.var_oc_params_set = pyo.Set(
        initialize=["a", "b"],
        doc="Set containing the necessary parameters for calculating the variable operating costs for all continuous options.",
    )
    m.operating_params.opt_var_oc_params = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        m.operating_params.var_oc_params_set,
        initialize=opt_var_oc_params_initialization,
        doc="Holds all the variable operating costs parameter values for all continuous options.",
    )
    m.operating_params.workers_per_discr_unit = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=workers_per_discr_unit,
        doc="The number of operators needed per discrete unit for discrete options.",
    )
    m.operating_params.yearly_cost_per_unit = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=yearly_cost_per_unit,
        doc="The operating costs per discrete unit for discrete options.",
    )
    m.operating_params.capital_cost_per_unit = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=capital_cost_per_unit,
        doc="The capital cost per discrete unit for discrete options.",
    )
    m.operating_params.processing_rate = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=processing_rate,
        doc="The processing rate per discrete unit for discrete options.",
    )
    m.operating_params.num_operators = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=num_operators,
        doc="The number of operators per continuous options.",
    )
    m.operating_params.labor_rate = pyo.Param(
        initialize=labor_rate, doc="The yearly wage per operator."
    )
    # This is a hardcoded value for now. Setting to 100, as having this many units would likely not be reasonable anyways.
    m.operating_params.max_dis_workers = pyo.Param(
        initialize=100, doc="The max number of disassembly units possible."
    )
    # This is a hardcoded value for now. Setting to 100, as having this many operators would likely not be reasonable anyways.
    m.operating_params.max_operators = pyo.Param(
        initialize=100,
        doc="The max number of operators possible for the entire process.",
    )


###################################################################################################
### Costing Parameters
def check_costing_params(m, discretized_capex):
    """
    This function checks that all the costing parameters are feasible.

    Args:
        m: pyomo model.
        discretized_capex: (dict) Discretized cost by flows entering for each continuous option.
    """
    ### Define parameters necessary for tests
    # Define a set of all the keys in the discretized_capex dict.
    discretized_capex_opts_set = set(discretized_capex.keys())

    ### Define necessary data structures for tests.
    # Define a list of continuous options which are missing discretized capex data.
    missing_continuous_opts = []
    # Define a list for tracking the discrete options for which discretized capex data is defined for.
    discrete_opts = []
    # Define a list for tracking the options in which inconsistent discretized capex is define.
    inconsistent_data_point_opts = []

    ### Run tests
    ## Check that discretized capex provided for all continuous options
    # and check that discretized capex not provided for options that utilize discrete units.
    for opt in m.supe_form_params.continuous_opts_set:
        # Keep track of the continuous opts for which discretized data is not defined for.
        if opt not in discretized_capex_opts_set:
            missing_continuous_opts.append(opt)
    # Check that discretized capex not provided for discrete opts.
    for opt in m.supe_form_params.discrete_opts_set:
        # Keep track of the discrete options for which discretized capex data is defined for.
        if opt in discretized_capex_opts_set:
            discrete_opts.append(opt)
    # Raise error if discretized capex is not provided for all continuous options.
    if missing_continuous_opts:
        raise ValueError(
            f"discretized_capex is missing values for the following continuous options: {missing_continuous_opts}. "
        )
    # Raise error if discretized capex provided for any discrete options.
    if discrete_opts:
        raise ValueError(
            f"discretized_capex contains values for the following discrete options: {discrete_opts}. "
        )
    ## Check that all options have the same number of discretized data points for flows entering and costs.
    for opt in m.supe_form_params.continuous_opts_set:
        # Keep track of the number of flowrate data points defined for the option.
        opt_num_flow_data_points = len(discretized_capex[opt]["Flowrates"])
        # Keep track of the number of cost data points defined for the option.
        opt_num_cost_data_points = len(discretized_capex[opt]["Costs"])
        # Keep track of the options for which the number of flowrate and cost data points are not the same.
        if opt_num_flow_data_points != opt_num_cost_data_points:
            inconsistent_data_point_opts.append(opt)
    # Raise error if there are options with an inconsistent number of data points within discretized_capex.
    if inconsistent_data_point_opts:
        raise ValueError(
            f"Inconsistent number of data points for Flowrates and Costs within discretized_capex for the following options: {inconsistent_data_point_opts}. "
        )


def add_costing_params(m, discretized_capex):
    """
    This function adds all the costing parameters to a block.

    Args:
        m: pyomo model.
        discretized_capex: (dict) Discretized cost by flows entering for each continuous option
    """
    ### Define parameters from user input.
    # Define a dict to hold discretized flowrate data for each option.
    flowrates_data = {}
    for opt in m.supe_form_params.continuous_opts_set:
        flowrates_data[opt] = discretized_capex[opt]["Flowrates"]
    # Define a dict to hold discretized cost data for each option.
    costs_data = {}
    for opt in m.supe_form_params.continuous_opts_set:
        costs_data[opt] = discretized_capex[opt]["Costs"]

    ### Define necessary pyomo parameters.
    ## Create a block
    m.costing_params = pyo.Block(doc="Block to hold costing parameters.")
    ## Pyomo parameters
    m.costing_params.flowrates_data = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=flowrates_data,
        doc="Discretized flowrate data for all continuous options.",
    )
    m.costing_params.costs_data = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=costs_data,
        doc="Discretized costing data for all continuous options.",
    )


###################################################################################################
### Objective Function Parameters
def check_objective_function_params(m, obj_func):
    """
    This function checks that all the objective function params are feasible.

    Args:
        m: pyomo model.
        obj_func (str) Choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive.
    """
    ### Run tests
    if (obj_func != "NPV") and (obj_func != "COR"):
        raise ValueError(
            "Invalid choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive."
        )


def add_objective_function_params(m, obj_func):
    """
    This function adds all the objective function parameters to a block.

    Args:
        m: pyomo model.
        obj_func (str) Choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive.
    """
    ### Define necessary pyomo parameters.
    ## Create block
    m.obj_func_params = pyo.Block("Block to hold objective function parameters.")
    ## Pyomo parameters
    m.obj_func_params.obj = pyo.Param(
        initialize=obj_func, doc="Choice of objective function."
    )


def add_mass_balance_params(m):
    """
    This function adds all mass balance parameters to a block.

    Args:
        m: pyomo model.
    """
    ### Define necessary pyomo parameters.
    ## Create block
    m.mb_params = pyo.Block(doc="Block to hold all mass balance parameters.")
    ## Pyomo parameters
    # Flow entering each stage (except the last stage).
    m.mb_params.f_stages = pyo.RangeSet(
        1,
        m.supe_form_params.num_stages - 1,
        doc="Set of all stages except the last. Used to define the flow variable: 'f'.",
    )
    m.mb_params.flow_set = pyo.Set(
        initialize=m.supe_form_params.all_opts_set
        * m.feed_params.tracked_comps
        * m.plant_lifetime_params.operational_range,
        doc="Set of all options, tracked components, and operational years of the plant. Used to define the flow variables: 'f_in' and 'f_out'.",
    )


def add_mass_balance_vars(m):
    """
    This function adds all mass balance variables to a block.

    Args:
        m: pyomo model.
    """
    ### Define necessary pyomo variables.
    ## Create block
    m.mb_vars = pyo.Block(doc="Block to hold all mass balance variables.")
    ## Pyomo variables
    m.mb_vars.f = pyo.Var(
        m.mb_params.f_stages
        * m.feed_params.tracked_comps
        * m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each stage (except the last stage). See documentation for more details.",
    )
    m.mb_vars.f_in = pyo.Var(
        m.mb_params.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )
    m.mb_vars.f_out = pyo.Var(
        m.mb_params.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )


def add_mass_balance_cons(m):
    """
    This function adds all mass balance constraints to a block.

    Args:
        m: pyomo model.
    """
    ### Define necessary pyomo constraints
    ## Create block
    m.mb_cons = pyo.Block(doc="Block to hold all mass balance constraints.")

    ## Pyomo constraints
    @m.mb_cons.Constraint(
        m.supe_form_params.stages_set,
        m.feed_params.tracked_comps,
        m.plant_lifetime_params.operational_range,
        doc="Equation (1) from the documentation.",
    )
    def inlet_flow_cons(b, j, c, t):
        m = b.model()  # Get the main model
        if j == 1:
            return pyo.Constraint.Skip
        else:
            # Extract all the options available in stage 'j'
            num_options = range(1, m.supe_form_params.options_in_stage[j] + 1)
            return m.mb_vars.f[j - 1, c, t] == sum(
                m.mb_vars.f_in[j, k, c, t] for k in num_options
            )

    @m.mb_cons.Constraint(
        m.supe_form_params.stages_set,
        m.feed_params.tracked_comps,
        m.plant_lifetime_params.operational_range,
        doc="Equation (2) from the documentation.",
    )
    def init_flow_cons(b, j, c, t):
        m = b.model()  # Get the main model
        if j == 1:
            # Extract all the options available in stage 'j'
            num_options = range(1, m.supe_form_params.options_in_stage[j] + 1)
            return m.feed_params.feed_entering[t] * m.feed_params.prod_comp_mass[
                c
            ] == sum(m.mb_vars.f_in[j, k, c, t] for k in num_options)
        else:
            return pyo.Constraint.Skip

    @m.mb_cons.Constraint(
        m.mb_params.flow_set, doc="Equation (3) from the documentation."
    )
    def intermediate_flow_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        alpha = m.supe_form_params.option_eff[j, k, c]
        return m.mb_vars.f_in[j, k, c, t] * alpha == m.mb_vars.f_out[j, k, c, t]

    @m.mb_cons.Constraint(
        m.supe_form_params.stages_set,
        m.feed_params.tracked_comps,
        m.plant_lifetime_params.operational_range,
        doc="Equation (4) from the documentation.",
    )
    def outlet_flow_cons(b, j, c, t):
        m = b.model()  # Get the main model
        if j != m.supe_form_params.num_stages:
            # Extract all the options available in stage 'j'
            num_options = range(1, m.supe_form_params.options_in_stage[j] + 1)
            return (
                sum(m.mb_vars.f_out[j, k, c, t] for k in num_options)
                == m.mb_vars.f[j, c, t]
            )
        else:
            return pyo.Constraint.Skip


def add_logic_params(m):
    """
    This function adds all logic parameters to a block.

    Args:
        m: pyomo model.
    """
    ### Define parameters
    # Calculate Big-M values for each tracked component.
    m_val = {}
    for c in m.feed_params.tracked_comps:
        # Max value assumes 100% efficiency (no losses) over the plant's lifetime. Rounded up to ensure validity
        # in constraints.
        m_val[c] = math.ceil(
            m.feed_params.max_feed_entering * m.feed_params.prod_comp_mass[c]
        )

    ### Define necessary pyomo parameters.
    ## Create block
    m.logic_params = pyo.Block(doc="Block to hold logic parameters.")
    ## Pyomo parameters
    m.logic_params.big_m_val = pyo.Param(
        m.feed_params.tracked_comps,
        initialize=m_val,
        doc="Big-M parameters used in Equations (7) and (8) from the documentation.",
    )


def add_logic_vars(m):
    """
    This function adds all logic variables to a block.

    Args:
        m: pyomo model.
    """
    ### Define necessary pyomo variables.
    ## Create block
    m.logic_vars = pyo.Block(doc="Block to hold logic variables.")
    ## Pyomo variables
    m.logic_vars.option_binary_var = pyo.Var(
        m.supe_form_params.all_opts_set,
        domain=pyo.Binary,
        doc="Binary variables to indicate whether or not an option has been selected.",
    )


def add_logic_cons(m):
    """
    This function adds all logic constraints to a block.

    Args:
        m: pyomo model.
    """
    ### Define necessary pyomo constraints
    ## Create block
    m.logic_cons = pyo.Block(doc="Block to hold logic constraints.")

    ## Pyomo constraints
    @m.logic_cons.Constraint(
        m.supe_form_params.stages_set, doc="Equation (5) from the documentation."
    )
    def stage_binary_cons(b, j):
        m = b.model()  # Get the main model
        # Extract all the options available in stage 'j'
        num_options = range(1, m.supe_form_params.options_in_stage[j] + 1)
        return sum(m.logic_vars.option_binary_var[j, k] for k in num_options) == 1

    @m.logic_cons.Constraint(
        m.supe_form_params.all_opts_set, doc="Equation (6) from the documentation."
    )
    def connection_binary_cons(b, j, k):
        m = b.model()  # Get the main model
        if j != m.supe_form_params.num_stages:
            # Extract the set of options k' in stage j+1 connected to option k in stage j.
            opt_connections = m.supe_form_params.option_outlets[j, k]
            return m.logic_vars.option_binary_var[j, k] <= sum(
                m.logic_vars.option_binary_var[j + 1, kp] for kp in opt_connections
            )
        else:
            return pyo.Constraint.Skip

    # create constraints
    @m.logic_cons.Constraint(
        m.mb_params.flow_set, doc="Equation (7) from the documentation."
    )
    def f_in_big_m_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        return (
            m.mb_vars.f_in[j, k, c, t]
            <= m.logic_vars.option_binary_var[j, k] * m.logic_params.big_m_val[c]
        )

    @m.logic_cons.Constraint(
        m.mb_params.flow_set, doc="Equation (8) from the documentation."
    )
    def f_out_big_m_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        return (
            m.mb_vars.f_out[j, k, c, t]
            <= m.logic_vars.option_binary_var[j, k] * m.logic_params.big_m_val[c]
        )


def build_model(
    ###################################################################################################
    ### Plant Lifetime Parameters
    # The year that plant construction begins.
    plant_start: int,
    # The total lifetime of the plant, including plant construction. Must be at least three years.
    plant_lifetime: int,
    ###################################################################################################
    ###################################################################################################
    ### Feed parameters
    # Total feedstock available for recycling each year.
    available_feed: dict,
    # How much available feed is processed by the plant each year.
    collection_rate: float,
    # List of tracked components.
    tracked_comps: list,
    # Mass of tracked component per EOL product.
    prod_comp_mass: dict,
    ###################################################################################################
    ###################################################################################################
    ### Superstructure formulation parameters
    # Number of total stages.
    num_stages: int,
    # Number of options in each stage.
    options_in_stage: dict,
    # Set of options k' in stage j+1 connected to option k in stage j.
    option_outlets: dict,
    # List of options that utilize discrete units.
    discrete_opts: list,
    # Tracked component retention efficiency for each option.
    option_eff: dict,
    ###################################################################################################
    ###################################################################################################
    ### Operating Parameters
    # Profit per unit of product in terms of tracked components.
    profit: dict,
    # Holds the variable operating cost param for options that are continuous.
    # Variable operating costs assumed to be proportional to the feed entering the option.
    opt_var_oc_params: dict,
    # Number of workers needed per discrete unit for options that utilize discrete units.
    workers_per_discr_unit: dict,
    # Yearly operating costs per unit for options which utilize discrete units.
    yearly_cost_per_unit: dict,
    # Cost per unit for options which utilize discrete units.
    capital_cost_per_unit: dict,
    # Processing rate per unit for options that utilize discrete units.
    # In terms of units of incoming feed processed per year per unit.
    processing_rate: dict,
    # Number of operators needed for each option.
    num_operators: dict,
    # Yearly wage per operator.
    labor_rate: float,
    ###################################################################################################
    ###################################################################################################
    ### Costing Parameters
    # Define Python Dictionary with discretized cost by flows for each option.
    discretized_capex: dict,
    ###################################################################################################
    ###################################################################################################
    ### Objective Function Parameters
    # Choice of objective function. Options are 'NPV' or 'COR'. Selection is case-sensitive.
    obj_func: str,
    ###################################################################################################
):
    # Define model
    m = pyo.ConcreteModel()

    ### Plant lifetime parameters
    # Check that plant lifetime parameters are feasible.
    check_plant_lifetime_params(plant_lifetime)
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
        m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
    )
    # Create separate block to hold superstructure formulation parameters.
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
    )

    ### Operating parameters
    # Check that operating parameters are feasible.
    check_operating_params(
        m,
        profit,
        opt_var_oc_params,
        workers_per_discr_unit,
        yearly_cost_per_unit,
        capital_cost_per_unit,
        processing_rate,
        num_operators,
        labor_rate,
    )
    # Create separate block to hold operating parameters.
    add_operating_params(
        m,
        profit,
        opt_var_oc_params,
        workers_per_discr_unit,
        yearly_cost_per_unit,
        capital_cost_per_unit,
        processing_rate,
        num_operators,
        labor_rate,
    )

    ### Costing parameters
    # Check that costing parameters are feasible.
    check_costing_params(m, discretized_capex)
    # Create a separate block to hold costing parameters.
    add_costing_params(m, discretized_capex)
    ### Objective function parameters
    # Check that objective function parameters are feasible.
    check_objective_function_params(m, obj_func)
    # Create a separate block to hold costing parameters.
    add_objective_function_params(m, obj_func)

    ### Mass balances
    # Generate mass balance parameters.
    add_mass_balance_params(m)
    # Generate mass balance variables.
    add_mass_balance_vars(m)
    # Generate mass balance constraints.
    add_mass_balance_cons(m)

    ### Logic constraints
    # Generate logic parameters.
    add_logic_params(m)
    # Generate logic variables.
    add_logic_vars(m)
    # Generate logic constraints.
    add_logic_cons(m)

    # define disassembly works set
    j_dis = 1
    m.J_dis = pyo.RangeSet(j_dis)
    m.K_dis = pyo.RangeSet(Options_in_stage[j_dis])  # options in disassembly stage
    dis_workers_range = pyo.RangeSet(0, max_dis_workers)
    jk_dis = []
    jkw_dis = []  # for declaring bin vars
    for k_dis in m.K_dis:
        jk_dis.append((j_dis, k_dis))
        for w_dis in pyo.RangeSet(0, max_dis_by_option[j_dis, k_dis]):
            jkw_dis.append((j_dis, k_dis, w_dis))

    m.disOpts = pyo.Set(within=m.J_dis * m.K_dis, initialize=jk_dis)
    m.DisOptWorkersSet = pyo.Set(
        within=m.J_dis * m.K_dis * dis_workers_range, initialize=jkw_dis
    )
    m.DisOptWorkers = pyo.Var(m.DisOptWorkersSet, domain=pyo.Binary)

    ################################################ Cash Flow Constraints ################################################
    m.OC_var_cons = pyo.ConstraintList()
    m.DisWorkerCons = pyo.ConstraintList()  # for the disassembly stage
    # calculate profit
    m.final_opt_set = pyo.Set(initialize=final_opt_list)  # list of final node list
    m.profit_opt_cons = pyo.ConstraintList()

    m.Byprods = pyo.Set(initialize=byprods)  # set of byproducts

    # opex and profit only calculated once production starts
    for t in pyo.RangeSet(prod_start, plant_end):
        m.plantYear[t].OC_var = pyo.Var(m.OptSet, domain=pyo.Reals)

        # must be enough dis units to handle incoming feed each year.
        m.DisWorkerCons.add(
            expr=m.plantYear[t].P_entering
            <= sum(
                Dis_Rate[elem]
                * sum(
                    i * m.DisOptWorkers[elem + (i,)]
                    for i in pyo.RangeSet(0, max_dis_by_option[elem])
                )
                for elem in m.disOpts
            )
        )

        for elem in m.disOpts:
            # only 1 'amount' of workers can be chosen.
            m.DisWorkerCons.add(
                expr=sum(
                    m.DisOptWorkers[elem + (i,)]
                    for i in pyo.RangeSet(0, max_dis_by_option[elem])
                )
                == 1
            )

            # if a disassembly option is not chosen, then that 'amount' must be 0
            m.DisWorkerCons.add(expr=m.binOpt[elem] >= 1 - m.DisOptWorkers[elem + (0,)])

        # calculate variable costs for disassembly stage
        j = 1
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            option = (j, k)
            m.OC_var_cons.add(
                expr=m.plantYear[t].OC_var[option]
                == sum(
                    m.DisOptWorkers[option + (i,)] * i
                    for i in pyo.RangeSet(0, max_dis_by_option[elem])
                )
                * YCU[option]
            )

        # calculate variable costs for rest of stages
        for elem in m.OptSet:
            j = 0
            k = 0
            j, k = elem
            if j > 1:  # already calculated opex for discretized options
                m.OC_var_cons.add(
                    expr=m.plantYear[t].OC_var[elem]
                    == N_OC_var[elem]["a"]
                    * sum(m.plantYear[t].F_in[elem + (c,)] for c in m.KeyComps)
                    + N_OC_var[elem]["b"] * m.binOpt[elem]
                )

        # calculate total variable cost
        m.plantYear[t].OC_var_total = pyo.Var(domain=pyo.Reals)
        m.OC_var_cons.add(
            expr=m.plantYear[t].OC_var_total
            == sum(m.plantYear[t].OC_var[elem] for elem in m.OptSet)
        )

        m.plantYear[t].ProfitOpt = pyo.Var(m.final_opt_set)

        if obj_func == "NPV":
            for opt in m.final_opt_set:
                m.profit_opt_cons.add(
                    expr=m.plantYear[t].ProfitOpt[opt]
                    == sum(
                        m.plantYear[t].F_out[opt + (c,)] * Profit[opt][c]
                        for c in m.KeyComps
                    )
                )
        elif obj_func == "COR":
            for opt in m.final_opt_set:
                m.profit_opt_cons.add(
                    expr=m.plantYear[t].ProfitOpt[opt]
                    == sum(
                        m.plantYear[t].F_out[opt + (c,)] * REE_to_REO_Conversion[c]
                        for c in m.KeyComps
                    )
                    * m.COR
                )
        else:
            sys.exit("Neither COR nor NPV specified as objective function.")

        m.plantYear[t].Profit = pyo.Var(domain=pyo.NonNegativeReals)
        m.plantYear[t].profit_con = pyo.Constraint(
            expr=m.plantYear[t].Profit
            == sum(m.plantYear[t].ProfitOpt[opt] for opt in m.final_opt_set)
        )

    ################################################ OC_fixed Constraints ################################################
    # calculate the cost of labor
    m.workers = pyo.Var(m.OptSet, domain=pyo.Reals)
    m.worker_cons = pyo.ConstraintList()
    m.COL_cons = pyo.ConstraintList()

    # calculate the number of workers for each stage
    dis_workers_range = pyo.RangeSet(0, max_dis_workers)
    j = 0
    for j in pyo.RangeSet(1, numStages):
        k = 0
        elem = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            elem = (j, k)

            if j == 1:
                m.worker_cons.add(
                    expr=m.workers[j, k]
                    == sum(
                        i * m.DisOptWorkers[j, k, i]
                        for i in pyo.RangeSet(0, max_dis_by_option[elem])
                    )
                    * num_workers[elem]
                )

            else:
                m.worker_cons.add(
                    expr=m.workers[j, k] == m.binOpt[j, k] * num_workers[elem]
                )

    m.workers_range = pyo.RangeSet(0, max_workers)
    m.bin_workers = pyo.Var(m.workers_range, domain=pyo.Binary)
    m.total_workers = pyo.Var(domain=pyo.Reals)

    m.worker_cons.add(
        expr=sum(m.workers[elem] for elem in m.OptSet)
        <= sum(i * m.bin_workers[i] for i in m.workers_range)
    )
    m.worker_cons.add(expr=sum(m.bin_workers[i] for i in m.workers_range) == 1)

    m.worker_cons.add(
        expr=m.total_workers == sum(i * m.bin_workers[i] for i in m.workers_range)
    )

    # calculate total COL
    m.COL_Total = pyo.Var(domain=pyo.NonNegativeReals)
    m.COL_Total_con = pyo.Constraint(expr=m.COL_Total == m.total_workers * labor_rate)

    ################################################ CAPEX Constraints ################################################
    t_max = maxFeedEnteringYear

    # set of non-disassembly options

    # TPC = BEC * LF
    m.BEC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.BEC_cons = pyo.ConstraintList()

    # make a var for the max flow entering each option. Used to calculate capex
    maxFlowUB = maxFeedEntering * sum(Prod_comp_mass[c] for c in Tracked_comps)
    m.BEC_max_flow = pyo.Var(m.OptSet_no_dis, bounds=(0, maxFlowUB))
    m.BEC_max_flow_cons = pyo.ConstraintList()

    # Calculate disassembly BEC (disassembly stage is always the first stage)
    k = 0
    for k in pyo.RangeSet(Options_in_stage[1]):
        m.BEC_cons.add(
            expr=m.BEC[1, k]
            == sum(
                i * m.DisOptWorkers[1, k, i]
                for i in pyo.RangeSet(0, max_dis_by_option[1, k])
            )
            * CU[1, k]
        )

    # calculate BEC for rest of stages
    for elem in m.OptSet:
        j, k = elem
        if j > 1:  # already calculated capex for disassembly stage
            # calculate max flow stream for capex
            # find the max flow (flow must be >= all flows entering it)
            for t in pyo.RangeSet(prod_start, plant_end):
                m.BEC_max_flow_cons.add(
                    expr=m.BEC_max_flow[elem]
                    >= sum(m.plantYear[t].F_in[elem + (c,)] for c in m.KeyComps)
                )

    j = 0
    k = 0
    for j in pyo.RangeSet(2, numStages):
        for k in pyo.RangeSet(Options_in_stage[j]):
            # get x and y data
            flowData = list(Discretized_CAPEX[str((j, k))]["Flowrates"].values())
            CAPEXData = list(Discretized_CAPEX[str((j, k))]["Costs"].values())

            # use m.add_component to generate all piecewise functions
            # piecewise = Piecewise(yval, xval, *kwargs)
            piecewise = pyo.Piecewise(
                m.BEC[j, k],
                m.BEC_max_flow[j, k],
                pw_pts=flowData,
                pw_constr_type="EQ",
                f_rule=CAPEXData,
                pw_repn="SOS2",
            )
            optName = (j, k)
            print("Piecewise_Node" + str(optName))
            m.add_component("Piecewise_Node" + str(optName), piecewise)

    m.TPC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.TPC_cons = pyo.ConstraintList()
    # TPC for disassembly (no installation factor for disassembly)
    j = 1  # disassembly takes place in first stage always
    k = 0
    for k in pyo.RangeSet(Options_in_stage[j]):
        m.TPC_cons.add(expr=m.TPC[j, k] == m.BEC[j, k])

    # multiply all TPCs by lang factor
    j = 0
    for j in pyo.RangeSet(2, numStages):
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            m.TPC_cons.add(expr=m.TPC[j, k] == m.BEC[j, k] * LF)

    # calculate total TPC
    m.Total_TPC = pyo.Var(domain=pyo.NonNegativeReals)
    m.Total_TPC_con = pyo.Constraint(
        expr=m.Total_TPC == sum(m.TPC[elem] for elem in m.OptSet)
    )

    # calculate TOC
    m.TOC = pyo.Var(domain=pyo.NonNegativeReals)
    m.TOC_con = pyo.Constraint(expr=m.TOC == m.Total_TPC * TOC_factor)

    # calculate node TOCs
    m.node_TOC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.node_TOC_cons = pyo.ConstraintList()

    j = 0
    for j in pyo.RangeSet(1, numStages):
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            elem = (j, k)
            m.node_TOC_cons.add(
                expr=m.node_TOC[elem]
                == m.TPC[elem] + 0.027 * m.TPC[elem] + 0.15 * m.TPC[elem]
            )

    #### Construct cash flows
    m.CF = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.CF_cons = pyo.ConstraintList()
    ## CAPEX cons
    # yearly TOC expenditure
    m.TOC_exp = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.TOC_exp_cons = pyo.ConstraintList()

    ## OPEX cons
    # net earnings
    # gross earnings
    m.GE = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.GE_cons = pyo.ConstraintList()
    # revenue
    m.Rev = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.Rev_cons = pyo.ConstraintList()
    # fixed operation costs
    m.OC_fixed = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.OC_fixed_cons = pyo.ConstraintList()
    # variable operating costs
    m.OC_var = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    # plant overhead costs
    m.OH = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.OH_cons = pyo.ConstraintList()

    t = 0
    for t in plant_life_range:
        # capital expenditure only in first three years
        if t <= (plant_start + 2):
            m.TOC_exp_cons.add(
                expr=m.TOC_exp[t]
                == ((1 + i_CAP_esc) ** (t - plant_start))
                * f_exp[t - plant_start]
                * m.TOC
            )
        else:
            # no capital escalation after third year
            m.TOC_exp_cons.add(expr=m.TOC_exp[t] == 0)

        # opex and revenue only once production begins (starting second year)
        if prod_start <= t:
            # fixed operating costs
            m.OC_fixed_cons.add(
                expr=m.OC_fixed[t]
                == m.COL_Total
                + 0.02 * m.Total_TPC
                + 0.1 * m.COL_Total
                + 0.01 * m.plantYear[t].Profit
                + 0.2 * m.COL_Total
                + 0.25 * m.COL_Total
                + 0.01 * m.Total_TPC
            )

            # variable costs
            m.OC_var_cons.add(expr=m.OC_var[t] == m.plantYear[t].OC_var_total)

            ## revenue
            # consider byproduct valorization if specified by user
            if consider_byprod_val == True:
                m.Rev_cons.add(
                    expr=m.Rev[t]
                    == m.plantYear[t].Profit + m.plantYear[t].Byprod_Profit
                )
            else:  # otherwise don't
                m.Rev_cons.add(expr=m.Rev[t] == m.plantYear[t].Profit)

        else:
            m.OC_fixed_cons.add(expr=m.OC_fixed[t] == 0)
            m.OC_var_cons.add(expr=m.OC_var[t] == 0)
            m.Rev_cons.add(expr=m.Rev[t] == 0)

        # calculate plant overhead
        m.OH_cons.add(expr=m.OH[t] == 0.2 * (m.OC_fixed[t] + m.OC_var[t]))

        # calculate the gross earnings
        m.GE_cons.add(
            expr=m.GE[t]
            == (m.Rev[t] - m.OC_fixed[t] - m.OC_var[t] - m.OH[t])
            * (1 + i_OC_esc) ** (t - 2024)
        )

        # calculate the cash flow
        m.CF_cons.add(expr=m.CF[t] == m.GE[t] - m.TOC_exp[t])

    ### testing
    m.bin_test_cons = pyo.ConstraintList()

    if obj_func == "NPV":

        def obj_rule(m):
            return sum(
                m.CF[t] / (1 + ATWACC) ** (t - plant_start) for t in plant_life_range
            )

        m.obj = pyo.Objective(rule=obj_rule, sense=pyo.maximize)

    elif obj_func == "COR":
        m.NPV = pyo.Var(domain=pyo.Reals)
        m.NPV_con1 = pyo.Constraint(
            expr=m.NPV
            == sum(
                m.CF[t] / (1 + ATWACC) ** (t - plant_start) for t in plant_life_range
            )
        )
        m.NPV_con2 = pyo.Constraint(expr=m.NPV == 0)
        m.obj = pyo.Objective(expr=m.COR, sense=pyo.minimize)

    return m
