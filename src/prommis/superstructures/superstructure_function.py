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
        plant_lifetime: (int) the total lifetime of the plant, including plant construction. Must be at least three years.
    """

    ## Check that plant lifetime is at least three years
    if plant_lifetime < 3:
        raise ValueError("Plant lifetime must be a minimum of three years.")


def add_plant_lifetime_params_block(m, plant_start, plant_lifetime):
    """
    This function builds the rest of the plant lifetime parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model
        plant_start: (int) the year that plant construction begins.
        plant_lifetime: (int) the total lifetime of the plant, including plant construction. Must be at least three years.
    """

    # Create a block to hold plant lifetime parameters
    m.plant_lifetime_params = pyo.Block()

    m.plant_lifetime_params.plant_start = pyo.Param(
        initialize=plant_start, doc="The year that plant construction begins."
    )
    m.plant_lifetime_params.plant_lifetime = pyo.Param(
        initialize=plant_lifetime,
        doc="The total lifetime of the plant, including plant construction. Must be at least three years.",
    )

    ## Calculate other necessary params from user input
    # first year is construction
    prod_start = plant_start + 1
    m.plant_lifetime_params.prod_start = pyo.Param(
        initialize=prod_start, doc="The first year of plant production."
    )
    # final year plant is in production
    plant_end = plant_start + plant_lifetime - 1
    m.plant_lifetime_params.plant_end = pyo.Param(
        initialize=plant_end, doc="The final year of plant production."
    )
    # lifetime of the plant
    m.plant_lifetime_params.plant_life_range = pyo.RangeSet(
        plant_start, plant_end, doc="Lifetime of the plant."
    )
    # operational lifetime of the plant
    m.plant_lifetime_params.operational_range = pyo.RangeSet(
        prod_start, plant_end, doc="Operational lifetime of the plant."
    )

    return m


###################################################################################################
### Feed Parameters
def check_feed_params(
    m, available_feed, collection_rate, tracked_comps, prod_comp_mass
):
    """
    The function checks that the feed parameter inputs are feasible.

    Args:
        m: pyomo model
        available_feed: (dict) total feedstock available for recycling each year
        collection_rate: (float) how much available feed is processed by the plant each year
        tracked_comps: (list) list of tracked components
        prod_comp_mass: (dict) mass of tracked components per EOL product
    """

    ## Check that available feed is provided for each year of plant operation.
    # Extract the relevant sets
    feed_years = set(available_feed.keys())
    operational_years = set(m.plant_lifetime_params.operational_range.data())
    lifetime_years = set(m.plant_lifetime_params.plant_life_range.data())

    if feed_years != operational_years:
        raise ValueError(
            "Years of available_feed do not match the plant's operational period. "
            f"Expected years: {sorted(operational_years)}, "
            f"but got: {feed_years}"
        )

    ## Check that none of the available feeds passed are negative
    if any(v < 0 for v in available_feed.values()):
        negative_years = [year for year, value in available_feed.items() if value < 0]
        raise ValueError(
            f"available_feed contains negative values for years: {negative_years}. "
            "Feedstock availability cannot be negative."
        )

    ## Check that available feed is not all zero
    if not any(value != 0 for value in available_feed.values()):
        raise ValueError(
            "All values in available_feed are zero. At least one year must have non-zero feedstock available."
        )

    ## Check that collection rate is positive value
    if collection_rate <= 0:
        raise ValueError("Collection rate must be a positive value.")

    ## Check that at least one tracked component is specified
    if not tracked_comps:
        raise ValueError(
            "tracked_comps list is empty. At least one component must be tracked."
        )

    ## Check an amount contained within an EOL product is specified for each tracked component
    prod_comp_mass_keys = set(prod_comp_mass.keys())
    tracked_comps_set = set(tracked_comps)
    if prod_comp_mass_keys != tracked_comps_set:
        raise ValueError(
            f"prod_comp_mass keys don't match up with the set of tracked components."
        )

    ## Check that amounts contained within EOL product for each tracked component is non-negative
    if any(v < 0 for v in prod_comp_mass.values()):
        negative_prods = [prod for prod, amount in prod_comp_mass.items() if amount < 0]
        raise ValueError(
            f"prod_comp_mass contains negative values for the tracked components: {negative_prods}. "
            "Amounts cannot be negative"
        )

    ## Raise warning if amounts contained within EOL product for a tracked component is zero
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
        m: pyomo model
        available_feed: (dict) total feedstock available for recycling each year
        collection_rate: (float) collection rate for how much available feed is processed by the plant each year
        tracked_comps: (list) list of tracked components
        prod_comp_mass: (dict) mass of tracked components per EOL product
    """

    # Create a block to hold plant lifetime parameters
    m.feed_params = pyo.Block()

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

    ## Calculate other necessary params from user input
    # calculate feed entering parameter based on yearly available feedstock and collection rate
    feed_entering = copy.deepcopy(available_feed)
    for key in feed_entering:
        feed_entering[key] = available_feed[key] * collection_rate
    m.feed_params.feed_entering = pyo.Param(
        m.plant_lifetime_params.operational_range,
        initialize=feed_entering,
        doc="The amount of feed entering the plant each year.",
    )

    # max feed entering plant over production period
    max_feed_entering = max(feed_entering.values())
    m.feed_params.max_feed_entering = pyo.Param(
        initialize=max_feed_entering,
        doc="The max yearly feed that enters the plant over the production period.",
    )

    # year in which max feed enters plant
    max_feed_entering_year = max(feed_entering, key=feed_entering.get)
    m.feed_params.max_feed_entering_year = pyo.Param(
        initialize=max_feed_entering_year,
        doc="The year that the max feed enters the plant.",
    )

    return m


###################################################################################################
### Superstructure Formulation Parameters
def check_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
):
    """
    The function checks that the superstructure formulation parameter inputs are feasible.

    Args:
        m: pyomo model
        num_stages: (int) number of total stages
        options_in_stage: (dict) number of options in each stage
        options_outlets: (dict) set of options k' in stage j+1 connected to option k in stage j
        discrete_opts: (list) list of options that utilize discrete units
        option_eff: (dict) tracked component retention efficiency for each option
    """

    ## Check that there is at least 2 stage
    if num_stages < 2:
        raise ValueError("There must be at least 2 processing stages.")

    ## The stages must be numbered starting at 1 and count up.
    # Extract the relevant sets
    num_stages_set = set(pyo.RangeSet(num_stages).data())
    options_in_stage_keys_set = set(options_in_stage.keys())

    if num_stages_set != options_in_stage_keys_set:
        raise ValueError(
            "Stages must start at 1 and count up. Each stage must contain at least 1 option. options_in_stage does not follow this convention. "
            f"Expected keys: {num_stages_set}, "
            f"but got: {options_in_stage_keys_set}"
        )

    ## Check that connections between options in superstructure are feasible
    # Check that each option in stage j is connected to an option in stage j+1
    missing_values = [key for key, value in option_outlets.items() if value is None]
    if missing_values:
        raise ValueError(
            f"Options {missing_values} are missing connections in the next stage, as defined by option_outlets."
        )

    # Check that each option in stage j is connected to an option in the preceding stage, stage j-1
    disconnected_options = []

    # Iterate over stages starting from 2 (since stage 1 has no predecessors)
    for current_stage in range(2, num_stages + 1):
        previous_stage = current_stage - 1

        # Get number of options in current and previous stages
        num_current_options = options_in_stage[current_stage]
        num_previous_options = options_in_stage[previous_stage]

        # Check each option in the current stage
        for current_option in range(1, num_current_options + 1):
            connected = False

            # Check all options in the previous stage
            for previous_option in range(1, num_previous_options + 1):
                # Get outlets for (previous_stage, previous_option)
                outlets = option_outlets.get((previous_stage, previous_option), [])
                if current_option in outlets:
                    connected = True
                    break  # No need to check further

            if not connected:
                disconnected_options.append((current_stage, current_option))

    # Raise error if disconnected options exist
    if disconnected_options:
        error_msg = "The following options are not connected from the previous stage:\n"
        error_msg += "\n".join(
            f"  - Stage {stage}, Option {option}"
            for stage, option in disconnected_options
        )
        raise ValueError(error_msg)
    
    ## Check that all discrete options listed are feasible
    all_opts_set = {
        (j, k)
        for j in range(1, num_stages + 1)
        for k in range(1, options_in_stage[j] + 1)
    }
    discr_opts_set = set(discrete_opts)

    if not (discr_opts_set <= all_opts_set):
        raise ValueError("Discrete options listed are not feasible.")

    ## Check that an option efficiency is defined for each option and is nonnegative
    # Extract the relevant sets
    tracked_comps = set(m.feed_params.tracked_comps.data())

    # Lists to track issues
    missing_efficiencies = []
    negative_effs = []

    for j in range(1, num_stages + 1):
        for k in range(1, options_in_stage[j] + 1):

            # Check for missing efficiencies
            option_eff_tracked_comp_keys = set(option_eff[(j, k)].keys())
            if option_eff_tracked_comp_keys != tracked_comps:
                missing = tracked_comps - option_eff_tracked_comp_keys
                missing_efficiencies.append((j, k, missing))

            # Check for negative efficiencies
            negative_comps = [c for c, eff in option_eff[j, k].items() if eff < 0]
            if negative_comps:
                negative_effs.append((j, k, negative_comps))

    # Raise error for missing efficiencies
    if missing_efficiencies:
        msg = "Efficiencies not specified for all tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing components: {missing}"
            for j, k, missing in missing_efficiencies
        )
        raise ValueError(msg)

    # Raise error for negative efficiencies
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
        m: pyomo model
        num_stages: (int) number of total stages
        options_in_stage: (dict) number of options in each stage
        option_outlets: (dict) set of options k' in stage j+1 connected to option k in stage j
        discrete_opts: (list) list of options that utilize discrete units
        option_eff: (dict) tracked component retention efficiency for each option
    """

    m.supe_form_params = pyo.Block()

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

    # Define a parameter for the max number of options in any of the stages
    max_options = max(options_in_stage.values())
    m.supe_form_params.max_options = pyo.Param(
        initialize=max_options, doc="The max number of options in any of the stages."
    )
    m.supe_form_params.max_options_set = pyo.RangeSet(
        1, max_options, doc="Set containing max number of options in any of the stages."
    )

    ## Build a set of all the options in the superstructure
    # Create a set to hold all options in the superstructure. Indexed by stages 'j', and option 'k' in stage 'j'
    m.supe_form_params.all_opts_set = pyo.Set(
        initialize=(
            (j, k)
            for j in m.supe_form_params.stages_set
            for k in range(1, options_in_stage[j] + 1)
        ),
        doc="Set containing all options in the superstructure."
    )

    ## Build a set containing all discrete options
    # Define necessary set
    discrete_opts_set = set(discrete_opts)
    m.supe_form_params.discrete_opts_set = pyo.Set(
        initialize=(
            (opt)
            for opt in discrete_opts_set
        ),
        doc="Set containing all the options which utilize discrete units."
    )

    ## Build a set containing all continuous options in the superstructure
    # Define necessary sets
    all_opts_set = set(m.supe_form_params.all_opts_set.data())
    continuous_opts_set = all_opts_set - discrete_opts_set

    m.supe_form_params.continuous_opts_set = pyo.Set(
        initialize=(
            (opt)
            for opt in continuous_opts_set
        ),
        doc="Set containing all the continuous options."
    )


    m.supe_form_params.option_outlets = pyo.Param(
        m.supe_form_params.all_opts_set,
        initialize=option_outlets,
        doc="Defines the set of options k' in stage j+1 connected to option k in stage j.",
    )

    ## Define a parameter to hold the tracked component efficiencies for each option in the superstructure
    # Efficiency defined for each tracked component 'c' for each option 'k' in each stage 'j'
    def option_eff_initialize(m, j, k, c):
        return option_eff[(j, k)][c]

    m.supe_form_params.option_eff = pyo.Param(
        m.supe_form_params.all_opts_set,
        m.feed_params.tracked_comps,
        initialize=option_eff_initialize,
    )

    # Define a set containing all of the options in the final stage
    final_opts_list = [
        (num_stages, k) for k in range(1, options_in_stage[num_stages] + 1)
    ]
    m.supe_form_params.final_opts_set = pyo.Set(
        initialize=final_opts_list,
        doc="Set containing all of the options in the final stage.",
    )

    return m


###################################################################################################
### Operating Parameters
def check_operating_params(
    m,
    profit,
    opt_var_oc_params,
    workers_per_discr_unit,
    yearly_cost_per_unit,
    cost_per_unit,
    processing_rate,
    num_workers,
    labor_rate,
):
    """
    This function checks that all the operating parameters are feasible.

    Args:
        m: pyomo model.
        profit: (dict) profit per unit of product in terms of tracked components.
        opt_var_oc_params: (dict) holds the variable operating cost param for options that are continuous. Variable operating costs assumed to be proportional to the feed entering the option.
        workers_per_discr_unit: (dict) number of workers needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) yearly operating costs per unit for options which utilize discrete units.
        cost_per_unit: (dict) cost per unit for options which utilize discrete units.
        processing_rate: (dict) disassembly rate per unit for options that utilize discrete units. In terms of units of incoming feed processed per year per unit.
        num_workers: (dict) number of workers needed for each option.
        labor_rate: (float) yearly wage per worker.
    """

    ## Check that profit per product is defined for all options in the final stage
    profit_opt_keys = set(profit.keys())
    final_opts_set = set(m.supe_form_params.final_opts_set)

    if profit_opt_keys != final_opts_set:
        raise ValueError(
            "Must include profit per unit of product for all options in the final stage."
        )

    ## Check that profit per product is in terms of the tracked components for all options in the final stage and are all nonnegative
    # create relevant sets
    tracked_comps = set(m.feed_params.tracked_comps.data())
    # list to track issues
    missing_profit_comps = []
    negative_profit_comps = []

    for opt in m.supe_form_params.final_opts_set:

        # Check for missing profits
        profit_tracked_comps = set(profit[opt].keys())
        if profit_tracked_comps != tracked_comps:
            missing = tracked_comps - profit_tracked_comps
            missing_profit_comps.append(opt + (missing,))

        # Check for negative profits
        negative_comps = [c for c, profit in profit[opt].items() if profit < 0]
        if negative_comps:
            negative_profit_comps.append(opt + (negative_comps,))

    if missing_profit_comps:
        msg = "Profits not specified for all tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing components: {missing}"
            for j, k, missing in missing_profit_comps
        )
        raise ValueError(msg)

    if negative_profit_comps:
        msg = "Profits some tracked components are listed as negative in the following options\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) negative components: {negative}"
            for j, k, negative in negative_profit_comps
        )
        raise ValueError(msg)

    ## Check that variable operating cost params are defined for all continuous options
    opt_var_oc_params_keys = set(opt_var_oc_params.keys())
    discr_opts_set = set(m.supe_form_params.discrete_opts_set.data())
    continuous_opts_set = set(m.supe_form_params.continuous_opts_set.data())

    if opt_var_oc_params_keys != continuous_opts_set:
        raise ValueError(
            "Variable operating cost params not defined for all continuous options."
        )

    ## Check that both necessary variable operating cost parameters ('a' and 'b') are defined for all continuous options
    var_oc_params = set(["a", "b"])
    # list to track issues
    missing_var_oc_params = []

    for opt in continuous_opts_set:

        # check for missing variable operating cost parameters
        params = set(opt_var_oc_params[opt].keys())
        if params != var_oc_params:
            missing = var_oc_params - params
            missing_var_oc_params.append(opt + (missing,))

    if missing_var_oc_params:
        msg = "not all variable operating cost parameters defined in the following options:\n"
        msg += "\n".join(
            f"  Option (stage={j}, option={k}) missing parameters: {missing}"
            for j, k, missing in missing_var_oc_params
        )
        raise ValueError(msg)

    ## Check that workers per discrete unit is defined for all options that utilize discrete units
    workers_per_discr_unit_keys = set(workers_per_discr_unit.keys())

    if workers_per_discr_unit_keys != discr_opts_set:
        raise ValueError("workers_per_discr_unit not defined for all discrete options.")

    ## Check that workers per discrete unit are all defined to be non-negative
    negative_workers_per_discr_units = [
        opt for opt, w in workers_per_discr_unit.items() if w < 0
    ]
    if negative_workers_per_discr_units:
        raise ValueError("Workers per discrete unit must all be non-negative.")

    ## Check that yearly cost per unit is defined for all discrete options
    yearly_cost_per_unit_keys = set(yearly_cost_per_unit.keys())
    if yearly_cost_per_unit_keys != discr_opts_set:
        raise ValueError(
            "yearly_cost_per_unit must be defined for all discrete options."
        )

    ## Check that the yearly cost per unit values are all non-negative
    negative_yearly_cost_per_unit = [
        opt for opt, cost in yearly_cost_per_unit.items() if cost < 0
    ]
    if negative_yearly_cost_per_unit:
        raise ValueError("yearly_cost_per_unit values must all be non-negative.")

    ## Check that cost per unit is defined for all discrete options
    cost_per_unit_keys = set(cost_per_unit.keys())
    if cost_per_unit_keys != discr_opts_set:
        raise ValueError("cost_per_unit must be defined for all discrete options.")

    ## Check that the cost per unit values are all non-negative
    negative_cost_per_unit = [opt for opt, cost in cost_per_unit.items() if cost < 0]
    if negative_cost_per_unit:
        raise ValueError("cost_per_unit values must all be non-negative.")

    ## Check that prosessing rate is defined for all discrete options
    processing_rate_keys = set(processing_rate.keys())
    if processing_rate_keys != discr_opts_set:
        raise ValueError("processing_rate must be defined for all discrete options.")

    ## Check that all processing rates are positive
    nonpositive_processing_rate = [
        opt for opt, rate in processing_rate.items() if rate <= 0
    ]
    if nonpositive_processing_rate:
        raise ValueError("processing rates must all be positive.")

    ## Check that num_workers is defined for all continuous options
    num_workers_keys = set(num_workers.keys())
    if num_workers_keys != continuous_opts_set:
        raise ValueError(
            "num_workers must be defined for all continuous options, and must not be defined for discrete options."
        )

    ## Check that num_workers values are all non-negative
    negative_num_workers = [opt for opt, workers in num_workers.items() if workers < 0]
    if negative_num_workers:
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
    cost_per_unit,
    processing_rate,
    num_workers,
    labor_rate,
):
    """
    This function builds the rest of the operating parameters from the ones provided by the user, and adds them all to a
    block.

    Args:
        m: pyomo model
        profit: (dict) profit per unit of product in terms of tracked components
        opt_var_oc_params: (dict) holds the variable operating cost params for options that don't utilize discrete units.
        workers_per_discr_unit: (dict) number of workers needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) yearly operating costs per unit for options which utilize discrete units.
        cost_per_unit: (dict) cost per unit for options which utilize discrete units.
        processing_rate: (dict) disassembly rate per unit for options that utilize discrete units. In terms of units of incoming feed processed per year per unit.
        num_workers: (dict) number of workers needed by option for options that don't utilize discrete units.
        labor_rate: (float) yearly wage per worker
    """

    m.operating_params = pyo.Block()

    ## Define a parameter to hold the profit for all options in the final stage in terms of the tracked components
    def profit_initialization(m, j, k, c):
        return profit[(j, k)][c]
    
    m.operating_params.profit = pyo.Param(
        m.supe_form_params.final_opts_set,
        m.feed_params.tracked_comps,
        initialize=profit_initialization,
        doc="Holds the profit for all options in the final stage in terms of the tracked components."
    )
    
    # Define set containing the variable operating cost parameters
    m.operating_params.var_oc_params_set = pyo.Set(initialize=["a", "b"], doc="Set containing the necessary parameters for calculating the variable operating costs for all continuous options.")

    ## Define a parameter to hold the variable operating cost params for continuous options
    def opt_var_oc_params_initialization(m, j, k, var_oc_param):
        return opt_var_oc_params[(j, k)][var_oc_param]
    
    m.operating_params.opt_var_oc_params = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        m.operating_params.var_oc_params_set,
        initialize=opt_var_oc_params_initialization,
        doc="Holds all the variable operating costs parameter values for all continuous options."
    )
    m.operating_params.opt_var_oc_params.display()

    ## Define a parameter to hold the number of workers needed per discrete unit for options that utilize discrete units
    # m.operating_params.

    # # calculate max disassembly units possible for each option
    # max_dis_by_option = copy.deepcopy(Dis_Rate)
    # for key in max_dis_by_option.keys():
    #     max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
    # # max_dis_workers = max(max_dis_by_option.values()) + 10
    # max_dis_workers = max(max_dis_by_option.values())
    # # calculate max possible workers for process
    # max_workers = max_dis_workers + numStages * math.ceil(max(num_workers.values()))



def build_model(
    ###################################################################################################
    ### Plant Lifetime Parameters
    # first year of plant lifetime
    plant_start: int,
    # lifetime of plant
    plant_lifetime: int,
    ###################################################################################################
    ###################################################################################################
    ### Feed parameters
    # Total feedstock available for recycling each year
    # of the form -> year: amount, ...
    available_feed: dict,
    # collection rate for how much of the available feed is processed by the plant each year
    collection_rate: float,
    # list of tracked components (this is assumed to be the elemental rare earths and other elemental contaminants.
    # Refer to examples)
    tracked_comps: list,
    # mass of tracked component per EOL Product (kg component / EOL product)
    # of the form -> component: amount, ...
    prod_comp_mass: dict,
    ###################################################################################################
    ###################################################################################################
    ### Superstructure formulation parameters
    # number of total stages
    num_stages: int,
    # number of options in each stage. Of the form -> stage number: number of options
    options_in_stage: dict,
    # set of options k' in stage j+1 connected to option k in stage j
    # of the form -> (j, k): [1, 2, 3...], ...
    option_outlets: dict,
    # dictionary of tracked component retention efficiency for each option
    # of the form -> (j, k): {"comp1": eff1, "comp2": eff2, ...}, ...
    option_eff: dict,
    ###################################################################################################
    ###################################################################################################
    ### Operating Parameters
    # profit per kg of product in terms of tracked components
    # of the form -> (j, k): {"comp1": price1, "comp2": price2, ...}, ...
    profit: dict,
    # conversion of kg REE/Fe to kg REO/Fe2O3
    # of the form -> "REE": conv_factor1, ...
    # see examples for more detail
    ree_to_reo_conversion: dict,
    # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it.
    # OPEX = a*F_in + b*y
    # of the form -> (j, k): {"a": val1, "b": val2}, ...
    opt_operating_cost_params: dict,
    # list of the options that utilize discrete units (typically disassembly options)
    # of the form -> [(1, 1), (1, 2), ...]
    dis_opts: list,
    # number of workers needed per discrete unit for options that utilize discrete units
    workers_per_discrete_unit: dict,
    # number of workers needed by option for options that do not utilize discrete units
    # of the form -> (j, k): val1, ...
    num_workers: dict,
    labor_rate: float,  # yearly wage per worker
    # yearly operating costs per unit ($/unit*yr) for options which utilize discrete units
    # of the form -> (j,k): cost1, ...
    yearly_cost_per_unit: dict,
    # installation cost per unit for options which utilize discrete units
    # of the form -> (j,k): cost1, ...
    cost_per_unit: dict,
    # disassembly rate for each disassembly option (in terms of EOL products disassembled per year per unit)
    # of the form -> (j,k): rate1, ...
    Dis_Rate: dict,
    ###################################################################################################
    ###################################################################################################
    ### Costing Parameters
    # Define Python Dictionary with discretized cost by flows for each option.
    # see example for form
    Discretized_CAPEX: dict,
    ###################################################################################################
    ###################################################################################################
    ### Choice of objective function. Options are 'NPV' or 'COR'
    obj_func: str,
    ###################################################################################################
    ###################################################################################################
    ### Consideration of environmental impacts parameters
    # boolean to decide whether or not to consider environmental impacts
    consider_environ_impacts: bool,
    # environmental impacts matrix (kg CO2e per incoming flowrate)
    # of the form -> (j, k): factor1, ...
    environ_impacts: dict,
    epsilon: float,  # epsilon factor for generating Pareto front
    ###################################################################################################
    ###################################################################################################
    ### Byproduct valorization
    # boolean to decide whether or not to consider the valorization of byproducts
    consider_byprod_val: bool,
    # list of byproducts.
    # of the form -> [byprod1, byprod2, ...]
    byprods: list,
    # dictionary of values for each byproduct ($/kg). Negative value indicates it cost money to dispose of the byproduct
    # of the form -> byprod1: value1, ...
    byprod_vals: dict,
    # dictionary keeping track of which tracked component produces which byproduct
    # of the form -> byprod1: tracked_comp1, ...
    tracked_comp_for_byprod: dict,
    # dictionary tracking which options produce a given byproduct
    # of the form -> byprod1: [(j1, k1), (j2, k2), ...], ...
    byprod_options: dict,
    # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
    # of the form -> byprod1: {(j1, k1): eff1, ...}, ...
    byprod_options_eff: dict,
    # Conversion factors of tracked component to byproduct (kg byproduct / kg iron)
    # of the form -> byprod1: conv_factor1, ...
    TC_to_byproduct: dict,
):
    ### Plant lifetime parameters
    # Check that plant lifetime parameters are feasible
    check_plant_lifetime_params(plant_lifetime)
    # Create separate block to hold plant lifetime parameters
    m = add_plant_lifetime_params_block(m, plant_start, plant_lifetime)

    ### Feed parameters
    # Check that feed parameters are feasible
    check_feed_parameters(
        m, available_feed, collection_rate, tracked_comps, prod_comp_mass
    )
    # Create separate block to hold feed parameters
    m = add_feed_params_block(
        m, available_feed, collection_rate, tracked_comps, prod_comp_mass
    )

    ### Superstructure formulation parameters
    #

    ###################################################################################################
    ###################################################################################################
    ### Calculate other superstructure formulation parameters
    maxOptions = max(Options_in_stage.values())  # max options in any of the stages
    # make a list of the options in the final stage
    final_opt_list = [(numStages, j) for j in pyo.RangeSet(Options_in_stage[numStages])]
    ###################################################################################################
    ###################################################################################################
    ### Calculate other operating parameters
    # calculate max disassembly units possible for each option
    max_dis_by_option = copy.deepcopy(Dis_Rate)
    for key in max_dis_by_option.keys():
        max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
    # max_dis_workers = max(max_dis_by_option.values()) + 10
    max_dis_workers = max(max_dis_by_option.values())
    # calculate max possible workers for process
    max_workers = max_dis_workers + numStages * math.ceil(max(num_workers.values()))
    ###################################################################################################

    m = pyo.ConcreteModel()

    if obj_func == "COR":
        m.COR = pyo.Var(domain=pyo.NonNegativeReals)

    def yearly_mass_flow_block_rule(b, t):

        # These sets will be used to define the flow vars
        b.J = pyo.RangeSet(numStages)  # number of stages
        b.K = pyo.RangeSet(maxOptions)  # max options in a stage
        b.KeyComps = pyo.Set(initialize=Tracked_comps)  # key components
        jk = []  # for declaring bin vars
        jkc = []  # for declaring flow vars
        for j in b.J:
            for k in pyo.RangeSet(Options_in_stage[j]):
                jk.append((j, k))
                for c in b.KeyComps:
                    jkc.append((j, k, c))

        # mass balance constraints and flows only present after construction is complete (1 year)
        # amount of EOL products chosen to be recycled
        b.P_entering = pyo.Var(domain=pyo.NonNegativeReals)

        # F is stream entering each stage (except last stage). indexed by j, c, t.
        # Thus,
        b.F_stages = pyo.RangeSet(1, numStages - 1)

        b.F = pyo.Var(b.F_stages * b.KeyComps, domain=pyo.NonNegativeReals)

        # F^in is stream entering each option
        # F^out is stream leaving each option
        b.FlowSet = pyo.Set(within=b.J * b.K * b.KeyComps, initialize=jkc)
        b.F_in = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)
        b.F_out = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)

        # set of options
        b.OptSet = pyo.Set(within=b.J * b.K, initialize=jk)

        ### Mass Balances
        b.init_flow_cons = pyo.ConstraintList()  # eqn. 2
        b.inlet_flow_cons = pyo.ConstraintList()  # eqn. 1
        b.intermediate_flow_cons = pyo.ConstraintList()  # eqn. 3
        b.outlet_flow_cons = pyo.ConstraintList()  # eqn. 4

        b.inlet_flow_cons.add(expr=m.plantYear[t].P_entering == Feed_entering[t])

        for j in b.J:
            num_options = pyo.RangeSet(Options_in_stage[j])

            if j == 1:
                # eqn. 2
                for c in b.KeyComps:
                    b.init_flow_cons.add(
                        expr=b.P_entering * Prod_comp_mass[c]
                        == sum(b.F_in[j, k, c] for k in num_options)
                    )

            else:
                # eqn. 1
                for c in b.KeyComps:
                    b.inlet_flow_cons.add(
                        expr=b.F[j - 1, c] == sum(b.F_in[j, k, c] for k in num_options)
                    )

            if j != numStages:
                # eqn. 4
                for c in b.KeyComps:
                    b.outlet_flow_cons.add(
                        expr=sum(b.F_out[j, k, c] for k in num_options) == b.F[j, c]
                    )

        # eqn. 3
        for j in b.J:
            num_options = pyo.RangeSet(Options_in_stage[j])
            for k in num_options:
                for c in b.KeyComps:
                    a = Option_Eff[(j, k)][c]
                    b.intermediate_flow_cons.add(
                        expr=b.F_in[j, k, c] * a == b.F_out[j, k, c]
                    )

        ### Calculation of yearly byproducts produced
        # only if user specifies this
        if consider_byprod_val == True:
            b.Byprods = pyo.Set(initialize=byprods)  # byproducts
            b.total_yearly_byprod = pyo.Var(
                b.Byprods
            )  # yearly amounts of a byproduct produced
            b.yearly_byprod_cons = pyo.ConstraintList()
            for byprod in b.Byprods:
                byprod_val = byprod_vals[byprod]
                c = tracked_comp_for_byprod[byprod]

                b.yearly_byprod_cons.add(
                    expr=b.total_yearly_byprod[byprod]
                    == sum(
                        b.F_in[byprod_option + (c,)]
                        * byprod_options_eff[byprod][byprod_option]
                        * TC_to_byproduct[byprod]
                        for byprod_option in byprod_options[byprod]
                    )
                )

        ### consideration of yearly environmental impacts
        # only if user species this
        if consider_environ_impacts == True:
            b.total_yearly_GWP = pyo.Var(domain=pyo.NonNegativeReals)
            b.yearly_GWP = pyo.Var(b.OptSet, domain=pyo.NonNegativeReals)
            b.GWP_cons = pyo.ConstraintList()

            # calculate GWP associated with each technology option
            for j in b.J:
                num_options = pyo.RangeSet(Options_in_stage[j])
                for k in num_options:
                    b.GWP_cons.add(
                        expr=b.yearly_GWP[j, k]
                        == sum(b.F_in[j, k, c] for c in b.KeyComps)
                        * environ_impacts[j, k]
                    )

            # calculate the total GWP for the year
            b.GWP_cons.add(
                expr=b.total_yearly_GWP == sum(b.yearly_GWP[elem] for elem in b.OptSet)
            )

    # build mass balance constraints for each year the plant is in operation.
    m.plant_year = pyo.RangeSet(prod_start, plant_end)
    m.plantYear = pyo.Block(m.plant_year, rule=yearly_mass_flow_block_rule)

    # These sets will be used to define the binary and flow vars
    m.J = pyo.RangeSet(numStages)  # number of stages
    m.K = pyo.RangeSet(maxOptions)  # max options in a stage
    m.KeyComps = pyo.Set(initialize=Tracked_comps)  # key components
    jk_dis = []  # set of options
    jk_no_dis = []  # set of options (excluding disassembly)
    for j in m.J:
        for k in pyo.RangeSet(Options_in_stage[j]):
            jk_dis.append((j, k))
            if j != 1:
                jk_no_dis.append((j, k))
    m.OptSet = pyo.Set(within=m.J * m.K, initialize=jk_dis)  # set of options
    m.OptSet_no_dis = pyo.Set(within=m.J * m.K, initialize=jk_no_dis)  # set of options

    ### declare binary variables
    m.binOpt = pyo.Var(m.OptSet, domain=pyo.Binary)

    ### Logical Constraints
    # eqn. 5
    m.stage_bin_cons = pyo.ConstraintList()
    for j in m.J:
        num_options = pyo.RangeSet(Options_in_stage[j])
        m.stage_bin_cons.add(expr=sum(m.binOpt[j, k] for k in num_options) == 1)

    # eqn. 6
    m.connect_bin_cons = pyo.ConstraintList()
    j = 0
    for j in pyo.RangeSet(1, numStages - 1):
        num_options = pyo.RangeSet(Options_in_stage[j])

        k = 0
        for k in num_options:
            opt_connects = Option_outlets[(j, k)]
            m.connect_bin_cons.add(
                expr=1
                - m.binOpt[j, k]
                + sum(m.binOpt[j + 1, kp] for kp in opt_connects)
                >= 1
            )

    ## big-M constraints
    m.big_M_cons = pyo.ConstraintList()
    M = {"Nd": 0, "Dy": 0, "Fe": 0}
    val = 0
    for c in m.KeyComps:
        val = math.ceil(maxFeedEntering * Prod_comp_mass[c])
        M[c] = val

    for t in pyo.RangeSet(prod_start, plant_end):
        for j in m.J:
            num_options = pyo.RangeSet(Options_in_stage[j])

            for k in num_options:
                for c in m.KeyComps:
                    m.big_M_cons.add(
                        expr=m.plantYear[t].F_in[j, k, c] <= m.binOpt[j, k] * M[c]
                    )  # eqn. 7
                    m.big_M_cons.add(
                        expr=m.plantYear[t].F_out[j, k, c] <= m.binOpt[j, k] * M[c]
                    )  # eqn. 8

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

        ## Calculate profit generated from byproduct valorization
        # only if specified by user
        if consider_byprod_val == True:
            m.plantYear[t].Byprod_Profit = pyo.Var(domain=pyo.Reals)
            m.plantYear[t].byprod_profit_con = pyo.Constraint(
                expr=m.plantYear[t].Byprod_Profit
                == sum(
                    m.plantYear[t].total_yearly_byprod[byprod] * byprod_vals[byprod]
                    for byprod in m.Byprods
                )
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

    ### environmental impacts
    # only consider if option is enabled by user
    if consider_environ_impacts == True:
        m.GWP = pyo.Var(domain=pyo.NonNegativeReals)
        m.GWP_cons = pyo.ConstraintList()

        m.GWP_cons.add(
            expr=m.GWP
            == sum(
                m.plantYear[t].total_yearly_GWP
                for t in pyo.RangeSet(prod_start, plant_end)
            )
        )
        # add in epsilon constraint
        m.epsilon_con = pyo.Constraint(expr=m.GWP <= epsilon)

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
