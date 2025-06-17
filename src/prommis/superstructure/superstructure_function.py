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
import warnings

import pyomo.environ as pyo


###################################################################################################
### Plant Lifetime Parameters
def check_plant_lifetime_params(plant_start, plant_lifetime):
    """
    This function checks that the lifetime parameters are feasible.

    Args:
        plant_lifetime: (int) The total lifetime of the plant, including plant construction. Must be at least three years.
    """
    ### Check types and structure.
    ## Check that plant_start is of type int.
    if not isinstance(plant_start, int):
        raise TypeError("plant_start is not of type int.")

    ## Check that plant_lifetime is of type int.
    if not isinstance(plant_lifetime, int):
        raise TypeError("plant_lifetime is not of type int.")

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
    ### Check types and structure.
    ## Check that available_feed is not an empty dict.
    if not available_feed:
        raise TypeError("available_feed dict is empty.")

    ## Check that structure of available_feed is correct.
    for key, value in available_feed.items():
        if not isinstance(key, int):
            raise TypeError(f"key {key} in available_feed dict is not of type int.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} in available_feed dict is not of type int or float."
            )

    ## Check that collection rate is of type int or float.
    if not isinstance(collection_rate, (int, float)):
        raise TypeError(f"collection_rate is not of type int or float.")

    ## Check that tracked_comps is a list.
    if not isinstance(tracked_comps, list):
        raise TypeError(f"tracked_comps is not of type list.")

    ## Check that all components in tracked_comps are of type str.
    for c in tracked_comps:
        if not isinstance(c, str):
            raise TypeError(f"Value {c} in tracked_comps is not of type str.")

    ## Check that prod_comp_mass is of type dict.
    if not isinstance(prod_comp_mass, dict):
        raise TypeError(f"prod_comp_mass is not of type dict.")

    ## Check that structure of prod_comp_mass is correct.
    for key, value in prod_comp_mass.items():
        if not isinstance(key, str):
            raise TypeError(f"key {key} in prod_comp_mass dict is not of type str.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} in prod_comp_mass dict is not of type int or float."
            )
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
    m, num_stages, options_in_stage, option_outlets, option_efficiencies
):
    """
    The function checks that the superstructure formulation parameter inputs are feasible.

    Args:
        m: pyomo model.
        num_stages: (int) Number of total stages.
        options_in_stage: (dict) Number of options in each stage.
        option_outlets: (dict) Set of options k' in stage j+1 connected to option k in stage j.
        option_efficiencies: (dict) Tracked component retention efficiency for each option.
    """
    ### Check types and structure.
    ## Check that num_stages is of type int.
    if not isinstance(num_stages, int):
        raise TypeError("num_stages is not of type int.")

    ## Check that options_in_stage is of type dict.
    if not isinstance(options_in_stage, dict):
        raise TypeError("options_in_stage is not of type dict.")

    ## Check that structure of options_in_stage is correct.
    for key, val in options_in_stage.items():
        if not isinstance(key, int):
            raise TypeError(f"key {key} in options_in_stage is not of type int.")
        if not isinstance(val, int):
            raise TypeError(f"value {val} in options_in_stage is not of type int.")
        if val <= 0:
            raise ValueError(
                f"number of options in stage {key} is infeasible. Number of options must be a positive int."
            )

    ## Check that option_outlets is of type dict.
    if not isinstance(option_outlets, dict):
        raise TypeError("option_outlets is not of type dict.")

    ## Check that structure of option_outlets is correct.
    for key, val in option_outlets.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in option_outlets is not of type tuple.")
        if not isinstance(val, list):
            raise TypeError(f"value {val} in option_outlets is not of type list.")
        if not val:
            raise ValueError(
                f"value for option {key} in option_outlets is an empty list."
            )
        if not all(type(x) is int and x > 0 for x in val):
            raise TypeError(
                f"value for option {key} in option_outlets incorrectly defined. All elements must be positive integers."
            )

    ## Check that option_efficiencies is of type dict.
    if not isinstance(option_efficiencies, dict):
        raise TypeError("option_efficiencies is not of type dict.")

    ## Check that structure of option_efficiencies is correct.
    for key, inner_dict in option_efficiencies.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in option_efficiencies is not of type tuple.")
        if not isinstance(inner_dict, dict):
            raise TypeError(
                f"value for option_efficiencies[{key}] is not of type dict."
            )
        for c, eff in inner_dict.items():
            if not isinstance(c, str):
                raise TypeError(
                    f"key {c} for inner dict of option_efficiency[{key}] is not of type str."
                )
            if not isinstance(eff, (int, float)):
                raise TypeError(
                    f"value for option_efficiencies[{key}][{c}] is not of type int or float."
                )
            if not ((0 <= eff) and (eff <= 1)):
                raise ValueError(
                    f"Efficiency for option {key}, component {c} is not in [0, 1]."
                )

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
    # Define a set for all the tracked components.
    tracked_comps = set(m.feed_params.tracked_comps.data())
    # Define a set for the keys of option_outlets.
    option_outlets_keys = set(option_outlets.keys())
    # Define a set for all the options in the superstructure (except for the last stage).
    initial_and_intermediate_opts_set = {
        (j, k) for j in range(1, num_stages) for k in range(1, options_in_stage[j] + 1)
    }
    # Define a set of all the keys of option_efficiencies.
    option_efficiencies_keys = set(option_efficiencies.keys())

    ### Define necessary data structures for tests.
    # Define a list to track the options in stage j which are not connected to an option in
    # stage j+1.
    missing_values = []
    # Define a list to track the options which don't have a connection from the previous
    # stage.
    disconnected_options = []
    # Define a list to track the missing efficiencies in the option_efficiencies dict.
    missing_effs = []
    # Set of options that don't have defined outlets (not in option_outlets).
    undefined_opt_outlets = []
    # Set of options whose outlets are specified, but aren't in the superstructure.
    undefined_opts = []
    # Create a list to track options for which efficiencies are not specified.
    opt_missing_efficiencies = []
    # Define a list to track the infeasible opts from option_efficiencies.
    infeasible_opts = []

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

    ## Check that outlets are defined for all options in the superstructure (except the last stage). i.e. check that each option in stage j is
    # connected to an option in stage j+1.
    # Find the options whose outlets are not defined.
    undefined_opt_outlets = initial_and_intermediate_opts_set - option_outlets_keys
    # Raise error if any are found.
    if undefined_opt_outlets:
        raise ValueError(
            f"The following options don't have defined outlets in option_outlets: {undefined_opt_outlets}"
        )

    ## Check that all options whose outlets are specified exist in the superstructure (as defined by options_in_stage).
    # Find the options that don't exist in the superstructure (as defined by options_in_stage).
    undefined_opts = option_outlets_keys - initial_and_intermediate_opts_set
    if undefined_opts:
        raise ValueError(
            f"The following options have defined outlets, but don't exist in the superstructure, as defined by options_in_stage: {undefined_opts}"
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
        raise ValueError(
            f"The following options are not connected from the previous stage: {disconnected_options}"
        )

    ## Check that all option outlets are feasible.
    for option, outlets in option_outlets.items():
        # Extract the stage number and option number for the current option.
        stage, option_number = option
        for outlet in outlets:
            # Construct the option in stage j+1 that the current option is connected to.
            next_opt = (stage + 1, outlet)
            # Check that this option is feasible. Raise an error if not.
            if next_opt not in all_opts_set:
                raise ValueError(
                    f"Option {next_opt} is not in the superstructure as defined by options_in_stage."
                )

    ## Check that option efficiencies are specified for all options in the superstructure.
    for opt in all_opts_set:
        if opt not in option_efficiencies_keys:
            opt_missing_efficiencies.append(opt)

    # Raise an error if there are any options for which no efficiencies are defined.
    if opt_missing_efficiencies:
        raise ValueError(
            f"Efficiences not specified for the following options in option_efficiencies: {opt_missing_efficiencies}."
        )

    ## Check that there are no options for which efficiencies are defined that are not in the superstructure.
    for opt in option_efficiencies_keys:
        if opt not in all_opts_set:
            infeasible_opts.append(opt)

    # Raise an error if there are any infeasible options in option_efficiencies.
    if infeasible_opts:
        raise ValueError(
            f"The following options for which efficiencies are specified in option_efficiencies do not exist in the superstructure: {infeasible_opts}"
        )

    ## Check that an option efficiency is defined for each tracked component in each option
    for j in range(1, num_stages + 1):
        for k in range(1, options_in_stage[j] + 1):
            # Check for missing efficiencies.
            option_efficiencies_tracked_comp_keys = set(
                option_efficiencies[(j, k)].keys()
            )
            # Keep track of the components for which an efficiency is not defined.
            missing = tracked_comps - option_efficiencies_tracked_comp_keys
            if missing:
                missing_effs.append((j, k, missing))
    # Raise an error if there are missing efficiencies.
    if missing_effs:
        msg = "Efficiencies not specified for all tracked components in the following options:\n"
        msg += "\n".join(
            f"  Option ({j}, {k}) missing components: {missing}"
            for j, k, missing in missing_effs
        )
        raise ValueError(msg)

    ## Check that all components in option_efficiencies are defined.
    for opt, inner_dict in option_efficiencies.items():
        for comp in inner_dict.keys():
            if comp not in tracked_comps:
                raise ValueError(
                    f"component {comp} for option {opt} in option_efficiencies is not defined as a tracked component."
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
        doc="Set containing all the options which utilize discrete units (discrete options). These are the options in the first stage, "
        "which utilize discrete units (or operators) to disassemble the incoming end-of-life products.",
    )
    m.supe_form_params.continuous_opts_set = pyo.Set(
        initialize=((opt) for opt in continuous_opts_set),
        doc="Set containing all the options in the stage which don't utilize discrete units (continuous options). "
        "All options after the first stage are assumed to be continuous. The sizing of these options is calculated "
        "based on the sum of the tracked components entering the option.",
    )
    m.supe_form_params.option_outlets = pyo.Param(
        m.supe_form_params.all_opts_set,
        initialize=option_outlets,
        doc="Defines the set of options k' in stage j+1 connected to option k in stage j.",
    )
    m.supe_form_params.option_efficiencies = pyo.Param(
        m.supe_form_params.all_opts_set,
        m.feed_params.tracked_comps,
        initialize=option_efficiencies_initialize,
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
    operators_per_discrete_unit,
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
        operators_per_discrete_unit: (dict) Number of operators needed per discrete unit for options that utilize discrete units.
        yearly_cost_per_unit: (dict) Yearly operating costs per unit for options which utilize discrete units.
        capital_cost_per_unit: (dict) Cost per unit for options which utilize discrete units.
        processing_rate: (dict) Processing rate per unit for options that utilize discrete units. In terms of units of incoming feed processed per year per unit.
        num_operators: (dict) Number of operators needed for each option.
        labor_rate: (float) Yearly wage per operator.
    """
    ### Check types and structure.
    ## Check that profit is of type dict.
    if not isinstance(profit, dict):
        raise TypeError("profit is not of type dict.")

    ## Check that structure of profit is correct.
    for key, inner_dict in profit.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in profit is not of type tuple.")
        if not isinstance(inner_dict, dict):
            raise TypeError(f"value for profit[{key}] is not of type dict.")
        for comp, value in inner_dict.items():
            if not isinstance(comp, str):
                raise TypeError(
                    f"key {comp} in inner dict of profit[{key}] is not of type str."
                )
            if not isinstance(value, (int, float)):
                raise TypeError(
                    f"value of profit[{key}][{comp}] is not of type int or float."
                )

    ## Check that structure of opt_var_oc_params is correct.
    for key, inner_dict in opt_var_oc_params.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in opt_var_oc_params is not of type tuple.")
        if not isinstance(inner_dict, dict):
            raise TypeError(f"value of opt_var_oc_params[{key}] is not of type dict.")
        for oc_param, value in inner_dict.items():
            if not isinstance(oc_param, str):
                raise TypeError(
                    f"key {oc_param} in opt_var_oc_params[{key}] is not of type str."
                )
            if not isinstance(value, (int, float)):
                raise TypeError(
                    f"value of opt_var_oc_params[{key}][{oc_param}] is not of type int or float."
                )

    ## Check that structure of operators_per_discrete_unit is correct.
    for key, value in operators_per_discrete_unit.items():
        if not isinstance(key, tuple):
            raise TypeError(
                f"key {key} in operators_per_discrete_unit is not of type tuple."
            )
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} of operators_per_discrete_unit[{key}] is not of type int or float."
            )

    ## Check that structure of yearly_cost_per_unit is correct.
    for key, value in yearly_cost_per_unit.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in yearly_cost_per_unit is not of type tuple.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} of yearly_cost_per_unit[{key}] is not of type int or float."
            )

    ## Check that structure of capital_cost_per_unit is correct.
    for key, value in capital_cost_per_unit.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in capital_cost_per_unit is not of type tuple.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} of capital_cost_per_unit[{key}] is not of type int or float."
            )

    ## Check that structure of processing_rate is correct.
    for key, value in processing_rate.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in processing_rate is not of type tuple.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} of processing_rate[{key}] is not of type int or float."
            )

    ## Check that structure of num_operators is correct.
    for key, value in num_operators.items():
        if not isinstance(key, tuple):
            raise TypeError(f"key {key} in num_operators is not of type tuple.")
        if not isinstance(value, (int, float)):
            raise TypeError(
                f"value {value} of num_operators[{key}] is not of type int or float."
            )

    ## Check that labor_rate is of type int or float.
    if not isinstance(labor_rate, (int, float)):
        raise TypeError("labor_rate is not of type int or float.")

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
    discrete_opts_set = set(m.supe_form_params.discrete_opts_set.data())
    # Define a set containing all the continuous options.
    continuous_opts_set = set(m.supe_form_params.continuous_opts_set.data())
    # Define a set containing the necessary variable operating parameters for all continuous options.
    var_oc_params = set(["a", "b"])
    # Define a set containing the keys of the operators_per_discrete_unit dict.
    operators_per_discrete_unit_keys = set(operators_per_discrete_unit.keys())
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
    # Define a list for tracking the options which define negative operators per discrete unit.
    negative_operators_per_discrete_units = []
    # Define a list for tracking the options which define a negative yearly cost per unit
    negative_yearly_cost_per_unit = []
    # Define a list for tracking the discrete options which define a negative capital cost per unit.
    negative_capital_cost_per_unit = []
    # Define a list for tracking the discrete options which define a non-positive processing rate for a unit.
    nonpositive_processing_rate = []
    # Define a list for tracking the continuous options which define a negative number of operators.
    negative_num_operators = []
    # Define a list to track the options for which this is not defined.
    undefined_profit_opts = []
    # Define a list to track the infeasible options.
    infeasible_profit_opts = []
    # Define a list to track the continuous options for which variable operating costs are not
    # defined.
    opts_missing_oc_params = []
    # Define a list to track the options for which variable operating costs are defined but aren't feasible.
    infeasible_opts_oc_params = []

    ### Run tests
    ## Check that profit per product is defined for all options in the final stage.
    undefined_profit_opts = final_opts_set - profit_opt_keys
    if undefined_profit_opts:
        raise ValueError(
            f"Profit not defined for the following options in the final stage: {undefined_profit_opts}"
        )

    ## Check that all options for which profits are defined are feasible (exist in the final stage of the superstructure).
    infeasible_profit_opts = profit_opt_keys - final_opts_set
    if infeasible_profit_opts:
        raise ValueError(
            f"The following options for which profits are defined do not exist in the final stage of the superstructure: {infeasible_profit_opts}."
        )

    ## Check that profit per product is in terms of the tracked components for all options in the final stage and are all nonnegative.
    for opt in m.supe_form_params.final_opts_set:
        # Check for missing profits.
        profit_tracked_comps = set(profit[opt].keys())
        # Keep track of the tracked components for which no profit per component is defined.
        missing = tracked_comps - profit_tracked_comps
        if missing:
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
            f"  Option ({j}, {k}) missing components: {missing}"
            for j, k, missing in missing_profit_comps
        )
        raise ValueError(msg)
    # Raise error if there are options with negative profits per tracked component defined.
    if negative_profit_comps:
        msg = "Profits some tracked components are listed as negative in the following options\n"
        msg += "\n".join(
            f"  Option ({j}, {k}) negative components: {negative}"
            for j, k, negative in negative_profit_comps
        )
        raise ValueError(msg)

    ## Check that all components in profit are defined tracked components.
    for opt, inner_dict in profit.items():
        for comp in inner_dict.keys():
            if comp not in tracked_comps:
                raise ValueError(
                    f"component: {comp} as defined in profit for option {opt} is not a tracked component."
                )

    ## Check that variable operating cost params are defined for all continuous options.
    opts_missing_oc_params = continuous_opts_set - opt_var_oc_params_keys
    # Raise an error if there are opts missing oc params.
    if opts_missing_oc_params:
        raise ValueError(
            f"Variable operating cost params not defined for the following continuous options: {opts_missing_oc_params}."
        )

    ## Check that all all options for which variable operating costs params are defined are feasible.
    infeasible_opts_oc_params = opt_var_oc_params_keys - continuous_opts_set
    # Raise an error if there are infeasible options found.
    if infeasible_opts_oc_params:
        raise ValueError(
            f"The following options for which variable operating costs params were defined do not exist as continuous options in the superstructure: {infeasible_opts_oc_params}."
        )

    ## Check that both necessary variable operating cost parameters ('a' and 'b') are defined for all continuous options.
    for opt in continuous_opts_set:
        opt_params_set = set(opt_var_oc_params[opt].keys())
        if opt_params_set != var_oc_params:
            raise ValueError(
                f"Variable operating cost params incorrectly defined in opt_var_oc_params for option: {opt}."
            )

    ## Check that operators per discrete unit is defined for all options that utilize discrete units.
    if operators_per_discrete_unit_keys != discrete_opts_set:
        raise ValueError(
            "operators_per_discrete_unit must be defined, and only defined, for all discrete options."
        )

    ## Check that operators per discrete unit are all defined to be non-negative.
    negative_operators_per_discrete_units = [
        opt for opt, w in operators_per_discrete_unit.items() if w < 0
    ]
    # If there are negative operators defined per discrete unit for any discrete options, raise an error.
    if negative_operators_per_discrete_units:
        raise ValueError("Operators per discrete unit must all be non-negative.")

    ## Check that yearly cost per unit is defined for all discrete options.
    if yearly_cost_per_unit_keys != discrete_opts_set:
        raise ValueError(
            "yearly_cost_per_unit must be defined, and only defined for, all discrete options."
        )

    ## Check that the yearly cost per unit values are all non-negative.
    negative_yearly_cost_per_unit = [
        opt for opt, cost in yearly_cost_per_unit.items() if cost < 0
    ]
    # If there any discrete options that define a negative yearly cost per unit, raise an error.
    if negative_yearly_cost_per_unit:
        raise ValueError("yearly_cost_per_unit values must all be non-negative.")

    ## Check that cost per unit is defined for all discrete options.
    if capital_cost_per_unit_set != discrete_opts_set:
        raise ValueError(
            "capital_cost_per_unit must be defined, and only defined for, all discrete options."
        )

    ## Check that the cost per unit values are all non-negative.
    negative_capital_cost_per_unit = [
        opt for opt, cost in capital_cost_per_unit.items() if cost < 0
    ]
    # If there are any discrete options which define a negative capital cost per unit, raise an error.
    if negative_capital_cost_per_unit:
        raise ValueError("capital_cost_per_unit values must all be non-negative.")

    ## Check that prosessing rate is defined for all discrete options.
    if processing_rate_keys != discrete_opts_set:
        raise ValueError(
            "processing_rate must be defined, and only defined for, all discrete options."
        )

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
            "num_operators must be defined, and only defined for, all continuous options."
        )

    ## Check that num_operators values are all non-negative
    negative_num_operators = [
        opt for opt, operators in num_operators.items() if operators < 0
    ]
    if negative_num_operators:
        raise ValueError("number of operators must all be non-negative.")

    ## Check that labor rate is non-negative
    if labor_rate < 0:
        raise ValueError("labor rate must be non-negative.")


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
            m.feed_params.max_feed_entering / processing_rate[key]
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
    m.operating_params.operators_per_discrete_unit = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=operators_per_discrete_unit,
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
        doc="The processing rate per discrete unit for discrete options. In terms of number of end-of-life products "
        "disassembled per discrete unit per year.",
    )
    m.operating_params.num_operators = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=num_operators,
        doc="The number of operators per continuous options.",
    )
    m.operating_params.labor_rate = pyo.Param(
        initialize=labor_rate, doc="The yearly wage per operator."
    )
    m.operating_params.discrete_units_per_option = pyo.Param(
        m.supe_form_params.discrete_opts_set,
        initialize=discrete_units_per_option,
        doc="The number of discrete units per option needed to disassemble all the incoming end-of-life products over the operational "
        "lifetime of the plant.",
    )


###################################################################################################
### Costing Parameters
def check_discretized_costing_params(m, discretized_purchased_equipment_cost):
    """
    This function checks that all the costing parameters are feasible.

    Args:
        m: pyomo model.
        discretized_purchased_equipment_cost: (dict) Discretized cost by flows entering for each continuous option.
    """
    ### Check types and structure.
    ## Check that discretized_purchased_equipment_cost is of type dict.
    if not isinstance(discretized_purchased_equipment_cost, dict):
        raise TypeError("discretized_purchased_equipment_cost is not of type dict.")

    ## Check structure of discretized_purchased_equipment_cost.
    for opt, inner_dict in discretized_purchased_equipment_cost.items():
        if not isinstance(opt, tuple):
            raise TypeError(
                f"opt {opt} in discretized_purchased_equipment_cost is not of type tuple."
            )
        if not isinstance(inner_dict, dict):
            raise TypeError(
                f"discretized_purchased_equipment_cost[{opt}] is not of type dict."
            )
        if not inner_dict:
            raise ValueError(
                f"discretized_purchased_equipment_cost[{opt}] is an empty dict."
            )
        # For each continuous option, discretized flowrates ("Flowrates"), followed by discretized costs ("Costs") must be provided.
        if [key for key in inner_dict.keys()] != ["Flowrates", "Costs"]:
            msg = "For each option in discretized_purchased_equipment_cost, discretized flowrates, followed by discretized costs must be provided. "
            msg += "This must be of the form: 'Flowrates': [list of data points], 'Costs': [list of data points]. "
            msg += f"This was not the case for the option: {opt}."
            raise ValueError(msg)
        for data_type, data in inner_dict.items():
            if not isinstance(data, list):
                raise TypeError(
                    f"Discretized data for option: {opt}, data type: {data_type} is not of type list."
                )
            if not data:
                raise ValueError(
                    f"Empy list passed for option: {opt}, data type: {data_type}."
                )

    ### Define parameters necessary for tests
    # Define a set of all the keys in the discretized_purchased_equipment_cost dict.
    discretized_purchased_equipment_cost_opts_set = set(
        discretized_purchased_equipment_cost.keys()
    )

    ### Define necessary data structures for tests.
    # Define a list of continuous options which are missing discretized capex data.
    missing_continuous_opts = []
    # Define a list for tracking the discrete options for which discretized capex data is defined for.
    discrete_opts = []
    # Define a list for tracking the options in which inconsistent discretized capex is defined.
    inconsistent_data_point_opts = []
    # Define a list for tracking the infeasible opts for which discretized capex is defined.
    infeasible_opts = []

    ### Run tests
    ## Check that discretized capex provided for all continuous options
    # and check that discretized capex not provided for options that utilize discrete units.
    for opt in m.supe_form_params.continuous_opts_set:
        # Keep track of the continuous opts for which discretized data is not defined for.
        if opt not in discretized_purchased_equipment_cost_opts_set:
            missing_continuous_opts.append(opt)
    # Check that discretized capex not provided for discrete opts.
    for opt in m.supe_form_params.discrete_opts_set:
        # Keep track of the discrete options for which discretized capex data is defined for.
        if opt in discretized_purchased_equipment_cost_opts_set:
            discrete_opts.append(opt)
    # Check all options in discretized capex are feasible.
    for opt in discretized_purchased_equipment_cost_opts_set:
        # Keep track of the infeasible opts
        if opt not in m.supe_form_params.continuous_opts_set:
            infeasible_opts.append(opt)
    # Raise error if discretized capex is not provided for all continuous options.
    if missing_continuous_opts:
        raise ValueError(
            f"discretized_purchased_equipment_cost is missing values for the following continuous options: {missing_continuous_opts}. "
        )
    # Raise error if discretized capex provided for any discrete options.
    if discrete_opts:
        raise ValueError(
            f"discretized_purchased_equipment_cost contains values for the following discrete options: {discrete_opts}. "
        )
    # Raise error if discretized capex provided for any infeasible options.
    if infeasible_opts:
        raise ValueError(
            f"discretized_purchased_equipment_cost contains values for the following options which don't exist as continuous options in the superstructure: {infeasible_opts}."
        )

    ## Check that all options have the same number of discretized data points for flows entering and costs.
    for opt in m.supe_form_params.continuous_opts_set:
        # Keep track of the number of flowrate data points defined for the option.
        opt_num_flow_data_points = len(
            discretized_purchased_equipment_cost[opt]["Flowrates"]
        )
        # Keep track of the number of cost data points defined for the option.
        opt_num_cost_data_points = len(
            discretized_purchased_equipment_cost[opt]["Costs"]
        )
        # Keep track of the options for which the number of flowrate and cost data points are not the same.
        if opt_num_flow_data_points != opt_num_cost_data_points:
            inconsistent_data_point_opts.append(opt)
    # Raise error if there are options with an inconsistent number of data points within discretized_purchased_equipment_cost.
    if inconsistent_data_point_opts:
        raise ValueError(
            f"Inconsistent number of data points for Flowrates and Costs within discretized_purchased_equipment_cost for the following options: {inconsistent_data_point_opts}. "
        )


def add_discretized_costing_params(m, discretized_purchased_equipment_cost):
    """
    This function adds all the costing parameters to a block.

    Args:
        m: pyomo model.
        discretized_purchased_equipment_cost: (dict) Discretized cost by flows entering for each continuous option
    """
    ### Define parameters from user input.
    # Define a dict to hold discretized flowrate data for each option.
    flowrates_data = {}
    for opt in m.supe_form_params.continuous_opts_set:
        flowrates_data[opt] = discretized_purchased_equipment_cost[opt]["Flowrates"]
    # Define a dict to hold discretized cost data for each option.
    costs_data = {}
    for opt in m.supe_form_params.continuous_opts_set:
        costs_data[opt] = discretized_purchased_equipment_cost[opt]["Costs"]

    ### Define necessary pyomo parameters.
    ## Create a block
    m.discretized_costing_params = pyo.Block(doc="Block to hold costing parameters.")
    ## Pyomo parameters
    m.discretized_costing_params.flowrates_data = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=flowrates_data,
        doc="Discretized flowrate data for all continuous options.",
    )
    m.discretized_costing_params.costs_data = pyo.Param(
        m.supe_form_params.continuous_opts_set,
        initialize=costs_data,
        doc="Discretized costing data for all continuous options.",
    )


###################################################################################################
### Mass Balances
def add_mass_balance_params(m):
    """
    This function builds the mass balance parameters.

    Args:
        m: pyomo model.
    """
    ## Define parameters
    # Calculate Big-M values for each tracked component.
    m_val = {}
    for c in m.feed_params.tracked_comps:
        # Max value assumes 100% efficiency (no losses) over the plant's lifetime. Rounded up to ensure validity
        # in constraints.
        m_val[c] = math.ceil(
            m.feed_params.max_feed_entering * m.feed_params.prod_comp_mass[c]
        )

    ## Create block
    m.mass_balances = pyo.Block(doc="Block to hold mass balances.")

    ## Pyomo parameters
    # Flow entering each stage (except the last stage).
    m.mass_balances.f_stages = pyo.RangeSet(
        1,
        m.supe_form_params.num_stages - 1,
        doc="Set of all stages except the last. Used to define the flow variable: 'f'.",
    )
    m.mass_balances.flow_set = pyo.Set(
        initialize=m.supe_form_params.all_opts_set
        * m.feed_params.tracked_comps
        * m.plant_lifetime_params.operational_range,
        doc="Set of all options, tracked components, and operational years of the plant. Used to define the flow variables: 'f_in' and 'f_out'.",
    )
    m.mass_balances.big_m_val = pyo.Param(
        m.feed_params.tracked_comps,
        initialize=m_val,
        doc="Big-M parameters used in Equations (7) and (8) from the documentation.",
    )


def add_mass_balance_vars(m):
    """
    This function builds the mass balance variables.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.mass_balances.f = pyo.Var(
        m.mass_balances.f_stages
        * m.feed_params.tracked_comps
        * m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each stage (except the last stage). See documentation for more details.",
    )
    m.mass_balances.f_in = pyo.Var(
        m.mass_balances.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )
    m.mass_balances.f_out = pyo.Var(
        m.mass_balances.flow_set,
        domain=pyo.NonNegativeReals,
        doc="Flow entering each option. See documentation for more details.",
    )
    m.mass_balances.option_binary_var = pyo.Var(
        m.supe_form_params.all_opts_set,
        domain=pyo.Binary,
        doc="Binary variables to indicate whether or not an option has been selected.",
    )


def add_mass_balance_cons(m):
    """
    This function builds the mass balance constraints.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    @m.mass_balances.Constraint(
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
            return m.mass_balances.f[j - 1, c, t] == sum(
                m.mass_balances.f_in[j, k, c, t] for k in num_options
            )

    @m.mass_balances.Constraint(
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
            ] == sum(m.mass_balances.f_in[j, k, c, t] for k in num_options)
        else:
            return pyo.Constraint.Skip

    @m.mass_balances.Constraint(
        m.mass_balances.flow_set, doc="Equation (3) from the documentation."
    )
    def intermediate_flow_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        alpha = m.supe_form_params.option_efficiencies[j, k, c]
        return (
            m.mass_balances.f_in[j, k, c, t] * alpha
            == m.mass_balances.f_out[j, k, c, t]
        )

    @m.mass_balances.Constraint(
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
                sum(m.mass_balances.f_out[j, k, c, t] for k in num_options)
                == m.mass_balances.f[j, c, t]
            )
        else:
            return pyo.Constraint.Skip

    @m.mass_balances.Constraint(
        m.supe_form_params.stages_set, doc="Equation (5) from the documentation."
    )
    def stage_binary_cons(b, j):
        m = b.model()  # Get the main model
        # Extract all the options available in stage 'j'
        num_options = range(1, m.supe_form_params.options_in_stage[j] + 1)
        return sum(m.mass_balances.option_binary_var[j, k] for k in num_options) == 1

    @m.mass_balances.Constraint(
        m.supe_form_params.all_opts_set, doc="Equation (6) from the documentation."
    )
    def connection_binary_cons(b, j, k):
        m = b.model()  # Get the main model
        if j != m.supe_form_params.num_stages:
            # Extract the set of options k' in stage j+1 connected to option k in stage j.
            opt_connections = m.supe_form_params.option_outlets[j, k]
            return m.mass_balances.option_binary_var[j, k] <= sum(
                m.mass_balances.option_binary_var[j + 1, kp] for kp in opt_connections
            )
        else:
            return pyo.Constraint.Skip

    @m.mass_balances.Constraint(
        m.mass_balances.flow_set, doc="Equation (7) from the documentation."
    )
    def f_in_big_m_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        return (
            m.mass_balances.f_in[j, k, c, t]
            <= m.mass_balances.option_binary_var[j, k] * m.mass_balances.big_m_val[c]
        )

    @m.mass_balances.Constraint(
        m.mass_balances.flow_set, doc="Equation (8) from the documentation."
    )
    def f_out_big_m_cons(b, j, k, c, t):
        m = b.model()  # Get the main model
        return (
            m.mass_balances.f_out[j, k, c, t]
            <= m.mass_balances.option_binary_var[j, k] * m.mass_balances.big_m_val[c]
        )


###################################################################################################
### Costing
def add_costing_params(m):
    """
    This function builds the costing parameters.

    Args:
        m: pyomo model.
    """
    ## Define parameters
    # Create an upper bound for the total flow entering a continuous option.
    max_flow_upper_bound = sum(
        m.mass_balances.big_m_val[c] for c in m.feed_params.tracked_comps
    )

    ## Create blocks
    m.costing = pyo.Block(doc="Block to hold costing.")
    m.net_present_value = pyo.Block(doc="Block for holding net present value.")
    m.cost_of_recovery = pyo.Block(doc="Block for holding cost of recovery.")

    ## Pyomo parameters
    m.costing.max_flow_upper_bound = pyo.Param(
        initialize=max_flow_upper_bound,
        doc="Upper bound for the total flow entering a continuous option.",
    )
    m.costing.lang_factor = pyo.Param(initialize=2.97)
    m.costing.i_operating_expense_escalation = pyo.Param(
        initialize=0.03, doc="Operating expenses escalation rate."
    )
    m.costing.i_capital_escalation = pyo.Param(
        initialize=0.036, doc="Capital expenses escalation rate."
    )
    m.costing.discount_factor = pyo.Param(
        initialize=0.0577, doc="Discount factor for calculate the net present value."
    )


def add_costing_vars(m):
    """
    This function builds the costing variables.

    Args:
        m: pyomo model.
    """
    ## Pyomo variables
    m.costing.opt_profit = pyo.Var(
        m.supe_form_params.final_opts_set,
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The profit generated by each option in the final processing stage each year.",
    )
    m.costing.total_profit = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total profit generated by the plant each year.",
    )
    m.costing.total_byproduct_profit = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.Reals,
        doc="The total profit generated by the plant each year from the valorization of byproducts.",
    )
    m.costing.npv = pyo.Var(domain=pyo.Reals, doc="The net present value.")
    m.cost_of_recovery.cor = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The cost of recovery."
    )
    m.costing.purchased_equipment_cost = pyo.Var(
        m.supe_form_params.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The cost of purchased equipment.",
    )
    m.costing.flow_entering = pyo.Var(
        m.supe_form_params.continuous_opts_set,
        bounds=(0, m.costing.max_flow_upper_bound),
        doc="The max total flow that enters each continuous option over the lifetime"
        "of the plant.",
    )
    m.costing.total_plant_cost = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The total plant cost."
    )
    m.costing.financing = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The total financing cost of the plant."
    )
    m.costing.other_costs = pyo.Var(
        domain=pyo.NonNegativeReals, doc="'Other costs' associated with the plant."
    )
    m.costing.total_overnight_cost = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The total overnight cost of the plant."
    )
    m.costing.opt_var_operating_expense = pyo.Var(
        m.supe_form_params.all_opts_set * m.plant_lifetime_params.operational_range,
        domain=pyo.Reals,
        doc="Yearly variable operating expense for each option.",
    )
    m.costing.total_var_operating_expense = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.Reals,
        doc="Total yearly variable operating expense.",
    )
    m.costing.opt_operators = pyo.Var(
        m.supe_form_params.all_opts_set,
        domain=pyo.NonNegativeReals,
        doc="The number of operators needed for each option.",
    )
    m.costing.total_operators = pyo.Var(
        domain=pyo.NonNegativeIntegers,
        doc="The total number of operators needed for the process. Must be an integer value.",
    )
    m.costing.cost_of_labor = pyo.Var(
        domain=pyo.NonNegativeReals, doc="The cost of labor for the process."
    )
    m.costing.m_and_sm = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Maintenance & Supply Materials (M&SM)."
    )
    m.costing.sa_and_qa_qc = pyo.Var(
        domain=pyo.NonNegativeReals,
        doc="Sample Analysis & Quality Assurance/Quality Control (SA&QA/QC).",
    )
    m.costing.s_ip_r_and_d = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Sales, Intellectual Property, and Research & Development (S,IP,R&D).",
    )
    m.costing.a_and_sl = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Administrative & Supporting Labor (A&SL)."
    )
    m.costing.fb = pyo.Var(domain=pyo.NonNegativeReals, doc="Fringe Benefits (FB).")
    m.costing.pt_and_i = pyo.Var(
        domain=pyo.NonNegativeReals, doc="Property Taxes & Insurance (PT&I)."
    )
    m.costing.total_fixed_operating_expense = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="Total yearly fixed operating expense.",
    )
    m.costing.total_overnight_cost_expended = pyo.Var(
        m.plant_lifetime_params.plant_life_range,
        domain=pyo.NonNegativeReals,
        doc="The total overnight cost expended each year.",
    )
    m.costing.plant_overhead = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The yearly plant overhead.",
    )
    m.costing.total_operating_expense = pyo.Var(
        m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The total operating expense each year.",
    )
    m.costing.cash_flow = pyo.Var(
        m.plant_lifetime_params.plant_life_range,
        domain=pyo.Reals,
        doc="The yearly cash flow.",
    )


def add_costing_cons(m):
    """
    This function builds the costing constraints.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    @m.net_present_value.Constraint(
        m.supe_form_params.final_opts_set,
        m.plant_lifetime_params.operational_range,
        doc="Calculates the profit from each option each year in the final stage when the net present value is selected as "
        "the objective function.",
    )
    def calculate_final_opts_profit_cons(b, j, k, t):
        return m.costing.opt_profit[j, k, t] == sum(
            m.mass_balances.f_out[j, k, c, t] * m.operating_params.profit[j, k, c]
            for c in m.feed_params.tracked_comps
        )

    @m.net_present_value.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total yearly profit when the net present value is selected as the objective function.",
    )
    def calculate_total_profit_cons(b, t):
        return (
            m.costing.total_profit[t]
            == sum(
                m.costing.opt_profit[opt, t]
                for opt in m.supe_form_params.final_opts_set
            )
            + m.costing.total_byproduct_profit[t]
        )

    @m.cost_of_recovery.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total yearly profit when the cost of recovery is selected as the objective function.",
    )
    def calculate_total_profit_cons(b, t):
        return (
            m.costing.total_profit[t]
            == m.cost_of_recovery.cor
            * sum(
                sum(
                    m.mass_balances.f_out[opt, c, t]
                    for c in m.feed_params.tracked_comps
                )
                for opt in m.supe_form_params.final_opts_set
            )
            + m.costing.total_byproduct_profit[t]
        )

    @m.costing.Constraint(
        m.supe_form_params.discrete_opts_set,
        doc="Calculates the purchased cost of equipment for all discrete options. "
        "Done by multiplying the number of discrete units by the capital cost per unit.",
    )
    def discrete_opts_purchased_equipment_cost_cons(b, j, k):
        m = b.model()  # Get the main model
        return (
            m.costing.purchased_equipment_cost[j, k]
            == m.operating_params.discrete_units_per_option[j, k]
            * m.operating_params.capital_cost_per_unit[j, k]
        )

    @m.costing.Constraint(
        m.supe_form_params.continuous_opts_set,
        m.plant_lifetime_params.operational_range,
        doc="Constraint to determine the max flow entering each continuous option over the lifetime of the plant.",
    )
    def max_flow_entering_cons(b, j, k, t):
        m = b.model()  # Get the main model
        return m.costing.flow_entering[j, k] >= sum(
            m.mass_balances.f_in[j, k, c, t] for c in m.feed_params.tracked_comps
        )

    def piecewise_rule(b, j, k):
        m = b.model()  # Get the main model

        flow_data = m.discretized_costing_params.flowrates_data[j, k]
        purchased_equipment_cost_data = m.discretized_costing_params.costs_data[j, k]

        # use m.add_component to generate all piecewise functions
        # piecewise = Piecewise(yval, xval, *kwargs)
        piecewise = pyo.Piecewise(
            m.costing.purchased_equipment_cost[j, k],
            m.costing.flow_entering[j, k],
            pw_pts=flow_data,
            pw_constr_type="EQ",
            f_rule=purchased_equipment_cost_data,
            pw_repn="SOS2",
        )
        b.add_component("Option_" + str((j, k)) + "_Piecewise_Constraint", piecewise)

    m.costing.piecewise_cons = pyo.Block(
        m.supe_form_params.continuous_opts_set,
        rule=piecewise_rule,
        doc="This block holds all the piecewise constraints for calculating the "
        "purchased equipment costs for all the continuous option.",
    )

    @m.costing.Constraint(
        doc="Calculates the total plant cost of the plant. See Equation (x) in the documentation."
    )
    def calculate_total_plant_cost_con(b):
        m = b.model()  # Get the main model
        return m.costing.total_plant_cost == sum(
            m.costing.purchased_equipment_cost[opt]
            for opt in m.supe_form_params.discrete_opts_set
        ) + m.costing.lang_factor * sum(
            m.costing.purchased_equipment_cost[opt]
            for opt in m.supe_form_params.continuous_opts_set
        )

    @m.costing.Constraint(
        doc="Calculates the total financing cost of the plant. Assumed to be 2.7% of the total plant cost."
    )
    def calculate_financing_cost_con(b):
        m = b.model()  # Get the main model
        return m.costing.financing == 0.027 * m.costing.total_plant_cost

    @m.costing.Constraint(
        doc="Calculates the 'other costs' of the plant. Assumed to be 15% of the total plant cost."
    )
    def calculate_other_costs_con(b):
        m = b.model()  # Get the main model
        return m.costing.other_costs == 0.15 * m.costing.total_plant_cost

    @m.costing.Constraint(
        doc="Calculates the total overnight cost of the plant. "
        "Equal to the total plant cost + financing + 'other costs'."
    )
    def calculate_total_overnight_cost_con(b):
        m = b.model()  # Get the main model
        return (
            m.costing.total_overnight_cost
            == m.costing.total_plant_cost + m.costing.financing + m.costing.other_costs
        )

    @m.costing.Constraint(
        m.supe_form_params.all_opts_set,
        m.plant_lifetime_params.operational_range,
        doc="Calculates yearly operating expense for each option.",
    )
    def calculate_opt_yearly_variable_expense_cons(b, j, k, t):
        if (j, k) in m.supe_form_params.discrete_opts_set:
            return (
                m.costing.opt_var_operating_expense[j, k, t]
                == m.operating_params.discrete_units_per_option[j, k]
                * m.operating_params.yearly_cost_per_unit[j, k]
                * m.mass_balances.option_binary_var[j, k]
            )
        else:
            return (
                m.costing.opt_var_operating_expense[j, k, t]
                == m.operating_params.opt_var_oc_params[j, k, "a"]
                * sum(
                    m.mass_balances.f_in[j, k, c, t]
                    for c in m.feed_params.tracked_comps
                )
                + m.operating_params.opt_var_oc_params[j, k, "b"]
                * m.mass_balances.option_binary_var[j, k]
            )

    @m.costing.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total yearly variable operating cost.",
    )
    def calculate_total_yearly_variable_operating_costs_cons(b, t):
        return m.costing.total_var_operating_expense[t] == sum(
            m.costing.opt_var_operating_expense[opt, t]
            for opt in m.supe_form_params.all_opts_set
        )

    @m.costing.Constraint(
        m.supe_form_params.all_opts_set,
        doc="Calculate the number of operators needed for each option.",
    )
    def calculate_operators_per_option_cons(b, j, k):
        if (j, k) in m.supe_form_params.discrete_opts_set:
            return (
                m.costing.opt_operators[j, k]
                == m.operating_params.discrete_units_per_option[j, k]
                * m.operating_params.operators_per_discrete_unit[j, k]
                * m.mass_balances.option_binary_var[j, k]
            )
        else:
            return (
                m.costing.opt_operators[j, k]
                == m.operating_params.num_operators[j, k]
                * m.mass_balances.option_binary_var[j, k]
            )

    @m.costing.Constraint(
        doc="Calculate the number of operators needed for the entire process."
    )
    def calculate_total_operators_cons(b):
        return m.costing.total_operators >= sum(
            m.costing.opt_operators[opt] for opt in m.supe_form_params.all_opts_set
        )

    @m.costing.Constraint(doc="Calculates the cost of labor for the process.")
    def calculate_cost_of_labor_con(b):
        return (
            m.costing.cost_of_labor
            == m.costing.total_operators * m.operating_params.labor_rate
        )

    @m.costing.Constraint(
        doc="Calculate M&SM. Assumed to be 2% of the total plant cost."
    )
    def calculate_m_and_sm_con(b):
        return m.costing.m_and_sm == 0.02 * m.costing.total_plant_cost

    @m.costing.Constraint(
        doc="Calculate SA&QA/QC. Assumed to be 10% of the cost of labor."
    )
    def calculate_sa_and_qa_qc_con(b):
        return m.costing.sa_and_qa_qc == 0.1 * m.costing.cost_of_labor

    @m.costing.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculate S,IP,R&D. Assumed to be 1% of the revenue each year.",
    )
    def calculate_s_ip_r_and_d_con(b, t):
        return m.costing.s_ip_r_and_d[t] == 0.01 * m.costing.total_profit[t]

    @m.costing.Constraint(doc="Calculate A&SL. Assumed to be 20% of the cost of labor.")
    def calculate_a_and_sl_con(b):
        return m.costing.a_and_sl == 0.2 * m.costing.cost_of_labor

    @m.costing.Constraint(doc="Calculate FB. Assumed to be 25% of the cost of labor.")
    def calculate_fb_con(b):
        return m.costing.fb == 0.25 * m.costing.cost_of_labor

    @m.costing.Constraint(
        doc="Calculate PT&I. Assumed to be 1% of the total plant cost."
    )
    def calculate_pt_and_i_con(b):
        return m.costing.pt_and_i == 0.01 * m.costing.total_plant_cost

    @m.costing.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total yearly fixed operating cost.",
    )
    def calculate_total_yearly_fixed_operating__costs_cons(b, t):
        return (
            m.costing.total_fixed_operating_expense[t]
            == m.costing.cost_of_labor
            + m.costing.m_and_sm
            + m.costing.sa_and_qa_qc
            + m.costing.s_ip_r_and_d[t]
            + m.costing.a_and_sl
            + m.costing.fb
            + m.costing.pt_and_i
        )

    @m.costing.Constraint(
        m.plant_lifetime_params.plant_life_range,
        doc="Calculates the total overnight cost expended in a given year. It is assumed that 10% is expended "
        "in the first year, 60% in the second year, and 30% in the third year.",
    )
    def calculate_total_overnight_cost_expended_cons(b, t):
        if t == m.plant_lifetime_params.plant_start:
            return (
                m.costing.total_overnight_cost_expended[t]
                == 0.1 * m.costing.total_overnight_cost
            )
        elif t == (m.plant_lifetime_params.plant_start + 1):
            return (
                m.costing.total_overnight_cost_expended[t]
                == 0.6 * m.costing.total_overnight_cost
            )
        elif t == (m.plant_lifetime_params.plant_start + 2):
            return (
                m.costing.total_overnight_cost_expended[t]
                == 0.3 * m.costing.total_overnight_cost
            )
        else:
            return m.costing.total_overnight_cost_expended[t] == 0

    @m.costing.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the plant overhead each year. Assumed to be 20% of the fixed and variable "
        "operating costs.",
    )
    def calculate_plant_overhead_cons(b, t):
        return m.costing.plant_overhead[t] == 0.2 * (
            m.costing.total_var_operating_expense[t]
            + m.costing.total_fixed_operating_expense[t]
        )

    @m.costing.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total operating expenses each year. Includes variable, fixed, and plant overhead.",
    )
    def calculate_total_operating_expense_cons(b, t):
        return (
            m.costing.total_operating_expense[t]
            == m.costing.total_var_operating_expense[t]
            + m.costing.total_fixed_operating_expense[t]
            + m.costing.plant_overhead[t]
        )

    @m.costing.Constraint(
        m.plant_lifetime_params.plant_life_range,
        doc="Calculate the cash flow for each year of the plant's lifetime.",
    )
    def calculate_yearly_costing(b, t):
        if t == m.plant_lifetime_params.plant_start:
            return m.costing.cash_flow[t] == -m.costing.total_overnight_cost_expended[t]
        else:
            return m.costing.cash_flow[t] == (
                m.costing.total_profit[t] - m.costing.total_operating_expense[t]
            ) * (1 + m.costing.i_operating_expense_escalation) ** (
                t - m.plant_lifetime_params.plant_start
            ) - m.costing.total_overnight_cost_expended[
                t
            ] * (
                1 + m.costing.i_capital_escalation
            ) ** (
                t - m.plant_lifetime_params.plant_start
            )

    @m.costing.Constraint(doc="Calculates the net present value of the plant.")
    def calculate_net_present_value_con(b):
        return m.costing.npv == sum(
            m.costing.cash_flow[t]
            / (
                (1 + m.costing.discount_factor)
                ** (t - m.plant_lifetime_params.plant_start)
            )
            for t in m.plant_lifetime_params.plant_life_range
        )

    @m.cost_of_recovery.Constraint(
        doc="Sets the net present value to zero when the cost of recovery is chosen as an "
        "objective function."
    )
    def set_net_present_value_to_zero_con(b):
        return m.costing.npv == 0

    ## Set objective function
    m.net_present_value.obj = pyo.Objective(
        expr=m.costing.npv,
        sense=pyo.maximize,
        doc="Objective function is maximizing the net present value when the NPV is chosen.",
    )
    m.cost_of_recovery.obj = pyo.Objective(
        expr=m.cost_of_recovery.cor,
        sense=pyo.minimize,
        doc="Objective function is minimizing the cost of recovery when the COR is chosen.",
    )


###################################################################################################
### Environmental impacts
def check_environmental_impact_params(
    m, consider_environmental_impacts, options_environmental_impacts, epsilon
):
    """
    This function checks that all the environmental impact parameters are feasible.

    Args:
        m: pyomo model.
        consider_environmental_impacts: (bool) Choice of whether or not to consider environmental impacts.
        options_environmental_impacts: (dict) Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.
        epsilon: (float) Epsilon factor for generating the Pareto front.
    """
    ## Check that consider_environmental_impacts is of type bool. Must be provided regardles on whether or not environmental impacts are considered.
    if not isinstance(consider_environmental_impacts, bool):
        raise TypeError("consider_environmental_impacts is not of type bool.")
    ### Only need to check feasibility of the rest of the parameters if user wants to consider environmental impacts
    # (consider_environmental_impacts is True).
    if consider_environmental_impacts:
        ### Check types and structure.
        ## Check that options_environmental_impacts is of type dict.
        if not isinstance(options_environmental_impacts, dict):
            raise TypeError("options_environmental_impacts is not of type dict.")

        ## Check that structure of options_environmental_impacts is correct.
        for opt, val in options_environmental_impacts.items():
            if not isinstance(opt, tuple):
                raise TypeError(
                    f"key {opt} in options_environmental_impacts is not of type tuple."
                )
            if not isinstance(val, (int, float)):
                raise TypeError(
                    f"value of options_environmental_impacts[{opt}] is not of type int or float."
                )

        ## Check that epsilon is of type int or float.
        if not isinstance(epsilon, (int, float)):
            raise TypeError("epsilon is not of type int or float.")

        ### Define parameters necessary for tests.
        # Define a set for all the options in the superstructure.
        all_opts_set = set(m.supe_form_params.all_opts_set.data())
        # Define a set for all the keys in the consider_environmental_impacts dict.
        options_environmental_impacts_keys_set = set(
            options_environmental_impacts.keys()
        )

        ### Run tests
        ## Check that environmental impacts matrix contains entry for all options in the superstructure.
        if all_opts_set != options_environmental_impacts_keys_set:
            msg = "options_environmental_impacts must provide values for all options in the superstructure. "
            msg += "Please check that environmental impacts are defined for all options, and that there aren't any defined for infeasible options (options that don't exist within th superstructure)."
            raise ValueError(msg)


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
    m.environmental_impacts = pyo.Block(doc="Block to hold environmental impacts.")

    ## Pyomo parameters
    m.environmental_impacts.consider_environmental_impacts = pyo.Param(
        initialize=consider_environmental_impacts,
        within=pyo.Boolean,
        doc="Choice of whether or not to consider environmental impacts.",
    )
    m.environmental_impacts.options_environmental_impacts = pyo.Param(
        m.supe_form_params.all_opts_set,
        initialize=options_environmental_impacts,
        doc="Environmental impacts matrix. Unit chosen indicator per unit of incoming flowrate.",
    )
    m.environmental_impacts.epsilon = pyo.Param(
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
    m.environmental_impacts.option_yearly_impacts = pyo.Var(
        m.supe_form_params.all_opts_set,
        m.plant_lifetime_params.operational_range,
        doc="Yearly environmental impacts generated by each option.",
    )
    m.environmental_impacts.total_yearly_impacts = pyo.Var(
        m.plant_lifetime_params.operational_range,
        doc="Total yearly environmental impacts generated by the process.",
    )
    m.environmental_impacts.total_impacts = pyo.Var(
        doc="The total environmental impacts generated by the process over its entire lifetime."
    )


def add_environmental_impact_cons(m):
    """
    This function builds the environmental impact constraints.

    Args:
        m: pyomo model.
    """

    ## Pyomo constraints
    @m.environmental_impacts.Constraint(
        m.supe_form_params.all_opts_set,
        m.plant_lifetime_params.operational_range,
        doc="Calculate the yearly environmental impacts generated by each option.",
    )
    def calculate_opt_yearly_impacts_con(b, j, k, t):
        return (
            m.environmental_impacts.option_yearly_impacts[j, k, t]
            == sum(
                m.mass_balances.f_in[j, k, c, t] for c in m.feed_params.tracked_comps
            )
            * m.environmental_impacts.options_environmental_impacts[j, k]
        )

    @m.environmental_impacts.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculate the total yearly environmental impacts generated by the process.",
    )
    def calculate_yearly_impacts_con(b, t):
        return m.environmental_impacts.total_yearly_impacts[t] == sum(
            m.environmental_impacts.option_yearly_impacts[opt, t]
            for opt in m.supe_form_params.all_opts_set
        )

    @m.environmental_impacts.Constraint(
        doc="Calculate the total environmental impacts generated by the process throughout it's lifetime."
    )
    def calculate_total_impacts_con(b):
        return m.environmental_impacts.total_impacts == sum(
            m.environmental_impacts.total_yearly_impacts[t]
            for t in m.plant_lifetime_params.operational_range
        )

    @m.environmental_impacts.Constraint(
        doc="Epsilon constraint. Total process impacts must be less than epsilon."
    )
    def epsilon_con(b):
        return m.environmental_impacts.total_impacts <= m.environmental_impacts.epsilon


###################################################################################################
### Byproduct Valorization
def check_byproduct_valorization_params(
    m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
):
    """
    This function checks that all the byproduct valorization parameters are feasible.

    Args:
        m: pyomo model.
        consider_byproduct_valorization: (bool) Decide whether or not to consider the valorization of byproducts.
        byproduct_values: (dict) Byproducts considered, and their value ($/kg).
        byproduct_opt_conversions: (dict) Defines the conversion factors for different byproducts for different options.
    """
    ### Only need to check feasibility of parameters if user wants to consider environmental impacts (consider_byprod_val is True).
    if consider_byproduct_valorization:
        ### Check typos and structure.
        ## Check that consider_byproduct_valorization is of type bool.
        if not isinstance(consider_byproduct_valorization, bool):
            raise TypeError("consider_byproduct_valorization not of type bool.")

        ## Check that structure of byproduct_values is correct.
        for key, val in byproduct_values.items():
            if not isinstance(key, str):
                raise TypeError(f"key {key} in byproduct_values is not of type str.")
            if not isinstance(val, (int, float)):
                raise TypeError(
                    f"value {val} in byproduct_values is not of type int or float."
                )

        ## Check that structure of byproduct_opt_conversions is correct
        for opt, inner_dict in byproduct_opt_conversions.items():
            # Check value is a dictionary
            if not isinstance(inner_dict, dict):
                raise TypeError(f"Value for option {opt} is not a dictionary.")
            # Check that all keys are strings, values are int/float
            for byprod, conversion in inner_dict.items():
                if not isinstance(byprod, str):
                    raise TypeError(
                        f"Key in inner dictionary of byproduct_opt_conversions[{opt}] is not of type string."
                    )
                if not (isinstance(conversion, (int, float))):
                    raise TypeError(
                        f"Value in inner dictionary of byproduct_opt_conversions[{opt}] is not of type int or float."
                    )

        ## Check that inner dictionary of byproduct_opt_conversions is not empty.
        for opt, inner_dict in byproduct_opt_conversions.items():
            # Track the options for which the inner dict is empty.
            if not inner_dict:
                # If some options have empty dicts in byproduct_opt_conversions, raise an error.
                raise TypeError(f"Empty dict passed for the option: {opt}.")

        ### Define parameters necessary for tests.
        # Create a set of the byproducts considered
        byproducts_considered = set(byproduct_values.keys())
        # Create a set of the options that produce byproducts.
        byproduct_producing_options = set(byproduct_opt_conversions.keys())
        # Create a list to track the missing byproduct values.
        incorrect_byproduct_values = []
        # Create a list to track the infeasible options.
        infeasible_options = []
        # Create a list to track byproducts that are not defined.
        undefined_byproducts = []
        # Create a list to track the byproducts being produced by options.
        produced_byproducts = []

        ### Run tests
        ## Check that a value of type int or float is defined for each byproduct.
        for byprod, val in byproduct_values.items():
            # Track the byproducts for which incorrect value is defined.
            if not isinstance(val, (int, float)):
                incorrect_byproduct_values.append(byprod)

        # Raise an error if there are byproducts for which values are incorrectly defined.
        if incorrect_byproduct_values:
            msg = "Values incorrectly specified for the following byproducts (should be of type int or float):\n"
            msg += "\n".join(
                f"  {incorrect}" for incorrect in incorrect_byproduct_values
            )
            raise ValueError(msg)

        ## Check that the set of options that produce byproducts is feasible (i.e. they exist within the superstructure).
        for opt in byproduct_producing_options:
            # Keep track of the infeasible options.
            if opt not in m.supe_form_params.all_opts_set:
                infeasible_options.append(opt)

        # Raise an error if there are infeasible options defined.
        if infeasible_options:
            msg = "Some options defined are infeasible (do not exist in the superstructure):\n"
            msg += "\n".join(f"  {infeasible}" for infeasible in infeasible_options)
            raise ValueError(msg)

        ## Check that all byproducts are defined.
        for opt, inner_dict in byproduct_opt_conversions.items():
            # Keep track of undefined byproducts.
            undefined = [
                byprod
                for byprod in inner_dict.keys()
                if byprod not in byproducts_considered
            ]
            if undefined:
                undefined_byproducts.append((opt, undefined))

        # Raise an error if there are undefined byproducts.
        if undefined_byproducts:
            msg = "The following byproducts are not defined:\n"
            msg += "\n".join(
                f"  Option: {opt} produces the undefined byproduct(s): {undefined}"
                for opt, undefined in undefined_byproducts
            )
            raise ValueError(msg)

        ## Check that all byproducts considered have at least one option that is producing it.
        for opt, inner_dict in byproduct_opt_conversions.items():
            for byprod in inner_dict.keys():
                if byprod not in produced_byproducts:
                    produced_byproducts.append(byprod)

        # turn list of produced byproducts into a set
        produced_byproducts_set = set(produced_byproducts)

        # Raise an error if there are considered byproducts that aren't being produced.
        if byproducts_considered != produced_byproducts_set:
            missing = byproducts_considered - produced_byproducts_set
            msg = "The following considered byproducts are not being produced by any options:\n"
            msg += "\n".join(f"{byprod}" for byprod in missing)
            raise ValueError(msg)


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
    ## Create blocks
    m.byproduct_valorization_params = pyo.Block()

    ## Pyomo parameters
    m.byproduct_valorization_params.consider_byproduct_valorization = pyo.Param(
        initialize=consider_byproduct_valorization,
        within=pyo.Boolean,
        doc="Choice of whether or not to consider byproduct valorization.",
    )
    m.byproduct_valorization_params.byproducts_set = pyo.Set(
        initialize=byproducts, doc="Set of byproducts considered."
    )
    m.byproduct_valorization_params.byproduct_values = pyo.Param(
        m.byproduct_valorization_params.byproducts_set,
        initialize=byproduct_values,
        doc="Defines the value of each byproduct ($/kg).",
    )
    m.byproduct_valorization_params.byproduct_opts_set = pyo.Set(
        initialize=byproduct_opt_conversions.keys(),
        doc="Set of options that produce byproducts.",
    )
    m.byproduct_valorization_params.opt_byproduct_set = pyo.Set(
        initialize=(
            (opt, byproduct)
            for opt in byproduct_opt_conversions.keys()
            for byproduct in byproduct_opt_conversions[opt].keys()
        ),
        doc="Set of options that produce byproducts, and the byproducts that they produce.",
    )

    m.byproduct_valorization_params.byproduct_opt_conversion = pyo.Param(
        m.byproduct_valorization_params.opt_byproduct_set,
        initialize=byproduct_opt_conversion_initialize,
        doc="Defines the conversion factors for all byproducts for all options that produce them.",
    )
    m.byproduct_valorization_params.byproduct_producing_opts = pyo.Param(
        m.byproduct_valorization_params.byproducts_set,
        initialize=byproduct_producing_opts,
        doc="Pyomo parameter which defines what options produce each byproduct.",
    )


def add_byproduct_valorization_vars(m):
    """
    This function builds the byproduct valorization variables.

    Args:
        m: pyomo model.
    """
    ## Create blocks
    m.byproduct_valorization = pyo.Block()

    ## Pyomo variables
    m.byproduct_valorization.byproduct_produced = pyo.Var(
        m.byproduct_valorization_params.byproducts_set
        * m.plant_lifetime_params.operational_range,
        domain=pyo.NonNegativeReals,
        doc="The amount of each byproduct produced each year.",
    )
    m.byproduct_valorization.byproduct_profit = pyo.Var(
        m.byproduct_valorization_params.byproducts_set
        * m.plant_lifetime_params.operational_range,
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
    m.no_byproduct_valorization = pyo.Block()

    ## Pyomo constraints
    @m.no_byproduct_valorization.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Sets the profit from the valorization of byproducts to zero.",
    )
    def calculate_opt_byprod_val_cons(b, t):
        return m.costing.total_byproduct_profit[t] == 0

    @m.byproduct_valorization.Constraint(
        m.byproduct_valorization_params.byproducts_set,
        m.plant_lifetime_params.operational_range,
        doc="Calculates the amount of each byproduct produced each year.",
    )
    def calculate_byproduct_produced_cons(b, byprod, t):
        return m.byproduct_valorization.byproduct_produced[byprod, t] == sum(
            sum(m.mass_balances.f_in[opt, c, t] for c in m.feed_params.tracked_comps)
            * m.byproduct_valorization_params.byproduct_opt_conversion[opt, byprod]
            for opt in m.byproduct_valorization_params.byproduct_producing_opts[byprod]
        )

    @m.byproduct_valorization.Constraint(
        m.byproduct_valorization_params.byproducts_set,
        m.plant_lifetime_params.operational_range,
        doc="Calculate the yearly profit from each byproduct.",
    )
    def calculate_byproduct_profit_cons(b, byprod, t):
        return (
            m.byproduct_valorization.byproduct_profit[byprod, t]
            == m.byproduct_valorization.byproduct_produced[byprod, t]
            * m.byproduct_valorization_params.byproduct_values[byprod]
        )

    @m.byproduct_valorization.Constraint(
        m.plant_lifetime_params.operational_range,
        doc="Calculates the total yearly profit from the valorization of byproducts.",
    )
    def calculate_opt_byprod_val_cons(b, t):
        return m.costing.total_byproduct_profit[t] == sum(
            m.byproduct_valorization.byproduct_profit[byprod, t]
            for byprod in m.byproduct_valorization_params.byproducts_set
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
        m.cost_of_recovery.deactivate()
    else:
        # deactivate the net present value constraints and objective function if the cost of recovery objective function is chosen.
        m.net_present_value.deactivate()

    # Check if environmental impacts are considered.
    if pyo.value(m.environmental_impacts.consider_environmental_impacts) == False:
        # deactivate environmental impact constraints if they are not considered.
        m.environmental_impacts.deactivate()

    # Check if byproduct valorization is considered.  Different constraints are considered depending on if they are considered or not.
    if (
        pyo.value(m.byproduct_valorization_params.consider_byproduct_valorization)
        == True
    ):
        # deactivate the constraints that're associated with no byproduct valorization if it is considered.
        m.no_byproduct_valorization.deactivate()
    else:
        # deactivate the constraints that're associated with byproduct valorization if it is not considered.
        m.byproduct_valorization.deactivate()
