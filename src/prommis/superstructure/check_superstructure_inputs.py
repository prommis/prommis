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

import warnings

import pyomo.environ as pyo


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
    operational_years = set(m.fs.operational_range.data())
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

    ## Check that option_outlets is not empty
    if not option_outlets:
        raise TypeError("option_outlets is an empty dict.")

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

    ## Check that option_efficiencies is not an empty dict.
    if not option_efficiencies:
        raise TypeError("option_efficiencies is an empty dict.")

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
    tracked_comps = set(m.fs.tracked_comps.data())
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


### Operating Parameters
profit = {
    (5, 1): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
    (5, 2): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
    (5, 3): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
    (5, 4): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
    (5, 5): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
}
opt_var_oc_params = {
    # level 2
    (2, 1): {"a": 0.0053, "b": 7929.7},
    (2, 2): {"a": 0.0015, "b": 2233.16},
    (2, 3): {"a": 0.0034, "b": 0},
    (2, 4): {"a": 0.0117, "b": 0},
    # level 3
    (3, 1): {"a": 15.594, "b": 4e6},
    (3, 2): {"a": 35.58463, "b": 4e6},
    (3, 3): {"a": 1.8359, "b": 0},
    (3, 4): {"a": 3.7414, "b": 2378.6},
    (3, 5): {"a": 10.35427, "b": 2378.6},
    (3, 6): {"a": 1.58, "b": 0},
    # level 4
    (4, 1): {"a": 0, "b": 0},
    (4, 2): {"a": 111.09, "b": 254606},
    (4, 3): {"a": 0, "b": 0},
    (4, 4): {"a": 0, "b": 0},
    # level 5
    (5, 1): {"a": 0.4997, "b": 89832},
    (5, 2): {"a": 9.8127, "b": 964921},
    (5, 3): {"a": 9.8127, "b": 964921},
    (5, 4): {"a": 2.17, "b": 0},
    (5, 5): {"a": 6.7063559004, "b": 0},
}
operators_per_discrete_unit = {
    (1, 1): 1,
    (1, 2): 0,
}
yearly_cost_per_unit = {
    (1, 1): 0,
    (1, 2): 280,
}
capital_cost_per_unit = {
    (1, 1): 0,
    (1, 2): 200000,
}
processing_rate = {
    (1, 1): 7868,
    (1, 2): 52453,
}
num_operators = {
    (2, 1): 0.65,
    (2, 2): 0.65,
    (2, 3): 0.65,
    (2, 4): 0.65,
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
    (3, 4): 0.45,
    (3, 5): 0.45,
    (3, 6): 1.15,
    (4, 1): 0,
    (4, 2): 1.3,
    (4, 3): 0,
    (4, 4): 0,
    (5, 1): 1.05,
    (5, 2): 0.75,
    (5, 3): 0.75,
    (5, 4): 1.15,
    (5, 5): 1.15,
}
labor_rate = 8000 * 38.20


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
    final_opts_set = set(m.fs.final_opts_set)
    # Define a set of all the tracked components.
    tracked_comps = set(m.fs.tracked_comps.data())
    # Define a set containing the keys in the opt_var_oc_params dict.
    opt_var_oc_params_keys = set(opt_var_oc_params.keys())
    # Define a set containing all the discrete options.
    discrete_opts_set = set(m.fs.discrete_opts_set.data())
    # Define a set containing all the continuous options.
    continuous_opts_set = set(m.fs.continuous_opts_set.data())
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
    for opt in m.fs.final_opts_set:
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
                    f"Empty list passed for option: {opt}, data type: {data_type}."
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
    for opt in m.fs.continuous_opts_set:
        # Keep track of the continuous opts for which discretized data is not defined for.
        if opt not in discretized_purchased_equipment_cost_opts_set:
            missing_continuous_opts.append(opt)
    # Check that discretized capex not provided for discrete opts.
    for opt in m.fs.discrete_opts_set:
        # Keep track of the discrete options for which discretized capex data is defined for.
        if opt in discretized_purchased_equipment_cost_opts_set:
            discrete_opts.append(opt)
    # Check all options in discretized capex are feasible.
    for opt in discretized_purchased_equipment_cost_opts_set:
        # Keep track of the infeasible opts
        if opt not in m.fs.continuous_opts_set:
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
    for opt in m.fs.continuous_opts_set:
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
    ## Check that consider_environmental_impacts is of type bool. Must be provided regardless on whether or not environmental impacts are considered.
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
        all_opts_set = set(m.fs.all_opts_set.data())
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
            if opt not in m.fs.all_opts_set:
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

def check_objective_function_choice(obj_func):
    """
    This function checks that the choice of objective function is feasible.

    Args:
        obj_func: choice of objective function. Options are 'NPV' or 'COR'. Case sensitive.
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