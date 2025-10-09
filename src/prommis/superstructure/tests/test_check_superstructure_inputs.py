import pyomo.environ as pyo

import pytest

from prommis.superstructure.add_superstructure_blocks import (
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
    add_operating_params,
    add_plant_lifetime_params,
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
from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice
from prommis.superstructure.superstructure_function import define_custom_units

### Define correct parameters
# define custom units
define_custom_units()
obj_func = ObjectiveFunctionChoice.NET_PRESENT_VALUE

plant_start = 2024
plant_lifetime = 5

available_feed = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
}
collection_rate = 0.1
tracked_comps = ["Nd", "Fe"]
prod_comp_mass = {
    "Nd": 0.206 * 3,
    "Fe": 0.691 * 3,
}

num_stages = 3
options_in_stage = {1: 1, 2: 2, 3: 3}
option_outlets = {
    (1, 1): [1, 2],
    (2, 1): [1],
    (2, 2): [2, 3],
}

option_efficiencies = {
    (1, 1): {"Nd": 1, "Fe": 1},
    (2, 1): {"Nd": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Fe": 1},
    (3, 1): {"Nd": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Fe": 0},
    (3, 3): {"Nd": 0.925, "Fe": 0},
}

profit = {
    (3, 1): {"Nd": 45.4272, "Fe": 0},
    (3, 2): {"Nd": 69.888, "Fe": 0},
    (3, 3): {"Nd": 45.4272, "Fe": 0},
}
opt_var_oc_params = {
    (2, 1): {"a": 0.0053, "b": 7929.7},
    (2, 2): {"a": 0.0015, "b": 2233.16},
    (3, 1): {"a": 15.594, "b": 4e6},
    (3, 2): {"a": 35.58463, "b": 4e6},
    (3, 3): {"a": 1.8359, "b": 0},
}
operators_per_discrete_unit = {(1, 1): 1}
yearly_cost_per_unit = {(1, 1): 0}
capital_cost_per_unit = {(1, 1): 0}
processing_rate = {(1, 1): 7868}
num_operators = {
    (2, 1): 0.65,
    (2, 2): 0.65,
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
}
labor_rate = 8000 * 38.20

discretized_equipment_cost = {
    (2, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            10130.08515,
            31353.21173,
            48788.84678,
            60305.81927,
            77063.4884,
            117214.7546,
            151018.0699,
            195698.5419,
        ],
    },
    (2, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            11702.08515,
            39023.21173,
            62134.84678,
            77539.81927,
            100113.4884,
            154792.7546,
            201326.0699,
            263374.5419,
        ],
    },
    (3, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            343228.652,
            482425.4684,
            618182.0594,
            743750.2902,
            844443.0443,
            978479.5225,
            1183834.522,
            1440660.587,
        ],
    },
    (3, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            643228.652,
            782425.4684,
            918182.0594,
            1043750.2902,
            1144443.0443,
            1278479.5225,
            1483834.522,
            1740660.587,
        ],
    },
    (3, 3): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            423074.7216,
            3042779.121,
            5348359.01,
            6921261.68,
            9251002.61,
            14933803.33,
            19762044.75,
            26151302.79,
        ],
    },
}

consider_environmental_impacts = True
options_environmental_impacts = {
    (1, 1): 0,
    (2, 1): 0,
    (2, 2): 1000,
    (3, 1): 600,
    (3, 2): 0,
    (3, 3): 600,
}
epsilon = 1

consider_byproduct_valorization = True
byproduct_values = {
    "Jarosite": -0.17,
    "Iron oxide": 10,
    "Residue": -0.17,
}
byproduct_opt_conversions = {
    (3, 1): {"Jarosite": 0.75},
    (3, 2): {"Iron oxide": 1},
    (3, 3): {"Residue": 0.25},
}


## test check_objective_function_choice
def test_check_objective_function_choice():
    ## test correct input
    check_objective_function_choice(obj_func)

    ## test incorrect inputs - objective function not specified as ObjectiveFunctionChoice enum
    incorrect_obj_choice = "net_present_value"
    with pytest.raises(TypeError):
        check_objective_function_choice(incorrect_obj_choice)


## test check_plant_lifetime_params
def test_check_plant_lifetime_params():
    ## Define necessary parameters
    # incorrect plant_start
    incorrect_plant_start = 2020.0

    ## test correct input
    check_plant_lifetime_params(plant_start, plant_lifetime)

    ## test incorrect inputs 1 - plant_lifetime not int
    incorrect_plant_lifetime1 = 15.0
    with pytest.raises(TypeError):
        check_plant_lifetime_params(plant_start, incorrect_plant_lifetime1)

    ## test incorrect inputs 2 - plant_lifetime not at least 3
    incorrect_plant_lifetime2 = 2
    with pytest.raises(ValueError):
        check_plant_lifetime_params(plant_start, incorrect_plant_lifetime2)

    ## test incorrect inputs 3 - plant_start is not of type int
    incorrect_plant_start3 = 2024.0
    with pytest.raises(TypeError):
        check_plant_lifetime_params(incorrect_plant_start3, plant_lifetime)


def test_check_feed_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)

    ## test correct inputs
    check_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)

    ## test incorrect inputs 1 - empty available feed dict
    incorrect_available_feed1 = {}
    with pytest.raises(TypeError):
        check_feed_params(
            m, incorrect_available_feed1, collection_rate, tracked_comps, prod_comp_mass
        )

    ## test incorrect inputs 2 - key in available_feed dict is not of type int
    incorrect_available_feed2 = {
        2025.0: 290273,
        2026: 274648,
        2027: 286512,
        2028: 487819,
    }
    with pytest.raises(TypeError):
        check_feed_params(
            m, incorrect_available_feed2, collection_rate, tracked_comps, prod_comp_mass
        )

    ## test incorrect inputs 3 - value in available_feed dict is not of type int or float
    incorrect_available_feed3 = {
        2025: "290273",
        2026: 274648,
        2027: 286512,
        2028: 487819,
    }
    with pytest.raises(TypeError):
        check_feed_params(
            m, incorrect_available_feed3, collection_rate, tracked_comps, prod_comp_mass
        )

    ## test incorrect inputs 4 - collection rate not of type int or float
    incorrect_collection_rate4 = "0.1"
    with pytest.raises(TypeError):
        check_feed_params(
            m, available_feed, incorrect_collection_rate4, tracked_comps, prod_comp_mass
        )

    ## test incorrect inputs 5 - tracked_comps is not a list
    incorrect_tracked_comps5 = 1
    with pytest.raises(TypeError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            incorrect_tracked_comps5,
            prod_comp_mass,
        )

    ## test incorrect inputs 6 - not all components in tracked_comps are of type str
    incorrect_tracked_comps6 = ["Nd", 8]
    with pytest.raises(TypeError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            incorrect_tracked_comps6,
            prod_comp_mass,
        )

    ## test incorrect inputs 7 - prod_comp_mass not a dict
    incorrect_prod_comp_mass7 = 7
    with pytest.raises(TypeError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass7,
        )

    ## test incorrect inputs 8 - key of prod_comp_mass is not a str
    incorrect_prod_comp_mass8 = {
        8: 0.206 * 3,
        "Fe": 0.691 * 3,
    }
    with pytest.raises(TypeError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass8,
        )

    ## test incorrect inputs 9 - value in prod_comp_mass dict is not int or float
    incorrect_prod_comp_mass8 = {
        "Nd": "0.206 * 3",
        "Fe": 0.691 * 3,
    }
    with pytest.raises(TypeError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass8,
        )

    ## test incorrect inputs 10 - available feed not provided for each year of plant operation
    incorrect_available_feed10 = {
        2025: 290273,
        2026: 274648,
        2027: 286512,
    }
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            incorrect_available_feed10,
            collection_rate,
            tracked_comps,
            prod_comp_mass,
        )

    ## test incorrect inputs 11 - available feed negative for some years
    incorrect_available_feed11 = {
        2025: 290273,
        2026: 274648,
        2027: 286512,
        2028: -487819,
    }
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            incorrect_available_feed11,
            collection_rate,
            tracked_comps,
            prod_comp_mass,
        )

    ## test incorrect inputs 12 - available feed is all zero
    incorrect_available_feed12 = {
        2025: 0,
        2026: 0,
        2027: 0,
        2028: 0,
    }
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            incorrect_available_feed12,
            collection_rate,
            tracked_comps,
            prod_comp_mass,
        )

    ## test incorrect inputs 13 - collection rate is negative
    incorrect_collection_rate13 = -0.1
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            available_feed,
            incorrect_collection_rate13,
            tracked_comps,
            prod_comp_mass,
        )

    ## test incorrect inputs 14 - no tracked comps specified
    incorrect_tracked_comps14 = []
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            incorrect_tracked_comps14,
            prod_comp_mass,
        )

    ## test incorrect inputs 15 - amount contained within EOL product is not specified for each tracked component
    incorrect_prod_comp_mass15 = {
        "Nd": 0.206 * 3,
    }
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass15,
        )

    ## test incorrect inputs 16 - amount contained within EOL product is negative for a tracked component
    incorrect_prod_comp_mass16 = {
        "Nd": -0.206 * 3,
        "Fe": 0.691 * 3,
    }
    with pytest.raises(ValueError):
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass16,
        )

    ## test incorrect inputs 17 - check that warning is raised if amounts contained within EOL product for tracked component is 0
    incorrect_prod_comp_mass17 = {
        "Nd": 0,
        "Fe": 0.691 * 3,
    }
    with pytest.warns():
        check_feed_params(
            m,
            available_feed,
            collection_rate,
            tracked_comps,
            incorrect_prod_comp_mass17,
        )


def test_check_supe_formulation_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)
    # add feed parameters
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)

    ## test correct inputs
    check_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )

    ## test incorrect inputs 1 - num_stages not int
    incorrect_num_stages1 = 3.0
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            incorrect_num_stages1,
            options_in_stage,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 2 - options_in_stage is not a dict
    incorrect_options_in_stage2 = "test"
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage2,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 3 - options_in_stage is empty
    incorrect_options_in_stage3 = {}
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage3,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 4 - key in options_in_stage is not of type int
    incorrect_options_in_stage4 = {1.0: 1, 2: 2, 3: 3}
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage4,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 5 - value in options_in_stage is not of type int
    incorrect_options_in_stage5 = {1: 1.0, 2: 2, 3: 3}
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage5,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 6 - value in options_in_stage is negative
    incorrect_options_in_stage6 = {1: -1, 2: 2, 3: 3}
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage6,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 7 - option_outlets is not a dict
    incorrect_option_outlets7 = 7
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets7,
            option_efficiencies,
        )

    ## test incorrect inputs 8 - option_outlets is an empty dict
    incorrect_option_outlets8 = {}
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets8,
            option_efficiencies,
        )

    ## test incorrect inputs 9 - option_efficiencies is not a dict
    incorrect_option_efficiencies9 = 9
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies9,
        )

    ## test incorrect inputs 10 - option_efficiencies is an empty dict
    incorrect_option_efficiencies10 = {}
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies10,
        )

    ## test incorrect inputs 11 - key in option_efficiencies is not a tuple
    incorrect_option_efficiencies11 = {
        1: {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies11,
        )

    ## test incorrect inputs 12 - incorrect structure of option_efficiencies. No inner dict
    incorrect_option_efficiencies12 = {
        (1, 1): "Nd",
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies12,
        )

    ## test incorrect inputs 13 - incorrect structure of option_efficiencies. key for inner dict is not of type str
    incorrect_option_efficiencies13 = {
        (1, 1): {8: 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies13,
        )

    ## test incorrect inputs 14 - incorrect structure of option_efficiencies. value for inner dict is not of type str
    incorrect_option_efficiencies14 = {
        (1, 1): {"Nd": "1", "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies14,
        )

    ## test incorrect inputs 15 - Efficiency value is infeasible
    incorrect_option_efficiencies15 = {
        (1, 1): {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 2},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies15,
        )

    ## test incorrect inputs 16 - options_in_stage not provided for all stages
    incorrect_options_in_stage16 = {1: 1, 2: 2}
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage16,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 17 - options_in_stage defined for stages that don't exist
    incorrect_options_in_stage17 = {1: 1, 2: 2, 3: 3, 4: 4}
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            incorrect_options_in_stage17,
            option_outlets,
            option_efficiencies,
        )

    ## test incorrect inputs 18 - there are <2 stages
    incorrect_num_stages18 = 1
    incorrect_options_in_stage18 = {1: 1}
    incorrect_option_outlets18 = {
        (1, 1): [1, 2],
    }
    incorrect_option_efficiencies18 = {
        (1, 1): {"Nd": 1, "Fe": 1},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            incorrect_num_stages18,
            incorrect_options_in_stage18,
            incorrect_option_outlets18,
            incorrect_option_efficiencies18,
        )

    ## test incorrect inputs 19 - there are disconnected options
    incorrect_num_stages19 = 3
    incorrect_options_in_stage19 = {1: 1, 2: 2, 3: 3}
    incorrect_option_outlets19 = {
        (1, 1): [1, 2],
        (2, 1): [2],
        (2, 2): [2, 3],
    }
    incorrect_option_efficiencies19 = {
        (1, 1): {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            incorrect_num_stages19,
            incorrect_options_in_stage19,
            incorrect_option_outlets19,
            incorrect_option_efficiencies19,
        )

    ## test incorrect inputs 20 - infeasible outlet in option_outlets
    incorrect_num_stages20 = 3
    incorrect_options_in_stage20 = {1: 1, 2: 2, 3: 3}
    incorrect_option_outlets20 = {
        (1, 1): [1, 2],
        (2, 1): [1, 4],
        (2, 2): [2, 3],
    }
    incorrect_option_efficiencies20 = {
        (1, 1): {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            incorrect_num_stages20,
            incorrect_options_in_stage20,
            incorrect_option_outlets20,
            incorrect_option_efficiencies20,
        )

    ## test incorrect inputs 21 - option efficiencies not specified for all options in superstructure
    incorrect_option_efficiencies21 = {
        (1, 1): {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies21,
        )

    ## test incorrect inputs 22 - option efficiencies defined for options that are not in the superstructure
    incorrect_option_efficiencies22 = {
        (1, 1): {"Nd": 1, "Fe": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (2, 3): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies22,
        )

    ## test incorrect inputs 23 - option efficiencies not defined for all tracked components in all options
    incorrect_option_efficiencies23 = {
        (1, 1): {"Nd": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies23,
        )

    ## test incorrect inputs 24 - some components in option_efficiencies not defined as a tracked component
    incorrect_option_efficiencies24 = {
        (1, 1): {"Nd": 1, "Fe": 1, "Ce": 1},
        (2, 1): {"Nd": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Fe": 1},
        (3, 1): {"Nd": 0.985, "Fe": 0},
        (3, 2): {"Nd": 0.985, "Fe": 0},
        (3, 3): {"Nd": 0.925, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            option_outlets,
            incorrect_option_efficiencies24,
        )

    ## test incorrect inputs 25 - key in option_outlets is not of type tuple
    incorrect_option_outlets25 = {
        1: [1, 2],
        (2, 1): [1],
        (2, 2): [2, 3],
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets25,
            option_efficiencies,
        )

    ## test incorrect inputs 26 - value in option_outlets is not of type list
    incorrect_option_outlets26 = {
        (1, 1): {1, 2},
        (2, 1): [1],
        (2, 2): [2, 3],
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets26,
            option_efficiencies,
        )

    ## test incorrect inputs 27 - value in option_outlets is an empty list
    incorrect_option_outlets27 = {
        (1, 1): [],
        (2, 1): [1],
        (2, 2): [2, 3],
    }
    with pytest.raises(ValueError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets27,
            option_efficiencies,
        )

    ## test incorrect inputs 28 - option_outlets list doesn't contain all positive integers
    incorrect_option_outlets28 = {
        (1, 1): [1, -2],
        (2, 1): [1],
        (2, 2): [2, 3],
    }
    with pytest.raises(TypeError):
        check_supe_formulation_params(
            m,
            num_stages,
            options_in_stage,
            incorrect_option_outlets28,
            option_efficiencies,
        )


def test_check_operating_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)
    # add feed parameters
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # add superstructure formulations parameters
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )

    ## test correct inputs
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

    ## test incorrect inputs 1 - profit not of type dict
    incorrect_profit1 = []
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit1,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 2 - profit is empty
    incorrect_profit2 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit2,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 3 - key in profit is not of type tuple
    incorrect_profit3 = {
        3: {"Nd": 45.4272, "Fe": 0},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit3,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 4 - incorrect structure of profit. No inner dict
    incorrect_profit4 = {
        (3, 1): "Nd",
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit4,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 5 - key in inner dict of profit is not of type str
    incorrect_profit5 = {
        (3, 1): {6: 45.4272, "Fe": 0},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit5,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 6 - value of in inner dict of profit is not of type int or float
    incorrect_profit6 = {
        (3, 1): {"Nd": "45.4272", "Fe": 0},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            incorrect_profit6,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 7 - opt_var_oc_params is empty
    incorrect_opt_var_oc_params7 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params7,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 8 - key in opt_var_oc_params is not of type tuple
    incorrect_opt_var_oc_params8 = {
        2: {"a": 0.0053, "b": 7929.7},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params8,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 9 - opt_var_oc_params has incorrect structure. no inner dict.
    incorrect_opt_var_oc_params9 = {
        (2, 1): 1,
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params9,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 10 - opt_var_oc_params inner dict key is not of type str
    incorrect_opt_var_oc_params10 = {
        (2, 1): {1: 0.0053, "b": 7929.7},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params10,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 11 - opt_var_oc_params inner dict value is not of type int or float
    incorrect_opt_var_oc_params11 = {
        (2, 1): {"a": "0.0053", "b": 7929.7},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params10,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 12 - operators_per_discrete_unit is empty
    incorrect_operators_per_discrete_unit12 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            incorrect_operators_per_discrete_unit12,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 13 - key in operators_per_discrete_unit is not of type tuple
    incorrect_operators_per_discrete_unit13 = {1: 1}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            incorrect_operators_per_discrete_unit13,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 14 - value in operators_per_discrete_unit is not of type int or float
    incorrect_operators_per_discrete_unit14 = {(1, 1): "1"}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            incorrect_operators_per_discrete_unit14,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 15 - yearly_cost_per_unit is empty
    incorrect_yearly_cost_per_unit15 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            incorrect_yearly_cost_per_unit15,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 16 - key in yearly_cost_per_unit is not of type tuple
    incorrect_yearly_cost_per_unit16 = {1: 0}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            incorrect_yearly_cost_per_unit16,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 17 - value of yearly_cost_per_unit is not of type int or float
    incorrect_yearly_cost_per_unit17 = {(1, 1): "0"}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            incorrect_yearly_cost_per_unit17,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 18 - capital_cost_per_unit is empty
    incorrect_capital_cost_per_unit18 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            incorrect_capital_cost_per_unit18,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 19 - capital_cost_per_unit key is not of type tuple
    incorrect_capital_cost_per_unit19 = {1: 0}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            incorrect_capital_cost_per_unit19,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 20 - capital_cost_per_unit value is not of type int or float
    incorrect_capital_cost_per_unit20 = {(1, 1): "0"}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            incorrect_capital_cost_per_unit20,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 21 - processing_rate is empty
    incorrect_processing_rate21 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            incorrect_processing_rate21,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 22 - processing_rate key is not of type tuple
    incorrect_processing_rate22 = {1: 7868}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            incorrect_processing_rate22,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 23 - processing_rate value is not of type int or float
    incorrect_processing_rate23 = {(1, 1): "7868"}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            incorrect_processing_rate23,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 24 - num_operators is empty
    incorrect_num_operators24 = {}
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            incorrect_num_operators24,
            labor_rate,
        )

    ## test incorrect inputs 25 - key in num_operators is not of type tuple
    incorrect_num_operators25 = {
        2: 0.65,
        (2, 2): 0.65,
        (3, 1): 1.6,
        (3, 2): 1.6,
        (3, 3): 1.3,
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            incorrect_num_operators25,
            labor_rate,
        )

    ## test incorrect inputs 26 - value for num_operators is not of type int or float
    incorrect_num_operators26 = {
        (2, 1): "0.65",
        (2, 2): 0.65,
        (3, 1): 1.6,
        (3, 2): 1.6,
        (3, 3): 1.3,
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            incorrect_num_operators26,
            labor_rate,
        )

    ## test incorrect inputs 27 - labor_rate is not of type int or float
    incorrect_labor_rate27 = "a"
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            incorrect_labor_rate27,
        )

    ## test incorrect inputs 28 - profit per product not defined for all options in the final stage
    incorrect_profit28 = {
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            incorrect_profit28,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 29 - profit per product defined for options not in final stage
    incorrect_profit29 = {
        (2, 1): {"Nd": 45.4272, "Fe": 0},
        (3, 1): {"Nd": 45.4272, "Fe": 0},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            incorrect_profit29,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 30 - profit per product defined for all key components for all options in final stage
    incorrect_profit30 = {
        (3, 1): {"Nd": 45.4272},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            incorrect_profit30,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 31 - profit per product defined as negative for some key components
    incorrect_profit31 = {
        (3, 1): {"Nd": -45.4272, "Fe": 0},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            incorrect_profit31,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 32 - not all components in profit are tracked components
    incorrect_profit32 = {
        (3, 1): {"Nd": 45.4272, "Fe": 0, "a": 1},
        (3, 2): {"Nd": 69.888, "Fe": 0},
        (3, 3): {"Nd": 45.4272, "Fe": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            incorrect_profit32,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 33 - variable operating cost params are not defined for all continuous options
    incorrect_opt_var_oc_params33 = {
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params33,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 34 - variable operating costs defined for infeasible options
    incorrect_opt_var_oc_params34 = {
        (1, 1): {"a": 0.0053, "b": 7929.7},
        (2, 1): {"a": 0.0053, "b": 7929.7},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params34,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 35 - both necessary variable operating cost parameters ('a' and 'b') not defined for all continuous options.
    incorrect_opt_var_oc_params35 = {
        (2, 1): {"a": 0.0053},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params35,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 36 - operators per discrete unit not defined for all options that utilize discrete units
    incorrect_operators_per_discrete_unit36 = {(1, 3): 1}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            incorrect_operators_per_discrete_unit36,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 37 - operators per discrete unit is negative for some options
    incorrect_operators_per_discrete_unit37 = {(1, 1): -1}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            incorrect_operators_per_discrete_unit37,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 38 - yearly_cost_per_unit not defined for all discrete options
    incorrect_yearly_cost_per_unit38 = {(1, 3): 0}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            incorrect_yearly_cost_per_unit38,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 39 - yearly_cost_per_unit values negative for some options
    incorrect_yearly_cost_per_unit39 = {(1, 1): -1}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            incorrect_yearly_cost_per_unit39,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 40 - capital_cost_per_unit not defined for all discrete options
    incorrect_capital_cost_per_unit40 = {(1, 3): 0}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            incorrect_capital_cost_per_unit40,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 41 - negative capital cost defined for some discrete options
    incorrect_capital_cost_per_unit41 = {(1, 1): -1}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            incorrect_capital_cost_per_unit41,
            processing_rate,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 42 - processing rate not defined for all discrete options
    incorrect_processing_rate42 = {(1, 3): 7868}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            incorrect_processing_rate42,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 43 - negative processing rate defined for some discrete options
    incorrect_processing_rate43 = {(1, 1): -7868}
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            incorrect_processing_rate43,
            num_operators,
            labor_rate,
        )

    ## test incorrect inputs 44 - num_operators not defined for all continuous options
    incorrect_num_operators44 = {
        (2, 2): 0.65,
        (3, 1): 1.6,
        (3, 2): 1.6,
        (3, 3): 1.3,
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            incorrect_num_operators44,
            labor_rate,
        )

    ## test incorrect inputs 45 - num_operators defined as negative for some continuous options
    incorrect_num_operators45 = {
        (2, 1): -0.65,
        (2, 2): 0.65,
        (3, 1): 1.6,
        (3, 2): 1.6,
        (3, 3): 1.3,
    }
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            incorrect_num_operators45,
            labor_rate,
        )

    ## test incorrect inputs 46 - labor rate is defined as negative
    incorrect_labor_rate46 = -8000 * 38.20
    with pytest.raises(ValueError):
        check_operating_params(
            m,
            profit,
            opt_var_oc_params,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            incorrect_labor_rate46,
        )

    ## test incorrect inputs 47 - value for inner dict of opt_var_oc_params not of type int or float.
    incorrect_opt_var_oc_params47 = {
        (2, 1): {"a": "0.0053", "b": 7929.7},
        (2, 2): {"a": 0.0015, "b": 2233.16},
        (3, 1): {"a": 15.594, "b": 4e6},
        (3, 2): {"a": 35.58463, "b": 4e6},
        (3, 3): {"a": 1.8359, "b": 0},
    }
    with pytest.raises(TypeError):
        check_operating_params(
            m,
            profit,
            incorrect_opt_var_oc_params47,
            operators_per_discrete_unit,
            yearly_cost_per_unit,
            capital_cost_per_unit,
            processing_rate,
            num_operators,
            labor_rate,
        )


def test_check_discretized_costing_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)
    # add feed parameters
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # add superstructure formulations parameters
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )
    # Add operating parameters
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

    ## test correct inputs
    check_discretized_costing_params(m, discretized_equipment_cost)

    ## test incorrect inputs 1 - discretized_equipment_cost is not of type dict
    incorrect_discretized_equipment_cost1 = 5
    with pytest.raises(TypeError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost1)

    ## test incorrect inputs 2 - discretized_equipment_cost is an empty dict
    incorrect_discretized_equipment_cost2 = {}
    with pytest.raises(TypeError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost2)

    ## test incorrect inputs 3 - key in discretized_equipment_cost is not of type tuple
    incorrect_discretized_equipment_cost3 = {
        2: {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        }
    }
    with pytest.raises(TypeError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost3)

    ## test incorrect inputs 4 - discretized_equipment_cost has incorrect structure. no inner dict
    incorrect_discretized_equipment_cost4 = {
        (2, 1): [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ]
    }
    with pytest.raises(TypeError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost4)

    ## test incorrect inputs 5 - discretized_equipment_cost has incorrect structure. inner dict is empty
    incorrect_discretized_equipment_cost5 = {(2, 1): {}}
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost5)

    ## test incorrect inputs 6 - discretized_equipment_cost has incorrect structure. value for inner dict not of type list
    incorrect_discretized_equipment_cost6 = {
        (2, 1): {
            "Flowrates": {
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            },
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        },
    }
    with pytest.raises(TypeError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost6)

    ## test incorrect inputs 7 - discretized_equipment_cost has incorrect structure. value for inner dict is empty list
    incorrect_discretized_equipment_cost7 = {
        (2, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [],
        },
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost7)

    ## test incorrect inputs 8 - discretized_equipment_cost missing values for some continuous options
    incorrect_discretized_equipment_cost8 = {
        (2, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        }
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost8)

    ## test incorrect inputs 9 - discretized_equipment_cost contains values for some discrete options
    incorrect_discretized_equipment_cost9 = {
        (1, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        },
        (2, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        },
        (2, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                11702.08515,
                39023.21173,
                62134.84678,
                77539.81927,
                100113.4884,
                154792.7546,
                201326.0699,
                263374.5419,
            ],
        },
        (3, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                343228.652,
                482425.4684,
                618182.0594,
                743750.2902,
                844443.0443,
                978479.5225,
                1183834.522,
                1440660.587,
            ],
        },
        (3, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                643228.652,
                782425.4684,
                918182.0594,
                1043750.2902,
                1144443.0443,
                1278479.5225,
                1483834.522,
                1740660.587,
            ],
        },
        (3, 3): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                423074.7216,
                3042779.121,
                5348359.01,
                6921261.68,
                9251002.61,
                14933803.33,
                19762044.75,
                26151302.79,
            ],
        },
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost9)

    ## test incorrect inputs 10 - discretized_equipment_cost contains values for non-existent continuous options
    incorrect_discretized_equipment_cost10 = {
        (2, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        },
        (2, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                11702.08515,
                39023.21173,
                62134.84678,
                77539.81927,
                100113.4884,
                154792.7546,
                201326.0699,
                263374.5419,
            ],
        },
        (3, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                343228.652,
                482425.4684,
                618182.0594,
                743750.2902,
                844443.0443,
                978479.5225,
                1183834.522,
                1440660.587,
            ],
        },
        (3, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                643228.652,
                782425.4684,
                918182.0594,
                1043750.2902,
                1144443.0443,
                1278479.5225,
                1483834.522,
                1740660.587,
            ],
        },
        (3, 3): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                423074.7216,
                3042779.121,
                5348359.01,
                6921261.68,
                9251002.61,
                14933803.33,
                19762044.75,
                26151302.79,
            ],
        },
        (3, 4): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                423074.7216,
                3042779.121,
                5348359.01,
                6921261.68,
                9251002.61,
                14933803.33,
                19762044.75,
                26151302.79,
            ],
        },
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost10)

    ## test incorrect inputs 11 - discretized_equipment_cost contains different number of discretized data points between flowrates and costs for some options
    incorrect_discretized_equipment_cost11 = {
        (2, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
        },
        (2, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                11702.08515,
                39023.21173,
                62134.84678,
                77539.81927,
                100113.4884,
                154792.7546,
                201326.0699,
                263374.5419,
            ],
        },
        (3, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                343228.652,
                482425.4684,
                618182.0594,
                743750.2902,
                844443.0443,
                978479.5225,
                1183834.522,
                1440660.587,
            ],
        },
        (3, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                643228.652,
                782425.4684,
                918182.0594,
                1043750.2902,
                1144443.0443,
                1278479.5225,
                1483834.522,
                1740660.587,
            ],
        },
        (3, 3): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                423074.7216,
                3042779.121,
                5348359.01,
                6921261.68,
                9251002.61,
                14933803.33,
                19762044.75,
            ],
        },
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost11)

    ## test incorrect inputs 12 - discretized_equipment_cost contains values for some discrete options
    incorrect_discretized_equipment_cost12 = {
        (2, 1): {
            "Costs": [
                0.0,
                10130.08515,
                31353.21173,
                48788.84678,
                60305.81927,
                77063.4884,
                117214.7546,
                151018.0699,
                195698.5419,
            ],
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
        },
        (2, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                11702.08515,
                39023.21173,
                62134.84678,
                77539.81927,
                100113.4884,
                154792.7546,
                201326.0699,
                263374.5419,
            ],
        },
        (3, 1): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                343228.652,
                482425.4684,
                618182.0594,
                743750.2902,
                844443.0443,
                978479.5225,
                1183834.522,
                1440660.587,
            ],
        },
        (3, 2): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                643228.652,
                782425.4684,
                918182.0594,
                1043750.2902,
                1144443.0443,
                1278479.5225,
                1483834.522,
                1740660.587,
            ],
        },
        (3, 3): {
            "Flowrates": [
                0.0,
                36480.0,
                634240.0,
                1434800.0,
                2083760.0,
                3171200.0,
                6342400.0,
                9513600.0,
                14270400.0,
            ],
            "Costs": [
                0.0,
                423074.7216,
                3042779.121,
                5348359.01,
                6921261.68,
                9251002.61,
                14933803.33,
                19762044.75,
                26151302.79,
            ],
        },
    }
    with pytest.raises(ValueError):
        check_discretized_costing_params(m, incorrect_discretized_equipment_cost12)


def test_check_environmental_impact_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)
    # add feed parameters
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # add superstructure formulations parameters
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )
    # Add operating parameters
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
    # Add discretized costing parameters.
    add_discretized_costing_params(m, discretized_equipment_cost)
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

    ## test correct inputs
    check_environmental_impact_params(
        m, consider_environmental_impacts, options_environmental_impacts, epsilon
    )

    ## test incorrect inputs 1 - consider_environmental_impacts is not of type bool
    incorrect_consider_environmental_impacts1 = "True"
    with pytest.raises(TypeError):
        check_environmental_impact_params(
            m,
            incorrect_consider_environmental_impacts1,
            options_environmental_impacts,
            epsilon,
        )

    ## test incorrect inputs 2 - options_environmental_impacts is not of type dict
    incorrect_options_environmental_impacts2 = []
    with pytest.raises(TypeError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            incorrect_options_environmental_impacts2,
            epsilon,
        )

    ## test incorrect inputs 3 - options_environmental_impcts is an empty dict
    incorrect_options_environmental_impacts3 = {}
    with pytest.raises(ValueError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            incorrect_options_environmental_impacts3,
            epsilon,
        )

    ## test incorrect inputs 4 - key in options_environmental_impcts is not of type tuple
    incorrect_options_environmental_impacts4 = {
        1: 0,
        (2, 1): 0,
        (2, 2): 1000,
        (3, 1): 600,
        (3, 2): 0,
        (3, 3): 600,
    }
    with pytest.raises(TypeError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            incorrect_options_environmental_impacts4,
            epsilon,
        )

    ## test incorrect inputs 5 - value in options_environmental_impcts is not of type int or float
    incorrect_options_environmental_impacts5 = {
        (1, 1): "0",
        (2, 1): 0,
        (2, 2): 1000,
        (3, 1): 600,
        (3, 2): 0,
        (3, 3): 600,
    }
    with pytest.raises(TypeError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            incorrect_options_environmental_impacts5,
            epsilon,
        )

    ## test incorrect inputs 6 - epsilon is not of type int or float
    incorrect_epsilon6 = "1"
    with pytest.raises(TypeError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            options_environmental_impacts,
            incorrect_epsilon6,
        )

    ## test incorrect inputs 7 - environmental impacts matrix doesn't contain entries for all options in the superstructure
    incorrect_options_environmental_impacts7 = {
        (2, 1): 0,
        (2, 2): 1000,
        (3, 1): 600,
        (3, 2): 0,
        (3, 3): 600,
    }
    with pytest.raises(ValueError):
        check_environmental_impact_params(
            m,
            consider_environmental_impacts,
            incorrect_options_environmental_impacts7,
            epsilon,
        )


def test_check_byproduct_valorization_params():
    # create new model
    m = pyo.ConcreteModel()
    # add objective function
    add_objective_function_choice_param(m, obj_func)
    # add plant lifetime parameters
    add_plant_lifetime_params(m, plant_start, plant_lifetime)
    # add feed parameters
    add_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
    # add superstructure formulations parameters
    add_supe_formulation_params(
        m, num_stages, options_in_stage, option_outlets, option_efficiencies
    )
    # Add operating parameters
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
    # Add discretized costing parameters.
    add_discretized_costing_params(m, discretized_equipment_cost)
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

    # add params, vars, and cons for environmental impacts.
    add_environmental_impact_params(m, options_environmental_impacts, epsilon)
    add_environmental_impact_vars(m)
    add_environmental_impact_cons(m)

    ## test correct inputs
    check_byproduct_valorization_params(
        m, consider_byproduct_valorization, byproduct_values, byproduct_opt_conversions
    )

    ## test incorrect inputs 1 - consider_byproduct_valorization is not of type bool.
    incorrect_consider_byproduct_valorization1 = "True"
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            incorrect_consider_byproduct_valorization1,
            byproduct_values,
            byproduct_opt_conversions,
        )

    ## test incorrect inputs 2 - byproduct_values is not of type dict
    incorrect_byproduct_values2 = []
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            incorrect_byproduct_values2,
            byproduct_opt_conversions,
        )

    ## test incorrect inputs 3 - byproduct_values is an empty dict
    incorrect_byproduct_values3 = {}
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            incorrect_byproduct_values3,
            byproduct_opt_conversions,
        )

    ## test incorrect inputs 4 - key in byproduct_values is not of type str
    incorrect_byproduct_values4 = {
        1: 0.75,
        "Iron oxide": 1,
        "Residue": 0.25,
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            incorrect_byproduct_values4,
            byproduct_opt_conversions,
        )

    ## test incorrect inputs 5 - value in byproduct_values is not of type int or float
    incorrect_byproduct_values5 = {
        "Jarosite": "0.75",
        "Iron oxide": 1,
        "Residue": 0.25,
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            incorrect_byproduct_values5,
            byproduct_opt_conversions,
        )

    ## test incorrect inputs 6 - byproduct_opt_conversions is not of type dict
    incorrect_byproduct_opt_conversions6 = []
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions6,
        )

    ## test incorrect inputs 7 - byproduct_opt_conversions is an empty dict
    incorrect_byproduct_opt_conversions7 = {}
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions7,
        )

    ## test incorrect inputs 8 - inner dict for byproduct_opt_conversions is empty
    incorrect_byproduct_opt_conversions8 = {
        (3, 1): {},
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"Residue": 0.25},
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions8,
        )

    ## test incorrect inputs 9 - incorrect structure for byproduct_opt_conversions. no inner dict
    incorrect_byproduct_opt_conversions9 = {
        (3, 1): [1, 2],
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"Residue": 0.25},
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions9,
        )

    ## test incorrect inputs 10 - incorrect structure for byproduct_opt_conversions. key in inner dict is not of type str.
    incorrect_byproduct_opt_conversions10 = {
        (3, 1): {1: 1},
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"Residue": 0.25},
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions10,
        )

    ## test incorrect inputs 11 - incorrect structure for byproduct_opt_conversions. value in inner dict is not of type int or float.
    incorrect_byproduct_opt_conversions11 = {
        (3, 1): {"Jarosite": "0.75"},
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"Residue": 0.25},
    }
    with pytest.raises(TypeError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions11,
        )

    ## test incorrect inputs 12 - byproduct_opt_conversions contains infeasible options (options that don't exist in superstructure).
    incorrect_byproduct_opt_conversions12 = {
        (3, 1): {"Jarosite": 0.75},
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"Residue": 0.25},
        (3, 4): {"Residue": 0.25},
    }
    with pytest.raises(ValueError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions12,
        )

    ## test incorrect inputs 13 - byproduct_opt_conversions contains undefined byproducts (not defined in byproduct_values).
    incorrect_byproduct_opt_conversions13 = {
        (3, 1): {"Jarosite": 0.75},
        (3, 2): {"Iron oxide": 1},
        (3, 3): {"acid": 0.25},
    }
    with pytest.raises(ValueError):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            byproduct_values,
            incorrect_byproduct_opt_conversions13,
        )

    ## test incorrect inputs 14 - some byproducts defined in byproduct_values are not produced by any options
    incorrect_byproduct_values14 = {
        "Jarosite": -0.17,
        "Iron oxide": 10,
        "Residue": -0.17,
        "acid": -0.2,
    }
    with pytest.warns(UserWarning):
        check_byproduct_valorization_params(
            m,
            consider_byproduct_valorization,
            incorrect_byproduct_values14,
            byproduct_opt_conversions,
        )
