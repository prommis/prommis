import pyomo.environ as pyo

from prommis.superstructures.superstructure_function import (
    add_feed_params_block,
    add_plant_lifetime_params_block,
    add_supe_formulation_params,
    add_operating_params,
    check_feed_params,
    check_operating_params,
    check_plant_lifetime_params,
    check_supe_formulation_params,
)

### Build model
m = pyo.ConcreteModel()


### Plant Lifetime Params
plant_start = 2024
plant_lifetime = 15

check_plant_lifetime_params(plant_lifetime)
add_plant_lifetime_params_block(m, plant_start, plant_lifetime)

### Feed Params
available_feed = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
    2029: 592637,
    2030: 571054,
    2031: 498472,
    2032: 506565,
    2033: 566355,
    2034: 669094,
    2035: 719057,
    2036: 762656,
    2037: 1434637,
    2038: 1697805,
}
collection_rate = 0.1
tracked_comps = ["Nd", "Dy", "Fe"]
prod_comp_mass = {
    "Nd": 0.206 * 3,
    "Dy": 0.103 * 3,
    "Fe": 0.691 * 3,
}

check_feed_params(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)
add_feed_params_block(m, available_feed, collection_rate, tracked_comps, prod_comp_mass)

### Superstructure Formulation Parameters
num_stages = 5
options_in_stage = {
    1: 2,
    2: 4,
    3: 6,
    4: 4,
    5: 5,
}
option_outlets = {
    # level 1
    (1, 1): [1, 2, 3, 4],
    (1, 2): [1, 2, 3, 4],
    # level 2
    (2, 1): [1, 2, 3, 6],
    (2, 2): [1, 2, 3, 6],
    (2, 3): [1, 2, 3, 6],
    (2, 4): [4, 5],
    # level 3
    (3, 1): [1],
    (3, 2): [1],
    (3, 3): [2, 3],
    (3, 4): [2, 3],
    (3, 5): [2, 3],
    (3, 6): [4],
    # level 4
    (4, 1): [1],
    (4, 2): [2],
    (4, 3): [3],
    (4, 4): [4, 5],
}
discrete_opts = [(1, 1), (1, 2)]
option_eff = {
    # Level 1 yields
    (1, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (1, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
    # level 2 yields
    (2, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    # level 3 yields
    (3, 1): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
    (3, 3): {"Nd": 0.925, "Dy": 0.98, "Fe": 0},
    (3, 4): {"Nd": 1, "Dy": 1, "Fe": 0},
    (3, 5): {"Nd": 1, "Dy": 1, "Fe": 0},
    (3, 6): {"Nd": 1, "Dy": 1, "Fe": 0.403},
    # level 4 yields
    (4, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 2): {"Nd": 1, "Dy": 0.899, "Fe": 0},
    (4, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    # level 5 yields
    (5, 1): {
        "Nd": 1,
        "Dy": 1,
        "Fe": 0,
    },
    (5, 2): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 3): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 4): {
        "Nd": 0.98,
        "Dy": 0.98,
        "Fe": 0,
    },
    (5, 5): {
        "Nd": 0.98,
        "Dy": 0.98,
        "Fe": 0,
    },
}

check_supe_formulation_params(
    m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff
)
add_supe_formulation_params(m, num_stages, options_in_stage, option_outlets, discrete_opts, option_eff)

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
workers_per_discr_unit = {
    (1, 1): 1,
    (1, 2): 0,
}
yearly_cost_per_unit = {
    (1, 1): 0,
    (1, 2): 280,
}
cost_per_unit = {
    (1, 1): 0,
    (1, 2): 200000,
}
processing_rate = {
    (1, 1): 7868,
    (1, 2): 52453,
}
num_workers = {
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

check_operating_params(
    m,
    profit,
    opt_var_oc_params,
    workers_per_discr_unit,
    yearly_cost_per_unit,
    cost_per_unit,
    processing_rate,
    num_workers,
    labor_rate,
)
add_operating_params(
    m,
    profit,
    opt_var_oc_params,
    workers_per_discr_unit,
    yearly_cost_per_unit,
    cost_per_unit,
    processing_rate,
    num_workers,
    labor_rate,
)
