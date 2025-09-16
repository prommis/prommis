# import pyomo.environ as pyo
# from pyomo.util.check_units import assert_units_consistent

# from idaes.core.scaling.util import report_scaling_factors
from idaes.core.solvers import get_solver

from prommis.superstructure.superstructure_function import (
    ObjectiveFunction,
    SuperstructureScaler,
    build_model,
)

# from idaes.core.util import DiagnosticsToolbox


#################################################################################################
### Choice of objectie function
obj_func = ObjectiveFunction.cost_of_recovery

### Plant Lifetime Params
plant_start = 2024
plant_lifetime = 6

### Feed Params
available_feed = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
    2029: 592637,
}
collection_rate = 0.1
tracked_comps = ["Nd", "Fe"]
prod_comp_mass = {
    "Nd": 0.309 * 3,
    "Fe": 0.691 * 3,
}

### Superstructure Formulation Parameters
num_stages = 3
options_in_stage = {
    1: 2,
    2: 3,
    3: 4,
}
option_outlets = {
    # level 1
    (1, 1): [1, 2, 3],
    (1, 2): [1, 2, 3],
    # level 2
    (2, 1): [2, 3, 4],
    (2, 2): [1, 2, 3],
    (2, 3): [1, 2],
}

option_efficiencies = {
    # Level 1 yields
    (1, 1): {"Nd": 1, "Fe": 1},
    (1, 2): {"Nd": 1, "Fe": 1},
    # level 2 yields
    (2, 1): {"Nd": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Fe": 1},
    (2, 3): {"Nd": 1, "Fe": 1},
    # level 3 yields
    (3, 1): {"Nd": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Fe": 0},
    (3, 3): {"Nd": 0.925, "Fe": 0},
    (3, 4): {"Nd": 1, "Fe": 0},
}

### Operating Parameters
profit = {
    (3, 1): {"Nd": 45.4272, "Fe": 0},
    (3, 2): {"Nd": 69.888, "Fe": 0},
    (3, 3): {"Nd": 45.4272, "Fe": 0},
    (3, 4): {"Nd": 100, "Fe": 0},
}
opt_var_oc_params = {
    # level 2
    (2, 1): {"a": 0.0053, "b": 7929.7},
    (2, 2): {"a": 0.0015, "b": 2233.16},
    (2, 3): {"a": 0.0034, "b": 0},
    # level 3
    (3, 1): {"a": 15.594, "b": 4e6},
    (3, 2): {"a": 35.58463, "b": 4e6},
    (3, 3): {"a": 1.8359, "b": 0},
    (3, 4): {"a": 3.7414, "b": 2378.6},
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
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
    (3, 4): 0.45,
}
labor_rate = 8000 * 38.20

### Discretized costing Parameters
discretized_purchased_equipment_cost = {
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
            101,
            313,
            487,
            603,
            770,
            1172,
            1510,
            1956,
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
    (2, 3): {
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
            67799.0,
            272739.0,
            406588.0,
            488109.0,
            599667.0,
            842665.0,
            1028482.0,
            1255541.0,
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
            226.0,
            446.0,
            713.0,
            1270.0,
            1541.0,
            2920.0,
            3652.0,
            5323.0,
        ],
    },
}

### Objective Function Parameters
# obj_func = "NPV"
# obj_func = "COR"

### Environmnetal Impact Parameters
consider_environmental_impacts = True
options_environmental_impacts = {
    (1, 1): 0,
    (1, 2): 1000,
    (2, 1): 0,
    (2, 2): 1000,
    (2, 3): 600,
    (3, 1): 600,
    (3, 2): 0,
    (3, 3): 600,
    (3, 4): 800,
}
epsilon = 1e12
# epsilon = 1

### Byproduct Valorization Parameters
consider_byproduct_valorization = True
byproduct_values = {
    "Jarosite": -0.17,
    "Iron oxide": 10,
    "Residue": -0.17,
    "Iron hydroxide": -0.17,
}
byproduct_opt_conversions = {
    (3, 1): {"Jarosite": 0.75},
    (3, 2): {"Iron oxide": 1},
    (3, 3): {"Residue": 0.25},
    (3, 4): {"Iron hydroxide": 0.5},
}

#################################################################################################
### Build model
m = build_model(
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
)

scaler = SuperstructureScaler()
scaler.scale_model(m)

solver = get_solver(solver="gurobi")
solver.options["NumericFocus"] = 2
results = solver.solve(m, tee="True")


# report_scaling_factors(m.fs)
# report_scaling_factors(m.fs.costing)
# report_scaling_factors(m.fs.byproduct_valorization)
# report_scaling_factors(m.fs.environmental_impacts)

# m.fs.costing.obj.display()
# m.fs.option_binary_var.display()

# for idx, v in m.fs.f.items():
#     comp = idx[1]
#     print('idx', idx)
#     print('comp', comp)
#     print('v', v)
