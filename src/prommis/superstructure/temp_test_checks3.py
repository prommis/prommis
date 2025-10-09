# import pyomo.environ as pyo
# from pyomo.util.check_units import assert_units_consistent

# from idaes.core.scaling.util import report_scaling_factors
from idaes.core.solvers import get_solver

from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice
from prommis.superstructure.report_superstructure_results import (  # report_superstructure_costing,; report_superstructure_environmental_impacts,; report_superstructure_streams,
    report_superstructure_costing,
    report_superstructure_results_overview,
)
from prommis.superstructure.superstructure_function import (
    SuperstructureScaler,
    build_model,
)

#################################################################################################
### Choice of objectie function
obj_func = ObjectiveFunctionChoice.NET_PRESENT_VALUE
# obj_func = ObjectiveFunctionChoice.COST_OF_RECOVERY

### Plant Lifetime Params
plant_start = 2024
plant_lifetime = 15

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

### Superstructure Formulation Parameters
num_stages = 5
options_in_stage = {
    1: 2,
    2: 4,
    3: 4,
    4: 4,
    5: 4,
}
option_outlets = {
    # level 1
    (1, 1): [1, 2, 3, 4],
    (1, 2): [1, 2, 3, 4],
    # level 2
    (2, 1): [1, 2, 4],
    (2, 2): [1, 2, 4],
    (2, 3): [1, 2, 4],
    (2, 4): [3],
    # level 3
    (3, 1): [1],
    (3, 2): [2, 3],
    (3, 3): [2, 3],
    (3, 4): [4],
    # level 4
    (4, 1): [1],
    (4, 2): [2],
    (4, 3): [3],
    (4, 4): [4],
}
option_efficiencies = {
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
    (3, 2): {"Nd": 0.925, "Dy": 0.98, "Fe": 0},
    (3, 3): {"Nd": 1, "Dy": 1, "Fe": 0},
    (3, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    # level 4 yields
    (4, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 2): {"Nd": 1, "Dy": 0.899, "Fe": 0},
    (4, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    # level 5 yields
    (5, 1): {"Nd": 1, "Dy": 1, "Fe": 0},  # grouped together with hydrom. extract.
    (5, 2): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 3): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 4): {
        "Nd": 0.98,
        "Dy": 0.98,
        "Fe": 0,
    },  # grouped together with Cu(NO3)2 diss.
}

### Operating Parameters
profit = {
    (5, 1): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
    (5, 2): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
    (5, 3): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
    (5, 4): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
}
opt_var_oc_params = {
    # level 2
    (2, 1): {"a": 0.0061, "b": 0},
    (2, 2): {"a": 0.0017, "b": 0},
    (2, 3): {"a": 0.0205, "b": 0},
    (2, 4): {"a": 54192, "b": 0},
    # level 3
    (3, 1): {"a": 15.987, "b": 0},
    (3, 2): {"a": 1.8359, "b": 0},
    (3, 3): {"a": 3.7416, "b": 0},
    (3, 4): {"a": 1.58, "b": 0},
    # level 4
    (4, 1): {"a": 0, "b": 0},
    (4, 2): {"a": 111.11, "b": 0},
    (4, 3): {"a": 0, "b": 0},
    (4, 4): {"a": 0, "b": 0},
    # level 5
    (5, 1): {"a": 0.5093, "b": 0},
    (5, 2): {"a": 8.6505, "b": 0},
    (5, 3): {"a": 8.6505, "b": 0},
    (5, 4): {"a": 2.17, "b": 0},
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
    (3, 2): 1.3,
    (3, 3): 0.45,
    (3, 4): 1.15,
    (4, 1): 0,
    (4, 2): 1.3,
    (4, 3): 0,
    (4, 4): 0,
    (5, 1): 1.05,
    (5, 2): 0.75,
    (5, 3): 0.75,
    (5, 4): 1.15,
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
    (2, 4): {
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
            121352.0,
            490387.0,
            732185.0,
            879652.0,
            1081651.0,
            1522296.0,
            1859740.0,
            2272553.0,
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
            226790.0,
            446435.0,
            713714.0,
            1270105.0,
            1541353.0,
            2920751.0,
            3652064.0,
            5323087.0,
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
            169491.0,
            940300.0,
            1534578.0,
            1919653.0,
            2469724.0,
            3743401.0,
            4774426.0,
            6089420.0,
        ],
    },
    (4, 1): {
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
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    },
    (4, 2): {
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
            450073.0,
            1623424.0,
            2349169.0,
            2778687.0,
            3349410.0,
            4585349.0,
            5503177.0,
            6590434.0,
        ],
    },
    (4, 3): {
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
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    },
    (4, 4): {
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
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    },
    (5, 1): {
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
            197813.1853,
            324151.478,
            468747.756,
            547368.2774,
            655614.4213,
            800184.1752,
            974415.7068,
            1114534.98,
        ],
    },
    (5, 2): {
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
            222906.0,
            354009.0,
            490597.0,
            562047.0,
            679397.0,
            912244.0,
            1097498.0,
            1297052.0,
        ],
    },
    (5, 3): {
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
            222906.0,
            354009.0,
            490597.0,
            562047.0,
            679397.0,
            912244.0,
            1097498.0,
            1297052.0,
        ],
    },
    (5, 4): {
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
            154685.0,
            858157.0,
            1400520.0,
            1751956.0,
            2253973.0,
            3416384.0,
            4357340.0,
            5557458.0,
        ],
    },
}

### Environmnetal Impact Parameters
consider_environmental_impacts = False
options_environmental_impacts = {
    (1, 1): 0,
    (1, 2): 1000,
    (2, 1): 0,
    (2, 2): 1000,
    (2, 3): 600,
    (2, 4): 800,
    (3, 1): 600,
    (3, 2): 0,
    (3, 3): 600,
    (3, 4): 800,
    (4, 1): 0,
    (4, 2): 800,
    (4, 3): 600,
    (4, 4): 1000,
    (5, 1): 0,
    (5, 2): 800,
    (5, 3): 600,
    (5, 4): 800,
}
epsilon = 1e16
# epsilon = 1

### Byproduct Valorization Parameters
consider_byproduct_valorization = False
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
    (5, 4): {"Iron hydroxide": 0.5},
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

## scale model
scaler = SuperstructureScaler()
scaler.scale_model(m)


## Solve model
solver = get_solver(solver="gurobi")
solver.options["NumericFocus"] = 3
results = solver.solve(m, tee=True)

report_superstructure_results_overview(m, results)
report_superstructure_costing(m, results)
# report_superstructure_streams(m, results)
# report_superstructure_environmental_impacts(m, results)

m.fs.costing.discrete_units_per_option.display()

# m.fs.costing.obj.display()

### Print out results
# m.fs.option_binary_var.display()

# assert_units_consistent(m.fs)
# m.fs.option_efficiencies.display()
# m.fs.byproduct_valorization.opt_byproduct_set.display()
# m.fs.byproduct_valorization.byproduct_opt_conversion.display()
# m.fs.costing.byproduct_values.display()

# this is wat is causing issues
# m.fs.byproduct_valorization.byproduct_producing_opts.display()

# check blocks
# for block in m.component_objects(pyo.Block, descend_into=True):
#     try:
#         assert_units_consistent(block)
#         print(f"Units consistent for block: {block.name}")
#     except TypeError as e:
#         print(f"Skipped block {block.name}: {e}")

# check constraints
# for constr in m.component_objects(pyo.Constraint, descend_into=True):
#     try:
#         assert_units_consistent(constr)
#         print(f"Units consistent for constraint: {constr.name}")
#     except TypeError as e:
#         print(f"Skipped constraint {constr.name}: {e}")
#     except Exception as e:
#         print(f"Units inconsistent in constraint {constr.name}: {e}")
# for constr in m.component_objects(pyo.Constraint, descend_into=True):
#     assert_units_consistent(constr)

# # check expressions
# for expr in m.component_objects(pyo.Expression, descend_into=True):
#     try:
#         assert_units_consistent(expr)
#         print(f"Units consistent for expression: {expr.name}")
#     except TypeError as e:
#         print(f"Skipped expression {expr.name}: {e}")
#     except Exception as e:
#         print(f"Units inconsistent in expression {expr.name}: {e}")
# for expr in m.component_objects(pyo.Expression, descend_into=True):
#     assert_units_consistent(expr)

# m.fs.available_feed.display()
# m.fs.feed_entering.display()

# m.fs.num_stages.display()
# m.fs.stages_set.display()
# m.fs.options_in_stage.display()
# m.fs.all_opts_set.display()
# m.fs.discrete_opts_set.display()
# m.fs.continuous_opts_set.display()
# m.fs.option_outlets.display()
# m.fs.option_efficiencies.display()
# m.fs.final_opts_set.display()

# Print out sets to check for redundancies
# add plant lifetime params
# m.fs.plant_life_range.display()
# m.fs.operational_range.display()

# add feed params
# m.fs.tracked_comps.display()

# add superstructure formulation params
# m.fs.max_options.display()
# m.fs.stages_set.display()
# m.fs.all_opts_set.display()
# m.fs.discrete_opts_set.display()
# m.fs.continuous_opts_set.display()
# m.fs.option_outlet_pairs.display()
# m.fs.final_opts_set.display()

# add operating params

# add discretized costing params
# m.fs.continuous_opts_discretized_costing_data_points.display()

# add mass balance params
# m.fs.f_stages.display()

# m.fs.option_efficiencies.display()

# dt = DiagnosticsToolbox(m)

# dt.display_components_with_inconsistent_units()

# from pyomo.util.check_units import assert_units_consistent

# # Call the scaling report
# from idaes.core.scaling.util import report_scaling_factors

# assert_units_consistent(m)

# for block in m.fs.costing.piecewise_cons.values():
#     for c in block.component_data_objects(pyo.Constraint, descend_into=True):
#         c.deactivate()


# report_scaling_factors(m.fs)
# report_scaling_factors(m.fs.costing)
# report_scaling_factors(m.fs.byproduct_valorization)
# report_scaling_factors(m.fs.environmental_impacts)

# import logging
# from io import StringIO
# from pyomo.environ import Constraint
# from pyomo.common.log import LoggingIntercept
# from pyomo.util.infeasible import (
# log_active_constraints,
# log_close_to_bounds,
# log_infeasible_bounds,
# log_infeasible_constraints,)
# output = StringIO()
# with LoggingIntercept(output, "pyomo.util.infeasible", logging.INFO):
#     log_infeasible_constraints(m)
# # print(output.getvalue().splitlines())
# for line in output.getvalue().splitlines():
#     print(line)

# m.display()
