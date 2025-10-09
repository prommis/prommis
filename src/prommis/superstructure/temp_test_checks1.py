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
plant_lifetime = 11

### Feed Params
available_feed = {
    2025: 20000,
    2026: 20000,
    2027: 20000,
    2028: 20000,
    2029: 20000,
    2030: 20000,
    2031: 20000,
    2032: 20000,
    2033: 20000,
    2034: 20000,
}
collection_rate = 0.6
tracked_comps = ["Nd", "Fe"]
prod_comp_mass = {
    "Nd": 7.5e-4,
    "Fe": 0.00175,
}

### Superstructure Formulation Parameters
num_stages = 4
options_in_stage = {
    1: 2,
    2: 3,
    3: 3,
    4: 4,
}
option_outlets = {
    # level 1
    (1, 1): [1, 2, 3],
    (1, 2): [1, 2, 3],
    # level 2
    (2, 1): [2, 3],
    (2, 2): [2, 3],
    (2, 3): [1],
    # level 3
    (3, 1): [1],
    (3, 2): [2],
    (3, 3): [3, 4],
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
    (3, 3): {"Nd": 1, "Fe": 1},
    # level 4 yields
    (4, 1): {"Nd": 1, "Fe": 0},
    (4, 2): {"Nd": 0.98, "Fe": 0},
    (4, 3): {"Nd": 0.98, "Fe": 0},
    (4, 4): {"Nd": 0.98, "Fe": 0},
}

### Operating Parameters
profit = {
    (4, 1): {"Nd": 70, "Fe": 0},
    (4, 2): {"Nd": 70, "Fe": 0},
    (4, 3): {"Nd": 70, "Fe": 0},
    (4, 4): {"Nd": 70, "Fe": 0},
}
opt_var_oc_params = {
    # level 2
    (2, 1): {"a": 0.0053, "b": 7929.7 / 1000},
    (2, 2): {"a": 0.0015, "b": 2233.16 / 1000},
    (2, 3): {"a": 0.0034, "b": 0},
    # level 3
    (3, 1): {"a": 15.594, "b": 4e6 / 1000},
    (3, 2): {"a": 35.845, "b": 4e6 / 1000},
    (3, 3): {"a": 0, "b": 0},
    # level 4
    (4, 1): {"a": 0.4997, "b": 898320 / 1000},
    (4, 2): {"a": 9.8352, "b": 677751 / 1000},
    (4, 3): {"a": 2.17, "b": 0},
    (4, 4): {"a": 6.808, "b": 0},
}
operators_per_discrete_unit = {
    (1, 1): 1,
    (1, 2): 0,
}
yearly_cost_per_unit = {
    (1, 1): 0,
    (1, 2): 280 / 1000,
}
capital_cost_per_unit = {
    (1, 1): 0,
    (1, 2): 200000 / 1000,
}
processing_rate = {
    (1, 1): 181132 / 1000,
    (1, 2): 523636 / 1000,
}
num_operators = {
    (2, 1): 0.65,
    (2, 2): 0.65,
    (2, 3): 0.65,
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
    (4, 1): 0,
    (4, 2): 0.75,
    (4, 3): 1.15,
    (4, 4): 1.15,
}
labor_rate = 8000 * 38.20 / 1000

### Discretized costing Parameters
discretized_purchased_equipment_cost = {
    (2, 1): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            10.13,
            31.353,
            48.789,
            60.306,
        ],
    },
    (2, 2): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            11.702,
            39.023,
            62.135,
            77.54,
        ],
    },
    (2, 3): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            67.799,
            272.739,
            406.588,
            488.109,
        ],
    },
    (3, 1): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            343.229,
            482.425,
            618.182,
            743.75,
        ],
    },
    (3, 2): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            643.229,
            782.425,
            918.182,
            1043.75,
        ],
    },
    (3, 3): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    },
    (4, 1): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            197.813,
            324.151,
            468.748,
            547.368,
        ],
    },
    (4, 2): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            222.906,
            354.009,
            490.597,
            562.047,
        ],
    },
    (4, 3): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            154.6847193,
            354.009,
            490.597,
            562.047,
        ],
    },
    (4, 4): {
        "Flowrates": [
            0.0,
            36.48,
            634.24,
            1434.8,
            2083.76,
        ],
        "Costs": [
            0.0,
            404.6847193,
            604.009,
            740.597,
            812.047,
        ],
    },
}

### Environmnetal Impact Parameters
consider_environmental_impacts = False

### Byproduct Valorization Parameters
consider_byproduct_valorization = False

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
    None,
    None,
    ### Byproduct valorization parameters
    consider_byproduct_valorization,
    None,
    None,
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

# m.fs.costing.discrete_units_per_option.display()

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
