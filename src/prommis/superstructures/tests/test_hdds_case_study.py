#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import copy
import math

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Objective,
    RangeSet,
    SolverFactory,
    Var,
    assert_optimal_termination,
    value,
)

from idaes.core.solvers import get_solver

# model statistics
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

import pytest

from prommis.superstructures.superstructure_function import build_model

solver_available = SolverFactory("gurobi").available()
if solver_available:
    solver = get_solver(solver="gurobi")
else:
    solver = None


@pytest.fixture(scope="module")
def get_common_params():
    return {
        ###################################################################################################
        ### Plant Lifetime Parameters
        "plant_start": 2024,  # start of plant production
        "plant_lifetime": 15,  # lifetime of plant
        ###################################################################################################
        ###################################################################################################
        ### Feed parameters
        # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
        "Available_feed": {
            2025: 19482.463,
            2026: 19362.314,
            2027: 17893.125,
            2028: 17651.492,
            2029: 16370.492,
            2030: 13854.916,
            2031: 13284.074,
            2032: 11991.115,
            2033: 10870.423,
            2034: 9787.703,
            2035: 8743.312,
            2036: 7734.329,
            2037: 6758.895,
            2038: 5813.304,
        },
        # collection rate for how much of the available feed is processed by the plant each year
        "CR": 0.6,
        "Tracked_comps": ["Nd", "Fe"],  # tracked components
        # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
        "Prod_comp_mass": {
            "Nd": 7.5e-4,
            "Fe": 0.00175,
        },
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        "numStages": 4,  # number of total stages
        # number of options in each stage
        "Options_in_stage": {
            1: 4,
            2: 4,
            3: 6,
            4: 4,
        },
        # set of options k' in stage j+1 connected to option k in stage j
        "Option_outlets": {
            # level 1
            (1, 1): [1, 2, 3, 4],
            (1, 2): [1, 2, 3, 4],
            (1, 3): [1, 2],
            (1, 4): [1, 2, 3, 4],
            # level 2
            (2, 1): [1, 2, 3, 6],
            (2, 2): [1, 2, 3, 6],
            (2, 3): [1, 2, 3, 6],
            (2, 4): [4, 5],
            # level 3
            (3, 1): [1],
            (3, 2): [1],
            (3, 3): [2],
            (3, 4): [2],
            (3, 5): [2],
            (3, 6): [3, 4],
        },
        # dictionary of tracked component retention efficiency for each option
        "Option_Eff": {
            # Level 1 yields
            (1, 1): {"Nd": 1, "Fe": 1},
            (1, 2): {"Nd": 1, "Fe": 1},
            (1, 3): {"Nd": 1, "Fe": 1},
            (1, 4): {"Nd": 1, "Fe": 1},
            # level 2 yields
            (2, 1): {"Nd": 1, "Fe": 1},
            (2, 2): {"Nd": 1, "Fe": 1},
            (2, 3): {"Nd": 1, "Fe": 1},
            (2, 4): {"Nd": 1, "Fe": 1},
            # level 3 yields
            (3, 1): {"Nd": 0.985, "Fe": 0},
            (3, 2): {"Nd": 0.985, "Fe": 0},
            (3, 3): {"Nd": 0.925, "Fe": 0},
            (3, 4): {"Nd": 1, "Fe": 0},
            (3, 5): {"Nd": 1, "Fe": 0},
            (3, 6): {"Nd": 1, "Fe": 0.403},
            # level 4 yields
            (4, 1): {"Nd": 1, "Fe": 0},
            (4, 2): {"Nd": 0.98, "Fe": 0},
            (4, 3): {"Nd": 0.98, "Fe": 0},
            (4, 4): {"Nd": 0.98, "Fe": 0},
        },
        ###################################################################################################
        ###################################################################################################
        ### Operating Parameters
        # profit (k$) per metric tonnes of product in terms of tracked components
        "Profit": {
            (4, 1): {"Nd": 69.888, "Fe": 0},
            (4, 2): {"Nd": 69.888, "Fe": 0},
            (4, 3): {"Nd": 69.888, "Fe": 0},
            (4, 4): {"Nd": 69.888, "Fe": 0},
        },
        # conversion of kg REE/Fe to kg REO/Fe2O3
        "REE_to_REO_Conversion": {"Nd": 1.664, "Fe": 1.43},
        # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes)
        # OPEX = a*F_in + b*y
        "N_OC_var": {
            # level 2
            (2, 1): {"a": 0.0053, "b": 7929.7 / 1000},
            (2, 2): {"a": 0.0015, "b": 2233.16 / 1000},
            (2, 3): {"a": 0.0034, "b": 0},
            (2, 4): {"a": 0.0117, "b": 0},
            # level 3
            (3, 1): {"a": 15.594, "b": 4e6 / 1000},
            (3, 2): {"a": 35.845, "b": 4e6 / 1000},
            (3, 3): {"a": 1.8359, "b": 0},
            (3, 4): {"a": 3.7414, "b": 3476.7 / 1000},
            (3, 5): {"a": 10.4404, "b": 3476.7 / 1000},
            (3, 6): {"a": 1.58, "b": 0},
            # level 4
            (4, 1): {"a": 0.4997, "b": 898320 / 1000},
            (4, 2): {"a": 9.8352, "b": 677751 / 1000},
            (4, 3): {"a": 2.17, "b": 0},
            (4, 4): {"a": 6.807857755993, "b": 0},
        },
        # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
        "num_workers": {
            ##############################################
            # for disassembly stage, its operators per unit
            (1, 1): 1,
            (1, 2): 0,
            (1, 3): 1,
            (1, 4): 0.3,
            ##############################################
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
            (4, 2): 0.75,
            (4, 3): 1.15,
            (4, 4): 1.15,
        },
        "labor_rate": 8000
        * 38.20
        / 1000,  # yearly wage per type of labor (k$/operator/yr)
        # yearly operating costs per unit (k$/unit/yr)
        "YCU": {
            (1, 1): 0,
            (1, 2): 280 / 1000,
            (1, 3): 0,
            (1, 4): 10367.608 / 1000,
        },
        # cost per disassembly stage unit for each disassembly option (k$/unit)
        "CU": {
            (1, 1): 0,
            (1, 2): 200000 / 1000,
            (1, 3): 0,
            (1, 4): 50000 / 1000,
        },
        # disassembly rate (thousands of EOL HDDs/yr/unit)
        "Dis_Rate": {
            (1, 1): 181132 / 1000,
            (1, 2): 523636 / 1000,
            (1, 3): 960000 / 1000,
            (1, 4): 2.16e7 / 1000,
        },
        ###################################################################################################
        ###################################################################################################
        ### Costing Parameters
        "LF": 2.97,  # Lang Factor
        "TOC_factor": 1.177,  # Overnight costs factor
        "ATWACC": 0.0577,  # discount rate. (default of 5.77%)
        "i_OC_esc": 0.03,  # opex, revenue (default of 3%)
        "i_CAP_esc": 0.036,  # capex escalation rate (default of 3.6%)
        "f_exp": [  # capital expenditure schedule (default of 10%, 60%, 30%)
            0.1,
            0.6,
            0.3,
        ],
        # Define Python Dictionary with discretized cost by flows for each option
        "Discretized_CAPEX": {
            "(2, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 10.13,
                    "2": 31.353,
                    "3": 48.789,
                    "4": 60.306,
                    "5": 77.063,
                    "6": 117.215,
                    "7": 151.018,
                    "8": 195.699,
                },
            },
            "(2, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 11.702,
                    "2": 39.023,
                    "3": 62.135,
                    "4": 77.54,
                    "5": 100.113,
                    "6": 154.793,
                    "7": 201.326,
                    "8": 263.375,
                },
            },
            "(2, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 67.799,
                    "2": 272.739,
                    "3": 406.588,
                    "4": 488.109,
                    "5": 599.667,
                    "6": 842.665,
                    "7": 1028.482,
                    "8": 1255.541,
                },
            },
            "(2, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 121.352,
                    "2": 490.387,
                    "3": 732.185,
                    "4": 879.652,
                    "5": 1081.651,
                    "6": 1522.296,
                    "7": 1859.74,
                    "8": 2272.553,
                },
            },
            "(3, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 343.229,
                    "2": 482.425,
                    "3": 618.182,
                    "4": 743.75,
                    "5": 844.443,
                    "6": 978.48,
                    "7": 1183.835,
                    "8": 1440.661,
                },
            },
            "(3, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 643.229,
                    "2": 782.425,
                    "3": 918.182,
                    "4": 1043.75,
                    "5": 1144.443,
                    "6": 1278.48,
                    "7": 1483.835,
                    "8": 1740.661,
                },
            },
            "(3, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 423.075,
                    "2": 3042.779,
                    "3": 5348.359,
                    "4": 6921.262,
                    "5": 9251.003,
                    "6": 14933.803,
                    "7": 19762.045,
                    "8": 26151.303,
                },
            },
            "(3, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 226.79,
                    "2": 446.435,
                    "3": 713.714,
                    "4": 1270.105,
                    "5": 1541.353,
                    "6": 2920.751,
                    "7": 3652.064,
                    "8": 5323.087,
                },
            },
            "(3, 5)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 476.79,
                    "2": 696.435,
                    "3": 963.714,
                    "4": 1520.105,
                    "5": 1791.353,
                    "6": 3170.751,
                    "7": 3902.064,
                    "8": 5573.087,
                },
            },
            "(3, 6)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 169.491,
                    "2": 940.3,
                    "3": 1534.578,
                    "4": 1919.653,
                    "5": 2469.724,
                    "6": 3743.401,
                    "7": 4774.426,
                    "8": 6089.42,
                },
            },
            "(4, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 197.813,
                    "2": 324.151,
                    "3": 468.748,
                    "4": 547.368,
                    "5": 655.614,
                    "6": 800.184,
                    "7": 974.416,
                    "8": 1114.535,
                },
            },
            "(4, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 222.906,
                    "2": 354.009,
                    "3": 490.597,
                    "4": 562.047,
                    "5": 679.397,
                    "6": 912.244,
                    "7": 1097.498,
                    "8": 1297.052,
                },
            },
            "(4, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 154.6847193,
                    "2": 354.009,
                    "3": 490.597,
                    "4": 562.047,
                    "5": 679.397,
                    "6": 912.244,
                    "7": 1097.498,
                    "8": 1297.052,
                },
            },
            "(4, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36.48,
                    "2": 634.24,
                    "3": 1434.8,
                    "4": 2083.76,
                    "5": 3171.2,
                    "6": 6342.4,
                    "7": 9513.6,
                    "8": 14270.4,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 404.6847193,
                    "2": 604.009,
                    "3": 740.597,
                    "4": 812.047,
                    "5": 929.397,
                    "6": 1162.244,
                    "7": 1347.498,
                    "8": 1547.052,
                },
            },
        },
        ###################################################################################################
        ###################################################################################################
        # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
        "environ_impacts": {
            (1, 1): 0,
            (1, 2): 600,
            (1, 3): 800,
            (1, 4): 1000,
            (2, 1): 0,
            (2, 2): 1000,
            (2, 3): 800,
            (2, 4): 600,
            (3, 1): 1000,
            (3, 2): 0,
            (3, 3): 800,
            (3, 4): 600,
            (3, 5): 600,
            (3, 6): 1200,
            (4, 1): 0,
            (4, 2): 800,
            (4, 3): 1200,
            (4, 4): 1000,
        },
        "epsilon": 1e16,  # epsilon factor for generating Pareto front
        ###################################################################################################
        ###################################################################################################
        ### Byproduct valorization
        # list of byproducts
        "byprods": ["Jarosite", "Iron oxide", "Residue", "Iron hydroxide"],
        # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
        "byprod_vals": {
            "Jarosite": -0.17,
            "Iron oxide": 10,
            "Residue": -0.17,
            "Iron hydroxide": -0.17,
        },
        # dictionary keeping track of which tracked component produces which byproduct
        "tracked_comp_for_byprod": {
            "Jarosite": "Fe",
            "Iron oxide": "Fe",
            "Residue": "Fe",
            "Iron hydroxide": "Fe",
        },
        # dictionary tracking which options produce a given byproduct
        "byprod_options": {
            "Jarosite": [(3, 1)],
            "Iron oxide": [(3, 2), (3, 5), (3, 6), (4, 4)],
            "Residue": [(3, 3)],
            "Iron hydroxide": [(3, 4), (4, 3)],
        },
        # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
        "byprod_options_eff": {
            "Jarosite": {(3, 1): 1},
            "Iron oxide": {(3, 2): 1, (3, 5): 1, (3, 6): 1, (4, 4): 1},
            "Residue": {(3, 3): 1},
            "Iron hydroxide": {(3, 4): 0.597, (4, 3): 1},  # means 40.3% of Fe remains
        },
        # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
        "TC_to_byproduct": {
            "Jarosite": 2.893,
            "Iron oxide": 1.430,
            "Residue": 1,
            "Iron hydroxide": 1.914,
        },
    }


class TestNPV(object):
    @pytest.fixture(scope="class")
    def NPV_model(self, get_common_params):
        common_params = get_common_params
        m = build_model(
            ###################################################################################################
            ### Plant Lifetime Parameters
            plant_start=common_params["plant_start"],  # start of plant production
            plant_lifetime=common_params["plant_lifetime"],  # lifetime of plant
            ###################################################################################################
            ###################################################################################################
            ### Feed parameters
            # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
            Available_feed=common_params["Available_feed"],
            # collection rate for how much of the available feed is processed by the plant each year
            CR=common_params["CR"],
            Tracked_comps=common_params["Tracked_comps"],  # tracked components
            # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
            Prod_comp_mass=common_params["Prod_comp_mass"],
            ###################################################################################################
            ###################################################################################################
            ### Superstructure formulation parameters
            numStages=common_params["numStages"],  # number of total stages
            Options_in_stage=common_params[  # number of options in each stage
                "Options_in_stage"
            ],
            # set of options k' in stage j+1 connected to option k in stage j
            Option_outlets=common_params[  # set of options k' in stage j+1 connected to option k in stage j
                "Option_outlets"
            ],
            # dictionary of tracked component retention efficiency for each option
            Option_Eff=common_params["Option_Eff"],
            ###################################################################################################
            ###################################################################################################
            ### Operating Parameters
            # profit (k$) per metric tonnes of product in terms of tracked components
            Profit=common_params["Profit"],
            # conversion of kg REE/Fe to kg REO/Fe2O3
            REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
            # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes)
            # OPEX = a*F_in + b*y
            N_OC_var=common_params["N_OC_var"],
            # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
            num_workers=common_params["num_workers"],
            labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
            # yearly operating costs per unit (k$/unit/yr)
            YCU=common_params["YCU"],
            # cost per disassembly stage unit for each disassembly option (k$/unit)
            CU=common_params["CU"],
            # disassembly rate (thousands of EOL HDDs/yr/unit)
            Dis_Rate=common_params["Dis_Rate"],
            ###################################################################################################
            ###################################################################################################
            ### Costing Parameters
            LF=common_params["LF"],  # Lang Factor
            TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
            ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
            i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
            i_CAP_esc=common_params[  # capex escalation rate (default of 3.6%)
                "i_CAP_esc"
            ],
            f_exp=common_params[  # capital expenditure schedule (default of 10%, 60%, 30%)
                "f_exp"
            ],
            # Define Python Dictionary with discretized cost by flows for each option.
            Discretized_CAPEX=common_params["Discretized_CAPEX"],
            ###################################################################################################
            ###################################################################################################
            ### Choice of objective function. Options are 'NPV' or 'COR'.capitalize
            obj_func="NPV",
            ###################################################################################################
            ###################################################################################################
            ### Consideration of environmental impacts parameters
            # boolean to decide whether or not to consider environmental impacts
            consider_environ_impacts=False,
            # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
            environ_impacts=common_params["environ_impacts"],
            epsilon=common_params[  # epsilon factor for generating Pareto front
                "epsilon"
            ],
            ###################################################################################################
            ###################################################################################################
            ### Byproduct valorization
            # boolean to decide whether or not to consider the valorization of byproducts
            consider_byprod_val=False,
            # list of byproducts
            byprods=common_params["byprods"],
            # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
            byprod_vals=common_params["byprod_vals"],
            # dictionary keeping track of which tracked component produces which byproduct
            tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
            # dictionary tracking which options produce a given byproduct
            byprod_options=common_params["byprod_options"],
            # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
            byprod_options_eff=common_params["byprod_options_eff"],
            # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
            TC_to_byproduct=common_params["TC_to_byproduct"],
        )

        return m

    def test_build(self, NPV_model, get_common_params):
        # start of plant production
        prod_start = get_common_params["plant_start"] + 1
        # final year plant is in production
        plant_end = (
            get_common_params["plant_start"] + get_common_params["plant_lifetime"] - 1
        )
        # plant operational period
        operational_range = RangeSet(prod_start, plant_end)

        assert isinstance(NPV_model, ConcreteModel)

        assert isinstance(NPV_model.plantYear, Block)

        for t in operational_range:
            assert isinstance(NPV_model.plantYear[t].P_entering, Var)
            assert isinstance(NPV_model.plantYear[t].F_in, Var)
            assert isinstance(NPV_model.plantYear[t].F_out, Var)

            assert isinstance(NPV_model.plantYear[t].init_flow_cons, Constraint)
            assert isinstance(NPV_model.plantYear[t].inlet_flow_cons, Constraint)
            assert isinstance(NPV_model.plantYear[t].intermediate_flow_cons, Constraint)
            assert isinstance(NPV_model.plantYear[t].outlet_flow_cons, Constraint)

        assert isinstance(NPV_model.binOpt, Var)

        assert isinstance(NPV_model.stage_bin_cons, Constraint)
        assert isinstance(NPV_model.connect_bin_cons, Constraint)
        assert isinstance(NPV_model.big_M_cons, Constraint)

        assert isinstance(NPV_model.DisOptWorkers, Var)

        assert isinstance(NPV_model.OC_var_cons, Constraint)
        assert isinstance(NPV_model.DisWorkerCons, Constraint)
        assert isinstance(NPV_model.profit_opt_cons, Constraint)

        for t in operational_range:
            assert isinstance(NPV_model.plantYear[t].OC_var, Var)
            assert isinstance(NPV_model.plantYear[t].OC_var_total, Var)
            assert isinstance(NPV_model.plantYear[t].ProfitOpt, Var)
            assert isinstance(NPV_model.plantYear[t].Profit, Var)

            assert isinstance(NPV_model.plantYear[t].profit_con, Constraint)

        assert isinstance(NPV_model.workers, Var)

        assert isinstance(NPV_model.worker_cons, Constraint)
        assert isinstance(NPV_model.COL_cons, Constraint)

        assert isinstance(NPV_model.bin_workers, Var)
        assert isinstance(NPV_model.total_workers, Var)
        assert isinstance(NPV_model.COL_Total, Var)

        assert isinstance(NPV_model.COL_Total_con, Constraint)

        assert isinstance(NPV_model.BEC, Var)

        assert isinstance(NPV_model.BEC_cons, Constraint)

        assert isinstance(NPV_model.BEC_max_flow, Var)

        assert isinstance(NPV_model.BEC_max_flow_cons, Constraint)

        assert isinstance(NPV_model.TPC, Var)

        assert isinstance(NPV_model.TPC_cons, Constraint)

        assert isinstance(NPV_model.Total_TPC, Var)

        assert isinstance(NPV_model.Total_TPC_con, Constraint)

        assert isinstance(NPV_model.TOC, Var)

        assert isinstance(NPV_model.TOC_con, Constraint)

        assert isinstance(NPV_model.node_TOC, Var)

        assert isinstance(NPV_model.node_TOC_cons, Constraint)

        assert isinstance(NPV_model.CF, Var)

        assert isinstance(NPV_model.CF_cons, Constraint)

        assert isinstance(NPV_model.TOC_exp, Var)

        assert isinstance(NPV_model.TOC_exp_cons, Constraint)

        assert isinstance(NPV_model.GE, Var)

        assert isinstance(NPV_model.GE_cons, Constraint)

        assert isinstance(NPV_model.Rev, Var)

        assert isinstance(NPV_model.Rev_cons, Constraint)

        assert isinstance(NPV_model.OC_fixed, Var)

        assert isinstance(NPV_model.OC_fixed_cons, Constraint)

        assert isinstance(NPV_model.OC_var, Var)

        assert isinstance(NPV_model.OH, Var)

        assert isinstance(NPV_model.OH_cons, Constraint)

        assert isinstance(NPV_model.bin_test_cons, Constraint)

        assert isinstance(NPV_model.obj, Objective)

        assert number_variables(NPV_model) == 1961
        assert number_total_constraints(NPV_model) == 2609
        assert number_unused_variables(NPV_model) == 0
        assert degrees_of_freedom(NPV_model) == 641

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, NPV_model):
        solver.options["NumericFocus"] = 3

        results = solver.solve(NPV_model)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, NPV_model, get_common_params):
        ##### Define Parameters for tests
        ###################################################################################################
        ### Plant Lifetime Parameters
        plant_start = get_common_params["plant_start"]  # start of plant production
        plant_lifetime = get_common_params["plant_lifetime"]  # lifetime of plant
        # first year is construction
        prod_start = plant_start + 1
        plant_end = (  # final year plant is in production
            plant_start + plant_lifetime - 1
        )
        # total plant lifetime
        plant_life_range = RangeSet(plant_start, plant_end)
        # operational lifetime of the plant
        operational_range = RangeSet(prod_start, plant_end)
        ###################################################################################################
        ###################################################################################################
        ### Feed parameters
        # Total feedstock available for recycling each year
        Available_feed = get_common_params["Available_feed"]
        # collection rate for how much of the available feed is processed by the plant each year (fraction)
        CR = get_common_params["CR"]
        # list of tracked components
        Tracked_comps = get_common_params["Tracked_comps"]
        # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
        Prod_comp_mass = get_common_params["Prod_comp_mass"]
        # calculate feed entering parameter based on yearly available feedstock and collection rate
        Feed_entering = copy.deepcopy(Available_feed)
        for key in Feed_entering:
            Feed_entering[key] = Available_feed[key] * CR
        # calculate max feed that can enter the plant
        maxFeedEntering = max(  # max feed entering plant over production period
            Feed_entering.values()
        )
        maxFeedEnteringYear = max(  # year in which max feed enters plant
            Feed_entering, key=Feed_entering.get
        )
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        numStages = get_common_params["numStages"]  # no. of total stages
        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]
        # set of options k' in stage j+1 connected to option k in stage j
        Option_outlets = get_common_params["Option_outlets"]
        # dictionary of tracked component retention efficiency for each option
        Option_Eff = get_common_params["Option_Eff"]
        maxOptions = max(Options_in_stage.values())  # max options in any of the stages
        # make a list of the options in the final stage
        final_opt_list = [(numStages, j) for j in RangeSet(Options_in_stage[numStages])]
        ###################################################################################################
        ###################################################################################################
        ### Operating Parameters
        # profit (k$) per metric tonnes of product in terms of tracked components
        Profit = get_common_params["Profit"]
        # conversion of kg REE/Fe to kg REO/Fe2O3
        REE_to_REO_Conversion = get_common_params["REE_to_REO_Conversion"]
        # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes).
        # OPEX = a*F_in + b*y
        N_OC_var = get_common_params["N_OC_var"]
        # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
        num_workers = get_common_params["num_workers"]
        labor_rate = get_common_params["labor_rate"]  # yearly wage per type of labor
        # yearly operating costs per unit (k$/unit/yr)
        YCU = get_common_params["YCU"]
        # cost per disassembly stage unit for each disassembly option (k$/unit)
        CU = get_common_params["CU"]
        # disassembly rate (thousands of EOL HDDs/yr/unit)
        Dis_Rate = get_common_params["Dis_Rate"]
        # calculate max disassembly units possible for each option
        max_dis_by_option = copy.deepcopy(Dis_Rate)
        for key in max_dis_by_option.keys():
            max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
        max_dis_workers = max(max_dis_by_option.values())
        # calculate max possible workers for the process
        max_workers = max_dis_workers + numStages * math.ceil(max(num_workers.values()))
        ###################################################################################################
        ###################################################################################################
        ### Costing Parameters
        LF = get_common_params["LF"]  # Lang Factor
        TOC_factor = get_common_params["TOC_factor"]  # Overnight costs factor
        ATWACC = get_common_params["ATWACC"]  # discount rate
        i_OC_esc = get_common_params["i_OC_esc"]  # opex, revenue escalation rate
        i_CAP_esc = get_common_params["i_CAP_esc"]  # capex escalation rate
        f_exp = get_common_params["f_exp"]  # capital expenditure schedule
        # Define Python Dictionary with discretized cost by flows for each option.
        Discretized_CAPEX = get_common_params["Discretized_CAPEX"]
        ###################################################################################################
        ###################################################################################################
        ### Other Parameters needed for tests
        # make a list of all options
        all_opts_list = []
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                all_opts_list.append((j, k))
        # list of stages that should be chosen for optimal process
        opt_stages = [(1, 4), (2, 2), (3, 6), (4, 3)]
        # number of manual disassembly workers
        man_dis_workers = 0
        # number of automatic disassembly units
        auto_dis_units = 0
        # number of mechanical fracturing units
        mech_frac_units = 0
        # number of shredder units
        shred_units = 1
        # list of disassembly workers/units tuple (stage, option, no. dis. workers/units)
        dis_stage_workers = [
            (1, 1, man_dis_workers),
            (1, 2, auto_dis_units),
            (1, 3, mech_frac_units),
            (1, 4, shred_units),
        ]
        # define number of workers for overall process
        total_workers = 4
        # create dictionary to hold yearly flow values
        yearly_F_vals = {
            key1: {
                key2: {key3: None for key3 in Tracked_comps}
                for key2 in RangeSet(numStages - 1)
            }
            for key1 in operational_range
        }
        yearly_F_in_vals = {
            key1: {
                key2: {key3: None for key3 in Tracked_comps} for key2 in all_opts_list
            }
            for key1 in operational_range
        }
        yearly_F_out_vals = {
            key1: {
                key2: {key3: None for key3 in Tracked_comps} for key2 in all_opts_list
            }
            for key1 in operational_range
        }
        # create dictionary to track yearly total operating costs
        OC_var = {key: None for key in plant_life_range}
        OC_var[plant_start] = 0  # operation doesn't start until second year
        # create dictionary to track yearly total profit
        revenue = {key: 0 for key in plant_life_range}
        revenue[plant_start] = 0  # operation doesn't start until second year
        # cost of labor
        COL_total = total_workers * labor_rate
        # BEC max flows for options chosen in optimal pathway (except for disassembly stage)
        BEC_max_flows = {
            (2, 2): 29.2236945,
            (3, 6): 29.2236945,
            (4, 3): 17.01111256845,
        }
        # BEC for options chosen in optimal pathway
        BECs = {
            (1, 4): 50,
            (2, 2): 9.374333142516448,
            (3, 6): 135.77722600053457,
            (4, 3): 72.13155626730786,
        }
        # TPC of each element
        elem_TPC = {
            (1, 4): BECs[(1, 4)],
            (2, 2): LF * BECs[(2, 2)],
            (3, 6): LF * BECs[(3, 6)],
            (4, 3): LF * BECs[(4, 3)],
        }
        # TOC of each element
        elem_TOC = {
            (1, 4): TOC_factor * elem_TPC[(1, 4)],
            (2, 2): TOC_factor * elem_TPC[(2, 2)],
            (3, 6): TOC_factor * elem_TPC[(3, 6)],
            (4, 3): TOC_factor * elem_TPC[(4, 3)],
        }
        # TPC
        Total_TPC = sum(elem_TPC.values())
        # TOC
        TOC = Total_TPC * TOC_factor
        # create dictionary to track yearly TOC expenditure
        TOC_exp = {key: None for key in plant_life_range}
        # create dictionary to track yearly fixed operatinc costs
        OC_fixed = {key: None for key in plant_life_range}
        # create dictionary to track yearly plant overhead
        OH = {key: None for key in plant_life_range}
        # create dictionary to track yearly gross earnings
        GE = {key: None for key in plant_life_range}
        # create dictionary to track yearly cash flows
        CF = {key: None for key in plant_life_range}
        ###################################################################################################

        for t in operational_range:
            # test P_entering for each year
            assert value(NPV_model.plantYear[t].P_entering) == pytest.approx(
                Feed_entering[t], rel=1e-8
            )

            # test all F for each year
            for j in RangeSet(numStages - 1):
                for c in Tracked_comps:
                    F_val = (
                        Feed_entering[t]
                        * Prod_comp_mass[c]
                        * math.prod(
                            Option_Eff[opt_stages[stage]][c]
                            for stage in RangeSet(0, j - 1)
                        )
                    )
                    assert value(NPV_model.plantYear[t].F[j, c]) == pytest.approx(
                        F_val,
                        rel=1e-8,
                    )
                    yearly_F_vals[t][j][c] = F_val

            # test all F_in for each year
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        F_in_val = 0
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            F_in_val = 0
                            assert value(
                                NPV_model.plantYear[t].F_in[j, k, c]
                            ) == pytest.approx(F_in_val, rel=1e-8)
                            yearly_F_in_vals[t][(j, k)][c] = F_in_val

                        else:
                            if j == 1:
                                F_in_val = Feed_entering[t] * Prod_comp_mass[c]
                                assert value(
                                    NPV_model.plantYear[t].F_in[j, k, c]
                                ) == pytest.approx(F_in_val, rel=1e-8)
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val

                            else:
                                F_in_val = (
                                    Feed_entering[t]
                                    * Prod_comp_mass[c]
                                    * math.prod(
                                        Option_Eff[opt_stages[stage]][c]
                                        for stage in RangeSet(0, j - 2)
                                    )
                                )
                                assert value(
                                    NPV_model.plantYear[t].F_in[j, k, c]
                                ) == pytest.approx(F_in_val, rel=1e-8)
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val

            # test all F_out for each year
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        F_out_val = 0
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            F_out_val = 0
                            assert value(
                                NPV_model.plantYear[t].F_out[j, k, c]
                            ) == pytest.approx(F_out_val, rel=1e-8)
                            yearly_F_out_vals[t][(j, k)][c] = F_out_val

                        else:
                            F_out_val = (
                                Feed_entering[t]
                                * Prod_comp_mass[c]
                                * math.prod(
                                    Option_Eff[opt_stages[stage]][c]
                                    for stage in RangeSet(0, j - 1)
                                )
                            )
                            assert value(
                                NPV_model.plantYear[t].F_out[j, k, c]
                            ) == pytest.approx(F_out_val, rel=1e-8)
                            yearly_F_out_vals[t][(j, k)][c] = F_out_val

        # test binary variables
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                if (j, k) in opt_stages:
                    assert value(NPV_model.binOpt[j, k]) == pytest.approx(1, rel=1e-8)
                else:
                    assert value(NPV_model.binOpt[j, k]) == pytest.approx(0, rel=1e-8)

        # test number of disassembly workers
        j_dis = 1
        for k_dis in RangeSet(Options_in_stage[j_dis]):
            for num_worker in RangeSet(max_dis_by_option[j_dis, k_dis]):
                if (j_dis, k_dis, num_worker) in dis_stage_workers:
                    assert value(
                        NPV_model.DisOptWorkers[j_dis, k_dis, num_worker]
                    ) == pytest.approx(1, rel=1e-8)

                else:
                    assert value(
                        NPV_model.DisOptWorkers[j_dis, k_dis, num_worker]
                    ) == pytest.approx(0, rel=1e-8)

        # test variable cost for each option
        for t in operational_range:
            tot_OC_var = 0
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    elem_OC_var = 0
                    # checking disassembly stages
                    if j == 1:
                        if k == 4:
                            elem_OC_var = shred_units * YCU[(j, k)]
                            assert value(
                                NPV_model.plantYear[t].OC_var[(j, k)]
                            ) == pytest.approx(elem_OC_var, rel=1e-8)
                            tot_OC_var += elem_OC_var
                        else:
                            elem_OC_var = 0
                            assert value(
                                NPV_model.plantYear[t].OC_var[(j, k)]
                            ) == pytest.approx(elem_OC_var, rel=1e-8)

                    # checking rest
                    if j != 1:
                        if (j, k) not in opt_stages:
                            elem_OC_var = 0
                            assert value(
                                NPV_model.plantYear[t].OC_var[(j, k)]
                            ) == pytest.approx(elem_OC_var, rel=1e-8)

                        else:
                            elem_OC_var = (
                                N_OC_var[(j, k)]["a"]
                                * sum(
                                    yearly_F_in_vals[t][(j, k)][c]
                                    for c in Tracked_comps
                                )
                                + N_OC_var[(j, k)]["b"]
                            )
                            assert value(
                                NPV_model.plantYear[t].OC_var[(j, k)]
                            ) == pytest.approx(elem_OC_var, rel=1e-8)
                            tot_OC_var += elem_OC_var
            OC_var[t] = tot_OC_var

        # test total yearly variable operating costs
        for t in operational_range:
            assert value(NPV_model.plantYear[t].OC_var_total) == pytest.approx(
                OC_var[t], rel=1e-8
            )

        # test profit from each opt in final stage
        for t in operational_range:
            j = 4
            for k in RangeSet(Options_in_stage[j]):
                total_profit = sum(
                    yearly_F_out_vals[t][(j, k)][c] * Profit[(j, k)][c]
                    for c in Tracked_comps
                )
                assert value(NPV_model.plantYear[t].ProfitOpt[(j, k)]) == pytest.approx(
                    total_profit, rel=1e-8
                )
                revenue[t] += total_profit

        # test total yearly profit
        for t in operational_range:
            assert value(NPV_model.plantYear[t].Profit) == pytest.approx(revenue[t])

        # test number of workers for each option
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                elem = (j, k)

                if elem not in opt_stages:
                    assert value(NPV_model.workers[elem]) == pytest.approx(0, rel=1e-8)
                else:
                    assert value(NPV_model.workers[elem]) == pytest.approx(
                        num_workers[elem], rel=1e-8
                    )

        # test number of workers
        for i in RangeSet(0, max_workers):
            if i == total_workers:
                assert value(NPV_model.bin_workers[i]) == pytest.approx(1, rel=1e-8)

            else:
                assert value(NPV_model.bin_workers[i]) == pytest.approx(0, rel=1e-8)

        # test total workers
        assert value(NPV_model.total_workers) == pytest.approx(total_workers, rel=1e-8)

        # test cost of labor (COL)
        assert value(NPV_model.COL_Total) == pytest.approx(COL_total, rel=1e-8)

        # test max flow vars and BECs
        for j in RangeSet(1, numStages):
            for k in RangeSet(Options_in_stage[j]):
                elem = (j, k)

                if j == 1:
                    if elem == (1, 4):
                        assert value(NPV_model.BEC[elem]) == pytest.approx(
                            BECs[elem], rel=1e-8
                        )
                    else:
                        assert value(NPV_model.BEC[elem]) == pytest.approx(0, rel=1e-8)

                else:
                    if elem not in opt_stages:
                        assert value(NPV_model.BEC_max_flow[elem]) == pytest.approx(
                            0, rel=1e-8
                        )
                        assert value(NPV_model.BEC[elem]) == pytest.approx(0, rel=1e-8)

                    elif elem == (2, 2):
                        assert value(NPV_model.BEC_max_flow[elem]) == pytest.approx(
                            BEC_max_flows[elem], rel=1e-8
                        )
                        assert value(NPV_model.BEC[elem]) == pytest.approx(
                            BECs[elem], rel=1e-8
                        )

                    elif elem == (3, 6):
                        assert value(NPV_model.BEC_max_flow[elem]) == pytest.approx(
                            BEC_max_flows[elem], rel=1e-8
                        )
                        assert value(NPV_model.BEC[elem]) == pytest.approx(
                            BECs[elem], rel=1e-8
                        )

                    elif elem == (4, 4):
                        assert value(NPV_model.BEC_max_flow[elem]) == pytest.approx(
                            BEC_max_flows[elem], rel=1e-8
                        )
                        assert value(NPV_model.BEC[elem]) == pytest.approx(
                            BECs[elem], rel=1e-8
                        )

        # test TPC
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                elem = (j, k)

                if elem not in opt_stages:
                    assert value(NPV_model.TPC[elem]) == pytest.approx(0, rel=1e-8)

                elif elem == (1, 4):
                    assert value(NPV_model.TPC[elem]) == pytest.approx(
                        elem_TPC[elem], rel=1e-8
                    )

                elif elem == (2, 2):
                    assert value(NPV_model.TPC[elem]) == pytest.approx(
                        elem_TPC[elem], rel=1e-8
                    )

                elif elem == (3, 6):
                    assert value(NPV_model.TPC[elem]) == pytest.approx(
                        elem_TPC[elem], rel=1e-8
                    )

                elif elem == (4, 3):
                    assert value(NPV_model.TPC[elem]) == pytest.approx(
                        elem_TPC[elem], rel=1e-8
                    )

        # test total TPC
        assert value(NPV_model.Total_TPC) == pytest.approx(Total_TPC, rel=1e-8)

        # test TOC
        assert value(NPV_model.TOC) == pytest.approx(TOC, rel=1e-8)

        # test node TOCs
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                elem = (j, k)

                if elem not in opt_stages:
                    assert value(NPV_model.node_TOC[elem]) == pytest.approx(0, rel=1e-8)

                elif elem == (1, 4):
                    assert value(NPV_model.node_TOC[elem]) == pytest.approx(
                        elem_TOC[elem], rel=1e-8
                    )

                elif elem == (2, 2):
                    assert value(NPV_model.node_TOC[elem]) == pytest.approx(
                        elem_TOC[elem], rel=1e-8
                    )

                elif elem == (3, 6):
                    assert value(NPV_model.node_TOC[elem]) == pytest.approx(
                        elem_TOC[elem], rel=1e-8
                    )

                elif elem == (4, 3):
                    assert value(NPV_model.node_TOC[elem]) == pytest.approx(
                        elem_TOC[elem], rel=1e-8
                    )

        # test TOC expenditure
        for t in plant_life_range:
            if t < plant_start + 3:  # capital expended over first three years
                TOC_exp[t] = (
                    ((1 + i_CAP_esc) ** (t - plant_start))
                    * f_exp[t - plant_start]
                    * TOC
                )
            else:
                TOC_exp[t] = 0

            assert value(NPV_model.TOC_exp[t]) == pytest.approx(TOC_exp[t], rel=1e-8)

        # test fixed and variable operating costs, revenue
        for t in plant_life_range:
            if t == plant_start:
                OC_fixed[t] = 0
                assert value(NPV_model.OC_fixed[t]) == pytest.approx(
                    OC_fixed[t], rel=1e-8
                )
                assert value(NPV_model.OC_var[t]) == pytest.approx(OC_var[t], rel=1e-8)
                assert value(NPV_model.Rev[t]) == pytest.approx(revenue[t], rel=1e-8)
            else:
                OC_fixed[t] = 1.55 * COL_total + 0.03 * Total_TPC + 0.01 * revenue[t]
                assert value(NPV_model.OC_fixed[t]) == pytest.approx(
                    OC_fixed[t], rel=1e-8
                )
                assert value(NPV_model.OC_var[t]) == pytest.approx(OC_var[t], rel=1e-8)
                assert value(NPV_model.Rev[t]) == pytest.approx(revenue[t], rel=1e-8)

        # test plant overhead, gross earnings, and cash flows
        for t in plant_life_range:
            OH[t] = 0.2 * (OC_fixed[t] + OC_var[t])
            GE[t] = (revenue[t] - OC_fixed[t] - OC_var[t] - OH[t]) * (1 + i_OC_esc) ** (
                t - 2024
            )
            CF[t] = GE[t] - TOC_exp[t]

            assert value(NPV_model.OH[t]) == pytest.approx(OH[t], rel=1e-8)
            assert value(NPV_model.GE[t]) == pytest.approx(GE[t], rel=1e-8)
            assert value(NPV_model.CF[t]) == pytest.approx(CF[t], rel=1e-8)

        # test NPV objective
        NPV = sum(CF[t] / ((1 + ATWACC) ** (t - plant_start)) for t in plant_life_range)
        assert value(NPV_model.obj) == pytest.approx(NPV, rel=1e-8)


class TestCOR(object):
    @pytest.fixture(scope="class")
    def COR_model(self, get_common_params):
        common_params = get_common_params
        m = build_model(
            ###################################################################################################
            ### Plant Lifetime Parameters
            plant_start=common_params["plant_start"],  # start of plant production
            plant_lifetime=common_params["plant_lifetime"],  # lifetime of plant
            ###################################################################################################
            ###################################################################################################
            ### Feed parameters
            # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
            Available_feed=common_params["Available_feed"],
            # collection rate for how much of the available feed is processed by the plant each year
            CR=common_params["CR"],
            Tracked_comps=common_params["Tracked_comps"],  # tracked components
            # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
            Prod_comp_mass=common_params["Prod_comp_mass"],
            ###################################################################################################
            ###################################################################################################
            ### Superstructure formulation parameters
            numStages=common_params["numStages"],  # number of total stages
            Options_in_stage=common_params[  # number of options in each stage
                "Options_in_stage"
            ],
            # set of options k' in stage j+1 connected to option k in stage j
            Option_outlets=common_params[  # set of options k' in stage j+1 connected to option k in stage j
                "Option_outlets"
            ],
            # dictionary of tracked component retention efficiency for each option
            Option_Eff=common_params["Option_Eff"],
            ###################################################################################################
            ###################################################################################################
            ### Operating Parameters
            # profit (k$) per metric tonnes of product in terms of tracked components
            Profit=common_params["Profit"],
            # conversion of kg REE/Fe to kg REO/Fe2O3
            REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
            # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes).
            # OPEX = a*F_in + b*y
            N_OC_var=common_params["N_OC_var"],
            # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
            num_workers=common_params["num_workers"],
            labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
            # yearly operating costs per unit (k$/unit/yr)
            YCU=common_params["YCU"],
            # cost per disassembly stage unit for each disassembly option (k$/unit)
            CU=common_params["CU"],
            # disassembly rate (thousands of EOL HDDs/yr/unit)
            Dis_Rate=common_params["Dis_Rate"],
            ###################################################################################################
            ###################################################################################################
            ### Costing Parameters
            LF=common_params["LF"],  # Lang Factor
            TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
            ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
            i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
            i_CAP_esc=common_params[  # capex escalation rate (default of 3.6%)
                "i_CAP_esc"
            ],
            f_exp=common_params[  # capital expenditure schedule (default of 10%, 60%, 30%)
                "f_exp"
            ],
            # Define Python Dictionary with discretized cost by flows for each option.
            Discretized_CAPEX=common_params["Discretized_CAPEX"],
            ###################################################################################################
            ###################################################################################################
            ### Choice of objective function. Options are 'NPV' or 'COR'.capitalize
            obj_func="COR",
            ###################################################################################################
            ###################################################################################################
            ### Consideration of environmental impacts parameters
            # boolean to decide whether or not to consider environmental impacts
            consider_environ_impacts=False,
            # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
            environ_impacts=common_params["environ_impacts"],
            epsilon=common_params[  # epsilon factor for generating Pareto front
                "epsilon"
            ],
            ###################################################################################################
            ###################################################################################################
            ### Byproduct valorization
            # boolean to decide whether or not to consider the valorization of byproducts
            consider_byprod_val=False,
            # list of byproducts
            byprods=common_params["byprods"],
            # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
            byprod_vals=common_params["byprod_vals"],
            # dictionary keeping track of which tracked component produces which byproduct
            tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
            # dictionary tracking which options produce a given byproduct
            byprod_options=common_params["byprod_options"],
            # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
            byprod_options_eff=common_params["byprod_options_eff"],
            # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
            TC_to_byproduct=common_params["TC_to_byproduct"],
        )

        return m

    def test_build(self, COR_model, get_common_params):
        # start of plant production
        prod_start = get_common_params["plant_start"] + 1
        # final year plant is in production
        plant_end = (
            get_common_params["plant_start"] + get_common_params["plant_lifetime"] - 1
        )
        # plant operational period
        operational_range = RangeSet(prod_start, plant_end)

        assert isinstance(COR_model.COR, Var)
        assert isinstance(COR_model.NPV, Var)
        assert isinstance(COR_model.NPV_con1, Constraint)
        assert isinstance(COR_model.NPV_con2, Constraint)

        assert number_variables(COR_model) == 1963
        assert number_total_constraints(COR_model) == 2611
        assert number_unused_variables(COR_model) == 0
        assert degrees_of_freedom(COR_model) == 641

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, COR_model):
        solver.options["NumericFocus"] = 3

        results = solver.solve(COR_model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, COR_model, get_common_params):
        ##### Define Parameters for tests
        ###################################################################################################
        ### Plant Lifetime Parameters
        plant_start = get_common_params["plant_start"]  # start of plant production
        plant_lifetime = get_common_params["plant_lifetime"]  # lifetime of plant
        # first year is construction
        prod_start = plant_start + 1
        plant_end = (  # final year plant is in production
            plant_start + plant_lifetime - 1
        )
        # total plant lifetime
        plant_life_range = RangeSet(plant_start, plant_end)
        # operational lifetime of the plant
        operational_range = RangeSet(prod_start, plant_end)
        ###################################################################################################
        ###################################################################################################
        ### Feed parameters
        # Total feedstock available for recycling each year
        Available_feed = get_common_params["Available_feed"]
        # collection rate for how much of the available feed is processed by the plant each year (fraction)
        CR = get_common_params["CR"]
        # list of tracked components
        Tracked_comps = get_common_params["Tracked_comps"]
        # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
        Prod_comp_mass = get_common_params["Prod_comp_mass"]
        # calculate feed entering parameter based on yearly available feedstock and collection rate
        Feed_entering = copy.deepcopy(Available_feed)
        for key in Feed_entering:
            Feed_entering[key] = Available_feed[key] * CR
        # calculate max feed that can enter the plant
        maxFeedEntering = max(  # max feed entering plant over production period
            Feed_entering.values()
        )
        maxFeedEnteringYear = max(  # year in which max feed enters plant
            Feed_entering, key=Feed_entering.get
        )
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        numStages = get_common_params["numStages"]  # no. of total stages
        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]
        # set of options k' in stage j+1 connected to option k in stage j
        Option_outlets = get_common_params["Option_outlets"]
        # dictionary of tracked component retention efficiency for each option
        Option_Eff = get_common_params["Option_Eff"]
        maxOptions = max(Options_in_stage.values())  # max options in any of the stages
        # make a list of the options in the final stage
        final_opt_list = [(numStages, j) for j in RangeSet(Options_in_stage[numStages])]
        ###################################################################################################
        ###################################################################################################
        ### Operating Parameters
        # profit (k$) per metric tonnes of product in terms of tracked components
        Profit = get_common_params["Profit"]
        # conversion of kg REE/Fe to kg REO/Fe2O3
        REE_to_REO_Conversion = get_common_params["REE_to_REO_Conversion"]
        # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes).
        # OPEX = a*F_in + b*y
        N_OC_var = get_common_params["N_OC_var"]
        # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
        num_workers = get_common_params["num_workers"]
        labor_rate = get_common_params["labor_rate"]  # yearly wage per type of labor
        # yearly operating costs per unit (k$/unit/yr)
        YCU = get_common_params["YCU"]
        # cost per disassembly stage unit for each disassembly option (k$/unit)
        CU = get_common_params["CU"]
        # disassembly rate (thousands of EOL HDDs/yr/unit)
        Dis_Rate = get_common_params["Dis_Rate"]
        # calculate max disassembly units possible for each option
        max_dis_by_option = copy.deepcopy(Dis_Rate)
        for key in max_dis_by_option.keys():
            max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
        max_dis_workers = max(max_dis_by_option.values())
        # calculate max possible workers for the process
        max_workers = max_dis_workers + numStages * math.ceil(max(num_workers.values()))
        ###################################################################################################
        ###################################################################################################
        ### Costing Parameters
        LF = get_common_params["LF"]  # Lang Factor
        TOC_factor = get_common_params["TOC_factor"]  # Overnight costs factor
        ATWACC = get_common_params["ATWACC"]  # discount rate
        i_OC_esc = get_common_params["i_OC_esc"]  # opex, revenue escalation rate
        i_CAP_esc = get_common_params["i_CAP_esc"]  # capex escalation rate
        f_exp = get_common_params["f_exp"]  # capital expenditure schedule
        # Define Python Dictionary with discretized cost by flows for each option.
        Discretized_CAPEX = get_common_params["Discretized_CAPEX"]
        ###################################################################################################
        ###################################################################################################
        ### Other Parameters needed for tests
        # cost of recovery
        COR = 253.93284856
        # profit for the final option in the optimal pathway by year (k$ per year)
        profit = {
            2025: 3630.401484,
            2026: 3608.012676,
            2027: 3334.241033,
            2028: 3289.214652,
            2029: 3050.510526,
            2030: 2581.753016,
            2031: 2475.381165,
            2032: 2234.448575,
            2033: 2025.616566,
            2034: 1823.860335,
            2035: 1629.246408,
            2036: 1441.230479,
            2037: 1259.466138,
            2038: 1083.262802,
        }

        # test COR
        assert value(COR_model.COR) == pytest.approx(COR, rel=1e-8)

        # test the profit from each option in the final stage
        for t in operational_range:
            j = 4
            for k in RangeSet(Options_in_stage[j]):
                if (j, k) == (4, 3):
                    assert value(
                        COR_model.plantYear[t].ProfitOpt[(j, k)]
                    ) == pytest.approx(profit[t], rel=1e-8)
                else:
                    assert value(
                        COR_model.plantYear[t].ProfitOpt[(j, k)]
                    ) == pytest.approx(0, rel=1e-8)

        # test the NPV
        assert value(COR_model.NPV) == pytest.approx(0, rel=1e-8)

        # test the objective function
        assert value(COR_model.obj) == pytest.approx(COR, rel=1e-8)


class TestEnvironmentalImpacts(object):
    @pytest.fixture(scope="class")
    def EI_model(self, get_common_params):
        common_params = get_common_params
        m = build_model(
            ###################################################################################################
            ### Plant Lifetime Parameters
            plant_start=common_params["plant_start"],  # start of plant production
            plant_lifetime=common_params["plant_lifetime"],  # lifetime of plant
            ###################################################################################################
            ###################################################################################################
            ### Feed parameters
            # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
            Available_feed=common_params["Available_feed"],
            # collection rate for how much of the available feed is processed by the plant each year
            CR=common_params["CR"],
            Tracked_comps=common_params["Tracked_comps"],  # tracked components
            # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
            Prod_comp_mass=common_params["Prod_comp_mass"],
            ###################################################################################################
            ###################################################################################################
            ### Superstructure formulation parameters
            numStages=common_params["numStages"],  # number of total stages
            Options_in_stage=common_params[  # number of options in each stage
                "Options_in_stage"
            ],
            # set of options k' in stage j+1 connected to option k in stage j
            Option_outlets=common_params[  # set of options k' in stage j+1 connected to option k in stage j
                "Option_outlets"
            ],
            # dictionary of tracked component retention efficiency for each option
            Option_Eff=common_params["Option_Eff"],
            ###################################################################################################
            ###################################################################################################
            ### Operating Parameters
            # profit (k$) per metric tonnes of product in terms of tracked components
            Profit=common_params["Profit"],
            # conversion of kg REE/Fe to kg REO/Fe2O3
            REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
            # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes).
            # OPEX = a*F_in + b*y
            N_OC_var=common_params["N_OC_var"],
            # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
            num_workers=common_params["num_workers"],
            labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
            # yearly operating costs per unit (k$/unit/yr)
            YCU=common_params["YCU"],
            # cost per disassembly stage unit for each disassembly option
            CU=common_params["CU"],
            # disassembly rate (thousands of EOL HDDs/yr/unit)
            Dis_Rate=common_params["Dis_Rate"],
            ###################################################################################################
            ###################################################################################################
            ### Costing Parameters
            LF=common_params["LF"],  # Lang Factor
            TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
            ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
            i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
            i_CAP_esc=common_params[  # capex escalation rate (default of 3.6%)
                "i_CAP_esc"
            ],
            f_exp=common_params[  # capital expenditure schedule (default of 10%, 60%, 30%)
                "f_exp"
            ],
            # Define Python Dictionary with discretized cost by flows for each option.
            Discretized_CAPEX=common_params["Discretized_CAPEX"],
            ###################################################################################################
            ###################################################################################################
            ### Choice of objective function. Options are 'NPV' or 'COR'.capitalize
            obj_func="NPV",
            ###################################################################################################
            ###################################################################################################
            ### Consideration of environmental impacts parameters
            # boolean to decide whether or not to consider environmental impacts
            consider_environ_impacts=True,
            # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
            environ_impacts=common_params["environ_impacts"],
            epsilon=common_params[  # epsilon factor for generating Pareto front
                "epsilon"
            ],
            ###################################################################################################
            ###################################################################################################
            ### Byproduct valorization
            # boolean to decide whether or not to consider the valorization of byproducts
            consider_byprod_val=False,
            # list of byproducts
            byprods=common_params["byprods"],
            # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
            byprod_vals=common_params["byprod_vals"],
            # dictionary keeping track of which tracked component produces which byproduct
            tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
            # dictionary tracking which options produce a given byproduct
            byprod_options=common_params["byprod_options"],
            # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
            byprod_options_eff=common_params["byprod_options_eff"],
            # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
            TC_to_byproduct=common_params["TC_to_byproduct"],
        )

        return m

    def test_build(self, EI_model, get_common_params):
        # start of plant production
        prod_start = get_common_params["plant_start"] + 1
        # final year plant is in production
        plant_end = (
            get_common_params["plant_start"] + get_common_params["plant_lifetime"] - 1
        )
        # plant operational period
        operational_range = RangeSet(prod_start, plant_end)
        # no. of total stages
        numStages = get_common_params["numStages"]
        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]

        for t in operational_range:
            assert isinstance(EI_model.plantYear[t].total_yearly_GWP, Var)
            assert isinstance(EI_model.GWP_cons, Constraint)

            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    print((j, k))
                    assert isinstance(EI_model.plantYear[t].yearly_GWP, Var)

        assert isinstance(EI_model.GWP, Var)
        assert isinstance(EI_model.GWP_cons, Constraint)
        assert isinstance(EI_model.epsilon_con, Constraint)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, EI_model):
        solver.options["NumericFocus"] = 3

        results = solver.solve(EI_model)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, EI_model, get_common_params):
        ##### Define Parameters for tests
        ###################################################################################################
        ### Plant Lifetime Parameters
        plant_start = get_common_params["plant_start"]  # start of plant production
        plant_lifetime = get_common_params["plant_lifetime"]  # lifetime of plant
        # first year is construction
        prod_start = plant_start + 1
        plant_end = (  # final year plant is in production
            plant_start + plant_lifetime - 1
        )
        # total plant lifetime
        plant_life_range = RangeSet(plant_start, plant_end)
        # operational lifetime of the plant
        operational_range = RangeSet(prod_start, plant_end)
        ###################################################################################################
        ###################################################################################################
        ### Feed parameters
        # Total feedstock available for recycling each year
        Available_feed = get_common_params["Available_feed"]
        # collection rate for how much of the available feed is processed by the plant each year (fraction)
        CR = get_common_params["CR"]
        # list of tracked components
        Tracked_comps = get_common_params["Tracked_comps"]
        # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
        Prod_comp_mass = get_common_params["Prod_comp_mass"]
        # calculate feed entering parameter based on yearly available feedstock and collection rate
        Feed_entering = copy.deepcopy(Available_feed)
        for key in Feed_entering:
            Feed_entering[key] = Available_feed[key] * CR
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        numStages = get_common_params["numStages"]  # no. of total stages
        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]
        # dictionary of tracked component retention efficiency for each option
        Option_Eff = get_common_params["Option_Eff"]
        ###################################################################################################
        ###################################################################################################
        ### Consideration of environmental impacts parameters
        environ_impacts = get_common_params["environ_impacts"]
        ###################################################################################################
        ###################################################################################################
        ### Other Parameters needed for tests
        # make a list of all options
        all_opts_list = []
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                all_opts_list.append((j, k))
        # list of stages that should be chosen for optimal process
        opt_stages = [(1, 4), (2, 2), (3, 6), (4, 3)]
        # create dictionary to hold yearly F_in values
        yearly_F_in_vals = {
            key1: {
                key2: {key3: None for key3 in Tracked_comps} for key2 in all_opts_list
            }
            for key1 in operational_range
        }
        # calculate F_in vals and store in dictionary
        for t in operational_range:
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        F_in_val = 0
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            F_in_val = 0
                            yearly_F_in_vals[t][(j, k)][c] = F_in_val
                        else:
                            if j == 1:
                                F_in_val = Feed_entering[t] * Prod_comp_mass[c]
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val

                            else:
                                F_in_val = (
                                    Feed_entering[t]
                                    * Prod_comp_mass[c]
                                    * math.prod(
                                        Option_Eff[opt_stages[stage]][c]
                                        for stage in RangeSet(0, j - 2)
                                    )
                                )
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val
        # create a dict to store yearly GWP
        yearly_GWP = {
            key1: {key2: None for key2 in all_opts_list} for key1 in operational_range
        }
        # calculate yearly GWP vals and store in dictionary
        for t in operational_range:
            for opt in all_opts_list:
                yearly_GWP[t][opt] = (
                    sum(yearly_F_in_vals[t][opt][c] for c in Tracked_comps)
                    * environ_impacts[opt]
                )
        # create dict to store total yearly GWP
        total_yearly_GWP = {key1: None for key1 in operational_range}
        # calculate total yearly GWP vals and store in dictionary
        for t in operational_range:
            total_yearly_GWP[t] = sum(yearly_GWP[t][opt] for opt in all_opts_list)
        # Calculate GWP
        GWP = sum(total_yearly_GWP[t] for t in operational_range)
        ###################################################################################################

        for t in operational_range:
            for opt in all_opts_list:
                # test yearly GWP
                assert value(EI_model.plantYear[t].yearly_GWP[opt]) == pytest.approx(
                    yearly_GWP[t][opt], rel=1e-8
                )

            # test total yearly GWP
            assert value(EI_model.plantYear[t].total_yearly_GWP) == pytest.approx(
                total_yearly_GWP[t], rel=1e-8
            )

        assert value(EI_model.GWP) == pytest.approx(GWP, rel=1e-8)


class TestByprodVal(object):
    @pytest.fixture(scope="class")
    def BV_model(self, get_common_params):
        common_params = get_common_params
        m = build_model(
            ###################################################################################################
            ### Plant Lifetime Parameters
            plant_start=common_params["plant_start"],  # start of plant production
            plant_lifetime=common_params["plant_lifetime"],  # lifetime of plant
            ###################################################################################################
            ###################################################################################################
            ### Feed parameters
            # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
            Available_feed=common_params["Available_feed"],
            # collection rate for how much of the available feed is processed by the plant each year
            CR=common_params["CR"],
            Tracked_comps=common_params["Tracked_comps"],  # tracked components
            # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
            Prod_comp_mass=common_params["Prod_comp_mass"],
            ###################################################################################################
            ###################################################################################################
            ### Superstructure formulation parameters
            numStages=common_params["numStages"],  # number of total stages
            Options_in_stage=common_params[  # number of options in each stage
                "Options_in_stage"
            ],
            # set of options k' in stage j+1 connected to option k in stage j
            Option_outlets=common_params[  # set of options k' in stage j+1 connected to option k in stage j
                "Option_outlets"
            ],
            # dictionary of tracked component retention efficiency for each option
            Option_Eff=common_params["Option_Eff"],
            ###################################################################################################
            ###################################################################################################
            ### Operating Parameters
            # profit (k$) per metric tonnes of product in terms of tracked components
            Profit=common_params["Profit"],
            # conversion of kg REE/Fe to kg REO/Fe2O3
            REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
            # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes)
            # OPEX = a*F_in + b*y
            N_OC_var=common_params["N_OC_var"],
            # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
            num_workers=common_params["num_workers"],
            labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
            # yearly operating costs per unit (k$/unit/yr)
            YCU=common_params["YCU"],
            # cost per disassembly stage unit for each disassembly option (k$/unit)
            CU=common_params["CU"],
            # disassembly rate (thousands of EOL HDDs/yr/unit)
            Dis_Rate=common_params["Dis_Rate"],
            ###################################################################################################
            ###################################################################################################
            ### Costing Parameters
            LF=common_params["LF"],  # Lang Factor
            TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
            ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
            i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
            i_CAP_esc=common_params[  # capex escalation rate (default of 3.6%)
                "i_CAP_esc"
            ],
            f_exp=common_params[  # capital expenditure schedule (default of 10%, 60%, 30%)
                "f_exp"
            ],
            # Define Python Dictionary with discretized cost by flows for each option.
            Discretized_CAPEX=common_params["Discretized_CAPEX"],
            ###################################################################################################
            ###################################################################################################
            ### Choice of objective function. Options are 'NPV' or 'COR'.capitalize
            obj_func="NPV",
            ###################################################################################################
            ###################################################################################################
            ### Consideration of environmental impacts parameters
            # boolean to decide whether or not to consider environmental impacts
            consider_environ_impacts=False,
            # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
            environ_impacts=common_params["environ_impacts"],
            epsilon=common_params[  # epsilon factor for generating Pareto front
                "epsilon"
            ],
            ###################################################################################################
            ###################################################################################################
            ### Byproduct valorization
            # boolean to decide whether or not to consider the valorization of byproducts
            consider_byprod_val=True,
            # list of byproducts
            byprods=common_params["byprods"],
            # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
            byprod_vals=common_params["byprod_vals"],
            # dictionary keeping track of which tracked component produces which byproduct
            tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
            # dictionary tracking which options produce a given byproduct
            byprod_options=common_params["byprod_options"],
            # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
            byprod_options_eff=common_params["byprod_options_eff"],
            # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
            TC_to_byproduct=common_params["TC_to_byproduct"],
        )

        return m

    def test_build(self, BV_model, get_common_params):
        # start of plant production
        prod_start = get_common_params["plant_start"] + 1
        # final year plant is in production
        plant_end = (
            get_common_params["plant_start"] + get_common_params["plant_lifetime"] - 1
        )
        # plant operational period
        operational_range = RangeSet(prod_start, plant_end)

        for t in operational_range:
            assert isinstance(BV_model.plantYear[t].total_yearly_byprod, Var)
            assert isinstance(BV_model.plantYear[t].yearly_byprod_cons, Constraint)

            assert isinstance(BV_model.plantYear[t].Byprod_Profit, Var)
            assert isinstance(BV_model.plantYear[t].byprod_profit_con, Constraint)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, BV_model):
        solver.options["NumericFocus"] = 3

        results = solver.solve(BV_model)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, BV_model, get_common_params):
        ##### Define Parameters for tests
        ###################################################################################################
        ### Plant Lifetime Parameters
        plant_start = get_common_params["plant_start"]  # start of plant production
        plant_lifetime = get_common_params["plant_lifetime"]  # lifetime of plant
        # first year is construction
        prod_start = plant_start + 1
        plant_end = (  # final year plant is in production
            plant_start + plant_lifetime - 1
        )
        # total plant lifetime
        plant_life_range = RangeSet(plant_start, plant_end)
        # operational lifetime of the plant
        operational_range = RangeSet(prod_start, plant_end)
        ###################################################################################################
        ###################################################################################################
        ### Feed parameters
        # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
        Available_feed = get_common_params["Available_feed"]
        # collection rate for how much of the available feed is processed by the plant each year (fraction)
        CR = get_common_params["CR"]
        # list of tracked components
        Tracked_comps = get_common_params["Tracked_comps"]  # tracked components
        # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
        Prod_comp_mass = get_common_params["Prod_comp_mass"]
        # calculate feed entering parameter based on yearly available feedstock and collection rate
        Feed_entering = copy.deepcopy(Available_feed)
        for key in Feed_entering:
            Feed_entering[key] = Available_feed[key] * CR
        # calculate max feed that can enter the plant
        maxFeedEntering = max(  # max feed entering plant over production period
            Feed_entering.values()
        )
        maxFeedEnteringYear = max(  # year in which max feed enters plant
            Feed_entering, key=Feed_entering.get
        )
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        numStages = get_common_params["numStages"]  # no. of total stages
        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]
        # set of options k' in stage j+1 connected to option k in stage j
        Option_outlets = get_common_params["Option_outlets"]
        # dictionary of tracked component retention efficiency for each option
        Option_Eff = get_common_params["Option_Eff"]
        maxOptions = max(Options_in_stage.values())  # max options in any of the stages
        # make a list of the options in the final stage
        final_opt_list = [(numStages, j) for j in RangeSet(Options_in_stage[numStages])]
        ###################################################################################################
        ###################################################################################################
        ### Operating Parameters
        # profit (k$) per metric tonnes of product in terms of tracked components
        Profit = get_common_params["Profit"]
        # conversion of kg REE/Fe to kg REO/Fe2O3
        REE_to_REO_Conversion = get_common_params["REE_to_REO_Conversion"]
        # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes)
        # OPEX = a*F_in + b*y
        N_OC_var = get_common_params["N_OC_var"]
        # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
        num_workers = get_common_params["num_workers"]
        labor_rate = get_common_params["labor_rate"]  # yearly wage per type of labor
        # yearly operating costs per unit (k$/unit/yr)
        YCU = get_common_params["YCU"]
        # cost per disassembly stage unit for each disassembly option (k$/unit)
        CU = get_common_params["CU"]
        # disassembly rate (thousands of EOL HDDs/yr/unit)
        Dis_Rate = get_common_params["Dis_Rate"]
        # calculate max disassembly units possible for each option
        max_dis_by_option = copy.deepcopy(Dis_Rate)
        for key in max_dis_by_option.keys():
            max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
        max_dis_workers = max(max_dis_by_option.values())
        # calculate max possible workers for the process
        max_workers = max_dis_workers + numStages * math.ceil(max(num_workers.values()))
        ###################################################################################################
        ###################################################################################################
        ### Costing Parameters
        LF = get_common_params["LF"]  # Lang Factor
        TOC_factor = get_common_params["TOC_factor"]  # Overnight costs factor
        ATWACC = get_common_params["ATWACC"]  # discount rate
        i_OC_esc = get_common_params["i_OC_esc"]  # opex, revenue escalation rate
        i_CAP_esc = get_common_params["i_CAP_esc"]  # capex escalation rate
        f_exp = get_common_params["f_exp"]  # capital expenditure schedule
        # Define Python Dictionary with discretized cost by flows for each option.
        Discretized_CAPEX = get_common_params["Discretized_CAPEX"]
        ###################################################################################################
        ###################################################################################################
        ### Byproduct valorization
        # list of byproducts
        byprods = get_common_params["byprods"]
        # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
        byprod_vals = get_common_params["byprod_vals"]
        # dictionary keeping track of which tracked component produces which byproduct
        tracked_comp_for_byprod = get_common_params["tracked_comp_for_byprod"]
        # dictionary tracking which options produce a given byproduct
        byprod_options = get_common_params["byprod_options"]
        # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
        byprod_options_eff = get_common_params["byprod_options_eff"]
        # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
        TC_to_byproduct = get_common_params["TC_to_byproduct"]
        ###################################################################################################
        ###################################################################################################
        ### Other Parameters needed for tests
        # make a list of all options
        all_opts_list = []
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                all_opts_list.append((j, k))
        # list of stages that should be chosen for optimal process
        opt_stages = [(1, 2), (2, 2), (3, 6), (4, 3)]
        # create dictionary to hold yearly F_in values
        yearly_F_in_vals = {
            key1: {
                key2: {key3: None for key3 in Tracked_comps} for key2 in all_opts_list
            }
            for key1 in operational_range
        }
        # calculate F_in vals and store in dictionary
        for t in operational_range:
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        F_in_val = 0
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            F_in_val = 0
                            yearly_F_in_vals[t][(j, k)][c] = F_in_val
                        else:
                            if j == 1:
                                F_in_val = Feed_entering[t] * Prod_comp_mass[c]
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val

                            else:
                                F_in_val = (
                                    Feed_entering[t]
                                    * Prod_comp_mass[c]
                                    * math.prod(
                                        Option_Eff[opt_stages[stage]][c]
                                        for stage in RangeSet(0, j - 2)
                                    )
                                )
                                yearly_F_in_vals[t][(j, k)][c] = F_in_val
        # create a dictionary to hold the yearly amounts of each byproduct produced
        total_yearly_byprod = {
            key1: {key2: None for key2 in byprods} for key1 in operational_range
        }
        # calculate yearly amounts of byproducts produced and store in dict
        for t in operational_range:
            for byprod in byprods:
                byprod_val = byprod_vals[byprod]
                c = tracked_comp_for_byprod[byprod]

                total_yearly_byprod[t][byprod] = sum(
                    yearly_F_in_vals[t][byprod_option][c] * TC_to_byproduct[byprod]
                    for byprod_option in byprod_options[byprod]
                )
        # create a dictionary to hold the yearly profit generated from byproduct valorization
        Byprod_Profit = {key1: None for key1 in operational_range}
        # calculate yearly amount of profit generated from byproducts and store in dict
        for t in operational_range:
            Byprod_Profit[t] = sum(
                total_yearly_byprod[t][byprod] * byprod_vals[byprod]
                for byprod in byprods
            )
        ###################################################################################################

        # test yearly byproducts produced
        for t in operational_range:
            for byprod in byprods:
                assert value(
                    BV_model.plantYear[t].total_yearly_byprod[byprod]
                ) == pytest.approx(total_yearly_byprod[t][byprod], rel=1e-8)

        # test profit generated from byproduct valorization
        for t in operational_range:
            assert value(BV_model.plantYear[t].Byprod_Profit) == pytest.approx(
                Byprod_Profit[t], rel=1e-8
            )


# Check if function exits if an incorrect objective function is chosen.
class TestNoObjectiveFunction(object):
    @pytest.fixture(scope="class")
    def No_Obj_Model(self, get_common_params):
        common_params = get_common_params

        with pytest.raises(SystemExit) as excinfo:
            m = build_model(
                ###################################################################################################
                ### Plant Lifetime Parameters
                plant_start=common_params["plant_start"],  # start of plant production
                plant_lifetime=common_params["plant_lifetime"],  # lifetime of plant
                ###################################################################################################
                ###################################################################################################
                ### Feed parameters
                # Total feedstock available for recycling each year (in terms of thousands of HDDs per year)
                Available_feed=common_params["Available_feed"],
                # collection rate for how much of the available feed is processed by the plant each year
                CR=common_params["CR"],
                Tracked_comps=common_params["Tracked_comps"],  # tracked components
                # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
                Prod_comp_mass=common_params["Prod_comp_mass"],
                ###################################################################################################
                ###################################################################################################
                ### Superstructure formulation parameters
                numStages=common_params["numStages"],  # number of total stages
                Options_in_stage=common_params[  # number of options in each stage
                    "Options_in_stage"
                ],
                # set of options k' in stage j+1 connected to option k in stage j
                Option_outlets=common_params[  # set of options k' in stage j+1 connected to option k in stage j
                    "Option_outlets"
                ],
                # dictionary of tracked component retention efficiency for each option
                Option_Eff=common_params["Option_Eff"],
                ###################################################################################################
                ###################################################################################################
                ### Operating Parameters
                # profit (k$) per metric tonnes of product in terms of tracked components
                Profit=common_params["Profit"],
                # conversion of kg REE/Fe to kg REO/Fe2O3
                REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
                # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes)
                # OPEX = a*F_in + b*y
                N_OC_var=common_params["N_OC_var"],
                # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
                num_workers=common_params["num_workers"],
                labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
                # yearly operating costs per unit (k$/unit/yr)
                YCU=common_params["YCU"],
                # cost per disassembly stage unit for each disassembly option (k$/unit)
                CU=common_params["CU"],
                # disassembly rate (thousands of EOL HDDs/yr/unit)
                Dis_Rate=common_params["Dis_Rate"],
                ###################################################################################################
                ###################################################################################################
                ### Costing Parameters
                LF=common_params["LF"],  # Lang Factor
                TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
                ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
                i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
                i_CAP_esc=common_params[  # capex escalation rate (default of 3.6%)
                    "i_CAP_esc"
                ],
                f_exp=common_params[  # capital expenditure schedule (default of 10%, 60%, 30%)
                    "f_exp"
                ],
                # Define Python Dictionary with discretized cost by flows for each option.
                Discretized_CAPEX=common_params["Discretized_CAPEX"],
                ###################################################################################################
                ###################################################################################################
                ### Choice of objective function. Options are 'NPV' or 'COR'.capitalize
                obj_func="npv", # NPV passed as lower-case. should trigger error
                ###################################################################################################
                ###################################################################################################
                ### Consideration of environmental impacts parameters
                # boolean to decide whether or not to consider environmental impacts
                consider_environ_impacts=False,
                # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
                environ_impacts=common_params["environ_impacts"],
                epsilon=common_params[  # epsilon factor for generating Pareto front
                    "epsilon"
                ],
                ###################################################################################################
                ###################################################################################################
                ### Byproduct valorization
                # boolean to decide whether or not to consider the valorization of byproducts
                consider_byprod_val=False,
                # list of byproducts
                byprods=common_params["byprods"],
                # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
                byprod_vals=common_params["byprod_vals"],
                # dictionary keeping track of which tracked component produces which byproduct
                tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
                # dictionary tracking which options produce a given byproduct
                byprod_options=common_params["byprod_options"],
                # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
                byprod_options_eff=common_params["byprod_options_eff"],
                # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes tracked component)
                TC_to_byproduct=common_params["TC_to_byproduct"],
            )
        # Check the exit code and message
        assert excinfo.value.code == "Neither COR nor NPV specified as objective function."