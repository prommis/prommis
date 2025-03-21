#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import copy
import math
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    RangeSet,
    value,
    SolverFactory,
    Var,
    Block,
    Constraint,
    Objective,
    # Set,
)

# import sys

# model statistics
# from idaes.core.util.model_statistics import (
#     # number_total_constraints,
#     # number_variables,
#     # number_unused_variables,
#     # degrees_of_freedom,
# )

from idaes.core.solvers import get_solver

from prommis.superstructures.version2.superstructure_v2 import build_model

# try:
#     solver = SolverFactory("test")
#     gurobi_available = True
# except KeyError:
#     gurobi_available = False

solver_available = SolverFactory("conopt").available()
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
        # Total feedstock available for recycling each year
        "Available_feed": {
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
        },
        # collection rate for how much of the available feed is processed by the plant each year
        "CR": 0.1,
        "Tracked_comps": ["Nd", "Dy", "Fe"],  # tracked components
        # mass of tracked component per EOL Product (kg component / EOL product)
        "Prod_comp_mass": {
            "Nd": 0.206 * 3,
            "Dy": 0.103 * 3,
            "Fe": 0.691 * 3,
        },
        ###################################################################################################
        ###################################################################################################
        ### Superstructure formulation parameters
        "numStages": 5,  # number of total stages
        # number of options in each stage
        "Options_in_stage": {
            1: 2,
            2: 4,
            3: 6,
            4: 4,
            5: 5,
        },
        # set of options k' in stage j+1 connected to option k in stage j
        "Option_outlets": {
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
        },
        # dictionary of tracked component retention efficiency for each option
        "Option_Eff": {
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
        },
        ###################################################################################################
        ###################################################################################################
        ### Operating Parameters
        # profit per kg of product in terms of tracked components
        "Profit": {
            (5, 1): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 2): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
            (5, 3): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 4): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 5): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
        },
        # conversion of kg REE/Fe to kg REO/Fe2O3
        "REE_to_REO_Conversion": {"Nd": 1.664, "Dy": 1.147, "Fe": 1.43},
        # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it.
        # OPEX = a*F_in + b*y
        "N_OC_var": {
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
        },
        # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
        "num_workers": {
            (1, 1): 1,
            (1, 2): 0,
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
        },
        "labor_rate": 8000 * 38.20,  # yearly wage per type of labor
        # yearly operating costs per unit ($/unit*yr)
        "YCU": {
            (1, 1): 0,
            (1, 2): 280,
        },
        # cost per disassembly stage unit for each disassembly option
        "CU": {
            (1, 1): 0,
            (1, 2): 200000,
        },
        # disassembly rate for each disassembly option (in terms of EOL products disassembled per year per unit)
        "Dis_Rate": {
            (1, 1): 7868,
            (1, 2): 52453,
        },
        ###################################################################################################
        ###################################################################################################
        ### Costing Parameters
        "LF": 2.97,  # Lang Factor
        "TOC_factor": 1.177,  # Overnight costs factor
        "ATWACC": 0.0577,  # discount rate. (default of 5.77%)
        "i_OC_esc": 0.03,  # opex, revenue (default of 3%)
        "i_CAP_esc": 0.036,  # capex escalation rate (default of 3.6%)
        "f_exp": [
            0.1,
            0.6,
            0.3,
        ],  # capital expenditure schedule (default of 10%, 60%, 30%)
        # Define Python Dictionary with discretized cost by flows for each option.
        "Discretized_CAPEX": {
            "(2, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 10130.08515,
                    "2": 31353.21173,
                    "3": 48788.84678,
                    "4": 60305.81927,
                    "5": 77063.4884,
                    "6": 117214.7546,
                    "7": 151018.0699,
                    "8": 195698.5419,
                },
            },
            "(2, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 11702.08515,
                    "2": 39023.21173,
                    "3": 62134.84678,
                    "4": 77539.81927,
                    "5": 100113.4884,
                    "6": 154792.7546,
                    "7": 201326.0699,
                    "8": 263374.5419,
                },
            },
            "(2, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 67799.0,
                    "2": 272739.0,
                    "3": 406588.0,
                    "4": 488109.0,
                    "5": 599667.0,
                    "6": 842665.0,
                    "7": 1028482.0,
                    "8": 1255541.0,
                },
            },
            "(2, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 121352.0,
                    "2": 490387.0,
                    "3": 732185.0,
                    "4": 879652.0,
                    "5": 1081651.0,
                    "6": 1522296.0,
                    "7": 1859740.0,
                    "8": 2272553.0,
                },
            },
            "(3, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 343228.652,
                    "2": 482425.4684,
                    "3": 618182.0594,
                    "4": 743750.2902,
                    "5": 844443.0443,
                    "6": 978479.5225,
                    "7": 1183834.522,
                    "8": 1440660.587,
                },
            },
            "(3, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 643228.652,
                    "2": 782425.4684,
                    "3": 918182.0594,
                    "4": 1043750.2902,
                    "5": 1144443.0443,
                    "6": 1278479.5225,
                    "7": 1483834.522,
                    "8": 1740660.587,
                },
            },
            "(3, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 423074.7216,
                    "2": 3042779.121,
                    "3": 5348359.01,
                    "4": 6921261.68,
                    "5": 9251002.61,
                    "6": 14933803.33,
                    "7": 19762044.75,
                    "8": 26151302.79,
                },
            },
            "(3, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 226790.0,
                    "2": 446435.0,
                    "3": 713714.0,
                    "4": 1270105.0,
                    "5": 1541353.0,
                    "6": 2920751.0,
                    "7": 3652064.0,
                    "8": 5323087.0,
                },
            },
            "(3, 5)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 476790.0,
                    "2": 696435.0,
                    "3": 963714.0,
                    "4": 1520105.0,
                    "5": 1791353.0,
                    "6": 3170751.0,
                    "7": 3902064.0,
                    "8": 5573087.0,
                },
            },
            "(3, 6)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 169491.0,
                    "2": 940300.0,
                    "3": 1534578.0,
                    "4": 1919653.0,
                    "5": 2469724.0,
                    "6": 3743401.0,
                    "7": 4774426.0,
                    "8": 6089420.0,
                },
            },
            "(4, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 0.0,
                    "2": 0.0,
                    "3": 0.0,
                    "4": 0.0,
                    "5": 0.0,
                    "6": 0.0,
                    "7": 0.0,
                    "8": 0.0,
                },
            },
            "(4, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 450073.0,
                    "2": 1623424.0,
                    "3": 2349169.0,
                    "4": 2778687.0,
                    "5": 3349410.0,
                    "6": 4585349.0,
                    "7": 5503177.0,
                    "8": 6590434.0,
                },
            },
            "(4, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 0.0,
                    "2": 0.0,
                    "3": 0.0,
                    "4": 0.0,
                    "5": 0.0,
                    "6": 0.0,
                    "7": 0.0,
                    "8": 0.0,
                },
            },
            "(4, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 0.0,
                    "2": 0.0,
                    "3": 0.0,
                    "4": 0.0,
                    "5": 0.0,
                    "6": 0.0,
                    "7": 0.0,
                    "8": 0.0,
                },
            },
            "(5, 1)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 197813.1853,
                    "2": 324151.478,
                    "3": 468747.756,
                    "4": 547368.2774,
                    "5": 655614.4213,
                    "6": 800184.1752,
                    "7": 974415.7068,
                    "8": 1114534.98,
                },
            },
            "(5, 2)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 222906.0,
                    "2": 354009.0,
                    "3": 490597.0,
                    "4": 562047.0,
                    "5": 679397.0,
                    "6": 912244.0,
                    "7": 1097498.0,
                    "8": 1297052.0,
                },
            },
            "(5, 3)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 222906.0,
                    "2": 354009.0,
                    "3": 490597.0,
                    "4": 562047.0,
                    "5": 679397.0,
                    "6": 912244.0,
                    "7": 1097498.0,
                    "8": 1297052.0,
                },
            },
            "(5, 4)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 154685.0,
                    "2": 858157.0,
                    "3": 1400520.0,
                    "4": 1751956.0,
                    "5": 2253973.0,
                    "6": 3416384.0,
                    "7": 4357340.0,
                    "8": 5557458.0,
                },
            },
            "(5, 5)": {
                "Flowrates": {
                    "0": 0.0,
                    "1": 36480.0,
                    "2": 634240.0,
                    "3": 1434800.0,
                    "4": 2083760.0,
                    "5": 3171200.0,
                    "6": 6342400.0,
                    "7": 9513600.0,
                    "8": 14270400.0,
                },
                "Costs": {
                    "0": 0.0,
                    "1": 404685.0,
                    "2": 1108157.0,
                    "3": 1650520.0,
                    "4": 2001956.0,
                    "5": 2503973.0,
                    "6": 3666384.0,
                    "7": 4607340.0,
                    "8": 5807458.0,
                },
            },
        },
        ###################################################################################################
        ###################################################################################################
        # environmental impacts matrix (kg CO2e per kg incoming flowrate)
        "environ_impacts": {
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
            (3, 5): 800,
            (3, 6): 1000,
            (4, 1): 0,
            (4, 2): 800,
            (4, 3): 600,
            (4, 4): 1000,
            (5, 1): 0,
            (5, 2): 800,
            (5, 3): 600,
            (5, 4): 800,
            (5, 5): 1000,
        },
        "epsilon": 1,  # epsilon factor for generating Pareto front
        ###################################################################################################
        ###################################################################################################
        ### Byproduct valorization
        # list of byproducts
        "byprods": ["Jarosite", "Iron oxide", "Residue", "Iron hydroxide"],
        # dictionary of values for each byproduct ($/kg). Negative value indicates it cost money to dispose of the byproduct
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
            "Iron oxide": [(3, 2), (3, 5), (3, 6), (5, 5)],
            "Residue": [(3, 3)],
            "Iron hydroxide": [(3, 4), (5, 4)],
        },
        # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
        "byprod_options_eff": {
            "Jarosite": {(3, 1): 1},
            "Iron oxide": {(3, 2): 1, (3, 5): 1, (3, 6): 1, (5, 5): 1},
            "Residue": {(3, 3): 1},
            "Iron hydroxide": {(3, 4): 0.597, (5, 4): 1},  # means 40.3% of Fe remains
        },
        # Conversion factors of tracked component to byproduct (kg byproduct / kg iron)
        "Fe_to_byproduct": {
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
            # Total feedstock available for recycling each year
            Available_feed=common_params["Available_feed"],
            # collection rate for how much of the available feed is processed by the plant each year
            CR=common_params["CR"],
            Tracked_comps=common_params["Tracked_comps"],  # tracked components
            # mass of tracked component per EOL Product (kg component / EOL product)
            Prod_comp_mass=common_params["Prod_comp_mass"],
            ###################################################################################################
            ###################################################################################################
            ### Superstructure formulation parameters
            numStages=common_params["numStages"],  # number of total stages
            Options_in_stage=common_params[
                "Options_in_stage"
            ],  # number of options in each stage
            # set of options k' in stage j+1 connected to option k in stage j
            Option_outlets=common_params[
                "Option_outlets"
            ],  # set of options k' in stage j+1 connected to option k in stage j
            # dictionary of tracked component retention efficiency for each option
            Option_Eff=common_params["Option_Eff"],
            ###################################################################################################
            ###################################################################################################
            ### Operating Parameters
            # profit per kg of product in terms of tracked components
            Profit=common_params["Profit"],
            # conversion of kg REE/Fe to kg REO/Fe2O3
            REE_to_REO_Conversion=common_params["REE_to_REO_Conversion"],
            # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it.
            # OPEX = a*F_in + b*y
            N_OC_var=common_params["N_OC_var"],
            # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
            num_workers=common_params["num_workers"],
            labor_rate=common_params["labor_rate"],  # yearly wage per type of labor
            # yearly operating costs per unit ($/unit*yr)
            YCU=common_params["YCU"],
            # cost per disassembly stage unit for each disassembly option
            CU=common_params["CU"],
            # disassembly rate for each disassembly option (in terms of EOL products disassembled per year per unit)
            Dis_Rate=common_params["Dis_Rate"],
            ###################################################################################################
            ###################################################################################################
            ### Costing Parameters
            LF=common_params["LF"],  # Lang Factor
            TOC_factor=common_params["TOC_factor"],  # Overnight costs factor
            ATWACC=common_params["ATWACC"],  # discount rate. (default of 5.77%)
            i_OC_esc=common_params["i_OC_esc"],  # opex, revenue (default of 3%)
            i_CAP_esc=common_params[
                "i_CAP_esc"
            ],  # capex escalation rate (default of 3.6%)
            f_exp=common_params[
                "f_exp"
            ],  # capital expenditure schedule (default of 10%, 60%, 30%)
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
            # environmental impacts matrix (kg CO2e per kg incoming flowrate)
            environ_impacts=common_params["environ_impacts"],
            epsilon=common_params[
                "epsilon"
            ],  # epsilon factor for generating Pareto front
            ###################################################################################################
            ###################################################################################################
            ### Byproduct valorization
            # boolean to decide whether or not to consider the valorization of byproducts
            consider_byprod_val=False,
            # list of byproducts
            byprods=common_params["byprods"],
            # dictionary of values for each byproduct ($/kg). Negative value indicates it cost money to dispose of the byproduct
            byprod_vals=common_params["byprod_vals"],
            # dictionary keeping track of which tracked component produces which byproduct
            tracked_comp_for_byprod=common_params["tracked_comp_for_byprod"],
            # dictionary tracking which options produce a given byproduct
            byprod_options=common_params["byprod_options"],
            # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
            byprod_options_eff=common_params["byprod_options_eff"],
            # Conversion factors of tracked component to byproduct (kg byproduct / kg iron)
            Fe_to_byproduct=common_params["Fe_to_byproduct"],
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

        # assert number_variables(NPV_model) == 2821
        # assert number_total_constraints(NPV_model) == 4020
        # assert number_unused_variables(NPV_model) == 0
        # assert degrees_of_freedom(NPV_model) == 890

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solve(self, NPV_model):
        solver.options["NumericFocus"] = 2

        results = solver.solve(NPV_model)
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    def test_solution(self, NPV_model, get_common_params):
        # start of plant production
        prod_start = get_common_params["plant_start"] + 1
        # final year plant is in production
        plant_end = (
            get_common_params["plant_start"] + get_common_params["plant_lifetime"] - 1
        )
        # plant operational period
        operational_range = RangeSet(prod_start, plant_end)

        # feed available for recycling each year
        Available_feed = get_common_params["Available_feed"]
        # collection rate of available feed processed by plant
        CR = get_common_params["CR"]

        # total number of stages in superstructure
        numStages = get_common_params["numStages"]

        # tracked components in superstructure
        Tracked_comps = get_common_params["Tracked_comps"]

        # mass of tracked component per EOL Product (kg component / EOL product)
        Prod_comp_mass = get_common_params["Prod_comp_mass"]

        # number of options in each stage
        Options_in_stage = get_common_params["Options_in_stage"]

        # dictionary of tracked component retention efficiency for each option
        Option_Eff = get_common_params["Option_Eff"]

        # list of stages that should be chosen for optimal process
        opt_stages = [(1, 2), (2, 2), (3, 6), (4, 4), (5, 4)]

        # calculate feed entering parameter based on yearly available feedstock and collection rate
        Feed_entering = copy.deepcopy(Available_feed)
        for key in Feed_entering:
            Feed_entering[key] = Available_feed[key] * CR

        # calculate max feed that can enter the plant
        maxFeedEntering = max(Feed_entering.values())

        # disassembly rate for each disassembly option (in terms of EOL products disassembled per year per unit)
        Dis_Rate = get_common_params["Dis_Rate"]

        # calculate max disassembly units possible for each option
        max_dis_by_option = copy.deepcopy(Dis_Rate)
        for key in max_dis_by_option.keys():
            max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])

        # # disassembly stage is the first stage
        # J_dis = RangeSet(1)
        # # set of options in the disassembly stage
        # K_dis = RangeSet(Options_in_stage[J_dis])
        # # set of possible disassembly workers
        # dis_workers_range = pyo.RangeSet()

        for t in operational_range:
            # test P_entering for each year
            assert pytest.approx(Feed_entering[t], abs=1e-8) == value(
                NPV_model.plantYear[t].P_entering
            )

            # test all F for each year
            for j in RangeSet(numStages - 1):
                for c in Tracked_comps:
                    assert pytest.approx(
                        Feed_entering[t]
                        * Prod_comp_mass[c]
                        * math.prod(
                            Option_Eff[opt_stages[stage]][c]
                            for stage in RangeSet(0, j - 1)
                        ),
                        abs=1e-8,
                    ) == value(NPV_model.plantYear[t].F[j, c])

            # test all F_in for each year
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            assert (
                                pytest.approx(
                                    value(NPV_model.plantYear[t].F_in[j, k, c]),
                                    abs=1e-8,
                                )
                                == 0
                            )

                        else:
                            if j == 1:
                                assert pytest.approx(
                                    Feed_entering[t] * Prod_comp_mass[c],
                                    abs=1e-8,
                                ) == value(NPV_model.plantYear[t].F_in[j, k, c])

                            else:
                                assert pytest.approx(
                                    Feed_entering[t]
                                    * Prod_comp_mass[c]
                                    * math.prod(
                                        Option_Eff[opt_stages[stage]][c]
                                        for stage in RangeSet(0, j - 2)
                                    ),
                                    abs=1e-8,
                                ) == value(NPV_model.plantYear[t].F_in[j, k, c])

            # test all F_out for each year
            for j in RangeSet(numStages):
                for k in RangeSet(Options_in_stage[j]):
                    for c in Tracked_comps:
                        # flow is zero for all stages that aren't part of the optimal process
                        if (j, k) not in opt_stages:
                            assert (
                                pytest.approx(
                                    value(NPV_model.plantYear[t].F_out[j, k, c]),
                                    abs=1e-8,
                                )
                                == 0
                            )

                        else:
                            assert pytest.approx(
                                Feed_entering[t]
                                * Prod_comp_mass[c]
                                * math.prod(
                                    Option_Eff[opt_stages[stage]][c]
                                    for stage in RangeSet(0, j - 1)
                                ),
                                abs=1e-8,
                            ) == value(NPV_model.plantYear[t].F_out[j, k, c])

        # test binary variables
        for j in RangeSet(numStages):
            for k in RangeSet(Options_in_stage[j]):
                if (j, k) in opt_stages:
                    assert pytest.approx(value(NPV_model.binOpt[j, k]), abs=1e-8) == 1
                else:
                    assert pytest.approx(value(NPV_model.binOpt[j, k]), abs=1e-8) == 0
