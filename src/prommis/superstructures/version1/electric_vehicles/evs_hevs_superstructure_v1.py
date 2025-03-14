from prommis.superstructures.version1.superstructure_v1 import run_model
import pyomo.environ as pyo

Discretized_CAPEX = {
    "(2, 1)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 23865,
            "2": 116745,
            "3": 161895,
            "4": 196725,
            "5": 225750,
            "6": 251550,
        },
    },
    "(2, 2)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 132654.25,
            "2": 775881.66,
            "3": 1121298.83,
            "4": 1395914.64,
            "5": 1631618.42,
            "6": 1842766.45,
        },
    },
    "(2, 3)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 132654.25,
            "2": 775881.66,
            "3": 1121298.83,
            "4": 1395915.64,
            "5": 1631618.42,
            "6": 1842766.45,
        },
    },
    "(3, 1)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 278640,
            "2": 559215,
            "3": 737235,
            "4": 837855,
            "5": 855270,
            "6": 968790,
        },
    },
    "(3, 2)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0},
    },
    "(3, 3)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 1523490,
            "2": 3989325,
            "3": 5877240,
            "4": 8275995,
            "5": 9511170,
            "6": 10881795,
        },
    },
    "(4, 1)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 3354645,
            "2": 4658835,
            "3": 5365110,
            "4": 5976570,
            "5": 6386790,
            "6": 6940845,
        },
    },
    "(4, 2)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 2887727.13,
            "2": 5256412.76,
            "3": 6171274.73,
            "4": 7082790.87,
            "5": 8080977.5,
            "6": 9052114.7,
        },
    },
    "(4, 3)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 1660937.13,
            "2": 2458402.76,
            "3": 3168154.73,
            "4": 3712020.87,
            "5": 4234197.5,
            "6": 4628704.7,
        },
    },
    "(4, 4)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 1262604.16,
            "2": 8814757.97,
            "3": 13202885.52,
            "4": 16771844,
            "5": 19891533.53,
            "6": 22713708.63,
        },
    },
    "(4, 5)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 2800009.41,
            "2": 4629165,
            "3": 5227080,
            "4": 5879820,
            "5": 6651240,
            "6": 7416855,
        },
    },
    "(4, 6)": {
        "Flowrates": {
            "0": 0,
            "1": 8226,
            "2": 209786,
            "3": 411346,
            "4": 612906,
            "5": 814466,
            "6": 1016029,
        },
        "Costs": {
            "0": 0,
            "1": 1573219.41,
            "2": 1831155,
            "3": 2223960,
            "4": 2509050,
            "5": 2804460,
            "6": 2993445,
        },
    },
}

# flowData = Discretized_CAPEX[str((i, j))]["Flowrates"]
# print(flowData)

model = run_model(
    ###################################################################################################
    ### Plant lifetime paramters
    # define start of construction for plant
    plant_start=2024,
    # define lifetime of the plant
    plant_lifetime=15,
    ###################################################################################################
    ###################################################################################################
    ### Feed parameters
    # Total feedstock available for recycling
    feedstock_availability={
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
    # define collection rate
    CR=0.1,
    # define components to track in the superstructure.
    tracked_comps=["Nd", "Dy", "Fe"],
    # mass (kg) of each tracked component per EV/HEV.
    # composition is
    prod_comp_mass={"Nd": 0.206 * 3, "Dy": 0.103 * 3, "Fe": 0.691 * 3},
    ###################################################################################################
    ###################################################################################################
    ### Superstructure formulation parameters
    numStages=4,
    Options_in_stage={1: 2, 2: 3, 3: 3, 4: 6},
    # set of options k' in stage j+1 connected to option k in stage j
    Option_outlets={
        # level 1
        (1, 1): [1, 2, 3],
        (1, 2): [1, 2, 3],
        # level 2
        (2, 1): [1, 2],
        (2, 2): [1, 2],
        (2, 3): [3],
        # level 3
        (3, 1): [1],
        (3, 2): [2, 3, 4],
        (3, 3): [5, 6],
    },
    # dictionary of option and efficiency for each component.
    Option_Eff={
        # Level 1 yields
        (1, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
        (1, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
        # level 2 yields
        (2, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
        (2, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
        (2, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
        # level 3 yields
        (3, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
        (3, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
        (3, 3): {"Nd": 1, "Dy": 1, "Fe": 0},
        # level 4 yields
        (4, 1): {"Nd": 0.95, "Dy": 0.95, "Fe": 0},
        (4, 2): {"Nd": 0.94, "Dy": 0.98, "Fe": 0},
        (4, 3): {"Nd": 0.94, "Dy": 0.98, "Fe": 0},
        (4, 4): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
        (4, 5): {"Nd": 0.989, "Dy": 0.989, "Fe": 0},
        (4, 6): {"Nd": 0.989, "Dy": 0.989, "Fe": 0},
    },
    ###################################################################################################
    ###################################################################################################
    ### Operating Parameters
    ## disassembly stage is the first stage of the superstructure. It is cost differently because discrete units must be utilized.
    # cost per disassembly stage unit for each disassembly option
    CU={
        (1, 1): 0,
        (1, 2): 200000,
    },
    # define the yearly operating costs per unit ($/unit*yr)
    YCU={
        (1, 1): 190960,
        (1, 2): 339.2,
    },
    # define yearly disassembly rates per unit (EOL products/yr*unit)
    yearly_dis_rate={
        (1, 1): 7868,
        (1, 2): 52800,
    },
    ## for all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it.
    N_OPEX={
        (2, 1): {"a": 0.0047, "b": 0},  # done!
        (2, 2): {"a": 0.0104, "b": 0},  # done!
        (2, 3): {"a": 0.0341, "b": 0},  # done!
        (3, 1): {"a": 67.809, "b": 0},  # done!
        (3, 2): {"a": 0, "b": 0},  # done!
        (3, 3): {"a": 20.923, "b": 0},  # done!
        (4, 1): {"a": 137.57, "b": 703756},  # done!
        (4, 2): {"a": 320.7, "b": 283591},  # done!
        (4, 3): {"a": 15.893, "b": 283288},  # done
        (4, 4): {"a": 4.7714, "b": 0},  # done!
        (4, 5): {"a": 315.3, "b": 283591},  # done!
        (4, 6): {"a": 10.496, "b": 283288},  # done!
    },
    # it is assumed that the options in the final stage of the superstructure produce the final products.
    # therefore, the sales price of final product ($/kg) in terms of the tracked components for the options in the final stage are defined below:
    N_prod_price={
        (4, 1): {"Nd": 52.416, "Dy": 197.8575, "Fe": 0},
        (4, 2): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
        (4, 3): {"Nd": 52.416, "Dy": 197.8575, "Fe": 0},
        (4, 4): {"Nd": 52.416, "Dy": 197.8575, "Fe": 0},
        (4, 5): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
        (4, 6): {"Nd": 52.416, "Dy": 197.8575, "Fe": 0},
    },
    # Define Python Dictionary with discretized cost by flows for each option.
    Discretized_CAPEX=Discretized_CAPEX,
    # define interest rate
    i_rate=0.1,
)

### Run model
model.results = pyo.SolverFactory("gurobi").solve(model, tee=True)
model.display()
print("The NPV is $", pyo.value(model.obj))
import math

for t in pyo.RangeSet(2024, 2038):
    print(
        "Cash flow in year " + str(t) + " is " + str(pyo.value(model.plantYear[t].CF))
    )
