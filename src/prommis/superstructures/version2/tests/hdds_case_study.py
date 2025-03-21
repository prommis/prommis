import copy

# import math
# import sys
# import pyomo.environ as pyo

from prommis.superstructures.version2.superstructure_v2 import run_model

### cost of purchased equipment for iron valorization processes (k$).
jarosite_capex = 300000 / 1000  # considered.
iron_hydroxide_capex = 250000 / 1000  # considered.

### OPEX parameters for iron valorization
# 10 k$/metric tonnes of jarosite processed
jaro_opex_param = 0.7 * 2.893 * 10

# 5 k$/metric tonnes of iron hydroxide processed
AFDE_iron_hydrox_opex_param = 0.4846246349 * 1.914 * 5  # AFDE Process
Sel_Leach_iron_hydrox_opex_param = 0.7 * 1.914 * 5  # Selective Leaching Process


# costs in terms of metric tons of incoming flow vs. thousands of US dollars
Discretized_CAPEX = {
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
            "1": 343.229 + jarosite_capex,
            "2": 482.425 + jarosite_capex,
            "3": 618.182 + jarosite_capex,
            "4": 743.75 + jarosite_capex,
            "5": 844.443 + jarosite_capex,
            "6": 978.48 + jarosite_capex,
            "7": 1183.835 + jarosite_capex,
            "8": 1440.661 + jarosite_capex,
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
            "1": 226.79 + iron_hydroxide_capex,
            "2": 446.435 + iron_hydroxide_capex,
            "3": 713.714 + iron_hydroxide_capex,
            "4": 1270.105 + iron_hydroxide_capex,
            "5": 1541.353 + iron_hydroxide_capex,
            "6": 2920.751 + iron_hydroxide_capex,
            "7": 3652.064 + iron_hydroxide_capex,
            "8": 5323.087 + iron_hydroxide_capex,
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
            "1": 154.6847193 + iron_hydroxide_capex,
            "2": 354.009 + iron_hydroxide_capex,
            "3": 490.597 + iron_hydroxide_capex,
            "4": 562.047 + iron_hydroxide_capex,
            "5": 679.397 + iron_hydroxide_capex,
            "6": 912.244 + iron_hydroxide_capex,
            "7": 1097.498 + iron_hydroxide_capex,
            "8": 1297.052 + iron_hydroxide_capex,
        },
    },
}

# number of EOL Desktops with HDDs
EOLDesktops = {
    2025: 3880104,
    2026: 3712231,
    2027: 3285000,
    2028: 3087313,
    2029: 2788261,
    2030: 2332477,
    2031: 2205356,
    2032: 1963099,
    2033: 1747667,
    2034: 1548760,
    2035: 1359222,
    2036: 1181324,
    2037: 1015165,
    2038: 858384,
}

# number of EOL laptops with HDDs
EOLLaptops = {
    2025: 15602359,
    2026: 15650083,
    2027: 14608125,
    2028: 14564179,
    2029: 13582609,
    2030: 11522439,
    2031: 11078718,
    2032: 10028016,
    2033: 9122756,
    2034: 8238943,
    2035: 7384090,
    2036: 6553005,
    2037: 5743730,
    2038: 4954920,
}

# sum up hdds from both potential sources
HDD_input_flow = copy.deepcopy(EOLDesktops)
for key in HDD_input_flow:
    # flow in terms of thousands of EOL products
    HDD_input_flow[key] = (EOLDesktops[key] + EOLLaptops[key]) / 1000

m = run_model(
    ###################################################################################################
    ### Plant Lifetime Parameters
    plant_start=2024,  # start of plant production
    plant_lifetime=15,  # lifetime of plant
    ###################################################################################################
    ###################################################################################################
    ### Feed parameters
    # Total feedstock available for recycling each year
    Available_feed=HDD_input_flow,
    # collection rate for how much of the available feed is processed by the plant each year
    CR=0.6,
    Tracked_comps=["Nd", "Fe"],  # tracked components
    # mass of tracked component per EOL Product (metric tonnes of component / 1000 EOL product)
    Prod_comp_mass={"Nd": 7.5e-4, "Fe": 0.00175},
    ###################################################################################################
    ###################################################################################################
    ### Superstructure formulation parameters
    numStages=4,  # number of total stages
    Options_in_stage={1: 4, 2: 4, 3: 6, 4: 4},  # number of options in each stage
    # set of options k' in stage j+1 connected to option k in stage j
    Option_outlets={
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
    Option_Eff={
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
    # profit per kg of product in terms of tracked components
    Profit={
        (4, 1): {"Nd": 69.888, "Fe": 0},
        (4, 2): {"Nd": 69.888, "Fe": 0},
        (4, 3): {"Nd": 69.888, "Fe": 0},
        (4, 4): {"Nd": 69.888, "Fe": 0},
    },
    # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it (metric tonnes).
    # OPEX (k$) = a*F_in + b*y
    N_OC_var={
        # level 2
        (2, 1): {"a": 0.0053, "b": 7929.7 / 1000},
        (2, 2): {"a": 0.0015, "b": 2233.16 / 1000},
        (2, 3): {"a": 0.0034, "b": 0},
        (2, 4): {"a": 0.0117, "b": 0},
        # level 3
        (3, 1): {"a": 15.594, "b": 4e6 / 1000},
        (3, 2): {"a": 15.594 + jaro_opex_param, "b": 4e6 / 1000},
        (3, 3): {"a": 1.8359, "b": 0},
        (3, 4): {"a": 3.7414, "b": 3476.7 / 1000},
        (3, 5): {"a": 3.7414 + Sel_Leach_iron_hydrox_opex_param, "b": 3476.7 / 1000},
        (3, 6): {"a": 1.58, "b": 0},
        # level 4
        (4, 1): {"a": 0.4997, "b": 898320 / 1000},
        (4, 2): {"a": 9.8352, "b": 677751 / 1000},
        (4, 3): {"a": 2.17, "b": 0},
        (4, 4): {"a": 2.17 + AFDE_iron_hydrox_opex_param, "b": 0},
    },
    # number of workers, and type, needed by option (for disassembly stage, its operators per unit)
    num_workers={
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
    labor_rate=8000 * 38.20 / 1000,  # yearly wage per type of labor (k$/operator/yr)
    # yearly operating costs per unit (k$/unit/yr)
    YCU={
        (1, 1): 0,
        (1, 2): 280 / 1000,
        (1, 3): 0,
        (1, 4): 10367.608 / 1000,
    },
    # cost per disassembly stage unit for each disassembly option (k$/unit)
    CU={
        (1, 1): 0,
        (1, 2): 200000 / 1000,
        (1, 3): 0,
        (1, 4): 50000 / 1000,
    },
    # disassembly rate (thousands of EOL HDDs/yr/unit)
    # for each disassembly option (in terms of EOL products disassembled per year per unit)
    Dis_Rate={
        (1, 1): 181132 / 1000,
        (1, 2): 523636 / 1000,
        (1, 3): 960000 / 1000,
        (1, 4): 2.16e7 / 1000,
    },
    ###################################################################################################
    ###################################################################################################
    ### Costing Parameters
    LF=2.97,  # Lang Factor
    TOC_factor=1.177,  # Overnight costs factor
    ATWACC=0.0577,  # discount rate. (default of 5.77%)
    i_OC_esc=0.03,  # opex, revenue (default of 3%)
    i_CAP_esc=0.036,  # capex escalation rate (default of 3.6%)
    f_exp=[0.1, 0.6, 0.3],  # capital expenditure schedule (default of 10%, 60%, 30%)
    # Define Python Dictionary with discretized cost by flows for each option.
    Discretized_CAPEX=Discretized_CAPEX,
    ###################################################################################################
    ###################################################################################################
    # Choice of objective function. Options are 'NPV' or 'COR'.capitalize
    obj_func="NPV",
    # conversion of metric tonnes REE/Fe to metric tonnes REO/Fe2O3
    REE_to_REO_Conversion={"Nd": 1.664, "Fe": 1.43},
    ###################################################################################################
    ###################################################################################################
    ### Consideration of environmental impacts parameters
    # boolean to decide whether or not to consider environmental impacts
    consider_environ_impacts=True,
    # environmental impacts matrix (kg CO2e per metric tonne of incoming flowrate)
    # Environmental Impacts Matrix
    environ_impacts={
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
    epsilon=1e16,  # epsilon factor for generating Pareto front
    ###################################################################################################
    ###################################################################################################
    ### Byproduct valorization
    # boolean to decide whether or not to consider the valorization of byproducts
    consider_byprod_val=True,
    # list of byproducts
    byprods=["Jarosite", "Iron oxide", "Residue", "Iron hydroxide"],
    # dictionary of values for each byproduct (k$/metric tonnes). Negative value indicates it cost money to dispose of the byproduct
    byprod_vals={
        "Jarosite": -0.17,
        "Iron oxide": 10,
        "Residue": -0.17,
        "Iron hydroxide": -0.17,
    },
    # dictionary keeping track of which tracked component produces which byproduct
    tracked_comp_for_byprod={
        "Jarosite": "Fe",
        "Iron oxide": "Fe",
        "Residue": "Fe",
        "Iron hydroxide": "Fe",
    },
    # dictionary tracking which options produce a given byproduct
    byprod_options={
        "Jarosite": [(3, 1)],
        "Iron oxide": [(3, 2), (3, 5), (3, 6), (4, 4)],
        "Residue": [(3, 3)],
        "Iron hydroxide": [(3, 4), (4, 3)],
    },
    # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
    byprod_options_eff={
        "Jarosite": {(3, 1): 1},
        "Iron oxide": {(3, 2): 1, (3, 5): 1, (3, 6): 1, (4, 4): 1},
        "Residue": {(3, 3): 1},
        "Iron hydroxide": {(3, 4): 0.597, (4, 3): 1},  # means 40.3% of Fe remains
    },
    # Conversion factors of tracked component to byproduct (metric tonnes byproduct / metric tonnes iron)
    Fe_to_byproduct={
        "Jarosite": 2.893,
        "Iron oxide": 1.430,
        "Residue": 1,
        "Iron hydroxide": 1.914,
    },
)

m.obj.display()
m.binOpt.display()
m.GWP.display()
