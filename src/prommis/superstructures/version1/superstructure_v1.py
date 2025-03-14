import pyomo.environ as pyo
from math import ceil
import copy


def run_model(
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
    Discretized_CAPEX=None,
    # define interest rate
    i_rate=0.1,
):

    ### calculate other plant lifetime parameters
    # first year is construction
    prod_start = plant_start + 1

    # final year plant is in production
    plant_end = plant_start + plant_lifetime - 1

    # lifetime of the plant
    plant_life_range = pyo.RangeSet(plant_start, plant_end)

    # operational lifetime of the plant
    operational_range = pyo.RangeSet(prod_start, plant_end)

    ### Calculate other feed parameters
    # Calculate flow of EOL Products being processed by the plant given the availability of feedstock and the collection rate
    feed_entering = copy.deepcopy(feedstock_availability)
    for key in feed_entering.keys():
        feed_entering[key] = feedstock_availability[key] * CR

    # max feed entering the plant. Used to bound the flow constraints (max possible)
    maxFeed = max(feed_entering.values())

    # get the year the max flow enters the plant. Equipment will be sized to handle this flow.
    # Used to calculate the CAPEX.
    maxFeedYear = max(feed_entering, key=feed_entering.get)

    ### Calculate rest of superstructure formulation parameters
    # the maximum number of options in a stage. Used to build flowrates later
    maxOptions = max(Options_in_stage.values())

    ### Calculate rest of operating parameters
    # make a list of the options in the final stage
    final_opt_list = [(numStages, j) for j in pyo.RangeSet(Options_in_stage[numStages])]

    # calculate the number of disassembly units needed for each disassembly option
    Dis_Units = copy.deepcopy(yearly_dis_rate)
    for key in yearly_dis_rate.keys():
        Dis_Units[key] = ceil(maxFeed / yearly_dis_rate[key])

    m = pyo.ConcreteModel()

    b = pyo.ConcreteModel()

    def cash_flow_block_rule(b, t):

        # These sets will be used to define the binary and flow vars
        b.J = pyo.RangeSet(numStages)  # number of stages
        b.K = pyo.RangeSet(maxOptions)  # max options in a stage
        b.TrackedComps = pyo.Set(initialize=tracked_comps)  # key components
        jk = []  # for declaring bin vars
        jkc = []  # for declaring flow vars
        for j in b.J:
            for k in pyo.RangeSet(Options_in_stage[j]):
                jk.append((j, k))
                for c in b.TrackedComps:
                    jkc.append((j, k, c))

        ### declare binary variables
        b.OptSet = pyo.Set(within=b.J * b.K, initialize=jk)
        b.binOpt = pyo.Var(b.OptSet, domain=pyo.Binary)

        ### Logical Constraints
        # eqn. 5
        b.stage_bin_cons = pyo.ConstraintList()
        for j in b.J:
            num_options = pyo.RangeSet(Options_in_stage[j])
            b.stage_bin_cons.add(expr=sum(b.binOpt[j, k] for k in num_options) == 1)

        # eqn. 6
        b.connect_bin_cons = pyo.ConstraintList()
        for j in pyo.RangeSet(1, numStages - 1):
            num_options = pyo.RangeSet(Options_in_stage[j])

            for k in num_options:
                opt_connects = Option_outlets[(j, k)]
                b.connect_bin_cons.add(
                    expr=1
                    - b.binOpt[j, k]
                    + sum(b.binOpt[j + 1, kp] for kp in opt_connects)
                    >= 1
                )

        # mass balance constraints and flows only present after construction is complete (1 year)
        if t >= prod_start:
            ### declare stream vars
            # phi is user-defined stream entering first stage. Indexed by c, t.
            # Thus,
            key_comps_feed = (
                []
            )  # this is a list of the flowrates of the key comps in year t.
            for comp in b.TrackedComps:
                key_comps_feed.append(prod_comp_mass[comp] * feed_entering[t])

            key_comps_feed_dict = dict(
                zip(tracked_comps, key_comps_feed)
            )  # dict mapping key comp to feed in year t
            b.Phi = pyo.Param(b.TrackedComps, initialize=key_comps_feed_dict)

            # F is stream entering each stage (except last stage). indexed by j, c, t.
            # Thus,
            b.F_stages = pyo.RangeSet(1, numStages - 1)

            b.F = pyo.Var(b.F_stages * b.TrackedComps, domain=pyo.NonNegativeReals)

            # F^in is stream entering each option
            # F^out is stream leaving each option
            b.FlowSet = pyo.Set(within=b.J * b.K * b.TrackedComps, initialize=jkc)
            b.F_in = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)
            b.F_out = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)

            ### Mass Balances
            b.init_flow_cons = pyo.ConstraintList()  # eqn. 2
            b.inlet_flow_cons = pyo.ConstraintList()  # eqn. 1
            b.interm_flow_cons = pyo.ConstraintList()  # eqn. 3
            b.outlet_flow_cons = pyo.ConstraintList()  # eqn. 4
            for j in b.J:
                num_options = pyo.RangeSet(Options_in_stage[j])

                if j == 1:
                    # eqn. 2
                    for c in b.TrackedComps:
                        b.init_flow_cons.add(
                            expr=b.Phi[c] == sum(b.F_in[j, k, c] for k in num_options)
                        )

                else:
                    # eqn. 1
                    for c in b.TrackedComps:
                        b.inlet_flow_cons.add(
                            expr=b.F[j - 1, c]
                            == sum(b.F_in[j, k, c] for k in num_options)
                        )

                if j != numStages:
                    # eqn. 4
                    for c in b.TrackedComps:
                        b.outlet_flow_cons.add(
                            expr=sum(b.F_out[j, k, c] for k in num_options) == b.F[j, c]
                        )

            # eqn. 3
            for j in b.J:
                num_options = pyo.RangeSet(Options_in_stage[j])
                for k in num_options:
                    for c in b.TrackedComps:
                        a = Option_Eff[(j, k)][c]
                        b.interm_flow_cons.add(
                            expr=b.F_in[j, k, c] * a == b.F_out[j, k, c]
                        )

            ## big-M constraints
            b.big_M_cons = pyo.ConstraintList()
            M = maxFeed * sum(prod_comp_mass[c] for c in b.TrackedComps)

            for j in b.J:
                num_options = pyo.RangeSet(Options_in_stage[j])

                for k in num_options:
                    for c in b.TrackedComps:
                        b.big_M_cons.add(
                            expr=b.F_in[j, k, c] <= b.binOpt[j, k] * M
                        )  # eqn. 7
                        b.big_M_cons.add(
                            expr=b.F_out[j, k, c] <= b.binOpt[j, k] * M
                        )  # eqn. 8

        ################################################ Cash Flow Constraints ################################################
        # opex and profit only calculated once production starts
        if t >= prod_start:

            b.Opt_OPEX = pyo.Var(b.OptSet, domain=pyo.Reals)
            b.opex_cons = pyo.ConstraintList()

            # calculate disassembly OPEX (disassembly stage is always the first stage)
            for k in pyo.RangeSet(Options_in_stage[1]):
                b.opex_cons.add(
                    expr=b.Opt_OPEX[1, k]
                    == Dis_Units[1, k] * YCU[1, k] * b.binOpt[1, k]
                )

            # calculate OPEX for rest of stages
            for elem in b.OptSet:
                j, k = elem
                if j > 1:  # already calculated opex for disassembly stage
                    b.opex_cons.add(
                        expr=b.Opt_OPEX[elem]
                        == N_OPEX[elem]["a"]
                        * sum(b.F_in[elem + (c,)] for c in b.TrackedComps)
                        + N_OPEX[elem]["b"] * b.binOpt[elem]
                    )

            # calculate total opex
            b.OPEX = pyo.Var(domain=pyo.Reals)
            b.OPEX_con = pyo.Constraint(
                expr=b.OPEX == sum(b.Opt_OPEX[elem] for elem in b.OptSet)
            )

            # calculate profit
            b.final_opt_set = pyo.Set(
                initialize=final_opt_list
            )  # list of final node list
            b.ProfitOpt = pyo.Var(b.final_opt_set)
            b.profit_opt_cons = pyo.ConstraintList()
            for opt in b.final_opt_set:
                b.profit_opt_cons.add(
                    expr=b.ProfitOpt[opt]
                    == sum(
                        b.F_out[opt + (c,)] * N_prod_price[opt][c]
                        for c in b.TrackedComps
                    )
                )

            b.Profit = pyo.Var(domain=pyo.NonNegativeReals)
            b.profit_con = pyo.Constraint(
                expr=b.Profit == sum(b.ProfitOpt[opt] for opt in b.final_opt_set)
            )

            # build the cash flow
            b.CF = pyo.Var(domain=pyo.Reals)

            b.CF_con = pyo.Constraint(expr=b.CF == b.Profit - b.OPEX)

    # build logic and mass balance constraints for each year the plant is in operation.
    m.plant_year = pyo.RangeSet(plant_start, plant_end)
    m.plantYear = pyo.Block(m.plant_year, rule=cash_flow_block_rule)

    ### add constraints to ensure binaries don't change from year to year
    m.J = pyo.RangeSet(numStages)
    m.K = pyo.RangeSet(maxOptions)
    jk = []
    for j in m.J:
        numOptions = pyo.RangeSet(Options_in_stage[j])
        for k in numOptions:
            jk.append((j, k))

    m.binOptSet = pyo.Set(within=m.J * m.K, initialize=jk)

    def binary_rule(m, t, j, k):
        if t > (plant_start):
            return m.plantYear[t].binOpt[j, k] == m.plantYear[t - 1].binOpt[j, k]

    m.time_bin_cons = pyo.Constraint(
        pyo.RangeSet(prod_start, plant_end), m.binOptSet, rule=binary_rule
    )

    t_max = maxFeedYear
    t_capex = plant_start

    # CAPEX only in first year
    m.plantYear[t_capex].Opt_CAPEX = pyo.Var(
        m.plantYear[t_capex].OptSet, domain=pyo.NonNegativeReals
    )
    m.plantYear[t_capex].capex_cons = pyo.ConstraintList()

    # make a var for the max flow entering each option. Used to calculate capex
    maxFlowUB = maxFeed * sum(prod_comp_mass[c] for c in tracked_comps)
    m.plantYear[t_capex].capex_max_flow = pyo.Var(
        m.plantYear[t_capex].OptSet, bounds=(0, maxFlowUB)
    )
    m.plantYear[t_capex].capex_max_flow_cons = pyo.ConstraintList()

    # Calculate disassembly CAPEX (disassembly stage is always the first stage)
    for k in pyo.RangeSet(Options_in_stage[1]):
        m.plantYear[t_capex].capex_cons.add(
            expr=m.plantYear[t_capex].Opt_CAPEX[1, k]
            == Dis_Units[1, k] * CU[1, k] * m.plantYear[t_capex].binOpt[1, k]
        )

    # calculate CAPEX for rest of stages
    for elem in m.plantYear[t_capex].OptSet:
        j, k = elem
        if j > 1:  # already calculated capex for disassembly stage
            # calculate max flow stream for capex
            m.plantYear[t_capex].capex_max_flow_cons.add(
                expr=m.plantYear[t_capex].capex_max_flow[elem]
                == sum(
                    m.plantYear[t_max].F_in[elem + (c,)]
                    for c in m.plantYear[t_capex].TrackedComps
                )
            )

    ### Create piecewise-linear constraints for CAPEX vs maxFlow
    # from CapexParams import CAPEXParams

    for i in pyo.RangeSet(2, numStages):
        for j in pyo.RangeSet(Options_in_stage[i]):
            # get x and y data
            # flowData = 0 # x data
            # CAPEXData = 0 # y data
            flowData = list(Discretized_CAPEX[str((i, j))]["Flowrates"].values())
            CAPEXData = list(Discretized_CAPEX[str((i, j))]["Costs"].values())
            # flowData, CAPEXData = CAPEXParams((i,j))

            # use m.add_component to generate all piecewise functions
            # piecewise = Piecewise(yval, xval, *kwargs)
            piecewise = pyo.Piecewise(
                m.plantYear[t_capex].Opt_CAPEX[i, j],
                m.plantYear[t_capex].capex_max_flow[i, j],
                pw_pts=flowData,
                pw_constr_type="EQ",
                f_rule=CAPEXData,
                pw_repn="SOS2",
            )
            optName = (i, j)
            print("Piecewise_Node" + str(optName))
            m.add_component("Piecewise_Node" + str(optName), piecewise)

    # calculate total capex
    m.plantYear[t_capex].Total_CAPEX = pyo.Var(domain=pyo.NonNegativeReals)
    m.plantYear[t_capex].capex_cons.add(
        expr=m.plantYear[t_capex].Total_CAPEX
        == sum(
            m.plantYear[t_capex].Opt_CAPEX[elem] for elem in m.plantYear[t_capex].OptSet
        )
    )

    # calculate cash flow for first year
    m.plantYear[t_capex].CF = pyo.Var(domain=pyo.Reals)
    m.plantYear[t_capex].CF_con = pyo.Constraint(
        expr=m.plantYear[t_capex].CF == -m.plantYear[t_capex].Total_CAPEX
    )

    # Maximize NPV
    def obj_rule(m):
        return sum(
            m.plantYear[t].CF / (1 + i_rate) ** (t - (plant_start))
            for t in pyo.RangeSet(plant_start, plant_end)
        )

    m.obj = pyo.Objective(rule=obj_rule, sense=pyo.maximize)

    return m
