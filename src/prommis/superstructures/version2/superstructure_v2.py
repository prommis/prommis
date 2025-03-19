#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Superstructure Code version 2
=============================

Author: Chris Laliwala
"""


import pyomo.environ as pyo
import sys
import copy
import math


def build_model(
    ###################################################################################################
    ### Plant Lifetime Parameters
    plant_start: int,  # start of plant production
    plant_lifetime: int,  # lifetime of plant
    ###################################################################################################
    ###################################################################################################
    ### Feed parameters
    # Total feedstock available for recycling each year
    # of the form -> year: amount, ...
    Available_feed: dict,
    # collection rate for how much of the available feed is processed by the plant each year (fraction)
    CR: float,
    # list of tracked components (this is assumed to be the elemental rare earths and other elemental contaminants.
    # Refer to examples)
    Tracked_comps: list,
    # mass of tracked component per EOL Product (kg component / EOL product)
    # of the form -> component: amount, ...
    Prod_comp_mass: dict,
    ###################################################################################################
    ###################################################################################################
    ### Superstructure formulation parameters
    numStages: int,  # number of total stages
    # number of options in each stage. Of the form -> stage number: number of options
    Options_in_stage: dict,
    # set of options k' in stage j+1 connected to option k in stage j
    # of the form -> (j, k): [1, 2, 3...], ...
    Option_outlets: dict,  # set of options k' in stage j+1 connected to option k in stage j
    # dictionary of tracked component retention efficiency for each option
    # of the form -> (j, k): {"comp1": eff1, "comp2": eff2, ...}, ...
    Option_Eff: dict,
    ###################################################################################################
    ###################################################################################################
    ### Operating Parameters
    # profit per kg of product in terms of tracked components
    # of the form -> (j, k): {"comp1": price1, "comp2": price2, ...}, ...
    Profit: dict,
    # conversion of kg REE/Fe to kg REO/Fe2O3
    # of the form -> "REE": conv_factor1, ...
    # see examples for more detail
    REE_to_REO_Conversion: dict,
    # For all options excluding the disassembly stage, the OPEX costs are linearly related to the flow entering it.
    # OPEX = a*F_in + b*y
    # of the form -> (j, k): {"a": val1, "b": val2}, ...
    N_OC_var: dict,
    # number of workers needed by option (for disassembly stage, its operators per unit)
    # of the form -> (j, k): val1, ...
    num_workers: dict,
    labor_rate: float,  # yearly wage per worker
    # yearly operating costs per unit ($/unit*yr) for disassembly stage
    # of the form -> (j,k): cost1, ...
    YCU: dict,
    # cost per disassembly stage unit for each disassembly option
    # of the form -> (j,k): cost1, ...
    CU: dict,
    # disassembly rate for each disassembly option (in terms of EOL products disassembled per year per unit)
    # of the form -> (j,k): rate1, ...
    Dis_Rate: dict,
    ###################################################################################################
    ###################################################################################################
    ### Costing Parameters
    LF: float,  # Lang Factor
    TOC_factor: float,  # Overnight costs factor
    ATWACC: float,  # discount rate
    i_OC_esc: float,  # opex, revenue (default of 3%)
    i_CAP_esc: float,  # capex escalation rate (default of 3.6%)
    # capital expenditure schedule
    # user should provide a list of the breakdown of % capex expenditure over the first three years
    # of the form -> [percentage1, percentage2, percentage3]
    f_exp: list,
    # Define Python Dictionary with discretized cost by flows for each option.
    # see example for form
    Discretized_CAPEX: dict,
    ###################################################################################################
    ###################################################################################################
    ### Choice of objective function. Options are 'NPV' or 'COR'
    obj_func: str,
    ###################################################################################################
    ###################################################################################################
    ### Consideration of environmental impacts parameters
    # boolean to decide whether or not to consider environmental impacts
    consider_environ_impacts: bool,
    # environmental impacts matrix (kg CO2e per incoming flowrate)
    # of the form -> (j, k): factor1, ...
    environ_impacts: dict,
    epsilon: float,  # epsilon factor for generating Pareto front
    ###################################################################################################
    ###################################################################################################
    ### Byproduct valorization
    # boolean to decide whether or not to consider the valorization of byproducts
    consider_byprod_val: bool,
    # list of byproducts.
    # of the form -> [byprod1, byprod2, ...]
    byprods: list,
    # dictionary of values for each byproduct ($/kg). Negative value indicates it cost money to dispose of the byproduct
    # of the form -> byprod1: value1, ...
    byprod_vals: dict,
    # dictionary keeping track of which tracked component produces which byproduct
    # of the form -> byprod1: tracked_comp1, ...
    tracked_comp_for_byprod: dict,
    # dictionary tracking which options produce a given byproduct
    # of the form -> byprod1: [(j1, k1), (j2, k2), ...], ...
    byprod_options: dict,
    # dictionary tracking byproduct recovery efficiency for each option (in terms of tracked component)
    # of the form -> byprod1: {(j1, k1): eff1, ...}, ...
    byprod_options_eff: dict,
    # Conversion factors of tracked component to byproduct (kg byproduct / kg iron)
    # of the form -> byprod1: conv_factor1, ...
    Fe_to_byproduct: dict,
):
    ###################################################################################################
    ### Calculate other plant lifetime parameters
    # first year is construction
    prod_start = plant_start + 1
    plant_end = plant_start + plant_lifetime - 1  # final year plant is in production
    plant_life_range = pyo.RangeSet(plant_start, plant_end)  # lifetime of the plant
    # operational lifetime of the plant
    operational_range = pyo.RangeSet(prod_start, plant_end)
    ###################################################################################################
    ###################################################################################################
    ### Calculate other feed parameters
    # calculate feed entering parameter based on yearly available feedstock and collection rate
    Feed_entering = copy.deepcopy(Available_feed)
    for key in Feed_entering:
        Feed_entering[key] = Available_feed[key] * CR
    maxFeedEntering = max(
        Feed_entering.values()
    )  # max feed entering plant over production period
    maxFeedEnteringYear = max(
        Feed_entering, key=Feed_entering.get
    )  # year in which max feed enters plant
    ###################################################################################################
    ###################################################################################################
    ### Calculate other superstructure formulation parameters
    maxOptions = max(Options_in_stage.values())  # max options in any of the stages
    # make a list of the options in the final stage
    final_opt_list = [(numStages, j) for j in pyo.RangeSet(Options_in_stage[numStages])]
    ###################################################################################################
    ###################################################################################################
    ### Calculate other operating parameters
    # calculate max disassembly units possible for each option
    max_dis_by_option = copy.deepcopy(Dis_Rate)
    for key in max_dis_by_option.keys():
        max_dis_by_option[key] = math.ceil(maxFeedEntering / Dis_Rate[key])
    max_dis_workers = max(max_dis_by_option.values()) + 10
    # calculate max possible workers for process
    max_workers = max_dis_workers + 10
    ###################################################################################################

    m = pyo.ConcreteModel()
    b = pyo.ConcreteModel()

    if obj_func == "COR":
        m.COR = pyo.Var(domain=pyo.NonNegativeReals)

    def cash_flow_block_rule(b, t):

        # These sets will be used to define the flow vars
        b.J = pyo.RangeSet(numStages)  # number of stages
        b.K = pyo.RangeSet(maxOptions)  # max options in a stage
        b.KeyComps = pyo.Set(initialize=Tracked_comps)  # key components
        jk = []  # for declaring bin vars
        jkc = []  # for declaring flow vars
        for j in b.J:
            for k in pyo.RangeSet(Options_in_stage[j]):
                jk.append((j, k))
                for c in b.KeyComps:
                    jkc.append((j, k, c))

        # mass balance constraints and flows only present after construction is complete (1 year)
        # amount of EOL products chosen to be recycled
        b.P_entering = pyo.Var(domain=pyo.NonNegativeReals)

        # F is stream entering each stage (except last stage). indexed by j, c, t.
        # Thus,
        b.F_stages = pyo.RangeSet(1, numStages - 1)

        b.F = pyo.Var(b.F_stages * b.KeyComps, domain=pyo.NonNegativeReals)

        # F^in is stream entering each option
        # F^out is stream leaving each option
        b.FlowSet = pyo.Set(within=b.J * b.K * b.KeyComps, initialize=jkc)
        b.F_in = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)
        b.F_out = pyo.Var(b.FlowSet, domain=pyo.NonNegativeReals)

        # set of options
        b.OptSet = pyo.Set(within=b.J * b.K, initialize=jk)

        ### Mass Balances
        b.init_flow_cons = pyo.ConstraintList()  # eqn. 2
        b.inlet_flow_cons = pyo.ConstraintList()  # eqn. 1
        b.interm_flow_cons = pyo.ConstraintList()  # eqn. 3
        b.outlet_flow_cons = pyo.ConstraintList()  # eqn. 4

        b.inlet_flow_cons.add(expr=m.plantYear[t].P_entering == Feed_entering[t])

        for j in b.J:
            num_options = pyo.RangeSet(Options_in_stage[j])

            if j == 1:
                # eqn. 2
                for c in b.KeyComps:
                    b.init_flow_cons.add(
                        expr=b.P_entering * Prod_comp_mass[c]
                        == sum(b.F_in[j, k, c] for k in num_options)
                    )

            else:
                # eqn. 1
                for c in b.KeyComps:
                    b.inlet_flow_cons.add(
                        expr=b.F[j - 1, c] == sum(b.F_in[j, k, c] for k in num_options)
                    )

            if j != numStages:
                # eqn. 4
                for c in b.KeyComps:
                    b.outlet_flow_cons.add(
                        expr=sum(b.F_out[j, k, c] for k in num_options) == b.F[j, c]
                    )

        # eqn. 3
        for j in b.J:
            num_options = pyo.RangeSet(Options_in_stage[j])
            for k in num_options:
                for c in b.KeyComps:
                    a = Option_Eff[(j, k)][c]
                    b.interm_flow_cons.add(expr=b.F_in[j, k, c] * a == b.F_out[j, k, c])

        ### Calculation of yearly byproducts produced
        # only if user specifies this
        if consider_byprod_val == True:
            b.Byprods = pyo.Set(initialize=byprods)  # byproducts
            b.total_yearly_byprod = pyo.Var(
                b.Byprods
            )  # yearly amounts of a byproduct produced
            b.yearly_byprod_cons = pyo.ConstraintList()
            for byprod in b.Byprods:
                byprod_val = byprod_vals[byprod]
                c = tracked_comp_for_byprod[byprod]

                b.yearly_byprod_cons.add(
                    expr=b.total_yearly_byprod[byprod]
                    == sum(
                        b.F_in[byprod_option + (c,)]
                        * byprod_options_eff[byprod][byprod_option]
                        * Fe_to_byproduct[byprod]
                        for byprod_option in byprod_options[byprod]
                    )
                )

        ### consideration of yearly environmental impacts
        # only if user species this
        if consider_environ_impacts == True:
            b.total_yearly_GWP = pyo.Var(domain=pyo.NonNegativeReals)
            b.yearly_GWP = pyo.Var(b.OptSet, domain=pyo.NonNegativeReals)
            b.GWP_cons = pyo.ConstraintList()

            # calculate GWP associated with each technology option
            for j in b.J:
                num_options = pyo.RangeSet(Options_in_stage[j])
                for k in num_options:
                    b.GWP_cons.add(
                        expr=b.yearly_GWP[j, k]
                        == sum(b.F_in[j, k, c] for c in b.KeyComps)
                        * environ_impacts[j, k]
                    )

            # calculate the total GWP for the year
            b.GWP_cons.add(
                expr=b.total_yearly_GWP == sum(b.yearly_GWP[elem] for elem in b.OptSet)
            )

    # build mass balance constraints for each year the plant is in operation.
    m.plant_year = pyo.RangeSet(prod_start, plant_end)
    m.plantYear = pyo.Block(m.plant_year, rule=cash_flow_block_rule)

    # These sets will be used to define the binary and flow vars
    m.J = pyo.RangeSet(numStages)  # number of stages
    m.K = pyo.RangeSet(maxOptions)  # max options in a stage
    m.KeyComps = pyo.Set(initialize=Tracked_comps)  # key components
    jk_dis = []  # set of options
    jk_no_dis = []  # set of options (excluding disassembly)
    for j in m.J:
        for k in pyo.RangeSet(Options_in_stage[j]):
            jk_dis.append((j, k))
            if j != 1:
                jk_no_dis.append((j, k))
    m.OptSet = pyo.Set(within=m.J * m.K, initialize=jk_dis)  # set of options
    m.OptSet_no_dis = pyo.Set(within=m.J * m.K, initialize=jk_no_dis)  # set of options

    ### declare binary variables
    m.binOpt = pyo.Var(m.OptSet, domain=pyo.Binary)

    ### Logical Constraints
    # eqn. 5
    m.stage_bin_cons = pyo.ConstraintList()
    for j in m.J:
        num_options = pyo.RangeSet(Options_in_stage[j])
        m.stage_bin_cons.add(expr=sum(m.binOpt[j, k] for k in num_options) == 1)

    # eqn. 6
    m.connect_bin_cons = pyo.ConstraintList()
    j = 0
    for j in pyo.RangeSet(1, numStages - 1):
        num_options = pyo.RangeSet(Options_in_stage[j])

        k = 0
        for k in num_options:
            opt_connects = Option_outlets[(j, k)]
            m.connect_bin_cons.add(
                expr=1
                - m.binOpt[j, k]
                + sum(m.binOpt[j + 1, kp] for kp in opt_connects)
                >= 1
            )

    ## big-M constraints
    m.big_M_cons = pyo.ConstraintList()
    M = {"Nd": 0, "Dy": 0, "Fe": 0}
    val = 0
    for c in m.KeyComps:
        val = maxFeedEntering * Prod_comp_mass[c]
        M[c] = val * 10

    for t in pyo.RangeSet(prod_start, plant_end):
        for j in m.J:
            num_options = pyo.RangeSet(Options_in_stage[j])

            for k in num_options:
                for c in m.KeyComps:
                    m.big_M_cons.add(
                        expr=m.plantYear[t].F_in[j, k, c] <= m.binOpt[j, k] * M[c]
                    )  # eqn. 7
                    m.big_M_cons.add(
                        expr=m.plantYear[t].F_out[j, k, c] <= m.binOpt[j, k] * M[c]
                    )  # eqn. 8

    # define disassembly works set
    m.J_dis = pyo.RangeSet(1)  # number of stages
    m.K_dis = pyo.RangeSet(maxOptions)  # max options in a stage
    dis_workers_range = pyo.RangeSet(0, max_dis_workers)
    jk_dis = []
    jkw_dis = []  # for declaring bin vars
    for j_dis in m.J_dis:
        for k_dis in pyo.RangeSet(Options_in_stage[j_dis]):
            jk_dis.append((j_dis, k_dis))
            for w_dis in dis_workers_range:
                jkw_dis.append((j_dis, k_dis, w_dis))

    m.disOpts = pyo.Set(within=m.J_dis * m.K_dis, initialize=jk_dis)
    m.DisOptWorkersSet = pyo.Set(
        within=m.J_dis * m.K_dis * dis_workers_range, initialize=jkw_dis
    )
    m.DisOptWorkers = pyo.Var(m.DisOptWorkersSet, domain=pyo.Binary)

    ################################################ Cash Flow Constraints ################################################
    m.OC_var_cons = pyo.ConstraintList()
    m.DisWorkerCons = pyo.ConstraintList()  # for the disassembly stage
    # calculate profit
    m.final_opt_set = pyo.Set(initialize=final_opt_list)  # list of final node list
    m.profit_opt_cons = pyo.ConstraintList()

    m.Byprods = pyo.Set(initialize=byprods)  # set of byproducts

    # opex and profit only calculated once production starts
    for t in pyo.RangeSet(prod_start, plant_end):
        m.plantYear[t].OC_var = pyo.Var(m.OptSet, domain=pyo.Reals)

        # must be enough dis units to handle incoming feed each year.
        m.DisWorkerCons.add(
            expr=m.plantYear[t].P_entering
            <= sum(
                Dis_Rate[elem]
                * sum(i * m.DisOptWorkers[elem + (i,)] for i in dis_workers_range)
                for elem in m.disOpts
            )
        )

        for elem in m.disOpts:
            # only 1 'amount' of workers can be chosen.
            m.DisWorkerCons.add(
                expr=sum(m.DisOptWorkers[elem + (i,)] for i in dis_workers_range) == 1
            )

            # if a disassembly option is not chosen, then that 'amount' must be 0
            m.DisWorkerCons.add(expr=m.binOpt[elem] >= 1 - m.DisOptWorkers[elem + (0,)])

        # calculate variable costs for disassembly stage
        j = 1
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            option = (j, k)
            m.OC_var_cons.add(
                expr=m.plantYear[t].OC_var[option]
                == sum(m.DisOptWorkers[option + (i,)] * i for i in dis_workers_range)
                * YCU[option]
            )

        # calculate variable costs for rest of stages
        for elem in m.OptSet:
            j = 0
            k = 0
            j, k = elem
            if j > 1:  # already calculated opex for discretized options
                m.OC_var_cons.add(
                    expr=m.plantYear[t].OC_var[elem]
                    == N_OC_var[elem]["a"]
                    * sum(m.plantYear[t].F_in[elem + (c,)] for c in m.KeyComps)
                    + N_OC_var[elem]["b"] * m.binOpt[elem]
                )

        # calculate total variable cost
        m.plantYear[t].OC_var_total = pyo.Var(domain=pyo.Reals)
        m.OC_var_cons.add(
            expr=m.plantYear[t].OC_var_total
            == sum(m.plantYear[t].OC_var[elem] for elem in m.OptSet)
        )

        m.plantYear[t].ProfitOpt = pyo.Var(m.final_opt_set)

        if obj_func == "NPV":
            for opt in m.final_opt_set:
                m.profit_opt_cons.add(
                    expr=m.plantYear[t].ProfitOpt[opt]
                    == sum(
                        m.plantYear[t].F_out[opt + (c,)] * Profit[opt][c]
                        for c in m.KeyComps
                    )
                )
        elif obj_func == "COR":
            for opt in m.final_opt_set:
                m.profit_opt_cons.add(
                    expr=m.plantYear[t].ProfitOpt[opt]
                    == sum(
                        m.plantYear[t].F_out[opt + (c,)] * REE_to_REO_Conversion[c]
                        for c in m.KeyComps
                    )
                    * m.COR
                )

        else:
            sys.exit("No objective function chosen.")

        m.plantYear[t].Profit = pyo.Var(domain=pyo.NonNegativeReals)
        m.plantYear[t].profit_con = pyo.Constraint(
            expr=m.plantYear[t].Profit
            == sum(m.plantYear[t].ProfitOpt[opt] for opt in m.final_opt_set)
        )

        ## Calculate profit generated from byproduct valorization
        # only if specified by user
        if consider_byprod_val == True:
            m.plantYear[t].Byprod_Profit = pyo.Var(domain=pyo.Reals)
            m.plantYear[t].byprod_profit_con = pyo.Constraint(
                expr=m.plantYear[t].Byprod_Profit
                == sum(
                    m.plantYear[t].total_yearly_byprod[byprod] * byprod_vals[byprod]
                    for byprod in m.Byprods
                )
            )

    ################################################ OC_fixed Constraints ################################################
    # calculate the cost of labor
    m.workers = pyo.Var(m.OptSet, domain=pyo.Reals)
    m.worker_cons = pyo.ConstraintList()
    m.COL_cons = pyo.ConstraintList()

    # calculate the number of workers for each stage
    dis_workers_range = pyo.RangeSet(0, max_dis_workers)
    j = 0
    for j in pyo.RangeSet(1, numStages):
        k = 0
        elem = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            elem = (j, k)

            if j == 1:
                m.worker_cons.add(
                    expr=m.workers[j, k]
                    == sum(i * m.DisOptWorkers[j, k, i] for i in dis_workers_range)
                    * num_workers[elem]
                )

            else:
                m.worker_cons.add(
                    expr=m.workers[j, k] == m.binOpt[j, k] * num_workers[elem]
                )

    m.workers_range = pyo.RangeSet(0, max_workers)
    m.bin_workers = pyo.Var(m.workers_range, domain=pyo.Binary)
    m.total_workers = pyo.Var(domain=pyo.Reals)

    m.worker_cons.add(
        expr=sum(m.workers[elem] for elem in m.OptSet)
        <= sum(i * m.bin_workers[i] for i in m.workers_range)
    )
    m.worker_cons.add(expr=sum(m.bin_workers[i] for i in m.workers_range) == 1)

    m.worker_cons.add(
        expr=m.total_workers == sum(i * m.bin_workers[i] for i in m.workers_range)
    )

    # calculate total COL
    m.COL_Total = pyo.Var(domain=pyo.NonNegativeReals)
    m.COL_Total_con = pyo.Constraint(expr=m.COL_Total == m.total_workers * labor_rate)

    ################################################ CAPEX Constraints ################################################
    t_max = maxFeedEnteringYear

    # set of non-disassembly options

    # TPC = BEC * LF
    m.BEC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.BEC_cons = pyo.ConstraintList()

    # make a var for the max flow entering each option. Used to calculate capex
    maxFlowUB = maxFeedEntering * sum(Prod_comp_mass[c] for c in Tracked_comps)
    m.BEC_max_flow = pyo.Var(m.OptSet_no_dis, bounds=(0, maxFlowUB))
    m.BEC_max_flow_cons = pyo.ConstraintList()

    # Calculate disassembly BEC (disassembly stage is always the first stage)
    k = 0
    for k in pyo.RangeSet(Options_in_stage[1]):
        m.BEC_cons.add(
            expr=m.BEC[1, k]
            == sum(i * m.DisOptWorkers[1, k, i] for i in dis_workers_range) * CU[1, k]
        )

    # calculate BEC for rest of stages
    for elem in m.OptSet:
        j, k = elem
        if j > 1:  # already calculated capex for disassembly stage
            # calculate max flow stream for capex
            # find the max flow (flow must be >= all flows entering it)
            for t in pyo.RangeSet(prod_start, plant_end):
                m.BEC_max_flow_cons.add(
                    expr=m.BEC_max_flow[elem]
                    >= sum(m.plantYear[t].F_in[elem + (c,)] for c in m.KeyComps)
                )

    j = 0
    k = 0
    for j in pyo.RangeSet(2, numStages):
        for k in pyo.RangeSet(Options_in_stage[j]):
            # get x and y data
            flowData = list(Discretized_CAPEX[str((j, k))]["Flowrates"].values())
            CAPEXData = list(Discretized_CAPEX[str((j, k))]["Costs"].values())

            # use m.add_component to generate all piecewise functions
            # piecewise = Piecewise(yval, xval, *kwargs)
            piecewise = pyo.Piecewise(
                m.BEC[j, k],
                m.BEC_max_flow[j, k],
                pw_pts=flowData,
                pw_constr_type="EQ",
                f_rule=CAPEXData,
                pw_repn="SOS2",
            )
            optName = (j, k)
            print("Piecewise_Node" + str(optName))
            m.add_component("Piecewise_Node" + str(optName), piecewise)

    m.TPC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.TPC_cons = pyo.ConstraintList()
    # TPC for disassembly (no installation factor for disassembly)
    j = 1  # disassembly takes place in first stage always
    k = 0
    for k in pyo.RangeSet(Options_in_stage[j]):
        m.TPC_cons.add(expr=m.TPC[j, k] == m.BEC[j, k])

    # multiply all TPCs by lang factor
    j = 0
    for j in pyo.RangeSet(2, numStages):
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            m.TPC_cons.add(expr=m.TPC[j, k] == m.BEC[j, k] * LF)

    # calculate total TPC
    m.Total_TPC = pyo.Var(domain=pyo.NonNegativeReals)
    m.Total_TPC_con = pyo.Constraint(
        expr=m.Total_TPC == sum(m.TPC[elem] for elem in m.OptSet)
    )

    # calculate TOC
    m.TOC = pyo.Var(domain=pyo.NonNegativeReals)
    m.TOC_con = pyo.Constraint(expr=m.TOC == m.Total_TPC * TOC_factor)

    # calculate node TOCs
    m.node_TOC = pyo.Var(m.OptSet, domain=pyo.NonNegativeReals)
    m.node_TOC_cons = pyo.ConstraintList()

    j = 0
    for j in pyo.RangeSet(1, numStages):
        k = 0
        for k in pyo.RangeSet(Options_in_stage[j]):
            elem = (j, k)
            m.node_TOC_cons.add(
                expr=m.node_TOC[elem]
                == m.TPC[elem] + 0.027 * m.TPC[elem] + 0.15 * m.TPC[elem]
            )

    #### Construct cash flows
    m.CF = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.CF_cons = pyo.ConstraintList()
    ## CAPEX cons
    # yearly TOC expenditure
    m.TOC_exp = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.TOC_exp_cons = pyo.ConstraintList()

    ## OPEX cons
    # net earnings
    # gross earnings
    m.GE = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.GE_cons = pyo.ConstraintList()
    # revenue
    m.Rev = pyo.Var(plant_life_range, domain=pyo.Reals)
    m.Rev_cons = pyo.ConstraintList()
    # fixed operation costs
    m.OC_fixed = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.OC_fixed_cons = pyo.ConstraintList()
    # variable operating costs
    m.OC_var = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    # plant overhead costs
    m.OH = pyo.Var(plant_life_range, domain=pyo.NonNegativeReals)
    m.OH_cons = pyo.ConstraintList()

    t = 0
    for t in plant_life_range:
        # capital expenditure only in first three years
        if t <= (plant_start + 2):
            m.TOC_exp_cons.add(
                expr=m.TOC_exp[t]
                == ((1 + i_CAP_esc) ** (t - plant_start))
                * f_exp[t - plant_start]
                * m.TOC
            )
        else:
            # no capital escalation after third year
            m.TOC_exp_cons.add(expr=m.TOC_exp[t] == 0)

        # opex and revenue only once production begins (starting second year)
        if prod_start <= t:
            # fixed operating costs
            m.OC_fixed_cons.add(
                expr=m.OC_fixed[t]
                == m.COL_Total
                + 0.02 * m.Total_TPC
                + 0.1 * m.COL_Total
                + 0.01 * m.plantYear[t].Profit
                + 0.2 * m.COL_Total
                + 0.25 * m.COL_Total
                + 0.01 * m.Total_TPC
            )

            # variable costs
            m.OC_var_cons.add(expr=m.OC_var[t] == m.plantYear[t].OC_var_total)

            ## revenue
            # consider byproduct valorization if specified by user
            if consider_byprod_val == True:
                m.Rev_cons.add(
                    expr=m.Rev[t]
                    == m.plantYear[t].Profit + m.plantYear[t].Byprod_Profit
                )
            else:  # otherwise don't
                m.Rev_cons.add(expr=m.Rev[t] == m.plantYear[t].Profit)

        else:
            m.OC_fixed_cons.add(expr=m.OC_fixed[t] == 0)
            m.OC_var_cons.add(expr=m.OC_var[t] == 0)
            m.Rev_cons.add(expr=m.Rev[t] == 0)

        # calculate plant overhead
        m.OH_cons.add(expr=m.OH[t] == 0.2 * (m.OC_fixed[t] + m.OC_var[t]))

        # calculate the gross earnings
        m.GE_cons.add(
            expr=m.GE[t]
            == (m.Rev[t] - m.OC_fixed[t] - m.OC_var[t] - m.OH[t])
            * (1 + i_OC_esc) ** (t - 2023)
        )

        # calculate the cash flow
        m.CF_cons.add(expr=m.CF[t] == m.GE[t] - m.TOC_exp[t])

    ### testing
    m.bin_test_cons = pyo.ConstraintList()

    ### environmental impacts
    # only consider if option is enabled by user
    if consider_environ_impacts == True:
        m.GWP = pyo.Var(domain=pyo.NonNegativeReals)
        m.GWP_cons = pyo.ConstraintList()

        m.GWP_cons.add(
            expr=m.GWP
            == sum(
                m.plantYear[t].total_yearly_GWP
                for t in pyo.RangeSet(prod_start, plant_end)
            )
        )
        # add in epsilon constraint
        m.epsilon_con = pyo.Constraint(expr=m.GWP <= epsilon)

    if obj_func == "NPV":

        def obj_rule(m):
            return sum(
                m.CF[t] / (1 + ATWACC) ** (t - plant_start) for t in plant_life_range
            )

        m.obj = pyo.Objective(rule=obj_rule, sense=pyo.maximize)

    elif obj_func == "COR":
        m.NPV = pyo.Var(domain=pyo.Reals)
        m.NPV_con1 = pyo.Constraint(
            expr=m.NPV
            == sum(
                m.CF[t] / (1 + ATWACC) ** (t - plant_start) for t in plant_life_range
            )
        )
        m.NPV_con2 = pyo.Constraint(expr=m.NPV == 0)
        m.obj = pyo.Objective(expr=m.COR, sense=pyo.minimize)
    else:
        sys.exit("no objective function was specified")

    # solver = pyo.SolverFactory("gurobi")
    # solver.options["NumericFocus"] = 2

    # # Enable the solution pool
    # solver.options["PoolSearchMode"] = 2  # final all solutions within the gap
    # solver.options["PoolSolutions"] = 10  # store up to 10 solutions in the pool
    # solver.options["PoolGap"] = 0  # look for multiple solutions

    # m.results = pyo.SolverFactory("gurobi").solve(m, tee=True)

    return m


def solve_model(m):
    solver = pyo.SolverFactory("gurobi")
    solver.options["NumericFocus"] = 2

    res = pyo.SolverFactory("gurobi").solve(m, tee=False)

    return res
