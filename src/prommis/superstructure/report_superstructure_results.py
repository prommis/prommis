#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Report Superstructure Results Code
==================================

Author: Chris Laliwala
"""
import math
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
from pyomo.environ import units as pyunits
from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice


def report_superstructure_results_overview(m, results=None):
    """
    Args:
        m: Pyomo model
        results: Solver results object
    """
    ## Check solver status from results object
    if results is not None:
        if not (
            (results.solver.status == SolverStatus.ok)
            and (results.solver.termination_condition == TerminationCondition.optimal)
        ):
            print(f"Model was not solved optimally.")
            print(f"Status: {results.solver.status}")
            print(f"Termination: {results.solver.termination_condition}")
            return None
    else:
        print("No solver results provided. Cannot verify if model was solved.")
        return None

    ## Define variables necessary for reporting results.
    # List of options that were chosen by the superstructure.
    chosen_opts = []
    # Final value of NPV.
    NPV_value = 0
    # Final value of COR.
    COR_value = 0
    # Total impacts generated over the lifetime of the plant
    total_impacts = 0

    print("\n\nSuperstructure Results Overview:")
    ## Report chosen process
    # try to get value of first variable first
    # List of options that were chosen.
    for opt in m.fs.all_opts_set:
        if hasattr(m.fs, "option_binary_var") and opt in m.fs.option_binary_var:
            try:
                if pyo.value(m.fs.option_binary_var[opt]) == 1:
                    chosen_opts.append(opt)
            except ValueError:
                # Model hasn't been solved yet
                pass

    # Create the formatted string representation of the chosen options.
    chosen_opts_str = " -> ".join(str(opt) for opt in chosen_opts)
    print(f"\nChosen Pathway: {chosen_opts_str}")

    ## Report NPV if it was chosen as the objective function.
    if (
        pyo.value(m.fs.objective_function_choice)
        == ObjectiveFunctionChoice.NET_PRESENT_VALUE.value
    ):
        NPV_value = pyo.value(m.fs.costing.net_present_value)
        print(f"\nNPV: {NPV_value:,.2f} USD")

    ## Report COR if it was chosen.
    else:
        COR_value = pyo.value(m.fs.costing.cost_of_recovery)
        print(f"\nCost of Recovery: {COR_value:,.2f} USD/kg")

    # Report total environmental impacts if it was chosen to track them.
    if hasattr(m.fs, "environmental_impacts"):
        total_impacts = pyo.value(m.fs.environmental_impacts.total_impacts)
        print(f"\nTotal Impacts: {total_impacts:,.2f}")


def report_superstructure_costing(m, results=None):
    """
    Args:
        m: Pyomo model
        results: Solver results object
    """
    ## Check solver status from results object
    if results is not None:
        if not (
            (results.solver.status == SolverStatus.ok)
            and (results.solver.termination_condition == TerminationCondition.optimal)
        ):
            print(f"Model was not solved optimally.")
            print(f"Status: {results.solver.status}")
            print(f"Termination: {results.solver.termination_condition}")
            return None
    else:
        print("No solver results provided. Cannot verify if model was solved.")
        return None

    ### Define variables necessary for reporting results.
    ## Variables for reporting objective function value
    # Final value of NPV.
    NPV_value = 0
    # Final value of COR.
    COR_value = 0
    # interest rate
    discount_factor = 0
    ## Variables for reporting cash flows
    # units for cash flow
    cash_flow_units_str = None
    # cash flow
    cash_flow = 0
    # sales and opex escalation rate
    i_esc_sales_opex = 0
    # capital expenses escalation rate
    i_esc_capex = 0
    # interest rate units
    interest_rate_units = 0
    ## Variables for reporting revenue
    # revenue from plant main product
    main_prod_rev = 0
    # revenue from plant byproducts
    byprod_rev = 0
    # total revenue
    total_rev = 0
    # units for revenue
    rev_units_str = None
    ## Variables for reporting capital expenses
    # cost of equipment for option
    equip_cost = 0
    # total equipment cost
    total_equipment_cost = 0
    # units for equipment cost
    equip_cost_units_str = None
    # total plant cost
    TPC_value = 0
    # lang factor
    LF_value = 0
    # total overnight cost
    TOC_value = 0
    # financing factor
    financing_factor_value = 0
    # other costs factor
    other_costs_factor_value = 0
    # total overnight costs expended
    TOC_expended_value = 0
    # units for total overnight costs expended
    TOC_expended_units_str = None
    ## Variables for reporting operating expenses
    # number of operators
    num_operators = 0
    # units for yearly operating expenses
    yearly_opex_units_str = None
    # yearly variable operating expenses
    yearly_variable_opex = 0
    # yearly fixed operating expenses
    yearly_fixed_opex = 0
    # yearly operating expenses
    yearly_total_opex = 0
    # cost of labor
    cost_of_labor = 0
    # Maintenance & Supply Materials (M&SM)
    m_and_sm = 0
    # Sample Analysis & Quality Assurance/Quality Control (SA&QA/QC)
    sa_and_qa_qc = 0
    # Sales, Intellectual Property, and Research & Development (S,IP,R&D)
    s_ip_r_and_d = 0
    # Administrative & Supporting Labor (A&SL)
    a_and_sl = 0
    # Fringe Benefits (FB)
    fb = 0
    # Property Taxes & Insurance (PT&I)
    pt_and_i = 0

    print("\n\nSuperstructure Costing:")
    ### Report Objective Function Value
    ## Report NPV if it was chosen as the objective function.
    if (
        pyo.value(m.fs.objective_function_choice)
        == ObjectiveFunctionChoice.NET_PRESENT_VALUE.value
    ):
        NPV_value = pyo.value(m.fs.costing.net_present_value)
        print(f"\nNPV: {NPV_value:,.2f} USD")

    ## Report COR if it was chosen.
    else:
        COR_value = pyo.value(m.fs.costing.cost_of_recovery)
        print(f"\nCost of Recovery: {COR_value:,.2f} USD/kg")

    discount_factor = pyo.value(m.fs.costing.discount_factor)
    print(f"\n\tdiscount factor: {discount_factor:.2%}")

    ### Report Cash Flows
    print("\nCash Flows:")
    print("")
    print(f"    {'Year':<6} : {'Cash Flow':<15} : {'Units':<8}")
    for t in m.fs.plant_life_range:
        if cash_flow_units_str is None:  # Capture units from first non-zero entry
            cash_flow_units_str = str(m.fs.costing.cash_flow[t].get_units())

        cash_flow = pyo.value(m.fs.costing.cash_flow[t])

        print(f"    {str(t):<6} : {cash_flow:<15,.2f} {cash_flow_units_str:<8}")

    print(f"\n\tCash Flow Interest Rates:")
    print("")
    print(f"    {'Key':<22} : {'Value':<6} : {'Units':<8}")
    i_esc_sales_opex = pyo.value(m.fs.costing.i_operating_expense_escalation)
    interest_rate_units = str(m.fs.costing.i_operating_expense_escalation.get_units())
    print(
        f"    {'OPEX escalation rate':<22} : {i_esc_sales_opex:<6.2%} : {interest_rate_units:<8}"
    )
    i_esc_capex = pyo.value(m.fs.costing.i_capital_escalation)
    print(
        f"    {'CAPEX escalation rate':<22} : {i_esc_capex:<6.2%} : {interest_rate_units:<8}"
    )

    ### Report Revenue
    print("\nRevenue:")
    print("")
    print(
        f"    {'Year':<6} : {'Main Product Revenue':<20} : {'Byproduct Revenue':<20} : {'Total Revenue':<20} : {'Units':<8}"
    )

    for t in m.fs.operational_range:
        if rev_units_str is None:  # Capture units from first non-zero entry
            rev_units_str = str(m.fs.costing.total_profit[t].get_units())

        main_prod_rev = pyo.value(m.fs.costing.main_product_profit[t])
        total_rev = pyo.value(m.fs.costing.total_profit[t])

        # Check if user selected byproduct valorization. Otherwise, set to 0.
        if hasattr(m.fs.costing, "total_byproduct_profit"):
            byprod_rev = pyo.value(m.fs.costing.total_byproduct_profit[t])
        else:
            byprod_rev = 0

        print(
            f"    {str(t):<6} : {main_prod_rev:<20,.2f} : {byprod_rev:<20,.2f} : {total_rev:<20,.2f} : {rev_units_str:<8}"
        )

    ### Report Capital Expenses
    print("\nCapital Expenses:")

    print("\tEquipment Cost Breakdown:")
    print("")
    print(f"    {'Opt':<10} : {'Equip. Cost':<12} : {'Units':<8}")

    for opt in m.fs.all_opts_set:
        equip_cost = pyo.value(m.fs.costing.equipment_cost[opt])

        if not math.isclose(equip_cost, 0, abs_tol=1e-6):
            if equip_cost_units_str is None:  # Capture units from first non-zero entry
                equip_cost_units_str = str(m.fs.costing.equipment_cost[opt].get_units())

            total_equipment_cost += equip_cost  # Add to total
            print(
                f"    {str(opt):<10} : {equip_cost:>12,.2f} : {equip_cost_units_str:<8}"
            )

    # Print total at bottom
    print(f"    {'-'*10} : {'-'*12} : {'-'*8}")
    print(
        f"    {'TOTAL':<10} : {total_equipment_cost:>12,.2f} : {equip_cost_units_str:<8}"
    )

    TPC_value = pyo.value(m.fs.costing.total_plant_cost)
    print(f"\n\tTotal Plant Cost: {TPC_value:,.2f} USD ")
    LF_value = pyo.value(m.fs.costing.lang_factor)
    print(f"\n\t\tLang Factor: {LF_value:,.2f}")

    TOC_value = pyo.value(m.fs.costing.total_overnight_cost)
    print(f"\n\tTotal Overnight Cost: {TOC_value:,.2f} USD")
    financing_factor_value = pyo.value(m.fs.costing.financing_factor)
    print(f"\n\t\tfinancing factor: {financing_factor_value:.2%}")
    other_costs_factor_value = pyo.value(m.fs.costing.other_costs_factor)
    print(f"\n\t\tother costs factor: {other_costs_factor_value:.2%}")

    print(f"\n\tTotal Overnight Cost Expended:")
    print("")
    print(f"    {'Year':<10} : {'Total Overnight Cost Expended':<30} : {'Units':<8}")
    for t in m.fs.plant_life_range:
        if TOC_expended_units_str is None:  # Capture units from first non-zero entry
            TOC_expended_units_str = str(
                m.fs.costing.total_overnight_cost_expended[t].get_units()
            )

        TOC_expended_value = pyo.value(m.fs.costing.total_overnight_cost_expended[t])
        print(
            f"    {str(t):<10} : {TOC_expended_value:<30,.2f} : {TOC_expended_units_str:<8}"
        )

    ### Report Operating Expenses
    print(f"\nOperating Expenses:")

    num_operators = pyo.value(m.fs.costing.total_operators)
    print(f"\n\tTotal Number of Operators: {num_operators:.0f}")

    print(f"\n\tBreakdown of Yearly Operating Expenses:")
    print("")
    print(
        f"    {'Year':<10} : {'Variable Operating Expenses':<30} : {'Fixed Operating Expenses':<30} : {'Total Operating Expenses':<30} : {'Units':<8}"
    )
    for t in m.fs.operational_range:
        if yearly_opex_units_str is None:  # Capture units from first non-zero entry
            yearly_opex_units_str = str(
                m.fs.costing.total_operating_expense[t].get_units()
            )

        yearly_variable_opex = pyo.value(
            m.fs.costing.aggregate_variable_operating_cost[t]
        )
        yearly_fixed_opex = pyo.value(m.fs.costing.aggregate_fixed_operating_cost[t])
        yearly_total_opex = pyo.value(m.fs.costing.total_operating_expense[t])

        print(
            f"    {str(t):<10} : {yearly_variable_opex:<30,.2f} : {yearly_fixed_opex:<30,.2f} : {yearly_total_opex:<30,.2f} : {yearly_opex_units_str:<8}"
        )

    print(f"\n\t\tFurther Breakdown of Yearly Fixed Operating Expenses:")
    print("")
    print(
        f"    {'Year':<6} : {'Cost of Labor':<15} : {'M & SM':<15} : {'SA & QA/QC':<15} : {'S, IP, R & D':<15} : {'A & SL':<15} : {'FB':<15} : {'PT & I':<15} : {'Units':<8}"
    )
    for t in m.fs.operational_range:
        cost_of_labor = pyo.value(m.fs.costing.cost_of_labor)
        m_and_sm = pyo.value(m.fs.costing.m_and_sm)
        sa_and_qa_qc = pyo.value(m.fs.costing.sa_and_qa_qc)
        s_ip_r_and_d = pyo.value(m.fs.costing.s_ip_r_and_d[t])
        a_and_sl = pyo.value(m.fs.costing.a_and_sl)
        fb = pyo.value(m.fs.costing.fb)
        pt_and_i = pyo.value(m.fs.costing.pt_and_i)

        print(
            f"    {str(t):<6} : {cost_of_labor:<15,.2f} : {m_and_sm:<15,.2f} : {sa_and_qa_qc:<15,.2f} : {s_ip_r_and_d:<15,.2f} : {a_and_sl:<15,.2f} : {fb:<15,.2f} : {pt_and_i:<15,.2f} : {yearly_opex_units_str:<8}"
        )


def report_superstructure_streams(m, results=None):
    """
    Args:
        m: Pyomo model
        results: Solver results object
    """
    ## Check solver status from results object
    if results is not None:
        if not (
            (results.solver.status == SolverStatus.ok)
            and (results.solver.termination_condition == TerminationCondition.optimal)
        ):
            print(f"Model was not solved optimally.")
            print(f"Status: {results.solver.status}")
            print(f"Termination: {results.solver.termination_condition}")
            return None
    else:
        print("No solver results provided. Cannot verify if model was solved.")
        return None

    print("\n\nSuperstructure Streams:")
    ### Define variables necessary for reporting results
    ## Variables for reporting flow entering each stage
    # units of flow entering each stage variable
    f_units_str = None
    # flow entering each stage value
    f_value = 0
    ## Variables for reporting flow entering each option
    # units of flow entering each option
    f_in_units_str = None
    # flow entering each option
    f_in_value = 0
    ## Variables for reporting flow exiting each option
    # units for flow exiting each option
    f_out_units_str = None
    # flow exiting each option
    f_out_value = 0
    ## Variables for reporting byproducts produced
    # units for byproducts produced
    byprod_produced_units_str = None
    # byproducts produced
    byprod_produced_value = 0

    print("\nFlow Entering Stage Variables (f):")
    print("")

    # Get units from the variable
    f_units_str = str(m.fs.f[next(iter(m.fs.f))].get_units())

    # Create header
    print(
        f"    {'Year':<6} : {'Stage':<8} : {'Component':<12} : {'Flow':<15} : {'Units':<8}"
    )
    print(f"    {'-'*6} : {'-'*8} : {'-'*12} : {'-'*15} : {'-'*8}")

    for t in m.fs.operational_range:
        for stage in m.fs.f_stages:
            for comp in m.fs.tracked_comps:
                f_value = pyo.value(m.fs.f[stage, comp, t])

                # Only display non-zero flows to keep table manageable
                if not math.isclose(f_value, 0, abs_tol=1e-6):
                    print(
                        f"    {str(t):<6} : {str(stage):<8} : {str(comp):<12} : {f_value:<15,.2f} : {f_units_str:<8}"
                    )

    print("\nFlow In Option Variables (f_in):")
    print("")

    # Get units from the variable
    f_in_units_str = str(m.fs.f_in[next(iter(m.fs.f_in))].get_units())

    # Create header
    print(
        f"    {'Year':<6} : {'Option':<12} : {'Component':<12} : {'Flow In':<15} : {'Units':<8}"
    )
    print(f"    {'-'*6} : {'-'*12} : {'-'*12} : {'-'*15} : {'-'*8}")

    for t in m.fs.operational_range:
        for opt in m.fs.all_opts_set:
            for comp in m.fs.tracked_comps:
                f_in_value = pyo.value(m.fs.f_in[opt, comp, t])

                # Only display non-zero flows to keep table manageable
                if not math.isclose(f_in_value, 0, abs_tol=1e-6):
                    print(
                        f"    {str(t):<6} : {str(opt):<12} : {str(comp):<12} : {f_in_value:<15,.2f} : {f_in_units_str:<8}"
                    )

    print("\nFlow Out Option Variables (f_out):")
    print("")

    # Get units from the variable
    f_out_units_str = str(m.fs.f_out[next(iter(m.fs.f_out))].get_units())

    # Create header
    print(
        f"    {'Year':<6} : {'Option':<12} : {'Component':<12} : {'Flow Out':<15} : {'Units':<8}"
    )
    print(f"    {'-'*6} : {'-'*12} : {'-'*12} : {'-'*15} : {'-'*8}")

    for t in m.fs.operational_range:
        for opt in m.fs.all_opts_set:
            for comp in m.fs.tracked_comps:
                f_out_value = pyo.value(m.fs.f_out[opt, comp, t])

                # Only display non-zero flows to keep table manageable
                if not math.isclose(f_out_value, 0, abs_tol=1e-6):
                    print(
                        f"    {str(t):<6} : {str(opt):<12} : {str(comp):<12} : {f_out_value:<15,.2f} : {f_out_units_str:<8}"
                    )

    # only report byproduct valorization if user specified model to calculate this
    if hasattr(m.fs, "byproduct_valorization"):
        print("\nByproduct Produced Variables:")
        print("")

        # Get units from the variable
        byprod_produced_units_str = str(
            m.fs.byproduct_valorization.byproduct_produced[
                next(iter(m.fs.byproduct_valorization.byproduct_produced))
            ].get_units()
        )

        # Create header
        print(
            f"    {'Year':<6} : {'Byproduct':<15} : {'Amount Produced':<18} : {'Units':<8}"
        )
        print(f"    {'-'*6} : {'-'*15} : {'-'*18} : {'-'*8}")

        for t in m.fs.operational_range:
            for byproduct in m.fs.byproduct_valorization.byproducts_set:
                byprod_produced_value = pyo.value(
                    m.fs.byproduct_valorization.byproduct_produced[byproduct, t]
                )

                # Only display non-zero production to keep table manageable
                if not math.isclose(byprod_produced_value, 0, abs_tol=1e-6):
                    print(
                        f"    {str(t):<6} : {str(byproduct):<15} : {byprod_produced_value:<18,.2f} : {byprod_produced_units_str:<8}"
                    )


def report_superstructure_environmental_impacts(m, results=None):
    """
    Args:
        m: Pyomo model
        results: Solver results object
    """
    ## Check solver status from results object
    if results is not None:
        if not (
            (results.solver.status == SolverStatus.ok)
            and (results.solver.termination_condition == TerminationCondition.optimal)
        ):
            print(f"Model was not solved optimally.")
            print(f"Status: {results.solver.status}")
            print(f"Termination: {results.solver.termination_condition}")
            return None
    else:
        print("No solver results provided. Cannot verify if model was solved.")
        return None

    ### Define variables necessary for reporting results
    ## Variables for reporting superstructure environmental impacts
    # total impacts produced over the lifetime of the plant
    total_impacts = 0
    # Epsilon factor for generating the Pareto front
    epsilon_value = 0
    # total yearly impacts
    total_yearly_impacts = 0
    # units for total yearly impacts
    total_yearly_impacts_units_str = None
    # total yearly impacts produced by each option
    option_yearly_impacts = 0
    # units for total yearly impacts produced by each option
    option_yearly_impacts_units_str = None

    print("\n\nSuperstructure Environmental Impacts:")
    # Report total environmental impacts if it was chosen to track them.
    if hasattr(m.fs, "environmental_impacts"):
        total_impacts = pyo.value(m.fs.environmental_impacts.total_impacts)
        print(f"\nTotal Impacts: {total_impacts:,.2f}")

        epsilon_value = pyo.value(m.fs.environmental_impacts.epsilon)
        print(f"\n\tEpsilon factor: {epsilon_value:,.2f}")

        print(f"\n\tImpacts broken down by year:")
        print("")
        # Get units from the variable
        total_yearly_impacts_units_str = str(
            m.fs.environmental_impacts.total_yearly_impacts[
                next(iter(m.fs.environmental_impacts.total_yearly_impacts))
            ].get_units()
        )

        # Create header
        print(f"    {'Year':<6} : {'Total Environmental Impact':<25} : {'Units':<8}")
        print(f"    {'-'*6} : {'-'*25} : {'-'*8}")

        for t in m.fs.operational_range:
            total_yearly_impacts = pyo.value(m.fs.environmental_impacts.total_yearly_impacts[t])

            # Only display non-zero impacts to keep table manageable
            if not math.isclose(total_yearly_impacts, 0, abs_tol=1e-6):
                print(
                    f"    {str(t):<6} : {total_yearly_impacts:<25,.4f} : {total_yearly_impacts_units_str:<8}"
                )

        print(f"\n\tImpacts broken down by year and option:")
        print("")

        # Get units from the variable
        option_yearly_impacts_units_str = str(
            m.fs.environmental_impacts.option_yearly_impacts[
                next(iter(m.fs.environmental_impacts.option_yearly_impacts))
            ].get_units()
        )

        # Create header
        print(
            f"    {'Year':<6} : {'Option':<12} : {'Environmental Impact':<20} : {'Units':<8}"
        )
        print(f"    {'-'*6} : {'-'*12} : {'-'*20} : {'-'*8}")

        for t in m.fs.operational_range:
            for opt in m.fs.all_opts_set:
                option_yearly_impacts = pyo.value(
                    m.fs.environmental_impacts.option_yearly_impacts[opt, t]
                )

                # Only display non-zero impacts to keep table manageable
                if not math.isclose(option_yearly_impacts, 0, abs_tol=1e-6):
                    print(
                        f"    {str(t):<6} : {str(opt):<12} : {option_yearly_impacts:<20,.4f} : {option_yearly_impacts_units_str:<8}"
                    )

    else:
        print("\n\tUser has not specified model to track environmental impacts.")
