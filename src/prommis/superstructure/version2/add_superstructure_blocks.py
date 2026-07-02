#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Element-flow superstructure model blocks
========================================

Modular ``add_*`` builders for the bilinear stream-mixing (element-flow) superstructure,
organized in the style of the PrOMMiS fixed-yield superstructure package. Each concern
(plant lifetime, feed, superstructure formulation, operating, costing, mass balance,
environmental impacts) has its own ``add_*_params`` / ``add_*_vars`` / ``add_*_cons``
function. ``build_model`` in :mod:`superstructure_function` calls them in order.

The function is **generic**: it makes no assumptions about a particular feedstock. All
feed, topology, disassembly, price, and equipment data is supplied by the caller; only the
generic QGESS economic factors carry defaults.

All model components live on the flowsheet block ``m.fs`` (costing on ``m.fs.costing``,
environmental impacts on ``m.fs.environmental_impacts``), mirroring PrOMMiS, and use
lowercase snake_case names. Python-side data needed across builders is stashed on the
private bag ``m._build_data``.

The mathematical formulation (element flows, composition-preserving split, fixed per-element
efficiencies / split-fractions) is identical to
``new_formulation_test/build_superstructure_model_with_costing.py``; only the code
organization and names differ.

Author: Chris Laliwala
"""

import pyomo.environ as pyo

from .objective_function_enums import ObjectiveFunctionChoice


def _flatten(stage_to_items):
    """Flatten ``{stage: [items]}`` into a list of ``(stage, item)`` pairs."""
    return [(stage, item) for stage, items in stage_to_items.items() for item in items]


def _bag(m):
    """Return (creating if needed) the private Python-data bag on the model."""
    if not hasattr(m, "_build_data"):
        m._build_data = {}
    return m._build_data


# ---------------------------------------------------------------------------------------
# Objective-function choice
# ---------------------------------------------------------------------------------------
def add_objective_function_choice_param(m, objective_function_choice):
    """Create the flowsheet block and record the objective-function choice.

    Args:
        m (ConcreteModel): Pyomo model.
        objective_function_choice (ObjectiveFunctionChoice): Objective selection.
    """
    if not isinstance(objective_function_choice, ObjectiveFunctionChoice):
        raise TypeError("objective_function_choice must be an ObjectiveFunctionChoice member.")
    m.fs = pyo.Block()
    _bag(m)["objective_function_choice"] = objective_function_choice


# ---------------------------------------------------------------------------------------
# Plant-lifetime parameters
# ---------------------------------------------------------------------------------------
def add_plant_lifetime_params(m, plant_start_year, plant_end_year):
    """Add the planning-horizon range sets and the start/end-year parameters.

    The first year (``plant_start_year``) is the construction year; production happens in
    ``operating_years`` = ``plant_start_year + 1 .. plant_end_year``.

    Args:
        m (ConcreteModel): Pyomo model.
        plant_start_year (int): Year construction begins.
        plant_end_year (int): Final year of the planning horizon.
    """
    bag = _bag(m)
    bag["plant_start_year"] = plant_start_year
    bag["plant_end_year"] = plant_end_year

    m.fs.plant_start_year = pyo.Param(initialize=plant_start_year, doc="Construction start year.")
    m.fs.plant_end_year = pyo.Param(initialize=plant_end_year, doc="Final planning-horizon year.")
    m.fs.plant_years = pyo.RangeSet(plant_start_year, plant_end_year, doc="All planning-horizon years.")
    m.fs.operating_years = pyo.RangeSet(
        plant_start_year + 1, plant_end_year, doc="Years in which the plant produces."
    )


# ---------------------------------------------------------------------------------------
# Feed parameters
# ---------------------------------------------------------------------------------------
def add_feed_params(m, feed_products, tracked_elements, feed_composition, feed_availability, flow_upper_bound):
    """Add feed sets, the global flow bound, and feed composition / availability parameters.

    Args:
        m (ConcreteModel): Pyomo model.
        feed_products (list): Feed (end-of-life product) streams.
        tracked_elements (list): Tracked elements, e.g. ``["Nd", "Dy", "Fe"]``.
        feed_composition (dict): ``{feed: {element: mass_fraction}}``.
        feed_availability (dict): ``{feed: {year: amount_available}}``.
        flow_upper_bound (float): Global upper bound on flow variables.
    """
    bag = _bag(m)
    bag["flow_upper_bound"] = flow_upper_bound
    fs = m.fs

    fs.feed_products_set = pyo.Set(initialize=feed_products, doc="Feed (EOL product) streams.")
    fs.tracked_elements_set = pyo.Set(initialize=tracked_elements, doc="Tracked elements.")

    fs.feed_composition = pyo.Param(
        fs.feed_products_set, fs.tracked_elements_set,
        initialize={(f, e): feed_composition[f][e] for f in feed_products for e in tracked_elements},
        mutable=True, doc="Mass fraction of each element in each feed.",
    )
    fs.feed_availability = pyo.Param(
        fs.feed_products_set, fs.operating_years,
        initialize={(f, t): feed_availability[f][t] for f in feed_products for t in fs.operating_years},
        mutable=True, doc="Amount of each feed available for recycling each operating year.",
    )


# ---------------------------------------------------------------------------------------
# Superstructure-formulation parameters
# ---------------------------------------------------------------------------------------
def add_supe_formulation_params(
    m,
    stages,
    discrete_stages,
    stage_inlet_pools,
    stage_outlet_pools,
    produced_pools,
    final_product_pools,
    stage_efficiencies,
    stage_split_fractions,
    discrete_unit_rate,
    workers_per_discrete_unit,
    operators_per_stage,
    stage_capacity_big_m,
):
    """Add stage/pool sets, connectivity, per-element efficiency/split parameters, and sizing.

    Args:
        m (ConcreteModel): Pyomo model.
        stages (list): All stage identifiers.
        discrete_stages (list): Stages that use discrete disassembly units (the
            remaining stages are continuous "processing" stages).
        stage_inlet_pools (dict): ``{stage: [pools entering the stage]}``.
        stage_outlet_pools (dict): ``{stage: [pools produced by the stage]}``.
        produced_pools (list): Non-feed pools produced within the process (intermediates +
            products).
        final_product_pools (list): Pools that are saleable final products.
        stage_efficiencies (dict): ``{stage: {element: retention_efficiency}}``.
        stage_split_fractions (dict): ``{stage: {pool: {element: split_fraction}}}``.
        discrete_unit_rate (dict): ``{discrete_stage: EOL products processed / unit / year}``.
        workers_per_discrete_unit (dict): ``{discrete_stage: workers per unit}``.
        operators_per_stage (dict): ``{processing_stage: number of operators}``.
        stage_capacity_big_m (dict): ``{processing_stage: big-M capacity bound}``.
    """
    bag = _bag(m)
    fs = m.fs

    processing_stages = [j for j in stages if j not in discrete_stages]
    bag["discrete_unit_rate"] = dict(discrete_unit_rate)
    bag["workers_per_discrete_unit"] = dict(workers_per_discrete_unit)
    bag["operators_per_stage"] = dict(operators_per_stage)
    bag["stage_capacity_big_m"] = dict(stage_capacity_big_m)

    # --- sets ---
    fs.stages_set = pyo.Set(initialize=stages, doc="All process stages.")
    fs.discrete_stages_set = pyo.Set(initialize=discrete_stages, doc="Discrete disassembly stages.")
    fs.processing_stages_set = pyo.Set(initialize=processing_stages, doc="Continuous processing stages.")
    fs.produced_pools_set = pyo.Set(initialize=produced_pools, doc="Pools produced within the process.")
    fs.all_pools_set = pyo.Set(
        initialize=fs.feed_products_set | fs.produced_pools_set, doc="All pools (feeds + produced)."
    )
    fs.final_product_pools_set = pyo.Set(initialize=final_product_pools, doc="Saleable final-product pools.")
    fs.stage_inlet_pairs = pyo.Set(
        initialize=_flatten(stage_inlet_pools), dimen=2, doc="(stage, pool) inlet connectivity."
    )
    fs.stage_outlet_pairs = pyo.Set(
        initialize=_flatten(stage_outlet_pools), dimen=2, doc="(stage, pool) outlet connectivity."
    )

    # --- per-element efficiency / split-fraction parameters ---
    fs.stage_efficiency = pyo.Param(
        fs.stages_set, fs.tracked_elements_set,
        initialize={(j, e): stage_efficiencies[j][e] for j in stages for e in fs.tracked_elements_set},
        mutable=True, doc="Per-element retention efficiency of each stage.",
    )
    fs.split_fraction = pyo.Param(
        fs.stage_outlet_pairs, fs.tracked_elements_set,
        initialize={
            (j, p, e): stage_split_fractions[j][p][e]
            for (j, p) in fs.stage_outlet_pairs for e in fs.tracked_elements_set
        },
        mutable=True, doc="Fraction of each recovered element routed to each product pool.",
    )

    # --- disassembly / processing sizing parameters ---
    fs.discrete_unit_rate = pyo.Param(
        fs.discrete_stages_set, initialize=discrete_unit_rate, mutable=True,
        doc="EOL products disassembled per unit per year.",
    )
    fs.workers_per_discrete_unit = pyo.Param(
        fs.discrete_stages_set, initialize=workers_per_discrete_unit, mutable=True,
        doc="Workers required per disassembly unit.",
    )
    fs.operators_per_stage = pyo.Param(
        fs.processing_stages_set, initialize=operators_per_stage, mutable=True,
        doc="Operators required by each continuous processing stage.",
    )
    fs.stage_capacity_big_m = pyo.Param(
        fs.processing_stages_set, initialize=stage_capacity_big_m, mutable=True,
        doc="Big-M capacity bound used by the stage-selection constraint.",
    )


# ---------------------------------------------------------------------------------------
# Operating parameters
# ---------------------------------------------------------------------------------------
def add_operating_params(
    m,
    utilities,
    stage_utilities,
    utility_consumption,
    reagents,
    stage_reagents,
    reagent_consumption,
    byproducts,
    stage_byproducts,
    byproduct_generation,
    product_values,
    reagent_values,
    utility_values,
    intermediate_disposal_costs,
    byproduct_values,
    labor_rate,
):
    """Add utility/reagent/byproduct sets + consumption rates and all price data.

    Args:
        m (ConcreteModel): Pyomo model.
        utilities (list): Utility types.
        stage_utilities (dict): ``{stage: [utilities consumed]}``.
        utility_consumption (dict): ``{stage: {utility: consumption per unit inlet}}``.
        reagents (list): Reagent types.
        stage_reagents (dict): ``{stage: [reagents consumed]}``.
        reagent_consumption (dict): ``{stage: {reagent: consumption per unit inlet}}``.
        byproducts (list): Byproduct types.
        stage_byproducts (dict): ``{stage: [byproducts produced]}``.
        byproduct_generation (dict): ``{stage: {byproduct: generation per unit inlet}}``.
        product_values (dict): ``{pool: {element: price ($/kg), "DF": discount_factor}}``.
        reagent_values (dict): ``{reagent: price ($/kg or $/m^3)}``.
        utility_values (dict): ``{utility: price ($/MMBTU)}``.
        intermediate_disposal_costs (dict): ``{pool: disposal cost ($/kg)}``.
        byproduct_values (dict): ``{byproduct: value ($/kg; negative = disposal cost)}``.
        labor_rate (float): Yearly wage per operator ($/year).
    """
    bag = _bag(m)
    bag["product_values"] = product_values
    bag["reagent_values"] = reagent_values
    bag["utility_values"] = utility_values
    bag["intermediate_disposal_costs"] = intermediate_disposal_costs
    bag["byproduct_values"] = byproduct_values
    bag["labor_rate"] = labor_rate
    fs = m.fs

    fs.utilities_set = pyo.Set(initialize=utilities, doc="Utility types.")
    fs.reagents_set = pyo.Set(initialize=reagents, doc="Reagent types.")
    fs.byproducts_set = pyo.Set(initialize=byproducts, doc="Byproduct types.")
    fs.stage_utility_pairs = pyo.Set(initialize=_flatten(stage_utilities), dimen=2, doc="(stage, utility) pairs.")
    fs.stage_reagent_pairs = pyo.Set(initialize=_flatten(stage_reagents), dimen=2, doc="(stage, reagent) pairs.")
    fs.stage_byproduct_pairs = pyo.Set(
        initialize=_flatten(stage_byproducts), dimen=2, doc="(stage, byproduct) pairs."
    )

    fs.utility_consumption_rate = pyo.Param(
        fs.stage_utility_pairs,
        initialize={(j, u): utility_consumption[j][u] for (j, u) in fs.stage_utility_pairs},
        mutable=True, doc="Utility consumed per unit of total inlet flow to the stage.",
    )
    fs.reagent_consumption_rate = pyo.Param(
        fs.stage_reagent_pairs,
        initialize={(j, r): reagent_consumption[j][r] for (j, r) in fs.stage_reagent_pairs},
        mutable=True, doc="Reagent consumed per unit of total inlet flow to the stage.",
    )
    fs.byproduct_generation_rate = pyo.Param(
        fs.stage_byproduct_pairs,
        initialize={(j, b): byproduct_generation[j][b] for (j, b) in fs.stage_byproduct_pairs},
        mutable=True, doc="Byproduct generated per unit of total inlet flow to the stage.",
    )


# ---------------------------------------------------------------------------------------
# Costing parameters (CAPEX data + generic QGESS economic factors with defaults)
# ---------------------------------------------------------------------------------------
def add_costing_params(
    m,
    capital_cost_per_discrete_unit,
    discretized_equipment_cost,
    lang_factor,
    operating_expense_escalation,
    capital_escalation,
    discount_rate,
    total_overnight_cost_factor,
    capex_expenditure_schedule,
    min_production_fraction,
    capacity_expansion_lower_bound,
    capacity_expansion_upper_bound,
    maintenance_and_supply_factor,
    sample_analysis_and_qa_qc_factor,
    sales_ip_and_rnd_factor,
    admin_and_support_labor_factor,
    fringe_benefits_factor,
    property_tax_and_insurance_factor,
    plant_overhead_factor,
    stockpile_holding_cost,
):
    """Create the costing block and add CAPEX data + generic (mutable) economic factors.

    The economic factors follow NETL QGESS / Keim et al. and carry defaults; they are NOT
    feed-specific. ``capital_cost_per_discrete_unit`` and ``discretized_equipment_cost``
    are case-specific equipment data supplied by the caller. ``stockpile_holding_cost`` is the
    per-feed cost of holding a unit of feed in the stockpile for a year.
    """
    bag = _bag(m)
    bag["discretized_equipment_cost"] = discretized_equipment_cost
    bag["capex_expenditure_schedule"] = list(capex_expenditure_schedule)

    fs = m.fs
    fs.costing = pyo.Block()
    c = fs.costing

    c.capital_cost_per_discrete_unit = pyo.Param(
        fs.discrete_stages_set, initialize=capital_cost_per_discrete_unit, mutable=True,
        doc="Purchased cost per disassembly unit ($).",
    )

    # Holding cost charged on the year-end feed stockpile (variable, linear; 0 by default so
    # stockpiling is free unless a cost is supplied). Mutable so it can be swept for sensitivity
    # analysis without rebuilding the model.
    c.stockpile_holding_cost = pyo.Param(
        fs.feed_products_set, initialize=(stockpile_holding_cost or {}), default=0.0, mutable=True,
        doc="Cost of holding one unit of a feed in the stockpile for one year ($/feed-unit/yr).",
    )

    c.lang_factor = pyo.Param(initialize=lang_factor, mutable=True, doc="Lang factor (TPC from equipment cost).")
    c.operating_expense_escalation = pyo.Param(
        initialize=operating_expense_escalation, mutable=True, doc="Operating-expense escalation rate."
    )
    c.capital_escalation = pyo.Param(
        initialize=capital_escalation, mutable=True, doc="Capital-expense escalation rate."
    )
    c.discount_rate = pyo.Param(initialize=discount_rate, mutable=True, doc="Discount rate (ATWACC) for NPV.")
    c.total_overnight_cost_factor = pyo.Param(
        initialize=total_overnight_cost_factor, mutable=True,
        doc="TOC = factor * TPC (bundles financing + other owner's costs).",
    )
    c.min_production_fraction = pyo.Param(
        initialize=min_production_fraction, mutable=True, doc="Minimum production as a fraction of capacity."
    )
    c.capacity_expansion_lower_bound = pyo.Param(
        initialize=capacity_expansion_lower_bound, mutable=True, doc="Lower bound on capacity expansion/decommission."
    )
    c.capacity_expansion_upper_bound = pyo.Param(
        initialize=capacity_expansion_upper_bound, mutable=True, doc="Upper bound on capacity expansion/decommission."
    )
    c.maintenance_and_supply_factor = pyo.Param(
        initialize=maintenance_and_supply_factor, mutable=True, doc="Maintenance & supply materials (fraction of TPC)."
    )
    c.sample_analysis_and_qa_qc_factor = pyo.Param(
        initialize=sample_analysis_and_qa_qc_factor, mutable=True, doc="Sample analysis & QA/QC (fraction of COL)."
    )
    c.sales_ip_and_rnd_factor = pyo.Param(
        initialize=sales_ip_and_rnd_factor, mutable=True, doc="Sales, IP, R&D (fraction of revenue)."
    )
    c.admin_and_support_labor_factor = pyo.Param(
        initialize=admin_and_support_labor_factor, mutable=True, doc="Administrative & supporting labor (fraction of COL)."
    )
    c.fringe_benefits_factor = pyo.Param(
        initialize=fringe_benefits_factor, mutable=True, doc="Fringe benefits (fraction of COL)."
    )
    c.property_tax_and_insurance_factor = pyo.Param(
        initialize=property_tax_and_insurance_factor, mutable=True, doc="Property taxes & insurance (fraction of TPC)."
    )
    c.plant_overhead_factor = pyo.Param(
        initialize=plant_overhead_factor, mutable=True, doc="Plant overhead (fraction of total operating cost)."
    )
    c.labor_rate = pyo.Param(initialize=bag["labor_rate"], mutable=True, doc="Yearly wage per operator ($/year).")


# ---------------------------------------------------------------------------------------
# Mass-balance variables
# ---------------------------------------------------------------------------------------
def add_mass_balance_vars(m):
    """Add the element-flow core variables and derived lumped-mass / consumption variables."""
    bag = _bag(m)
    fub = bag["flow_upper_bound"]
    fs = m.fs

    # feed-side lumped masses
    fs.feed_purchased = pyo.Var(
        fs.feed_products_set, fs.plant_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each feed purchased each year.",
    )
    fs.feed_processed = pyo.Var(
        fs.feed_products_set, fs.plant_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each feed processed each year.",
    )
    fs.feed_stockpiled = pyo.Var(
        fs.feed_products_set, fs.plant_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each feed held in stockpile at end of year.",
    )

    fs.pool_mass_entering_stage = pyo.Var(
        fs.stage_inlet_pairs, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total mass of a pool entering a stage.",
    )
    fs.pool_mass_produced = pyo.Var(
        fs.produced_pools_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total mass of each produced pool.",
    )
    fs.pool_mass_disposed = pyo.Var(
        fs.produced_pools_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each produced pool sent to disposal.",
    )
    fs.pool_mass_processed = pyo.Var(
        fs.produced_pools_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each produced pool processed further.",
    )

    # === ELEMENT-FLOW core: *_flow[..., element, ...] = mass of that element ===
    fs.pool_element_flow = pyo.Var(
        fs.all_pools_set, fs.tracked_elements_set, fs.operating_years,
        domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each element present in each pool.",
    )
    fs.routing_fraction = pyo.Var(
        fs.stage_inlet_pairs, fs.operating_years, bounds=(0, 1),
        doc="Fraction of a pool routed to a stage (composition-preserving split).",
    )
    fs.element_entering_stage = pyo.Var(
        fs.stage_inlet_pairs, fs.tracked_elements_set, fs.operating_years,
        domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each element entering a stage from a given pool.",
    )
    fs.stage_inlet_flow = pyo.Var(
        fs.stages_set, fs.tracked_elements_set, fs.operating_years,
        domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each element at a stage inlet (mixed over all incoming pools).",
    )
    fs.stage_outlet_flow = pyo.Var(
        fs.stages_set, fs.tracked_elements_set, fs.operating_years,
        domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each element recovered at a stage outlet.",
    )
    fs.stage_product_flow = pyo.Var(
        fs.stage_outlet_pairs, fs.tracked_elements_set, fs.operating_years,
        domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Mass of each element in each product stream of a stage.",
    )

    # lumped masses kept as derived quantities (costing/capacity blocks reference them)
    fs.stage_inlet_mass = pyo.Var(
        fs.stages_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total mass entering each stage (sum over elements).",
    )
    fs.stage_outlet_mass = pyo.Var(
        fs.stages_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total mass recovered from each stage.",
    )
    fs.stage_product_mass = pyo.Var(
        fs.stage_outlet_pairs, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total mass of each product stream of a stage.",
    )

    # utilities / reagents / byproducts
    fs.total_utility_consumed = pyo.Var(
        fs.utilities_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total consumption of each utility.",
    )
    fs.stage_utility_consumed = pyo.Var(
        fs.stage_utility_pairs, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Utility consumed at each stage.",
    )
    fs.total_reagent_consumed = pyo.Var(
        fs.reagents_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total consumption of each reagent.",
    )
    fs.stage_reagent_consumed = pyo.Var(
        fs.stage_reagent_pairs, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Reagent consumed at each stage.",
    )
    fs.total_byproduct_produced = pyo.Var(
        fs.byproducts_set, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Total production of each byproduct.",
    )
    fs.stage_byproduct_produced = pyo.Var(
        fs.stage_byproduct_pairs, fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, fub),
        doc="Byproduct produced at each stage.",
    )


# ---------------------------------------------------------------------------------------
# Mass-balance constraints
# ---------------------------------------------------------------------------------------
def add_mass_balance_cons(m, purity_gated_stages, initial_state):
    """Add the element-flow mass balances and the (only) bilinear split.

    The single bilinear block is ``split_con`` (``element_entering_stage =
    routing_fraction * pool_element_flow``). The ``purity_gated_stages`` inlet gate belongs to
    the generic formulation (any feedstock) and is added here. Global-convergence cuts are not
    added here -- see :mod:`new_superstructure_function.superstructure.convergence_cuts`.

    Args:
        m (ConcreteModel): Pyomo model.
        purity_gated_stages (dict | None): ``{stage: allowed_element}`` — gated stages may
            carry only the allowed element at the inlet (others forced to zero).
        initial_state (dict | None): Inherited boundary state ``{"feed_stockpiled": {...},
            "stage_capacity": {...}, "discrete_units": {...}}`` at the start year.
    """
    bag = _bag(m)
    fs = m.fs
    t_start = bag["plant_start_year"]
    t_end = bag["plant_end_year"]
    final_products = set(fs.final_product_pools_set)

    @fs.Constraint(fs.feed_products_set, fs.operating_years, doc="Feed availability (Eq 1).")
    def feed_availability_con(fs, i, t):
        return fs.feed_purchased[i, t] <= fs.feed_availability[i, t]

    @fs.Constraint(fs.feed_products_set, fs.operating_years, doc="Processed <= purchased + stockpile (Eqs 2-3).")
    def feed_processing_limit_con(fs, i, t):
        if t == (t_start + 1):
            sp0 = initial_state["feed_stockpiled"][i] if initial_state else 0
            return fs.feed_processed[i, t] <= fs.feed_purchased[i, t] + sp0
        return fs.feed_processed[i, t] <= fs.feed_purchased[i, t] + fs.feed_stockpiled[i, t - 1]

    @fs.Constraint(fs.feed_products_set, fs.operating_years, doc="Stockpile balance (Eqs 4-5).")
    def stockpile_balance_con(fs, i, t):
        if t == (t_start + 1):
            sp0 = initial_state["feed_stockpiled"][i] if initial_state else 0
            return fs.feed_stockpiled[i, t] == sp0 + fs.feed_purchased[i, t] - fs.feed_processed[i, t]
        return fs.feed_stockpiled[i, t] == fs.feed_stockpiled[i, t - 1] + fs.feed_purchased[i, t] - fs.feed_processed[i, t]

    @fs.Constraint(fs.feed_products_set, fs.operating_years, doc="Empty stockpile at horizon end (Eq 6).")
    def final_stockpile_empty_con(fs, i, t):
        if t == t_end:
            return fs.feed_stockpiled[i, t] == 0
        return pyo.Constraint.Skip

    @fs.Constraint(fs.feed_products_set, fs.operating_years, doc="Feed routed to stages (Eq 7).")
    def feed_routing_con(fs, i, t):
        return fs.feed_processed[i, t] == sum(
            fs.pool_mass_entering_stage[j, i, t] for (j, pool) in fs.stage_inlet_pairs if i == pool
        )

    @fs.Constraint(fs.produced_pools_set, fs.operating_years, doc="Produced - disposed = processed (Eq 8).")
    def pool_processing_balance_con(fs, i, t):
        if i not in final_products:
            return fs.pool_mass_processed[i, t] == fs.pool_mass_produced[i, t] - fs.pool_mass_disposed[i, t]
        return pyo.Constraint.Skip

    for i in final_products:
        for t in fs.operating_years:
            fs.pool_mass_processed[i, t].fix(0)
            fs.pool_mass_disposed[i, t].fix(0)

    @fs.Constraint(fs.produced_pools_set, fs.operating_years, doc="Processed pool routed to stages (Eq 9).")
    def pool_routing_con(fs, i, t):
        return fs.pool_mass_processed[i, t] == sum(
            fs.pool_mass_entering_stage[j, i, t] for (j, pool) in fs.stage_inlet_pairs if i == pool
        )

    @fs.Constraint(fs.stages_set, fs.operating_years, doc="Total stage inlet mass (Eq 10).")
    def stage_inlet_mass_con(fs, j, t):
        return fs.stage_inlet_mass[j, t] == sum(
            fs.pool_mass_entering_stage[stage, i, t] for (stage, i) in fs.stage_inlet_pairs if j == stage
        )

    @fs.Constraint(fs.feed_products_set, fs.tracked_elements_set, fs.operating_years, doc="Fixed feed composition.")
    def feed_composition_con(fs, i, e, t):
        return fs.pool_element_flow[i, e, t] == fs.feed_composition[i, e] * fs.feed_processed[i, t]

    # --- COMPOSITION-PRESERVING SPLIT: the only bilinear block (Eq 11a) ---
    @fs.Constraint(fs.stage_inlet_pairs, fs.tracked_elements_set, fs.operating_years, doc="Split (Eq 11a, bilinear).")
    def split_con(fs, j, i, e, t):
        return fs.element_entering_stage[j, i, e, t] == fs.routing_fraction[j, i, t] * fs.pool_element_flow[i, e, t]

    @fs.Constraint(fs.stage_inlet_pairs, fs.operating_years, doc="Entering-mass link (Eq 11b).")
    def entering_mass_link_con(fs, j, i, t):
        return fs.pool_mass_entering_stage[j, i, t] == sum(
            fs.element_entering_stage[j, i, e, t] for e in fs.tracked_elements_set
        )

    # --- MIX at stage inlet (linear, Eq 11) ---
    @fs.Constraint(fs.stages_set, fs.tracked_elements_set, fs.operating_years, doc="Inlet mixing (Eq 11).")
    def stage_inlet_mixing_con(fs, j, e, t):
        return fs.stage_inlet_flow[j, e, t] == sum(
            fs.element_entering_stage[stage, i, e, t] for (stage, i) in fs.stage_inlet_pairs if j == stage
        )

    # PURITY GATE: gated stages may carry only their allowed element at the inlet.
    if purity_gated_stages:
        fs.purity_gate_cons = pyo.ConstraintList()
        for j, allowed in purity_gated_stages.items():
            allowed_set = {allowed} if isinstance(allowed, str) else set(allowed)
            for e in fs.tracked_elements_set:
                if e not in allowed_set:
                    for t in fs.operating_years:
                        fs.purity_gate_cons.add(expr=fs.stage_inlet_flow[j, e, t] == 0)

    # --- RECOVER (linear, Eq 12) ---
    @fs.Constraint(fs.stages_set, fs.tracked_elements_set, fs.operating_years, doc="Per-element recovery (Eq 12).")
    def recovery_con(fs, j, e, t):
        return fs.stage_outlet_flow[j, e, t] == fs.stage_inlet_flow[j, e, t] * fs.stage_efficiency[j, e]

    @fs.Constraint(fs.stages_set, fs.operating_years, doc="Outlet-mass link.")
    def outlet_mass_link_con(fs, j, t):
        return fs.stage_outlet_mass[j, t] == sum(fs.stage_outlet_flow[j, e, t] for e in fs.tracked_elements_set)

    # --- DISTRIBUTE to product streams (linear, Eq 13) ---
    @fs.Constraint(fs.stage_outlet_pairs, fs.tracked_elements_set, fs.operating_years, doc="Distribution (Eq 13).")
    def distribution_con(fs, j, i, e, t):
        return fs.stage_product_flow[j, i, e, t] == fs.stage_outlet_flow[j, e, t] * fs.split_fraction[j, i, e]

    @fs.Constraint(fs.stage_outlet_pairs, fs.operating_years, doc="Product-mass link.")
    def product_mass_link_con(fs, j, i, t):
        return fs.stage_product_mass[j, i, t] == sum(
            fs.stage_product_flow[j, i, e, t] for e in fs.tracked_elements_set
        )

    # --- RE-MERGE produced streams into pools (linear, Eqs 14-15) ---
    @fs.Constraint(fs.produced_pools_set, fs.operating_years, doc="Pool re-merge, mass (Eq 14).")
    def pool_remerge_mass_con(fs, i, t):
        return fs.pool_mass_produced[i, t] == sum(
            fs.stage_product_mass[stage, i, t] for (stage, pool) in fs.stage_outlet_pairs if i == pool
        )

    @fs.Constraint(fs.produced_pools_set, fs.tracked_elements_set, fs.operating_years, doc="Pool re-merge, element (Eq 15).")
    def pool_remerge_element_con(fs, i, e, t):
        return fs.pool_element_flow[i, e, t] == sum(
            fs.stage_product_flow[stage, i, e, t] for (stage, pool) in fs.stage_outlet_pairs if i == pool
        )

    # --- utilities / reagents / byproducts (linear in stage_inlet_mass) ---
    @fs.Constraint(fs.utilities_set, fs.operating_years, doc="Total utility consumption (Eq 16).")
    def total_utility_con(fs, u, t):
        return fs.total_utility_consumed[u, t] == sum(
            fs.stage_utility_consumed[j, u, t] for (j, util) in fs.stage_utility_pairs if u == util
        )

    @fs.Constraint(fs.stage_utility_pairs, fs.operating_years, doc="Stage utility consumption (Eq 17).")
    def stage_utility_con(fs, j, u, t):
        return fs.stage_utility_consumed[j, u, t] == fs.stage_inlet_mass[j, t] * fs.utility_consumption_rate[j, u]

    @fs.Constraint(fs.reagents_set, fs.operating_years, doc="Total reagent consumption (Eq 18).")
    def total_reagent_con(fs, r, t):
        return fs.total_reagent_consumed[r, t] == sum(
            fs.stage_reagent_consumed[j, r, t] for (j, reagent) in fs.stage_reagent_pairs if r == reagent
        )

    @fs.Constraint(fs.stage_reagent_pairs, fs.operating_years, doc="Stage reagent consumption (Eq 19).")
    def stage_reagent_con(fs, j, r, t):
        return fs.stage_reagent_consumed[j, r, t] == fs.stage_inlet_mass[j, t] * fs.reagent_consumption_rate[j, r]

    @fs.Constraint(fs.byproducts_set, fs.operating_years, doc="Total byproduct production (Eq 23).")
    def total_byproduct_con(fs, b, t):
        return fs.total_byproduct_produced[b, t] == sum(
            fs.stage_byproduct_produced[j, b, t] for (j, byp) in fs.stage_byproduct_pairs if b == byp
        )

    @fs.Constraint(fs.stage_byproduct_pairs, fs.operating_years, doc="Stage byproduct production (Eq 24).")
    def stage_byproduct_con(fs, j, b, t):
        return fs.stage_byproduct_produced[j, b, t] == fs.stage_inlet_mass[j, t] * fs.byproduct_generation_rate[j, b]


# ---------------------------------------------------------------------------------------
# Capacity / workforce variables and constraints
# ---------------------------------------------------------------------------------------
def add_capacity_workforce_vars(m):
    """Add design binaries, capacity-expansion variables, and disassembly-unit counts."""
    fs = m.fs
    fs.stage_selected = pyo.Var(fs.stages_set, domain=pyo.Binary, doc="1 if a stage is built/selected.")
    fs.stage_capacity = pyo.Var(
        fs.processing_stages_set, fs.plant_years, domain=pyo.NonNegativeReals,
        doc="Cumulative capacity of each processing stage.",
    )
    fs.capacity_expansion = pyo.Var(
        fs.processing_stages_set, fs.plant_years, domain=pyo.NonNegativeReals, bounds=(0, 1e3),
        doc="Capacity added in a year.",
    )
    fs.capacity_decommission = pyo.Var(
        fs.processing_stages_set, fs.plant_years, domain=pyo.NonNegativeReals,
        doc="Capacity removed in a year.",
    )
    fs.capacity_expansion_indicator = pyo.Var(
        fs.processing_stages_set, fs.plant_years, domain=pyo.Binary, doc="1 if expansion occurs."
    )
    fs.capacity_decommission_indicator = pyo.Var(
        fs.processing_stages_set, fs.plant_years, domain=pyo.Binary, doc="1 if decommission occurs."
    )
    fs.discrete_units = pyo.Var(
        fs.discrete_stages_set, fs.plant_years, domain=pyo.NonNegativeIntegers, bounds=(0, 100),
        doc="Number of disassembly units in service.",
    )
    fs.discrete_units_added = pyo.Var(
        fs.discrete_stages_set, fs.plant_years, domain=pyo.NonNegativeIntegers, bounds=(0, 100),
        doc="Number of disassembly units added in a year.",
    )


def add_capacity_workforce_cons(m, initial_state):
    """Add disassembly-unit, stage-selection, and capacity-expansion constraints."""
    bag = _bag(m)
    fs = m.fs
    t_start, t_end = bag["plant_start_year"], bag["plant_end_year"]

    # disassembly unit counts / throughput
    fs.discrete_unit_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        if t == t_start:
            for j in fs.discrete_stages_set:
                w0 = initial_state["discrete_units"][j] if initial_state else 0
                fs.discrete_unit_cons.add(expr=fs.discrete_units[j, t] == w0)
        else:
            for j in fs.discrete_stages_set:
                fs.discrete_unit_cons.add(
                    expr=fs.discrete_units[j, t - 1] + fs.discrete_units_added[j, t - 1]
                    == fs.discrete_units[j, t]
                )
            for j in fs.discrete_stages_set:
                fs.discrete_unit_cons.add(
                    expr=fs.stage_inlet_mass[j, t] <= fs.discrete_unit_rate[j] * fs.discrete_units[j, t]
                )
        if t == t_end:
            for j in fs.discrete_stages_set:
                fs.discrete_unit_cons.add(expr=fs.discrete_units_added[j, t] == 0)

    # Tie the design binary to discrete-stage usage so stage_selected is a faithful
    # "stage is built/used" indicator for discrete stages, mirroring its role for
    # continuous stages. (i) units -- and hence flow -- only if selected; (ii) selected
    # only if units are ever deployed. NPV-neutral: stage_selected is otherwise unused for
    # discrete stages, so pinning it does not change the optimum.
    unit_ub = 100  # upper bound on discrete_units (see add_capacity_workforce_vars)
    fs.discrete_selection_cons = pyo.ConstraintList()
    for j in fs.discrete_stages_set:
        for t in fs.plant_years:
            fs.discrete_selection_cons.add(expr=fs.discrete_units[j, t] <= unit_ub * fs.stage_selected[j])
        fs.discrete_selection_cons.add(
            expr=fs.stage_selected[j] <= sum(fs.discrete_units[j, t] for t in fs.plant_years)
        )

    # stage selection (loose big-M) for continuous stages
    fs.stage_selection_cons = pyo.ConstraintList()
    for t in fs.operating_years:
        for j in fs.processing_stages_set:
            fs.stage_selection_cons.add(
                expr=fs.stage_inlet_mass[j, t] <= fs.stage_capacity_big_m[j] * fs.stage_selected[j]
            )

    # capacity expansion / decommission / production limits
    lb = fs.costing.capacity_expansion_lower_bound
    ub = fs.costing.capacity_expansion_upper_bound
    mp = fs.costing.min_production_fraction
    fs.capacity_expansion_bound_cons = pyo.ConstraintList()
    fs.capacity_balance_cons = pyo.ConstraintList()
    fs.production_capacity_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        for j in fs.processing_stages_set:
            fs.capacity_expansion_bound_cons.add(expr=lb * fs.capacity_expansion_indicator[j, t] <= fs.capacity_expansion[j, t])
            fs.capacity_expansion_bound_cons.add(expr=fs.capacity_expansion[j, t] <= ub * fs.capacity_expansion_indicator[j, t])
            fs.capacity_expansion_bound_cons.add(expr=lb * fs.capacity_decommission_indicator[j, t] <= fs.capacity_decommission[j, t])
            fs.capacity_expansion_bound_cons.add(expr=fs.capacity_decommission[j, t] <= ub * fs.capacity_decommission_indicator[j, t])

        if t == t_start:
            for j in fs.processing_stages_set:
                q0 = initial_state["stage_capacity"][j] if initial_state else 0
                fs.capacity_balance_cons.add(fs.stage_capacity[j, t] == q0)
        else:
            for j in fs.processing_stages_set:
                fs.capacity_balance_cons.add(
                    expr=fs.stage_capacity[j, t - 1] + fs.capacity_expansion[j, t - 1] - fs.capacity_decommission[j, t - 1]
                    == fs.stage_capacity[j, t]
                )
                fs.production_capacity_cons.add(expr=mp * fs.stage_capacity[j, t] <= fs.stage_inlet_mass[j, t])
                fs.production_capacity_cons.add(expr=fs.stage_inlet_mass[j, t] <= fs.stage_capacity[j, t])


# ---------------------------------------------------------------------------------------
# Costing variables
# ---------------------------------------------------------------------------------------
def add_costing_vars(m):
    """Add revenue, OPEX, CAPEX, cash-flow, and NPV variables on ``m.fs.costing``."""
    fs = m.fs
    c = m.fs.costing
    c.revenue = pyo.Var(fs.plant_years, domain=pyo.NonNegativeReals, doc="Revenue each year.")
    c.variable_operating_cost = pyo.Var(fs.plant_years, domain=pyo.NonNegativeReals, doc="Variable OPEX each year.")
    c.fixed_operating_cost = pyo.Var(fs.plant_years, domain=pyo.NonNegativeReals, doc="Fixed OPEX each year.")
    c.cost_of_labor = pyo.Var(fs.operating_years, domain=pyo.NonNegativeReals, doc="Cost of labor each operating year.")
    c.cash_flow = pyo.Var(fs.plant_years, domain=pyo.Reals, doc="Cash flow each year.")
    c.net_present_value = pyo.Var(domain=pyo.Reals, doc="Net present value.")
    c.total_plant_cost = pyo.Var(domain=pyo.NonNegativeReals, doc="Total plant cost (first-year build).")
    c.total_overnight_cost = pyo.Var(domain=pyo.NonNegativeReals, doc="Total overnight cost (first-year build).")
    c.total_overnight_cost_expended = pyo.Var(fs.plant_years, domain=pyo.NonNegativeReals, doc="TOC expended each year.")
    c.plant_overhead = pyo.Var(fs.plant_years, domain=pyo.NonNegativeReals, doc="Plant overhead each year.")
    c.total_plant_cost_expansion = pyo.Var(fs.operating_years, domain=pyo.NonNegativeReals, doc="TPC of expansions.")
    c.total_overnight_cost_expansion = pyo.Var(fs.operating_years, domain=pyo.NonNegativeReals, doc="TOC of expansions.")
    c.bare_equipment_cost = pyo.Var(fs.stages_set, fs.plant_years, domain=pyo.NonNegativeReals, doc="Bare equipment cost.")
    c.total_operators = pyo.Var(fs.operating_years, domain=pyo.NonNegativeIntegers, bounds=(0, 100), doc="Total operators (integer).")
    c.total_operators_unrounded = pyo.Var(fs.operating_years, domain=pyo.NonNegativeReals, bounds=(0, 100), doc="Total operators before rounding up.")


# ---------------------------------------------------------------------------------------
# Capital-cost constraints
# ---------------------------------------------------------------------------------------
def add_capital_cost_cons(m):
    """Add bare-equipment cost (discrete + piecewise), TPC, and TOC constraints."""
    bag = _bag(m)
    fs = m.fs
    c = m.fs.costing
    t_start = bag["plant_start_year"]
    discretized_equipment_cost = bag["discretized_equipment_cost"]

    c.bare_equipment_cost_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        for j in fs.discrete_stages_set:
            c.bare_equipment_cost_cons.add(
                expr=c.bare_equipment_cost[j, t] == fs.discrete_units_added[j, t] * c.capital_cost_per_discrete_unit[j]
            )

    for t in fs.plant_years:
        for j in fs.processing_stages_set:
            capacity_points = list(discretized_equipment_cost[str(j)]["Capacity"].values())
            cost_points = list(discretized_equipment_cost[str(j)]["Costs"].values())
            piecewise = pyo.Piecewise(
                c.bare_equipment_cost[j, t],
                fs.capacity_expansion[j, t],
                pw_pts=capacity_points,
                pw_constr_type="EQ",
                f_rule=cost_points,
                pw_repn="SOS2",
            )
            c.add_component("piecewise_equipment_cost_" + str(j) + "_" + str(t), piecewise)

    c.total_plant_cost_cons = pyo.ConstraintList()
    c.total_overnight_cost_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        if t == t_start:
            c.total_plant_cost_cons.add(
                expr=c.total_plant_cost
                == sum(c.bare_equipment_cost[j, t] for j in fs.discrete_stages_set)
                + sum(c.bare_equipment_cost[j, t] * c.lang_factor for j in fs.processing_stages_set)
            )
            c.total_overnight_cost_cons.add(expr=c.total_overnight_cost == c.total_plant_cost * c.total_overnight_cost_factor)
        else:
            c.total_plant_cost_cons.add(
                expr=c.total_plant_cost_expansion[t]
                == sum(c.bare_equipment_cost[j, t] for j in fs.discrete_stages_set)
                + sum(c.bare_equipment_cost[j, t] * c.lang_factor for j in fs.processing_stages_set)
            )
            c.total_overnight_cost_cons.add(
                expr=c.total_overnight_cost_expansion[t] == c.total_plant_cost_expansion[t] * c.total_overnight_cost_factor
            )


# ---------------------------------------------------------------------------------------
# Operating-cost constraints (labor, revenue, fixed/variable OPEX, overhead)
# ---------------------------------------------------------------------------------------
def add_operating_cost_cons(m):
    """Add operator-count, cost-of-labor, revenue, fixed/variable OPEX, and overhead constraints.

    Revenue uses the element-flow form: any element in any product pool may carry a price,
    so byproducts (e.g. an Fe-oxide pool) need no separate set or block.
    """
    bag = _bag(m)
    fs = m.fs
    c = m.fs.costing
    t_start = bag["plant_start_year"]
    product_values = bag["product_values"]
    byproduct_values = bag["byproduct_values"]
    reagent_values = bag["reagent_values"]
    utility_values = bag["utility_values"]
    intermediate_disposal_costs = bag["intermediate_disposal_costs"]
    schedule = bag["capex_expenditure_schedule"]

    # total operators (rounded up to integer)
    c.operator_count_cons = pyo.ConstraintList()
    for t in fs.operating_years:
        c.operator_count_cons.add(
            expr=sum(fs.discrete_units[j, t] * fs.workers_per_discrete_unit[j] for j in fs.discrete_stages_set)
            + sum(fs.operators_per_stage[j] * fs.stage_selected[j] for j in fs.processing_stages_set)
            == c.total_operators_unrounded[t]
        )
        c.operator_count_cons.add(expr=c.total_operators_unrounded[t] <= c.total_operators[t])

    # cost of labor
    c.cost_of_labor_cons = pyo.ConstraintList()
    for t in fs.operating_years:
        c.cost_of_labor_cons.add(expr=c.cost_of_labor[t] == c.total_operators[t] * c.labor_rate)

    # TOC expended, revenue, fixed/variable OPEX, overhead
    c.overnight_cost_expended_cons = pyo.ConstraintList()
    c.revenue_cons = pyo.ConstraintList()
    c.variable_operating_cost_cons = pyo.ConstraintList()
    c.fixed_operating_cost_cons = pyo.ConstraintList()
    c.plant_overhead_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        if t <= (t_start + 2):
            c.overnight_cost_expended_cons.add(expr=c.total_overnight_cost_expended[t] == c.total_overnight_cost * schedule[t - t_start])
        else:
            c.overnight_cost_expended_cons.add(expr=c.total_overnight_cost_expended[t] == 0)

        if (t_start + 1) <= t:
            c.revenue_cons.add(
                expr=c.revenue[t]
                == sum(
                    product_values[i]["DF"]
                    * sum(fs.pool_element_flow[i, e, t] * product_values[i][e] for e in fs.tracked_elements_set)
                    for i in fs.produced_pools_set
                    if i in product_values
                )
                + sum(fs.total_byproduct_produced[b, t] * byproduct_values[b] for b in fs.byproducts_set)
            )
            c.fixed_operating_cost_cons.add(
                expr=c.fixed_operating_cost[t]
                == c.cost_of_labor[t]
                + c.maintenance_and_supply_factor * c.total_plant_cost
                + c.sample_analysis_and_qa_qc_factor * c.cost_of_labor[t]
                + c.sales_ip_and_rnd_factor * c.revenue[t]
                + c.admin_and_support_labor_factor * c.cost_of_labor[t]
                + c.fringe_benefits_factor * c.cost_of_labor[t]
                + c.property_tax_and_insurance_factor * c.total_plant_cost
            )
            c.variable_operating_cost_cons.add(
                expr=c.variable_operating_cost[t]
                == sum(fs.total_reagent_consumed[r, t] * reagent_values[r] for r in fs.reagents_set)
                + sum(fs.total_utility_consumed[u, t] * utility_values[u] for u in fs.utilities_set)
                + sum(fs.pool_mass_disposed[i, t] * intermediate_disposal_costs[i]
                      for i in fs.produced_pools_set if i in intermediate_disposal_costs)
                + sum(c.stockpile_holding_cost[i] * fs.feed_stockpiled[i, t] for i in fs.feed_products_set)
            )
        else:
            c.revenue_cons.add(expr=c.revenue[t] == 0)
            c.variable_operating_cost_cons.add(expr=c.variable_operating_cost[t] == 0)
            c.fixed_operating_cost_cons.add(expr=c.fixed_operating_cost[t] == 0)

        c.plant_overhead_cons.add(
            expr=c.plant_overhead[t] == c.plant_overhead_factor * (c.fixed_operating_cost[t] + c.variable_operating_cost[t])
        )


# ---------------------------------------------------------------------------------------
# Cash-flow constraints and NPV
# ---------------------------------------------------------------------------------------
def add_cash_flow_cons(m):
    """Add the discounted cash-flow constraints and the NPV definition."""
    bag = _bag(m)
    fs = m.fs
    c = m.fs.costing
    t_start = bag["plant_start_year"]
    op_esc = c.operating_expense_escalation
    cap_esc = c.capital_escalation
    discount = c.discount_rate

    c.cash_flow_cons = pyo.ConstraintList()
    for t in fs.plant_years:
        if t == t_start:
            c.cash_flow_cons.add(
                expr=c.cash_flow[t]
                == (c.revenue[t] - c.variable_operating_cost[t] - c.fixed_operating_cost[t] - c.plant_overhead[t])
                * ((1 + op_esc) ** (t - t_start))
                - c.total_overnight_cost_expended[t] * ((1 + cap_esc) ** (t - t_start))
            )
        else:
            c.cash_flow_cons.add(
                expr=c.cash_flow[t]
                == (c.revenue[t] - c.variable_operating_cost[t] - c.fixed_operating_cost[t] - c.plant_overhead[t])
                * ((1 + op_esc) ** (t - t_start))
                - (c.total_overnight_cost_expended[t] + c.total_overnight_cost_expansion[t])
                * ((1 + cap_esc) ** (t - t_start))
            )

    c.net_present_value_con = pyo.Constraint(
        expr=c.net_present_value == sum(c.cash_flow[t] / ((1 + discount) ** (t - t_start)) for t in fs.plant_years),
        doc="Net present value (sum of discounted cash flows).",
    )


# ---------------------------------------------------------------------------------------
# Objective function
# ---------------------------------------------------------------------------------------
def add_costing_objective_functions(m, objective_function_choice):
    """Add the objective. Only NPV maximization is implemented.

    Raises:
        NotImplementedError: If ``COST_OF_RECOVERY`` is selected.
    """
    if objective_function_choice == ObjectiveFunctionChoice.NET_PRESENT_VALUE:
        m.fs.costing.objective = pyo.Objective(expr=m.fs.costing.net_present_value, sense=pyo.maximize)
    elif objective_function_choice == ObjectiveFunctionChoice.COST_OF_RECOVERY:
        raise NotImplementedError(
            "Cost-of-recovery is not implemented for the element-flow (bilinear) "
            "superstructure; only NET_PRESENT_VALUE is supported."
        )
    else:
        raise ValueError(f"Unsupported objective function: {objective_function_choice}")


# ---------------------------------------------------------------------------------------
# Environmental impacts (optional block, keyed on the total stage inlet mass)
# ---------------------------------------------------------------------------------------
def add_environmental_impact_params(m, stage_environmental_impacts, epsilon,
                                    augmecon_delta=0.0, augmecon_impact_range=1.0):
    """Add the environmental-impacts sub-block, per-stage intensity, and epsilon bound.

    Args:
        m (ConcreteModel): Pyomo model.
        stage_environmental_impacts (dict): ``{stage: impact intensity per unit total inlet
            mass}``.
        epsilon (float): Upper bound on total lifetime impacts (Pareto constraint).
        augmecon_delta (float): AUGMECON augmentation weight ``delta`` on the normalized surplus
            in the objective. ``0`` (default) recovers the plain epsilon-constraint method; a
            small positive value (Mavrotas 2009) rewards leaving impact budget unused, which
            excludes weakly-Pareto-optimal solutions.
        augmecon_impact_range (float): Impact range ``r = EI_max - EI_min`` used to normalize the
            surplus so a single ``delta`` works across impact categories of different magnitude.
    """
    fs = m.fs
    fs.environmental_impacts = pyo.Block()
    e = fs.environmental_impacts
    e.stage_environmental_impact_intensity = pyo.Param(
        fs.stages_set, initialize=stage_environmental_impacts, mutable=True,
        doc="Environmental-impact intensity per unit of total inlet mass to each stage.",
    )
    e.epsilon = pyo.Param(
        initialize=float(epsilon), mutable=True, domain=pyo.NonNegativeReals,
        doc="Upper bound on total lifetime environmental impacts.",
    )
    e.augmecon_delta = pyo.Param(
        initialize=float(augmecon_delta), mutable=True, domain=pyo.NonNegativeReals,
        doc="AUGMECON augmentation weight on the normalized surplus (0 => plain epsilon-constraint).",
    )
    e.augmecon_impact_range = pyo.Param(
        initialize=float(augmecon_impact_range), mutable=True, domain=pyo.PositiveReals,
        doc="Impact range r = EI_max - EI_min used to normalize the AUGMECON surplus.",
    )


def add_environmental_impact_vars(m):
    """Add per-stage-yearly, yearly-total, and lifetime-total impact variables."""
    fs = m.fs
    e = fs.environmental_impacts
    e.stage_yearly_impact = pyo.Var(fs.stages_set, fs.operating_years, domain=pyo.Reals, doc="Per-stage yearly impact.")
    e.yearly_impact = pyo.Var(fs.operating_years, domain=pyo.Reals, doc="Total yearly impact.")
    e.total_impact = pyo.Var(domain=pyo.Reals, doc="Total lifetime impact.")
    e.surplus = pyo.Var(domain=pyo.NonNegativeReals, doc="AUGMECON surplus: unused impact budget (epsilon - total_impact).")


def add_environmental_impact_cons(m):
    """Add the environmental-impact accounting constraints and the epsilon-constraint.

    Keyed on the element-flow model's total stage inlet mass ``stage_inlet_mass[j, t]``
    (= sum over elements of ``stage_inlet_flow[j, e, t]``), mirroring PrOMMiS's per-option
    inlet-flow formulation.
    """
    fs = m.fs
    e = fs.environmental_impacts

    @e.Constraint(fs.stages_set, fs.operating_years, doc="Per-stage yearly environmental impact.")
    def stage_yearly_impact_con(e, j, t):
        return e.stage_yearly_impact[j, t] == e.stage_environmental_impact_intensity[j] * fs.stage_inlet_mass[j, t]

    @e.Constraint(fs.operating_years, doc="Total yearly environmental impact.")
    def yearly_impact_con(e, t):
        return e.yearly_impact[t] == sum(e.stage_yearly_impact[j, t] for j in fs.stages_set)

    @e.Constraint(doc="Total lifetime environmental impact.")
    def total_impact_con(e):
        return e.total_impact == sum(e.yearly_impact[t] for t in fs.operating_years)

    # AUGMECON epsilon-constraint: equality with a non-negative surplus absorbing the unused
    # budget (Mavrotas 2009). With augmecon_delta == 0 this is equivalent to total_impact <= epsilon.
    @e.Constraint(doc="Epsilon-constraint (AUGMECON equality form) on total lifetime impacts.")
    def epsilon_con(e):
        return e.total_impact + e.surplus == e.epsilon

    # Augment the (already-built) NPV objective with the normalized-surplus reward. The mutable
    # delta/range Params let the Pareto sweep switch augmentation on (and set the range) without
    # rebuilding; delta == 0 leaves the objective at plain max-NPV.
    m.fs.costing.objective.set_value(
        m.fs.costing.net_present_value
        + e.augmecon_delta * e.surplus / e.augmecon_impact_range
    )
