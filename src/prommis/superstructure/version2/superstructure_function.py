#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Element-Flow (Bilinear Stream-Mixing) Superstructure Function
=============================================================

``build_model`` constructs a multiperiod, non-convex MIQCP superstructure that identifies
optimal pathways for recovering elements from end-of-life products while **tracking variable
stream compositions through mixing and splitting nodes** (a generalized pooling/blending
problem). Composition is carried as per-element mass flows; the only bilinear block is the
composition-preserving split. The function is feedstock-agnostic: all feed, topology,
disassembly, price, and equipment data is supplied by the caller.

This mirrors the layout of the PrOMMiS fixed-yield superstructure
(``prommis.superstructure``) — a thin orchestrator over modular ``add_*`` builders, the
``ObjectiveFunctionChoice`` enum, custom Pyomo units, and a ``CustomScalerBase`` scaler —
but implements the element-flow formulation rather than the single-extent fixed-yield MILP.

Generic global-convergence cuts (:mod:`.convergence_cuts`) are applied at the end of the build
by default and can be toggled off via the ``add_*_cuts`` arguments. They are feedstock-agnostic
-- auto-detected from ``feed_composition`` and the topology -- so they remain valid for any
element set or feeds.

Author: Chris Laliwala
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from idaes.core.scaling import CustomScalerBase

from .objective_function_enums import ObjectiveFunctionChoice
from .convergence_cuts import (
    apply_component_invariance_cuts,
    apply_bounded_component_cuts,
    apply_routing_sum_cuts,
    apply_ratio_cuts,
)
from .add_superstructure_blocks import (
    add_capacity_workforce_cons,
    add_capacity_workforce_vars,
    add_capital_cost_cons,
    add_cash_flow_cons,
    add_costing_objective_functions,
    add_costing_params,
    add_costing_vars,
    add_environmental_impact_cons,
    add_environmental_impact_params,
    add_environmental_impact_vars,
    add_feed_params,
    add_mass_balance_cons,
    add_mass_balance_vars,
    add_objective_function_choice_param,
    add_operating_cost_cons,
    add_operating_params,
    add_plant_lifetime_params,
    add_supe_formulation_params,
)


def define_custom_units():
    """Define the custom Pyomo units used throughout the superstructure."""
    pyunits.load_definitions_from_strings(["EOL_Product = [item]"])
    pyunits.load_definitions_from_strings(["Operator = [item]"])
    pyunits.load_definitions_from_strings(["Disassembly_Unit = [item]"])
    pyunits.load_definitions_from_strings(["USD = [currency]"])


class SuperstructureScaler(CustomScalerBase):
    """Custom scaler for the element-flow superstructure (mirrors PrOMMiS's scaler).

    Provided for layout parity and optional IDAES-based scaling. It is **not** required for
    the documented solve recipe (Gurobi ``NonConvex=2`` is invoked directly on the model).
    """

    def variable_scaling_routine(self, model, overwrite: bool = False, submodel_scalers: dict = None):
        fs = model.fs
        for _, v in fs.pool_element_flow.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for _, v in fs.stage_inlet_flow.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for _, v in fs.stage_outlet_flow.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)
        for _, v in fs.stage_inlet_mass.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)

        self.set_variable_scaling_factor(fs.costing.net_present_value, 1e-6, overwrite=overwrite)
        for _, v in fs.costing.cash_flow.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for _, v in fs.costing.revenue.items():
            self.set_variable_scaling_factor(v, 1e-5, overwrite=overwrite)
        for _, v in fs.costing.bare_equipment_cost.items():
            self.set_variable_scaling_factor(v, 1e-4, overwrite=overwrite)

    def constraint_scaling_routine(self, model, overwrite: bool = False, submodel_scalers: dict = None):
        fs = model.fs
        for _, con in fs.stage_inlet_mass_con.items():
            self.scale_constraint_by_nominal_value(con, scheme="inverse_maximum", overwrite=overwrite)
        for _, con in fs.recovery_con.items():
            self.scale_constraint_by_nominal_value(con, scheme="inverse_maximum", overwrite=overwrite)
        for _, con in fs.costing.cash_flow_cons.items():
            self.scale_constraint_by_nominal_value(con, scheme="inverse_maximum", overwrite=overwrite)


def build_model(
    # objective
    objective_function_choice,
    # plant lifetime
    plant_start_year,
    plant_end_year,
    # feed (no feedstock-specific defaults)
    feed_products,
    tracked_elements,
    feed_composition,
    feed_availability,
    # superstructure topology
    stages,
    discrete_stages,
    stage_inlet_pools,
    stage_outlet_pools,
    produced_pools,
    final_product_pools,
    stage_efficiencies,
    stage_split_fractions,
    # disassembly / processing sizing
    discrete_unit_rate,
    workers_per_discrete_unit,
    operators_per_stage,
    stage_capacity_big_m,
    capital_cost_per_discrete_unit,
    discretized_equipment_cost,
    # utilities / reagents / byproducts
    utilities,
    stage_utilities,
    utility_consumption,
    reagents,
    stage_reagents,
    reagent_consumption,
    byproducts,
    stage_byproducts,
    byproduct_generation,
    # economics (prices)
    product_values,
    reagent_values,
    utility_values,
    intermediate_disposal_costs,
    byproduct_values,
    labor_rate,
    # numerics & options
    flow_upper_bound,
    stockpile_holding_cost=None,
    purity_gated_stages=None,
    initial_state=None,
    # global-convergence cuts (applied after the model is built; on by default)
    add_component_invariance_cuts=True,
    add_bounded_component_cuts=True,
    add_routing_sum_cuts=True,
    add_ratio_cuts=True,
    # environmental impacts (optional)
    consider_environmental_impacts=False,
    stage_environmental_impacts=None,
    epsilon=None,
    augmecon_delta=0.0,
    augmecon_impact_range=1.0,
    # generic QGESS economic factors (defaults; NOT feedstock-specific)
    lang_factor=2.97,
    operating_expense_escalation=0.03,
    capital_escalation=0.036,
    discount_rate=0.0577,
    total_overnight_cost_factor=1.177,
    capex_expenditure_schedule=(0.1, 0.6, 0.3),
    min_production_fraction=0.5,
    capacity_expansion_lower_bound=0.0,
    capacity_expansion_upper_bound=1e6,
    maintenance_and_supply_factor=0.02,
    sample_analysis_and_qa_qc_factor=0.1,
    sales_ip_and_rnd_factor=0.01,
    admin_and_support_labor_factor=0.2,
    fringe_benefits_factor=0.25,
    property_tax_and_insurance_factor=0.01,
    plant_overhead_factor=0.2,
):
    """Build the element-flow (bilinear stream-mixing) superstructure model.

    Args:
        objective_function_choice (ObjectiveFunctionChoice): ``NET_PRESENT_VALUE`` (only
            implemented option) or ``COST_OF_RECOVERY`` (raises ``NotImplementedError``).
        plant_start_year (int): Construction year; production runs ``plant_start_year + 1``
            through ``plant_end_year``.
        plant_end_year (int): Final planning-horizon year.
        feed_products (list): Feed (end-of-life product) streams.
        tracked_elements (list): Tracked elements (e.g. ``["Nd", "Dy", "Fe"]``).
        feed_composition (dict): ``{feed: {element: mass_fraction}}``.
        feed_availability (dict): ``{feed: {year: amount_available}}``.
        stages (list): All stage identifiers.
        discrete_stages (list): Stages using discrete disassembly units.
        stage_inlet_pools (dict): ``{stage: [pools entering]}``.
        stage_outlet_pools (dict): ``{stage: [pools produced]}``.
        produced_pools (list): Non-feed pools produced within the process.
        final_product_pools (list): Saleable final-product pools.
        stage_efficiencies (dict): ``{stage: {element: retention_efficiency}}``.
        stage_split_fractions (dict): ``{stage: {pool: {element: split_fraction}}}``.
        discrete_unit_rate (dict): ``{discrete_stage: EOL products / unit / year}``.
        workers_per_discrete_unit (dict): ``{discrete_stage: workers / unit}``.
        operators_per_stage (dict): ``{processing_stage: operators}``.
        stage_capacity_big_m (dict): ``{processing_stage: big-M capacity bound}``.
        capital_cost_per_discrete_unit (dict): ``{discrete_stage: $/unit}``.
        discretized_equipment_cost (dict): ``{stage: {"Capacity": {...}, "Costs": {...}}}``
            piecewise-linear equipment-cost breakpoints for processing stages.
        utilities, stage_utilities, utility_consumption: utility sets / connectivity / rates.
        reagents, stage_reagents, reagent_consumption: reagent sets / connectivity / rates.
        byproducts, stage_byproducts, byproduct_generation: byproduct sets / connectivity / rates.
        product_values (dict): ``{pool: {element: price, "DF": discount_factor}}``.
        reagent_values (dict): ``{reagent: price}``.
        utility_values (dict): ``{utility: price}``.
        intermediate_disposal_costs (dict): ``{pool: disposal cost}``.
        byproduct_values (dict): ``{byproduct: value (negative = disposal cost)}``.
        labor_rate (float): Yearly wage per operator.
        flow_upper_bound (float): Global upper bound on flow variables.
        purity_gated_stages (dict | None): ``{stage: allowed_element}`` inlet purity gate.
        initial_state (dict | None): Inherited boundary state for rolling-horizon runs.
        add_component_invariance_cuts (bool): Linearize the split of every element with a constant
            cross-feed fraction (default ``True``). See
            :mod:`new_superstructure_function.superstructure.convergence_cuts`.
        add_bounded_component_cuts (bool): Apply the two-sided bracket for every element whose
            cross-feed fraction varies within a narrow band (default ``True``); complementary to the
            component-invariance cut, which handles the exactly-constant elements.
        add_routing_sum_cuts (bool): Apply the ``sum_j routing_fraction <= 1`` cut (default ``True``).
        add_ratio_cuts (bool): Apply the pairwise composition-ratio cuts (default ``True``).
        consider_environmental_impacts (bool): Add the optional environmental-impacts block.
        stage_environmental_impacts (dict | None): ``{stage: impact intensity}`` (required if
            ``consider_environmental_impacts``).
        epsilon (float | None): Upper bound on total lifetime impacts (Pareto constraint).
        augmecon_delta (float): AUGMECON augmentation weight (Mavrotas 2009); 0 (default) gives
            the plain epsilon-constraint, a small positive value excludes weakly-Pareto solutions.
        augmecon_impact_range (float): Impact range r = EI_max - EI_min normalizing the surplus.
        lang_factor, ... plant_overhead_factor: generic QGESS economic factors (defaults).

    Returns:
        ConcreteModel: the Pyomo superstructure model (objective on ``m.fs.costing``).
    """
    define_custom_units()

    m = pyo.ConcreteModel()

    add_objective_function_choice_param(m, objective_function_choice)
    add_plant_lifetime_params(m, plant_start_year, plant_end_year)
    add_feed_params(m, feed_products, tracked_elements, feed_composition, feed_availability, flow_upper_bound)
    add_supe_formulation_params(
        m, stages, discrete_stages, stage_inlet_pools, stage_outlet_pools, produced_pools,
        final_product_pools, stage_efficiencies, stage_split_fractions, discrete_unit_rate,
        workers_per_discrete_unit, operators_per_stage, stage_capacity_big_m,
    )
    add_operating_params(
        m, utilities, stage_utilities, utility_consumption, reagents, stage_reagents,
        reagent_consumption, byproducts, stage_byproducts, byproduct_generation, product_values,
        reagent_values, utility_values, intermediate_disposal_costs, byproduct_values, labor_rate,
    )
    add_costing_params(
        m, capital_cost_per_discrete_unit, discretized_equipment_cost, lang_factor,
        operating_expense_escalation, capital_escalation, discount_rate, total_overnight_cost_factor,
        capex_expenditure_schedule, min_production_fraction, capacity_expansion_lower_bound,
        capacity_expansion_upper_bound, maintenance_and_supply_factor, sample_analysis_and_qa_qc_factor,
        sales_ip_and_rnd_factor, admin_and_support_labor_factor, fringe_benefits_factor,
        property_tax_and_insurance_factor, plant_overhead_factor, stockpile_holding_cost,
    )

    add_mass_balance_vars(m)
    add_mass_balance_cons(m, purity_gated_stages, initial_state)
    add_capacity_workforce_vars(m)
    add_capacity_workforce_cons(m, initial_state)

    add_costing_vars(m)
    add_capital_cost_cons(m)
    add_operating_cost_cons(m)
    add_cash_flow_cons(m)
    add_costing_objective_functions(m, objective_function_choice)

    if consider_environmental_impacts:
        if stage_environmental_impacts is None or epsilon is None:
            raise ValueError(
                "stage_environmental_impacts and epsilon must be provided when "
                "consider_environmental_impacts is True."
            )
        add_environmental_impact_params(
            m, stage_environmental_impacts, epsilon,
            augmecon_delta=augmecon_delta, augmecon_impact_range=augmecon_impact_range,
        )
        add_environmental_impact_vars(m)
        add_environmental_impact_cons(m)

    # Global-convergence cuts (feedstock-agnostic; on by default). Applied last so the full
    # model -- in particular the bilinear split -- already exists.
    if add_component_invariance_cuts:
        apply_component_invariance_cuts(m, feed_composition)
    if add_bounded_component_cuts:
        apply_bounded_component_cuts(m, feed_composition)
    if add_routing_sum_cuts:
        apply_routing_sum_cuts(m)
    if add_ratio_cuts:
        apply_ratio_cuts(m, feed_composition)

    return m
