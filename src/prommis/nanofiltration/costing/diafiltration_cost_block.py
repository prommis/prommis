"""
Class to build the costing model for diafiltration-precipitation flowsheet

Reference: watertap > watertap > costing > costing_base.py
"""

from idaes.core import FlowsheetCostingBlockData, declare_process_block_class
from pyomo.environ import Expression, Var, units
from watertap.core.util.misc import is_constant_up_to_units


@declare_process_block_class("DiafiltrationCostingBlock")
class DiafiltrationCostingBlockData(FlowsheetCostingBlockData):
    """
    The base class for costing diafiltration-precitipiation flowsheets

    methods that come directly from reference file:
        build
        _build_common_global_params
        _get_costing_method_for
        register_flow_type
    """

    def build(self):
        return super().build()

    # TODO: add any important parameters/constraints to the costing block
    # e.g., energy consumption, levelized cost

    def _build_common_process_costs(self):
        """
        Build the common process costs to WaterTAP (Diafiltration) Costing Packages.
        The currency units should already be registered.

        The derived class should add constraints for total_capital_cost
        and total_operating_cost
        """
        self.total_capital_cost = Var(
            initialize=0,
            doc="Total capital cost of the process",
            units=self.base_currency,
        )
        self.total_operating_cost = Var(
            initialize=0,
            doc="Total operating cost of process per operating period",
            units=self.base_currency / self.base_period,
        )

        # # initialize variables for individual capital costs
        # self.membrane_capital_cost = Var(
        #     initialize=0,
        #     units=self.base_currency,
        #     doc="Total membrane capital cost",
        # )
        # self.pump_capital_cost = Var(
        #     initialize=0,
        #     units=self.base_currency,
        #     doc="Total pump capital cost",
        # )

    def _build_common_global_params(self):
        """
        Build the global parameters common to WaterTAP (Diafiltration) Costing Packages.
        The currency units should already be registered.

        The derived class should define the capital_recovery_factor.
        """

        # TODO: add this variable to operating cost calculations
        self.utilization_factor = Var(
            initialize=0.9,
            doc="Plant capacity utilization [fraction of uptime]",
            units=units.dimensionless,
        )

        # TODO: update this for value used in unit_model/diafiltration.py
        # self.electricity_cost = Var(
        #     initialize=0.07,
        #     doc="Electricity cost",
        #     units=units.USD_2018 / units.kWh,
        # )
        # self.defined_flows["electricity"] = self.electricity_cost

        # self.electrical_carbon_intensity = Var(
        #     initialize=0.475,
        #     doc="Grid carbon intensity [kgCO2_eq/kWh]",
        #     units=units.kg / units.kWh,
        # )

        self.capital_recovery_factor = Expression(
            expr=0,  # TODO: verify
            doc="Capital annualization factor [fraction of investment cost/year]",
        )

        self.TPEC = Var(
            initialize=3.4,  # TODO: verify
            doc="Total Purchased Equipment Cost (TPEC)",
            units=units.dimensionless,
        )

        self.TIC = Var(
            initialize=1.65,  # TODO: verify
            doc="Total Installed Cost (TIC)",
            units=units.dimensionless,
        )

        self.fix_all_vars()

    def _get_costing_method_for(self, unit_model):
        # """
        # Allow the unit model to register its default costing method,
        # either through an attribute named "default_costing_method"
        # or by naming the default costing method "default_costing_method"
        # """
        # if hasattr(unit_model, "default_costing_method"):
        #     return unit_model.default_costing_method
        # return super()._get_costing_method_for(unit_model)
        unit_model.default_costing_method()

    # TODO: determine if this method is still needed
    def register_flow_type(self, flow_type, cost):
        """
        This method allows users to register new material and utility flows
        with the FlowsheetCostingBlock for use when costing flows.
        If `cost` is a constant (up to units), then this method creates a new
        `Var` on the FlowsheetCostingBlock named f`{flow_type}_cost`.
        Otherwise `cost` is a non-constant expression and this method will
        create a new `Expression` on the FlowsheetCostingBlock named
        f`{flow_type}_cost` whose value is fixed to `cost`.

        If a component named f`{flow_type}_cost` already exists on the
        FlowsheetCostingBlock, then an error is raised unless f`{flow_type}_cost`
        is `cost`. If f`{flow_type}_cost` is `cost`, no error is raised and
        the existing component f`{flow_type}_cost` is used to cost the flow.

        Args:
            flow_type: string name to represent flow type
            cost: a Pyomo expression with units representing the flow cost
        """

        flow_cost_name = flow_type + "_cost"
        current_flow_cost = self.component(flow_cost_name)
        if (current_flow_cost is None) and (not is_constant_up_to_units(cost)):
            cost_expr = Expression(expr=cost)
            self.add_component(flow_cost_name, cost_expr)
            super().register_flow_type(flow_type, cost_expr)
        else:
            # all other cases are handled in the base class
            super().register_flow_type(flow_type, cost)
