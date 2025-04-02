#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Class to build the costing model for diafiltration-precipitation flowsheet
"""

from pyomo.environ import Expression, Var, units

from idaes.core import FlowsheetCostingBlockData, declare_process_block_class


@declare_process_block_class("DiafiltrationCostingBlock")
class DiafiltrationCostingBlockData(FlowsheetCostingBlockData):
    """
    The base class for costing diafiltration-precipitation flowsheets
    """

    def build(self):
        return super().build()

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

    def _build_common_global_params(self):
        """
        Build the global parameters common to WaterTAP (Diafiltration) Costing Packages.
        The currency units should already be registered.

        The derived class should define the capital_recovery_factor.
        """
        self.utilization_factor = Var(
            initialize=0.9,
            doc="Plant capacity utilization [fraction of uptime]",
            units=units.dimensionless,
        )
        self.utilization_factor.fix()

        self.capital_recovery_factor = Expression(
            expr=0,
            doc="Capital annualization factor [fraction of investment cost/year]",
        )

        self.TPEC = Var(
            initialize=3.4,  # TODO: verify
            doc="Total Purchased Equipment Cost (TPEC)",
            units=units.dimensionless,
        )
        self.TPEC.fix()

        self.TIC = Var(
            initialize=1.65,  # TODO: verify
            doc="Total Installed Cost (TIC)",
            units=units.dimensionless,
        )
        self.TIC.fix()
