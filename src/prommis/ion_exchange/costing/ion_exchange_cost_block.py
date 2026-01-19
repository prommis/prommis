#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.

# This model was developed using the following model as a reference:
# https://github.com/prommis/prommis/blob/main/src/prommis/nanofiltration/costing/diafiltration_cost_block.py
#####################################################################################################

#################################################################################
# This model was derived from:
# https://github.com/watertap-org/watertap/blob/main/watertap/costing/watertap_costing_package.py
# and
# https://github.com/watertap-org/watertap/blob/723576a4a05596cc1f897376863c965f22efa99f/watertap/costing/util.py

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from idaes.core import FlowsheetCostingBlockData, declare_process_block_class


""" modified by: Soraya Rawlings """


@declare_process_block_class("IXCostingBlock")
class IXCostingBlockData(FlowsheetCostingBlockData):
    """
    Base class for costing the ion exchange model
    """

    def build(self):
        return super().build()

    def _build_common_process_costs(self):
        """Build the common process costs to WaterTAP Ion Exchange Costing
        Packages.  The currency units should already be registered.

        The derived class should add constraints for total_capital_cost
        and total_operating_cost

        """

        self.total_capital_cost = pyo.Var(
            initialize=0,
            doc="Total capital cost of the process",
            units=self.base_currency,
        )
        self.total_operating_cost = pyo.Var(
            initialize=0,
            doc="Total operating cost of process per operating period",
            units=self.base_currency / self.base_period,
        )
        self.electricity_cost = pyo.Var(
            initialize=0.07,
            doc="Electricity cost",
            units=self.base_currency / pyo.units.kWh,
        )
        self.electricity_cost.fix()

        # self.flow_type = pyo.Set()
        # self.flow_types.add("electricity")
        self.defined_flows["electricity"] = self.electricity_cost

    def _build_common_global_params(self):
        """Build the global parameters common to WaterTAP Ion Exchange
        Costing Packages.  The currency units should already be
        registered.

        The derived class should define the capital_recovery_factor.
        """

        self.utilization_factor = pyo.Var(
            initialize=1,
            doc="Plant capacity utilization [fraction of uptime]",
            units=pyo.units.dimensionless,
        )
        self.utilization_factor.fix()

        self.capital_recovery_factor = pyo.Expression(
            expr=0.1,
            doc="Capital annualization factor [fraction of investment cost/year]",
        )

        self.TPEC = pyo.Var(
            initialize=4.121212121212121,
            doc="Total Purchased Equipment Cost (TPEC)",
            units=pyo.units.dimensionless,
        )
        self.TPEC.fix()

        self.TIC = pyo.Var(
            initialize=2.0,
            doc="Total Installed Cost (TIC)",
            units=pyo.units.dimensionless,
        )
        self.TIC.fix()
