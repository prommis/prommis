#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module including power consumption for solid crushing; breakage probability and distribution.
"""

__author__ = "Lingyan Deng"
__version__ = "1.0.0"

from pyomo.environ import Var, Constraint, exp 
from pyomo.common.config import ConfigBlock, ConfigValue
from idaes.core import declare_process_block_class
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

@declare_process_block_class("CrushAndBreakageUnit")
class CrushAndBreakageUnitData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package"))

    def build(self):
        super(CrushAndBreakageUnitData, self).build()

        self.BWI = Var(initialize=12, doc="Bond Work Index, kWh/tonne")
        self.F80 = Var(initialize=200, doc="Feed Particle Size that allows 80% passing, micrometer")
        self.P80 = Var(initialize=80, doc="Product Particle Size that allows 80% passing, micrometer")
        self.Massflow = Var(initialize=2, doc="Mass flow rate, tonne/h")
        self.x = Var(initialize=60, doc="Particle size, micrometer")
        self.x50 = Var(initialize=80, doc="Median particle size, micrometer")
        self.n = Var(initialize=1.5, doc="Distribution width parameter")

        # CrushPower calculation as a constraint. This is the solid crushPower calculation as a constraint. 
        self.CrushPower = Constraint(
            expr=10 * self.Massflow * self.BWI * (1 / self.P80**0.5 - 1 / self.F80**0.5) == 10 * self.Massflow * self.BWI * (1 / self.P80**0.5 - 1 / self.F80**0.5))

        # BreakageDistribution calculation as a constraint. This is the equation for accumulative fraction of solid breakage probability distribution smaller than size x
        self.BreakageDistribution = Constraint(
            expr=1 - exp(-(self.x / self.x50) ** self.n) == self.soliddistribution)
