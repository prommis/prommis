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

from functools import partial
from pyomo.environ import Var, log, Constraint, units as pyunits
from pyomo.common.config import ConfigValue, ConfigDict, In
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CrushAndBreakageUnit")
class CrushAndBreakageUnitData(UnitModelBlockData):
    CONFIG = (
        ConfigDict()
    )  # or  CONFIG = UnitModelBlockData.CONFIG() (not sure what's the difference)
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Crush unit is steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Crush unit has no holdup.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigDict(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigDict with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state = True,
            **self.config.property_package_args,
        )
        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state = True,
            **self.config.property_package_args,
        )
       
        tref = self.flowsheet().time.first()
        statevars = self.properties_in[tref].define_state_vars()
        
        for k, v in statevars.items():
            if k not in ["particle_size_median", "particle_size_width"]:
                idx = v.index_set()
                c = Constraint(self.flowsheet().time, idx, doc = f"{k} constraint", 
                               rule = partial(_state_rule, state = k))
                self.add_component(k+"_constraint", c)



        # Add Ports
        self.add_port("inlet", self.properties_in)
        self.add_port("outlet", self.properties_out)

        self.work = Var(
            self.flowsheet().time,
            units=pyunits.W,
            initialize=3915.17,
            doc="Work required to increase crush the solid",
        )

        """ Breakage Distribution calculation as a constraint. 
        This is the equation for accumulative fraction of solid breakage 
        probability distribution smaller than size x=feed80size
        """
        
        sunit = self.properties_in[tref].particle_size_median.get_units()

        @self.Expression(self.flowsheet().time, doc="Feed P80 Size")
        def feed_p80(self, t):
            return (
                self.properties_in[t].particle_size_median
                / sunit
                * (-log(1 - 0.8))
                ** (self.properties_in[t].particle_size_width / 2)
            )

        @self.Expression(self.flowsheet().time, doc="product p80 size")
        def prod_p80(self, t):
            return (
                self.properties_out[t].particle_size_median
                / sunit
                * (-log(1 - 0.8))
                ** (self.properties_out[t].particle_size_width / 2)
            )

        @self.Constraint(self.flowsheet().time, doc="Crusher work constraint")
        def crush_work_eq(self, t):
            return self.work[t] == (
                10  # 10 is an empirical correlation, this should not be changed.
                * self.properties_in[t].flow_mass
                * self.config.property_package.bond_work_index
                * (1 / (self.prod_p80[t]) ** 0.5 - 1 / (self.feed_p80[t]) ** 0.5)
            )


def _state_rule(b, time, index, state):
    sin = b.properties_in[time].define_state_vars()[state]
    sout = b.properties_out[time].define_state_vars()[state]
    return sin[index] == sout[index]