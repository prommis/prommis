#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Crusher 
=======

Author: Lingyan Deng

The Crusher module includes power consumption for solid crushing. It is a function of particle size distribution, mass flow rate, and bond work index.

Degrees of Freedom
------------------

A Crusher module has two degrees of freedom, which are the output of "particle_size_median" and "particle_size_width". 

Model Structure
---------------

The Crusher model includes one inlet Port (inlet) and one outlet Port (outlet). The properties of the Crusher Unit model is mainly the particle size distribution. 

Additional Constraints
----------------------

Crusher adds one additional constraint to calculate the work required to crush the particles.

.. math:: work_{t} = 10 * m_{t, in} * BWI * \left(\frac{1}{\sqrt{P_{t, prod, 80}}} - \frac{1}{\sqrt{P_{t, feed, 80}}}\right) 

where :math:`work_{t}` is the work required to crush the particles at :math:`t` time, 10 is an empirical value and should not be changed, :math:`m_{t, in}` is the inlet mass flow rate
at :math:`t` time, :math:`BWI` is the bond work index of particles, :math:`P_{t, prod, 80}` is production particle size with 80% passing the mesh at :math:`t` time, :math:`P_{t, feed, 80}` is 
feed particle size with 80% passing the mesh at :math:`t` time.

Expressions
-----------

Crusher includes two expressions to calculate the size of particles that has 80% passing the mesh for both feed and product particles.

.. math:: P_{t, feed, 80} = \frac{S_{t, in, median}}{unit} * \left(-\log(1 - 0.8)\right)^{\frac{SW_{t, in}}{2}}

where :math:`P_{t, feed, 80}` is the feed particle size that has 80% passing the mesh at :math:`t` time, :math:`\frac{S_{t, in, median}}{unit}` is the median particle size of input at :math:`t` time and unitless. The 
default particle size is in micrometer. The :math:`SW_{t, in}` is the particle size width of input at :math:`t` time. 

.. math:: P_{t, prod, 80} = \frac{S_{t, out, median}}{unit} * \left(-\log(1 - 0.8)\right)^{\frac{SW_{t, out}}{2}}

where :math:`P_{t, prod, 80}` is the product particle size that has 80% passing the mesh at :math:`t` time, :math:`\frac{S_{t, out, median}}{unit}` is the median particle size of output at :math:`t` time and unitless. The 
default particle size is in micrometer. The :math:`SW_{t, in}` is the particle size width of output at :math:`t` time. 

Variables
---------

Crusher add the following additional variables beyond those created in property packages.

================ ====== ================================================================================================
Variable         Name   Notes
================ ====== ================================================================================================
:math:`work_{t}`  work
================ ====== ================================================================================================

"""

from functools import partial

from pyomo.common.config import ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Var, log
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Crusher")
class CrusherData(UnitModelBlockData):
    CONFIG = ConfigDict()
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

        # Build state blocks.
        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )
        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )

        tref = self.flowsheet().time.first()
        statevars = self.properties_in[tref].define_state_vars()

        # A constraint of In=Out for all state variables which are not related to particle size.
        for k, v in statevars.items():
            if k not in ["particle_size_median", "particle_size_width"]:
                idx = v.index_set()
                c = Constraint(
                    self.flowsheet().time,
                    idx,
                    doc=f"{k} constraint",
                    rule=partial(_state_rule, state=k),
                )
                self.add_component(k + "_constraint", c)

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
        probability distribution smaller than size x=feed_p80
        """

        sunit = self.properties_in[tref].particle_size_median.get_units()

        @self.Expression(self.flowsheet().time, doc="Feed P80 size")
        def feed_p80(self, t):
            return (
                self.properties_in[t].particle_size_median
                / sunit
                * (-log(1 - 0.8)) ** (self.properties_in[t].particle_size_width / 2)
            )

        @self.Expression(self.flowsheet().time, doc="Product P80 size")
        def prod_p80(self, t):
            return (
                self.properties_out[t].particle_size_median
                / sunit
                * (-log(1 - 0.8)) ** (self.properties_out[t].particle_size_width / 2)
            )

        @self.Constraint(self.flowsheet().time, doc="Crusher work constraint")
        def crush_work_eq(self, t):
            return self.work[t] == (
                10  # 10 is an empirical correlation, this should not be changed.
                * self.properties_in[t].flow_mass
                * self.config.property_package.bond_work_index
                * (1 / (self.prod_p80[t]) ** 0.5 - 1 / (self.feed_p80[t]) ** 0.5)
            )

    def _get_stream_table_contents(self, time_point=0):
        # Dictionary to hold data for all streams
        io_dict = {"Inlet": self.properties_in, "Outlet": self.properties_out}
        return create_stream_table_dataframe(io_dict, time_point=time_point)

    def _get_performance_contents(self, time_point=0):
        # Report
        var_dict = {
            "Work Required (W)": value(self.work[time_point]),
        }
        return {"vars": var_dict}


def _state_rule(b, time, index, state):
    sin = b.properties_in[time].define_state_vars()[state]
    sout = b.properties_out[time].define_state_vars()[state]
    return sin[index] == sout[index]
