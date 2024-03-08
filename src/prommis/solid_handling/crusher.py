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

from pyomo.environ import Var, Expression, Constraint, NonNegativeReals, log
from pyomo.common.config import ConfigValue, ConfigBlock, In
from pyomo.network import Port
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CrushAndBreakageUnit")
class CrushAndBreakageUnitData(UnitModelBlockData):
    CONFIG = (
        ConfigBlock()
    )  # or  CONFIG = UnitModelBlockData.CONFIG() (not sure what's the difference)
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Crush units do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Crush units do not have defined volume, thus
this must be False.""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
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
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
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
        super(CrushAndBreakageUnitData, self).build()

        # Build Control Volume
        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_phase_equilibrium=False
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=False,
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        solid_prop = self.config.property_package
        solid_prop.flow_mass  # mass low rate, property unit kg/hr, unit needed tonne/hr
        solid_prop.feed50size  # Feed median particle size, micrometer
        solid_prop.prod50size  # Product median particle size, micrometer

        units_meta = self.control_volume.config.property_package.get_metadata()
        self.feed80size = Var(
            self.flowsheet().time,
            initialize=200.0,
            doc="Feed Particle Size that has 80% passing, micrometer",
            units=units_meta.get_derived_units(
                "length"
            ),  # unit should be changed to micrometer
        )

        # BreakageDistribution calculation as a constraint. This is the equation for accumulative fraction of solid breakage probability distribution smaller than size x=feed80size
        @self.Constraint(self.flowsheet().time, doc="feed size constraint")
        def feed_size_eq(self, t):
            return self.feed80size[t] == (
                (log(1 - self.control_volume.properties_in[t].probfeed80))
                ** (self.control_volume.properties_in[t].nfeed / 2)
                * self.control_volume.properties_in[t].feed50size
            )

        self.prod80size = Var(
            self.flowsheet().time,
            initialize=80.0,
            doc="Feed Particle Size that has 80% passing, micrometer",
            units=units_meta.get_derived_units(
                "length"
            ),  # unit should be changed to micrometer
        )

        @self.Constraint(self.flowsheet().time, doc="product size constraint")
        def prod_size_eq(self, t):
            return self.prod80size[t] == (
                (log(1 - self.control_volume.properties_in[t].probprod80))
                ** (self.control_volume.properties_in[t].nprod / 2)
                * self.control_volume.properties_in[t].prod50size
            )

        self.crushpower = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Work required to increase crush the solid",
            units=units_meta.get_derived_units("power"),
        )

        @self.Constraint(self.flowsheet().time, doc="Crusher work constraint")
        def crush_work_eq(self, t):
            return self.crushpower[t] == (
                10
                * self.control_volume.properties_in[t].flow_mass
                * self.control_volume.properties_in[t].bwi
                * (
                    1 / self.control_volume.properties_in[t].prod80size ** 0.5
                    - 1 / self.control_volume.properties_in[t].feed80size ** 0.5
                )
            )
