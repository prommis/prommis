##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2023
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
###############################################################################
"""
IDAES REE oxalate salt roaster unit model.

This is a heterogeneous reactor model for roasting rare earth oxalate salts
produced from the oxalate precipitation process to form rare earth oxides.

The reactions involved are listed below:
RE2(C2O4)3_xH2O + 1.5O2 -> RE2O3 + 6CO2(g) + xH2O(g)
Fe2(C2O4)3_2H2O + 1.5O2 -> Fe2O3 + 6CO2(g) + 2H2O(g)
Al2(C2O4)2_H2O + 1.5O2 -> Al2O3 + 6CO2(g) + H2O(g)
Note that Fe and Al are considered as impurity in the solid reactant.

Surface moisture is vaporized to become H2O in the gas phase

The feed stream is from oxalate acid precipitator and solid filter
which may contain surface moisture.
The heat for the reaction can be provided by either external heating, modeled
as a user-specified heat duty or by the combustion of fossil fuel with air to form
a hot O2-containing flue gas entering the reactor through a gas inlet port.
The O2 is used to oxidize the feed salt to form REE oxides.
The gas inlet should contain enough O2 for the oxidation reactions.
The gas outlet port represents the gas product stream leaving the reactor.
The solid product streams include a recovered final product stream and
a dust stream containing the fine particles which are not recovered.

Current model assumes that all hydrous oxalates are decomposted without
considering the kinetics or mass transfer. Therefore, it is a steady-state
model only.

Full material and energy balances are modeled. Note that the standard heats of formation of
many RE oxalates and RE oxides are unavailable in the literature and, therefore, the average
values of La, Ce, Pr, and Nd compounds are used in case of unavailable data.

The outlet gas product stream is assumed to have the same temperature as the two solid
product streams.
if the product temperature is specified as a user input, the heat duty added to the
reactor is calculately. If the heat duty is given, the product temperature will be calculated.

Currently, no ports for the solid inlet and outlet streams are used. The mass flow rate
and composition of the solid reactant are specified as input variables inside the model.
The mass flow rate and the composition of the product streams are also declared as model
variables.

"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, Bool

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault

from idaes.core.util.config import is_physical_parameter_block, DefaultBool
import idaes.logger as idaeslog

# Additional import for the unit operation
from pyomo.environ import (
    Var,
    Param,
    Set,
    units as pyunits,
)
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

__author__ = "Jinliang Ma"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------------------------------------
@declare_process_block_class("REEOxalateRoaster")
class REEOxalateRoasterData(UnitModelBlockData):
    """
    Simple 0D roaster model with mass and energy balance only
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=DefaultBool,
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
        "property_package_gas",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Gas property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_gas",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing gas property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "property_package_precipitate",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for precipitate control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_precipitate",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing precipitate property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "metal_list",
        ConfigValue(
            default=[
                "Al",
                "Fe",
                "Sc",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Gd",
                "Dy",
            ],
            domain=list,
            description="List of components in solid oxalate feed",
            doc="""A dict of the components of interest in the mixture.
        Keys are component names and values are configuration arguments to
        be passed to Component on construction.
        """,
        ),
    )

    def build(self):
        # Call TranslatorData build to setup dynamics
        super(REEOxalateRoasterData, self).build()

        # Build Holdup Block
        # gas phase inlet stream
        self.gas_in = self.config.property_package_gas.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args_gas,
        )
        # gas phase outlet stream
        self.gas_out = self.config.property_package_gas.build_state_block(
            self.flowsheet().time, **self.config.property_package_args_gas
        )
        # solid phase inlet stream from precipitator
        self.solid_in = self.config.property_package_precipitate.build_state_block(
            self.flowsheet().time, **self.config.property_package_args_precipitate
        )
        self.add_port("gas_inlet", self.gas_in)
        self.add_port("gas_outlet", self.gas_out)
        self.add_port("solid_inlet", self.solid_in)

        # Add Geometry
        if self.config.has_holdup is True:
            self.volume = Var(
                initialize=1, units=pyunits.m**3, doc="volume of the reactor"
            )
            self.volume.fix()

        # Add heat duty
        if self.config.has_heat_transfer is True:
            self.heat_duty = Var(
                self.flowsheet().time,
                initialize=0,
                units=pyunits.W,
                doc="heat duty added to the reactor",
            )

        if self.config.has_pressure_change is True:
            self.deltaP = Var(
                self.flowsheet().time,
                initialize=0,
                units=pyunits.Pa,
                doc="pressure drop from inlet to outlet",
            )

        # metal list is a user-specified list of metals contained in the solid feed stream including impurity metals
        self.metal_list = Set(
            initialize=self.config.metal_list, doc="List of metals in solid feed stream"
        )

        # Construct performance equations
        self._make_params()
        self._make_vars()
        self._make_mass_balance()
        self._make_energy_balance()
        self._make_momentum_balance()

    def _make_params(self):
        """This section is for parameters within this model."""
        # reference temperature at standard condition
        self.temp_ref = Param(initialize=298.15, units=pyunits.K)

        # atomic mass of elements involved in kg/mol
        self.metal_list_all = [
            "Al",
            "Fe",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "Th",
        ]
        self.am_metal_list_all = Param(
            self.metal_list_all, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.xH2O_oxalate_list_all = Param(
            self.metal_list_all, initialize=10, mutable=True
        )
        self.mw_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=0.1,
            mutable=True,
            units=pyunits.kg / pyunits.mol,
        )
        self.mw_oxide_list_all = Param(
            self.metal_list_all,
            initialize=0.1,
            mutable=True,
            units=pyunits.kg / pyunits.mol,
        )
        self.am_H = Param(initialize=0.0010078, units=pyunits.kg / pyunits.mol)
        self.am_C = Param(initialize=0.012011, units=pyunits.kg / pyunits.mol)
        self.am_N = Param(initialize=0.014007, units=pyunits.kg / pyunits.mol)
        self.am_O = Param(initialize=0.015999, units=pyunits.kg / pyunits.mol)
        self.am_metal_list_all["Al"] = 0.026982
        self.am_metal_list_all["Fe"] = 0.055845
        self.am_metal_list_all["Sc"] = 0.044956
        self.am_metal_list_all["Y"] = 0.088906
        self.am_metal_list_all["La"] = 0.13891
        self.am_metal_list_all["Ce"] = 0.14012
        self.am_metal_list_all["Pr"] = 0.14091
        self.am_metal_list_all["Nd"] = 0.14424
        self.am_metal_list_all["Pm"] = 0.145
        self.am_metal_list_all["Sm"] = 0.15036
        self.am_metal_list_all["Eu"] = 0.15196
        self.am_metal_list_all["Gd"] = 0.15725
        self.am_metal_list_all["Tb"] = 0.15893
        self.am_metal_list_all["Dy"] = 0.1625
        self.am_metal_list_all["Ho"] = 0.16493
        self.am_metal_list_all["Er"] = 0.16726
        self.am_metal_list_all["Tm"] = 0.16893
        self.am_metal_list_all["Yb"] = 0.17304
        self.am_metal_list_all["Lu"] = 0.17497
        self.am_metal_list_all["Th"] = 0.23204
        self.xH2O_oxalate_list_all["Al"] = 1
        self.xH2O_oxalate_list_all["Fe"] = 2

        self.mw_H2O = Param(
            initialize=self.am_H * 2 + self.am_O, units=pyunits.kg / pyunits.mol
        )

        for i in self.metal_list_all:
            self.mw_oxalate_list_all[i] = (
                self.am_metal_list_all[i] * 2
                + 3 * (self.am_C * 2 + self.am_O * 4)
                + self.xH2O_oxalate_list_all[i] * (self.am_H + self.am_O)
            )
        for i in self.metal_list_all:
            self.mw_oxide_list_all[i] = self.am_metal_list_all[i] * 2 + 3 * self.am_O

        # molar standard enthalpy at 298.15 K initialized to the average of La, Ce, Pr, and Nd
        self.enth0_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=-6350044,
            mutable=True,
            units=pyunits.J / pyunits.mol,
        )
        # molar standard enthalpy at 298.15 K initialized to the average of La, Ce, Pr, and Nd
        self.enth0_oxide_list_all = Param(
            self.metal_list_all,
            initialize=-1801865.75,
            mutable=True,
            units=pyunits.J / pyunits.mol,
        )
        # molar heat capacity of oxalate initialized to the average of La, Ce, Pr, and Nd
        self.cp0_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=45.751,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp1_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=0,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        # molar heat capacity of oxide initialized to the average of La, Ce, Pr, and Nd
        self.cp0_oxide_list_all = Param(
            self.metal_list_all,
            initialize=110.3625,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp1_oxide_list_all = Param(
            self.metal_list_all,
            initialize=0.033887,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )

        self.enth0_oxalate_list_all["Al"] = -3397000
        self.enth0_oxalate_list_all["Fe"] = -6782000
        self.enth0_oxalate_list_all["La"] = -5916176
        self.enth0_oxalate_list_all["Ce"] = -6782000
        self.enth0_oxalate_list_all["Pr"] = -5920000
        self.enth0_oxalate_list_all["Nd"] = -6782000

        self.enth0_oxide_list_all["Fe"] = -825500
        self.enth0_oxide_list_all["Al"] = -1675000
        self.enth0_oxide_list_all["La"] = -1793702
        self.enth0_oxide_list_all["Ce"] = -1796191
        self.enth0_oxide_list_all["Pr"] = -1809664
        self.enth0_oxide_list_all["Nd"] = -1807906
        self.enth0_oxide_list_all["Sc"] = -1908866
        self.enth0_oxide_list_all["Y"] = -1932800

        self.cp0_oxalate_list_all["La"] = 45.751
        self.cp0_oxalate_list_all["Ce"] = 45.751
        self.cp0_oxalate_list_all["Pr"] = 45.751
        self.cp0_oxalate_list_all["Nd"] = 45.751

        self.cp1_oxalate_list_all["La"] = 0
        self.cp1_oxalate_list_all["Ce"] = 0
        self.cp1_oxalate_list_all["Pr"] = 0
        self.cp1_oxalate_list_all["Nd"] = 0

        self.cp0_oxide_list_all["La"] = 107.72
        self.cp0_oxide_list_all["Ce"] = 115.78
        self.cp0_oxide_list_all["Pr"] = 112.82
        self.cp0_oxide_list_all["Nd"] = 105.13

        self.cp1_oxide_list_all["La"] = 0.026114
        self.cp1_oxide_list_all["Ce"] = 0.03477
        self.cp1_oxide_list_all["Pr"] = 0.034364
        self.cp1_oxide_list_all["Nd"] = 0.0403

        # unit constants used for the expressions of liquid water enthalpy
        self.enth_mol_const = Param(
            initialize=1,
            units=pyunits.J / pyunits.mol,
            doc="1 unit of molar enthalpy in J/mol",
        )
        self.cp_mas_const = Param(
            initialize=1,
            units=pyunits.J / pyunits.kg / pyunits.K,
            doc="1 unit of mass heat capacity in J/kg-K",
        )

    def _make_vars(self):
        """This section declares variables within this model."""
        self.flow_mol_comp_feed = Var(
            self.flowsheet().config.time,
            self.metal_list,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of oxalate in solid feed stream",
        )

        self.flow_mol_moist_feed = Var(
            self.flowsheet().config.time,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of liquid water in solid feed stream",
        )

        self.frac_comp_recovery = Var(
            self.flowsheet().config.time,
            self.metal_list,
            initialize=0.95,
            doc="fraction of oxide recovery",
        )

        self.flow_mol_comp_product = Var(
            self.flowsheet().config.time,
            self.metal_list,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of oxide in product stream",
        )

        self.flow_mol_comp_dust = Var(
            self.flowsheet().config.time,
            self.metal_list,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of oxide in dust stream",
        )

        self.enth_mol_comp_feed = Var(
            self.flowsheet().config.time,
            self.metal_list,
            units=pyunits.J / pyunits.mol,
            doc="molar enthalpy of individual species of solid oxalate feed",
        )

        self.enth_mol_comp_product = Var(
            self.flowsheet().config.time,
            self.metal_list,
            units=pyunits.J / pyunits.mol,
            doc="molar enthalpy of individual species of solid oxide product",
        )

    def _make_mass_balance(self):
        """This section contains equations for mass balance within this model."""

        # Currently the solid feed port contains anhydrous REE oxalate
        # Convert solid inlet port anhydrous oxalate mol flow to flow_mol_comp_feed
        # Since currently only Ce is in the precipiate property package, set all others to zero
        reversed_react = dict(
            map(reversed, self.config.property_package_precipitate.react.items())
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.metal_list,
            doc="component flow of oxalates",
        )
        def flow_mol_comp_feed_eqn(b, t, i):
            return b.flow_mol_comp_feed[t, i] == pyunits.convert(
                b.solid_in[t].flow_mol_comp[reversed_react[i]],
                to_units=pyunits.mol / pyunits.second,
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package_gas.component_list,
            doc="component flow of outlet gas stream",
        )
        def flow_mol_outlet_eqn(b, t, i):
            if i == "H2O":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + sum(
                        b.flow_mol_comp_feed[t, j] * b.xH2O_oxalate_list_all[j]
                        for j in b.metal_list
                    )
                    + b.flow_mol_moist_feed[t]
                )
            elif i == "O2":
                return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[
                    i
                ] - 1.5 * sum(b.flow_mol_comp_feed[t, j] for j in b.metal_list)
            elif i == "CO2":
                return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[
                    i
                ] + 6 * sum(b.flow_mol_comp_feed[t, j] for j in b.metal_list)
            else:
                return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[i]

        @self.Constraint(
            self.flowsheet().config.time,
            self.metal_list,
            doc="mole flow rate of product",
        )
        def flow_mol_product(b, t, i):
            return (
                b.flow_mol_comp_product[t, i]
                == b.flow_mol_comp_feed[t, i] * b.frac_comp_recovery[t, i]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.metal_list,
            doc="mole flow rate of product",
        )
        def flow_mol_dust(b, t, i):
            return b.flow_mol_comp_dust[t, i] == b.flow_mol_comp_feed[t, i] * (
                1 - b.frac_comp_recovery[t, i]
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="total mass flow rate of recovered product",
        )
        def flow_mas_product(b, t):
            return sum(
                b.flow_mol_comp_product[t, i] * b.mw_oxide_list_all[i]
                for i in b.metal_list
            )

        @self.Expression(
            self.flowsheet().config.time, doc="total mass flow rate of dust product"
        )
        def flow_mas_dust(b, t):
            return sum(
                b.flow_mol_comp_dust[t, i] * b.mw_oxide_list_all[i]
                for i in b.metal_list
            )

        @self.Expression(
            self.flowsheet().config.time,
            self.metal_list,
            doc="component mass fraction of metal oxide in recovered product",
        )
        def mass_frac_comp_product(b, t, i):
            return (
                b.flow_mol_comp_product[t, i]
                * b.mw_oxide_list_all[i]
                / b.flow_mas_product[t]
            )

    def _make_energy_balance(self):
        # molar enthalpy of feed solid
        @self.Constraint(
            self.flowsheet().config.time,
            self.metal_list,
            doc="molar enthalpy of individual oxalate",
        )
        def enth_mol_comp_feed_eqn(b, t, i):
            temp = b.solid_in[t].temperature
            return b.enth_mol_comp_feed[t, i] == b.enth0_oxalate_list_all[
                i
            ] + b.cp0_oxalate_list_all[i] * (
                temp - b.temp_ref
            ) + 0.5 * b.cp1_oxalate_list_all[
                i
            ] * (
                temp * temp - b.temp_ref * b.temp_ref
            )

        # molar enthalpy of moisture in feed solid, liquid water at 25 C is 68.3174 kcal/mol
        @self.Expression(
            self.flowsheet().config.time,
            doc="molar enthalpy of moisture in solid feed stream",
        )
        def enth_mol_moist_feed(b, t):
            return (
                -b.enth_mol_const * 68317.4 * 4.184
                + b.cp_mas_const
                * 4182
                * b.mw_H2O
                * (b.solid_in[t].temperature - b.temp_ref)
            )

        # molar enthalpy of product solid
        @self.Constraint(
            self.flowsheet().config.time,
            self.metal_list,
            doc="molar enthalpy of individual oxide product",
        )
        def enth_mol_comp_product_eqn(b, t, i):
            temp_prod = b.gas_out[t].temperature
            return b.enth_mol_comp_product[t, i] == b.enth0_oxide_list_all[
                i
            ] + b.cp0_oxide_list_all[i] * (
                temp_prod - b.temp_ref
            ) + 0.5 * b.cp1_oxide_list_all[
                i
            ] * (
                temp_prod * temp_prod - b.temp_ref * b.temp_ref
            )

        # enthalpy in + heat in == enthalpy out
        @self.Constraint(
            self.flowsheet().config.time,
            doc="enthalpy balance equation for both phases",
        )
        def energy_balance_eqn(b, t):
            if self.config.has_heat_transfer is True:
                heat = b.heat_duty[t]
            else:
                heat = 0
            return (
                sum(
                    b.enth_mol_comp_feed[t, i] * b.flow_mol_comp_feed[t, i]
                    for i in b.metal_list
                )
                + b.flow_mol_moist_feed[t] * b.enth_mol_moist_feed[t]
                + b.gas_in[t].flow_mol * b.gas_in[t].enth_mol
                + heat
                == sum(
                    b.enth_mol_comp_product[t, i]
                    * (b.flow_mol_comp_product[t, i] + b.flow_mol_comp_dust[t, i])
                    for i in b.metal_list
                )
                + b.gas_out[t].flow_mol * b.gas_out[t].enth_mol
            )

    def _make_momentum_balance(self):
        @self.Constraint(self.flowsheet().config.time, doc="momentum balance equation")
        def momentum_balance_eqn(b, t):
            if self.config.has_pressure_change is True:
                return b.gas_out[t].pressure == b.gas_in[t].pressure
            else:
                return b.gas_out[t].pressure == b.gas_in[t].pressure + b.deltaP[t]

    def set_initial_condition(self):
        pass

    def initialize_build(
        blk,
        state_args_gas_in=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine.
        1.- initialize state blocks, using an initial guess for inlet
        gas inlet.
        2.- guess gas outlet component molar flowrates,
        Temperature, and Pressure. Initialize flue gas state block.
        3.- Then, solve complete model.

        Keyword Arguments:
            state_args_gas_in : a dict of arguments to be passed to the property
                           package(s) for the inlet gas state block to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize inlet property blocks
        blk.gas_in.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_gas_in
        )
        init_log.info_high("Initialization Step 1 Complete.")

        # initialize outlet gas property block
        blk.gas_out.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_gas_in
        )
        init_log.info_high("Initialization Step 2 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # set a default platen heat scaling factor
        if self.config.has_heat_transfer is True:
            for v in self.heat_duty.values():
                if iscale.get_scaling_factor(v, warning=True) is None:
                    iscale.set_scaling_factor(v, 1e-6)

        # set energy balance constraint scaling factor
        if self.config.has_heat_transfer is True:
            for t, c in self.energy_balance_eqn.items():
                sf = iscale.get_scaling_factor(
                    self.heat_duty[t], default=1e-6, warning=True
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)
