#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
IDAES REE Oxalate Roaster Unit Model
====================================

REE Oxalate Roaster is a unit operation to convert oxalates of rare earth and gangue metal
elements to metal oxides through thermal decomposition and oxidation. There are 18 rare earth elements
including `Sc`, `Y`, `La`, `Ce`, `Pr`, `Nd`, `Pm`, `Sm`, `Eu`, `Gd`, `Tb`, `Dy`, `Ho`, `Er`, `Tm`, `Yb`, `Lu`,
and `Th` in the model. The user can specify a subset of those 18 elements through configuration when creating
the unit model. The 3 gangue elements considered in the model are `Fe`, `Al`, and `Ca`.

The feed oxalate mixture stream is typically from a precipitator in an upstream process. It is assumed
that all oxalates are in their hydrate forms. In case the anhydrous oxalate feed flow rates are specified in
the property package of the solid feed stream, the molar flow rates are converted to the corresponding hydrate flow rates. 
The molecular formula of an oxalate hydrate can be expressed in a general form as :ce:`RE2(C2O4)3 \\cdot xH2O` where RE is
one of the 18 rare earth elements and x is the number of water molecules associated with the hyrate. The three 
gangue oxalate hydrates considered in the model are :ce:`Fe2(C2O4)3 \\cdot 2H2O`, :ce:`Al2(C2O4)3 \\cdot H2O`, and
:ce:`CaC2O4 \\cdot H2O`, for `Fe`, `Al`, and `Ca` elements, respectively.

The feed stream contains surface moisture. The amount of surface moisture entering the reactor is specified by
a liquid inlet that contains a fraction of the liquid outlet of the precipitator.


Physical Changes
----------------

Moisture in the solid feed stream from the liquid inlet is vaporized.


Reactions
---------

The thermal decomposition reactions modeled for rare earth oxalates and three gaugue oxalates are as
listed below:

- :ce:`RE2(C2O4)3 \\cdot xH2O(s) + 1.5O2(g) -> RE2O3(s) + 6CO2(g) + xH2O(g)`
- :ce:`Fe2(C2O4)3 \\cdot 2H2O(s) + 1.5O2(g) -> Fe2O3(s) + 6CO2(g) + 2H2O(g)`
- :ce:`Al2(C2O4)3 \\cdot H2O(s) + 1.5O2(g) -> Al2O3(s) + 6CO2(g) + H2O(g)`
- :ce:`CaC2O4 \\cdot H2O(s) + 0.5O2(g) -> CaO(s) + 2CO2(g) + H2O(g)`

In the first reaction, RE represent any of 18 rare earth elements and x is the number of water molecules
associated with the oxalate. Typically x=10 for most rare earth elements.
It is assumed in the current model that the reaction is carried out at a high enough temperature and
long enough residence time that the reactants are completely decomposed. It is also assumed that a gas
feed stream provides enough :ce:`O2` reactant to carry out the reactions. Since the kinetics of the reactions
are not modeled, the current unit model is valid for steady-state simulations only.


Thermal Properties
------------------

The standard heats of formation and heat capacities of solid components involved are defined as parameters in this model.
The default values of those parameters are obtained from three sources as listed below:

1. NIST Chemistry WebBook
2. Wagman, D.D., W.H. Evans, V.B. Parker, R.H.Schumm, I. Halow, S.M. Bailey, K.L. Churney,
   R.L. Nuttall, "The NBS tables of chemical thermodynamic properties-Selected values for
   inorganic and C1 and C2 organic substances in SI units," Journal of Physical and Chemical
   Reference Data, 11(2), 1982
3. Kotz, J.C., P.M. Treichel, J. Townsend, D. Treichel, "Chemistry and Chemical Reactivity,"
   9th Edition, Cengage Learning, 2014

The solid heat capacity model is simplified as a linear function of temperature. Since the data for the rare earth
components are very limited, default parameters based on the average values of :ce:`La`, :ce:`Ce`, :ce:`Pr`, and :ce:`Nd`
are used for the species if no reported data are found.

The gas phase properties are calculated based on user configured property package.


Mass Balance
------------

The content of the surface moisture specified by the liquid inlet stream is vaporized and enters the gas phase.
The other species in the liquid solution, including metal elements in the liquid inlet stream, are ignored.

The species mass balance is based on complete conversion of solid reactants such that the molar flow rates of
individual metals (rare earth and gangue elements) are conserved. For the species in the gas phase, the :ce:`O2`
is consumed while :ce:`CO2` and :ce:`H2O` are produced. For any other species in the gas feed stream that does not
participate in any reactions, its molar flow rate in the gas product stream is the same as that in the
inlet stream. Note that the user needs to make sure that the gas feed stream contains enough :ce:`O2` to avoid
negative flow rate of :ce:`O2` in the gas product stream.

Two solid product streams are modeled, one representing the fine solid product particles carried out by the gas
product stream and the other are the remaining solid oxides leaving the reactor that is recovered as the final
metal oxide product. The fractions of recovery for individual metal oxides are specified by the user as model inputs.


Energy Balance
--------------

The model considers the heat required to vaporize the moisture content of the solid feed based on the enthalpy
increase from the liquid water in the solid feed stream to water vapor in the gas outlet stream.
The heats of reactions are also considered along with the heat capacities of the reactants and product species.
The total enthalpy including the standard heat of formation and sensible heat of an individual component (either
a gas or a solid species) is used for energy balance calculation.
The outlet gas product stream is assumed to have the same temperature as the two solid product streams.


Heat Source
-----------

The heat to the reactor can be provided either by external heating as a user input or by the combustion of a fossil fuel
with air to form a hot :ce:`O2`-containing flue gas. The gas inlet stream is an :ce:`O2`-containing hot flue gas.


Streams
-------

- **Gas Inlet Stream**: :ce:`O2`-containing hot flue gas.
- **Gas Outlet Stream**: Gas product leaving the reactor.
- **Solid Inlet Stream**: Solid feed of oxalate mixture entering the reactor.
- **Liquid Inlet Stream**: Surface moisture in the solid feed as a liquid inlet entering the reactor.

Note that the two solid product streams (fine solid carried by the gas stream and the recovered solid product stream)
are not defined as outlet streams and their flow rates are defined as model variables.

"""

# Import Pyomo libraries
from pyomo.common.config import Bool, ConfigBlock, ConfigValue

# Additional import for the unit operation
from pyomo.environ import Param, Set, Var
from pyomo.environ import units as pyunits

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# Import IDAES cores
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.solvers import get_solver
from idaes.core.util.config import DefaultBool, is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe

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
        "property_package_precipitate_solid",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Solid precipitate property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_precipitate_solid",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing solid precipitate property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "property_package_precipitate_liquid",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Liquid precipitate property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_precipitate_liquid",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing liquid precipitate property packages",
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
                "Ca",
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
            doc="""A list of metal elements contained in solid oxalate feed
        """,
        ),
    )

    def build(self):
        # Call TranslatorData build to setup dynamics
        super(REEOxalateRoasterData, self).build()

        # Attributed for storing contents of reporting output
        self._stream_table_dict = {}

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
        self.solid_in = (
            self.config.property_package_precipitate_solid.build_state_block(
                self.flowsheet().time,
                **self.config.property_package_args_precipitate_solid,
            )
        )
        # liquid phase inlet stream from precipitator
        self.liquid_in = (
            self.config.property_package_precipitate_liquid.build_state_block(
                self.flowsheet().time,
                **self.config.property_package_args_precipitate_liquid,
            )
        )
        self.add_port("gas_inlet", self.gas_in)
        self.add_port("gas_outlet", self.gas_out)
        self.add_port("solid_inlet", self.solid_in)
        self.add_port("liquid_inlet", self.liquid_in)

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

        # List of all possible elements considered in this unit model
        self.metal_list_all = [
            "Al",
            "Fe",
            "Ca",
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
        # atomic mass of each metal element
        self.am_metal_list_all = Param(
            self.metal_list_all, mutable=True, units=pyunits.kg / pyunits.mol
        )
        # count of C2O4 in each oxalate molecule
        self.xC2O4_oxalate_list_all = Param(
            self.metal_list_all, initialize=3, mutable=True
        )
        # count of associated H2O in each hydrate
        self.xH2O_oxalate_list_all = Param(
            self.metal_list_all, initialize=10, mutable=True
        )
        # molecular weight of each oxalate molecule
        self.mw_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=0.1,
            mutable=True,
            units=pyunits.kg / pyunits.mol,
        )
        # molecular weight of each oxide molecule
        self.mw_oxide_list_all = Param(
            self.metal_list_all,
            initialize=0.1,
            mutable=True,
            units=pyunits.kg / pyunits.mol,
        )
        # atomic masses of all elements involved
        self.am_H = Param(initialize=0.0010078, units=pyunits.kg / pyunits.mol)
        self.am_C = Param(initialize=0.012011, units=pyunits.kg / pyunits.mol)
        self.am_N = Param(initialize=0.014007, units=pyunits.kg / pyunits.mol)
        self.am_O = Param(initialize=0.015999, units=pyunits.kg / pyunits.mol)
        self.am_metal_list_all["Al"] = 0.026982
        self.am_metal_list_all["Fe"] = 0.055845
        self.am_metal_list_all["Ca"] = 0.040078
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
        # The count of C2O4 unit in oxalate is 3 for all oxalate except CaC2O4_H2O
        self.xC2O4_oxalate_list_all["Ca"] = 1
        # The count of H2O in the oxalate hydrate is 10 except three listed below
        self.xH2O_oxalate_list_all["Al"] = 1
        self.xH2O_oxalate_list_all["Fe"] = 2
        self.xH2O_oxalate_list_all["Ca"] = 1

        self.mw_H2O = Param(
            initialize=self.am_H * 2 + self.am_O, units=pyunits.kg / pyunits.mol
        )

        for i in self.metal_list_all:
            if i == "Ca":
                self.mw_oxalate_list_all[i] = (
                    self.am_metal_list_all[i]
                    + self.xC2O4_oxalate_list_all[i] * (self.am_C * 2 + self.am_O * 4)
                    + self.xH2O_oxalate_list_all[i] * (self.am_H * 2 + self.am_O)
                )
            else:
                self.mw_oxalate_list_all[i] = (
                    self.am_metal_list_all[i] * 2
                    + self.xC2O4_oxalate_list_all[i] * (self.am_C * 2 + self.am_O * 4)
                    + self.xH2O_oxalate_list_all[i] * (self.am_H * 2 + self.am_O)
                )

        for i in self.metal_list_all:
            if i == "Ca":
                self.mw_oxide_list_all[i] = self.am_metal_list_all[i] + self.am_O
            else:
                self.mw_oxide_list_all[i] = (
                    self.am_metal_list_all[i] * 2 + 3 * self.am_O
                )

        # molar standard enthalpy of oxalate at 298.15 K initialized to the average of La, Ce, Pr, and Nd
        # Data from Wagman et al (1982)
        self.enth0_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=-6350044,
            mutable=True,
            units=pyunits.J / pyunits.mol,
        )
        # molar standard enthalpy of oxide at 298.15 K initialized to the average of La, Ce, Pr, and Nd
        self.enth0_oxide_list_all = Param(
            self.metal_list_all,
            initialize=-1801865.75,
            mutable=True,
            units=pyunits.J / pyunits.mol,
        )
        # heat capacity Cp = Cp0 + Cp1*T (a linear function of temperature in K)
        # constant term of molar heat capacity of oxalate initialized to the average of La, Ce, Pr, and Nd
        self.cp0_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=45.751,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        # linear coefficient of molar heat capacity of oxalate initialized to the average of La, Ce, Pr, and Nd
        self.cp1_oxalate_list_all = Param(
            self.metal_list_all,
            initialize=0,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        # constant molar heat capacity of oxide initialized to the average of La, Ce, Pr, and Nd
        self.cp0_oxide_list_all = Param(
            self.metal_list_all,
            initialize=110.3625,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        # linear coefficient of molar heat capacity of oxide initialized to the average of La, Ce, Pr, and Nd
        self.cp1_oxide_list_all = Param(
            self.metal_list_all,
            initialize=0.033887,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        # Oxalate standard enthalpy available in literature
        self.enth0_oxalate_list_all["Al"] = (
            -3397000
        )  # Kotz et al (2014), anhydrous data
        self.enth0_oxalate_list_all["Fe"] = (
            -2572300
        )  # Wagman et al (1982), anhydrous data
        self.enth0_oxalate_list_all["Ca"] = -1674860  # Wagman et al (1982)
        self.enth0_oxalate_list_all["La"] = -5916176  # Not in Wagman et al
        self.enth0_oxalate_list_all["Ce"] = -6782000  # Wagman et al (1982)
        self.enth0_oxalate_list_all["Pr"] = -5920000  # Wagman et al (1982)
        self.enth0_oxalate_list_all["Nd"] = -6782000  # Wagman et al (1982)
        # Oxide standard enthalpy available in literature
        self.enth0_oxide_list_all["Fe"] = -825500  # NIST WebBook
        self.enth0_oxide_list_all["Al"] = -1675700  # NIST WebBook
        self.enth0_oxide_list_all["Ca"] = -635090  # NIST WebBook
        self.enth0_oxide_list_all["La"] = -1793702  # Wagman et al (1982)
        self.enth0_oxide_list_all["Ce"] = -1796191  # Wagman et al (1982)
        self.enth0_oxide_list_all["Pr"] = -1809664  # Wagman et al (1982)
        self.enth0_oxide_list_all["Nd"] = -1807906  # Wagman et al (1982)
        self.enth0_oxide_list_all["Sc"] = -1908820  # Wagman et al (1982)
        self.enth0_oxide_list_all["Y"] = -1905310  # Wagman et al (1982)
        # Heat capacity of most oxalates except Ca are unavailable, use the default value
        self.cp0_oxalate_list_all["Ca"] = 152.8  # Wagman et al (1982)

        # Heat capacity of oxide available in literature
        self.cp0_oxide_list_all["La"] = 107.72  # revised based on Wagman et al (1982)
        self.cp0_oxide_list_all["Ce"] = 115.78  # revised based on Wagman et al (1982)
        self.cp0_oxide_list_all["Pr"] = 112.82  # revised based on Wagman et al (1982)
        self.cp0_oxide_list_all["Nd"] = 105.13  # revised based on Wagman et al (1982)
        self.cp0_oxide_list_all["Al"] = 28.039  # NIST WebBook
        self.cp0_oxide_list_all["Fe"] = 80.623  # NIST WebBook
        self.cp0_oxide_list_all["Ca"] = 47.2  # NIST WebBook

        self.cp1_oxide_list_all["La"] = 0.026114  # revised based on Wagman et al (1982)
        self.cp1_oxide_list_all["Ce"] = 0.03477  # revised based on Wagman et al (1982)
        self.cp1_oxide_list_all["Pr"] = 0.034364  # revised based on Wagman et al (1982)
        self.cp1_oxide_list_all["Nd"] = 0.0403  # revised based on Wagman et al (1982)
        self.cp1_oxide_list_all["Al"] = 0.17156  # NIST WebBook
        self.cp1_oxide_list_all["Fe"] = 0.09936  # NIST WebBook
        self.cp1_oxide_list_all["Ca"] = 0.00299  # NIST WebBook

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
            doc="mole flow rate of liquid water as surface moisture from liquid inlet",
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

        self.flow_mass_product = Var(
            self.flowsheet().config.time,
            bounds=(0.0, None),
            units=pyunits.kg / pyunits.s,
            doc="mass flow rate of total oxides in product stream",
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
        reversed_react = dict(
            map(reversed, self.config.property_package_precipitate_solid.react.items())
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
            doc="surface moisture molar flow rate",
        )
        def flow_mol_moist_feed_eqn(b, t):
            return b.flow_mol_moist_feed[t] == pyunits.convert(
                b.liquid_in[t].flow_mol_comp["H2O"],
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
                ] - sum(
                    b.xC2O4_oxalate_list_all[j] / 2 * b.flow_mol_comp_feed[t, j]
                    for j in b.metal_list
                )
            elif i == "CO2":
                return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[
                    i
                ] + sum(
                    b.xC2O4_oxalate_list_all[j] * 2 * b.flow_mol_comp_feed[t, j]
                    for j in b.metal_list
                )
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

        @self.Constraint(
            self.flowsheet().config.time,
            doc="total mass flow rate of recovered product",
        )
        def flow_mass_product_eqn(b, t):
            return b.flow_mass_product[t] == sum(
                b.flow_mol_comp_product[t, i] * b.mw_oxide_list_all[i]
                for i in b.metal_list
            )

        @self.Expression(
            self.flowsheet().config.time, doc="total mass flow rate of dust product"
        )
        def flow_mass_dust(b, t):
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
                / b.flow_mass_product[t]
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
        # since the liquid inlet does not contain temperature, use the solid inlet temperature
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

        # molar enthalpy of solid product
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

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            self._stream_table_dict, time_point=time_point
        )

    def _get_performance_contents(self, time_point=0):
        exprs = {}

        for j in self.config.metal_list:
            exprs[f"Product {j} Mass Fraction"] = self.mass_frac_comp_product[
                time_point, j
            ]

        return {"exprs": exprs}
