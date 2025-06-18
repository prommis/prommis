#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
IDAES REE Feed Roaster Unit Model
==================================

This model represents a roaster/calcination unit for Rare Earth Element (REE) feedstock, which includes rare earth minerals, gangue/impurity minerals, moisture, and combustible organic materials.

Reactions
---------

The reactions of impurities involved are listed below:

1. Kaolinite calcination:
   
   .. ce::
      Al2O3 \cdot 2SiO2 \cdot 2H2O -> Al2O3 + 2 SiO2 + 2 H2O (g)

   Typically occurs above 400°C, assuming complete conversion.

2. Limestone calcination:

   .. ce::
      CaCO3 -> CaO + CO2 (g)

   Typically occurs above 850°C; conversion is a user input.

3. Pyrite combustion:

   .. ce::
      FeS2 + 2.75 O2 -> 0.5 Fe2O3 + 2 SO2 (g)

   Typically occurs above 600°C, assuming complete conversion.

Combustion of organic elements is modeled as follows:

- :ce:`C + O2 -> CO2`
- :ce:`H + 0.25 O2 -> 0.5 H2O`
- :ce:`O -> 0.5 O2`
- :ce:`N -> 0.5 N2`
- :ce:`S + O2 -> SO2`

Physical Changes
----------------

Moisture in the feed stream is vaporized.

Composition
-----------

Impurity minerals are assumed to be a mixture of :ce:`Un2O3`, :ce:`CaCO3`, :ce:`SiO2`, :ce:`Al2O3`, kaolinite, and pyrite, where :ce:`Un` is an unknown element with the same atomic mass as :ce:`Al`.

The feed stream also contains organic material including :ce:`C`, :ce:`H`, :ce:`O`, :ce:`N`, :ce:`S` elements. The composition of the organic material is specified by the user.

Heat Source
-----------

The heat to the reactor can be provided either by external heating as a user input or by the combustion of a fossil fuel with air to form a hot :ce:`O2`-containing flue gas. The gas inlet stream is an :ce:`O2`-containing hot flue gas.

Streams
-------

- **Gas Inlet Stream**: :ce:`O2`-containing hot flue gas.
- **Gas Outlet Stream**: Gas product leaving the reactor.
- **Solid Outlet Stream**: Recovered solid product leaving the reactor.


Thermal Properties
------------------

The standard heats of formation and heat capacities of solid components involved are defined as parameters in this model. The default values of those parameters are obtained from two sources as listed below:

1. NIST Chemistry WebBook
2. Wagman, D.D., W.H. Evans, V.B. Parker, R.H.Schumm, I. Halow, S.M. Bailey, K.L. Churney,
   R.L. Nuttall, "The NBS tables of chemical thermodynamic properties-Selected values for
   inorganic and C1 and C2 organic substances in SI units," Journal of Physical and Chemical
   Reference Data, 11(2), 1982
3. Merrick, D., "Mathematical models of the thermal decomposition of coal, 2. Specific heats and heats of reaction,"
   Fuel, 62, pp540-546, 1983

The NIST WebBook data are used for the properties of :ce:`Al2O3`, :ce:`SiO2`, :ce:`CaO`, :ce:`Fe2O3`, and `pyrite`. Note that the heat capacity model is simplified as a linear function of temperature.
The data of Wagman et al are used for the properties of :ce:`CaCO3` and `kaolinite`.
The gas phase properties are calculated based on user configured property package.
The heat capacity of organic part of the feed is usually a function of temperature and elemental composition of C, H, O, N, and S elements according to Merrick (1983).
For simplicity, a constant heat capacity of 1260 J/kg-K in the range reported by Merrick is used in this model.

Assumptions
-----------

- No kinetics or mass transfer is considered for the calcination of impurity minerals.
- User specifies the conversion of limestone.
- Calcination of kaolinite and combustion of pyrite are assumed to be complete.
- Conversion of insoluble REE mineral to dissolvable mineral for each element is a user input.
- Final solid product is split to a recovered product stream and a dust stream with user-specified recovery fractions for the impurity and individual RE elements.
- Rare earth minerals, being in ppm level, are ignored in energy balance.
- Material balance of the REE is considered for the element only; the forms of the RE compounds (in salt or oxide forms) are not considered.
- If the product temperature is specified as a user input, the heat duty will be calculated. If the heat duty is given, the product temperature will be calculated.
- Temperatures of solid and gas products are assumed to be the same.
- No port for the solid inlet stream is used. The mass flow rate and composition of the solid reactant are specified as input variables inside the model. The mass flow rate and the composition of the solid product and dust streams are also declared as model variables. When mapping the solid products to the ``solid_outlet`` port, only the components defined in the :mod:`prommis.leaching.leach_solids_properties` module are mapped. The other species are discarded.

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
from idaes.core.util.constants import Constants as const

__author__ = "Jinliang Ma"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------------------------------------
@declare_process_block_class("REEFeedRoaster")
class REEFeedRoasterData(UnitModelBlockData):
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
        "gas_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "gas_property_package_args",
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
    CONFIG.declare(
        "solid_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "solid_property_package_args",
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
        "ree_list",
        ConfigValue(
            domain=list,
            description="List of RE elements in solid feed stream",
            doc="""A dict of the components of interest in the mixture.
        Keys are component names and values are configuration arguments to
        be passed to Component on construction.
        """,
        ),
    )

    def build(self):
        # Call TranslatorData build to setup dynamics
        super(REEFeedRoasterData, self).build()

        # Build Holdup Block
        # gas phase inlet stream
        self.gas_in = self.config.gas_property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.gas_property_package_args
        )
        # gas phase outlet stream
        self.gas_out = self.config.gas_property_package.build_state_block(
            self.flowsheet().time, **self.config.gas_property_package_args
        )
        # solid phase product stream
        self.solid_out = self.config.solid_property_package.build_state_block(
            self.flowsheet().time, **self.config.solid_property_package_args
        )

        self.add_port("gas_inlet", self.gas_in)
        self.add_port("gas_outlet", self.gas_out)
        self.add_port("solid_outlet", self.solid_out)

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

        # list of all possible REE
        self.ree_list_all = [
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

        # list of considered impurity mineral species in solid feed stream
        self.impurity_list = ["Un2O3", "CaCO3", "SiO2", "Al2O3", "Kaolinite", "Pyrite"]

        # list of organic elements
        self.organic_ele_list = ["C", "H", "O", "N", "S"]

        # list of solid product species
        self.product_list = ["Un2O3", "CaCO3", "CaO", "SiO2", "Al2O3", "Fe2O3"]

        # list of impurity elements
        self.impurity_ele_list = ["Un", "Ca", "Si", "Al", "Fe"]

        # user-specified list of rare earth elements
        self.ree_list = Set(
            initialize=self.config.ree_list,
            doc="List of rare earth elements in solid feed stream",
        )

        # atomic mass of elements involved in reactions in kg/mol
        self.am_H = Param(initialize=0.0010078, units=pyunits.kg / pyunits.mol)
        self.am_C = Param(initialize=0.012011, units=pyunits.kg / pyunits.mol)
        self.am_N = Param(initialize=0.014007, units=pyunits.kg / pyunits.mol)
        self.am_O = Param(initialize=0.015999, units=pyunits.kg / pyunits.mol)
        self.am_S = Param(initialize=0.032065, units=pyunits.kg / pyunits.mol)
        self.am_Un = Param(initialize=0.026982, units=pyunits.kg / pyunits.mol)
        self.am_Al = Param(initialize=0.026982, units=pyunits.kg / pyunits.mol)
        self.am_Si = Param(initialize=0.028086, units=pyunits.kg / pyunits.mol)
        self.am_Ca = Param(initialize=0.040078, units=pyunits.kg / pyunits.mol)
        self.am_Fe = Param(initialize=0.055845, units=pyunits.kg / pyunits.mol)

        # atomic mass of organic element
        self.am_comp_organic_ele = Param(
            self.organic_ele_list, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.am_comp_organic_ele["C"] = self.am_C
        self.am_comp_organic_ele["H"] = self.am_H
        self.am_comp_organic_ele["O"] = self.am_O
        self.am_comp_organic_ele["N"] = self.am_N
        self.am_comp_organic_ele["S"] = self.am_S

        # atomic mass of REE
        self.am_ree_list_all = Param(
            self.ree_list_all, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.am_ree_list_all["Sc"] = 0.044956
        self.am_ree_list_all["Y"] = 0.088906
        self.am_ree_list_all["La"] = 0.13891
        self.am_ree_list_all["Ce"] = 0.14012
        self.am_ree_list_all["Pr"] = 0.14091
        self.am_ree_list_all["Nd"] = 0.14424
        self.am_ree_list_all["Pm"] = 0.145
        self.am_ree_list_all["Sm"] = 0.15036
        self.am_ree_list_all["Eu"] = 0.15196
        self.am_ree_list_all["Gd"] = 0.15725
        self.am_ree_list_all["Tb"] = 0.15893
        self.am_ree_list_all["Dy"] = 0.1625
        self.am_ree_list_all["Ho"] = 0.16493
        self.am_ree_list_all["Er"] = 0.16726
        self.am_ree_list_all["Tm"] = 0.16893
        self.am_ree_list_all["Yb"] = 0.17304
        self.am_ree_list_all["Lu"] = 0.17497
        self.am_ree_list_all["Th"] = 0.23204

        # molecular weights of compounds
        self.mw_N2 = Param(initialize=self.am_N * 2, units=pyunits.kg / pyunits.mol)
        self.mw_O2 = Param(initialize=self.am_O * 2, units=pyunits.kg / pyunits.mol)
        self.mw_H2O = Param(
            initialize=self.am_H * 2 + self.am_O, units=pyunits.kg / pyunits.mol
        )
        self.mw_CO2 = Param(
            initialize=self.am_C + self.am_O * 2, units=pyunits.kg / pyunits.mol
        )
        self.mw_SO2 = Param(
            initialize=self.am_S + self.am_O * 2, units=pyunits.kg / pyunits.mol
        )
        self.mw_Un2O3 = Param(
            initialize=self.am_Un * 2 + self.am_O * 3, units=pyunits.kg / pyunits.mol
        )
        self.mw_CaO = Param(
            initialize=self.am_Ca + self.am_O, units=pyunits.kg / pyunits.mol
        )
        self.mw_SiO2 = Param(
            initialize=self.am_Si + self.am_O * 2, units=pyunits.kg / pyunits.mol
        )
        self.mw_Al2O3 = Param(
            initialize=self.am_Al * 2 + self.am_O * 3, units=pyunits.kg / pyunits.mol
        )
        self.mw_Fe2O3 = Param(
            initialize=self.am_Fe * 2 + self.am_O * 3, units=pyunits.kg / pyunits.mol
        )
        self.mw_CaCO3 = Param(
            initialize=self.am_Ca + self.am_C + self.am_O * 3,
            units=pyunits.kg / pyunits.mol,
        )
        self.mw_Kaolinite = Param(
            initialize=self.mw_Al2O3 + self.mw_SiO2 * 2 + self.mw_H2O * 2,
            units=pyunits.kg / pyunits.mol,
        )
        self.mw_Pyrite = Param(
            initialize=self.am_Fe + self.am_S * 2, units=pyunits.kg / pyunits.mol
        )

        # molecular weights of impurity mineral species in solid feed
        self.mw_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.mw_comp_impurity["Un2O3"] = self.mw_Un2O3
        self.mw_comp_impurity["CaCO3"] = self.mw_CaCO3
        self.mw_comp_impurity["SiO2"] = self.mw_SiO2
        self.mw_comp_impurity["Al2O3"] = self.mw_Al2O3
        self.mw_comp_impurity["Kaolinite"] = self.mw_Kaolinite
        self.mw_comp_impurity["Pyrite"] = self.mw_Pyrite

        # molecular weights of solid product species
        self.mw_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.mw_comp_product["Un2O3"] = self.mw_Un2O3
        self.mw_comp_product["CaCO3"] = self.mw_CaCO3
        self.mw_comp_product["CaO"] = self.mw_CaO
        self.mw_comp_product["Al2O3"] = self.mw_Al2O3
        self.mw_comp_product["SiO2"] = self.mw_SiO2
        self.mw_comp_product["Fe2O3"] = self.mw_Fe2O3

        # standard enthalpy of formation compound at 298.15 K
        self.enth0_Un2O3 = Param(initialize=-1675700, units=pyunits.J / pyunits.mol)
        self.enth0_CaO = Param(initialize=-635090, units=pyunits.J / pyunits.mol)
        self.enth0_SiO2 = Param(initialize=-910700, units=pyunits.J / pyunits.mol)
        self.enth0_Al2O3 = Param(initialize=-1675700, units=pyunits.J / pyunits.mol)
        self.enth0_Fe2O3 = Param(initialize=-825500, units=pyunits.J / pyunits.mol)
        self.enth0_CaCO3 = Param(initialize=-1206920, units=pyunits.J / pyunits.mol)
        self.enth0_Pyrite = Param(initialize=-171544, units=pyunits.J / pyunits.mol)
        self.enth0_Kaolinite = Param(initialize=-4119600, units=pyunits.J / pyunits.mol)

        # heat capacity modeled as linear function of temperature T in K.
        # Cp = Cp0 + Cp1*T
        self.cp0_Un2O3 = Param(
            initialize=28.039, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_CaO = Param(
            initialize=47.563, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_SiO2 = Param(
            initialize=31.414, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_Al2O3 = Param(
            initialize=28.039, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_Fe2O3 = Param(
            initialize=80.623, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_CaCO3 = Param(
            initialize=76.009, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_Pyrite = Param(
            initialize=59.09, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_Kaolinite = Param(
            initialize=246.14, units=pyunits.J / pyunits.mol / pyunits.K
        )

        self.cp1_Un2O3 = Param(
            initialize=0.17156, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_CaO = Param(
            initialize=0.0059258, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_SiO2 = Param(
            initialize=0.054286, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_Al2O3 = Param(
            initialize=0.17156, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_Fe2O3 = Param(
            initialize=0.09936, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_CaCO3 = Param(
            initialize=0.046296, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_Pyrite = Param(
            initialize=0.024232, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_Kaolinite = Param(
            initialize=0.0, units=pyunits.J / pyunits.mol / pyunits.K**2
        )

        # standard heat of formation of impurity species
        self.enth0_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.J / pyunits.mol
        )
        self.enth0_comp_impurity["Un2O3"] = self.enth0_Un2O3
        self.enth0_comp_impurity["CaCO3"] = self.enth0_CaCO3
        self.enth0_comp_impurity["SiO2"] = self.enth0_SiO2
        self.enth0_comp_impurity["Al2O3"] = self.enth0_Al2O3
        self.enth0_comp_impurity["Kaolinite"] = self.enth0_Kaolinite
        self.enth0_comp_impurity["Pyrite"] = self.enth0_Pyrite

        # standard heat of formation of solid product species
        self.enth0_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.J / pyunits.mol
        )
        self.enth0_comp_product["Un2O3"] = self.enth0_Un2O3
        self.enth0_comp_product["CaCO3"] = self.enth0_CaCO3
        self.enth0_comp_product["CaO"] = self.enth0_CaO
        self.enth0_comp_product["Al2O3"] = self.enth0_Al2O3
        self.enth0_comp_product["SiO2"] = self.enth0_SiO2
        self.enth0_comp_product["Fe2O3"] = self.enth0_Fe2O3

        # 1st heat capacity coefficient of impurity species
        self.cp0_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_comp_impurity["Un2O3"] = self.cp0_Un2O3
        self.cp0_comp_impurity["CaCO3"] = self.cp0_CaCO3
        self.cp0_comp_impurity["SiO2"] = self.cp0_SiO2
        self.cp0_comp_impurity["Al2O3"] = self.cp0_Al2O3
        self.cp0_comp_impurity["Kaolinite"] = self.cp0_Kaolinite
        self.cp0_comp_impurity["Pyrite"] = self.cp0_Pyrite

        # 1st heat capacity coefficient of solid product species
        self.cp0_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_comp_product["Un2O3"] = self.cp0_Un2O3
        self.cp0_comp_product["CaCO3"] = self.cp0_CaCO3
        self.cp0_comp_product["CaO"] = self.cp0_CaO
        self.cp0_comp_product["Al2O3"] = self.cp0_Al2O3
        self.cp0_comp_product["SiO2"] = self.cp0_SiO2
        self.cp0_comp_product["Fe2O3"] = self.cp0_Fe2O3

        # 2nd heat capacity coefficient of impurity species
        self.cp1_comp_impurity = Param(
            self.impurity_list,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        self.cp1_comp_impurity["Un2O3"] = self.cp1_Un2O3
        self.cp1_comp_impurity["CaCO3"] = self.cp1_CaCO3
        self.cp1_comp_impurity["SiO2"] = self.cp1_SiO2
        self.cp1_comp_impurity["Al2O3"] = self.cp1_Al2O3
        self.cp1_comp_impurity["Kaolinite"] = self.cp1_Kaolinite
        self.cp1_comp_impurity["Pyrite"] = self.cp1_Pyrite

        # 2nd heat capacity coefficient of solid product species
        self.cp1_comp_product = Param(
            self.product_list,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        self.cp1_comp_product["Un2O3"] = self.cp1_Un2O3
        self.cp1_comp_product["CaCO3"] = self.cp1_CaCO3
        self.cp1_comp_product["CaO"] = self.cp1_CaO
        self.cp1_comp_product["Al2O3"] = self.cp1_Al2O3
        self.cp1_comp_product["SiO2"] = self.cp1_SiO2
        self.cp1_comp_product["Fe2O3"] = self.cp1_Fe2O3

        # organic heat capacity, currently assuming constant based on Merrick (1983)
        self.cp_organic = Param(
            initialize=1260, units=pyunits.J / pyunits.kg / pyunits.K
        )

        # unit constants used for the expressions of HHV and heat of formation of organic material and liquid water enthalpy
        self.enth_mas_const = Param(
            initialize=1,
            units=pyunits.J / pyunits.kg,
            doc="1 unit of mass enthalpy in J/kg",
        )
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

        self.flow_mass_feed = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=pyunits.kg / pyunits.s,
            doc="total mass flow rate of solid feed including surface moisture and combustible organic material",
        )

        self.mass_frac_moist_feed = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            doc="mass fraction of moisture in solid feed",
        )

        self.mass_frac_organic_feed = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            doc="mass fraction of total organic materials in solid feed",
        )

        self.mass_frac_comp_impurity_feed = Var(
            self.flowsheet().config.time,
            self.impurity_list,
            initialize=0.2,
            doc="mass fraction of impurity minerals in solid feed stream on dry organic free basis",
        )

        self.mass_frac_comp_organic_feed = Var(
            self.flowsheet().config.time,
            self.organic_ele_list,
            initialize=0.25,
            doc="mass fraction of organic elements in solid feed stream on dry ash free basis",
        )

        self.ppm_comp_ree_ins_feed = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=10,
            doc="fraction of insoluble rare earth element in ppm per unit mass of impurity minerals in solid feed",
        )

        self.ppm_comp_ree_dis_feed = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=10,
            doc="mass fraction of dissovable rare earth element in ppm per unit mass of impurity minerals in solid feed",
        )

        self.xconv_comp_ins = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=0.9,
            doc="fraction of conversion of insoluble REE to dissovable REE",
        )

        self.frac_comp_ree_recovery = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=0.95,
            doc="fraction of recovery of each REE in the product stream",
        )

        self.frac_impurity_recovery = Var(
            self.flowsheet().config.time,
            initialize=0.95,
            doc="fractional of recovery of impurity product in the product stream",
        )

        # solid feed temperature
        self.temp_feed = Var(
            self.flowsheet().config.time,
            initialize=298.15,
            units=pyunits.K,
            doc="solid feed temperature",
        )

        self.xconv_caco3 = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            doc="fractional conversion of CaCO3 calcination reaction",
        )

        self.flow_mol_comp_product_total = Var(
            self.flowsheet().config.time,
            self.product_list,
            initialize=10,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of total solid product from calcination of impurity minerals",
        )

        self.flow_mol_comp_product_recovered = Var(
            self.flowsheet().config.time,
            self.product_list,
            initialize=10,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of recovered solid product from calcination of impurity minerals",
        )

        self.flow_mol_comp_product_dust = Var(
            self.flowsheet().config.time,
            self.product_list,
            initialize=1,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of solid product in dust stream from calcination of impurity minerals",
        )

        self.flow_mass_product_recovered = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=pyunits.kg / pyunits.s,
            doc="mass flow rate of recovered solid product from calcination of impurity minerals",
        )

        self.flow_mass_product_dust = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            units=pyunits.kg / pyunits.s,
            doc="mass flow rate of solid product in dust stream from calcination of impurity minerals",
        )

        self.mass_frac_comp_impurity_ele_product = Var(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            initialize=0.1,
            doc="mass fraction of impurity element in solid product stream",
        )

        self.ppm_comp_ree_ins_product = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=10,
            doc="fraction of insoluble rare earth element in product stream",
        )

        self.ppm_comp_ree_dis_product = Var(
            self.flowsheet().config.time,
            self.ree_list,
            initialize=10,
            doc="fraction of dissovable rare earth element in product stream",
        )

    def _make_mass_balance(self):
        """This section contains equations for mass balance within the reactor model."""

        # Vaporized moisture mole flow rate
        @self.Expression(
            self.flowsheet().config.time, doc="vaporized moisture mass flowrate [kg/s]"
        )
        def flow_mol_moist_feed(b, t):
            return b.flow_mass_feed[t] * b.mass_frac_moist_feed[t] / b.mw_H2O

        # Organic material mass flow rate
        @self.Expression(
            self.flowsheet().config.time, doc="Organic material mass flowrate [kg/s]"
        )
        def flow_mass_organic_feed(b, t):
            return b.flow_mass_feed[t] * b.mass_frac_organic_feed[t]

        # Organic element mole flow rate
        @self.Expression(
            self.flowsheet().config.time,
            self.organic_ele_list,
            doc="Organic material mass flowrate [kg/s]",
        )
        def flow_mol_comp_organic_feed(b, t, i):
            return (
                b.flow_mass_organic_feed[t]
                * b.mass_frac_comp_organic_feed[t, i]
                / b.am_comp_organic_ele[i]
            )

        # impurity mineral mass flow rate
        @self.Expression(
            self.flowsheet().config.time,
            doc="Impurity mineral material mass flowrate [kg/s]",
        )
        def flow_mass_impurity_feed(b, t):
            return b.flow_mass_feed[t] * (
                1 - b.mass_frac_moist_feed[t] - b.mass_frac_organic_feed[t]
            )

        # mass flow rates of individual mineral impurities
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_list,
            doc="mass flow rate of individual impurity minerals [kg/s]",
        )
        def flow_mass_comp_impurity_feed(b, t, i):
            return b.flow_mass_impurity_feed[t] * b.mass_frac_comp_impurity_feed[t, i]

        # mole flow rates of individual mineral impurity
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_list,
            doc="mole flow rate of individual impurity minerals [mol/s]",
        )
        def flow_mol_comp_impurity_feed(b, t, i):
            return b.flow_mass_comp_impurity_feed[t, i] / b.mw_comp_impurity[i]

        # the gas product should contain at least N2, O2, H2O, CO2, and SO2
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.gas_property_package.component_list,
            doc="component flow of outlet gas stream",
        )
        def flow_mol_outlet_eqn(b, t, i):
            if i == "N2":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + b.flow_mol_comp_organic_feed[t, "N"] / 2
                )
            elif i == "O2":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + b.flow_mol_comp_organic_feed[t, "O"] / 2
                    - b.flow_mol_comp_organic_feed[t, "C"]
                    - b.flow_mol_comp_organic_feed[t, "S"]
                    - b.flow_mol_comp_organic_feed[t, "H"] / 4
                    - b.flow_mol_comp_impurity_feed[t, "Pyrite"] * 2.75
                )
            elif i == "H2O":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + b.flow_mol_comp_organic_feed[t, "H"] / 2
                    + b.flow_mol_comp_impurity_feed[t, "Kaolinite"] * 2
                    + b.flow_mol_moist_feed[t]
                )
            elif i == "CO2":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + b.flow_mol_comp_organic_feed[t, "C"]
                    + b.flow_mol_comp_impurity_feed[t, "CaCO3"] * b.xconv_caco3[t]
                )
            elif i == "SO2":
                return (
                    b.gas_out[t].flow_mol_comp[i]
                    == b.gas_in[t].flow_mol_comp[i]
                    + b.flow_mol_comp_organic_feed[t, "S"]
                    + b.flow_mol_comp_impurity_feed[t, "Pyrite"] * 2
                )
            else:
                return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[i]

        # total impurity solid product based on calcination reactions
        @self.Constraint(
            self.flowsheet().config.time,
            self.product_list,
            doc="mole flow rate of total product before dust carry over",
        )
        def flow_mol_comp_product_total_eqn(b, t, i):
            if i == "CaCO3":
                return b.flow_mol_comp_product_total[
                    t, i
                ] == b.flow_mol_comp_impurity_feed[t, "CaCO3"] * (1 - b.xconv_caco3[t])
            elif i == "CaO":
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, "CaCO3"] * b.xconv_caco3[t]
                )
            elif i == "SiO2":
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, "SiO2"]
                    + b.flow_mol_comp_impurity_feed[t, "Kaolinite"] * 2
                )
            elif i == "Al2O3":
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, "Al2O3"]
                    + b.flow_mol_comp_impurity_feed[t, "Kaolinite"]
                )
            elif i == "Fe2O3":
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, "Pyrite"] / 2
                )
            else:  # i == "Un2O3"
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, "Un2O3"]
                )

        # mole flow rate of recovered solid product stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.product_list,
            doc="mole flow rate of recovered product stream",
        )
        def flow_mol_comp_product_recovered_eqn(b, t, i):
            return (
                b.flow_mol_comp_product_recovered[t, i]
                == b.flow_mol_comp_product_total[t, i] * b.frac_impurity_recovery[t]
            )

        # mole flow rate of solid product in dust stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.product_list,
            doc="mole flow rate of recovered product stream",
        )
        def flow_mol_comp_product_dust_eqn(b, t, i):
            return b.flow_mol_comp_product_dust[t, i] == b.flow_mol_comp_product_total[
                t, i
            ] * (1 - b.frac_impurity_recovery[t])

        # mass flow rate of solid product recovered
        @self.Constraint(self.flowsheet().config.time, doc="mass flow rate of product")
        def flow_mass_product_recovered_eqn(b, t):
            return b.flow_mass_product_recovered[t] == sum(
                b.flow_mol_comp_product_recovered[t, i] * b.mw_comp_product[i]
                for i in b.product_list
            )

        # mass flow rate of solid product in dust stream
        @self.Constraint(self.flowsheet().config.time, doc="mass flow rate of product")
        def flow_mass_product_dust_eqn(b, t):
            return b.flow_mass_product_dust[t] == sum(
                b.flow_mol_comp_product_dust[t, i] * b.mw_comp_product[i]
                for i in b.product_list
            )

        # mass fraction of impurity metal in the solid feed stream on dry organic free basis
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            doc="mass fraction of impurity element in dry organic free solid feed stream",
        )
        def mass_frac_comp_impurity_ele_feed(b, t, i):
            if i == "Ca":
                return (
                    b.flow_mol_comp_impurity_feed[t, "CaCO3"]
                    * b.am_Ca
                    / b.flow_mass_impurity_feed[t]
                )
            elif i == "Si":
                return (
                    (
                        b.flow_mol_comp_impurity_feed[t, "SiO2"]
                        + b.flow_mol_comp_impurity_feed[t, "Kaolinite"] * 2
                    )
                    * b.am_Si
                    / b.flow_mass_impurity_feed[t]
                )
            elif i == "Al":
                return (
                    (
                        b.flow_mol_comp_impurity_feed[t, "Al2O3"]
                        + b.flow_mol_comp_impurity_feed[t, "Kaolinite"]
                    )
                    * b.am_Al
                    * 2
                    / b.flow_mass_impurity_feed[t]
                )
            elif i == "Fe":
                return (
                    b.flow_mol_comp_impurity_feed[t, "Pyrite"]
                    * b.am_Fe
                    / b.flow_mass_impurity_feed[t]
                )
            else:  # i==Un
                return (
                    b.flow_mol_comp_impurity_feed[t, "Un2O3"]
                    * b.am_Un
                    * 2
                    / b.flow_mass_impurity_feed[t]
                )

        # mass fraction of impurity metal in the product stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            doc="mass fraction of impurity element in solid product",
        )
        def mass_frac_comp_impurity_ele_product_eqn(b, t, i):
            if i == "Ca":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == (
                        b.flow_mol_comp_product_recovered[t, "CaCO3"]
                        + b.flow_mol_comp_product_recovered[t, "CaO"]
                    )
                    * b.am_Ca
                )
            elif i == "Si":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "SiO2"] * b.am_Si
                )
            elif i == "Al":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Al2O3"] * 2 * b.am_Al
                )
            elif i == "Fe":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Fe2O3"] * 2 * b.am_Fe
                )
            else:  # i=='Un':
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Un2O3"] * 2 * b.am_Un
                )

        # insoluble solid ppm in the product stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.ree_list,
            doc="insoluble solid ppm in solid product",
        )
        def ppm_comp_ree_ins_product_eqn(b, t, i):
            return (
                b.ppm_comp_ree_ins_product[t, i] * b.flow_mass_product_recovered[t]
                == b.flow_mass_impurity_feed[t]
                * b.ppm_comp_ree_ins_feed[t, i]
                * (1 - b.xconv_comp_ins[t, i])
                * b.frac_comp_ree_recovery[t, i]
            )

        # dissovable solid ppm in the product stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.ree_list,
            doc="dissovable solid ppm in solid product",
        )
        def ppm_comp_ree_dis_product_eqn_eqn(b, t, i):
            return (
                b.ppm_comp_ree_dis_product[t, i] * b.flow_mass_product_recovered[t]
                == b.flow_mass_impurity_feed[t]
                * (
                    b.ppm_comp_ree_ins_feed[t, i] * b.xconv_comp_ins[t, i]
                    + b.ppm_comp_ree_dis_feed[t, i]
                )
                * b.frac_comp_ree_recovery[t, i]
            )

        # solid product with 13 species defined in CoalRefuseParameters
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.solid_property_package.component_list,
            doc="component flow of outlet solid stream",
        )
        def solid_outlet_comp_eqn(b, t, i):
            if i == "Al2O3" or i == "CaO" or i == "Fe2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    * pyunits.convert(
                        b.solid_out[t].flow_mass, to_units=pyunits.kg / pyunits.second
                    )
                    == b.flow_mol_comp_product_recovered[t, i]
                    * b.config.solid_property_package.mw[i]
                )
            elif i == "inerts":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    * pyunits.convert(
                        b.solid_out[t].flow_mass, to_units=pyunits.kg / pyunits.second
                    )
                    == b.flow_mol_comp_product_recovered[t, "Un2O3"] * b.mw_Un2O3
                    + b.flow_mol_comp_product_recovered[t, "SiO2"] * b.mw_SiO2
                    + b.flow_mol_comp_product_recovered[t, "CaCO3"] * b.mw_CaCO3
                )
            elif i == "Sc2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Sc"]
                        + b.ppm_comp_ree_dis_product[t, "Sc"]
                    )
                    / 1e6
                )
            elif i == "Y2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Y"]
                        + b.ppm_comp_ree_dis_product[t, "Y"]
                    )
                    / 1e6
                )
            elif i == "La2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "La"]
                        + b.ppm_comp_ree_dis_product[t, "La"]
                    )
                    / 1e6
                )
            elif i == "Ce2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Ce"]
                        + b.ppm_comp_ree_dis_product[t, "Ce"]
                    )
                    / 1e6
                )
            elif i == "Pr2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Pr"]
                        + b.ppm_comp_ree_dis_product[t, "Pr"]
                    )
                    / 1e6
                )
            elif i == "Nd2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Nd"]
                        + b.ppm_comp_ree_dis_product[t, "Nd"]
                    )
                    / 1e6
                )
            elif i == "Sm2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Sm"]
                        + b.ppm_comp_ree_dis_product[t, "Sm"]
                    )
                    / 1e6
                )
            elif i == "Gd2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Gd"]
                        + b.ppm_comp_ree_dis_product[t, "Gd"]
                    )
                    / 1e6
                )
            elif i == "Dy2O3":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    == (
                        b.ppm_comp_ree_ins_product[t, "Dy"]
                        + b.ppm_comp_ree_dis_product[t, "Dy"]
                    )
                    / 1e6
                )
            else:
                pass

    def _make_energy_balance(self):
        # HHV based on Dulong's formula
        @self.Expression(
            self.flowsheet().config.time, doc="HHV based on Dulong's formula [J/kg]"
        )
        def HHV_Dulong(b, t):
            return b.enth_mas_const * (
                232367.224
                * (
                    145.44 * b.mass_frac_comp_organic_feed[t, "C"]
                    + 620
                    * (
                        b.mass_frac_comp_organic_feed[t, "H"]
                        - b.mass_frac_comp_organic_feed[t, "O"] / 8
                    )
                    + 41 * b.mass_frac_comp_organic_feed[t, "S"]
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Standard heat of formation of organic materials [J/kg]",
        )
        def Hf_organic(b, t):
            dhcoal = -b.HHV_Dulong[t] + const.gas_constant * b.temp_ref / 2 * (
                -b.mass_frac_comp_organic_feed[t, "H"] / 2 / b.am_H
                + b.mass_frac_comp_organic_feed[t, "O"] / b.am_O
                + b.mass_frac_comp_organic_feed[t, "N"] / b.am_N
            )
            dhc = (
                -b.enth_mol_const
                * 94052
                * 4.184
                * b.mass_frac_comp_organic_feed[t, "C"]
                / b.am_C
            )
            dhh = (
                -b.enth_mol_const
                * 68317.4
                * 4.184
                * b.mass_frac_comp_organic_feed[t, "H"]
                / b.am_H
                / 2
            )
            dhs = (
                -b.enth_mol_const
                * 70940
                * 4.184
                * b.mass_frac_comp_organic_feed[t, "S"]
                / b.am_S
            )
            return dhc + dhh + dhs - dhcoal

        # molar enthalpy of impurity minerals in solid feed
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_list,
            doc="molar enthalpy of individual impurity minerals",
        )
        def enth_mol_comp_impurity_feed(b, t, i):
            return (
                b.enth0_comp_impurity[i]
                + b.cp0_comp_impurity[i] * (b.temp_feed[t] - b.temp_ref)
                + 0.5
                * b.cp1_comp_impurity[i]
                * (b.temp_feed[t] * b.temp_feed[t] - b.temp_ref * b.temp_ref)
            )

        # molar enthalpy of minerals in solid product
        @self.Expression(
            self.flowsheet().config.time,
            self.product_list,
            doc="molar enthalpy of individual product minerals",
        )
        def enth_mol_comp_product(b, t, i):
            temp_prod = b.gas_out[t].temperature
            return (
                b.enth0_comp_product[i]
                + b.cp0_comp_product[i] * (temp_prod - b.temp_ref)
                + 0.5
                * b.cp1_comp_product[i]
                * (temp_prod * temp_prod - b.temp_ref * b.temp_ref)
            )

        # molar enthalpy of moisture in feed solid, liquid water at 25 C is 68.3174 kcal/mol
        @self.Expression(
            self.flowsheet().config.time,
            doc="molar enthalpy of moisture in solid feed stream",
        )
        def enth_mol_moist_feed(b, t):
            return (
                -b.enth_mol_const * 68317.4 * 4.184
                + b.cp_mas_const * 4182 * b.mw_H2O * (b.temp_feed[t] - b.temp_ref)
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
                    b.enth_mol_comp_impurity_feed[t, i]
                    * b.flow_mol_comp_impurity_feed[t, i]
                    for i in b.impurity_list
                )
                + b.flow_mol_moist_feed[t] * b.enth_mol_moist_feed[t]
                + b.flow_mass_organic_feed[t]
                * (b.Hf_organic[t] + b.cp_organic * (b.temp_feed[t] - b.temp_ref))
                + b.gas_in[t].flow_mol * b.gas_in[t].enth_mol
                + heat
                == sum(
                    b.enth_mol_comp_product[t, i] * b.flow_mol_comp_product_total[t, i]
                    for i in b.product_list
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
        # Initialize inlet gas property block
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
