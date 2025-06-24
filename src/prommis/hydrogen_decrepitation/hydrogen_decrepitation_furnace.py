#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Hydrogen Decrepitation Unit Model
==================================

This model represents a gas-fired or electrically-heated hydrogen decrepitation unit for Rare Earth Permanent Magnet (REPM) feedstocks, specifically recycling hard drive disks (HDDs) containing neodymium compounds.

The process and reactions modeled are described by

Alireza Habibzadeh, Mehmet Ali Kucuker, and Mertol Gökelma
Review on the Parameters of Recycling NdFeB Magnets via a Hydrogenation Process
ACS Omega 2023 8 (20), 17431-17445
DOI: 10.1021/acsomega.3c00299

The hydrogen decrepitation and dehydrogenation process consists of the two reactions below:

1. Hydrogenation of the Nd2Fe14B matrix:

   .. ce::
      Nd2Fe14B (s) + H2 (g) -> NdFe14B \cdot H2 (powder)

   Typically occurs between 25-400°C and 0.5-15000 kPa.


2. Dehydrogenation of the powder product:

   .. ce::
      NdFe14B \cdot H2 (powder) -> Nd2Fe14B (powder) + H2 (g)

   Occurs during downstream processing at high temperatures, e.g. heating, leaching, sintering. This model assumes that this reaction occurs inside the hydrogen decrepitation unit due to heat released during the process.

Therefore, the net reaction in this process is:
   
   .. ce::
      Nd2Fe14B (s) + H2 (g) + \Delta H -> Nd2Fe14B (powder) + H2 (g)

Two additional reactions can occur if an Nd-rich phase exists:

3. Hydrogenation of Nd-rich phase:
   
   .. ce::
      Nd (s) + H2 (g) -> Nd \cdot H2 (powder)

   Typically occurs between 25-400°C and 0.5-15000 kPa.


4. Dehydrogenation of the Nd-rich phase:

   .. ce::
      Nd \cdot H2 (powder) -> Nd (powder) + H2 (g)

   Occurs during downstream processing at high temperatures, e.g. heating, leaching, sintering.



Heat Source
-----------

The heat to the reactor is provided by external heating, either gas-fired or electric.

Streams
-------

- **Gas Inlet Stream**: :ce:`H2` gas entering the reactor.
- **Gas Outlet Stream**: `H2` gas leaving the reactor.
- **Solid Inlet Stream**: REPM solid sample entering the reactor.
- **Solid Outlet Stream**: REPM powder leaving the reactor.


Thermal Properties
------------------

The properties of components involved are defined as parameters in this model. The default values of those parameters are obtained from two sources as listed below:

1. NIST Chemistry WebBook
2. Neodymium (Nd) (chemicalaid.com) https://periodictable.chemicalaid.com/element.php/Nd
3. Properties of the rare earth metals and compounds. Gibson, John A; Harvey, Gifford S. Ohio State University Research Foundation Columbus. https://apps.dtic.mil/sti/trecms/pdf/AD0479391.pdf
4. Masao Morishita, Taichi Abe, Ai Nozaki, Ikuo Ohnuma, Keisuke Kamon, Calorimetric study of Nd2Fe14B: Heat capacity, standard Gibbs energy of formation and magnetic entropy, Thermochimica Acta, Volume 690, 2020, 178672, ISSN 0040-6031, https://doi.org/10.1016/j.tca.2020.178672.

The data from NIST WebBook provides the atomic masses of elements present in the process.
The data from Chemical Aid provides the standard enthalpy of formation of Nd.
The data from Gibson and Harvey provides molar heat capacities of Nd and Nd2Fe14B as linear functions of temperature.
The data from Morishita, et al. provides the standard enthalpy of formation of Nd2Fe14B.
The data from Biegański provides the molar heat capacity of NdH2.
The gas phase properties are calculated based on user configured property package.


"""

# Import Pyomo libraries
from pyomo.common.config import Bool, ConfigBlock, ConfigValue

# Additional import for the unit operation
from pyomo.environ import Param, PositiveReals, Set, Var, exp
from pyomo.environ import units as pyunits

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# Import IDAES cores
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.solvers import get_solver
from idaes.core.util.config import DefaultBool, is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.util.math import smooth_max
from idaes.core.util.tables import create_stream_table_dataframe

__author__ = "Brandon Paul, Akintomiwa Ojo"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------------------------------------
@declare_process_block_class("REPMHydrogenDecrepitationFurnace")
class REPMHydrogenDecrepitationFurnaceData(UnitModelBlockData):
    """
    Simple 0D hydrogen decrepitation furnace model with mass and energy balance only
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
    CONFIG.declare(
        "number_of_units",
        ConfigValue(
            default=1,
            domain=int,
            description="Integer value of number of furnaces required",
            doc="""Applies to several model equations to account for multiple units.
        """,
        ),
    )

    def build(self):
        # Call TranslatorData build to setup dynamics
        super(REPMHydrogenDecrepitationFurnaceData, self).build()

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
        # solid phase inlet stream
        self.solid_in = self.config.solid_property_package.build_state_block(
            self.flowsheet().time, **self.config.solid_property_package_args
        )
        # solid phase product stream
        self.solid_out = self.config.solid_property_package.build_state_block(
            self.flowsheet().time, **self.config.solid_property_package_args
        )

        self.add_port("gas_inlet", self.gas_in)
        self.add_port("gas_outlet", self.gas_out)
        self.add_port("solid_inlet", self.solid_in)
        self.add_port("solid_outlet", self.solid_out)

        # Add Geometry
        if self.config.has_holdup is True:
            self.volume = Var(
                initialize=1, units=pyunits.m**3, doc="volume of the reactor"
            )
            self.volume.fix()

        # Add heat duty
        if self.config.has_heat_transfer is True:
            self.supplied_heat_duty = Var(
                self.flowsheet().time,
                initialize=0,
                units=pyunits.W,
                doc="heat duty added to the reactor in addition to that required by decrepitation",
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
        self._make_geometry()
        self._make_mass_balance()
        self._make_energy_balance()
        self._make_momentum_balance()

    def _make_params(self):
        """This section is for parameters within this model."""

        # reference temperature at standard condition
        self.temp_ref = Param(initialize=298.15, units=pyunits.K)

        # list of all possible REE
        self.ree_list_all = [
            "Nd",
        ]

        # list of considered impurity mineral species in solid feed stream
        self.impurity_list = [
            "Nd",
            "Nd2Fe14B",
        ]

        # list of solid product species
        self.product_list = [
            "Nd",
            "Nd2Fe14B",
        ]

        # list of impurity elements
        self.impurity_ele_list = [
            "Nd",
            "Fe",
            "B",
        ]

        # user-specified list of rare earth elements
        self.ree_list = Set(
            initialize=self.config.ree_list,
            doc="List of rare earth elements in solid feed stream",
        )

        # atomic mass of elements involved in reactions in kg/mol [1]
        self.am_H = Param(initialize=0.0010078, units=pyunits.kg / pyunits.mol)
        self.am_Nd = Param(initialize=0.14424, units=pyunits.kg / pyunits.mol)
        self.am_Fe = Param(initialize=0.055845, units=pyunits.kg / pyunits.mol)
        self.am_B = Param(initialize=0.010811, units=pyunits.kg / pyunits.mol)

        # atomic mass of REE
        self.am_ree_list_all = Param(
            self.ree_list_all, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.am_ree_list_all["Nd"] = self.am_Nd

        # molecular weights of compounds
        self.mw_Nd2Fe14B = Param(
            initialize=self.am_Nd * 2 + self.am_Fe * 14 + self.am_B,
            units=pyunits.kg / pyunits.mol,
        )

        # molecular weights of impurity mineral species in solid feed
        self.mw_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.mw_comp_impurity["Nd"] = self.am_Nd
        self.mw_comp_impurity["Nd2Fe14B"] = self.mw_Nd2Fe14B

        # molecular weights of solid product species
        self.mw_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.kg / pyunits.mol
        )
        self.mw_comp_product["Nd"] = self.am_Nd
        self.mw_comp_product["Nd2Fe14B"] = self.mw_Nd2Fe14B

        # standard enthalpy of formation compound at 298.15 K
        self.enth0_Nd = Param(initialize=326.9e3, units=pyunits.J / pyunits.mol)  # [2]
        self.enth0_Nd2Fe14B = Param(
            initialize=-93.39e3, units=pyunits.J / pyunits.mol
        )  # [4]

        # heat capacity modeled as linear function of temperature T in K.
        # Cp = Cp0 + Cp1*T
        self.cp0_Nd = Param(
            initialize=22.44, units=pyunits.J / pyunits.mol / pyunits.K
        )  # [3]
        self.cp0_Nd2Fe14B = Param(
            initialize=272.59, units=pyunits.J / pyunits.mol / pyunits.K
        )

        self.cp1_Nd = Param(
            initialize=0.0162, units=pyunits.J / pyunits.mol / pyunits.K**2
        )
        self.cp1_Nd2Fe14B = Param(
            initialize=0.6487, units=pyunits.J / pyunits.mol / pyunits.K**2
        )

        # standard heat of formation of impurity species
        self.enth0_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.J / pyunits.mol
        )
        self.enth0_comp_impurity["Nd"] = self.enth0_Nd
        self.enth0_comp_impurity["Nd2Fe14B"] = self.enth0_Nd2Fe14B

        # standard heat of formation of solid product species
        self.enth0_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.J / pyunits.mol
        )
        self.enth0_comp_product["Nd"] = self.enth0_Nd
        self.enth0_comp_product["Nd2Fe14B"] = self.enth0_Nd2Fe14B

        # 1st heat capacity coefficient of impurity species
        self.cp0_comp_impurity = Param(
            self.impurity_list, mutable=True, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_comp_impurity["Nd"] = self.cp0_Nd
        self.cp0_comp_impurity["Nd2Fe14B"] = self.cp0_Nd2Fe14B

        # 1st heat capacity coefficient of solid product species
        self.cp0_comp_product = Param(
            self.product_list, mutable=True, units=pyunits.J / pyunits.mol / pyunits.K
        )
        self.cp0_comp_product["Nd"] = self.cp0_Nd
        self.cp0_comp_product["Nd2Fe14B"] = self.cp0_Nd2Fe14B

        # 2nd heat capacity coefficient of impurity species
        self.cp1_comp_impurity = Param(
            self.impurity_list,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        self.cp1_comp_impurity["Nd"] = self.cp1_Nd
        self.cp1_comp_impurity["Nd2Fe14B"] = self.cp1_Nd2Fe14B

        # 2nd heat capacity coefficient of solid product species
        self.cp1_comp_product = Param(
            self.product_list,
            mutable=True,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )
        self.cp1_comp_product["Nd"] = self.cp1_Nd
        self.cp1_comp_product["Nd2Fe14B"] = self.cp1_Nd2Fe14B

        # operating parameters
        self.ramp_up_time = Param(
            initialize=300,
            units=pyunits.s,
            mutable=True,
            doc="Time required to reach operating temperature",
        )
        self.decrepitation_duration = Param(
            initialize=10800,
            units=pyunits.s,
            mutable=True,
            doc="Amount of time the furnace operates at its operating temperature",
        )
        self.preparation_time = Param(
            initialize=3600,
            units=pyunits.s,
            mutable=True,
            doc="Time required to setup the furnace",
        )
        self.cool_down_time = Param(
            initialize=3600,
            units=pyunits.s,
            mutable=True,
            doc="Time required for the decrepitated sample to cool down to room temperature",
        )
        self.sample_density = Param(
            initialize=7500,
            units=pyunits.kg / pyunits.m**3,
            doc="Mass density of the sample",
        )
        self.chamber_to_sample_ratio = Param(
            initialize=2,
            units=pyunits.dimensionless,
            mutable=True,
            doc="ratio of the volume of furnace chamber to the volume of the sample",
        )

        # vendor material properties for first insulation, ceramide fiber
        # https://www.grainger.com/product/23AR50?gucid=N:N:PS:Paid:GGL:CSM-2295:4P7A1P:20501231&gad_source=1&gclid=CjwKCAjwnqK1BhBvEiwAi7o0Xxz1EhsrD9nk0yk9SgHwN7fzjC8X249MA5NulaF_JSsTjC290z-HChoCSyYQAvD_BwE&gclsrc=aw.ds
        self.length_insulation1 = Param(
            initialize=7.62,
            units=pyunits.m,
            doc="The length of a unit of the first insulation material",
        )
        self.width_insulation1 = Param(
            initialize=0.6096,
            units=pyunits.m,
            doc="The width of a unit of the first insulation material",
        )
        self.thickness_insulation1 = Param(
            initialize=0.0254,
            units=pyunits.m,
            doc="The thickness of a unit of the first insulation material",
        )
        self.weight_insulation1 = Param(
            initialize=15.42,
            units=pyunits.kg,
            doc="The weight of a unit of the first insulation material",
        )
        self.thermal_cond_insulation_material1 = Param(
            initialize=0.33,
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="The thermal conductivity of a unit of the first insulation material",
        )
        self.specific_heat_capacity_insulation1 = Param(
            initialize=1.08,
            units=pyunits.kJ / pyunits.kg / pyunits.K,
            doc="The specific heat capacity of a unit of the first insulation material",
        )

        # vendor material properties for first metal material, SS304
        self.density_metal1 = Param(
            initialize=7473.57,
            units=pyunits.kg / pyunits.m**3,
            doc="The density of the first metal material",
        )
        self.thermal_cond_metal_material1 = Param(
            initialize=13.53,
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="The thermal conductivity of the first metal material",
        )
        self.specific_heat_capacity_metal1 = Param(
            initialize=0.468,
            units=pyunits.kJ / pyunits.kg / pyunits.K,
            doc="The specific heat capacity of the first metal material",
        )

        # vendor material properties for second insulation, fiberglass
        self.length_insulation2 = Param(
            initialize=1.19,
            units=pyunits.m,
            doc="The length of a unit of the second insulation material",
        )
        self.width_insulation2 = Param(
            initialize=0.381,
            units=pyunits.m,
            doc="The width of a unit of the second insulation material",
        )
        self.thickness_insulation2 = Param(
            initialize=0.0889,
            units=pyunits.m,
            doc="The thickness of a unit of the second insulation material",
        )
        self.weight_insulation2 = Param(
            initialize=10.43,
            units=pyunits.kg,
            doc="The weight of a unit of the second insulation material",
        )
        self.thermal_cond_insulation_material2 = Param(
            initialize=0.069,
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="The thermal conductivity of a unit of the second insulation material",
        )
        self.specific_heat_capacity_insulation2 = Param(
            initialize=0.9,
            units=pyunits.kJ / pyunits.kg / pyunits.K,
            doc="The specific heat capacity of a unit of the second insulation material",
        )

        # vendor material properties for second metal material, carbon steel
        self.density_metal2 = Param(
            initialize=7861.09,
            units=pyunits.kg / pyunits.m**3,
            doc="The density of the second metal material",
        )
        self.thermal_cond_metal_material2 = Param(
            initialize=45,
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="The thermal conductivity of the second metal material",
        )
        self.specific_heat_capacity_metal2 = Param(
            initialize=0.502416,
            units=pyunits.kJ / pyunits.kg / pyunits.K,
            doc="The specific heat capacity of the second metal material",
        )

        # heat transfer parameters
        self.max_temperature = Param(
            initialize=1173.15,
            units=pyunits.K,
            doc="Maximum allowable operating temperature",
        )
        self.temperature_metal_material1 = Param(
            initialize=950.13,
            units=pyunits.K,
            doc="Temperature of metal material 1 under steady operating conditions",
        )
        self.temperature_insulation_material2 = Param(
            initialize=950.00,
            units=pyunits.K,
            doc="Temperature of insulation material 2 under steady operating conditions",
        )
        self.temperature_metal_material2 = Param(
            initialize=333.18,
            units=pyunits.K,
            doc="Temperature of metal material 2 under steady operating conditions",
        )
        self.temperature_furnace_ext_surface = Param(
            initialize=333.15,
            units=pyunits.K,
            doc="Temperature of external surface of the furnace under steady operating conditions",
        )
        self.air_heat_transfer_coeff = Param(
            initialize=5,
            units=pyunits.W / (pyunits.m**2 * pyunits.K),
            doc="Heat transfer coefficient of air (natural/free convection) at ambient temperature",
        )
        self.ref_temp = Param(
            initialize=298.15,
            units=pyunits.K,
            doc="Reference temperature - 25 degree celsius",
        )
        self.min_quantity_insulation = Param(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Smooth max parameter for relative volume of insulation materials",
        )
        self.insulation1_max_delta_temp = Param(
            initialize=875,
            units=pyunits.K,
            mutable=True,
            doc="Temperature change of insulation material 1 from reference temperature to steady operating condition",
        )
        self.metal1_max_delta_temp = Param(
            initialize=651.98,
            units=pyunits.K,
            mutable=True,
            doc="Temperature change of metal 1 from reference temperature to steady operating condition",
        )
        self.insulation2_max_delta_temp = Param(
            initialize=651.85,
            units=pyunits.K,
            mutable=True,
            doc="Temperature change of insulation material 2 from reference temperature to steady operating condition",
        )
        self.metal2_max_delta_temp = Param(
            initialize=35.03,
            units=pyunits.K,
            mutable=True,
            doc="Temperature change of metal 2 from reference temperature to steady operating condition",
        )

    def _make_vars(self):
        """This section declares variables within this model."""

        self.flow_vol_feed = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=pyunits.m**3 / pyunits.s,
            doc="solid feed volumetric flow",
        )

        self.temp_feed = Var(
            self.flowsheet().config.time,
            initialize=298.15,
            units=pyunits.K,
            doc="solid feed temperature",
        )

        self.temp_prod = Var(
            self.flowsheet().config.time,
            initialize=298.15,
            units=pyunits.K,
            doc="solid product temperature after cooling",
        )

        self.operating_temperature = Var(
            self.flowsheet().config.time,
            initialize=443.15,
            units=pyunits.K,
            doc="operating temperature of furnace",
        )

        self.flow_mol_comp_product_total = Var(
            self.flowsheet().config.time,
            self.product_list,
            initialize=10,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of total solid product from decrepitation of impurity minerals",
        )

        self.flow_mol_comp_product_recovered = Var(
            self.flowsheet().config.time,
            self.product_list,
            initialize=10,
            units=pyunits.mol / pyunits.s,
            doc="mole flow rate of recovered solid product from decrepitation of impurity minerals",
        )

        self.flow_mass_product_recovered = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=pyunits.kg / pyunits.s,
            doc="mass flow rate of recovered solid product from calcination of impurity minerals",
        )

        self.mass_frac_comp_impurity_ele_product = Var(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            initialize=0.1,
            doc="mass fraction of impurity element in solid product stream",
        )

        # operating variables
        self.sample_heat_capacity = Var(
            self.flowsheet().config.time,
            initialize=0.44,
            units=pyunits.kJ / (pyunits.kg * pyunits.K),
            doc="heat capacity of the sample",
        )

        self.sample_mass = Var(
            self.flowsheet().config.time,
            initialize=62.025,
            units=pyunits.kg,
            doc="mass of the sample",
        )

        self.sample_volume = Var(
            self.flowsheet().config.time,
            initialize=0.002,
            units=pyunits.m**3,
            doc="volume of the sample",
        )

        self.furnace_chamber_volume = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=1,
            bounds=(1e-12, 1e12),
            units=pyunits.m**3,
            doc="internal volume of the furnace",
        )

        self.heat_loss = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=1e3,
            units=pyunits.W,
            bounds=(1e-6, None),
            doc="Heat loss from the external surface of the furnace",
        )

        self.temperature_insulation_material1 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=1050,
            units=pyunits.K,
            bounds=(950.13, 1173.15),
            doc="Temperature of insulation material 1 under steady operating conditions",
        )

        self.thickness_insulation_material1 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=0.12,
            units=pyunits.m,
            bounds=(0, None),
            doc="Thickness of insulation material 1",
        )

        self.relative_thickness_ratio_insulation_material1 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of insulation material 1",
        )

        self.thickness_metal_material1 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=3 * 1e-3,
            units=pyunits.m,
            bounds=(0, None),
            doc="Thickness of metal material 1",
        )

        self.relative_thickness_ratio_metal_material1 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=0.05,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of metal material 1",
        )

        self.thickness_insulation_material2 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=0.15,
            units=pyunits.m,
            bounds=(0, None),
            doc="Thickness of insulation material 2",
        )

        self.relative_thickness_ratio_insulation_material2 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=1,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of insulation material 2",
        )

        self.thickness_metal_material2 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=4 * 1e-3,
            units=pyunits.m,
            bounds=(0, None),
            doc="Thickness of metal material 2",
        )

        self.relative_thickness_ratio_metal_material2 = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=0.01,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of metal material 2",
        )

        self.total_heat_duty = Var(
            self.flowsheet().config.time,
            within=PositiveReals,
            initialize=3500,
            units=pyunits.W,
            bounds=(0, None),
            doc="Total heat duty required for decrepitation",
        )

    def _make_geometry(self):
        """This section contains equations related to the geometry of the reactor model."""

        # feed properties
        @self.Constraint(
            self.flowsheet().config.time,
            doc="volumetric flow of solid feed",
        )
        def flow_vol_feed_eqn(b, t):
            return b.flow_vol_feed[t] == pyunits.convert(
                b.solid_in[t].flow_mass / b.sample_density,
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="mass of solid feed",
        )
        def sample_mass_eqn(b, t):
            return b.sample_mass[t] == pyunits.convert(
                b.solid_in[t].flow_mass * b.decrepitation_duration, to_units=pyunits.kg
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="volume of solid feed",
        )
        def sample_volume_eqn(b, t):
            return b.sample_volume[t] == pyunits.convert(
                b.flow_vol_feed[t] * b.decrepitation_duration, to_units=pyunits.m**3
            )

        # furnace chamber dimensions
        @self.Constraint(
            self.flowsheet().config.time,
            doc="internal volume of the furnace",
        )
        def furnace_chamber_volume_eqn(b, t):
            return b.furnace_chamber_volume[t] == pyunits.convert(
                b.chamber_to_sample_ratio * b.sample_volume[t], to_units=pyunits.m**3
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Radius of the furnace chamber"
        )
        def radius_chamber(b, t):
            return (b.furnace_chamber_volume[t] / (6 * Constants.pi)) ** (1 / 3)

        @self.Expression(
            self.flowsheet().config.time, doc="Length of the furnace chamber"
        )
        def length_chamber(b, t):
            # 3 times the diameter
            return 6 * b.radius_chamber[t]

        @self.Expression(
            self.flowsheet().config.time, doc="Volume of insulation material 1"
        )
        def volume_insulation1(b, t):
            return (
                Constants.pi
                * b.length_chamber[t]
                * (
                    ((b.radius_chamber[t] + b.thickness_insulation_material1[t]) ** 2)
                    - (b.radius_chamber[t] ** 2)
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Required quantity of insulation material 1",
        )
        def quantity_insulation1(b, t):
            return smooth_max(
                b.min_quantity_insulation,
                (
                    b.volume_insulation1[t]
                    / (
                        b.length_insulation1
                        * b.width_insulation1
                        * b.thickness_insulation1
                    )
                ),
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Total weight of insulation material 1"
        )
        def total_weight_insulation1(b, t):
            return b.quantity_insulation1[t] * b.weight_insulation1

        @self.Expression(
            self.flowsheet().config.time,
            doc="Internal diameter of the structure formed after attaching metal 1",
        )
        def internal_diameter_metal1(b, t):
            return 2 * (b.radius_chamber[t] + b.thickness_insulation_material1[t])

        @self.Expression(
            self.flowsheet().config.time,
            doc="External diameter of the structure formed after attaching metal 1",
        )
        def external_diameter_metal1(b, t):
            return 2 * (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Required volume of metal material 1"
        )
        def volume_metal1(b, t):
            return pyunits.convert(
                (
                    Constants.pi
                    * (
                        (b.external_diameter_metal1[t] ** 2)
                        - (b.internal_diameter_metal1[t] ** 2)
                    )
                    * b.length_chamber[t]
                )
                / 4,
                to_units=pyunits.inches**3,
            )

        @self.Expression(self.flowsheet().config.time, doc="Weight of metal material 1")
        def weight_metal1(b, t):
            return b.volume_metal1[t] * pyunits.convert(
                b.density_metal1, to_units=pyunits.pound / (pyunits.inches**3)
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Volume of insulation material 2"
        )
        def volume_insulation2(b, t):
            return (
                Constants.pi
                * b.length_chamber[t]
                * (
                    (
                        (
                            b.radius_chamber[t]
                            + b.thickness_insulation_material1[t]
                            + b.thickness_metal_material1[t]
                            + b.thickness_insulation_material2[t]
                        )
                        ** 2
                    )
                    - (
                        (
                            b.radius_chamber[t]
                            + b.thickness_insulation_material1[t]
                            + b.thickness_metal_material1[t]
                        )
                        ** 2
                    )
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Required quantity of insulation material 2",
        )
        def quantity_insulation2(b, t):
            return smooth_max(
                b.min_quantity_insulation,
                (
                    b.volume_insulation2[t]
                    / (
                        b.length_insulation2
                        * b.width_insulation2
                        * b.thickness_insulation2
                    )
                ),
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Total weight of insulation material 2"
        )
        def total_weight_insulation2(b, t):
            return b.quantity_insulation2[t] * b.weight_insulation2

        @self.Expression(
            self.flowsheet().config.time,
            doc="Internal diameter of the structure formed after attaching metal 2",
        )
        def internal_diameter_metal2(b, t):
            return 2 * (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
                + b.thickness_insulation_material2[t]
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="External diameter of the structure formed after attaching metal 2",
        )
        def external_diameter_metal2(b, t):
            return 2 * (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
                + b.thickness_insulation_material2[t]
                + b.thickness_metal_material2[t]
            )

        @self.Expression(
            self.flowsheet().config.time, doc="Required volume of metal material 2"
        )
        def volume_metal2(b, t):
            return pyunits.convert(
                (
                    Constants.pi
                    * (
                        (b.external_diameter_metal2[t] ** 2)
                        - (b.internal_diameter_metal2[t] ** 2)
                    )
                    * b.length_chamber[t]
                )
                / 4,
                to_units=pyunits.inches**3,
            )

        @self.Expression(self.flowsheet().config.time, doc="Weight of metal material 2")
        def weight_metal2(b, t):
            return b.volume_metal2[t] * pyunits.convert(
                b.density_metal2, to_units=pyunits.pound / (pyunits.inches**3)
            )

        @self.Expression(
            self.flowsheet().config.time, doc="External surface area of the furnace"
        )
        def furnace_external_surface_area(b, t):
            return (
                Constants.pi * b.external_diameter_metal2[t] * b.length_chamber[t]
            ) + ((Constants.pi * (b.external_diameter_metal2[t] ** 2)) / 2)

        @self.Expression(
            self.flowsheet().config.time,
            doc="Energy required to raise the temperature of furnace material from room temperature to temperature at steady state",
        )
        def heat_furnace_material(b, t):
            return (
                (
                    b.total_weight_insulation1[t]
                    * b.specific_heat_capacity_insulation1
                    * b.insulation1_max_delta_temp
                )
                + (
                    (pyunits.convert(b.weight_metal1[t], to_units=pyunits.kg))
                    * b.specific_heat_capacity_metal1
                    * b.metal1_max_delta_temp
                )
                + (
                    b.total_weight_insulation2[t]
                    * b.specific_heat_capacity_insulation2
                    * b.insulation2_max_delta_temp
                )
                + (
                    (pyunits.convert(b.weight_metal2[t], to_units=pyunits.kg))
                    * b.specific_heat_capacity_metal2
                    * b.metal2_max_delta_temp
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Energy lost in the decrepitation duration",
        )
        def energy_consumption(b, t):
            return pyunits.convert(
                b.heat_loss[t] * b.decrepitation_duration,
                to_units=pyunits.kJ,
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Energy required to raise the temperature of sample to operating temperature",
        )
        def heat_sample_material(b, t):
            return (
                b.sample_mass[t]
                * b.sample_heat_capacity[t]
                * (b.operating_temperature[t] - b.ref_temp)
            )

        @self.Expression(self.flowsheet().config.time, doc="Total heat duty needed")
        def total_heat(b, t):
            return (
                b.heat_furnace_material[t]
                + b.energy_consumption[t]
                + b.heat_sample_material[t]
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Total duration of the decrepitation process",
        )
        def processing_time(b, t):
            return (
                (pyunits.convert(b.preparation_time, to_units=pyunits.hr))
                + (pyunits.convert(b.ramp_up_time, to_units=pyunits.hr))
                + (pyunits.convert(b.decrepitation_duration, to_units=pyunits.hr))
                + (pyunits.convert(b.cool_down_time, to_units=pyunits.hr))
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Total heat duty required for decrepitation",
        )
        def total_heat_duty_constraint(b, t):
            return b.total_heat_duty[t] == pyunits.convert(
                b.total_heat[t] / b.processing_time[t], to_units=pyunits.W
            )

    def _make_mass_balance(self):
        """This section contains equations for mass balance within the reactor model."""

        # impurity mineral mass flow rate
        @self.Expression(
            self.flowsheet().config.time,
            doc="Impurity mineral material mass flowrate [kg/s]",
        )
        def flow_mass_impurity_feed(b, t):
            return sum(
                b.solid_in[t].flow_mass * b.solid_in[t].mass_frac_comp[i]
                for i in b.impurity_list
            )

        # mass flow rates of individual mineral impurities
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_list,
            doc="mass flow rate of individual impurity minerals [kg/s]",
        )
        def flow_mass_comp_impurity_feed(b, t, i):
            return b.solid_in[t].flow_mass * b.solid_in[t].mass_frac_comp[i]

        # mole flow rates of individual mineral impurity
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_list,
            doc="mole flow rate of individual impurity minerals [mol/s]",
        )
        def flow_mol_comp_impurity_feed(b, t, i):
            return b.flow_mass_comp_impurity_feed[t, i] / b.mw_comp_impurity[i]

        # the gas product is H2, no other species are carried out with it
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.gas_property_package.component_list,
            doc="component flow of outlet gas stream",
        )
        def flow_mol_outlet_eqn(b, t, i):
            return b.gas_out[t].flow_mol_comp[i] == b.gas_in[t].flow_mol_comp[i]

        # total impurity solid product based on decrepitation reactions
        @self.Constraint(
            self.flowsheet().config.time,
            self.product_list,
            doc="mole flow rate of total product",
        )
        def flow_mol_comp_product_total_eqn(b, t, i):
            if i == "Nd" or i == "Nd2Fe14B":
                return (
                    b.flow_mol_comp_product_total[t, i]
                    == b.flow_mol_comp_impurity_feed[t, i]
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
                == b.flow_mol_comp_product_total[t, i]
            )

        # mass flow rate of solid product recovered
        @self.Constraint(self.flowsheet().config.time, doc="mass flow rate of product")
        def flow_mass_product_recovered_eqn(b, t):
            return b.flow_mass_product_recovered[t] == sum(
                b.flow_mol_comp_product_recovered[t, i] * b.mw_comp_product[i]
                for i in b.product_list
            )

        # mass fraction of impurity metal in the solid feed stream
        @self.Expression(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            doc="mass fraction of impurity element in solid feed stream",
        )
        def mass_frac_comp_impurity_ele_feed(b, t, i):
            if i == "Fe":
                return (
                    b.flow_mol_comp_impurity_feed[t, "Nd2Fe14B"]
                    * b.am_Fe
                    * 14
                    / b.flow_mass_impurity_feed[t]
                )

        # mass fraction of impurity metal in the product stream
        @self.Constraint(
            self.flowsheet().config.time,
            self.impurity_ele_list,
            doc="mass fraction of impurity element in solid product",
        )
        def mass_frac_comp_impurity_ele_product_eqn(b, t, i):
            if i == "Nd":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Nd2Fe14B"] * b.am_Nd * 2
                )
            elif i == "Fe":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Nd2Fe14B"] * b.am_Fe * 14
                )
            elif i == "B":
                return (
                    b.mass_frac_comp_impurity_ele_product[t, i]
                    * b.flow_mass_product_recovered[t]
                    == b.flow_mol_comp_product_recovered[t, "Nd2Fe14B"] * b.am_B
                )

        # solid product with 13 species defined in CoalRefuseParameters
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.solid_property_package.component_list,
            doc="component flow of outlet solid stream",
        )
        def solid_outlet_comp_eqn(b, t, i):
            if i == "Nd" or i == "Nd2Fe14B":
                return (
                    b.solid_out[t].mass_frac_comp[i]
                    * pyunits.convert(
                        b.solid_out[t].flow_mass, to_units=pyunits.kg / pyunits.second
                    )
                    == b.flow_mol_comp_product_recovered[t, i]
                    * b.config.solid_property_package.mw[i]
                )
            else:
                pass

    def _make_energy_balance(self):
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
            return (
                b.enth0_comp_product[i]
                + b.cp0_comp_product[i] * (b.temp_prod[t] - b.temp_ref)
                + 0.5
                * b.cp1_comp_product[i]
                * (b.temp_prod[t] * b.temp_prod[t] - b.temp_ref * b.temp_ref)
            )

        # enthalpy in + heat in == enthalpy out
        @self.Constraint(
            self.flowsheet().config.time,
            doc="enthalpy balance equation for both phases",
        )
        def energy_balance_eqn(b, t):
            if self.config.has_heat_transfer is True:
                heat = b.supplied_heat_duty[t]
            else:
                heat = 0
            return (
                sum(
                    b.enth_mol_comp_impurity_feed[t, i]
                    * b.flow_mol_comp_impurity_feed[t, i]
                    for i in b.impurity_list
                )
                + b.gas_in[t].flow_mol * b.gas_in[t].enth_mol
                + heat
                == sum(
                    b.enth_mol_comp_product[t, i] * b.flow_mol_comp_product_total[t, i]
                    for i in b.product_list
                )
                + b.gas_out[t].flow_mol * b.gas_out[t].enth_mol
            )

        # sample heat capacity
        @self.Constraint(
            self.flowsheet().config.time,
            doc="solid sample heat capacity",
        )
        def sample_heat_capacity_eqn(b, t):
            # lever rule, Cp = sum(mass_frac_i * Cp_mol_i * mw_i for all components i)
            return b.sample_heat_capacity[t] == pyunits.convert(
                b.solid_in[t].mass_frac_comp["Nd2Fe14B"]
                * (b.cp0_Nd2Fe14B + b.cp1_Nd2Fe14B * b.temp_feed[t])
                / b.mw_Nd2Fe14B
                + b.solid_in[t].mass_frac_comp["Nd"]
                * (b.cp0_Nd + b.cp1_Nd * b.temp_feed[t])
                / b.am_Nd,
                to_units=pyunits.J / (pyunits.g * pyunits.K),
            )

        # heat loss from radiation
        @self.Constraint(self.flowsheet().config.time, doc="heat loss from radiation")
        def heat_loss_radiation(b, t):
            return b.heat_loss[t] == 0.70 * 2 * Constants.pi * b.radius_chamber[
                t
            ] * b.length_chamber[t] * Constants.stefan_constant * (
                (b.max_temperature**4) - (b.temperature_insulation_material1[t] ** 4)
            )

        # relative thickness and heat loss of insulation material 1
        @self.Constraint(
            self.flowsheet().config.time,
            doc="relative thickness of insulation material 1",
        )
        def relative_thickness_insulation_material1(b, t):
            return (
                exp(b.relative_thickness_ratio_insulation_material1[t])
                == (b.radius_chamber[t] + b.thickness_insulation_material1[t])
                / b.radius_chamber[t]
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="heat loss of insulation material 1"
        )
        def heat_loss_insulation_material1(b, t):
            return b.heat_loss[t] * b.relative_thickness_ratio_insulation_material1[
                t
            ] == (
                2
                * Constants.pi
                * b.thermal_cond_insulation_material1
                * b.length_chamber[t]
                * (
                    b.temperature_insulation_material1[t]
                    - b.temperature_metal_material1
                )
            )

        # relative thickness and heat loss of metal material 1
        @self.Constraint(
            self.flowsheet().config.time, doc="relative thickness of metal material 1"
        )
        def relative_thickness_metal_material1(b, t):
            return exp(b.relative_thickness_ratio_metal_material1[t]) == (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
            ) / (b.radius_chamber[t] + b.thickness_insulation_material1[t])

        @self.Constraint(
            self.flowsheet().config.time, doc="heat loss of metal material 1"
        )
        def heat_loss_metal_material1(b, t):
            return b.heat_loss[t] * b.relative_thickness_ratio_metal_material1[t] == (
                2
                * Constants.pi
                * b.thermal_cond_metal_material1
                * b.length_chamber[t]
                * (b.temperature_metal_material1 - b.temperature_insulation_material2)
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="relative thickness of insulation material 2",
        )
        def relative_thickness_insulation_material2(b, t):
            return exp(b.relative_thickness_ratio_insulation_material2[t]) == (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
                + b.thickness_insulation_material2[t]
            ) / (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="heat loss of insulation material 2"
        )
        def heat_loss_insulation_material2(b, t):
            return b.heat_loss[t] * b.relative_thickness_ratio_insulation_material2[
                t
            ] == (
                (
                    2
                    * Constants.pi
                    * b.thermal_cond_insulation_material2
                    * b.length_chamber[t]
                )
                * (b.temperature_insulation_material2 - b.temperature_metal_material2)
            )

        # relative thickness and heat loss of metal material 2
        @self.Constraint(
            self.flowsheet().config.time, doc="relative thickness of metal material 2"
        )
        def relative_thickness_metal_material2(b, t):
            return exp(b.relative_thickness_ratio_metal_material2[t]) == (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
                + b.thickness_insulation_material2[t]
                + b.thickness_metal_material2[t]
            ) / (
                b.radius_chamber[t]
                + b.thickness_insulation_material1[t]
                + b.thickness_metal_material1[t]
                + b.thickness_insulation_material2[t]
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="heat loss of metal material 2"
        )
        def heat_loss_metal_material2(b, t):
            return b.heat_loss[t] * b.relative_thickness_ratio_metal_material2[t] == (
                2
                * Constants.pi
                * b.thermal_cond_metal_material2
                * b.length_chamber[t]
                * (b.temperature_metal_material2 - b.temperature_furnace_ext_surface)
            )

        # heat loss from external air convection outside the furnace
        @self.Constraint(
            self.flowsheet().config.time,
            doc="heat loss from external air convection outside the furnace",
        )
        def heat_loss_external_air_convection(b, t):
            return b.heat_loss[t] == (
                b.air_heat_transfer_coeff
                * 2
                * Constants.pi
                * b.length_chamber[t]
                * (b.temperature_furnace_ext_surface - b.ref_temp)
                * (
                    b.radius_chamber[t]
                    + b.thickness_insulation_material1[t]
                    + b.thickness_metal_material1[t]
                    + b.thickness_insulation_material2[t]
                    + b.thickness_metal_material2[t]
                )
            )

    def _make_momentum_balance(self):
        @self.Constraint(self.flowsheet().config.time, doc="momentum balance equation")
        def momentum_balance_eqn(b, t):
            if self.config.has_pressure_change is True:
                return b.gas_out[t].pressure == b.gas_in[t].pressure + b.deltaP[t]
            else:
                return b.gas_out[t].pressure == b.gas_in[t].pressure

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

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Solid Inlet": self.solid_inlet,
                "Solid Outlet": self.solid_outlet,
                "Gas Inlet": self.gas_inlet,
                "Gas Outlet": self.gas_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        exprs = {}
        exprs["Volumetric Feed Flow"] = self.flow_vol_feed[0]
        exprs["Heat Duty For Decrepitation"] = self.total_heat_duty[0]
        exprs["Heat Duty Supplemental To Decrepitation"] = self.supplied_heat_duty[0]
        exprs["Sample Mass"] = self.sample_mass[0]
        exprs["Sample Volume"] = self.sample_volume[0]
        exprs["Furnace Chamber Volume"] = self.furnace_chamber_volume[0]
        exprs["Furnace Chamber Radius"] = self.radius_chamber[0]
        exprs["Furnace Chamber Length"] = self.length_chamber[0]
        exprs["Volume Insulation 1"] = self.volume_insulation1[0]
        exprs["Weight Insulation 1"] = self.total_weight_insulation1[0]
        exprs["Internal Diameter Metal 1"] = self.internal_diameter_metal1[0]
        exprs["External Diameter Metal 1"] = self.external_diameter_metal1[0]
        exprs["Volume Metal 1"] = self.volume_metal1[0]
        exprs["Weight Metal 1"] = self.weight_metal1[0]
        exprs["Volume Insulation 2"] = self.volume_insulation2[0]
        exprs["Weight Insulation 2"] = self.total_weight_insulation2[0]
        exprs["Internal Diameter Metal 2"] = self.internal_diameter_metal2[0]
        exprs["External Diameter Metal 2"] = self.external_diameter_metal2[0]
        exprs["Volume Metal 2"] = self.volume_metal2[0]
        exprs["Weight Metal 2"] = self.weight_metal2[0]
        exprs["Furnace External Surface Area"] = self.furnace_external_surface_area[0]
        exprs["Decrepitation Duration"] = self.decrepitation_duration
        exprs["Total Processing Time"] = self.processing_time[0]

        return {"exprs": exprs}
