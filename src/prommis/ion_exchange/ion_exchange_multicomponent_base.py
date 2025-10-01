#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This model is a modified version of the model from
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/ion_exchange_base.py

# WaterTAP License Agreement

# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import json
from copy import deepcopy

# Import Pyomo
import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES libraries and components
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import InitializationError, ConfigurationError

# Import WaterTAP libraries
from watertap.core import ControlVolume0DBlock, InitializationMixin
from watertap.core.solvers import get_solver
from watertap.core.util.initialization import interval_initializer

# Import modified version of IX costing model
from prommis.ion_exchange.ion_exchange_costing import cost_ion_exchange


""" This is the Ion Exchange (IX) base model that contains relevant
general equations of the IX unit. This model calls a specific
equilibiurm model that uses the Clark model with Freundlich
equilibrium parameters. 

REFERENCES

[1] LeVan, M. D., Carta, G., & Yon, C. M. (2019).  Section 16:
Adsorption and Ion Exchange.  Perry's Chemical Engineers' Handbook,
9th Edition.

[2] Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., &
Tchobanoglous, G. (2012).  Chapter 16: Ion Exchange.  MWH's Water
Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

[3] DOWEX Ion Exchange Resins Water Conditioning Manual
https://www.lenntech.com/Data-sheets/Dowex-Ion-Exchange-Resins-Water-Conditioning-Manual-L.pdf

[4] Inamuddin, & Luqman, M. (2012).  Ion Exchange Technology I: Theory
and Materials.

[5] Vassilis J. Inglezakis and Stavros G. Poulopoulos Adsorption, Ion
Exchange and Catalysis: Design of Operations and Environmental
Applications (2006).  doi.org/10.1016/B978-0-444-52783-7.X5000-9

[6] Michaud, C.F. (2013) Hydrodynamic Design, Part 8: Flow Through Ion
Exchange Beds Water Conditioning & Purification Magazine (WC&P)
https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/

[7] EPA-WBS, https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models

Assumption: We are always below the resin_max_capacity

modified by: E. Soraya Rawlings

"""

__author__ = "Kurban Sitterley"


_log = idaeslog.getLogger(__name__)


class IonExchangeType(StrEnum):
    anion = "anion"
    cation = "cation"
    demineralize = "demineralize"


class RegenerantChem(StrEnum):
    HCl = "HCl"
    NaOH = "NaOH"
    H2SO4 = "H2SO4"
    NaCl = "NaCl"
    MeOH = "MeOH"
    single_use = "single_use"


# ESR: Add resin options
class Resin(StrEnum):
    A850 = "A850"
    XAD7HP = "XAD7HP"
    D50WX8 = "D50WX8"
    S950 = "S950"
    HPR1100Na = "HPR1100Na"
    HPR1200H = "HPR1200H"
    C160 = "C160"


class IsothermType(StrEnum):
    langmuir = "langmuir"
    freundlich = "freundlich"


class IonExchangeBaseData(InitializationMixin, UnitModelBlockData):
    """
    Base for zero-order ion exchange model.
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
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
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
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
            default=EnergyBalanceType.none,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.none.
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
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
        **default** - MomentumBalanceType.pressureTotal.
        **Valid values:** {
        **MomentumBalanceType.none** - exclude momentum balances,
        **MomentumBalanceType.pressureTotal** - single pressure balance for material,
        **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
        **MomentumBalanceType.momentumTotal** - single momentum balance for material,
        **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "target_component",
        ConfigValue(
            default="",
            domain=str,
            description="Designates targeted species for removal",
        ),
    )

    CONFIG.declare(
        "number_traps",
        ConfigValue(
            default=5,
            domain=int,
            description="Designates number of trapezoids",
        ),
    )

    CONFIG.declare(
        "c_trap_min",
        ConfigValue(
            default=1 - 3,
            domain=float,
            description="Minimum relative breakthrough concentration for estimating area under curve",
        ),
    )

    # NOTE: If resin is not "single_use", add regeneration
    CONFIG.declare(
        "regenerant",
        ConfigValue(
            default="single_use",
            domain=In(RegenerantChem),
            description="Chemical used for regeneration of fixed bed",
        ),
    )

    # Add resin name from a predefined list of names. [TODO: Improve
    # the names to make them easy to identify]
    CONFIG.declare(
        "resin",
        ConfigValue(
            default=Resin.S950,
            domain=In(Resin),
            description="Resin used for the calculation of the bed expasion and pressure drop",
        ),
    )

    # Add path to resin data as a .json file
    CONFIG.declare(
        "resin_data_path",
        ConfigValue(
            default="resin_data.json",
            domain=str,
            description="Resin data json file path ",
        ),
    )

    CONFIG.declare(
        "hazardous_waste",
        ConfigValue(
            default=False,
            domain=bool,
            description="Designates if resin and residuals contain hazardous material",
        ),
    )

    CONFIG.declare(
        "isotherm",
        ConfigValue(
            default=IsothermType.langmuir,
            domain=In(IsothermType),
            description="Designates the isotherm type to use for equilibrium calculations",
        ),
    )

    CONFIG.declare(
        "add_steady_state_approximation",
        ConfigValue(
            default=False,
            domain=bool,
            description="Designates whether to add the variables and constraints necessary to estimate steady-state effluent concentration with trapezoid method.",
        ),
    )

    def build(self):
        super().build()

        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        if target_component != "":
            if self.config.property_package.charge_comp[target_component].value > 0:
                self.ion_exchange_type = IonExchangeType.cation
            elif self.config.property_package.charge_comp[target_component].value < 0:
                self.ion_exchange_type = IonExchangeType.anion
            else:
                assert target_component in ["TDS", "Alkalinity"]
                self.ion_exchange_type = IonExchangeType.demineralize
                # raise ConfigurationError("Target ion must have non-zero charge.")

        self.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

        self.process_flow = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.process_flow.add_state_blocks(has_phase_equilibrium=False)
        self.process_flow.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.process_flow.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=False
        )
        self.process_flow.add_isothermal_assumption()
        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        prop_in = self.process_flow.properties_in[0]

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        if self.config.regenerant != RegenerantChem.single_use:

            # Add regeneration stream only when the resin is not of
            # "single_use"
            self.regeneration_stream = self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of regeneration stream",
                **tmp_dict,
            )
            regen = self.regeneration_stream[0]
            self.add_outlet_port(name="regen", block=self.regeneration_stream)

        # ==========PARAMETERS==========

        # Add resin parameters
        self._add_parameters_resin(
            resin=self.config.resin, resin_data_path=self.config.resin_data_path
        )

        self.underdrain_h = pyo.Param(
            initialize=0.5,
            mutable=True,
            units=pyo.units.m,
            doc="Underdrain height",  # from ref[1] Perry's
        )

        self.distributor_h = pyo.Param(
            initialize=0.5,
            mutable=True,
            units=pyo.units.m,
            doc="Distributor height",  # from ref[1] Perry's
        )

        # Particle Peclet number correlation Eq. 4.100 in
        # ref[4]. NOTE: From ref[5], this is the "L" term and it has a
        # value of 0.050 for downflow and 0.52 for upflow.
        self.Pe_p_A = pyo.Param(
            initialize=0.05,
            units=pyo.units.dimensionless,
            doc="Peclet particle equation A parameter",
        )

        # From ref[5], this is the "k" term and it has a value of
        # 0.48 for downflow and -0.65 for upflow.
        self.Pe_p_exp = pyo.Param(
            initialize=0.48,
            units=pyo.units.dimensionless,
            doc="Peclet particle equation exponent",
        )

        # Sherwood number as a function of Reynolds and Schmidt number
        # Table 16-9 in ref[1] Perry's
        # Wilson and Geankoplis, Ind. Eng. Chem. Fundam., 5, 9 (1966)
        self.Sh_A = pyo.Param(
            initialize=2.4,
            units=pyo.units.dimensionless,
            doc="Sherwood equation A parameter",
        )
        self.Sh_exp_A = pyo.Param(
            initialize=0.66,
            units=pyo.units.dimensionless,
            doc="Sherwood equation exponent A",
        )
        self.Sh_exp_B = pyo.Param(
            initialize=0.34,
            units=pyo.units.dimensionless,
            doc="Sherwood equation exponent B",
        )
        self.Sh_exp_C = pyo.Param(
            initialize=0.33,
            units=pyo.units.dimensionless,
            doc="Sherwood equation exponent C",
        )

        # ==========IX parameters==========

        self.pump_efficiency = pyo.Param(
            initialize=0.8,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Pump efficiency",
        )

        # Add regeneration-related variables only when regeneration is
        # True (or not "single_use")
        if self.config.regenerant != RegenerantChem.single_use:

            self.regeneration_time = pyo.Param(
                initialize=1800,
                mutable=True,
                units=pyo.units.s,
                doc="Regeneration time",
            )

            self.service_to_regen_flow_ratio = pyo.Param(
                initialize=3,
                mutable=True,
                units=pyo.units.dimensionless,
                doc="Ratio of service flow rate to regeneration flow rate",
            )

            # [TODO: Refine the following parameters once we include
            # the regeneration stage.]

            # Rinse, Regen, Backwashing params
            self.rinse_bed_volumes = pyo.Param(
                initialize=5,
                mutable=True,
                doc="Number of bed volumes for rinse step",
            )

            @self.Expression(doc="Rinse time")
            def rinse_time(b):
                return b.ebct * b.rinse_bed_volumes

            @self.Expression(doc="Waste time")
            def waste_time(b):
                return b.regeneration_time + b.backwash_time + b.rinse_time

            @self.Expression(doc="Cycle time")
            def cycle_time(b):
                return b.target_breakthrough_time + b.waste_time

        self.backwashing_rate = pyo.Param(
            initialize=5,
            mutable=True,
            units=pyo.units.m / pyo.units.hour,
            doc="Backwash loading rate [m/hr]",
        )
        self.backwash_time = pyo.Param(
            initialize=600,
            mutable=True,
            units=pyo.units.s,
            doc="Backwash time",
        )
        self.redundant_column_freq = pyo.Param(
            initialize=4,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Frequency for redundant columns",
        )

        # ==========VARIABLES==========

        self.resin_diam = pyo.Var(
            initialize=7e-4,
            bounds=(0, 5e-3),  # add bounds valid for resins for REEs recovery
            # domain=NonNegativeReals,
            units=pyo.units.m,
            doc="Resin bead diameter",
        )
        self.resin_density = pyo.Var(
            initialize=700,
            bounds=(100, 2000),  # add bounds valid for resins for REEs recovery
            # domain=NonNegativeReals,
            units=pyo.units.kg / pyo.units.m**3,
            doc="Resin bulk density",
        )
        self.bed_volume = pyo.Var(
            initialize=2,
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.m**3,
            doc="Bed volume per column",
        )
        self.bed_volume_total = pyo.Var(
            initialize=2,
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.m**3,
            doc="Total bed volume",
        )
        self.bed_depth = pyo.Var(
            initialize=1,
            bounds=(0.75, 2),  # from ref[7] EPA-WBS guidance
            # domain=NonNegativeReals,
            units=pyo.units.m,
            doc="Bed depth",
        )
        self.bed_porosity = pyo.Var(
            initialize=0.4,
            bounds=(0.3, 0.8),
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Bed porosity",
        )

        self.column_height = pyo.Var(
            initialize=2,
            bounds=(0, 4.26),  # from ref[7] EPA-WBS guidance
            # domain=NonNegativeReals,
            units=pyo.units.m,
            doc="Column height",
        )
        self.bed_diameter = pyo.Var(
            initialize=1,
            bounds=(0.75, 4.26),  # from ref[7] EPA-WBS guidance
            # domain=NonNegativeReals,
            units=pyo.units.m,
            doc="Column diameter",
        )

        self.number_columns = pyo.Var(
            initialize=2,
            bounds=(1, None),
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Number of operational columns for ion exchange process",
        )
        self.number_columns_redundant = pyo.Var(
            initialize=1,
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Number of redundant columns for ion exchange process",
        )

        # For the multicomponent model, add a new variable to save the
        # breakthrough time for the target component and use it
        # instead of the original breakthrough_time. [TODO: Figure out
        # a way to calculate the maximum breakthrough. For now we use
        # the target component.]
        self.target_breakthrough_time = pyo.Var(
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.s,
            doc="(Max) Breakthrough time of target component",
        )
        self.ebct = pyo.Var(
            initialize=520,
            bounds=(90, None),
            # domain=NonNegativeReals,
            units=pyo.units.s,
            doc="Empty bed contact time",
        )

        # ====== Hydrodynamic variables ====== #

        self.loading_rate = pyo.Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # ref[1], [2], and [7] (MWH, Perry's, and EPA-WBS)
            # domain=NonNegativeReals,
            units=pyo.units.m / pyo.units.s,
            doc="Superficial velocity through bed",
        )
        self.service_flow_rate = pyo.Var(
            initialize=10,
            bounds=(1, 40),
            # domain=NonNegativeReals,
            units=pyo.units.hr**-1,
            doc="Service flow rate [BV/hr]",
        )

        # ====== Dimensionless variables ====== #

        self.N_Re = pyo.Var(
            initialize=4.3,
            bounds=(0.001, 60),  # ref[1] Perry's - bounds relevant to N_Sh regression
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Reynolds number",
        )
        self.N_Pe_particle = pyo.Var(
            initialize=0.1,
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Peclet particle number",
        )
        self.N_Pe_bed = pyo.Var(
            initialize=1000,  # ref[4] Inamuddin/Luqman - N_Pe_bed > ~100 considered to be plug flow
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Peclet bed number",
        )

        # ==========EXPRESSIONS==========

        @self.Expression(doc="Flow per column")
        def flow_per_column(b):
            return prop_in.flow_vol_phase["Liq"] / b.number_columns

        @self.Expression(doc="Cross-sectional area of one column")
        def bed_area(b):
            return pyo.units.convert(
                Constants.pi * (b.bed_diameter / 2) ** 2, to_units=pyo.units.m**2
            )

        if self.config.regenerant != RegenerantChem.single_use:

            # Add regeneration variables/expressions when the resin is
            # not "single_use" and regeneration is needed. [TODO: Check
            # if all these expressions/variables should be included in
            # this "if" when including the regeneration stage.]
            @self.Expression(doc="Rinse time")
            def rinse_time(b):
                return b.ebct * b.rinse_bed_volumes

            @self.Expression(doc="Waste time")
            def waste_time(b):
                return b.regeneration_time + b.backwash_time + b.rinse_time

        @self.Expression(doc="Cycle time")
        def cycle_time(b):
            if self.config.regenerant == RegenerantChem.single_use:
                return b.target_breakthrough_time
            else:
                return b.target_breakthrough_time + b.waste_time

        # Add resin-specific properties
        self.add_properties_resin()

        if self.config.regenerant != RegenerantChem.single_use:

            # If resin is not single use, add regeneration
            # units. [TODO: Confirm these expressions are correct once
            # we add the regeneration stage.]

            @self.Expression(doc="Regen pump power")
            def regen_pump_power(b):
                return pyo.units.convert(
                    (
                        b.pressure_drop
                        * (
                            prop_in.flow_vol_phase["Liq"]
                            / b.service_to_regen_flow_ratio
                        )
                    )
                    / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.regeneration_time / b.cycle_time)

            @self.Expression(doc="Regen tank volume")
            def regen_tank_vol(b):
                return (
                    prop_in.flow_vol_phase["Liq"] / b.service_to_regen_flow_ratio
                ) * b.regeneration_time

            @self.Constraint(
                doc="Isothermal assumption for regen stream",
            )
            def eq_isothermal_regen_stream(b):
                return prop_in.temperature == regen.temperature

            @self.Constraint(
                doc="Isobaric assumption for regen stream",
            )
            def eq_isobaric_regen_stream(b):
                return prop_in.pressure == regen.pressure

        @self.Expression(doc="Backwashing flow rate")
        def bw_flow(b):
            return (
                pyo.units.convert(
                    b.backwashing_rate, to_units=pyo.units.m / pyo.units.s
                )
                * b.bed_area
                * b.number_columns
            )

        @self.Expression(doc="Rinse flow rate")
        def rinse_flow(b):
            return b.loading_rate * b.bed_area * b.number_columns

        # Add regen terms when regeneration is needed (or not
        # "single_use")
        if self.config.regenerant != "single_use":

            @self.Expression(doc="Backwash pump power")
            def bw_pump_power(b):
                return pyo.units.convert(
                    (b.pressure_drop * b.bw_flow) / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.backwash_time / b.cycle_time)

            @self.Expression(doc="Rinse pump power")
            def rinse_pump_power(b):
                time = b.rinse_time / b.cycle_time
                return pyo.units.convert(
                    (b.pressure_drop * b.rinse_flow) / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.rinse_time / b.cycle_time)

            @self.Expression(doc="Main pump power")
            def main_pump_power(b):
                return pyo.units.convert(
                    (b.pressure_drop * prop_in.flow_vol_phase["Liq"])
                    / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.target_breakthrough_time / b.cycle_time)

        @self.Expression(doc="Bed expansion from backwashing")
        def bed_expansion_h(b):
            return b.bed_expansion_frac * b.bed_depth

        @self.Expression(doc="Free board needed")
        def free_board(b):
            return b.distributor_h + b.underdrain_h + b.bed_expansion_h

        @self.Expression(doc="Volume per column")
        def column_volume(b):
            return b.column_height * b.bed_area

        @self.Expression(doc="Total column volume required")
        def column_volume_total(b):
            return b.number_columns * b.column_volume

        @self.Expression(doc="Contact time")
        def t_contact(b):
            return b.ebct * b.bed_porosity

        @self.Expression(doc="Interstitial velocity")
        def vel_inter(b):
            return b.loading_rate / b.bed_porosity

        @self.Expression(doc="Total number of columns")
        def number_columns_total(b):
            return b.number_columns + b.number_columns_redundant

        @self.Constraint(doc="Reynolds number")
        def eq_Re(b):  # Eq. 3.358, ref[5] Inglezakis + Poulopoulos
            return b.N_Re * prop_in.visc_k_phase["Liq"] == (
                b.loading_rate * b.resin_diam
            )

        @self.Constraint(doc="Bed Peclet number")
        def eq_Pe_bed(b):
            return b.N_Pe_bed * b.resin_diam == b.N_Pe_particle * b.bed_depth

        @self.Constraint(doc="Particle Peclet number")
        def eq_Pe_p(b):  # Eq. 3.313, Inglezakis + Poulopoulos, for downflow
            return b.N_Pe_particle == b.Pe_p_A * b.N_Re**b.Pe_p_exp

        # =========== RESIN & COLUMN ===========

        @self.Constraint(doc="Empty bed contact time")
        def eq_ebct(b):
            # NOTE: Replace original equation, since units were not
            # consistent
            return b.ebct == b.bed_volume / b.flow_per_column

        @self.Constraint(doc="Loading rate")
        def eq_loading_rate(b):
            return b.loading_rate * b.bed_area == prop_in.flow_vol_phase["Liq"]

        @self.Constraint(doc="Service flow rate")
        def eq_service_flow_rate(b):
            return b.service_flow_rate * b.bed_volume_total == pyo.units.convert(
                prop_in.flow_vol_phase["Liq"],
                to_units=pyo.units.m**3 / pyo.units.hr,
            )

        @self.Constraint(doc="Bed volume per operational column")
        def eq_bed_volume(b):
            return b.bed_volume == b.bed_area * b.bed_depth

        @self.Constraint(doc="Total bed volume")
        def eq_bed_design(b):
            return b.bed_volume_total == b.bed_volume * b.number_columns

        @self.Constraint(doc="Column height")
        def eq_column_height(b):
            return b.column_height == b.bed_depth + b.free_board

    def add_ss_approximation(self, ix_model_type=None):

        prop_in = self.process_flow.properties_in[0]
        self.num_traps = self.config.number_traps
        self.trap_disc = range(self.num_traps + 1)
        self.trap_index = self.trap_disc[1:]

        # Read c_trap_min from its configuration value. The c_trap_min
        # is the minimum relative breakthrough concentration for
        # estimating area under curve. [TODO: Check if this is enough
        # for the case of multiple components.]
        self.c_trap_min = float(self.config.c_trap_min)

        self.c_traps = pyo.Var(
            self.trap_disc,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
        )
        self.tb_traps = pyo.Var(
            self.trap_disc,
            initialize=1e6,
            bounds=(0, None),
            units=pyo.units.second,
            doc="Breakthrough times for estimating area under breakthrough curve",
        )
        self.traps = pyo.Var(
            self.trap_index,
            initialize=0.01,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Trapezoid areas for estimating area under breakthrough curve",
        )
        self.c_norm_avg = pyo.Var(
            self.target_component_set,
            initialize=0.25,
            bounds=(0, 2),
            units=pyo.units.dimensionless,
            doc="Sum of trapezoid areas",
        )

        # Fix conditions at initial point in the trapezoids
        self.c_traps[0].fix(0)
        self.tb_traps[0].fix(0)

        @self.Constraint(
            self.target_component_set,
            self.trap_index,
            doc="Evenly spaced c_norm for trapezoids",
        )
        def eq_c_traps(b, j, k):
            return b.c_traps[k] * (b.num_traps - 1) == (
                b.c_trap_min * (b.num_traps - 1)
                + (b.trap_disc[k] - 1) * (b.c_norm[j] - b.c_trap_min)
            )

        @self.Constraint(
            self.trap_index,
            doc="Breakthru time calc for trapezoids",
        )
        def eq_tb_traps(b, k):
            if ix_model_type == "clark":
                bv_traps = (b.tb_traps[k] * b.loading_rate) / b.bed_depth
                left_side = pyo.exp(
                    (
                        (b.mass_transfer_coeff * b.bed_depth * (b.freundlich_n - 1))
                        / (b.bv_50 * b.loading_rate)
                    )
                    * (b.bv_50 - bv_traps)
                )
                right_side = ((1 / b.c_traps[k]) ** (b.freundlich_n - 1) - 1) / (
                    2 ** (b.freundlich_n - 1) - 1
                )
                return left_side - right_side == 0
            else:
                return Constraint.Skip

        @self.Constraint(self.trap_index, doc="Area of trapezoids")
        def eq_traps(b, k):
            return b.traps[k] * b.tb_traps[self.num_traps] == (
                (b.tb_traps[k] - b.tb_traps[k - 1])
                * ((b.c_traps[k] + b.c_traps[k - 1]) / 2)
            )

        @self.Constraint(
            self.target_component_set, doc="Average relative effluent concentration"
        )
        def eq_c_norm_avg(b, j):
            return b.c_norm_avg[j] == sum(b.traps[k] for k in b.trap_index)

        @self.Constraint(
            self.target_component_set,
            doc="CV mass transfer term",
        )
        def eq_mass_transfer_term(b, j):
            return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
                "Liq", j
            ) == -b.process_flow.mass_transfer_term[0, "Liq", j]

    # Add resin parameters and expressions in separate functions that
    # read data from a .json file
    def _add_parameters_resin(self, resin=None, resin_data_path=None):
        """Method to add resin-specific parameters needed to calculate the
        pressure drop and bed expansion properties for resins.

        This method adds parameters for the selected resin.

        """

        # Read resin parameters from .json file
        with open(resin_data_path) as data_file:
            resin_params = json.load(data_file)

        self.resin_params = resin_params

        self.bed_expansion_param_a = resin_params[resin]["bed_expansion"]["param_a"]
        self.bed_expansion_param_b = resin_params[resin]["bed_expansion"]["param_b"]
        self.bed_expansion_param_c = resin_params[resin]["bed_expansion"]["param_c"]

        self.pressure_drop_param_a = resin_params[resin]["pressure_drop"]["param_a"]
        self.pressure_drop_param_b = resin_params[resin]["pressure_drop"]["param_b"]
        self.pressure_drop_param_c = resin_params[resin]["pressure_drop"]["param_c"]

        # Parameters to calculate the pressure drop
        self.pressure_drop_param_A = pyo.Param(
            initialize=self.pressure_drop_param_a,
            mutable=True,
            units=pyo.units.psi / pyo.units.m,
            doc="Parameter a for pressure drop calculation",
        )
        self.pressure_drop_param_B = pyo.Param(
            initialize=self.pressure_drop_param_b,
            mutable=True,
            units=(pyo.units.psi * pyo.units.hr) / pyo.units.m**2,
            doc="Parameter b for pressure drop calculation",
        )
        self.pressure_drop_param_C = pyo.Param(
            initialize=self.pressure_drop_param_c,
            mutable=True,
            units=(pyo.units.psi * pyo.units.hr**2) / pyo.units.m**3,
            doc="Parameter c for pressure drop calculation",
        )

        # Parameters to calculate the bed expansion fraction
        self.bed_expansion_frac_param_A = pyo.Param(
            initialize=self.bed_expansion_param_a,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Parameter a for bed expansion fraction calculation",
        )
        self.bed_expansion_frac_param_B = pyo.Param(
            initialize=self.bed_expansion_param_b,
            mutable=True,
            units=pyo.units.hr / pyo.units.m,
            doc="Parameter b for bed expansion fraction calculation",
        )
        self.bed_expansion_frac_param_C = pyo.Param(
            initialize=self.bed_expansion_param_c,
            mutable=True,
            units=pyo.units.hr**2 / pyo.units.m**2,
            doc="Parameter c for bed expansion fraction calculation",
        )

    def add_properties_resin(self):
        """Method to add resin-specific properties for the selected resin.

        NOTE: When using this method, keep in mind the temperature for
        the data for each property given in the .json file

        """

        # Pressure drop (psi/m of resin bed depth) is a function of
        # loading rate (loading_rate) in m/hr following the equation:
        # p_drop (psi/m) = p_drop_A + p_drop_B * loading_rate + p_drop_C * loading_rate**2
        @self.Expression(doc="Pressure drop")
        def pressure_drop(b):
            loading_rate_conv = pyo.units.convert(
                b.loading_rate, to_units=pyo.units.m / pyo.units.hr
            )
            return (
                b.pressure_drop_param_A
                + b.pressure_drop_param_B * loading_rate_conv
                + b.pressure_drop_param_C * loading_rate_conv**2
            ) * b.bed_depth

        # Bed expansion is calculated as a fraction of the bed_depth
        # These coefficients are used to calculate that fraction
        # (bed_expansion_frac) as a function of backwash rate
        # (bw_rate, m/hr):
        # bed_expansion_frac = bed_expansion_A bed_expansion_B *
        #                      bw_rate bed_expansion_C * bw_rate**2
        @self.Expression(doc="Bed expansion fraction from backwashing")
        def bed_expansion_frac(b):
            return (
                b.bed_expansion_frac_param_A
                + b.bed_expansion_frac_param_B * b.backwashing_rate
                + b.bed_expansion_frac_param_C * b.backwashing_rate**2
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.process_flow.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")

        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
            state_dict = self.process_flow.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args_out = deepcopy(state_args)

        state_args_regen = deepcopy(state_args)

        # Add regen terms when regeneration is needed
        if self.config.regenerant != "single_use":
            self.regeneration_stream.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_regen,
            )

        init_log.info("Initialization Step 1c Complete.")
        # interval_initializer(self)

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not pyo.check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.target_breakthrough_time) is None:
            iscale.set_scaling_factor(self.target_breakthrough_time, 1e-6)

        if iscale.get_scaling_factor(self.N_Re) is None:
            iscale.set_scaling_factor(self.N_Re, 1)

        if iscale.get_scaling_factor(self.N_Pe_particle) is None:
            iscale.set_scaling_factor(self.N_Pe_particle, 1e2)

        if iscale.get_scaling_factor(self.N_Pe_bed) is None:
            iscale.set_scaling_factor(self.N_Pe_bed, 1e-3)

        if iscale.get_scaling_factor(self.number_columns) is None:
            iscale.set_scaling_factor(self.number_columns, 1)

        if iscale.get_scaling_factor(self.resin_diam) is None:
            iscale.set_scaling_factor(self.resin_diam, 1e4)

        if iscale.get_scaling_factor(self.resin_density) is None:
            iscale.set_scaling_factor(self.resin_density, 1e-3)

        if iscale.get_scaling_factor(self.bed_volume_total) is None:
            iscale.set_scaling_factor(self.bed_volume_total, 0.1)

        if iscale.get_scaling_factor(self.bed_depth) is None:
            iscale.set_scaling_factor(self.bed_depth, 1)

        if iscale.get_scaling_factor(self.bed_porosity) is None:
            iscale.set_scaling_factor(self.bed_porosity, 10)

        if iscale.get_scaling_factor(self.column_height) is None:
            iscale.set_scaling_factor(self.column_height, 1)

        if iscale.get_scaling_factor(self.bed_diameter) is None:
            iscale.set_scaling_factor(self.bed_diameter, 1)

        if iscale.get_scaling_factor(self.service_flow_rate) is None:
            iscale.set_scaling_factor(self.service_flow_rate, 0.1)

        if iscale.get_scaling_factor(self.ebct) is None:
            iscale.set_scaling_factor(self.ebct, 1e-2)

        if iscale.get_scaling_factor(self.loading_rate) is None:
            iscale.set_scaling_factor(self.loading_rate, 1e3)

        if iscale.get_scaling_factor(self.bv) is None:
            iscale.set_scaling_factor(self.bv, 1e-4)

    def _get_stream_table_contents(self, time_point=0):

        if self.config.regenerant != "single_use":
            return create_stream_table_dataframe(
                {
                    "Feed Inlet": self.inlet,
                    "Liquid Outlet": self.outlet,
                    "Regen Outlet": self.regen,
                },
                time_point=time_point,
            )
        else:
            return create_stream_table_dataframe(
                {
                    "Feed Inlet": self.inlet,
                    "Liquid Outlet": self.outlet,
                },
                time_point=time_point,
            )

    def _get_performance_contents(self, time_point=0):

        # Add relevant Params, Expressions, differences for CONFIGs
        target_component = self.config.target_component
        var_dict = {}
        var_dict["Max Breakthrough Time"] = self.target_breakthrough_time
        var_dict["EBCT"] = self.ebct
        var_dict["Number Columns"] = self.number_columns
        var_dict["Bed Volume Total"] = self.bed_volume_total
        var_dict["Bed Depth"] = self.bed_depth
        # var_dict["Col. Height to Diam. Ratio"] = self.col_height_to_diam_ratio
        var_dict["Bed Porosity"] = self.bed_porosity
        var_dict["Service Flow Rate [BV/hr]"] = self.service_flow_rate
        var_dict["Bed Velocity"] = self.loading_rate
        var_dict["Resin Particle Diameter"] = self.resin_diam
        var_dict["Resin Bulk Density"] = self.resin_density
        var_dict["Reynolds Number"] = self.N_Re
        var_dict["Peclet Number (bed)"] = self.N_Pe_bed
        var_dict["Peclet Number (particle)"] = self.N_Pe_particle

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_ion_exchange
