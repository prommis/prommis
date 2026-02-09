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

r"""

Ion Exchange Multicomponent (IXMC) Base Model
=============================================

The Ion Exchange Multicomponent (IXMC) base model addresses the
following specific aspects of the ion exchange process:

- Provides the foundational structure for the IXMC model, ensuring
  consistency across all calculations and enabling users to define
  global key variables and parameters.

- Contains all global variables and parameters required for modeling
  the ion exchange column, including terms related to column geometry
  and resin properties. Most of these variables, parameters, and
  expressions are included under the **Model Components** Table in the
  `WaterTAP IX documentation
  <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/ion_exchange_0D.html>`_.

"""

import json
from copy import deepcopy

# Import Pyomo
import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES libraries and components
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core import (
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
from idaes.core.util.exceptions import ConfigurationError

# Import WaterTAP libraries
from watertap.core import ControlVolume0DBlock, InitializationMixin
from watertap.core.solvers import get_solver


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
    cation = "cation"
    # [ESR WIP: The current ion exchange model is used for the
    # separation of rare earth elements (REEs), which implies a cation
    # exchange type. We commented other types for now. We will
    # consider if we need to integrate the alternative types in future
    # iterations.]
    # anion = "anion"
    # demineralize = "demineralize"


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
            description="Resin used for the calculation of the bed expansion and pressure drop",
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

    # [ESR WIP: The ion_exchange_multicomponent model (Clark)
    # incorporates the trapezoidal rule within the model. For now, we
    # will comment this out and evaluate whether to keep it in the
    # base model in future iterations.]
    # CONFIG.declare(
    #     "add_steady_state_approximation",
    #     ConfigValue(
    #         default=False,
    #         domain=bool,
    #         description="Designates whether to add the variables and constraints necessary to estimate steady-state effluent concentration with trapezoid method.",
    #     ),
    # )

    def build(self):
        super().build()

        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        if target_component != "":
            if self.config.property_package.charge_comp[target_component].value > 0:
                self.ion_exchange_type = IonExchangeType.cation
            else:
                # [ESR WIP: The current ion exchange model is used for
                # the separation of rare earth elements (REEs), which
                # implies a cation exchange type. We commented other
                # types for now. We will consider if we need to
                # integrate the alternative types in future
                # iterations.]
                raise ConfigurationError(
                    "The current ion exchange model is limited to cation exchange methods and alternative techniques are not addressed at this time."
                )

            # elif self.config.property_package.charge_comp[target_component].value < 0:
            #     self.ion_exchange_type = IonExchangeType.anion
            # else:
            #     assert target_component in ["TDS", "Alkalinity"]
            #     self.ion_exchange_type = IonExchangeType.demineralize
            #     # raise ConfigurationError("Target ion must have non-zero charge.")

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

        # ==========PARAMETERS==========

        self.pump_efficiency = pyo.Param(
            initialize=0.8,
            mutable=True,
            units=pyo.units.dimensionless,
            doc="Pump efficiency",
        )

        # Add regeneration-related parameters only when regeneration
        # is True (or not "single_use"). [ESR TODO: Refine these
        # parameters once we include the regeneration stage.]
        if self.config.regenerant != RegenerantChem.single_use:

            self.regen_dose = pyo.Param(
                initialize=300,
                units=pyo.units.kg / pyo.units.m**3,  # kg regenerant / m3 resin
                mutable=True,
                doc="Regenerant dose required for regeneration per volume of resin",
            )

            self.regen_soln_conc = pyo.Param(
                initialize=107,
                units=pyo.units.kg / pyo.units.m**3,  # kg regenerant / m3 solution
                mutable=True,
                doc="Concentration of regeneration solution",
            )

            # For NaCl, saturation concentration is 372 kg/m3
            self.regen_soln_conc_sat = pyo.Param(
                initialize=372,  # default is for NaCl
                units=pyo.units.kg / pyo.units.m**3,  # kg regenerant / m3 solution
                mutable=True,
                doc="Concentration of regeneration solution at saturation",
            )

            # For NaCl solutions, density for 8-14% wt is ~1050-1100
            # see: https://www.handymath.com/cgi-bin/nacltble.cgi?submit=Entry
            self.regen_soln_dens = pyo.Param(
                initialize=1080,
                units=pyo.units.kg / pyo.units.m**3,  # kg solution / m3 solution
                mutable=True,
                doc="Density of regeneration solution",
            )

            self.regen_recycle = pyo.Param(
                initialize=1,
                units=pyo.units.dimensionless,
                mutable=True,
                doc="Number of cycles the regenerant can be reused before disposal",
            )

            self.regen_rate = pyo.Param(
                initialize=2.4,  # EPA-WBS: 0.25-0.4 gpm/ft3 resin; 2.0-3.2 (m3/hr)/m3 resin
                mutable=True,
                units=pyo.units.hr**-1,
                doc="Regeneration flow rate per resin volume",
            )

            self.num_regen_columns = pyo.Param(
                initialize=1,  # next integer above waste_time / breakthrough_time
                units=pyo.units.dimensionless,
                mutable=True,
                doc="Number of columns regenerated at a time",
            )

            self.rinse_bed_volumes = pyo.Param(
                initialize=5,
                mutable=True,
                doc="Number of bed volumes for rinse step",
            )

        self.backwash_loading_rate = pyo.Param(
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
            domain=pyo.NonNegativeReals,
            units=pyo.units.m,
            doc="Resin bead diameter",
        )
        self.resin_density = pyo.Var(
            initialize=700,
            bounds=(100, 2000),  # add bounds valid for resins for REEs recovery
            domain=pyo.NonNegativeReals,
            units=pyo.units.kg / pyo.units.m**3,
            doc="Resin bulk density",
        )
        self.bed_volume = pyo.Var(
            initialize=2,
            bounds=(0, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.m**3,
            doc="Bed volume per column",
        )
        self.bed_volume_total = pyo.Var(
            initialize=2,
            bounds=(0, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.m**3,
            doc="Total bed volume",
        )
        self.bed_depth = pyo.Var(
            initialize=1,
            bounds=(0.75, 2),  # from ref[7] EPA-WBS guidance
            domain=pyo.NonNegativeReals,
            units=pyo.units.m,
            doc="Bed depth",
        )
        self.bed_porosity = pyo.Var(
            initialize=0.4,
            bounds=(0.3, 0.8),
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Bed porosity",
        )

        self.column_height = pyo.Var(
            initialize=2,
            bounds=(0, 4.26),  # from ref[7] EPA-WBS guidance
            domain=pyo.NonNegativeReals,
            units=pyo.units.m,
            doc="Column height",
        )
        self.bed_diameter = pyo.Var(
            initialize=1,
            bounds=(0.75, 4.26),  # from ref[7] EPA-WBS guidance
            domain=pyo.NonNegativeReals,
            units=pyo.units.m,
            doc="Column diameter",
        )

        self.number_columns = pyo.Var(
            initialize=2,
            bounds=(1, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Number of operational columns for ion exchange process",
        )
        self.number_columns_redundant = pyo.Var(
            initialize=1,
            bounds=(0, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Number of redundant columns for ion exchange process",
        )

        # When multicomponent, define a new variable to save the
        # breakthrough time for the target component and use it
        # instead of the original breakthrough_time. [ESR TODO: Figure
        # out a way to calculate the maximum breakthrough. For now we
        # use the target component.]
        self.target_breakthrough_time = pyo.Var(
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.s,
            doc="(Max) Breakthrough time of target component",
        )
        self.ebct = pyo.Var(
            initialize=520,
            bounds=(90, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.s,
            doc="Empty bed contact time",
        )

        # ====== Hydrodynamic variables ====== #

        self.loading_rate = pyo.Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # ref[1], [2], and [7] (MWH, Perry's, and EPA-WBS)
            domain=pyo.NonNegativeReals,
            units=pyo.units.m / pyo.units.s,
            doc="Superficial velocity through bed",
        )
        self.service_flow_rate = pyo.Var(
            initialize=10,
            bounds=(1, 40),
            domain=pyo.NonNegativeReals,
            units=pyo.units.hr**-1,
            doc="Service flow rate [BV/hr]",
        )

        # ====== Dimensionless variables ====== #

        self.N_Re = pyo.Var(
            initialize=4.3,
            bounds=(0.001, 60),  # ref[1] Perry's - bounds relevant to N_Sh regression
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Reynolds number",
        )
        self.N_Pe_particle = pyo.Var(
            initialize=0.1,
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Peclet particle number",
        )
        self.N_Pe_bed = pyo.Var(
            initialize=1000,  # ref[4] Inamuddin/Luqman - N_Pe_bed > ~100 considered to be plug flow
            domain=pyo.NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Peclet bed number",
        )

        # ==========EXPRESSIONS==========

        # Add resin-specific properties
        self.add_properties_resin()

        @self.Expression(doc="Flow per column")
        def flow_per_column(b):
            return prop_in.flow_vol_phase["Liq"] / b.number_columns

        @self.Expression(doc="Cross-sectional area of one column")
        def bed_area(b):
            return pyo.units.convert(
                Constants.pi * (b.bed_diameter / 2) ** 2, to_units=pyo.units.m**2
            )

        # [ESR WIP: Define waste time expression that depends only on
        # variables when "single_use" is selected.]
        waste_time_expr = self.backwash_time

        if self.config.regenerant != RegenerantChem.single_use:

            # Add regeneration expressions when regeneration is needed
            # (or resin is not "single_use"). [ESR TODO: Confirm these
            # expressions should be included in this "if" when
            # including the regeneration stage.]
            @self.Expression(doc="Rinse time")
            def rinse_time(b):
                return b.ebct * b.rinse_bed_volumes

            @self.Expression(doc="Regen time")
            def regeneration_time(b):
                return pyo.units.convert(
                    b.regen_dose / b.regen_soln_conc / b.regen_rate,
                    to_units=pyo.units.s,
                )

            waste_time_expr += self.regeneration_time + self.rinse_time

            @self.Expression(doc="Rinse flow rate")
            def rinse_flow_rate(b):
                return b.loading_rate * b.bed_area * b.number_columns

            @self.Expression(doc="Regen flow rate")
            def regen_flow_rate(b):
                return pyo.units.convert(
                    b.bed_volume * b.number_columns * b.regen_rate,
                    to_units=pyo.units.m**3 / pyo.units.s,
                )

        @self.Expression(doc="Waste time")
        def waste_time(b):
            return waste_time_expr

        @self.Expression(doc="Cycle time")
        def cycle_time(b):
            return b.target_breakthrough_time + b.waste_time

        @self.Expression(doc="Fraction of cycle time for regen + backwash + rinse")
        def frac_waste_time(b):
            return b.waste_time / b.cycle_time

        @self.Expression(doc="Backwashing flow rate")
        def backwash_flow_rate(b):
            return (
                pyo.units.convert(
                    b.backwash_loading_rate, to_units=pyo.units.m / pyo.units.s
                )
                * b.bed_area
                * b.number_columns
            )

        @self.Expression(doc="Backwash pump power")
        def backwash_pump_power(b):
            return pyo.units.convert(
                (b.pressure_drop * b.backwash_flow_rate) / b.pump_efficiency,
                to_units=pyo.units.kilowatts,
            ) * (b.backwash_time / b.cycle_time)

        @self.Expression(doc="Main pump power")
        def main_pump_power(b):
            return pyo.units.convert(
                (b.pressure_drop * prop_in.flow_vol_phase["Liq"]) / b.pump_efficiency,
                to_units=pyo.units.kilowatts,
            ) * (b.target_breakthrough_time / b.cycle_time)

        # Add regen terms when regeneration is needed (or not
        # "single_use")
        if self.config.regenerant != "single_use":

            @self.Expression(doc="Rinse pump power")
            def rinse_pump_power(b):
                return pyo.units.convert(
                    (b.pressure_drop * b.rinse_flow_rate) / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.rinse_time / b.cycle_time)

            @self.Expression(doc="Regen pump power")
            def regen_pump_power(b):
                return pyo.units.convert(
                    (b.pressure_drop * b.regen_flow_rate * b.regen_recycle)
                    / b.pump_efficiency,
                    to_units=pyo.units.kilowatts,
                ) * (b.regeneration_time / b.cycle_time)

            @self.Expression(doc="Regen tank volume")
            def regen_tank_vol(b):
                # Large enough to hold saturated regen solution for
                # one regeneration cycle
                return pyo.units.convert(
                    b.num_regen_columns
                    * ((b.bed_volume * b.regen_dose) / b.regen_soln_conc_sat),
                    to_units=pyo.units.m**3,
                )

            @self.Expression(doc="Total residuals volume: backwash + regen + rinse")
            def total_residuals_vol(b):
                return pyo.units.convert(
                    b.backwash_flow_rate * b.backwash_time
                    + b.regen_flow_rate * b.regeneration_time
                    + b.rinse_flow_rate * b.rinse_time,
                    to_units=pyo.units.m**3,
                )

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

        # =========== CONSTRAINTS ===========

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

        # =========== CONSTRAINTS RESIN & COLUMN ===========

        @self.Constraint(doc="Empty bed contact time")
        def eq_ebct(b):
            # [ESR NOTE: Replace original equation since units were
            # not consistent.]
            return b.ebct * b.flow_per_column == b.bed_volume

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

        # [ESR NOTE: Make sure the pressure ratio has positive values
        # even at low loading rate values.]
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
                + b.bed_expansion_frac_param_B * b.backwash_loading_rate
                + b.bed_expansion_frac_param_C * b.backwash_loading_rate**2
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

        iscale.set_scaling_factor(self.target_breakthrough_time, 1e-6)
        iscale.set_scaling_factor(self.N_Re, 1)
        iscale.set_scaling_factor(self.N_Pe_particle, 1e2)
        iscale.set_scaling_factor(self.N_Pe_bed, 1e-3)
        iscale.set_scaling_factor(self.number_columns, 1)
        iscale.set_scaling_factor(self.resin_diam, 1e4)
        iscale.set_scaling_factor(self.resin_density, 1e-3)
        iscale.set_scaling_factor(self.bed_volume_total, 0.1)
        iscale.set_scaling_factor(self.bed_depth, 1)
        iscale.set_scaling_factor(self.bed_porosity, 10)
        iscale.set_scaling_factor(self.column_height, 1)
        iscale.set_scaling_factor(self.bed_diameter, 1)
        iscale.set_scaling_factor(self.service_flow_rate, 0.1)
        iscale.set_scaling_factor(self.ebct, 1e-2)
        iscale.set_scaling_factor(self.loading_rate, 1e3)

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


# [ESR WIP: The ion_exchange_multicomponent model (Clark) incorporates
# the trapezoidal rule within the model. For now, we will comment this
# out and evaluate whether to retain it in the base model.]

# def add_ss_approximation(blk, ix_model_type=None):

#     prop_in = blk.process_flow.properties_in[0]
#     blk.num_traps = slf.config.number_traps
#     blk.trap_disc = range(blk.num_traps + 1)
#     blk.trap_index = blk.trap_disc[1:]

#     blk.c_trap_min = Param(  # TODO: make CONFIG option
#         initialize=0.01,
#         default=0.01,
#         mutable=True,
#         units=pyunits.dimensionless,
#         doc="Minimum relative breakthrough concentration for estimating area under curve",
#     )

#     blk.c_traps = Var(
#         blk.reactive_component_set,
#         blk.trap_disc,
#         initialize=0.5,
#         bounds=(0, 1),
#         units=pyunits.dimensionless,
#         doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
#     )

#     blk.tb_traps = Var(
#         blk.reactive_component_set,
#         blk.trap_disc,
#         initialize=1e6,
#         bounds=(0, None),
#         units=pyunits.second,
#         doc="Breakthrough times for estimating area under breakthrough curve",
#     )

#     blk.traps = Var(
#         blk.reactive_component_set,
#         blk.trap_index,
#         initialize=0.01,
#         bounds=(0, 1),
#         units=pyunits.dimensionless,
#         doc="Trapezoid areas for estimating area under breakthrough curve",
#     )

#     for c in blk.reactive_component_set:
#         blk.c_traps[(c, 0)].fix(0)
#         blk.tb_traps[(c, 0)].fix(0)

#     blk.c_norm_avg = Var(
#         blk.reactive_component_set,
#         initialize=0.25,
#         bounds=(0, 2),
#         units=pyunits.dimensionless,
#         doc="Sum of trapezoid areas",
#     )

#     if ix_model_type == "cphsdm":

#         blk.min_tb_traps = Var(
#             blk.trap_index,
#             initialize=1e8,
#             bounds=(0, None),
#             units=pyunits.s,
#             doc="Minimum operational time of the bed by discrete element",
#         )

#         blk.throughput_traps = Var(
#             blk.trap_index,
#             initialize=1,
#             bounds=(0, None),
#             units=pyunits.dimensionless,
#             doc="Specific throughput of the bed by discrete element",
#         )

#         @blk.Constraint(
#             blk.trap_index,
#             doc="Minimum operational time of the bed from fresh to achieve a constant pattern solution by discrete element",
#         )
#         def eq_min_tb_traps(b, k):
#             return (
#                 b.min_tb_traps[k]
#                 == b.min_t_contact * (b.solute_dist_param + 1) * b.throughput_traps[k]
#             )

#         if blk.config.cphsdm_calculation_method == "surrogate":

#             blk.throughput_surrogate_traps = SurrogateBlock(
#                 blk.trap_index, concrete=True
#             )
#             for j in blk.reactive_component_set:
#                 for tr in blk.trap_index:
#                     blk.throughput_surrogate_traps[tr].build_model(
#                         blk.throughput_surrogate,
#                         input_vars=[
#                             blk.freundlich_ninv,
#                             blk.N_Bi,
#                             blk.c_traps[j, tr],
#                         ],
#                         output_vars=[blk.throughput_traps[tr]],
#                     )

#         elif blk.config.cphsdm_calculation_method == "input":

#             @blk.Constraint(
#                 blk.reactive_component_set,
#                 blk.trap_index,
#                 doc="Throughput constraint based on empirical 5-parameter regression by discretized element",
#             )
#             def eq_throughput_traps(b, j, k):
#                 return b.throughput_traps[k] == b.b0 + b.b1 * (
#                     b.c_traps[j, k] ** b.b2
#                 ) + b.b3 / (1.01 - b.c_traps[j, k] ** b.b4)

#     @blk.Constraint(
#         blk.reactive_component_set,
#         blk.trap_index,
#         doc="Evenly spaced c_norm for trapezoids",
#     )
#     def eq_c_traps(b, j, k):
#         return b.c_traps[j, k] == smooth_max(
#             1e-5,
#             b.c_trap_min
#             + (b.trap_disc[k] - 1) * ((b.c_norm[j] - b.c_trap_min) / (b.num_traps - 1)),
#         )

#     @blk.Constraint(
#         blk.reactive_component_set,
#         blk.trap_index,
#         doc="Breakthrough time calc for trapezoids",
#     )
#     def eq_tb_traps(b, j, k):
#         if ix_model_type == "clark":
#             bv_traps = (b.tb_traps[j, k] * b.loading_rate) / b.bed_depth
#             left_side = (
#                 (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
#                 / (b.bv_50[j] * b.loading_rate)
#             ) * (b.bv_50[j] - bv_traps)

#             right_side = log(
#                 ((1 / b.c_traps[j, k]) ** (b.freundlich_n[j] - 1) - 1)
#                 / (2 ** (b.freundlich_n[j] - 1) - 1)
#             )
#             return left_side - right_side == 0
#         elif ix_model_type == "cphsdm":
#             return b.tb_traps[j, k] == b.min_tb_traps[k] + (
#                 b.t_contact - b.min_t_contact
#             ) * (b.solute_dist_param + 1)
#         else:
#             return Constraint.Skip

#     @blk.Constraint(
#         blk.reactive_component_set, blk.trap_index, doc="Area of trapezoids"
#     )
#     def eq_traps(b, j, k):
#         return b.traps[j, k] == (b.tb_traps[j, k] - b.tb_traps[j, k - 1]) / b.tb_traps[
#             j, b.num_traps
#         ] * ((b.c_traps[j, k] + b.c_traps[j, k - 1]) / 2)

#     @blk.Constraint(
#         blk.reactive_component_set, doc="Average relative effluent concentration"
#     )
#     def eq_c_norm_avg(b, j):
#         return b.c_norm_avg[j] == sum(b.traps[j, k] for k in b.trap_index)

#     @blk.Constraint(
#         blk.reactive_component_set,
#         doc="CV mass transfer term",
#     )
#     def eq_mass_transfer_term(b, j):
#         return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
#             "Liq", j
#         ) == -b.process_flow.mass_transfer_term[0, "Liq", j]

#     @blk.Constraint(blk.reactive_component_set, doc="Regeneration stream mass flow")
#     def eq_mass_transfer_regen(b, j):
#         return (
#             b.regeneration_stream[0].get_material_flow_terms("Liq", j)
#             == -b.process_flow.mass_transfer_term[0, "Liq", j]
#         )

#     @blk.Constraint(doc="Regeneration stream flow rate")
#     def eq_regen_flow_rate(b):
#         return b.regeneration_stream[0].flow_vol_phase["Liq"] == pyunits.convert(
#             b.rinse_flow_rate * (b.rinse_time / b.cycle_time)
#             + b.backwash_flow_rate * (b.backwash_time / b.cycle_time)
#             + b.regen_flow_rate * (b.regeneration_time / b.cycle_time),
#             to_units=pyunits.m**3 / pyunits.s,
#         )
