#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This model was derived from:
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/ion_exchange_multi_component.py

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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

Ion Exchange Multicomponent (IXMC) Model
========================================

Author: E. Soraya Rawlings

The Ion Exchange Multicomponent (IXMC) model represents a fixed-bed
ion exchange process for the recovery of critical minerals and rare
earth elements (REEs) from diverse sources, including waste streams
generated during mining and industrial processes. The IXMC introduces
new features to address the complexities of multicomponent ion
exchange, enabling the prediction of process equilibrium,
hydrodynamics, bed and column geometry, and capital and operating
costs. The model extends the capabilities of the `WaterTAP Ion Exchange unit model <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/ion_exchange_0D.html>`_.

Configuration Arguments
-----------------------

The IXMC unit model requires the use of the `Multi-component aqueous
solution (MCAS) property package
<https://watertap.readthedocs.io/en/latest/technical_reference/property_models/mc_aq_sol.html>`_
. Users must specify the solute ions and provide molecular weight,
diffusivity, and charge data for each ion in the system.

The primary configuration arguments are:

#. ``property_package``: PhysicalParameterBlock for stream properties.
#. ``property_package_args``: Arguments for constructing the property package.
#. ``target_component``: The primary ion to be removed via ion exchange.
#. ``regenerant``: Regenerant configuration (e.g., ``single#.use``, ``NaCl``).
#. ``resin``: Resin name (e.g., ``S950``)
#. ``resin_data_path``: Path to the JSON file containing resin properties.
#. ``reactive_ions``: List of additional ions to be separated alongside the target ion.
#. ``number_of_trapezoids``: Number of trapezoids for trapezoidal rule calculations.
#. ``minimum_concentration_trapezoids``: Minimum concentration for ions to be used during trapezoidal rule.

Degrees of Freedom
------------------

In addition to specifying inlet feed state variables (temperature,
pressure, and molar flowrate for each ion :math:`j`), users must
provide additional degrees of freedom for both single-component and
multicomponent configurations to fully define the model. Common
variables to fix include:

#. ``resin_diam``: Diameter of the resin bead.
#. ``resin_density``: Bulk density of the resin.
#. ``bed_depth``: Height of the resin bed within column.
#. ``bed_diameter``: Diameter of the column.
#. ``bed_porosity``: Porosity of the resin bed.
#. ``number_columns``: Number of columns in operation.
#. ``loading_rate`` or ``ebct``: Superficial velocity through the resin bed or empty bed contact time.

Equilibrium-related variables to fix include:

#. ``freundlich_n_j``: Freundlich isotherm exponent for each ion :math:`j`
#. ``bv_50``: Bed volumes at fifty percent breakthrough for each ion :math:`j`
#. ``conc_comp_norm_breakthrough_j``: Normalized breakthrough concentration for each ion :math:`j`
#. ``mass_transfer_coeff_j``: Mass transfer coefficient for each ion :math:`j`

Model Structure
---------------

The IXMC model is organized into three distinct models, each designed
to address specific aspects of the ion exchange process:

1. **ion_exchange_multicomponent_base**: The base model that provides the
   foundational structure for the IXMC model, defining all global
   key variables and parameters, including column geometry and resin
   properties.

2. **ion_exchange_multicomponent**: The equilibrium model that
   incorporates the Freundlich multicomponent equations, using the
   Clark equation to describe the adsorption behavior of ions on the
   resin surface. It also applies all variables and equations for
   the trapezoidal rule calculations to determine breakthrough curve
   behavior and estimate effluent concentrations.

3. **ion_exchange_cost_block** and **ion_exchange_cost_model**: The
   costing models that implement costing variables and parameters
   using standard PrOMMiS costing features, ensuring compatibility
   with existing costing methodologies.

Assumptions
-----------

The model is implemented as a 0D control volume with the following
assumptions:

#. Single liquid phase
#. Steady state operation
#. Single and multiple solutes and one solvent (water)
#. Plug flow and isothermal conditions
#. Freundlich isotherm for equilibrium modeling
#. If ``single-use`` regenerant configuration is used, backwashing,
  regeneration, and rinsing steps are not modeled and associated costs
  are zero.
#. For multicomponent systems, competitive effects between ions are not
  currently modeled (see Limitations).

Sets
----

The model defines the following sets:

#. ``Components``: :math:`j` (e.g., ``['H2O', 'Cation_+', 'Anion_-', 'Inert']``)
#. ``Target Component``: :math:`j` (e.g., ``['Cation_+']``)
#. ``Reactive Ions``: :math:`j` (list of ions to be separated)
#. ``Inert``: :math:`j` (list of ions that are not separated)

Parameters
----------

The IXMC model has the following parameters:

.. csv-table:: Key operating and resin-specific parameters
   :header: "Parameter", "Description"

   ":math:`underdrain_h`", "Height of the column underdrain."
   ":math:`distributor_h`", "Height of the column distributor."
   ":math:`Pe_p_A`", "Peclet particle equation parameter A."
   ":math:`Pe_p_exp`", "Peclet particle equation exponent."
   ":math:`Sh_A`", "Sherwood equation parameter A."
   ":math:`Sh_exp_A`", "Sherwood equation exponent A."
   ":math:`Sh_exp_B`", "Sherwood equation exponent B."
   ":math:`Sh_exp_C`", "Sherwood equation exponent C."
   ":math:`pump_efficiency`", "Pump efficiency."
   ":math:`backwash_loading_rate`", "Backwash loading rate."
   ":math:`backwash_time`", "Backwash time."
   ":math:`pressure_drop_param_{A}`", "Pressure drop equation parameter A (resin-specific)."
   ":math:`pressure_drop_param_{B}`", "Pressure drop equation parameter B (resin-specific)."
   ":math:`pressure_drop_param_{C}`", "Pressure drop equation parameter C (resin-specific)."
   ":math:`bed_expansion_param_{A}`", "Bed expansion fraction equation parameter A (resin-specific)."
   ":math:`bed\_expansion\_frac\_B`", "Bed expansion fraction equation parameter B (resin-specific)."
   ":math:`bed\_expansion\_frac\_C`", "Bed expansion fraction equation parameter C (resin-specific)."
   ":math:`regen_dose`", "Regenerant dose required per volume of resin (if regenerant is not ``single_use``)."
   ":math:`regen_soln_conc`", "Concentration of regenerant solution (if regenerant is not ``single_use``)."
   ":math:`regen_soln_conc_sat`", "Saturation concentration of regenerant solution (if regenerant is not ``single_use``)."
   ":math:`regen_soln_dens`", "Density of regenerant solution (if regenerant is not ``single_use``)."
   ":math:`regen_cycle`", "Number of cycles the regenerant can be reused (if regenerant is not ``single_use``)."
   ":math:`num_regen_columns`", "Number of regeneration columns (if regenerant is not ``single_use``)."
   ":math:`rinse_bed_volumes`", "Number of bed volumes for the rinse step (if regenerant is not ``single_use``)."

Parameters from Resin Data
^^^^^^^^^^^^^^^^^^^^^^^^^^

Resin-specific parameters are unique to each resin and are typically
obtained from manufacturer data sheeta, which provide physical
properties and experimental plots of bed expansion and pressure drop
as functions of flow rate, temperature, and other operating
conditions. Curve-fitting techniques are commonly used to derive these
coefficients from the data. These parameters are stored in a JSON
file. The tables below summarize the available resin types and the key
properties included for each resin.

.. csv-table:: Resin names and their types included in the resin JSON file.
   :header: "Resin Name", "Resin Type"
   :name: resin_table

   "A850", "Strong-base Type I Acrylic Anion Exchange :math:`[1]`"
   "S950", "Chelating Polystyrene-Divinylbenzene :math:`[2]`"

.. csv-table:: Overview of the key properties included for each resin in the JSON file and their units.
   :header: "Property", "Units"

   "Functional Group", "N/A"
   "Bed Expansion Parameters (A, B, C)", ":math:`dimensionless`"
   "Pressure Drop Parameters (A, B, C)", ":math:`psi/m`"
   "Diameter", ":math:`m`"
   "Bulk Density", ":math:`kg/m^3`"
   "Porosity (takes a value from 0 to 1)", ":math:`dimensionless` "
   "Reference", "N/A"

Variables
---------

The IXMC model adds the following variables:

.. csv-table:: Key Parameters Needed in the IXMC Model
   :header: "Name", "Symbol", "Description"

    "``Inlet Temperature``", ":math:`T_{in}`", "Inlet temperature of the feed stream (:math:`\text{K}`)."
    "``Inlet Pressure``", ":math:`P_{in}`", "Inlet pressure of the feed stream (:math:`\text{Pa}`)."
    "``Component Molar Flow Rate``", ":math:`\dot{n}_{j}`", "Molar flow rate of each component (:math:`\text{mol/s}`)."
    "``Bed Linear Velocity``", ":math:`u_{bed}`", "Linear velocity through the resin bed (:math:`\text{m/s}`)."
    "``Interstitial Velocity``", ":math:`u_{inter}`", "Interstitial velocity in the bed (:math:`\text{m/s}`)."
    "``Regenerant Recycle Cycles``", ":math:`regen_{recycle}`", "Number of cycles before regenerant disposal."

    "``Resin Diameter``", ":math:`resin_diam`", "Diameter of individual resin beads (:math:`\text{m}`)."
    "``Resin Density``", ":math:`resin_density`", "Bulk density of the resin (:math:`\text{kg/L}`)."
    "``Bed Volume``", ":math:`bed_volume`", "Volume of resin bed in column (:math:`\text{m3}`)."
    "``Total Bed Volume``", ":math:`V_{bed,tot}`", "Total volume of the resin beds (:math:`\text{m}^3`)."
    "``Bed Depth``", ":math:`Z`", "Depth of the resin bed (:math:`\text{m}`)."
    "``Bed Porosity``", ":math:`\varepsilon`", "Porosity of the resin bed (dimensionless).    "``Column Height``", ":math:`H_{col}`", "Height of the column (:math:`\text{m}`)."
    "``Column Diameter``", ":math:`D_{col}`", "Diameter of the column (:math:`\text{m}`)."
"
    "``Number of Columns``", ":math:`n_{col}`", "Number of operational columns in the system."
    "``Target Component Breakthrough Time``", ":math:`target_breakthrough_time}`", "The time the target ion takes to appear at the outlet at the desired final concentration (:math:`\text{s}`)."
    "``Empty Bed Contact Time``", ":math:`ebct`", "The time the liquid remains in contact with the resin (:math:`\text{s}`)."
    "``Loading Rate``", ":math:`loading_rate`", "Superficial velocity through the resin bed (:math:`\text{m/s}`)."
    "``Service Flow Rate``", ":math:`service_flow_rate`", "Service flow rate in the column (:math:`1 / \text{hr}`)."
    "``Reynolds Number``", ":math:`N_Re`", "Reynolds number (dimensionless)."
    "``Peclet Particle Number``", ":math:`N_Pe_particle`", "Peclet particle number (dimensionless)."
    "``Peclet Bed Number``", ":math:`N_Pe_bed`", "Peclet bed number (dimensionless)."
    "``Bed Area``", ":math:`bed_area`", "The cross-sectional area of the resin bed, calculated based on the column dimensions."
    "``Breakthrough Time``", ":math:`target_breakthrough_time}`", "The time each reactive ion :math:`j` takes to appear at the outlet at the desired final concentration (:math:`\text{s}`)."
    "``Sherwood Number``", ":math:`N_Sh`", "Sherwood number (dimensionless)."
    "``Normalized Concentration at Final Breakthrough``", ":math:`conc_comp_norm_breakthrough`", "Normalized final breakthrough concentration for each ion :math:`j` (dimensionless)."
    "``Normalized Concentrations for Trapezoids``", ":math:`conc_comp_norm_trapezoids_j`", "Normalized breakthrough concentrations for each ion :math:`j` for each trapezoid (dimensionless)."
    "``Breakthrough Time for Trapezoids``", ":math:`breakthrough_time_trapezoids}`", "The time each reactive ion :math:`j` takes to appear at each trapezoid (:math:`\text{s}`)."
    "``Trapezoids``", ":math:`trapezoids}`", "The trapezoid areas for each ion :math:`j` for estimatinf the area under the breakthrough curve (dimensionless)."
    "``Freundlich Coefficient n``", ":math:`{freundlich_n_j}`", "The Freundlich isotherm coefficient :math:`n` that characterizes the adsorption capacity of each ion :math:`j` in the solution."
    "``Mass Transfer Coefficient``", ":math:`k_{T_j}`", "The coefficient that quantifies the rate of mass transfer for each ion in the solution."
    "``Bed Volumes at 50% Influent Concentration``", ":math:`BV_{{50}_j}`", "The volume of influent required to reach 50% of the initial concentration for each ion :math:`j` in the solution."

These variables are calculated using constraints and expressions
defined within the models (see Model Structure above).

Limitations of the IXMC Model
-----------------------------

1. **Steady-State Approximation**: The IXMC model uses a steady-state
   continuous approximation to represent the ion exchange batch
   process. This simplified model will not capture transient dynamics
   or time-dependent behavior observed in real-world ion exchange
   operations.

2. **Ion Competitiveness**: The IXMC multicomponent model currently
   does not account for ion competitiveness for the resin sites during
   the adsorption stage. These competitive effects include
   interactions between ions and the resin.

References
----------

| [1] Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). MWH's Water Treatment, Chapter 16: Ion Exchange, Figures 16-14 and 16-15. John Wiley & Sons, Inc.
| [2] Resin Product data sheet: URL: https://www.purolite.com/product/mts9500
| [3] LeVan, M. D., Carta, G., & Yon, C. M. (2019). Section 16: Adsorption and Ion Exchange. Perry's Chemical Engineers' Handbook, 9th Edition.
| [4] Inamuddin, & Luqman, M. (2012). Ion Exchange Technology I: Theory and Materials.
| [5] Inglezakis, V. J., & Poulopoulos, S. G. (2006). Adsorption, Ion Exchange and Catalysis: Design of Operations and Environmental Applications. doi:10.1016/B978-0-444-52783-7.X5000-9
| [6] Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle, R., Oppenheimer, J., & Jacangelo, J. G. (2023). Fundamental kinetic constants for breakthrough of per| and polyfluoroalkyl substances at varying empty bed contact times: Theoretical analysis and pilot scale demonstration. Chemical Engineering Journal, 464. doi:10.1016/j.cej.2023.142587
| [7] United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Models. https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models

"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue
from pyomo.common.config import ListOf

# Import IDAES cores
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError

from prommis.ion_exchange.ion_exchange_multicomponent_base import (
    IonExchangeBaseData,
    IonExchangeType,
)

""" This model contains the Clark model equations using Freundlich
equilibrium equations. This is a modified version of multicomponent
model from the WaterTAP developers:
https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/ion_exchange_multi_component.py

REFERENCES

[1] Vassilis J. Inglezakis and Stavros G. Poulopoulos Adsorption, Ion
Exchange and Catalysis: Design of Operations and Environmental
Applications (2006).  doi.org/10.1016/B978-0-444-52783-7.X5000-9

[2] Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle,
R., Oppenheimer, J., & Jacangelo, J. G. (2023).  Fundamental kinetic
constants for breakthrough of per- and polyfluoroalkyl substances at
varying empty bed contact times: Theoretical analysis and pilot scale
demonstration.  Chemical Engineering Journal,
464. doi:10.1016/j.cej.2023.142587

[3] Clark, R. M. (1987).  Evaluating the cost and performance of
field-scale granular activated carbon systems.  Environ Sci Technol,
21(6), 573-580. doi:10.1021/es00160a008


"""

__author__ = "Kurban Sitterley, Soraya Rawlings"


@declare_process_block_class("IonExchangeMultiComp")
class IonExchangeMultiCompData(IonExchangeBaseData):
    """
    Ion exchange multi component Clark model
    """

    CONFIG = IonExchangeBaseData.CONFIG()

    CONFIG.declare(
        "reactive_ions",
        ConfigValue(
            default=list(),
            domain=ListOf(str),
            description="Designates other reactive species",
        ),
    )

    CONFIG.declare(
        "number_of_trapezoids",
        ConfigValue(
            default=5,
            domain=int,
            description="Designates number of trapezoids",
        ),
    )

    CONFIG.declare(
        "minimum_concentration_trapezoids",
        ConfigValue(
            default=1e-3,
            domain=float,
            description="Minimum relative breakthrough concentration for estimating area under curve",
        ),
    )

    def build(self):

        # Validation that the flowsheet has one single time point
        t0 = self.flowsheet().time.first()

        super().build()

        prop_in = self.process_flow.properties_in[t0]
        p0 = list(prop_in.phase_list)[0]

        # [ESR WIP: Add regen terms when regeneration is needed]
        if self.config.regenerant != "single_use":
            regen = self.regeneration_stream[t0]

        comps = self.config.property_package.component_list
        target_component = self.config.target_component
        reactive_ions = self.config.reactive_ions

        assert target_component in reactive_ions
        self.target_component_set = pyo.Set(initialize=[target_component])

        # [ESR WIP: Define if inerts are needed in future
        # cases. Commented for now since it is not used in our current
        # case.]
        self.reactive_ion_set = pyo.Set(initialize=reactive_ions)
        self.inert_set = pyo.Set(initialize=comps - self.reactive_ion_set)

        if len(self.target_component_set) > 1:
            raise ConfigurationError(
                f"IonExchangeMultiComp can only accept a single target ion but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        else:
            # [ESR WIP: The current ion exchange model is used for the
            # separation of rare earth elements (REEs), which implies
            # a cation exchange type. We removed other types for
            # now. We will consider if we need to integrate the
            # alternative types in future iterations.]
            raise ConfigurationError(
                f"The current ion exchange model is limited to cation exchange methods, but the target component {target_component} has a charge of {self.config.property_package.charge_comp[target_component].value}."
            )

        # [ESR update: Change the inerts for the ones in the set.]
        for j in self.inert_set:
            self.process_flow.mass_transfer_term[:, p0, j].fix(0)

            if self.config.regenerant != "single_use":
                regen.get_material_flow_terms(p0, j).fix(0)

        # [ESR WIP: Define terms for trapezoidal rule. NOTE: the
        # trap_disc is a discretization index/parameter that defines
        # how the range between the minimum normalized concentration
        # (min_conc_comp_norm_trapezoids) and the final normalized
        # breakthrough concentration (conc_comp_norm_breakthrough) is
        # broken up.]
        self.number_of_trapezoids = self.config.number_of_trapezoids
        self.trap_disc = range(self.number_of_trapezoids + 1)
        self.trap_index = self.trap_disc[1:]

        min_conc_comp_norm_trapezoids_init = {}
        for i in self.reactive_ion_set:
            min_conc_comp_norm_trapezoids_init[i] = (
                self.config.minimum_concentration_trapezoids
            )
        self.min_conc_comp_norm_trapezoids = pyo.Param(
            self.reactive_ion_set,
            mutable=True,
            initialize=min_conc_comp_norm_trapezoids_init,
            doc="Minimum normalized concentration for each reactive ion species",
        )

        # [ESR WIP: Bring breakthrough time and bv here from base model
        # since these variables depend on the set of reactive
        # ions. Also, add a constraint to calculate the maximum
        # breakthrough time. For now, this maximum time is for the
        # target component. Think if there is a better way to do
        # this.]
        self.breakthrough_time = pyo.Var(
            self.reactive_ion_set,
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            domain=pyo.NonNegativeReals,
            units=pyo.units.s,
            doc="Breakthrough time at final breakthrough normalized concentration",
        )

        @self.Constraint()
        def eq_target_breakthru_time(b):
            return (
                self.target_breakthrough_time
                == self.breakthrough_time[target_component]
            )

        # Add variables for dimensionless number
        self.N_Sh = pyo.Var(
            self.reactive_ion_set,
            initialize=30,
            units=pyo.units.dimensionless,
            doc="Sherwood number",
        )

        # Add relevant variables for IX column and trapezoidal rule
        self.conc_comp_norm_breakthrough = pyo.Var(
            self.reactive_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Normalized final breakthrough concentration [Ct/C0] of reactive ions",
        )

        self.conc_comp_norm_trapezoids = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Normalized breakthrough concentrations for each trapezoid for estimating area under breakthrough curve",
        )

        self.breakthrough_time_trapezoids = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=1,
            bounds=(0, None),
            units=pyo.units.second,
            doc="Breakthrough times for each trapezoid for estimating area under breakthrough curve",
        )

        self.trapezoids = pyo.Var(
            self.reactive_ion_set,
            self.trap_index,
            initialize=0.01,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Trapezoid areas for estimating area under breakthrough curve",
        )

        for c in self.reactive_ion_set:
            self.conc_comp_norm_trapezoids[(c, t0)].fix(0)
            self.breakthrough_time_trapezoids[(c, t0)].fix(0)

        self.conc_comp_norm_breakthrough_trapezoids = pyo.Var(
            self.reactive_ion_set,
            initialize=0.25,
            bounds=(0, 2),
            units=pyo.units.dimensionless,
            doc="Summation of normalized concentrations of all trapezoid areas for each reactive ion",
        )

        self.freundlich_n = pyo.Var(
            self.reactive_ion_set,
            initialize=1.5,
            bounds=(0, None),
            units=pyo.units.dimensionless,
            doc="Freundlich isotherm exponent",
        )

        self.mass_transfer_coeff = pyo.Var(  # k_T
            self.reactive_ion_set,
            initialize=0.001,
            units=pyo.units.s**-1,
            bounds=(0, None),
            doc="Mass transfer coefficient for Clark model (kT)",
        )

        self.bv_50 = pyo.Var(  # BV_50
            self.reactive_ion_set,
            initialize=2e5,
            bounds=(0, None),
            units=pyo.units.dimensionless,
            doc="Bed volumes of feed at 50 percent breakthrough",
        )

        @self.Expression(self.reactive_ion_set, doc="Breakthrough concentration")
        def breakthrough_conc_mass_comp(b, j):
            return (
                b.conc_comp_norm_breakthrough[j] * prop_in.conc_mass_phase_comp[p0, j]
            )

        # Add Expression/Constraint to calculate dimensionless numbers
        @self.Expression(self.reactive_ion_set, doc="Schmidt number")
        def N_Sc(b, j):  # Eq. 3.359, ref[1] Inglezakis + Poulopoulos
            return prop_in.visc_k_phase[p0] / prop_in.diffus_phase_comp[p0, j]

        @self.Constraint(self.reactive_ion_set, doc="Sherwood number")
        def eq_Sh(b, j):  # Eq. 3.346, ref[1] Inglezakis + Poulopoulos
            return (
                b.N_Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.N_Re**b.Sh_exp_B
                * b.N_Sc[j] ** b.Sh_exp_C
            )

        # [ESR updates: Add an if statement to consider this
        # constraint ONLY for the target component since it was adding
        # 1 degree of freedom when solving for multiple ions. [TODO:
        # Review this change once we add the regeneration stage.]
        if self.config.regenerant != "single_use":

            @self.Constraint(
                self.reactive_ion_set, doc="Mass transfer for regeneration stream"
            )
            def eq_mass_transfer_regen(b, j):
                return (
                    regen.get_material_flow_terms(p0, j)
                    == -b.process_flow.mass_transfer_term[t0, p0, j]
                )

        # [ESR updates: Add multicomponent set for breakthrough_time
        # and bv]
        @self.Expression(self.reactive_ion_set, doc="Bed volumes at breakthrough")
        def bv(b, r):
            return b.breakthrough_time[r] * b.loading_rate / b.bed_depth

        @self.Constraint(
            self.reactive_ion_set, doc="Clark equation with fundamental constants"
        )  # ref[2], Croll et al (2023), Eq.9
        def eq_clark(b, j):
            left_side = pyo.exp(
                (
                    (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
                    / (b.bv_50[j] * b.loading_rate)
                )
                * (b.bv_50[j] - b.bv[j])
            )
            right_side = (
                (1 / b.conc_comp_norm_breakthrough[j]) ** (b.freundlich_n[j] - 1) - 1
            ) / (2 ** (b.freundlich_n[j] - 1) - 1)
            return left_side - right_side == 0

        # [ESR notes: This is the concentration value at the kth
        # “trapezoid” boundary for ion j.]
        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Evenly spaced conc_comp_norm_breakthrough for trapezoids",
        )
        def eq_conc_comp_norm_trapezoids(b, j, k):
            if k == max(b.trap_index):
                return (
                    b.conc_comp_norm_trapezoids[j, k]
                    == b.conc_comp_norm_breakthrough[j]
                )
            else:
                # [ESR updates: Reformulated to avoid denominators.]
                return b.conc_comp_norm_trapezoids[j, k] * (
                    b.number_of_trapezoids - 1
                ) == (
                    b.min_conc_comp_norm_trapezoids[j] * (b.number_of_trapezoids - 1)
                    + (b.trap_disc[k] - 1)
                    * (
                        b.conc_comp_norm_breakthrough[j]
                        - b.min_conc_comp_norm_trapezoids[j]
                    )
                )

        # [ESR updates: Add bv_trapezoids as an expression. This term
        # corresponds to the bed volume corresponding to the kth
        # discretization point for ion j.]
        @self.Expression(self.reactive_ion_set, self.trap_disc, doc="BV for trapezoids")
        def bv_trapezoids(b, j, k):
            return (b.breakthrough_time_trapezoids[j, k] * b.loading_rate) / b.bed_depth

        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Breakthru time calc for trapezoids",
        )
        def eq_breakthrough_time_trapezoids(b, j, k):
            if k == max(self.trap_index):
                return b.breakthrough_time_trapezoids[j, k] == b.breakthrough_time[j]
            else:
                left_side = pyo.exp(
                    (
                        (
                            b.mass_transfer_coeff[j]
                            * b.bed_depth
                            * (b.freundlich_n[j] - 1)
                        )
                        / (b.bv_50[j] * b.loading_rate)
                    )
                    * (b.bv_50[j] - b.bv_trapezoids[j, k])
                )
                right_side = (
                    (1 / b.conc_comp_norm_trapezoids[j, k]) ** (b.freundlich_n[j] - 1)
                    - 1
                ) / (2 ** (b.freundlich_n[j] - 1) - 1)
                return left_side - right_side == 0

        # [ESR WIP: Try to make changes in the normalization term in
        # eq_trapezoids to see if that improves the breakthrough time
        # calculation. This normalization helps to express each
        # trapezoid’s width as a fraction of the total interval. This
        # normalization helps ensure that the integrated output (when
        # all trapezoidal areas are summed) gives a proper,
        # dimensionless measure of the overall breakthrough
        # behavior. NOTE: Re-write to avoid denominators]
        @self.Constraint(
            self.reactive_ion_set, self.trap_index, doc="Area of trapezoids"
        )
        def eq_trapezoids(b, j, k):
            # Norm term is the normalization term

            # Original normalization value
            norm_term = b.breakthrough_time_trapezoids[j, self.number_of_trapezoids]

            # # Use total interval value for normalization
            # norm_term = b.breakthrough_time_trapezoids[j, self.number_of_trapezoids] - b.breakthrough_time_trapezoids[j, 1]

            # No normalization
            # norm_term = 1

            return b.trapezoids[j, k] * norm_term == (
                # width of the trapezoid
                b.breakthrough_time_trapezoids[j, k]
                - b.breakthrough_time_trapezoids[j, k - 1]
            ) * (
                # average height between the two points
                (
                    b.conc_comp_norm_trapezoids[j, k]
                    + b.conc_comp_norm_trapezoids[j, k - 1]
                )
                / 2
            )

        @self.Constraint(
            self.reactive_ion_set,
            doc="Summation of relative effluent concentration for each trapezoid area",
        )
        def eq_conc_comp_norm_breakthrough_trapezoids(b, j):
            return b.conc_comp_norm_breakthrough_trapezoids[j] == sum(
                b.trapezoids[j, k] for k in b.trap_index
            )

        # [ESR WIP: Define total_mass_removed as an Expression.]
        @self.Expression(
            self.reactive_ion_set, doc="Total mass of ion removed from the liquid phase"
        )
        def total_mass_removed(b, j):
            return (
                1 - b.conc_comp_norm_breakthrough_trapezoids[j]
            ) * prop_in.get_material_flow_terms(p0, j)

        @self.Constraint(
            self.reactive_ion_set,
            doc="CV mass transfer term",
        )
        def eq_mass_transfer_reactive_ions(b, j):
            return (
                b.total_mass_removed[j] == -b.process_flow.mass_transfer_term[t0, p0, j]
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # [ESR updates: Add new variables from base model.]
        iscale.set_scaling_factor(self.breakthrough_time, 1e-4)
        iscale.set_scaling_factor(self.bv, 1e-1)
        iscale.set_scaling_factor(self.freundlich_n, 1)
        iscale.set_scaling_factor(self.mass_transfer_coeff, 1e2)
        iscale.set_scaling_factor(self.bv_50, 1e-3)

        sf = iscale.get_scaling_factor(self.breakthrough_time)
        iscale.set_scaling_factor(self.breakthrough_time_trapezoids, sf)

        iscale.set_scaling_factor(self.conc_comp_norm_trapezoids, 1)
        iscale.set_scaling_factor(self.trapezoids, 1e3)
        iscale.set_scaling_factor(self.conc_comp_norm_breakthrough_trapezoids, 1e2)

        for ind, c in self.eq_clark.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e-2)

        for ind, c in self.eq_trapezoids.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e2)
