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

The Ion Exchange Multicomponent (IXMC) model extends the WaterTAP Ion
Exchange (IX) unit model to improve the recovery of critical minerals
and rare earth elements (REEs) from multicomponent systems, including
waste streams generated during mining and industrial processes. The
IXMC model introduces new features to address the challenges of
multicomponent ion exchange processes while building upon the basic
assumptions and methodologies of the original IX model. Refer to the
`WaterTAP IX documentation <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/ion_exchange_0D.html>`_ for more information on the original IX model, including its
assumptions, equations, and implementation.


Key Features
------------
The IXMC model introduces new key features, including:

1. **Multicomponent separation**:

   - Supports the separation of multiple components by introducing a
     list of "reactive" ions that can be separated alongside the
     target ion. This allows for individual calculations of
     breakthrough times, composition profiles, and adsorption dynamics
     using the trapezoidal rule for each ion in the system.
   - Increases flexibility by providing individual control over final
     concentrations, making the model suited for systems where
     selective recovery is critical.

2. **Resin Manipulation**:

   - Allows users to select, define, and manipulate ion exchange
     resins, including specifying resin properties such as bead
     diameter, bulk density, and surface area per volume specific for
     each resin of interest.
   - Allows the modification of the resin file to incorporate new
     resins properties tailored to specific REEs or metal separations
     without altering the core equations of the IXMC model.

3. **Extended Freundlich Isotherm Implementation**:

   - Implements the Freundlich equilibrium model, which is well-suited
     for heterogeneous surfaces and non-ideal adsorption processes, to
     allow users to fit experimental breakthrough data and derive
     equilibrium parameters for multicomponent systems.

4. **New Configuration Options**:

   - Introduces new configuration arguments to the unit model to
     provide greater control over operational parameters, including
     initial concentrations, number of trapezoids for trapezoidal rule
     calculations, and the file path for resin properties.

5. **Regeneration Stream and Stage Separation**:

   - Distinguishes the exchange stage (adsorption) from the
     regeneration stage (desorption) using the configuration argument
     `regenerant`, allowing for greater flexibility in modeling
     the regeneration process.

6. **Redefined Dimensionless Numbers**:

   - Redefined dimensionless numbers and other key design terms as
     Pyomo `Expression` components, improving model flexibility,
     scalability, and convergence.


Introduction
------------

The IXMC model builds upon the original capabilities of the IX model
and retains the four primary operational steps of ion exchange
processes:

(1) Service
(2) Backwashing
(3) Regeneration
(4) Rinsing

Note that, even though the regeneration is added as a stream in the
IXMC model, it is currently only used for costing purposes. However,
as mentioned in the **Key Features** section (point 5), the IXMC model
separates all the terms required during this stage using the
configuration argument `regenerant`. When it is set to `single_use`,
the IXMC model excludes all equations related to the regeneration
stage. Conversely, when selecting a regenerant type, such as `NaCl`,
the regeneration stream (implemented as a Pyomo `port` component) and
associated regeneration equations are added to the model. This
capability provides greater flexibility in modeling ion exchange
process since it distinguishes between adsorption and desorption
processes. The integration of equilibrium equations to accurately
model both stages is part of the ongoing and future work.

Additionally, the IXMC model enhances multicomponent separation
capabilities by introducing a set of reactive ions (as mentioned in
the **Key Features** section in point 1), which can be separated
alongside the selected target ion. This feature enables precise
predictions of column performance, resin bed design characteristics,
and adsorption dynamics for multicomponent systems.


Model Structure
---------------

The IXMC model is separated into three distinct models, each
designed to address specific aspects of the ion exchange process:

1. **Base Model**:

   - Provides the foundational structure for the IXMC model, ensuring
     consistency across all calculations and enabling users to define
     global key variables and parameters.
   - Contains all global variables and parameters required for
     modeling the ion exchange column, including terms related to
     column geometry and resin properties. Most of these variables,
     parameters, and expressions are included under the **Model
     Components** Table in the `WaterTAP IX documentation
     <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/ion_exchange_0D.html>`_.

2. **Equilibrium Model**:

   - Adds the Freundlich multicomponent equations, which describe the
     adsorption behavior of ions on the resin surface.
   - Includes the equations for the trapezoidal rule to calculate
     compositions along the breakthrough curve. This numerical
     approach allows for accurate estimation of effluent
     concentrations and adsorption dynamics over time.

3. **Costing Block and Model**:

   - Modifies the original `IX costing package
     <https://watertap.readthedocs.io/en/latest/technical_reference/costing/ion_exchange.html>`_
     from WaterTAP. These models implement all costing variables and
     parameters using the standard `PrOMMiS costing features
     <https://prommis.readthedocs.io/en/latest/tutorials/costing_basic_features.html>`_
     to ensure compatibility with existing costing methodologies.

How to Use the IXMC Model
-------------------------

The IXMC model is designed to handle multicomponent systems, with a
focus on REE recovery. To solve an example, we follow three sequential steps:

1. **Resin-Specific Step**:

   - Ensure that the resin selected for the ion exchange process meets
     their specific requirements and that all resin-specific
     parameters and equations are accurately defined.
   - If the resin of interest is not available, analyze the resin data
     product sheet provided by the manufacturer to extract polynomial
     coefficients for bed expansion and pressure drop equations. These
     coefficients are essential for accurately modeling the
     operational performance of the column. For details on the current
     available resins and how we determine these coefficients, refer
     to the **Resin-Specific Information** section below.

2. **Parameter Estimation Step**:

   - Estimate equilibrium parameters using a parameter estimation
     model based on known breakthrough data. For our specific model
     under the Freundlich equilibrium, the estimated parameters are
     the Freundlich isotherm exponent (:math:`n`), mass transfer
     coefficient (:math:`k_T`), and bed volumes at 50% of the influent
     concentration (:math:`BV_{50}`).
   - Refer to the `Parmest documentation
     <https://pyomo.readthedocs.io/en/6.8.0/contributed_packages/parmest/index.html>`_,
     for a detailed guidance on parameter estimation.

3. **Data Implementation and Example Construction Step**:

   - Integrate all calculated and known parameters into the IXMC model
     to solve the system.
   - Customize variables, bounds, and separation conditions to match
     their specific case requirements. Additionally, adjust
     configuration options to define multicomponent property packages
     and specify target ions, ensuring the model accurately reflects
     the desired operational conditions. For assistance in defining
     your degrees of freedom, refer to the **Degrees of Freedom**
     section under Freundlich in the `WaterTAP IX documentation
     <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/ion_exchange_0D.html>`_.

This procedure mirrors the approach used in the original IX model,
with additional enhancements for multicomponent systems.


Resin-Specific Parameters and Equations
----------------------------------------

The IXMC model uses resin-specific parameters and equations to
calculate the pressure drop (:math:`p_{drop}`) and the bed expansion
fraction (:math:`H_{expan}`), considering the physical and chemical
properties of the resin. The IXMC model reads the resin parameters,
represented as polynomial coefficients :math:`p_{drop, A}, p_{drop,B},
p_{drop, C}, H_{expan, A}, H_{expan, B}` and :math:`H_{expan, C}`,
from a JSON file and calculates :math:`p_{drop}` and :math:`H_{expan}`
values using the equations in Table 1. Consider :math:`u_{bed}` and
:math:`u_{bw}` represent the linear velocity through the bed and the
backwashing loading rate, respectively. The equations in Table 1 allow
the IXMC model to accurately predict hydrodynamic behavior and
operational performance, ensuring compatibility with multiple resin
types.

.. csv-table:: Resin-Specific Equations for Pressure Drop and Bed Expansion.
   :header: "Description", "Equation"

   "Pressure drop (psi/m)", ":math:`p_{drop} = p_{drop, A} + p_{drop,B}u_{bed} + p_{drop, C}u_{bed}^{2}`"
   "Bed expansion fraction (dimensionless)", ":math:`H_{expan} = H_{expan, A} + H_{expan, B} u_{bw} + H_{expan, C} u_{bw}^{2}`"

The polynomial coefficients are specific to the resin being used and
are typically obtained from the resin data production sheet, which
include experimental plots of bed expansion and pressure drop as
functions of flow rate, temperature, and other operating
conditions. Curve-fitting techniques are often applied to derive these
coefficients from the data. Note that, to ensure accuracy, the output
of the equations should be validated against the data provided in the
resin product sheet. This validation step is essential for confirming
the reliability of the model predictions.

Resin Data in JSON File
^^^^^^^^^^^^^^^^^^^^^^^

Tables 2 and 3 summarize the key properties and information available
for the available resins. Note that most of the key values for these
properties can be found in the resin data sheets provided by the
manufacturer.

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


Modified Model Components
=========================

The IXMC model modifies variables, parameters, and expressions to
support multicomponent systems. Table 4 summarizes some of the key
terms that were modified in the model.

.. csv-table::  Key Parameters Needed in the IXMC Model
   :header: "Type", "Name", "Symbol", "Description"

   Variable, "``Empty Bed Contact Time``", ":math:`ebct`", "The time the liquid remains in contact with the resin, calculated as :math:`\frac{V_{res}}{flow_{col}}` where :math:`V_{res}` is the bed volume and :math:`flow_{col}` is the volumetric flow per column."
   Variable, "``Bed Area``", ":math:`bed_{area}`", "The cross-sectional area of the resin bed, calculated based on the column dimensions."
   Variable, "``Breakthrough Time``", ":math:`t_{{break}_j}`", "The time required for each ion :math:`j` to appear at the outlet at the desired final concentration."
   Variable, "``Freundlich Coefficient n``", ":math:`{n_j}`", "The Freundlich isotherm coefficient :math:`n` that characterizes the adsorption capacity of each ion :math:`j` in the solution."
   Variable, "``Bed Volumes at 50% Influent Concentration``", ":math:`BV_{{50}_j}`", "The volume of influent required to reach 50% of the initial concentration for each ion :math:`j` in the solution."
   Variable, "``Mass Transfer Coefficient``", ":math:`k_{T_j}`", "The coefficient that quantifies the rate of mass transfer for each ion in the solution."

Limitations of the IXMC Model
=============================

The IXMC model introduces new advancements for multicomponent ion
exchange processes but it is important to acknowledge its current
limitations and areas for future development. These limitations
highlight opportunities for further refinement and extension of the
model to improve its accuracy and applicability.

1. **Steady-State Approximation**: The IXMC model uses a steady-state continuous approximation to represent the ion exchange batch process. This simplification may not fully capture transient dynamics or time-dependent behavior observed in real-world ion exchange operations.

2. **Ion Competitiveness**: The IXMC multicomponent model currently does not account for ion competitiveness for the resin sites during the adsorption stage. These competitive effects include interactions between ions and the resin.

References
----------

| [1] Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012). MWH's Water Treatment, Chapter 16: Ion Exchange, Figures 16-14 and 16-15. John Wiley & Sons, Inc.

| [2] Resin Product data sheet: URL: https://www.purolite.com/product/mts9500


"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue

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


modified by: Soraya Rawlings

"""

__author__ = "Kurban Sitterley"


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
            domain=list,
            description="Designates other reactive species",
        ),
    )

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]

        # [ESR WIP: Add regen terms when regeneration is needed]
        if self.config.regenerant != "single_use":
            regen = self.regeneration_stream[0]

        comps = self.config.property_package.component_list
        target_component = self.config.target_component
        reactive_ions = self.config.reactive_ions

        assert target_component in reactive_ions
        self.target_component_set = pyo.Set(initialize=[target_component])

        # [ESR WIP: Define if inerts are needed in future
        # cases. Commented for now since it is not used in our current
        # case.]
        # inerts = comps - self.target_component_set
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
            # a cation exchange type. We commented other types for
            # now. We will consider if we need to integrate the
            # alternative types in future iterations.]
            raise ConfigurationError(
                "The current ion exchange model is limited to cation exchange methods and alternative techniques are not addressed at this time."
            )

        # elif self.config.property_package.charge_comp[target_component].value < 0:
        #     self.ion_exchange_type = IonExchangeType.anion
        # else:
        #     raise ConfigurationError("Target ion must have non-zero charge.")

        # [ESR update: Change the inerts for the ones in the set.]
        for j in self.inert_set:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)

            if self.config.regenerant != "single_use":
                regen.get_material_flow_terms("Liq", j).fix(0)

        # [ESR WIP: Define terms for trapezoidal rule. NOTE: the
        # trap_disc is a discretization index/parameter that defines
        # how the range between c_trap_min and c_norm is broken up.]
        self.num_traps = int(self.config.number_traps)
        self.trap_disc = range(self.num_traps + 1)
        self.trap_index = self.trap_disc[1:]

        self.eps = pyo.Param(initialize=1e-4, mutable=True)

        self.c_trap_min = {}
        for i in self.reactive_ion_set:
            self.c_trap_min[i] = float(self.config.c_trap_min)

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
            doc="Breakthrough time",
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
        self.c_norm = pyo.Var(
            self.reactive_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        )

        self.c_traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
        )

        self.tb_traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=1,
            bounds=(0, None),
            units=pyo.units.second,
            doc="Breakthrough times for estimating area under breakthrough curve",
        )

        self.traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_index,
            initialize=0.01,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Trapezoid areas for estimating area under breakthrough curve",
        )

        for c in self.reactive_ion_set:
            self.c_traps[(c, 0)].fix(0)
            self.tb_traps[(c, 0)].fix(0)

        self.c_norm_avg = pyo.Var(
            self.reactive_ion_set,
            initialize=0.25,
            bounds=(0, 2),
            units=pyo.units.dimensionless,
            doc="Sum of trapezoid areas",
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
        def c_breakthru(b, j):
            return b.c_norm[j] * prop_in.conc_mass_phase_comp["Liq", j]

        # Add Expression/Constraint to calculate dimensionless numbers
        @self.Expression(self.reactive_ion_set, doc="Schmidt number")
        def N_Sc(b, j):  # Eq. 3.359, ref[1] Inglezakis + Poulopoulos
            return prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]

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
                    regen.get_material_flow_terms("Liq", j)
                    == -b.process_flow.mass_transfer_term[0, "Liq", j]
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
            right_side = ((1 / b.c_norm[j]) ** (b.freundlich_n[j] - 1) - 1) / (
                2 ** (b.freundlich_n[j] - 1) - 1
            )
            return left_side - right_side == 0

        # [ESR notes: This is the concentration value at the kth
        # “trapezoid” boundary for ion j.]
        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Evenly spaced c_norm for trapezoids",
        )
        def eq_c_traps(b, j, k):
            if k == max(b.trap_index):
                return b.c_traps[j, k] == b.c_norm[j]
            else:
                # [ESR updates: Reformulated to avoid denominators.]
                return b.c_traps[j, k] * (b.num_traps - 1) == (
                    b.c_trap_min[j] * (b.num_traps - 1)
                    + (b.trap_disc[k] - 1) * (b.c_norm[j] - b.c_trap_min[j])
                )

        # [ESR updates: Add bv_traps as an expression. This term
        # corresponds to the bed volume corresponding to the kth
        # discretization point for ion j.]
        @self.Expression(self.reactive_ion_set, self.trap_disc, doc="BV for trapezoids")
        def bv_traps(b, j, k):
            return (b.tb_traps[j, k] * b.loading_rate) / b.bed_depth

        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Breakthru time calc for trapezoids",
        )
        def eq_tb_traps(b, j, k):
            if k == max(self.trap_index):
                return b.tb_traps[j, k] == b.breakthrough_time[j]
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
                    * (b.bv_50[j] - b.bv_traps[j, k])
                )
                right_side = ((1 / b.c_traps[j, k]) ** (b.freundlich_n[j] - 1) - 1) / (
                    2 ** (b.freundlich_n[j] - 1) - 1
                )
                return left_side - right_side == 0

        # [ESR WIP: Try to make changes in the normalization term in
        # eq_traps to see if that improves the breakthrough time
        # calculation. This normalization helps to express each
        # trapezoid’s width as a fraction of the total interval. This
        # normalization helps ensure that the integrated output (when
        # all trapezoidal areas are summed) gives a proper,
        # dimensionless measure of the overall breakthrough
        # behavior. NOTE: Re-write to avoid denominators]
        @self.Constraint(
            self.reactive_ion_set, self.trap_index, doc="Area of trapezoids"
        )
        def eq_traps(b, j, k):
            # Norm term is the normalization term

            # Original normalization value
            norm_term = b.tb_traps[j, self.num_traps]

            # # Use total interval value for normalization
            # norm_term = b.tb_traps[j, self.num_traps] - b.tb_traps[j, 1]

            # No normalization
            # norm_term = 1

            return b.traps[j, k] * norm_term == (
                # width of the trapezoid
                b.tb_traps[j, k]
                - b.tb_traps[j, k - 1]
            ) * (
                # average height between the two points
                (b.c_traps[j, k] + b.c_traps[j, k - 1])
                / 2
            )

        @self.Constraint(
            self.reactive_ion_set, doc="Average relative effluent concentration"
        )
        def eq_c_norm_avg(b, j):
            return b.c_norm_avg[j] == sum(b.traps[j, k] for k in b.trap_index)

        # [ESR WIP: Define total_mass_removed as an Expression.]
        @self.Expression(
            self.reactive_ion_set, doc="Total mass of ion removed from the liquid phase"
        )
        def total_mass_removed(b, j):
            return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms("Liq", j)

        @self.Constraint(
            # self.target_component_set,
            self.reactive_ion_set,
            doc="CV mass transfer term",
        )
        def eq_mass_transfer_reactive_ions(b, j):
            # return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
            #     "Liq", j
            # ) == -b.process_flow.mass_transfer_term[0, "Liq", j]
            return (
                b.total_mass_removed[j]
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
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
        iscale.set_scaling_factor(self.tb_traps, sf)

        iscale.set_scaling_factor(self.c_traps, 1)
        iscale.set_scaling_factor(self.traps, 1e3)
        iscale.set_scaling_factor(self.c_norm_avg, 1e2)

        for ind, c in self.eq_clark.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e-2)

        for ind, c in self.eq_traps.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e2)
