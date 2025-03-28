#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Two Salt Diafiltration Unit Model
=================================

Author: Molly Dougher

This multi-component model for the diafiltration membrane model is a two salt system with a common
anion. The membrane is designed for use in the diafiltration cascade, i.e., is a spiral-wound
membrane module.

Configuration Arguments
-----------------------

The two-salt diafiltration model requires a property package that includes the moles of dissociated
ions in solution, as well as the valency, molar mass, and reflection coefficient of each ion in 
solution. 

Additionally, there are two required arguments, ``NFEx`` and ``NFEz``, to specfiy the desired number 
of finite elements across the width and thickness of the membrane, respectively.

Model Structure
---------------

There are three phases in the two-salt diafiltration model: the retentate, the membrane, and the
permeate. The retentate and the permeate are only discretized with respect to ``x``, while the
membrane is discretized with respect to both ``x`` and ``z``. The resulting system of partial
differential algebraic equations is solved by discretizing with the forward finite element method.

Assumptions
-----------

The partition coefficients of all ions in solution are equal to 1, meaning the concentration of an 
ion ``i`` just outside the membrane is equal to the concentration of ion ``i`` just inside the 
membrane. We also assume that the membrane has no surface charged, which is a valid assumption for
NF90 membranes.
    
The formation of a boundary layer at the membrane surface due to concentration polarization is 
neglected for mathematical simplicity.

The dominating transport mechanism within the bulk/retentate solution is convection in the 
``x``-direction (parallel to the membrane surface). The dominating transport mechanism within the 
permeate solution is convection in the ``z``-direction (perpendicular to the membrane surface).

The transport mechanisms modeled within the membrane are convection, diffusion, and electromigration.
Diffusion within the membrane that is normal to the pore walls is ignored, meaning the concentration 
gradient of ion ``i`` within the membrane only has a ``z``-component (perpendicular to the membrane 
surface).

Sets
----

The two-salt diafiltration model defines the following discrete sets for ions in the system.

.. math:: \mathcal{I}=\{\mathrm{Li^+,Co^{2+},Cl^-}\}

There are 2 continuous sets for each of length dimension: the ``x``-direction parallel to the membrane
surface and the ``z``-direction perpendicular to the membrane surface. ``x`` and ``z`` are non-
dimensionalized with the membrane width (``w``) and thickness (``l``), respectively, to improve numerics.

.. math:: \bar{x} \in \mathbb{R} \; | \; 0 \leq \bar{x} \geq 1

.. math:: \bar{z} \in \mathbb{R} \; | \; 0 \leq \bar{z} \geq 1

Parameters
----------

============================== ============================================= ================================= ============== =========================
Parameter                      Description                                    Name                             Default Value  Units
============================== ============================================= ================================= ============== =========================
:math:`c_{\mathrm{Co^{2+}},d}` concentration of cobalt in the diafiltrate    ``diafiltrate_conc_mass_cobalt``   0.2           :math:`kg m^{-3}` 
:math:`c_{\mathrm{Li^+},d}`    concentration of lithium in the diafiltrate   ``diafiltrate_conc_mass_lithium``  0.1           :math:`kg m^{-3}` 
:math:`c_{\mathrm{Co^{2+}},f}` concentration of cobalt in the feed           ``feed_conc_mass_cobalt``          17            :math:`kg m^{-3}` 
:math:`c_{\mathrm{Li^+},f}`    concentration of lithium in the feed          ``feed_conc_mass_lithium``         1.7           :math:`kg m^{-3}` 
:math:`l`                      thickness of the membrane                     ``membrane_thickness``             1e-07         :math:`m`
:math:`L`                      length of the membrane                        ``membrane_length``                100           :math:`m`
:math:`L_p`                    hydraulic permeability of the membrane        ``membrane_permeability``          0.01          :math:`m h^{-1} bar^{-1}`
:math:`\delta P`               applied pressure to the membrane              ``applied_pressure``               10            :math:`bar`
:math:`q_d`                    volumetic flow rate of the diafiltrate        ``diafiltrate_flow_volume``        30            :math:`m^3 h^{-1}`
:math:`q_f`                    volumetic flow rate of the feed               ``feed_flow_volume``               100           :math:`m^3 h^{-1}`
:math:`T`                      temperature of the system                     ``temperature``                    298           :math:`K`
:math:`w`                      width of the membrane                         ``membrane_width``                 1             :math:`m`
============================== ============================================= ================================= ============== =========================

Variables
---------

==================================== =========================================== =============================== ======================== ================================
Variable                             Description                                 Name                            Units                    Notes
==================================== =========================================== =============================== ======================== ================================
:math:`c_{\mathrm{Cl^-},m}`          concentration of chlorine in the membrane   ``membrane_conc_mass_chlorine`` :math:`kg m^{-3}`        discretized over ``x`` and ``z``
:math:`c_{\mathrm{Co^{2+}},m}`       concentration of cobalt in the membrane     ``membrane_conc_mass_cobalt``   :math:`kg m^{-3}`        discretized over ``x`` and ``z``
:math:`c_{\mathrm{Li^+},m}`          concentration of lithium in the membrane    ``membrane_conc_mass_lithium``  :math:`kg m^{-3}`        discretized over ``x`` and ``z``
:math:`c_{\mathrm{Cl^-},p}`          concentration of chlorine in the membrane   ``permeate_conc_mass_chlorine`` :math:`kg m^{-3}`        discretized over ``x``
:math:`c_{\mathrm{Co^{2+}},p}`       concentration of cobalt in the permeate     ``permeate_conc_mass_cobalt``   :math:`kg m^{-3}`        discretized over ``x``
:math:`c_{\mathrm{Li^+},p}`          concentration of lithium in the permeate    ``permeate_conc_mass_lithium``  :math:`kg m^{-3}`        discretized over ``x``
:math:`c_{\mathrm{Cl^-},r}`          concentration of chlorine in the retentate  ``retentate_conc_mass_chlorine``:math:`kg m^{-3}`        discretized over ``x``
:math:`c_{\mathrm{Co^{2+}},r}`       concentration of cobalt in the retentate    ``retentate_conc_mass_cobalt``  :math:`kg m^{-3}`        discretized over ``x``
:math:`c_{\mathrm{Li^+},r}`          concentration of lithium in the retentate   ``retentate_conc_mass_lithium`` :math:`kg m^{-3}`        discretized over ``x``
:math:`D_{\mathrm{Li^+,Li^+}}`       cross-diffusion coefficient (Li,Li)         ``D_lithium_lithium``           :math:`m^2 h^{-1}`       discretized over ``x`` and ``z``
:math:`D_{\mathrm{Li^+,Co^{2+}}}`    cross-diffusion coefficient (Li,Co)         ``D_lithium_cobalt``            :math:`m^2 h^{-1}`       discretized over ``x`` and ``z``
:math:`D_{\mathrm{Co^{2+},Li^+}}`    cross-diffusion coefficient (Co,Li)         ``D_cobalt_lithium``            :math:`m^2 h^{-1}`       discretized over ``x`` and ``z``
:math:`D_{\mathrm{Co^{2+},Co^{2+}}}` cross-diffusion coefficient (Co,Co)         ``D_cobalt_cobalt``             :math:`m^2 h^{-1}`       discretized over ``x`` and ``z``
:math:`j_{\mathrm{Cl^-}}`            mass flux of chlorine across the membrane   ``mass_flux_chlorine``          :math:`kg m^{-2} h^{-1}` discretized over ``x``
:math:`j_{\mathrm{Co^{2+}}}`         mass flux of cobalt across the membrane     ``mass_flux_cobalt``            :math:`kg m^{-2} h^{-1}` discretized over ``x``
:math:`j_{\mathrm{Li^+}}`            mass flux of lithium across the membrane    ``mass_flux_lithium``           :math:`kg m^{-2} h^{-1}` discretized over ``x``
:math:`J_w`                          water flux across the membrane `            ``volume_flux_water``           :math:`m^3 m^{-2} h^{-1}`discretized over ``x``
:math:`\delta \pi`                   osmotic pressure of feed-side fluid         ``osmotic_pressure``            :math:`bar`              discretized over ``x``
:math:`q_p`                          volumetic flow rate of the permeate         ``permeate_flow_volume``        :math:`m^3 h^{-1}`       discretized over ``x``
:math:`q_r`                          volumetic flow rate of the retentate        ``retentate_flow_volume``       :math:`m^3 h^{-1}`       discretized over ``x``
==================================== =========================================== =============================== ======================== ================================

Derivative Variables
--------------------

================================================================== =========================================== ==================================== ================== ================================
Variable                                                           Description                                 Name                                 Units              Notes
================================================================== =========================================== ==================================== ================== ================================
:math:`\frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}}{\mathrm{d}\bar{x}}` cobalt concentration gradient in retentate  ``d_retentate_conc_mass_cobalt_dx``  :math:`kg m^{-3}`  discretized over ``x``
:math:`\frac{\mathrm{d}c_{\mathrm{Li^+},r}}{\mathrm{d}\bar{x}}`    lithium concentration gradient in retentate ``d_retentate_conc_mass_lithium_dx`` :math:`kg m^{-3}`  discretized over ``x``
:math:`\frac{\mathrm{d}q_r}{\mathrm{d}\bar{x}}`                    retentate flow rate gradient                ``d_retentate_flow_volume_dx``       :math:`m^3 h^{-1}` discretized over ``x``
:math:`\frac{\partial c_{\mathrm{Co^{2+}},m}}{\partial \bar{z}}`   cobalt concentration gradient in membrane   ``d_membrane_conc_mass_cobalt_dz``   :math:`kg m^{-3}`  discretized over ``x`` and ``z``
:math:`\frac{\partial c_{\mathrm{Li^+},m}}{\partial \bar{z}}`      lithium concentration gradient in membrane  ``d_membrane_conc_mass_lithium_dz``  :math:`kg m^{-3}`  discretized over ``x`` and ``z``
================================================================== =========================================== ==================================== ================== ================================

Constraints
-----------

Will be added pending finalization of model and developing automated workflow of converting latex
document to sphinx format

"""

from pyomo.common.config import ConfigBlock, ConfigValue
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    Suffix,
    TransformationFactory,
    Var,
    units,
)

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants


@declare_process_block_class("TwoSaltDiafiltration")
class TwoSaltDiafiltrationData(UnitModelBlockData):
    """
    Two Salt Diafiltration Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for membrane system",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}
""",
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
see property package for documentation}
""",
        ),
    )
    CONFIG.declare(
        "NFEx",
        ConfigValue(
            doc="Number of discretization points in the x-direction",
        ),
    )
    CONFIG.declare(
        "NFEz",
        ConfigValue(
            doc="Number of discretization points in the z-direction",
        ),
    )

    def build(self):
        """
        Build method for the two salt diafiltration unit model
        """
        super().build()

        # TODO: generalize to any 2 cations and 1 anion

        self.add_mutable_parameters()
        self.add_variables()
        self.add_constraints()
        self.discretize_model()
        self.add_scaling_factors()

    def add_mutable_parameters(self):
        """
        Adds parameters for the two salt diafiltration unit model

        Assigns default values that can be changed by the user during implementation
        """
        self.membrane_thickness = Param(
            initialize=1e-7,
            mutable=True,
            units=units.m,
            doc="Thickness of membrane (z-direction)",
        )
        self.membrane_width = Param(
            initialize=1,
            mutable=True,
            units=units.m,
            doc="Width of the membrane (x-direction)",
        )
        self.membrane_length = Param(
            initialize=100,
            mutable=True,
            units=units.m,
            doc="Length of the membrane, wound radially",
        )
        self.applied_pressure = Param(
            initialize=10,
            mutable=True,
            units=units.bar,
            doc="Pressure applied to membrane",
        )
        self.membrane_permeability = Param(
            initialize=0.01,
            mutable=True,
            units=units.m / units.h / units.bar,
            doc="Hydraulic permeability coefficient",
        )
        self.feed_flow_volume = Param(
            initialize=100,
            mutable=True,
            units=units.m**3 / units.h,
            doc="Volumetric flow rate of the feed",
        )
        self.feed_conc_mass_lithium = Param(
            initialize=1.7,
            mutable=True,
            units=units.kg / units.m**3,
            doc="Mass concentration of lithium in the feed",
        )
        self.feed_conc_mass_cobalt = Param(
            initialize=17,
            mutable=True,
            units=units.kg / units.m**3,
            doc="Mass concentration of cobalt in the feed",
        )
        self.diafiltrate_flow_volume = Param(
            initialize=30,
            mutable=True,
            units=units.m**3 / units.h,
            doc="Volumetric flow rate of the diafiltrate",
        )
        self.diafiltrate_conc_mass_lithium = Param(
            initialize=0.1,
            mutable=True,
            units=units.kg / units.m**3,
            doc="Mass concentration of lithium in the diafiltrate",
        )
        self.diafiltrate_conc_mass_cobalt = Param(
            initialize=0.2,
            mutable=True,
            units=units.kg / units.m**3,
            doc="Mass concentration of cobalt in the diafiltrate",
        )
        self.temperature = Param(
            initialize=298,
            mutable=True,
            units=units.K,
            doc="System temperature",
        )

    def add_variables(self):
        """
        Adds variables for the two salt diafiltration unit model
        """
        # define length scales
        self.x_bar = ContinuousSet(bounds=(0, 1))
        self.z_bar = ContinuousSet(bounds=(0, 1))

        ## dependent on x_hat
        self.volume_flux_water = Var(
            self.x_bar,
            initialize=0.03,
            units=units.m**3 / units.m**2 / units.h,
            domain=NonNegativeReals,
            doc="Volumetric water flux of water across the membrane",
        )
        self.mass_flux_lithium = Var(
            self.x_bar,
            initialize=0.05,
            units=units.kg / units.m**2 / units.h,
            domain=NonNegativeReals,
            doc="Mass flux of lithium across the membrane (z-direction, x-dependent)",
        )
        self.mass_flux_cobalt = Var(
            self.x_bar,
            initialize=0.05,
            units=units.kg / units.m**2 / units.h,
            domain=NonNegativeReals,
            doc="Mass flux of cobalt across the membrane (z-direction, x-dependent)",
        )
        self.mass_flux_chlorine = Var(
            self.x_bar,
            initialize=0.05,
            units=units.kg / units.m**2 / units.h,
            domain=NonNegativeReals,
            doc="Mass flux of chlorine across the membrane (z-direction, x-dependent)",
        )
        self.retentate_flow_volume = Var(
            self.x_bar,
            initialize=130,
            units=units.m**3 / units.h,
            domain=NonNegativeReals,
            doc="Volumetric flow rate of the retentate, x-dependent",
        )
        self.retentate_conc_mass_lithium = Var(
            self.x_bar,
            initialize=1.33,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of lithium in the retentate, x-dependent",
        )
        self.retentate_conc_mass_cobalt = Var(
            self.x_bar,
            initialize=13.1,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of cobalt in the retentate, x-dependent",
        )
        self.retentate_conc_mass_chlorine = Var(
            self.x_bar,
            initialize=22.5,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of chlorine in the retentate, x-dependent",
        )
        self.permeate_flow_volume = Var(
            self.x_bar,
            initialize=1e-8,
            units=units.m**3 / units.h,
            domain=NonNegativeReals,
            doc="Volumetric flow rate of the permeate, x-dependent",
        )
        self.permeate_conc_mass_lithium = Var(
            self.x_bar,
            initialize=1e-10,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of lithium in the permeate, x-dependent",
        )
        self.permeate_conc_mass_cobalt = Var(
            self.x_bar,
            initialize=1e-10,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of cobalt in the permeate, x-dependent",
        )
        self.permeate_conc_mass_chlorine = Var(
            self.x_bar,
            initialize=1e-10,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of chlorine in the retentate, x-dependent",
        )
        self.osmotic_pressure = Var(
            self.x_bar,
            initialize=5,
            units=units.bar,
            domain=NonNegativeReals,
            doc="Osmostic pressure of the feed-side fluid",
        )

        ## dependent on z_hat and x_hat
        self.membrane_conc_mass_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=1.7,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of lithium in the membrane, x- and z-dependent",
        )
        self.membrane_conc_mass_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=13.1,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of cobalt in the membrane, x- and z-dependent",
        )
        self.membrane_conc_mass_chlorine = Var(
            self.x_bar,
            self.z_bar,
            initialize=22.5,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of chlorine in the membrane, x- and z-dependent",
        )
        self.D_lithium_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=-1e-9,
            units=units.m**2 / units.h,
            doc="Linearized cross diffusion coefficient for lithium-lithium",
        )
        self.D_lithium_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=-1e-11,
            units=units.m**2 / units.h,
            doc="Linearized cross diffusion coefficient for lithium-cobalt",
        )
        self.D_cobalt_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=-1e-10,
            units=units.m**2 / units.h,
            doc="Linearized cross diffusion coefficient for cobalt-lithium",
        )
        self.D_cobalt_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=-1e-9,
            units=units.m**2 / units.h,
            doc="Linearized cross diffusion coefficient for cobalt-cobalt",
        )

        # define the (partial) derivative variables
        self.d_retentate_conc_mass_lithium_dx = DerivativeVar(
            self.retentate_conc_mass_lithium,
            wrt=self.x_bar,
            units=units.kg / units.m**3,
        )
        self.d_retentate_conc_mass_cobalt_dx = DerivativeVar(
            self.retentate_conc_mass_cobalt,
            wrt=self.x_bar,
            units=units.kg / units.m**3,
        )
        self.d_retentate_flow_volume_dx = DerivativeVar(
            self.retentate_flow_volume,
            wrt=self.x_bar,
            units=units.m**3 / units.h,
        )
        self.d_membrane_conc_mass_lithium_dz = DerivativeVar(
            self.membrane_conc_mass_lithium,
            wrt=self.z_bar,
            units=units.kg / units.m**3,
        )
        self.d_membrane_conc_mass_cobalt_dz = DerivativeVar(
            self.membrane_conc_mass_cobalt,
            wrt=self.z_bar,
            units=units.kg / units.m**3,
        )

    def add_constraints(self):
        # mass balance constraints
        def _overall_mass_balance(self, x):
            return self.d_retentate_flow_volume_dx[x] == (
                -self.volume_flux_water[x] * self.membrane_length * self.membrane_width
            )

        self.overall_mass_balance = Constraint(self.x_bar, rule=_overall_mass_balance)

        def _lithium_mass_balance(self, x):
            return (
                self.retentate_flow_volume[x] * self.d_retentate_conc_mass_lithium_dx[x]
            ) == (
                (
                    self.volume_flux_water[x] * self.retentate_conc_mass_lithium[x]
                    - self.mass_flux_lithium[x]
                )
                * self.membrane_length
                * self.membrane_width
            )

        self.lithium_mass_balance = Constraint(self.x_bar, rule=_lithium_mass_balance)

        def _cobalt_mass_balance(self, x):
            return (
                self.retentate_flow_volume[x] * self.d_retentate_conc_mass_cobalt_dx[x]
            ) == (
                (
                    self.volume_flux_water[x] * self.retentate_conc_mass_cobalt[x]
                    - self.mass_flux_cobalt[x]
                )
                * self.membrane_length
                * self.membrane_width
            )

        self.cobalt_mass_balance = Constraint(self.x_bar, rule=_cobalt_mass_balance)

        def _general_mass_balance_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_lithium[x] * self.retentate_flow_volume[x]
                + self.permeate_conc_mass_lithium[x] * self.permeate_flow_volume[x]
            ) == (
                self.feed_flow_volume * self.feed_conc_mass_lithium
                + self.diafiltrate_flow_volume * self.diafiltrate_conc_mass_lithium
            )

        self.general_mass_balance_lithium = Constraint(
            self.x_bar, rule=_general_mass_balance_lithium
        )

        def _general_mass_balance_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_cobalt[x] * self.retentate_flow_volume[x]
                + self.permeate_conc_mass_cobalt[x] * self.permeate_flow_volume[x]
            ) == (
                self.feed_flow_volume * self.feed_conc_mass_cobalt
                + self.diafiltrate_flow_volume * self.diafiltrate_conc_mass_cobalt
            )

        self.general_mass_balance_cobalt = Constraint(
            self.x_bar, rule=_general_mass_balance_cobalt
        )

        # transport constraints (geometric)
        def _geometric_flux_equation_overall(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.permeate_flow_volume[x]
                == self.volume_flux_water[x]
                * x
                * self.membrane_length
                * self.membrane_width
            )

        self.geometric_flux_equation_overall = Constraint(
            self.x_bar, rule=_geometric_flux_equation_overall
        )

        def _geometric_flux_equation_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return self.mass_flux_lithium[x] == (
                self.permeate_conc_mass_lithium[x] * self.volume_flux_water[x]
            )

        self.geometric_flux_equation_lithium = Constraint(
            self.x_bar, rule=_geometric_flux_equation_lithium
        )

        def _geometric_flux_equation_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return self.mass_flux_cobalt[x] == (
                self.permeate_conc_mass_cobalt[x] * self.volume_flux_water[x]
            )

        self.geometric_flux_equation_cobalt = Constraint(
            self.x_bar, rule=_geometric_flux_equation_cobalt
        )

        # transport constraints (first principles)
        def _lumped_water_flux(self, x):
            if x == 0:
                return Constraint.Skip
            return self.volume_flux_water[x] == (
                self.membrane_permeability
                * (self.applied_pressure - self.osmotic_pressure[x])
            )

        self.lumped_water_flux = Constraint(self.x_bar, rule=_lumped_water_flux)

        def _D_lithium_lithium_calculation(self, x, z):
            return self.D_lithium_lithium[x, z] == (
                (-3.87e-6 * units.m**2 / units.h)
                + (
                    (-6.56e-8 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_lithium[x, z])
                )
                + (
                    (2.58e-8 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_cobalt[x, z])
                )
            )

        self.D_lithium_lithium_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_lithium_lithium_calculation
        )

        def _D_lithium_cobalt_calculation(self, x, z):
            return self.D_lithium_cobalt[x, z] == (
                (-4.50e-7 * units.m**2 / units.h)
                + (
                    (-1.70e-7 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_lithium[x, z])
                )
                + (
                    (6.67e-8 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_cobalt[x, z])
                )
            )

        self.D_lithium_cobalt_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_lithium_cobalt_calculation
        )

        def _D_cobalt_lithium_calculation(self, x, z):
            return self.D_cobalt_lithium[x, z] == (
                (-6.47e-7 * units.m**2 / units.h)
                + (
                    (4.10e-8 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_lithium[x, z])
                )
                + (
                    (-1.61e-8 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_cobalt[x, z])
                )
            )

        self.D_cobalt_lithium_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_cobalt_lithium_calculation
        )

        def _D_cobalt_cobalt_calculation(self, x, z):
            return self.D_cobalt_cobalt[x, z] == (
                (-3.56e-6 * units.m**2 / units.h)
                + (
                    (3.91e-7 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_lithium[x, z])
                )
                + (
                    (-1.53e-7 * units.m**5 / units.kg / units.h)
                    * (self.membrane_conc_mass_cobalt[x, z])
                )
            )

        self.D_cobalt_cobalt_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_cobalt_cobalt_calculation
        )

        def _lithium_flux_membrane(self, x, z):
            if z == 0:
                return Constraint.Skip

            return self.mass_flux_lithium[x] == (
                self.membrane_conc_mass_lithium[x, z] * self.volume_flux_water[x]
                + (
                    self.D_lithium_lithium[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mass_lithium_dz[x, z]
                )
                + (
                    self.D_lithium_cobalt[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mass_cobalt_dz[x, z]
                )
            )

        self.lithium_flux_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_lithium_flux_membrane
        )

        def _cobalt_flux_membrane(self, x, z):
            if z == 0:
                return Constraint.Skip

            return self.mass_flux_cobalt[x] == (
                self.membrane_conc_mass_cobalt[x, z] * self.volume_flux_water[x]
                + (
                    self.D_cobalt_lithium[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mass_lithium_dz[x, z]
                )
                + (
                    self.D_cobalt_cobalt[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mass_cobalt_dz[x, z]
                )
            )

        self.cobalt_flux_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_cobalt_flux_membrane
        )

        def _chlorine_flux_membrane(self, x):
            return self.mass_flux_chlorine[x] == -(
                (
                    self.config.property_package.charge["Li"]
                    / self.config.property_package.charge["Cl"]
                )
                * (
                    self.config.property_package.molar_mass["Cl"]
                    / self.config.property_package.molar_mass["Li"]
                )
                * self.mass_flux_lithium[x]
            ) - (
                (
                    self.config.property_package.charge["Co"]
                    / self.config.property_package.charge["Cl"]
                )
                * (
                    self.config.property_package.molar_mass["Cl"]
                    / self.config.property_package.molar_mass["Co"]
                )
                * self.mass_flux_cobalt[x]
            )

        self.chlorine_flux_membrane = Constraint(
            self.x_bar, rule=_chlorine_flux_membrane
        )

        # other physical constraints
        def _osmotic_pressure_calculation(self, x):
            return self.osmotic_pressure[x] == units.convert(
                (
                    (
                        self.config.property_package.num_solutes
                        * Constants.gas_constant  # J / mol / K
                        * self.temperature
                    )
                    * (
                        self.config.property_package.sigma["Li"]
                        / self.config.property_package.molar_mass["Li"]
                        * (
                            self.retentate_conc_mass_lithium[x]
                            - self.permeate_conc_mass_lithium[x]
                        )
                        + self.config.property_package.sigma["Co"]
                        / self.config.property_package.molar_mass["Co"]
                        * (
                            self.retentate_conc_mass_cobalt[x]
                            - self.permeate_conc_mass_cobalt[x]
                        )
                        + self.config.property_package.sigma["Cl"]
                        / self.config.property_package.molar_mass["Cl"]
                        * (
                            self.retentate_conc_mass_chlorine[x]
                            - self.permeate_conc_mass_chlorine[x]
                        )
                    )
                ),
                to_units=units.bar,
            )

        self.osmotic_pressure_calculation = Constraint(
            self.x_bar, rule=_osmotic_pressure_calculation
        )

        def _electroneutrality_retentate(self, x):
            return 0 == (
                self.config.property_package.charge["Li"]
                * self.retentate_conc_mass_lithium[x]
                / self.config.property_package.molar_mass["Li"]
                + self.config.property_package.charge["Co"]
                * self.retentate_conc_mass_cobalt[x]
                / self.config.property_package.molar_mass["Co"]
                + self.config.property_package.charge["Cl"]
                * self.retentate_conc_mass_chlorine[x]
                / self.config.property_package.molar_mass["Cl"]
            )

        self.electroneutrality_retentate = Constraint(
            self.x_bar, rule=_electroneutrality_retentate
        )

        def _electroneutrality_membrane(self, x, z):
            if z == 0:
                return Constraint.Skip
            return 0 == (
                self.config.property_package.charge["Li"]
                * self.membrane_conc_mass_lithium[x, z]
                / self.config.property_package.molar_mass["Li"]
                + self.config.property_package.charge["Co"]
                * self.membrane_conc_mass_cobalt[x, z]
                / self.config.property_package.molar_mass["Co"]
                + self.config.property_package.charge["Cl"]
                * self.membrane_conc_mass_chlorine[x, z]
                / self.config.property_package.molar_mass["Cl"]
            )

        self.electroneutrality_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_electroneutrality_membrane
        )

        # boundary conditions
        def _retentate_membrane_interface_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_lithium[x]
                == self.membrane_conc_mass_lithium[x, 0]
            )

        self.retentate_membrane_interface_lithium = Constraint(
            self.x_bar, rule=_retentate_membrane_interface_lithium
        )

        def _retentate_membrane_interface_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_cobalt[x]
                == self.membrane_conc_mass_cobalt[x, 0]
            )

        self.retentate_membrane_interface_cobalt = Constraint(
            self.x_bar, rule=_retentate_membrane_interface_cobalt
        )

        def _retentate_membrane_interface_chlorine(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_chlorine[x]
                == self.membrane_conc_mass_chlorine[x, 0]
            )

        self.retentate_membrane_interface_chlorine = Constraint(
            self.x_bar, rule=_retentate_membrane_interface_chlorine
        )

        def _membrane_permeate_interface_lithium(self, x):
            return (
                self.permeate_conc_mass_lithium[x]
                == self.membrane_conc_mass_lithium[x, 1]
            )

        self.membrane_permeate_interface_lithium = Constraint(
            self.x_bar, rule=_membrane_permeate_interface_lithium
        )

        def _membrane_permeate_interface_cobalt(self, x):
            return (
                self.permeate_conc_mass_cobalt[x]
                == self.membrane_conc_mass_cobalt[x, 1]
            )

        self.membrane_permeate_interface_cobalt = Constraint(
            self.x_bar, rule=_membrane_permeate_interface_cobalt
        )

        def _membrane_permeate_interface_chlorine(self, x):
            return (
                self.permeate_conc_mass_chlorine[x]
                == self.membrane_conc_mass_chlorine[x, 1]
            )

        self.membrane_permeate_interface_chlorine = Constraint(
            self.x_bar, rule=_membrane_permeate_interface_chlorine
        )

        # initial conditions
        def _initial_retentate_flow_volume(self):
            return self.retentate_flow_volume[0] == (
                self.feed_flow_volume + self.diafiltrate_flow_volume
            )

        self.initial_retentate_flow_volume = Constraint(
            rule=_initial_retentate_flow_volume
        )

        def _initial_retentate_conc_mass_lithium(self):
            return self.retentate_conc_mass_lithium[0] == (
                (
                    self.feed_flow_volume * self.feed_conc_mass_lithium
                    + self.diafiltrate_flow_volume * self.diafiltrate_conc_mass_lithium
                )
                / (self.feed_flow_volume + self.diafiltrate_flow_volume)
            )

        self.initial_retentate_conc_mass_lithium = Constraint(
            rule=_initial_retentate_conc_mass_lithium
        )

        def _initial_retentate_conc_mass_cobalt(self):
            return self.retentate_conc_mass_cobalt[0] == (
                (
                    self.feed_flow_volume * self.feed_conc_mass_cobalt
                    + self.diafiltrate_flow_volume * self.diafiltrate_conc_mass_cobalt
                )
                / (self.feed_flow_volume + self.diafiltrate_flow_volume)
            )

        self.initial_retentate_conc_mass_cobalt = Constraint(
            rule=_initial_retentate_conc_mass_cobalt
        )

        # set "zero" initial values to a sufficiently small value
        # flow rates: 1e-8
        # concentrations: 1e-10
        # derivatives: 1e-15
        def _initial_permeate_flow_volume(self):
            return self.permeate_flow_volume[0] == (1e-8 * units.m**3 / units.h)

        self.initial_permeate_flow_volume = Constraint(
            rule=_initial_permeate_flow_volume
        )

        def _initial_membrane_interface_lithium(self):
            return self.membrane_conc_mass_lithium[0, 0] == (
                1e-10 * units.kg / units.m**3
            )

        self.initial_membrane_interface_lithium = Constraint(
            rule=_initial_membrane_interface_lithium
        )

        def _initial_membrane_interface_cobalt(self):
            return self.membrane_conc_mass_cobalt[0, 0] == (
                1e-10 * units.kg / units.m**3
            )

        self.initial_membrane_interface_cobalt = Constraint(
            rule=_initial_membrane_interface_cobalt
        )

        def _initial_membrane_interface_chlorine(self):
            return self.membrane_conc_mass_chlorine[0, 0] == (
                1e-10 * units.kg / units.m**3
            )

        self.initial_membrane_interface_chlorine = Constraint(
            rule=_initial_membrane_interface_chlorine
        )

        def _initial_permeate_conc_mass_lithium(self):
            return self.permeate_conc_mass_lithium[0] == (1e-10 * units.kg / units.m**3)

        self.initial_permeate_conc_mass_lithium = Constraint(
            rule=_initial_permeate_conc_mass_lithium
        )

        def _initial_permeate_conc_mass_cobalt(self):
            return self.permeate_conc_mass_cobalt[0] == (1e-10 * units.kg / units.m**3)

        self.initial_permeate_conc_mass_cobalt = Constraint(
            rule=_initial_permeate_conc_mass_cobalt
        )

        def _initial_d_retentate_conc_mass_lithium_dx(self):
            return self.d_retentate_conc_mass_lithium_dx[0] == (
                1e-15 * units.kg / units.m**3
            )

        self.initial_d_retentate_conc_mass_lithium_dx = Constraint(
            rule=_initial_d_retentate_conc_mass_lithium_dx
        )

        def _initial_d_retentate_conc_mass_cobalt_dx(self):
            return self.d_retentate_conc_mass_cobalt_dx[0] == (
                1e-15 * units.kg / units.m**3
            )

        self.initial_d_retentate_conc_mass_cobalt_dx = Constraint(
            rule=_initial_d_retentate_conc_mass_cobalt_dx
        )

        def _initial_d_retentate_flow_volume_dx(self):
            return self.d_retentate_flow_volume_dx[0] == (1e-15 * units.m**3 / units.h)

        self._initial_d_retentate_flow_volume_dx = Constraint(
            rule=_initial_d_retentate_flow_volume_dx
        )

    def discretize_model(self):
        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self, wrt=self.x_bar, nfe=self.config.NFEx, scheme="FORWARD"
        )
        discretizer.apply_to(
            self, wrt=self.z_bar, nfe=self.config.NFEz, scheme="FORWARD"
        )

    def add_scaling_factors(self):
        """
        Apply scaling factors to certain constraints to improve solver performance

        Args:
            m: Pyomo model
        """
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add scaling factors for poorly scaled variables
        for x in self.x_bar:
            self.scaling_factor[self.retentate_flow_volume[x]] = 1e-2

            for z in self.z_bar:
                self.scaling_factor[self.D_lithium_lithium[x, z]] = 1e8
                self.scaling_factor[self.D_lithium_cobalt[x, z]] = 1e8
                self.scaling_factor[self.D_cobalt_lithium[x, z]] = 1e8
                self.scaling_factor[self.D_cobalt_cobalt[x, z]] = 1e8

        # Add scaling factors for poorly scaled constraints
        for x in self.x_bar:

            for z in self.z_bar:
                self.scaling_factor[self.D_lithium_lithium_calculation[x, z]] = 1e15
                self.scaling_factor[self.D_lithium_cobalt_calculation[x, z]] = 1e15
                self.scaling_factor[self.D_cobalt_lithium_calculation[x, z]] = 1e15
                self.scaling_factor[self.D_cobalt_cobalt_calculation[x, z]] = 1e15

                if z != 0:
                    self.scaling_factor[self.lithium_flux_membrane[x, z]] = 1e10
                    self.scaling_factor[self.cobalt_flux_membrane[x, z]] = 1e10

    # TODO: add ports
