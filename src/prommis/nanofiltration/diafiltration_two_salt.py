#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Two-Salt Diafiltration Unit Model
=================================

Author: Molly Dougher

This multi-component model for the diafiltration membrane model is a two-salt system with a common anion. The membrane is designed for use in the diafiltration cascade, i.e., is a spiral-wound membrane module.

Configuration Arguments
-----------------------

The Two-Salt Diafiltration model requires a property package that includes the moles of dissociated ions in solution, as well as the valency, molar mass, and reflection coefficient of each ion in solution. 

There are three arguments with default values: the membrane width (``membrane_width``), the membrane length (``membrane_length``), and the applied pressure (`applied_pressure``). To run a simulation (with zero degrees of freedom), these variables must be fixed. Default values of :math:`1 m`, :math:`100 m`, and :math:`10 bar`, respectively, are provided if not specified by the user.

Additionally, there are two required arguments, ``NFEx`` and ``NFEz``, to specfiy the desired number of finite elements across the width and thickness of the membrane, respectively.

Degrees of Freedom
------------------

The Two-Salt Diafiltration unit model has three degrees of freedom: the width of the membrane (i.e., the length of the membrane module), the length of the membrane, and the pressure applied to the membrane system.

Model Structure
---------------

There are three phases in the Two-Salt Diafiltration model: the retentate, the membrane, and the permeate. The retentate and the permeate are only discretized with respect to :math:`x`, while the membrane is discretized with respect to both :math:`x` and :math:`z`. The resulting system of partial differential algebraic equations is solved by discretizing with the forward finite element method.

Assumptions
-----------

The partition coefficients of all ions in solution are equal to 1, meaning the concentration of an ion :math:`i` just outside the membrane is equal to the concentration of ion :math:`i` just inside the membrane. We also assume that the membrane has no surface charge, which is a valid assumption for NF90 membranes.

The formation of a boundary layer at the membrane surface due to concentration polarization is neglected for mathematical simplicity.

The dominating transport mechanism within the bulk/retentate solution is convection in the :math:`x`-direction (parallel to the membrane surface). The dominating transport mechanism within the permeate solution is convection in the :math:`z`-direction (perpendicular to the membrane surface).

The transport mechanisms modeled within the membrane are convection, diffusion, and electromigration. Diffusion within the membrane that is normal to the pore walls is ignored, meaning the concentration gradient of ion :math:`i` within the membrane only has a :math:`z`-component (perpendicular to the membrane surface).

Sets
----

The Two-Salt Diafiltration model defines the following discrete sets for ions in the system.

.. math:: \mathcal{I}=\{\mathrm{Li^+,Co^{2+},Cl^-}\}

There are 2 continuous sets for each length dimension: the :math:`x`-direction parallel to the membrane surface (width) and the :math:`z`-direction perpendicular to the membrane surface (thickness). :math:`x` and :math:`z` are non-dimensionalized with the membrane width (:math:`w`) and thickness (:math:`l`), respectively, to improve numerics.

.. math:: \bar{x} \in \mathbb{R} \| 0 \leq \bar{x} \leq 1
.. math:: \bar{z} \in \mathbb{R} \| 0 \leq \bar{z} \leq 1

Default Model Parameters
------------------------

The Two-Salt Diafiltration model has the following parameters.

============================== =============================================== ================================= ============= =========================
Parameter                      Description                                     Name                              Default Value Units
============================== =============================================== ================================= ============= =========================
:math:`c_{\mathrm{Co^{2+}},d}` concentration of cobalt ion in the diafiltrate  ``diafiltrate_conc_mass_cobalt``  0.2           :math:`kg m^{-3}` 
:math:`c_{\mathrm{Li^+},d}`    concentration of lithium ion in the diafiltrate ``diafiltrate_conc_mass_lithium`` 0.1           :math:`kg m^{-3}` 
:math:`c_{\mathrm{Co^{2+}},f}` concentration of cobalt ion in the feed         ``feed_conc_mass_cobalt``         17            :math:`kg m^{-3}` 
:math:`c_{\mathrm{Li^+},f}`    concentration of lithium ion in the feed        ``feed_conc_mass_lithium``        1.7           :math:`kg m^{-3}` 
:math:`l`                      thickness of the membrane                       ``membrane_thickness``            1e-07         :math:`m`
:math:`L_p`                    hydraulic permeability of the membrane          ``membrane_permeability``         0.01          :math:`m h^{-1} bar^{-1}`
:math:`q_d`                    volumetic flow rate of the diafiltrate          ``diafiltrate_flow_volume``       30            :math:`m^3 h^{-1}`
:math:`q_f`                    volumetic flow rate of the feed                 ``feed_flow_volume``              100           :math:`m^3 h^{-1}`
:math:`T`                      temperature of the system                       ``temperature``                   298           :math:`K`
============================== =============================================== ================================= ============= =========================

Variables
---------

The Two-Salt Diafiltration model adds the following variables.

==================================== ============================================== ================================ ========================= ========================================
Variable                             Description                                    Name                             Units                     Notes
==================================== ============================================== ================================ ========================= ========================================
:math:`c_{\mathrm{Cl^-},m}`          concentration of chlorine ion in the membrane  ``membrane_conc_mass_chlorine``  :math:`kg m^{-3}`         discretized over :math:`x` and :math:`z`
:math:`c_{\mathrm{Co^{2+}},m}`       concentration of cobalt ion in the membrane    ``membrane_conc_mass_cobalt``    :math:`kg m^{-3}`         discretized over :math:`x` and :math:`z`
:math:`c_{\mathrm{Li^+},m}`          concentration of lithium ion in the membrane   ``membrane_conc_mass_lithium``   :math:`kg m^{-3}`         discretized over :math:`x` and :math:`z`
:math:`c_{\mathrm{Cl^-},p}`          concentration of chlorine ion in the membrane  ``permeate_conc_mass_chlorine``  :math:`kg m^{-3}`         discretized over :math:`x`
:math:`c_{\mathrm{Co^{2+}},p}`       concentration of cobalt ion in the permeate    ``permeate_conc_mass_cobalt``    :math:`kg m^{-3}`         discretized over :math:`x`
:math:`c_{\mathrm{Li^+},p}`          concentration of lithium ion in the permeate   ``permeate_conc_mass_lithium``   :math:`kg m^{-3}`         discretized over :math:`x`
:math:`c_{\mathrm{Cl^-},r}`          concentration of chlorine ion in the retentate ``retentate_conc_mass_chlorine`` :math:`kg m^{-3}`         discretized over :math:`x`
:math:`c_{\mathrm{Co^{2+}},r}`       concentration of cobalt ion in the retentate   ``retentate_conc_mass_cobalt``   :math:`kg m^{-3}`         discretized over :math:`x`
:math:`c_{\mathrm{Li^+},r}`          concentration of lithium ion in the retentate  ``retentate_conc_mass_lithium``  :math:`kg m^{-3}`         discretized over :math:`x`
:math:`D_{\mathrm{Li^+,Li^+}}`       cross-diffusion coefficient (Li,Li)            ``D_lithium_lithium``            :math:`m^2 h^{-1}`        discretized over :math:`x` and :math:`z`
:math:`D_{\mathrm{Li^+,Co^{2+}}}`    cross-diffusion coefficient (Li,Co)            ``D_lithium_cobalt``             :math:`m^2 h^{-1}`        discretized over :math:`x` and :math:`z`
:math:`D_{\mathrm{Co^{2+},Li^+}}`    cross-diffusion coefficient (Co,Li)            ``D_cobalt_lithium``             :math:`m^2 h^{-1}`        discretized over :math:`x` and :math:`z`
:math:`D_{\mathrm{Co^{2+},Co^{2+}}}` cross-diffusion coefficient (Co,Co)            ``D_cobalt_cobalt``              :math:`m^2 h^{-1}`        discretized over :math:`x` and :math:`z`
:math:`j_{\mathrm{Cl^-}}`            mass flux of chlorine ion across the membrane  ``mass_flux_chlorine``           :math:`kg m^{-2} h^{-1}`  discretized over :math:`x`
:math:`j_{\mathrm{Co^{2+}}}`         mass flux of cobalt ion across the membrane    ``mass_flux_cobalt``             :math:`kg m^{-2} h^{-1}`  discretized over :math:`x`
:math:`j_{\mathrm{Li^+}}`            mass flux of lithium ion across the membrane   ``mass_flux_lithium``            :math:`kg m^{-2} h^{-1}`  discretized over :math:`x`
:math:`J_w`                          water flux across the membrane                 ``volume_flux_water``            :math:`m^3 m^{-2} h^{-1}` discretized over :math:`x`
:math:`L`                            length of the membrane                         ``membrane_length``              :math:`m`
:math:`\Delta \pi`                   osmotic pressure of feed-side fluid            ``osmotic_pressure``             :math:`bar`               discretized over :math:`x`
:math:`\Delta P`                     applied pressure to the membrane               ``applied_pressure``             :math:`bar`
:math:`q_p`                          volumetic flow rate of the permeate            ``permeate_flow_volume``         :math:`m^3 h^{-1}`        discretized over :math:`x`
:math:`q_r`                          volumetic flow rate of the retentate           ``retentate_flow_volume``        :math:`m^3 h^{-1}`        discretized over :math:`x`
:math:`w`                            width of the membrane                          ``membrane_width``               :math:`m`
==================================== ============================================== ================================ ========================= ========================================

Derivative Variables
--------------------

The Two-Salt Diafiltration model adds the following derivative variables.

================================================================== =============================================== ==================================== ================== ========================================
Variable                                                           Description                                     Name                                 Units              Notes
================================================================== =============================================== ==================================== ================== ========================================
:math:`\frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}}{\mathrm{d}\bar{x}}` cobalt ion concentration gradient in retentate  ``d_retentate_conc_mass_cobalt_dx``  :math:`kg m^{-3}`  discretized over :math:`x`
:math:`\frac{\mathrm{d}c_{\mathrm{Li^+},r}}{\mathrm{d}\bar{x}}`    lithium ion concentration gradient in retentate ``d_retentate_conc_mass_lithium_dx`` :math:`kg m^{-3}`  discretized over :math:`x`
:math:`\frac{\mathrm{d}q_r}{\mathrm{d}\bar{x}}`                    retentate flow rate gradient                    ``d_retentate_flow_volume_dx``       :math:`m^3 h^{-1}` discretized over :math:`x`
:math:`\frac{\partial c_{\mathrm{Co^{2+}},m}}{\partial \bar{z}}`   cobalt ion concentration gradient in membrane   ``d_membrane_conc_mass_cobalt_dz``   :math:`kg m^{-3}`  discretized over :math:`x` and :math:`z`
:math:`\frac{\partial c_{\mathrm{Li^+},m}}{\partial \bar{z}}`      lithium ion concentration gradient in membrane  ``d_membrane_conc_mass_lithium_dz``  :math:`kg m^{-3}`  discretized over :math:`x` and :math:`z`
================================================================== =============================================== ==================================== ================== ========================================

Constraints
-----------

Differential mass balances in the retentate:

.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}} = - J_w(\bar{x}) wL
.. math:: \frac{\mathrm{d}c_{\mathrm{Li^+},r}(\bar{x})}{\mathrm{d}\bar{x}} = \left( J_w(\bar{x}) c_{\mathrm{Li^+},r}(\bar{x}) - j_{\mathrm{Li^+}}(\bar{x}) \right) \frac{Lw}{q_r(\bar{x})}
.. math:: \frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}(\bar{x})}{\mathrm{d}\bar{x}} = \left( J_w(\bar{x}) c_{\mathrm{Co^{2+}},r}(\bar{x}) - j_{\mathrm{Co^{2+}}}(\bar{x}) \right) \frac{Lw}{q_r(\bar{x})}

Overall component mass balances:

.. math:: q_r(\bar{x})c_{\mathrm{Li^+},r}(\bar{x}) + q_p(\bar{x})c_{\mathrm{Li^+},p}(\bar{x}) = q_f c_{\mathrm{Li^+},f}+q_d c_{\mathrm{Li^+},d} \forall \bar{x} \neq 0
.. math:: q_r(\bar{x})c_{\mathrm{Co^{2+}},r}(\bar{x}) + q_p(\bar{x})c_{\mathrm{Co^{2+}},p}(\bar{x}) = q_f c_{\mathrm{Li^+},f}+q_d c_{\mathrm{Li^+},d} \forall \bar{x} \neq 0

Electroneutrality in the retentate:

.. math:: c_{\mathrm{Cl^-},r}(\bar{x}) = - \frac{z_{\mathrm{Li^+}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Li^+}}} c_{\mathrm{Li^+},r}(\bar{x}) - \frac{z_{\mathrm{Co^{2+}}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Co^{2+}}}} c_{\mathrm{Co^{2+}},r}(\bar{x})

Overall water flux through the membrane:

.. math:: J_w (\bar{x}) = L_p (\Delta P - \Delta \pi (\bar{x})) \forall \bar{x} \neq 0

Osmotic pressure:

.. math:: \Delta \pi (\bar{x}) =  n \mathrm{R} \mathrm{T} \sum_{i \in \mathcal{I}} \frac{\sigma_i}{MW_i} (c_{i,r}(\bar{x})-c_{i,p}(\bar{x}))

Nernst-Plank equations for the ion flux through the membrane:

.. math:: j_{\mathrm{Li^+}}(\bar{x}) = c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) J_w(\bar{x}) + \frac{D_{\mathrm{Li^+,Li^+}}}{l} \frac{\partial c_{\mathrm{Li^+},m}(\bar{x},\bar{z})}{\partial \bar{z}} + \frac{D_{\mathrm{Li^+,Co^{2+}}}}{l} \frac{\partial c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})}{\partial \bar{z}} \forall \bar{z} \neq 0
.. math:: j_{\mathrm{Co^{2+}}}(\bar{x}) = c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}) J_w(\bar{x}) + \frac{D_{\mathrm{Co^{2+},Li^+}}}{l} \frac{\partial c_{\mathrm{Li^+},m}(\bar{x},\bar{z})}{\partial \bar{z}} + \frac{D_{\mathrm{Co^{2+},Co^{2+}}}}{l} \frac{\partial c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})}{\partial \bar{z}} \forall \bar{z} \neq 0

with linearized cross-diffusion coefficients:

.. math:: D_{\mathrm{Li^+,Li^+}} = \beta_0 + \beta_1 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_2 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: D_{\mathrm{Li^+,Co^{2+}}} = \beta_3 + \beta_4 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_5 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: D_{\mathrm{Co^{2+},Li^+}} = \beta_6 + \beta_7 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_8 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: D_{\mathrm{Co^{2+},Co^{2+}}} = \beta_9 + \beta_{10} c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_{11} c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})

that have the following regressed parameter values:

================== ========= ================
Parameter          Value     Units
================== ========= ================
:math:`\beta_0`    -3.87e-06 :math:`m^2/h`
:math:`\beta_1`    -6.56e-08 :math:`m^5/kg/h`
:math:`\beta_2`    2.58e-08  :math:`m^5/kg/h`
:math:`\beta_3`    -4.50e-07 :math:`m^2/h`
:math:`\beta_4`    -1.70e-07 :math:`m^5/kg/h`
:math:`\beta_5`    6.67e-08  :math:`m^5/kg/h`
:math:`\beta_6`    -6.47e-07 :math:`m^2/h`
:math:`\beta_7`    4.10e-08  :math:`m^5/kg/h`
:math:`\beta_8`    -1.61e-08 :math:`m^5/kg/h`
:math:`\beta_9`    -3.56e-06 :math:`m^2/h`
:math:`\beta_{10}` 3.91e-07  :math:`m^5/kg/h`
:math:`\beta_{11}` -1.53e-07 :math:`m^5/kg/h`
================== ========= ================

No applied potential on the system:

.. math:: j_{\mathrm{Cl^-}}(\bar{x}) = - \frac{z_{\mathrm{Li^+}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Li^+}}} j_{\mathrm{Li^+}}(\bar{x}) - \frac{z_{\mathrm{Co^{2+}}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Co^{2+}}}} j_{\mathrm{Co^{2+}}}(\bar{x})

Electroneutrality in the membrane:

.. math:: c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}) = -\frac{z_{\mathrm{Li^+}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Li^+}}} c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) -\frac{z_{\mathrm{Co^{2+}}}}{z_{\mathrm{Cl^-}}} \frac{MW_{\mathrm{Cl^-}}}{MW_{\mathrm{Co^{2+}}}} c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}) \forall \bar{z} \neq 0

Mass balance (via convection) on the permeate outlet:

.. math:: q_p(\bar{x}) = \bar{x} wL J_w(\bar{x}) \forall \bar{x} \neq 0
.. math:: c_{\mathrm{Li^+},p}(\bar{x}) = \frac{j_{\mathrm{Li^+}}(\bar{x})}{J_w(\bar{x})} \forall \bar{x} \neq 0
.. math:: c_{\mathrm{Co^{2+}},p}(\bar{x}) = \frac{j_{\mathrm{Co^{2+}}}(\bar{x})}{J_w(\bar{x})} \forall \bar{x} \neq 0

Relate the concentrations at each phase interface:

.. math:: c_{\mathrm{Li^+},r}(\bar{x})=c_{\mathrm{Li^+},m}(\bar{x},\bar{z}=0)
.. math:: c_{\mathrm{Co^{2+}},r}(\bar{x})=c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}=0)
.. math:: c_{\mathrm{Cl^-},r}(\bar{x})=c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=0)
.. math:: c_{\mathrm{Li^+},m}(\bar{x},\bar{z}=1)=c_{\mathrm{Li^+},p}(\bar{x})
.. math:: c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}=1)=c_{\mathrm{Co^{2+}},p}(\bar{x})
.. math:: c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=1)=c_{\mathrm{Cl^-},p}(\bar{x})

Initial and boundary conditions:

.. math:: q_r(\bar{x}=0)=q_f+q_d
.. math:: q_r(\bar{x}=0) c_{\mathrm{Li^+},r}(\bar{x}=0)=q_f c_{\mathrm{Li^+},f} + q_d c_{\mathrm{Li^+},d}
.. math:: q_r(\bar{x}=0) c_{\mathrm{Co^{2+}},r}(\bar{x}=0)=q_f c_{\mathrm{Co^{2+}},f} + q_d c_{\mathrm{Co^{2+}},d}
.. math:: q_p(\bar{x}=0) = 0
.. math:: c_{\mathrm{Li^+},p} (\bar{x}=0) = 0
.. math:: c_{\mathrm{Co^{2+}},p} (\bar{x}=0) = 0
.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=0
.. math:: \frac{\mathrm{d}c_{\mathrm{Li^+},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=0
.. math:: \frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=0

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
    Two-Salt Diafiltration Unit Model Class.
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
**Valid values:** {see property package for documentation}
""",
        ),
    )
    CONFIG.declare(
        "membrane_width",
        ConfigValue(
            default=1,
            doc="Width of the membrane (x-direction) in meters",
        ),
    )
    CONFIG.declare(
        "membrane_length",
        ConfigValue(
            default=100,
            doc="Length of the membrane, wound radially in meters",
        ),
    )
    CONFIG.declare(
        "applied_pressure",
        ConfigValue(
            default=10,
            doc="Pressure applied to membrane in bar",
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
        Adds parameters for the two salt diafiltration unit model.

        Assigns default values that can be changed by the user during implementation.
        """
        self.membrane_thickness = Param(
            initialize=1e-7,
            mutable=True,
            units=units.m,
            doc="Thickness of membrane (z-direction)",
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
        Adds variables for the two salt diafiltration unit model.
        """
        # define length scales
        self.x_bar = ContinuousSet(bounds=(0, 1))
        self.z_bar = ContinuousSet(bounds=(0, 1))

        self.membrane_width = Var(
            initialize=self.config.membrane_width,
            units=units.m,
            domain=NonNegativeReals,
            doc="Width of the membrane (x-direction)",
        )
        self.membrane_length = Var(
            initialize=self.config.membrane_length,
            units=units.m,
            domain=NonNegativeReals,
            doc="Length of the membrane, wound radially",
        )
        self.applied_pressure = Var(
            initialize=self.config.applied_pressure,
            units=units.bar,
            domain=NonNegativeReals,
            doc="Pressure applied to membrane",
        )

        ## dependent on x_bar
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

        ## dependent on bar and x_bar
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
            doc="Lithium concentration gradient in the retentate",
        )
        self.d_retentate_conc_mass_cobalt_dx = DerivativeVar(
            self.retentate_conc_mass_cobalt,
            wrt=self.x_bar,
            units=units.kg / units.m**3,
            doc="Cobalt concentration gradient in the retentate",
        )
        self.d_retentate_flow_volume_dx = DerivativeVar(
            self.retentate_flow_volume,
            wrt=self.x_bar,
            units=units.m**3 / units.h,
            doc="Volume flow gradient in the retentate",
        )
        self.d_membrane_conc_mass_lithium_dz = DerivativeVar(
            self.membrane_conc_mass_lithium,
            wrt=self.z_bar,
            units=units.kg / units.m**3,
            doc="Lithium concentration gradient wrt z in the membrane",
        )
        self.d_membrane_conc_mass_cobalt_dz = DerivativeVar(
            self.membrane_conc_mass_cobalt,
            wrt=self.z_bar,
            units=units.kg / units.m**3,
            doc="Cobalt concentration gradient wrt z in the membrane",
        )

    def add_constraints(self):
        """
        Adds model constraints for the two salt diafiltration unit model.
        """

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
        Assigns scaling factors to certain variables and constraints to
        improve solver performance.
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
