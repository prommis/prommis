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

This multi-component model for the diafiltration membrane model a two salt system with a common
anion. The membrane is designed for use in the diafiltration cascade, i.e., is a spiral-wound
membrane module.

Configuration Arguments
-----------------------

Model Structure
---------------

Additional Constraints
----------------------

Assumptions
-----------

Variables
---------
"""

from pyomo.common.config import ConfigBlock, ConfigValue

from pyomo.dae import (
    ContinuousSet,
    DerivativeVar,
)
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    TransformationFactory,
    units,
    Var,
)
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
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
        "membrane_length",
        ConfigValue(
            doc="Length of the membrane, wound radially (m)",
        ),
    )
    CONFIG.declare(
        "membrane_width",
        ConfigValue(
            doc="Width of the membrane, parallel to the surface (m)",
        ),
    )
    CONFIG.declare(
        "membrane_thickness",
        ConfigValue(
            doc="Thickness of membrane (m)",
        ),
    )
    CONFIG.declare(
        "membrane_permeability",
        ConfigValue(
            doc="Hydraulic permeability coefficient (m/h/bar)",
        ),
    )
    CONFIG.declare(
        "applied_pressure",
        ConfigValue(
            doc="Pressure applied to membrane (bar)",
        ),
    )
    CONFIG.declare(
        "feed_flow_volume",
        ConfigValue(
            doc="Volumetric flow rate of the feed (m3/h)",
        ),
    )
    CONFIG.declare(
        "feed_conc_mass_lithium",
        ConfigValue(
            doc="Mass concentration of lithium in the feed (kg/m3)",
        ),
    )
    CONFIG.declare(
        "feed_conc_mass_cobalt",
        ConfigValue(
            doc="Mass concentration of cobalt in the feed (kg/m3)",
        ),
    )
    CONFIG.declare(
        "diafiltrate_flow_volume",
        ConfigValue(
            doc="Volumetric flow rate of the diafiltrate (m3/h)",
        ),
    )
    CONFIG.declare(
        "diafiltrate_conc_mass_lithium",
        ConfigValue(
            doc="Mass concentration of lithium in the diafiltrate (kg/m3)",
        ),
    )
    CONFIG.declare(
        "diafiltrate_conc_mass_cobalt",
        ConfigValue(
            doc="Mass concentration of cobalt in the diafiltrate (kg/m3)",
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

        self.add_variables()
        self.add_constraints()
        self.discretize_model()

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
            initialize=0,
            units=units.m**3 / units.h,
            domain=NonNegativeReals,
            doc="Volumetric flow rate of the permeate, x-dependent",
        )
        self.permeate_conc_mass_lithium = Var(
            self.x_bar,
            initialize=0,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of lithium in the permeate, x-dependent",
        )
        self.permeate_conc_mass_cobalt = Var(
            self.x_bar,
            initialize=0,
            units=units.kg / units.m**3,
            domain=NonNegativeReals,
            doc="Mass concentration of cobalt in the permeate, x-dependent",
        )
        self.permeate_conc_mass_chlorine = Var(
            self.x_bar,
            initialize=0,
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
                -self.volume_flux_water[x]
                * self.config.membrane_length
                * self.config.membrane_width
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
                * self.config.membrane_length
                * self.config.membrane_width
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
                * self.config.membrane_length
                * self.config.membrane_width
            )

        self.cobalt_mass_balance = Constraint(self.x_bar, rule=_cobalt_mass_balance)

        def _general_mass_balance_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_conc_mass_lithium[x] * self.retentate_flow_volume[x]
                + self.permeate_conc_mass_lithium[x] * self.permeate_flow_volume[x]
            ) == (
                self.config.feed_flow_volume * self.config.feed_conc_mass_lithium
                + self.config.diafiltrate_flow_volume
                * self.config.diafiltrate_conc_mass_lithium
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
                self.config.feed_flow_volume * self.config.feed_conc_mass_cobalt
                + self.config.diafiltrate_flow_volume
                * self.config.diafiltrate_conc_mass_cobalt
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
                * self.config.membrane_length
                * self.config.membrane_width
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
                self.config.membrane_permeability
                * (self.config.applied_pressure - self.osmotic_pressure[x])
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
                    / self.config.membrane_thickness
                    * self.d_membrane_conc_mass_lithium_dz[x, z]
                )
                + (
                    self.D_lithium_cobalt[x, z]
                    / self.config.membrane_thickness
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
                    / self.config.membrane_thickness
                    * self.d_membrane_conc_mass_lithium_dz[x, z]
                )
                + (
                    self.D_cobalt_cobalt[x, z]
                    / self.config.membrane_thickness
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
                        * 298
                        * units.K  # assume room temp
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
                self.config.feed_flow_volume + self.config.diafiltrate_flow_volume
            )

        self.initial_retentate_flow_volume = Constraint(
            rule=_initial_retentate_flow_volume
        )

        def _initial_permeate_flow_volume(self):
            return self.permeate_flow_volume[0] == (0 * units.m**3 / units.h)

        self.initial_permeate_flow_volume = Constraint(
            rule=_initial_permeate_flow_volume
        )

        def _initial_retentate_conc_mass_lithium(self):
            return self.retentate_conc_mass_lithium[0] == (
                (
                    self.config.feed_flow_volume * self.config.feed_conc_mass_lithium
                    + self.config.diafiltrate_flow_volume
                    * self.config.diafiltrate_conc_mass_lithium
                )
                / (self.config.feed_flow_volume + self.config.diafiltrate_flow_volume)
            )

        self.initial_retentate_conc_mass_lithium = Constraint(
            rule=_initial_retentate_conc_mass_lithium
        )

        def _initial_retentate_conc_mass_cobalt(self):
            return self.retentate_conc_mass_cobalt[0] == (
                (
                    self.config.feed_flow_volume * self.config.feed_conc_mass_cobalt
                    + self.config.diafiltrate_flow_volume
                    * self.config.diafiltrate_conc_mass_cobalt
                )
                / (self.config.feed_flow_volume + self.config.diafiltrate_flow_volume)
            )

        self.initial_retentate_conc_mass_cobalt = Constraint(
            rule=_initial_retentate_conc_mass_cobalt
        )

        def _initial_membrane_interface_lithium(self):
            return self.membrane_conc_mass_lithium[0, 0] == (0 * units.kg / units.m**3)

        self.initial_membrane_interface_lithium = Constraint(
            rule=_initial_membrane_interface_lithium
        )

        def _initial_membrane_interface_cobalt(self):
            return self.membrane_conc_mass_cobalt[0, 0] == (0 * units.kg / units.m**3)

        self.initial_membrane_interface_cobalt = Constraint(
            rule=_initial_membrane_interface_cobalt
        )

        def _initial_membrane_interface_chlorine(self):
            return self.membrane_conc_mass_chlorine[0, 0] == (0 * units.kg / units.m**3)

        self.initial_membrane_interface_chlorine = Constraint(
            rule=_initial_membrane_interface_chlorine
        )

        def _initial_permeate_conc_mass_lithium(self):
            return self.permeate_conc_mass_lithium[0] == (0 * units.kg / units.m**3)

        self.initial_permeate_conc_mass_lithium = Constraint(
            rule=_initial_permeate_conc_mass_lithium
        )

        def _initial_permeate_conc_mass_cobalt(self):
            return self.permeate_conc_mass_cobalt[0] == (0 * units.kg / units.m**3)

        self.initial_permeate_conc_mass_cobalt = Constraint(
            rule=_initial_permeate_conc_mass_cobalt
        )

        def _initial_d_retentate_conc_mass_lithium_dx(self):
            return self.d_retentate_conc_mass_lithium_dx[0] == (
                0 * units.kg / units.m**3
            )

        self.initial_d_retentate_conc_mass_lithium_dx = Constraint(
            rule=_initial_d_retentate_conc_mass_lithium_dx
        )

        def _initial_d_retentate_conc_mass_cobalt_dx(self):
            return self.d_retentate_conc_mass_cobalt_dx[0] == (
                0 * units.kg / units.m**3
            )

        self.initial_d_retentate_conc_mass_cobalt_dx = Constraint(
            rule=_initial_d_retentate_conc_mass_cobalt_dx
        )

        def _initial_d_retentate_flow_volume_dx(self):
            return self.d_retentate_flow_volume_dx[0] == (0 * units.m**3 / units.h)

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

    # TODO: add ports
