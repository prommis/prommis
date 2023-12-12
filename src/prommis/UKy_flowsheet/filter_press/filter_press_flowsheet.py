#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Initial flowsheet for filter press

Authors: Marcus Holly
"""


from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TransformationFactory,
    units,
    Var,
)
from pyomo.network import Arc

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.separator import (
    Separator,
    SplittingType,
    SeparatorInitializer,
)
from idaes.models.unit_models import Product
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed

from prommis.UKy_flowsheet.filter_press.filter_press_solids_properties import FilterPressSolidsParameters
from prommis.UKy_flowsheet.filter_press.filter_press_solution_properties import FilterPressParameters
from idaes.core.util.exceptions import ConfigurationError


# Build flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.liq_prop = FilterPressParameters()
m.fs.solid_prop = FilterPressSolidsParameters()

m.fs.liq_feed = Feed(property_package=m.fs.liq_prop)
m.fs.sol_feed = Feed(property_package=m.fs.solid_prop)

# Separates filtrate from residual liquid in the filter cake
m.fs.sep = Separator(
    property_package=m.fs.liq_prop,
    outlet_list=["filtrate", "residual_liquid"],
    split_basis=SplittingType.totalFlow,
)

# Liquid outlets
m.fs.filtrate = Product(property_package=m.fs.liq_prop)
m.fs.residual_liquid = Product(property_package=m.fs.liq_prop)

# Solid outlet
m.fs.filter_cake = Product(property_package=m.fs.solid_prop)

# Treat solids as a pass-through
m.fs.s01 = Arc(source=m.fs.sol_feed.outlet, destination=m.fs.filter_cake.inlet)

# Separate liquids into filtrate and residual liquid
m.fs.s02 = Arc(source=m.fs.liq_feed.outlet, destination=m.fs.sep.inlet)
m.fs.s03 = Arc(source=m.fs.sep.filtrate, destination=m.fs.filtrate.inlet)
m.fs.s04 = Arc(source=m.fs.sep.residual_liquid, destination=m.fs.residual_liquid.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

# Liquid inlet feed conditions
m.fs.liq_feed.flow_vol.fix(100 * units.L / units.hour)
m.fs.liq_feed.conc_mole[0, "TDS"].fix(0.1 * units.mol / units.L)
m.fs.liq_feed.temperature.fix(298 * units.K)
m.fs.liq_feed.pressure.fix(1e5 * units.Pa)

# Solid inlet feed conditions
m.fs.sol_feed.flow_mass.fix(25 * units.kg / units.hour)
m.fs.sol_feed.mass_frac_comp[0, "solid1"].fix(0.5 * units.kg / units.kg)
m.fs.sol_feed.mass_frac_comp[0, "solid2"].fix(0.5 * units.kg / units.kg)
m.fs.sol_feed.temperature.fix(298 * units.K)
m.fs.sol_feed.pressure.fix(1e5 * units.Pa)

# Define variables - these would go in a unit model

m.fs.sep.cake_area = Var(
    initialize=0.1,
    units=units.meter ** 2,
    doc="Cross-sectional area of filter cake",
)
# Filter cake is assumed to be incompressible, so this is constant
m.fs.sep.cake_voidage = Var(
    initialize=0.01,
    units=units.dimensionless,
    doc="Voidage of filter cake",
)
m.fs.sep.solids_density = Var(
    initialize=2000,
    units=units.kg / units.meter ** 3,
    doc="Density of the solids",
)
m.fs.sep.filtrate_density = Var(
    initialize=1000,
    units=units.kg / units.meter ** 3,
    doc="Density of filtrate",
)
m.fs.sep.filtrate_viscosity = Var(
    initialize=0.01,
    units=units.Pa * units.hr,
    doc="Density of filtrate",
)
m.fs.sep.specific_surface = Var(
    initialize=1,
    units=units.meter ** -1,
    doc="Specific surface of the particles",
)
# TODO: Not sure exactly what is meant by this term - is it constant?
m.fs.sep.solids_mass_frac = Var(
    initialize=0.2,
    units=units.dimensionless,
    doc="Mass fraction of solids in the original suspension",
)
m.fs.sep.pressure_drop = Var(
    initialize=0.1e5,
    units=units.Pa,
    doc="Constant pressure difference",
)

# Equation 7.3
specific_resistance = (
        5 * (1 - m.fs.sep.cake_voidage) ** 2
        * m.fs.sep.specific_surface ** 2
        / m.fs.sep.cake_voidage ** 3
)

# Specific resistance is constant for incompressible cakes
m.fs.sep.specific_resistance = Var(
    initialize=specific_resistance,
    units=units.meter ** -2,
    doc="Specific resistance of the filter cake",
)

# Equation 7.7
cake_volume_frac = (
        m.fs.sep.solids_mass_frac
        * m.fs.sep.filtrate_density
        / (
                (1 - m.fs.sep.solids_mass_frac)
                * (1 - m.fs.sep.cake_voidage) * m.fs.sep.solids_density
                - m.fs.sep.solids_mass_frac * m.fs.sep.cake_voidage * m.fs.sep.filtrate_density
        )
)

m.fs.sep.cake_volume_frac = Var(
    initialize=cake_volume_frac,
    units=units.dimensionless,
    doc="Volume of cake deposited by unit volume of filtrate",
)

residual_liquid_split_frac = cake_volume_frac / (1 + cake_volume_frac)

m.fs.sep.split_fraction[0, "residual_liquid"].fix(residual_liquid_split_frac)

print(f"Degrees of freedom before separator initialization are: {degrees_of_freedom(m)}")

initializer = SeparatorInitializer()
initializer.initialize(m.fs.sep)

# Separator outlet conditions
m.fs.sep.filtrate.temperature.fix(298 * units.K)
# m.fs.sep.filtrate.pressure.fix(1e5 * units.Pa)

m.fs.sep.residual_liquid.temperature.fix(298 * units.K)
m.fs.sep.residual_liquid.pressure.fix(1e5 * units.Pa)

print(f"Degrees of freedom after separator initialization are: {degrees_of_freedom(m)}")

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

m.fs.sep.filtrate.display()
m.fs.sep.residual_liquid.display()
m.fs.filter_cake.display()

#TODO: Not sure how to handle pressure drop or "volume of filtrate which has passed in time t" term

# if set_P_constant == True and set_Q_constant == True:
#     raise ConfigurationError(
#         "Pressure and flow rate cannot both be held constant. Change one of these settings to False."
#     )

# Eq. 7.9
# if set_Q_constant == True:
#     pressure_drop = (
#             m.fs.sep.filtrate.flow_vol
#             * m.fs.sep.specific_resistance
#             * m.fs.sep.filtrate_viscosity
#             * "volume of filtrate which has passed in time t"
#             * m.fs.sep.cake_volume_frac
#             / m.fs.sep.cake_area**2
#     )

# Eq. 7.12
# if set_P_constant == True:
#     m.fs.sep.pressure_drop.fix()
#     filtrate_flow = (
#         2 * m.fs.sep.cake_area**2
#         * (-m.fs.sep.pressure_drop)
#         / (
#             m.fs.sep.specific_resistance
#             * m.fs.sep.filtrate_viscosity
#             * m.fs.sep.cake_volume_frac
#             "volume of filtrate which has passed in time t"
#         )
#     )

#
