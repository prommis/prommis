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
Initial flowsheet for UKy leaching process

Authors: Andrew Lee
"""


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    units,
    Var,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.mscontactor import (
    MSContactor,
    MSContactorInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from leach_solution_properties import LeachSolutionParameters
from leach_solids_properties import CoalRefuseParameters
from leach_reactions import CoalRefuseLeachingReactions


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.leach_soln = LeachSolutionParameters()
m.fs.coal = CoalRefuseParameters()
m.fs.leach_rxns = CoalRefuseLeachingReactions()

m.fs.leach = MSContactor(
    number_of_finite_elements=1,
    streams={
        "liquid": {
            "property_package": m.fs.leach_soln,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        "solid": {
            "property_package": m.fs.coal,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
    },
    heterogeneous_reactions=m.fs.leach_rxns,
)

# Liquid feed state
m.fs.leach.liquid_inlet.flow_vol.fix(224.3 * units.L / units.hour)
m.fs.leach.liquid_inlet.conc_mass_metals.fix(1e-10 * units.mg / units.L)
m.fs.leach.liquid_inlet.conc_mole_acid[0, "H"].fix(2 * 0.05 * units.mol / units.L)
m.fs.leach.liquid_inlet.conc_mole_acid[0, "HSO4"].fix(1e-8 * units.mol / units.L)
m.fs.leach.liquid_inlet.conc_mole_acid[0, "SO4"].fix(0.05 * units.mol / units.L)

# Solid feed state
m.fs.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
m.fs.leach.solid_inlet.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Sc2O3"].fix(2.77966E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Y2O3"].fix(3.28653E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "La2O3"].fix(6.77769E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Ce2O3"].fix(0.000156161 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Pr2O3"].fix(1.71438E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Nd2O3"].fix(6.76618E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Sm2O3"].fix(1.47926E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Gd2O3"].fix(1.0405E-05 * units.kg / units.kg)
m.fs.leach.solid_inlet.mass_frac_comp[0, "Dy2O3"].fix(7.54827E-06 * units.kg / units.kg)

# Reactor volume
m.fs.leach.volume = Var(
    m.fs.time,
    m.fs.leach.elements,
    initialize=1,
    units=units.litre,
    doc="Volume of each finite element."
)
m.fs.leach.volume.fix(100 * units.gallon)

def rule_heterogeneous_reaction_extent(b, t, s, r):
    return b.heterogeneous_reaction_extent[t, s, r] == b.heterogeneous_reactions[t, s].reaction_rate[r]*b.volume[t,s]

m.fs.leach.heterogeneous_reaction_extent_constraint = Constraint(
    m.fs.time,
    m.fs.leach.elements,
    m.fs.leach_rxns.reaction_idx,
    rule=rule_heterogeneous_reaction_extent,
)

print(degrees_of_freedom(m))

# TODO: Energy transfer term is not needed but still gets built
# m.fs.leach.del_component(m.fs.leach.energy_transfer_term)

initializer = MSContactorInitializer()
initializer.initialize(m.fs.leach)

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

m.fs.leach.liquid_outlet.display()
m.fs.leach.solid_outlet.display()

m.fs.leach.heterogeneous_reactions[0, 1].c_max.display()

for j in m.fs.coal.component_list:
    f_in = m.fs.leach.solid_inlet.flow_mass[0]
    f_out = m.fs.leach.solid_outlet.flow_mass[0]
    x_in = m.fs.leach.solid_inlet.mass_frac_comp[0, j]
    x_out = m.fs.leach.solid_outlet.mass_frac_comp[0, j]

    r = value(1 - f_out*x_out/(f_in*x_in))*100

    print(f"Recovery {j}: {r}")

from math import log10
print(f"pH in {-log10(value(m.fs.leach.liquid_inlet.conc_mole_acid[0, 'H']))}")
print(f"pH out {-log10(value(m.fs.leach.liquid_outlet.conc_mole_acid[0, 'H']))}")

# # Conservation
# h_in = (
#     m.fs.leach.liquid_inlet.flow_vol[0]
#     * (m.fs.leach.liquid_inlet.conc_mole_acid[0, 'HSO4']
#        + m.fs.leach.liquid_inlet.conc_mole_acid[0, 'H']
#        + 2/18e-3))
# h_out = (
#     m.fs.leach.liquid_outlet.flow_vol[0]
#     * (m.fs.leach.liquid_outlet.conc_mole_acid[0, 'HSO4']
#        + m.fs.leach.liquid_outlet.conc_mole_acid[0, 'H']
#        + 2/18e-3))
# print(f"H: {value(h_in-h_out)}")
# o_in_l = (
#     m.fs.leach.liquid_inlet.flow_vol[0]
#     * (4*m.fs.leach.liquid_inlet.conc_mole_acid[0, 'HSO4']
#        + 4*m.fs.leach.liquid_inlet.conc_mole_acid[0, 'SO4']
#        + 1/18e-3))
# o_out_l = (
#     m.fs.leach.liquid_outlet.flow_vol[0]
#     * (4*m.fs.leach.liquid_outlet.conc_mole_acid[0, 'HSO4']
#        + 4*m.fs.leach.liquid_outlet.conc_mole_acid[0, 'SO4']
#        + 1/18e-3))
# biatomic= ["Sc2O3", "Y2O3", "La2O3", "Ce2O3", "Pr2O3", "Nd2O3", "Sm2O3", "Gd2O3", "Dy2O3", "Al2O3", "Fe2O3"]
# o_in_s = (
#     m.fs.leach.solid_inlet.flow_mass[0]
#     * (3*sum(m.fs.leach.solid_inlet.mass_frac_comp[0, j]/m.fs.coal.mw[j] for j in biatomic)
#        + m.fs.leach.solid_inlet.mass_frac_comp[0, "CaO"]/m.fs.coal.mw["CaO"])
# )
# o_out_s = (
#     m.fs.leach.solid_outlet.flow_mass[0]
#     * (3*sum(m.fs.leach.solid_outlet.mass_frac_comp[0, j]/m.fs.coal.mw[j] for j in biatomic)
#        + m.fs.leach.solid_outlet.mass_frac_comp[0, "CaO"]/m.fs.coal.mw["CaO"])
# )
# print(f"O: {value(o_in_l-o_out_l+o_in_s-o_out_s)} ({value(o_in_l-o_out_l)}, {value(o_in_s-o_out_s)})")
# print("S " + str(value(
#         m.fs.leach.liquid_inlet.flow_vol[0]
#         * (m.fs.leach.liquid_inlet.conc_mole_acid[0, 'HSO4']+m.fs.leach.liquid_inlet.conc_mole_acid[0, 'SO4'])
#         - m.fs.leach.liquid_outlet.flow_vol[0]
#         * (m.fs.leach.liquid_outlet.conc_mole_acid[0, 'HSO4']+m.fs.leach.liquid_outlet.conc_mole_acid[0, 'SO4'])
#     )))
#
# for j in m.fs.leach_soln.dissolved_metals_set:
#     if j == "Ca":
#         k = "CaO"
#         n = 1
#     else:
#         k = f"{j}2O3"
#         n = 2
#
#     l_in = m.fs.leach.liquid_inlet.flow_vol[0] * m.fs.leach.liquid_inlet.conc_mass_metals[0, j]
#     l_out = m.fs.leach.liquid_outlet.flow_vol[0] * m.fs.leach.liquid_outlet.conc_mass_metals[0, j]
#     s_in = m.fs.leach.solid_inlet.flow_mass[0]*m.fs.leach.solid_inlet.mass_frac_comp[0, k]
#     s_out = m.fs.leach.solid_outlet.flow_mass[0] * m.fs.leach.solid_outlet.mass_frac_comp[0, k]
#
#     l_side = value((l_in-l_out)/m.fs.leach_soln.mw[j]*1e-3)
#     s_side = value((s_in - s_out) / m.fs.coal.mw[k] * n*1e3)
#
#     print(f"{j}: {l_side+s_side} ({l_side} {s_side})")

# m.fs.leach.heterogeneous_reaction_extent.display()
m.fs.leach.heterogeneous_reactions[0, 1].reaction_rate.display()