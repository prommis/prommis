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

from math import log10

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
    Var,
    value,
)

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.mscontactor import (
    MSContactor,
    MSContactorInitializer,
)
from idaes.core.util import DiagnosticsToolbox

from prommis_workspace.leaching.leach_solution_properties import LeachSolutionParameters
from prommis_workspace.leaching.leach_solids_properties import CoalRefuseParameters
from prommis_workspace.leaching.leach_reactions import CoalRefuseLeachingReactions


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
m.fs.leach.liquid_inlet.conc_mass_comp.fix(1e-10 * units.mg / units.L)

m.fs.leach.liquid_inlet.conc_mass_comp[0, "H"].fix(2 * 0.05 * 1e3 * units.mg / units.L)
m.fs.leach.liquid_inlet.conc_mass_comp[0, "HSO4"].fix(1e-8 * units.mg / units.L)
m.fs.leach.liquid_inlet.conc_mass_comp[0, "SO4"].fix(0.05 * 96e3 * units.mg / units.L)

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
    return (
        b.heterogeneous_reaction_extent[t, s, r]
        == b.heterogeneous_reactions[t, s].reaction_rate[r]*b.volume[t,s]
    )

m.fs.leach.heterogeneous_reaction_extent_constraint = Constraint(
    m.fs.time,
    m.fs.leach.elements,
    m.fs.leach_rxns.reaction_idx,
    rule=rule_heterogeneous_reaction_extent,
)

dt = DiagnosticsToolbox(m)
dt.assert_no_structural_warnings()

# -------------------------------------------------------------------------------------
# Scaling
m.scaling_factor = Suffix(direction=Suffix.EXPORT)

for j in m.fs.coal.component_list:
    if j not in ["Al2O3", "Fe2O3", "CaO", "inerts"]:
        m.scaling_factor[m.fs.leach.solid[0.0, 1].mass_frac_comp[j]] = 1e5
        m.scaling_factor[m.fs.leach.solid_inlet_state[0.0].mass_frac_comp[j]] = 1e5
        m.scaling_factor[m.fs.leach.heterogeneous_reactions[0.0, 1].reaction_rate[j]] = 1e5
        m.scaling_factor[m.fs.leach.solid[0.0, 1].conversion_eq[j]] = 1e3
        m.scaling_factor[m.fs.leach.solid_inlet_state[0.0].conversion_eq[j]] = 1e3
        m.scaling_factor[m.fs.leach.heterogeneous_reactions[0.0, 1].reaction_rate_eq[j]] = 1e5

scaling = TransformationFactory('core.scale_model')
scaled_model = scaling.create_using(m, rename=False)

# -------------------------------------------------------------------------------------

initializer = MSContactorInitializer()
try:
    initializer.initialize(scaled_model.fs.leach)
except:
    pass

solver = SolverFactory("ipopt")
solver.solve(scaled_model, tee=True)

scaling.propagate_solution(scaled_model, m)

m.fs.leach.liquid_outlet.display()
m.fs.leach.solid_outlet.display()

for j in m.fs.coal.component_list:
    f_in = m.fs.leach.solid_inlet.flow_mass[0]
    f_out = m.fs.leach.solid_outlet.flow_mass[0]
    x_in = m.fs.leach.solid_inlet.mass_frac_comp[0, j]
    x_out = m.fs.leach.solid_outlet.mass_frac_comp[0, j]

    r = value(1 - f_out*x_out/(f_in*x_in))*100

    print(f"Recovery {j}: {r}")

print(f"pH in {-log10(value(m.fs.leach.liquid_inlet_state[0].conc_mol_comp['H']))}")
print(f"pH out {-log10(value(m.fs.leach.liquid[0, 1].conc_mol_comp['H']))}")

m.fs.leach.solid[0, 1].conversion.display()
m.fs.leach.liquid[0, 1].dens_mol.display()
