from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.util import DiagnosticsToolbox

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.ree_aq_distribution import REESolExAqParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.prop_a = REESolExAqParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 3

m.fs.solex = SolventExtraction(
    number_of_finite_elements=number_of_stages,
    dynamic=False,
    aqueous_stream={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
)

m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sc"] = 99.9 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "La"] = 75.2 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ce"] = 95.7 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Pr"] = 96.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Nd"] = 99.2 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sm"] = 99.9 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Gd"] = 98.6 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Dy"] = 99.9 / 100

m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(820)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(270)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)

m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol.fix(4.4)

m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(7.54e-10)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(4.955e-9)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(1.491e-7)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(321.34)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(5.67e-6)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(1.78e-05)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(4.019e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(6.73e-6)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(1.82e-5)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(3.285e-6)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(1.55e-6)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(9e-7)

m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

print(dof(m))

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()

# Initializing of the model

initializer = BlockTriangularizationInitializer()
initializer.initialize(m.fs.solex)
assert initializer.summary[m.fs.solex]["status"] == InitializationStatus.Ok

# Solving of the model

solver = SolverFactory("ipopt")
solver.options["bound_push"] = 1e-8
solver.options["mu_init"] = 1e-8
solver.solve(m, tee=True)

# Final organic outlet display
m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp.display()
m.fs.solex.mscontactor.organic[0, 1].conc_mol_comp.display()

# Final aqueous outlets display
m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp.display()
m.fs.solex.mscontactor.aqueous[0, 3].conc_mol_comp.display()
