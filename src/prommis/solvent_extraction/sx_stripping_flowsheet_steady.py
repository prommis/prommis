from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom as dof

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.ree_aq_distribution import REESolExAqParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.prop_a = REESolExAqParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 1

m.fs.solex = SolventExtraction(
    number_of_finite_elements=number_of_stages,
    dynamic=False,
    aqueous_stream={
        "property_package": m.fs.prop_a,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    aqueous_to_organic=False,
)

m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Al"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ca"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Fe"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sc"] = 1 - 98.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Y"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "La"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ce"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Pr"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Nd"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sm"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Gd"] = 1 - 0.5 / 100
m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Dy"] = 1 - 0.5 / 100

m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e-9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(128.94)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(878.9)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(44.723)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(288.01)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(108.94)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(287.49)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(731.14)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(127.25)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(352.60)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(64.22)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(30.39)
m.fs.solex.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(17.63)

m.fs.solex.mscontactor.aqueous_inlet_state[0].flow_vol.fix(0.13)

m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(52.249)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(356.139)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(18.122)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(39.3)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(44.14)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(116.49)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(296.26)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(51.56)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(142.88)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(26.02)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(12.31)
m.fs.solex.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(7.14)

m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

print(dof(m))

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

# Final aqueous outlets display
m.fs.solex.mscontactor.aqueous[0, 1].conc_mass_comp.display()
