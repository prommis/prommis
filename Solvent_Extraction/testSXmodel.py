import numpy as np

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.initialization.block_triangularization import BlockTriangularizationInitializer
from idaes.core.initialization import InitializationStatus

from REESXmodel import REESX

from REEAqdistribution import REESolExAqParameters
from REEOgdistribution import REESolExOgParameters

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_a = REESolExAqParameters()
m.fs.prop_o = REESolExOgParameters()

m.fs.solex = REESX(number_of_finite_elements=1, dynamic=False,
                       aqueous_streams = {"Acidsoln":{"property_package":m.fs.prop_a, "flow_direction":1}},
                       organic_streams = {"Orgacid":{"property_package":m.fs.prop_o, "flow_direction":2}})


m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Al"].fix(820)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Fe"].fix(270)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Si"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Pm"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Eu"].fix(25.81)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Tb"].fix(75.28)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ho"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Er"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Tm"].fix(41.60)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Yb"].fix(65.65)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Lu"].fix(31.71)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Th"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["U"].fix(0.01)

m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

m.fs.solex.Acidsoln[0,1].aqueous_vol.fix(1000)

m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Al"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ca"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Fe"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Si"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sc"].fix(321.34)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Y"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["La"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ce"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Pr"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Nd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Pm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Eu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Gd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Tb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Dy"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ho"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Er"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Tm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Yb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Lu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Th"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["U"].fix(0)

m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01)

m.fs.solex.Orgacid[0,1].organic_vol.fix(1000)

print(dof(m))

# BTInitializer not working for steady state 

#initializer = BlockTriangularizationInitializer()
#initializer.initialize(m.fs.solex)
#assert initializer.summary[m.fs.solex]["status"] == InitializationStatus.Ok

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

# Final organic outlet display
m.fs.solex.Orgacid[0,1].conc_mass_comp.display()

# Final aqueous outlets display
m.fs.solex.Acidsoln[0,1].conc_mass_comp.display()
