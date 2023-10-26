from pyomo.environ import *

import numpy as np

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.core.initialization.block_triangularization import BlockTriangularizationInitializer
from idaes.core.initialization import InitializationStatus

from REESXmodel import REESX

from REEAqdistribution import REESolExAqParameters
from REEOgdistribution import REESolExOgParameters

#t_start = 0
#t_end = 1 * 60
#nfe = 1

#time_set = np.linspace(t_start, t_end, nfe+1)

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_a = REESolExAqParameters()
m.fs.prop_o = REESolExOgParameters()

m.fs.solex = REESX(number_of_finite_elements=1, dynamic=False,
                       aqueous_streams = {"Acidsoln":{"property_package":m.fs.prop_a, "flow_direction":1}},
                       organic_streams = {"Orgacid":{"property_package":m.fs.prop_o, "flow_direction":2}})

#TransformationFactory('dae.finite_difference').apply_to(m, nfe=nfe, wrt=m.fs.time, scheme='BACKWARD')


#Aqueous feed fixing
m.fs.solex.Acidsoln_inlet_state[0].concentration["Al"].fix(820)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Ca"].fix(5230)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Fe"].fix(270)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Si"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Sc"].fix(209.31)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Y"].fix(637.74)
m.fs.solex.Acidsoln_inlet_state[0].concentration["La"].fix(2032.77)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Ce"].fix(4516.13)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Pr"].fix(756.64)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Nd"].fix(2047.85)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Pm"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Sm"].fix(369.1)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Eu"].fix(25.81)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Gd"].fix(174.38)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Tb"].fix(75.28)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Dy"].fix(101.12)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Ho"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Er"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Tm"].fix(41.60)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Yb"].fix(65.65)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Lu"].fix(31.71)
m.fs.solex.Acidsoln_inlet_state[0].concentration["Th"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].concentration["U"].fix(0.01)

m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

m.fs.solex.Acidsoln[0,1].aqueous_vol.fix(1000)

#Organic feed fixing
m.fs.solex.Orgacid_inlet_state[0].concentration["Al"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Ca"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Fe"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Si"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Sc"].fix(321.34)
m.fs.solex.Orgacid_inlet_state[0].concentration["Y"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["La"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Ce"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Pr"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Nd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Pm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Sm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Eu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Gd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Tb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Dy"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Ho"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Er"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Tm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Yb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Lu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["Th"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].concentration["U"].fix(0)

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
m.fs.solex.Orgacid[0,1].concentration.display()

# Final aqueous outlets display
m.fs.solex.Acidsoln[0,1].concentration.display()
