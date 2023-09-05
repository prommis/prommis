from pyomo.environ import (
    ConcreteModel,
    SolverFactory, 
)

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom as dof

from REESXmodel import REESX

from REEAqdistribution import REESolExAqParameters
from REEOgdistribution import REESolExOgParameters

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_a = REESolExAqParameters()
m.fs.prop_o = REESolExOgParameters()

m.fs.solex = REESX(number_of_finite_elements=3,
                       aqueous_streams = {"Acidsoln":{"property_package":m.fs.prop_a, "flow_direction":1}},
                       organic_streams = {"Orgacid":{"property_package":m.fs.prop_o, "flow_direction":2}})


#Aqueous feed fixing
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Al"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ca"].fix(0.02)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Fe"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Si"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sc"].fix(0.92)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Y"].fix(2.82)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["La"].fix(8.98)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ce"].fix(19.94)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Pr"].fix(3.34)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Nd"].fix(9.04)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Pm"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sm"].fix(1.63)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Eu"].fix(0.11)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Gd"].fix(0.77)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Tb"].fix(0.33)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Dy"].fix(0.45)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ho"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Er"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Tm"].fix(0.18)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Yb"].fix(0.29)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Lu"].fix(0.14)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Th"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["U"].fix(0)

m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

#Organic feed fixing
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Al"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ca"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Fe"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Si"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sc"].fix(19.93)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Y"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["La"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ce"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pr"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Nd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Eu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Gd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Tb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Dy"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ho"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Er"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Tm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Yb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Lu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Th"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["U"].fix(0)

m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01)

print(dof(m))

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

# All organic outlets display
m.fs.solex.Orgacid[0,1].flow_mass.display()
m.fs.solex.Orgacid[0,2].flow_mass.display()
m.fs.solex.Orgacid[0,3].flow_mass.display()

# All aqueous outlets display
m.fs.solex.Acidsoln[0,1].flow_mass.display()
m.fs.solex.Acidsoln[0,2].flow_mass.display()
m.fs.solex.Acidsoln[0,3].flow_mass.display()
