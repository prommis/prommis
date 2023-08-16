from pyomo.environ import (
    ConcreteModel,
    Constraint,
    SolverFactory,
    units,
    Var,
    value,
)

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom as dof

from REESXmodel import REESX

from REEdistribution import REESolExParameters

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop = REESolExParameters()

m.fs.solex = REESX(number_of_finite_elements=1,
                       aqstreams = {"Acidsoln":{"property_package":m.fs.prop}},
                       ogstreams = {"Orgacid":{"property_package":m.fs.prop}})


#Aqueous feed fixing
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Al"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ca"].fix(0.02)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Fe"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Si"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Sc"].fix(0.92)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Y"].fix(2.82)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["La"].fix(8.98)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ce"].fix(19.94)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Pr"].fix(3.34)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Nd"].fix(9.04)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Pm"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Sm"].fix(1.63)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Eu"].fix(0.11)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Gd"].fix(0.77)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Tb"].fix(0.33)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Dy"].fix(0.45)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ho"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Er"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Tm"].fix(0.18)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Yb"].fix(0.29)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Lu"].fix(0.14)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Th"].fix(0)
m.fs.solex.Acidsoln_inlet_state[0].mass_flow["U"].fix(0)

m.fs.solex.Acidsoln_inlet_state[0].volumetric_flow.fix(4.4)

#Organic feed fixing
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Al"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ca"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Fe"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Si"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Sc"].fix(19.93)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Y"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["La"].fix(2.09)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ce"].fix(0.86)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Pr"].fix(0.12)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Nd"].fix(0.07)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Pm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Sm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Eu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Gd"].fix(0.01)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Tb"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Dy"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ho"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Er"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Tm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Yb"].fix(0.05)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Lu"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["Th"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].mass_flow["U"].fix(0)

m.fs.solex.Orgacid_inlet_state[0].volumetric_flow.fix(62.01)

print(dof(m))

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

m.fs.solex.Orgacid[0,1].mass_flow.display()
m.fs.solex.Acidsoln[0,1].mass_flow.display()


