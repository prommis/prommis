from pyomo.environ import (
    ConcreteModel,
    SolverFactory, 
)

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom as dof

from workspace.UKy_flowsheet.Solvent_Extraction.REESXmodel import REESX
from workspace.UKy_flowsheet.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from workspace.UKy_flowsheet.Solvent_Extraction.REEOgdistribution import REESolExOgParameters

from pyomo.util.check_units import assert_units_consistent

from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)

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
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sc"].fix(0.92)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Y"].fix(2.82)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["La"].fix(8.98)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ce"].fix(19.94)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Pr"].fix(3.34)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Nd"].fix(9.04)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sm"].fix(1.63)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Gd"].fix(0.77)
m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Dy"].fix(0.45)

m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

#Organic feed fixing
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Al"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ca"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Fe"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sc"].fix(19.93)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Y"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["La"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ce"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pr"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Nd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sm"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Gd"].fix(0)
m.fs.solex.Orgacid_inlet_state[0].flow_mass["Dy"].fix(0)

m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01)

print(dof(m))

assert_units_consistent(m)

initializer = BlockTriangularizationInitializer()
initializer.initialize(m.fs.solex)
assert initializer.summary[m.fs.solex]["status"] == InitializationStatus.Ok

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
