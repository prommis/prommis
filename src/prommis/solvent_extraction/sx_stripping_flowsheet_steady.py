#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
""" 
Demonstration flowsheet for steady state solvent extraction stripping process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import ConcreteModel

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
)

"""
Method of building a solvent extraction model with a specified number of stages
and with two separate property packages for the two inlet streams.
This is a stripping operation, so an additional argument regarding the direction 
of mass transfer flow has to be specified.
"""

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()

number_of_stages = 1

m.fs.solex = SolventExtraction(
    number_of_finite_elements=number_of_stages,
    dynamic=False,
    aqueous_stream={
        "property_package": m.fs.leach_soln,
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

"""
Specification of the values of the partition coefficients of the elements
based on the values provided in the REESim file. 

"""

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

"""
Fixing the inlet conditions of the two feed streams to the solvent extraction model,
based on a case study of University of Kentucky pilot plant.

"""

m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(1e-9)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-9)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-9)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(128.94)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(878.9)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(44.723)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(288.01)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(108.94)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(287.49)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(731.14)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(127.25)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(352.60)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(64.22)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(30.39)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(17.63)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-8)

m.fs.solex.aqueous_inlet.flow_vol.fix(62.01)

m.fs.solex.organic_inlet.conc_mass_comp[0, "Al"].fix(52.249)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Ca"].fix(356.139)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Fe"].fix(18.122)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Sc"].fix(39.3)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Y"].fix(44.14)
m.fs.solex.organic_inlet.conc_mass_comp[0, "La"].fix(116.49)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Ce"].fix(296.26)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Pr"].fix(51.56)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Nd"].fix(142.88)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Sm"].fix(26.02)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Gd"].fix(12.31)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Dy"].fix(7.14)

m.fs.solex.organic_inlet.flow_vol.fix(62.01)

"""
Initialization of the model, which gives a good starting point.

"""

initializer = SolventExtractionInitializer()
initializer.initialize(m.fs.solex)

"""
Solution of the model and display of the final results.

"""

solver = get_solver("ipopt")
solver.solve(m, tee=True)

# Final organic outlet display
m.fs.solex.mscontactor.organic_outlet.conc_mass_comp.display()

# Final aqueous outlets display
m.fs.solex.mscontactor.aqueous_outlet.conc_mass_comp.display()
