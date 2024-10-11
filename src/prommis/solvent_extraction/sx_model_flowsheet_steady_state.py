#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
""" 
Demonstration flowsheet for steady state solvent extraction loading process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import ConcreteModel

import numpy as np

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
This is a loading operation, so no additional argument has to be specified.

"""

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
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

"""
Specification of the values of the partition coefficients of the elements
based on the values provided in the REESim file. 

"""

stage_number = np.arange(1, number_of_stages + 1)

for s in stage_number:
    if s == 1:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 5.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 24.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 99.1 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = 99.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 32.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = 58.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 58.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = 87.6 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = 99.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = 69.8 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = 96.6 / 100
    else:
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 4.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 12.3 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 6.4 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 16.7 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = 99.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 23.2 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = 24.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 15.1 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = 99.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = 99.9 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = 7.6 / 100
        m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = 5 / 100

"""
Fixing the inlet conditions of the two feed streams to the solvent extraction model,
based on a case study of University of Kentucky pilot plant.

"""
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(1.755)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(3999.818)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(693.459)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)
m.fs.solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-8)

m.fs.solex.aqueous_inlet.flow_vol.fix(62.01)

m.fs.solex.organic_inlet.conc_mass_comp[0, "Al"].fix(1.267e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Ca"].fix(2.684e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Fe"].fix(2.873e-6)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Sc"].fix(1.734)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Y"].fix(2.179e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "La"].fix(0.000105)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Ce"].fix(0.00031)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Pr"].fix(3.711e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Nd"].fix(0.000165)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Sm"].fix(1.701e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Gd"].fix(3.357e-5)
m.fs.solex.organic_inlet.conc_mass_comp[0, "Dy"].fix(8.008e-6)

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
m.fs.solex.organic_outlet.conc_mass_comp.display()

# Final aqueous outlets display
m.fs.solex.aqueous_outlet.conc_mass_comp.display()
