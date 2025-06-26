from pyomo.environ import ConcreteModel, units, value, TransformationFactory

import numpy as np

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json
from idaes.core.util.scaling import set_scaling_factor

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom as dof

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.compound_solvent_extraction import (
    CompoundSolventExtraction,
    CompoundSolventExtractionInitializer,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.reaxn = SolventExtractionReactions()

dosage = 5
m.fs.reaxn.extractant_dosage = dosage

m.fs.compound_solex = CompoundSolventExtraction(
    number_of_finite_elements=3,
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
    heterogeneous_reaction_package=m.fs.reaxn,
    has_holdup=True,
    settler_transformation_method="dae.finite_difference",
    settler_transformation_scheme="BACKWARD",
    settler_finite_elements=4,
)

m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(10.75)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(100)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e4)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)

m.fs.compound_solex.aqueous_inlet.flow_vol.fix(62.01)

m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(975.8e3 * dosage / 100)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1.267e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(2.684e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(2.873e-6)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1.734)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(2.179e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "La_o"].fix(0.000105)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.00031)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(3.711e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(0.000165)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1.701e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(3.357e-5)
m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(8.008e-6)

m.fs.compound_solex.organic_inlet.flow_vol.fix(62.01)

m.fs.compound_solex.solvent_extraction_mixer[:].unit.mscontactor.aqueous[
    :, :
].temperature.fix(305.15 * units.K)
m.fs.compound_solex.solvent_extraction_mixer[:].unit.mscontactor.organic[
    :, :
].temperature.fix(305.15 * units.K)

m.fs.compound_solex.solvent_extraction_mixer[:].unit.mscontactor.volume[:].fix(
    0.4 * units.m**3
)
m.fs.compound_solex.solvent_extraction_mixer[:].unit.mscontactor.volume_frac_stream[
    :, :, "organic"
].fix(0.5)

m.fs.compound_solex.solvent_extraction_organic_settler[:].unit.area.fix(1)
m.fs.compound_solex.solvent_extraction_aqueous_settler[:].unit.area.fix(1)
m.fs.compound_solex.solvent_extraction_aqueous_settler[:].unit.length.fix(1)
m.fs.compound_solex.solvent_extraction_organic_settler[:].unit.length.fix(1)

for s in m.fs.compound_solex.elements:
    for x in m.fs.compound_solex.solvent_extraction_aqueous_settler[
        s
    ].unit.length_domain:
        if x != 0:
            m.fs.compound_solex.solvent_extraction_aqueous_settler[s].unit.properties[
                :, x
            ].temperature.fix(305.15 * units.K)
            P = value(
                m.fs.compound_solex.solvent_extraction_aqueous_settler[s]
                .unit.properties[0, 0]
                .pressure
            )
            m.fs.compound_solex.solvent_extraction_aqueous_settler[s].unit.properties[
                :, x
            ].pressure.fix(P * units.Pa)
    for x in m.fs.compound_solex.solvent_extraction_organic_settler[
        s
    ].unit.length_domain:
        if x != 0:
            m.fs.compound_solex.solvent_extraction_organic_settler[s].unit.properties[
                :, x
            ].temperature.fix(305.15 * units.K)
            P = value(
                m.fs.compound_solex.solvent_extraction_organic_settler[s]
                .unit.properties[0, 0]
                .pressure
            )
            m.fs.compound_solex.solvent_extraction_organic_settler[s].unit.properties[
                :, x
            ].pressure.fix(P * units.Pa)

# for s in m.fs.compound_solex.elements:
#     set_scaling_factor(
#         m.fs.compound_solex.solvent_extraction_mixer[s]
#         .unit.mscontactor.aqueous[0.0, 1]
#         .hso4_dissociation,
#         10,
#     )

print("Degrees of freedom:", dof(m.fs.compound_solex))
m.fs.compound_solex.fix_initialization_states()
print("DoF after fix:", dof(m.fs.compound_solex))

initializer = CompoundSolventExtractionInitializer()
initializer.initialize(m.fs.compound_solex)

# scaling = TransformationFactory("core.scale_model")
# scaled_model = scaling.create_using(m, rename=False)

# print(dof(scaled_model))

solver = get_solver("ipopt_v2")
solver.options["max_iter"] = 10000
# # solver.options["tol"] = 1e-5
results = solver.solve(m, tee=True)

# scaling.propagate_solution(scaled_model, m)

# to_json(m, fname="compound_solvent_extraction.json", human_read=True)
