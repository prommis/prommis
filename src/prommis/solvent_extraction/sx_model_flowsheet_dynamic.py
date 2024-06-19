""" 
Demonstration flowsheet for dynamic state solvent extraction loading process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import ConcreteModel, SolverFactory, units, TransformationFactory

import numpy as np

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
)


def build_model():
    """
    Method of building a dynamic solvent extraction model with a specified number of
    stages and with two separate property packages for the two inlet streams.
    This is a loading operation, so no additional argument has to be specified.

    """

    m = ConcreteModel()

    time_duration = 24

    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )

    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()

    number_of_stages = 3

    m.fs.solex = SolventExtraction(
        number_of_finite_elements=number_of_stages,
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

    return m


def set_inputs(m):
    """
    Discretization of the time domain, and specification of the partition coefficients,
    volume, volume fractions, and the initial conditions of state variables for the components
    for all the stages.

    """

    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=5, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")

    """
    Specifications of the partition coefficients, volume and volume fractions for all
    the stages.

    """

    m.fs.solex.mscontactor.volume[:].fix(0.4)

    m.fs.solex.mscontactor.volume_frac_stream[:, :, "organic"].fix(0.4)

    number_of_stages = 3
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
    Fixation of the inlet conditions and the initial state values for all the components.

    """

    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(1.755)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(3999.818)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(693.459)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(422.375)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(109.542)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(688.266)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(0.032)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(0.124)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(0.986)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(2.277)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(0.303)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(0.946)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(0.097)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(0.2584)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(1e-7)

    m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0, :, "Ka2"].fix(0)

    m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(1.267e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(2.684e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(2.873e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(1.734)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(2.179e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(0.000105)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(0.00031)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(3.711e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(0.000165)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(1.701e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(3.357e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(8.008e-6)

    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Y"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ce"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Nd"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sm"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Gd"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Dy"].fix(1e-7)

    m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)


if __name__ == "__main__":

    time_duration = 24
    number_of_stages = 3

    # Call build model function

    m = build_model()
    set_inputs(m)

    """
    Initialization of the model, which gives a good starting point.

    """
    initializer = BlockTriangularizationInitializer(constraint_tolerance=1e-4)
    # initializer = SolventExtractionInitializer()
    initializer.initialize(m.fs.solex)

    """
    Solution of the model and display of the final results.

    """

    solver = SolverFactory("ipopt")
    solver.options["bound_push"] = 1e-8
    solver.options["mu_init"] = 1e-8
    solver.solve(m, tee=True)

    # Final organic outlet display
    m.fs.solex.mscontactor.organic[time_duration, 1].conc_mass_comp.display()
    m.fs.solex.mscontactor.organic[time_duration, 1].conc_mol_comp.display()

    # Final aqueous outlets display
    m.fs.solex.mscontactor.aqueous[
        time_duration, number_of_stages
    ].conc_mass_comp.display()
    m.fs.solex.mscontactor.aqueous[
        time_duration, number_of_stages
    ].conc_mol_comp.display()
