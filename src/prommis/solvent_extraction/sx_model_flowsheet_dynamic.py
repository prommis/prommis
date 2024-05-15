""" 
Demonstration flowsheet for dynamic state solvent extraction loading process
using parameters and data derived from West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import ConcreteModel, SolverFactory, units, TransformationFactory

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

from idaes.core.util.model_diagnostics import DiagnosticsToolbox


def build_model():

    """
    Method of building a dynamic solvent extraction model with a specified number of 
    stages and with two separate property packages for the two inlet streams.
    This is a loading operation, so no additional argument has to be specified.

    """

    m = ConcreteModel()

    time_duration = 60

    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, time_duration], time_units=units.hour)

    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()

    number_of_stages = 1

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

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=2, wrt=m.fs.time, scheme="FORWARD")
    #m.discretizer.apply_to(m, nfe=2, wrt=m.fs.time, scheme="BACKWARD")

    """
    Specifications of the partition coefficients, volume and volume fractions for all
    the stages.

    """

    m.fs.solex.mscontactor.volume[:].fix(0.4)

    m.fs.solex.mscontactor.volume_frac_stream[:, :, "organic"].fix(0.4)

    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sc"] = 99.9 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "La"] = 75.2 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Ce"] = 95.7 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Pr"] = 96.5 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Nd"] = 99.2 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Gd"] = 98.6 / 100
    m.fs.solex.partition_coefficient[:, "aqueous", "organic", "Dy"] = 99.9 / 100

    """
    Fixation of the inlet conditions and the initial state values for all the components.

    """

    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(820)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(5230)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(270)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(209.31)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(637.74)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(2032.77)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(4516.13)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(756.64)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(2047.85)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(369.1)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(174.38)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(101.12)

    # m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["H2O"].fix(1e6)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(1e-4)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(1e-4)
    #m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["H"].fix(1e-4)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(1e-9)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(1e-9)

    #m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0, :, "Ka2"].fix(0)

    m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(4.4)

    #m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(4.4)

    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(7.54e-10)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(4.955e-9)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(1.491e-7)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(321.34)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(5.67e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(1.78e-05)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(4.019e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(6.73e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(1.82e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(3.285e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(1.55e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(9e-7)

    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Y"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ce"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Nd"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sm"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Gd"].fix(1e-9)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Dy"].fix(1e-9)

    m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

    #m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)


if __name__ == "__main__":

    time_duration = 60
    number_of_stages = 1

    # Call build model function

    m = build_model()
    set_inputs(m)

    dt = DiagnosticsToolbox(m)
    dt.report_structural_issues()

    """
    Initialization of the model, which gives a good starting point.

    """

    initializer = BlockTriangularizationInitializer()
    initializer.initialize(m.fs.solex)
    assert initializer.summary[m.fs.solex]["status"] == InitializationStatus.Ok

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
    m.fs.solex.mscontactor.aqueous[time_duration, number_of_stages].conc_mass_comp.display()
    m.fs.solex.mscontactor.aqueous[time_duration, number_of_stages].conc_mol_comp.display()
