#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Demonstration of dynamic flowsheet for LeachTrain unit model using
parameters and data for West Kentucky No. 13 coal refuse.

Authors: Arkoprabho Dasgupta, Akintomiwa Ojo
"""
import matplotlib.pyplot as plt
from pyomo.environ import ConcreteModel, TransformationFactory, units, Var, Constraint
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import FlowsheetBlock
from idaes.core.util import from_json
from idaes.core.solvers import get_solver

from prommis.leaching.leach_train import LeachingTrain
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters


def build_model(time_duration, number_of_tanks):
    """
    Method to build a single stage leaching system using data for
    West Kentucky No. 13 coal refuse.
    Args:
        time_duration: Duration of the simulation in hours.
        number_of_tanks: Number of tanks in the leaching train.
    Returns:
        m: ConcreteModel object with the leaching system.
    """
    m = ConcreteModel()

    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )

    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    m.fs.leach = LeachingTrain(
        number_of_tanks=number_of_tanks,
        liquid_phase={
            "property_package": m.fs.leach_soln,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        solid_phase={
            "property_package": m.fs.coal,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        reaction_package=m.fs.leach_rxns,
    )

    return m


def discretization(m):
    """
    Discretization of the time domain
    """

    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=6, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    """
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


def set_inputs(m, perturb_time):
    """
    Set inlet conditions to leach reactor based on one case study from
    University of Kentucky pilot plant study. The values of the time discrete
    variables at initial time are fixed to the steady state values.
    Args:
        m: ConcreteModel object with the leaching system.
        perturb_time: Time at which the perturbation is applied.
    Returns:
        None
    """

    # Liquid feed state
    for t in m.fs.time:
        if t <= perturb_time:
            m.fs.leach.liquid_inlet.flow_vol[t].fix(224.3 * units.L / units.hour)
        else:
            m.fs.leach.liquid_inlet.flow_vol[t].fix(224.3 * units.L / units.hour)

    m.fs.leach.liquid_inlet.conc_mass_comp.fix(1e-10 * units.mg / units.L)

    m.fs.leach.liquid_inlet.conc_mass_comp[:, "H"].fix(
        2 * 0.05 * 1e3 * units.mg / units.L
    )
    m.fs.leach.liquid_inlet.conc_mass_comp[:, "HSO4"].fix(1e-8 * units.mg / units.L)
    m.fs.leach.liquid_inlet.conc_mass_comp[:, "SO4"].fix(
        0.05 * 96e3 * units.mg / units.L
    )

    # Solid feed state
    m.fs.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach.solid_inlet.mass_frac_comp[:, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[:, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Sc2O3"].fix(
        2.77966e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Y2O3"].fix(
        3.28653e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "La2O3"].fix(
        6.77769e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Ce2O3"].fix(
        0.000156161 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Pr2O3"].fix(
        1.71438e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Nd2O3"].fix(
        6.76618e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Sm2O3"].fix(
        1.47926e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Gd2O3"].fix(
        1.0405e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[:, "Dy2O3"].fix(
        7.54827e-06 * units.kg / units.kg
    )

    # Fixing the volume and the volumetric phases of the leach reactor
    m.fs.leach.volume.fix(100 * units.gallon)
    m.fs.leach.mscontactor.volume.fix(100 * units.gallon)
    m.fs.leach.mscontactor.volume_frac_stream[0, :, "liquid"].fix(0.5)

    @m.Constraint(m.fs.time, m.fs.leach.mscontactor.elements)
    def volume_fraction_rule(m, t, s):
        if t == m.fs.time.first():
            return Constraint.Skip
        else:
            theta_s = m.fs.leach.mscontactor.volume_frac_stream[t, s, "solid"]
            theta_l = m.fs.leach.mscontactor.volume_frac_stream[t, s, "liquid"]
            v_l = m.fs.leach.mscontactor.liquid[t, s].flow_vol
            solid_dens_mass = m.fs.leach.config.solid_phase[
                "property_package"
            ].dens_mass
            v_s = m.fs.leach.mscontactor.solid[t, s].flow_mass / solid_dens_mass
            return v_l * theta_s == v_s * theta_l

    # Fixing the variable values at t=0
    m.fs.leach.mscontactor.liquid[0, :].flow_vol.fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["H"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["HSO4"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Cl"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Sc"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Y"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["La"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Ce"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Pr"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Nd"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Sm"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Gd"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Dy"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Al"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Ca"].fix()
    m.fs.leach.mscontactor.liquid[0, :].conc_mass_comp["Fe"].fix()

    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["inerts"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Y2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["La2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Ce2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Pr2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Nd2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Sm2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Gd2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Dy2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Al2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["CaO"].fix()
    m.fs.leach.mscontactor.solid[0, :].mass_frac_comp["Fe2O3"].fix()
    m.fs.leach.mscontactor.solid[0, :].flow_mass.fix()

    m.fs.leach.mscontactor.liquid_inherent_reaction_extent[0.0, 1, "Ka2"].fix()


if __name__ == "__main__":

    time_duration = 24
    perturb_time = 12
    number_of_tanks = 1

    # Call the build_model function to create the model
    m = build_model(time_duration, number_of_tanks)

    # Discretize the model
    discretization(m)

    # Import steady state values from JSON file
    from_json(m, fname="leaching.json")

    # Initialize the model at steady state values
    copy_first_steady_state(m)

    # Set the inputs for the model
    set_inputs(m, perturb_time)

    # Solve the model
    solver = get_solver("ipopt")
    solver.solve(m, tee=True)

    # Solid stream outlet values at final time
    m.fs.leach.mscontactor.solid[time_duration, number_of_tanks].flow_mass.pprint()
    m.fs.leach.mscontactor.solid[time_duration, number_of_tanks].mass_frac_comp.pprint()

    # Liquid stream outlet values at final time
    m.fs.leach.mscontactor.liquid[time_duration, number_of_tanks].flow_vol.pprint()
    m.fs.leach.mscontactor.liquid[
        time_duration, number_of_tanks
    ].conc_mass_comp.pprint()

    # Plotting the results for REE oxides
    REE_set = m.fs.coal.component_list - ["inerts", "Al2O3", "Fe2O3", "CaO"]
    for e in REE_set:
        plt.plot(m.fs.time, m.fs.leach.recovery[:, e]())
    plt.legend(REE_set)
    plt.xlabel("Time (h)")
    plt.ylabel("Recovery %")
    plt.axvline(
        x=perturb_time,
        color="r",
        linestyle="--",
        label="Perturbation at t=12h",
    )
    plt.title("REE oxide recovery variation wrt time, with perturbation at t=12 ")
    plt.figure()
    plt.plot(m.fs.time, m.fs.leach.mscontactor.solid[:, :].flow_mass())
    plt.title("Solid mass flow rate variation wrt time")
    plt.xlabel("Time (h)")
    plt.ylabel("Solid mass flow rate (kg/h)")
