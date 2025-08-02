#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    ConcreteModel,
    units,
    TransformationFactory,
    Var,
)
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import from_json, StoreSpec

from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.mixer_settler_extraction import (
    MixerSettlerExtraction,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)

time_duration = 12


def build_model(dosage, number_of_stages, time_duration):
    """
    Method to build a dynamic flowsheet for mixer settler solvent extraction.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )
    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.mixer_settler_ex = MixerSettlerExtraction(
        number_of_stages=number_of_stages,
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

    return m


def discretization_scheme(m):
    """
    Discretize the mixer settler solvent extraction model
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=4, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


# Fixing inlet conditions


def set_inputs(m, dosage, perturb_time):
    """
    Set inlet conditions to the mixer settler solvent extraction model and fixing
    the parameters of the model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        None

    """

    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H2O"].fix(1e6)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "H"].fix(10.75)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "SO4"].fix(100)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "HSO4"].fix(1e4)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Al"].fix(422.375)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ca"].fix(109.542)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Cl"].fix(1e-7)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Fe"].fix(688.266)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sc"].fix(0.032)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Y"].fix(0.124)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "La"].fix(0.986)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Ce"].fix(2.277)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Pr"].fix(0.303)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Nd"].fix(0.946)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Sm"].fix(0.097)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Gd"].fix(0.2584)
    m.fs.mixer_settler_ex.aqueous_inlet.conc_mass_comp[:, "Dy"].fix(0.047)

    for t in m.fs.time:
        if t <= perturb_time:
            m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(62.01)
        else:
            m.fs.mixer_settler_ex.aqueous_inlet.flow_vol[t].fix(72.01)

    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Kerosene"].fix(820e3)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Al_o"].fix(1.267e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ca_o"].fix(2.684e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Fe_o"].fix(2.873e-6)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sc_o"].fix(1.734)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Y_o"].fix(2.179e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "La_o"].fix(0.000105)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Ce_o"].fix(0.00031)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Pr_o"].fix(3.711e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Nd_o"].fix(0.000165)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Sm_o"].fix(1.701e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Gd_o"].fix(3.357e-5)
    m.fs.mixer_settler_ex.organic_inlet.conc_mass_comp[:, "Dy_o"].fix(8.008e-6)

    m.fs.mixer_settler_ex.organic_inlet.flow_vol.fix(62.01)

    # Fixing mixer parameters

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    # Fixing settler parameters

    m.fs.mixer_settler_ex.organic_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.area.fix(1)
    m.fs.mixer_settler_ex.aqueous_settler[:].unit.length.fix(0.4)
    m.fs.mixer_settler_ex.organic_settler[:].unit.length.fix(0.4)


def set_initial_conditions(m):
    """
    Set initial conditions at time=0 for the mixer-settler solvent extraction model
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None

    """

    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4"]:
            m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[
                0, :
            ].conc_mass_comp[e].fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.volume_frac_stream[
        0, :, "aqueous"
    ].fix()
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous[0, :].flow_vol.fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.aqueous_inherent_reaction_extent[
        0.0, :, "Ka2"
    ].fix()

    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].flow_vol.fix()
    m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.organic[0, :].conc_mass_comp[
        "DEHPA"
    ].fix()

    for e in m.fs.reaxn.element_list:
        m.fs.mixer_settler_ex.mixer[:].unit.mscontactor.heterogeneous_reaction_extent[
            0.0, :, f"{e}_mass_transfer"
        ].fix()

    # set variable values in the settler at t=0

    for s in m.fs.mixer_settler_ex.elements:
        for x in m.fs.mixer_settler_ex.aqueous_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.leach_soln.component_list:
                    if e not in ["H2O", "HSO4"]:
                        m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.mixer_settler_ex.aqueous_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()
                m.fs.mixer_settler_ex.aqueous_settler[s].unit.inherent_reaction_extent[
                    0, x, "Ka2"
                ].fix()

        for x in m.fs.mixer_settler_ex.organic_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.prop_o.component_list:
                    if e not in ["Kerosene"]:
                        m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.mixer_settler_ex.organic_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()


def build_model_and_discretize(dosage, number_of_stages, time_duration):
    """
    Method to build a dynamic model for mixer settler solvent extraction and discretize
    the model.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    """

    m = build_model(dosage, number_of_stages, time_duration)
    discretization_scheme(m)

    return m


def import_steady_value(m, path_name):
    """
    A function to import the steady state values of the mixer-settler solvent extraction
    model to the dynamic model for initializing it.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        path_name: name of the path of the json file
    Returns:
        None
    """
    from_json(m, fname=path_name, wts=StoreSpec.value())


def initialize_set_input_and_initial_conditions(m, dosage, perturb_time):
    """
    Function to initialize, set inlet values and give initial conditions to the dynamic
    mixer settler solvent extraction model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
        dosage: percentage dosage of extractant to the system.
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        None

    """

    copy_first_steady_state(m)
    set_inputs(m, dosage, perturb_time)
    set_initial_conditions(m)


def solve_model(m):
    """
    A function to solve the initialized mixer-settler solvent extraction model.
    Args:
        m: ConcreteModel object with the mixer-settler solvent extraction system.
    Returns:
        None
    """
    solver = get_solver("ipopt_v2")
    results = solver.solve(m, tee=True)
    return results


def main(dosage, number_of_stages, time_duration, perturb_time, path_name):
    """
    Function to build a dynamic model, discretize it, initialize it, set input values
    and give initial conditions, then solve the model.
    Args:
        dosage: percentage dosage of extractant to the system.
        number_of_stages: number of stages in the mixer settler model.
        time_duration = total time of operation of the model
        path_name: name of the path of the json file
        perturb_time : time at which a perturbation is added in the flowsheet, should
        be lesser than the time of operation.
    Returns:
        m: ConcreteModel object with the mixer-settler solvent extraction system.

    """
    m = build_model_and_discretize(dosage, number_of_stages, time_duration)
    import_steady_value(m, path_name)
    initialize_set_input_and_initial_conditions(m, dosage, perturb_time)
    results = solve_model(m)

    return m, results


dosage = 5
number_of_stages = 3
time_duration = 12
perturb_time = 4

if __name__ == "__main__":
    m, results = main(
        dosage,
        number_of_stages,
        time_duration,
        perturb_time,
        path_name="mixer_settler_extraction.json",
    )
