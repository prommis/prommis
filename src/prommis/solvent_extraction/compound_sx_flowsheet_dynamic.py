from pyomo.environ import (
    ConcreteModel,
    units,
    TransformationFactory,
    Var,
    value,
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
from prommis.solvent_extraction.compound_solvent_extraction import (
    CompoundSolventExtraction,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)

time_duration = 12


def build_model(dosage, number_of_stages, time_duration):
    """
    Build model
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )
    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.compound_solex = CompoundSolventExtraction(
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=True,
        settler_transformation_method="dae.finite_difference",
        settler_transformation_scheme="BACKWARD",
        settler_finite_elements=4,
    )

    return m


def discretization_scheme(m):
    """
    Discretize the model
    """
    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=4, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")


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


# Fixing inlet conditions


def set_inputs(m, dosage, perturb_time):

    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "H2O"].fix(1e6)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "H"].fix(10.75)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "SO4"].fix(100)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "HSO4"].fix(1e4)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Al"].fix(422.375)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Ca"].fix(109.542)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Cl"].fix(1e-7)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Fe"].fix(688.266)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Sc"].fix(0.032)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Y"].fix(0.124)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "La"].fix(0.986)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Ce"].fix(2.277)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Pr"].fix(0.303)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Nd"].fix(0.946)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Sm"].fix(0.097)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Gd"].fix(0.2584)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[:, "Dy"].fix(0.047)

    for t in m.fs.time:
        if t <= perturb_time:
            m.fs.compound_solex.aqueous_inlet.flow_vol[t].fix(62.01)
        else:
            m.fs.compound_solex.aqueous_inlet.flow_vol[t].fix(72.01)

    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Kerosene"].fix(820e3)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Al_o"].fix(1.267e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Ca_o"].fix(2.684e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Fe_o"].fix(2.873e-6)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Sc_o"].fix(1.734)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Y_o"].fix(2.179e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "La_o"].fix(0.000105)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Ce_o"].fix(0.00031)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Pr_o"].fix(3.711e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Nd_o"].fix(0.000165)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Sm_o"].fix(1.701e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Gd_o"].fix(3.357e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[:, "Dy_o"].fix(8.008e-6)

    m.fs.compound_solex.organic_inlet.flow_vol.fix(62.01)

    # Fixing mixer parameters

    m.fs.compound_solex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs.compound_solex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.compound_solex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    # Fixing settler parameters

    m.fs.compound_solex.organic_settler[:].unit.area.fix(1)
    m.fs.compound_solex.aqueous_settler[:].unit.area.fix(1)
    m.fs.compound_solex.aqueous_settler[:].unit.length.fix(0.4)
    m.fs.compound_solex.organic_settler[:].unit.length.fix(0.4)


def set_initial_guess(m):

    for e in m.fs.leach_soln.component_list:
        if e not in ["H2O", "HSO4"]:
            m.fs.compound_solex.mixer[:].unit.mscontactor.aqueous[0, :].conc_mass_comp[
                e
            ].fix()

    m.fs.compound_solex.mixer[:].unit.mscontactor.volume_frac_stream[
        0, :, "aqueous"
    ].fix()
    m.fs.compound_solex.mixer[:].unit.mscontactor.aqueous[0, :].flow_vol.fix()

    m.fs.compound_solex.mixer[:].unit.mscontactor.aqueous_inherent_reaction_extent[
        0.0, :, "Ka2"
    ].fix()

    m.fs.compound_solex.mixer[:].unit.mscontactor.organic[0, :].flow_vol.fix()
    m.fs.compound_solex.mixer[:].unit.mscontactor.organic[0, :].conc_mass_comp[
        "DEHPA"
    ].fix()

    for e in m.fs.reaxn.element_list:
        m.fs.compound_solex.mixer[:].unit.mscontactor.heterogeneous_reaction_extent[
            0.0, :, f"{e}_mass_transfer"
        ].fix()

    for s in m.fs.compound_solex.elements:
        for x in m.fs.compound_solex.aqueous_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.leach_soln.component_list:
                    if e not in ["H2O", "HSO4"]:
                        m.fs.compound_solex.aqueous_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.compound_solex.aqueous_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()
                m.fs.compound_solex.aqueous_settler[s].unit.inherent_reaction_extent[
                    0, x, "Ka2"
                ].fix()

        for x in m.fs.compound_solex.organic_settler[s].unit.length_domain:
            if x != 0:
                for e in m.fs.prop_o.component_list:
                    if e != "Kerosene":
                        m.fs.compound_solex.organic_settler[s].unit.properties[
                            0, x
                        ].conc_mass_comp[e].fix()
                m.fs.compound_solex.organic_settler[s].unit.properties[
                    0, x
                ].flow_vol.fix()


if __name__ == "__main__":

    dosage = 5
    number_of_stages = 3
    time_duration = 12
    perturb_time = 4

    m = build_model(dosage, number_of_stages, time_duration)
    discretization_scheme(m)
    from_json(m, fname="compound_solvent_extraction.json", wts=StoreSpec.value())
    copy_first_steady_state(m)
    set_inputs(m, dosage, perturb_time)
    set_initial_guess(m)

    solver = get_solver("ipopt_v2")
    solver.options["max_iter"] = 10000
    results = solver.solve(m, tee=True)

    percentage_recovery = {}
    REE_set = ["Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]

    for e in REE_set:
        percentage_recovery[e] = [
            (
                1
                - (
                    m.fs.compound_solex.aqueous_outlet.conc_mass_comp[t, e]()
                    * m.fs.compound_solex.aqueous_outlet.flow_vol[t]()
                )
                / (
                    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[t, e]()
                    * m.fs.compound_solex.aqueous_inlet.flow_vol[t]()
                )
            )
            * 100
            for t in m.fs.time
        ]
