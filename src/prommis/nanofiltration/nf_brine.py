####################################################################
# This code has been adapted from a test in the WaterTAP repo to
# simulate a simple separation of lithium and magnesium
#
# watertap > unit_models > tests > test_nanofiltration_DSPMDE_0D.py
# test defined as test_pressure_recovery_step_2_ions()
#
# also used the following flowsheet as a reference
# watertap > examples > flowsheets > nf_dspmde > nf.py
#
# https://github.com/watertap-org/watertap/blob/main/tutorials/nawi_spring_meeting2023.ipynb
####################################################################

"""
Nanofiltration flowsheet for Donnan steric pore model with dielectric exclusion
"""

# import statements
from math import floor, log

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    TransformationFactory,
    assert_optimal_termination,
    maximize,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    ActivityCoefficientModel,
    DensityCalculation,
    MCASParameterBlock,
)
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
from watertap.unit_models.pressure_changer import Pump


def main():
    """
    Builds and solves the NF flowsheet
    """
    solver = get_solver()
    m = build()

    initialize(m, solver)
    print("init_okay")

    assert degrees_of_freedom(m) == 0
    optimize(m, solver)
    print("solved box problem")
    m.fs.unit.report()

    unfix_opt_vars(m)
    add_obj(m)
    add_con(m, pressure_limit=7e6)
    optimize(m, solver)
    m.fs.unit.report()
    print("Optimal NF feed pressure (Bar)", m.fs.pump.outlet.pressure[0].value / 1e5)
    print("Optimal area (m2)", m.fs.unit.area.value)
    print(
        "Optimal NF vol recovery (%)",
        m.fs.unit.recovery_vol_phase[0.0, "Liq"].value * 100,
    )
    print(
        "Optimal Li rejection (%)",
        m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Li_+"].value * 100,
    )
    print(
        "Optimal Mg rejection (%)",
        m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Mg_2+"].value * 100,
    )
    print(
        "Feed Mg:Li ratio (mass)",
        (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
        / (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
    )
    print(
        "Permeate Mg:Li ratio (mass)",
        (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
        / (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069),
    )

    return m


def set_default_feed(m, solver):
    """
    Fixes the concentrations used to initialize the feed

    Approximates the concentration of Salar de Atacama (kg/m3 = g/L)

    Cl- concentration will get overridden to enforce electroneutrality
    """
    conc_mass_phase_comp = {"Li_+": 1.19, "Mg_2+": 7.31, "Cl_-": 143.72}
    set_nf_feed(
        blk=m.fs,
        solver=solver,
        flow_mass_h2o=1,  # arbitrary for now
        conc_mass_phase_comp=conc_mass_phase_comp,
    )


def define_feed_comp():
    """
    Defines the ion properties needed for the DSPM-DE property package

    Ions include lithium, magnesium, and chloride, assuming LiCl and MgCl2 salts

    diffusivity:
    - https://www.aqion.de/site/diffusion-coefficients
    - very confident

    molecular weights:
    - very confident

    Stokes radius:
    - average values from https://www.sciencedirect.com/science/article/pii/S138358661100637X
    - medium confident (averaged values from multiple studies)
    - reasonable orders of magnitude

    ion charge:
    - very confident

    The activity coefficient options are ideal or davies
    """
    default = {
        "solute_list": ["Li_+", "Mg_2+", "Cl_-"],
        "diffusivity_data": {
            ("Liq", "Li_+"): 1.03e-09,
            ("Liq", "Mg_2+"): 0.705e-09,
            ("Liq", "Cl_-"): 2.03e-09,
        },
        "mw_data": {"H2O": 0.018, "Li_+": 0.0069, "Mg_2+": 0.024, "Cl_-": 0.035},
        "stokes_radius_data": {
            "Li_+": 3.61e-10,
            # "Mg_2+": 4.07e-10,
            # "Cl_-": 3.28e-10
            # adjusted Cl and Mg to values from nf.py'
            "Cl_-": 0.121e-9,
            "Mg_2+": 0.347e-9,
        },
        "charge": {"Li_+": 1, "Mg_2+": 2, "Cl_-": -1},
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
    }
    return default


def build():
    """
    Builds the NF flowsheet
    """
    # create the model
    m = ConcreteModel()

    # create the flowsheet
    m.fs = FlowsheetBlock(dynamic=False)

    # define the property model
    default = define_feed_comp()
    m.fs.properties = MCASParameterBlock(**default)

    # add the feed and product streams
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.permeate = Product(property_package=m.fs.properties)
    m.fs.retentate = Product(property_package=m.fs.properties)

    # define unit models
    m.fs.pump = Pump(property_package=m.fs.properties)
    m.fs.unit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)

    # connect the streams and blocks
    m.fs.feed_to_pump = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.pump_to_nf = Arc(source=m.fs.pump.outlet, destination=m.fs.unit.inlet)
    m.fs.nf_to_permeate = Arc(
        source=m.fs.unit.permeate, destination=m.fs.permeate.inlet
    )
    m.fs.nf_to_retentate = Arc(
        source=m.fs.unit.retentate, destination=m.fs.retentate.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def fix_init_vars(m):
    """
    Fixes the initial variables needed to create 0 DOF
    """

    # pump variables
    m.fs.pump.efficiency_pump[0].fix(0.75)
    m.fs.pump.control_volume.properties_in[0].temperature.fix(298.15)
    m.fs.pump.control_volume.properties_in[0].pressure.fix(101325)
    m.fs.pump.outlet.pressure[0].fix(2e5)
    iscale.set_scaling_factor(m.fs.pump.control_volume.work, 1e-4)

    # membrane operation
    m.fs.unit.recovery_vol_phase[0, "Liq"].setub(0.95)
    m.fs.unit.spacer_porosity.fix(0.85)
    m.fs.unit.channel_height.fix(5e-4)
    m.fs.unit.velocity[0, 0].fix(0.1)
    m.fs.unit.area.fix(100)
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)

    # variables for calculating mass transfer coefficient with spiral wound correlation
    m.fs.unit.spacer_mixing_efficiency.fix()
    m.fs.unit.spacer_mixing_length.fix()

    # membrane properties
    m.fs.unit.radius_pore.fix(0.5e-9)
    m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
    m.fs.unit.membrane_charge_density.fix(-60)
    m.fs.unit.dielectric_constant_pore.fix(41.3)
    iscale.calculate_scaling_factors(m)


def unfix_opt_vars(m):
    """
    Unfixes select variables to enable optimization with DOF>0
    """
    m.fs.pump.outlet.pressure[0].unfix()
    m.fs.unit.area.unfix()


def add_obj(m):
    """
    Adds objectives to the pyomo model
    """
    # limit Li loss
    # m.fs.obj = Objective(
    #     expr = m.fs.retentate.flow_mol_phase_comp[0,"Liq", "Li_+"],
    #     # sense = maximize
    # )

    # maxmize the reduction in Mg:Li ratio
    m.fs.obj = Objective(
        expr=(
            (
                (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
                / (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069)
            )
            - (
                (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
                / (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069)
            )
        ),
        sense=maximize,
    )


def add_con(m, pressure_limit):
    """
    Adds constraints to the pyomo model
    """
    # # limit the Li rejection
    # m.fs.li_rejection_con = Constraint(
    #     expr=m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Li_+"] >= 0.2
    # )

    # bound the feed pressure to a reasonable value for nanofiltration
    # choose an upper limit of 70 bar (https://doi.org/10.1021/acs.est.2c08584)
    m.fs.pressure_con = Constraint(expr=m.fs.pump.outlet.pressure[0] <= pressure_limit)


def optimize(m, solver):
    """
    Optimizes the flowsheet
    """
    print(f"Optimizing with {format(degrees_of_freedom(m))} DOFs")
    simulation_results = solver.solve(m, tee=True)
    assert_optimal_termination(simulation_results)
    return simulation_results


def initialize(m, solver):
    """
    Initializes the flowsheet units
    """
    set_default_feed(m, solver)
    fix_init_vars(m)

    m.fs.feed.initialize(optarg=solver.options)
    propagate_state(m.fs.feed_to_pump)

    m.fs.pump.initialize(optarg=solver.options)
    propagate_state(m.fs.pump_to_nf)

    m.fs.unit.initialize(optarg=solver.options)
    propagate_state(m.fs.nf_to_permeate)
    propagate_state(m.fs.nf_to_retentate)

    m.fs.permeate.initialize(optarg=solver.options)
    m.fs.retentate.initialize(optarg=solver.options)


def set_nf_feed(blk, solver, flow_mass_h2o, conc_mass_phase_comp):  # kg/m3
    """
    Calculates the concentration of the feed solution in molar flow rate
    """
    if solver is None:
        solver = get_solver()

    # fix the inlet flow to the block as water flowrate
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(flow_mass_h2o)

    # fix the ion cncentrations and unfix ion flows
    for ion, x in conc_mass_phase_comp.items():
        blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].fix(x)
        blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].unfix()
    # solve for the new flow rates
    solver.solve(blk.feed)
    # fix new water concentration
    blk.feed.properties[0].conc_mass_phase_comp["Liq", "H2O"].fix()
    # unfix ion concentrations and fix flows
    for ion, x in conc_mass_phase_comp.items():
        blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].unfix()
        blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix()
        blk.feed.properties[0].flow_mass_phase_comp["Liq", ion].unfix()
    blk.feed.properties[0].conc_mass_phase_comp["Liq", "H2O"].unfix()
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    blk.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix()

    set_nf_feed_scaling(blk)

    # assert electroneutrality
    blk.feed.properties[0].assert_electroneutrality(
        defined_state=True, adjust_by_ion="Cl_-", get_property="flow_mol_phase_comp"
    )

    # switching to concentration for ease of adjusting in UI
    # addresses error in fixing flow_mol_phase_comp
    for ion, x in conc_mass_phase_comp.items():
        blk.feed.properties[0].conc_mass_phase_comp["Liq", ion].unfix()
        blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix()


def calc_scale(value):
    """
    Calculates a default scaling value
    """
    return -1 * floor(log(value, 10))


def set_nf_feed_scaling(blk):
    """
    Calculates the default scaling for the feed solution
    """
    _add = 0
    for i in blk.feed.properties[0].flow_mol_phase_comp:
        scale = calc_scale(blk.feed.properties[0].flow_mol_phase_comp[i].value)
        print(f"{i} flow_mol_phase_comp scaling factor = {10**(scale+_add)}")
        blk.properties.set_default_scaling(
            "flow_mol_phase_comp", 10 ** (scale + _add), index=i
        )


if __name__ == "__main__":
    main()
