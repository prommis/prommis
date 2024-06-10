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
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    TransformationFactory,
    floor,
    log10,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed, Product

import matplotlib.pyplot as plt
import numpy as np
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    ActivityCoefficientModel,
    DensityCalculation,
    MCASParameterBlock,
)
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
from watertap.unit_models.pressure_changer import Pump

_log = idaeslog.getLogger(__name__)


def main():
    """
    Builds and solves the NF flowsheet

    Returns:
        m: pyomo model
    """

    # initialize lists to store sensitivity data
    (
        area,
        li_rejection,
        mg_rejection,
        mg_li_ratio,
        feed_ratio,
        feed_pressure,
        recovery_vals,
    ) = initialize_sensitivity()

    # solve the system at each point
    for recovery in recovery_vals:
        solver = get_solver()
        m = build()

        initialize(m, solver)
        _log.info("Initialization Okay")

        if degrees_of_freedom(m) != 0:
            raise ValueError("Degrees of freedom were not equal to zero")
        solve_model(m, solver)
        _log.info("Solved Box Problem")

        unfix_optimization_variables(m)
        add_objective(m)
        add_pressure_constraint(m, pressure_limit=None)
        add_recovery_constraint(m, recovery_limit=recovery)
        solve_model(m, solver)
        collect_plot_data(
            m, area, li_rejection, mg_rejection, mg_li_ratio, feed_ratio, feed_pressure
        )
        print_information(m)

    # create the sensitvity analysis plots using the data collected above
    plot(
        area,
        li_rejection,
        mg_rejection,
        mg_li_ratio,
        feed_ratio,
        feed_pressure,
        recovery_vals,
    )

    return m


def set_default_feed(m, solver):
    """
    Fixes the concentrations used to initialize the feed using the
    concentration of the Salar de Atacama (kg/m3 = g/L)
    Note: Cl- concentration will get overridden to enforce electroneutrality

    Args:
        m: pyomo model
        solver: optimization solver
    """
    conc_mass_phase_comp = {"Li_+": 1.19, "Mg_2+": 7.31, "Cl_-": 143.72}
    set_nf_feed(
        blk=m.fs,
        solver=solver,
        flow_mass_h2o=1,  # arbitrary for now
        conc_mass_phase_comp=conc_mass_phase_comp,
    )


def define_feed_composition():
    """
    Returns the ion properties needed for the DSPM-DE property package

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

    Returns:
        m: pyomo model
    """
    # create the model
    m = ConcreteModel()

    # create the flowsheet
    m.fs = FlowsheetBlock(dynamic=False)

    # define the property model
    default = define_feed_composition()
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


def fix_initial_variables(m):
    """
    Fixes the initial variables needed to create 0 DOF

    Args:
        m: pyomo model
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


def unfix_optimization_variables(m):
    """
    Unfixes select variables to enable optimization with DOF>0

    Args:
        m: pyomo model
    """
    m.fs.pump.outlet.pressure[0].unfix()
    m.fs.unit.area.unfix()


def add_objective(m):
    """
    Adds objective to the pyomo model

    Args:
        m: pyomo model
    """
    # limit Li loss
    m.fs.objective = Objective(
        expr=m.fs.retentate.flow_mol_phase_comp[0, "Liq", "Li_+"]
    )


def add_pressure_constraint(m, pressure_limit):
    """
    Adds feed pressure constraint to the pyomo model

    Args:
        m: pyomo model
        pressure_limit: upper bound on the outlet pump pressure
    """
    if pressure_limit is None:
        pressure_limit = 7e6
    # bound the feed pressure to a reasonable value for nanofiltration
    # choose an upper limit of 70 bar (https://doi.org/10.1021/acs.est.2c08584)
    m.fs.pressure_constraint = Constraint(
        expr=m.fs.pump.outlet.pressure[0] <= pressure_limit
    )


def add_recovery_constraint(m, recovery_limit):
    """
    Adds recovery constraint to the pyomo model

    Args:
        m: pyomo model
        recovery_limit: upper bound on the volume recovery
    """
    if recovery_limit is None:
        recovery_limit = 0.8
    # limit the NF recovery
    m.fs.recovery_constraint = Constraint(
        expr=m.fs.unit.recovery_vol_phase[0.0, "Liq"] <= recovery_limit
    )


def solve_model(m, solver):
    """
    Optimizes the flowsheet

    Args:
        m: pyomo model
        solver: optimization solver
    """
    _log.info(f"Optimizing with {format(degrees_of_freedom(m))} DOFs")
    simulation_results = solver.solve(m, tee=True)
    if simulation_results.solver.termination_condition != "optimal":
        raise ValueError("The solver did not return optimal termination")
    return simulation_results


def initialize_sensitivity():
    """
    Makes plots to perform a sensitivity analysis on the nanofiltration flowsheet

    Returns:
        area: list to store optimal membrane area (m2)
        li_rejection: list to store lithium rejection
        mg_rejection: list to store magnesium rejection
        mg_li_ratio: list to store Mg:Li mass ratio of the permeate
        feed_ratio: list to store Mg:Li mass ratio of the feed
        feed_pressure: list to store the optimal feed pressure (bar)
        recovery_vals: list that holds the volume recovery values to test
    """
    # initialize lists to store data
    area = []  # m2
    li_rejection = []
    mg_rejection = []
    mg_li_ratio = []
    feed_ratio = []
    feed_pressure = []  # bar

    # provide values to constrain
    recovery_vals = np.arange(0.2, 1, 0.1)
    return (
        area,
        li_rejection,
        mg_rejection,
        mg_li_ratio,
        feed_ratio,
        feed_pressure,
        recovery_vals,
    )


def collect_plot_data(
    m, area, li_rejection, mg_rejection, mg_li_ratio, feed_ratio, feed_pressure
):
    """
    Stores the relevant information after each flowsheet solve to prepare plots

    Args:
        m: pyomo model
        area: list to store optimal membrane area (m2)
        li_rejection: list to store lithium rejection
        mg_rejection: list to store magnesium rejection
        mg_li_ratio: list to store Mg:Li mass ratio of the permeate
        feed_ratio: list to store Mg:Li mass ratio of the feed
        feed_pressure: list to store the optimal feed pressure (bar)
    """
    area.append(m.fs.unit.area.value)
    li_rejection.append(
        m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Li_+"].value
    )
    mg_rejection.append(
        m.fs.unit.rejection_intrinsic_phase_comp[0, "Liq", "Mg_2+"].value
    )
    mg_li_ratio.append(
        (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
        / (m.fs.permeate.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069)
    )
    feed_ratio.append(
        (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Mg_2+"].value / 0.024)
        / (m.fs.feed.flow_mol_phase_comp[0, "Liq", "Li_+"].value / 0.0069)
    )
    feed_pressure.append(m.fs.pump.outlet.pressure[0].value / 1e5)


def print_information(m):
    """
    Prints relevant information about the system
    """
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


def plot(
    area,
    li_rejection,
    mg_rejection,
    mg_li_ratio,
    feed_ratio,
    feed_pressure,
    recovery_vals,
):
    """
    Creates four subplots of the nanofiltration system, reporting
    ion rejection, Mg:Li ratio, membrane area, and feed pressure
    as the volume recovery of the membrane changes

    Args:
        area: list to store optimal membrane area (m2)
        li_rejection: list to store lithium rejection
        mg_rejection: list to store magnesium rejection
        mg_li_ratio: list to store Mg:Li mass ratio of the permeate
        feed_ratio: list to store Mg:Li mass ratio of the feed
        feed_pressure: list to store the optimal feed pressure (bar)
        recovery_vals: list that holds the volume recovery values to test
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    ax1.plot(recovery_vals, li_rejection, "-o")
    ax1.plot(recovery_vals, mg_rejection, "-o")
    ax1.legend(["Li+", "Mg2+"])
    ax1.set_title("Ion Rejection vs Volume Recovery")
    ax1.set_ylabel("Rejection")

    ax2.plot(recovery_vals, area, "-o")
    ax2.set_title("Membrane Area vs Volume Recovery")
    ax2.set_ylabel("Area (m2)")

    ax3.plot(recovery_vals, mg_li_ratio, "-o")
    ax3.plot(recovery_vals, feed_ratio, "r")
    ax3.legend(["Permeate", "Feed"])
    ax3.set_title("Mg:Li Mass Ratio vs Volume Recovery")
    ax3.set_xlabel("Recovery")
    ax3.set_ylabel("Mg:Li")

    ax4.plot(recovery_vals, feed_pressure, "-o")
    ax4.set_title("Feed Pressure vs Volume Recovery")
    ax4.set_xlabel("Recovery")
    ax4.set_ylabel("Feed Pressure (bar)")

    plt.show()


def initialize(m, solver):
    """
    Initializes the flowsheet units

    Args:
        m: pyomo model
        solver: optimization solver
    """
    set_default_feed(m, solver)
    fix_initial_variables(m)

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

    Args:
        blk: flowsheet block
        solver: optimization solver
        flow_mass_h2o: inlet water flow rate (feed)
        conc_mass_phase_conc: mass concentration (feed)
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


def calculate_scale(value):
    """
    Calculates a default scaling value
    """
    return -1 * floor(log10(value))


def set_nf_feed_scaling(blk):
    """
    Calculates the default scaling for the feed solution
    """
    _add = 0
    for i in blk.feed.properties[0].flow_mol_phase_comp:
        scale = calculate_scale(blk.feed.properties[0].flow_mol_phase_comp[i].value)
        print(f"{i} flow_mol_phase_comp scaling factor = {10**(scale+_add)}")
        blk.properties.set_default_scaling(
            "flow_mol_phase_comp", 10 ** (scale + _add), index=i
        )


if __name__ == "__main__":
    main()
