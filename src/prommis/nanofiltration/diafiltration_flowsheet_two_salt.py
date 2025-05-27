#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Sample flowsheet for the diafiltration cascade.

Author: Molly Dougher
"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TransformationFactory,
    assert_optimal_termination,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.unit_models import Feed, Product

import matplotlib.pyplot as plt
from pandas import DataFrame

from prommis.nanofiltration.diafiltration_stream_properties import (
    DiafiltrationStreamParameter,
)
from prommis.nanofiltration.diafiltration_solute_properties import SoluteParameter
from prommis.nanofiltration.diafiltration_two_salt import TwoSaltDiafiltration


def main():
    """
    Builds and solves flowsheet with two-salt diafiltration unit model.
    """
    # build flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.stream_properties = DiafiltrationStreamParameter()
    m.fs.properties = SoluteParameter()

    # update parameter inputs if desired
    build_membrane_parameters(m)

    # add feed blocks for feed and diafiltrate
    m.fs.feed_block = Feed(property_package=m.fs.stream_properties)
    m.fs.diafiltrate_block = Feed(property_package=m.fs.stream_properties)

    # add the membrane unit model
    m.fs.membrane = TwoSaltDiafiltration(
        property_package=m.fs.properties,
        NFEx=10,
        NFEz=5,
    )

    # add product blocks for retentate and permeate
    m.fs.retentate_block = Product(property_package=m.fs.stream_properties)
    m.fs.permeate_block = Product(property_package=m.fs.stream_properties)

    # fix the degrees of freedom to their default values
    fix_variables(m)

    # add and connect flowsheet streams
    add_and_connect_streams(m)

    # check structural warnings
    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    # solve model
    solve_model(m)

    # check numerical warnings
    dt.assert_no_numerical_warnings()

    # visualize the results
    # currently no partitioning in the membrane
    plot_results(m)
    plot_membrane_results(m)


def build_membrane_parameters(m):
    """
    Updates parameters needed in two salt diafiltration unit model if desired

    Args:
        m: Pyomo model
    """
    pass


def fix_variables(m):
    # fix degrees of freedom in the membrane
    m.fs.membrane.membrane_width.fix()
    m.fs.membrane.membrane_length.fix()
    m.fs.membrane.applied_pressure.fix()

    # fix degrees of freedom in the flowsheet
    m.fs.membrane.feed_flow_volume.fix()
    m.fs.membrane.feed_conc_mass_comp[0, "Li"].fix()
    m.fs.membrane.feed_conc_mass_comp[0, "Co"].fix()
    m.fs.membrane.feed_conc_mass_comp[0, "Cl"].fix()

    m.fs.membrane.diafiltrate_flow_volume.fix()
    m.fs.membrane.diafiltrate_conc_mass_comp[0, "Li"].fix()
    m.fs.membrane.diafiltrate_conc_mass_comp[0, "Co"].fix()
    m.fs.membrane.diafiltrate_conc_mass_comp[0, "Cl"].fix()


def add_and_connect_streams(m):
    m.fs.feed_stream = Arc(
        source=m.fs.feed_block.outlet,
        destination=m.fs.membrane.feed_inlet,
    )
    m.fs.diafiltrate_stream = Arc(
        source=m.fs.diafiltrate_block.outlet,
        destination=m.fs.membrane.diafiltrate_inlet,
    )
    m.fs.retentate_stream = Arc(
        source=m.fs.membrane.retentate_outlet,
        destination=m.fs.retentate_block.inlet,
    )
    m.fs.permeate_stream = Arc(
        source=m.fs.membrane.permeate_outlet,
        destination=m.fs.permeate_block.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def solve_model(m):
    """
    Solves scaled model.

    Args:
        m: Pyomo model
    """
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=True)
    assert_optimal_termination(results)

    scaling.propagate_solution(scaled_model, m)


def plot_results(m):
    """
    Plots concentration and flux variables across the width of the membrane module.

    Args:
        m: Pyomo model
    """
    # store values for x-coordinate
    x_axis_values = []

    # store values for concentration of lithium in the retentate
    conc_ret_lith = []
    # store values for concentration of lithium in the permeate
    conc_perm_lith = []
    # store values for concentration of cobalt in the retentate
    conc_ret_cob = []
    # store values for concentration of cobalt in the permeate
    conc_perm_cob = []

    # store values for water flux across membrane
    water_flux = []
    # store values for mass flux of lithium across membrane
    lithium_flux = []

    # store values for percent recovery
    percent_recovery = []

    # store values for lithium rejection
    lithium_rejection = []
    # store values for lithium solute passage
    lithium_sieving = []
    # store values for cobalt rejection
    cobalt_rejection = []
    # store values for cobalt solute passage
    cobalt_sieving = []

    for x_val in m.fs.membrane.x_bar:
        x_axis_values.append(x_val * value(m.fs.membrane.membrane_width))
        conc_ret_lith.append(
            value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Li"])
        )
        conc_perm_lith.append(
            value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Li"])
        )
        conc_ret_cob.append(
            value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Co"])
        )
        conc_perm_cob.append(
            value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Co"])
        )

        water_flux.append(value(m.fs.membrane.volume_flux_water[x_val]))
        lithium_flux.append(value(m.fs.membrane.mass_flux_lithium[x_val]))

        lithium_rejection.append(
            (
                1
                - (
                    value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Li"])
                    / value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Li"])
                )
            )
            * 100
        )
        lithium_sieving.append(
            (
                value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Li"])
                / value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Li"])
            )
        )
        cobalt_rejection.append(
            (
                1
                - (
                    value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Co"])
                    / value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Co"])
                )
            )
            * 100
        )
        cobalt_sieving.append(
            (
                value(m.fs.membrane.permeate_conc_mass_comp[0, x_val, "Co"])
                / value(m.fs.membrane.retentate_conc_mass_comp[0, x_val, "Co"])
            )
        )

        percent_recovery.append(
            (
                value(m.fs.membrane.permeate_flow_volume[0, x_val])
                / value(m.fs.membrane.feed_flow_volume[0])
                * 100
            )
        )

    fig1, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        3, 2, dpi=125, figsize=(12, 10)
    )

    ax1.plot(x_axis_values, conc_ret_lith, linewidth=2, label="retentate")
    ax1.plot(x_axis_values, conc_perm_lith, linewidth=2, label="permeate")
    ax1.set_ylabel(
        "Lithium Concentration (kg/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.tick_params(direction="in", labelsize=10)
    ax1.legend()

    ax2.plot(x_axis_values, conc_ret_cob, linewidth=2, label="retentate")
    ax2.plot(x_axis_values, conc_perm_cob, linewidth=2, label="permeate")
    ax2.set_ylabel(
        "Cobalt Concentration (kg/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax2.tick_params(direction="in", labelsize=10)
    ax2.legend()

    ax3.plot(x_axis_values, water_flux, linewidth=2)
    ax3.set_xlabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    ax3.set_ylabel("Water Flux (m$^3$/m$^2$/h)", fontsize=10, fontweight="bold")
    ax3.tick_params(direction="in", labelsize=10)

    ax4.plot(x_axis_values, lithium_flux, linewidth=2)
    ax4.set_xlabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    ax4.set_ylabel(
        "Mass Flux of Lithium\n (kg/m$^2$/h)", fontsize=10, fontweight="bold"
    )
    ax4.tick_params(direction="in", labelsize=10)

    ax5.plot(x_axis_values, lithium_rejection, linewidth=2, label="lithium")
    ax5.plot(x_axis_values, cobalt_rejection, linewidth=2, label="cobalt")
    ax5.set_xlabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    ax5.set_ylabel("Solute Rejection (%)", fontsize=10, fontweight="bold")
    ax5.tick_params(direction="in", labelsize=10)
    ax5.legend()

    # ax6.plot(x_axis_values, lithium_sieving, linewidth=2, label="lithium")
    # ax6.plot(x_axis_values, cobalt_sieving, linewidth=2, label="cobalt")
    # ax6.set_xlabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    # ax6.set_ylabel("Solute Passage", fontsize=10, fontweight="bold")
    # ax6.tick_params(direction="in", labelsize=10)
    # ax6.legend()

    ax6.plot(x_axis_values, percent_recovery, linewidth=2)
    ax6.set_xlabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    ax6.set_ylabel("Percent Recovery (%)", fontsize=10, fontweight="bold")
    ax6.tick_params(direction="in", labelsize=10)

    plt.show()


def plot_membrane_results(m):
    """
    Plots concentrations within the membrane.

    Args:
        m: Pyomo model
    """
    x_axis_values = []
    z_axis_values = []

    for x_val in m.fs.membrane.x_bar:
        x_axis_values.append(x_val * value(m.fs.membrane.membrane_width))
    for z_val in m.fs.membrane.z_bar:
        z_axis_values.append(z_val * value(m.fs.membrane.membrane_thickness) * 1e9)
    # store values for concentration of lithium in the membrane
    conc_mem_lith = []
    conc_mem_lith_dict = {}
    # store values for concentration of cobalt in the membrane
    conc_mem_cob = []
    conc_mem_cob_dict = {}
    # store values for concentration of chlorine in the membrane
    conc_mem_chl = []
    conc_mem_chl_dict = {}

    for z_val in m.fs.membrane.z_bar:
        for x_val in m.fs.membrane.x_bar:
            conc_mem_lith.append(
                value(m.fs.membrane.membrane_conc_mass_lithium[x_val, z_val])
            )
            conc_mem_cob.append(
                value(m.fs.membrane.membrane_conc_mass_cobalt[x_val, z_val])
            )
            conc_mem_chl.append(
                value(m.fs.membrane.membrane_conc_mass_chlorine[x_val, z_val])
            )

        conc_mem_lith_dict[f"{z_val}"] = conc_mem_lith
        conc_mem_cob_dict[f"{z_val}"] = conc_mem_cob
        conc_mem_chl_dict[f"{z_val}"] = conc_mem_chl
        conc_mem_lith = []
        conc_mem_cob = []
        conc_mem_chl = []

    conc_mem_lith_df = DataFrame(index=x_axis_values, data=conc_mem_lith_dict)
    conc_mem_cob_df = DataFrame(index=x_axis_values, data=conc_mem_cob_dict)
    conc_mem_chl_df = DataFrame(index=x_axis_values, data=conc_mem_chl_dict)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=125, figsize=(15, 7))
    lithium_plot = ax1.pcolor(
        z_axis_values, x_axis_values, conc_mem_lith_df, cmap="Blues"
    )
    ax1.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Membrane Width (m)", fontsize=10, fontweight="bold")
    ax1.set_title(
        "Lithium Concentration\n in Membrane (kg/m$^3$)", fontsize=10, fontweight="bold"
    )
    ax1.tick_params(direction="in", labelsize=10)
    fig.colorbar(lithium_plot, ax=ax1)

    cobalt_plot = ax2.pcolor(
        z_axis_values, x_axis_values, conc_mem_cob_df, cmap="Blues"
    )
    ax2.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax2.set_title(
        "Cobalt Concentration\n in Membrane (kg/m$^3$)", fontsize=10, fontweight="bold"
    )
    ax2.tick_params(direction="in", labelsize=10)
    fig.colorbar(cobalt_plot, ax=ax2)

    chlorine_plot = ax3.pcolor(
        z_axis_values, x_axis_values, conc_mem_chl_df, cmap="Blues"
    )
    ax3.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax3.set_title(
        "Chlorine Concentration\n in Membrane (kg/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax3.tick_params(direction="in", labelsize=10)
    fig.colorbar(chlorine_plot, ax=ax3)

    plt.show()


if __name__ == "__main__":
    main()
