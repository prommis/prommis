#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Sample flowsheet for the multi-component diafiltration cascade.

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

from prommis.nanofiltration.multi_component_diafiltration_stream_properties import (
    MultiComponentDiafiltrationStreamParameter,
)
from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)
from prommis.nanofiltration.multi_component_diafiltration import (
    MultiComponentDiafiltration,
)


def main():
    """
    Builds and solves flowsheet with multi-component diafiltration unit model
    for a two-salt (LiCl + CoCl2) solution.
    """
    # build flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # specify the feed
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Cl": 821},
        "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
    }

    m.fs.stream_properties = MultiComponentDiafiltrationStreamParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    # add feed blocks for feed and diafiltrate
    m.fs.feed_block = Feed(property_package=m.fs.stream_properties)
    m.fs.diafiltrate_block = Feed(property_package=m.fs.stream_properties)

    # add the membrane unit model
    m.fs.membrane = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    # update parameter inputs if desired
    update_membrane_parameters(m)

    # add product blocks for retentate and permeate
    m.fs.retentate_block = Product(property_package=m.fs.stream_properties)
    m.fs.permeate_block = Product(property_package=m.fs.stream_properties)

    # fix the degrees of freedom
    fix_variables(m, inlet_flow_volume, inlet_concentration)

    # initialize membrane model
    initialized_membrane_model = m.fs.membrane.default_initializer()
    initialized_membrane_model.initialize(m.fs.membrane)

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
    overall_results_plot = plot_results(m)
    membrane_results_plot = plot_membrane_results(m)

    print(overall_results_plot)

    return (m, overall_results_plot, membrane_results_plot)


def update_membrane_parameters(m):
    """
    Updates parameters needed in multi-component diafiltration unit model if desired.

    Args:
        m: Pyomo model
    """
    pass


def fix_variables(m, inlet_flow_volume, inlet_concentration):
    # fix the nine degrees of freedom in the membrane
    m.fs.membrane.total_module_length.fix()
    m.fs.membrane.total_membrane_length.fix()
    m.fs.membrane.applied_pressure.fix()

    m.fs.membrane.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.membrane.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.membrane.time:
        for j in m.fs.membrane.solutes:
            m.fs.membrane.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.membrane.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )


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
    Plots concentration and flux variables across the length of the membrane module.

    Args:
        m: Pyomo model
    """
    # store values for x-coordinate
    x_axis_values = []

    # store values for concentration of Li in the retentate
    conc_ret_lith = []
    # store values for concentration of Li in the permeate
    conc_perm_lith = []
    # store values for concentration of Co in the retentate
    conc_ret_cob = []
    # store values for concentration of Co in the permeate
    conc_perm_cob = []

    # store values for water flux across membrane
    water_flux = []
    # store values for mol flux of Li across membrane
    Li_flux = []

    # store values for percent recovery
    percent_recovery = []

    # store values for Li rejection
    Li_rejection = []
    # store values for Li solute passage
    Li_sieving = []
    # store values for Co rejection
    Co_rejection = []
    # store values for Co solute passage
    Co_sieving = []

    for x_val in m.fs.membrane.dimensionless_module_length:
        if x_val != 0:
            x_axis_values.append(x_val * value(m.fs.membrane.total_module_length))
            conc_ret_lith.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Li"])
            )
            conc_perm_lith.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Li"])
            )
            conc_ret_cob.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Co"])
            )
            conc_perm_cob.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Co"])
            )

            water_flux.append(value(m.fs.membrane.volume_flux_water[0, x_val]))
            Li_flux.append(value(m.fs.membrane.molar_ion_flux[0, x_val, "Li"]))

            Li_rejection.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Li"])
                        / value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Li"])
                    )
                )
                * 100
            )
            Li_sieving.append(
                (
                    value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Li"])
                    / value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Li"])
                )
            )
            Co_rejection.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Co"])
                        / value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Co"])
                    )
                )
                * 100
            )
            Co_sieving.append(
                (
                    value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "Co"])
                    / value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "Co"])
                )
            )

            percent_recovery.append(
                (
                    value(m.fs.membrane.permeate_flow_volume[0, x_val])
                    / (
                        value(m.fs.membrane.feed_flow_volume[0])
                        + value(m.fs.membrane.diafiltrate_flow_volume[0])
                    )
                    * 100
                )
            )

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        3, 2, dpi=100, figsize=(12, 10)
    )

    ax1.plot(x_axis_values, conc_ret_lith, linewidth=2, label="retentate")
    ax1.plot(x_axis_values, conc_perm_lith, linewidth=2, label="permeate")
    ax1.set_ylabel(
        "Lithium Concentration (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.tick_params(direction="in", labelsize=10)
    ax1.legend()

    ax2.plot(x_axis_values, conc_ret_cob, linewidth=2, label="retentate")
    ax2.plot(x_axis_values, conc_perm_cob, linewidth=2, label="permeate")
    ax2.set_ylabel(
        "Cobalt Concentration (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax2.tick_params(direction="in", labelsize=10)
    ax2.legend()

    ax3.plot(x_axis_values, water_flux, linewidth=2)
    ax3.set_xlabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax3.set_ylabel("Water Flux (m$^3$/m$^2$/h)", fontsize=10, fontweight="bold")
    ax3.tick_params(direction="in", labelsize=10)

    ax4.plot(x_axis_values, Li_flux, linewidth=2)
    ax4.set_xlabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax4.set_ylabel("Lithium Molar Flux (mol/m$^2$/h)", fontsize=10, fontweight="bold")
    ax4.tick_params(direction="in", labelsize=10)

    ax5.plot(x_axis_values, Li_rejection, linewidth=2, label="Li")
    ax5.plot(x_axis_values, Co_rejection, linewidth=2, label="Co")
    ax5.set_xlabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax5.set_ylabel("Solute Rejection (%)", fontsize=10, fontweight="bold")
    ax5.tick_params(direction="in", labelsize=10)
    ax5.legend()

    ax6.plot(x_axis_values, percent_recovery, linewidth=2)
    ax6.set_xlabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax6.set_ylabel("Percent Recovery (%)", fontsize=10, fontweight="bold")
    ax6.tick_params(direction="in", labelsize=10)

    plt.show()

    return fig


def plot_membrane_results(m):
    """
    Plots concentrations within the membrane.

    Args:
        m: Pyomo model
    """
    x_axis_values = []
    z_axis_values = []

    for x_val in m.fs.membrane.dimensionless_module_length:
        if x_val != 0:
            x_axis_values.append(x_val * value(m.fs.membrane.total_module_length))
    for z_val in m.fs.membrane.dimensionless_membrane_thickness:
        z_axis_values.append(
            z_val * value(m.fs.membrane.total_membrane_thickness) * 1e9
        )
    # store values for concentration of Li in the membrane
    conc_mem_lith = []
    conc_mem_lith_dict = {}
    # store values for concentration of Co in the membrane
    conc_mem_cob = []
    conc_mem_cob_dict = {}
    # store values for concentration of Cl in the membrane
    conc_mem_chl = []
    conc_mem_chl_dict = {}

    for z_val in m.fs.membrane.dimensionless_membrane_thickness:
        for x_val in m.fs.membrane.dimensionless_module_length:
            if x_val != 0:
                conc_mem_lith.append(
                    value(m.fs.membrane.membrane_conc_mol_comp[0, x_val, z_val, "Li"])
                )
                conc_mem_cob.append(
                    value(m.fs.membrane.membrane_conc_mol_comp[0, x_val, z_val, "Co"])
                )
                conc_mem_chl.append(
                    value(m.fs.membrane.membrane_conc_mol_comp[0, x_val, z_val, "Cl"])
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
    Li_plot = ax1.pcolor(z_axis_values, x_axis_values, conc_mem_lith_df, cmap="Greens")
    ax1.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax1.set_title(
        "Lithium Concentration\n in Membrane (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.tick_params(direction="in", labelsize=10)
    fig.colorbar(Li_plot, ax=ax1)

    Co_plot = ax2.pcolor(z_axis_values, x_axis_values, conc_mem_cob_df, cmap="Blues")
    ax2.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax2.set_title(
        "Cobalt Concentration\n in Membrane (mol/m$^3$)", fontsize=10, fontweight="bold"
    )
    ax2.tick_params(direction="in", labelsize=10)
    fig.colorbar(Co_plot, ax=ax2)

    Cl_plot = ax3.pcolor(z_axis_values, x_axis_values, conc_mem_chl_df, cmap="Oranges")
    ax3.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax3.set_title(
        "Chloride Concentration\n in Membrane (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax3.tick_params(direction="in", labelsize=10)
    fig.colorbar(Cl_plot, ax=ax3)

    plt.show()

    return fig


if __name__ == "__main__":
    main()
