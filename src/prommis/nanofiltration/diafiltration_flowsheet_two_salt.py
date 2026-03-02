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
    for a two-salt (lithium chloride and cobalt chloride) solution.
    """
    # build flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # specify the feed
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"lithium": 245, "cobalt": 288, "chloride": 821},
        "diafiltrate": {"lithium": 14, "cobalt": 3, "chloride": 20},
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
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
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
    overall_results_plot = plot_results_by_length(m)
    boundary_layer_results_plot = plot_results_by_thickness(m, phase="Boundary Layer")
    membrane_results_plot = plot_results_by_thickness(m, phase="Membrane")
    rejection_plot = plot_rejection_versus_concentration(m)

    return (
        m,
        overall_results_plot,
        boundary_layer_results_plot,
        membrane_results_plot,
        rejection_plot,
    )


def update_membrane_parameters(m):
    """
    Updates parameters needed in multi-component diafiltration unit model if desired.

    Args:
        m: Pyomo model
    """
    pass


def fix_variables(m, inlet_flow_volume, inlet_concentration):
    # fix degrees of freedom in the membrane
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


def plot_results_by_length(m):
    """
    Plots concentration and flux variables across the length of the membrane module.

    Args:
        m: Pyomo model
    """
    # store values for x-coordinate
    x_axis_values = []

    # store values for concentration of lithium in the retentate
    conc_ret_lith = []
    # store values for concentration of lithium at interface of BL and M
    conc_int_lith = []
    # store values for concentration of lithium in the permeate
    conc_perm_lith = []
    # store values for concentration of cobalt in the retentate
    conc_ret_cob = []
    # store values for concentration of lithium at interface of BL and M
    conc_int_cob = []
    # store values for concentration of cobalt in the permeate
    conc_perm_cob = []

    # store values for water flux across membrane
    water_flux = []
    # store values for mol flux of lithium across membrane
    lithium_flux = []
    # store values for mol flux of cobalt across membrane
    cobalt_flux = []

    # store values for percent recovery
    percent_recovery = []

    # store values for lithium rejection (observed)
    lithium_rejection_obs = []
    # store values for lithium rejection (actual)
    lithium_rejection_act = []
    # store values for cobalt rejection (observed)
    cobalt_rejection_obs = []
    # store values for cobalt rejection (actual)
    cobalt_rejection_act = []

    for x_val in m.fs.membrane.dimensionless_module_length:
        if x_val != 0:
            x_axis_values.append(x_val * value(m.fs.membrane.total_module_length))
            conc_ret_lith.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "lithium"])
            )
            conc_int_lith.append(
                value(
                    m.fs.membrane.boundary_layer_conc_mol_comp[0, x_val, 1, "lithium"]
                )
            )
            conc_perm_lith.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
            )
            conc_ret_cob.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "cobalt"])
            )
            conc_int_cob.append(
                value(m.fs.membrane.boundary_layer_conc_mol_comp[0, x_val, 1, "cobalt"])
            )
            conc_perm_cob.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
            )

            water_flux.append(value(m.fs.membrane.volume_flux_water[0, x_val]))
            lithium_flux.append(
                value(m.fs.membrane.molar_ion_flux[0, x_val, "lithium"])
            )
            cobalt_flux.append(value(m.fs.membrane.molar_ion_flux[0, x_val, "cobalt"]))

            lithium_rejection_obs.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
                        / value(
                            m.fs.membrane.retentate_conc_mol_comp[0, x_val, "lithium"]
                        )
                    )
                )
                * 100
            )
            lithium_rejection_act.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
                        / value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, 1, "lithium"
                            ]
                        )
                    )
                )
                * 100
            )
            cobalt_rejection_obs.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
                        / value(
                            m.fs.membrane.retentate_conc_mol_comp[0, x_val, "cobalt"]
                        )
                    )
                )
                * 100
            )
            cobalt_rejection_act.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
                        / value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, 1, "cobalt"
                            ]
                        )
                    )
                )
                * 100
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
    ax1.plot(x_axis_values, conc_int_lith, linewidth=2, label="interface")
    ax1.plot(x_axis_values, conc_perm_lith, linewidth=2, label="permeate")
    ax1.set_ylabel(
        "Lithium Concentration (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.tick_params(direction="in", labelsize=10)
    ax1.legend()

    ax2.plot(x_axis_values, conc_ret_cob, linewidth=2, label="retentate")
    ax2.plot(x_axis_values, conc_int_cob, linewidth=2, label="interface")
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

    ax4.plot(x_axis_values, lithium_flux, linewidth=2, label="lithium")
    ax4.plot(x_axis_values, cobalt_flux, linewidth=2, label="cobalt")
    ax4.set_xlabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax4.set_ylabel("Solute Molar Flux (mol/m$^2$/h)", fontsize=10, fontweight="bold")
    ax4.tick_params(direction="in", labelsize=10)
    ax4.legend()

    ax5.plot(
        x_axis_values, lithium_rejection_obs, linewidth=2, label="lithium (observed)"
    )
    ax5.plot(
        x_axis_values,
        lithium_rejection_act,
        "--",
        linewidth=2,
        label="lithium (actual)",
    )
    ax5.plot(
        x_axis_values, cobalt_rejection_obs, linewidth=2, label="cobalt (observed)"
    )
    ax5.plot(
        x_axis_values, cobalt_rejection_act, "--", linewidth=2, label="cobalt (actual)"
    )
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


def plot_results_by_thickness(m, phase):
    """
    Plots concentrations within the boundary layer or membrane.

    Args:
        m: Pyomo model
    """
    x_axis_values = []
    z_axis_values = []

    # store values for concentration of lithium
    conc_lith = []
    conc_lith_dict = {}
    # store values for concentration of cobalt
    conc_cob = []
    conc_cob_dict = {}
    # store values for concentration of chloride
    conc_chl = []
    conc_chl_dict = {}

    for x_val in m.fs.membrane.dimensionless_module_length:
        if x_val != 0:
            x_axis_values.append(x_val * value(m.fs.membrane.total_module_length))

    if phase == "Boundary Layer":
        for z_val in m.fs.membrane.dimensionless_boundary_layer_thickness:
            z_axis_values.append(
                z_val * value(m.fs.membrane.total_boundary_layer_thickness) * 1e6
            )
            for x_val in m.fs.membrane.dimensionless_module_length:
                if x_val != 0:
                    conc_lith.append(
                        value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, z_val, "lithium"
                            ]
                        )
                    )
                    conc_cob.append(
                        value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, z_val, "cobalt"
                            ]
                        )
                    )
                    conc_chl.append(
                        value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, z_val, "chloride"
                            ]
                        )
                    )

            conc_lith_dict[f"{z_val}"] = conc_lith
            conc_cob_dict[f"{z_val}"] = conc_cob
            conc_chl_dict[f"{z_val}"] = conc_chl
            conc_lith = []
            conc_cob = []
            conc_chl = []

    elif phase == "Membrane":
        for z_val in m.fs.membrane.dimensionless_membrane_thickness:
            z_axis_values.append(
                z_val * value(m.fs.membrane.total_membrane_thickness) * 1e9
            )
            for x_val in m.fs.membrane.dimensionless_module_length:
                if x_val != 0:
                    conc_lith.append(
                        value(
                            m.fs.membrane.membrane_conc_mol_comp[
                                0, x_val, z_val, "lithium"
                            ]
                        )
                    )
                    conc_cob.append(
                        value(
                            m.fs.membrane.membrane_conc_mol_comp[
                                0, x_val, z_val, "cobalt"
                            ]
                        )
                    )
                    conc_chl.append(
                        value(
                            m.fs.membrane.membrane_conc_mol_comp[
                                0, x_val, z_val, "chloride"
                            ]
                        )
                    )

            conc_lith_dict[f"{z_val}"] = conc_lith
            conc_cob_dict[f"{z_val}"] = conc_cob
            conc_chl_dict[f"{z_val}"] = conc_chl
            conc_lith = []
            conc_cob = []
            conc_chl = []

    conc_lith_df = DataFrame(index=x_axis_values, data=conc_lith_dict)
    conc_cob_df = DataFrame(index=x_axis_values, data=conc_cob_dict)
    conc_chl_df = DataFrame(index=x_axis_values, data=conc_chl_dict)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=125, figsize=(15, 7))
    lithium_plot = ax1.pcolor(z_axis_values, x_axis_values, conc_lith_df, cmap="Greens")
    if phase == "Boundary Layer":
        ax1.set_xlabel("Boundary Layer Thickness (um)", fontsize=10, fontweight="bold")
    elif phase == "Membrane":
        ax1.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Module Length (m)", fontsize=10, fontweight="bold")
    ax1.set_title(
        f"Lithium Concentration\n in {phase} (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.tick_params(direction="in", labelsize=10)
    fig.colorbar(lithium_plot, ax=ax1)

    cobalt_plot = ax2.pcolor(z_axis_values, x_axis_values, conc_cob_df, cmap="Blues")
    if phase == "Boundary Layer":
        ax2.set_xlabel("Boundary Layer Thickness (um)", fontsize=10, fontweight="bold")
    elif phase == "Membrane":
        ax2.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax2.set_title(
        f"Cobalt Concentration\n in {phase} (mol/m$^3$)", fontsize=10, fontweight="bold"
    )
    ax2.tick_params(direction="in", labelsize=10)
    fig.colorbar(cobalt_plot, ax=ax2)

    chloride_plot = ax3.pcolor(
        z_axis_values, x_axis_values, conc_chl_df, cmap="Oranges"
    )
    if phase == "Boundary Layer":
        ax3.set_xlabel("Boundary Layer Thickness (um)", fontsize=10, fontweight="bold")
    elif phase == "Membrane":
        ax3.set_xlabel("Membrane Thickness (nm)", fontsize=10, fontweight="bold")
    ax3.set_title(
        f"Chloride Concentration\n in {phase} (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax3.tick_params(direction="in", labelsize=10)
    fig.colorbar(chloride_plot, ax=ax3)

    plt.show()

    return fig


def plot_rejection_versus_concentration(m):
    """
    Plots rejection versus retentate-side concentration.

    Args:
        m: Pyomo model
    """
    # store values for concentration of lithium in the retentate
    conc_ret_lith = []
    # store values for concentration of lithium at interface of BL and M
    conc_int_lith = []
    # store values for concentration of lithium in the permeate
    conc_perm_lith = []
    # store values for concentration of cobalt in the retentate
    conc_ret_cob = []
    # store values for concentration of lithium at interface of BL and M
    conc_int_cob = []
    # store values for concentration of cobalt in the permeate
    conc_perm_cob = []

    # store values for lithium rejection (observed)
    lithium_rejection_obs = []
    # store values for lithium rejection (actual)
    lithium_rejection_act = []
    # store values for cobalt rejection (observed)
    cobalt_rejection_obs = []
    # store values for cobalt rejection (actual)
    cobalt_rejection_act = []

    for x_val in m.fs.membrane.dimensionless_module_length:
        if x_val != 0:
            conc_ret_lith.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "lithium"])
            )
            conc_int_lith.append(
                value(
                    m.fs.membrane.boundary_layer_conc_mol_comp[0, x_val, 1, "lithium"]
                )
            )
            conc_perm_lith.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
            )
            conc_ret_cob.append(
                value(m.fs.membrane.retentate_conc_mol_comp[0, x_val, "cobalt"])
            )
            conc_int_cob.append(
                value(m.fs.membrane.boundary_layer_conc_mol_comp[0, x_val, 1, "cobalt"])
            )
            conc_perm_cob.append(
                value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
            )

            lithium_rejection_obs.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
                        / value(
                            m.fs.membrane.retentate_conc_mol_comp[0, x_val, "lithium"]
                        )
                    )
                )
                * 100
            )
            lithium_rejection_act.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "lithium"])
                        / value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, 1, "lithium"
                            ]
                        )
                    )
                )
                * 100
            )
            cobalt_rejection_obs.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
                        / value(
                            m.fs.membrane.retentate_conc_mol_comp[0, x_val, "cobalt"]
                        )
                    )
                )
                * 100
            )
            cobalt_rejection_act.append(
                (
                    1
                    - (
                        value(m.fs.membrane.permeate_conc_mol_comp[0, x_val, "cobalt"])
                        / value(
                            m.fs.membrane.boundary_layer_conc_mol_comp[
                                0, x_val, 1, "cobalt"
                            ]
                        )
                    )
                )
                * 100
            )

    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=100, figsize=(10, 5))

    ax1.plot(conc_ret_lith, lithium_rejection_obs, linewidth=2, label="observed")
    ax1.plot(conc_ret_lith, lithium_rejection_act, linewidth=2, label="actual")
    ax1.set_xlabel(
        "Lithium Concentration (Feed-Side) (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax1.set_ylabel("Percent Rejection (%)", fontsize=10, fontweight="bold")
    ax1.tick_params(direction="in", labelsize=10)
    ax1.legend()

    ax2.plot(conc_ret_cob, cobalt_rejection_obs, linewidth=2, label="observed")
    ax2.plot(conc_ret_cob, cobalt_rejection_act, linewidth=2, label="actual")
    ax2.set_xlabel(
        "Cobalt Concentration (Feed-Side) (mol/m$^3$)",
        fontsize=10,
        fontweight="bold",
    )
    ax2.set_ylabel("Percent Rejection (%)", fontsize=10, fontweight="bold")
    ax2.tick_params(direction="in", labelsize=10)
    ax2.legend()

    plt.show()

    return fig


if __name__ == "__main__":
    main()
