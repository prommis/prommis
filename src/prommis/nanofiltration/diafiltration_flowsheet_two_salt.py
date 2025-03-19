from pyomo.environ import (
    ConcreteModel,
    Param,
    Set,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
    value,
)

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import matplotlib.pyplot as plt
import seaborn as sns
from pandas import DataFrame


from diafiltration_two_salt import TwoSaltDiafiltration
from diafiltration_solute_properties import SoluteParameter


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteParameter()

    # parameter inputs
    build_membrane_parameters(m)

    m.fs.membrane = TwoSaltDiafiltration(
        property_package=m.fs.properties,
        membrane_length=m.membrane_length,
        membrane_width=m.membrane_width,
        membrane_thickness=m.membrane_thickness,
        membrane_permeability=m.Lp,
        applied_pressure=m.dP,
        feed_flow_volume=m.feed_flow_volume,
        feed_conc_mass_lithium=m.feed_conc_mass_lithium,
        feed_conc_mass_cobalt=m.feed_conc_mass_cobalt,
        diafiltrate_flow_volume=m.diafiltrate_flow_volume,
        diafiltrate_conc_mass_lithium=m.diafiltrate_conc_mass_lithium,
        diafiltrate_conc_mass_cobalt=m.diafiltrate_conc_mass_cobalt,
        NFEx=5,
        NFEz=5,
    )

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    solve_model(m)  # TODO: debug numerical scaling
    dt.report_numerical_issues()

    plot_results(m)
    plot_membrane_results(m)


def build_membrane_parameters(m):

    m.membrane_thickness = Param(
        initialize=1e-7,
        units=units.m,
        doc="Thickness of membrane (z-direction)",
    )
    m.membrane_width = Param(
        initialize=1,
        units=units.m,
        doc="Width of the membrane (x-direction)",
    )
    m.membrane_length = Param(
        initialize=100,
        units=units.m,
        doc="Length of the membrane, wound radially",
    )
    m.dP = Param(
        initialize=10,
        units=units.bar,
        doc="Pressure applied to membrane",
    )
    m.Lp = Param(
        initialize=0.01,
        units=units.m / units.h / units.bar,
        doc="Hydraulic permeability coefficient",
    )
    m.feed_flow_volume = Param(
        initialize=100,
        units=units.m**3 / units.h,
        doc="Volumetric flow rate of the feed",
    )
    m.feed_conc_mass_lithium = Param(
        initialize=1.7,
        units=units.kg / units.m**3,
        doc="Mass concentration of lithium in the feed",
    )
    m.feed_conc_mass_cobalt = Param(
        initialize=17,
        units=units.kg / units.m**3,
        doc="Mass concentration of cobalt in the feed",
    )
    m.diafiltrate_flow_volume = Param(
        initialize=30,
        units=units.m**3 / units.h,
        doc="Volumetric flow rate of the diafiltrate",
    )
    m.diafiltrate_conc_mass_lithium = Param(
        initialize=0.1,
        units=units.kg / units.m**3,
        doc="Mass concentration of lithium in the diafiltrate",
    )
    m.diafiltrate_conc_mass_cobalt = Param(
        initialize=0.2,
        units=units.kg / units.m**3,
        doc="Mass concentration of cobalt in the diafiltrate",
    )


def solve_model(m):
    set_scaling(m)
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    solver = SolverFactory("ipopt")
    solver.solve(scaled_model, tee=True)

    scaling.propagate_solution(scaled_model, m)


def set_scaling(m):
    """
    Apply scaling factors to certain constraints to improve solver performance

    Args:
        m: Pyomo model
    """
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Add scaling factors for poorly scaled variables
    for x in m.fs.membrane.x_bar:
        m.scaling_factor[m.fs.membrane.retentate_flow_volume[x]] = 1e-2
        m.scaling_factor[m.fs.membrane.retentate_conc_mass_cobalt[x]] = 1e-1
        m.scaling_factor[m.fs.membrane.retentate_conc_mass_chlorine[x]] = 1e-1
        for z in m.fs.membrane.z_bar:
            m.scaling_factor[m.fs.membrane.D_lithium_lithium[x, z]] = 1e8
            m.scaling_factor[m.fs.membrane.D_lithium_cobalt[x, z]] = 1e8
            m.scaling_factor[m.fs.membrane.D_cobalt_lithium[x, z]] = 1e8
            m.scaling_factor[m.fs.membrane.D_cobalt_cobalt[x, z]] = 1e8

            m.scaling_factor[m.fs.membrane.volume_flux_water[x]] = 1e2
            m.scaling_factor[m.fs.membrane.mass_flux_lithium[x]] = 1e2
            m.scaling_factor[m.fs.membrane.mass_flux_cobalt[x]] = 1e2
            m.scaling_factor[m.fs.membrane.mass_flux_chlorine[x]] = 1e2

    # Add scaling factors for poorly scaled constraints
    for x in m.fs.membrane.x_bar:
        for z in m.fs.membrane.z_bar:
            m.scaling_factor[m.fs.membrane.D_lithium_lithium_calculation[x, z]] = 1e12
            m.scaling_factor[m.fs.membrane.D_lithium_cobalt_calculation[x, z]] = 1e12
            m.scaling_factor[m.fs.membrane.D_cobalt_lithium_calculation[x, z]] = 1e12
            m.scaling_factor[m.fs.membrane.D_cobalt_cobalt_calculation[x, z]] = 1e12

            if z != 0:
                m.scaling_factor[m.fs.membrane.lithium_flux_membrane[x, z]] = 1e8
                m.scaling_factor[m.fs.membrane.cobalt_flux_membrane[x, z]] = 1e8


def plot_results(m):
    x_plot = []
    conc_ret_lith = []
    conc_perm_lith = []
    conc_ret_cob = []
    conc_perm_cob = []

    water_flux = []
    lithium_flux = []

    for x_val in m.fs.membrane.x_bar:
        x_plot.append(x_val * value(m.fs.membrane.config.membrane_width))
        conc_ret_lith.append(value(m.fs.membrane.retentate_conc_mass_lithium[x_val]))
        conc_perm_lith.append(value(m.fs.membrane.permeate_conc_mass_lithium[x_val]))
        conc_ret_cob.append(value(m.fs.membrane.retentate_conc_mass_cobalt[x_val]))
        conc_perm_cob.append(value(m.fs.membrane.permeate_conc_mass_cobalt[x_val]))

        water_flux.append(value(m.fs.membrane.volume_flux_water[x_val]))
        lithium_flux.append(value(m.fs.membrane.mass_flux_lithium[x_val]))

    fig1, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        3, 2, dpi=125, figsize=(10, 7)
    )

    ax1.plot(x_plot, conc_ret_lith, linewidth=2)
    ax1.set_ylim(1.3, 1.4)
    ax1.set_ylabel(
        "Retentate-side Lithium\n Concentration (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax1.tick_params(direction="in", labelsize=10)

    ax2.plot(x_plot, conc_perm_lith, linewidth=2)
    ax2.set_ylabel(
        "Permeate-side Lithium\n Concentration (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax2.tick_params(direction="in", labelsize=10)

    ax3.plot(x_plot, conc_ret_cob, linewidth=2)
    ax3.set_ylim(13, 13.5)
    ax3.set_ylabel(
        "Retentate-side Cobalt\n Concentration (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax3.tick_params(direction="in", labelsize=10)

    ax4.plot(x_plot, conc_perm_cob, linewidth=2)
    ax4.set_ylabel(
        "Permeate-side Cobalt\n Concentration (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax4.tick_params(direction="in", labelsize=10)

    ax5.plot(x_plot, water_flux, linewidth=2)
    ax5.set_xlabel("Membrane Length (m)", fontsize=10, fontweight="bold")
    ax5.set_ylabel("Water Flux (m3/m2/h)", fontsize=10, fontweight="bold")
    ax5.tick_params(direction="in", labelsize=10)

    ax6.plot(x_plot, lithium_flux, linewidth=2)
    ax6.set_xlabel("Membrane Length (m)", fontsize=10, fontweight="bold")
    ax6.set_ylabel("Mass Flux of Lithium\n (kg/m2/h)", fontsize=10, fontweight="bold")
    ax6.tick_params(direction="in", labelsize=10)

    plt.show()


def plot_membrane_results(m):
    x_vals = []
    z_vals = []

    for x_val in m.fs.membrane.x_bar:
        x_vals.append(x_val)
    for z_val in m.fs.membrane.z_bar:
        z_vals.append(z_val)

    c_lith_mem = []
    c_cob_mem = []
    c_chl_mem = []

    c_lith_mem_dict = {}
    c_cob_mem_dict = {}
    c_chl_mem_dict = {}

    for z_val in m.fs.membrane.z_bar:
        for x_val in m.fs.membrane.x_bar:
            c_lith_mem.append(
                value(m.fs.membrane.membrane_conc_mass_lithium[x_val, z_val])
            )
            c_cob_mem.append(
                value(m.fs.membrane.membrane_conc_mass_cobalt[x_val, z_val])
            )
            c_chl_mem.append(
                value(m.fs.membrane.membrane_conc_mass_chlorine[x_val, z_val])
            )

        c_lith_mem_dict[f"{z_val}"] = c_lith_mem
        c_cob_mem_dict[f"{z_val}"] = c_cob_mem
        c_chl_mem_dict[f"{z_val}"] = c_chl_mem
        c_lith_mem = []
        c_cob_mem = []
        c_chl_mem = []

    c_lith_mem_df = DataFrame(index=x_vals, data=c_lith_mem_dict)
    c_cob_mem_df = DataFrame(index=x_vals, data=c_cob_mem_dict)
    c_chl_mem_df = DataFrame(index=x_vals, data=c_chl_mem_dict)

    figs, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=125, figsize=(15, 7))
    sns.heatmap(
        ax=ax1,
        data=c_lith_mem_df,
        cmap="mako",
    )
    ax1.tick_params(axis="x", labelrotation=45)
    ax1.set_xlabel("z (dimensionless)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("x (dimensionless)", fontsize=10, fontweight="bold")
    ax1.invert_yaxis()
    ax1.set_title(
        "Lithium Concentration\n in Membrane (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax1.tick_params(direction="in", labelsize=10)

    sns.heatmap(
        ax=ax2,
        data=c_cob_mem_df,
        cmap="mako",
    )
    ax2.tick_params(axis="x", labelrotation=45)
    ax2.set_xlabel("z (dimensionless)", fontsize=10, fontweight="bold")
    # ax2.set_ylabel("x (dimensionless)", fontsize=10, fontweight="bold")
    ax2.invert_yaxis()
    ax2.set_title(
        "Cobalt Concentration\n in Membrane (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax2.tick_params(direction="in", labelsize=10)

    sns.heatmap(
        ax=ax3,
        data=c_chl_mem_df,
        cmap="mako",
    )
    ax3.tick_params(axis="x", labelrotation=45)
    ax3.set_xlabel("z (dimensionless)", fontsize=10, fontweight="bold")
    # ax3.set_ylabel("x (dimensionless)", fontsize=10, fontweight="bold")
    ax3.invert_yaxis()
    ax3.set_title(
        "Chlorine Concentration\n in Membrane (kg/m3)", fontsize=10, fontweight="bold"
    )
    ax3.tick_params(direction="in", labelsize=10)

    plt.show()


if __name__ == "__main__":
    main()
