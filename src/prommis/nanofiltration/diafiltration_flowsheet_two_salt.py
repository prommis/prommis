from pyomo.environ import (
    ConcreteModel,
    Param,
    Set,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
)

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

from diafiltration_two_salt import TwoSaltDiafiltration
from diafiltration_solute_properties import SoluteParameter


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # m.fs.solutes = Set(["lithium", "cobalt", "chlorine"])
    m.fs.properties = (
        SoluteParameter()
    )  # TODO: call relevent parameters from property package

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
        NFEx=10,
        NFEz=8,
    )

    # m.fs.membrane.display()

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    solve_model(m)  # TODO: debug numerical scaling
    dt.report_numerical_issues()


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
        for z in m.fs.membrane.z_bar:
            m.scaling_factor[m.fs.membrane.D_lithium_lithium[x, z]] = 1e6
            m.scaling_factor[m.fs.membrane.D_lithium_cobalt[x, z]] = 1e6
            m.scaling_factor[m.fs.membrane.D_cobalt_lithium[x, z]] = 1e6
            m.scaling_factor[m.fs.membrane.D_cobalt_cobalt[x, z]] = 1e6

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
                m.scaling_factor[m.fs.membrane.lithium_flux_membrane[x, z]] = 1e2
                m.scaling_factor[m.fs.membrane.cobalt_flux_membrane[x, z]] = 1e2


if __name__ == "__main__":
    main()
