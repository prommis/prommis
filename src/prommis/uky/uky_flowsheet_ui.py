#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Interface for UKy plant model in WaterTAP UI

The `export_to_ui()` function provides all the functions that the UI
will invoke to export the inputs and outputs, build the flowsheet,
and solve the flowsheet.

Authors: Dan Gunter (LBNL), Marcus Holly (KeyLogic)
"""
__author__ = "Dan Gunter"

# third party
import pyomo.environ as pyo

from idaes import logger as idaeslog

from idaes_flowsheet_processor.api import FlowsheetCategory, FlowsheetInterface

# package
from prommis.uky.uky_flowsheet import (
    build,
    initialize_system,
    set_operating_conditions,
    set_scaling,
    solve_system,
)

_log = idaeslog.getLogger(__name__)


def export_to_ui():
    """Hook called by the UI to get the interface to the flowsheet."""
    return FlowsheetInterface(
        name="UKy",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        get_diagram=get_diagram,
        requires_idaes_solver=True,
        category=FlowsheetCategory.wastewater,
        build_options={},
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    """Export input and output variables for the UKy flowsheet."""
    _log.info(f"begin/setup-UI-exports build_options={build_options}")

    kwargs.get("", None)  # eliminate not-used warning

    # Chemical components
    comp = {
        "Al",
        "Ca",
        "Ce",
        "Dy",
        "Fe",
        "Gd",
        "La",
        "Nd",
        "Pr",
        "Sc",
        "Sm",
        "Y",
    }

    # Chemical components - oxides
    comp_ox = {
        "Al2O3",
        "CaO",
        "Ce2O3",
        "Dy2O3",
        "Fe2O3",
        "Gd2O3",
        "La2O3",
        "Nd2O3",
        "Pr2O3",
        "Sc2O3",
        "Sm2O3",
        "Y2O3",
    }

    # Organic components
    comp_org = {
        "Al_o",
        "Ca_o",
        "Ce_o",
        "Dy_o",
        "Fe_o",
        "Gd_o",
        "La_o",
        "Nd_o",
        "Pr_o",
        "Sc_o",
        "Sm_o",
        "Y_o",
    }

    # Liquid chemical components
    comp_liq = {"H", "H2O", "HSO4", "SO4"}

    # Chemical components - precipitates
    comp_precip = {
        "Al2(C2O4)3(s)",
        "Ce2(C2O4)3(s)",
        "Dy2(C2O4)3(s)",
        "Fe2(C2O4)3(s)",
        "Gd2(C2O4)3(s)",
        "La2(C2O4)3(s)",
        "Nd2(C2O4)3(s)",
        "Pr2(C2O4)3(s)",
        "Sc2(C2O4)3(s)",
        "Sm2(C2O4)3(s)",
        "Y2(C2O4)3(s)",
    }

    # Export the leach liquid feed and its mass components, as inputs
    llf = flowsheet.leach_liquid_feed
    exports.add(
        obj=llf.flow_vol[0],
        name="Leach liquid feed rate",
        ui_units=pyo.units.l / pyo.units.hour,
        display_units="L/h",
        rounding=2,
        description="Leach liquid feed volumetric flow rate",
        is_input=True,
        is_output=False,
        input_category="Leaching liquid feed",
    )
    # Mass comp
    for compound, compound_name in (
        ("H", "hydrogen"),
        ("SO4", "SO4"),
        ("HSO4", "HSO4"),
    ):
        exports.add(
            obj=llf.conc_mass_comp[0, compound],
            name=f"Leach liquid feed {compound_name}",
            description=f"Leach liquid feed {compound_name} mass composition",
            ui_units=pyo.units.mg / pyo.units.l,
            display_units="mg/L",
            rounding=3,
            is_input=True,
            is_output=False,
            input_category="Leaching liquid feed",
        )

    # Export the leaching solid feed, and its mass flow components, as inputs
    lsf = flowsheet.leach_solid_feed
    category = "Leaching solid feed"
    comp_solid_in = {"inerts"}.union(comp_ox)
    for compound in comp_solid_in:
        exports.add(
            obj=lsf.mass_frac_comp[0, compound],
            name=f"Leach solid feed {compound}",
            description=f"Leach solid feed {compound} fractional composition",
            rounding=3,
            is_input=True,
            is_output=False,
            input_category=category,
        )
    exports.add(
        obj=lsf.flow_mass[0],
        name="Leach solid feed mass flow",
        rounding=3,
        ui_units=pyo.units.kg / pyo.units.hour,
        display_units="kg/hr",
        is_input=True,
        is_output=False,
        input_category=category,
    )

    # Export the roaster inlet as an input
    category = "Roaster"
    rst = flowsheet.roaster
    exports.add(
        obj=rst.gas_inlet.temperature[0],
        name="Gas inlet temperature",
        rounding=3,
        ui_units=pyo.units.K,
        display_units="K",
        is_input=True,
        is_output=False,
        input_category=category,
    )
    exports.add(
        obj=rst.gas_inlet.pressure[0],
        name="Gas inlet temperature",
        rounding=2,
        ui_units=pyo.units.Pa,
        display_units="Pa",
        is_input=True,
        is_output=False,
        input_category=category,
    )

    # Export the roaster outlet as an output
    exports.add(
        obj=rst.gas_outlet.temperature[0],
        name="Gas outlet temperature",
        rounding=3,
        ui_units=pyo.units.K,
        display_units="K",
        is_input=False,
        is_output=True,
        input_category=category,
    )
    exports.add(
        obj=rst.gas_outlet.pressure[0],
        name="Gas pressure",
        rounding=2,
        ui_units=pyo.units.Pa,
        display_units="Pa",
        is_input=False,
        is_output=True,
        input_category=category,
    )
    exports.add(
        obj=rst.gas_out[0].flow_mol,
        name="Gas molar flow",
        ui_units=pyo.units.mol / pyo.units.s,
        display_units="mol/s",
        rounding=3,
        is_input=False,
        is_output=True,
        input_category=category,
    )

    # Export the leach solid outputs, which includes overall mass flow,
    # and mass fraction of oxides and inerts.
    category = "solids"
    leach = flowsheet.leach
    exports.add(
        obj=leach.solid_outlet.flow_mass[0],
        name=f"solid flow mass",
        rounding=4,
        ui_units=pyo.units.kg / pyo.units.hour,
        display_units="kg/hr",
        description=f"solid flow mass",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    for c in comp_ox:
        name = f"leaching solid mass fraction of {c}"
        obj = leach.solid_outlet.mass_frac_comp[0, c]
        exports.add(
            obj=obj,
            name=name,
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )
    exports.add(
        obj=leach.solid_outlet.mass_frac_comp[0, "inerts"],
        name=f"leaching solid mass fraction of inerts",
        rounding=4,
        description=f"leaching solid mass fraction of inert components outlet",
        is_input=False,
        is_output=True,
        output_category="solid outlet",
    )

    # Export leach liquid outputs, which includes the liquid flow
    # and liquid mass compositions
    category = "leaching"
    exports.add(
        obj=leach.liquid_outlet.flow_vol[0],
        name=f"liquid flow volume",
        ui_units=pyo.units.l / pyo.units.hour,
        display_units="l/h",
        rounding=4,
        description=f"liquid flow volume",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    for c in comp.union(comp_liq):
        name = f"leaching liquid mass composition fraction of {c}"
        obj_name = leach.liquid_outlet.conc_mass_comp[0, c]
        exports.add(
            obj=obj_name,
            name=name,
            ui_units=pyo.units.mg / pyo.units.l,
            display_units="mg/l",
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )

    # Export the outputs for the solex rougher and cleaner
    for stype in {
        "rougher_load",
        "rougher_scrub",
        "rougher_strip",
        "cleaner_load",
        "cleaner_strip",
    }:
        category = f"solex {stype}"
        block = getattr(flowsheet, f"solex_{stype}")
        for ltype in {"organic", "aqueous"}:
            if stype == "rougher" and ltype == "aqueous":
                # add aqueous components for the aqueous rougher
                complist = comp.union(comp_liq)
            elif ltype == "organic":
                complist = comp_org
            else:
                complist = comp
            # export the output for each component
            for c in complist:
                name = f"solex {stype} {ltype} liquid mass composition fraction {c}"
                _log.debug(f"export: {name}")
                obj = getattr(block.mscontactor, f"{ltype}_outlet").conc_mass_comp[0, c]
                exports.add(
                    obj=obj,
                    name=name,
                    ui_units=pyo.units.mg / pyo.units.l,
                    display_units="mg/l",
                    rounding=4,
                    description=f"{name} outlet",
                    is_input=False,
                    is_output=True,
                    output_category=category,
                )

    # Export the outputs for the precipitator, including overall flow
    # as well as concentration mass composition for chemical components
    # and precipitate components.
    category = "precipitator"
    precipitator = flowsheet.precipitator
    exports.add(
        obj=precipitator.cv_aqueous.properties_out[0].flow_vol,
        name=f"precipitator aqueous out",
        ui_units=pyo.units.l / pyo.units.hour,
        display_units="liters/hour",
        rounding=4,
        description=f"precipitator aqueous properties out",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    for c in comp:
        name = f"precipitator aqueous concentration mass composition of {c}"
        obj = precipitator.cv_aqueous.properties_out[0].conc_mass_comp[c]
        exports.add(
            obj=obj,
            name=name,
            ui_units=pyo.units.mg / pyo.units.l,
            display_units="mg/l",
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )
    exports.add(
        obj=precipitator.precipitate_outlet.temperature[0],
        name="precipitator outlet temperature",
        rounding=4,
        ui_units=pyo.units.K,
        display_units="K",
        description="temperature of the precipitator's precipitate outlet",
        is_input=False,
        is_output=True,
        output_category=category,
    )
    for c in comp_precip:
        name = f"precipitate molar flow composition of {c}"
        obj = precipitator.precipitate_outlet.flow_mol_comp[0, c]
        exports.add(
            obj=obj,
            name=name,
            ui_units=pyo.units.mol / pyo.units.hour,
            display_units="moles/hour",
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )
    _log.debug(f"exports:\n{exports.model_dump_json()}")
    _log.info(f"end/setup-UI-exports build_options={build_options}")


def build_flowsheet(build_options=None, **kwargs):
    """Called by the UI to build the flowsheet.
    Does not solve the flowsheet, but does set operating conditions, scaling, and
    initialize the system.
    """
    _log.info(f"begin/build-flowsheet build_options={build_options}")
    m = build()
    set_operating_conditions(m)
    set_scaling(m)
    scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)
    initialize_system(scaled_model)
    _log.info(f"end/build-flowsheet build_options={build_options}")
    return scaled_model


def get_diagram(build_options):
    """Return a diagram to be shown in the UI for this flowsheet."""
    return "uky_flowsheet_ui.png"


def solve_flowsheet(flowsheet=None):
    """Solve a built/initialized flowsheet."""

    m = build()

    set_operating_conditions(m)

    set_scaling(m)

    scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    initialize_system(scaled_model)

    results = solve_system(scaled_model)

    return results
