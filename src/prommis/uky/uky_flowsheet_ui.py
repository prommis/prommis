#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Interface for UKy plant model in WaterTAP UI

Authors: Dan Gunter (LBNL)
"""
__author__ = "Dan Gunter"

# stdlib
# third party
import pyomo.environ as pyo
from idaes import logger as idaeslog
from watertap.ui.fsapi import FlowsheetInterface
from watertap.ui.fsapi import FlowsheetCategory as FC

# package
from prommis.uky.uky_flowsheet import (
    build,
    set_operating_conditions,
    set_scaling,
    solve,
    initialize_system,
)

_log = idaeslog.getLogger(__name__)


def export_to_ui():
    return FlowsheetInterface(
        name="UKy",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        get_diagram=get_diagram,
        requires_idaes_solver=True,
        category=FC.wastewater,
        build_options={},
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    # wlog = logging.getLogger("idaes.watertap.ui.fsapi")
    # wlog.setLevel(logging.DEBUG)

    exports.from_csv(file="uky_flowsheet_ui.csv", flowsheet=flowsheet)

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
    comp_liq = {"H", "H2O", "HSO4", "SO4"}
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

    category = "solids"
    exports.add(
        deferred_obj=f"leach.solid_outlet.flow_mass[0]",
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
        obj_name = f"leach.solid_outlet.mass_frac_comp[0, '{c}']"
        exports.add(
            deferred_obj=obj_name,
            name=name,
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )
    exports.add(
        deferred_obj=f"leach.solid_outlet.mass_frac_comp[0, 'inerts']",
        name=f"leaching solid mass fraction of inerts",
        rounding=4,
        description=f"leaching solid mass fraction of inert components outlet",
        is_input=False,
        is_output=True,
        output_category="solid outlet",
    )

    category = "leaching"
    exports.add(
        deferred_obj="leach.liquid_outlet.flow_vol[0]",
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
        obj_name = f"leach.liquid_outlet.conc_mass_comp[0, '{c}']"
        exports.add(
            deferred_obj=obj_name,
            name=name,
            ui_units=pyo.units.mg / pyo.units.l,
            display_units="mg/l",
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )

    for stype in {"rougher", "cleaner"}:
        category = f"solex {stype}"
        # organic
        for c in comp:
            name = f"solex {stype} organic liquid mass composition fraction {c}"
            obj_name = (
                f"solex_{stype}.mscontactor.organic_outlet"
                f".conc_mass_comp[0, '{c}']"
            )
            exports.add(
                deferred_obj=obj_name,
                name=name,
                ui_units=pyo.units.mg / pyo.units.l,
                display_units="mg/l",
                rounding=4,
                description=f"{name} outlet",
                is_input=False,
                is_output=True,
                output_category=category,
            )
        # aqueous
        complist = comp.union(comp_liq) if stype == "rougher" else comp
        for c in complist:
            name = f"solex {stype} aqueous liquid mass composition fraction {c}"
            obj_name = (
                f"solex_{stype}.mscontactor.aqueous_outlet"
                f".conc_mass_comp[0, '{c}']"
            )
            exports.add(
                deferred_obj=obj_name,
                name=name,
                ui_units=pyo.units.mg / pyo.units.l,
                display_units="mg/l",
                rounding=4,
                description=f"{name} outlet",
                is_input=False,
                is_output=True,
                output_category=category,
            )

    category = "precipitator"
    exports.add(
        deferred_obj="precipitator.cv_aqueous.properties_out[0].flow_vol",
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
        obj_name = f"precipitator.cv_aqueous.properties_out[0].conc_mass_comp['{c}']"
        exports.add(
            deferred_obj=obj_name,
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
        deferred_obj="precipitator.precipitate_outlet.temperature[0]",
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
        obj_name = f"precipitator.precipitate_outlet.flow_mol_comp[0, '{c}']"
        exports.add(
            deferred_obj=obj_name,
            name=name,
            ui_units=pyo.units.mol / pyo.units.hour,
            display_units="moles/hour",
            rounding=4,
            description=f"{name} outlet",
            is_input=False,
            is_output=True,
            output_category=category,
        )
    _log.info(f"exports:\n{exports.json(indent=2)}")


def build_flowsheet(build_options=None, **kwargs):
    m = build()
    set_operating_conditions(m)
    scaled_model = set_scaling(m)
    # assert_units_consistent(scaled_model)
    # assert degrees_of_freedom(scaled_model) == 0
    initialize_system(scaled_model)
    return scaled_model


def get_diagram(build_options):
    return "uky_ui.png"


def solve_flowsheet(flowsheet=None):
    return solve(flowsheet)
