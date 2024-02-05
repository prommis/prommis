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
import logging
# third party
from idaes import logger as idaeslog
from watertap.ui.fsapi import FlowsheetInterface
from watertap.ui.fsapi import FlowsheetCategory as FC

# package
from prommis.uky.uky_flowsheet import (
    build,
    set_operating_conditions,
    set_scaling,
    solve,
    initialize_system
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
    wlog = logging.getLogger("idaes.watertap.ui.fsapi")
    wlog.setLevel(logging.DEBUG)

    exports.from_csv(file="uky_flowsheet_ui.csv", flowsheet=flowsheet)
    fs = flowsheet

    # TODO: Fix compositions to match test_uky_flowsheet exactly

    compositions = {
        "liquid1": ("'Al', 'Ca', 'Ce', 'Dy', 'Fe', 'Gd', 'H', 'H2O', 'HSO4', 'SO4', 'La', 'Nd', 'Pr',
        'Sc', 'Sm', 'Y'),
        "solex"
    }
    units = {
        "conc_mass_comp": {
            "outlets": [
                ("leach liquid", "fs.leach.liquid_outlet", liquid1),
            "solex rougher mscontactor organic":
                "fs.solex_rougher.mscontactor.organic_outlet",
            "solex rougher mscontactor aqueous":
                "fs.solex_rougher.mscontactor.aqueous_outlet",
            "solex cleaner mscontactor organic":
                "fs.solex_cleaner.mscontactor.organic_outlet",
            "solex cleaner mscontactor aqueous":
                "fs.solex_cleaner.mscontactor.aqueous_outlet",
            "precipitator cv aqueous":
                "fs.precipitator.cv_aqueous.properties_out[0]",
            ],
            "value_name": "mass composition concentration"
        },
        "flow_mol_comp": {
            "value_name": "flow molar composition",
            "compositions": ("Al2(C2O4)3(s)", "Dy2(C2O4)3(s)", "Fe2(C2O4)3(s)",
                             "Gd2(C2O4)3(s)", "La2(C2O4)3(s)", "Nd2(C2O4)3(s)",
                             "Pr2(C2O4)3(s)", "Sc2(C2O4)3(s)", "Sm2(C2O4)3(s)",
                             "Y2(C2O4)3(s)"),
            "outlets": {
                "precipitator precipitate": fs.precipitator.precipitate_outlet
            }
        }
    }
    # add exports for: {type of value} X {unit} X {composition}
    for value_type, data in units.items():
        for name, outlet in data["outlets"].items():
            for comp in data["compositions"]:
                obj_name = f"{outlet}.{value_type}[0,'{comp}']"
                _log.debug(f"start: add UI export for obj={obj_name}")
                exports.add(
                    deferred_obj=obj_name,
                    name=name,
                    rounding=2,
                    description=f"{name} outlet",
                    is_input=False,
                    is_output=True,
                    output_category=data["value_name"]
                )
                _log.info(f"end: add UI export for obj={obj_name}")


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


# for terminal debugging
if __name__ == "__main__":
    import argparse
    from idaes import logger
    import logging
    import sys

    p = argparse.ArgumentParser()
    p.add_argument("-v", "--verbose", action="count", default=0)
    a = p.parse_args()

    if a.verbose > 0:
        log = logger.getLogger("watertap.ui.fsapi")
        if a.verbose > 1:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.INFO)

    iface = export_to_ui()
    iface.build()

    iface.fs_exp.to_csv(output=sys.stdout)
