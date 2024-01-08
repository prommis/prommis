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
from watertap.ui.fsapi import FlowsheetInterface  # pylint: disable=import-error

# package
from prommis.uky.uky_flowsheet import (
    build,
    set_operating_conditions,
    set_scaling,
    solve,
)


def export_to_ui():
    return FlowsheetInterface(
        name="UKy",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
        get_diagram=get_diagram,
        requires_idaes_solver=True,
        category="PrOMMiS",
        build_options={},
    )


def export_variables(flowsheet=None, exports=None, build_options=None, **kwargs):
    exports.from_csv(file="uky_flowsheet_ui.csv", flowsheet=flowsheet)


def build_flowsheet(build_options=None, **kwargs):
    m = build()
    set_operating_conditions(m)
    scaled_model = set_scaling(m)
    # initialize_system(scaled_model) XXX: Currently not working
    return m


def get_diagram(build_options):
    return "uky_ui.png"


def solve_flowsheet(flowsheet=None):
    return solve(flowsheet)  # XXX: also not working


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
