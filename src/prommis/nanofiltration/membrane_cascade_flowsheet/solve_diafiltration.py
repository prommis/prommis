#!/usr/bin/env python
#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Executable file for generating and solving diafiltration model."""

import sys

from pyomo.environ import SolverFactory

from idaes.core.util.model_statistics import report_statistics

from prommis.nanofiltration.membrane_cascade_flowsheet import utils
from prommis.nanofiltration.membrane_cascade_flowsheet.diafiltration_flowsheet_model import (
    DiafiltrationModel,
)


def main(args):
    """Driver for creating diafiltration model."""
    # collect arguments
    # check if arguments are given. Use default if not
    if len(args) == 1:
        print("No args provided")
        print('Using default "stage" 3 10')
        mix_style = "stage"
        num_s = 3
        num_t = 10
    else:
        if args != 3:
            raise ValueError(
                'Must provide args for "mixing" "number stages" "number tubes" '
                'e.g. "stage" 3 10'
            )
        args = args[1:]
        mix_style = args[0]
        num_s = int(args[1])
        num_t = int(args[2])

    # set relevant parameter values
    NS = num_s  # number of stages
    NT = num_t  # number of tubes
    solutes = ["Li", "Co"]
    flux = 0.1
    sieving_coefficient = {"Li": 1.3, "Co": 0.5}
    feed = {
        "solvent": 100,  # m^3/hr of water
        "Li": 1.7 * 100,  # kg/hr
        "Co": 17 * 100,  # kg/hr
    }
    diaf = {
        "solvent": 30,  # m^3/hr of water
        "Li": 0.1 * 30,  # kg/hr
        "Co": 0.2 * 30,  # kg/hr
    }
    prec = True

    df = DiafiltrationModel(
        NS=NS,
        NT=NT,
        solutes=solutes,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        feed=feed,
        diafiltrate=diaf,
        precipitate=prec,
        precipitate_yield={
            "permeate": {"Li": 0.81, "Co": 0.01},
            "retentate": {"Li": 0.20, "Co": 0.89},
        },
    )

    mixing = mix_style
    m = df.build_flowsheet(mixing=mixing)
    df.initialize(m, mixing=mixing, precipitate=prec)
    df.unfix_dof(m, mixing=mixing, precipitate=prec)
    m.fs.precipitator["retentate"].volume.fix(500)
    m.fs.precipitator["permeate"].volume.fix(500)
    report_statistics(m)

    # solve model
    # R is used for the Li LB constraint.
    # This can be changed to any desired LB.
    m.recovery_li = 0.8
    solver = SolverFactory("ipopt")
    result = solver.solve(m, tee=True)

    # NOTE These percent recoveries are for precipitators
    m.prec_perc_co.display()
    m.prec_perc_li.display()

    # Print all relevant flow information
    utils.report_values(m)


if __name__ == "__main__":
    main(sys.argv)
