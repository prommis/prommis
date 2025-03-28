#!/usr/bin/env python
#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Executable file for generating and solving diafiltration model."""

from pyomo.environ import SolverFactory, Suffix, TransformationFactory

import idaes.logger as idaeslog
from idaes.core.util.model_statistics import report_statistics

from prommis.nanofiltration.membrane_cascade_flowsheet import utils
from prommis.nanofiltration.membrane_cascade_flowsheet.diafiltration_flowsheet_model import (
    DiafiltrationModel,
)

_log = idaeslog.getLogger(__name__)


def main():
    """Driver for creating diafiltration model."""
    # collect arguments
    mix_style, num_s, num_t = set_arguments()

    solutes = ["Li", "Co"]
    flux = 0.1  # m3 / m2 / h
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
    precipitate = True
    simple_costing = False

    # setup for diafiltration model
    df = DiafiltrationModel(
        NS=num_s,
        NT=num_t,
        solutes=solutes,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        feed=feed,
        diafiltrate=diaf,
        precipitate=precipitate,
        precipitate_yield={
            "permeate": {"Li": 0.81, "Co": 0.01},
            "retentate": {"Li": 0.20, "Co": 0.89},
        },
        simple_costing=simple_costing,
    )

    # model initialization
    m = df.build_flowsheet(mixing=mix_style)
    df.initialize(m, mixing=mix_style)
    df.unfix_dof(m, mixing=mix_style)
    report_statistics(m)

    df.add_costing(m)
    df.add_costing_objectives(m)

    report_statistics(m)

    # set recovery lower bounds
    lithium_recovery = 0.8
    cobalt_recovery = 0.8

    m.fs.split_diafiltrate.inlet.flow_vol.setub(10000)

    solve_scaled_model(
        m,
        L=lithium_recovery,
        C=cobalt_recovery,
        precipitator_simple_costing=simple_costing,
    )

    # # solve model
    # try:
    #     solve_model(m, L=lithium_recovery, C=cobalt_recovery)
    # except Exception:
    #     _log.info("The solver did not return an optimal solution, trying scaling...")
    #     try:
    #         solve_scaled_model(m, L=lithium_recovery, C=cobalt_recovery)
    #     except Exception:
    #         _log.info("Failure to solve scaled model")

    # NOTE These percent recoveries are for precipitators
    m.prec_perc_co.display()
    m.prec_perc_li.display()

    m.fs.costing.total_annualized_cost.display()

    # Print all relevant flow information
    vals = utils.report_values(m)
    utils.visualize_flows(
        num_boxes=num_s, num_sub_boxes=num_t, conf=mix_style, model=vals
    )


def set_arguments(
    mixing="tube",
    num_stages=3,
    num_tubes=10,
):
    mix_style = mixing  # specify 'tube' or 'stage'
    num_s = num_stages  # specify number of stages
    num_t = num_tubes  # specify number of elements/tubes

    return (mix_style, num_s, num_t)


def set_scaling(m, precipitator_simple_costing):
    """
    Apply scaling factors to certain constraints to improve solver performance

    Args:
        m: Pyomo model
    """
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Add scaling factors for poorly scaled variables
    if precipitator_simple_costing == False:
        m.scaling_factor[
            m.fs.precipitator["retentate"].costing.precipitator_diameter
        ] = 1e2
        m.scaling_factor[
            m.fs.precipitator["permeate"].costing.precipitator_diameter
        ] = 1e2
    m.scaling_factor[m.fs.cascade.costing.SEC] = 1e3


def solve_model(m, L, C):
    m.R = L
    m.Rco = C

    solver = SolverFactory("ipopt")
    result = solver.solve(m, tee=True)

    if result.solver.termination_condition == "infeasible":
        raise Exception("The solver returned an infeasible solution")

    if result.solver.termination_condition == "maxIterations":
        raise Exception("The maximum iteration limit was reached")

    return result


def solve_scaled_model(m, L, C, precipitator_simple_costing):
    m.R = L
    m.Rco = C

    scaling = TransformationFactory("core.scale_model")
    solver = SolverFactory("ipopt")

    set_scaling(m, precipitator_simple_costing=precipitator_simple_costing)
    scaled_model = scaling.create_using(m, rename=False)
    result = solver.solve(scaled_model, tee=True)
    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, m)

    return result


if __name__ == "__main__":
    main()
    # Print all relevant flow information
    utils.report_values(m)
