#!/usr/bin/env python
#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Executable file for generating and solving diafiltration model."""

import sys

from pyomo.environ import (
    SolverFactory,
    Suffix,
    TransformationFactory,
    assert_optimal_termination,
    value,
)

from idaes.core.util import to_json, from_json
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import report_statistics
from idaes.core.util.scaling import constraint_autoscale_large_jac

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
    )

    # model initialization
    m = df.build_flowsheet(mixing=mix_style)

    saved_initialization = False
    if saved_initialization:
        from_json(m, fname="initialized_model_stage_3_10")
    else:
        df.initialize(m, mixing=mix_style, precipitate=precipitate)
        to_json(m, fname="initialized_model_stage_3_10")

    df.unfix_dof(m, mixing=mix_style, precipitate=precipitate)
    m.fs.split_diafiltrate.inlet.flow_vol.setub(200)
    report_statistics(m)

    costing = True
    atmospheric_pressure = 101325  # ambient pressure, Pa
    operating_pressure = 145  # nanofiltration operating pressure, psi
    simple_costing = False
    if costing:
        df.add_costing(
            m,
            NS=num_s,
            flux=flux,
            feed=feed,
            diaf=diaf,
            precipitate=precipitate,
            atmospheric_pressure=atmospheric_pressure,
            operating_pressure=operating_pressure,
            simple_costing=simple_costing,
        )
        df.add_costing_objectives(m)

    # set recovery lower bounds
    lithium_recovery = 0.8
    cobalt_recovery = 0.8

    solve_scaled_model(
        m,
        L=lithium_recovery,
        C=cobalt_recovery,
        NS=num_s,
        costing=costing,
        simple_costing=simple_costing,
    )

    dt = DiagnosticsToolbox(m)
    # some flows are at their bounds of zero
    dt.report_numerical_issues()

    if costing:
        # Verify the feed pump operating pressure workaround is valid
        # assume this additional cost is less than half a cent
        if value(m.fs.feed_pump.costing.variable_operating_cost) >= 0.005:
            raise ValueError(
                "The variable  operating cost of the feed pump as calculated in the feed"
                "pump costing block is not negligible. This operating cost is already"
                "accounted for via the membrane's pressure drop specific energy consumption."
            )

    # NOTE These percent recoveries are for precipitators
    m.prec_perc_co.display()
    m.prec_perc_li.display()

    m.fs.costing.total_annualized_cost.display()

    # Print all relevant flow information
    vals = utils.report_values(m)
    utils.visualize_flows(
        num_boxes=num_s, num_sub_boxes=num_t, conf=mix_style, model=vals
    )


def set_scaling(m, NS, costing, simple_costing):
    """
    Apply scaling factors to certain constraints to improve solver performance

    Args:
        m: Pyomo model
    """
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Add scaling factors for poorly scaled variables
    if costing:
        m.scaling_factor[m.fs.costing.aggregate_capital_cost] = 1e-5
        m.scaling_factor[m.fs.costing.aggregate_fixed_operating_cost] = 1e-4
        m.scaling_factor[m.fs.costing.aggregate_variable_operating_cost] = 1e-5
        m.scaling_factor[m.fs.costing.total_capital_cost] = 1e-5
        m.scaling_factor[m.fs.costing.total_operating_cost] = 1e-5
        m.scaling_factor[m.fs.costing.maintenance_labor_chemical_operating_cost] = 1e-4
        m.scaling_factor[m.fs.costing.total_annualized_cost] = 1e-5

        for n in range(1, NS + 1):
            m.scaling_factor[m.fs.stage[n].costing.capital_cost] = 1e-5
            m.scaling_factor[m.fs.stage[n].costing.fixed_operating_cost] = 1e-4
            m.scaling_factor[m.fs.stage[n].costing.membrane_area] = 1e-3

        m.scaling_factor[m.fs.feed_pump.costing.capital_cost] = 1e-4
        m.scaling_factor[m.fs.diafiltrate_pump.costing.capital_cost] = 1e-3
        m.scaling_factor[m.fs.diafiltrate_pump.costing.variable_operating_cost] = 1e-3

        if simple_costing:
            m.scaling_factor[m.fs.feed_pump.costing.pump_power_factor_simple] = 1e-2
            m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e-4
            m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_power_factor_simple] = (
                1e-2
            )

        if not simple_costing:
            m.scaling_factor[m.fs.cascade.costing.variable_operating_cost] = 1e-5
            m.scaling_factor[m.fs.cascade.costing.pressure_drop] = 1e-2
            m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e-4
            m.scaling_factor[m.fs.feed_pump.costing.pump_head] = 1e6
            m.scaling_factor[m.fs.feed_pump.costing.pump_power] = 1e2
            m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_head] = 1e-2
            m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_power] = 1e-5

        for prod in ["retentate", "permeate"]:
            m.scaling_factor[m.fs.precipitator[prod].costing.capital_cost] = 1e-5
            if not simple_costing:
                m.scaling_factor[m.fs.precipitator[prod].costing.base_cost_per_unit] = (
                    1e-4
                )

    # Add scaling factors for poorly scaled constraints
    constraint_autoscale_large_jac(m)


def solve_scaled_model(m, L, C, NS, costing, simple_costing):
    m.recovery_li = L
    m.recovery_co = C

    scaling = TransformationFactory("core.scale_model")
    solver = SolverFactory("ipopt")

    set_scaling(m, NS, costing=costing, simple_costing=simple_costing)
    scaled_model = scaling.create_using(m, rename=False)
    result = solver.solve(scaled_model, tee=True)
    assert_optimal_termination(result)
    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, m)

    return result


if __name__ == "__main__":
    main(sys.argv)
