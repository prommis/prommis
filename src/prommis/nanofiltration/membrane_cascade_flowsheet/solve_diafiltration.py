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
    Objective,
    RangeSet,
    SolverFactory,
    Suffix,
    TransformationFactory,
    assert_optimal_termination,
    units,
)

from idaes.core import (
    UnitModelBlock,
    UnitModelCostingBlock,
)
from idaes.core.util import to_json, from_json
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import report_statistics
from idaes.core.util.scaling import constraint_autoscale_large_jac

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCosting,
    DiafiltrationCostingData,
)
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
    atmospheric_pressure = 101325  # ambient pressure, Pa
    operating_pressure = 145  # nanofiltration operating pressure, psi
    simple_costing = True

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
    m.fs.precipitator["retentate"].volume.fix(500)
    m.fs.precipitator["permeate"].volume.fix(500)
    report_statistics(m)

    costing = True
    if costing:
        add_costing(
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
        add_costing_objectives(m)

    report_statistics(m)

    # set recovery lower bounds
    lithium_recovery = 0.8
    # cobalt_recovery = 0.8

    solve_scaled_model(
        m,
        L=lithium_recovery,
        # C=cobalt_recovery,
        NS=num_s,
        costing=costing,
        simple_costing=simple_costing,
    )

    # TODO: debug solver convergence failures
    dt = DiagnosticsToolbox(m)
    dt.report_numerical_issues()
    # dt.display_variables_at_or_outside_bounds()

    # m.fs.cascade.costing.display()
    # m.fs.feed_pump.costing.display()
    # m.fs.diafiltrate_pump.costing.display()
    # m.fs.precipitator["retentate"].costing.display()
    # m.fs.precipitator["permeate"].costing.display()

    # NOTE These percent recoveries are for precipitators
    # m.prec_perc_co.display()
    # m.prec_perc_li.display()

    # m.fs.costing.total_annualized_cost.display()

    # Print all relevant flow information
    # vals = utils.report_values(m)
    # utils.visualize_flows(
    #     num_boxes=num_s, num_sub_boxes=num_t, conf=mix_style, model=vals
    # )


def add_costing(
    m,
    NS,
    flux,
    feed,
    diaf,
    precipitate,
    atmospheric_pressure,
    operating_pressure,
    simple_costing,
):
    """
    Adds custom costing block to the flowsheet
    """
    m.fs.costing = DiafiltrationCosting()

    # Create dummy variables to store the UnitModelCostingBlocks
    # These are needed because the sieving coefficient model does not account for pressure
    m.fs.cascade = UnitModelBlock()  # to cost the pressure drop
    m.fs.feed_pump = UnitModelBlock()  # to cost feed pump
    m.fs.diafiltrate_pump = UnitModelBlock()  # to cost diafiltrate pump

    # TODO: remove for simple costing
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop,
        costing_method_arguments={
            "water_flux": flux * units.m**3 / units.m**2 / units.h,
            "vol_flow_feed": feed["solvent"] * units.m**3 / units.h,  # cascade feed
            "vol_flow_perm": sum(
                m.fs.split_permeate[i].product.flow_vol[0] for i in RangeSet(NS)
            ),  # cascade permeate
        },
    )
    # TODO: add UnitModelCostingBlock for feed_pump for simple costing
    if not simple_costing:
        m.fs.feed_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_pump,
            costing_method_arguments={
                "inlet_pressure": atmospheric_pressure * units.Pa,  # 14.7 psia
                "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
                * units.psi,  # this should make m.fs.feed_pump.costing.variable_operating_cost ~0
                "inlet_vol_flow": feed["solvent"] * units.m**3 / units.h,  # feed
                "simple_costing": simple_costing,
            },
        )
    m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": atmospheric_pressure * units.Pa,  # 14.7 psia
            "outlet_pressure": operating_pressure * units.psi,
            "inlet_vol_flow": diaf["solvent"] * units.m**3 / units.h,  # diafiltrate
            "simple_costing": simple_costing,
        },
    )
    # membrane stage cost blocks
    for n in range(1, NS + 1):
        m.fs.stage[n].costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_membranes,
            costing_method_arguments={
                "membrane_length": m.fs.stage[n].length,
                "membrane_width": m.fs.stage[n].width,
            },
        )

    if precipitate:
        for prod in ["retentate", "permeate"]:
            m.fs.precipitator[prod].costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method=DiafiltrationCostingData.cost_precipitator,
                costing_method_arguments={
                    "precip_volume": m.fs.precipitator[prod].volume,
                    "simple_costing": simple_costing,
                },
            )

    m.fs.costing.cost_process()


def add_costing_objectives(m):
    """
    Method to add cost objective to flowsheet for performing optimization

    Args:
        m: Pyomo model
    """
    m.co_obj.deactivate()
    m.li_lb.deactivate()
    m.prec_co_obj.deactivate()
    m.prec_co_lb.activate()

    def cost_obj(m):
        return m.fs.costing.total_annualized_cost

    m.cost_objecticve = Objective(rule=cost_obj)


def set_scaling(m, NS, costing, simple_costing):
    """
    Apply scaling factors to certain constraints to improve solver performance

    Args:
        m: Pyomo model
    """
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # TODO: determine appropriate scaling factors

    # Add scaling factors for poorly scaled variables
    if costing:
        m.scaling_factor[m.fs.costing.aggregate_capital_cost] = 1e-6
        m.scaling_factor[m.fs.costing.aggregate_fixed_operating_cost] = 1e-4
        m.scaling_factor[m.fs.costing.aggregate_variable_operating_cost] = 1e-5
        m.scaling_factor[m.fs.costing.total_capital_cost] = 1e-6
        m.scaling_factor[m.fs.costing.total_operating_cost] = 1e-5
        m.scaling_factor[m.fs.costing.maintenance_labor_chemical_operating_cost] = 1e-5

        for n in range(1, NS + 1):
            m.scaling_factor[m.fs.stage[n].costing.capital_cost] = 1e-4
            m.scaling_factor[m.fs.stage[n].costing.fixed_operating_cost] = 1e-4
            m.scaling_factor[m.fs.stage[n].costing.membrane_area] = 1e-3

        m.scaling_factor[m.fs.cascade.costing.variable_operating_cost] = 1e-5
        m.scaling_factor[m.fs.cascade.costing.pressure_drop] = 1e-2
        m.scaling_factor[m.fs.cascade.costing.SEC] = 1e1

        # m.scaling_factor[m.fs.feed_pump.costing.capital_cost] = 1e-4
        # m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e3

        m.scaling_factor[m.fs.diafiltrate_pump.costing.capital_cost] = 1e-4
        m.scaling_factor[m.fs.diafiltrate_pump.costing.variable_operating_cost] = 1e-4

        if simple_costing:
            m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_power_factor_simple] = (
                1e-2
            )

        if not simple_costing:
            m.scaling_factor[m.fs.feed_pump.costing.capital_cost] = 1e-4
            m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e3
            m.scaling_factor[m.fs.feed_pump.costing.pump_head] = 1e6
            # m.scaling_factor[m.fs.feed_pump.costing.pump_power] = 1e2
            m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_head] = 1e-2
            # m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_power] = 1e-5

        for prod in ["retentate", "permeate"]:
            m.scaling_factor[m.fs.precipitator[prod].costing.capital_cost] = 1e-5
            if simple_costing == False:
                m.scaling_factor[
                    m.fs.precipitator[prod].costing.precipitator_diameter
                ] = 1e1

    # Add scaling factors for poorly scaled constraints
    constraint_autoscale_large_jac(m)


def solve_scaled_model(m, L, NS, costing, simple_costing):
    m.recovery_li = L
    # m.Rco = C

    scaling = TransformationFactory("core.scale_model")
    solver = SolverFactory("ipopt")

    set_scaling(m, NS, costing=costing, simple_costing=simple_costing)
    scaled_model = scaling.create_using(m, rename=False)
    result = solver.solve(scaled_model, tee=True)
    # assert_optimal_termination(result)
    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, m)

    return result


if __name__ == "__main__":
    main(sys.argv)
