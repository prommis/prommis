#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Example of optimizing a flowsheet using the REE costing library.
"""

__author__ = (
    "Costing Team (B. Paul, A. Fritz, A. Ojo, L. Deng, A. Dasgupta, and M. Zamarripa)"
)
__version__ = "1.0.0"

from pyomo.environ import (
    Constraint,
    Expression,
    Objective,
    Param,
    Suffix,
    TransformationFactory,
    minimize,
)
from pyomo.environ import units as pyunits
from pyomo.environ import value

from idaes.core import UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox, degrees_of_freedom
from idaes.core.util.scaling import constraint_autoscale_large_jac

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCostingData,
)
from prommis.nanofiltration.diafiltration import (
    add_product_constraints,
    build_model,
    initialize_model,
    solve_model,
    unfix_opt_variables,
)
from prommis.uky.costing.ree_plant_capcost import QGESSCosting


# --- I/O reporting helpers (keeps the script's simple print style) ---
def _v(x):
    try:
        return float(value(x))
    except Exception:
        return None


def _purity_mass(stream, li_key="Li", co_key="Co"):
    """
    Mass-based purity among Li and Co using conc_mass_solute if available.
    Returns (purity_Li, purity_Co) or (None, None) if not computable.
    """
    try:
        li = _v(stream.conc_mass_solute[0, li_key])
        co = _v(stream.conc_mass_solute[0, co_key])
    except Exception:
        return None, None
    if li is None or co is None:
        return None, None
    tot = li + co
    if tot and tot > 0:
        return li / tot, co / tot
    return None, None


def print_io_snap(fs, tag="STATE"):
    """
    Prints:
      - Initial FEED flow + Li/Co concentrations (stage3 retentate side-stream @ element 10)
      - Initial DIAFILTRATE flow + Li/Co concentrations (mix2.inlet_1)
      - PRODUCT PERMEATE (stage3.permeate_outlet): flow + Li/Co purity + Li recovery fraction
      - PRODUCT RETENTATE (stage1.retentate_outlet): flow + Li/Co purity + Co recovery fraction
    """
    root = fs.parent_block()

    print("\n" + "=" * 72)
    print(f"I/O SNAPSHOT: {tag}")
    print("=" * 72)

    # -------- Initial feeds (fixed, not inside any stage loop) --------
    # FEED: side-stream into stage 3 at element 10
    feed_state = fs.stage3.retentate_side_stream_state[0, 10]
    print("[FEED  (initial; stage3.retentate_side_stream_state[0,10])]")
    print(f"  flow_vol (m\u00b3/hr): {_v(feed_state.flow_vol)}")
    print(f"  conc_Li (kg/m\u00b3):  {_v(feed_state.conc_mass_solute['Li'])}")
    print(f"  conc_Co (kg/m\u00b3):  {_v(feed_state.conc_mass_solute['Co'])}")

    # DIAFILTRATE: mixer 2 inlet_1
    diaf = fs.mix2.inlet_1
    print("\n[DIAFILTRATE (initial; mix2.inlet_1)]")
    print(f"  flow_vol (m\u00b3/hr): {_v(diaf.flow_vol[0])}")
    print(f"  conc_Li (kg/m\u00b3):  {_v(diaf.conc_mass_solute[0,'Li'])}")
    print(f"  conc_Co (kg/m\u00b3):  {_v(diaf.conc_mass_solute[0,'Co'])}")

    # -------- Product streams (final, not intermediate elements) --------
    perm = fs.stage3.permeate_outlet  # product permeate (final)
    ret = fs.stage1.retentate_outlet  # product retentate (final)

    print("\n[PRODUCT PERMEATE (stage3.permeate_outlet)]")
    print(f"  flow_vol (m\u00b3/hr): {_v(perm.flow_vol[0])}")
    print(f"  Li concentration (kg/m\u00b3): {_v(perm.conc_mass_solute[0, 'Li'])}")
    print(f"  Co concentration (kg/m\u00b3): {_v(perm.conc_mass_solute[0, 'Co'])}")
    print(f"  Li recovery fraction: {_v(fs.Li_recovery)}")

    li_p, co_p = _purity_mass(perm)
    print(f"  purity_Li: {li_p if li_p is not None else 'N/A'}")
    print(f"  purity_Co: {co_p if co_p is not None else 'N/A'}")

    print("\n[PRODUCT RETENTATE (stage1.retentate_outlet)]")
    print(f"  flow_vol (m\u00b3/hr): {_v(ret.flow_vol[0])}")
    print(f"  Li concentration (kg/m\u00b3): {_v(ret.conc_mass_solute[0, 'Li'])}")
    print(f"  Co concentration (kg/m\u00b3): {_v(ret.conc_mass_solute[0, 'Co'])}")
    print(f"  Co recovery fraction: {_v(fs.Co_recovery)}")
    li_r, co_r = _purity_mass(ret)
    print(f"  purity_Li: {li_r if li_r is not None else 'N/A'}")
    print(f"  purity_Co: {co_r if co_r is not None else 'N/A'}")

    # --------- Selectivity_coefficient, membrane width, operating pressure ----------
    sel_Li = (
        _v(fs.sieving_coefficient["Li"]) if hasattr(fs, "sieving_coefficient") else None
    )
    sel_Co = (
        _v(fs.sieving_coefficient["Co"]) if hasattr(fs, "sieving_coefficient") else None
    )
    print("\n[PARAMETERS]")
    print(f"  sieving_coefficient_Li : {sel_Li if sel_Li is not None else 'N/A'}")
    print(f"  sieving_coefficient_Co: {sel_Co if sel_Co is not None else 'N/A'}")
    print(
        f"  membrane_width (m.w): {_v(getattr(root, 'w', None)) if hasattr(root, 'w') else 'N/A'}"
    )
    print(
        f"  operating_pressure (psi): {_v(getattr(root, 'operating_pressure', None)) if hasattr(root,'operating_pressure') else 'N/A'}"
    )

    print("=" * 72 + "\n")


def build_costing(m):

    # Create dummy variables to store the UnitModelCostingBlocks
    # These are needed because the sieving coefficient model does not account for pressure
    m.fs.cascade = UnitModelBlock()  # to cost the pressure drop
    m.fs.feed_pump = UnitModelBlock()  # to cost feed pump
    m.fs.diafiltrate_pump = UnitModelBlock()  # to cost diafiltrate pump

    m.fs.costing = QGESSCosting()
    m.fs.stage1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage1.length,
            "membrane_width": m.w,
        },
    )
    m.fs.stage2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage2.length,
            "membrane_width": m.w,
        },
    )
    m.fs.stage3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.stage3.length,
            "membrane_width": m.w,
        },
    )
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop,
        costing_method_arguments={
            "water_flux": m.Jw,
            "vol_flow_feed": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # cascade feed
            "vol_flow_perm": m.fs.stage3.permeate_outlet.flow_vol[
                0
            ],  # cascade permeate
        },
    )
    m.fs.feed_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,  # 14.7 psia
            "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
            * pyunits.psi,  # this should make m.fs.feed_pump.costing.variable_operating_cost ~0
            "inlet_vol_flow": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # feed
        },
    )
    m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,  # 14.7 psia
            "outlet_pressure": m.operating_pressure,
            "inlet_vol_flow": m.fs.stage3.retentate_inlet.flow_vol[0],  # diafiltrate
        },
    )

    m.fs.Li_product = Expression(
        expr=pyunits.convert(
            m.fs.stage3.permeate_outlet.flow_vol[0]
            * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"],
            to_units=pyunits.kg / pyunits.h,
        )
    )

    m.fs.Co_product = Expression(
        expr=pyunits.convert(
            m.fs.stage1.retentate_outlet.flow_vol[0]
            * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"],
            to_units=pyunits.kg / pyunits.h,
        )
    )

    m.fs.Li_feed = Expression(
        expr=pyunits.convert(
            m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"],
            to_units=pyunits.kg/pyunits.h,
        )
    )

    m.fs.Co_feed = Expression(
        expr=pyunits.convert(
            m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"],
            to_units=pyunits.kg/pyunits.h,
        )
    )

    m.fs.Li_recovery = Expression(expr=m.fs.Li_product / m.fs.Li_feed)
    m.fs.Co_recovery = Expression(expr=m.fs.Co_product / m.fs.Co_feed)


    # Operation parameters to use later
    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    m.fs.annual_operating_hours = Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=True,
        units=pyunits.hours / pyunits.a,
    )

    # Define the recovery rate

    m.fs.recovery_rate_per_year = Expression(
        expr=pyunits.convert(
            (m.fs.Li_product + m.fs.Co_product) * m.fs.annual_operating_hours,
            to_units=pyunits.kg / pyunits.year,
        )
    )

    m.fs.costing.build_process_costs(
        Lang_factor=2,  # includes installation, material, construction
        labor_types=[],  # labor costs already included in maintenance, admin
        fixed_OM=True,
        variable_OM=True,
        resources=[],  # List of strings of resources.
        rates=[],  # Resource consumption rate.
        prices={},  # Resource prices is not considered.
        recovery_rate_per_year=m.fs.recovery_rate_per_year,
        CE_index_year="2021",
    )


def apply_design_limits(
    m,
    Lmin=0.1,
    Lmax=10000.0,
    Li_max=20 * pyunits.kg / pyunits.m**3,
    Co_max=200 * pyunits.kg / pyunits.m**3,
    sc_min=0.01,
    sc_max=0.99,
):  # assume a boundary
    # 1) Length bounds for each stage
    for st in (m.fs.stage1, m.fs.stage2, m.fs.stage3):
        st.length.setlb(Lmin * pyunits.m)
        st.length.setub(Lmax * pyunits.m)

    # 2) Product concentration variable upper bounds
    #    Permeate product: stage3.permeate_outlet
    #    Retentate product: stage1.retentate_outlet
    perm = m.fs.stage3.permeate_outlet
    ret = m.fs.stage1.retentate_outlet

    # Permeate specs (UB only, keep LB whatever the model already has)
    perm.conc_mass_solute[0, "Li"].setub(Li_max)
    perm.conc_mass_solute[0, "Co"].setub(Co_max)

    # Retentate specs
    ret.conc_mass_solute[0, "Li"].setub(Li_max)
    ret.conc_mass_solute[0, "Co"].setub(Co_max)

    # 3) Stage-cut (sc) bounds according to : Q_perm / Q_feed in [sc_min, sc_max]
    stages = [m.fs.stage1, m.fs.stage2, m.fs.stage3]
    for i, st in enumerate(stages, start=1):
        Q_perm = st.permeate_outlet.flow_vol[0]
        if i < 3:
            Q_feed = st.retentate_inlet.flow_vol[0]
        else:
            Q_feed = (
                st.retentate_inlet.flow_vol[0]
                + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            )
        m.add_component(f"stage{i}_cut_lb", Constraint(expr=Q_perm >= sc_min * Q_feed))
        m.add_component(f"stage{i}_cut_ub", Constraint(expr=Q_perm <= sc_max * Q_feed))


def compute_stage_cuts(m, eps=1e-12):
    # Stage 1
    Qf1 = m.fs.stage1.retentate_inlet.flow_vol[0]
    Qp1 = m.fs.stage1.permeate_outlet.flow_vol[0]
    sc1 = value(Qp1 / Qf1) if abs(value(Qf1)) > eps else float("nan")

    # Stage 2
    Qf2 = m.fs.stage2.retentate_inlet.flow_vol[0]
    Qp2 = m.fs.stage2.permeate_outlet.flow_vol[0]
    sc2 = value(Qp2 / Qf2) if abs(value(Qf2)) > eps else float("nan")

    # Stage 3 (add fresh side-stream @ element 10)
    Qf3 = (
        m.fs.stage3.retentate_inlet.flow_vol[0]
        + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
    )
    Qp3 = m.fs.stage3.permeate_outlet.flow_vol[0]
    sc3 = value(Qp3 / Qf3) if abs(value(Qf3)) > eps else float("nan")

    return [sc1, sc2, sc3]


def print_stage_cuts(m, label="STAGE CUTS"):
    sc = compute_stage_cuts(m)
    print("\n" + "=" * 60)
    print(f"{label}")
    print("=" * 60)
    for i, s in enumerate(sc, start=1):
        print(f"Stage {i} cut (Q_perm/Q_feed): {s:.6f}")
    print("=" * 60 + "\n")


def apply_sieving_bounds_and_unfix(
    m, li_bounds=(0.75, 1.3), co_bounds=(0.05, 0.5), li_start=1.3, co_start=0.5
):
    sc = m.fs.sieving_coefficient
    # Li
    sc["Li"].setlb(li_bounds[0])
    sc["Li"].setub(li_bounds[1])
    if li_start is not None:
        sc["Li"].set_value(li_start)
    sc["Li"].unfix()
    # Co
    sc["Co"].setlb(co_bounds[0])
    sc["Co"].setub(co_bounds[1])
    if co_start is not None:
        sc["Co"].set_value(co_start)
    sc["Co"].unfix()


def build_optimization(m):
    apply_design_limits(
        m,
        Lmin=0.1,
        Lmax=10000.0,
        Li_max=20 * pyunits.kg / pyunits.m**3,
        Co_max=200 * pyunits.kg / pyunits.m**3,
        sc_min=0.01,
        sc_max=0.99,
    )
    apply_sieving_bounds_and_unfix(
        m, li_bounds=(0.75, 1.3), co_bounds=(0.05, 0.5), li_start=1.3, co_start=0.5
    )

    def cost_obj(m):
        return m.fs.costing.cost_of_recovery

    m.cost_objective = Objective(rule=cost_obj, sense=minimize)

    unfix_opt_variables(m)
    add_product_constraints(m, Li_recovery_bound=0.945, Co_recovery_bound=0.635)


def scale_and_solve_model(m):

    # set scaling
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Add scaling factors for poorly scaled constraints
    constraint_autoscale_large_jac(m)

    # Add scaling factors for poorly scaled variables
    m.scaling_factor[m.fs.cascade.costing.variable_operating_cost] = 1e-5
    m.scaling_factor[m.fs.feed_pump.costing.capital_cost] = 1e-5
    m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e3
    m.scaling_factor[m.fs.feed_pump.costing.pump_head] = 1e5
    m.scaling_factor[m.fs.feed_pump.costing.pump_power] = 1e2
    m.scaling_factor[m.fs.diafiltrate_pump.costing.variable_operating_cost] = 1e-5
    m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_head] = 1e-1
    m.scaling_factor[m.fs.diafiltrate_pump.costing.pump_power] = 1e-5

    m.scaling_factor[m.fs.costing.annualized_cost] = 1e-5
    m.scaling_factor[m.fs.costing.total_variable_OM_cost] = 1e-5
    m.scaling_factor[m.fs.costing.total_fixed_OM_cost] = 1e-4

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)
    solve_model(scaled_model, tee=False)
    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, m)


if __name__ == "__main__":
    m = build_model()

    dt = DiagnosticsToolbox(m)
    assert degrees_of_freedom(m) == 0

    initialize_model(m)
    solve_model(m, tee=False)
    dt.assert_no_numerical_warnings()

    build_costing(m)

    solve_model(m, tee=False)
    dt.assert_no_numerical_warnings()

    # fix some initial long stage lengths so optimization starts from a feasible place with recovery bounds
    # chose the initial guess based on Multi-Stream Contactor Tutorial (Solution).
    m.fs.stage1.length.fix(754)
    m.fs.stage2.length.fix(758)
    m.fs.stage3.length.fix(756)

    solve_model(m, tee=False)
    dt.assert_no_numerical_warnings()
    print_io_snap(m.fs, tag="BEFORE OPTIMIZATION")
    print_stage_cuts(m, label="STAGE CUTS — BEFORE OPTIMIZATION")

    print(
        "\nStage lengths prior to optimization: ",
        [m.fs.stage1.length.value, m.fs.stage2.length.value, m.fs.stage3.length.value],
        " ",
        pyunits.get_units(m.fs.stage1.length),
    )
    print(
        "\nTotal length: ",
        value(m.fs.stage1.length + m.fs.stage2.length + m.fs.stage3.length),
        " ",
        pyunits.get_units(m.fs.stage1.length),
    )
    print("\nCost results prior to optimization:")
    m.fs.costing.report()

    build_optimization(m)

    scale_and_solve_model(m)
    dt.assert_no_numerical_warnings()

    print_io_snap(m.fs, tag="AFTER OPTIMIZATION")
    print_stage_cuts(m, label="STAGE CUTS — AFTER OPTIMIZATION")
    print(
        "\nStage lengths after optimization: ",
        [m.fs.stage1.length.value, m.fs.stage2.length.value, m.fs.stage3.length.value],
        " ",
        pyunits.get_units(m.fs.stage1.length),
    )
    print(
        "\nTotal length: ",
        value(m.fs.stage1.length + m.fs.stage2.length + m.fs.stage3.length),
        " ",
        pyunits.get_units(m.fs.stage1.length),
    )
    print("\nCost results after optimization:")
    m.fs.costing.report()
