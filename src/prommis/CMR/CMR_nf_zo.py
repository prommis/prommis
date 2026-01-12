#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    Set,
    Suffix,
    TransformationFactory,
    value,
    Var,
)

from idaes.core import FlowsheetBlock
from idaes.core.scaling import set_scaling_factor
from idaes.core.solvers import get_solver

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.nanofiltration_ZO import NanofiltrationZO
from idaes.core import UnitModelCostingBlock
from watertap.costing.unit_models.nanofiltration import cost_nanofiltration
from watertap.costing import WaterTAPCosting
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.misc import StrEnum
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from pyomo.network import Arc

__author__ = "Chenyu Wang"


class Case(StrEnum):
    case1 = "case1"
    case2 = "case2"


def CMR_nf_case(case=Case.case1, simplified_routine=False):
    m = build(case=case, simplified_routine=simplified_routine)
    set_scaling(m)
    initialize_system(m)
    m, results = solve(m)
    m.fs.NF.area.unfix()
    m.fs.NF.flux_vol_solvent[0, "H2O"].fix(
        10 * pyo.units.L / (pyo.units.m**2 * pyo.units.hr)
    )
    m, results = solve(m)
    display_performance_metrics(m)
    display_costing(m)

    return m, results


def build(case, simplified_routine=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.properties = MCASParameterBlock(
        solute_list=[
            "Co_2+",
            "Ca_2+",
            "Cu_2+",
            "Fe_3+",
            "Nd_3+",
            "Ni_2+",
            "Pr_3+",
            "Na_+",
            "Cr_6+",
            "Sn_2+",
            "Zn_2+",
            "Pb_2+",
            "Dy_3+",
            "B_3+",
            "Gd_3+",
            "Mn_2+",
            "Si_4+",
            "SO4_2-",
        ],
        # diffusivity_data={
        #     ("Liq", "Ca_2+"): 9.2e-10,
        #     ("Liq", "Mg_2+"): 7.06e-10,
        #     ("Liq", "Na_+"): 1.33e-09,
        #     ("Liq", "SO4_2-"): 2.03e-09,
        # },
        # mw_data={
        #     "H2O": 0.018,
        #     "Ca_2+": 0.04,
        #     "Mg_2+": 0.024,
        #     "Na_+": 0.023,
        #     "SO4_2-": 0.035,
        # },
        # stokes_radius_data={
        #     "Ca_2+": 3.09e-10,
        #     "Mg_2+": 3.47e-10,
        #     "SO4_2-": 1.21e-10,
        #     "Na_+": 1.84e-10,
        # },
        charge={
            "Co_2+": 2,
            "Ca_2+": 2,
            "Cu_2+": 2,
            "Fe_3+": 3,
            "Nd_3+": 3,
            "Ni_2+": 2,
            "Pr_3+": 3,
            "Na_+": 1,
            "Cr_6+": 6,
            "Sn_2+": 2,
            "Zn_2+": 2,
            "Pb_2+": 2,
            "Dy_3+": 3,
            "B_3+": 3,
            "Gd_3+": 3,
            "Mn_2+": 3,
            "Si_4+": 4,
            "SO4_2-": -2,
        },
    )

    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.NF = NanofiltrationZO(
        property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.NF.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing, costing_method=cost_nanofiltration
    )
    membrane_cost = 50
    m.fs.NF.costing.membrane_cost.set_value(membrane_cost)

    m.fs.s01 = Arc(source=m.fs.P1.outlet, destination=m.fs.NF.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Fix other inlet state variables
    # fully specify system
    # Case 1
    if case is Case.case1:
        m.fs.P1.control_volume.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.0000281653,
                ("conc_mass_phase_comp", ("Liq", "Co_2+")): value(
                    105.732e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Ca_2+")): value(
                    11.3742e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Cu_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Fe_3+")): value(
                    17454.9915e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Nd_3+")): value(
                    521.0505e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Ni_2+")): value(
                    15.9399e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Pr_3+")): value(
                    134.4078e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Na_+")): value(
                    6.3279e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Cr_6+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Sn_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Zn_2+")): value(
                    668.3544e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Pb_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Dy_3+")): value(
                    0.801e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "B_3+")): value(
                    72.0099e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Gd_3+")): value(
                    2.6433e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Mn_2+")): value(
                    40.2903e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Si_4+")): value(
                    20.1852e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "SO4_2-")): value(
                    10000e-3
                ),  # feed mass concentration
            },  # volumetric feed flowrate [-]
            hold_state=True,  # fixes the calculated component mass flow rates
        )
    elif case is Case.case2:
        m.fs.P1.control_volume.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.0000281653,
                ("conc_mass_phase_comp", ("Liq", "Co_2+")): value(
                    7.503e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Ca_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Cu_2+")): value(
                    72.2871e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Fe_3+")): value(
                    1248.983e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Nd_3+")): value(
                    1.9065e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Ni_2+")): value(
                    1.7753e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Pr_3+")): value(
                    0.4264e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Na_+")): value(
                    3.3907e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Cr_6+")): value(
                    0.0287e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Sn_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Zn_2+")): value(
                    223.2409e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Pb_2+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Dy_3+")): value(
                    0.0738e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "B_3+")): value(
                    0
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Gd_3+")): value(
                    0.2214e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Mn_2+")): value(
                    0.3895e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "Si_4+")): value(
                    1.1808e-3
                ),  # feed mass concentration
                ("conc_mass_phase_comp", ("Liq", "SO4_2-")): value(
                    10000e-3
                ),  # feed mass concentration
            },  # volumetric feed flowrate [-]
            hold_state=True,  # fixes the calculated component mass flow rates
        )

    # if case is Case.case1:
    #     m.fs.NF.feed_side.properties_in.calculate_state(
    #         var_args={
    #             ("flow_vol_phase", "Liq"): 1e-2,
    #             ("conc_mass_phase_comp", ("Liq", "Co_2+")): value(
    #                 105.732e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Ca_2+")): value(
    #                 11.3742e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Cu_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Fe_3+")): value(
    #                 17454.9915e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Nd_3+")): value(
    #                 521.0505e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Ni_2+")): value(
    #                 15.9399e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Pr_3+")): value(
    #                 134.4078e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Na_+")): value(
    #                 6.3279e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Cr_6+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Sn_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Zn_2+")): value(
    #                 668.3544e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Pb_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Dy_3+")): value(
    #                 0.801e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "B_3+")): value(
    #                 72.0099e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Gd_3+")): value(
    #                 2.6433e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Mn_2+")): value(
    #                 40.2903e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Si_4+")): value(
    #                 20.1852e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "SO4_2-")): value(
    #                 10000e-3
    #             ),  # feed mass concentration
    #         },  # volumetric feed flowrate [-]
    #         hold_state=True,  # fixes the calculated component mass flow rates
    #     )
    # elif case is Case.case2:
    #     m.fs.NF.feed_side.properties_in.calculate_state(
    #         var_args={
    #             ("flow_vol_phase", "Liq"): 1e-2,
    #             ("conc_mass_phase_comp", ("Liq", "Co_2+")): value(
    #                 7.503e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Ca_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Cu_2+")): value(
    #                 72.2871e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Fe_3+")): value(
    #                 1248.983e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Nd_3+")): value(
    #                 1.9065e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Ni_2+")): value(
    #                 1.7753e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Pr_3+")): value(
    #                 0.4264e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Na_+")): value(
    #                 3.3907e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Cr_6+")): value(
    #                 0.0287e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Sn_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Zn_2+")): value(
    #                 223.2409e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Pb_2+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Dy_3+")): value(
    #                 0.0738e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "B_3+")): value(
    #                 0
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Gd_3+")): value(
    #                 0.2214e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Mn_2+")): value(
    #                 0.3895e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "Si_4+")): value(
    #                 1.1808e-3
    #             ),  # feed mass concentration
    #             ("conc_mass_phase_comp", ("Liq", "SO4_2-")): value(
    #                 10000e-3
    #             ),  # feed mass concentration
    #         },  # volumetric feed flowrate [-]
    #         hold_state=True,  # fixes the calculated component mass flow rates
    #     )

    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(1e3)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Co_2+"].fix(0.01)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Ca_2+"].fix(0.019)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Cu_2+"].fix(0.001)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Fe_3+"].fix(1.190)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Nd_3+"].fix(0.020)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Ni_2+"].fix(0.012)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Pr_3+"].fix(0.006)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Na_+"].fix(0.005)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Cr_6+"].fix(0.001)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Sn_2+"].fix(0.001)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Zn_2+"].fix(0.002)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Pb_2+"].fix(0.001)
    # m.fs.NF.feed_side.properties_in[0].flow_mol_phase_comp["Liq", "Dy_3+"].fix(0.001)

    # Pump
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(10 * pyo.units.bar)

    # Nanofiltration
    m.fs.NF.feed_side.properties_in[0].temperature.fix(298.15)
    # m.fs.NF.feed_side.properties_in[0].pressure.fix(10*pyo.units.bar)

    m.fs.NF.recovery_vol_phase[0, "Liq"].fix(0.95)
    # m.fs.NF.flux_vol_solvent.fix(1.67e-6)
    m.fs.NF.rejection_phase_comp.fix(1e-10)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Co_2+"].fix(0.98)
    # m.fs.NF.rejection_phase_comp[0, "Liq", "Ca_2+"].fix(0.92)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Cu_2+"].fix(0.98)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Fe_3+"].fix(0.98)
    # m.fs.NF.rejection_phase_comp[0, "Liq", "Nd_3+"].fix(0.9945)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Ni_2+"].fix(0.98)
    # m.fs.NF.rejection_phase_comp[0, "Liq", "Pr_3+"].fix(0.95)  # Pr 59; Nd 60
    # m.fs.NF.rejection_phase_comp[0, "Liq", "Na_+"].fix(0.796)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Cr_6+"].fix(0.93)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Sn_2+"].fix(0.95)  # Sn 50; Nd 60; Zn 30
    m.fs.NF.rejection_phase_comp[0, "Liq", "Zn_2+"].fix(0.98)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Pb_2+"].fix(0.99)
    # m.fs.NF.rejection_phase_comp[0, "Liq", "Dy_3+"].fix(0.95)  # Dy 66; Nd 60
    m.fs.NF.rejection_phase_comp[0, "Liq", "Mn_2+"].fix(0.89)

    if simplified_routine is False:
        m.fs.NF.rejection_phase_comp[0, "Liq", "Ca_2+"].fix(0.92)
        m.fs.NF.rejection_phase_comp[0, "Liq", "Nd_3+"].fix(0.9945)
        m.fs.NF.rejection_phase_comp[0, "Liq", "Pr_3+"].fix(0.95)  # Pr 59; Nd 60
        m.fs.NF.rejection_phase_comp[0, "Liq", "Na_+"].fix(0.796)
        m.fs.NF.rejection_phase_comp[0, "Liq", "Dy_3+"].fix(0.95)  # Dy 66; Nd 60
        m.fs.NF.rejection_phase_comp[0, "Liq", "Gd_3+"].fix(0.95)

    m.fs.NF.area.fix(10)
    # m.fs.NF.flux_vol_solvent[0, "H2O"].fix(9.51 * pyo.units.L / (pyo.units.m**2 * pyo.units.hr))
    # m.fs.NF.deltaP.fix(-9*pyo.units.bar)
    m.fs.NF.feed_side.properties_out[0.0].pressure.fix(1 * pyo.units.bar)
    # m.fs.NF.feed_side.properties_in[0].assert_electroneutrality(
    #     defined_state=True, adjust_by_ion="SO4_2-"
    # )

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.NF.properties_permeate[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.NF.properties_permeate[0].flow_vol)

    return m


def set_scaling(m):
    # Scale model
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "H2O"], 1e5)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Co_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Cu_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Fe_3+"], 1e0)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Nd_3+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Ni_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Pr_3+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Cr_6+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Sn_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Zn_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Pb_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Dy_3+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "B_3+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Gd_3+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Mn_2+"], 1e2)
    set_scaling_factor(m.fs.NF.inlet.flow_mol_phase_comp[0, "Liq", "Si_4+"], 1e2)
    calculate_scaling_factors(m)


def initialize_system(m):
    # Initialize system
    # m.fs.NF.initialize()
    try:
        m.fs.NF.initialize()
    except:
        pass
    m.fs.costing.initialize()


def solve(m):
    # Solve
    solver = get_solver(
        "ipopt_v2", writer_config={"scale_model": True, "linear_presolve": True}
    )
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m, results


def display_performance_metrics(m):
    print("\n---- Feed Metrics ----")
    f_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].flow_vol,
        to_units=pyo.units.gal / pyo.units.day,
    )
    print(f"Influent flow: " f"{pyo.value(f_in):.3g}" f"{pyo.units.get_units(f_in)}")
    f_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )
    Co_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Co_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Co2+ feed mass concentration: "
        f"{pyo.value(Co_in):.3g}"
        f"{pyo.units.get_units(Co_in)}"
    )
    Ca_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Ca_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ca2+ feed mass concentration: "
        f"{pyo.value(Ca_in):.3g}"
        f"{pyo.units.get_units(Ca_in)}"
    )
    Cu_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Cu_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Cu2+ feed mass concentration: "
        f"{pyo.value(Cu_in):.3g}"
        f"{pyo.units.get_units(Cu_in)}"
    )
    Fe_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Fe_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Fe3+ feed mass concentration: "
        f"{pyo.value(Fe_in):.3g}"
        f"{pyo.units.get_units(Fe_in)}"
    )
    Nd_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Nd_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Nd3+ feed mass concentration: "
        f"{pyo.value(Nd_in):.3g}"
        f"{pyo.units.get_units(Nd_in)}"
    )
    Ni_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Ni_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ni2+ feed mass concentration: "
        f"{pyo.value(Ni_in):.3g}"
        f"{pyo.units.get_units(Ni_in)}"
    )
    Pr_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Pr_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Pr3+ feed mass concentration: "
        f"{pyo.value(Pr_in):.3g}"
        f"{pyo.units.get_units(Pr_in)}"
    )
    Na_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Na_+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Na+ feed mass concentration: "
        f"{pyo.value(Na_in):.3g}"
        f"{pyo.units.get_units(Na_in)}"
    )
    Cr_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Cr_6+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Cr6+ feed mass concentration: "
        f"{pyo.value(Cr_in):.3g}"
        f"{pyo.units.get_units(Cr_in)}"
    )
    Sn_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Sn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Sn2+ feed mass concentration: "
        f"{pyo.value(Sn_in):.3g}"
        f"{pyo.units.get_units(Sn_in)}"
    )
    Zn_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Zn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Zn2+ feed mass concentration: "
        f"{pyo.value(Zn_in):.3g}"
        f"{pyo.units.get_units(Zn_in)}"
    )
    # Pb_in = pyo.units.convert(
    #     m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Pb_2+"],
    #     to_units=pyo.units.mg / pyo.units.L,
    # )
    # print(
    #     f"Pb2+ feed mass concentration: "
    #     f"{pyo.value(Pb_in):.3g}"
    #     f"{pyo.units.get_units(Pb_in)}"
    # )
    Dy_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Dy_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Dy3+ feed mass concentration: "
        f"{pyo.value(Dy_in):.3g}"
        f"{pyo.units.get_units(Dy_in)}"
    )
    B_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "B_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"B3+ feed mass concentration: "
        f"{pyo.value(B_in):.3g}"
        f"{pyo.units.get_units(B_in)}"
    )
    Gd_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Gd_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Gd3+ feed mass concentration: "
        f"{pyo.value(Gd_in):.3g}"
        f"{pyo.units.get_units(Gd_in)}"
    )
    Mn_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Mn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mn2+ feed mass concentration: "
        f"{pyo.value(Mn_in):.3g}"
        f"{pyo.units.get_units(Mn_in)}"
    )
    Si_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "Si_4+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Si4+ feed mass concentration: "
        f"{pyo.value(Si_in):.3g}"
        f"{pyo.units.get_units(Si_in)}"
    )

    print("\n---- Permeate Metrics ----")
    f_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].flow_vol,
        to_units=pyo.units.gal / pyo.units.day,
    )
    print(
        f"Permeate flow: "
        f"{pyo.value(f_permeate):.3g}"
        f"{pyo.units.get_units(f_permeate)}"
    )
    Co_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Co_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Co2+ permeate mass concentration: "
        f"{pyo.value(Co_permeate):.3g}"
        f"{pyo.units.get_units(Co_permeate)}"
    )
    Ca_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Ca_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ca2+ permeate mass concentration: "
        f"{pyo.value(Ca_permeate):.3g}"
        f"{pyo.units.get_units(Ca_permeate)}"
    )
    Cu_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Cu_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Cu2+ permeate mass concentration: "
        f"{pyo.value(Cu_permeate):.3g}"
        f"{pyo.units.get_units(Cu_permeate)}"
    )
    Fe_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Fe_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Fe3+ permeate mass concentration: "
        f"{pyo.value(Fe_permeate):.3g}"
        f"{pyo.units.get_units(Fe_permeate)}"
    )
    Nd_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Nd_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Nd3+ permeate mass concentration: "
        f"{pyo.value(Nd_permeate):.3g}"
        f"{pyo.units.get_units(Nd_permeate)}"
    )
    Ni_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Ni_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ni2+ permeate mass concentration: "
        f"{pyo.value(Ni_permeate):.3g}"
        f"{pyo.units.get_units(Ni_permeate)}"
    )
    Pr_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Pr_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Pr3+ permeate mass concentration: "
        f"{pyo.value(Pr_permeate):.3g}"
        f"{pyo.units.get_units(Pr_permeate)}"
    )
    Na_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Na_+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Na+ permeate mass concentration: "
        f"{pyo.value(Na_permeate):.3g}"
        f"{pyo.units.get_units(Na_permeate)}"
    )
    Cr_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Cr_6+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Cr6+ permeate mass concentration: "
        f"{pyo.value(Cr_permeate):.3g}"
        f"{pyo.units.get_units(Cr_permeate)}"
    )
    Sn_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Sn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Sn2+ permeate mass concentration: "
        f"{pyo.value(Sn_permeate):.3g}"
        f"{pyo.units.get_units(Sn_permeate)}"
    )
    Zn_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Zn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Zn2+ permeate mass concentration: "
        f"{pyo.value(Zn_permeate):.3g}"
        f"{pyo.units.get_units(Zn_permeate)}"
    )
    # Pb_permeate = pyo.units.convert(
    #     m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Pb_2+"],
    #     to_units=pyo.units.mg / pyo.units.L,
    # )
    # print(
    #     f"Pb2+ permeate mass concentration: "
    #     f"{pyo.value(Pb_permeate):.3g}"
    #     f"{pyo.units.get_units(Pb_permeate)}"
    # )
    Dy_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Dy_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Dy3+ permeate mass concentration: "
        f"{pyo.value(Dy_permeate):.3g}"
        f"{pyo.units.get_units(Dy_permeate)}"
    )
    B_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "B_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"B3+ permeate mass concentration: "
        f"{pyo.value(B_permeate):.3g}"
        f"{pyo.units.get_units(B_permeate)}"
    )
    Gd_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Gd_3+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Gd3+ permeate mass concentration: "
        f"{pyo.value(Gd_permeate):.3g}"
        f"{pyo.units.get_units(Gd_permeate)}"
    )
    Mn_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Mn_2+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mn2+ permeate mass concentration: "
        f"{pyo.value(Mn_permeate):.3g}"
        f"{pyo.units.get_units(Mn_permeate)}"
    )
    Si_permeate = pyo.units.convert(
        m.fs.NF.properties_permeate[0].conc_mass_phase_comp["Liq", "Si_4+"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Si4+ permeate mass concentration: "
        f"{pyo.value(Si_permeate):.3g}"
        f"{pyo.units.get_units(Si_permeate)}"
    )

    print("\n---- System Performance Metrics ----")
    f_in = pyo.units.convert(
        m.fs.NF.feed_side.properties_in[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )
    # print(f"Influent flow: " f"{pyo.value(f_in):.3g}" f"{pyo.units.get_units(f_in)}")
    # f_permeate = pyo.units.convert(
    #     m.fs.NF.properties_permeate[0].flow_vol,
    #     to_units=pyo.units.m**3 / pyo.units.hr,
    # )
    # print(
    #     f"Permeate flow: "
    #     f"{pyo.value(f_permeate):.3g}"
    #     f"{pyo.units.get_units(f_permeate)}"
    # )
    # f_retentate = pyo.units.convert(
    #     m.fs.NF.properties_retentate[0].flow_vol,
    #     to_units=pyo.units.m**3 / pyo.units.hr,
    # )
    # print(
    #     f"Retentate flow: "
    #     f"{pyo.value(f_retentate):.3g}"
    #     f"{pyo.units.get_units(f_retentate)}"
    # )
    water_recovery = m.fs.NF.recovery_vol_phase[0, "Liq"]
    print(f"Volumetric-based recovery: " f"{pyo.value(water_recovery):.3g}")
    Co_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Co_2+"]
    print(f"Co2+ mass rejection: " f"{pyo.value(1 - Co_recovery):.3g}")
    Ca_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Ca_2+"]
    print(f"Ca2+ mass rejection: " f"{pyo.value(1 - Ca_recovery):.3g}")
    Cu_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Cu_2+"]
    print(f"Cu2+ mass rejection: " f"{pyo.value(1 - Cu_recovery):.3g}")
    Fe_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Fe_3+"]
    print(f"Fe3+ mass rejection: " f"{pyo.value(1 - Fe_recovery):.3g}")
    Nd_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Nd_3+"]
    print(f"Nd3+ mass rejection: " f"{pyo.value(1 - Nd_recovery):.3g}")
    Ni_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Ni_2+"]
    print(f"Ni2+ mass rejection: " f"{pyo.value(1 - Ni_recovery):.3g}")
    Pr_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Pr_3+"]
    print(f"Pr3+ mass rejection: " f"{pyo.value(1 - Pr_recovery):.3g}")
    Na_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Na_+"]
    print(f"Na+ mass rejection: " f"{pyo.value(1 - Na_recovery):.3g}")
    Cr_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Cr_6+"]
    print(f"Cr6+ mass rejection: " f"{pyo.value(1 - Cr_recovery):.3g}")
    Sn_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Sn_2+"]
    print(f"Sn2+ mass rejection: " f"{pyo.value(1 - Sn_recovery):.3g}")
    Zn_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Zn_2+"]
    print(f"Zn2+ mass rejection: " f"{pyo.value(1 - Zn_recovery):.3g}")
    # Pb_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Pb_2+"]
    # print(f"Pb2+ mass rejection: " f"{pyo.value(1 - Pb_recovery):.3g}")
    Dy_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Dy_3+"]
    print(f"Dy3+ mass rejection: " f"{pyo.value(1 - Dy_recovery):.3g}")
    B_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "B_3+"]
    print(f"B3+ mass rejection: " f"{pyo.value(1 - B_recovery):.3g}")
    Gd_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Gd_3+"]
    print(f"Gd3+ mass rejection: " f"{pyo.value(1 - Gd_recovery):.3g}")
    Mn_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Mn_2+"]
    print(f"Mn3+ mass rejection: " f"{pyo.value(1 - Mn_recovery):.3g}")
    Si_recovery = m.fs.NF.recovery_mass_phase_comp[0, "Liq", "Si_4+"]
    print(f"Si4+ mass rejection: " f"{pyo.value(1 - Si_recovery):.3g}")

    area = pyo.units.convert(m.fs.NF.area, to_units=pyo.units.m**2)
    print(f"Area: " f"{pyo.value(area):.3g}" f"{pyo.units.get_units(area)}")

    flux_vol_solvent = pyo.units.convert(
        m.fs.NF.flux_vol_solvent[0, "H2O"],
        to_units=pyo.units.L / (pyo.units.m**2 * pyo.units.hr),
    )
    print(
        f"Solvent volumetric flux: "
        f"{pyo.value(flux_vol_solvent):.3g}"
        f"{pyo.units.get_units(flux_vol_solvent)}"
    )

    deltaP = pyo.units.convert(m.fs.NF.deltaP[0], to_units=pyo.units.bar)
    print(
        f"Pressure drop: " f"{pyo.value(deltaP):.3g}" f"{pyo.units.get_units(deltaP)}"
    )


def display_costing(m):
    print("\n---- Cost Metrics ----")
    print("Levelized cost of water: %.3g $/m3" % pyo.value(m.fs.costing.LCOW))

    print(
        "Total operating cost: %.3g $/yr" % pyo.value(m.fs.costing.total_operating_cost)
    )
    print("Total capital cost: %.3g $" % pyo.value(m.fs.costing.total_capital_cost))

    print(
        "Total annualized cost: %.3g $/yr"
        % pyo.value(m.fs.costing.total_annualized_cost)
    )
    print("Capital cost pump: %.3g $" % pyo.value(m.fs.P1.costing.capital_cost))
    print(
        "Capital cost nanofiltration: %.3g $" % pyo.value(m.fs.NF.costing.capital_cost)
    )


if __name__ == "__main__":
    m, results = CMR_nf_case(case=Case.case1, simplified_routine=False)
