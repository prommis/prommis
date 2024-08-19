#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diafiltration cascade flowsheet for separating lithium and cobalt.
"""

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Set,
    TransformationFactory,
    Var,
    log,
    units,
    value,
)
from pyomo.network import Arc

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    UnitModelBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import (
    Mixer,
    MixerInitializer,
    MixingType,
    MomentumMixingType,
    MSContactor,
    MSContactorInitializer,
)

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCosting,
    DiafiltrationCostingData,
)
from prommis.nanofiltration.diafiltration_properties import LiCoParameters

# Global constants
Jw = 0.1 * units.m / units.hour
w = 1.5 * units.m
diafiltrate_inlet_pressure = 101325 * units.Pa
diafiltrate_outlet_pressure = 145 * units.psi
Q_feed = 100 * units.m**3 / units.h
C_Li_feed = 1.7 * units.kg / units.m**3
C_Co_feed = 17 * units.kg / units.m**3
Q_diaf = 30 * units.m**3 / units.h
C_Li_diaf = 0.1 * units.kg / units.m**3
C_Co_diaf = 0.2 * units.kg / units.m**3


def main():
    m = build_and_init_model()
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building model"
        )

    # dt = DiagnosticsToolbox(m)
    # dt.report_structural_issues()
    # dt.display_components_with_inconsistent_units()
    # dt.display_potential_evaluation_errors()
    # dt.report_numerical_issues()

    add_costing(m)
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building cost block"
        )

    # dt.report_structural_issues()
    # dt.display_components_with_inconsistent_units()
    # dt.display_potential_evaluation_errors()
    # dt.report_numerical_issues()

    unfix_variables(m)
    add_constraints(m)
    add_objective(m)
    print(f"DOF = {degrees_of_freedom(m)}")
    solve_model(m)
    print_information(m)


def build_and_init_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = LiCoParameters()
    m.fs.stage1, m.fs.stage2, m.fs.stage3 = build_stages(m)
    m.fs.mix1, m.fs.mix2 = build_mixers(m)
    add_streams(m)
    check_model(m)
    add_transport_constraints(m)
    initialize(m)
    assert degrees_of_freedom(m) == 0
    add_expressions(m)
    return m


def build_stages(m):
    m.fs.stage1 = MSContactor(
        number_of_finite_elements=10,
        streams={
            "retentate": {
                "property_package": m.fs.properties,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            "permeate": {
                "property_package": m.fs.properties,
                "has_feed": False,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        },
    )
    m.fs.stage1.length = Var(units=units.m)
    m.fs.stage1.length.fix(10)
    m.fs.stage1.retentate_inlet.flow_vol[0].fix(Q_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].fix(C_Co_feed)

    m.fs.stage2 = MSContactor(
        number_of_finite_elements=10,
        streams={
            "retentate": {
                "property_package": m.fs.properties,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            "permeate": {
                "property_package": m.fs.properties,
                "has_feed": False,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        },
    )
    m.fs.stage2.length = Var(units=units.m)
    m.fs.stage2.length.fix(10)

    m.fs.stage3 = MSContactor(
        number_of_finite_elements=10,
        streams={
            "retentate": {
                "property_package": m.fs.properties,
                "side_streams": [10],
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            "permeate": {
                "property_package": m.fs.properties,
                "has_feed": False,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        },
    )
    m.fs.stage3.length = Var(units=units.m)
    m.fs.stage3.length.fix(10)

    return m.fs.stage1, m.fs.stage2, m.fs.stage3


def build_mixers(m):
    m.fs.mix1 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )
    m.fs.mix2 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )
    return m.fs.mix1, m.fs.mix2


def add_streams(m):
    m.fs.stream1 = Arc(
        source=m.fs.stage1.permeate_outlet,
        destination=m.fs.mix1.inlet_2,
    )
    m.fs.stream2 = Arc(
        source=m.fs.mix1.outlet,
        destination=m.fs.stage2.retentate_inlet,
    )
    m.fs.stream3 = Arc(
        source=m.fs.stage2.retentate_outlet,
        destination=m.fs.stage1.retentate_inlet,
    )
    m.fs.stage1.retentate_inlet.flow_vol[0].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

    m.fs.mix1.inlet_1.flow_vol[0].fix(Q_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].fix(C_Co_feed)

    m.fs.stream4 = Arc(
        source=m.fs.stage2.permeate_outlet,
        destination=m.fs.mix2.inlet_2,
    )
    m.fs.stream5 = Arc(
        source=m.fs.mix2.outlet,
        destination=m.fs.stage3.retentate_inlet,
    )
    m.fs.stream6 = Arc(
        source=m.fs.stage3.retentate_outlet,
        destination=m.fs.mix1.inlet_1,
    )
    # Unfix Feed to Mixer 1 inlet 1
    m.fs.mix1.inlet_1.flow_vol[0].unfix()
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].unfix()
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].unfix()

    # Fix Diafiltrate feed to Mixer 2 inlet 1
    m.fs.mix2.inlet_1.flow_vol[0].fix(30)
    m.fs.mix2.inlet_1.conc_mass_solute[0, "Li"].fix(0.1)
    m.fs.mix2.inlet_1.conc_mass_solute[0, "Co"].fix(0.2)

    # Fix Feed stream at element 10 of stage 3
    m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol.fix(Q_feed)
    m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"].fix(C_Li_feed)
    m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"].fix(C_Co_feed)

    TransformationFactory("network.expand_arcs").apply_to(m)


# TODO: move these checks to test file
def check_model(m):
    assert isinstance(
        m.fs.stage3, MSContactor
    )  # check that stage3 exists and is an MSContactor

    # Retentate side checks
    assert hasattr(
        m.fs.stage3, "retentate_inlet_state"
    )  # check that there is a retentate feed
    assert hasattr(
        m.fs.stage3, "retentate_side_stream_state"
    )  # check that a side stream exists
    for k in m.fs.stage3.retentate_side_stream_state:
        assert k == (0, 10)  # check that the side stream only exists at element 10
    assert not hasattr(
        m.fs.stage3, "retentate_energy_balance"
    )  # check that there are no energy balances
    assert not hasattr(
        m.fs.stage3, "retentate_pressure_balance"
    )  # check that there are no pressure balances

    # Permeate side checks
    assert not hasattr(
        m.fs.stage3, "permeate_inlet_state"
    )  # check that there is no permeate feed
    assert not hasattr(
        m.fs.stage3, "permeate_side_stream_state"
    )  # check that there are no side streams on permeate side
    assert not hasattr(
        m.fs.stage3, "permeate_energy_balance"
    )  # check that there are no energy balances
    assert not hasattr(
        m.fs.stage3, "permeate_pressure_balance"
    )  # check that there are no pressure balances


def initialize(m):
    initializer = MSContactorInitializer()
    initializer.initialize(m.fs.stage1)

    propagate_state(
        destination=m.fs.mix1.inlet_2,
        source=m.fs.stage1.permeate_outlet,
    )
    mix_initializer = MixerInitializer()
    mix_initializer.initialize(m.fs.mix1)

    propagate_state(
        destination=m.fs.stage2.retentate_inlet,
        source=m.fs.mix1.outlet,
    )
    initializer.initialize(m.fs.stage2)

    propagate_state(
        destination=m.fs.mix2.inlet_2,
        source=m.fs.stage2.permeate_outlet,
    )
    mix_initializer.initialize(m.fs.mix2)

    propagate_state(
        destination=m.fs.stage3.retentate_inlet,
        source=m.fs.mix2.outlet,
    )
    initializer.initialize(m.fs.stage3)

    # Lets try again
    # Initial feed guess is pure diafiltrate (no recycle)
    # Remember to only set the value, not fix the inlet
    m.fs.stage3.retentate_inlet.flow_vol[0].set_value(Q_diaf)
    m.fs.stage3.retentate_inlet.conc_mass_solute[0, "Li"].set_value(C_Li_diaf)
    m.fs.stage3.retentate_inlet.conc_mass_solute[0, "Co"].set_value(C_Co_diaf)

    initializer.initialize(m.fs.stage3)


def add_transport_constraints(m):
    m.fs.solutes = Set(initialize=["Li", "Co"])

    m.fs.sieving_coefficient = Var(
        m.fs.solutes,
        units=units.dimensionless,
    )
    m.fs.sieving_coefficient["Li"].fix(1.3)
    m.fs.sieving_coefficient["Co"].fix(0.5)

    # TODO: update to @blk.Constraint() format
    def solvent_rule(b, s):
        return (
            b.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == Jw * b.length * w * m.fs.properties.dens_H2O / 10
        )

    def solute_rule(b, s, j):
        if s == 1:
            in_state = b.retentate_inlet_state[0]
        else:
            sp = b.elements.prev(s)
            in_state = b.retentate[0, sp]

        return log(b.retentate[0, s].conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * log(in_state.flow_vol) == log(in_state.conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * log(
            b.retentate[0, s].flow_vol
        )

    def stage3_solute_rule(b, s, j):
        if s == 1:
            q_in = b.retentate_inlet_state[0].flow_vol
            c_in = b.retentate_inlet_state[0].conc_mass_solute[j]
        elif s == 10:
            sp = b.elements.prev(s)
            q_in = (
                b.retentate[0, sp].flow_vol
                + b.retentate_side_stream_state[0, 10].flow_vol
            )
            c_in = (
                b.retentate[0, sp].conc_mass_solute[j] * b.retentate[0, sp].flow_vol
                + b.retentate_side_stream_state[0, 10].conc_mass_solute[j]
                * b.retentate_side_stream_state[0, 10].flow_vol
            ) / q_in
        else:
            sp = b.elements.prev(s)
            q_in = b.retentate[0, sp].flow_vol
            c_in = b.retentate[0, sp].conc_mass_solute[j]

        return log(b.retentate[0, s].conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * log(q_in) == log(c_in) + (m.fs.sieving_coefficient[j] - 1) * log(
            b.retentate[0, s].flow_vol
        )

    m.fs.stage1.solvent_flux = Constraint(
        m.fs.stage1.elements,
        rule=solvent_rule,
    )
    m.fs.stage1.solute_sieving = Constraint(
        m.fs.stage1.elements,
        m.fs.solutes,
        rule=solute_rule,
    )
    m.fs.stage2.solvent_flux = Constraint(
        m.fs.stage2.elements,
        rule=solvent_rule,
    )
    m.fs.stage2.solute_sieving = Constraint(
        m.fs.stage2.elements,
        m.fs.solutes,
        rule=solute_rule,
    )
    m.fs.stage3.solvent_flux = Constraint(
        m.fs.stage3.elements,
        rule=solvent_rule,
    )
    m.fs.stage3.solute_sieving = Constraint(
        m.fs.stage3.elements,
        m.fs.solutes,
        rule=stage3_solute_rule,
    )


def solve_model(m):
    solver = get_solver()
    solver.solve(m, tee=True)


def add_expressions(m):
    m.Li_recovery = Expression(
        expr=m.fs.stage3.permeate_outlet.flow_vol[0]
        * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
        / (
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Li"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"]
        )
    )
    m.Li_purity = Expression(
        expr=m.fs.stage3.permeate_outlet.flow_vol[0]
        * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
        / (
            (
                m.fs.stage3.permeate_outlet.flow_vol[0]
                * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
            )
            + (
                m.fs.stage3.permeate_outlet.flow_vol[0]
                * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Co"]
            )
        )
    )
    m.Co_recovery = Expression(
        expr=m.fs.stage1.retentate_outlet.flow_vol[0]
        * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
        / (
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Co"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"]
        )
    )
    m.Co_purity = Expression(
        expr=m.fs.stage1.retentate_outlet.flow_vol[0]
        * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
        / (
            (
                m.fs.stage1.retentate_outlet.flow_vol[0]
                * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
            )
            + (
                m.fs.stage1.retentate_outlet.flow_vol[0]
                * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Li"]
            )
        )
    )


def get_expression_values(m):
    R_Li = round(value(m.Li_recovery) * 100, 2)
    R_Co = round(value(m.Co_recovery) * 100, 2)

    P_Li = round(value(m.Li_purity) * 100, 2)
    P_Co = round(value(m.Co_purity) * 100, 2)

    return R_Li, R_Co, P_Li, P_Co


def get_model_values(m):
    membrane_length = (
        value(m.fs.stage1.length)
        + value(m.fs.stage2.length)
        + value(m.fs.stage3.length)
    ) * units.m

    # the feed is the stream entering stage 3 in the last element
    feed_flow_rate = (
        value(m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol)
        * units.m**3
        / units.hr
    )
    # the diafiltrate is the stream entering stage 3 in the first element
    diafiltrate_flow_rate = (
        value(m.fs.stage3.retentate_inlet.flow_vol[0]) * units.m**3 / units.hr
    )
    # the permeate is the stream exiting stage 3 in the last element
    permeate_flow_rate = (
        value(m.fs.stage3.permeate_outlet.flow_vol[0]) * units.m**3 / units.hr
    )
    return membrane_length, feed_flow_rate, diafiltrate_flow_rate, permeate_flow_rate


def add_costing(
    m,
    membrane_width=w,
    water_flux=Jw,
    inlet_pressure=diafiltrate_inlet_pressure,
    outlet_pressure=diafiltrate_outlet_pressure,
):
    (
        membrane_length,
        feed_flow_rate,
        diafiltrate_flow_rate,
        permeate_flow_rate,
    ) = get_model_values(m)

    # creating dummy variables to store the UnitModelCostingBlocks
    m.fs.membrane = UnitModelBlock()
    m.fs.pump = UnitModelBlock()

    m.fs.costing = DiafiltrationCosting()
    m.fs.membrane.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": membrane_length,
            "membrane_width": membrane_width,
            "water_flux": water_flux,
            "vol_flow_feed": feed_flow_rate,
            "vol_flow_perm": permeate_flow_rate,
        },
    )
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": inlet_pressure,
            "outlet_pressure": outlet_pressure,
            "inlet_vol_flow": diafiltrate_flow_rate,
        },
    )
    m.fs.costing.cost_process()


def unfix_variables(m):
    m.fs.stage1.length.unfix()
    m.fs.stage2.length.unfix()
    m.fs.stage3.length.unfix()


def add_constraints(m):
    m.Co_recovery_constraint = Constraint(expr=m.Co_recovery >= 0.4)
    m.Li_recovery_constraint = Constraint(expr=m.Li_recovery >= 0.95)

    # m.Li_purity_constraint = pyo.Constraint(expr=m.Li_purity >= 0.8)
    # m.Co_purity_constraint = pyo.Constraint(expr=m.Co_purity >= 0.5)


def add_objective(m):
    def cost_obj(m):
        return m.fs.costing.total_annualized_cost

    m.cost_objecticve = Objective(rule=cost_obj)


# TODO: update functionality to specify ann output stream
def print_information(m):
    R_Li, R_Co, P_Li, P_Co = get_expression_values(m)

    print(f"The lithium recovery is {R_Li}% at purity {P_Li}")
    print(f"The cobalt recovery is {R_Co}% at purity {P_Co}")

    print("\nmembrane area")
    print(f"{value(m.fs.membrane.costing.membrane_area)} m2")

    print("\nmembrane capital cost")
    print(f"${value(m.fs.membrane.costing.capital_cost)}")

    print("\npump capital cost")
    print(f"${value(m.fs.pump.costing.capital_cost)}")

    print("\naggregate capital cost")
    print(f"${value(m.fs.costing.aggregate_capital_cost)}")

    print("\ntotal capital cost")
    print(f"${value(m.fs.costing.total_capital_cost)}")

    print("\ntotal operating cost")
    print(f"${value(m.fs.costing.total_operating_cost)}")

    print("\ntotal annualized cost")
    print(f"${value(m.fs.costing.total_annualized_cost)}")


if __name__ == "__main__":
    main()
