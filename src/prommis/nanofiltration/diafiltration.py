#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diafiltration cascade flowsheet for separating lithium and cobalt. Serves as an example of
custom cost model.
"""

from pyomo.environ import (
    ConcreteModel,
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
    """
    Builds and solves the diafiltration flowsheet with cost

    Returns:
        m: Pyomo model
    """
    m = build_model()
    initialize_model(m)
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building model"
        )

    dt = DiagnosticsToolbox(m)
    dt.report_structural_issues()

    add_costing(m)
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building cost block"
        )

    dt.report_numerical_issues()

    unfix_opt_variables(m)
    add_product_constraints(m)
    add_objective(m)
    print(f"DOF = {degrees_of_freedom(m)}")
    solve_model(m)
    print_information(m)

    return m


def build_model():
    """
    Builds the diafiltration flowsheet

    Returns:
        m: Pyomo model
    """
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = LiCoParameters()
    m.fs.solutes = Set(initialize=["Li", "Co"])
    m.fs.sieving_coefficient = Var(
        m.fs.solutes,
        units=units.dimensionless,
    )
    m.fs.sieving_coefficient["Li"].fix(1.3)
    m.fs.sieving_coefficient["Co"].fix(0.5)

    m.fs.stage1, m.fs.stage2, m.fs.stage3, m.fs.mix1, m.fs.mix2 = build_unit_models(m)

    add_streams(m)
    fix_values(m)
    add_useful_expressions(m)  # adds recovery, purity, and membrane length expressions

    return m


def build_unit_models(m):
    """
    Adds the membrane stages and mixers to the flowsheet

    Args:
        m: Pyomo model

    Returns:
        m.fs.stage1: first membrane stage
        m.fs.stage2: second membrane stage
        m.fs.stage3: third membrane stage
        m.fs.mix1: mixer for stage 1 permeate and recycle into stage 2
        m.fs.mix2: mixer for stage 2 permeate and diafiltrate ito stage 3
    """
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
    m.fs.stage1.length.fix(10)  # initialize
    add_stage1_constraints(blk=m.fs.stage1, model=m)

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
    m.fs.stage2.length.fix(10)  # initialize
    add_stage2_constraints(blk=m.fs.stage2, model=m)

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
    m.fs.stage3.length.fix(10)  # initialize
    add_stage3_constraints(blk=m.fs.stage3, model=m)

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
    return m.fs.stage1, m.fs.stage2, m.fs.stage3, m.fs.mix1, m.fs.mix2


def add_stage1_constraints(blk, model):
    @blk.Constraint(blk.elements)
    def stage_1_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == Jw * blk.length * w * model.fs.properties.dens_H2O / 10
        )

    @blk.Constraint(blk.elements, model.fs.solutes)
    def stage_1_solute_sieving(blk, s, j):
        if s == 1:
            in_state = blk.retentate_inlet_state[0]
        else:
            sp = blk.elements.prev(s)
            in_state = blk.retentate[0, sp]

        return log(blk.retentate[0, s].conc_mass_solute[j]) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(in_state.flow_vol) == log(in_state.conc_mass_solute[j]) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(
            blk.retentate[0, s].flow_vol
        )


def add_stage2_constraints(blk, model):
    @blk.Constraint(blk.elements)
    def stage_2_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == Jw * blk.length * w * model.fs.properties.dens_H2O / 10
        )

    @blk.Constraint(blk.elements, model.fs.solutes)
    def stage_2_solute_sieving(blk, s, j):
        if s == 1:
            in_state = blk.retentate_inlet_state[0]
        else:
            sp = blk.elements.prev(s)
            in_state = blk.retentate[0, sp]

        return log(blk.retentate[0, s].conc_mass_solute[j]) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(in_state.flow_vol) == log(in_state.conc_mass_solute[j]) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(
            blk.retentate[0, s].flow_vol
        )


def add_stage3_constraints(blk, model):
    @blk.Constraint(blk.elements)
    def stage_3_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == Jw * blk.length * w * model.fs.properties.dens_H2O / 10
        )

    @blk.Constraint(blk.elements, model.fs.solutes)
    def stage_3_solute_sieving(blk, s, j):
        if s == 1:
            q_in = blk.retentate_inlet_state[0].flow_vol
            c_in = blk.retentate_inlet_state[0].conc_mass_solute[j]
        elif s == 10:
            sp = blk.elements.prev(s)
            q_in = (
                blk.retentate[0, sp].flow_vol
                + blk.retentate_side_stream_state[0, 10].flow_vol
            )
            c_in = (
                blk.retentate[0, sp].conc_mass_solute[j] * blk.retentate[0, sp].flow_vol
                + blk.retentate_side_stream_state[0, 10].conc_mass_solute[j]
                * blk.retentate_side_stream_state[0, 10].flow_vol
            ) / q_in
        else:
            sp = blk.elements.prev(s)
            q_in = blk.retentate[0, sp].flow_vol
            c_in = blk.retentate[0, sp].conc_mass_solute[j]

        return log(blk.retentate[0, s].conc_mass_solute[j]) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(q_in) == log(c_in) + (model.fs.sieving_coefficient[j] - 1) * log(
            blk.retentate[0, s].flow_vol
        )


def add_streams(m):
    """
    Defines and connects streams in the flowsheet

    Args:
        m: Pyomo model
    """
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
    TransformationFactory("network.expand_arcs").apply_to(m)


def fix_values(m):
    """
    Fixes the volumetric flow and concentration of streams

    Args:
        m: Pyomo model
    """
    m.fs.stage1.retentate_inlet.flow_vol[0].fix(Q_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].fix(C_Co_feed)

    m.fs.stage1.retentate_inlet.flow_vol[0].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

    m.fs.mix1.inlet_1.flow_vol[0].fix(Q_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].fix(C_Co_feed)

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


def initialize_model(m):
    """
    Method to initialize the diafiltration flowhseet

    Args:
        m: Pyomo model
    """
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


def add_useful_expressions(m):
    """
    Method to add recovery, purity, and membrane length expressions for convenience

    Args:
        m: Pyomo model
    """
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

    m.membrane_length = Expression(
        expr=m.fs.stage1.length + m.fs.stage2.length + m.fs.stage3.length
    )


def solve_model(m):
    """
    Method to solve the diafiltration flowsheet

    Args:
        m: Pyomo model
    """
    solver = get_solver()
    solver.solve(m, tee=True)


def add_costing(
    m,
    membrane_width=w,
    water_flux=Jw,
    inlet_pressure=diafiltrate_inlet_pressure,
    outlet_pressure=diafiltrate_outlet_pressure,
):
    """
    Method to add costing block to the flowsheet

    Args:
        m: Pyomo model
        membrane_width: width of membranes in cascade (tube length) (m)
        water_flux: flux of water across the membrane (m/h)
        inlet_pressure: inlet pressure of the diafiltrate to pump (Pa)
        outlet_pressure: outlet pressure of difiltrate from pump (psi)
    """
    # creating dummy variables to store the UnitModelCostingBlocks
    m.fs.membrane = UnitModelBlock()
    m.fs.pump = UnitModelBlock()

    m.fs.costing = DiafiltrationCosting()
    m.fs.membrane.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.membrane_length,  # total membrane length
            "membrane_width": membrane_width,
            "water_flux": water_flux,
            "vol_flow_feed": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # feed
            "vol_flow_perm": m.fs.stage3.permeate_outlet.flow_vol[0],  # permeate
        },
    )
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": inlet_pressure,
            "outlet_pressure": outlet_pressure,
            "inlet_vol_flow": m.fs.stage3.retentate_inlet.flow_vol[0],  # diafiltrate
        },
    )
    m.fs.costing.cost_process()


def unfix_opt_variables(m):
    """
    Method to unfix variables for performing optimization with DOF>0

    Args:
        m: Pyomo model
    """
    m.fs.stage1.length.unfix()
    m.fs.stage2.length.unfix()
    m.fs.stage3.length.unfix()


def add_product_constraints(m):
    """
    Method to add recovery/purity constraints to the flowsheet for performing optimization

    Args:
        m: Pyomo model
    """

    @m.Constraint()
    def Co_recovery_constraint(m):
        return m.Co_recovery >= 0.4

    @m.Constraint()
    def Li_recovery_constraint(m):
        return m.Li_recovery >= 0.95

    # @m.Constraint()
    # def Li_purity_constraint(m):
    #     return m.Li_purity >= 0.8

    # @m.Constraint()
    # def Co_purity_constraint(m):
    #     return m.Co_purity >= 0.5


def add_objective(m):
    """
    Method to add cost objective to flowsheet for performing optimization

    Args:
        m: Pyomo model
    """

    def cost_obj(m):
        return m.fs.costing.total_annualized_cost

    m.cost_objecticve = Objective(rule=cost_obj)


# TODO: update functionality to specify an output stream
def print_information(m):
    """
    Prints relevant information after solving the model

    Args:
        m: Pyomo model
    """
    print(
        f"The lithium recovery is {value(m.Li_recovery) * 100}% at purity {value(m.Li_purity) * 100}"
    )
    print(
        f"The cobalt recovery is {value(m.Co_recovery) * 100}% at purity {value(m.Co_purity) * 100}"
    )

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
