#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diafiltration cascade flowsheet case study example for separating lithium and cobalt. This flowsheet
models a specific system from literature (cascade III in Figure 2 [1]) and serves as an example of
implementing a custom cost model.

[1] https://pubs.acs.org/doi/full/10.1021/acssuschemeng.2c02862
"""

from pyomo.environ import (
    ConcreteModel,
    Expression,
    Objective,
    Param,
    Set,
    Suffix,
    TransformationFactory,
    Var,
    log,
    units,
    value,
)
from pyomo.network import Arc

import idaes.logger as idaeslog
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

_log = idaeslog.getLogger(__name__)


def main():
    """
    Builds and solves the diafiltration flowsheet with cost

    Returns:
        m: Pyomo model
    """
    m = build_model()
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building model"
        )

    add_costing(m)
    if degrees_of_freedom(m) != 0:
        raise ValueError(
            "Degrees of freedom were not equal to zero after building cost block"
        )

    initialize_model(m)
    _log.info("Initialization Okay")

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    solve_model(m)
    _log.info("Solved Square Problem")

    dt.assert_no_numerical_warnings()

    unfix_opt_variables(m)
    add_product_constraints(m, Li_recovery_bound=0.95, Co_recovery_bound=0.635)
    add_objective(m)

    # Create a scaled version of the model to solve
    set_scaling(m)
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)
    solve_model(scaled_model)
    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, m)

    # TODO: add Boolean variable to calculate pump OPEX
    # Verify the feed pump operating pressure workaround is valid
    # assume this additional cost is less than half a cent
    if value(m.fs.feed_pump.costing.variable_operating_cost) >= 0.005:
        raise ValueError(
            "The variable  operating cost of the feed pump as calculated in the feed"
            "pump costing block is not negligible. This operating cost is already"
            "accounted for via the membrane's pressure drop specific energy consumption."
        )

    dt.assert_no_numerical_warnings()
    print_information(m)

    return m


def build_model():
    """
    Builds the diafiltration flowsheet

    Returns:
        m: Pyomo model
    """
    m = ConcreteModel()
    # These parameters will change for different membranes or cascade configurations
    add_global_flowsheet_parameters(m)

    m.fs = FlowsheetBlock(dynamic=False)

    # This membrane systems uses a constant sieving coefficient model
    m.fs.properties = LiCoParameters()
    m.fs.solutes = Set(initialize=["Li", "Co"])
    m.fs.sieving_coefficient = Var(
        m.fs.solutes,
        units=units.dimensionless,
    )
    m.fs.sieving_coefficient["Li"].fix(1.3)
    m.fs.sieving_coefficient["Co"].fix(0.5)

    m.fs.stage1, m.fs.stage2, m.fs.stage3, m.fs.mix1, m.fs.mix2 = build_unit_models(m)

    # Connect the unit models
    add_streams(m)
    # Set the desired values for this system
    fix_stream_values(m)
    # Add recovery, purity, and membrane length expressions
    add_useful_expressions(m)

    return m


def add_global_flowsheet_parameters(m):
    """
    Adds global parameters to the Pyomo model. These values are system-dependent.

    Args:
        m: Pyomo model
    """
    m.Jw = Param(
        initialize=0.1,
        mutable=True,
        doc="Water flux",  # assumed constant in this model
        units=units.m**3 / units.m**2 / units.h,
    )
    m.w = Param(
        initialize=1.5,
        mutable=True,
        doc="Membrane module length",
        units=units.m,
    )
    m.atmospheric_pressure = Param(
        initialize=101325,
        doc="Atmospheric pressure in Pascal",
        units=units.Pa,
    )
    m.operating_pressure = Param(
        initialize=145,
        mutable=True,
        doc="Membrane operating pressure",
        units=units.psi,
    )
    m.Q_feed = Param(
        initialize=100,
        mutable=True,
        doc="Cascade feed flow",
        units=units.m**3 / units.h,
    )
    m.C_Li_feed = Param(
        initialize=1.7,
        mutable=True,
        doc="Lithium concentration in cascade feed",
        units=units.kg / units.m**3,
    )
    m.C_Co_feed = Param(
        initialize=17,
        mutable=True,
        doc="Cobalt concentration in cascade feed",
        units=units.kg / units.m**3,
    )
    m.Q_diaf = Param(
        initialize=30,
        mutable=True,
        doc="Cascade diafiltrate flow",
        units=units.m**3 / units.h,
    )
    m.C_Li_diaf = Param(
        initialize=0.1,
        mutable=True,
        doc="Lithium concentration in diafiltrate",
        units=units.kg / units.m**3,
    )
    m.C_Co_diaf = Param(
        initialize=0.2,
        mutable=True,
        doc="Cobalt concentration in diafiltrate",
        units=units.kg / units.m**3,
    )


def build_unit_models(m):
    """
    Adds the membrane stages and mixers to the flowsheet to build the cascade of interest

    Args:
        m: Pyomo model

    Returns:
        m.fs.stage1: first membrane stage
        m.fs.stage2: second membrane stage
        m.fs.stage3: third membrane stage
        m.fs.mix1: mixer for stage 1 permeate and recycle into stage 2
        m.fs.mix2: mixer for stage 2 permeate and diafiltrate into stage 3
    """
    m.fs.stage1 = MSContactor(
        number_of_finite_elements=10,  # assume 10 tubes per stage
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
        number_of_finite_elements=10,  # assume 10 tubes per stage
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
        number_of_finite_elements=10,  # assume 10 tubes per stage
        streams={
            "retentate": {
                "property_package": m.fs.properties,
                "side_streams": [
                    10
                ],  # the feed enters in the last element of this stage
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

    # Add a mixer for the stage 1 permeate and recycle into stage 2
    m.fs.mix1 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )
    # Add a mixer for the stage 2 permeate and diafiltrate into stage 3
    m.fs.mix2 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )
    return m.fs.stage1, m.fs.stage2, m.fs.stage3, m.fs.mix1, m.fs.mix2


def add_stage1_constraints(blk, model):
    """
    Adds the constraints for the sieving coefficient model to stage 1

    Args:
        blk: the MSContactor block (stage 1)
        model: the Pyomo model (m)
    """

    @blk.Constraint(blk.elements)
    def stage_1_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == model.Jw * blk.length * model.w * model.fs.properties.dens_H2O / 10
        )

    @blk.Constraint(blk.elements, model.fs.solutes)
    def stage_1_solute_sieving(blk, s, j):
        if s == 1:
            in_state = blk.retentate_inlet_state[0]
        else:
            sp = blk.elements.prev(s)
            in_state = blk.retentate[0, sp]

        return log(
            blk.retentate[0, s].conc_mass_solute[j] * units.m**3 / units.kg
        ) + (model.fs.sieving_coefficient[j] - 1) * log(
            in_state.flow_vol * units.h / units.m**3
        ) == log(
            in_state.conc_mass_solute[j] * units.m**3 / units.kg
        ) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(
            blk.retentate[0, s].flow_vol * units.h / units.m**3
        )


def add_stage2_constraints(blk, model):
    """
    Adds the constraints for the sieving coefficient model to stage 2

    Args:
        blk: the MSContactor block (stage 2)
        model: the Pyomo model (m)
    """

    @blk.Constraint(blk.elements)
    def stage_2_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == model.Jw * blk.length * model.w * model.fs.properties.dens_H2O / 10
        )

    @blk.Constraint(blk.elements, model.fs.solutes)
    def stage_2_solute_sieving(blk, s, j):
        if s == 1:
            in_state = blk.retentate_inlet_state[0]
        else:
            sp = blk.elements.prev(s)
            in_state = blk.retentate[0, sp]

        return log(
            blk.retentate[0, s].conc_mass_solute[j] * units.m**3 / units.kg
        ) + (model.fs.sieving_coefficient[j] - 1) * log(
            in_state.flow_vol * units.h / units.m**3
        ) == log(
            in_state.conc_mass_solute[j] * units.m**3 / units.kg
        ) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(
            blk.retentate[0, s].flow_vol * units.h / units.m**3
        )


def add_stage3_constraints(blk, model):
    """
    Adds the constraints for the sieving coefficient model to stage 3

    Args:
        blk: the MSContactor block (stage 3)
        model: the Pyomo model (m)
    """

    @blk.Constraint(blk.elements)
    def stage_3_solvent_flux(blk, s):
        return (
            blk.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
            == model.Jw * blk.length * model.w * model.fs.properties.dens_H2O / 10
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

        return log(
            blk.retentate[0, s].conc_mass_solute[j] * units.m**3 / units.kg
        ) + (model.fs.sieving_coefficient[j] - 1) * log(
            q_in * units.h / units.m**3
        ) == log(
            c_in * units.m**3 / units.kg
        ) + (
            model.fs.sieving_coefficient[j] - 1
        ) * log(
            blk.retentate[0, s].flow_vol * units.h / units.m**3
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


def fix_stream_values(m):
    """
    Fixes the volumetric flow and concentration of streams

    Args:
        m: Pyomo model
    """
    m.fs.stage1.retentate_inlet.flow_vol[0].fix(m.Q_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].fix(m.C_Li_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].fix(m.C_Co_feed)

    m.fs.stage1.retentate_inlet.flow_vol[0].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

    m.fs.mix1.inlet_1.flow_vol[0].fix(m.Q_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].fix(m.C_Li_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].fix(m.C_Co_feed)

    # Unfix Feed to Mixer 1 inlet 1
    m.fs.mix1.inlet_1.flow_vol[0].unfix()
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].unfix()
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].unfix()

    # Fix Diafiltrate feed to Mixer 2 inlet 1
    m.fs.mix2.inlet_1.flow_vol[0].fix(30)
    m.fs.mix2.inlet_1.conc_mass_solute[0, "Li"].fix(0.1)
    m.fs.mix2.inlet_1.conc_mass_solute[0, "Co"].fix(0.2)

    # Fix Feed stream at element 10 of stage 3
    m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol.fix(m.Q_feed)
    m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"].fix(
        m.C_Li_feed
    )
    m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"].fix(
        m.C_Co_feed
    )


def initialize_model(m):
    """
    Method to initialize the diafiltration flowsheet

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
    m.fs.stage3.retentate_inlet.flow_vol[0].set_value(m.Q_diaf)
    m.fs.stage3.retentate_inlet.conc_mass_solute[0, "Li"].set_value(m.C_Li_diaf)
    m.fs.stage3.retentate_inlet.conc_mass_solute[0, "Co"].set_value(m.C_Co_diaf)

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
    solver.options = {"max_iter": 3000}
    results = solver.solve(m, tee=True)
    if results.solver.termination_condition != "optimal":
        raise ValueError("The solver did not return optimal termination")

    return results


def add_costing(m):
    """
    Method to add costing block to the flowsheet

    Args:
        m: Pyomo model
    """
    # Create dummy variables to store the UnitModelCostingBlocks
    # These are needed because the sieving coefficient model does not account for pressure
    m.fs.cascade = UnitModelBlock()  # to cost the pressure drop
    m.fs.feed_pump = UnitModelBlock()  # to cost feed pump
    m.fs.diafiltrate_pump = UnitModelBlock()  # to cost diafiltrate pump

    m.fs.costing = DiafiltrationCosting()
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
            "inlet_pressure": m.atmospheric_pressure
            + units.convert(m.fs.cascade.costing.pressure_drop, to_units=units.Pa),
            "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
            * units.psi,  # this should make m.fs.feed_pump.costing.fixed_operating_cost ~0
            "inlet_vol_flow": m.fs.stage3.retentate_side_stream_state[
                0, 10
            ].flow_vol,  # feed
        },
    )
    m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,
            "outlet_pressure": m.operating_pressure,
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
    # For this example we will optimize over the membrane area
    m.fs.stage1.length.unfix()
    m.fs.stage2.length.unfix()
    m.fs.stage3.length.unfix()


def add_product_constraints(
    m,
    Li_recovery_bound,
    Co_recovery_bound,
    recovery=True,
    Li_purity_bound=None,
    Co_purity_bound=None,
    purity=False,
):
    """
    Method to add recovery/purity constraints to the flowsheet for performing optimization

    Args:
        m: Pyomo model
    """
    if recovery:

        @m.Constraint()
        def Li_recovery_constraint(m):
            return m.Li_recovery >= Li_recovery_bound

        @m.Constraint()
        def Co_recovery_constraint(m):
            return m.Co_recovery >= Co_recovery_bound

    # By default, these constraints are not added for this example
    # Toggle the Boolean and add bound arguments to apply to flowsheet
    if purity:
        if Li_purity_bound == None:
            raise ValueError("A lithium product purity bound was not provided")

        if Co_purity_bound == None:
            raise ValueError("A cobalt product purity bound was not provided")

        @m.Constraint()
        def Li_purity_constraint(m):
            return m.Li_purity >= Li_purity_bound

        @m.Constraint()
        def Co_purity_constraint(m):
            return m.Co_purity >= Co_purity_bound


def add_objective(m):
    """
    Method to add cost objective to flowsheet for performing optimization

    Args:
        m: Pyomo model
    """

    def cost_obj(m):
        return m.fs.costing.total_annualized_cost

    m.cost_objecticve = Objective(rule=cost_obj)


def set_scaling(m):
    """
    Apply scaling factors to certain constraints to improve solver performance

    Args:
        m: Pyomo model
    """
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Add scaling factors for poorly scaled constraints
    m.scaling_factor[m.fs.cascade.costing.variable_operating_cost_constraint] = 1e-4
    m.scaling_factor[m.fs.feed_pump.costing.capital_cost_constraint] = 1e-5
    m.scaling_factor[m.fs.diafiltrate_pump.costing.capital_cost_constraint] = 1e-4
    m.scaling_factor[m.fs.feed_pump.costing.pump_head_equation] = 1e6
    m.scaling_factor[m.fs.feed_pump.costing.pump_power_equation] = 1e9
    m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost_constraint] = 1e8
    m.scaling_factor[m.fs.costing.aggregate_capital_cost_constraint] = 1e-5
    m.scaling_factor[m.fs.costing.aggregate_variable_operating_cost_constraint] = 1e-5
    m.scaling_factor[m.fs.costing.total_capital_cost_constraint] = 1e-5
    m.scaling_factor[
        m.fs.costing.maintenance_labor_chemical_operating_cost_constraint
    ] = 1e-4
    m.scaling_factor[m.fs.costing.total_operating_cost_constraint] = 1e-5

    # Add scaling factors for poorly scaled variables
    m.scaling_factor[m.fs.cascade.costing.variable_operating_cost] = 1e-5
    m.scaling_factor[m.fs.feed_pump.costing.capital_cost] = 1e-5
    m.scaling_factor[m.fs.feed_pump.costing.variable_operating_cost] = 1e8
    m.scaling_factor[m.fs.feed_pump.costing.pump_head] = 1e6
    m.scaling_factor[m.fs.feed_pump.costing.pump_power] = 1e9
    m.scaling_factor[m.fs.diafiltrate_pump.costing.capital_cost] = 1e-4
    m.scaling_factor[m.fs.costing.aggregate_capital_cost] = 1e-5
    m.scaling_factor[m.fs.costing.aggregate_variable_operating_cost] = 1e-5
    m.scaling_factor[m.fs.costing.total_capital_cost] = 1e-5
    m.scaling_factor[m.fs.costing.total_operating_cost] = 1e-5
    m.scaling_factor[m.fs.costing.maintenance_labor_chemical_operating_cost] = 1e-4


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

    print("\ntotal membrane area")
    print(f"{value(m.membrane_length)*value(m.w)} m2")

    print("\ntotal annualized cost")
    print(f"${value(m.fs.costing.total_annualized_cost)}")


if __name__ == "__main__":
    main()
