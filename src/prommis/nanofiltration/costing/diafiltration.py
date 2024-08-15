#########################################################
# MSContactor flowsheet code adapted from Andrew Lee.
#########################################################


import pyomo.environ as pyo
from pyomo.network import Arc

from idaes.core.util.model_diagnostics import DiagnosticsToolbox

from idaes.core import (
    Component,
    declare_process_block_class,
    FlowsheetBlock,
    MaterialBalanceType,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    UnitModelBlock,
    UnitModelCostingBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, propagate_state
from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCosting,
    DiafiltrationCostingData,
)

from idaes.models.unit_models import (
    MSContactor,
    MSContactorInitializer,
    Mixer,
    MixingType,
    MomentumMixingType,
    MixerInitializer,
)


class _StateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        fix_state_vars(self)

    def initialization_routine(self):
        pass


@declare_process_block_class("LiCoStateBlock", block_class=_StateBlock)
class LiCoStateBlock1Data(StateBlockData):
    def build(self):
        super().build()

        self.flow_vol = pyo.Var(
            units=pyo.units.m**3 / pyo.units.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_solute = pyo.Var(
            ["Li", "Co"],
            units=pyo.units.kg / pyo.units.m**3,
            bounds=(1e-8, None),
        )

    def get_material_flow_terms(self, p, j):
        if j == "solvent":
            # Assume constant density of pure water
            return self.flow_vol * self.params.dens_H2O
        else:
            return self.flow_vol * self.conc_mass_solute[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_solute": self.conc_mass_solute,
        }


@declare_process_block_class("LiCoParameters")
class LiCoParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent = Component()
        self.Li = Component()
        self.Co = Component()

        self.dens_H2O = pyo.Param(
            default=1000,
            units=pyo.units.kg / pyo.units.m**3,
        )

        self._state_block_class = LiCoStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyo.units.hour,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.mol,
                "temperature": pyo.units.K,
            }
        )


# Global constants
Jw = 0.1 * pyo.units.m / pyo.units.hour
w = 1.5 * pyo.units.m
diff_Co = 2.64e-06  # m2/h   https://www.aqion.de/site/diffusion-coefficients
diafiltrate_inlet_pressure = 101325 * pyo.units.Pa
diafiltrate_outlet_pressure = 145 * pyo.units.psi
Q_feed = 100 * pyo.units.m**3 / pyo.units.h
C_Li_feed = 1.7 * pyo.units.kg / pyo.units.m**3
C_Co_feed = 17 * pyo.units.kg / pyo.units.m**3
Q_diaf = 30 * pyo.units.m**3 / pyo.units.h
C_Li_diaf = 0.1 * pyo.units.kg / pyo.units.m**3
C_Co_diaf = 0.2 * pyo.units.kg / pyo.units.m**3


def main():
    m = build_and_init_model()

    (
        R_Li,
        R_Co,
        P_Li,
        P_Co,
        membrane_length,
        feed_flow_rate,
        diafiltrate_flow_rate,
        permeate_flow_rate,
    ) = get_values(m)

    add_costing(
        m,
        membrane_length,
        membrane_width=w,
        water_flux=Jw,
        feed_flow=feed_flow_rate,
        perm_flow=permeate_flow_rate,
        inlet_pressure=diafiltrate_inlet_pressure,
        outlet_pressure=diafiltrate_outlet_pressure,
        diaf_flow=diafiltrate_flow_rate,
    )

    # print(f"DOF = {degrees_of_freedom(m)}")
    # dt = DiagnosticsToolbox(m)
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
    # dt.report_structural_issues()
    # dt.report_numerical_issues()


def build_and_init_model():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = LiCoParameters()

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

    # Add mass transfer variables and constraints
    m.fs.solutes = pyo.Set(initialize=["Li", "Co"])

    m.fs.sieving_coefficient = pyo.Var(
        m.fs.solutes,
        units=pyo.units.dimensionless,
    )
    m.fs.sieving_coefficient["Li"].fix(1.3)
    m.fs.sieving_coefficient["Co"].fix(0.5)

    m.fs.stage1.length = pyo.Var(units=pyo.units.m)

    # Start by initializing with a short length
    # Too long and we deplete solvent due to lack of recycles
    m.fs.stage1.length.fix(10)

    solvent_rule, solute_rule, stage3_solute_rule = create_rules(m)

    m.fs.stage1.solvent_flux = pyo.Constraint(
        m.fs.stage1.elements,
        rule=solvent_rule,
    )

    m.fs.stage1.solute_sieving = pyo.Constraint(
        m.fs.stage1.elements,
        m.fs.solutes,
        rule=solute_rule,
    )

    m.fs.stage1.retentate_inlet.flow_vol[0].fix(Q_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].fix(C_Co_feed)

    assert degrees_of_freedom(m.fs.stage1) == 0

    initializer = MSContactorInitializer()
    initializer.initialize(m.fs.stage1)

    # Create an instance of the MSContactor model for the second stage here
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

    # Also add the necessary variables and constraints for the material flow terms
    # Do not forget to fix any design variables you might need to add
    m.fs.stage2.length = pyo.Var(units=pyo.units.m)
    m.fs.stage2.length.fix(10)

    m.fs.stage2.solvent_flux = pyo.Constraint(
        m.fs.stage2.elements,
        rule=solvent_rule,
    )

    m.fs.stage2.solute_sieving = pyo.Constraint(
        m.fs.stage2.elements,
        m.fs.solutes,
        rule=solute_rule,
    )

    m.fs.mix1 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.stream1 = Arc(
        source=m.fs.stage1.permeate_outlet,
        destination=m.fs.mix1.inlet_2,
    )
    m.fs.stream2 = Arc(
        source=m.fs.mix1.outlet,
        destination=m.fs.stage2.retentate_inlet,
    )

    # Add an Arc connecting the retentate of the 2nd stage to the retentate feed of the 1st stage
    m.fs.stream3 = Arc(
        source=m.fs.stage2.retentate_outlet,
        destination=m.fs.stage1.retentate_inlet,
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.stage1.retentate_inlet.flow_vol[0].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
    m.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

    m.fs.mix1.inlet_1.flow_vol[0].fix(Q_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Li"].fix(C_Li_feed)
    m.fs.mix1.inlet_1.conc_mass_solute[0, "Co"].fix(C_Co_feed)

    assert degrees_of_freedom(m) == 0

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

    # Add third stage, including side stream feed at element 10
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

    m.fs.stage3.length = pyo.Var(units=pyo.units.m)
    m.fs.stage3.length.fix(10)

    m.fs.stage3.solvent_flux = pyo.Constraint(
        m.fs.stage3.elements,
        rule=solvent_rule,
    )
    m.fs.stage3.solute_sieving = pyo.Constraint(
        m.fs.stage3.elements,
        m.fs.solutes,
        rule=stage3_solute_rule,
    )

    # Add second feed mixer
    m.fs.mix2 = Mixer(
        num_inlets=2,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    # Add 3 Arcs to connect third stage to remainder of flowsheet
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

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

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

    assert degrees_of_freedom(m) == 0

    # Initialize Mixer 2
    propagate_state(
        destination=m.fs.mix2.inlet_2,
        source=m.fs.stage2.permeate_outlet,
    )

    mix_initializer.initialize(m.fs.mix2)

    # Try to initialize stage 3
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

    add_expressions(m)

    return m


def create_rules(m):
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

        return pyo.log(b.retentate[0, s].conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * pyo.log(in_state.flow_vol) == pyo.log(in_state.conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * pyo.log(
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

        return pyo.log(b.retentate[0, s].conc_mass_solute[j]) + (
            m.fs.sieving_coefficient[j] - 1
        ) * pyo.log(q_in) == pyo.log(c_in) + (
            m.fs.sieving_coefficient[j] - 1
        ) * pyo.log(
            b.retentate[0, s].flow_vol
        )

    return (solvent_rule, solute_rule, stage3_solute_rule)


def solve_model(m):
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=True)


def add_expressions(m):
    m.Li_recovery = pyo.Expression(
        expr=m.fs.stage3.permeate_outlet.flow_vol[0]
        * m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
        / (
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Li"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"]
        )
    )

    m.Li_purity = pyo.Expression(
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

    m.Co_recovery = pyo.Expression(
        expr=m.fs.stage1.retentate_outlet.flow_vol[0]
        * m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
        / (
            m.fs.mix2.inlet_1.flow_vol[0] * m.fs.mix2.inlet_1.conc_mass_solute[0, "Co"]
            + m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * m.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"]
        )
    )

    m.Co_purity = pyo.Expression(
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


def get_values(m):
    R_Li = round(pyo.value(m.Li_recovery) * 100, 2)
    R_Co = round(pyo.value(m.Co_recovery) * 100, 2)

    P_Li = round(pyo.value(m.Li_purity) * 100, 2)
    P_Co = round(pyo.value(m.Co_purity) * 100, 2)

    membrane_length = (
        pyo.value(m.fs.stage1.length)
        + pyo.value(m.fs.stage2.length)
        + pyo.value(m.fs.stage3.length)
    ) * pyo.units.m

    # the feed is the stream entering stage 3 in the last element
    feed_flow_rate = (
        pyo.value(m.fs.stage3.retentate_side_stream_state[0, 10].flow_vol)
        * pyo.units.m**3
        / pyo.units.hr
    )

    # the diafiltrate is the stream entering stage 3 in the first element
    diafiltrate_flow_rate = (
        pyo.value(m.fs.stage3.retentate_inlet.flow_vol[0])
        * pyo.units.m**3
        / pyo.units.hr
    )

    # the permeate is the stream exiting stage 3 in the last element
    permeate_flow_rate = (
        pyo.value(m.fs.stage3.permeate_outlet.flow_vol[0])
        * pyo.units.m**3
        / pyo.units.hr
    )

    return (
        R_Li,
        R_Co,
        P_Li,
        P_Co,
        membrane_length,
        feed_flow_rate,
        diafiltrate_flow_rate,
        permeate_flow_rate,
    )


def add_costing(
    m,
    membrame_length,
    membrane_width,
    water_flux,
    feed_flow,
    perm_flow,
    inlet_pressure,
    outlet_pressure,
    diaf_flow,
):
    # creating dummy variables to store the UnitModelCostingBlocks
    # TODO: assign costing blocks to each stage
    m.fs.membrane = UnitModelBlock()
    m.fs.pump = UnitModelBlock()

    m.fs.costing = DiafiltrationCosting()
    m.fs.membrane.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": membrame_length,
            "membrane_width": membrane_width,
            "water_flux": water_flux,
            "vol_flow_feed": feed_flow,
            "vol_flow_perm": perm_flow,
        },
    )

    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": inlet_pressure,
            "outlet_pressure": outlet_pressure,
            "inlet_vol_flow": diaf_flow,
        },
    )
    m.fs.costing.cost_process()


def unfix_variables(m):
    m.fs.stage1.length.unfix()
    m.fs.stage2.length.unfix()
    m.fs.stage3.length.unfix()


def add_constraints(m):
    m.Co_recovery_constraint = pyo.Constraint(expr=m.Co_recovery >= 0.4)
    m.Li_recovery_constraint = pyo.Constraint(expr=m.Li_recovery >= 0.95)

    # m.Li_purity_constraint = pyo.Constraint(expr=m.Li_purity >= 0.8)
    # m.Co_purity_constraint = pyo.Constraint(expr=m.Co_purity >= 0.5)


def add_objective(m):
    def cost_obj(m):
        return m.fs.costing.total_annualized_cost

    m.cost_objecticve = pyo.Objective(rule=cost_obj)


def print_information(m):
    (
        R_Li,
        R_Co,
        P_Li,
        P_Co,
        membrane_length,
        feed_flow_rate,
        diafiltrate_flow_rate,
        permeate_flow_rate,
    ) = get_values(m)

    print(f"The lithium recovery is {R_Li}% at purity {P_Li}")
    print(f"The cobalt recovery is {R_Co}% at purity {P_Co}")

    print("\nmembrane area")
    print(f"{pyo.value(m.fs.membrane.costing.membrane_area)} m2")

    print("\nmembrane capital cost")
    print(f"${pyo.value(m.fs.membrane.costing.capital_cost)}")

    print("\npump capital cost")
    print(f"${pyo.value(m.fs.pump.costing.capital_cost)}")

    print("\naggregate capital cost")
    print(f"${pyo.value(m.fs.costing.aggregate_capital_cost)}")

    print("\ntotal capital cost")
    print(f"${pyo.value(m.fs.costing.total_capital_cost)}")

    print("\ntotal operating cost")
    print(f"${pyo.value(m.fs.costing.total_operating_cost)}")

    print("\ntotal annualized cost")
    print(f"${pyo.value(m.fs.costing.total_annualized_cost)}")


if __name__ == "__main__":
    main()
