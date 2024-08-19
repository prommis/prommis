import pytest
import pyomo.environ as pyo
from pyomo.environ import units as pyunits, value, assert_optimal_termination
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from prommis.uky.costing.models.hydrogen_decrepitation import (
    REEEquipmentCostingData,
    REEEquipmentCosting,
)

_log = idaeslog.getLogger(__name__)

solver = get_solver("ipopt")


def base_model():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_units=pyunits.s)
    m.fs.costing = REEEquipmentCosting()
    m.fs.hydrogen_decrepitation_furnace = UnitModelBlock()

    return m


class TestREEEquipmentCosting:
    @pytest.fixture(scope="class")
    def model(self):
        model = base_model()
        return model

    @pytest.mark.component
    def test_base_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    # costing test
    @pytest.mark.component
    def test_hd_costing(self, model):
        model.fs.hydrogen_decrepitation_furnace.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=REEEquipmentCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "ramp_up_time": 300 * pyunits.s,
                "operating_temperature": 443.15 * pyunits.K,
                "decrepitation_duration": 10800 * pyunits.s,
                "preparation_time": 3600 * pyunits.s,
                "cool_down_time": 3600 * pyunits.s,
                "sample_heat_capacity": 0.44 * pyunits.kJ / (pyunits.kg * pyunits.K),
                "sample_mass": None,
                "sample_volume": 0.002 * ((pyunits.m) ** 3),
                "sample_density": 7500 * pyunits.kg / (pyunits.m**3),
                "chamber_to_sample_ratio": 2 * pyunits.dimensionless,
                "length_insulation1": 7.62 * pyunits.m,
                "width_insulation1": 0.6096 * pyunits.m,
                "thickness_insulation1": 0.0254 * pyunits.m,
                "price_insulation1": 183.81 * pyunits.USD_Jan_2024,
                "weight_insulation1": 15.42 * pyunits.kg,
                "thermal_cond_insulation_material1": 0.33
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "price_metal1": 3.14 * pyunits.USD_Jan_2024 / pyunits.kg,
                "density_metal1": 7473.57 * pyunits.kg / (pyunits.m**3),
                "thermal_cond_metal_material1": 13.53
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "length_insulation2": 1.19 * pyunits.m,
                "width_insulation2": 0.381 * pyunits.m,
                "thickness_insulation2": 0.0889 * pyunits.m,
                "price_insulation2": 47.00 * pyunits.USD_Jan_2024,
                "weight_insulation2": 10.43 * pyunits.kg,
                "thermal_cond_insulation_material2": 0.069
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "price_metal2": 1.50 * pyunits.USD_Jan_2024 / pyunits.kg,
                "density_metal2": 7861.09 * pyunits.kg / (pyunits.m**3),
                "thermal_cond_metal_material2": 45 * pyunits.W / pyunits.m / pyunits.K,
                "hours_per_shift": 8 * pyunits.hr,
                "shifts_per_day": 3 * (pyunits.day) ** (-1),
                "specific_heat_capacity_insulation1": 1.08
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_metal1": 0.468
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_insulation2": 0.9
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_metal2": 0.502416
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "operating_days_per_year": 336 * pyunits.day / pyunits.year,
                "efficiency": 0.95 * pyunits.dimensionless,
                "electricity_rate": 0.081
                * pyunits.USD_Jan_2024
                / (pyunits.kW * pyunits.hr),
                "labor_rate": 75 * pyunits.USD_Jan_2024 / pyunits.hr,
                "temperature_controller_price": 129.00 * pyunits.USD_Jan_2024,
                "number_of_units": 1 * pyunits.dimensionless,
                "engineering_and_drafting": 1000 * pyunits.USD_Jan_2024,
            },
        )

        assert (
            model.fs.hydrogen_decrepitation_furnace.costing.costing_package.base_currency
            == pyunits.USD_Jan_2024
        )
        assert (
            model.fs.hydrogen_decrepitation_furnace.costing.costing_package.base_period
            == pyunits.year
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.furnace_chamber_volume,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio3,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio4,
            pyo.Var,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.costing.OPEX, pyo.Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_volume, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_mass, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.stef_bolt_constant,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.max_temperature, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_insulation_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_metal_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_metal_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_insulation_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_metal_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_metal_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_furnace_ext_surface,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.air_heat_transfer_coeff,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.ref_temp, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.length_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.width_insulation1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.min_quantity, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.weight_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.labor_rate, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.density_metal1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.length_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.width_insulation2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.weight_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.density_metal2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_metal1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_metal2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.ramp_up_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.efficiency, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.decrepitation_duration,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_temperature,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_heat_capacity,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.preparation_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.cool_down_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.hours_per_shift, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.shifts_per_day, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_days_per_year,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.electricity_rate, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_controller_price,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.engineering_and_drafting,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.number_of_units, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint1,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint2,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint3,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint4,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint5,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint6,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint7,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint8,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint9,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint10,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_cost_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_eq,
            pyo.Constraint,
        )

    @pytest.mark.component
    def test_solve1(self, model):
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solved_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_result1(self, model):
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.furnace_chamber_volume
        ) == pytest.approx(0.004, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss
        ) == pytest.approx(145.35, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material1
        ) == pytest.approx(0.12258, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material1
        ) == pytest.approx(1168.89, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material1
        ) == pytest.approx(0.0050267, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material2
        ) == pytest.approx(0.17447, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material2
        ) == pytest.approx(0.0076342, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
        ) == pytest.approx(5353.04, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit
        ) == pytest.approx(5353.03, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit
        ) == pytest.approx(1397.77, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost
        ) == pytest.approx(1397.77, rel=1e-4)

    # Test with sample_mass known and sample volume is None
    @pytest.mark.component
    def test_sample_mass_known(self, model):
        model.fs.hydrogen_decrepitation_furnace.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=REEEquipmentCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "ramp_up_time": 300 * pyunits.s,
                "operating_temperature": 443.15 * pyunits.K,
                "decrepitation_duration": 10800 * pyunits.s,
                "preparation_time": 3600 * pyunits.s,
                "cool_down_time": 3600 * pyunits.s,
                "sample_heat_capacity": 0.44 * pyunits.kJ / (pyunits.kg * pyunits.K),
                "sample_mass": 62.025 * pyunits.kg,  # kg
                "sample_volume": None,  # m**3
                "sample_density": 7500 * pyunits.kg / (pyunits.m**3),
                "chamber_to_sample_ratio": 2 * pyunits.dimensionless,
                "length_insulation1": 7.62 * pyunits.m,
                "width_insulation1": 0.6096 * pyunits.m,
                "thickness_insulation1": 0.0254 * pyunits.m,
                "price_insulation1": 183.81 * pyunits.USD_Jan_2024,
                "weight_insulation1": 15.42 * pyunits.kg,
                "thermal_cond_insulation_material1": 0.33
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "price_metal1": 3.14 * pyunits.USD_Jan_2024 / pyunits.kg,
                "density_metal1": 7473.57 * pyunits.kg / (pyunits.m**3),
                "thermal_cond_metal_material1": 13.53
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "length_insulation2": 1.19 * pyunits.m,
                "width_insulation2": 0.381 * pyunits.m,
                "thickness_insulation2": 0.0889 * pyunits.m,
                "price_insulation2": 47.00 * pyunits.USD_Jan_2024,
                "weight_insulation2": 10.43 * pyunits.kg,
                "thermal_cond_insulation_material2": 0.069
                * pyunits.W
                / pyunits.m
                / pyunits.K,
                "price_metal2": 1.50 * pyunits.USD_Jan_2024 / pyunits.kg,
                "density_metal2": 7861.09 * pyunits.kg / (pyunits.m**3),
                "thermal_cond_metal_material2": 45 * pyunits.W / pyunits.m / pyunits.K,
                "hours_per_shift": 8 * pyunits.hr,
                "shifts_per_day": 3 * (pyunits.day) ** (-1),
                "specific_heat_capacity_insulation1": 1.08
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_metal1": 0.468
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_insulation2": 0.9
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "specific_heat_capacity_metal2": 0.502416
                * pyunits.kJ
                / (pyunits.kg * pyunits.K),
                "operating_days_per_year": 336 * pyunits.day / pyunits.year,
                "efficiency": 0.95 * pyunits.dimensionless,
                "electricity_rate": 0.081
                * pyunits.USD_Jan_2024
                / (pyunits.kW * pyunits.hr),
                "labor_rate": 75 * pyunits.USD_Jan_2024 / pyunits.hr,
                "temperature_controller_price": 129.00 * pyunits.USD_Jan_2024,
                "number_of_units": 1 * pyunits.dimensionless,
                "engineering_and_drafting": 1000 * pyunits.USD_Jan_2024,
            },
        )

        assert (
            model.fs.hydrogen_decrepitation_furnace.costing.costing_package.base_currency
            == pyunits.USD_Jan_2024
        )
        assert (
            model.fs.hydrogen_decrepitation_furnace.costing.costing_package.base_period
            == pyunits.year
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.furnace_chamber_volume,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material1,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio3,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material2,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.relative_thickness_ratio4,
            pyo.Var,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.costing.OPEX, pyo.Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost, pyo.Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost,
            pyo.Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_volume, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_mass, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.stef_bolt_constant,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.max_temperature, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_insulation_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_metal_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_metal_material1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_insulation_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_metal_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thermal_cond_metal_material2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_furnace_ext_surface,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.air_heat_transfer_coeff,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.ref_temp, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.length_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.width_insulation1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.min_quantity, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.weight_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.labor_rate, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.density_metal1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal1, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.length_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.width_insulation2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.weight_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.density_metal2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal2, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_insulation1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_metal1,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_insulation2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.specific_heat_capacity_metal2,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.ramp_up_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.efficiency, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.decrepitation_duration,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_temperature,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.sample_heat_capacity,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.preparation_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.cool_down_time, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.hours_per_shift, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.shifts_per_day, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_days_per_year,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.electricity_rate, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_controller_price,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.engineering_and_drafting,
            pyo.Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.number_of_units, pyo.Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint1,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint2,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint3,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint4,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint5,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint6,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint7,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint8,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint9,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss_constraint10,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_cost_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_eq,
            pyo.Constraint,
        )

    @pytest.mark.component
    def test_solve2(self, model):
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solved_model_diagnostics2(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_result2(self, model):
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.furnace_chamber_volume
        ) == pytest.approx(0.01654, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.heat_loss
        ) == pytest.approx(279.02, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material1
        ) == pytest.approx(0.14897, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_insulation_material1
        ) == pytest.approx(1169.98, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material1
        ) == pytest.approx(0.0056312, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_insulation_material2
        ) == pytest.approx(0.18378, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.thickness_metal_material2
        ) == pytest.approx(0.0076473, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
        ) == pytest.approx(8245.66, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit
        ) == pytest.approx(8245.66, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit
        ) == pytest.approx(2473.80, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost
        ) == pytest.approx(2473.80, rel=1e-4)
