import pytest
import pyomo.environ as pyo
from pyomo.environ import units as pyunits, value, assert_optimal_termination
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from prommis.uky.costing.models.HD_CE import (
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
                "ramp_up_time": 300,  # (in seconds)
                "operating_temperature": 443.15,  # in Kelvin
                "decrepitation_duration": 10800,  # (in seconds)
                "preparation_time": 3600,  # in seconds
                "cool_down_time": 3600,  # in seconds
                "sample_heat_capacity": 0.44,  # in kJ/(kg*K)
                "sample_mass": None,  # kg
                "sample_volume": 0.002,  # m**3
                "sample_density": 7500,  # in kg/m**3
                "chamber_to_sample_ratio": 2,
                "length_insulation1": 7.62,  # in m
                "width_insulation1": 0.6096,  # in m
                "thickness_insulation1": 0.0254,  # in m
                "price_insulation1": 183.81,  # in USD
                "weight_insulation1": 15.42,  # in kg
                "thermal_cond_insulation_material1": 0.33,  # W/(m*K)
                "price_metal1": 3.14,  # USD/kg
                "density_metal1": 7473.57,  # kg/m**3
                "thermal_cond_metal_material1": 13.53,  # W/(m*K)
                "length_insulation2": 1.19,  # in m
                "width_insulation2": 0.381,  # in m
                "thickness_insulation2": 0.0889,  # in m
                "price_insulation2": 47.00,  # in USD
                "weight_insulation2": 10.43,  # in kg
                "thermal_cond_insulation_material2": 0.069,  # W/(m*K)
                "price_metal2": 1.50,  # USD/kg
                "density_metal2": 7861.09,  # kg/m**3
                "thermal_cond_metal_material2": 45,  # W/(m*K)
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "specific_heat_capacity_insulation1": 1.08,  # KJ/(kg*K)
                "specific_heat_capacity_metal1": 0.468,  # KJ/(kg*K)
                "specific_heat_capacity_insulation2": 0.9,  # KJ/(kg*K)
                "specific_heat_capacity_metal2": 0.502416,  # KJ/(kg*K)
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "electricity_rate": 0.081,  # USD/kWhr
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "number_of_units": 1,
                "engineering_and_drafting": 1000,  # USD
            },
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
                "ramp_up_time": 300,  # (in seconds)
                "operating_temperature": 443.15,  # in Kelvin
                "decrepitation_duration": 10800,  # (in seconds)
                "preparation_time": 3600,  # in seconds
                "cool_down_time": 3600,  # in seconds
                "sample_heat_capacity": 0.44,  # in kJ/(kg*K)
                "sample_mass": 62.025,  # kg
                "sample_volume": None,  # m**3
                "sample_density": 7500,  # in kg/m**3
                "chamber_to_sample_ratio": 2,
                "length_insulation1": 7.62,  # in m
                "width_insulation1": 0.6096,  # in m
                "thickness_insulation1": 0.0254,  # in m
                "price_insulation1": 183.81,  # in USD
                "weight_insulation1": 15.42,  # in kg
                "thermal_cond_insulation_material1": 0.33,  # W/(m*K)
                "price_metal1": 3.14,  # USD/kg
                "density_metal1": 7473.57,  # kg/m**3
                "thermal_cond_metal_material1": 13.53,  # W/(m*K)
                "length_insulation2": 1.19,  # in m
                "width_insulation2": 0.381,  # in m
                "thickness_insulation2": 0.0889,  # in m
                "price_insulation2": 47.00,  # in USD
                "weight_insulation2": 10.43,  # in kg
                "thermal_cond_insulation_material2": 0.069,  # W/(m*K)
                "price_metal2": 1.50,  # USD/kg
                "density_metal2": 7861.09,  # kg/m**3
                "thermal_cond_metal_material2": 45,  # W/(m*K)
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "specific_heat_capacity_insulation1": 1.08,  # KJ/(kg*K)
                "specific_heat_capacity_metal1": 0.468,  # KJ/(kg*K)
                "specific_heat_capacity_insulation2": 0.9,  # KJ/(kg*K)
                "specific_heat_capacity_metal2": 0.502416,  # KJ/(kg*K)
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "electricity_rate": 0.081,  # USD/kWhr
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "number_of_units": 1,
                "engineering_and_drafting": 1000,  # USD
            },
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
