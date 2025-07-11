from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    Var,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

import pytest

from prommis.hydrogen_decrepitation.costing.cost_hydrogen_decrepitation_furnace import (
    HydrogenDecrepitationCosting,
    HydrogenDecrepitationCostingData,
)
from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)
from prommis.hydrogen_decrepitation.repm_solids_properties import REPMParameters

_log = idaeslog.getLogger(__name__)

solver = get_solver("ipopt")


def base_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {
        "H2",
    }
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_solid = REPMParameters(
        doc="solid property",
    )

    m.fs.prop_gas.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.prop_gas.set_default_scaling("pressure", 1e-5)
    m.fs.prop_gas.set_default_scaling("temperature", 1e-2)
    m.fs.prop_gas.set_default_scaling("flow_mol", 1e1)
    m.fs.prop_gas.set_default_scaling("flow_mol_phase", 1e1)
    m.fs.prop_gas.set_default_scaling("_energy_density_term", 1e-4)
    m.fs.prop_gas.set_default_scaling("phase_frac", 1)

    _mf_scale = {
        "H2": 1,
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )

    m.fs.hydrogen_decrepitation_furnace = REPMHydrogenDecrepitationFurnace(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Nd",
        ],
    )

    m.fs.hydrogen_decrepitation_furnace.solid_in[0].flow_mass.fix(
        0.0057367 * pyunits.kg / pyunits.s
    )
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd2Fe14B"].fix(0.99)
    m.fs.hydrogen_decrepitation_furnace.solid_in[0].mass_frac_comp["Nd"] = 0.01

    # don't fix, already have mole frac balance so just need initial value
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.mole_frac_comp[0, "H2"].fix(1)

    # inlet flue gas mole flow rate, stoichiometric on molar basis with REPM
    @m.fs.hydrogen_decrepitation_furnace.Constraint(m.fs.time)
    def flow_mol_gas_constraint(b, t):
        return b.gas_inlet.flow_mol[t] == (
            sum(
                b.flow_mol_comp_impurity_feed[t, c]
                for c in m.fs.prop_solid.component_list
            )
        )

    # operating parameters
    m.fs.hydrogen_decrepitation_furnace.deltaP.fix(0)
    m.fs.hydrogen_decrepitation_furnace.operating_temperature.fix(443.15)
    m.fs.hydrogen_decrepitation_furnace.gas_inlet.pressure.fix(101325)

    m.fs.hydrogen_decrepitation_furnace.decrepitation_duration.set_value(
        10800 * pyunits.s
    )
    m.fs.hydrogen_decrepitation_furnace.sample_density.set_value(
        m.fs.prop_solid.dens_mass
    )  # 7500 kg/m3
    m.fs.hydrogen_decrepitation_furnace.chamber_to_sample_ratio.set_value(2)

    # solid temperature, cools back to inlet temperature during shutdown
    m.fs.hydrogen_decrepitation_furnace.temp_feed.fix(298.15)
    m.fs.hydrogen_decrepitation_furnace.temp_prod.fix(298.15)

    # gas temperature, assume comes in at operating temperature
    m.fs.hydrogen_decrepitation_furnace.gas_in[0].temperature.fix(443.15)

    # no additional heat is supplied other than what's required for decrepitation
    m.fs.hydrogen_decrepitation_furnace.supplied_heat_duty.fix(0)

    return m


class TestHydrogenDecrepitationCostingGasFired:
    @pytest.fixture(scope="class")
    def model(self):
        model = base_model()
        return model

    @pytest.mark.component
    def test_base_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_base_model_attributes(self, model):
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.furnace_chamber_volume,
            Var,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.heat_loss, Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_metal_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_metal_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_insulation_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_metal_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_metal_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_volume,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_mass,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.max_temperature, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_insulation_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_metal_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_metal_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_insulation_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_insulation_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_metal_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_metal_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_furnace_ext_surface,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.air_heat_transfer_coeff,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.ref_temp, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.length_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.width_insulation1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.min_quantity_insulation, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.weight_insulation1,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.density_metal1, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.length_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.width_insulation2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.weight_insulation2,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.density_metal2, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_metal1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_metal2,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.ramp_up_time, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.decrepitation_duration,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.operating_temperature,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_heat_capacity,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.preparation_time, Param
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.cool_down_time, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_radiation,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_insulation_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_insulation_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_metal_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_metal_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_insulation_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_insulation_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_metal_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_metal_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_external_air_convection,
            Constraint,
        )

        assert model.fs.hydrogen_decrepitation_furnace.config.number_of_units == 1

    @pytest.mark.component
    def test_hd_costing(self, model):
        model.fs.costing = HydrogenDecrepitationCosting()

        model.fs.hydrogen_decrepitation_furnace.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=HydrogenDecrepitationCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "price_insulation1": 183.81,  # in USD
                "price_metal1": 3.14,  # USD/kg
                "price_insulation2": 47.00,  # in USD
                "price_metal2": 1.50,  # USD/kg
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "utility_rate": 0.081,  # USD/kWhr
                "heating_mode": 0,
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "engineering_and_drafting": 1000,  # USD
            },
        )

        assert isinstance(model.fs.hydrogen_decrepitation_furnace.costing.OPEX, Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit, Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost, Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.labor_rate, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.efficiency, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.hours_per_shift, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.shifts_per_day, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_days_per_year,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.utility_rate, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_controller_price,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.engineering_and_drafting,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_cost_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_eq,
            Constraint,
        )

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solved_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_results(self, model):

        assert (
            pyunits.get_units(
                model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
            )
            == pyunits.USD_Jan_2024
        )

        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
        ) == pytest.approx(6581.66, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit
        ) == pytest.approx(6581.66, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit
        ) == pytest.approx(2592.89, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost
        ) == pytest.approx(2592.89, rel=1e-4)


class TestHydrogenDecrepitationCostingElectric:
    @pytest.fixture(scope="class")
    def model(self):
        model = base_model()
        return model

    @pytest.mark.component
    def test_base_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_base_model_attributes(self, model):
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.furnace_chamber_volume,
            Var,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.heat_loss, Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_insulation_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_metal_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_metal_material1,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_insulation_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_metal_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_ratio_metal_material2,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_volume,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_mass,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.max_temperature, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_insulation_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_metal_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_metal_material1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_insulation_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_insulation_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_metal_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thermal_cond_metal_material2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.temperature_furnace_ext_surface,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.air_heat_transfer_coeff,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.ref_temp, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.length_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.width_insulation1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.min_quantity_insulation, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.weight_insulation1,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.density_metal1, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.length_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.width_insulation2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.thickness_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.weight_insulation2,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.density_metal2, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_insulation1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_metal1,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_insulation2,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.specific_heat_capacity_metal2,
            Param,
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.ramp_up_time, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.decrepitation_duration,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.operating_temperature,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.sample_heat_capacity,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.preparation_time, Param
        )
        assert isinstance(model.fs.hydrogen_decrepitation_furnace.cool_down_time, Param)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_radiation,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_insulation_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_insulation_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_metal_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_metal_material1,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_insulation_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_insulation_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.relative_thickness_metal_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_metal_material2,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.heat_loss_external_air_convection,
            Constraint,
        )

        assert model.fs.hydrogen_decrepitation_furnace.config.number_of_units == 1

    @pytest.mark.component
    def test_hd_costing(self, model):
        model.fs.costing = HydrogenDecrepitationCosting()

        model.fs.hydrogen_decrepitation_furnace.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=HydrogenDecrepitationCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "price_insulation1": 183.81,  # in USD
                "price_metal1": 3.14,  # USD/kg
                "price_insulation2": 47.00,  # in USD
                "price_metal2": 1.50,  # USD/kg
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "utility_rate": 0.081,  # USD/kWhr
                "heating_mode": 1,
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "engineering_and_drafting": 1000,  # USD
            },
        )

        assert isinstance(model.fs.hydrogen_decrepitation_furnace.costing.OPEX, Var)
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit, Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost, Var
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost,
            Var,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.labor_rate, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal1, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_insulation2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.price_metal2, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.efficiency, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.hours_per_shift, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.shifts_per_day, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_days_per_year,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.utility_rate, Param
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.temperature_controller_price,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.engineering_and_drafting,
            Param,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.operating_cost_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit_eq,
            Constraint,
        )
        assert isinstance(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_eq,
            Constraint,
        )

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_solved_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_results(self, model):

        assert (
            pyunits.get_units(
                model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
            )
            == pyunits.USD_Jan_2024
        )

        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
        ) == pytest.approx(7014.91, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit
        ) == pytest.approx(7014.91, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit
        ) == pytest.approx(2592.89, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost
        ) == pytest.approx(2592.89, rel=1e-4)


@pytest.mark.unit
def test_invalid_parent_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = HydrogenDecrepitationCosting()
    m.fs.unit = UnitModelBlock()

    with pytest.raises(TypeError) as e:
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=HydrogenDecrepitationCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "price_insulation1": 183.81,  # in USD
                "price_metal1": 3.14,  # USD/kg
                "price_insulation2": 47.00,  # in USD
                "price_metal2": 1.50,  # USD/kg
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "utility_rate": 0.081,  # USD/kWhr
                "heating_mode": 0,
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "engineering_and_drafting": 1000,  # USD
            },
        )

    # split up string match checks since block identifiers will change on each run
    assert "Parent block is of type " in str(e.value)
    assert "ScalarUnitModelBlock" in str(e.value)
    assert " and should be of type " in str(e.value)
    assert "REPMHydrogenDecrepitationFurnace" in str(e.value)
    assert " to use costing model." in str(e.value)


@pytest.mark.unit
def test_invalid_heating_mode():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = HydrogenDecrepitationCosting()

    gas_species = {
        "H2",
    }
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_solid = REPMParameters(
        doc="solid property",
    )

    m.fs.unit = REPMHydrogenDecrepitationFurnace(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Nd",
        ],
    )

    with pytest.raises(
        TypeError,
        match="Valid heating modes are either 0: electric-fired or 1: gas-fired.",
    ):
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=HydrogenDecrepitationCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "price_insulation1": 183.81,  # in USD
                "price_metal1": 3.14,  # USD/kg
                "price_insulation2": 47.00,  # in USD
                "price_metal2": 1.50,  # USD/kg
                "hours_per_shift": 8,  # hr
                "shifts_per_day": 3,
                "operating_days_per_year": 336,  # days
                "efficiency": 0.95,
                "utility_rate": 0.081,  # USD/kWhr
                "heating_mode": 2,
                "labor_rate": 75,  # USD/hr
                "temperature_controller_price": 129.00,  # USD
                "engineering_and_drafting": 1000,  # USD
            },
        )
