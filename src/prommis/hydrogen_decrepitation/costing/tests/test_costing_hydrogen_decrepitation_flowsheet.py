import pytest
from pyomo.environ import(
    ConcreteModel,
    Var,
    Param,
    Constraint,
    Expression,
    )
from pyomo.environ import units as pyunits, value, assert_optimal_termination
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
import idaes.logger as idaeslog
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

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)
from prommis.hydrogen_decrepitation.hydrogen_decrepitation_flowsheet import (
    initialize_and_solve,
    main,
)
from prommis.hydrogen_decrepitation.repm_solids_properties import REPMParameters

from prommis.hydrogen_decrepitation.costing.cost_hydrogen_decrepitation_furnace import (
    HydrogenDecrepitationCostingData,
    HydrogenDecrepitationCosting,
)
from prommis.uky.costing.ree_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
    custom_REE_plant_currency_units,
)

_log = idaeslog.getLogger(__name__)

solver = get_solver("ipopt")


class TestHydrogenDecrepitationQGESS:
    @pytest.fixture(scope="class")
    def model(self):
        model = main()
        return model

    @pytest.mark.component
    def test_flowsheet_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_build_unit_model_costing(self, model):
        model.fs.costing = QGESSCosting()

        CE_index_year = "Jan_2024"

        # update the base currency for the costing block
        model.fs.costing.base_currency = pyunits.USD_Jan_2024

        # Source 2, 2.1 is HDD shredder
        # this is a constant-cost unit, where n_equip is the scaling parameter
        HDD_Recycling_shredder_accounts = ["2.1"]
        model.fs.shredder.n_equip = Var(
            initialize=1, units=pyunits.dimensionless
        )  # 1 shredder for 2700 HDDs per hour, 600 2.5 inch and 2100 3.5 inch
        model.fs.shredder.n_equip.fix()

        model.fs.shredder.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": HDD_Recycling_shredder_accounts,
                "scaled_param": model.fs.shredder.n_equip,  # 1 shredder
                "source": 2,
                # no. units is the scaling parameter for constant-cost units,
                # so use n_equip below to specify the number of shredders
                "n_equip": 1,
                "scale_down_parallel_equip": False,
                "CE_index_year": CE_index_year,
            },
        )

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

        assert isinstance(model.fs.shredder.costing.bare_erected_cost, Var)

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
    def test_build_process_costing(self, model):
        CE_index_year = "Jan_2024"

        model.fs.costing.build_process_costs(
            labor_types=["unskilled",],
            labor_rate=[38.20,],
            operators_per_shift=[2,],
            hours_per_shift=8,
            shifts_per_day=3,
            operating_days_per_year=336,
            # powder product is a mixed basket that would require additional processing to separate
            pure_product_output_rates={
                "Nd2Fe14B": 0.0000 * pyunits.kg / pyunits.hr,
                "Nd": 0.0000 * pyunits.kg / pyunits.hr,
                },
            mixed_product_output_rates={
                "Nd2Fe14B":
                    model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass *
                    model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd2Fe14B"],
                "Nd":
                    model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass *
                    model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd"],
                },
                # https://www.msesupplies.com/products/mse-pro-99-9-neodymium-iron-boron-magnetic-powder-ndfeb-100g?variant=40921835307066
            sale_prices={"Nd2Fe14B": 3899.90 * pyunits.USD_2023 / pyunits.kg},
            CE_index_year=CE_index_year,
            )

    @pytest.mark.component
    def test_costing_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

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

        assert pyunits.get_units(model.fs.shredder.costing.bare_erected_cost["2.1"]) == pyunits.MUSD_Jan_2024

        assert value(
            model.fs.shredder.costing.bare_erected_cost[
                "2.1"
            ]
        ) == pytest.approx(0.065465, rel=1e-4)


        assert pyunits.get_units(model.fs.hydrogen_decrepitation_furnace.costing.capital_cost) == pyunits.USD_Jan_2024

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


        assert pyunits.get_units(model.fs.costing.total_plant_cost) == pyunits.MUSD_Jan_2024
        model.fs.costing.report()

        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            0.21398, rel=1e-4
        )
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.072047, rel=1e-4)
        assert value(model.fs.costing.total_installation_cost) == pytest.approx(
            0.14193, rel=1e-4
        )
        assert value(model.fs.costing.other_plant_costs) == pytest.approx(
            1.0000e-12, rel=1e-4
        )
        assert value(model.fs.costing.annual_operating_labor_cost) == pytest.approx(
            0.77011, rel=1e-4
        )
        assert value(model.fs.costing.annual_technical_labor_cost) == pytest.approx(
            1.0000e-12, rel=1e-4
        )
        assert value(model.fs.costing.annual_labor_cost) == pytest.approx(0.77011, rel=1e-4)
        assert value(model.fs.costing.maintenance_and_material_cost) == pytest.approx(
            0.0042796, rel=1e-4
        )
        assert value(model.fs.costing.quality_assurance_and_control_cost) == pytest.approx(
            0.077011, rel=1e-4
        )
        assert value(model.fs.costing.sales_patenting_and_research_cost) == pytest.approx(
            2.0855, rel=1e-4
        )
        assert value(model.fs.costing.admin_and_support_labor_cost) == pytest.approx(
            0.15402, rel=1e-4
        )
        assert value(model.fs.costing.property_taxes_and_insurance_cost) == pytest.approx(
            0.0021398, rel=1e-4
        )
        assert value(model.fs.costing.other_fixed_costs) == pytest.approx(1.0000e-12, rel=1e-4)
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(3.0930, rel=1e-4)
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(417.09, rel=1e-4)
