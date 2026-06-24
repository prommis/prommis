from pyomo.environ import Constraint, Param, Var, assert_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import pytest

from prommis.hydrogen_decrepitation.costing.cost_hydrogen_decrepitation_furnace import (
    HydrogenDecrepitationCostingData,
)
from prommis.hydrogen_decrepitation.hydrogen_decrepitation_flowsheet import main
from prommis.uky.costing.ree_plant_capcost import QGESSCosting, QGESSCostingData

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

        # HD-specific params that QGESSCosting doesn't define
        model.fs.costing.price_insulation1 = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024,
        )
        model.fs.costing.price_insulation1.set_value(183.81 * pyunits.USD_Jan_2024)

        model.fs.costing.labor_rate = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024 / pyunits.hr,
        )
        model.fs.costing.labor_rate.set_value(75 * pyunits.USD_Jan_2024 / pyunits.hr)

        model.fs.costing.price_metal1 = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024 / pyunits.kg,
        )
        model.fs.costing.price_metal1.set_value(
            3.14 * pyunits.USD_Jan_2024 / pyunits.kg
        )

        model.fs.costing.price_insulation2 = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024,
        )
        model.fs.costing.price_insulation2.set_value(47.00 * pyunits.USD_Jan_2024)

        model.fs.costing.price_metal2 = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024 / pyunits.kg,
        )
        model.fs.costing.price_metal2.set_value(
            1.50 * pyunits.USD_Jan_2024 / pyunits.kg
        )

        model.fs.costing.efficiency = Param(
            mutable=True,
            units=pyunits.dimensionless,
        )
        model.fs.costing.efficiency.set_value(0.95 * pyunits.dimensionless)

        model.fs.costing.capacity_factor = Param(
            mutable=True,
            units=pyunits.dimensionless,
        )
        model.fs.costing.capacity_factor.set_value(0.92)

        model.fs.costing.utility_rate = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024 / (pyunits.kW * pyunits.hr),
        )
        model.fs.costing.utility_rate.set_value(
            0.081 * pyunits.USD_Jan_2024 / (pyunits.kW * pyunits.hr)
        )

        model.fs.costing.temperature_controller_price = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024,
        )
        model.fs.costing.temperature_controller_price.set_value(
            129.00 * pyunits.USD_Jan_2024
        )

        model.fs.costing.engineering_and_drafting = Param(
            mutable=True,
            units=pyunits.USD_Jan_2024,
        )
        model.fs.costing.engineering_and_drafting.set_value(1000 * pyunits.USD_Jan_2024)

        model.fs.costing.CE_index_year = Param(
            initialize="Jan_2024",
            mutable=True,
        )

        model.fs.hydrogen_decrepitation_furnace.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=HydrogenDecrepitationCostingData.cost_hydrogen_decrepitation_furnace,
            costing_method_arguments={
                "heating_mode": 0,
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
        assert isinstance(model.fs.costing.price_insulation1, Param)
        assert isinstance(model.fs.costing.labor_rate, Param)
        assert isinstance(model.fs.costing.price_metal1, Param)
        assert isinstance(model.fs.costing.price_insulation2, Param)
        assert isinstance(model.fs.costing.price_metal2, Param)
        assert isinstance(model.fs.costing.efficiency, Param)
        assert isinstance(model.fs.costing.capacity_factor, Param)
        assert isinstance(model.fs.costing.utility_rate, Param)
        assert isinstance(
            model.fs.costing.temperature_controller_price,
            Param,
        )
        assert isinstance(
            model.fs.costing.engineering_and_drafting,
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
            labor_types=[
                "unskilled",
            ],
            labor_rate=[
                38.20,
            ],
            operators_per_shift=[
                2,
            ],
            capacity_factor=0.92,
            # powder product is a mixed basket that would require additional processing to separate
            pure_product_output_rates={
                "Nd2Fe14B": 0.0000 * pyunits.kg / pyunits.hr,
                "Nd": 0.0000 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Nd2Fe14B": model.fs.hydrogen_decrepitation_furnace.solid_out[
                    0
                ].flow_mass
                * model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp[
                    "Nd2Fe14B"
                ],
                "Nd": model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass
                * model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp[
                    "Nd"
                ],
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
        assert (
            pyunits.get_units(model.fs.shredder.costing.bare_erected_cost["2.1"])
            == pyunits.MUSD_Jan_2024
        )

        assert value(
            model.fs.shredder.costing.bare_erected_cost["2.1"]
        ) == pytest.approx(0.065465, rel=1e-4)

        assert (
            pyunits.get_units(
                model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
            )
            == pyunits.USD_Jan_2024
        )

        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.capital_cost
        ) == pytest.approx(6567.89, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.base_cost_per_unit
        ) == pytest.approx(6567.89, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost_per_unit
        ) == pytest.approx(2586.43, rel=1e-4)
        assert value(
            model.fs.hydrogen_decrepitation_furnace.costing.variable_operating_cost
        ) == pytest.approx(2586.43, rel=1e-4)

        assert (
            pyunits.get_units(model.fs.costing.total_plant_cost)
            == pyunits.MUSD_Jan_2024
        )

        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            0.21394, rel=1e-4
        )
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.072033, rel=1e-4)
        assert value(model.fs.costing.total_installation_cost) == pytest.approx(
            0.14190, rel=1e-4
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
        assert value(model.fs.costing.annual_labor_cost) == pytest.approx(
            0.77011, rel=1e-4
        )
        assert value(model.fs.costing.maintenance_and_material_cost) == pytest.approx(
            0.0042788, rel=1e-4
        )
        assert value(
            model.fs.costing.quality_assurance_and_control_cost
        ) == pytest.approx(0.077011, rel=1e-4)
        assert value(
            model.fs.costing.sales_patenting_and_research_cost
        ) == pytest.approx(2.0855, rel=1e-4)
        assert value(model.fs.costing.admin_and_support_labor_cost) == pytest.approx(
            0.15404, rel=1e-4
        )
        assert value(
            model.fs.costing.property_taxes_and_insurance_cost
        ) == pytest.approx(0.0021394, rel=1e-4)
        assert value(model.fs.costing.other_fixed_costs) == pytest.approx(
            1.0000e-12, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            3.0930, rel=1e-4
        )
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(
            417.09, rel=1e-4
        )
