from pyomo.environ import (
    Param,
    Var,
    units as pyunits,
    PositiveReals,
)
from idaes.core import (
    declare_process_block_class,
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)
from idaes.core.util.math import smooth_max
import idaes.logger as idaeslog

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)

_log = idaeslog.getLogger(__name__)


def custom_REE_plant_currency_units():
    """
    Define conversion rates for US Dollars based on CE Index.
    """
    register_idaes_currency_units()
    if "USD_Jan_2024" in pyunits._pint_registry:  # pylint: disable=protected-access
        # Assume that custom REE plant units have already been registered
        # Log a message and end
        _log.debug(
            "Custom REE plant currency units (USD_Jan_2024) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                # from https://toweringskills.com/financial-analysis/cost-indices/ as of 01/02/2024
                "USD_Jan_2024 = 500/795.4 * USD_CE500",
            ]
        )


@declare_process_block_class("HydrogenDecrepitationCosting")
class HydrogenDecrepitationCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index
    custom_REE_plant_currency_units()

    def build_global_params(self):
        # Set the base year for all costs
        self.base_currency = pyunits.USD_Jan_2024
        # Set a base period for all operating costs
        self.base_period = pyunits.year

    def cost_hydrogen_decrepitation_furnace(
        blk,
        # ramp_up_time=300 * pyunits.s,
        # operating_temperature=443.15 * pyunits.K,
        # decrepitation_duration=10800 * pyunits.s,
        # preparation_time=3600 * pyunits.s,
        # cool_down_time=3600 * pyunits.s,
        # sample_heat_capacity=0.44 * pyunits.kJ / (pyunits.kg * pyunits.K),
        # sample_mass=None,
        # sample_volume=None,
        # sample_density=7500 * pyunits.kg / (pyunits.m**3),
        # chamber_to_sample_ratio=2 * pyunits.dimensionless,
        # length_insulation1=7.62 * pyunits.m,
        # width_insulation1=0.6096 * pyunits.m,
        # thickness_insulation1=0.0254 * pyunits.m,
        price_insulation1=183.81 * pyunits.USD_Jan_2024,
        # weight_insulation1=15.42 * pyunits.kg,
        # thermal_cond_insulation_material1=0.33 * pyunits.W / pyunits.m / pyunits.K,
        price_metal1=3.14 * pyunits.USD_Jan_2024 / pyunits.kg,
        # density_metal1=7473.57 * pyunits.kg / (pyunits.m**3),
        # thermal_cond_metal_material1=13.53 * pyunits.W / pyunits.m / pyunits.K,
        # length_insulation2=1.19 * pyunits.m,
        # width_insulation2=0.381 * pyunits.m,
        # thickness_insulation2=0.0889 * pyunits.m,
        price_insulation2=47.00 * pyunits.USD_Jan_2024,
        # weight_insulation2=10.43 * pyunits.kg,
        # thermal_cond_insulation_material2=0.069 * pyunits.W / pyunits.m / pyunits.K,
        price_metal2=1.50 * pyunits.USD_Jan_2024 / pyunits.kg,
        # density_metal2=7861.09 * pyunits.kg / (pyunits.m**3),
        # thermal_cond_metal_material2=45 * pyunits.W / pyunits.m / pyunits.K,
        hours_per_shift=8 * pyunits.hr,
        shifts_per_day=3 * (pyunits.day) ** (-1),
        # specific_heat_capacity_insulation1=1.08 * pyunits.kJ / (pyunits.kg * pyunits.K),
        # specific_heat_capacity_metal1=0.468 * pyunits.kJ / (pyunits.kg * pyunits.K),
        # specific_heat_capacity_insulation2=0.9 * pyunits.kJ / (pyunits.kg * pyunits.K),
        # specific_heat_capacity_metal2=0.502416 * pyunits.kJ / (pyunits.kg * pyunits.K),
        operating_days_per_year=336 * pyunits.day / pyunits.year,
        efficiency=0.95 * pyunits.dimensionless,
        utility_rate=0.081 * pyunits.USD_Jan_2024 / (pyunits.kW * pyunits.hr),
        heating_mode=0,
        labor_rate=75 * pyunits.USD_Jan_2024 / pyunits.hr,
        temperature_controller_price=129.00 * pyunits.USD_Jan_2024,
        # number_of_units=1 * pyunits.dimensionless,
        engineering_and_drafting=1000 * pyunits.USD_Jan_2024,
        CE_index_year="Jan_2024",
    ):
        """
        Cost Estimation of an Hydrogen Decrepitation Furnace

        This method calculates the capital cost of a hydrogen decreptitation furnace.

        Args:
        blk: A unit-level costing block where costing variables and
             constraints can be added to.
        ramp_up_time: Time required to reach operating temperature.
        operating_temperature: Optimum temperature at which decrepitation occurs.
        decrepitation_duration: Amount of time the furnace operates at its
                                operating temperature (in hr).
        preparation_time: Time required to setup the furnace.
        cool_down_time: Time required for the decrepitated sample to cool down
                        to room temperature.
        sample_heat_capacity: Heat capacity of the sample (default is sintered NdFeB magnet).
        sample_mass: Mass of the sample.
        sample_volume: Volume of the sample; if not set, will calculate from other inputs.
        sample_density: Density of the sample; default is the density of sintered NdFeB magnet
        chamber_to_sample_ratio: ratio of the volume of furnace chamber to the volume of the sample.
                            Default is set as 2
        length_insulation1: The length of a unit of the first insulation material
                            (default is ceramic fiber), as specified by the vendor.
        width_insulation1: The width of a unit of the first insulation material,
                           as specified by the vendor.
        thickness_insulation1: The thickness of a unit of the first insulation material,
                               as specified by the vendor.
        price_insulation1: The price of a unit of the first insulation material,
                           as specified by the vendor.
        weight_insulation1: The weight of a unit of the first insulation material,
                            as specified by the vendor.
        thermal_cond_insulation_material1: Thermal conductivity of insulation material 1 in W/(m*K)
        insulation material 1 specifications from https://www.grainger.com/product/23AR50?gucid=N:N:PS:Paid:GGL:CSM-2295:4P7A1P:20501231&gad_source=1&gclid=CjwKCAjwnqK1BhBvEiwAi7o0Xxz1EhsrD9nk0yk9SgHwN7fzjC8X249MA5NulaF_JSsTjC290z-HChoCSyYQAvD_BwE&gclsrc=aw.ds
        price_metal1: Price of metal material 1 in USD/kg (default price is for stainless steel 304)
                      from https://mepsinternational.com/gb/en/products/world-stainless-steel-prices as of 01/01/2024
        density_metal1: Density of metal material 1 in kg/m**3 (default density is for stainless steel 304)
        thermal_cond_metal_material1: Thermal conductivity of metal material 1 (in W/(m*K))
        length_insulation2: The length of a unit of the second insulation material
                            (default is fiberglass), as specified by the vendor.
        width_insulation2: The width of a unit of the second insulation material,
                           as specified by the vendor.
        thickness_insulation2: The thickness of a unit of the second insulation material,
                               as specified by the vendor.
        price_insulation2: The price of a unit of the second insulation material,
                           as specified by the vendor.
        weight_insulation2: The weight of a unit of the second insulation material,
                            as specified by the vendor.
        thermal_cond_insulation_material2: Thermal conductivity of insulation material 2 in W/(m*K)
        insulation material 2 specifications from https://www.menards.com/main/building-materials/insulation/insulation-rolls-batts/r13-thermafiber-reg-ultrabatt-trade-mineral-wool-insulation-3-5-x-15-x-47/1200064/p-37325794453617-c-5780.htm
        price_metal2: Price of metal material 2 in USD/kg (default price is for carbon steel)
                      from https://mepsinternational.com/gb/en/products/north-america-steel-prices as of 01/01/2024
        density_metal2: Density of metal material 2 in kg/m**3 (default density is for carbon steel)
        thermal_cond_metal_material2: Thermal conductivity of metal material 2 (in W/(m*K))
        hours_per_shift: Number of hours per shift
        shifts_per_day: Number of shifts per day
        specific_heat_capacity_insulation1: The specific heat capacity of the first insulation material
        specific_heat_capacity_metal1: The specific heat capacity of the first metallic material
        specific_heat_capacity_insulation2: The specific heat capacity of the second insulation material
        specific_heat_capacity_metal2: The specific heat capacity of the second metallic material
        operating_days_per_year: Number of operating days per year
        efficiency: Power usage efficiency
        utility_rate: Unit rate for energy consumption
                        from https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_03 as of 01/01/2024
        heating_mode: mode of operation of the heating system. 0 - electricity, 1 - gas-fired
        labor_rate: Hourly rate of plant operators in project dollar year;
        temperature_controller_price: Unit cost of the temperature controller
                        from https://www.iothrifty.com/products/n480d-low-cost-pid-temperature-controller?variant=42613548187885&gad_source=1&gclid=Cj0KCQjw8J6wBhDXARIsAPo7QA832gIq8ZE3oqLv_OZKHllF8qHXsZ8l48v-urT1a3OCoxQifOp6E1YaAjKJEALw_wcB
                        as of 01/01/2024
        number_of_units: Number of units of the hydrogen decrepitation furnace
                         to be estimated.
        engineering_and_drafting: A flat fee charged for design, drafting, analysis, consultation,
                        and project management.
        """
        if not isinstance(blk.parent_block(), REPMHydrogenDecrepitationFurnace):
            raise TypeError("Parent block is of type ", blk.parent_block(), " and should be of type ",
                            REPMHydrogenDecrepitationFurnace, " to use costing model.")

        if heating_mode != 0 and heating_mode != 1:
            raise TypeError("Valid heating modes are either 0: electric-fired or 1: gas-fired.")

        # Material costs

        blk.price_insulation1 = Param(
            units=blk.costing_package.base_currency, mutable=True
        )
        blk.price_insulation1.set_value(price_insulation1)

        @blk.Expression(doc="Total cost of insulation material 1")
        def material_cost_insulation1(b):
            return b.parent_block().quantity_insulation1[0] * b.price_insulation1

        blk.labor_rate = Param(
            units=blk.costing_package.base_currency / pyunits.hr, mutable=True
        )
        blk.labor_rate.set_value(labor_rate)

        blk.price_metal1 = Param(
            mutable=True,
            units=blk.costing_package.base_currency / pyunits.kg,
            doc="Price of stainless steel 304 (in $/kg)",
        )
        blk.price_metal1.set_value(price_metal1)

        @blk.Expression(doc="Cost of metal material 1")
        def material_cost_metal1(b):
            return (
                pyunits.convert(b.parent_block().weight_metal1[0], to_units=pyunits.kg) * b.price_metal1
            )

        blk.price_insulation2 = Param(
            units=blk.costing_package.base_currency, mutable=True
        )
        blk.price_insulation2.set_value(price_insulation2)

        @blk.Expression(doc="Material cost for insulation material 2")
        def material_cost_insulation2(b):
            return b.parent_block().quantity_insulation2[0] * b.price_insulation2

        blk.price_metal2 = Param(
            mutable=True,
            units=blk.costing_package.base_currency / pyunits.kg,
            doc="Price of carbon steel (in $/kg)",
        )
        blk.price_metal2.set_value(price_metal2)

        @blk.Expression(doc="Cost of metal material 2")
        def material_cost_metal2(b):
            return (
                pyunits.convert(b.parent_block().weight_metal2[0], to_units=pyunits.kg) * b.price_metal2
            )

        # Energy costs

        blk.efficiency = Param(units=pyunits.dimensionless, mutable=True)
        blk.efficiency.set_value(efficiency)

        @blk.Expression(doc="Power rating of the furnace")
        def furnace_power_rating(b):
            return (b.parent_block().heat_furnace_material[0] + b.parent_block().heat_sample_material[0]) / (
                b.parent_block().ramp_up_time * b.efficiency
            )

        blk.hours_per_shift = Param(units=pyunits.hr, mutable=True)
        blk.hours_per_shift.set_value(hours_per_shift)
        blk.shifts_per_day = Param(units=1 / pyunits.day, mutable=True)
        blk.shifts_per_day.set_value(shifts_per_day)
        blk.operating_days_per_year = Param(
            units=pyunits.day / pyunits.year,
            mutable=True,
        )
        blk.operating_days_per_year.set_value(operating_days_per_year)

        @blk.Expression(doc="Number of batches processed annually")
        def total_runs(b):
            return (
                b.hours_per_shift * b.shifts_per_day * b.operating_days_per_year
            ) / b.parent_block().processing_time[0]

        @blk.Expression(doc="Heat duty required annually")
        def annual_heat_duty(b):
            return b.parent_block().total_heat_duty[0] * b.total_runs

        blk.utility_rate = Param(
            units=blk.costing_package.base_currency / (pyunits.kW * pyunits.hr),
            mutable=True,
        )
        blk.utility_rate.set_value(utility_rate)

        blk.OPEX = Var(
            within=PositiveReals,
            initialize=2e6,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Operating expenditure (in USD)",
        )

        @blk.Constraint()
        def operating_cost_eq(b):
            return b.OPEX == pyunits.convert(
                (b.parent_block().total_heat_duty[0] / b.efficiency)
                * b.hours_per_shift * b.shifts_per_day * b.operating_days_per_year
                * b.utility_rate,
                to_units=blk.costing_package.base_currency / blk.costing_package.base_period
                )

        blk.eps = Param(
            initialize=1e-4,
            units=pyunits.USD_Jan_2024,
            doc="eps controls the smoothness of the smoothmax approximation.",
        )

        if heating_mode == 0:

            @blk.Expression(doc="Cost of heating device")
            # Approximate cost correlation based on heater coil price from
            # https://hvacdirect.com/goodman-5-kilowatt-16-200-btu-package-unit-heater-coil-hkp-05c.html
            # https://www.hvacpartsshop.com/l99h008-heating-element-5kw/
            # https://www.alpinehomeair.com/product/air-handlers-electric-furnaces/heater-coils/blueridge/bmahk20
            # https://kaleidoroasters.com/products/kaleido-m10-heating-element?variant=47908078584088
            # https://surpluscityliquidators.com/products/36-kw-electric-heat-kit-7510-ton-air-handlers-208230603.html?srsltid=AfmBOooYvdG8SXAoZufjtif3FUDxYv5JEP180WmuctSv8827yoOfNOWQ0EI
            # as of 01/01/2024
            def cost_heating_device(b):
                return pyunits.convert(
                    (
                        (
                            9.2596
                            * ((pyunits.USD_Jan_2024 * pyunits.s) / pyunits.kJ)
                            * b.furnace_power_rating
                        )
                        + 66.76 * pyunits.USD_Jan_2024
                    ),
                    to_units=b.costing_package.base_currency,
                )

        else:

            @blk.Expression(doc="Cost of heating device")
            # Approximate cost correlation based on gas burner price from
            # https://www.combustion-plus.com/product/maxon-industrial-burner-ople40sunss11dcs-2/
            # https://horizonpfm.com/midco-ec300-economite-gas-burner-300-000-btu-hr/ec300/?setCurrencyId=1&sku=EC300&gad_source=1&gclid=Cj0KCQjw6oi4BhD1ARIsAL6pox2As5QlWrze2Yvyp0BjCFlZ1ydT5OwNJHnM7vTYDpnW79R8nVFP3FIaAsTsEALw_wcB
            # https://burnerparts.com/eclipse-airheat-burner-ah0200-v2-in-stock-ready-to-ship.html?gad_source=4&gclid=Cj0KCQjw6oi4BhD1ARIsAL6pox25C7_IYEYxEoMxizuLW3mVaGdK8acUwlN8JOoeEjhi2zr7nh8BzA0aAvioEALw_wcB
            # https://combustionsupplies.com/products/eclipse-thermjet-tj0100-high-velocity-burner-new-1?variant=43165518463202&currency=USD&utm_medium=product_sync&utm_source=google&utm_content=sag_organic&utm_campaign=sag_organic&gad_source=1&gclid=Cj0KCQjw6oi4BhD1ARIsAL6pox0YMnAnjvsiEjkvVHJwGKSU02d10VVwEzjD-bkaJ7PZIdnmmfbXP7IaAvv9EALw_wcB
            # https://oswaldsupply.com/products/economite-ec300-gas-conversion-burner-min-90-000-max-300-000-btu-hr?variant=33944221450377&kw=EC300&c=Shopping&utm_source=google&utm_medium=cpc&utm_campaign=PerformanceMax-SmartShopping&kw=&ad=&matchtype=&adposition=&c=PerformanceMaxUS&gad_source=1&gclid=Cj0KCQjw6oi4BhD1ARIsAL6pox2c_z-E247Mledh5o9CpHrsJYBodWY8NB5eGzHcnytENYiLSLvCwvQaAtOjEALw_wcB
            def cost_heating_device(b):
                return pyunits.convert(
                    (
                        (
                            11.616
                            * ((pyunits.USD_Jan_2024 * pyunits.s) / pyunits.kJ)
                            * b.furnace_power_rating
                        )
                        - 117.89 * pyunits.USD_Jan_2024
                    ),
                    to_units=b.costing_package.base_currency,
                )

        blk.eps2 = Param(
            initialize=1e-4,
            units=blk.costing_package.base_currency,
            doc="eps2 controls the smoothness of the smoothmax approximation.",
        )

        @blk.Expression(doc="Total labor cost")
        def labor_cost(b):
            return smooth_max(
                3 * pyunits.hr * b.labor_rate,
                (
                    (
                        6.79
                        * (pyunits.hr / (pyunits.m**2))
                        * b.parent_block().furnace_external_surface_area[0]
                    )
                    * b.labor_rate
                ),
                eps=b.eps2,
            )

        blk.temperature_controller_price = Param(
            units=blk.costing_package.base_currency,
            mutable=True,
        )
        blk.temperature_controller_price.set_value(temperature_controller_price)

        blk.engineering_and_drafting = Param(
            units=blk.costing_package.base_currency, mutable=True
        )
        blk.engineering_and_drafting.set_value(engineering_and_drafting)

        @blk.Expression(doc="Overhead Cost")
        def overhead_cost(b):
            return 0.2 * (
                b.material_cost_insulation1
                + b.material_cost_metal1
                + b.material_cost_insulation2
                + b.material_cost_metal2
                + b.temperature_controller_price
                + b.cost_heating_device
                + b.labor_cost
                + b.engineering_and_drafting
            )

        blk.base_cost_per_unit = Var(
            within=PositiveReals,
            initialize=5 * 1e5,
            units=blk.costing_package.base_currency,
            doc="Base cost per unit",
        )

        @blk.Constraint()
        def base_cost_per_unit_eq(b):
            return b.base_cost_per_unit == (
                b.material_cost_insulation1
                + b.material_cost_metal1
                + b.material_cost_insulation2
                + b.material_cost_metal2
                + b.temperature_controller_price
                + b.cost_heating_device
                + b.labor_cost
                + b.engineering_and_drafting
                + b.overhead_cost
            )

        @blk.Expression(doc="Base cost for all installed units")
        def base_cost(b):
            return b.base_cost_per_unit * b.parent_block().config.number_of_units

        blk.capital_cost = Var(
            within=PositiveReals,
            initialize=5 * 1e4,
            units=blk.costing_package.base_currency,
            doc="Capital expenditure",
        )

        @blk.Constraint()
        def capital_cost_eq(b):
            return b.capital_cost == pyunits.convert(
                b.base_cost, to_units=b.costing_package.base_currency
            )

        blk.variable_operating_cost_per_unit = Var(
            initialize=1e5,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Variable operating cost for a furnace unit (in USD)",
        )

        @blk.Constraint()
        def variable_operating_cost_per_unit_eq(b):
            return b.variable_operating_cost_per_unit == pyunits.convert(
                b.OPEX,
                to_units=b.costing_package.base_currency
                / b.costing_package.base_period,
            )

        blk.variable_operating_cost = Var(
            initialize=1e6,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Total variable operating cost for all furnace units",
        )

        @blk.Constraint()
        def variable_operating_cost_eq(b):
            return (
                b.variable_operating_cost
                == b.variable_operating_cost_per_unit * b.parent_block().config.number_of_units
            )
