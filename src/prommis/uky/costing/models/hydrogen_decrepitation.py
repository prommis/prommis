from pyomo.environ import (
    Constraint,
    Param,
    Var,
    units,
    PositiveReals,
    exp,
    Any,
)
from idaes.core import (
    declare_process_block_class,
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)
from idaes.core.util.constants import Constants
from idaes.core.util.math import smooth_max
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def custom_REE_plant_currency_units():
    """
    Define conversion rates for US Dollars based on CE Index.
    """
    register_idaes_currency_units()
    if "USD_Jan_2024" in units._pint_registry:  # pylint: disable=protected-access
        # Assume that custom REE plant units have already been registered
        # Log a message and end
        _log.debug(
            "Custom REE plant currency units (USD_Jan_2024) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        units.load_definitions_from_strings(
            [
                # from https://toweringskills.com/financial-analysis/cost-indices/ as of 01/02/2024
                "USD_Jan_2024 = 500/795.4 * USD_CE500",
            ]
        )


@declare_process_block_class("REEEquipmentCosting")
class REEEquipmentCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index
    custom_REE_plant_currency_units()

    def build_global_params(self):
        # Set the base year for all costs
        self.base_currency = units.USD_Jan_2024
        # Set a base period for all operating costs
        self.base_period = units.year

    def cost_hydrogen_decrepitation_furnace(
        blk,
        ramp_up_time=300 * units.s,
        operating_temperature=443.15 * units.K,
        decrepitation_duration=10800 * units.s,
        preparation_time=3600 * units.s,
        cool_down_time=3600 * units.s,
        sample_heat_capacity=0.44 * units.kJ / (units.kg * units.K),
        sample_mass=None,
        sample_volume=None,
        sample_density=7500 * units.kg / (units.m**3),
        chamber_to_sample_ratio=2 * units.dimensionless,
        length_insulation1=7.62 * units.m,
        width_insulation1=0.6096 * units.m,
        thickness_insulation1=0.0254 * units.m,
        price_insulation1=183.81 * units.USD_Jan_2024,
        weight_insulation1=15.42 * units.kg,
        thermal_cond_insulation_material1=0.33 * units.W / units.m / units.K,
        price_metal1=3.14 * units.USD_Jan_2024 / units.kg,
        density_metal1=7473.57 * units.kg / (units.m**3),
        thermal_cond_metal_material1=13.53 * units.W / units.m / units.K,
        length_insulation2=1.19 * units.m,
        width_insulation2=0.381 * units.m,
        thickness_insulation2=0.0889 * units.m,
        price_insulation2=47.00 * units.USD_Jan_2024,
        weight_insulation2=10.43 * units.kg,
        thermal_cond_insulation_material2=0.069 * units.W / units.m / units.K,
        price_metal2=1.50 * units.USD_Jan_2024 / units.kg,
        density_metal2=7861.09 * units.kg / (units.m**3),
        thermal_cond_metal_material2=45 * units.W / units.m / units.K,
        hours_per_shift=8 * units.hr,
        shifts_per_day=3 * (units.day) ** (-1),
        specific_heat_capacity_insulation1=1.08 * units.kJ / (units.kg * units.K),
        specific_heat_capacity_metal1=0.468 * units.kJ / (units.kg * units.K),
        specific_heat_capacity_insulation2=0.9 * units.kJ / (units.kg * units.K),
        specific_heat_capacity_metal2=0.502416 * units.kJ / (units.kg * units.K),
        operating_days_per_year=336 * units.day / units.year,
        efficiency=0.95 * units.dimensionless,
        electricity_rate=0.081 * units.USD_Jan_2024 / (units.kW * units.hr),
        labor_rate=75 * units.USD_Jan_2024 / units.hr,
        temperature_controller_price=129.00 * units.USD_Jan_2024,
        number_of_units=1 * units.dimensionless,
        engineering_and_drafting=1000 * units.USD_Jan_2024,
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
        electricity_rate: Unit rate for energy consumption
                        from https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_03 as of 01/01/2024
        labor_rate: Hourly rate of plant operators in project dollar year;
        temperature_controller_price: Unit cost of the temperature controller
                        from https://www.iothrifty.com/products/n480d-low-cost-pid-temperature-controller?variant=42613548187885&gad_source=1&gclid=Cj0KCQjw8J6wBhDXARIsAPo7QA832gIq8ZE3oqLv_OZKHllF8qHXsZ8l48v-urT1a3OCoxQifOp6E1YaAjKJEALw_wcB
                        as of 01/01/2024
        number_of_units: Number of units of the hydrogen decrepitation furnace
                         to be estimated.
        engineering_and_drafting: A flat fee charged for design, drafting, analysis, consultation,
                        and project management.
        """

        if sample_mass is None and sample_volume is None:
            raise TypeError("Sample volume and/or sample mass is required")

        blk.furnace_chamber_volume = Var(
            initialize=0.0164,
            units=units.m**3,
            doc="Internal volume of the furnace in cubic meter",
        )
        blk.sample_volume = Param(units=units.m**3, within=Any, mutable=True)
        blk.sample_volume.set_value(sample_volume)
        blk.sample_density = Param(units=units.kg / (units.m**3), mutable=True)
        blk.sample_density.set_value(sample_density)
        blk.chamber_to_sample_ratio = Param(units=units.dimensionless, mutable=True)
        blk.chamber_to_sample_ratio.set_value(chamber_to_sample_ratio)

        if sample_mass is None:
            blk.furnace_chamber_volume.fix(
                blk.chamber_to_sample_ratio * blk.sample_volume
            )
            sample_mass = blk.sample_density * blk.sample_volume
            blk.sample_mass = Param(units=units.kg, mutable=True)
            blk.sample_mass.set_value(sample_mass)

        else:
            blk.sample_mass = Param(units=units.kg, mutable=True)
            blk.sample_mass.set_value(sample_mass)
            blk.furnace_chamber_volume.fix(
                blk.chamber_to_sample_ratio * (blk.sample_mass / blk.sample_density)
            )

        @blk.Expression(doc="Radius of the furnace chamber")
        def radius_chamber(blk):
            return (blk.furnace_chamber_volume / (6 * Constants.pi)) ** (1 / 3)

        @blk.Expression(doc="Length of the furnace chamber")
        def length_chamber(blk):
            # 3 times the diameter
            return 6 * blk.radius_chamber

        blk.heat_loss = Var(
            within=PositiveReals,
            initialize=1e2,
            units=units.W,
            bounds=(1e-6, None),
            doc="Heat loss from the external surface of the furnace (in Watts)",
        )

        blk.stef_bolt_constant = Param(
            initialize=5.67 * 1e-8,
            units=units.W / ((units.m**2) * (units.K**4)),
            doc="Stefan-Boltzmann constant",
        )
        blk.max_temperature = Param(
            initialize=1173.15,
            units=units.K,
            doc="Maximum allowable operating temperature (in Kelvin)",
        )
        blk.temperature_insulation_material1 = Var(
            within=PositiveReals,
            initialize=1050,
            units=units.K,
            bounds=(950.13, 1173.15),
            doc="Temperature of insulation material 1 under steady operating conditions (in Kelvin)",
        )
        blk.heat_loss_constraint1 = Constraint(
            expr=blk.heat_loss
            == 0.70
            * 2
            * Constants.pi
            * (blk.radius_chamber)
            * blk.length_chamber
            * blk.stef_bolt_constant
            * ((blk.max_temperature**4) - (blk.temperature_insulation_material1**4))
        )

        blk.thermal_cond_insulation_material1 = Param(
            mutable=True, units=units.W / units.m / units.K
        )
        blk.thermal_cond_insulation_material1.set_value(
            thermal_cond_insulation_material1
        )
        blk.temperature_metal_material1 = Param(
            initialize=950.13,
            units=units.K,
            within=PositiveReals,
            doc="Temperature of metal 1 under steady operating conditions (in Kelvin)",
        )
        blk.thickness_insulation_material1 = Var(
            within=PositiveReals,
            initialize=0.12,
            units=units.m,
            bounds=(0, None),
            doc="Thickness of insulation material 1 (in meter)",
        )

        blk.relative_thickness_ratio1 = Var(
            within=PositiveReals,
            initialize=2.33,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of insulation material 1",
        )

        blk.heat_loss_constraint2 = Constraint(
            expr=exp(blk.relative_thickness_ratio1)
            == (blk.radius_chamber + blk.thickness_insulation_material1)
            / blk.radius_chamber
        )
        blk.heat_loss_constraint3 = Constraint(
            expr=blk.heat_loss * blk.relative_thickness_ratio1
            == (
                2
                * Constants.pi
                * blk.thermal_cond_insulation_material1
                * blk.length_chamber
                * (
                    blk.temperature_insulation_material1
                    - blk.temperature_metal_material1
                )
            )
        )

        blk.thermal_cond_metal_material1 = Param(
            mutable=True,
            units=units.W / units.m / units.K,
        )
        blk.thermal_cond_metal_material1.set_value(thermal_cond_metal_material1)
        blk.temperature_insulation_material2 = Param(
            initialize=950.00,
            units=units.K,
            doc="Temperature of insulation material 2 under steady operating conditions (in Kelvin)",
        )
        blk.thickness_metal_material1 = Var(
            within=PositiveReals,
            initialize=3 * 1e-3,
            units=units.m,
            bounds=(0, None),
            doc="Thickness of metal material 1 (in meter)",
        )

        blk.relative_thickness_ratio2 = Var(
            within=PositiveReals,
            initialize=1.05,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of metal material 1",
        )

        blk.heat_loss_constraint4 = Constraint(
            expr=exp(blk.relative_thickness_ratio2)
            == (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
            )
            / (blk.radius_chamber + blk.thickness_insulation_material1)
        )
        blk.heat_loss_constraint5 = Constraint(
            expr=blk.heat_loss * blk.relative_thickness_ratio2
            == (
                (
                    (
                        2
                        * Constants.pi
                        * blk.thermal_cond_metal_material1
                        * blk.length_chamber
                    )
                    * (
                        blk.temperature_metal_material1
                        - blk.temperature_insulation_material2
                    )
                )
            )
        )

        blk.thermal_cond_insulation_material2 = Param(
            mutable=True, units=units.W / units.m / units.K
        )
        blk.thermal_cond_insulation_material2.set_value(
            thermal_cond_insulation_material2
        )
        blk.temperature_metal_material2 = Param(
            initialize=333.18,
            units=units.K,
            doc="Temperature of metal 2 under steady operating conditions (in Kelvin)",
        )
        blk.thickness_insulation_material2 = Var(
            within=PositiveReals,
            initialize=0.15,
            units=units.m,
            bounds=(0, None),
            doc="Thickness of insulation material 2 (in meter)",
        )

        blk.relative_thickness_ratio3 = Var(
            within=PositiveReals,
            initialize=1.5,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of insulation material 2",
        )

        blk.heat_loss_constraint6 = Constraint(
            expr=exp(blk.relative_thickness_ratio3)
            == (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
            )
            / (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
            )
        )
        blk.heat_loss_constraint7 = Constraint(
            expr=blk.heat_loss * blk.relative_thickness_ratio3
            == (
                (
                    2
                    * Constants.pi
                    * blk.thermal_cond_insulation_material2
                    * blk.length_chamber
                )
                * (
                    blk.temperature_insulation_material2
                    - blk.temperature_metal_material2
                )
            )
        )

        blk.thermal_cond_metal_material2 = Param(
            mutable=True, units=units.W / units.m / units.K
        )
        blk.thermal_cond_metal_material2.set_value(thermal_cond_metal_material2)
        blk.thickness_metal_material2 = Var(
            within=PositiveReals,
            initialize=4 * 1e-3,
            units=units.m,
            bounds=(0, None),
            doc="Thickness of metal material 2 (in meter)",
        )
        blk.temperature_furnace_ext_surface = Param(
            initialize=333.15,
            units=units.K,
            doc="Temperature of external surface of the furnace under steady operating conditions (in Kelvin)",
        )

        blk.relative_thickness_ratio4 = Var(
            within=PositiveReals,
            initialize=1.005,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Relative thickness ratio of metal material 2",
        )

        blk.heat_loss_constraint8 = Constraint(
            expr=exp(blk.relative_thickness_ratio4)
            == (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
                + blk.thickness_metal_material2
            )
            / (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
            )
        )
        blk.heat_loss_constraint9 = Constraint(
            expr=blk.heat_loss * blk.relative_thickness_ratio4
            == (
                (
                    2
                    * Constants.pi
                    * blk.thermal_cond_metal_material2
                    * blk.length_chamber
                )
                * (
                    blk.temperature_metal_material2
                    - blk.temperature_furnace_ext_surface
                )
            )
        )

        blk.air_heat_transfer_coeff = Param(
            initialize=5,
            units=units.W / ((units.m**2) * units.K),
            doc="Heat transfer coefficient of air (natural/free convection) at ambient temperature",
        )
        blk.ref_temp = Param(
            initialize=298.15,
            units=units.K,
            doc="Reference temperature - 25 degree celsius ",
        )
        blk.heat_loss_constraint10 = Constraint(
            expr=blk.heat_loss
            == blk.air_heat_transfer_coeff
            * 2
            * Constants.pi
            * blk.length_chamber
            * (blk.temperature_furnace_ext_surface - blk.ref_temp)
            * (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
                + blk.thickness_metal_material2
            )
        )

        @blk.Expression(doc="Volume of insulation material 1")
        def volume_insulation1(blk):
            return (
                Constants.pi
                * blk.length_chamber
                * (
                    ((blk.radius_chamber + blk.thickness_insulation_material1) ** 2)
                    - (blk.radius_chamber**2)
                )
            )

        blk.length_insulation1 = Param(units=units.m, mutable=True)
        blk.length_insulation1.set_value(length_insulation1)
        blk.width_insulation1 = Param(units=units.m, mutable=True)
        blk.width_insulation1.set_value(width_insulation1)
        blk.thickness_insulation1 = Param(units=units.m, mutable=True)
        blk.thickness_insulation1.set_value(thickness_insulation1)
        blk.min_quantity = Param(
            initialize=1,
            units=units.dimensionless,
            doc="Lower bound: 1",
        )

        @blk.Expression(doc="Required quantity of insulation material 1")
        def quantity_insulation1(blk):
            return smooth_max(
                blk.min_quantity,
                (
                    blk.volume_insulation1
                    / (
                        blk.length_insulation1
                        * blk.width_insulation1
                        * blk.thickness_insulation1
                    )
                ),
            )

        blk.weight_insulation1 = Param(units=units.kg, mutable=True)
        blk.weight_insulation1.set_value(weight_insulation1)

        @blk.Expression(doc="Total weight of insulation material 1")
        def total_weight_insulation1(blk):
            return blk.quantity_insulation1 * blk.weight_insulation1

        blk.price_insulation1 = Param(
            units=blk.costing_package.base_currency, mutable=True
        )
        blk.price_insulation1.set_value(price_insulation1)

        @blk.Expression(doc="Total cost of insulation material 1")
        def material_cost_insulation1(blk):
            return blk.quantity_insulation1 * blk.price_insulation1

        blk.labor_rate = Param(
            units=blk.costing_package.base_currency / units.hr, mutable=True
        )
        blk.labor_rate.set_value(labor_rate)

        @blk.Expression(
            doc="Internal diameter of the structure formed after attaching metal 1"
        )
        def internal_diameter_metal1(blk):
            return 2 * (blk.radius_chamber + blk.thickness_insulation_material1)

        @blk.Expression(
            doc="External diameter of the structure formed after attaching metal 1"
        )
        def external_diameter_metal1(blk):
            return 2 * (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
            )

        @blk.Expression(doc="Required volume of metal material 1")
        def volume_metal1(blk):
            vol = (
                Constants.pi
                * (
                    (blk.external_diameter_metal1**2)
                    - (blk.internal_diameter_metal1**2)
                )
                * blk.length_chamber
            ) / 4
            return units.convert(vol, to_units=units.inches**3)

        blk.density_metal1 = Param(
            mutable=True,
            units=units.kg / (units.m**3),
            doc="Density of stainless steel 304",
        )
        blk.density_metal1.set_value(density_metal1)

        @blk.Expression(doc="Weight of metal material 1")
        def weight_metal1(blk):
            return blk.volume_metal1 * units.convert(
                blk.density_metal1, to_units=units.pound / (units.inches**3)
            )

        blk.price_metal1 = Param(
            mutable=True,
            units=blk.costing_package.base_currency / units.kg,
            doc="Price of stainless steel 304 (in $/kg)",
        )
        blk.price_metal1.set_value(price_metal1)

        @blk.Expression(doc="Cost of metal material 1")
        def material_cost_metal1(blk):
            return (
                units.convert(blk.weight_metal1, to_units=units.kg) * blk.price_metal1
            )

        # Insulation 2
        @blk.Expression(doc="Volume of insulation material 2")
        def volume_insulation2(blk):
            return (
                Constants.pi
                * blk.length_chamber
                * (
                    (
                        (
                            blk.radius_chamber
                            + blk.thickness_insulation_material1
                            + blk.thickness_metal_material1
                            + blk.thickness_insulation_material2
                        )
                        ** 2
                    )
                    - (
                        (
                            blk.radius_chamber
                            + blk.thickness_insulation_material1
                            + blk.thickness_metal_material1
                        )
                        ** 2
                    )
                )
            )

        blk.length_insulation2 = Param(units=units.m, mutable=True)
        blk.length_insulation2.set_value(length_insulation2)
        blk.width_insulation2 = Param(units=units.m, mutable=True)
        blk.width_insulation2.set_value(width_insulation2)
        blk.thickness_insulation2 = Param(units=units.m, mutable=True)
        blk.thickness_insulation2.set_value(thickness_insulation2)

        @blk.Expression(doc="Required quantity of insulation material 2")
        def quantity_insulation2(blk):
            return smooth_max(
                blk.min_quantity,
                (
                    blk.volume_insulation2
                    / (
                        blk.length_insulation2
                        * blk.width_insulation2
                        * blk.thickness_insulation2
                    )
                ),
            )

        blk.weight_insulation2 = Param(units=units.kg, mutable=True)
        blk.weight_insulation2.set_value(weight_insulation2)

        @blk.Expression(doc="Total weight of insulation material 2")
        def total_weight_insulation2(blk):
            return blk.quantity_insulation2 * blk.weight_insulation2

        blk.price_insulation2 = Param(
            units=blk.costing_package.base_currency, mutable=True
        )
        blk.price_insulation2.set_value(price_insulation2)

        @blk.Expression(doc="Material cost for insulation material 2")
        def material_cost_insulation2(blk):
            return blk.quantity_insulation2 * blk.price_insulation2

        # Metal 2
        @blk.Expression(
            doc="Internal diameter of the structure formed after attaching metal 2"
        )
        def internal_diameter_metal2(blk):
            return 2 * (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
            )

        @blk.Expression(
            doc="Internal diameter of the structure formed after attaching metal 2"
        )
        def external_diameter_metal2(blk):
            return 2 * (
                blk.radius_chamber
                + blk.thickness_insulation_material1
                + blk.thickness_metal_material1
                + blk.thickness_insulation_material2
                + blk.thickness_metal_material2
            )

        @blk.Expression(doc="Required volume of metal material 2")
        def volume_metal2(blk):
            vol = (
                Constants.pi
                * (
                    (blk.external_diameter_metal2**2)
                    - (blk.internal_diameter_metal2**2)
                )
                * blk.length_chamber
            ) / 4
            return units.convert(vol, to_units=units.inches**3)

        blk.density_metal2 = Param(
            mutable=True,
            units=units.kg / (units.m**3),
            doc="Density of carbon steel",
        )
        blk.density_metal2.set_value(density_metal2)

        @blk.Expression(doc="Weight of metal material 2")
        def weight_metal2(blk):
            return blk.volume_metal2 * units.convert(
                blk.density_metal2, to_units=units.pound / (units.inches**3)
            )

        blk.price_metal2 = Param(
            mutable=True,
            units=blk.costing_package.base_currency / units.kg,
            doc="Price of carbon steel (in $/kg)",
        )
        blk.price_metal2.set_value(price_metal2)

        @blk.Expression(doc="Cost of metal material 2")
        def material_cost_metal2(blk):
            return (
                units.convert(blk.weight_metal2, to_units=units.kg) * blk.price_metal2
            )

        @blk.Expression(doc="External surface area of the furnace")
        def furnace_external_surface_area(blk):
            return (
                Constants.pi * blk.external_diameter_metal2 * blk.length_chamber
            ) + ((Constants.pi * (blk.external_diameter_metal2**2)) / 2)

        # energy requirements

        blk.insulation1_max_delta_temp = Param(
            initialize=875,
            units=units.K,
            mutable=True,
            doc="Temperature change of insulation material 1 from reference temperature to steady operating condition",
        )
        blk.metal1_max_delta_temp = Param(
            initialize=651.98,
            units=units.K,
            mutable=True,
            doc="Temperature change of metal 1 from reference temperature to steady operating condition",
        )
        blk.insulation2_max_delta_temp = Param(
            initialize=651.85,
            units=units.K,
            mutable=True,
            doc="Temperature change of insulation material 2 from reference temperature to steady operating condition",
        )
        blk.metal2_max_delta_temp = Param(
            initialize=35.03,
            units=units.K,
            mutable=True,
            doc="Temperature change of metal 2 from reference temperature to steady operating condition",
        )
        blk.specific_heat_capacity_insulation1 = Param(
            units=units.kJ / (units.kg * units.K),
            mutable=True,
        )
        blk.specific_heat_capacity_insulation1.set_value(
            specific_heat_capacity_insulation1
        )
        blk.specific_heat_capacity_metal1 = Param(
            units=units.kJ / (units.kg * units.K),
            mutable=True,
        )
        blk.specific_heat_capacity_metal1.set_value(specific_heat_capacity_metal1)
        blk.specific_heat_capacity_insulation2 = Param(
            units=units.kJ / (units.kg * units.K),
            mutable=True,
        )
        blk.specific_heat_capacity_insulation2.set_value(
            specific_heat_capacity_insulation2
        )
        blk.specific_heat_capacity_metal2 = Param(
            units=units.kJ / (units.kg * units.K),
            mutable=True,
        )
        blk.specific_heat_capacity_metal2.set_value(specific_heat_capacity_metal2)

        @blk.Expression(
            doc="Energy required to raise the temperature of furnace material from room temperature to temperature at steady state"
        )
        def heat_duty1(blk):
            return (
                (
                    blk.total_weight_insulation1
                    * blk.specific_heat_capacity_insulation1
                    * blk.insulation1_max_delta_temp
                )
                + (
                    (units.convert(blk.weight_metal1, to_units=units.kg))
                    * blk.specific_heat_capacity_metal1
                    * blk.metal1_max_delta_temp
                )
                + (
                    blk.total_weight_insulation2
                    * blk.specific_heat_capacity_insulation2
                    * blk.insulation2_max_delta_temp
                )
                + (
                    (units.convert(blk.weight_metal2, to_units=units.kg))
                    * blk.specific_heat_capacity_metal2
                    * blk.metal2_max_delta_temp
                )
            )

        blk.ramp_up_time = Param(units=units.s, mutable=True)
        blk.ramp_up_time.set_value(ramp_up_time)
        blk.efficiency = Param(units=units.dimensionless, mutable=True)
        blk.efficiency.set_value(efficiency)
        blk.decrepitation_duration = Param(units=units.s, mutable=True)
        blk.decrepitation_duration.set_value(decrepitation_duration)

        @blk.Expression(doc="Energy lost in the decrepitation duration")
        def energy_consumption(blk):
            return units.convert(
                (blk.heat_loss * blk.decrepitation_duration) / blk.efficiency,
                to_units=units.kJ,
            )

        blk.operating_temperature = Param(units=units.K, mutable=True)
        blk.operating_temperature.set_value(operating_temperature)
        blk.sample_heat_capacity = Param(
            units=units.kJ / (units.kg * units.K),
            mutable=True,
        )
        blk.sample_heat_capacity.set_value(sample_heat_capacity)

        @blk.Expression(
            doc="Energy required to raise the temperature of sample to operating temperature"
        )
        def heat_duty3(blk):
            return (
                blk.sample_mass
                * blk.sample_heat_capacity
                * (blk.operating_temperature - blk.ref_temp)
            )

        @blk.Expression(doc="Power rating of the furnace")
        def furnace_power_rating(blk):
            return (blk.heat_duty1 + blk.heat_duty3) / (
                blk.ramp_up_time * blk.efficiency
            )

        @blk.Expression(doc="Total heat duty needed")
        def total_heat_duty(blk):
            return blk.heat_duty1 + blk.energy_consumption + blk.heat_duty3

        blk.preparation_time = Param(units=units.s, mutable=True)
        blk.preparation_time.set_value(preparation_time)
        blk.cool_down_time = Param(units=units.s, mutable=True)
        blk.cool_down_time.set_value(cool_down_time)

        @blk.Expression(doc="Total duration of the decrepitation process")
        def processing_time(blk):
            return (
                (units.convert(blk.preparation_time, to_units=units.hr))
                + (units.convert(blk.ramp_up_time, to_units=units.hr))
                + (units.convert(blk.decrepitation_duration, to_units=units.hr))
                + (units.convert(blk.cool_down_time, to_units=units.hr))
            )

        blk.hours_per_shift = Param(units=units.hr, mutable=True)
        blk.hours_per_shift.set_value(hours_per_shift)
        blk.shifts_per_day = Param(units=1 / units.day, mutable=True)
        blk.shifts_per_day.set_value(shifts_per_day)
        blk.operating_days_per_year = Param(
            units=units.day / units.year,
            mutable=True,
        )
        blk.operating_days_per_year.set_value(operating_days_per_year)

        @blk.Expression(doc="Number of batches processed annually")
        def total_runs(blk):
            return (
                blk.hours_per_shift * blk.shifts_per_day * blk.operating_days_per_year
            ) / blk.processing_time

        @blk.Expression(doc="Heat duty required annually")
        def annual_heat_duty(blk):
            return blk.total_heat_duty * blk.total_runs

        blk.electricity_rate = Param(
            units=blk.costing_package.base_currency / (units.kW * units.hr),
            mutable=True,
        )
        blk.electricity_rate.set_value(electricity_rate)

        blk.OPEX = Var(
            within=PositiveReals,
            initialize=2e4,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Operating expenditure (in USD)",
        )

        @blk.Constraint()
        def operating_cost_eq(blk):
            return blk.OPEX == blk.annual_heat_duty * units.convert(
                blk.electricity_rate,
                to_units=blk.costing_package.base_currency / units.kJ,
            )

        blk.eps = Param(
            initialize=1e-4,
            units=units.USD_Jan_2024,
            doc="eps controls the smoothness of the smoothmax approximation.",
        )

        @blk.Expression(doc="Cost of electric heating coil")
        # based on heater coil price from
        # https://hvacdirect.com/goodman-5-kilowatt-16-200-btu-package-unit-heater-coil-hkp-05c.html
        # as of 01/01/2024
        def cost_heating_coil(blk):
            return units.convert(
                smooth_max(
                    85 * units.USD_Jan_2024,
                    (
                        15.50
                        * ((units.USD_Jan_2024 * units.s) / units.kJ)
                        * blk.furnace_power_rating
                    ),
                    eps=blk.eps,
                ),
                to_units=blk.costing_package.base_currency,
            )

        @blk.Expression(doc="Total weight of furnace")
        def furnace_weight(blk):
            return (
                units.convert(blk.total_weight_insulation1, to_units=units.pound)
                + blk.weight_metal1
                + units.convert(blk.total_weight_insulation2, to_units=units.pound)
                + blk.weight_metal2
            )

        @blk.Expression(doc="Volume of the furnace")
        def furnace_volume(blk):
            vol = (
                Constants.pi * (blk.external_diameter_metal2**2) * blk.length_chamber
            ) / 4
            return units.convert(vol, to_units=units.ft**3)

        blk.eps2 = Param(
            initialize=1e-4,
            units=blk.costing_package.base_currency,
            doc="eps2 controls the smoothness of the smoothmax approximation.",
        )

        @blk.Expression(doc="Total labor cost")
        def labor_cost(blk):
            return smooth_max(
                3 * units.hr * blk.labor_rate,
                (
                    (
                        6.79
                        * (units.hr / (units.m**2))
                        * blk.furnace_external_surface_area
                    )
                    * blk.labor_rate
                ),
                eps=blk.eps2,
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
        def overhead_cost(blk):
            return 0.2 * (
                blk.material_cost_insulation1
                + blk.material_cost_metal1
                + blk.material_cost_insulation2
                + blk.material_cost_metal2
                + blk.temperature_controller_price
                + blk.cost_heating_coil
                + blk.labor_cost
                + blk.engineering_and_drafting
            )

        blk.base_cost_per_unit = Var(
            within=PositiveReals,
            initialize=5 * 1e4,
            units=blk.costing_package.base_currency,
            doc="Base cost per unit",
        )
        blk.number_of_units = Param(units=units.dimensionless, mutable=True)
        blk.number_of_units.set_value(number_of_units)

        @blk.Constraint()
        def base_cost_per_unit_eq(blk):
            return blk.base_cost_per_unit == (
                blk.material_cost_insulation1
                + blk.material_cost_metal1
                + blk.material_cost_insulation2
                + blk.material_cost_metal2
                + blk.temperature_controller_price
                + blk.cost_heating_coil
                + blk.labor_cost
                + blk.engineering_and_drafting
                + blk.overhead_cost
            )

        @blk.Expression(doc="Base cost for all installed units")
        def base_cost(blk):
            return blk.base_cost_per_unit * blk.number_of_units

        blk.capital_cost = Var(
            within=PositiveReals,
            initialize=5 * 1e4,
            units=blk.costing_package.base_currency,
            doc="Capital expenditure",
        )

        @blk.Constraint()
        def capital_cost_eq(blk):
            return blk.capital_cost == units.convert(
                blk.base_cost, to_units=blk.costing_package.base_currency
            )

        blk.variable_operating_cost_per_unit = Var(
            initialize=1e5,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Variable operating cost for a furnace unit (in USD)",
        )

        @blk.Constraint()
        def variable_operating_cost_per_unit_eq(blk):
            return blk.variable_operating_cost_per_unit == units.convert(
                blk.OPEX,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )

        blk.variable_operating_cost = Var(
            initialize=1e6,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Total variable operating cost for all furnace units",
        )

        @blk.Constraint()
        def variable_operating_cost_eq(blk):
            return (
                blk.variable_operating_cost
                == blk.variable_operating_cost_per_unit * blk.number_of_units
            )
