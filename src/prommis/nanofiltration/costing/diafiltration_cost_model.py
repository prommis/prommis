#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Flowsheet costing block for diafiltration flowsheet model
"""

from pyomo.environ import (
    Constraint,
    Expression,
    NonNegativeReals,
    Param,
    Var,
    units,
    value,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class, register_idaes_currency_units
from idaes.core.util.constants import Constants
from idaes.models.costing.SSLW import SSLWCostingData

from prommis.nanofiltration.costing.diafiltration_cost_block import (
    DiafiltrationCostingBlockData,
)


@declare_process_block_class("DiafiltrationCosting")
class DiafiltrationCostingData(DiafiltrationCostingBlockData):
    """
    Costing block for the diafiltration flowsheet
    """

    def build_global_params(self):
        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        # Set the base year for all costs
        self.base_currency = units.USD_2021

        # Set a base period for all operating costs
        self.base_period = units.year

    def build_process_costs(
        self,
    ):
        """
        Builds the process-wide costing
        """
        # initialize the common global parameters
        self._build_common_global_params()

        # add total_capital_cost and total_operating_cost
        self._build_common_process_costs()

        self.factor_total_investment = Param(
            initialize=2,
            domain=NonNegativeReals,
            doc="Total investment factor [investment cost/equipment cost]",
            units=units.dimensionless,
        )
        self.factor_maintenance_labor_chemical = Param(
            initialize=0.03,
            domain=NonNegativeReals,
            doc="Maintenance-labor-chemical factor [fraction of investment cost/year]",
            units=units.year**-1,
        )
        self.factor_capital_annualization = Param(
            initialize=0.1,
            domain=NonNegativeReals,
            doc="Capital annualization factor [fraction of investment cost/year]",
            units=units.year**-1,
        )
        self.capital_recovery_factor.expr = self.factor_capital_annualization

        self.maintenance_labor_chemical_operating_cost = Var(
            initialize=1e3,
            domain=NonNegativeReals,
            doc="Maintenance-labor-chemical operating cost",
            units=self.base_currency / self.base_period,
        )
        self.total_capital_cost_constraint = Constraint(
            expr=self.total_capital_cost
            == self.factor_total_investment * self.aggregate_capital_cost
        )
        self.maintenance_labor_chemical_operating_cost_constraint = Constraint(
            expr=self.maintenance_labor_chemical_operating_cost
            == self.factor_maintenance_labor_chemical * self.total_capital_cost
        )
        self.total_fixed_operating_cost = Expression(
            expr=self.aggregate_fixed_operating_cost
            + self.maintenance_labor_chemical_operating_cost,
            doc="Total fixed operating costs",
        )
        self.total_variable_operating_cost = Expression(
            expr=(
                (
                    self.aggregate_variable_operating_cost
                    + sum(self.aggregate_flow_costs[f] for f in self.used_flows)
                    * self.utilization_factor
                )
                if self.used_flows
                else self.aggregate_variable_operating_cost
            ),
            doc="Total variable operating cost of process per operating period",
        )
        self.total_operating_cost_constraint = Constraint(
            expr=self.total_operating_cost
            == (self.total_fixed_operating_cost + self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period",
        )
        self.total_annualized_cost = Expression(
            expr=(
                self.total_capital_cost * self.capital_recovery_factor
                + self.total_operating_cost
            ),
            doc="Total annualized cost of operation",
        )

    @staticmethod
    def initialize_build(self):
        calculate_variable_from_constraint(
            self.total_capital_cost, self.total_capital_cost_constraint
        )
        calculate_variable_from_constraint(
            self.maintenance_labor_chemical_operating_cost,
            self.maintenance_labor_chemical_operating_cost_constraint,
        )
        calculate_variable_from_constraint(
            self.total_operating_cost, self.total_operating_cost_constraint
        )

    def cost_membranes(
        blk,
        membrane_length,
        membrane_width,
    ):
        """
        Costing method for membranes

        References:
            https://doi.org/10.1016/j.ijggc.2019.03.018

        Args:
            membrane_length: total membrane length (m)
            membrane_width: membrane width (m)
        """

        blk.factor_membrane_replacement = Param(
            initialize=0.2,
            doc="Membrane replacement factor [fraction of membrane replaced/year]",
            units=units.year**-1,
        )
        blk.membrane_cost = Param(
            initialize=50,
            doc="Membrane cost",
            units=units.USD_2021 / (units.meter**2),  # TODO: validate reference year
        )

        # create the capital and operating cost variables
        blk.capital_cost = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit capital cost",
        )
        blk.fixed_operating_cost = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit fixed operating cost",
        )

        # calculate membrane area
        blk.membrane_area = Var(
            initialize=2700,
            domain=NonNegativeReals,
            doc="Membrane area in square meters",
            units=units.m**2,
        )

        @blk.Constraint()
        def membrane_area_equation(blk):
            return blk.membrane_area == units.convert(
                (membrane_length * membrane_width), to_units=units.m**2
            )

        @blk.Constraint()
        def capital_cost_constraint(blk):
            return blk.capital_cost == units.convert(
                (blk.membrane_cost * blk.membrane_area),
                to_units=blk.costing_package.base_currency,
            )

        @blk.Constraint()
        def fixed_operating_cost_constraint(blk):
            return blk.fixed_operating_cost == units.convert(
                (
                    blk.factor_membrane_replacement
                    * blk.membrane_cost
                    * blk.membrane_area
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )

    def cost_membrane_pressure_drop(
        blk,
        water_flux,
        vol_flow_feed,
        vol_flow_perm,
    ):
        """
        Costing method for membrane pressure drop

        Args:
            water_flux: water flux through membrane (m/h)
            vol_flow_feed: volumetric flow rate of feed (m3/h)
            vol_flow_perm: volumetric flow rate of permeate (m3/h)
        """
        blk.hydraulic_permeability = Param(
            initialize=3,
            mutable=True,
            doc="Hydraulic permeability (Lp) of the membrane",
            units=units.L / units.m**2 / units.hr / units.bar,
        )

        if not (hasattr(blk.costing_package, "electricity_cost")):
            blk.electricity_cost = Var(
                initialize=0.141,
                domain=NonNegativeReals,
                doc="Unit cost of electricity",
                units=units.USD_2021 / units.kWh,
            )

        else:
            blk.electricity_cost = Var(
                initialize=value(blk.costing_package.electricity_cost),
                domain=NonNegativeReals,
                doc="Unit cost of electricity",
                units=blk.costing_package.electricity_cost.units,
            )
        blk.electricity_cost.fix()

        blk.variable_operating_cost = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit variable operating cost",
        )

        # calculate pressure drop
        blk.pressure_drop = Var(
            initialize=483,
            domain=NonNegativeReals,  # we expect a positive value for pressure drop (Pin-Pout, where Pout<Pin)
            doc="Pressure drop over the membrane",
            units=units.psi,
        )

        Lp = units.convert(
            blk.hydraulic_permeability,
            to_units=units.m**3 / units.m**2 / units.hr / units.bar,
        )

        @blk.Constraint()
        def pressure_drop_equation(blk):
            return blk.pressure_drop == units.convert(
                (water_flux / Lp), to_units=units.psi
            )

        # calculate specific energy consumption
        blk.SEC = Var(
            initialize=3,
            domain=NonNegativeReals,
            doc="Specific energy consumption of feed pump",
            units=units.kWh / units.m**3,
        )

        dP = units.convert(blk.pressure_drop, to_units=units.Pa)

        @blk.Constraint()
        def SEC_equation(blk):
            return blk.SEC == units.convert(
                (vol_flow_feed * dP / vol_flow_perm), to_units=units.kWh / units.m**3
            )

        @blk.Constraint()
        def variable_operating_cost_constraint(blk):
            return blk.variable_operating_cost == units.convert(
                (blk.SEC * vol_flow_perm * blk.electricity_cost),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )

    def cost_pump(blk, inlet_pressure, outlet_pressure, inlet_vol_flow):
        """
        Costing method for pumps.

        References:
            https://doi.org/10.1016/j.memsci.2015.04.065
            Volk, Michael. Pump characteristics and applications. CRC Press, 2013.
            Moran, Seán. "Pump Sizing: Bridging the Gap Between Theory and Practice."
                The Best of Equipment Series (2016): 3.
            https://www.bls.gov/regions/midwest/data/averageenergyprices_selectedareas_table.htm

        Args:
            inlet_pressure: pressure of inlet stream to pump (Pa)
            outlet_pressure: pressure of outlet stream from pump (psi)
            inlet_vol_flow: volumetric flow rate of inlet stream to pump (m3/h)
        """
        blk.density = Param(
            initialize=1000,
            doc="Operating fluid density",
            units=units.kg / units.m**3,
        )
        blk.specific_gravity = Param(
            initialize=1,
            doc="Operating fluid specific gravity",
            units=units.dimensionless,
        )
        blk.pump_correlation_factor = Param(
            initialize=622.59,
            doc="Pump correlation factor (constant)",
            units=units.USD_1996 / (units.kPa * units.m**3 / units.hr) ** 0.39,
        )
        blk.pump_exponential_factor = Param(
            initialize=0.39,
            doc="Pump correlation factor (exponential)",
            units=units.dimensionless,
        )
        blk.pump_head_factor = Param(
            initialize=2.31,
            doc="Pump head factor",
            units=units.ft / units.psi,
        )
        blk.pump_power_factor = Param(
            initialize=3.6 * 10 ** (6),
            doc="Pump power factor",
            units=units.dimensionless,
        )
        blk.pump_efficiency = Param(
            initialize=0.7,
            doc="Pump efficiency",
            units=units.dimensionless,
        )
        if not (hasattr(blk.costing_package, "electricity_cost")):
            blk.electricity_cost = Var(
                initialize=0.141,
                domain=NonNegativeReals,
                doc="Unit cost of electricity",
                units=units.USD_2021 / units.kWh,
            )

        else:
            blk.electricity_cost = Var(
                initialize=value(blk.costing_package.electricity_cost),
                domain=NonNegativeReals,
                doc="Unit cost of electricity",
                units=blk.costing_package.electricity_cost.units,
            )
        blk.electricity_cost.fix()

        # create the capital and operating cost variables
        blk.capital_cost = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit capital cost",
        )
        blk.variable_operating_cost = Var(
            initialize=1e5,
            domain=NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit variable operating cost",
        )

        @blk.Constraint()
        def capital_cost_constraint(blk):
            return blk.capital_cost == units.convert(
                blk.pump_correlation_factor
                * (inlet_vol_flow * inlet_pressure) ** blk.pump_exponential_factor,
                to_units=blk.costing_package.base_currency,
            )

        # calculate the pump head: pump Ref [1] Eqn 1.1
        blk.pump_head = Var(
            initialize=10,
            domain=NonNegativeReals,
            doc="Pump head in meters",
            units=units.m,
        )

        @blk.Constraint()
        def pump_head_equation(blk):
            return blk.pump_head == units.convert(
                (outlet_pressure * blk.pump_head_factor / blk.specific_gravity),
                to_units=units.m,
            )

        # calculate the pump power: pump Ref [2] Eqn 7
        blk.pump_power = Var(
            initialize=10,
            domain=NonNegativeReals,
            doc="Pump power in kWh required for the operational period",
            units=units.kWh,
        )

        grav_constant = units.convert(
            Constants.acceleration_gravity, to_units=units.m / units.hr**2
        )

        @blk.Constraint()
        def pump_power_equation(blk):
            return blk.pump_power == units.convert(
                (
                    units.convert(
                        (
                            inlet_vol_flow
                            * blk.density
                            * grav_constant
                            * blk.pump_head
                            / blk.pump_power_factor
                            / blk.pump_efficiency
                        ),
                        to_units=units.kW,
                    )
                    * blk.costing_package.base_period  # per one year
                ),
                to_units=units.kWh,
            )

        @blk.Constraint()
        def variable_operating_cost_constraint(blk):
            return blk.variable_operating_cost == units.convert(
                blk.pump_power
                * blk.electricity_cost
                / blk.costing_package.base_period,  # per one year
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )

    def cost_precipitator(
        blk,
        precip_volume,
        precip_headspace=1.2,
    ):
        """
        Costing method for precipitator unit. Assumes these are horizontal vessels
        made from 1.25 in thick carbon steel, includes platforms and ladders, and each
        instance is one unit (default args).

        References:
            residence time:
                https://www.sciencedirect.com/science/article/pii/S0304386X19309806
                https://onlinelibrary.wiley.com/doi/full/10.1002/ceat.201700667
                https://pubs.acs.org/doi/full/10.1021/acs.iecr.1c04876
                https://www.sciencedirect.com/science/article/pii/S0304386X01002134
            vessel dimensions:
                https://www.accessengineeringlibrary.com/content/book/9781260455410/back-matter/appendix1?implicit-login=true

        Args:
            precip_volume: volume of the precipitator as calculated by the unit model (m3)
            precip_headaspace: precipitator headspace percentage; default value is 20%
        """

        # calculate the volume needed
        blk.volume_capacity = Var(
            initialize=120,
            domain=NonNegativeReals,
            doc="Volume requirement of precipitator vessel",
            units=units.m**3,
        )

        @blk.Constraint()
        def volume_capacity_equation(blk):
            return blk.volume_capacity == units.convert(
                (precip_headspace * precip_volume), to_units=units.m**3
            )

        # include a length and diameter constraint
        blk.precipitator_diameter = Var(
            initialize=6,
            bounds=(1, 12),
            doc="Diameter of the precipitator vessel",
            units=units.ft,
        )
        blk.precipitator_diameter.fix()
        blk.precipitator_length = Var(
            initialize=8,
            bounds=(3, 10),
            doc="Length of the precipitator vessel",
            units=units.ft,
        )

        @blk.Constraint()
        def diameter_length_ratio_equation(blk):
            """
            Coefficients come from literature source noted in above docstring
            """
            return units.convert(blk.precipitator_length, to_units=units.inch) == (
                units.convert(
                    (
                        blk.volume_capacity
                        - units.convert(
                            (
                                2
                                * (0.954 * units.gal / units.ft**3)
                                * (
                                    units.convert(
                                        blk.precipitator_diameter, to_units=units.inch
                                    )
                                    / (12 * units.inch / units.ft)
                                )
                                ** 3
                            ),
                            to_units=units.m**3,
                        )
                    )
                    / units.convert(
                        (
                            (0.0034 * units.gal / units.inch**3)
                            * units.convert(
                                blk.precipitator_diameter, to_units=units.inch
                            )
                            ** 2
                        ),
                        to_units=units.m**2,
                    ),
                    to_units=units.inch,
                )
            )

        SSLWCostingData.cost_vessel(
            blk,
            vertical=False,  # horizontal vessel
            vessel_diameter=blk.precipitator_diameter,
            vessel_length=blk.precipitator_length,
        )
