#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This model was derived from:
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/costing/unit_models/ion_exchange.py
# and
# https://github.com/watertap-org/watertap/blob/main/watertap/costing/watertap_costing_package.py

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

"""
Flowsheet costing block for diafiltration flowsheet model
"""
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class, register_idaes_currency_units
from idaes.core.util.constants import Constants
from idaes.models.costing.SSLW import SSLWCostingData

from prommis.ion_exchange.costing.ion_exchange_cost_block import (
    IXCostingBlockData,
)

"""
References:
[1] https://www.intratec.us/chemical-markets/sulfuric-acid-price
"""


def build_regenerant_cost_param_block(blk, regenerant):

    # Define cost and purity based on the regenerant type
    if regenerant == "NaCl":
        cost_value = 0.09
        purity_value = 1.0
        year_value = pyo.units.USD_2020
    elif regenerant == "MeOH":
        cost_value = 3.395  # for 100% purity - ICIS
        purity_value = 1
        year_value = pyo.units.USD_2008
    elif regenerant == "NaOH":
        cost_value = 0.59  # for 30% sol'n - iDST
        purity_value = 0.30
        year_value = pyo.units.USD_2020
    elif regenerant == "H2SO4":
        cost_value = 0.119  # from ref[1]
        purity_value = 1
        year_value = pyo.units.USD_2019
    elif regenerant == "HCl":
        cost_value = 0.17  # for 37% sol'n - CatCost v 1.0.4
        purity_value = 0.37
        year_value = pyo.units.USD_2020
    else:
        raise ValueError(f"Unsupported regenerant type: {regenerant}")

    # Register cost and purity parameters
    blk.cost = pyo.Param(
        mutable=True,
        initialize=cost_value,
        doc=f"{regenerant} cost",
        units=year_value / pyo.units.kg,
    )

    blk.purity = pyo.Param(
        mutable=True,
        initialize=purity_value,
        doc=f"{regenerant} purity",
        units=pyo.units.dimensionless,
    )

    costing = blk.costing_package
    costing.register_flow_type(regenerant, blk.cost / blk.purity)


@declare_process_block_class("IXCosting")
class IXCostingData(IXCostingBlockData):
    """
    Costing block for the ion exchange model
    """

    def build_global_params(self):

        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2021

        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Build the common global parameters and process cost at the
        # start. NOTE: The common process costs adds the
        # total_capital_cost and total_operating_cost
        self._build_common_global_params()

        self._build_common_process_costs()

    def build_process_costs(
        self,
    ):
        """This method builds the process-wide costing"""

        self.total_investment_factor = pyo.Param(
            initialize=1,
            domain=pyo.NonNegativeReals,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.maintenance_labor_chemical_factor = pyo.Param(
            initialize=0.03,
            domain=pyo.NonNegativeReals,
            doc="Maintenance-labor-chemical factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )
        self.capital_annualization_factor = pyo.Param(
            initialize=0.1,
            domain=pyo.NonNegativeReals,
            doc="Capital annualization factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )
        self.capital_recovery_factor.expr = self.capital_annualization_factor

        self.maintenance_labor_chemical_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Maintenance-labor-chemical operating cost",
            units=self.base_currency / self.base_period,
        )

        self.total_capital_cost_constraint = pyo.Constraint(
            expr=self.total_capital_cost
            == self.total_investment_factor * self.aggregate_capital_cost
        )

        self.maintenance_labor_chemical_operating_cost_constraint = pyo.Constraint(
            expr=self.maintenance_labor_chemical_operating_cost
            == self.maintenance_labor_chemical_factor * self.aggregate_capital_cost,
            doc="Maintenance-labor-chemical operating cost",
        )
        self.total_fixed_operating_cost = pyo.Expression(
            expr=self.aggregate_fixed_operating_cost
            + self.maintenance_labor_chemical_operating_cost,
            doc="Total fixed operating costs",
        )

        self.total_variable_operating_cost = pyo.Expression(
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
        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == (self.total_fixed_operating_cost + self.total_variable_operating_cost),
            doc="Total operating cost of process per operating period",
        )
        self.total_annualized_cost = pyo.Expression(
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

    def aggregate_costs(self):
        """
        This method aggregates costs from all the unit models and flows
        registered with this FlowsheetCostingBlock and creates aggregate
        variables for these on the FlowsheetCostingBlock that can be used for
        further process-wide costing calculations.

        The following costing variables are aggregated from all the registered
        UnitModelCostingBlocks (if they exist):

        * capital_cost,
        * direct_capital_cost,
        * fixed_operating_cost, and
        * variable_operating_cost

        Additionally, aggregate flow variables are created for all registered
        flow types along with aggregate costs associated with each of these.

        Args:
            None
        """

        super().aggregate_costs()
        c_units = self.base_currency

        @self.Expression(doc="Aggregation Expression for direct capital cost")
        def aggregate_direct_capital_cost(blk):
            e = 0
            for u in self._registered_unit_costing:
                # Allow for units that might only have a subset of cost Vars
                if hasattr(u, "direct_capital_cost"):
                    e += pyo.units.convert(u.direct_capital_cost, to_units=c_units)
                elif hasattr(u, "capital_cost"):
                    raise RuntimeError(
                        f"Models with a capital_cost must also supply a direct_capital_cost. Found unit {u.unit_model} with `capital_cost` but no `direct_capital_cost`."
                    )

            return e

    def cost_ion_exchange(blk):
        """
        Volume-based capital cost for Ion Exchange
        """

        # Define coefficients for the ion exchange pressure vessels,
        # backwash/rinse tank, and regeneration solution tank costed
        # with power equation. Consider col_vol and tank_vol in
        # gallons. [ESR NOTE: Theoriginal terms were Var() but now I
        # converted them to Param()]. Below find the power equation
        # for each equipment:
        # pressure_vessel_cost = A * col_vol ** b
        # bw_tank_cost = A * tank_vol ** b
        # regen_tank_cost = A * tank_vol ** b
        blk.vessel_A_coeff = pyo.Param(
            initialize=1596.499333,
            units=pyo.units.USD_2020,
            doc="Ion exchange pressure vessel cost equation - A coeff., Carbon steel w/ stainless steel internals",
        )
        blk.vessel_b_coeff = pyo.Param(
            initialize=0.459496809,
            units=pyo.units.dimensionless,
            doc="Ion exchange pressure vessel cost equation - b coeff., Carbon steel w/ stainless steel internals",
        )

        blk.backwash_tank_A_coeff = pyo.Param(
            initialize=308.9371309,
            units=pyo.units.USD_2020,
            doc="Ion exchange backwash tank cost equation - A coeff., Steel tank",
        )
        blk.backwash_tank_b_coeff = pyo.Param(
            initialize=0.501467571,
            units=pyo.units.dimensionless,
            doc="Ion exchange backwash tank cost equation - b coeff., Steel tank",
        )

        blk.regen_tank_A_coeff = pyo.Param(
            initialize=57.02158923,
            units=pyo.units.USD_2020,
            doc="Ion exchange regen tank cost equation - A coeff. Stainless steel",
        )

        blk.regen_tank_b_coeff = pyo.Param(
            initialize=0.729325391,
            units=pyo.units.dimensionless,
            doc="Ion exchange regen tank cost equation - b coeff. Stainless steel",
        )

        blk.annual_resin_replacement_factor = pyo.Param(
            initialize=0.05,
            units=pyo.units.year**-1,
            doc="Fraction of ion exchange resin replaced per year, 4-5% of bed volume - EPA",
        )

        blk.hazardous_min_cost = pyo.Param(
            initialize=3240,
            units=pyo.units.USD_2020 / pyo.units.year,
            doc="Min cost per hazardous waste shipment - EPA",
        )

        blk.hazardous_resin_disposal = pyo.Param(
            initialize=347.10,
            units=pyo.units.USD_2020 * pyo.units.ton**-1,
            doc="Hazardous resin disposal cost - EPA",
        )

        blk.hazardous_regen_disposal = pyo.Param(
            initialize=3.64,
            units=pyo.units.USD_2020 * pyo.units.gal**-1,
            doc="Hazardous liquid disposal cost - EPA",
        )

        blk.regen_recycle = pyo.Var(
            initialize=1,
            units=pyo.units.dimensionless,
            doc="Number of cycles the regenerant can be reused before disposal",
        )
        blk.regen_recycle.fix()

        # Add capital and operating costs
        blk.capital_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit capital cost",
        )
        blk.fixed_operating_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit fixed operating cost",
        )

        # Define terms and unit conversions to match cost equations in
        # references
        ix_type = blk.unit_model.ion_exchange_type
        tot_num_col = (
            blk.unit_model.number_columns + blk.unit_model.number_columns_redundant
        )
        col_vol_gal = pyo.units.convert(
            blk.unit_model.column_volume, to_units=pyo.units.gal
        )
        bed_vol_ft3 = pyo.units.convert(
            blk.unit_model.bed_volume, to_units=pyo.units.ft**3
        )

        # Declare variables
        blk.capital_cost_vessel = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for one vessel",
        )
        blk.capital_cost_resin = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for resin for one vessel",
        )
        blk.capital_cost_backwash_tank = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for backwash + rinse solution tank",
        )
        blk.operating_cost_hazardous = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Operating cost for hazardous waste disposal",
        )
        blk.total_pumping_power = pyo.Var(
            initialize=1,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kilowatt,
            doc="Total pumping power required",
        )

        if blk.unit_model.config.regenerant != "single_use":

            # [ESR WIP: Add regenerant cost. For now call the function
            # but eventually use this function to construct a block.]
            build_regenerant_cost_param_block(blk, blk.unit_model.config.regenerant)

            # Add regenerant properties
            blk.regen_soln_dens = pyo.Param(
                initialize=1000,
                units=pyo.units.kg / pyo.units.m**3,
                mutable=True,
                doc="Density of regeneration solution",
            )
            blk.regen_dose = pyo.Param(
                initialize=300,
                units=pyo.units.kg / pyo.units.m**3,
                mutable=True,
                doc="Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]",
            )
            blk.capital_cost_regen_tank = pyo.Var(
                initialize=1e5,
                domain=pyo.NonNegativeReals,
                units=blk.costing_package.base_currency,
                doc="Capital cost for regeneration solution tank",
            )
            blk.flow_mass_regen_soln = pyo.Var(
                initialize=1,
                domain=pyo.NonNegativeReals,
                units=pyo.units.kg / pyo.units.year,
                doc="Regeneration solution flow",
            )

        # Declare dictionary with known costs values for anion and
        # cation resin type. When using a resin, always make sure to
        # select the cost for the right type. TODO: Add resin type in
        # configuration dictionary to select cost automatically.
        resin_cost = {
            "anion": {
                "typeA": {
                    "name": "Strong base polystyrenic gel-type Type II",
                    "cost": 244.00,
                    "units": pyo.units.USD_2021 / pyo.units.ft**3,
                    "reference": "EPA-WBS cost model: https://www.epa.gov/system/files/other-files/2022-03/anion-exchange-ae-.xlsm.xlsm",
                },
            },
            "cation": {
                "typeA": {
                    "name": "Strong acid polystyrenic gel-type",
                    "cost": 209.00,
                    "units": pyo.units.USD_2021 / pyo.units.ft**3,
                    "reference": "EPA-WBS cost model: https://www.epa.gov/system/files/other-files/2022-03/cation-exchange-ce-.xlsm.xlsm",
                },
                "typeB": {
                    "name": "Strong acid polystyrenic macroporous",
                    "cost": 231.00,
                    "units": pyo.units.USD_2021 / pyo.units.ft**3,
                    "reference": "EPA-WBS cost model: https://www.epa.gov/system/files/other-files/2022-03/cation-exchange-ce-.xlsm.xlsm",
                },
                "typeC": {
                    "name": "Chelating Amino Phosphonic",
                    "cost": 198.21,
                    "units": pyo.units.USD_2021
                    / pyo.units.ft**3,  # WIP: originally 2021
                    "reference": "https://www.alibaba.com/product-detail/Chelating-amino-phosphonic-resins-for-second_1600071099994.html",
                },
            },
        }

        if ix_type == "cation":

            # [ESR NOTE: Original costing package had this as a
            # variable but I changed it to be a parameter.]
            resin_type = "typeC"
            resin_cost = blk.cation_exchange_resin_cost = pyo.Param(
                initialize=resin_cost["cation"][resin_type]["cost"],
                units=resin_cost["cation"][resin_type]["units"],
                doc="Cation exchange resin cost per cubic ft",
            )

        else:

            # [ESR WIP: The current ion exchange model is used for the
            # separation of rare earth elements (REEs), which implies
            # a cation exchange type. We commented other types for
            # now. We will consider if we need to integrate the
            # alternative types in future iterations.]

            # resin_cost = blk.anion_exchange_resin_cost = pyo.Var(
            #     initialize=resin_cost["anion"]["typeA"]["cost"],
            #     units=resin_cost["anion"]["typeA"]["units"],
            #     doc="Anion exchange resin cost per cubic ft",
            # )

            raise ConfigurationError(
                "The current ion exchange model is limited to cation exchange methods and alternative techniques are not addressed at this time."
            )

        @blk.Constraint()
        def capital_cost_vessel_constraint(blk):
            return blk.capital_cost_vessel == pyo.units.convert(
                (
                    blk.vessel_A_coeff
                    * (col_vol_gal / pyo.units.gallon) ** blk.vessel_b_coeff
                ),
                to_units=blk.costing_package.base_currency,
            )

        @blk.Constraint()
        def capital_cost_resin_constraint(blk):
            return blk.capital_cost_resin == pyo.units.convert(
                resin_cost * bed_vol_ft3, to_units=blk.costing_package.base_currency
            )

        backwash_tank_vol_expr = (
            blk.unit_model.backwash_flow_rate * blk.unit_model.backwash_time
        )

        if blk.unit_model.config.regenerant == "single_use":

            blk.flow_vol_resin = pyo.Var(
                initialize=1e5,
                bounds=(0, None),
                units=pyo.units.m**3 / blk.costing_package.base_period,
                doc="Volumetric flow of resin per cycle",  # assumes you are only replacing the operational columns, t_cycle = t_breakthru
            )
            blk.single_use_resin_replacement_cost = pyo.Var(
                initialize=1e5,
                bounds=(0, None),
                units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
                doc="Operating cost for using single-use resin (i.e., no regeneration)",
            )

            # [ESR updates: Replace the original "breakthrough_time"
            # variable for the target_breakthrough_time. This new variable
            # is the one included in the multicomponent model.]
            @blk.Constraint()
            def flow_vol_resin_constraint(blk):
                return blk.flow_vol_resin == pyo.units.convert(
                    blk.unit_model.bed_volume_total
                    / blk.unit_model.target_breakthrough_time,
                    to_units=pyo.units.m**3 / blk.costing_package.base_period,
                )

            blk.mass_flow_resin = pyo.units.convert(
                blk.flow_vol_resin * blk.unit_model.resin_density,
                to_units=pyo.units.ton / blk.costing_package.base_period,
            )

        else:

            @blk.Expression()
            def regeneration_tank_vol(blk):
                return pyo.units.convert(
                    blk.unit_model.regen_tank_vol, to_units=pyo.units.gal
                )

            @blk.Constraint()
            def capital_cost_regen_tank_constraint(blk):
                return blk.capital_cost_regen_tank == pyo.units.convert(
                    blk.regen_tank_A_coeff
                    * (blk.regeneration_tank_vol / pyo.units.gallon)
                    ** blk.regen_tank_b_coeff,
                    to_units=blk.costing_package.base_currency,
                )

            backwash_tank_vol_expr += (
                blk.unit_model.rinse_flow_rate * blk.unit_model.rinse_time
            )

        @blk.Expression()
        def backwash_tank_vol(blk):
            return pyo.units.convert(backwash_tank_vol_expr, to_units=pyo.units.gal)

        @blk.Constraint()
        def capital_cost_backwash_tank_constraint(blk):
            return blk.capital_cost_backwash_tank == pyo.units.convert(
                blk.backwash_tank_A_coeff
                * (blk.backwash_tank_vol / pyo.units.gallon)
                ** blk.backwash_tank_b_coeff,
                to_units=blk.costing_package.base_currency,
            )

        # Add cost_factor as a variable and make it equal to
        # TIC. [NOTE: This only covers the case when TIC is
        # available. If TPEC is also available, it can be used here
        # too. Check watertap costing package add_cost_factor for more
        # details.]
        blk.cost_factor = pyo.Var(
            initialize=1.0,
            doc="Cost factor",
        )

        @blk.Constraint()
        def cost_factor_link(blk):
            return blk.cost_factor == blk.costing_package.TIC

        # Add direct capital cost. NOTE: This is needed to add the
        # aggregate costs in function aggregate_costs()
        @blk.Expression()
        def direct_capital_cost(b):
            return blk.capital_cost / blk.cost_factor

        if blk.unit_model.config.regenerant == "single_use":

            @blk.Constraint()
            def capital_cost_constraint(blk):
                return blk.capital_cost == blk.cost_factor * pyo.units.convert(
                    (blk.capital_cost_vessel + blk.capital_cost_resin) * tot_num_col,
                    to_units=blk.costing_package.base_currency,
                )

        else:

            @blk.Constraint()
            def capital_cost_constraint(blk):
                return blk.capital_cost == blk.cost_factor * pyo.units.convert(
                    (
                        (
                            (blk.capital_cost_vessel + blk.capital_cost_resin)
                            * tot_num_col
                        )
                        + blk.capital_cost_backwash_tank
                        + blk.capital_cost_regen_tank
                    ),
                    to_units=blk.costing_package.base_currency,
                )

        if blk.unit_model.config.hazardous_waste:

            if blk.unit_model.config.regenerant == "single_use":

                @blk.Constraint()
                def operating_cost_hazardous_constraint(blk):
                    return blk.operating_cost_hazardous == pyo.units.convert(
                        (blk.mass_flow_resin * blk.hazardous_resin_disposal)
                        + blk.hazardous_min_cost,
                        to_units=blk.costing_package.base_currency
                        / blk.costing_package.base_period,
                    )

            else:

                bed_mass_ton = pyo.units.convert(
                    blk.unit_model.bed_volume * blk.unit_model.resin_density,
                    to_units=pyo.units.ton,
                )

                @blk.Constraint()
                def operating_cost_hazardous_constraint(blk):
                    return blk.operating_cost_hazardous == pyo.units.convert(
                        (bed_mass_ton * tot_num_col * blk.hazardous_resin_disposal)
                        * blk.annual_resin_replacement_factor
                        + pyo.units.convert(
                            blk.flow_mass_regen_soln / blk.regen_soln_dens,
                            to_units=pyo.units.gal / pyo.units.year,
                        )
                        * blk.hazardous_regen_disposal
                        + blk.hazardous_min_cost,
                        to_units=blk.costing_package.base_currency
                        / blk.costing_package.base_period,
                    )

        else:

            blk.operating_cost_hazardous.fix(0)

        if blk.unit_model.config.regenerant == "single_use":

            @blk.Constraint()
            def single_use_resin_replacement_cost_constraint(blk):
                return blk.single_use_resin_replacement_cost == pyo.units.convert(
                    blk.flow_vol_resin * resin_cost,
                    to_units=blk.costing_package.base_currency
                    / blk.costing_package.base_period,
                )

            @blk.Constraint()
            def fixed_operating_cost_constraint(blk):
                return blk.fixed_operating_cost == (
                    blk.single_use_resin_replacement_cost + blk.operating_cost_hazardous
                )

        else:

            @blk.Constraint()
            def fixed_operating_cost_constraint(blk):
                return blk.fixed_operating_cost == (
                    blk.operating_cost_hazardous
                    + pyo.units.convert(
                        (
                            (
                                bed_vol_ft3
                                * tot_num_col
                                * blk.annual_resin_replacement_factor
                                * resin_cost
                            )
                        ),
                        to_units=blk.costing_package.base_currency
                        / blk.costing_package.base_period,
                    )
                )

            @blk.Constraint()
            def flow_mass_regen_soln_constraint(blk):
                return blk.flow_mass_regen_soln == pyo.units.convert(
                    (
                        (blk.regen_dose * blk.unit_model.bed_volume * tot_num_col)
                        / (blk.unit_model.cycle_time)
                    )
                    / blk.regen_recycle,
                    to_units=pyo.units.kg / pyo.units.year,
                )

            blk.costing_package.cost_flow(
                blk.flow_mass_regen_soln, blk.unit_model.config.regenerant
            )

        # [ESR updates: Remove regeneration terms and add them inside "if"
        # statement.]
        power_expr = blk.unit_model.main_pump_power + blk.unit_model.backwash_pump_power

        if blk.unit_model.config.regenerant != "single_use":

            power_expr += (
                blk.unit_model.rinse_pump_power + blk.unit_model.regen_pump_power
            )

        @blk.Constraint()
        def total_pumping_power_constr(blk):
            return blk.total_pumping_power == power_expr

        blk.costing_package.cost_flow(blk.total_pumping_power, "electricity")
