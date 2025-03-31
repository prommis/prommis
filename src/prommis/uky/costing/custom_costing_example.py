#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Example of a custom cost model compatible with the REE costing library.
"""

__author__ = "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, and M. Zamarripa)"
__version__ = "1.0.0"

from pyomo.environ import NonNegativeReals, Param, Var
from pyomo.environ import units as pyunits

from idaes.core import FlowsheetCostingBlockData, declare_process_block_class


@declare_process_block_class("CustomCosting")
class CustomCostingData(FlowsheetCostingBlockData):
    # Register custom currency units used in the custom costing model
    pyunits.load_definitions_from_strings(
        [
            "USD_custom = 500/700 * USD_CE500",
        ]
    )

    # define global parameters - base cost year and period
    def build_global_params(self):

        # Set the base year for all costs
        self.base_currency = pyunits.USD_2022
        # Set a base period for all operating costs
        self.base_period = pyunits.year

    # custom costing method
    def cost_custom_vessel(
        blk,  # when the costing block is built, blk will be the costing block itself
        volume_per_unit=1000 * pyunits.m**3,
        material="carbonsteel",
        water_injection_rate_per_unit=1 * pyunits.m**3 / pyunits.s,
        number_of_units=1,
    ):

        # make parameter for number of units
        blk.number_of_units = Param(initialize=number_of_units, mutable=False)

        # define the bare erected cost per unit

        material_factor_dict = (
            {  # material factors for each valid choice of shell material
                "carbonsteel": 1,
                "stainlessteel": 1.5,
                "aluminum": 2,
            }
        )

        blk.capital_cost_per_unit = Var(
            initialize=1000,
            units=blk.costing_package.base_currency,
            domain=NonNegativeReals,
            bounds=(0, None),
        )

        @blk.Constraint()
        def capital_cost_per_unit_eq(blk):
            # cost equation CAPITAL_COST = REF_COST * (VOLUME / REF_VOLUME)**0.6
            ref_cost = (
                10000 * pyunits.USD_custom
            )  # reference cost is in reference currency units
            ref_volume = 1000 * pyunits.m**3
            return blk.capital_cost_per_unit == pyunits.convert(
                material_factor_dict[material]
                * ref_cost
                * (volume_per_unit / ref_volume) ** 0.6,
                to_units=blk.costing_package.base_currency,  # convert to costing block base currency
            )

        # create a variable capital_cost that the REE Costing Framework can look for
        blk.capital_cost = Var(
            initialize=1000,
            units=blk.costing_package.base_currency,  # define in costing block base currency
            domain=NonNegativeReals,
            bounds=(0, None),
        )

        @blk.Constraint()
        def capital_cost_constraint(blk):
            return blk.capital_cost == blk.capital_cost_per_unit * blk.number_of_units

        # define the fixed costs

        # create a variable fixed_operating_cost that the REE Costing Framework can look for
        blk.fixed_operating_cost = Var(
            initialize=1000,
            units=blk.costing_package.base_currency
            / blk.costing_package.base_period,  # define in costing block base currency
            domain=NonNegativeReals,
            bounds=(0, None),
        )

        # set fixed OPEX = 5% of CAPEX
        # this is an arbitrary choice for this example, and not a general rule
        # in practice users should create Param and Var as needed for their model
        blk.fixed_OPEX_coefficient = Param(initialize=0.05, mutable=False)

        @blk.Constraint()
        def fixed_operating_cost_constraint(blk):
            return blk.fixed_operating_cost == pyunits.convert(
                blk.fixed_OPEX_coefficient * blk.capital_cost / pyunits.year,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )

        # define the variable costs

        blk.variable_operating_cost_per_unit = Var(
            initialize=1000,
            units=blk.costing_package.base_currency
            / blk.costing_package.base_period,  # define in costing block base currency / time
            domain=NonNegativeReals,
            bounds=(0, None),
        )

        # set variable OPEX = $0.0019 per gallon of water injected using the reference year (USD_custom)
        # this is an arbitrary choice for this example, and not a general rule
        blk.variable_opex_price = Param(
            initialize=0.00190,
            units=pyunits.USD_custom
            / pyunits.gal,  # define in reference currency units
            mutable=False,
        )

        @blk.Constraint()
        def variable_operating_cost_per_unit_eq(blk):
            return blk.variable_operating_cost_per_unit == pyunits.convert(
                blk.variable_opex_price * water_injection_rate_per_unit,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,  # define in costing block base currency / time
            )

        # create a variable variable_operating_cost that the REE Costing Framework can look for
        blk.variable_operating_cost = Var(
            initialize=1000,
            units=blk.costing_package.base_currency
            / blk.costing_package.base_period,  # define in costing block base currency / time
            domain=NonNegativeReals,
            bounds=(0, None),
        )

        @blk.Constraint()
        def variable_operating_cost_constraint(blk):
            return (
                blk.variable_operating_cost
                == blk.variable_operating_cost_per_unit * blk.number_of_units
            )
