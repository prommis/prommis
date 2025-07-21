#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
REE costing library
This method leverages NETL costing capabilities.

Other methods:

    - get_fixed_OM_costs() to cost fixed O&M costs
    - get_variable_OM_costs() to cost variable O&M costs
    - costing_initialization() to initialize costing blocks
    - display_total_plant_costs() to display total plant cost (TPC)
    - display_bare_erected_costs() to display BEC costs
    - get_total_BEC() to display the total BEC of the entire flowsheet
    - display_flowsheet_cost() to display flowsheet cost
    - calculate_REE_costing_bounds() to provide an estimate of costing bounds
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta,L. Deng, and M. Zamarripa)"
)
__version__ = "1.0.0"

import textwrap
from sys import stdout

from pyomo.common.config import ConfigValue, ListOf
from pyomo.common.dependencies import attempt_import
from pyomo.core.base.component import Component
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.units_container import InconsistentUnitsError, UnitsError
from pyomo.environ import ConcreteModel, Expression, Param, Reference, Var, log10
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import (
    FlowsheetBlock,
    FlowsheetCostingBlockData,
    UnitModelCostingBlock,
    declare_process_block_class,
    register_idaes_currency_units,
)
from idaes.core.util.math import smooth_max
from idaes.core.util.tables import stream_table_dataframe_to_string

from pandas import DataFrame

from prommis.uky.costing.costing_dictionaries import (
    load_default_resource_prices,
    load_default_sale_prices,
    load_REE_costing_dictionary,
)

_, watertap_costing_available = attempt_import("watertap.costing")
if watertap_costing_available:
    from watertap.core import ZeroOrderBaseData
    from watertap.costing import WaterTAPCosting
    from watertap.costing.zero_order_costing import ZeroOrderCosting

_log = idaeslog.getLogger(__name__)

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------


# example of adding custom units; taken from pp costing, can change later
def custom_REE_plant_currency_units():
    """
    Define conversion rates for US Dollars based on CE Index.
    """
    register_idaes_currency_units()
    if (
        "USD_2022" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_2025" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_UKy_2019" in pyunits._pint_registry  # pylint: disable=protected-access
    ):
        # Assume that custom REE plant units have already been registered
        # Log a message and end
        _log.info(
            "Custom REE plant currency units (USD_2022, USD_2025, USD_UKy_2019) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                # from UKy 2019 report
                "USD_UKy_2019 = 500/609.495 * USD_CE500",
                # from https://toweringskills.com/financial-analysis/cost-indices/ as of 9/26/2023
                "USD_2022 = 500/816.0 * USD_CE500",
                # from UKy 2023 report
                "USD_2025 = 500/815.59 * USD_CE500",
            ]
        )


@declare_process_block_class("QGESSCosting")
class QGESSCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index

    custom_REE_plant_currency_units()

    # set CONFIG

    CONFIG = FlowsheetCostingBlockData.CONFIG()

    CONFIG.declare(
        "discount_percentage",
        ConfigValue(
            default=None,
            domain=float,
            description="Rate of return used to discount future cash flows "
            "back to their present value. The value should be a percentage, "
            "for example 10 for a 10% discount. The NETL QGESS recommends "
            "setting the discount rate as the calculated after-tax weighted "
            "average cost of capital (ATWACC).",
        ),
    )
    CONFIG.declare(
        "plant_lifetime",
        ConfigValue(
            default=None,
            domain=float,
            description="Length of operating period in years.",
        ),
    )
    CONFIG.declare(
        "total_capital_cost",
        ConfigValue(
            default=None,
            description="Value for total capital cost (including equipment, "
            "installation, and other plant costs); ignored if no value is "
            "passed. Can be a Var, Param, or Expression with currency units, "
            "or can specify a cost year in the cost_year argument.",
        ),
    )
    CONFIG.declare(
        "annual_operating_cost",
        ConfigValue(
            default=None,
            description="Value for total operating cost; ignored if no value "
            "is passed. If a Var, Param, or Expression, must have the same "
            "units as the Var, Param, or Expression provided for total_capital_cost.",
        ),
    )
    CONFIG.declare(
        "annual_revenue",
        ConfigValue(
            default=None,
            description="Value for total revenue; ignored if no value is passed. "
            "If a Var, Param, or Expression, must have the same units "
            "as the Var, Param, or Expression provided for total_capital_cost.",
        ),
    )
    CONFIG.declare(
        "cost_year",
        ConfigValue(
            default=None,
            domain=str,
            description="Assumed project start year for costs, which is the "
            "basis for NPV results; ignored if no value is passed.",
        ),
    )
    CONFIG.declare(
        "has_capital_expenditure_period",
        ConfigValue(
            default=False,
            domain=bool,
            description="True/false flag whether a capital expenditure period "
            "occurs.",
        ),
    )
    CONFIG.declare(
        "capital_expenditure_percentages",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="A list of values that sum to 100 representing how "
            "capital costs are spread over a capital expenditure period; for "
            "example, an input of [10, 60, 30] is parsed as a 3-year period "
            "where capital costs are spread as 10% in year 1, 60% in year 2, "
            "and 30% in year 3. The capital period precedes the operating "
            "period, for example an input of [100] means that 100% of capital "
            "expenses occur in the year preceding the operating period. Set "
            "to None to indicate no expenditure period, which means that all "
            "capital expenses occur at the start of the plant lifetime t=0.",
        ),
    )
    CONFIG.declare(
        "capital_escalation_percentage",
        ConfigValue(
            default=3.6,
            domain=float,
            description="Rate at which capital costs escalate during the "
            "capital expenditure period. The value should be a percentage, "
            "for example 10 for a 10% escalation rate. Set to 0 to indicate "
            "there is no cost escalation in the expenditure period.",
        ),
    )
    CONFIG.declare(
        "capital_loan_interest_percentage",
        ConfigValue(
            default=6,
            domain=float,
            description="Interest rate for capital equipment loan repayment."
            "The value should be a percentage, for example 10 for a 10% "
            "interest rate.",
        ),
    )
    CONFIG.declare(
        "capital_loan_repayment_period",
        ConfigValue(
            default=10,
            domain=float,
            description="Length of loan repayment period in years.",
        ),
    )
    CONFIG.declare(
        "debt_percentage_of_CAPEX",
        ConfigValue(
            default=50,
            domain=float,
            description="Percentage of CAPEX financed by debt; ignored if "
            "debt_expression is not None. The value should be a percentage, "
            "for example 10 for a debt corresponding to 10% of the CAPEX. Set "
            "to zero to indicate no loans are taken out on capital.",
        ),
    )
    CONFIG.declare(
        "debt_expression",
        ConfigValue(
            default=None,
            description="Set the value or expression to calculate total debt. "
            "Enter a Pyomo expression of flowsheet or cost model variables. "
            "The expression should have currency units.",
        ),
    )
    CONFIG.declare(
        "operating_inflation_percentage",
        ConfigValue(
            default=3,
            domain=float,
            description="Inflation rate for operating costs during the "
            "operating period. The value should be a percentage, for example "
            "10 for a 10% inflation rate. Set to 0 to indicate no inflation.",
        ),
    )
    CONFIG.declare(
        "revenue_inflation_percentage",
        ConfigValue(
            default=3,
            domain=float,
            description="Inflation rate for revenue during the operating "
            "period. The value should be a percentage, for example 10 for a "
            "10% inflation rate. Set to 0 to indicate no inflation.",
        ),
    )

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyunits.USD_2021
        # Set a base period for all operating costs
        self.base_period = pyunits.year

    def build_process_costs(
        self,
        # arguments related to installation costs
        total_purchase_cost=None,
        Lang_factor=None,  # default percentages are effective Lang_factor of 2.97
        piping_materials_and_labor_percentage=20,
        electrical_materials_and_labor_percentage=20,
        instrumentation_percentage=8,
        plants_services_percentage=10,
        process_buildings_percentage=40,
        auxiliary_buildings_percentage=15,
        site_improvements_percentage=10,
        equipment_installation_percentage=17,
        field_expenses_percentage=12,
        project_management_and_construction_percentage=30,
        process_contingency_percentage=15,
        # arguments related to Fixed OM costs
        labor_types=None,
        labor_rate=None,
        labor_burden=25,
        operators_per_shift=None,
        hours_per_shift=8,
        shifts_per_day=3,
        operating_days_per_year=336,
        pure_product_output_rates=None,
        mixed_product_output_rates=None,
        mixed_product_sale_price_realization_factor=0.65,
        sale_prices=None,
        # arguments related to total owners costs
        land_cost=None,
        resources=None,
        rates=None,
        prices=None,
        fixed_OM=True,
        variable_OM=False,
        feed_input=None,
        efficiency=0.85,
        chemicals=None,
        additional_chemicals_cost=None,
        waste=None,
        additional_waste_cost=None,
        transport_cost_per_ton_product=None,
        recovery_rate_per_year=None,
        consider_taxes=False,
        income_tax_percentage=26,
        mineral_depletion_percentage=14,
        production_incentive_percentage=10,
        royalty_charge_percentage_of_revenue=6.5,
        CE_index_year="2021",
        watertap_blocks=None,
        calculate_NPV=False,
    ):
        """
        This method builds process-wide costing, including fixed and variable
        operating & maintenance costs, costs of production, cost of
        electricity and cost of capture.

        If individual percentages are provided (defaults to values above), this
        method creates constraints for the following plantwide costs:
        1. Total ancillary
        2. Piping materials and labor ancillary
        3. Electrical materials and labor ancillary
        4. Instrumentation ancillary
        5. Plant services ancillary
        6. Total buildings
        7. Process buildings
        8. Auxiliary buildings
        9. Site improvements buildings
        10. Total engineering procurement and construction management (EPCM)
        11. Equipment installation EPCM
        12. Field expenses EPCM
        13. Project management and construction EPCM
        14. Total contingency
        15. Process contingency
        16. (space for more costs, including contingency costs, in the future)

        These costs apply to the project as a whole and are scaled based on the
        total TPC.

        Args:
            total_purchase_cost: user-defined value for the total equipment
                purchase cost (not including installation or other plant costs).
                To use as the total plant cost, including installation, set the
                Lang_factor to 1.
            Lang_factor: single multiplicative factor to estimate installation
                costs; defaults to None and method will use percentages. The
                default percentages yield an effective Lang factor of 2.97.
            piping_materials_and_labor_percentage: Piping, materials and labor
                costs as a percentage of the total plant cost. The value
                should be a percentage, for example 10 for 10%. If Lang_factor
                is not None, this value will not be used.
            electrical_materials_and_labor_percentage: Electrical, materials
                and labor costs as a percentage of the total plant cost. The
                value should be a percentage, for example 10 for 10%. If
                Lang_factor is not None, this value will not be used.
            instrumentation_percentage: Instrumentation costs as a percentage
                of the total plant cost. The value should be a percentage, for
                example 10 for 10%. If Lang_factor is not None, this value
                will not be used.
            plants_services_percentage: Plant services costs as a percentage
                of the total plant cost. The value should be a percentage, for
                example 10 for 10%. If Lang_factor is not None, this value
                will not be used.
            process_buildings_percentage: Process buildings costs as a
                percentage of the total plant cost. The value should be a
                percentage, for example 10 for 10%. If Lang_factor is not None,
                this value will not be used.
            auxiliary_buildings_percentage: Auxiliary buildings costs as a
                percentage of the total plant cost. The value should be a
                percentage, for example 10 for 10%. If Lang_factor is not None,
                this value will not be used.
            site_improvements_percentage: Site improvements costs as a
                percentage of the total plant cost. The value should be a
                percentage, for example 10 for 10%. If Lang_factor is not None,
                this value will not be used.
            equipment_installation_percentage: Equipment installation costs as
                a percentage of the total plant cost. The value should be a
                percentage, for example 10 for 10%. If Lang_factor is not
                None, this value will not be used.
            field_expenses_percentage: Field expenses costs as a percentage of
                the total plant cost. The value should be a percentage, for
                example 10 for 10%. If Lang_factor is not None, this value
                will not be used.
            project_management_and_construction_percentage: Project management
                and construction costs as a percentage of the total plant cost.
                The value should be a percentage, for example 10 for 10%. If
                Lang_factor is not None, this value will not be used.
            process_contingency_percentage: Process contingency costs as a
                percentage of the total plant cost. The value should be a
                percentage, for example 10 for 10%. If Lang_factor is not None,
                this value will not be used.
            total_purchase_cost: The BEC that will be used to determine
                installation and fixed O&M costs. If the value is None, the
                function will try to use the BEC calculated from the individual
                units. This quantity should be a Pyomo Var or Param that will
                contain the BEC value.
            labor_type: list of types of operators present in plant; assumed to
                correspond with labor rate and operator per shift lists
            labor_rate: hourly rate of plant operators in project dollar year;
                defined as list corresponding to different operator types
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses; assumed constant for all operator types. The
                value should be a percentage, for example 10 for 10%.
            operators_per_shift: number of operators per shift; defined as list
                of operators per shift for each operator type
            hours_per_shift: number of hours per shift
            shifts_per_day: number of shifts per day
            operating_days_per_year: number of operating days per year
            feed_input: rate of feedstock input
            pure_product_output_rates: dict of production rates of each REE pure product
            mixed_product_output_rates: dict of production rates of each REE in the mixed product
            mixed_product_sale_price_realization_factor: multiplicative factor for selling impure products
            sale_prices: list setting sale prices of products
            land_cost: Expression, Var or Param to calculate land costs
            resources: list setting resources to cost
            rates: list setting flow rates of resources
            prices: list setting prices of resources
            fixed_OM: True/False flag for calculating fixed O&M costs
            variable_OM: True/False flag for calculating variable O&M costs
            efficiency: power usage efficiency, or fixed motor/distribution efficiency
            chemicals: string setting chemicals type for chemicals costs
            additional_chemicals_cost: Expression, Var or Param to calculate additional chemical costs.
            waste: string setting waste type for waste costs
            additional_waste_cost: Expression, Var or Param to calculate additional waste disposal costs.
            recovery_rate_per_year: Var or value to use for rate of REE recovered, in units
                of mass/year
            consider_taxes: True/False flag for calculating net tax owed. Defaults to False.
            income_tax_percentage: combined federal and state income tax percentage,
                usually between 26 - 40%. Here, it defaults to 26%.
            mineral_depletion_percentage: fixed tax deduction percentage for mineral depletion based on
                the type of mineral recovered, defaults to 14% of gross income excluding royalties
                as reported in the UKy report.
            production_incentive_percentage: tax deduction percentage for producing critical minerals,
                defaults to 10% of total production cost (excludes cost of feedstock).
            royalty_charge_percentage_of_revenue: Percentage of revenue charged as royalties;
                defaults to 6.5% as reported in the UKy report.
            transport_cost_per_ton_product: Expression, Var or Param to use for transport costs
                per ton of product (note, this is not part of the TOC)
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
            watertap_blocks: list of unit model blocks corresponding to watertap models
            calculate_NPV: True/false flag for calculating net present value (NPV).
        """

        # define costing library
        if hasattr(self, "library") and self.library == "REE":  # costing already exists
            raise RuntimeError(
                f"Costing for the block {self} already exists. Please ensure that "
                f"the costing build method is not called twice on the same "
                f"model."
            )
        self.library = "REE"

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2020."
            )

        if (
            fixed_OM is False and calculate_NPV is True
        ):  # NPV method invoked with fixed inputs
            self.calculate_NPV(fixed_OM, variable_OM)

        else:  # continue on with building the plant costs

            self.BEC_list = []
            self.watertap_fixed_costs_list = []
            # TODO commented as no WaterTAP models currently use this, may change in the future
            # self.watertap_variable_costs_list = []
            self.custom_fixed_costs_list = []
            self.custom_variable_costs_list = []

            if total_purchase_cost is None:
                self.get_total_BEC(CE_index_year, watertap_blocks)
            else:
                self.total_BEC = Var(
                    initialize=total_purchase_cost,
                    units=CE_index_units,
                )
                self.total_BEC.fix()

            # define variables
            if Lang_factor is None:
                # initialize parameters from specified percentages
                self.piping_materials_and_labor_percentage = Param(
                    mutable=True,
                    initialize=piping_materials_and_labor_percentage,
                    doc="Percentage of BEC used to estimate piping, materials and labor installation costs",
                    units=pyunits.percent,
                )

                self.electrical_materials_and_labor_percentage = Param(
                    mutable=True,
                    initialize=electrical_materials_and_labor_percentage,
                    doc="Percentage of BEC used to estimate electrical, materials and labor installation costs",
                    units=pyunits.percent,
                )

                self.instrumentation_percentage = Param(
                    mutable=True,
                    initialize=instrumentation_percentage,
                    doc="Percentage of BEC used to estimate instrumentation installation costs",
                    units=pyunits.percent,
                )

                self.plant_services_percentage = Param(
                    mutable=True,
                    initialize=plants_services_percentage,
                    doc="Percentage of BEC used to estimate plant services installation costs",
                    units=pyunits.percent,
                )

                self.process_buildings_percentage = Param(
                    mutable=True,
                    initialize=process_buildings_percentage,
                    doc="Percentage of BEC used to estimate process buildings installation costs",
                    units=pyunits.percent,
                )

                self.auxiliary_buildings_percentage = Param(
                    mutable=True,
                    initialize=auxiliary_buildings_percentage,
                    doc="Percentage of BEC used to estimate auxiliary buildings installation costs",
                    units=pyunits.percent,
                )

                self.site_improvements_percentage = Param(
                    mutable=True,
                    initialize=site_improvements_percentage,
                    doc="Percentage of BEC used to estimate site improvements installation costs",
                    units=pyunits.percent,
                )

                self.equipment_installation_percentage = Param(
                    mutable=True,
                    initialize=equipment_installation_percentage,
                    doc="Percentage of BEC used to estimate equipment installation costs",
                    units=pyunits.percent,
                )

                self.field_expenses_percentage = Param(
                    mutable=True,
                    initialize=field_expenses_percentage,
                    doc="Percentage of BEC used to estimate field expenses installation costs",
                    units=pyunits.percent,
                )

                self.project_management_and_construction_percentage = Param(
                    mutable=True,
                    initialize=project_management_and_construction_percentage,
                    doc="Percentage of BEC used to estimate project management and construction installation costs",
                    units=pyunits.percent,
                )

                self.process_contingency_percentage = Param(
                    mutable=True,
                    initialize=process_contingency_percentage,
                    doc="Percentage of BEC used to estimate process contingency installation costs",
                    units=pyunits.percent,
                )

                # ancillary cost variables
                self.ancillary_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Ancillary cost",
                    units=CE_index_units,
                )

                self.piping_materials_and_labor_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Piping, materials and labor ancillary cost",
                    units=CE_index_units,
                )

                self.electrical_materials_and_labor_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Electrical, materials and labor ancillary cost",
                    units=CE_index_units,
                )

                self.instrumentation_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Ancillary cost",
                    units=CE_index_units,
                )

                self.plant_services_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Ancillary cost",
                    units=CE_index_units,
                )

                # buildings cost variables
                self.buildings_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Buildings cost",
                    units=CE_index_units,
                )

                self.process_buildings_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Process buildings cost",
                    units=CE_index_units,
                )

                self.auxiliary_buildings_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Auxiliary buildings cost",
                    units=CE_index_units,
                )

                self.site_improvements_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Site improvements buildings cost",
                    units=CE_index_units,
                )

                # engineering, procurement and construction management cost variables
                self.epcm_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="EPCM cost",
                    units=CE_index_units,
                )

                self.equipment_installation_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Equipment installation EPCM cost",
                    units=CE_index_units,
                )

                self.field_expenses_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Field expenses EPCM cost",
                    units=CE_index_units,
                )

                self.project_management_and_construction_costs = Var(
                    initialize=self.total_BEC,
                    bounds=(0, None),
                    doc="Project management and construction EPCM cost",
                    units=CE_index_units,
                )

                # contingency cost variables - generic to support more contingency cost types in the future
                self.contingency_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Contingency cost",
                    units=CE_index_units,
                )

                self.process_contingency_costs = Var(
                    initialize=value(self.total_BEC),
                    bounds=(0, None),
                    doc="Contingency cost",
                    units=CE_index_units,
                )
            else:
                self.Lang_factor = Param(
                    initialize=Lang_factor,
                    mutable=True,
                    doc="Lang factor",
                    units=pyunits.dimensionless,
                )

            # total cost variables
            self.total_installation_cost = Var(
                initialize=self.total_BEC,
                bounds=(0, None),
                doc="Total installation cost",
                units=CE_index_units,
            )

            self.total_plant_cost = Var(
                initialize=self.total_BEC,
                bounds=(0, None),
                doc="Total plant cost",
                units=CE_index_units,
            )

            # add other plant costs to catch non-equipment capital costs, e.g. reagent fills
            self.other_plant_costs = Var(
                initialize=0,
                bounds=(0, None),
                doc="Additional plant costs",
                units=CE_index_units,
            )
            self.other_plant_costs.fix(1e-12)

            if Lang_factor is None:
                # constraints for calculating Ancillary costs
                @self.Constraint()
                def piping_materials_and_labor_cost_eq(c):
                    return c.piping_materials_and_labor_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.piping_materials_and_labor_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def electrical_materials_and_labor_cost_eq(c):
                    return c.electrical_materials_and_labor_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.electrical_materials_and_labor_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def instrumentation_cost_eq(c):
                    return c.instrumentation_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.instrumentation_percentage, to_units=pyunits.dimensionless
                        )
                    )

                @self.Constraint()
                def plant_services_cost_eq(c):
                    return c.plant_services_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.plant_services_percentage, to_units=pyunits.dimensionless
                        )
                    )

                @self.Constraint()
                def ancillary_cost_eq(c):
                    return c.ancillary_costs == (
                        c.piping_materials_and_labor_costs
                        + c.electrical_materials_and_labor_costs
                        + c.instrumentation_costs
                        + c.plant_services_costs
                    )

                # constraints for calculating Buildings costs
                @self.Constraint()
                def process_buildings_cost_eq(c):
                    return c.process_buildings_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.process_buildings_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def auxiliary_buildings_cost_eq(c):
                    return c.auxiliary_buildings_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.auxiliary_buildings_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def site_improvements_cost_eq(c):
                    return c.site_improvements_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.site_improvements_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def buildings_cost_eq(c):
                    return c.buildings_costs == (
                        c.process_buildings_costs
                        + c.auxiliary_buildings_costs
                        + c.site_improvements_costs
                    )

                # constraints for calculating Engineering, Procurement and Construction Management costs
                @self.Constraint()
                def equipment_installation_cost_eq(c):
                    return c.equipment_installation_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.equipment_installation_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def field_expenses_cost_eq(c):
                    return c.field_expenses_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.field_expenses_percentage, to_units=pyunits.dimensionless
                        )
                    )

                @self.Constraint()
                def project_management_and_construction_cost_eq(c):
                    return c.project_management_and_construction_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.project_management_and_construction_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def epcm_cost_eq(c):
                    return c.epcm_costs == (
                        c.equipment_installation_costs
                        + c.field_expenses_costs
                        + c.project_management_and_construction_costs
                    )

                # constraints for calculating Contingency costs
                @self.Constraint()
                def process_contingency_cost_eq(c):
                    return c.contingency_costs == (
                        c.total_BEC
                        * pyunits.convert(
                            c.process_contingency_percentage,
                            to_units=pyunits.dimensionless,
                        )
                    )

                @self.Constraint()
                def contingency_cost_eq(c):
                    return c.contingency_costs == (c.process_contingency_costs)

                @self.Constraint()
                def total_installation_cost_eq(c):
                    return c.total_installation_cost == (
                        c.ancillary_costs
                        + c.buildings_costs
                        + c.epcm_costs
                        + c.contingency_costs
                    )

            else:

                @self.Constraint()
                def total_installation_cost_eq(c):
                    return c.total_installation_cost == c.total_BEC * (
                        c.Lang_factor - 1
                    )

            # constraint for calculating TPC
            @self.Constraint()
            def total_plant_cost_eq(c):
                return c.total_plant_cost == (
                    c.total_BEC + c.total_installation_cost + c.other_plant_costs
                )

            # define land cost
            if land_cost is not None:
                if type(land_cost) in [Expression, ScalarExpression]:
                    if pyunits.get_units(land_cost) == pyunits.dimensionless:
                        self.land_cost = Expression(
                            expr=land_cost.expr * CE_index_units
                        )
                    else:
                        self.land_cost = Expression(
                            expr=pyunits.convert(
                                land_cost.expr, to_units=CE_index_units
                            )
                        )
                else:
                    if pyunits.get_units(land_cost) == pyunits.dimensionless:
                        self.land_cost = Expression(expr=land_cost * CE_index_units)
                    else:
                        self.land_cost = Expression(
                            expr=pyunits.convert(land_cost, to_units=CE_index_units)
                        )
            else:
                self.land_cost = Expression(expr=0 * CE_index_units)

            # define feed input, if passed
            if feed_input is not None:
                if (
                    pyunits.get_units(feed_input) == pyunits.dimensionless
                ):  # require units
                    raise UnitsError(
                        "The argument feed_input was passed as a dimensionless "
                        "quantity with no units. Please ensure that the feed "
                        "rate is passed in units of mass / time."
                    )
                else:
                    feed_input_rate = pyunits.convert(
                        feed_input, to_units=pyunits.ton / pyunits.hr
                    )
            else:
                feed_input_rate = None

            # build operating & maintenance costs
            if chemicals is None:
                self.chemicals_list = []
            else:
                self.chemicals_list = chemicals

            # define additional chemicals cost
            if additional_chemicals_cost is not None:
                if type(additional_chemicals_cost) in [Expression, ScalarExpression]:
                    if (
                        pyunits.get_units(additional_chemicals_cost)
                        == pyunits.dimensionless
                    ):
                        self.additional_chemicals_cost = Expression(
                            expr=additional_chemicals_cost.expr * CE_index_units
                        )
                    else:
                        self.additional_chemicals_cost = Expression(
                            expr=pyunits.convert(
                                additional_chemicals_cost.expr, to_units=CE_index_units
                            )
                        )
                else:
                    if (
                        pyunits.get_units(additional_chemicals_cost)
                        == pyunits.dimensionless
                    ):
                        self.additional_chemicals_cost = Expression(
                            expr=additional_chemicals_cost * CE_index_units
                        )
                    else:
                        self.additional_chemicals_cost = Expression(
                            expr=pyunits.convert(
                                additional_chemicals_cost, to_units=CE_index_units
                            )
                        )
            else:
                self.additional_chemicals_cost = Expression(expr=0 * CE_index_units)

            if waste is None:
                self.waste_list = []
            else:
                self.waste_list = waste

            # define waste cost
            if additional_waste_cost is not None:
                if type(additional_waste_cost) in [Expression, ScalarExpression]:
                    if (
                        pyunits.get_units(additional_waste_cost)
                        == pyunits.dimensionless
                    ):
                        self.additional_waste_cost = Expression(
                            expr=additional_waste_cost.expr * CE_index_units
                        )
                    else:
                        self.additional_waste_cost = Expression(
                            expr=pyunits.convert(
                                additional_waste_cost.expr, to_units=CE_index_units
                            )
                        )
                else:
                    if (
                        pyunits.get_units(additional_waste_cost)
                        == pyunits.dimensionless
                    ):
                        self.additional_waste_cost = Expression(
                            expr=additional_waste_cost * CE_index_units
                        )
                    else:
                        self.additional_waste_cost = Expression(
                            expr=pyunits.convert(
                                additional_waste_cost, to_units=CE_index_units
                            )
                        )
            else:
                self.additional_waste_cost = Expression(expr=0 * CE_index_units)

            if fixed_OM:
                self.get_fixed_OM_costs(
                    labor_types=labor_types,
                    labor_rate=labor_rate,
                    labor_burden=labor_burden,
                    operators_per_shift=operators_per_shift,
                    hours_per_shift=hours_per_shift,
                    shifts_per_day=shifts_per_day,
                    operating_days_per_year=operating_days_per_year,
                    pure_product_output_rates=pure_product_output_rates,
                    mixed_product_output_rates=mixed_product_output_rates,
                    mixed_product_sale_price_realization_factor=mixed_product_sale_price_realization_factor,
                    sale_prices=sale_prices,
                    CE_index_year=CE_index_year,
                )

            if variable_OM:
                self.get_variable_OM_costs(
                    efficiency=efficiency,
                    resources=resources,
                    rates=rates,
                    prices=prices,
                    feed_input_rate=feed_input_rate,
                    CE_index_year=CE_index_year,
                )

            # build system costs (owner's, total overnight costs, annualized costs,
            # and cost of recovery)

            self.total_overnight_capital = Expression(expr=self.total_plant_cost)

            self.tasc_toc_factor = Param(
                initialize=1.144,
                mutable=True,
                doc="TASC/TOC factor calculated from UKy report using 3 year "
                "expenditure period with 10/60/30 % expenditure at 3.6% "
                "escalation at 2.94% debt interest rate with 7.84% return on "
                "equity, 26% combined federal/state tax, and 50/50 % debt and "
                "equity financed.",
            )

            self.total_as_spent_cost = Expression(
                expr=self.total_overnight_capital * self.tasc_toc_factor
            )

            self.fixed_charge_factor = Param(
                initialize=0.1002,
                mutable=True,
                doc="Fixed charge rate calculated from UKy report using a 26% "
                "effective tax rate, a tax depreciation fraction of 2.231 over "
                "21 years of depreciation, a nominal capital recovery factor of "
                "0.0856, an after-tax weighted average cost of capital of 5.77%, "
                "= Present value of tax depreciation expense of 0.237",
            )
            self.annualized_cost = Expression(
                expr=self.fixed_charge_factor * self.total_as_spent_cost
            )

            if fixed_OM and variable_OM:
                if consider_taxes:
                    self.income_tax_percentage = Param(
                        initialize=income_tax_percentage,
                        mutable=True,
                        doc="Combined federal and state income tax percentage"
                        "usually between 26 - 40%",
                        units=pyunits.percent,
                    )
                    self.mineral_depletion_percentage = Param(
                        initialize=mineral_depletion_percentage,
                        mutable=True,
                        doc="tax deduction percentage for mineral depletion."
                        "default value of 14% is used, as reported in the UKy report",
                        units=pyunits.percent,
                    )
                    self.production_incentive_percentage = Param(
                        initialize=production_incentive_percentage,
                        mutable=True,
                        doc="tax deduction percentage for producing critical minerals"
                        "default value of 10% of total production cost",
                        units=pyunits.percent,
                    )
                    self.royalty_charge_percentage_of_revenue = Param(
                        initialize=royalty_charge_percentage_of_revenue,
                        mutable=True,
                        doc="Percentage of revenue charged as royalties",
                        units=pyunits.percent,
                    )
                    self.min_net_tax_owed = Param(
                        initialize=0,
                        doc="Minimum net tax owed in millions USD",
                        units=CE_index_units / pyunits.year,
                    )
                    self.eps = Param(
                        initialize=1e-4,
                        units=CE_index_units / pyunits.year,
                    )
                    self.net_tax_owed = Var(
                        initialize=0.40 * self.total_sales_revenue,
                        doc="Net tax owed in millions USD",
                        units=CE_index_units / pyunits.year,
                    )
                    self.income_tax = Var(
                        initialize=0.26 * self.total_sales_revenue,
                        doc="Income tax in millions USD",
                        units=CE_index_units / pyunits.year,
                    )
                    self.additional_tax_credit = Var(
                        initialize=0,
                        doc="Additional tax credit",
                        units=CE_index_units / pyunits.year,
                    )
                    self.additional_tax_credit.fix(1e-12)

                    self.additional_tax_owed = Var(
                        initialize=0,
                        doc="Additional tax owed",
                        units=CE_index_units / pyunits.year,
                    )
                    self.additional_tax_owed.fix(1e-12)

                    self.royalty_charge = Expression(
                        expr=pyunits.convert(
                            self.royalty_charge_percentage_of_revenue,
                            to_units=pyunits.dimensionless,
                        )
                        * self.total_sales_revenue
                    )
                    self.mineral_depletion_charge = Expression(
                        expr=pyunits.convert(
                            self.mineral_depletion_percentage,
                            to_units=pyunits.dimensionless,
                        )
                        * (self.total_sales_revenue - self.royalty_charge)
                    )
                    self.production_incentive_charge = Expression(
                        expr=pyunits.convert(
                            self.production_incentive_percentage,
                            to_units=pyunits.dimensionless,
                        )
                        * (
                            self.total_variable_OM_cost[0]
                            + self.total_fixed_OM_cost
                            + self.annualized_cost / pyunits.year
                        )
                    )

                    @self.Constraint()
                    def income_tax_eq(c):
                        return c.income_tax == pyunits.convert(
                            self.income_tax_percentage, to_units=pyunits.dimensionless
                        ) * (
                            self.total_sales_revenue
                            - (
                                self.total_variable_OM_cost[0]
                                + self.total_fixed_OM_cost
                                + self.annualized_cost / pyunits.year
                            )
                        )

                    @self.Constraint()
                    def net_tax_owed_eq(c):
                        return c.net_tax_owed == smooth_max(
                            self.min_net_tax_owed,
                            (
                                (
                                    self.income_tax
                                    + self.royalty_charge
                                    + self.additional_tax_owed
                                )
                                - (
                                    self.mineral_depletion_charge
                                    + self.production_incentive_charge
                                    + self.additional_tax_credit
                                )
                            ),
                            eps=self.eps,
                        )

                # build cost of recovery (COR)
                if recovery_rate_per_year is not None:
                    self.additional_cost_of_recovery = Var(
                        initialize=0,
                        doc="Additional cost to be added to the COR calculations"
                        + " in millions",
                        units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
                    )
                    self.additional_cost_of_recovery.fix()

                    if (
                        pyunits.get_units(recovery_rate_per_year)
                        == pyunits.dimensionless
                    ):
                        raise UnitsError(
                            "The argument recovery_rate_per_year was passed as a dimensionless "
                            "quantity with no units. Please ensure that the feed "
                            "rate is passed in units of mass / time."
                        )

                    if not hasattr(self, "recovery_rate_per_year"):
                        self.recovery_rate_per_year = Expression()

                    rec_rate_units = pyunits.get_units(recovery_rate_per_year)

                    # check that units are compatible
                    try:
                        conversion = (
                            value(
                                pyunits.convert(
                                    rec_rate_units,
                                    to_units=pyunits.kg / pyunits.year,
                                )
                            )
                            * pyunits.kg
                            / pyunits.year
                            / rec_rate_units
                        )
                    except InconsistentUnitsError:
                        raise UnitsError(
                            f"The argument recovery_rate_per_year was passed with units of "
                            f"{rec_rate_units} which cannot be converted to units of mass per year. "
                            f"Please ensure that recovery_rate_per_year is passed with rate units "
                            f"of mass per year (mass/a)."
                        )

                    self.recovery_rate_per_year.expr = (
                        recovery_rate_per_year * conversion
                    )

                    self.cost_of_recovery = Expression(
                        expr=(
                            pyunits.convert(
                                (
                                    self.annualized_cost / pyunits.year
                                    + self.total_fixed_OM_cost
                                    + self.total_variable_OM_cost[0]
                                    + (
                                        self.net_tax_owed
                                        if consider_taxes
                                        else 0 * CE_index_units / pyunits.year
                                    )
                                )
                                / (self.recovery_rate_per_year),
                                to_units=getattr(pyunits, "USD_" + CE_index_year)
                                / pyunits.kg,
                            )
                            + self.additional_cost_of_recovery
                        )
                    )

                    if transport_cost_per_ton_product is not None:
                        if (
                            isinstance(
                                transport_cost_per_ton_product,
                                (Expression, ScalarExpression, Param, Var),
                            )
                            and pyunits.get_units(transport_cost_per_ton_product)
                            == pyunits.dimensionless
                        ) or isinstance(transport_cost_per_ton_product, (int, float)):
                            self.transport_cost = (
                                transport_cost_per_ton_product
                                * 1e-6
                                * CE_index_units
                                / pyunits.ton
                                * pyunits.convert(
                                    self.recovery_rate_per_year,
                                    to_units=pyunits.ton / pyunits.year,
                                )
                            )
                        else:
                            self.transport_cost = pyunits.convert(
                                transport_cost_per_ton_product
                                * self.recovery_rate_per_year,
                                to_units=CE_index_units / pyunits.year,
                            )

                else:  # except the case where transport_cost_per_ton_product is passed but recovery_rate_per_year is not passed
                    if transport_cost_per_ton_product is not None:
                        raise AttributeError(
                            "If transport_cost_per_ton_product is not None, "
                            "recovery_rate_per_year cannot be None."
                        )

            if calculate_NPV:
                self.calculate_NPV(fixed_OM, variable_OM, consider_taxes)

    @staticmethod
    def initialize_build(*args, **kwargs):
        """
        Here we can add initialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now,  no additional process level costs to initialize

    def report(self, export=False):
        var_dict = {}
        var_dict["Plant Cost Units"] = str(pyunits.get_units(self.total_plant_cost))

        if hasattr(self, "total_plant_cost"):
            var_dict["Total Plant Cost"] = value(self.total_plant_cost)

        if hasattr(self, "total_BEC"):
            var_dict["Total Bare Erected Cost"] = value(self.total_BEC)

        if hasattr(self, "total_installation_cost"):
            var_dict["Total Installation Cost"] = value(self.total_installation_cost)

        if hasattr(self, "other_plant_costs"):
            var_dict["Total Other Plant Costs"] = value(self.other_plant_costs)

        if hasattr(self, "ancillary_costs"):
            var_dict["Summation of Ancillary Installation Costs"] = value(
                self.ancillary_costs
            )

        if hasattr(self, "piping_materials_and_labor_costs"):
            var_dict[
                "Total Ancillary Piping, Materials and Labor Installation Cost"
            ] = value(self.piping_materials_and_labor_costs)

        if hasattr(self, "electrical_materials_and_labor_costs"):
            var_dict[
                "Total Ancillary Electrical, Materials and Labor Installation Cost"
            ] = value(self.electrical_materials_and_labor_costs)

        if hasattr(self, "instrumentation_costs"):
            var_dict["Total Ancillary Instrumentation Installation Cost"] = value(
                self.instrumentation_costs
            )

        if hasattr(self, "plant_services_costs"):
            var_dict["Total Ancillary Plant Services Installation Cost"] = value(
                self.plant_services_costs
            )

        if hasattr(self, "buildings_costs"):
            var_dict["Summation of Buildings Installation Costs"] = value(
                self.buildings_costs
            )

        if hasattr(self, "process_buildings_costs"):
            var_dict["Total Process Buildings Installation Cost"] = value(
                self.process_buildings_costs
            )

        if hasattr(self, "auxiliary_buildings_costs"):
            var_dict["Total Auxiliary Buildings Installation Cost"] = value(
                self.auxiliary_buildings_costs
            )

        if hasattr(self, "site_improvements_costs"):
            var_dict["Total Site Improvements Buildings Installation Cost"] = value(
                self.site_improvements_costs
            )

        if hasattr(self, "epcm_costs"):
            var_dict["Summation of EPCM Installation Costs"] = value(self.epcm_costs)

        if hasattr(self, "equipment_installation_costs"):
            var_dict["Total Equipment Installation EPCM Installation Cost"] = value(
                self.equipment_installation_costs
            )

        if hasattr(self, "field_expenses_costs"):
            var_dict["Total Field Expenses EPCM Cost"] = value(
                self.field_expenses_costs
            )

        if hasattr(self, "project_management_and_construction_costs"):
            var_dict[
                "Total Project Management and Construction EPCM Installation Cost"
            ] = value(self.project_management_and_construction_costs)

        if hasattr(self, "process_contingency_costs"):
            var_dict["Total Process Contingency Installation Cost"] = value(
                self.process_contingency_costs
            )

        if hasattr(self, "contingency_costs"):
            var_dict["Summation of Contingency Installation Costs"] = value(
                self.contingency_costs
            )

        if hasattr(self, "total_fixed_OM_cost"):
            var_dict["Total Fixed Operating & Maintenance Cost"] = value(
                self.total_fixed_OM_cost
            )

        if hasattr(self, "annual_operating_labor_cost"):
            var_dict["Total Annual Operating Labor Cost"] = value(
                self.annual_operating_labor_cost
            )

            var_dict["Total Annual Technical Labor Cost"] = value(
                self.annual_technical_labor_cost
            )

            var_dict["Summation of Annual Labor Costs"] = value(self.annual_labor_cost)

        if hasattr(self, "maintenance_and_material_cost"):
            var_dict["Total Maintenance and Material Cost"] = value(
                self.maintenance_and_material_cost
            )

        if hasattr(self, "quality_assurance_and_control_cost"):
            var_dict["Total Quality Assurance and Control Cost"] = value(
                self.quality_assurance_and_control_cost
            )

        general_sales_and_admin = 0

        if hasattr(self, "sales_patenting_and_research_cost"):
            var_dict["Total Sales, Patenting and Research Cost"] = value(
                self.sales_patenting_and_research_cost
            )
            general_sales_and_admin += value(self.sales_patenting_and_research_cost)

            var_dict["Summation of Sales, Admin and Insurance Cost"] = value(
                general_sales_and_admin
            )

        if hasattr(self, "admin_and_support_labor_cost"):
            var_dict["Total Admin Support and Labor Cost"] = value(
                self.admin_and_support_labor_cost
            )
            general_sales_and_admin += value(self.admin_and_support_labor_cost)

            var_dict["Summation of Sales, Admin and Insurance Cost"] = value(
                general_sales_and_admin
            )

        if hasattr(self, "property_taxes_and_insurance_cost"):
            var_dict["Total Property Taxes and Insurance Cost"] = value(
                self.property_taxes_and_insurance_cost
            )
            general_sales_and_admin += value(self.property_taxes_and_insurance_cost)

            var_dict["Summation of Sales, Admin and Insurance Cost"] = value(
                general_sales_and_admin
            )

        if hasattr(self, "other_fixed_costs"):
            var_dict["Total Other Fixed Costs"] = value(self.other_fixed_costs)

        if hasattr(self, "variable_operating_costs"):
            if (0, "power") in self.variable_operating_costs.id_index_map().values():
                var_dict["Total Variable Power Cost"] = value(
                    self.variable_operating_costs[0, "power"]
                )

            if hasattr(self, "additional_waste_cost"):
                var_dict["Total Variable Waste Cost"] = value(
                    sum(
                        self.variable_operating_costs[0, waste]
                        for waste in self.waste_list
                    )
                    + self.additional_waste_cost
                )

            if hasattr(self, "additional_chemicals_cost"):
                var_dict["Total Variable Chemicals Cost"] = value(
                    sum(
                        self.variable_operating_costs[0, chemical]
                        for chemical in self.chemicals_list
                    )
                    + self.additional_chemicals_cost
                )

            var_dict["General Plant Overhead Cost"] = value(self.plant_overhead_cost[0])

            var_dict[
                "Total Plant Overhead Cost, Including Maintenance & Quality Assurance"
            ] = value(
                self.plant_overhead_cost[0]
                + self.maintenance_and_material_cost
                + self.quality_assurance_and_control_cost
            )

        if hasattr(self, "total_variable_OM_cost"):
            var_dict["Total Variable Operating & Maintenance Cost"] = value(
                self.total_variable_OM_cost[0]
            )

        if hasattr(self, "land_cost"):
            var_dict["Total Land Cost"] = value(self.land_cost)

        if hasattr(self, "transport_cost"):
            var_dict["Total Transport Cost"] = value(self.transport_cost)

        if hasattr(self, "total_sales_revenue"):
            var_dict["Total Sales Revenue Cost"] = value(self.total_sales_revenue)

        if hasattr(self, "npv"):
            var_dict["Net Present Value"] = value(self.npv)

        if hasattr(self, "royalty_charge"):
            var_dict["Royalty Charge"] = value(self.royalty_charge)

        if hasattr(self, "mineral_depletion_charge"):
            var_dict["Mineral Depletion Charge"] = value(self.mineral_depletion_charge)

        if hasattr(self, "production_incentive_charge"):
            var_dict["Production Incentive Charge"] = value(
                self.production_incentive_charge
            )

        if hasattr(self, "income_tax"):
            var_dict["Income Tax"] = value(self.income_tax)

        if hasattr(self, "net_tax_owed"):
            var_dict["Net Tax Owed"] = value(self.net_tax_owed)

        report_dir = {}
        report_dir["Value"] = {}
        report_dir["pos"] = {}

        count = 1
        for k, v in var_dict.items():
            report_dir["Value"][k] = value(v)
            report_dir["pos"][k] = count
            count += 1

        df = DataFrame.from_dict(report_dir, orient="columns")
        del df["pos"]
        if export:
            df.to_csv(f"{self.local_name}_report.csv")

        print("\n" + "=" * 84)
        print(f"{self.local_name}")
        print("-" * 84)
        stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
        print("\n" + "=" * 84 + "\n")

    # -----------------------------------------------------------------------------
    # REE Recovery Costing Library
    # -----------------------------------------------------------------------------
    def get_REE_costing(
        blk,
        cost_accounts,
        scaled_param,
        source,
        Lang_factor=None,
        n_equip=1,
        scale_down_parallel_equip=False,
        CE_index_year="2021",
        additional_costing_params=None,
        use_additional_costing_params=False,
    ):
        """
        The scaled cost is computed using reference values for different
        sources as listed below:
        1. University of Kentucky Fire Clay Seam (Hazard No. 4) Rejects

        Args:
            blk: A unit-level costing block where costing variables and
                constraints can be added to
            cost_accounts: A list of accounts to be included in the total cost
            scaled_param: the process parameter for the system(s) being costed;
                this is the total flow for all parallel trains of the system(s)
            source: integer representing the above categories
            Lang_factor: optional single Lang factor value used to calculate
                commercial-scale installation costs. If None, accounts must
                include necessary component factors.
            n_equip: Integer number of parallel equipment trains for unit
                operations; for example, enter '5' if a feed will be split
                among 5 identical units and then re-mixed
            scale_down_parallel_equip: Boolean flag whether to scale down
                parallel equipment trains, e.g. two trains scaled down will
                each be half the size/capacity of a single train, and two
                trains not scaled down will each be the same size as a single
                train (twice the capacity). If duplicating a fixed operation
                into multiple parallel trains, use the default value ('False').
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
            additional_costing_params: user-defined dictionary to append to
                existing cost accounts dictionary
            use_additional_costing_params: Boolean flag to use additional
                costing parameters when account names conflict with existing
                accounts data
        """
        # check to see if a costing block already exists
        if (
            blk.parent_block().name
            in blk.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                f"{blk.name} already has an attribute costing. "
                f"Check that you are not calling get_costing"
                f" twice on the same model"
            )

        # define costing library
        blk.library = "REE"

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2020."
            )

        # pull data for each account into dictionaries
        process_params = {}
        reference_units = {}
        account_names = {}
        exponents = {}
        reference_costs = {}
        reference_cost_units = {}
        reference_costs_init = {}
        reference_params = {}

        # load ree costing dictionary
        REE_costing_params = load_REE_costing_dictionary()

        # for compatibility with potential custom accounts, the loop handles
        # new sources, and new accounts for existing sources
        # Users should not be adding new entries for existing accounts

        costing_params = REE_costing_params  # initialize with baseline accounts
        if additional_costing_params is not None and additional_costing_params != {}:
            for new_costing_params in [
                additional_costing_params
            ]:  # merge new dictionaries sequentially
                # adding any provided custom params to the base dictionary
                # need to "freeze" dict so it is hashable for merging keys
                frozen_dict = {**costing_params}
                for sourcekey, sourceval in new_costing_params.items():
                    if (
                        sourcekey in frozen_dict.keys()
                    ):  # if sourcekey already exists, append any new accounts
                        for accountkey, accountval in new_costing_params[
                            sourcekey
                        ].items():
                            if (
                                accountkey in frozen_dict[sourcekey].keys()
                            ) and not use_additional_costing_params:
                                if accountkey not in cost_accounts:
                                    pass  # not the current account, don't fail here
                                else:  # this is not allowed
                                    raise ValueError(
                                        f"Data already exists for Account {accountkey} "
                                        f"using source {sourcekey}. "
                                        f"Please confirm that the custom "
                                        f"account dictionary is correct, or "
                                        f"add the new parameters as a new "
                                        f"account. To use the custom account "
                                        f"dictionary for all conflicts, please "
                                        f"pass the argument use_additional_costing_params "
                                        f"as True."
                                    )
                            else:  # conflict is the account passed, and overwrite it
                                frozen_dict[sourcekey][accountkey] = accountval
                    else:
                        frozen_dict[sourcekey] = sourceval
                costing_params = {k: frozen_dict[k] for k in sorted(frozen_dict)}

        for account in cost_accounts:
            try:  # look for data in json file info
                process_params[account] = costing_params[str(source)][account][
                    "Process Parameter"
                ]
                reference_units[account] = costing_params[str(source)][
                    cost_accounts[0]
                ]["Units"]
                account_names[account] = costing_params[str(source)][account][
                    "Account Name"
                ]
                exponents[account] = float(
                    costing_params[str(source)][account]["Exponent"]
                )
                reference_costs[account] = costing_params[str(source)][account]["BEC"]
                reference_cost_units[account] = costing_params[str(source)][account][
                    "BEC_units"
                ]
                reference_costs_init[account] = (
                    costing_params[str(source)][account]["BEC"] * 1e-3
                )

                if isinstance(process_params[account], list):
                    for i, processparam in enumerate(process_params[account]):
                        reference_params[account, processparam] = costing_params[
                            str(source)
                        ][account]["RP Value"][i]

                elif isinstance(process_params[account], str):
                    reference_params[account] = costing_params[str(source)][account][
                        "RP Value"
                    ]
            except KeyError:
                raise KeyError(
                    f"Account {account} could not be found in the dictionary for "
                    f"source {source}"
                )

        # check that all accounts use the same process parameter
        param_check = None
        for account in cost_accounts:
            param = process_params[account]
            if param_check is None:
                param_check = param
            elif param != param_check:
                raise ValueError(
                    f"{blk.name} cost accounts selected do not use "
                    f"the same process parameter"
                )

        # check that the user passed the correct units type and try to convert

        for account in cost_accounts:
            ref_units = reference_units[account]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                if "**" in ref_units[0]:
                    ref_units[0] = ref_units[0].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0][0]) ** int(
                            ref_units[0][1]
                        ) / getattr(pyunits, ref_units[1])
                    except AttributeError:
                        expected_units = str(
                            ref_units[0][0]
                            + "**"
                            + ref_units[0][1]
                            + "/"
                            + ref_units[1]
                        )
                        raise AttributeError(
                            f"Account {cost_accounts[0]} uses references units of "
                            f"{expected_units}. "
                            f"Cannot parse reference units as Pyomo unit containers. "
                            f"Check that source uses correct syntax for Pyomo "
                            f"unit containers, for example gpm should be "
                            f"gal/min, tpd should be ton/d and MMBtu should be "
                            f"MBtu (using Pyomo prefix)."
                        )
                elif "**" in ref_units[1]:
                    ref_units[1] = ref_units[1].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1][0]
                        ) ** int(ref_units[1][1])
                    except AttributeError:
                        expected_units = str(
                            ref_units[0]
                            + "/"
                            + ref_units[1][0]
                            + "**"
                            + ref_units[1][1]
                        )
                        raise AttributeError(
                            f"Account {cost_accounts[0]} uses references units of "
                            f"{expected_units}. "
                            f"Cannot parse reference units as Pyomo unit containers. "
                            f"Check that source uses correct syntax for Pyomo "
                            f"unit containers, for example gpm should be "
                            f"gal/min, tpd should be ton/d and MMBtu should be "
                            f"MBtu (using Pyomo prefix)."
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1]
                        )
                    except AttributeError:
                        expected_units = str(ref_units[0] + "/" + ref_units[1])
                        raise AttributeError(
                            f"Account {cost_accounts[0]} uses references units of "
                            f"{expected_units}. "
                            f"Cannot parse reference units as Pyomo unit containers. "
                            f"Check that source uses correct syntax for Pyomo "
                            f"unit containers, for example gpm should be "
                            f"gal/min, tpd should be ton/d and MMBtu should be "
                            f"MBtu (using Pyomo prefix)."
                        )

            else:
                if "**" in ref_units:
                    ref_units = ref_units.split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) ** int(ref_units[1])
                    except AttributeError:
                        expected_units = str(ref_units[0] + "/" + ref_units[1])
                        raise AttributeError(
                            f"Account {cost_accounts[0]} uses references units of "
                            f"{expected_units}. "
                            f"Cannot parse reference units as Pyomo unit containers. "
                            f"Check that source uses correct syntax for Pyomo "
                            f"unit containers, for example gpm should be "
                            f"gal/min, tpd should be ton/d and MMBtu should be "
                            f"MBtu (using Pyomo prefix)."
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units)
                    except AttributeError:
                        expected_units = str(ref_units[0] + "/" + ref_units[1])
                        raise AttributeError(
                            f"Account {cost_accounts[0]} uses references units of "
                            f"{expected_units}. "
                            f"Cannot parse reference units as Pyomo unit containers. "
                            f"Check that source uses correct syntax for Pyomo "
                            f"unit containers, for example gpm should be "
                            f"gal/min, tpd should be ton/d and MMBtu should be "
                            f"MBtu (using Pyomo prefix)."
                        )

            if isinstance(scaled_param, list):
                for sp in scaled_param:
                    if sp.get_units() is None:
                        raise ValueError(
                            f"Account {cost_accounts[0]} uses units of {ref_units}. "
                            f"Units of {sp.get_units()} were passed. "
                            f"Scaled_param must have units."
                        )
                    else:
                        try:
                            pyunits.convert(sp, ref_units)
                        except InconsistentUnitsError:
                            raise UnitsError(
                                f"Account {cost_accounts[0]} uses units of {ref_units}. "
                                f"Units of {sp.get_units()} were passed. "
                                f"Cannot convert unit containers."
                            )
            else:
                try:
                    if pyunits.get_units(scaled_param) is None:
                        raise UnitsError(
                            f"Account {cost_accounts[0]} uses units of {ref_units}. "
                            f"Units of {pyunits.get_units(scaled_param)} were passed. "
                            f"Scaled_param must have units."
                        )
                    else:
                        try:
                            pyunits.convert(scaled_param, ref_units)
                        except InconsistentUnitsError:
                            raise UnitsError(
                                f"Account {cost_accounts[0]} uses units of {ref_units}. "
                                f"Units of {pyunits.get_units(scaled_param)} were passed. "
                                f"Cannot convert unit containers."
                            )
                except InconsistentUnitsError:
                    raise UnitsError(
                        f"The expression {scaled_param.name} has inconsistent units."
                    )

        # Used by other functions for reporting results
        blk.account_names = account_names

        # define parameters
        blk.exp = Param(
            cost_accounts,
            mutable=True,
            initialize=exponents,
            doc="Exponential parameter for account",
        )

        blk.ref_cost = Param(
            cost_accounts,
            mutable=True,
            initialize=reference_costs,
            doc="Reference cost for account",
            # units not defined here, since every account could have different
            # currency units
        )

        if isinstance(process_params[cost_accounts[0]], list):
            if len(process_params[cost_accounts[0]]) > 1:
                blk.ref_param = Param(
                    cost_accounts,
                    process_params[cost_accounts[0]],
                    mutable=True,
                    initialize=reference_params,
                    doc="Reference parameter for account",
                )
        elif isinstance(process_params[cost_accounts[0]], str):
            blk.ref_param = Param(
                cost_accounts,
                mutable=True,
                initialize=reference_params,
                doc="Reference parameter for account",
            )

        # define variables
        blk.bare_erected_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, None),
            doc="Scaled bare erected cost",
            units=CE_index_units,
        )

        # constraint for scaling BEC
        @blk.Constraint(cost_accounts)
        def bare_erected_cost_eq(c, i):
            ref_units = reference_units[i]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                ref_units = getattr(pyunits, ref_units[0]) / getattr(
                    pyunits, ref_units[1]
                )
            else:
                ref_units = getattr(pyunits, ref_units)

            ref_cost_units = reference_cost_units[i]
            ref_cost_units = ref_cost_units.split("$")
            if ref_cost_units[0] == "":  # no prefix
                ref_cost_units = getattr(pyunits, "USD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "K":  # thousands of USD
                ref_cost_units = getattr(pyunits, "kUSD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "M":  # millions of USD
                ref_cost_units = getattr(pyunits, "MUSD_" + ref_cost_units[1])

            # determine reference parameter scaler based on train scaling
            if scale_down_parallel_equip:
                scaler = n_equip
            else:
                scaler = 1

            if isinstance(process_params[i], list):
                if len(process_params[i]) > 1:
                    return c.bare_erected_cost[i] == (
                        n_equip
                        * pyunits.convert(
                            c.ref_cost[i] * ref_cost_units, CE_index_units
                        )
                        * sum(
                            (
                                pyunits.convert(scaled_param[j], ref_units)
                                / (scaler * c.ref_param[i, p] * ref_units)
                            )
                            ** c.exp[i]
                            for j, p in enumerate(process_params[i])
                        )
                    )
            elif isinstance(process_params[i], str):
                return c.bare_erected_cost[i] == (
                    n_equip
                    * pyunits.convert(c.ref_cost[i] * ref_cost_units, CE_index_units)
                    * (
                        pyunits.convert(scaled_param, ref_units)
                        / (scaler * c.ref_param[i] * ref_units)
                    )
                    ** c.exp[i]
                )

        # add variable and constraint scaling
        for i in cost_accounts:
            iscale.set_scaling_factor(blk.bare_erected_cost[i], 1)
            iscale.constraint_scaling_transform(
                blk.bare_erected_cost_eq[i], 1e-3, overwrite=False
            )

    # -----------------------------------------------------------------------------

    # -----------------------------------------------------------------------------
    # Operation & Maintenance Costing Library
    # -----------------------------------------------------------------------------

    # pylint: disable-next=dangerous-default-value
    def get_fixed_OM_costs(
        b,
        labor_types=None,
        labor_rate=None,
        labor_burden=25,
        operators_per_shift=None,
        hours_per_shift=8,
        shifts_per_day=3,
        operating_days_per_year=336,
        pure_product_output_rates=None,
        mixed_product_output_rates=None,
        mixed_product_sale_price_realization_factor=0.65,
        sale_prices=None,
        CE_index_year="2021",
    ):
        """
        Args:
            b: costing block to add fixed cost variables and constraints to
            labor_types: list of types of operators present in plant; assumed to
                correspond with labor rate and operator per shift lists
            labor_rate: hourly rate of plant operators in project dollar year;
                defined as list corresponding to different operator types
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses; assumed constant for all operator types. The
                value should be a percentage, for example 10 for 10%.
            operators_per_shift: number of operators per shift; defined as list
                of operators per shift for each operator type
            hours_per_shift: number of hours per shift
            shifts_per_day: number of shifts per day
            operating_days_per_year: number of operating days per year
            pure_product_output_rates: dict of production rates of each REE pure product
            mixed_product_output_rates: dict of production rates of each REE in the mixed product
            mixed_product_sale_price_realization_factor: multiplicative factor for selling impure products
            sale_prices: dict of sale prices to be added to the premade dictionary
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars

        Returns:
            None
        """

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2020."
            )

        # check that required product arguments were passed
        if not isinstance(pure_product_output_rates, dict):
            raise TypeError("product_output_rates argument must be a dict")
        if not isinstance(mixed_product_output_rates, dict):
            raise TypeError("product_output_rates argument must be a dict")

        # dictionary of default sale prices
        # the currency units are of USD
        # Purity, purchase quantity, purchasing time, and location, all affect the cost.
        default_sale_prices = load_default_sale_prices()

        if sale_prices is None:
            sale_prices = {}

        # add entries from sale_prices to default_sale_prices
        if not isinstance(sale_prices, dict):
            raise TypeError("Dictionary of custom sale_prices must be a dict object.")
        else:
            for key in sale_prices.keys():
                default_sale_prices[key] = sale_prices[key]

        # raise error if the user included a product not in default_sale_prices
        if not set(pure_product_output_rates).issubset(default_sale_prices.keys()):
            raise AttributeError(
                f"A pure product was included that does not contain a "
                f"sale price. Sale prices exist for the following products: "
                f"{list(default_sale_prices.keys())}"
            )
        elif not set(mixed_product_output_rates).issubset(default_sale_prices.keys()):
            raise AttributeError(
                f"A mixed product was included that does not contain a "
                f"sale price. Sale prices exist for the following products: "
                f"{list(default_sale_prices.keys())}"
            )

        # set default values
        if labor_types is None:
            labor_types = [
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
            ]

        if labor_rate is None:
            labor_rate = [27.90, 23.26, 30.29, 24.06, 23.43, 46.82]

        if operators_per_shift is None:
            operators_per_shift = [2, 5, 2, 3, 1, 2]

        # make params
        b.labor_rate = Param(
            labor_types,
            initialize=dict(zip(labor_types, labor_rate)),
            mutable=True,
            units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.hr,
        )
        b.labor_burden = Param(
            initialize=labor_burden, mutable=True, units=pyunits.dimensionless
        )
        b.operators_per_shift = Param(
            labor_types,
            initialize=dict(zip(labor_types, operators_per_shift)),
            mutable=True,
            units=pyunits.dimensionless,
        )
        b.hours_per_shift = Param(
            initialize=hours_per_shift, mutable=True, units=pyunits.hr
        )
        b.shifts_per_day = Param(
            initialize=shifts_per_day, mutable=True, units=pyunits.day**-1
        )
        b.operating_days_per_year = Param(
            initialize=operating_days_per_year,
            mutable=True,
            units=pyunits.day / pyunits.year,
        )
        b.mixed_product_sale_price_realization_factor = Param(
            initialize=mixed_product_sale_price_realization_factor,
            mutable=True,
            units=pyunits.dimensionless,
        )

        # make vars
        b.annual_operating_labor_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Annual operating labor cost",
            units=CE_index_units / pyunits.year,
        )
        b.annual_technical_labor_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Annual technical labor cost",
            units=CE_index_units / pyunits.year,
        )
        b.annual_labor_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Annual labor cost",
            units=CE_index_units / pyunits.year,
        )
        b.maintenance_and_material_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Maintenance and material cost",
            units=CE_index_units / pyunits.year,
        )
        b.quality_assurance_and_control_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Quality assurance and control cost",
            units=CE_index_units / pyunits.year,
        )
        b.sales_patenting_and_research_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Sales, patenting and research cost",
            units=CE_index_units / pyunits.year,
        )
        b.admin_and_support_labor_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Admin and support labor cost",
            units=CE_index_units / pyunits.year,
        )
        b.property_taxes_and_insurance_cost = Var(
            initialize=1,
            bounds=(0, None),
            doc="Property taxes and insurance cost",
            units=CE_index_units / pyunits.year,
        )
        b.total_fixed_OM_cost = Var(
            initialize=4,
            bounds=(0, None),
            doc="Total fixed O&M costs",
            units=CE_index_units / pyunits.year,
        )
        b.total_sales_revenue = Var(
            initialize=4,
            bounds=(0, None),
            doc="Total sales revenue",
            units=CE_index_units / pyunits.year,
        )

        # variable for user to assign other fixed costs to,
        # fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Other fixed costs",
            units=CE_index_units / pyunits.year,
        )
        b.other_fixed_costs.fix(1e-12)

        # variable for user to assign watertap fixed costs to,
        # fixed to 0 by default
        b.watertap_fixed_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Watertap fixed costs",
            units=CE_index_units / pyunits.year,
        )

        # variable for user to assign custom fixed costs to,
        # constraint sets to sum of list, which is 0 for empty list
        b.custom_fixed_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Custom fixed costs",
            units=CE_index_units / pyunits.year,
        )

        # create constraints
        TPC = b.total_plant_cost  # quick reference to total_plant_cost
        operating_labor_types, technical_labor_types = [], []  # subset labor lists
        for i in labor_types:
            if i in ["skilled", "unskilled", "supervisor", "maintenance"]:
                operating_labor_types.append(i)
            elif i in ["technician", "engineer"]:
                technical_labor_types.append(i)
            else:
                raise ValueError(
                    f"Value {i} for labor_type is not allowed. "
                    f"Allowed labor types for operating labor include skilled,"
                    f"unskilled, supervisor and maintenance. Allowed labor types "
                    f"for direct labor include technician and engineer."
                )

        # calculated from labor rate, labor burden, and operators per shift
        @b.Constraint()
        def annual_operating_labor_cost_eq(c):
            return c.annual_operating_labor_cost == pyunits.convert(
                (
                    sum(
                        c.operators_per_shift[i] * c.labor_rate[i]
                        for i in operating_labor_types
                    )
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units / pyunits.year,
            )

        @b.Constraint()
        def annual_technical_labor_cost_eq(c):
            return c.annual_technical_labor_cost == pyunits.convert(
                (
                    sum(
                        c.operators_per_shift[i] * c.labor_rate[i]
                        for i in technical_labor_types
                    )
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units / pyunits.year,
            )

        @b.Constraint()
        def annual_labor_cost_eq(c):
            return c.annual_labor_cost == pyunits.convert(
                (c.annual_operating_labor_cost + c.annual_technical_labor_cost),
                CE_index_units / pyunits.year,
            )

        # maintenance cost is 2% of TPC
        @b.Constraint()
        def maintenance_and_material_cost_eq(c):
            return c.maintenance_and_material_cost == 0.02 * TPC / pyunits.year

        # quality assurance cost is 10% of operating labor
        @b.Constraint()
        def quality_assurance_and_control_cost_eq(c):
            return c.quality_assurance_and_control_cost == 0.10 * pyunits.convert(
                (c.annual_operating_labor_cost),
                CE_index_units / pyunits.year,
            )

        # sales cost is 0.5% of total revenue
        @b.Constraint()
        def sales_patenting_and_research_cost_eq(c):
            return c.sales_patenting_and_research_cost == 0.005 * pyunits.convert(
                (c.total_sales_revenue),
                CE_index_units / pyunits.year,
            )

        # admin cost is 20% of direct labor
        @b.Constraint()
        def admin_and_support_labor_cost_eq(c):
            return c.admin_and_support_labor_cost == 0.20 * pyunits.convert(
                (c.annual_operating_labor_cost),
                CE_index_units / pyunits.year,
            )

        # taxes are 1% of TPC
        @b.Constraint()
        def taxes_and_insurance_cost_eq(c):
            return c.property_taxes_and_insurance_cost == 0.01 * TPC / pyunits.year

        # sum of fixed O&M costs

        # sum of fixed operating costs of watertap units
        @b.Constraint()
        def sum_watertap_fixed_costs(c):
            if len(c.watertap_fixed_costs_list) == 0:
                return c.watertap_fixed_costs == 1e-12 * CE_index_units / pyunits.year
            else:
                return c.watertap_fixed_costs == sum(b.watertap_fixed_costs_list)

        # sum of fixed operating costs of custom units
        @b.Constraint()
        def sum_custom_fixed_costs(c):
            if len(c.custom_fixed_costs_list) == 0:
                return c.custom_fixed_costs == 1e-12 * CE_index_units / pyunits.year
            else:
                return c.custom_fixed_costs == sum(b.custom_fixed_costs_list)

        @b.Constraint()
        def total_fixed_OM_cost_eq(c):
            return c.total_fixed_OM_cost == (
                c.annual_labor_cost
                + c.maintenance_and_material_cost
                + c.quality_assurance_and_control_cost
                + c.admin_and_support_labor_cost
                + c.sales_patenting_and_research_cost
                + c.property_taxes_and_insurance_cost
                + c.other_fixed_costs
                + c.watertap_fixed_costs
                + c.custom_fixed_costs
            )

        @b.Constraint()
        def total_sales_revenue_eq(c):
            return c.total_sales_revenue == pyunits.convert(
                (
                    (
                        sum(
                            pure_product_output_rates[p] * default_sale_prices[p]
                            for p in pure_product_output_rates.keys()
                        )
                        + c.mixed_product_sale_price_realization_factor
                        * sum(
                            mixed_product_output_rates[p] * default_sale_prices[p]
                            for p in mixed_product_output_rates.keys()
                        )
                    )
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units / pyunits.year,
            )

    def get_variable_OM_costs(
        b,
        resources,
        rates,
        prices=None,
        feed_input_rate=None,
        CE_index_year="2021",
        efficiency=0.85,
    ):
        """
        Args:
            b: costing block to add fixed cost variables and constraints to
            resources: list of strings for the resources to be costed
            rates: list of pyomo vars for resource consumption rates
            prices: dict of resource prices to be added to the premade dictionary
            feed_input_rate: rate of feedstock input
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
            efficiency: power usage efficiency, or fixed motor/distribution efficiency

        Returns:
            None.

        """
        if feed_input_rate is not None:
            b.feed_input_rate = value(feed_input_rate) * pyunits.get_units(
                feed_input_rate
            )

        if prices is None:
            prices = {}

        if not hasattr(b.parent_block(), "time"):  # flowsheet is not dynamic
            b.parent_block().time = [0]
        if not hasattr(b.parent_block(), "time_units"):  # no time units set
            b.parent_block().time_units = pyunits.s
        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2021."
            )

        # assert arguments are correct types
        if not isinstance(resources, list):
            raise TypeError("resources argument must be a list")
        if not isinstance(rates, list):
            raise TypeError("rates argument must be a list")
        if not isinstance(prices, dict):
            raise TypeError("prices argument must be a dictionary")

        # assert lists are the same length
        if len(resources) != len(rates):
            raise AttributeError("resources and rates must be lists of the same length")

        # dictionary of default prices
        # the currency units are millions of USD, so all prices need a 1e-6 multiplier to get USD
        default_prices = load_default_resource_prices()

        # add entries from prices to default_prices
        for key in prices.keys():
            default_prices[key] = prices[key]

        # raise error if the user included a resource not in default_prices
        if not set(resources).issubset(default_prices.keys()):
            raise AttributeError(
                f"A resource was included that does not contain a "
                f"price. Prices exist for the following resources: "
                f"{list(default_prices.keys())}"
            )

        # create list of prices
        prices = [default_prices[r] for r in resources]

        # zip rates and prices into a dict accessible by resource
        resource_rates = dict(zip(resources, rates))
        resource_prices = dict(zip(resources, prices))

        # make vars
        b.variable_operating_costs = Var(
            b.parent_block().time,
            resources,
            initialize=2e-7,
            doc="Variable operating costs",
            units=CE_index_units / pyunits.year,
        )

        b.other_variable_costs = Var(
            b.parent_block().time,
            initialize=0,
            bounds=(0, None),
            doc="A variable to include non-standard O&M costs",
            units=CE_index_units / pyunits.year,
        )

        # assume the user is not using this
        b.other_variable_costs.fix(1e-12)

        # TODO commented as no WaterTAP models currently use this, may change in the future
        # variable for user to assign watertap variable costs to,
        # constraint sets to sum of list, which is 0 for empty list
        # b.watertap_variable_costs = Var(
        #     initialize=0,
        #     bounds=(0, None),
        #     doc="Watertap variable costs",
        #     units=CE_index_units / pyunits.year,
        # )

        # variable for user to assign custom variable costs to,
        # constraint sets to sum of list, which is 0 for empty list
        b.custom_variable_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Custom variable costs",
            units=CE_index_units / pyunits.year,
        )

        b.total_variable_OM_cost = Var(
            b.parent_block().time,
            initialize=4e-6,
            doc="Total variable operating and maintenance costs",
            units=CE_index_units / pyunits.year,
        )

        @b.Constraint(b.parent_block().time, resources)
        def variable_cost_eq(c, t, r):
            if r == "power":
                efficiency_factor = efficiency  # fixed motor efficiency
            else:
                efficiency_factor = (
                    1  # other costs don't have this, could add more later
                )
            return c.variable_operating_costs[t, r] == (
                pyunits.convert(
                    resource_prices[r]
                    * resource_rates[r][t]
                    / efficiency_factor
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year,
                    to_units=CE_index_units / pyunits.year,
                )
            )

        # TODO commented as no WaterTAP models currently use this, may change in the future
        # sum of variable operating costs of watertap units
        # @b.Constraint()
        # def sum_watertap_variable_costs(c):
        #     if len(c.watertap_variable_costs_list) == 0:
        #         return c.watertap_variable_costs == 1e-12 * CE_index_units / pyunits.year
        #     else:
        #         return c.watertap_variable_costs == sum(b.watertap_variable_costs_list)

        # sum of variable operating costs of custom units
        @b.Constraint()
        def sum_custom_variable_costs(c):
            if len(c.custom_variable_costs_list) == 0:
                return c.custom_variable_costs == 1e-12 * CE_index_units / pyunits.year
            else:
                return c.custom_variable_costs == sum(b.custom_variable_costs_list)

        if hasattr(b, "total_fixed_OM_cost"):
            # define overhead cost
            # plant overhead, 20% of direct costs - fixed OM, power, water, lease/land, chemicals, waste
            b.plant_overhead_cost = Var(
                b.parent_block().time,
                initialize=0,
                doc="Plant overhead costs",
                units=CE_index_units / pyunits.year,
            )

        if (0, "power") in b.variable_operating_costs.id_index_map().values():

            @b.Constraint(b.parent_block().time)
            def plant_overhead_cost_eq(c, t):
                return c.plant_overhead_cost[t] == 0.20 * (
                    c.total_fixed_OM_cost
                    + c.variable_operating_costs[0, "power"]
                    + c.land_cost / pyunits.year
                    + sum(
                        c.variable_operating_costs[0, chemical]
                        for chemical in c.chemicals_list
                    )
                    + sum(
                        c.variable_operating_costs[0, waste] for waste in c.waste_list
                    )
                    + c.additional_chemicals_cost / pyunits.year
                    + c.additional_waste_cost / pyunits.year
                    # TODO commented as no WaterTAP models currently use this, may change in the future
                    # + c.watertap_variable_costs
                    + c.custom_variable_costs
                )

        else:

            @b.Constraint(b.parent_block().time)
            def plant_overhead_cost_eq(c, t):
                return c.plant_overhead_cost[t] == 0.20 * (
                    c.total_fixed_OM_cost
                    + c.land_cost / pyunits.year
                    + sum(
                        c.variable_operating_costs[0, chemical]
                        for chemical in c.chemicals_list
                    )
                    + sum(
                        c.variable_operating_costs[0, waste] for waste in c.waste_list
                    )
                    + c.additional_chemicals_cost / pyunits.year
                    + c.additional_waste_cost / pyunits.year
                    # TODO commented as no WaterTAP models currently use this, may change in the future
                    # + c.watertap_variable_costs
                    + c.custom_variable_costs
                )

        @b.Constraint(b.parent_block().time)
        def total_variable_cost_eq(c, t):
            return (
                c.total_variable_OM_cost[t]
                == sum(c.variable_operating_costs[t, r] for r in resources)
                + c.other_variable_costs[t]
                + c.plant_overhead_cost[t]
                + c.land_cost / pyunits.year
                + c.additional_chemicals_cost / pyunits.year
                + c.additional_waste_cost / pyunits.year
                # TODO commented as no WaterTAP models currently use this, may change in the future
                # + c.watertap_variable_costs
                + c.custom_variable_costs
            )

    def initialize_fixed_OM_costs(b):
        # b is the flowsheet-level costing block
        if hasattr(b, "total_fixed_OM_cost"):
            calculate_variable_from_constraint(
                b.annual_operating_labor_cost, b.annual_labor_cost_eq
            )

            calculate_variable_from_constraint(
                b.maintenance_and_material_cost, b.maintenance_and_material_cost_eq
            )

            calculate_variable_from_constraint(
                b.quality_assurance_and_control_cost,
                b.quality_assurance_and_control_cost_eq,
            )

            calculate_variable_from_constraint(
                b.sales_patenting_and_research_cost,
                b.sales_patenting_and_research_cost_eq,
            )

            calculate_variable_from_constraint(
                b.admin_and_support_labor_cost,
                b.admin_and_support_labor_cost_eq,
            )

            calculate_variable_from_constraint(
                b.property_taxes_and_insurance_cost,
                b.taxes_and_insurance_cost_eq,
            )

            calculate_variable_from_constraint(
                b.total_fixed_OM_cost, b.total_fixed_OM_cost_eq
            )

    def initialize_variable_OM_costs(b):
        # b is the flowsheet-level costing block
        # initialization for power generation costs
        if hasattr(b, "variable_operating_costs"):
            for i in b.variable_operating_costs.keys():
                if hasattr(b, "variable_cost_eq"):
                    calculate_variable_from_constraint(
                        b.variable_operating_costs[i],
                        b.variable_cost_eq[i],
                    )

            for i in b.total_variable_OM_cost.keys():
                calculate_variable_from_constraint(
                    b.total_variable_OM_cost[i],
                    b.total_variable_cost_eq[i],
                )

    # -----------------------------------------------------------------------------
    # Costing Library Utility Functions
    # -----------------------------------------------------------------------------

    def costing_initialization(b):
        # b is the flowsheet-level costing block
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in b._registered_unit_costing:
                for key in o.costing.bare_erected_cost.keys():
                    calculate_variable_from_constraint(
                        o.bare_erected_cost[key],
                        o.bare_erected_cost_eq[key],
                    )
                    calculate_variable_from_constraint(
                        o.total_plant_cost[key],
                        o.total_plant_cost_eq[key],
                    )
            # make sure all installation cost variables are initialized
            installation_cost_list = [
                "total_plant_cost",
                "bare_erected_cost",
                "total_installation_cost",
                "ancillary_cost",
                "piping_materials_and_labor_cost",
                "electrical_materials_and_labor_cost",
                "instrumentation_cost",
                "plant_services_cost",
                "buildings_cost",
                "process_buildings_cost",
                "auxiliary_buildings_cost",
                "site_improvements_cost",
                "epcm_cost",
                "equipment_installation_cost",
                "field_expenses_cost",
                "project_management_and_construction_cost",
                "process_contingency_cost",
                "contingency_cost",
            ]
            for var in installation_cost_list:
                if hasattr(b, var):
                    calculate_variable_from_constraint(
                        getattr(b, var), getattr(b, var + "_eq")
                    )

    def display_total_plant_costs(b):
        print("-----Total Plant Costs-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "total_plant_cost"):
                print(
                    "%s: %.2f"
                    % (
                        value(o.name),
                        value(
                            sum(o.total_plant_cost[key] for key in o.total_plant_cost)
                        ),
                    )
                )

    def display_bare_erected_costs(b):
        print("-----Bare Erected Costs-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "bare_erected_cost"):
                print(
                    "%s: %.5f"
                    % (
                        value(o.name),
                        value(
                            sum(o.bare_erected_cost[key] for key in o.bare_erected_cost)
                        ),
                    )
                )

    def get_total_BEC(b, CE_index_year, watertap_blocks=None):
        # This method accepts a flowsheet-level costing block

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2020."
            )

        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [block.name for block in b._registered_unit_costing]:
                if hasattr(o, "bare_erected_cost"):  # added from cost accounts
                    for key in o.bare_erected_cost.keys():
                        b.BEC_list.append(
                            pyunits.convert(
                                o.bare_erected_cost[key], to_units=CE_index_units
                            )
                        )
                elif hasattr(o, "capital_cost"):  # added from custom model
                    b.BEC_list.append(
                        pyunits.convert(o.capital_cost, to_units=CE_index_units)
                    )
                    if hasattr(o, "fixed_operating_cost"):
                        b.custom_fixed_costs_list.append(
                            pyunits.convert(
                                o.fixed_operating_cost,
                                to_units=CE_index_units / pyunits.year,
                            )
                        )
                    if hasattr(o, "variable_operating_cost"):
                        b.custom_variable_costs_list.append(
                            pyunits.convert(
                                o.variable_operating_cost,
                                to_units=CE_index_units / pyunits.year,
                            )
                        )

        if watertap_blocks is not None:  # added from WaterTAP
            for w in watertap_blocks:
                m = ConcreteModel()
                m.fs = FlowsheetBlock(dynamic=False)
                if issubclass(w._ComponentDataClass, ZeroOrderBaseData):
                    m.fs.costing = ZeroOrderCosting()
                else:
                    m.fs.costing = WaterTAPCosting()
                w.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
                b.BEC_list.append(
                    pyunits.convert(w.costing.capital_cost, to_units=CE_index_units)
                )
                if hasattr(w.costing, "fixed_operating_cost"):
                    b.watertap_fixed_costs_list.append(
                        pyunits.convert(
                            w.costing.fixed_operating_cost,
                            to_units=CE_index_units / pyunits.year,
                        )
                    )
                # # TODO commented as no WaterTAP models currently use this, may change in the future
                # if hasattr(w.costing, "variable_operating_cost"):
                #     b.watertap_variable_costs_list.append(
                #         pyunits.convert(
                #             w.costing.variable_operating_cost,
                #             to_units=CE_index_units / pyunits.year,
                #         )
                #     )

        b.total_BEC = Var(
            initialize=100,
            bounds=(0, None),
            doc="Total TPC",
            # assume that total_plant_cost is in millions of
            # USD_year, where year is the CE_index_year users set
            units=CE_index_units,
        )

        @b.Constraint()
        def total_BEC_eq(c):
            return c.total_BEC == sum(b.BEC_list)

    def display_flowsheet_cost(b):
        # This method accepts a flowsheet-level costing block
        print("\n")
        print("Total bare erected cost: %.3f" % value(b.total_BEC))
        if hasattr(b, "total_overnight_capital"):
            print(
                "Total overnight (installed) equipment cost: %.3f"
                % value(b.total_overnight_capital)
            )
        if hasattr(b, "annualized_cost"):
            print("Total annualized capital cost: %.3f" % value(b.annualized_cost))
        print()
        if hasattr(b, "total_fixed_OM_cost"):
            print("Total annual fixed O&M cost: %.3f" % value(b.total_fixed_OM_cost))
        if hasattr(b, "total_variable_OM_cost"):
            print(
                "Total annual variable O&M cost: %.3f"
                % value(b.total_variable_OM_cost[0])
            )
        if hasattr(b, "total_fixed_OM_cost") and hasattr(b, "total_variable_OM_cost"):
            print(
                "Total annual O&M cost: %.3f"
                % value(b.total_fixed_OM_cost + b.total_variable_OM_cost[0])
            )
            if hasattr(b, "feed_input_rate"):
                print(
                    "Total annual O&M cost per ton feed processed: %.3f"
                    % value(
                        (b.total_fixed_OM_cost + b.total_variable_OM_cost[0])
                        * 1e6
                        / (
                            pyunits.convert(
                                b.feed_input_rate, to_units=pyunits.ton / pyunits.hr
                            )
                            * b.hours_per_shift
                            * b.shifts_per_day
                            * b.operating_days_per_year
                        )
                    )
                )
            if hasattr(b, "recovery_rate_per_year"):
                print(
                    "Total annual O&M cost per kg REE recovered: %.3f"
                    % value(
                        (b.total_fixed_OM_cost + b.total_variable_OM_cost[0])
                        * 1e6
                        / (
                            pyunits.convert(
                                b.recovery_rate_per_year,
                                to_units=pyunits.kg / pyunits.year,
                            )
                        )
                    )
                )
        print()
        if (
            hasattr(b, "annualized_cost")
            and hasattr(b, "total_fixed_OM_cost")
            and hasattr(b, "total_variable_OM_cost")
        ):
            print(
                "Total annualized plant cost: %.3f"
                % value(
                    b.annualized_cost
                    + (b.total_fixed_OM_cost + b.total_variable_OM_cost[0])
                    * pyunits.year
                )
            )
        if hasattr(b, "recovery_rate_per_year"):
            print(
                "Annual rate of recovery: %.3f kg/year REE recovered"
                % value(
                    pyunits.convert(
                        b.recovery_rate_per_year, to_units=pyunits.kg / pyunits.year
                    )
                )
            )
        if hasattr(b, "cost_of_recovery"):
            print(
                "Cost of recovery per kg REE recovered: %.3f"
                % value(b.cost_of_recovery)
            )
        print()

        if hasattr(b, "npv"):
            print("Net present value: %.3f" % value(b.npv))
        print("\n")

    def calculate_REE_costing_bounds(
        b, capacity, grade, CE_index_year, recalculate=False
    ):
        # adapted from https://doi.org/10.1038/s41893-023-01145-1
        # This method accepts a flowsheet-level costing block
        # capacity and grade should be variables with Pyomo units,
        # or values with Pyomo unit containers
        # CE_index_year should be a string currency unit, e.g. "2021"
        # recalculate tells method to delete and rebuild components

        if recalculate and hasattr(b, "components_already_built"):
            delattr(b, "capacity")
            delattr(b, "grade")
            delattr(b, "costing_lower_bound")
            delattr(b, "costing_upper_bound")
            delattr(b, "costing_lower_bound_eq")
            delattr(b, "costing_upper_bound_eq")

        if not hasattr(b, "capacity"):
            b.capacity = Var(
                initialize=value(pyunits.convert(capacity, to_units=pyunits.tonnes)),
                bounds=(0, None),
                doc="Feedstock capacity of site",
                units=pyunits.tonnes,
            )
            b.capacity.fix(capacity)
            _log.info(
                "New variable 'capacity' created as attribute of {}".format(b.name)
            )
        else:
            _log.info(
                "Flowsheet-level costing block {} already has attribute "
                "'capacity', moving on. Set 'recalculate' to True to delete "
                "old objects and recalculate for new inputs".format(b.name)
            )

        if not hasattr(b, "grade"):
            b.grade = Var(
                initialize=value(pyunits.convert(grade, to_units=pyunits.percent)),
                bounds=(0, None),
                doc="Grade percentage of site. The value should be a "
                "percentage, for example 10 for 10%.",
                units=pyunits.percent,
            )
            b.grade.fix(grade)
            _log.info("New variable 'grade' created as attribute of {}".format(b.name))
        else:
            _log.info(
                "Flowsheet-level costing block {} already has attribute "
                "'grade', moving on. Set 'recalculate' to True to delete "
                "old objects and recalculate for new inputs".format(b.name)
            )

        processes = {
            "Total Capital": [81, 1.4, -0.46, 0.063],
            "Total Operating": [27, 0.87, -0.087, 0.038],
            "Beneficiation": [2.7, 1.3, -0.15, 0.062],
            "Beneficiation, Chemical Extraction, Enrichment and Separation": [
                22,
                1.28,
                -0.059,
                0.046,
            ],
            "Chemical Extraction": [40, 2.9, -0.46, 0.14],
            "Chemical Extraction, Enrichment and Separation": [15, 15, -0.19, 0.28],
            "Enrichment and Separation": [6.7, 2.8, -0.16, 0.11],
            "Mining": [25, 2.5, -0.32, 0.095],
        }

        if not hasattr(b, "costing_lower_bound"):
            b.costing_lower_bound = Var(
                processes,
                initialize=1,
                bounds=(0, None),
                doc="Estimated lower bound on per unit production cost of site",
                units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            )
            _log.info(
                "New variable 'costing_lower_bound' created as attribute of {}".format(
                    b.name
                )
            )
        else:
            _log.info(
                "Flowsheet-level costing block {} already has attribute "
                "'costing_lower_bound', moving on. Set 'recalculate' to True to delete "
                "old objects and recalculate for new inputs".format(b.name)
            )

        if not hasattr(b, "costing_upper_bound"):
            b.costing_upper_bound = Var(
                processes,
                initialize=1,
                bounds=(0, None),
                doc="Estimated upper bound on per unit production cost of site",
                units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            )
            _log.info(
                "New variable 'costing_upper_bound' created as attribute of {}".format(
                    b.name
                )
            )
        else:
            _log.info(
                "Flowsheet-level costing block {} already has attribute "
                "'costing_upper_bound', moving on. Set 'recalculate' to True to delete "
                "old objects and recalculate for new inputs".format(b.name)
            )

        if not hasattr(b, "costing_lower_bound_eq"):

            @b.Constraint(processes)
            def costing_lower_bound_eq(c, p):
                return (
                    c.costing_lower_bound[p]
                    == pyunits.convert(
                        pyunits.USD_2022
                        * (processes[p][0] - processes[p][1])
                        * (
                            pyunits.convert(c.grade, to_units=pyunits.dimensionless)
                            * pyunits.convert(c.capacity, to_units=pyunits.tonnes)
                            / pyunits.tonnes
                        )
                        ** (processes[p][2] - processes[p][3]),
                        to_units=getattr(pyunits, "USD_" + CE_index_year),
                    )
                    / pyunits.kg
                )

            # assume model is already solved, so just calculate costing bounds here
            for i in b.costing_lower_bound.keys():
                calculate_variable_from_constraint(
                    b.costing_lower_bound[i],
                    b.costing_lower_bound_eq[i],
                )
            _log.info(
                "New constraint 'costing_lower_bounding_eq' created as attribute of {}".format(
                    b.name
                )
            )
        else:
            _log.info(
                "Flowsheet-level costing block {} already has indexed "
                "constraint 'costing_lower_bounding_eq', reporting existing results. "
                "Set 'recalculate' to True to delete old objects and recalculate "
                "for new inputs".format(b.name)
            )

        if not hasattr(b, "costing_upper_bound_eq"):

            @b.Constraint(processes)
            def costing_upper_bound_eq(c, p):
                return (
                    c.costing_upper_bound[p]
                    == pyunits.convert(
                        pyunits.USD_2022
                        * (processes[p][0] + processes[p][1])
                        * (
                            pyunits.convert(c.grade, to_units=pyunits.dimensionless)
                            * pyunits.convert(c.capacity, to_units=pyunits.tonnes)
                            / pyunits.tonnes
                        )
                        ** (processes[p][2] + processes[p][3]),
                        to_units=getattr(pyunits, "USD_" + CE_index_year),
                    )
                    / pyunits.kg
                )

            # assume model is already solved, so just calculate costing bounds here
            for i in b.costing_upper_bound.keys():
                calculate_variable_from_constraint(
                    b.costing_upper_bound[i],
                    b.costing_upper_bound_eq[i],
                )
            _log.info(
                "New constraint 'costing_upper_bounding_eq' created as attribute of {}".format(
                    b.name
                )
            )
        else:
            _log.info(
                "Flowsheet-level costing block {} already has indexed "
                "constraint 'costing_upper_bounding_eq', reporting existing results. "
                "Set 'recalculate' to True to delete old objects and recalculate "
                "for new inputs".format(b.name)
            )

        _log.info("\n\nPrinting calculated costing bounds for processes:")
        for p in processes:
            print(
                p,
                ": [",
                value(b.costing_lower_bound[p]),
                ", ",
                value(b.costing_upper_bound[p]),
                "]",
                getattr(pyunits, "USD_" + CE_index_year),
                "/kg",
            )

        # method has finished building components
        b.components_already_built = True

    def calculate_NPV(b, fixed_OM, variable_OM, consider_taxes=False):
        """
        Equations for cash flow expressions derived from the textbook
        Engineering Economy: Applying Theory to Practice, 3rd Ed. by Ted. G. Eschenbach.

        The net present value (NPV) is a representative measure of the "current
        day" value of a chemical plant over the total lifetime, including all
        cash flows.

        This method supports capital expenditure, loan repayment, inflation,
        and royalties. The NPV formulation assumes that negative cash flows
        consists capital and operating costs scaled to a constant present
        value. The general NPV formula with 100% of capital expenditure
        upfront at the start of the operating period and no capital or
        operating growth rate is

        NPV = [(REVENUE - OPEX - ROYALTIES) * P/A(r, N)] - CAPEX

        where P/A(r, N) is the series present worth factor; this factor scales
        a future cost to its present value from a known discount rate r and
        project lifetime N based on annuity growth over time. This factor is
        calculated as

        P/A(r, N) = [ 1 - (1+r)**(-N) ] / r

        where r is expressed as a decimal and N is expressed in years. In the
        NPV expression above, REVENUE is the constant annual revenue, OPEX is
        the constant annual operating cost, and CAPEX is the total capital cost.
        ROYALTIES are charged based on revenue, usually a fixed percentage.

        Operating costs and revenues are adjusted based on predicted annuity
        growth to obtain the present value. These expressions are implemented
        if there is no capital expenditure period or additional growth rate.

        Often, uniform series of cash flows do not start in the first period
        (t=1), but in some later period, T, after a known delay. In this case
        the delay is accounted for using

        PV_year1_cashflow = Cash_Flow_Value * P/A(r, N)

        PV_yearT_cashflow = Cash_Flow_Value * [P/A(r, N) - P/A(r, T-1)]

        ----------------------------------------------------------------------

        The general NPV formulation allows capital costs, operating costs,
        and revenues to escalate over time via geometric gradient growth, e.g.
        a constant proportional growth rate expressed as an escalation or
        inflation percentage. The general formulation is given by

        NPV = PV_Revenue - PV_Operating_Cost - PV_Royalties
              - PV_Capital_Cost - PV_Loan_Interest

        For costs escalating at a constant rate for the project lifetime, the
        series present worth factor is modified to account for escalation,
        yielding a modified formula

        P/A(r, g, N) = ( 1 - [ (1+g)**(N) ] * [(1+r)**(-N)] ) / (r - g)

        where r is the discount rate expressed as a decimal, N is the project
        lifetime, and g is the escalation rate (e.g. inflation) expressed as a
        decimal.

        The general formulation considers a capital escalation period followed
        by the operating period, and capital costs may be distributed across
        the capital escalation period rather than fully paid for upfront. For
        example, if the capital expenditures are distributed across a 3-year
        capital escalation period, the PV from the capital costs are given as

        PV_Capital_Cost = Y1_% * CAPEX * [P/A(r, gCap, 1) - P/A(r, gCap, 0)]
                          + Y2_% * CAPEX * [P/A(r, gCap, 2) - P/A(r, gCap, 1)]
                          + Y3_% * CAPEX * [P/A(r, gCap, 3) - P/A(r, gCap, 2)]

        where Y1_%, Y2_%, and Y3_% are the percentages of capital expenditure
        in each year expressed as decimals, CAPEX is the total capital cost from
        equipment purchasing, gCap is the capital escalation growth rate
        expressed as a decimal. The capital costs spent in each year are handled
        separately to properly account for the value growth over time. Loan
        repayment and interest owed are calculated as

        PV_Loan_Interest_Owed = Debt * [P/A(r, 0, Nloan) / P/A(iLoan_%, 0, Nloan) - 1]

        where Debt is the loan principal (typically a percentage of the CAPEX), Nloan
        is the loan repayment period, and iLoan_% is the capital equipment loan interest
        rate expressed as a decimal.

        Revenue, operating costs and royalties based on revenue escalate with
        standard inflation. Notably, these cash flows occur after any capital
        expenditure period, meaning that the annuity growth must be offset by
        the length of the capital expenditure period. This yields the expressions

        PV_Revenue = REVENUE * [ P/A(r, gRev, NOp+NCap) - P/A(r, gRev, NCap) ]

        PV_Operating_Cost = OPEX * [ P/A(r, gOp, NOp+NCap) - P/A(r, gOp, NCap) ]

        PV_Royalties = iRoy_% * REVENUE * [ P/A(r, gRev, NOp+NCap) - P/A(r, gRev, NCap) ]

        where REVENUE is the annual revenue, OPEX is the annual operating cost,
        gRev is the inflation or growth rate of revenue year-on-year expressed as
        a decimal, gOp is the inflation or growth rate of operating costs year-on-year
        expressed as a decimal, NOp is the length of the operating period or plant
        lifetime, NCap is the length of the capital expenditure period, and iRoy_%
        is the percentage of the revenue charged as royalties expressed as a decimal.
        The expressions above take the annuity growth during the entire analysis
        period (NOp+NCap) and subtract the capital expenditure period (NCap) as there
        is no operation or production during that time.

        Args:
            b: costing block to retrieve total plant cost (capital), total plant
                operating cost, and total revenue from, and add net present
                value (NPV) calculations to
            fixed_OM: True/False flag for calculating fixed O&M costs
            variable_OM: True/False flag for calculating variable O&M costs
        """

        # input verification

        input_list = [
            b.config.total_capital_cost,
            b.config.annual_operating_cost,
            b.config.annual_revenue,
            b.config.cost_year,
        ]

        if True in [
            i is not None for i in input_list
        ]:  # if one fixed input is passed, require all fixed inputs
            b.verify_calculate_from_inputs()
        elif not (
            fixed_OM and variable_OM
        ):  # if OM variables don't exist, fail and suggest ways to fix error
            raise AttributeError(
                "If capital, fixed O&M, or variable O&M costs are not calculated, "
                "then inputs for total_capital_cost, annual_operating_cost, and annual_revenue "
                "must be passed, and cost_year must be passed as a string, e.g. '2021'."
                "Alternatively, set fixed_OM and variable_OM to True to calculate O&M results."
            )
        else:  # otherwise use expected results from costing block
            b.verify_calculate_from_costing_block()

        # check required arguments
        QGESSCostingData.assert_config_argument_set(b, name="discount_percentage")
        QGESSCostingData.assert_config_argument_set(b, name="plant_lifetime")

        # check capital expenditure arguments
        if (
            b.config.has_capital_expenditure_period
            and b.config.capital_expenditure_percentages is None
        ):
            b.config.capital_expenditure_percentages = [10, 60, 30]
        elif not b.config.has_capital_expenditure_period:
            b.config.capital_expenditure_percentages = []

        if b.config.has_capital_expenditure_period:
            QGESSCostingData.verify_percentages_list(
                b.config, name="capital_expenditure_percentages"
            )

        # check optional expressions
        QGESSCostingData.assert_Pyomo_object(b.config, name="debt_expression")

        # build variables

        b.pv_capital_cost = Var(
            initialize=-b.CAPEX,
            bounds=(None, 0),
            doc="Present value of total lifetime capital costs; negative cash flow",
            units=b.cost_units,
        )

        b.loan_debt = Var(
            initialize=b.CAPEX,
            bounds=(0, 1e4),
            doc="total debt from loans in $MM",
            units=b.cost_units,
        )

        b.pv_loan_interest = Var(
            initialize=-b.CAPEX,
            bounds=(-1e4, 1e4),
            doc="present value of total lifetime loan interest in $MM; normally a negative cash flow, but can be positive depending on the discount and interest rates",
            units=b.cost_units,
        )

        b.pv_operating_cost = Var(
            initialize=-b.OPEX * b.config.plant_lifetime,
            bounds=(None, 0),
            doc="Present value of total lifetime operating costs; negative cash flow",
            units=b.cost_units,
        )

        b.pv_revenue = Var(
            initialize=b.REVENUE * b.config.plant_lifetime,
            bounds=(0, None),
            doc="Present value of total lifetime sales revenue; positive cash flow",
            units=b.cost_units,
        )

        if consider_taxes:
            b.pv_taxes = Var(
                initialize=-b.net_tax_owed * pyunits.year * b.config.plant_lifetime,
                bounds=(None, 0),
                doc="Present value of total lifetime tax owed; negative cash flow",
                units=b.cost_units,
            )

        b.npv = Var(
            initialize=(-b.CAPEX + (b.REVENUE - b.OPEX) * b.config.plant_lifetime),
            bounds=(None, None),
            doc="Present value of plant over entire capital and operation lifetime",
            units=b.cost_units,
        )

        # build parameters

        b.discount_percentage = Param(
            initialize=b.config.discount_percentage, units=pyunits.percent
        )
        b.plant_lifetime = Param(
            initialize=b.config.plant_lifetime, units=pyunits.years
        )

        if b.config.has_capital_expenditure_period:
            b.capital_expenditure_percentages = Param(
                range(len(b.config.capital_expenditure_percentages)),
                initialize=dict(
                    zip(
                        range(len(b.config.capital_expenditure_percentages)),
                        b.config.capital_expenditure_percentages,
                    )
                ),
            )

        b.capital_escalation_percentage = Param(
            initialize=b.config.capital_escalation_percentage, units=pyunits.percent
        )

        b.capital_loan_interest_percentage = Param(
            initialize=b.config.capital_loan_interest_percentage, units=pyunits.percent
        )

        b.capital_loan_repayment_period = Param(
            initialize=b.config.capital_loan_repayment_period, units=pyunits.years
        )

        b.debt_percentage_of_CAPEX = Param(
            initialize=b.config.debt_percentage_of_CAPEX, units=pyunits.percent
        )

        b.operating_inflation_percentage = Param(
            initialize=b.config.operating_inflation_percentage, units=pyunits.percent
        )

        b.revenue_inflation_percentage = Param(
            initialize=b.config.revenue_inflation_percentage, units=pyunits.percent
        )

        # define series present worth factor as an method so it can be called

        def series_present_worth_factor(r, g, N):
            """
            Returns expression for series present worth factor where r is the discount rate
            expressed as a decimal, N is the project lifetime, and g is the escalation rate
            (e.g. inflation) expressed as a decimal.
            """
            return (1 - ((1 + g) ** (N)) * ((1 + r) ** (-N))) / (r - g)

        # build constraints

        if b.config.has_capital_expenditure_period:

            @b.Constraint()
            def pv_capital_cost_constraint(c):
                # percentage of CAPEX is basis for each capital expenditure year
                # since the expenditure series restarts in each year, we need to split
                # the terms for each year out and subtract off the delayed years
                # PV_Capital_Cost = - (
                # %year1 * CAPEX * (P/A_year1 - P/A_year0)     change from year 1 only
                # + %year2 * CAPEX * (P/A_year2 - P/A_year1)   change from year 2 only
                # + %year3 * CAPEX * (P/A_year3 - P/A_year2)   change from year 2 only
                # + ...)
                # P/A_year0 = 0, which places each CAPEX expenditure at the end of each period

                return c.pv_capital_cost == -pyunits.convert(
                    sum(
                        pyunits.convert(
                            c.config.capital_expenditure_percentages[idx]
                            * pyunits.percent,
                            to_units=pyunits.dimensionless,
                        )
                        * c.CAPEX
                        * (  # P/A_year(i) - P/A_year(i-1))
                            series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.capital_escalation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx + 1,
                            )
                            - series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.capital_escalation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx,
                            )
                        )
                        for idx in range(len(c.config.capital_expenditure_percentages))
                    ),
                    to_units=c.cost_units,
                )

        else:

            @b.Constraint()
            def pv_capital_cost_constraint(c):
                # no expenditure period, so cash flow occurs at t=0 (project year)
                # PV_Capital_Cost = - CAPEX

                return c.pv_capital_cost == -pyunits.convert(
                    c.CAPEX,
                    to_units=c.cost_units,
                )

        if b.config.debt_expression is None:

            @b.Constraint()
            def loan_debt_constraint(c):
                # Debt  = %debt_charge_of_CAPEX * CAPEX

                return c.loan_debt == pyunits.convert(
                    pyunits.convert(
                        c.debt_percentage_of_CAPEX, to_units=pyunits.dimensionless
                    )
                    * c.CAPEX,
                    to_units=c.cost_units,
                )

        else:

            b.loan_debt = Reference(b.config.debt_expression)

        @b.Constraint()
        def pv_loan_interest_constraint(c):
            # PV_Loan_Interest_Owed = Debt * [P/A(r, 0, loan_length) / P/A(%interest, 0, loan_length) - 1]
            # when r > %interest, this is a negative value; the loan amount borrowed devalues faster than
            # the loan amount repaid
            # when r < %interest, this is a positive value; the loan amount borrowed devalues slower than
            # the loan amount repaid

            if c.config.debt_expression is None:

                return c.pv_loan_interest == pyunits.convert(
                    c.loan_debt
                    * (
                        series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            0,  # loan value does not grow over time
                            c.capital_loan_repayment_period / pyunits.year,
                        )
                        / series_present_worth_factor(
                            pyunits.convert(
                                c.capital_loan_interest_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            0,  # loan payments do not grow over time
                            c.capital_loan_repayment_period / pyunits.year,
                        )
                        - 1
                    ),
                    to_units=c.cost_units,
                )

            else:

                return c.pv_loan_interest == pyunits.convert(
                    c.loan_debt[None]
                    * (
                        series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            0,  # loan value does not grow over time
                            c.capital_loan_repayment_period / pyunits.year,
                        )
                        / series_present_worth_factor(
                            pyunits.convert(
                                c.capital_loan_interest_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            0,  # loan payments do not grow over time
                            c.capital_loan_repayment_period / pyunits.year,
                        )
                        - 1
                    ),
                    to_units=c.cost_units,
                )

        @b.Constraint()
        def pv_operating_cost_constraint(c):
            # OPEX starts after the capital expenditure period, so we need to account for a delay
            # PV_Operating_Cost = - OPEX * [ P/A(r, g, OPEX_end_year) - P/A(r, g, CAPEX_end_year) ]

            return c.pv_operating_cost == -pyunits.convert(
                c.OPEX
                * (
                    series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.operating_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        c.plant_lifetime / pyunits.year
                        + len(c.config.capital_expenditure_percentages),
                    )
                    - series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.operating_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        len(c.config.capital_expenditure_percentages),
                    )
                ),
                to_units=c.cost_units,
            )

        @b.Constraint()
        def pv_revenue_constraint(c):
            # Revenue starts after the capital expenditure period, so we need to account for a delay
            # PV_Revenue = - REVENUE * [ P/A(r, g, Revenue_end_year) - P/A(r, g, CAPEX_end_year) ]

            return c.pv_revenue == pyunits.convert(
                c.REVENUE
                * (
                    series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.revenue_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        c.plant_lifetime / pyunits.year
                        + len(c.config.capital_expenditure_percentages),
                    )
                    - series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.revenue_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        len(c.config.capital_expenditure_percentages),
                    )
                ),
                to_units=c.cost_units,
            )

        if consider_taxes:

            @b.Constraint()
            def pv_taxes_constraint(c):
                # Taxes start after the capital expenditure period, so we need to account for a delay
                # PV_taxes = net_tax_owed * [ P/A(r, 0, Operating_end_year) - P/A(r, 0, CAPEX_end_year) ]

                return c.pv_taxes == -pyunits.convert(
                    c.net_tax_owed
                    * pyunits.year
                    * (
                        series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage, to_units=pyunits.dimensionless
                            ),
                            pyunits.convert(
                                0,
                                to_units=pyunits.dimensionless,
                            ),
                            c.plant_lifetime / pyunits.year
                            + len(c.config.capital_expenditure_percentages),
                        )
                        - series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage, to_units=pyunits.dimensionless
                            ),
                            pyunits.convert(
                                0,
                                to_units=pyunits.dimensionless,
                            ),
                            len(c.config.capital_expenditure_percentages),
                        )
                    ),
                    to_units=c.cost_units,
                )

        @b.Constraint()
        def npv_constraint(c):

            return c.npv == pyunits.convert(
                c.pv_revenue
                + c.pv_capital_cost
                + c.pv_loan_interest
                + c.pv_operating_cost
                + (c.pv_taxes if consider_taxes else 0 * c.cost_units),
                to_units=c.cost_units,
            )

    def verify_calculate_from_costing_block(b):
        """
        Verify that parent block for NPV calculations has expected attributes.
        """
        try:
            b.CAPEX = b.total_BEC + b.total_installation_cost + b.other_plant_costs
            b.OPEX = (
                b.total_fixed_OM_cost * pyunits.year
                + b.total_variable_OM_cost[0] * pyunits.year
                + b.land_cost
            )
            b.REVENUE = Reference(b.total_sales_revenue)[None] * pyunits.year

            b.cost_units = pyunits.get_units(b.CAPEX)

        except AttributeError:
            raise AttributeError(
                "Expected FlowsheetCostingBlockData object "
                "with attributes total_BEC, total_installation_cost, "
                "total_fixed_OM_cost, total_variable_OM_cost, "
                "other_plant_costs, land_cost, and total_sales_revenue. "
                "Please confirm that b is a FlowsheetCostingBlockData object "
                "and that all expected attributes exist."
            )

    def verify_calculate_from_inputs(b):
        """
        Verify that expected inputs are set.
        """
        # if b is not a costing block, it must be a flowsheet block
        # variables and constraints will be added there
        for attr in [
            b.config.total_capital_cost,
            b.config.annual_operating_cost,
            b.config.annual_revenue,
            b.config.cost_year,
        ]:
            if attr is None:
                raise AttributeError(
                    "If capital, fixed O&M, or variable O&M costs are not calculated, "
                    "then inputs for total_capital_cost, annual_operating_cost, and annual_revenue "
                    "must be passed, and cost_year must be passed as a string, e.g. '2021'."
                )

        for attr in [
            "total_capital_cost",
            "annual_operating_cost",
            "annual_revenue",
        ]:
            # config value can be a number (int or float), or a Pyomo object
            if not type(getattr(b.config, attr)) in [int, float]:
                QGESSCostingData.assert_Pyomo_object(b.config, name=attr)

        else:
            # check if the cost arguments are variables or expressions with units and handle appropriately
            costs = {
                "CAPEX": b.config.total_capital_cost,
                "OPEX": b.config.annual_operating_cost,
                "REVENUE": b.config.annual_revenue,
            }
            b.cost_units = getattr(pyunits, "MUSD_" + b.config.cost_year)

            for key in costs.keys():
                # check if the object is a Reference
                if isinstance(costs[key], Component):
                    # it's a Pyomo object
                    if costs[key].is_reference():
                        costs[key] = costs[key][None]

                if type(costs[key]) in [Expression, ScalarExpression]:
                    if pyunits.get_units(costs[key]) == pyunits.dimensionless:
                        costs[key] = Expression(expr=costs[key].expr * b.cost_units)
                    else:
                        costs[key] = Expression(
                            expr=pyunits.convert(costs[key].expr, to_units=b.cost_units)
                        )
                else:
                    if pyunits.get_units(costs[key]) == pyunits.dimensionless:
                        costs[key] = Expression(expr=costs[key] * b.cost_units)
                    else:
                        costs[key] = Expression(
                            expr=pyunits.convert(costs[key], to_units=b.cost_units)
                        )

            # store for later use
            b.CAPEX = costs["CAPEX"]
            b.OPEX = costs["OPEX"]
            b.REVENUE = costs["REVENUE"]

    def assert_config_argument_set(b, name):
        """
        Verify that required arguments are set.
        """
        if getattr(b.config, name) is None:
            raise AttributeError(f"Required argument {name} not set")

    def verify_percentages_list(b, name):
        """
        Verify that percentage lists have expected properties.
        """
        percentages_list = getattr(b, name)

        if not isinstance(percentages_list, list):
            raise TypeError(
                f"{percentages_list} is not a list. "
                f"Argument {name} must be passed as a list."
            )
        if len(percentages_list) == 0:
            raise AttributeError(
                f"Argument {name} has a length of "
                "zero. List must have a nonzero length."
            )

        if not sum(percentages_list) == 100:
            raise AttributeError(
                f"Argument {name} has a sum of "
                f"{sum(percentages_list)}. List must sum to 100 percent."
            )

    def assert_Pyomo_object(b, name):
        """
        Verify that objects are of a supported type.
        """
        if getattr(b, name) is not None:
            obj = getattr(b, name)
            if not (
                isinstance(obj, Param)
                or isinstance(obj, Var)
                or isinstance(obj, Expression)
                or isinstance(obj, ScalarExpression)
            ):
                raise TypeError(
                    f"Argument {obj} of type {type(obj)} is not a supported object type. "
                    f"Ensure {name} is a Pyomo Param, Var, Expression, or ScalarExpression."
                )

    def economy_of_numbers(
        blk, cum_num_units, cost_FOAK, CE_index_year, learning_rate=0.04
    ):
        """
        Economy of Numbers (EoN) estimates the future profitability of novel/First-of-A-Kind (FOAK)
        equipment. This is because the cost of manufacturing a piece of equipment tends to decline
        as the cumulative production quantity rises, resulting from a consistent improvement in
        technical know-how.

        Y = A(X^-b)

        b = - log(1-R)/log(2)
        where Y is the cost of the Nth-of-A-Kind (NOAK) of the equipment, A is the cost of the FOAK,
        X is the cumulative number of units, b is the learning rate exponent, and R is the learning
        rate constant.

        The equations above are derived from Faber G, Ruttinger A, Strunge T, Langhorst T, Zimmermann A,
        van der Hulst M, Bensebaa F, Moni S and Tao L (2022) Adapting Technology Learning Curves for
        Prospective Techno-Economic and Life Cycle Assessments of Emerging Carbon Capture and Utilization
        Pathways. Front. Clim. 4:820261. doi: 10.3389/fclim.2022.820261

        Args:
            cum_num_units: The cumulative number of units.
            cost_FOAK: The cost of manufacturing the First-of-A-Kind equipment.
            CE_index_year: year for cost basis, e.g., "2021" to use 2021 dollars
            learning_rate: ranges between 0.01 - 0.1, depending on the level of maturity
                           (i.e., experimental, growing, proven, etc.)
                            as described in Rubin, E. S., Mantripragada, H., and Zhai, H.,
                            "An Assessment of the NETL Cost Estimation Methodology".
                            Department of Engineering and Public Policy, Carnegie Mellon University,
                            Pittsburgh, PA (2016). p. 31, Fig. 6-4.
        """

        blk.cum_num_units = Param(
            initialize=cum_num_units,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Cumulative number of units produced",
        )
        blk.learning_rate = Param(
            initialize=learning_rate,
            mutable=True,
            units=pyunits.dimensionless,
            doc="The learning factor reflects the level of maturity of the unit/technology",
        )

        blk.cost_NOAK = Var(
            initialize=1e5,
            bounds=(0, None),
            doc="Cost of the Nth-of-A-Kind of the unit",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        @blk.Expression(
            doc="This measures the rate at which the cost is reduced as cumulative units increases"
        )
        def learning_rate_exponent(b):

            return -log10(1 - b.learning_rate) / log10(2)

        @blk.Constraint()
        def cost_NOAK_eq(b):
            return b.cost_NOAK == pyunits.convert(
                cost_FOAK * ((b.cum_num_units) ** -(b.learning_rate_exponent)),
                to_units=getattr(pyunits, "MUSD_" + CE_index_year),
            )
