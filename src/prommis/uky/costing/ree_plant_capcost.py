#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
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

__author__ = "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, and M. Zamarripa)"
__version__ = "1.0.0"

import textwrap
from sys import stdout

from pyomo.common.dependencies import attempt_import
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.units_container import InconsistentUnitsError, UnitsError
from pyomo.environ import ConcreteModel, Constraint, Expression, Param, Var
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
from idaes.core.util.tables import stream_table_dataframe_to_string

from pandas import DataFrame

from prommis.uky.costing.costing_dictionaries import load_REE_costing_dictionary

_, watertap_costing_available = attempt_import("watertap.costing")
if watertap_costing_available:
    from watertap.costing import WaterTAPCosting

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
    # register_idaes_currency_units()
    custom_REE_plant_currency_units()

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

    # pylint: disable-next=dangerous-default-value
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
        nameplate_capacity=500,
        labor_types=[
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ],
        labor_rate=[27.90, 23.26, 30.29, 24.06, 23.43, 46.82],
        labor_burden=25,
        operators_per_shift=[2, 5, 2, 3, 1, 2],
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
        CE_index_year="2021",
        watertap_blocks=None,
    ):
        """
        This method builds process-wide costing, including fixed and variable
        operating & maintenance costs, costs of production, cost of
        electricity and cost of capture.

        If individual percentages are provided (defaults to values above), this
        method creates constraints for the following plantwide costs ($MM/yr):
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
                purchase cost. To use as the total plant cost, including
                installation, also set the Lang_factor to 1.
            Lang_factor: single multiplicative factor to estimate installation
                costs; defaults to None and method will use percentages. The
                default percentages yield an effective Lang factor of 2.97.
            piping_materials_and_labor_percentage: Piping, materials and labor
                costs as a percentage of the total plant cost. If Lang_factor
                is not None, this value will not be used.
            electrical_materials_and_labor_percentage: Electrical, materials
                and labor costs as a percentage of the total plant cost. If
                Lang_factor is not None, this value will not be used.
            instrumentation_percentage: Instrumentation costs as a percentage
                of the total plant cost. If Lang_factor is not None, this value
                will not be used.
            plants_services_percentage: Plant services costs as a percentage
                of the total plant cost. If Lang_factor is not None, this value
                will not be used.
            process_buildings_percentage: Process buildings costs as a
                percentage of the total plant cost. If Lang_factor is not None,
                this value will not be used.
            auxiliary_buildings_percentage: Auxiliary buildings costs as a
                percentage of the total plant cost. If Lang_factor is not None,
                this value will not be used.
            site_improvements_percentage: Site improvements costs as a
                percentage of the total plant cost. If Lang_factor is not None,
                this value will not be used.
            equipment_installation_percentage: Equipment installation costs as
                a percentage of the total plant cost. If Lang_factor is not
                None, this value will not be used.
            field_expenses_percentage: Field expenses costs as a percentage of
                the total plant cost. If Lang_factor is not None, this value
                will not be used.
            project_management_and_construction_percentage: Project management
                and construction costs as a percentage of the total plant cost.
                If Lang_factor is not None, this value will not be used.
            process_contingency_percentage: Process contingency costs as a
                percentage of the total plant cost. If Lang_factor is not None,
                this value will not be used.
            total_purchase_cost: The BEC in $MM that will be used to determine
                installation and fixed O&M costs. If the value is None, the
                function will try to use the BEC calculated from the individual
                units. This quantity should be a Pyomo Var or Param that will
                contain the BEC value.
            nameplate_capacity: rated plant output in short (US) ton/hr
            labor_type: list of types of operators present in plant; assumed to
                correspond with labor rate and operator per shift lists
            labor_rate: hourly rate of plant operators in project dollar year;
                defined as list corresponding to different operator types
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses; assumed constant for all operator types
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
            transport_cost_per_ton_product: Expression, Var or Param to use for transport costs
                per ton of product (note, this is not part of the TOC)
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
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

        if total_purchase_cost is None:
            self.get_total_BEC(CE_index_year, watertap_blocks)
        else:
            self.total_BEC = Var(
                initialize=total_purchase_cost,
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )
            self.total_BEC.fix()

        # define variables
        if Lang_factor is None:
            # initialize parameters from specified percentages
            self.piping_materials_and_labor_percentage = Param(
                mutable=True,
                initialize=piping_materials_and_labor_percentage / 100,
                doc="Piping, materials and labor",
            )

            self.electrical_materials_and_labor_percentage = Param(
                mutable=True,
                initialize=electrical_materials_and_labor_percentage / 100,
                doc="Electrical, materials and labor",
            )

            self.instrumentation_percentage = Param(
                mutable=True,
                initialize=instrumentation_percentage / 100,
                doc="Instrumentation",
            )

            self.plant_services_percentage = Param(
                mutable=True,
                initialize=plants_services_percentage / 100,
                doc="Plant services",
            )

            self.process_buildings_percentage = Param(
                mutable=True,
                initialize=process_buildings_percentage / 100,
                doc="Process buildings",
            )

            self.auxiliary_buildings_percentage = Param(
                mutable=True,
                initialize=auxiliary_buildings_percentage / 100,
                doc="Auxiliary buildings",
            )

            self.site_improvements_percentage = Param(
                mutable=True,
                initialize=site_improvements_percentage / 100,
                doc="Site improvements",
            )

            self.equipment_installation_percentage = Param(
                mutable=True,
                initialize=equipment_installation_percentage / 100,
                doc="Equipment installation",
            )

            self.field_expenses_percentage = Param(
                mutable=True,
                initialize=field_expenses_percentage / 100,
                doc="Field expenses",
            )

            self.project_management_and_construction_percentage = Param(
                mutable=True,
                initialize=project_management_and_construction_percentage / 100,
                doc="Project management and construction",
            )

            self.process_contingency_percentage = Param(
                mutable=True,
                initialize=process_contingency_percentage / 100,
                doc="Process contingency",
            )

            # ancillary cost variables
            self.ancillary_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.piping_materials_and_labor_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Piping, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.electrical_materials_and_labor_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Electrical, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.instrumentation_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.plant_services_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # buildings cost variables
            self.buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.process_buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Process buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.auxiliary_buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Auxiliary buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.site_improvements_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Site improvements buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # engineering, procurement and construction management cost variables
            self.epcm_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="EPCM cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.equipment_installation_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Equipment installation EPCM cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.field_expenses_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Field expenses EPCM cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.project_management_and_construction_costs = Var(
                initialize=self.total_BEC,
                bounds=(0, 1e4),
                doc="Project management and construction EPCM cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # contingency cost variables - generic to support more contingency cost types in the future
            self.contingency_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Contingency cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.process_contingency_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="Contingency cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
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
            bounds=(0, 1e4),
            doc="Total installation cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        self.total_plant_cost = Var(
            initialize=self.total_BEC,
            bounds=(0, 1e4),
            doc="Total plant cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        # add other plant costs to catch non-equipment capital costs, e.g. reagent fills
        self.other_plant_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Additional plant costs in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )
        self.other_plant_costs.fix(0)

        if Lang_factor is None:
            # rules for calculating Ancillary costs
            def piping_materials_and_labor_cost_rule(self):
                return self.piping_materials_and_labor_costs == (
                    self.total_BEC * self.piping_materials_and_labor_percentage
                )

            self.piping_materials_and_labor_cost_eq = Constraint(
                rule=piping_materials_and_labor_cost_rule
            )

            def electrical_materials_and_labor_cost_rule(self):
                return self.electrical_materials_and_labor_costs == (
                    self.total_BEC * self.electrical_materials_and_labor_percentage
                )

            self.electrical_materials_and_labor_cost_eq = Constraint(
                rule=electrical_materials_and_labor_cost_rule
            )

            def instrumentation_cost_rule(self):
                return self.instrumentation_costs == (
                    self.total_BEC * self.instrumentation_percentage
                )

            self.instrumentation_cost_eq = Constraint(rule=instrumentation_cost_rule)

            def plant_services_cost_rule(self, i):
                return self.plant_services_costs == (
                    self.total_BEC * self.plant_services_percentage
                )

            self.plant_services_cost_eq = Constraint(rule=plant_services_cost_rule)

            def ancillary_cost_rule(self):
                return self.ancillary_costs == (
                    self.piping_materials_and_labor_costs
                    + self.electrical_materials_and_labor_costs
                    + self.instrumentation_costs
                    + self.plant_services_costs
                )

            self.ancillary_cost_eq = Constraint(rule=ancillary_cost_rule)

            # rules for calculating Buildings costs
            def process_buildings_cost_rule(self):
                return self.process_buildings_costs == (
                    self.total_BEC * self.process_buildings_percentage
                )

            self.process_buildings_cost_eq = Constraint(
                rule=process_buildings_cost_rule
            )

            def auxiliary_buildings_cost_rule(self):
                return self.auxiliary_buildings_costs == (
                    self.total_BEC * self.auxiliary_buildings_percentage
                )

            self.auxiliary_buildings_cost_eq = Constraint(
                rule=auxiliary_buildings_cost_rule
            )

            def site_improvements_cost_rule(self):
                return self.site_improvements_costs == (
                    self.total_BEC * self.site_improvements_percentage
                )

            self.site_improvements_cost_eq = Constraint(
                rule=site_improvements_cost_rule
            )

            def buildings_cost_rule(self):
                return self.buildings_costs == (
                    self.process_buildings_costs
                    + self.auxiliary_buildings_costs
                    + self.site_improvements_costs
                )

            self.buildings_cost_eq = Constraint(rule=buildings_cost_rule)

            # rules for calculating Engineering, Procurement and Construction Management costs
            def equipment_installation_cost_rule(self):
                return self.equipment_installation_costs == (
                    self.total_BEC * self.equipment_installation_percentage
                )

            self.equipment_installation_cost_eq = Constraint(
                rule=equipment_installation_cost_rule
            )

            def field_expenses_cost_rule(self):
                return self.field_expenses_costs == (
                    self.total_BEC * self.field_expenses_percentage
                )

            self.field_expenses_cost_eq = Constraint(rule=field_expenses_cost_rule)

            def project_management_and_construction_cost_rule(self):
                return self.project_management_and_construction_costs == (
                    self.total_BEC * self.project_management_and_construction_percentage
                )

            self.project_management_and_construction_cost_eq = Constraint(
                rule=project_management_and_construction_cost_rule
            )

            def epcm_cost_rule(self):
                return self.epcm_costs == (
                    self.equipment_installation_costs
                    + self.field_expenses_costs
                    + self.project_management_and_construction_costs
                )

            self.epcm_cost_eq = Constraint(rule=epcm_cost_rule)

            # rules for calculating Contingency costs
            def process_contingency_cost_rule(self):
                return self.contingency_costs == (
                    self.total_BEC * self.process_contingency_percentage
                )

            self.process_contingency_cost_eq = Constraint(
                rule=process_contingency_cost_rule
            )

            def contingency_cost_rule(self):
                return self.contingency_costs == (self.process_contingency_costs)

            self.contingency_cost_eq = Constraint(rule=contingency_cost_rule)

            def total_installation_cost_rule(self):
                return self.total_installation_cost == (
                    self.ancillary_costs
                    + self.buildings_costs
                    + self.epcm_costs
                    + self.contingency_costs
                )

            self.total_installation_cost_eq = Constraint(
                rule=total_installation_cost_rule
            )
        else:

            def total_installation_cost_rule(self):
                return self.total_installation_cost == self.total_BEC * (
                    self.Lang_factor - 1
                )

            self.total_installation_cost_eq = Constraint(
                rule=total_installation_cost_rule
            )

        # rule for calculating TPC
        def total_plant_cost_rule(self):
            return self.total_plant_cost == (
                self.total_BEC + self.total_installation_cost + self.other_plant_costs
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # define land cost
        if land_cost is not None:
            if type(land_cost) in [Expression, ScalarExpression]:
                if pyunits.get_units(land_cost) == pyunits.dimensionless:
                    self.land_cost = Expression(expr=land_cost.expr * CE_index_units)
                else:
                    self.land_cost = Expression(
                        expr=pyunits.convert(land_cost.expr, to_units=CE_index_units)
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
            ):  # assume it's short ton per hour
                feed_input_rate = feed_input * pyunits.ton / pyunits.hr
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
                if pyunits.get_units(additional_waste_cost) == pyunits.dimensionless:
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
                if pyunits.get_units(additional_waste_cost) == pyunits.dimensionless:
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
            # build cost of recovery (COR)
            if recovery_rate_per_year is not None:
                self.additional_cost_of_recovery = Var(
                    initialize=0,
                    doc="Additional cost to be added to the COR calculations"
                    + " in millions",
                    units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
                )

                if not hasattr(self, "recovery_rate_per_year"):
                    if (
                        pyunits.get_units(recovery_rate_per_year)
                        == pyunits.dimensionless
                    ):  # assume it's in kg/year
                        self.recovery_rate_per_year = Param(
                            initialize=recovery_rate_per_year,
                            mutable=True,
                            units=pyunits.kg / pyunits.year,
                        )
                    else:  # use source units
                        self.recovery_rate_per_year = Param(
                            initialize=recovery_rate_per_year,
                            mutable=True,
                            units=pyunits.get_units(recovery_rate_per_year),
                        )
                    recovery_units_factor = 1

                rec_rate_units = pyunits.get_units(self.recovery_rate_per_year)

                # check that units are compatible
                try:
                    pyunits.convert(
                        self.recovery_rate_per_year * recovery_units_factor,
                        to_units=pyunits.kg / pyunits.year,
                    )
                except InconsistentUnitsError:
                    raise UnitsError(
                        f"The argument recovery_rate_per_year was passed with units of "
                        f"{rec_rate_units} which cannot be converted to units of mass per year. "
                        f"Please ensure that recovery_rate_per_year is passed with rate units "
                        f"of mass per year (mass/a) or dimensionless."
                    )

                # check that units are on an annual basis
                if str(rec_rate_units).split("/")[1] not in ["a", "year"]:
                    raise UnitsError(
                        f"The argument recovery_rate_per_year was passed with units of "
                        f"{rec_rate_units} and must be on an anuual basis. Please "
                        f"ensure that recovery_rate_per_year is passed with rate units "
                        f"of mass per year (mass/a) or dimensionless."
                    )

                self.cost_of_recovery = Expression(
                    expr=(
                        pyunits.convert(
                            (
                                self.annualized_cost / pyunits.year
                                + self.total_fixed_OM_cost / pyunits.year
                                + self.total_variable_OM_cost[0]
                            )
                            / (self.recovery_rate_per_year * recovery_units_factor),
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
                        # no units, assume $/ton
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

        if hasattr(self, "total_plant_cost"):
            var_dict["Total Plant Cost [$MM]"] = value(self.total_plant_cost)

        if hasattr(self, "total_BEC"):
            var_dict["Total Bare Erected Cost [$MM]"] = value(self.total_BEC)

        if hasattr(self, "total_installation_cost"):
            var_dict["Total Installation Cost [$MM]"] = value(
                self.total_installation_cost
            )

        if hasattr(self, "other_plant_costs"):
            var_dict["Total Other Plant Costs [$MM/year]"] = value(
                self.other_plant_costs
            )

        if hasattr(self, "ancillary_costs"):
            var_dict["Summation of Ancillary Installation Costs [$MM]"] = value(
                self.ancillary_costs
            )

        if hasattr(self, "piping_materials_and_labor_costs"):
            var_dict[
                "Total Ancillary Piping, Materials and Labor Installation Cost [$MM]"
            ] = value(self.piping_materials_and_labor_costs)

        if hasattr(self, "electrical_materials_and_labor_costs"):
            var_dict[
                "Total Ancillary Electrical, Materials and Labor Installation Cost [$MM]"
            ] = value(self.electrical_materials_and_labor_costs)

        if hasattr(self, "instrumentation_costs"):
            var_dict["Total Ancillary Instrumentation Installation Cost [$MM]"] = value(
                self.instrumentation_costs
            )

        if hasattr(self, "plant_services_costs"):
            var_dict["Total Ancillary Plant Services Installation Cost [$MM]"] = value(
                self.plant_services_costs
            )

        if hasattr(self, "buildings_costs"):
            var_dict["Summation of Buildings Installation Costs [$MM]"] = value(
                self.buildings_costs
            )

        if hasattr(self, "process_buildings_costs"):
            var_dict["Total Process Buildings Installation Cost [$MM]"] = value(
                self.process_buildings_costs
            )

        if hasattr(self, "auxiliary_buildings_costs"):
            var_dict["Total Auxiliary Buildings Installation Cost [$MM]"] = value(
                self.auxiliary_buildings_costs
            )

        if hasattr(self, "site_improvements_costs"):
            var_dict["Total Site Improvements Buildings Installation Cost [$MM]"] = (
                value(self.site_improvements_costs)
            )

        if hasattr(self, "epcm_costs"):
            var_dict["Summation of EPCM Installation Costs [$MM]"] = value(
                self.epcm_costs
            )

        if hasattr(self, "equipment_installation_costs"):
            var_dict["Total Equipment Installation EPCM Installation Cost [$MM]"] = (
                value(self.equipment_installation_costs)
            )

        if hasattr(self, "field_expenses_costs"):
            var_dict["Total Field Expenses EPCM Cost [$MM]"] = value(
                self.field_expenses_costs
            )

        if hasattr(self, "project_management_and_construction_costs"):
            var_dict[
                "Total Project Management and Construction EPCM Installation Cost [$MM]"
            ] = value(self.project_management_and_construction_costs)

        if hasattr(self, "process_contingency_costs"):
            var_dict["Total Process Contingency Installation Cost [$MM]"] = value(
                self.process_contingency_costs
            )

        if hasattr(self, "contingency_costs"):
            var_dict["Summation of Contingency Installation Costs [$MM]"] = value(
                self.contingency_costs
            )

        if hasattr(self, "total_fixed_OM_cost"):
            var_dict["Total Fixed Operating & Maintenance Cost [$MM/year]"] = value(
                self.total_fixed_OM_cost
            )

        if hasattr(self, "annual_operating_labor_cost"):
            var_dict["Total Annual Operating Labor Cost [$MM/year]"] = value(
                self.annual_operating_labor_cost
            )

            var_dict["Total Annual Technical Labor Cost [$MM/year]"] = value(
                self.annual_technical_labor_cost
            )

            var_dict["Summation of Annual Labor Costs [$MM/year]"] = value(
                self.annual_labor_cost
            )

        if hasattr(self, "maintenance_and_material_cost"):
            var_dict["Total Maintenance and Material Cost [$MM/year]"] = value(
                self.maintenance_and_material_cost
            )

        if hasattr(self, "quality_assurance_and_control_cost"):
            var_dict["Total Quality Assurance and Control Cost [$MM/year]"] = value(
                self.quality_assurance_and_control_cost
            )

        general_sales_and_admin = 0

        if hasattr(self, "sales_patenting_and_research_cost"):
            var_dict["Total Sales, Patenting and Research Cost [$MM/year]"] = value(
                self.sales_patenting_and_research_cost
            )
            general_sales_and_admin += value(self.sales_patenting_and_research_cost)

        if hasattr(self, "admin_and_support_labor_cost"):
            var_dict["Total Admin Support and Labor Cost [$MM/year]"] = value(
                self.admin_and_support_labor_cost
            )
            general_sales_and_admin += value(self.admin_and_support_labor_cost)

        if hasattr(self, "property_taxes_and_insurance_cost"):
            var_dict["Total Property Taxes and Insurance Cost [$MM/year]"] = value(
                self.property_taxes_and_insurance_cost
            )
            general_sales_and_admin += value(self.property_taxes_and_insurance_cost)

        var_dict["Summation of Sales, Admin and Insurance Cost [$MM/year]"] = value(
            general_sales_and_admin
        )

        if hasattr(self, "other_fixed_costs"):
            var_dict["Total Other Fixed Costs [$MM/year]"] = value(
                self.other_fixed_costs
            )

        if hasattr(self, "variable_operating_costs"):
            var_dict["Total Variable Power Cost [$MM/year]"] = value(
                self.variable_operating_costs[0, "power"]
            )

            if hasattr(self, "additional_waste_cost"):
                var_dict["Total Variable Waste Cost [$MM/year]"] = value(
                    sum(
                        self.variable_operating_costs[0, waste]
                        for waste in self.waste_list
                    )
                    + self.additional_waste_cost
                )

            if hasattr(self, "additional_chemicals_cost"):
                var_dict["Total Variable Chemicals Cost [$MM/year]"] = value(
                    sum(
                        self.variable_operating_costs[0, chemical]
                        for chemical in self.chemicals_list
                    )
                    + self.additional_chemicals_cost
                )

            var_dict["General Plant Overhead Cost [$MM/year]"] = value(
                self.plant_overhead_cost[0]
            )

            var_dict[
                "Total Plant Overhead Cost, Including Maintenance & Quality Assurance [$MM/year]"
            ] = value(
                self.plant_overhead_cost[0]
                + self.maintenance_and_material_cost
                + self.quality_assurance_and_control_cost
            )

        if hasattr(self, "total_variable_OM_cost"):
            var_dict["Total Variable Operating & Maintenance Cost [$MM/year]"] = value(
                self.total_variable_OM_cost[0]
            )

        if hasattr(self, "land_cost"):
            var_dict["Total Land Cost [$MM/year]"] = value(self.land_cost)

        if hasattr(self, "transport_cost"):
            var_dict["Total Transport Cost [$MM/year]"] = value(self.transport_cost)

        if hasattr(self, "total_sales_revenue"):
            var_dict["Total Sales Revenue Cost [$MM/year]"] = value(
                self.total_sales_revenue
            )

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


        Cost is in M$
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
            bounds=(0, 1e4),
            doc="Scaled bare erected cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        # rule for scaling BEC
        def bare_erected_cost_rule(costing, i):
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
            elif ref_cost_units[0] == "K":  # thousands of $
                ref_cost_units = getattr(pyunits, "kUSD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "M":  # millions of $
                ref_cost_units = getattr(pyunits, "MUSD_" + ref_cost_units[1])

            # determine reference parameter scaler based on train scaling
            if scale_down_parallel_equip:
                scaler = n_equip
            else:
                scaler = 1

            if isinstance(process_params[i], list):
                if len(process_params[i]) > 1:
                    return costing.bare_erected_cost[i] == (
                        n_equip
                        * pyunits.convert(
                            costing.ref_cost[i] * ref_cost_units, CE_index_units
                        )
                        * sum(
                            (
                                pyunits.convert(scaled_param[j], ref_units)
                                / (scaler * costing.ref_param[i, p] * ref_units)
                            )
                            ** costing.exp[i]
                            for j, p in enumerate(process_params[i])
                        )
                    )
            elif isinstance(process_params[i], str):
                return costing.bare_erected_cost[i] == (
                    n_equip
                    * pyunits.convert(
                        costing.ref_cost[i] * ref_cost_units, CE_index_units
                    )
                    * (
                        pyunits.convert(scaled_param, ref_units)
                        / (scaler * costing.ref_param[i] * ref_units)
                    )
                    ** costing.exp[i]
                )

        blk.bare_erected_cost_eq = Constraint(
            cost_accounts, rule=bare_erected_cost_rule
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
        labor_types=[
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ],
        labor_rate=[27.90, 23.26, 30.29, 24.06, 23.43, 46.82],
        labor_burden=25,
        operators_per_shift=[2, 5, 2, 3, 1, 2],
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
                labor expenses; assumed constant for all operator types
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
        # the currency units are millions of USD, so all prices need a 1e-6 multiplier to get USD
        default_sale_prices = {
            # pure elements
            "Sc": 6442 * 1e-6 * CE_index_units / pyunits.kg,
            "Y": 8 * 1e-6 * CE_index_units / pyunits.kg,
            "La": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "Ce": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "Pr": 63 * 1e-6 * CE_index_units / pyunits.kg,
            "Nd": 49 * 1e-6 * CE_index_units / pyunits.kg,
            "Sm": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "Eu": 174 * 1e-6 * CE_index_units / pyunits.kg,
            "Gd": 37 * 1e-6 * CE_index_units / pyunits.kg,
            "Tb": 471 * 1e-6 * CE_index_units / pyunits.kg,
            "Dy": 264 * 1e-6 * CE_index_units / pyunits.kg,
            "Ho": 61 * 1e-6 * CE_index_units / pyunits.kg,
            "Er": 39 * 1e-6 * CE_index_units / pyunits.kg,
            "Tm": 0 * 1e-6 * CE_index_units / pyunits.kg,  # price not available
            "Yb": 33 * 1e-6 * CE_index_units / pyunits.kg,
            "Lu": 906 * 1e-6 * CE_index_units / pyunits.kg,
            # oxides
            "Sc2O3": 4200 * 1e-6 * CE_index_units / pyunits.kg,
            "Y2O3": 6 * 1e-6 * CE_index_units / pyunits.kg,
            "La2O3": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "CeO2": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "Pr6O11": 52 * 1e-6 * CE_index_units / pyunits.kg,
            "Nd2O3": 42 * 1e-6 * CE_index_units / pyunits.kg,
            "Sm2O3": 2 * 1e-6 * CE_index_units / pyunits.kg,
            "Eu2O3": 150 * 1e-6 * CE_index_units / pyunits.kg,
            "Gd2O3": 32 * 1e-6 * CE_index_units / pyunits.kg,
            "Tb4O7": 400 * 1e-6 * CE_index_units / pyunits.kg,
            "Dy2O3": 230 * 1e-6 * CE_index_units / pyunits.kg,
            "Ho2O3": 53 * 1e-6 * CE_index_units / pyunits.kg,
            "Er2O3": 34 * 1e-6 * CE_index_units / pyunits.kg,
            "Tm2O3": 0 * 1e-6 * CE_index_units / pyunits.kg,  # price not available
            "Yb2O3": 29 * 1e-6 * CE_index_units / pyunits.kg,
            "Lu2O3": 797 * 1e-6 * CE_index_units / pyunits.kg,
        }

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
            initialize=operating_days_per_year, mutable=True, units=pyunits.day
        )
        b.mixed_product_sale_price_realization_factor = Param(
            initialize=mixed_product_sale_price_realization_factor,
            mutable=True,
            units=pyunits.dimensionless,
        )

        # make vars
        b.annual_operating_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Annual operating labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.annual_technical_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Annual technical labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.annual_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Annual labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.maintenance_and_material_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Maintenance and material cost in $MM/yr",
            units=CE_index_units,
        )
        b.quality_assurance_and_control_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Quality assurance and control cost in $MM/yr",
            units=CE_index_units,
        )
        b.sales_patenting_and_research_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Sales, patenting and research cost in $MM/yr",
            units=CE_index_units,
        )
        b.admin_and_support_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Admin and support labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.property_taxes_and_insurance_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="Property taxes and insurance cost in $MM/yr",
            units=CE_index_units,
        )
        b.total_fixed_OM_cost = Var(
            initialize=4,
            bounds=(0, 1e4),
            doc="Total fixed O&M costs in $MM/yr",
            units=CE_index_units,
        )
        b.total_sales_revenue = Var(
            initialize=4,
            bounds=(0, 1e4),
            doc="Total sales revenue in $MM/yr",
            units=CE_index_units,
        )

        # variable for user to assign other fixed costs to,
        # fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Other fixed costs in $MM/yr",
            units=CE_index_units,
        )
        b.other_fixed_costs.fix(0)

        # variable for user to assign watertap fixed costs to,
        # fixed to 0 by default
        b.watertap_fixed_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Watertap fixed costs in $MM/yr",
            units=CE_index_units,
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
        def annual_operating_labor_cost_rule(c):
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
                CE_index_units,
            )

        @b.Constraint()
        def annual_technical_labor_cost_rule(c):
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
                CE_index_units,
            )

        @b.Constraint()
        def annual_labor_cost_rule(c):
            return c.annual_labor_cost == pyunits.convert(
                (c.annual_operating_labor_cost + c.annual_technical_labor_cost),
                CE_index_units,
            )

        # maintenance cost is 2% of TPC
        @b.Constraint()
        def maintenance_and_material_cost_rule(c):
            return c.maintenance_and_material_cost == 0.02 * TPC

        # quality assurance cost is 10% of operating labor
        @b.Constraint()
        def quality_assurance_and_control_cost_rule(c):
            return c.quality_assurance_and_control_cost == 0.10 * pyunits.convert(
                (c.annual_operating_labor_cost),
                CE_index_units,
            )

        # sales cost is 0.5% of total revenue
        @b.Constraint()
        def sales_patenting_and_research_cost_rule(c):
            return c.sales_patenting_and_research_cost == 0.005 * pyunits.convert(
                (c.total_sales_revenue),
                CE_index_units,
            )

        # admin cost is 20% of direct labor
        @b.Constraint()
        def admin_and_support_labor_cost_rule(c):
            return c.admin_and_support_labor_cost == 0.20 * pyunits.convert(
                (c.annual_operating_labor_cost),
                CE_index_units,
            )

        # taxes are 1% of TPC
        @b.Constraint()
        def taxes_and_insurance_cost_rule(c):
            return c.property_taxes_and_insurance_cost == 0.01 * TPC

        # sum of fixed O&M costs

        # sum of fixed operating costs of membrane units
        @b.Constraint()
        def sum_watertap_fixed_cost(c):
            if not hasattr(c, "watertap_fixed_costs_list"):
                return c.watertap_fixed_costs == 0
            else:
                return c.watertap_fixed_costs == sum(b.watertap_fixed_costs_list)

        @b.Constraint()
        def total_fixed_OM_cost_rule(c):
            return c.total_fixed_OM_cost == (
                c.annual_labor_cost
                + c.maintenance_and_material_cost
                + c.quality_assurance_and_control_cost
                + c.admin_and_support_labor_cost
                + c.sales_patenting_and_research_cost
                + c.property_taxes_and_insurance_cost
                + c.other_fixed_costs
                + c.watertap_fixed_costs
            )

        @b.Constraint()
        def total_sales_revenue_rule(c):
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
                CE_index_units,
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
        if feed_input_rate is None:
            raise AttributeError(
                "No feed_input rate variable passed to main costing block."
            )
        else:
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
        default_prices = {
            "power": 0.07 * 1e-6 * CE_index_units / pyunits.kWh,
            "water": 1.90e-3 * 1e-6 * CE_index_units / pyunits.gallon,
            "diesel": 2 * 1e-6 * CE_index_units / pyunits.gal,
            "bioleaching_solution": 0.008 * 1e-6 * CE_index_units / pyunits.L,
            "H2SO4": 200 * 1e-6 * CE_index_units / pyunits.tonne,
            "natural_gas": 5.79 * 1e-3 * 1e-6 * CE_index_units / pyunits.ft**3,
            "polymer": 33.61 * 1e-6 * CE_index_units / pyunits.kg,
            "NAOH": 350.00 * 1e-6 * CE_index_units / pyunits.tonne,
            "CACO3": 80.00 * 1e-6 * CE_index_units / pyunits.tonne,
            "coal_calcite": 0.50 * 1e-6 * CE_index_units / pyunits.tonne,
            "HCL": 250.00 * 1e-6 * CE_index_units / pyunits.tonne,
            "oxalic_acid": 1.00 * 1e-6 * CE_index_units / pyunits.kg,
            "ascorbic_acid": 2.00 * 1e-6 * CE_index_units / pyunits.kg,
            "kerosene": 400.00 * 1e-6 * CE_index_units / pyunits.tonne,
            "D2EHPA": 15.00 * 1e-6 * CE_index_units / pyunits.kg,
            "NA2S": 360.00 * 1e-6 * CE_index_units / pyunits.tonne,
            "nonhazardous_solid_waste": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
            "nonhazardous_precipitate_waste": 5.00
            * 1e-6
            * CE_index_units
            / pyunits.ton,
            "dust_and_volatiles": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
        }

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
            doc="Variable operating costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        b.other_variable_costs = Var(
            b.parent_block().time,
            initialize=0,
            bounds=(0, 1e4),
            doc="A variable to include non-standard O&M costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        # assume the user is not using this
        b.other_variable_costs.fix(0)

        b.total_variable_OM_cost = Var(
            b.parent_block().time,
            initialize=4e-6,
            doc="Total variable operating and maintenance costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        @b.Constraint(b.parent_block().time, resources)
        def variable_cost_rule(c, t, r):
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
                    * c.operating_days_per_year
                    * pyunits.year**-1,
                    to_units=CE_index_units / pyunits.year,
                )
            )

        if hasattr(b, "total_fixed_OM_cost"):
            # define overhead cost
            # plant overhead, 20% of direct costs - fixed OM, power, water, lease/land, chemicals, waste
            b.plant_overhead_cost = Var(
                b.parent_block().time,
                initialize=0,
                doc="Plant overhead costs in $MM/year",
                units=CE_index_units / pyunits.year,
            )

        if (0, "power") in b.variable_operating_costs.id_index_map().values():

            @b.Constraint(b.parent_block().time)
            def plant_overhead_cost_rule(c, t):
                return c.plant_overhead_cost[t] == 0.20 * (
                    c.total_fixed_OM_cost / pyunits.year
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
                )

        else:

            @b.Constraint(b.parent_block().time)
            def plant_overhead_cost_rule(c, t):
                return c.plant_overhead_cost[t] == 0.20 * (
                    c.total_fixed_OM_cost / pyunits.year
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
                )

        @b.Constraint(b.parent_block().time)
        def total_variable_cost_rule(c, t):
            return (
                c.total_variable_OM_cost[t]
                == sum(c.variable_operating_costs[t, r] for r in resources)
                + c.other_variable_costs[t]
                + c.plant_overhead_cost[t]
                + c.land_cost / pyunits.year
                + c.additional_chemicals_cost / pyunits.year
                + c.additional_waste_cost / pyunits.year
            )

    def initialize_fixed_OM_costs(b):
        # b is the flowsheet-level costing block
        if hasattr(b, "total_fixed_OM_cost"):
            calculate_variable_from_constraint(
                b.annual_operating_labor_cost, b.annual_labor_cost_rule
            )

            calculate_variable_from_constraint(
                b.maintenance_and_material_cost, b.maintenance_and_material_cost_rule
            )

            calculate_variable_from_constraint(
                b.quality_assurance_and_control_cost,
                b.quality_assurance_and_control_cost_rule,
            )

            calculate_variable_from_constraint(
                b.sales_patenting_and_research_cost,
                b.sales_patenting_and_research_cost_rule,
            )

            calculate_variable_from_constraint(
                b.admin_and_support_labor_cost,
                b.admin_and_support_labor_cost_rule,
            )

            calculate_variable_from_constraint(
                b.property_taxes_and_insurance_cost,
                b.taxes_and_insurance_cost_rule,
            )

            calculate_variable_from_constraint(
                b.total_fixed_OM_cost, b.total_fixed_OM_cost_rule
            )

    def initialize_variable_OM_costs(b):
        # b is the flowsheet-level costing block
        # initialization for power generation costs
        if hasattr(b, "variable_operating_costs"):
            for i in b.variable_operating_costs.keys():
                if hasattr(b, "variable_cost_rule"):
                    calculate_variable_from_constraint(
                        b.variable_operating_costs[i],
                        b.variable_cost_rule[i],
                    )

            for i in b.total_variable_OM_cost.keys():
                calculate_variable_from_constraint(
                    b.total_variable_OM_cost[i],
                    b.total_variable_cost_rule[i],
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
                    "%s: $%.2f Million"
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
                    "%s: $%.5f Million"
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

        BEC_list = []
        b.watertap_fixed_costs_list = []

        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "bare_erected_cost"):
                for key in o.bare_erected_cost.keys():
                    BEC_list.append(o.bare_erected_cost[key])

        if watertap_blocks is not None:
            for w in watertap_blocks:
                m = ConcreteModel()
                m.fs = FlowsheetBlock(dynamic=False)
                m.fs.costing = WaterTAPCosting()
                w.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
                BEC_list.append(
                    pyunits.convert(w.costing.capital_cost, to_units=CE_index_units)
                )
                b.watertap_fixed_costs_list.append(
                    pyunits.convert(
                        w.costing.fixed_operating_cost * pyunits.year,
                        to_units=CE_index_units,
                    )
                )

        b.total_BEC = Var(
            initialize=100,
            bounds=(0, 1e4),
            doc="Total TPC in $MM",
            # assume that total_plant_cost is in millions of
            # USD_year, where year is the CE_index_year users set
            units=CE_index_units,
        )

        @b.Constraint()
        def total_BEC_eq(c):
            return c.total_BEC == sum(BEC_list)

    def display_flowsheet_cost(b):
        # This method accepts a flowsheet-level costing block
        print("\n")
        print("Total bare erected cost: $%.3f Million" % value(b.total_BEC))
        if hasattr(b, "total_overnight_capital"):
            print(
                "Total overnight (installed) equipment cost: $%.3f Million"
                % value(b.total_overnight_capital)
            )
        if hasattr(b, "annualized_cost"):
            print(
                "Total annualized capital cost: $%.3f Million"
                % value(b.annualized_cost)
            )
        print()
        if hasattr(b, "total_fixed_OM_cost"):
            print(
                "Total annual fixed O&M cost: $%.3f Million"
                % value(b.total_fixed_OM_cost)
            )
        if hasattr(b, "total_variable_OM_cost"):
            print(
                "Total annual variable O&M cost: $%.3f Million"
                % value(b.total_variable_OM_cost[0])
            )
        if hasattr(b, "total_fixed_OM_cost") and hasattr(b, "total_variable_OM_cost"):
            print(
                "Total annual O&M cost: $%.3f Million"
                % value(b.total_fixed_OM_cost + b.total_variable_OM_cost[0])
            )
            if hasattr(b, "feed_input_rate"):
                print(
                    "Total annual O&M cost: $%.3f per ton feed processed"
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
                    "Total annual O&M cost: $%.3f per kg REE recovered"
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
                "Total annualized plant cost: $%.3f Million"
                % value(
                    b.annualized_cost
                    + b.total_fixed_OM_cost
                    + b.total_variable_OM_cost[0]
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
                "Cost of recovery: $%.3f per kg REE recovered"
                % value(b.cost_of_recovery)
            )
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
                bounds=(0, 1e9),
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
                bounds=(0, 100),
                doc="Grade percentage of site",
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
                bounds=(0, 100),
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
                bounds=(0, 100),
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
                            / pyunits.dimensionless
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
                            / pyunits.dimensionless
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

        _log.info("\nPrinting calculated costing bounds for processes:")
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
