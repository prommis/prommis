#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
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
    "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, and M. Zamarripa)"
)
__version__ = "1.0.0"

from sys import stdout
import textwrap

from pandas import DataFrame

from pyomo.environ import (
    Param,
    Var,
    Constraint,
    Expression,
    value,
    units as pyunits,
)
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.units_container import InconsistentUnitsError, UnitsError
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from idaes.core import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)
from costing_dictionaries import (
    load_REE_costing_dictionary,
)

from idaes.core.util.tables import stream_table_dataframe_to_string
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class

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
        "USD_2008_Nov" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_2019_Sep" in pyunits._pint_registry  # pylint: disable=protected-access
    ):
        # Assume that custom REE plant units have already been registered
        # Log a message and end
        _log.debug(
            "Custom REE plant currency units (USD_2008_Nov, USD_2019_Sep) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                "USD_2008_Nov = 500/566.2 * USD_CE500",
                "USD_2019_Sep = 500/599.3 * USD_CE500",
                # from https://toweringskills.com/financial-analysis/cost-indices/ as of 9/26/2023
                "USD_2022 = 500/816.0 * USD_CE500",
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

    def build_process_costs(
        self,
        # arguments related to installation costs
        total_purchase_cost=None,
        Lang_factor=None, # default percentages are effective Lang_factor of 2.97
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
        capacity_factor=0.92,
        labor_types = ["skilled", "unskilled", "supervisor", "maintenance", "technician", "engineer"],
        labor_rate=[27.90, 23.26, 30.29, 24.06, 23.43, 46.82],
        labor_burden=25,
        operators_per_shift=[2, 5, 2, 3, 1, 2],
        hours_per_shift=8,
        shifts_per_day=3,
        operating_days_per_year=336,
        # arguments related to total owners costs
        land_cost=None,
        resources=None,
        rates=None,
        prices=None,
        fixed_OM=True,
        variable_OM=False,
        fuel=None,
        chemicals=None,
        chemicals_inventory=None,
        waste=None,
        transport_cost=None,
        tonne_REE_capture=None,
        CE_index_year="2021",
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
            nameplate_capacity: rated plant output in tonne/hr
            capacity_factor: multiplicative factor for normal operating capacity
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
            land_cost: Expression, Var or Param to calculate land costs
            resources: list setting resources to cost
            rates: list setting flow rates of resources
            prices: list setting prices of resources
            fixed_OM: True/False flag for calculating fixed O&M costs
            variable_OM: True/False flag for calculating variable O&M costs
            fuel: string setting fuel type for fuel costs
            chemicals: string setting chemicals type for chemicals costs
            chemicals_inventory: string setting chemicals type for inventory costs
            waste: string setting waste type for waste costs
            tonne_REE_capture: Var or value to use for tonnes of REE capture
                in one year, in units of tonnes (not tonnes/year)
            transport_cost: Expression, Var or Param to use for transport costs
                per ton of REE captured (note, this is not part of the TOC)
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
        """

        # define costing library
        self.library = "REE"

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        if total_purchase_cost is None:
            self.get_total_BEC(CE_index_year)
        else:
            self.total_BEC = Var(
                initialize=total_purchase_cost,
                units=getattr(pyunits, "MUSD_" + CE_index_year)
                )

        # define variables
        if Lang_factor is None:

            # initialize parameters from specified percentages
            self.piping_materials_and_labor_percentage = Param(
                mutable=True,
                initialize=piping_materials_and_labor_percentage/100,
                doc="piping, materials and labor",
            )
    
            self.electrical_materials_and_labor_percentage = Param(
                mutable=True,
                initialize=electrical_materials_and_labor_percentage/100,
                doc="electrical, materials and labor",
            )
    
            self.instrumentation_percentage = Param(
                mutable=True,
                initialize=instrumentation_percentage/100,
                doc="instrumentation",
            )
    
            self.plant_services_percentage = Param(
                mutable=True,
                initialize=plants_services_percentage/100,
                doc="plant services",
            )
    
            self.process_buildings_percentage = Param(
                mutable=True,
                initialize=process_buildings_percentage/100,
                doc="process buildings",
            )
    
            self.auxiliary_buildings_percentage = Param(
                mutable=True,
                initialize=auxiliary_buildings_percentage/100,
                doc="auxiliary buildings",
            )
    
            self.site_improvements_percentage = Param(
                mutable=True,
                initialize=site_improvements_percentage/100,
                doc="site improvements",
            )
    
            self.equipment_installation_percentage = Param(
                mutable=True,
                initialize=equipment_installation_percentage/100,
                doc="equipment installation",
            )
    
            self.field_expenses_percentage = Param(
                mutable=True,
                initialize=field_expenses_percentage/100,
                doc="field expenses",
            )
    
            self.project_management_and_construction_percentage = Param(
                mutable=True,
                initialize=project_management_and_construction_percentage/100,
                doc="project management and construction",
            )
    
            self.process_contingency_percentage = Param(
                mutable=True,
                initialize=process_contingency_percentage/100,
                doc="process contingency",
            )

            # ancillary cost variables
            self.ancillary_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.piping_materials_and_labor_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="piping, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.electrical_materials_and_labor_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="electrical, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.instrumentation_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.plant_services_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # buildings cost variables
            self.buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.process_buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="process buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.auxiliary_buildings_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="auxiliary buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.site_improvements_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="site improvements buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # engineering, procurement and construction management cost variables
            self.epcm_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.equipment_installation_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="equipment installation epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.field_expenses_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="field expenses epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.project_management_and_construction_costs = Var(
                initialize=self.total_BEC,
                bounds=(0, 1e4),
                doc="project management and construction epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # contingency cost variables - generic to support more contingency cost types in the future
            self.contingency_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="contingency cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            self.process_contingency_costs = Var(
                initialize=value(self.total_BEC),
                bounds=(0, 1e4),
                doc="contingency cost in $MM",
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
            doc="total installation cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        self.total_plant_cost = Var(
            initialize=self.total_BEC,
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        if Lang_factor is None:

            # rules for calculating Ancillary costs
            def piping_materials_and_labor_cost_rule(self):
                return self.piping_materials_and_labor_costs == (
                    self.total_BEC
                    * self.piping_materials_and_labor_percentage
                )

            self.piping_materials_and_labor_cost_eq = Constraint(
                rule=piping_materials_and_labor_cost_rule
            )

            def electrical_materials_and_labor_cost_rule(self):
                return self.electrical_materials_and_labor_costs == (
                    self.total_BEC
                    * self.electrical_materials_and_labor_percentage
                )

            self.electrical_materials_and_labor_cost_eq = Constraint(
                rule=electrical_materials_and_labor_cost_rule
            )

            def instrumentation_cost_rule(self):
                return self.instrumentation_costs == (
                    self.total_BEC * self.instrumentation_percentage
                )

            self.instrumentation_cost_eq = Constraint(
                rule=instrumentation_cost_rule
            )

            def plant_services_cost_rule(self, i):
                return self.plant_services_costs == (
                    self.total_BEC * self.plant_services_percentage
                )

            self.plant_services_cost_eq = Constraint(
                rule=plant_services_cost_rule
            )

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

            self.field_expenses_cost_eq = Constraint(
                rule=field_expenses_cost_rule
            )

            def project_management_and_construction_cost_rule(self):
                return self.project_management_and_construction_costs == (
                    self.total_BEC
                    * self.project_management_and_construction_percentage
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

            self.contingency_cost_eq = Constraint(
                rule=contingency_cost_rule
            )

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
                self.total_BEC + self.total_installation_cost
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # build operating & maintenance costs
        if fixed_OM:
            self.get_fixed_OM_costs(
                labor_rate=labor_rate,
                labor_burden=labor_burden,
                operators_per_shift=operators_per_shift,
                hours_per_shift=hours_per_shift,
                shifts_per_day=shifts_per_day,
                operating_days_per_year=operating_days_per_year,
                CE_index_year=CE_index_year,
            )

        # build system costs (owner's, total overnight costs, annualized costs,
        # and cost of recovery)

        # items retrieved from Battelle (2020) (https://www.osti.gov/servlets/purl/1631038)
        # are noted with [1]

        # items retrieved from UKy (2021) (https://www.osti.gov/servlets/purl/1798663)
        # are noted with [2]

        # items retrieved from NETL (2021)
        # (https://www.netl.doe.gov/projects/files/QGESSCostEstMethodforNETLAssessmentsofPowerPlantPerformance_022621.pdf)
        # are noted with [3]

        if fixed_OM and variable_OM:
            # total overnight cost requires fixed and variable OM costs
            self.pct_indirect_capital = Param(
                initialize=17 / 100, doc="Fixed percentage for indirect capital costs"
            )  # engineering & home office fees 7%, general facilities 10% [1]
               # project contigency 10%, process contingency 10% not included
               # because contingency cost is already part of installation cost
            self.pct_royalties = Param(
                initialize=0.5 / 100, doc="Fixed percentage for royalties cost"
            )  # [1]
            self.pct_construction_allowance = Param(
                initialize=(0.5 / 100) * (1 + self.pct_indirect_capital),
                doc="Fixed percentage for construction allowance cost"
            )  # [1]
            self.pct_inventory_capital = Param(
                initialize=(0.5 / 100) * (1 + self.pct_indirect_capital),
                doc="Fixed percentage for inventory capital cost"
            )  # [1]
            self.one_month_OM = Expression(
                expr=(
                    self.total_fixed_OM_cost
                    + self.total_variable_OM_cost[0] * self.capacity_factor * pyunits.year
                )
                / 12 # one month fixed and variable OM [1]
            )
            non_fuel_resources = resources  # duplicate resources list
            if fuel is not None:
                self.fuel_cost_OC = Expression(
                    expr=self.variable_operating_costs[0, fuel] / 12 * 2.25,
                    doc="Owner's costs - 2.25 months of fuel costs",  # [3]
                )
                non_fuel_resources.remove(fuel)  # remove fuel from the list

            if waste is not None:
                self.waste_costs_OC = Expression(
                    expr=(sum(self.variable_operating_costs[0, i] for i in waste) / 12)
                )  # [3]
            if chemicals is not None:
                self.chemical_costs_OC = Expression(
                    expr=(
                        sum(self.variable_operating_costs[0, i] for i in chemicals)
                        / 2  # six months of chemicals
                    )  # [3]
                )

            if chemicals_inventory is not None:
                self.chemical_inventory_costs_OC = Expression(
                    expr=(
                        (
                            sum(
                                self.variable_operating_costs[0, i]
                                for i in chemicals_inventory
                            )
                            / 6
                        )
                        * pyunits.year  # two months of chemicals inventory
                        + 0.005 * self.total_plant_cost
                    )
                )  # [3]

            self.non_fuel_and_waste_OC = Expression(
                expr=(
                    sum(self.variable_operating_costs[0, i] for i in non_fuel_resources)
                    / 12
                )
            )  # [3]

            if land_cost is not None:
                if isinstance(land_cost, (Expression, ScalarExpression)) or (
                    isinstance(land_cost, (Param, Var)) and land_cost.get_units is None
                ):
                    self.land_cost = land_cost * CE_index_units
                else:
                    self.land_cost = land_cost

            self.total_overnight_capital = Expression(
                expr=self.total_plant_cost
                # pre production costs
                + self.one_month_OM
                + (
                    self.chemical_inventory_costs_OC
                    if chemicals_inventory is not None
                    else 0 * CE_index_units
                )  # Initial Cost for Catalyst and Chemicals Inventory
                + (
                    self.maintenance_and_material_cost / 12 / pyunits.year # 1 month materials [3]
                    + self.non_fuel_and_waste_OC  # 1 month nonfuel consumables
                    + (
                        self.waste_costs_OC if waste is not None else 0 * CE_index_units
                    )  # 1 month waste
                    # inventory capital costs
                    + (
                        self.fuel_cost_OC if fuel is not None else 0 * CE_index_units
                    )  # 60 day fuel supply
                    # Other costs
                    + (
                        self.chemical_costs_OC
                        if chemicals is not None
                        else 0 * CE_index_units
                    )  # Initial Cost for Catalyst and Chemicals
                )
                * 1
                * pyunits.year  # variable costs for 1 year
                + (self.land_cost if land_cost is not None else 0 * CE_index_units)
                + self.total_plant_cost
                * (self.pct_indirect_capital
                   + self.pct_royalties
                   + self.pct_construction_allowance
                   + self.pct_inventory_capital)
            )

            self.tasc_toc_factor = Param(
                initialize=1.144,  # [2]
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
                initialize=0.1002,  # [2]
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

            self.additional_cost_of_recovery = Var(
                initialize=0,
                doc="additional cost to be added to the COR calculations"
                + " in millions",
                units=CE_index_units / pyunits.kg,
            )

            # build cost of recovery (COR)
            if tonne_REE_capture is not None:
                if not hasattr(self, "tonne_REE_capture"):
                    self.tonne_REE_capture = Param(
                        initialize=tonne_REE_capture, mutable=True, units=pyunits.tonne
                    )
                self.cost_of_recovery = Expression(
                    expr=(
                        (
                            self.annualized_cost / pyunits.year
                            + self.total_fixed_OM_cost / pyunits.year
                            + self.total_variable_OM_cost[0] * self.capacity_factor
                        )
                        / (
                            pyunits.convert(
                                value(self.tonne_REE_capture) * pyunits.tonne,
                                to_units=pyunits.kg)
                            / pyunits.year
                        )
                        + self.additional_cost_of_electricity
                    )
                )

                if transport_cost is not None:
                    if isinstance(transport_cost, (Expression, ScalarExpression)) or (
                        isinstance(transport_cost, (Param, Var))
                        and transport_cost.get_units is None
                    ):
                        self.transport_cost = (
                            transport_cost
                            * CE_index_units
                            / pyunits.tonne
                            * self.tonne_REE_capture
                        )
                    else:
                        self.transport_cost = transport_cost * self.tonne_REE_capture

            else:  # except the case where transport_cost is passed but tonne_REE_capture is not passed
                if transport_cost is not None:
                    raise Exception(
                        "If a transport_cost is not None, "
                        "tonne_REE_capture cannot be None."
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
            var_dict["Total Bare Erected Cost [$MM]"] = value(
                self.total_BEC
            )

        if hasattr(self, "total_installation_cost"):
            var_dict["Total Installation Cost [$MM]"] = value(
                self.total_installation_cost
            )

        if hasattr(self, "ancillary_costs"):
            var_dict["Total Ancillary Installation Cost [$MM]"] = value(
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
            var_dict[
                "Total Ancillary Instrumentation Installation Cost [$MM]"
            ] = value(self.instrumentation_costs)

        if hasattr(self, "plant_services_costs"):
            var_dict[
                "Total Ancillary Plant Services Installation Cost [$MM]"
            ] = value(self.plant_services_costs)

        if hasattr(self, "buildings_costs"):
            var_dict["Total Buildings Installation Cost [$MM]"] = value(
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
            var_dict[
                "Total Site Improvements Buildings Installation Cost [$MM]"
            ] = value(self.site_improvements_costs)

        if hasattr(self, "epcm_costs"):
            var_dict["Total EPCM Installation Cost [$MM]"] = value(self.epcm_costs)

        if hasattr(self, "equipment_installation_costs"):
            var_dict[
                "Total Equipment Installation EPCM Installation Cost [$MM]"
            ] = value(self.equipment_installation_costs)

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
            var_dict["Total Contingency Installation Cost [$MM]"] = value(
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

        if hasattr(self, "maintenance_and_material_cost"):
            var_dict["Total Maintenance and Material Cost [$MM/year]"] = value(
                self.maintenance_and_material_cost
            )

        if hasattr(self, "admin_and_support_labor_cost"):
            var_dict["Total Admin Support and Labor Cost [$MM/year]"] = value(
                self.admin_and_support_labor_cost
            )

        if hasattr(self, "sales_patenting_and_research_cost"):
            var_dict["Total Sales, Patenting and Research Cost [$MM/year]"] = value(
                self.sales_patenting_and_research_cost
            )

        if hasattr(self, "property_taxes_and_insurance"):
            var_dict["Total Property Taxes and Insurance Cost [$MM/year]"] = value(
                self.property_taxes_and_insurance
            )

        if hasattr(self, "other_fixed_costs"):
            var_dict["Total Other Fixed Costs [$MM/year]"] = value(
                self.other_fixed_costs
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
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(blk.name)
            )

        # define costing library
        blk.library = "REE"

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
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
            for (
                new_costing_params
            ) in additional_costing_params:  # merge new dictionaries sequentially
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
                                        "Data already exists for Account {} "
                                        "using source {}. "
                                        "Please confirm that the custom "
                                        "account dictionary is correct, or "
                                        "add the new parameters as a new "
                                        "account. To use the custom account "
                                        "dictionary for all conflicts, please "
                                        "pass the argument use_additional_costing_params "
                                        "as True.".format(accountkey, str(sourcekey))
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
                reference_units[account] = costing_params[str(source)][cost_accounts[0]][
                    "Units"
                ]
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
                print(
                    "KeyError: Account {} could not be found in the "
                    "dictionary for source {}".format(account, str(source))
                )

        # check that all accounts use the same process parameter
        param_check = None
        for account in cost_accounts:
            param = process_params[account]
            if param_check is None:
                param_check = param
            elif param != param_check:
                raise ValueError(
                    "{} cost accounts selected do not use "
                    "the same process parameter".format(blk.name)
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
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0][0]
                                + "**"
                                + ref_units[0][1]
                                + "/"
                                + ref_units[1],
                            )
                        )
                elif "**" in ref_units[1]:
                    ref_units[1] = ref_units[1].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1][0]
                        ) ** int(ref_units[1][1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0]
                                + "/"
                                + ref_units[1][0]
                                + "**"
                                + ref_units[1][1],
                            )
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1]
                        )
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            else:
                if "**" in ref_units:
                    ref_units = ref_units.split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) ** int(ref_units[1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units)
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            if isinstance(scaled_param, list):
                for sp in scaled_param:
                    if sp.get_units() is None:
                        raise ValueError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (cost_accounts[0], ref_units, sp.get_units())
                        )
                    else:
                        try:
                            pyunits.convert(sp, ref_units)
                        except InconsistentUnitsError:
                            raise Exception(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    sp.get_units(),
                                )
                            )
            else:
                try:
                    if pyunits.get_units(scaled_param) is None:
                        raise UnitsError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (
                                cost_accounts[0],
                                ref_units,
                                pyunits.get_units(scaled_param),
                            )
                        )
                    else:
                        try:
                            pyunits.convert(scaled_param, ref_units)
                        except InconsistentUnitsError:
                            raise UnitsError(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    pyunits.get_units(scaled_param),
                                )
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
            doc="exponential parameter for account",
        )

        blk.ref_cost = Param(
            cost_accounts,
            mutable=True,
            initialize=reference_costs,
            doc="reference cost for account",
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
                    doc="reference parameter for account",
                )
        elif isinstance(process_params[cost_accounts[0]], str):
            blk.ref_param = Param(
                cost_accounts,
                mutable=True,
                initialize=reference_params,
                doc="reference parameter for account",
            )

        # define variables
        blk.bare_erected_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
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

    def get_fixed_OM_costs(
        b,
        labor_types=["skilled", "unskilled", "supervisor", "maintenance", "technician", "engineer"],
        labor_rate=[27.90, 23.26, 30.29, 24.06, 23.43, 46.82],
        labor_burden=25,
        operators_per_shift=[2, 5, 2, 3, 1, 2],
        hours_per_shift=8,
        shifts_per_day=3,
        operating_days_per_year=336,
        CE_index_year="2021",
    ):
        """
        Args:
            b: costing block to add fixed cost variables and constraints to
            net_power: actual plant output in MW, only required if calculating
                variable costs
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
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars

        Returns:
            None
        """

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
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
            units=pyunits.dimensionless
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

        # make vars
        b.annual_operating_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="annual labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.maintenance_and_material_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="maintenance and material cost in $MM/yr",
            units=CE_index_units,
        )
        b.quality_assurance_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="quality assurance cost in $MM/yr",
            units=CE_index_units,
        )
        b.sales_patenting_and_research_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="sales, patenting and research cost in $MM/yr",
            units=CE_index_units,
        )
        b.admin_and_support_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="admin and support labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.property_taxes_and_insurance = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="property taxes and insurance cost in $MM/yr",
            units=CE_index_units,
        )
        b.total_fixed_OM_cost = Var(
            initialize=4,
            bounds=(0, 1e4),
            doc="total fixed O&M costs in $MM/yr",
            units=CE_index_units,
        )

        # variable for user to assign other fixed costs to,
        # fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="other fixed costs in $MM/yr",
            units=CE_index_units,
        )
        b.other_fixed_costs.fix(0)

        # create constraints
        TPC = b.total_plant_cost  # quick reference to total_plant_cost
        operating_labor_types, direct_labor_types = [], []  # subset labor lists
        for i in labor_types:
            if i in ["skilled", "unskilled", "supervisor", "maintenance"]:
                operating_labor_types.append(i)
            elif i in ["technician", "engineer"]:
                direct_labor_types.append(i)
            else:
                raise ValueError(
                    "Value {} for labor_type is not allowed. "
                    "Allowed labor types for operating labor include skilled,"
                    "unskilled, supervisor and maintenance. Allowed labor types "
                    "for direct labor include technician and engineer.".format(i)
                    )

        # calculated from labor rate, labor burden, and operators per shift
        @b.Constraint()
        def annual_labor_cost_rule(c):
            return c.annual_operating_labor_cost == pyunits.convert(
                (
                    sum(c.operators_per_shift[i] * c.labor_rate[i]
                        for i in labor_types)
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units,
            )

        # maintenance cost is 1% of TPC
        @b.Constraint()
        def maintenance_and_material_cost_rule(c):
            return c.maintenance_and_material_cost == 0.01 * TPC

        # quality assurance cost is 10% of operating labor
        @b.Constraint()
        def quality_assurance_cost_rule(c):
            return c.quality_assurance_cost == 0.10 * pyunits.convert(
                (
                    sum(c.operators_per_shift[i] * c.labor_rate[i]
                        for i in operating_labor_types)
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units,
            )

        # sales cost is 0.5% of total revenue
        # supposed to be based on revenue, we don't have that calculation
        # for now, make it 20% of direct labor costs
        @b.Constraint()
        def sales_patenting_and_research_cost_rule(c):
            return c.sales_patenting_and_research_cost == 0.20 * pyunits.convert(
                (
                    sum(c.operators_per_shift[i] * c.labor_rate[i]
                        for i in direct_labor_types)
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units,
            )

        # admin cost is 15% of direct labor
        @b.Constraint()
        def admin_and_support_labor_cost_rule(c):
            return c.admin_and_support_labor_cost == 0.15 * pyunits.convert(
                (
                    sum(c.operators_per_shift[i] * c.labor_rate[i]
                        for i in direct_labor_types)
                    * (1 + c.labor_burden / 100)
                    * c.hours_per_shift
                    * c.shifts_per_day
                    * c.operating_days_per_year
                ),
                CE_index_units,
            )

        # taxes are 2% of TPC
        @b.Constraint()
        def taxes_and_insurance_cost_rule(c):
            return c.property_taxes_and_insurance == 0.02 * TPC

        # sum of fixed O&M costs
        @b.Constraint()
        def total_fixed_OM_cost_rule(c):
            return c.total_fixed_OM_cost == (
                c.annual_operating_labor_cost
                + c.maintenance_and_material_cost
                + c.admin_and_support_labor_cost
                + c.sales_patenting_and_research_cost
                + c.property_taxes_and_insurance
                + c.other_fixed_costs
            )

    def get_variable_OM_costs(
        b,
        CE_index_year="2021",
    ):
        """
        Args:
            b: costing block to add fixed cost variables and constraints to
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars

        Returns:
            None.

        """
        pass  # coming soon

    def initialize_fixed_OM_costs(b):
        # b is the flowsheet-level costing block
        pass  # coming soon

    def initialize_variable_OM_costs(b):
        # b is the flowsheet-level costing block
        pass  # coming soon

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
                            sum(
                                o.bare_erected_cost[key]
                                for key in o.bare_erected_cost
                            )
                        ),
                    )
                )

    def get_total_BEC(b, CE_index_year):
        # This method accepts a flowsheet-level costing block

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        BEC_list = []

        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "bare_erected_cost"):
                for key in o.bare_erected_cost.keys():
                    BEC_list.append(o.bare_erected_cost[key])

        b.total_BEC = Var(
            initialize=100,
            bounds=(0, 1e4),
            doc="total TPC in $MM",
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
        print("Total installed equipment cost: $%.3f Million" % value(b.total_plant_cost))
        print("Total fixed O&M cost: $%.3f Million" % value(b.total_fixed_OM_cost))
        print("\n")

    def calculate_REE_costing_bounds(b, capacity, grade, CE_index_year, recalculate=False):
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
            delattr(b, "costing_lower_bound_index")
            delattr(b, "costing_upper_bound_index")
            delattr(b, "costing_lower_bound_eq")
            delattr(b, "costing_upper_bound_eq")
            delattr(b, "costing_lower_bound_eq_index")
            delattr(b, "costing_upper_bound_eq_index")

        if not hasattr(b, "capacity"):
            b.capacity = Var(
                initialize=value(pyunits.convert(capacity, to_units=pyunits.tonnes)),
                bounds=(0, 1E9),
                doc="feedstock capacity of site",
                units=pyunits.tonnes)
            print("New variable 'capacity' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has attribute "
                  "'capacity', moving on. Set 'recalculate' to True to delete "
                  "old objects and recalculate for new inputs".format(b.name))

        if not hasattr(b, "grade"):
            b.grade = Var(
                initialize=value(pyunits.convert(grade, to_units=pyunits.percent)),
                bounds=(0, 100),
                doc="grade percentage of site",
                units=pyunits.percent
                )
            print("New variable 'grade' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has attribute "
                  "'grade', moving on. Set 'recalculate' to True to delete "
                  "old objects and recalculate for new inputs".format(b.name))

        processes = {
            "Total Capital": [81, 1.4, -0.46, 0.063],
            "Total Operating": [27, 0.87, -0.087, 0.038],
            "Beneficiation": [2.7, 1.3, -0.15, 0.062],
            "Beneficiation, Chemical Extraction, Enrichment and Separation": [22, 1.28,  -0.059, 0.046],
            "Chemical Extraction": [40, 2.9, -0.46, 0.14],
            "Chemical Extraction, Enrichment and Separation": [15, 15, -0.19, 0.28],
            "Enrichment and Separation": [6.7, 2.8, -0.16, 0.11],
            "Mining": [25, 2.5, -0.32, 0.095],
            }


        if not hasattr(b, "costing_lower_bound"):
            b.costing_lower_bound = Var(
                processes,
                initialize=10,
                bounds=(0, 100),
                doc="estimated lower bound on per unit production cost of site",
                units=getattr(pyunits, "USD_" + CE_index_year)/pyunits.kg
                )
            print("New variable 'costing_lower_bound' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has attribute "
                  "'costing_lower_bound', moving on. Set 'recalculate' to True to delete "
                  "old objects and recalculate for new inputs".format(b.name))

        if not hasattr(b, "costing_upper_bound"):
            b.costing_upper_bound = Var(
                processes,
                initialize=10,
                bounds=(0, 100),
                doc="estimated upper bound on per unit production cost of site",
                units=getattr(pyunits, "USD_" + CE_index_year)/pyunits.kg
                )
            print("New variable 'costing_upper_bound' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has attribute "
                  "'costing_upper_bound', moving on. Set 'recalculate' to True to delete "
                  "old objects and recalculate for new inputs".format(b.name))

        if not hasattr(b, "costing_lower_bound_eq"):
            @b.Constraint(processes)
            def costing_lower_bound_eq(c, p):
                return (
                    c.costing_lower_bound[p] == pyunits.convert(
                        pyunits.USD_2022 * (processes[p][0] - processes[p][1]) * (
                            pyunits.convert(
                                c.grade,
                                to_units=pyunits.dimensionless)/pyunits.dimensionless *
                            pyunits.convert(
                                c.capacity,
                                to_units=pyunits.tonnes)/pyunits.tonnes
                            ) ** (processes[p][2] - processes[p][3]),
                    to_units=getattr(pyunits, "USD_" + CE_index_year)
                    ) / pyunits.kg
                    )

            # assume model is already solved, so just calculate costing bounds here
            for i in b.costing_lower_bound.keys():
                calculate_variable_from_constraint(
                    b.costing_lower_bound[i],
                    b.costing_lower_bound_eq[i],
                )
            print("New constraint 'costing_lower_bounding_eq' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has indexed "
                  "constraint 'costing_lower_bounding_eq', reporting existing results. "
                  "Set 'recalculate' to True to delete old objects and recalculate "
                  "for new inputs".format(b.name))

        if not hasattr(b, "costing_upper_bound_eq"):
            @b.Constraint(processes)
            def costing_upper_bound_eq(c, p):
                return (
                    c.costing_upper_bound[p] == pyunits.convert(
                        pyunits.USD_2022 * (processes[p][0] + processes[p][1]) * (
                            pyunits.convert(
                                c.grade,
                                to_units=pyunits.dimensionless)/pyunits.dimensionless *
                            pyunits.convert(
                                c.capacity,
                                to_units=pyunits.tonnes)/pyunits.tonnes
                            ) ** (processes[p][2] + processes[p][3]),
                    to_units=getattr(pyunits, "USD_" + CE_index_year)
                    ) / pyunits.kg
                    )

            # assume model is already solved, so just calculate costing bounds here
            for i in b.costing_upper_bound.keys():
                calculate_variable_from_constraint(
                    b.costing_upper_bound[i],
                    b.costing_upper_bound_eq[i],
                )
            print("New constraint 'costing_upper_bounding_eq' created as attribute of {}".format(b.name))
        else:
            print("Flowsheet-level costing block {} already has indexed "
                  "constraint 'costing_upper_bounding_eq', reporting existing results. "
                  "Set 'recalculate' to True to delete old objects and recalculate "
                  "for new inputs".format(b.name))

        print("\nPrinting calculated costing bounds for processes:")
        for p in processes:
            print(p, ": [",
                  value(b.costing_lower_bound[p]), ", ",
                  value(b.costing_upper_bound[p]), "]",
                  getattr(pyunits, "USD_" + CE_index_year), "/kg"
                  )

        # method has finished building components
        b.components_already_built = True
