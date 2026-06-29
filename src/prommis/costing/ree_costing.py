#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
REE costing library
This method leverages NETL costing capabilities.

- calculate_REE_costing_bounds() to provide an estimate of costing bounds

Other methods for calculating capital and equipment costs live in IDAES.
Methods for byproduct recovery and uncertainty quantification live in other
modules in this directory.
"""

# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, L. Deng, and M. Zamarripa)"
)
__version__ = "1.0.0"

from functools import partial

from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Expression, Var
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import UnitModelCostingBlock
from idaes.models.costing.QGESS import QGESSCostingData

# this is a duplicate reference to reporting method in IDAES, remove later
from idaes.models_extra.power_generation.costing.power_plant_costing_dictionaries import (
    report,
)

from prommis.costing.ree_costing_dictionaries import (
    load_default_resource_prices,
    load_default_sale_prices,
    load_REE_costing_dictionary,
    register_ree_currency_units,
)

_log = idaeslog.getLogger(__name__)

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------


def REEUnitModelCostingBlock(
    flowsheet_costing_block,
    costing_method,
    costing_method_arguments,
):
    """
    Wrapper method for constructing costing blocks for PrOMMiS REE unit models.

    Passes to IDAES UnitModelCostingBlock and sets some defaults
    """

    if "tech" not in costing_method_arguments:
        costing_method_arguments["tech"] = 10  # UKy

    if "ccs" not in costing_method_arguments:
        costing_method_arguments["ccs"] = "A"  # no ccs

    REE_costing_params = load_REE_costing_dictionary()  # UKy
    if "additional_costing_params" not in costing_method_arguments:
        costing_method_arguments["additional_costing_params"] = []
        costing_method_arguments["additional_costing_params"].append(REE_costing_params)
    else:
        # do some type checking here before calling the QGESS method

        if not (
            isinstance(costing_method_arguments["additional_costing_params"], list)
            and all(
                isinstance(item, dict)
                for item in costing_method_arguments["additional_costing_params"]
            )
        ):
            raise TypeError(
                "additional_costing_params must be a list of dicts, not a single dict, "
                "e.g. [{'1': data}, {'2': data},] or [{'1': data},] and not {'1': data}."
            )

        # create a list with the dictionaties
        # dictionaries_to_append = [REE_costing_params, costing_method_arguments["additional_costing_params"]]
        # costing_method_arguments["additional_costing_params"] = dictionaries_to_append
        costing_method_arguments["additional_costing_params"].insert(
            0, REE_costing_params
        )

    if "CEPCI_year" not in costing_method_arguments:
        costing_method_arguments["CEPCI_year"] = "2021"

    costing = UnitModelCostingBlock(
        flowsheet_costing_block=flowsheet_costing_block,
        costing_method=costing_method,
        costing_method_arguments=costing_method_arguments,
    )

    costing.costing_method_arguments = costing_method_arguments

    return costing


@declare_process_block_class("REECosting")
class REECostingData(QGESSCostingData):
    # Register currency and conversion rates based on CEPCI
    register_ree_currency_units()

    CONFIG = QGESSCostingData.CONFIG()

    # IDAES users must explicitly select a power plant tech
    # for PrOMMiS, if users don't select a tech value, assume they want 10 (UKy)
    # make 10 the default, but leave the domain as int so users can use tech 1-9 if they want to
    del CONFIG["tech"]

    CONFIG.declare(
        "tech",
        ConfigValue(
            default=10,
            domain=int,
            description="Integer corresponding to supported technology libraries, where 1-7 are various power plant types, 8-9 are specific case studies, and 10 is UKy REE.",
        ),
    )

    del CONFIG["CEPCI_year"]

    CONFIG.declare(
        "CEPCI_year",
        ConfigValue(
            default="2021",
            domain=str,
            description="Basis year for costing. Must be a supported value from 1990 to 2023 or a user-defined value. For details, see the IDAES 'costing_base.py' module",
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
        super().build_global_params()

        # Set the base year for all costs
        self.base_currency = getattr(pyunits, "USD_" + self.config.CEPCI_year)
        # Set a base period for all operating costs
        self.base_period = pyunits.year

        # TODO remove, mock some changes that need to patched into the IDAES QGESS module later
        if isinstance(self.Lang_factor, Expression):
            installation_components = {
                "piping_materials_and_labor_percentage": 20,
                "electrical_materials_and_labor_percentage": 20,
                "instrumentation_percentage": 8,
                "plants_services_percentage": 10,
                "process_buildings_percentage": 40,
                "auxiliary_buildings_percentage": 15,
                "site_improvements_percentage": 10,
                "equipment_installation_percentage": 17,
                "field_expenses_percentage": 12,
                "project_management_and_construction_percentage": 30,
                "process_contingency_percentage": 15,
            }

            delattr(self, "Lang_factor")
            self.Lang_factor = Expression(
                expr=pyunits.convert(
                    sum(
                        self.installation_components[k] for k in installation_components
                    ),
                    to_units=pyunits.dimensionless,
                )
                + 1  # added this
            )

    def build_REE_process_costs(self, **kwargs):
        """
        Wrapper method for building PrOMMiS REE process costs.

        Passes to QGESS build_process_costs() and sets some defaults
        """

        # TODO change to use additional dictionary arguments instead of looping through prices

        REE_resource_prices = load_default_resource_prices()

        if "resource_prices" not in kwargs:
            kwargs["resource_prices"] = REE_resource_prices
        else:
            # do some type checking before appending the REE default resource prices
            if not isinstance(kwargs["resource_prices"], dict):
                raise TypeError("Dictionary of resource prices must be a dict object.")
            for k in REE_resource_prices:
                kwargs["resource_prices"][k] = REE_resource_prices[k]

        REE_sale_prices = load_default_sale_prices()

        if "sale_prices" not in kwargs:
            kwargs["sale_prices"] = REE_sale_prices
        else:
            # do some type checking before appending the REE default sale prices
            if not isinstance(kwargs["sale_prices"], dict):
                raise TypeError(
                    "Dictionary of custom sale_prices must be a dict object."
                )
            for k in REE_sale_prices:
                kwargs["sale_prices"][k] = REE_sale_prices[k]

        QGESSCostingData.build_process_costs(self, **kwargs)

        # TODO remove, mock some changes that need to patched into the IDAES QGESS module later
        if hasattr(self, "total_TPC_eq"):
            delattr(self, "total_TPC_eq")

            @self.Constraint()
            def total_TPC_eq(c):
                # TPC = BEC that needs Lang factor applied to get TPC + TPC that is already calculated
                # apply location factor and economy of numbers here too
                return c.total_TPC == (
                    (
                        (
                            (sum(c.BEC_list) - sum(c.BEC_blocks_with_TPC_list))
                            * c.Lang_factor
                            + sum(c.TPC_blocks_with_TPC_list)
                        )
                        + c.other_plant_costs  # added this
                    )
                    # apply economy of numbers if enabled
                    * (
                        c.NOAK_factor if c.config.has_economy_of_numbers else 1
                    )  # applied to whole TPC
                )

        if hasattr(self, "total_variable_cost_eq"):
            delattr(self, "total_variable_cost_eq")

            @self.Constraint(self.parent_block().time)
            def total_variable_cost_eq(c, t):
                return (
                    c.total_variable_OM_cost[t]
                    == sum(
                        c.variable_operating_costs[t, r] for r in kwargs["resources"]
                    )
                    + (
                        c.plant_overhead_cost[t]
                        if hasattr(c, "plant_overhead_cost")
                        else 0 * c.CEPCI_units / pyunits.year
                    )
                    + (
                        c.land_cost
                        if c.land_cost_reoccurrence == "annual"
                        else 0 * c.CEPCI_units / pyunits.year
                    )
                    + (
                        c.additional_chemicals_cost
                        if c.additional_chemicals_cost_reoccurrence == "annual"
                        else 0 * c.CEPCI_units / pyunits.year
                    )
                    + (
                        c.additional_waste_cost
                        if c.additional_waste_cost_reoccurrence == "annual"
                        else 0 * c.CEPCI_units / pyunits.year
                    )
                    # for power plants, include maintenance material costs here if defined
                    + (
                        c.maintenance_material_cost
                        if hasattr(c, "maintenance_material_cost")
                        else 0 * c.CEPCI_units / pyunits.year
                    )
                    + c.other_variable_costs[t]  # added this
                    + c.custom_variable_costs  # added this
                )

        if hasattr(self, "capex"):  # this is used for the NPV calculation
            delattr(self, "capex")
            self.capex = Expression(
                expr=(
                    self.total_TPC
                    # removed double counting of other_plant_costs
                    + (
                        self.land_cost
                        if self.land_cost_reoccurrence == "one_time"
                        else 0 * self.CEPCI_units
                    )
                )
            )

    # TODO remove, mock reference to report method that needs to be patched into IDAES QGESS
    def report(self, export=False):
        """
        Call report method.
        """
        return report(self, export=export)

    def calculate_REE_costing_bounds(
        b,
        capacity,
        grade,
        report=False,
    ):
        # adapted from https://doi.org/10.1038/s41893-023-01145-1
        # This method accepts a flowsheet-level costing block
        # capacity and grade should be variables with Pyomo units,
        # or values with Pyomo unit containers

        b.capacity = Var(
            initialize=value(pyunits.convert(capacity, to_units=pyunits.tonnes)),
            bounds=(0, None),
            doc="Feedstock capacity of site",
            units=pyunits.tonnes,
        )
        b.capacity.fix(capacity)

        b.grade = Var(
            initialize=value(pyunits.convert(grade, to_units=pyunits.percent)),
            bounds=(0, None),
            doc="Grade percentage of site. The value should be a "
            "percentage, for example 10 for 10%.",
            units=pyunits.percent,
        )
        b.grade.fix(grade)

        b.processes = {
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

        b.costing_lower_bound = Var(
            b.processes,
            initialize=1,
            bounds=(0, None),
            doc="Estimated lower bound on per unit production cost of site",
            units=b.base_currency / pyunits.kg,
        )

        b.costing_upper_bound = Var(
            b.processes,
            initialize=1,
            bounds=(0, None),
            doc="Estimated upper bound on per unit production cost of site",
            units=b.base_currency / pyunits.kg,
        )

        def rule_costing_lower_bound_eq(c, p, ps):
            return (
                c.costing_lower_bound[p]
                == pyunits.convert(
                    pyunits.USD_2022
                    * (ps[p][0] - ps[p][1])
                    * (
                        pyunits.convert(c.grade, to_units=pyunits.dimensionless)
                        * pyunits.convert(c.capacity, to_units=pyunits.tonnes)
                        / pyunits.tonnes
                    )
                    ** (ps[p][2] - ps[p][3]),
                    to_units=b.base_currency,
                )
                / pyunits.kg
            )

        b.costing_lower_bound_eq = Constraint(
            [p for p in b.processes.keys()],
            rule=partial(rule_costing_lower_bound_eq, ps=b.processes),
        )

        # assume model is already solved, so just calculate costing bounds here
        for i in b.costing_lower_bound.keys():
            calculate_variable_from_constraint(
                b.costing_lower_bound[i],
                b.costing_lower_bound_eq[i],
            )

        def rule_costing_upper_bound_eq(c, p, ps):
            return (
                c.costing_upper_bound[p]
                == pyunits.convert(
                    pyunits.USD_2022
                    * (ps[p][0] + ps[p][1])
                    * (
                        pyunits.convert(c.grade, to_units=pyunits.dimensionless)
                        * pyunits.convert(c.capacity, to_units=pyunits.tonnes)
                        / pyunits.tonnes
                    )
                    ** (ps[p][2] + ps[p][3]),
                    to_units=b.base_currency,
                )
                / pyunits.kg
            )

        b.costing_upper_bound_eq = Constraint(
            [p for p in b.processes.keys()],
            rule=partial(rule_costing_upper_bound_eq, ps=b.processes),
        )

        # assume model is already solved, so just calculate costing bounds here
        for i in b.costing_upper_bound.keys():
            calculate_variable_from_constraint(
                b.costing_upper_bound[i],
                b.costing_upper_bound_eq[i],
            )

        if report:
            REECostingData.report_costing_bounds(b)

    def report_costing_bounds(b):
        """
        Display calculated costing bounds.
        """

        print("\n\nPrinting calculated costing bounds for processes:")
        for p in b.processes:
            print(
                p,
                ": [",
                value(b.costing_lower_bound[p]),
                ", ",
                value(b.costing_upper_bound[p]),
                "]",
                b.CEPCI_units,
                "/kg",
            )
