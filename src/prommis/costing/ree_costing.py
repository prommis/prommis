#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Power Plant costing library
This method leverages NETL costing capabilities. Two main methods have been
developed to calculate the capital cost of power generation plants:

1. Fossil fueled power plants (from SCPC to IGCC) (get_PP_costing)
2. Supercritical CO2 power cycles (direct and indirect) (get_sCO2_unit_cost)

Other methods:

    - check_sCO2_costing_bounds() to display a warning if costing model have
      been used outside the range that where designed for
    - get_ASU_cost() to cost air separation units
"""

# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (A. Noring, A. Deshpande, B. Paul, D. Caballero, and M. Zamarripa)"
)
__version__ = "1.0.0"

from pyomo.environ import units as pyunits
from pyomo.environ import value, Var
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.logger as idaeslog
from idaes.core import declare_process_block_class
from idaes.models.costing.QGESS import QGESSCostingData
from idaes.core.base.costing_base import UnitModelCostingBlock

from prommis.costing.ree_costing_dictionaries import (
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
        # append REE_costing_params to the existing list, but put it first so
        # that name conflict can check against it during loading
        costing_method_arguments["additional_costing_params"].insert(0, REE_costing_params)

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
    CONFIG["tech"] = 10

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        super().build_global_params()
        self.config.tech = 10

        # Set the base year for all costs
        self.base_currency = pyunits.USD_2021
        # Set a base period for all operating costs
        self.base_period = pyunits.year

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
