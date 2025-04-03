#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Author: Lingyan Deng

Figure 1 shows the byproduct recovery determination tree, referenced from Chapter 8 of book by Towler, 
Gavin, and Ray Sinnott. Chemical engineering design: principles, practice and economics of plant and 
process design. Butterworth-Heinemann, 2021. The waste disposal cost, estimated cost of 
conversion, and estimated added costs of by-product recovery should all be annualized cost. The
estimate potential revenue is also annually based. 
.. figure:: ../byproduct_recovery_determination_tree.png
    :width: 800
    :align: center


"""
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.environ import ConcreteModel, Expression, NonNegativeReals, Set, Var
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
    register_idaes_currency_units,
)

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ByproductRecovery")
class ByproductRecoveryData(UnitModelBlockData):
    CONFIG = ConfigDict()
    CONFIG.declare(
        "dynamic", ConfigValue(default=False, description="Steady-state model")
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(default=False, description="No holdup in this process"),
    )
    CONFIG.declare(
        "materials",
        ConfigValue(default=[], domain=list, description="List of materials to model"),
    )
    register_idaes_currency_units()

    def build(self):
        """
        Build the Byproduct Recovery Model supporting a flexible number of materials.
        """
        super().build()
        self.base_currency = pyunits.USD_2021

        # Set of materials (Dynamic Initialization)
        if not self.config.materials:
            raise ValueError(
                "⚠️ Material list cannot be empty! Provide at least one material."
            )

        self.materials = Set(initialize=self.config.materials)

        # Model Variables (Indexed by Material)
        self.material_production = Var(
            self.materials,
            initialize=100,
            units=pyunits.kg,
            doc="Amount of material produced",
        )
        self.market_value = Var(
            self.materials,
            initialize=10,
            units=pyunits.USD_2021 / pyunits.kg,
            doc="Market value per unit",
        )
        self.waste_disposal_cost = Var(
            self.materials,
            initialize=5,
            units=pyunits.USD_2021 / pyunits.kg,
            doc="Saved cots from avoiding waste disposal",
        )
        self.conversion_possible = Var(
            self.materials,
            initialize=0,
            domain=NonNegativeReals,
            doc="1 if conversion is needed, 0 otherwise",
        )
        self.conversion_cost = Var(
            self.materials,
            initialize=0,
            units=pyunits.USD_2021,
            doc="Cost of conversion if applicable",
        )
        self.added_process_steps = Var(
            self.materials,
            initialize=0,
            domain=NonNegativeReals,
            doc="1 if additional processing is required, 0 otherwise",
        )
        self.added_process_cost = Var(
            self.materials,
            initialize=0,
            units=pyunits.USD_2021,
            doc="Additional cost for purification",
        )

        # Revenue Calculation (Summing Over Materials)
        self.potential_revenue = Expression(
            expr=sum(
                self.material_production[m]
                * (self.market_value[m] + self.waste_disposal_cost[m])
                for m in self.materials
            ),
            doc="Total revenue from byproduct recovery across materials",
        )

        # Cost Calculation (Summing Over Materials)
        self.total_recovery_cost = Expression(
            expr=sum(
                self.conversion_cost[m] * self.conversion_possible[m]
                + self.added_process_cost[m] * self.added_process_steps[m]
                for m in self.materials
            ),
            doc="Total cost of byproduct recovery across materials",
        )

        # Net Benefit Calculation
        self.net_benefit = Expression(
            expr=self.potential_revenue - self.total_recovery_cost,
            doc="Net benefit of byproduct recovery across all materials",
        )

    def determine_financial_viability(self):
        """
        Evaluate whether byproduct recovery is financially viable.
        """
        net_benefit_value = value(self.net_benefit)
        if net_benefit_value > 0:
            return f"✅ Byproduct recovery is financially viable. Net Benefit: ${net_benefit_value:.2f}"
        else:
            return f"❌ Byproduct recovery is NOT financially viable. Loss: ${-net_benefit_value:.2f}"


def determine_example_usage():
    """
    Example Usage of ByproductRecovery within a Pyomo ConcreteModel.
    Returns the financial viability result.
    """
    # Define material list dynamically
    material_list = ["Aluminum", "Iron", "Copper", "Zinc"]

    # Create a Pyomo ConcreteModel
    model = ConcreteModel()

    # Add an instance of the ByproductRecovery process block
    model.recovery = ByproductRecovery(materials=material_list)

    # Define input values dynamically based on provided materials
    material_data = {
        "Aluminum": {
            "production": 100,
            "market_value": 10,
            "waste_disposal": 5,
            "conversion": 1,
            "conversion_cost": 200,
            "process_steps": 1,
            "process_cost": 50,
        },
        "Iron": {
            "production": 150,
            "market_value": 8,
            "waste_disposal": 4,
            "conversion": 0,
            "conversion_cost": 0,
            "process_steps": 0,
            "process_cost": 0,
        },
        "Copper": {
            "production": 50,
            "market_value": 15,
            "waste_disposal": 3,
            "conversion": 1,
            "conversion_cost": 100,
            "process_steps": 1,
            "process_cost": 30,
        },
        "Zinc": {
            "production": 80,
            "market_value": 12,
            "waste_disposal": 2,
            "conversion": 0,
            "conversion_cost": 0,
            "process_steps": 1,
            "process_cost": 20,
        },
    }

    # Set values dynamically based on the material list
    for m in material_list:
        data = material_data.get(m, {})  # Default to empty if material not in dict
        model.recovery.material_production[m].set_value(data.get("production", 0))
        model.recovery.market_value[m].set_value(data.get("market_value", 0))
        model.recovery.waste_disposal_cost[m].set_value(data.get("waste_disposal", 0))
        model.recovery.conversion_possible[m].set_value(data.get("conversion", 0))
        model.recovery.conversion_cost[m].set_value(data.get("conversion_cost", 0))
        model.recovery.added_process_steps[m].set_value(data.get("process_steps", 0))
        model.recovery.added_process_cost[m].set_value(data.get("process_cost", 0))

    # Compute net benefit
    net_benefit_value = value(model.recovery.net_benefit)

    # Evaluate the financial viability
    result = model.recovery.determine_financial_viability()

    return result, net_benefit_value


# Run the example only when executed as a script
if __name__ == "__main__":
    result = determine_example_usage()
    print("\n--- Byproduct Recovery Decision ---")
    print(result)
