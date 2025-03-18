#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Author: Lingyan Deng, Brandon Paul

This is a test file for testing the byproduct recovery framework. The example is a Li-Co diafiltration process to recover Li and Co as byproduct from battery recovery. Li-Co recovery
process is referenced from literature: Wamble, N.P., Eugene, E.A., Phillip, W.A., Dowling, A.W., 'Optimal Diafiltration Membrane Cascades Enable Green Recycling of Spent Lithium-Ion 
Batteries', ACS Sustainable Chem. Eng. 2022, 10, 12207−12225. 

The steps to test this framework including:
1. Import packages
2. Build the li-co diafiltration recovery flowsheet
3. add costing
4. Use the byproduct recovery framework to determine if li should be recovered 
5. Solve and display results


"""
# 1. Import packages
# Pyomo packages
from pyomo.environ import (
    ConcreteModel,
    Expression,
    Param,
    TransformationFactory,
    Var,
    value,
)

# IDAES packages
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

# PrOMMiS packages
from prommis.nanofiltration.diafiltration import (
    add_costing,
    add_objective,
    add_product_constraints,
    build_model,
    initialize_model,
    print_information,
    set_scaling,
    solve_model,
    unfix_opt_variables,
)
from prommis.uky.costing.costing_dictionaries import load_default_sale_prices
from prommis.uky.costing.determine_byproduct_recovery import (
    ByproductRecovery,
    determine_example_usage,
)

# Byproduct recovery determine package


class TestLiCoDiafiltration:
    # 2. Build the li-co diafiltration recovery flowsheet
    def build_LiCoDiafiltration_model(self):
        self.m = build_model()

        assert (
            degrees_of_freedom(self.m) == 0
        ), "Degrees of freedom should be zero after adding costing."
        add_costing(self.m)
        initialize_model(self.m)
        solve_model(self.m)

        unfix_opt_variables(self.m)
        add_product_constraints(self.m, Li_recovery_bound=0.95, Co_recovery_bound=0.635)
        add_objective(self.m)
        set_scaling(self.m)
        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(self.m, rename=False)
        solve_model(scaled_model)
        # Propagate results back to unscaled model
        scaling.propagate_solution(scaled_model, self.m)

        # Ensure feed pump OPEX is negligible
        assert (
            value(self.m.fs.feed_pump.costing.variable_operating_cost) < 0.005
        ), "Feed pump operating cost should be negligible."

        print_information(self.m)

        # Check basic model properties
        assert isinstance(
            self.m.fs.costing.total_annualized_cost, Expression
        ), "total_annualized_cost should be an Expression."
        assert isinstance(
            self.m.fs.stage3.permeate_outlet.flow_vol, (Var, Expression, Param)
        ), "stage3 permeate flow_vol should be an Var, Expression, or Param."
        assert isinstance(
            self.m.fs.stage1.retentate_outlet.flow_vol, (Var, Expression, Param)
        ), "stage1 retentate flow_vol should be an Var, Expression, or Param."

        # Store results for later tests
        self.total_annualized_cost = value(self.m.fs.costing.total_annualized_cost)

        # recovery mass flow rate kg/hr
        self.Li_recovery_mass = value(
            self.m.fs.stage3.permeate_outlet.flow_vol[0]
        ) * value(self.m.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"])
        self.Co_recovery_mass = value(
            self.m.fs.stage1.retentate_outlet.flow_vol[0]
        ) * value(self.m.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"])

    @pytest.mark.component
    def test_structural_issues(self):
        self.build_LiCoDiafiltration_model()
        dt = DiagnosticsToolbox(self.m)
        dt.report_structural_issues()
        dt.assert_no_numerical_warnings()
        dt.display_underconstrained_set()
        dt.display_potential_evaluation_errors()

    # 3. Access product price.
    def test_import_product_prices(self):
        """Test case for importing and verifying product prices."""

        self.build_LiCoDiafiltration_model()  # Ensure model is built first
        sale_prices = load_default_sale_prices()

        assert "Li" in sale_prices, "Lithium price missing in sale prices dictionary."
        assert "Co" in sale_prices, "Cobalt price missing in sale prices dictionary."

        self.Li_price = sale_prices["Li"]
        self.Co_price = sale_prices["Co"]

        # Check reasonable price values
        assert (
            value(self.Co_recovery_mass) > 0
        ), "Cobalt recovery mass should be positive."
        assert (
            value(self.Li_recovery_mass) > 0
        ), "Lithium recovery mass should be positive."

        print(f"Li price: {self.Li_price}")
        print(f"Co price: {self.Co_price}")

    # 4. Test framework to determine if the byproduct should be recovered
    def test_determine_byproduct_recovery(self):
        self.build_LiCoDiafiltration_model()  # Ensure model is built first
        self.test_import_product_prices()  # Ensure prices are loaded
        material_list = ["Lithium", "Cobalt"]  # Allow user input for flexible materials

        model = ConcreteModel()
        model.recovery_determine = ByproductRecovery(materials=material_list)

        # Define input values dynamically based on provided materials
        material_data = {
            "Lithium": {
                "production": self.Li_recovery_mass,
                "market_value": self.Li_price,
                "waste_disposal": 1,
                "conversion": 0,
                "conversion_cost": 0,
                "process_steps": 1,
                "process_cost": self.total_annualized_cost,
            },
            "Cobalt": {
                "production": self.Co_recovery_mass,
                "market_value": self.Co_price,
                "waste_disposal": 1,
                "conversion": 0,
                "conversion_cost": 0,
                "process_steps": 0,
                "process_cost": 0,
            },
        }

        # Set values dynamically based on material list
        for m in material_list:
            data = material_data.get(m, {})
            model.recovery_determine.material_production[m].set_value(
                data.get("production", 0)
            )
            model.recovery_determine.market_value[m].set_value(
                data.get("market_value", 0)
            )
            model.recovery_determine.waste_disposal_cost[m].set_value(
                data.get("waste_disposal", 0)
            )
            model.recovery_determine.conversion_possible[m].set_value(
                data.get("conversion", 0)
            )
            model.recovery_determine.conversion_cost[m].set_value(
                data.get("conversion_cost", 0)
            )
            model.recovery_determine.added_process_steps[m].set_value(
                data.get("process_steps", 0)
            )
            model.recovery_determine.added_process_cost[m].set_value(
                data.get("process_cost", 0)
            )

        potential_revenue = value(model.recovery_determine.potential_revenue)
        assert potential_revenue >= 0, "Potential revenue should be non-negative."

        determine_result = model.recovery_determine.determine_financial_viability()
        assert isinstance(
            determine_result, str
        ), f"Expected a string message, but got {type(determine_result)}"

        # Check the output string for financial viability
        net_benefit_value = value(model.recovery_determine.net_benefit)
        assert net_benefit_value == pytest.approx(-211162.615, rel=1e-4)

        if net_benefit_value > 0:
            expected_message = f"✅ Byproduct recovery is financially viable. Net Benefit: ${net_benefit_value:.2f}"
        else:
            expected_message = f"❌ Byproduct recovery is NOT financially viable. Loss: ${-net_benefit_value:.2f}"

        assert (
            determine_result == expected_message
        ), f"Unexpected output: {determine_result}"

        print("\n--- Byproduct Recovery Decision ---")
        print(determine_result)

    def test_example_usage(self):
        """
        Ensure that the example usage from script A runs correctly.
        """
        determine_result, net_benefit_value = determine_example_usage()
        assert net_benefit_value == pytest.approx(4920.00, rel=1e-4)
        # Ensure result is a string
        assert isinstance(
            determine_result, str
        ), f"Expected a string message, but got {type(determine_result)}"

        # Check that the result follows the expected pattern
        assert determine_result.startswith(
            "✅ Byproduct recovery is financially viable."
        ) or determine_result.startswith(
            "❌ Byproduct recovery is NOT financially viable."
        ), f"Unexpected output: {determine_result}"

    def test_results(self):
        """Check expected numerical and string outputs."""
        self.build_LiCoDiafiltration_model()  # Ensure model is built first
        print(f"Co recovery mass flow rate: {self.Co_recovery_mass} kg/hr")
        print(f"Li recovery mass flow rate: {self.Li_recovery_mass} kg/hr")

        # Ensure recovery masses are reasonable
        assert (
            value(self.Co_recovery_mass) > 0
        ), "Cobalt recovery mass should be positive."
        assert (
            value(self.Li_recovery_mass) > 0
        ), "Lithium recovery mass should be positive."

    def test_edge_case_empty_material_list(self):
        # Test that an empty material list raises the correct ValueError
        with pytest.raises(
            ValueError,
            match="⚠️ Material list cannot be empty! Provide at least one material.",
        ):
            empty_model = ConcreteModel()
            empty_model.recovery_determine = ByproductRecovery(materials=[])
