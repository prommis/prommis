#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Param,
    RangeSet,
    Set,
    SolverFactory,
    Var,
    assert_optimal_termination,
    value,
)

from idaes.core.scaling import get_scaling_factor
from idaes.core.solvers import get_solver

import pytest

from prommis.superstructure.objective_function_enums import ObjectiveFunctionChoice
from prommis.superstructure.report_superstructure_results import (
    report_superstructure_costing,
    report_superstructure_environmental_impacts,
    report_superstructure_results_overview,
    report_superstructure_streams,
)
from prommis.superstructure.superstructure_function import (
    SuperstructureScaler,
    build_model,
    define_custom_units,
)
from unittest.mock import patch

solver_available = SolverFactory("gurobi").available()
if solver_available:
    solver = get_solver(solver="gurobi")
else:
    solver = None

### Define common parameters
# define custom units
define_custom_units()
obj_func = ObjectiveFunctionChoice.NET_PRESENT_VALUE

plant_start = 2024
plant_lifetime = 5

available_feed = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
}
collection_rate = 0.1
tracked_comps = ["Nd", "Fe"]
prod_comp_mass = {
    "Nd": 0.206 * 3,
    "Fe": 0.691 * 3,
}

num_stages = 3
options_in_stage = {1: 1, 2: 2, 3: 3}
option_outlets = {
    (1, 1): [1, 2],
    (2, 1): [1],
    (2, 2): [2, 3],
}

option_efficiencies = {
    (1, 1): {"Nd": 1, "Fe": 1},
    (2, 1): {"Nd": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Fe": 1},
    (3, 1): {"Nd": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Fe": 0},
    (3, 3): {"Nd": 1, "Fe": 0},
}

profit = {
    (3, 1): {"Nd": 45.4272, "Fe": 0},
    (3, 2): {"Nd": 69.888, "Fe": 0},
    (3, 3): {"Nd": 45.4272, "Fe": 0},
}
opt_var_oc_params = {
    (2, 1): {"a": 0.0053, "b": 7929.7},
    (2, 2): {"a": 0.0015, "b": 2233.16},
    (3, 1): {"a": 15.594, "b": 4e6},
    (3, 2): {"a": 35.58463, "b": 4e6},
    (3, 3): {"a": 1.58, "b": 0},
}
operators_per_discrete_unit = {(1, 1): 1}
yearly_cost_per_unit = {(1, 1): 0}
capital_cost_per_unit = {(1, 1): 0}
processing_rate = {(1, 1): 7868}
num_operators = {
    (2, 1): 0.65,
    (2, 2): 0.65,
    (3, 1): 1.6,
    (3, 2): 1.6,
    (3, 3): 1.3,
}
labor_rate = 8000 * 38.20

discretized_purchased_equipment_cost = {
    (2, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            10130.08515,
            31353.21173,
            48788.84678,
            60305.81927,
            77063.4884,
            117214.7546,
            151018.0699,
            195698.5419,
        ],
    },
    (2, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            11702.08515,
            39023.21173,
            62134.84678,
            77539.81927,
            100113.4884,
            154792.7546,
            201326.0699,
            263374.5419,
        ],
    },
    (3, 1): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            343228.652,
            482425.4684,
            618182.0594,
            743750.2902,
            844443.0443,
            978479.5225,
            1183834.522,
            1440660.587,
        ],
    },
    (3, 2): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            643228.652,
            782425.4684,
            918182.0594,
            1043750.2902,
            1144443.0443,
            1278479.5225,
            1483834.522,
            1740660.587,
        ],
    },
    (3, 3): {
        "Flowrates": [
            0.0,
            36480.0,
            634240.0,
            1434800.0,
            2083760.0,
            3171200.0,
            6342400.0,
            9513600.0,
            14270400.0,
        ],
        "Costs": [
            0.0,
            4274.7216,
            30479.121,
            53459.01,
            69261.68,
            92510.61,
            143803.33,
            197644.75,
            261513.79,
        ],
    },
}

options_environmental_impacts = {
    (1, 1): 10,
    (2, 1): 1000,
    (2, 2): 10,
    (3, 1): 600,
    (3, 2): 600,
    (3, 3): 10,
}
epsilon = 1e16

byproduct_values = {
    "Jarosite": -0.17,
    "Iron oxide": 10,
    "Residue": -0.17,
}
byproduct_opt_conversions = {
    (3, 1): {"Jarosite": 0.75},
    (3, 2): {"Iron oxide": 1},
    (3, 3): {"Residue": 0.25},
}


class TestNPVPrintout(object):
    @pytest.fixture(scope="class")
    def model_and_results(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(m)

        results = solver.solve(m, tee=True)

        return m, results

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_results_overview(
        self, mock_print, model_and_results
    ):
        model, results = model_and_results
        report_superstructure_results_overview(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)
        assert "Superstructure Results Overview:" in printed_text
        assert "Chosen Pathway: (1, 1) -> (2, 2) -> (3, 3)" in printed_text
        assert "NPV: -16,440,490.52 USD" in printed_text

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @patch("builtins.print")
    def test_report_superstructure_costing(self, mock_print, model_and_results):
        model, results = model_and_results
        report_superstructure_costing(model, results)

        printed_output = "\n".join(call.args[0] for call in mock_print.call_args_list)

        # Check header
        assert "Superstructure Costing:" in printed_output

        # Check Net Present Value line (formatted number with comma and 2 decimals)
        assert "NPV: -16,440,490.52 USD" in printed_output

        # Check discount factor line
        assert "discount factor: 5.77%" in printed_output

        # Check Cash Flows section and example entry line
        assert "Cash Flows:" in printed_output
        assert "Year   : Cash Flow       : Units" in printed_output
        assert "2024   : -8,552.14       USD/a" in printed_output

        # Check Cash Flow Interest Rates header and example line
        assert "Cash Flow Interest Rates:" in printed_output
        assert "OPEX escalation rate   : 3.00%  : None" in printed_output
        assert "CAPEX escalation rate  : 3.60%  : None" in printed_output

        # Check Revenue section and example line
        assert "Revenue:" in printed_output
        assert (
            "Year   : Main Product Revenue : Byproduct Revenue    : Total Revenue        : Units"
            in printed_output
        )
        assert (
            "2025   : 814,912.70           : 0.00                 : 814,912.70"
            in printed_output
        )

        # Check Capital Expenses and Equipment Cost Breakdown sections
        assert "Capital Expenses:" in printed_output
        assert "Equipment Cost Breakdown:" in printed_output
        assert "Opt        : Equip. Cost  : Units" in printed_output
        assert "(2, 2)     :    16,034.64 : USD" in printed_output
        assert "(3, 3)     :     8,430.19 : USD" in printed_output
        assert "TOTAL      :    24,464.82 : USD" in printed_output
        assert "Total Plant Cost: 72,660.52 USD" in printed_output
        assert "Lang Factor: 2.97" in printed_output
        assert "Total Overnight Cost: 85,521.44 USD" in printed_output
        assert "financing factor: 2.70%" in printed_output
        assert "other costs factor: 15.00%" in printed_output

        # Check Total Overnight Cost Expended section and example line
        assert "Total Overnight Cost Expended:" in printed_output
        assert "2024       : 8,552.14                       : USD/a" in printed_output

        # Check Operating Expenses section and total operators
        assert "Operating Expenses:" in printed_output
        assert "Total Number of Operators: 9" in printed_output

        # Check Breakdown of Yearly Operating Expenses section and example line
        assert "Breakdown of Yearly Operating Expenses:" in printed_output
        assert (
            "2025       : 125,768.02                     : 4,273,448.94                   : 5,279,060.36                   : USD/a"
            in printed_output
        )

        # Check Further Breakdown of Yearly Fixed Operating Expenses section and example lines
        assert "Further Breakdown of Yearly Fixed Operating Expenses:" in printed_output
        assert (
            "Year   : Cost of Labor   : M & SM          : SA & QA/QC      : S, IP, R & D    : A & SL          : FB              : PT & I          : Units"
            in printed_output
        )
        assert (
            "2025   : 2,750,400.00    : 1,453.21        : 275,040.00      : 8,149.13        : 550,080.00      : 687,600.00      : 726.61          : USD/a"
            in printed_output
        )

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @patch("builtins.print")
    def test_report_superstructure_streams(self, mock_print, model_and_results):
        model, results = model_and_results
        report_superstructure_streams(model, results)

        printed_output = "\n".join(call.args[0] for call in mock_print.call_args_list)

        # Check main header
        assert "Superstructure Streams:" in printed_output

        # Flow Entering Stage Variables header and sample flow entries
        assert "Flow Entering Stage Variables (f):" in printed_output
        assert "Year   : Stage    : Component    : Flow" in printed_output
        # For example, check a known flow value
        assert "2025   : 1        : Nd           : 17,938.87" in printed_output
        assert "2025   : 1        : Fe           : 60,173.59" in printed_output

        # Flow In Option Variables header and sample lines
        assert "Flow In Option Variables (f_in):" in printed_output
        assert "Year   : Option       : Component    : Flow In" in printed_output
        assert "2025   : (1, 1)       : Nd           : 17,938.87" in printed_output
        assert "2025   : (1, 1)       : Fe           : 60,173.59" in printed_output

        # Flow Out Option Variables header and sample lines
        assert "Flow Out Option Variables (f_out):" in printed_output
        assert "Year   : Option       : Component    : Flow Out" in printed_output
        assert "2025   : (1, 1)       : Nd           : 17,938.87" in printed_output
        assert "2025   : (1, 1)       : Fe           : 60,173.59" in printed_output


class TestCORPrintout(object):
    @pytest.fixture(scope="class")
    def model_and_results(self):
        m = build_model(
            ### Choice of objective function
            obj_func=ObjectiveFunctionChoice.COST_OF_RECOVERY,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=False,
            byproduct_values=[],
            byproduct_opt_conversions=[],
        )

        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(m)

        results = solver.solve(m, tee=True)

        return m, results

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_results_overview(
        self, mock_print, model_and_results
    ):
        model, results = model_and_results
        report_superstructure_results_overview(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)
        assert "Cost of Recovery: 261.46 USD/kg" in printed_text

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_costing(self, mock_print, model_and_results):
        model, results = model_and_results
        report_superstructure_costing(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)
        assert "Cost of Recovery: 261.46 USD/kg" in printed_text


class TestByprodValPrintout(object):
    @pytest.fixture(scope="class")
    def model_and_results(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=False,
            options_environmental_impacts=[],
            epsilon=[],
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(m)

        results = solver.solve(m, tee=True)

        return m, results

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_streams(self, mock_print, model_and_results):
        model, results = model_and_results
        report_superstructure_streams(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)

        # Check header for byproduct section
        assert "Byproduct Produced Variables:" in printed_text

        # Check header columns in byproduct section
        assert "Year   : Byproduct       : Amount Produced    : Units" in printed_text

        # Check line separator in byproduct section
        assert (
            "------ : --------------- : ------------------ : --------" in printed_text
        )

        # Check example byproduct production lines
        assert "2025   : Residue         : 19,528.12          : kg/a" in printed_text
        assert "2026   : Residue         : 18,476.94          : kg/a" in printed_text
        assert "2027   : Residue         : 19,275.09          : kg/a" in printed_text
        assert "2028   : Residue         : 32,818.02          : kg/a" in printed_text


class TestEnvironImpactsPrintout(object):
    @pytest.fixture(scope="class")
    def model_and_results(self):
        m = build_model(
            ### Choice of objective function
            obj_func=obj_func,
            ### Plant lifetime parameters
            plant_start=plant_start,
            plant_lifetime=plant_lifetime,
            ### Feed parameters
            available_feed=available_feed,
            collection_rate=collection_rate,
            tracked_comps=tracked_comps,
            prod_comp_mass=prod_comp_mass,
            ### Superstructure formulation parameters
            num_stages=num_stages,
            options_in_stage=options_in_stage,
            option_outlets=option_outlets,
            option_efficiencies=option_efficiencies,
            ### Operating parameters
            profit=profit,
            opt_var_oc_params=opt_var_oc_params,
            operators_per_discrete_unit=operators_per_discrete_unit,
            yearly_cost_per_unit=yearly_cost_per_unit,
            capital_cost_per_unit=capital_cost_per_unit,
            processing_rate=processing_rate,
            num_operators=num_operators,
            labor_rate=labor_rate,
            ### Discretized costing parameters
            discretized_purchased_equipment_cost=discretized_purchased_equipment_cost,
            ### Environmental impacts parameters
            consider_environmental_impacts=True,
            options_environmental_impacts=options_environmental_impacts,
            epsilon=epsilon,
            ### Byproduct valorization parameters
            consider_byproduct_valorization=True,
            byproduct_values=byproduct_values,
            byproduct_opt_conversions=byproduct_opt_conversions,
        )

        # Set tolerance parameters
        solver.options["OptimalityTol"] = 1e-9  # Primal feasibility tolerance
        solver.options["FeasibilityTol"] = 1e-9  # Dual feasibility tolerance
        solver.options["NumericFocus"] = 3  # focus on getting correct solution

        # For MIP problems, you may also want:
        solver.options["MIPGap"] = 1e-9  # Relative MIP optimality gap
        solver.options["MIPGapAbs"] = 1e-9  # Absolute MIP optimality gap
        solver.options["IntFeasTol"] = 1e-9  # Integer feasibility tolerance

        # Create and apply the scaler
        scaler = SuperstructureScaler()
        scaler.scale_model(m)

        results = solver.solve(m, tee=True)

        return m, results

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_results_overview(
        self, mock_print, model_and_results
    ):
        model, results = model_and_results
        report_superstructure_results_overview(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)
        assert "Total Impacts: 10,811,781.40" in printed_text

    @pytest.mark.solver
    @pytest.mark.skipif(not solver_available, reason="Gurobi solver not available")
    @pytest.mark.component
    @patch("builtins.print")
    def test_report_superstructure_environmental_impacts(
        self, mock_print, model_and_results
    ):
        model, results = model_and_results
        report_superstructure_environmental_impacts(model, results)

        printed_text = "".join(call.args[0] for call in mock_print.call_args_list)

        assert "Superstructure Environmental Impacts:" in printed_text
        assert "Total Impacts: 10,811,781.40" in printed_text
        assert "Epsilon factor: 10,000,000,000,000,000.00" in printed_text

        assert "Impacts broken down by year:" in printed_text
        assert "Year   : Total Environmental Impact : Units" in printed_text
        assert "2025   : 2,343,373.9290            : 1/a" in printed_text
        assert "2026   : 2,217,233.3040            : 1/a" in printed_text
        assert "2027   : 2,313,011.3760            : 1/a" in printed_text
        assert "2028   : 3,938,162.7870            : 1/a" in printed_text

        assert "Impacts broken down by year and option:" in printed_text
        assert "Year   : Option       : Environmental Impact : Units" in printed_text
        # Check some known lines for different years and options
        assert "2025   : (1, 1)       : 781,124.6430         : 1/a" in printed_text
        assert "2025   : (2, 2)       : 781,124.6430         : 1/a" in printed_text
        assert "2025   : (3, 3)       : 781,124.6430         : 1/a" in printed_text
        assert "2026   : (1, 1)       : 739,077.7680         : 1/a" in printed_text
        assert "2026   : (2, 2)       : 739,077.7680         : 1/a" in printed_text
        assert "2026   : (3, 3)       : 739,077.7680         : 1/a" in printed_text
        assert "2027   : (1, 1)       : 771,003.7920         : 1/a" in printed_text
        assert "2027   : (2, 2)       : 771,003.7920         : 1/a" in printed_text
        assert "2027   : (3, 3)       : 771,003.7920         : 1/a" in printed_text
        assert "2028   : (1, 1)       : 1,312,720.9290       : 1/a" in printed_text
        assert "2028   : (2, 2)       : 1,312,720.9290       : 1/a" in printed_text
        assert "2028   : (3, 3)       : 1,312,720.9290       : 1/a" in printed_text
