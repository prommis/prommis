#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Author: Lingyan Deng, Brandon Paul

This is a test file for testing flowsheet optimization with PrOMMiS REE Costing. The example is a Li-Co diafiltration process to recover Li and Co as byproduct from battery recovery. Li-Co recovery
process is referenced from literature: Wamble, N.P., Eugene, E.A., Phillip, W.A., Dowling, A.W., 'Optimal Diafiltration Membrane Cascades Enable Green Recycling of Spent Lithium-Ion 
Batteries', ACS Sustainable Chem. Eng. 2022, 10, 12207−12225. 

The steps to test this framework include:
1. Import packages
2. Build the Li-Co diafiltration recovery flowsheet
3. Add costing
4. Solve and display results
5. Optimize stage lengths to minimize cost of recovery and display results


"""
from pyomo.environ import value

from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox, degrees_of_freedom

import pytest

from prommis.nanofiltration.diafiltration import (
    build_model,
    initialize_model,
    solve_model,
)
from prommis.uky.costing.diafiltration_flowsheet_optimization_example import (
    build_costing,
    build_optimization,
    print_io_snap,
    print_stage_cuts,
    scale_and_solve_model,
    main,
)
from prommis.uky.costing.ree_plant_capcost import QGESSCosting


class TestDiafiltrationOptimization:
    @pytest.fixture(scope="class")
    def model(self):
        model = build_model()

        return model

    @pytest.mark.component
    def test_build_LiCoDiafiltration_model(self, model):

        model = build_model()

        assert degrees_of_freedom(model) == 0
        dt = DiagnosticsToolbox(model)
        dt.report_structural_issues()
        dt.display_potential_evaluation_errors()
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.component
    def test_initialize_LiCoDiafiltration_model(self, model):
        initialize_model(model)

    @pytest.mark.component
    def test_solve_LiCoDiafiltration_model(self, model):
        solve_model(model, tee=False)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_cost_LiCoDiafiltration_model(self, model):
        build_costing(model)

        assert hasattr(model.fs, "costing")
        assert isinstance(model.fs.costing, QGESSCosting)
        assert hasattr(model.fs.stage1, "costing")
        assert isinstance(model.fs.stage1.costing, UnitModelCostingBlock)
        assert hasattr(model.fs.stage2, "costing")
        assert isinstance(model.fs.stage2.costing, UnitModelCostingBlock)
        assert hasattr(model.fs.stage3, "costing")
        assert isinstance(model.fs.stage3.costing, UnitModelCostingBlock)
        assert hasattr(model.fs.cascade, "costing")
        assert isinstance(model.fs.cascade.costing, UnitModelCostingBlock)
        assert hasattr(model.fs.feed_pump, "costing")
        assert isinstance(model.fs.feed_pump.costing, UnitModelCostingBlock)
        assert hasattr(model.fs.diafiltrate_pump, "costing")
        assert isinstance(model.fs.diafiltrate_pump.costing, UnitModelCostingBlock)

    @pytest.mark.component
    def test_solve_cost_LiCoDiafiltration_model(self, model):
        solve_model(model, tee=False)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # results in USD_2021
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(
            750, rel=1e-4
        )
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(
            150, rel=1e-4
        )
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(
            750, rel=1e-4
        )
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(
            150, rel=1e-4
        )
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(
            750, rel=1e-4
        )
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(
            150, rel=1e-4
        )
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(
            114445, rel=1e-4
        )
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(
            42145, rel=1e-4
        )
        assert value(
            model.fs.feed_pump.costing.variable_operating_cost
        ) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(
            26859, rel=1e-4
        )
        assert value(
            model.fs.diafiltrate_pump.costing.variable_operating_cost
        ) == pytest.approx(15468, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.071253, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            0.14251, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            0.0047252, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            0.15684, rel=1e-4
        )
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(
            0.0129372, rel=1e-4
        )

    @pytest.mark.component
    def test_presolve_LiCoDiafiltration_model(self, model):
        # fix some initial long stage lengths so optimization starts from a feasible place with recovery bounds
        model.fs.stage1.length.fix(754)
        model.fs.stage2.length.fix(758)
        model.fs.stage3.length.fix(756)

        solve_model(model, tee=False)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # results in USD_2021
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(
            56550, rel=1e-4
        )
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(
            11310, rel=1e-4
        )
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(
            56850, rel=1e-4
        )
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(
            11370, rel=1e-4
        )
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(
            56700, rel=1e-4
        )
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(
            11340, rel=1e-4
        )
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(
            114445, rel=1e-4
        )
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(
            42145, rel=1e-4
        )
        assert value(
            model.fs.feed_pump.costing.variable_operating_cost
        ) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(
            48545.7, rel=1e-4
        )
        assert value(
            model.fs.diafiltrate_pump.costing.variable_operating_cost
        ) == pytest.approx(70565.5, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.26079, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            0.521581, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            0.0496674, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            0.231946, rel=1e-4
        )
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(
            0.033803, rel=1e-4
        )

    @pytest.mark.component
    def test_optimize_LiCoDiafiltration_model(self, model):

        build_optimization(model)

        scale_and_solve_model(model)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # decision variables
        assert value(model.fs.stage1.length) == pytest.approx(1195.6, rel=1e-4)
        assert value(model.fs.stage2.length) == pytest.approx(584.4, rel=1e-4)
        assert value(model.fs.stage3.length) == pytest.approx(757.8, rel=1e-4)

        # results in USD_2021
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(
            89667.74, rel=1e-4
        )
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(
            17933.55, rel=1e-4
        )
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(
            43833.29, rel=1e-4
        )
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(
            8766.66, rel=1e-4
        )
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(
            56834.03, rel=1e-4
        )
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(
            11366.81, rel=1e-4
        )
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(
            114445, rel=1e-4
        )
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(
            42144.7, rel=1e-4
        )
        assert value(
            model.fs.feed_pump.costing.variable_operating_cost
        ) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(
            44905.28, rel=1e-4
        )
        assert value(
            model.fs.diafiltrate_pump.costing.variable_operating_cost
        ) == pytest.approx(57781.50, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.27739, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            0.55477, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            0.05471, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            0.21761, rel=1e-4
        )
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(
            0.032835, rel=1e-4
        )

    @pytest.mark.component
    def test_print_methods(self, model):
        print_io_snap(model.fs, tag="AFTER OPTIMIZATION")
        print_stage_cuts(model, label="STAGE CUTS — AFTER OPTIMIZATION")

    @pytest.mark.component
    def test_main(self, model):
        main()
