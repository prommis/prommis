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

The steps to test this framework including:
1. Import packages
2. Build the li-co diafiltration recovery flowsheet
3. add costing
4. Solve and display results
5. Optimize stage lengths to minimize cost of recovery and display results


"""
from pyomo.environ import (
    Objective,
    TransformationFactory,
    SolverFactory,
    value,
    minimize,
    units as pyunits,
    Constraint,
    Expression,
    Param,
    Var,
    Suffix,
    )

from idaes.core import (
    FlowsheetBlock,
    UnitModelBlock,
    UnitModelCostingBlock,
)
from idaes.core.util.model_diagnostics import DiagnosticsToolbox, degrees_of_freedom
from idaes.core.util.scaling import constraint_autoscale_large_jac

from prommis.nanofiltration.diafiltration import (
    add_costing,
    add_product_constraints,
    build_model,
    initialize_model,
    set_scaling,
    solve_model,
    unfix_opt_variables,
)

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCosting,
    DiafiltrationCostingData,
)

from prommis.uky.costing.ree_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)

from prommis.uky.costing.diafiltration_flowsheet_optimization_example import (
    build_costing,
    build_optimization,
    scale_and_solve_model,
)

import pytest

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
        # TODO address unbounded solute sieving expressions inside log functions
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
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(750, rel=1e-4)
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(150, rel=1e-4)
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(750, rel=1e-4)
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(150, rel=1e-4)
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(750, rel=1e-4)
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(150, rel=1e-4)
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(114445, rel=1e-4)
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(42145, rel=1e-4)
        assert value(model.fs.feed_pump.costing.variable_operating_cost) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(26859, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.variable_operating_cost) == pytest.approx(15468, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.071253, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(0.14251, rel=1e-4)
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(0.0047252, rel=1e-4)
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(0.15684, rel=1e-4)
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(11.56977, rel=1e-4)

    @pytest.mark.component
    def test_presolve_LiCoDiafiltration_model(self, model):
        # fix some initial long stage lengths so optimization starts from a feasible place with recovery bounds
        model.fs.stage1.length.fix(5000)
        model.fs.stage2.length.fix(1000)
        model.fs.stage3.length.fix(750)

        solve_model(model, tee=False)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # results in USD_2021
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(375000, rel=1e-4)
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(75000, rel=1e-4)
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(75000, rel=1e-4)
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(15000, rel=1e-4)
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(56250, rel=1e-4)
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(11250, rel=1e-4)
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(114445, rel=1e-4)
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(42145, rel=1e-4)
        assert value(model.fs.feed_pump.costing.variable_operating_cost) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(53003, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.variable_operating_cost) == pytest.approx(88391, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.60140, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(1.2028, rel=1e-4)
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(0.13733, rel=1e-4)
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(0.27087, rel=1e-4)
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(0.12937, rel=1e-4)

    @pytest.mark.component
    def test_optimize_LiCoDiafiltration_model(self, model):

        build_optimization(model)
    
        scale_and_solve_model(model)
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        # decision variables
        assert value(model.fs.stage1.length) == pytest.approx(1178.7, rel=1e-4)
        assert value(model.fs.stage2.length) == pytest.approx(588.18, rel=1e-4)
        assert value(model.fs.stage3.length) == pytest.approx(765.39, rel=1e-4)

        # results in USD_2021
        assert value(model.fs.stage1.costing.capital_cost) == pytest.approx(88406, rel=1e-4)
        assert value(model.fs.stage1.costing.fixed_operating_cost) == pytest.approx(17681, rel=1e-4)
        assert value(model.fs.stage2.costing.capital_cost) == pytest.approx(44114, rel=1e-4)
        assert value(model.fs.stage2.costing.fixed_operating_cost) == pytest.approx(8823, rel=1e-4)
        assert value(model.fs.stage3.costing.capital_cost) == pytest.approx(57404, rel=1e-4)
        assert value(model.fs.stage3.costing.fixed_operating_cost) == pytest.approx(11481, rel=1e-4)
        assert value(model.fs.cascade.costing.variable_operating_cost) == pytest.approx(114445, rel=1e-4)
        assert value(model.fs.feed_pump.costing.capital_cost) == pytest.approx(42145, rel=1e-4)
        assert value(model.fs.feed_pump.costing.variable_operating_cost) == pytest.approx(0.00338663, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.capital_cost) == pytest.approx(44989, rel=1e-4)
        assert value(model.fs.diafiltrate_pump.costing.variable_operating_cost) == pytest.approx(58057, rel=1e-4)

        # results in MUSD_2021
        assert value(model.fs.costing.total_BEC) == pytest.approx(0.27706, rel=1e-4)
        assert value(model.fs.costing.total_plant_cost) == pytest.approx(0.55411, rel=1e-4)
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(0.054608, rel=1e-4)
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(0.21792, rel=1e-4)
        assert value(model.fs.costing.cost_of_recovery) == pytest.approx(0.052949, rel=1e-4)
