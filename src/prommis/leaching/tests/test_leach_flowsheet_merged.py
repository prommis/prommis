#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    SolverFactory,
    value,
    units
)

from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.solvers import get_solver
from idaes.core.scaling.util import jacobian_cond
from idaes.core.util import DiagnosticsToolbox, from_json, StoreSpec

import pytest

from prommis.util import assert_solution_equivalent, copy_first_steady_state
from prommis.leaching.leach_flowsheet_merged import CocurrentSlurryLeachingFlowsheet

_expected_results_steady = {}
_rel = 1e-5
_abs = 1e-6

_expected_results_steady["recovery"] = {
    "inerts": (3.3306690738754696e-14, _rel, _abs),
    "Sc2O3": (5.486783341761903, _rel, _abs),
    "Y2O3": (15.092538416882373, _rel, _abs),
    "La2O3": (38.342494088335066, _rel, _abs),
    "Ce2O3": (42.074359313539645, _rel, _abs),
    "Pr2O3": (48.42509098259142, _rel, _abs),
    "Nd2O3": (45.79540997369025, _rel, _abs),
    "Sm2O3": (23.377407833768427, _rel, _abs),
    "Gd2O3": (61.06358012227954, _rel, _abs),
    "Dy2O3": (22.365403062815247, _rel, _abs),
    "Al2O3": (3.0037655469309033, _rel, _abs),
    "CaO": (51.45052459223044, _rel, _abs),
    "Fe2O3": (16.33850709753919, _rel, _abs),
}
_expected_results_steady["fs.leach.mscontactor.solid[0, 1].conversion_comp"] = {
    "Al2O3": (0.03003765546930494, _rel, _abs),
    "CaO": (0.5145052459223426, _rel, _abs),
    "Ce2O3": (0.4207435931354908, _rel, _abs),
    "Dy2O3": (0.22365403062817607, _rel, _abs),
    "Fe2O3": (0.16338507097541013, _rel, _abs),
    "Gd2O3": (0.6106358012228825, _rel, _abs),
    "La2O3": (0.38342494088341084, _rel, _abs),
    "Nd2O3": (0.45795409973702533, _rel, _abs),
    "Pr2O3": (0.4842509098265882, _rel, _abs),
    "Sc2O3": (0.05486783341762388, _rel, _abs),
    "Sm2O3": (0.23377407833771635, _rel, _abs),
    "Y2O3": (0.15092538416884166, _rel, _abs),
    "inerts": (0.0, _rel, _abs),
}

class TestSteadyStateModel:
    @pytest.fixture(scope="class")
    def model_steady_state(self):
        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet()
        m.fs.scale_model()

        return m


    @pytest.mark.unit
    def test_structural_issues(self, model_steady_state):
        dt = DiagnosticsToolbox(model_steady_state)
        dt.assert_no_structural_warnings()


    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model_steady_state):

        initializer = model_steady_state.fs.default_initializer()
        initializer.initialize(model_steady_state.fs)

        # Solve scaled model
        solver = get_solver("ipopt_v2")
        results = solver.solve(model_steady_state, tee=False)

        assert_optimal_termination(results)


    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model_steady_state):
        dt = DiagnosticsToolbox(model_steady_state)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model_steady_state, scaled=False) == pytest.approx(
            6.234582e12, rel=1e-3
        )
        assert jacobian_cond(model_steady_state, scaled=True) == pytest.approx(
            3827.42, rel=1e-3
        )


    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model_steady_state):
        m = model_steady_state

        @m.Expression(m.fs.coal.component_list)
        def recovery(b, j):
            f_in = b.fs.leach.solid_inlet.flow_mass[0]
            f_out = b.fs.leach.solid_outlet.flow_mass[0]
            x_in = b.fs.leach.solid_inlet.mass_frac_comp[0, j]
            x_out = b.fs.leach.solid_outlet.mass_frac_comp[0, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

        assert_solution_equivalent(m, expected_results=_expected_results_steady)

class TestSteadyStateModelWithHoldup:
    @pytest.fixture(scope="class")
    def model_holdup(self):
        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet(has_holdup=True)
        m.fs.scale_model()

        return m


    @pytest.mark.unit
    def test_structural_issues(self, model_holdup):
        dt = DiagnosticsToolbox(model_holdup)
        dt.assert_no_structural_warnings()


    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model_holdup):

        initializer = model_holdup.fs.default_initializer()
        initializer.initialize(model_holdup.fs)

        # Solve scaled model
        solver = get_solver("ipopt_v2")
        results = solver.solve(model_holdup, tee=False)

        assert_optimal_termination(results)


    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model_holdup):
        dt = DiagnosticsToolbox(model_holdup)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model_holdup, scaled=False) == pytest.approx(
            6.243792e12, rel=1e-3
        )
        assert jacobian_cond(model_holdup, scaled=True) == pytest.approx(59103.1, rel=1e-3)


    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model_holdup):
        m = model_holdup

        @m.Expression(m.fs.coal.component_list)
        def recovery(b, j):
            f_in = b.fs.leach.solid_inlet.flow_mass[0]
            f_out = b.fs.leach.solid_outlet.flow_mass[0]
            x_in = b.fs.leach.solid_inlet.mass_frac_comp[0, j]
            x_out = b.fs.leach.solid_outlet.mass_frac_comp[0, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

        assert_solution_equivalent(m, expected_results=_expected_results_steady)


class TestDynamicOneTank:
    @pytest.mark.fixture(scope="class")
    def model(self):
        perturb_time = 4

        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet(
            has_holdup=True,
            dynamic=True,
            time_set=range(25),
            time_units=units.hr,
        )
        
        m.fs.scale_model()
        from_json(m, fname="leaching_one_tank.json", wts=StoreSpec.value())

        m.fs.reduce_dae_index()
        m.fs.scale_dynamics(15)
        m.fs.fix_initial_conditions()
        copy_first_steady_state(m, m.fs.time)

        for t in m.fs.time:
            if t <= perturb_time:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(224.3 * units.L / units.hour)
            else:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(300 * units.L / units.hour)
            

        return m
    
    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):

        initializer = BlockTriangularizationInitializer()
        initializer.initialize(model.fs)

        # Solve scaled model
        solver = get_solver("ipopt_v2")
        results = solver.solve(model_holdup, tee=False)

        assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        expected_results = {
            "fs.leach.mscontactor.liquid[24.0,1].conc_mass_comp":{
                "Sc": (0.08080260729256393, None, None)
            }
        }