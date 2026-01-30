#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from copy import deepcopy

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    SolverFactory,
    value,
    Var,
    units
)
from pyomo.common.fileutils import this_file_dir

from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.solvers import get_solver
from idaes.core.scaling.util import jacobian_cond
from idaes.core.util import DiagnosticsToolbox, from_json, StoreSpec

import pytest

from prommis.util import assert_solution_equivalent, copy_first_steady_state
from prommis.leaching.leach_flowsheet_merged import CocurrentSlurryLeachingFlowsheet

_expected_results_steady_one_tank = {}
_rel = 1e-5
_abs = 1e-6

_expected_results_steady_one_tank["recovery"] = {
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
_expected_results_steady_one_tank["fs.leach.mscontactor.solid[0, 1].conversion_comp"] = {
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

_expected_results_steady_two_tanks = {}

_expected_results_steady_two_tanks[
    "fs.leach.mscontactor.solid[0, 1].conversion_comp"
] = _expected_results_steady_one_tank[
    "fs.leach.mscontactor.solid[0, 1].conversion_comp"
]

_expected_results_steady_two_tanks["fs.leach.mscontactor.solid[0, 2].conversion_comp"] = {
    "Al2O3": (3.117560e-02, _rel, _abs),
    "CaO": (7.481940e-01, _rel, _abs),
    "Ce2O3": (5.202264e-01, _rel, _abs),
    "Dy2O3": (3.049849e-01, _rel, _abs),
    "Fe2O3": (1.856473e-01, _rel, _abs),
    "Gd2O3": (6.646174e-01, _rel, _abs),
    "La2O3": (5.083381e-01, _rel, _abs),
    "Nd2O3": (5.454125e-01, _rel, _abs),
    "Pr2O3": (6.166401e-01, _rel, _abs),
    "Sc2O3": (6.510993e-02, _rel, _abs),
    "Sm2O3": (2.960138e-01, _rel, _abs),
    "Y2O3": (1.918858e-01, _rel, _abs),
    "inerts": (0.0, _rel, _abs),
}

_expected_results_steady_two_tanks["recovery"] = {
    "inerts": (0, _rel, _abs),
    "Sc2O3": (6.510992e+00, _rel, _abs),
    "Y2O3": (1.918858e+01, _rel, _abs),
    "La2O3": (5.083381e+01, _rel, _abs),
    "Ce2O3": (5.202264e+01, _rel, _abs),
    "Pr2O3": (6.166401e+01, _rel, _abs),
    "Nd2O3": (5.454125e+01, _rel, _abs),
    "Sm2O3": (2.960138e+01, _rel, _abs),
    "Gd2O3": (6.646174e+01, _rel, _abs),
    "Dy2O3": (3.049849e+01, _rel, _abs),
    "Al2O3": (3.117559e+00, _rel, _abs),
    "CaO": (7.481940e+01, _rel, _abs),
    "Fe2O3": (1.856473e+01, _rel, _abs),
}

class TestSteadyStateModel:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet()
        m.fs.scale_model()

        return m


    @pytest.mark.unit
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()


    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):

        initializer = model.fs.default_initializer()
        initializer.initialize(model.fs)

        solver = get_solver("ipopt_v2")
        results = solver.solve(model, tee=False)

        assert_optimal_termination(results)


    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            6.234582e12, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(
            3827.42, rel=1e-3
        )


    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        m = model

        @m.Expression(m.fs.coal.component_list)
        def recovery(b, j):
            f_in = b.fs.leach.solid_inlet.flow_mass[0]
            f_out = b.fs.leach.solid_outlet.flow_mass[0]
            x_in = b.fs.leach.solid_inlet.mass_frac_comp[0, j]
            x_out = b.fs.leach.solid_outlet.mass_frac_comp[0, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

        assert_solution_equivalent(m, expected_results=_expected_results_steady_one_tank)

class TestSteadyStateModelTwoTanks:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet(
            number_of_tanks = 2
        )
        m.fs.scale_model()

        return m


    @pytest.mark.unit
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()


    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):

        initializer = model.fs.default_initializer()
        initializer.initialize(model.fs)

        solver = get_solver("ipopt_v2")
        results = solver.solve(model, tee=False)

        assert_optimal_termination(results)


    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            1.14913e13, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(
            9657.39, rel=1e-3
        )


    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        m = model

        @m.Expression(m.fs.coal.component_list)
        def recovery(b, j):
            f_in = b.fs.leach.solid_inlet.flow_mass[0]
            f_out = b.fs.leach.solid_outlet.flow_mass[0]
            x_in = b.fs.leach.solid_inlet.mass_frac_comp[0, j]
            x_out = b.fs.leach.solid_outlet.mass_frac_comp[0, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

        assert_solution_equivalent(m, expected_results=_expected_results_steady_two_tanks)

class TestSteadyStateModelWithHoldup:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet(has_holdup=True)
        m.fs.scale_model()

        return m


    @pytest.mark.unit
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()


    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):

        initializer = model.fs.default_initializer()
        initializer.initialize(model.fs)

        solver = get_solver("ipopt_v2")
        results = solver.solve(model, tee=False)

        assert_optimal_termination(results)


    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            6.243792e12, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(59103.1, rel=1e-3)


    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        m = model

        @m.Expression(m.fs.coal.component_list)
        def recovery(b, j):
            f_in = b.fs.leach.solid_inlet.flow_mass[0]
            f_out = b.fs.leach.solid_outlet.flow_mass[0]
            x_in = b.fs.leach.solid_inlet.mass_frac_comp[0, j]
            x_out = b.fs.leach.solid_outlet.mass_frac_comp[0, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

        assert_solution_equivalent(m, expected_results=_expected_results_steady_one_tank)


class TestDynamicOneTank:
    @pytest.fixture(scope="class")
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

        # Set initial value for accumulation terms.
        for var in m.component_data_objects(ctype=Var, descend_into=True):
            if var.value is None:
                var.set_value(0)

        for t in m.fs.time:
            if t <= perturb_time:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(224.3 * units.L / units.hour)
            else:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(300 * units.L / units.hour)
            

        return m

    @pytest.mark.unit
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings(ignore_unit_consistency=True)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):
        solver = get_solver("ipopt_v2")
        results = solver.solve(model, tee=False)

        assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        # The solution takes chloride to it's lower bound, which
        # the diagnostics toolbox does not like. This isn't actually
        # a numerical issue, though, so just fudge the value away from 1e-20.
        for t in model.fs.time:
            for e in model.fs.leach.mscontactor.elements:
                if (
                model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value
                == model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].lb
                ):
                    model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value = (
                        1.01 * model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value
                    )

        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            3.4449e14, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(9.8495e6, rel=1e-3)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        # These results come from integrating the problem with
        # PETSc, so we expect only loose agreement.
        _abs = None
        _rel = 1e-2
        expected_results = {
            "fs.leach.mscontactor.liquid[24.0,1].conc_mass_comp":{
                "Sc": (8.0802607e-02, _rel, _abs),
                "Y": (2.9420954e-01, _rel, _abs),
                "La": (1.5865103e+00, _rel, _abs),
                "Ce": (4.2327277e+00, _rel, _abs),
                "Pr": (5.1833345e-01, _rel, _abs),
                "Nd": (2.0665492e+00, _rel, _abs),
                "Sm": (2.2415357e-01, _rel, _abs),
                "Gd": (4.6155189e-01, _rel, _abs),
                "Dy": (1.0447770e-01, _rel, _abs),
            },
            "fs.leach.mscontactor.solid[24.0,1].conversion_comp":{
                "Sc2O3": (5.8234784e-02, _rel, _abs),
                "Y2O3": (1.5362446e-01, _rel, _abs),
                "La2O3": (3.7899544e-01, _rel, _abs),
                "Ce2O3": (4.2832328e-01, _rel, _abs),
                "Pr2O3": (4.8365811e-01, _rel, _abs),
                "Nd2O3": (4.7454274e-01, _rel, _abs),
                "Sm2O3": (2.3762900e-01, _rel, _abs),
                "Gd2O3": (6.6387387e-01, _rel, _abs),
                "Dy2O3": (2.2042091e-01, _rel, _abs),
            }
        }
        assert_solution_equivalent(model, expected_results)

class TestDynamicTwoTanks:
    @pytest.fixture(scope="class")
    def model(self):
        perturb_time = 4

        m = ConcreteModel()
        m.fs = CocurrentSlurryLeachingFlowsheet(
            has_holdup=True,
            dynamic=True,
            time_set=range(25),
            time_units=units.hr,
            number_of_tanks=2,
        )   
        
        m.fs.scale_model()
        from_json(m, fname="leaching_two_tanks.json", wts=StoreSpec.value())

        m.fs.reduce_dae_index()
        m.fs.scale_dynamics(15)
        m.fs.fix_initial_conditions()
        copy_first_steady_state(m, m.fs.time)

        # Set initial value for accumulation terms.
        for var in m.component_data_objects(ctype=Var, descend_into=True):
            if var.value is None:
                var.set_value(0)

        for t in m.fs.time:
            if t <= perturb_time:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(224.3 * units.L / units.hour)
            else:
                m.fs.leach.liquid_inlet.flow_vol[t].fix(300 * units.L / units.hour)
            

        return m

    @pytest.mark.unit
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings(ignore_unit_consistency=True)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):
        solver = get_solver("ipopt_v2")
        results = solver.solve(model, tee=False)

        assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, model):
        # The solution takes chloride to it's lower bound, which
        # the diagnostics toolbox does not like. This isn't actually
        # a numerical issue, though, so just fudge the value away from 1e-20.
        for t in model.fs.time:
            for e in model.fs.leach.mscontactor.elements:
                if (
                model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value
                == model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].lb
                ):
                    model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value = (
                        1.01 * model.fs.leach.mscontactor.liquid[t, e].conc_mol_comp["Cl"].value
                    )

        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

        assert jacobian_cond(model, scaled=False) == pytest.approx(
            3.060857e15, rel=1e-3
        )
        assert jacobian_cond(model, scaled=True) == pytest.approx(3.00789e7, rel=1e-3)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        # These results come from integrating the problem with
        # PETSc, so we expect only loose agreement.
        _abs = None
        _rel = 1e-2
        expected_results = {
            "fs.leach.mscontactor.liquid[24.0,1].conc_mass_comp":{
                "Sc": (8.0802607e-02, _rel, _abs),
                "Y": (2.9420954e-01, _rel, _abs),
                "La": (1.5865103e+00, _rel, _abs),
                "Ce": (4.2327277e+00, _rel, _abs),
                "Pr": (5.1833345e-01, _rel, _abs),
                "Nd": (2.0665492e+00, _rel, _abs),
                "Sm": (2.2415357e-01, _rel, _abs),
                "Gd": (4.6155189e-01, _rel, _abs),
                "Dy": (1.0447770e-01, _rel, _abs),
            },
            "fs.leach.mscontactor.solid[24.0,1].conversion_comp":{
                "Sc2O3": (5.8234784e-02, _rel, _abs),
                "Y2O3": (1.5362446e-01, _rel, _abs),
                "La2O3": (3.7899544e-01, _rel, _abs),
                "Ce2O3": (4.2832328e-01, _rel, _abs),
                "Pr2O3": (4.8365811e-01, _rel, _abs),
                "Nd2O3": (4.7454274e-01, _rel, _abs),
                "Sm2O3": (2.3762900e-01, _rel, _abs),
                "Gd2O3": (6.6387387e-01, _rel, _abs),
                "Dy2O3": (2.2042091e-01, _rel, _abs),
            },
            "fs.leach.mscontactor.liquid[24.0,2].conc_mass_comp":{
                "Sc": (1.0301485e-01, _rel, _abs),
                "Y": (4.0293409e-01, _rel, _abs),
                "La": (2.2431474e+00, _rel, _abs),
                "Ce": (5.5868604e+00, _rel, _abs),
                "Pr": (7.0096265e-01, _rel, _abs),
                "Nd": (2.6234931e+00, _rel, _abs),
                "Sm": (3.0508031e-01, _rel, _abs),
                "Gd": (5.2969391e-01, _rel, _abs),
                "Dy": (1.5290345e-01, _rel, _abs),
            },
            "fs.leach.mscontactor.solid[24.0,2].conversion_comp":{
                "Sc2O3": (6.9370858e-02, _rel, _abs),
                "Y2O3": (1.9983959e-01, _rel, _abs),
                "La2O3": (5.1766375e-01, _rel, _abs),
                "Ce2O3": (5.4037591e-01, _rel, _abs),
                "Pr2O3": (6.3170014e-01, _rel, _abs),
                "Nd2O3": (5.7224648e-01, _rel, _abs),
                "Sm2O3": (3.0782506e-01, _rel, _abs),
                "Gd2O3": (7.1232093e-01, _rel, _abs),
                "Dy2O3": (3.1047217e-01, _rel, _abs),
            }
        }
        assert_solution_equivalent(model, expected_results)