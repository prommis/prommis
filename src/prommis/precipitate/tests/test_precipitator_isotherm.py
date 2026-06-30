#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import pytest
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    SolverFactory,
    TransformationFactory,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import InitializationStatus
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.testing import assert_solution_equivalent
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import unscaled_variables_generator

# from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
# from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
# from prommis.precipitate.precipitate_reactions import OxalatePrecipitationReactions
# from prommis.precipitate.precipitator_isotherm import (
#     OxalatePrecipitator,
#     OxalatePrecipitatorInitializer,
# )

from precipitate_liquid_properties_isotherm import AqueousParameter
from precipitate_solids_properties_isotherm import PrecipitateParameters
from precipitate_reactions import OxalatePrecipitationReactions
from precipitator_isotherm import (
    OxalatePrecipitator,
    OxalatePrecipitatorInitializer,
)

# from prommis.util import assert_solution_equivalent
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_aq = AqueousParameter()
    m.fs.properties_solid = PrecipitateParameters()
    m.fs.prec_rxns = OxalatePrecipitationReactions()

    m.fs.unit = OxalatePrecipitator(
        number_of_tanks=1,
        liquid_phase={
            "property_package": m.fs.properties_aq,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        solid_phase={
            "property_package": m.fs.properties_solid,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        reaction_package=m.fs.prec_rxns,
    )

    assert len(m.fs.unit.config) == 7

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.liquid_phase.property_package is m.fs.properties_aq
    assert m.fs.unit.config.solid_phase.property_package is m.fs.properties_solid
    assert m.fs.unit.config.reaction_package is m.fs.prec_rxns


# -----------------------------------------------------------------------------
class TestPrec(object):
    @pytest.fixture(scope="class")
    def prec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties_aq = AqueousParameter()
        m.fs.properties_solid = PrecipitateParameters()
        m.fs.prec_rxns = OxalatePrecipitationReactions()

        m.fs.unit = OxalatePrecipitator(
            number_of_tanks=1,
            liquid_phase={
                "property_package": m.fs.properties_aq,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            solid_phase={
                "property_package": m.fs.properties_solid,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            reaction_package=m.fs.prec_rxns,
        )

        m.fs.unit.aqueous_inlet.flow_vol[0].fix(100)

        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Al"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Y"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "La"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2C2O4"].fix(6400)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(100000)

        m.fs.unit.aqueous_inlet.temperature[0].fix(298.15)
        m.fs.unit.aqueous_inlet.pressure[0].fix(1e5)

        m.fs.unit.precipitate_outlet.temperature.fix(348.15)
        m.fs.unit.hydraulic_retention_time[0.0].fix(2)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, prec):
        assert hasattr(prec.fs.unit, "aqueous_inlet")
        assert len(prec.fs.unit.aqueous_inlet.vars) == 4
        assert hasattr(prec.fs.unit.aqueous_inlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_inlet, "conc_mass_comp")
        assert hasattr(prec.fs.unit.aqueous_inlet, "temperature")
        assert hasattr(prec.fs.unit.aqueous_inlet, "pressure")

        assert hasattr(prec.fs.unit, "aqueous_outlet")
        assert len(prec.fs.unit.aqueous_outlet.vars) == 4
        assert hasattr(prec.fs.unit.aqueous_outlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_outlet, "conc_mass_comp")
        assert hasattr(prec.fs.unit.aqueous_outlet, "temperature")
        assert hasattr(prec.fs.unit.aqueous_outlet, "pressure")

        assert hasattr(prec.fs.unit, "precipitate_outlet")
        assert len(prec.fs.unit.precipitate_outlet.vars) == 2
        assert hasattr(prec.fs.unit.precipitate_outlet, "flow_mol_comp")
        assert hasattr(prec.fs.unit.precipitate_outlet, "temperature")

        assert hasattr(prec.fs.unit, "init_solid_constraint")
        assert hasattr(prec.fs.unit, "temp_constraint")
        assert hasattr(prec.fs.unit, "heterogeneous_reaction_extent_constraint")
        assert hasattr(prec.fs.unit, "eq_hydraulic_retention")

        assert number_variables(prec.fs.unit) == 208
        assert number_total_constraints(prec.fs.unit) == 185
        assert number_unused_variables(prec.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, prec):
        assert_units_consistent(prec.fs.unit)

        dt = DiagnosticsToolbox(model=prec)

    @pytest.mark.unit
    def test_dof(self, prec):
        assert degrees_of_freedom(prec) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, prec):
        # scaling = TransformationFactory("core.scale_model")
        # scaled_model = scaling.create_using(prec, rename=False)
        # initializer = scaled_model.fs.unit.default_initializer()
        # initializer.initialize(scaled_model.fs.unit)

        # # Solve scaled model
        # solver = get_solver("ipopt_v2")
        # results = solver.solve(scaled_model, tee=False)

        # # Propagate results back to unscaled model
        # scaling.propagate_solution(scaled_model, prec)

        # assert_optimal_termination(results)

        initializer = OxalatePrecipitatorInitializer()
        initializer.initialize(prec.fs.unit)
        assert initializer.summary[prec.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.component
    def test_var_scaling(self, prec):
        unscaled_var_list = list(
            unscaled_variables_generator(prec.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, prec):
        dt = DiagnosticsToolbox(prec)
        dt.assert_no_numerical_warnings()
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        expected_results = {
            "unit.aqueous_outlet.flow_vol": {
                0: (100, None, 1e-0),
            },
            "unit.aqueous_outlet.conc_mass_comp": {
                (0, "Al"): (9.982, None, 1e-3),
                (0, "Ca"): (10, None, 1e-3),
                (0, "Ce"): (0.0991, None, 1e-3),
                (0, "Dy"): (2.3824, None, 1e-3),
                (0, "Fe"): (9.7864, None, 1e-3),
                (0, "Gd"): (0.4486, None, 1e-3),
                (0, "La"): (1.5715, 1e-3, None),
                (0, "Nd"): (0.1142, None, 1e-3),
                (0, "Pr"): (0.2118, None, 1e-3),
                (0, "Sc"): (6.3964, None, 1e-3),
                (0, "Sm"): (0.2183, None, 1e-3),
                (0, "Y"): (1.8401, None, 1e-3),
            },
            "unit.precipitate_outlet.flow_mol_comp": {
                (0, "Al2(C2O4)3(s)"): (3.20213e-05, None, 1e-6),
                (0, "Ca(C2O4)(s)"): (0.00, None, 1e-6),
                (0, "Ce2(C2O4)3(s)"): (0.003533, None, 1e-6),
                (0, "Dy2(C2O4)3(s)"): (0.002344, None, 1e-6),
                (0, "Fe2(C2O4)3(s)"): (0.000192, None, 1e-6),
                (0, "Gd2(C2O4)3(s)"): (0.003037, None, 1e-6),
                (0, "La2(C2O4)3(s)"): (0.003034, None, 1e-6),
                (0, "Nd2(C2O4)3(s)"): (0.003427, None, 1e-6),
                (0, "Pr2(C2O4)3(s)"): (0.0034732, None, 1e-6),
                (0, "Sc2(C2O4)3(s)"): (0.004009, None, 1e-6),
                (0, "Sm2(C2O4)3(s)"): (0.003253, None, 1e-6),
                (0, "Y2(C2O4)3(s)"): (0.004589, None, 1e-6),
            },
            "unit.precipitate_outlet.temperature": {
                0: (348.15, None, 1e-3),
            },
        }
        assert_solution_equivalent(prec.fs, expected_results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, prec):
        assert value(
            prec.fs.unit.aqueous_inlet.flow_vol[0] * prec.fs.properties_aq.dens_mass
        ) == pytest.approx(
            value(
                prec.fs.unit.aqueous_outlet.flow_vol[0]
                * prec.fs.properties_aq.dens_mass
            ),
            rel=1e-6,
            abs=1e-6,
        )

        stochiometry = {
            "Sc2(C2O4)3(s)": 2,
            "Y2(C2O4)3(s)": 2,
            "La2(C2O4)3(s)": 2,
            "Ce2(C2O4)3(s)": 2,
            "Pr2(C2O4)3(s)": 2,
            "Nd2(C2O4)3(s)": 2,
            "Sm2(C2O4)3(s)": 2,
            "Gd2(C2O4)3(s)": 2,
            "Dy2(C2O4)3(s)": 2,
            "Al2(C2O4)3(s)": 2,
            "Ca(C2O4)(s)": 1,
            "Fe2(C2O4)3(s)": 2,
        }

        reversed_react = dict(map(reversed, prec.fs.properties_solid.react.items()))
        pass_through_elements = ["Cl", "SO4", "H2O", "HSO4"]
        for j in prec.fs.properties_aq.dissolved_elements:
            if j in ["H", "H2C2O4"]:
                pass
            elif j in pass_through_elements:
                assert value(
                    prec.fs.unit.mscontactor.liquid_inlet_state[0.0].flow_mol_comp[j]
                ) == pytest.approx(
                    value(prec.fs.unit.mscontactor.liquid[0.0, 1].flow_mol_comp[j]),
                    rel=1e-5,
                    abs=1e-5,
                )
            else:
                assert value(
                    prec.fs.unit.mscontactor.liquid_inlet_state[0.0].flow_mol_comp[j]
                ) == pytest.approx(
                    value(
                        prec.fs.unit.mscontactor.liquid[0.0, 1].flow_mol_comp[j]
                        + (
                            prec.fs.unit.mscontactor.solid[0.0, 1].flow_mol_comp[
                                reversed_react[j]
                            ]
                            * stochiometry[reversed_react[j]]
                        )
                    ),
                    rel=1e-5,
                    abs=1e-5,
                )


# -----------------------------------------------------------------------------
class TestPrecRob(object):
    @pytest.fixture(scope="class")
    def prec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties_aq = AqueousParameter()
        m.fs.properties_solid = PrecipitateParameters()
        m.fs.prec_rxns = OxalatePrecipitationReactions()

        m.fs.unit = OxalatePrecipitator(
            number_of_tanks=1,
            liquid_phase={
                "property_package": m.fs.properties_aq,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            solid_phase={
                "property_package": m.fs.properties_solid,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            reaction_package=m.fs.prec_rxns,
        )

        m.fs.unit.aqueous_inlet.flow_vol[0].fix(100)

        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Al"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Y"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "La"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(10)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e-9)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2C2O4"].fix(12000)
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(100000)

        m.fs.unit.aqueous_inlet.temperature[0].fix(298.15)
        m.fs.unit.aqueous_inlet.pressure[0].fix(1e5)

        m.fs.unit.precipitate_outlet.temperature.fix(348.15)
        m.fs.unit.hydraulic_retention_time[0.0].fix(2)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, prec):
        assert hasattr(prec.fs.unit, "aqueous_inlet")
        assert len(prec.fs.unit.aqueous_inlet.vars) == 4
        assert hasattr(prec.fs.unit.aqueous_inlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_inlet, "conc_mass_comp")
        assert hasattr(prec.fs.unit.aqueous_inlet, "temperature")
        assert hasattr(prec.fs.unit.aqueous_inlet, "pressure")

        assert hasattr(prec.fs.unit, "aqueous_outlet")
        assert len(prec.fs.unit.aqueous_outlet.vars) == 4
        assert hasattr(prec.fs.unit.aqueous_outlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_outlet, "conc_mass_comp")
        assert hasattr(prec.fs.unit.aqueous_outlet, "temperature")
        assert hasattr(prec.fs.unit.aqueous_outlet, "pressure")

        assert hasattr(prec.fs.unit, "precipitate_outlet")
        assert len(prec.fs.unit.precipitate_outlet.vars) == 2
        assert hasattr(prec.fs.unit.precipitate_outlet, "flow_mol_comp")
        assert hasattr(prec.fs.unit.precipitate_outlet, "temperature")

        assert hasattr(prec.fs.unit, "init_solid_constraint")
        assert hasattr(prec.fs.unit, "temp_constraint")
        assert hasattr(prec.fs.unit, "heterogeneous_reaction_extent_constraint")
        assert hasattr(prec.fs.unit, "eq_hydraulic_retention")

        assert number_variables(prec.fs.unit) == 208
        assert number_total_constraints(prec.fs.unit) == 185
        assert number_unused_variables(prec.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, prec):
        assert_units_consistent(prec.fs.unit)

    @pytest.mark.unit
    def test_dof(self, prec):
        assert degrees_of_freedom(prec) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, prec):
        # scaling = TransformationFactory("core.scale_model")
        # scaled_model = scaling.create_using(prec, rename=False)
        # initializer = OxalatePrecipitatorInitializer()
        # initializer.initialize(prec.fs.unit)
        # # initializer = scaled_model.fs.unit.default_initializer()
        # # initializer.initialize(scaled_model.fs.unit)

        # # Solve scaled model
        # solver = get_solver("ipopt_v2")
        # results = solver.solve(scaled_model, tee=False)

        # # Propagate results back to unscaled model
        # scaling.propagate_solution(scaled_model, prec)

        # assert_optimal_termination(results)

        initializer = OxalatePrecipitatorInitializer()
        initializer.initialize(prec.fs.unit)
        assert initializer.summary[prec.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.component
    def test_var_scaling(self, prec):
        unscaled_var_list = list(
            unscaled_variables_generator(prec.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    @pytest.mark.solver
    def test_numerical_issues(self, prec):
        dt = DiagnosticsToolbox(prec)
        dt.assert_no_numerical_warnings()
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        expected_results = {
            "unit.aqueous_outlet.flow_vol": {
                0: (100, None, 1e-0),
            },
            "unit.aqueous_outlet.conc_mass_comp": {
                (0, "Al"): (9.730, None, 1e-3),
                (0, "Ca"): (9.9999, None, 1e-3),
                (0, "Ce"): (0.01782, None, 1e-3),
                (0, "Dy"): (0.1380, None, 1e-3),
                (0, "Fe"): (2.0863, None, 1e-3),
                (0, "Gd"): (0.0327, None, 1e-3),
                (0, "La"): (0.0924, 1e-3, None),
                (0, "Nd"): (0.0250, None, 1e-3),
                (0, "Pr"): (0.0254, None, 1e-3),
                (0, "Sc"): (0.1791, None, 1e-3),
                (0, "Sm"): (0.0212, None, 1e-3),
                (0, "Y"): (0.1080, None, 1e-3),
            },
            "unit.precipitate_outlet.flow_mol_comp": {
                (0, "Al2(C2O4)3(s)"): (0.0005, None, 1e-6),
                (0, "Ca(C2O4)(s)"): (1.3e-08, None, 1e-6),
                (0, "Ce2(C2O4)3(s)"): (0.003562, None, 1e-6),
                (0, "Dy2(C2O4)3(s)"): (0.003034, None, 1e-6),
                (0, "Fe2(C2O4)3(s)"): (0.007086, None, 1e-6),
                (0, "Gd2(C2O4)3(s)"): (0.003169, None, 1e-6),
                (0, "La2(C2O4)3(s)"): (0.003566, None, 1e-6),
                (0, "Nd2(C2O4)3(s)"): (0.003458, None, 1e-6),
                (0, "Pr2(C2O4)3(s)"): (0.003540, None, 1e-6),
                (0, "Sc2(C2O4)3(s)"): (0.010926, None, 1e-6),
                (0, "Sm2(C2O4)3(s)"): (0.003318, None, 1e-6),
                (0, "Y2(C2O4)3(s)"): (0.005564, None, 1e-6),
            },
            "unit.precipitate_outlet.temperature": {
                0: (348.15, None, 1e-3),
            },
        }
        assert_solution_equivalent(prec.fs, expected_results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, prec):
        assert value(
            prec.fs.unit.aqueous_inlet.flow_vol[0] * prec.fs.properties_aq.dens_mass
        ) == pytest.approx(
            value(
                prec.fs.unit.aqueous_outlet.flow_vol[0]
                * prec.fs.properties_aq.dens_mass
            ),
            rel=1e-6,
            abs=1e-6,
        )

        stochiometry = {
            "Sc2(C2O4)3(s)": 2,
            "Y2(C2O4)3(s)": 2,
            "La2(C2O4)3(s)": 2,
            "Ce2(C2O4)3(s)": 2,
            "Pr2(C2O4)3(s)": 2,
            "Nd2(C2O4)3(s)": 2,
            "Sm2(C2O4)3(s)": 2,
            "Gd2(C2O4)3(s)": 2,
            "Dy2(C2O4)3(s)": 2,
            "Al2(C2O4)3(s)": 2,
            "Ca(C2O4)(s)": 1,
            "Fe2(C2O4)3(s)": 2,
        }

        reversed_react = dict(map(reversed, prec.fs.properties_solid.react.items()))
        pass_through_elements = ["Cl", "SO4", "H2O", "HSO4"]
        for j in prec.fs.properties_aq.dissolved_elements:
            if j in ["H", "H2C2O4"]:
                pass
            elif j in pass_through_elements:
                assert value(
                    prec.fs.unit.mscontactor.liquid_inlet_state[0.0].flow_mol_comp[j]
                ) == pytest.approx(
                    value(prec.fs.unit.mscontactor.liquid[0.0, 1].flow_mol_comp[j]),
                    rel=1e-5,
                    abs=1e-5,
                )
            else:
                assert value(
                    prec.fs.unit.mscontactor.liquid_inlet_state[0.0].flow_mol_comp[j]
                ) == pytest.approx(
                    value(
                        prec.fs.unit.mscontactor.liquid[0.0, 1].flow_mol_comp[j]
                        + (
                            prec.fs.unit.mscontactor.solid[0.0, 1].flow_mol_comp[
                                reversed_react[j]
                            ]
                            * stochiometry[reversed_react[j]]
                        )
                    ),
                    rel=1e-5,
                    abs=1e-5,
                )
