"""
Demonstration flowsheet for EvaporationPond unit model.

Authors: Andrew Lee
"""
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    SolverFactory,
    units,
    value,
    Var,
)

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.core.initialization import (
    InitializationStatus,
)
from idaes.core.solvers import get_solver

from prommis.evaporation_pond.evaporation_pond import (
    EvaporationPond,
    EvaporationPondInitializer,
)
from prommis.evaporation_pond.tests.example_properties import BrineParameters
from prommis.evaporation_pond.tests.example_reactions import (
    BrineReactionParameters,
)


class TestExampleCase():
    @pytest.fixture(scope="class")
    def model(self):
        """
        Method to build a single stage evaporation pond model for testing
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.brine_props = BrineParameters()
        m.fs.brine_rxns = BrineReactionParameters(
            property_package=m.fs.brine_props,
        )

        m.fs.pond_1 = EvaporationPond(
            property_package=m.fs.brine_props,
            reaction_package=m.fs.brine_rxns,
        )

        m.fs.pond_1.inlet.flow_vol.fix(240 * units.L / units.s)
        m.fs.pond_1.inlet.conc_mass_comp[0, "Li"].fix(650 * units.mg / units.L)
        m.fs.pond_1.inlet.conc_mass_comp[0, "Na"].fix(650 * units.mg / units.L)
        m.fs.pond_1.inlet.conc_mass_comp[0, "Cl"].fix(4310 * units.mg / units.L)
        m.fs.pond_1.inlet.conc_mass_comp[0, "H2O"].fix(1 * units.kg / units.L)

        m.fs.pond_1.surface_area.fix(50000 * units.m**2)
        m.fs.pond_1.average_pond_depth.fix(0.5 * units.m)
        m.fs.pond_1.evaporation_rate.fix(4.75 * units.mm / units.day)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.pond_1.default_initializer is EvaporationPondInitializer

        assert isinstance(model.fs.pond_1.surface_area, Var)
        assert isinstance(model.fs.pond_1.average_pond_depth, Var)
        assert isinstance(model.fs.pond_1.volume, Var)
        assert isinstance(model.fs.pond_1.evaporation_rate, Var)
        assert isinstance(model.fs.pond_1.water_loss_rate, Var)
        assert isinstance(model.fs.pond_1.precipitation_rate, Var)
        assert isinstance(model.fs.pond_1.reaction_extent, Var)

        assert isinstance(model.fs.pond_1.eps, Param)
        assert isinstance(model.fs.pond_1.s_norm, Param)
        assert isinstance(model.fs.pond_1.s_scale, Param)

        assert isinstance(model.fs.pond_1.evaporation_constraint, Constraint)
        assert isinstance(model.fs.pond_1.volume_constraint, Constraint)
        assert isinstance(model.fs.pond_1.component_balances, Constraint)
        assert isinstance(model.fs.pond_1.stoichiometry_constraint, Constraint)
        assert isinstance(model.fs.pond_1.equilibrium_constraint, Constraint)

    @pytest.mark.unit
    def test_no_structural_warnings(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.pond_1._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Surface Area": model.fs.pond_1.surface_area[0],
                "Average Depth": model.fs.pond_1.average_pond_depth[0],
                "Volume": model.fs.pond_1.volume[0],
                "Evaporation Rate": model.fs.pond_1.evaporation_rate[0],
                "Water Loss Rate": model.fs.pond_1.water_loss_rate[0],
                "Precipitation Rate Cl": model.fs.pond_1.precipitation_rate[0, "Cl"],
                "Precipitation Rate H2O": model.fs.pond_1.precipitation_rate[0, "H2O"],
                "Precipitation Rate Li": model.fs.pond_1.precipitation_rate[0, "Li"],
                "Precipitation Rate Na": model.fs.pond_1.precipitation_rate[0, "Na"],
            }
        }

    @pytest.mark.component
    @pytest.mark.solver
    def test_initialization(self, model):
        initializer = EvaporationPondInitializer()
        initializer.initialize(model.fs.pond_1)

        assert initializer.summary[model.fs.pond_1]["status"] == InitializationStatus.Ok

    @pytest.mark.component
    @pytest.mark.solver
    def test_solve(self, model):
        solver = get_solver()
        res = solver.solve(model)

        assert_optimal_termination(res)

    @pytest.mark.unit
    def test_no_numerical_warnings(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model):
        assert value(model.fs.pond_1.outlet.flow_vol[0]) == pytest.approx(854104, rel=1e-5)
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "H2O"]) == pytest.approx(1e6, rel=1e-5)
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Cl"]) == pytest.approx(3933.85, rel=1e-5)
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Li"]) == pytest.approx(657.531, rel=1e-5)
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Na"]) == pytest.approx(381.203, rel=1e-5)

        assert value(model.fs.pond_1.volume[0]) == pytest.approx(25000, rel=1e-5)
        assert value(model.fs.pond_1.precipitation_rate[0, "Na"]) == pytest.approx(-10265.9, rel=1e-5)
        assert value(model.fs.pond_1.precipitation_rate[0, "Cl"]) == pytest.approx(-10265.9, rel=1e-5)
        assert value(model.fs.pond_1.precipitation_rate[0, "Li"]) == pytest.approx(0, abs=1e-5)
        assert value(model.fs.pond_1.precipitation_rate[0, "H2O"]) == pytest.approx(0, abs=1e-5)
        assert value(model.fs.pond_1.water_loss_rate[0]) == pytest.approx(-549768, rel=1e-5)

    @pytest.mark.component
    @pytest.mark.solver
    def test_conservation(self, model):
        assert value(
            model.fs.pond_1.inlet.flow_vol[0]
            * model.fs.pond_1.properties_in[0].conc_mole_comp["H2O"]
            + model.fs.pond_1.water_loss_rate[0]
        ) == pytest.approx(
            value(
                model.fs.pond_1.outlet.flow_vol[0]
                * model.fs.pond_1.properties_out[0].conc_mole_comp["H2O"]
            ),
            rel=1e-5,
        )

        assert value(
            model.fs.pond_1.inlet.flow_vol[0]
            * model.fs.pond_1.properties_in[0].conc_mole_comp["Li"]
        ) == pytest.approx(
            value(
                model.fs.pond_1.outlet.flow_vol[0]
                * model.fs.pond_1.properties_out[0].conc_mole_comp["Li"]
            ),
            rel=1e-5,
        )

        assert value(
            model.fs.pond_1.inlet.flow_vol[0]
            * model.fs.pond_1.properties_in[0].conc_mole_comp["Na"]
            + model.fs.pond_1.precipitation_rate[0, "Na"]
        ) == pytest.approx(
            value(
                model.fs.pond_1.outlet.flow_vol[0]
                * model.fs.pond_1.properties_out[0].conc_mole_comp["Na"]
            ),
            rel=1e-5,
        )

        assert value(
            model.fs.pond_1.inlet.flow_vol[0]
            * model.fs.pond_1.properties_in[0].conc_mole_comp["Cl"]
            + model.fs.pond_1.precipitation_rate[0, "Cl"]
        ) == pytest.approx(
            value(
                model.fs.pond_1.outlet.flow_vol[0]
                * model.fs.pond_1.properties_out[0].conc_mole_comp["Cl"]
            ),
            rel=1e-5,
        )
