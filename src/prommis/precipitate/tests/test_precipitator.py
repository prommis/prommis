from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import unscaled_variables_generator

import pytest

from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator

import idaes.logger as idaeslog

_log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

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

    m.fs.unit = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    assert len(m.fs.unit.config) == 11

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert not m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package_aqueous is m.fs.properties_aq
    assert m.fs.unit.config.property_package_precipitate is m.fs.properties_solid


# -----------------------------------------------------------------------------
class TestPrec(object):
    @pytest.fixture(scope="class")
    def prec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties_aq = AqueousParameter()
        m.fs.properties_solid = PrecipitateParameters()

        m.fs.unit = Precipitator(
            property_package_aqueous=m.fs.properties_aq,
            property_package_precipitate=m.fs.properties_solid,
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

        m.fs.unit.cv_precipitate.properties_in[0].temperature.fix(348.15)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, prec):
        assert hasattr(prec.fs.unit, "aqueous_inlet")
        assert len(prec.fs.unit.aqueous_inlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_inlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_inlet, "conc_mass_comp")

        assert hasattr(prec.fs.unit, "aqueous_outlet")
        assert len(prec.fs.unit.aqueous_outlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_outlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_outlet, "conc_mass_comp")

        assert hasattr(prec.fs.unit, "precipitate_outlet")
        assert len(prec.fs.unit.precipitate_outlet.vars) == 2
        assert hasattr(prec.fs.unit.precipitate_outlet, "flow_mol_comp")
        assert hasattr(prec.fs.unit.precipitate_outlet, "temperature")

        assert hasattr(prec.fs.unit, "generation")
        assert hasattr(prec.fs.unit, "mass_balance")
        assert hasattr(prec.fs.unit, "vol_balance")

        assert number_variables(prec.fs.unit) == 98
        assert number_total_constraints(prec.fs.unit) == 84
        assert number_unused_variables(prec.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, prec):
        assert_units_consistent(prec.fs.unit)

        dt = DiagnosticsToolbox(model=prec)
        dt.report_structural_issues()
        dt.display_underconstrained_set()
        dt.display_overconstrained_set()
        assert degrees_of_freedom(prec) == 0

    @pytest.mark.unit
    def test_dof(self, prec):
        assert degrees_of_freedom(prec) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, prec):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(prec)
        assert initializer.summary[prec]["status"] == InitializationStatus.Ok

    @pytest.mark.component
    def test_var_scaling(self, prec):
        unscaled_var_list = list(
            unscaled_variables_generator(prec.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, prec):
        solver = get_solver()
        results = solver.solve(prec)
        prec.fs.unit.aqueous_outlet.conc_mass_comp.display()
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        assert pytest.approx(100, abs=1e-0) == value(
            prec.fs.unit.aqueous_outlet.flow_vol[0]
        )
        assert pytest.approx(10.0, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Al"]
        )
        assert pytest.approx(10, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ca"]
        )
        assert pytest.approx(3.887, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ce"]
        )
        assert pytest.approx(2.161, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Dy"]
        )
        assert pytest.approx(8.897, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Fe"]
        )
        assert pytest.approx(2.019, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Gd"]
        )
        assert pytest.approx(5.570, rel=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "La"]
        )
        assert pytest.approx(2.2913, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Nd"]
        )
        assert pytest.approx(2.802, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Pr"]
        )
        assert pytest.approx(5.514, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sc"]
        )
        assert pytest.approx(2.039, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sm"]
        )
        assert pytest.approx(3.283, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Y"]
        )
        assert pytest.approx(1e-20, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Al2(C2O4)3(s)"]
        )
        # assert pytest.approx(0.0051150, abs=1e-6) == value(
        #     prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ca(C2O4)(s)"]
        # )
        assert pytest.approx(0.002181, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ce2(C2O4)3(s)"]
        )
        assert pytest.approx(0.002411, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Dy2(C2O4)3(s)"]
        )
        assert pytest.approx(0.000986, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Fe2(C2O4)3(s)"]
        )
        assert pytest.approx(0.002537, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Gd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.001594, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "La2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0026721, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Nd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0025540, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Pr2(C2O4)3(s)"]
        )
        assert pytest.approx(0.004989, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sc2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0026469, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sm2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003777, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Y2(C2O4)3(s)"]
        )
        assert pytest.approx(348.15, abs=1e-3) == value(
            prec.fs.unit.precipitate_outlet.temperature[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, prec):
        assert (
            abs(
                value(
                    prec.fs.unit.aqueous_inlet.flow_vol[0]
                    * prec.fs.properties_aq.dens_mass
                    - prec.fs.unit.aqueous_outlet.flow_vol[0]
                    * prec.fs.properties_aq.dens_mass
                )
            )
            <= 1e-6
        )

        reversed_react = dict(map(reversed, prec.fs.properties_solid.react.items()))
        for j in prec.fs.properties_aq.dissolved_elements:
            if j == "Ca":
                pass
            else:
                assert (
                    abs(
                        value(
                            prec.fs.unit.cv_aqueous.properties_in[0].flow_mol_comp[j]
                            - (
                                prec.fs.unit.cv_aqueous.properties_out[0].flow_mol_comp[
                                    j
                                ]
                                + (
                                    prec.fs.unit.precipitate_outlet.flow_mol_comp[
                                        0, reversed_react[j]
                                    ]
                                    * prec.fs.properties_solid.stoich[reversed_react[j]]
                                )
                            )
                        )
                    )
                    <= 1e-5
                )
