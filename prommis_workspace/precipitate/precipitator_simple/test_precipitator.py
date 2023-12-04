import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
    units,
)

from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from idaes.core.util.scaling import (
    unscaled_variables_generator,
)

from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from pyomo.util.check_units import assert_units_consistent

from precipitator_new import Precipitator
from precipitate_solids_properties import PrecipitateParameters
from precipitate_liquid_properties import AqueousParameter

from idaes.core.util.model_diagnostics import DiagnosticsToolbox

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

        assert number_variables(prec.fs.unit) == 74
        assert number_total_constraints(prec.fs.unit) == 60
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
        initializer.initialize(prec.fs.unit)
        assert initializer.summary[prec.fs.unit]["status"] == InitializationStatus.Ok

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        assert pytest.approx(100, abs=1e-0) == value(
            prec.fs.unit.aqueous_outlet.flow_vol[0]
        )
        assert pytest.approx(9.91, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Al"]
        )
        assert pytest.approx(10, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ca"]
        )
        assert pytest.approx(3.193, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ce"]
        )
        assert pytest.approx(1.284, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Dy"]
        )
        assert pytest.approx(9.756, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Fe"]
        )
        assert pytest.approx(1.199, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Gd"]
        )
        assert pytest.approx(4.849, rel=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "La"]
        )
        assert pytest.approx(1.845, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Nd"]
        )
        assert pytest.approx(2.2, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Pr"]
        )
        assert pytest.approx(6.839, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sc"]
        )
        assert pytest.approx(1.265, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sm"]
        )
        assert pytest.approx(2.554, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Y"]
        )
        assert pytest.approx(0.00016678, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Al2(C2O4)3(s)"]
        )
        # assert pytest.approx(0.0051150, abs=1e-6) == value(
        #     prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ca(C2O4)(s)"]
        # )
        assert pytest.approx(0.002429, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ce2(C2O4)3(s)"]
        )
        assert pytest.approx(0.002682, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Dy2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0002189, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Fe2(C2O4)3(s)"]
        )
        assert pytest.approx(0.002799, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Gd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.001854, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "La2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0028268, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Nd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0027677, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Pr2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003517, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sc2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0029047, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sm2(C2O4)3(s)"]
        )
        assert pytest.approx(0.004188, abs=1e-6) == value(
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
