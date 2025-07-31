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
    Suffix,
    SolverFactory,
    TransformationFactory,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import unscaled_variables_generator

from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitate_reactions import OxalatePrecipitationReactions
from prommis.precipitate.precipitator import (
    OxalatePrecipitator,
    OxalatePrecipitatorInitializer,
)

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
        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(prec, rename=False)
        initializer = prec.fs.unit.default_initializer()
        try:
            initializer.initialize(prec.fs.unit)
        except:
            pass

        # Solve scaled model
        solver = SolverFactory("ipopt")
        results = solver.solve(scaled_model, tee=False)

        # Propagate results back to unscaled model
        scaling.propagate_solution(scaled_model, prec)

        assert_optimal_termination(results)

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
        dt.report_numerical_issues()
        dt.display_variables_at_or_outside_bounds()
        dt.assert_no_numerical_warnings()
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        assert pytest.approx(100, abs=1e-0) == value(
            prec.fs.unit.aqueous_outlet.flow_vol[0]
        )
        assert pytest.approx(9.982, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Al"]
        )
        assert pytest.approx(10, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ca"]
        )
        assert pytest.approx(0.0991, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ce"]
        )
        assert pytest.approx(2.3824, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Dy"]
        )
        assert pytest.approx(9.7864, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Fe"]
        )
        assert pytest.approx(0.4486, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Gd"]
        )
        assert pytest.approx(1.5715, rel=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "La"]
        )
        assert pytest.approx(0.1142, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Nd"]
        )
        assert pytest.approx(0.2118, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Pr"]
        )
        assert pytest.approx(6.3964, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sc"]
        )
        assert pytest.approx(0.2183, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sm"]
        )
        assert pytest.approx(1.8401, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Y"]
        )
        assert pytest.approx(3.20213e-05, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Al2(C2O4)3(s)"]
        )
        assert pytest.approx(0.00, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ca(C2O4)(s)"]
        )
        assert pytest.approx(0.003533, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ce2(C2O4)3(s)"]
        )
        assert pytest.approx(0.002344, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Dy2(C2O4)3(s)"]
        )
        assert pytest.approx(0.000192, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Fe2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003037, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Gd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003034, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "La2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003427, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Nd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.0034732, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Pr2(C2O4)3(s)"]
        )
        assert pytest.approx(0.004009, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sc2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003253, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sm2(C2O4)3(s)"]
        )
        assert pytest.approx(0.004589, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Y2(C2O4)3(s)"]
        )
        assert pytest.approx(348.15, abs=1e-3) == value(
            prec.fs.unit.precipitate_outlet.temperature[0]
        )

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
        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(prec, rename=False)
        initializer = OxalatePrecipitatorInitializer()
        initializer.initialize(prec.fs.unit)

        # Solve scaled model
        solver = SolverFactory("ipopt")
        results = solver.solve(scaled_model, tee=False)

        # Propagate results back to unscaled model
        scaling.propagate_solution(scaled_model, prec)

        assert_optimal_termination(results)

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

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, prec):
        assert pytest.approx(100, abs=1e-0) == value(
            prec.fs.unit.aqueous_outlet.flow_vol[0]
        )
        assert pytest.approx(9.730, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Al"]
        )
        assert pytest.approx(9.9999, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ca"]
        )
        assert pytest.approx(0.01782, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ce"]
        )
        assert pytest.approx(0.1380, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Dy"]
        )
        assert pytest.approx(2.0863, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Fe"]
        )
        assert pytest.approx(0.0327, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Gd"]
        )
        assert pytest.approx(0.0924, rel=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "La"]
        )
        assert pytest.approx(0.0250, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Nd"]
        )
        assert pytest.approx(0.0254, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Pr"]
        )
        assert pytest.approx(0.1791, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sc"]
        )
        assert pytest.approx(0.0212, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Sm"]
        )
        assert pytest.approx(0.1080, abs=1e-3) == value(
            prec.fs.unit.aqueous_outlet.conc_mass_comp[0, "Y"]
        )
        assert pytest.approx(0.0005, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Al2(C2O4)3(s)"]
        )
        assert pytest.approx(1.3e-08, abs=1e-3) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ca(C2O4)(s)"]
        )
        assert pytest.approx(0.003562, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Ce2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003034, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Dy2(C2O4)3(s)"]
        )
        assert pytest.approx(0.007086, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Fe2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003169, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Gd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003566, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "La2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003458, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Nd2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003540, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Pr2(C2O4)3(s)"]
        )
        assert pytest.approx(0.010926, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sc2(C2O4)3(s)"]
        )
        assert pytest.approx(0.003318, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Sm2(C2O4)3(s)"]
        )
        assert pytest.approx(0.005564, abs=1e-6) == value(
            prec.fs.unit.precipitate_outlet.flow_mol_comp[0, "Y2(C2O4)3(s)"]
        )
        assert pytest.approx(348.15, abs=1e-3) == value(
            prec.fs.unit.precipitate_outlet.temperature[0]
        )

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
