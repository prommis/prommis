#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Demonstration flowsheet for EvaporationPond unit model.

Authors: Andrew Lee
"""

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    Set,
    Var,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    Component,
    FlowsheetBlock,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    ReactionBlockBase,
    ReactionBlockDataBase,
    ReactionParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.initialization import InitializationStatus
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError

import pytest

from prommis.evaporation_pond.evaporation_pond import (
    EvaporationPond,
    EvaporationPondInitializer,
)
from prommis.evaporation_pond.tests.example_properties import BrineParameters
from prommis.evaporation_pond.tests.example_reactions import BrineReactionParameters

# TODO: Conversion of bases?


@pytest.mark.unit
def test_invalid_solvent_id():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.brine_props = BrineParameters()
    m.fs.brine_rxns = BrineReactionParameters(
        property_package=m.fs.brine_props,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.pond_1 - foo was set as the solvent_id for "
        "EvaporationPond model, however this is not a valid component name "
        "in the property package provided.",
    ):
        m.fs.pond_1 = EvaporationPond(
            property_package=m.fs.brine_props,
            reaction_package=m.fs.brine_rxns,
            solvent_id="foo",
        )


@declare_process_block_class("DummyParameters")
class DummyParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()
        self.phase2 = Phase()

        self.comp1 = Component()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@pytest.mark.unit
def test_multiple_phases():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = DummyParameters()

    with pytest.raises(
        ConfigurationError,
        match="fs.pond_1 - EvaporationPond model only supports a single phase. However, "
        "the property package provided includes 2 phases.",
    ):
        m.fs.pond_1 = EvaporationPond(
            property_package=m.fs.props,
            solvent_id="comp1",
        )


@declare_process_block_class("DummyParameters2")
class DummyParameter2Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()
        self.comp1 = Component()

        self._state_block_class = DummyStateBlock2

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("DummyStateBlock2", block_class=StateBlock)
class DummyStateBlockData2(StateBlockData):
    """
    State block for lithium brine solutions.

    """

    def build(self):
        super().build()

        self.dummy_state = Var()

    def get_material_flow_basis(self):
        return MaterialFlowBasis.other

    def define_state_vars(self):
        return {
            "dummy_state": self.dummy_state,
        }


@pytest.mark.unit
def test_flow_basis_other():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = DummyParameters2()
    m.fs.brine_rxns = BrineReactionParameters(
        property_package=m.fs.props,
    )

    with pytest.raises(
        PropertyPackageError,
        match="fs.pond_1 - Property package uses a flow basis of 'other'. "
        "EvaporationPond model only supports mass or molar bases.",
    ):
        m.fs.pond_1 = EvaporationPond(
            property_package=m.fs.props,
            reaction_package=m.fs.brine_rxns,
            solvent_id="comp1",
        )


@declare_process_block_class("DummyReactionParameters")
class DummyReactionParametersData(ReactionParameterBlock):
    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = DummyReactionBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("DummyReactionBlock", block_class=ReactionBlockBase)
class DummyReactionData(ReactionBlockDataBase):
    def build(self):
        super().build()

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.other


@pytest.mark.unit
def test_reaction_basis_other():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = BrineParameters()
    m.fs.brine_rxns = DummyReactionParameters(
        property_package=m.fs.props,
    )

    with pytest.raises(
        PropertyPackageError,
        match="fs.pond_1 - Reaction package uses a reaction basis other than mass or "
        "mole. EvaporationPond model only supports mass or molar bases.",
    ):
        m.fs.pond_1 = EvaporationPond(
            property_package=m.fs.props,
            reaction_package=m.fs.brine_rxns,
        )


@declare_process_block_class("DummyReactionParameters2")
class DummyReactionParametersData2(ReactionParameterBlock):
    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = DummyReactionBlock2

        # Reaction Index
        self.equilibrium_reaction_idx = Set(initialize=["P1"])

        # Reaction Stoichiometry
        self.equilibrium_reaction_stoichiometry = {
            ("P1", "liquid", "Li"): 0,
            ("P1", "liquid", "Na"): -1,
            ("P1", "liquid", "Cl"): -1,
            ("P1", "liquid", "H2O"): 0,
        }

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("DummyReactionBlock2", block_class=ReactionBlockBase)
class DummyReactionData2(ReactionBlockDataBase):
    def build(self):
        super().build()

        self.solubility_product_P1 = Param(
            default=1,
            units=units.mass**2 / units.L**2,
            mutable=True,
            doc="Solubility constant for reaction P1",
        )

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.mass


@pytest.mark.unit
def test_reaction_basis_conversion_mass():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = BrineParameters()
    m.fs.brine_rxns = DummyReactionParameters2(
        property_package=m.fs.props,
    )

    m.fs.pond_1 = EvaporationPond(
        property_package=m.fs.props,
        reaction_package=m.fs.brine_rxns,
    )

    for j in ["Na", "Cl"]:
        assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", j].expr) == str(
            m.fs.pond_1.precipitation_rate[0, j]
            == -1 * m.fs.pond_1.reaction_extent[0, "P1"] / m.fs.props.mw[j]
        )
    for j in ["H2O", "Li"]:
        assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", j].expr) == str(
            m.fs.pond_1.precipitation_rate[0, j] == 0 * (units.mol * units.hour**-1)
        )

    assert_units_consistent(m)


@declare_process_block_class("DummyParameters3")
class DummyParameterData3(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Solutes
        self.Li = Component()
        self.Na = Component()
        self.Cl = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "Li": 6.94e-3,
                "Na": 22.99e-3,
                "Cl": 35.45e-3,
            },
        )

        self._state_block_class = DummyStateBlock3

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("DummyStateBlock3", block_class=StateBlock)
class DummyStateBlockData3(StateBlockData):
    def build(self):
        super().build()

        self.dummy_flow = Var(
            units=units.kg / units.hour,
        )

        self.conc_mole_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            bounds=(1e-20, None),
        )

    @property
    def mw(self):
        """Molecular weight of species"""
        return self.params.mw

    def get_material_flow_terms(self, _, j):
        return self.dummy_flow

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "conc_mole_comp": self.conc_mole_comp,
        }


@pytest.mark.unit
def test_reaction_basis_conversion_mole():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = DummyParameters3()
    m.fs.brine_rxns = BrineReactionParameters(
        property_package=m.fs.props,
    )

    m.fs.pond_1 = EvaporationPond(
        property_package=m.fs.props,
        reaction_package=m.fs.brine_rxns,
    )

    assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", "Na"].expr) == str(
        m.fs.pond_1.precipitation_rate[0, "Na"]
        == -m.fs.props.mw["Na"] * m.fs.pond_1.reaction_extent[0, "P1"]
    )
    assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", "Li"].expr) == str(
        m.fs.pond_1.precipitation_rate[0, "Li"]
        == -m.fs.props.mw["Li"] * m.fs.pond_1.reaction_extent[0, "P2"]
    )
    assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", "Cl"].expr) == str(
        m.fs.pond_1.precipitation_rate[0, "Cl"]
        == (
            -m.fs.pond_1.reaction_extent[0, "P1"] - m.fs.pond_1.reaction_extent[0, "P2"]
        )
        * m.fs.props.mw["Cl"]
    )
    assert str(m.fs.pond_1.stoichiometry_constraint[0, "liquid", "H2O"].expr) == str(
        m.fs.pond_1.precipitation_rate[0, "H2O"] == 0 * (units.kg * units.hour**-1)
    )

    assert_units_consistent(m)


class TestExampleCase:
    @pytest.fixture(scope="class")
    def model(self):
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
        assert value(model.fs.pond_1.outlet.flow_vol[0]) == pytest.approx(
            854104, rel=1e-5
        )
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "H2O"]) == pytest.approx(
            1e6, rel=1e-5
        )
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Cl"]) == pytest.approx(
            3933.85, rel=1e-5
        )
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Li"]) == pytest.approx(
            657.531, rel=1e-5
        )
        assert value(model.fs.pond_1.outlet.conc_mass_comp[0, "Na"]) == pytest.approx(
            381.203, rel=1e-5
        )

        assert value(model.fs.pond_1.volume[0]) == pytest.approx(25000, rel=1e-5)
        assert value(model.fs.pond_1.precipitation_rate[0, "Na"]) == pytest.approx(
            -10265.9, rel=1e-5
        )
        assert value(model.fs.pond_1.precipitation_rate[0, "Cl"]) == pytest.approx(
            -10265.9, rel=1e-5
        )
        assert value(model.fs.pond_1.precipitation_rate[0, "Li"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.pond_1.precipitation_rate[0, "H2O"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.pond_1.water_loss_rate[0]) == pytest.approx(
            -549768, rel=1e-5
        )

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
