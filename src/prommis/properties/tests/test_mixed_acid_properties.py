#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re

from pyomo.environ import ConcreteModel, Constraint, Param, Set, value, Var

from idaes.core import FlowsheetBlock
from idaes.core.scaling.util import (
    jacobian_cond,
    list_unscaled_constraints,
    list_unscaled_variables,
)
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.properties.mixed_acid_properties import (
    MixedAcidParameterBlock,
    MixedAcidPropertiesScaler,
    contaminant_list,
    ree_list,
)


class TestHClStripping(object):
    class TestDefinedStateTrue(object):
        @pytest.fixture(scope="class")
        def frame(self):
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            m.fs.prec_sol = MixedAcidParameterBlock()

            m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time, defined_state=True)

            return m

        @pytest.mark.unit
        def test_build(self, frame):
            m = frame

            assert len(m.fs.state) == 1

            assert isinstance(m.fs.state[0].flow_vol, Var)
            assert isinstance(m.fs.state[0].conc_mass_comp, Var)
            assert isinstance(m.fs.state[0].flow_mol_comp, Var)
            assert not m.fs.state[0].is_property_constructed("pH_phase")

            assert isinstance(m.fs.state[0].conc_mol_comp_eqn, Constraint)
            assert isinstance(m.fs.state[0].flow_mol_comp_eqn, Constraint)
            assert not m.fs.state[0].is_property_constructed("h2o_concentration_eqn")
            assert not m.fs.state[0].is_property_constructed("pH_phase_eqn")

            m.fs.state[0].flow_vol.set_value(10)
            for i in m.fs.prec_sol.component_list:
                m.fs.state[0].conc_mass_comp[i].set_value(0.1)
                m.fs.state[0].flow_mol_comp[i].set_value(0.5)

            m.fs.state[0].pH_phase["liquid"].set_value(1)
            assert isinstance(m.fs.state[0].pH_phase_eqn, Constraint)

            m.fs.state.fix_initialization_states()

            assert degrees_of_freedom(m.fs.state[0]) == 0

            assert m.fs.state[0].flow_vol.fixed
            for j in m.fs.prec_sol.component_list:
                assert m.fs.state[0].conc_mass_comp[j].fixed
                assert not m.fs.state[0].flow_mol_comp[j].fixed

                assert m.fs.state[0].flow_mol_comp_eqn[j].active

        @pytest.mark.unit
        def test_scaling(self, frame):
            m = frame
            assert m.fs.state[0].default_scaler is MixedAcidPropertiesScaler
            scaler_obj = m.fs.state[0].default_scaler()

            with pytest.raises(
                ValueError,
                match=re.escape(
                    "This scaler requires the user to provide a default scaling "
                    "factor for fs.state[0.0].flow_vol, but no default scaling "
                    "factor was set."
                ),
            ):
                scaler_obj.scale_model(m.fs.state[0])

            scaler_obj.default_scaling_factors["flow_vol"] = 0.1
            for ree in ree_list:
                scaler_obj.default_scaling_factors[f"conc_mass_comp[{ree}]"] = 0.1
            for contaminant in contaminant_list:
                scaler_obj.default_scaling_factors[f"conc_mass_comp[{contaminant}]"] = (
                    0.1
                )

            scaler_obj.scale_model(m.fs.state[0])

            assert len(list_unscaled_variables(m.fs.state[0])) == 0
            assert len(list_unscaled_constraints(m.fs.state[0])) == 0

            assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(
                2.401290e6
            )
            assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(33.24336)

            m.fs.state[0].conc_mass_comp.unfix()

            assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(
                2.401274e6
            )
            assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(36.685066)

    class TestDefinedStateFalse(object):
        @pytest.fixture(scope="class")
        def frame(self):
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            m.fs.prec_sol = MixedAcidParameterBlock()

            m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time, defined_state=False)

            return m

        @pytest.mark.unit
        def test_build(self, frame):
            m = frame

            assert len(m.fs.state) == 1

            assert isinstance(m.fs.state[0].flow_vol, Var)
            assert isinstance(m.fs.state[0].conc_mass_comp, Var)
            assert isinstance(m.fs.state[0].flow_mol_comp, Var)
            assert not m.fs.state[0].is_property_constructed("pH_phase")

            assert isinstance(m.fs.state[0].flow_mol_comp_eqn, Constraint)
            assert isinstance(m.fs.state[0].conc_mol_comp_eqn, Constraint)
            assert isinstance(m.fs.state[0].h2o_concentration_eqn, Constraint)
            assert m.fs.state[0].h2o_concentration_eqn.active
            assert not m.fs.state[0].is_property_constructed("pH_phase_eqn")

            m.fs.state[0].flow_vol.set_value(10)
            for i in m.fs.prec_sol.component_list:
                m.fs.state[0].conc_mass_comp[i].set_value(0.1)
                m.fs.state[0].flow_mol_comp[i].set_value(0.5)

            m.fs.state[0].pH_phase["liquid"].set_value(1)
            assert isinstance(m.fs.state[0].pH_phase_eqn, Constraint)

            m.fs.state.fix_initialization_states()

            assert degrees_of_freedom(m.fs.state[0]) == 0

            assert m.fs.state[0].flow_vol.fixed
            for j in m.fs.prec_sol.component_list:
                if j == "H2O":
                    assert not m.fs.state[0].conc_mass_comp[j].fixed
                else:
                    assert m.fs.state[0].conc_mass_comp[j].fixed
                assert not m.fs.state[0].flow_mol_comp[j].fixed

                assert m.fs.state[0].flow_mol_comp_eqn[j].active
            assert m.fs.state[0].h2o_concentration_eqn.active

        @pytest.mark.unit
        def test_scaling(self, frame):
            m = frame
            assert m.fs.state[0].default_scaler is MixedAcidPropertiesScaler
            scaler_obj = m.fs.state[0].default_scaler()

            with pytest.raises(
                ValueError,
                match=re.escape(
                    "This scaler requires the user to provide a default scaling "
                    "factor for fs.state[0.0].flow_vol, but no default scaling "
                    "factor was set."
                ),
            ):
                scaler_obj.scale_model(m.fs.state[0])

            scaler_obj.default_scaling_factors["flow_vol"] = 0.1
            for ree in ree_list:
                scaler_obj.default_scaling_factors[f"conc_mass_comp[{ree}]"] = 0.1
            for contaminant in contaminant_list:
                scaler_obj.default_scaling_factors[f"conc_mass_comp[{contaminant}]"] = (
                    0.1
                )

            scaler_obj.scale_model(m.fs.state[0])

            assert len(list_unscaled_variables(m.fs.state[0])) == 0
            assert len(list_unscaled_constraints(m.fs.state[0])) == 0

            assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(
                2.436490e6
            )
            assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(36.25460)

            m.fs.state[0].conc_mass_comp.unfix()

            assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(
                2.436474e6
            )
            assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(
                40.1640625
            )


class TestSulfuricAcidLeaching(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.leach_soln = MixedAcidParameterBlock(include_sulfates=True)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert len(model.fs.leach_soln.phase_list) == 1
        for k in model.fs.leach_soln.phase_list:
            assert k == "liquid"

        for k in model.fs.leach_soln.component_list:
            assert k in [
                "H2O",
                "H",
                "HSO4",
                "SO4",
                "Sc",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Gd",
                "Dy",
                "Al",
                "Ca",
                "Fe",
                "Cl",
            ]
            assert k in model.fs.leach_soln.mw

        assert model.fs.leach_soln._has_inherent_reactions

        assert isinstance(model.fs.leach_soln.inherent_reaction_idx, Set)

        assert model.fs.leach_soln.inherent_reaction_stoichiometry == {
            ("H2SO4_Ka2", "liquid", "H"): 1,
            ("H2SO4_Ka2", "liquid", "HSO4"): -1,
            ("H2SO4_Ka2", "liquid", "SO4"): 1,
            ("H2SO4_Ka2", "liquid", "H2O"): 0,
            ("H2SO4_Ka2", "liquid", "Sc"): 0,
            ("H2SO4_Ka2", "liquid", "Y"): 0,
            ("H2SO4_Ka2", "liquid", "La"): 0,
            ("H2SO4_Ka2", "liquid", "Ce"): 0,
            ("H2SO4_Ka2", "liquid", "Pr"): 0,
            ("H2SO4_Ka2", "liquid", "Nd"): 0,
            ("H2SO4_Ka2", "liquid", "Sm"): 0,
            ("H2SO4_Ka2", "liquid", "Gd"): 0,
            ("H2SO4_Ka2", "liquid", "Dy"): 0,
            ("H2SO4_Ka2", "liquid", "Al"): 0,
            ("H2SO4_Ka2", "liquid", "Ca"): 0,
            ("H2SO4_Ka2", "liquid", "Fe"): 0,
            ("H2SO4_Ka2", "liquid", "Cl"): 0,
        }

        assert isinstance(model.fs.leach_soln.k_eq, Param)
        assert len(model.fs.leach_soln.k_eq) == 1
        assert value(model.fs.leach_soln.k_eq["H2SO4_Ka2"]) == pytest.approx(
            10**-1.99, rel=1e-8
        )

        assert isinstance(model.fs.leach_soln.dens_mass, Param)
        assert value(model.fs.leach_soln.dens_mass) == pytest.approx(1, rel=1e-8)

    @pytest.mark.unit
    def test_build_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        assert len(model.fs.state) == 1

        assert isinstance(model.fs.state[0].flow_vol, Var)
        assert isinstance(model.fs.state[0].conc_mass_comp, Var)
        assert isinstance(model.fs.state[0].conc_mol_comp, Var)

        assert isinstance(model.fs.state[0].conc_mol_comp_eqn, Constraint)
        assert isinstance(model.fs.state[0].h2o_concentration_eqn, Constraint)
        assert isinstance(model.fs.state[0].inherent_equilibrium_eqn, Constraint)
        assert len(model.fs.state[0].inherent_equilibrium_eqn) == 1

        assert isinstance(model.fs.state[0].dens_mass, Param)

    @pytest.mark.unit
    def test_fix_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        model.fs.state[0].flow_vol.set_value(10)
        model.fs.state[0].conc_mass_comp[:].set_value(1)
        model.fs.state[0].conc_mol_comp[:].set_value(1)

        model.fs.state.fix_initialization_states()

        assert model.fs.state[0].flow_vol.fixed
        for j in model.fs.leach_soln.component_list:
            if j == "H2O" or j == "HSO4":
                assert not model.fs.state[0].conc_mass_comp[j].fixed
            else:
                assert model.fs.state[0].conc_mass_comp[j].fixed
            assert not model.fs.state[0].conc_mol_comp[j].fixed

            assert model.fs.state[0].conc_mol_comp_eqn[j].active

        assert model.fs.state[0].h2o_concentration_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn["H2SO4_Ka2"].active


class TestSulfuricAcidLeaching(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.leach_soln = MixedAcidParameterBlock(include_sulfates=True)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert len(model.fs.leach_soln.phase_list) == 1
        for k in model.fs.leach_soln.phase_list:
            assert k == "liquid"

        for k in model.fs.leach_soln.component_list:
            assert k in [
                "H2O",
                "H",
                "HSO4",
                "SO4",
                "Sc",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Gd",
                "Dy",
                "Al",
                "Ca",
                "Fe",
                "Cl",
            ]
            assert k in model.fs.leach_soln.mw

        assert model.fs.leach_soln._has_inherent_reactions

        assert isinstance(model.fs.leach_soln.inherent_reaction_idx, Set)

        assert model.fs.leach_soln.inherent_reaction_stoichiometry == {
            ("H2SO4_Ka2", "liquid", "H"): 1,
            ("H2SO4_Ka2", "liquid", "HSO4"): -1,
            ("H2SO4_Ka2", "liquid", "SO4"): 1,
            ("H2SO4_Ka2", "liquid", "H2O"): 0,
            ("H2SO4_Ka2", "liquid", "Sc"): 0,
            ("H2SO4_Ka2", "liquid", "Y"): 0,
            ("H2SO4_Ka2", "liquid", "La"): 0,
            ("H2SO4_Ka2", "liquid", "Ce"): 0,
            ("H2SO4_Ka2", "liquid", "Pr"): 0,
            ("H2SO4_Ka2", "liquid", "Nd"): 0,
            ("H2SO4_Ka2", "liquid", "Sm"): 0,
            ("H2SO4_Ka2", "liquid", "Gd"): 0,
            ("H2SO4_Ka2", "liquid", "Dy"): 0,
            ("H2SO4_Ka2", "liquid", "Al"): 0,
            ("H2SO4_Ka2", "liquid", "Ca"): 0,
            ("H2SO4_Ka2", "liquid", "Fe"): 0,
            ("H2SO4_Ka2", "liquid", "Cl"): 0,
        }

        assert isinstance(model.fs.leach_soln.k_eq, Param)
        assert len(model.fs.leach_soln.k_eq) == 1
        assert value(model.fs.leach_soln.k_eq["H2SO4_Ka2"]) == pytest.approx(
            10**-1.99, rel=1e-8
        )

        assert isinstance(model.fs.leach_soln.dens_mass, Param)
        assert value(model.fs.leach_soln.dens_mass) == pytest.approx(1, rel=1e-8)

    @pytest.mark.unit
    def test_build_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        assert len(model.fs.state) == 1

        assert isinstance(model.fs.state[0].flow_vol, Var)
        assert isinstance(model.fs.state[0].conc_mass_comp, Var)
        assert isinstance(model.fs.state[0].conc_mol_comp, Var)

        assert isinstance(model.fs.state[0].conc_mol_comp_eqn, Constraint)
        assert isinstance(model.fs.state[0].h2o_concentration_eqn, Constraint)
        assert isinstance(model.fs.state[0].inherent_equilibrium_eqn, Constraint)
        assert len(model.fs.state[0].inherent_equilibrium_eqn) == 1

        assert isinstance(model.fs.state[0].dens_mass, Param)

    @pytest.mark.unit
    def test_fix_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        model.fs.state[0].flow_vol.set_value(10)
        model.fs.state[0].conc_mass_comp[:].set_value(1)
        model.fs.state[0].conc_mol_comp[:].set_value(1)

        model.fs.state.fix_initialization_states()

        assert model.fs.state[0].flow_vol.fixed
        for j in model.fs.leach_soln.component_list:
            if j == "H2O" or j == "HSO4":
                assert not model.fs.state[0].conc_mass_comp[j].fixed
            else:
                assert model.fs.state[0].conc_mass_comp[j].fixed
            assert not model.fs.state[0].conc_mol_comp[j].fixed

            assert model.fs.state[0].conc_mol_comp_eqn[j].active

        assert model.fs.state[0].h2o_concentration_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn["H2SO4_Ka2"].active
        assert str(model.fs.state[0].inherent_equilibrium_eqn["H2SO4_Ka2"].body) == (
            "fs.leach_soln.k_eq[H2SO4_Ka2]*fs.state[0.0].conc_mol_comp[HSO4] "
            "- fs.state[0.0].conc_mol_comp[H]*fs.state[0.0].conc_mol_comp[SO4]"
        )


class TestOxalicAcid(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.leach_soln = MixedAcidParameterBlock(include_oxalates=True)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert len(model.fs.leach_soln.phase_list) == 1
        for k in model.fs.leach_soln.phase_list:
            assert k == "liquid"

        for k in model.fs.leach_soln.component_list:
            assert k in [
                "H2O",
                "H",
                "H2C2O4",
                "HC2O4",
                "C2O4",
                "Sc",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Gd",
                "Dy",
                "Al",
                "Ca",
                "Fe",
                "Cl",
            ]
            assert k in model.fs.leach_soln.mw

        assert model.fs.leach_soln._has_inherent_reactions

        assert isinstance(model.fs.leach_soln.inherent_reaction_idx, Set)

        rxn_set = {"H2C2O4_Ka1", "H2C2O4_Ka2"}
        for (
            r,
            _,
            j,
        ), nu_j in model.fs.leach_soln.inherent_reaction_stoichiometry.items():
            assert r in rxn_set
            if r == "H2C2O4_Ka1":
                if j == "H" or j == "HC2O4":
                    assert nu_j == 1
                elif j == "H2C2O4":
                    assert nu_j == -1
                else:
                    assert nu_j == 0
            else:
                if j == "H" or j == "C2O4":
                    assert nu_j == 1
                elif j == "HC2O4":
                    assert nu_j == -1
                else:
                    assert nu_j == 0

        assert isinstance(model.fs.leach_soln.k_eq, Param)
        assert len(model.fs.leach_soln.k_eq) == 2
        assert value(model.fs.leach_soln.k_eq["H2C2O4_Ka1"]) == pytest.approx(
            5.90e-2, rel=1e-8
        )
        assert value(model.fs.leach_soln.k_eq["H2C2O4_Ka2"]) == pytest.approx(
            6.40e-5, rel=1e-8
        )

        assert isinstance(model.fs.leach_soln.dens_mass, Param)
        assert value(model.fs.leach_soln.dens_mass) == pytest.approx(1, rel=1e-8)

    @pytest.mark.unit
    def test_build_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        assert len(model.fs.state) == 1

        assert isinstance(model.fs.state[0].flow_vol, Var)
        assert isinstance(model.fs.state[0].conc_mass_comp, Var)
        assert isinstance(model.fs.state[0].conc_mol_comp, Var)

        assert isinstance(model.fs.state[0].conc_mol_comp_eqn, Constraint)
        assert isinstance(model.fs.state[0].h2o_concentration_eqn, Constraint)
        assert isinstance(model.fs.state[0].inherent_equilibrium_eqn, Constraint)
        assert len(model.fs.state[0].inherent_equilibrium_eqn) == 2

        assert isinstance(model.fs.state[0].dens_mass, Param)

    @pytest.mark.unit
    def test_fix_state(self, model):
        model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

        model.fs.state[0].flow_vol.set_value(10)
        model.fs.state[0].conc_mass_comp[:].set_value(1)
        model.fs.state[0].conc_mol_comp[:].set_value(1)

        model.fs.state.fix_initialization_states()

        assert model.fs.state[0].flow_vol.fixed
        for j in model.fs.leach_soln.component_list:
            if j in {"H2O", "HC2O4", "C2O4"}:
                assert not model.fs.state[0].conc_mass_comp[j].fixed
            else:
                assert model.fs.state[0].conc_mass_comp[j].fixed
            assert not model.fs.state[0].conc_mol_comp[j].fixed

            assert model.fs.state[0].conc_mol_comp_eqn[j].active

        assert model.fs.state[0].h2o_concentration_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn.active
        assert model.fs.state[0].inherent_equilibrium_eqn["H2C2O4_Ka1"].active
        assert str(model.fs.state[0].inherent_equilibrium_eqn["H2C2O4_Ka1"].body) == (
            "fs.leach_soln.k_eq[H2C2O4_Ka1]*fs.state[0.0].conc_mol_comp[H2C2O4] "
            "- fs.state[0.0].conc_mol_comp[H]*fs.state[0.0].conc_mol_comp[HC2O4]"
        )
        assert model.fs.state[0].inherent_equilibrium_eqn["H2C2O4_Ka2"].active
        assert str(model.fs.state[0].inherent_equilibrium_eqn["H2C2O4_Ka2"].body) == (
            "fs.leach_soln.k_eq[H2C2O4_Ka2]*fs.state[0.0].conc_mol_comp[HC2O4] "
            "- fs.state[0.0].conc_mol_comp[H]*fs.state[0.0].conc_mol_comp[C2O4]"
        )
