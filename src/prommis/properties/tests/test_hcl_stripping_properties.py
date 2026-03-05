#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re

from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core import FlowsheetBlock
from idaes.core.scaling.util import (
    jacobian_cond,
    list_unscaled_constraints,
    list_unscaled_variables,
)
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.properties.hcl_stripping_properties import (
    HClStrippingParameterBlock,
    HClStrippingPropertiesScaler,
)


class TestDefinedStateTrue(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.prec_sol = HClStrippingParameterBlock()

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
        for i in m.fs.prec_sol.dissolved_elements:
            m.fs.state[0].conc_mass_comp[i].set_value(0.1)
            m.fs.state[0].flow_mol_comp[i].set_value(0.5)

        m.fs.state[0].pH_phase["liquid"].set_value(1)
        assert isinstance(m.fs.state[0].pH_phase_eqn, Constraint)

        m.fs.state.fix_initialization_states()

        assert degrees_of_freedom(m.fs.state[0]) == 0

        assert m.fs.state[0].flow_vol.fixed
        for j in m.fs.prec_sol.dissolved_elements:
            assert m.fs.state[0].conc_mass_comp[j].fixed
            assert not m.fs.state[0].flow_mol_comp[j].fixed

            assert m.fs.state[0].flow_mol_comp_eqn[j].active

    @pytest.mark.unit
    def test_scaling(self, frame):
        m = frame
        assert m.fs.state[0].default_scaler is HClStrippingPropertiesScaler
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
        scaler_obj.scale_model(m.fs.state[0])

        assert len(list_unscaled_variables(m.fs.state[0])) == 0
        assert len(list_unscaled_constraints(m.fs.state[0])) == 0

        assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(2.401286e6)
        assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(33.24336)

        m.fs.state[0].conc_mass_comp.unfix()

        assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(2.401271e6)
        assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(36.685066)


class TestDefinedStateFalse(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.prec_sol = HClStrippingParameterBlock()

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
        for i in m.fs.prec_sol.dissolved_elements:
            m.fs.state[0].conc_mass_comp[i].set_value(0.1)
            m.fs.state[0].flow_mol_comp[i].set_value(0.5)

        m.fs.state[0].pH_phase["liquid"].set_value(1)
        assert isinstance(m.fs.state[0].pH_phase_eqn, Constraint)

        m.fs.state.fix_initialization_states()

        assert degrees_of_freedom(m.fs.state[0]) == 0

        assert m.fs.state[0].flow_vol.fixed
        for j in m.fs.prec_sol.dissolved_elements:
            assert m.fs.state[0].conc_mass_comp[j].fixed
            assert not m.fs.state[0].flow_mol_comp[j].fixed

            assert m.fs.state[0].flow_mol_comp_eqn[j].active
        assert not m.fs.state[0].h2o_concentration_eqn.active

    @pytest.mark.unit
    def test_scaling(self, frame):
        m = frame
        assert m.fs.state[0].default_scaler is HClStrippingPropertiesScaler
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
        scaler_obj.scale_model(m.fs.state[0])

        assert len(list_unscaled_variables(m.fs.state[0])) == 0
        assert len(list_unscaled_constraints(m.fs.state[0])) == 0

        assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(2.401286e6)
        assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(33.24336)

        m.fs.state[0].conc_mass_comp["H2O"].unfix()
        m.fs.state[0].h2o_concentration_eqn.activate()

        assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(2.436486e6)
        assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(36.25460)

        m.fs.state[0].conc_mass_comp.unfix()

        assert jacobian_cond(m.fs.state[0], scaled=False) == pytest.approx(2.436470e6)
        assert jacobian_cond(m.fs.state[0], scaled=True) == pytest.approx(40.1640625)
