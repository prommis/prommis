#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re

from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.properties.hcl_stripping_properties import (
    HClStrippingParameterBlock,
    HClStrippingPropertiesScaler,
)

@pytest.fixture
def frame():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_sol = HClStrippingParameterBlock()

    m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time)
    
    return m

@pytest.mark.unit
def test_build(frame):
    m = frame

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_vol, Var)
    assert isinstance(m.fs.state[0].conc_mass_comp, Var)
    assert isinstance(m.fs.state[0].flow_mol_comp, Var)

    assert isinstance(m.fs.state[0].flow_mol_constraint, Constraint)

    m.fs.state[0].flow_vol.set_value(10)
    for i in m.fs.prec_sol.dissolved_elements:
        m.fs.state[0].conc_mass_comp[i].set_value(0.1)
        m.fs.state[0].flow_mol_comp[i].set_value(0.5)

    m.fs.state.fix_initialization_states()

    assert m.fs.state[0].flow_vol.fixed
    for j in m.fs.prec_sol.dissolved_elements:
        assert m.fs.state[0].conc_mass_comp[j].fixed
        assert not m.fs.state[0].flow_mol_comp[j].fixed

        assert m.fs.state[0].flow_mol_constraint[j].active

@pytest.mark.unit
def test_scaling(frame):
    m = frame
    assert m.fs.state[0].default_scaler is HClStrippingPropertiesScaler
    with pytest.raises(
        AssertionError,
        match=re.escape("foo")
    ):
        HClStrippingPropertiesScaler.scale_model(m)
