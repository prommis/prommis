from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_solid = PrecipitateParameters()

    m.fs.state = m.fs.prec_solid.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_mol_comp, Var)

    m.fs.state[0].flow_mol_comp[:].set_value(0.5)

    m.fs.state.fix_initialization_states()

    for j in m.fs.prec_solid.component_list:
        assert m.fs.state[0].flow_mol_comp[j].fixed
