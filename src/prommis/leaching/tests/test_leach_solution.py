from pyomo.environ import ConcreteModel, Constraint, Param, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_solution_properties import LeachSolutionParameters


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = LeachSolutionParameters()

    m.fs.state = m.fs.leach_soln.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_vol, Var)
    assert isinstance(m.fs.state[0].conc_mass_comp, Var)
    assert isinstance(m.fs.state[0].conc_mol_comp, Var)

    assert isinstance(m.fs.state[0].molar_concentration_constraint, Constraint)
    assert isinstance(m.fs.state[0].h2o_concentration, Constraint)
    assert isinstance(m.fs.state[0].hso4_dissociation, Constraint)

    assert isinstance(m.fs.state[0].dens_mol, Param)

    m.fs.state[0].flow_vol.set_value(10)
    m.fs.state[0].conc_mass_comp[:].set_value(1)
    m.fs.state[0].conc_mol_comp[:].set_value(1)

    m.fs.state.fix_initialization_states()

    assert m.fs.state[0].flow_vol.fixed
    for j in m.fs.leach_soln.component_list:
        assert m.fs.state[0].conc_mass_comp[j].fixed
        assert not m.fs.state[0].conc_mol_comp[j].fixed

        assert m.fs.state[0].molar_concentration_constraint[j].active

    assert not m.fs.state[0].h2o_concentration.active
    assert not m.fs.state[0].hso4_dissociation.active
