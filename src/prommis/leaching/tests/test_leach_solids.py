from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_solids_properties import CoalRefuseParameters


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_sol = CoalRefuseParameters()

    m.fs.state = m.fs.leach_sol.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_mass, Var)
    assert isinstance(m.fs.state[0].mass_frac_comp, Var)
    assert isinstance(m.fs.state[0].conversion, Var)

    assert isinstance(m.fs.state[0].conversion_eq, Constraint)
    assert isinstance(m.fs.state[0].sum_mass_frac, Constraint)

    m.fs.state[0].flow_mass.set_value(10)
    m.fs.state[0].mass_frac_comp[:].set_value(0.1)
    m.fs.state[0].conversion[:].set_value(0.5)

    m.fs.state.fix_initialization_states()

    assert m.fs.state[0].flow_mass.fixed
    for j in m.fs.leach_sol.component_list:
        assert m.fs.state[0].mass_frac_comp[j].fixed
        assert not m.fs.state[0].conversion[j].fixed

        assert m.fs.state[0].conversion_eq[j].active

        assert not m.fs.state[0].sum_mass_frac.active
