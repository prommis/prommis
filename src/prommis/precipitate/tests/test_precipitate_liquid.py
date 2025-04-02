#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.precipitate.precipitate_liquid_properties import AqueousParameter


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_sol = AqueousParameter()

    m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time)

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
