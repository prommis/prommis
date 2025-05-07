#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.nanofiltration.membrane_cascade_flowsheet.solute_property import (
    SoluteParameters,
)

# -----------------------------------------------------------------------------
# Test settings

# solutes
solutes = ["Li", "Co"]

# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = SoluteParameters(solutes=solutes)

    m.fs.state = m.fs.properties.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_vol, Var)
    assert isinstance(m.fs.state[0].flow_mass_solute, Var)

    m.fs.state[0].flow_vol.set_value(100)
    m.fs.state[0].flow_mass_solute[:].set_value(100)

    m.fs.state.fix_initialization_states()

    for j in m.fs.properties.component_list:
        if j == "solvent":
            assert m.fs.state[0].flow_vol.fixed
        else:
            assert m.fs.state[0].flow_mass_solute[j].fixed
