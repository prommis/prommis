#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.hydrogen_decrepitation.repm_solids_properties import REPMParameters


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.repm_sol = REPMParameters()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.repm_sol.phase_list) == 1
    for k in model.fs.repm_sol.phase_list:
        assert k == "solid"

    for k in model.fs.repm_sol.component_list:
        assert k in [
            "Nd",
            "Nd2Fe14B",
        ]
        assert k in model.fs.repm_sol.mw

    # assert isinstance(model.fs.repm_sol.dens_mass, Param)
    # assert value(model.fs.repm_sol.dens_mass) == pytest.approx(7500, rel=1e-8)


@pytest.mark.unit
def test_build_state(model):
    model.fs.state = model.fs.repm_sol.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_mass, Var)
    assert isinstance(model.fs.state[0].mass_frac_comp, Var)

    assert isinstance(model.fs.state[0].sum_mass_frac, Constraint)


@pytest.mark.unit
def test_fix_state(model):
    model.fs.state = model.fs.repm_sol.build_state_block(model.fs.time)

    model.fs.state[0].flow_mass.set_value(10)
    model.fs.state[0].mass_frac_comp[:].set_value(0.5)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_mass.fixed
    for j in model.fs.repm_sol.component_list:
        assert model.fs.state[0].mass_frac_comp[j].fixed

        assert not model.fs.state[0].sum_mass_frac.active
