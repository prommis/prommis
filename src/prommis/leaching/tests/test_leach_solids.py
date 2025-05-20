#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Param, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_solids_properties import CoalRefuseParameters


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_sol = CoalRefuseParameters()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.leach_sol.phase_list) == 1
    for k in model.fs.leach_sol.phase_list:
        assert k == "solid"

    for k in model.fs.leach_sol.component_list:
        assert k in [
            "inerts",
            "Sc2O3",
            "Y2O3",
            "La2O3",
            "Ce2O3",
            "Pr2O3",
            "Nd2O3",
            "Sm2O3",
            "Gd2O3",
            "Dy2O3",
            "Al2O3",
            "CaO",
            "Fe2O3",
        ]
        assert k in model.fs.leach_sol.mw
        assert k in model.fs.leach_sol.mass_frac_comp_initial

    assert isinstance(model.fs.leach_sol.dens_mass, Param)
    assert value(model.fs.leach_sol.dens_mass) == pytest.approx(2.4, rel=1e-8)


@pytest.mark.unit
def test_build_state(model):
    model.fs.state = model.fs.leach_sol.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_mass, Var)
    assert isinstance(model.fs.state[0].mass_frac_comp, Var)
    assert isinstance(model.fs.state[0].conversion, Var)

    assert isinstance(model.fs.state[0].conversion_eq, Constraint)
    assert isinstance(model.fs.state[0].sum_mass_frac, Constraint)


@pytest.mark.unit
def test_fix_state(model):
    model.fs.state = model.fs.leach_sol.build_state_block(model.fs.time)

    model.fs.state[0].flow_mass.set_value(10)
    model.fs.state[0].mass_frac_comp[:].set_value(0.1)
    model.fs.state[0].conversion[:].set_value(0.5)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_mass.fixed
    for j in model.fs.leach_sol.component_list:
        assert model.fs.state[0].mass_frac_comp[j].fixed
        assert not model.fs.state[0].conversion[j].fixed

        assert model.fs.state[0].conversion_eq[j].active

        assert not model.fs.state[0].sum_mass_frac.active
