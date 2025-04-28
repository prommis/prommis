#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Param, Set, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_solution_properties import LeachSolutionParameters


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = LeachSolutionParameters()

    return m


@pytest.mark.unit
def test_parameters(model):
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
        ("Ka2", "liquid", "H"): 1,
        ("Ka2", "liquid", "HSO4"): -1,
        ("Ka2", "liquid", "SO4"): 1,
        ("Ka2", "liquid", "H2O"): 0,
        ("Ka2", "liquid", "Sc"): 0,
        ("Ka2", "liquid", "Y"): 0,
        ("Ka2", "liquid", "La"): 0,
        ("Ka2", "liquid", "Ce"): 0,
        ("Ka2", "liquid", "Pr"): 0,
        ("Ka2", "liquid", "Nd"): 0,
        ("Ka2", "liquid", "Sm"): 0,
        ("Ka2", "liquid", "Gd"): 0,
        ("Ka2", "liquid", "Dy"): 0,
        ("Ka2", "liquid", "Al"): 0,
        ("Ka2", "liquid", "Ca"): 0,
        ("Ka2", "liquid", "Fe"): 0,
        ("Ka2", "liquid", "Cl"): 0,
    }

    assert isinstance(model.fs.leach_soln.Ka2, Param)
    assert value(model.fs.leach_soln.Ka2) == pytest.approx(10**-1.99, rel=1e-8)

    assert isinstance(model.fs.leach_soln.dens_mass, Param)
    assert value(model.fs.leach_soln.dens_mass) == pytest.approx(1, rel=1e-8)


@pytest.mark.unit
def test_build_state(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_vol, Var)
    assert isinstance(model.fs.state[0].conc_mass_comp, Var)
    assert isinstance(model.fs.state[0].conc_mol_comp, Var)

    assert isinstance(model.fs.state[0].molar_concentration_constraint, Constraint)
    assert isinstance(model.fs.state[0].h2o_concentration, Constraint)
    assert isinstance(model.fs.state[0].hso4_dissociation, Constraint)

    assert isinstance(model.fs.state[0].dens_mass, Param)


@pytest.mark.unit
def test_fix_state(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    model.fs.state[0].flow_vol.set_value(10)
    model.fs.state[0].conc_mass_comp[:].set_value(1)
    model.fs.state[0].conc_mol_comp[:].set_value(1)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_vol.fixed
    for j in model.fs.leach_soln.component_list:
        assert model.fs.state[0].conc_mass_comp[j].fixed
        assert not model.fs.state[0].conc_mol_comp[j].fixed

        assert model.fs.state[0].molar_concentration_constraint[j].active

    assert not model.fs.state[0].h2o_concentration.active
    assert not model.fs.state[0].hso4_dissociation.active
