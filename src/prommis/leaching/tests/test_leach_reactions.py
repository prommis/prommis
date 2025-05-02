#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import Block, ConcreteModel, Constraint, Var, units
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

import pytest

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions

RXN_LIST = [
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


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Dummy blocks for liquid and solid states
    m.fs.liquid = Block(m.fs.time)
    m.fs.liquid[0].flow_vol = Var(units=units.liter / units.hour)
    m.fs.liquid[0].conc_mol_comp = Var(["H"], units=units.mol / units.liter)

    m.fs.solid = Block(m.fs.time)
    m.fs.solid[0].flow_mass = Var(units=units.kg / units.hour)
    m.fs.solid[0].conversion = Var(RXN_LIST, units=units.dimensionless)

    m.fs.solid[0].params = Block()
    m.fs.solid[0].params.dens_mass = Var(units=units.kg / units.liter)

    # Leaching reaction parameters
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.leach_rxns.reaction_idx) == 12
    for k in model.fs.leach_rxns.reaction_idx:
        assert k in RXN_LIST
        assert k in model.fs.leach_rxns.A
        assert k in model.fs.leach_rxns.B

    assert isinstance(model.fs.leach_rxns.reaction_stoichiometry, dict)


@pytest.mark.unit
def test_build_reaction_block(model):
    model.fs.rxns = model.fs.leach_rxns.build_reaction_block(model.fs.time)

    assert len(model.fs.rxns) == 1

    assert isinstance(model.fs.rxns[0].reaction_rate, Var)
    assert len(model.fs.rxns[0].reaction_rate) == 12
    for k in model.fs.rxns[0].reaction_rate:
        assert k in RXN_LIST

    assert isinstance(model.fs.rxns[0].reaction_rate_eq, Constraint)
    assert len(model.fs.rxns[0].reaction_rate_eq) == 12
    for k in model.fs.rxns[0].reaction_rate_eq:
        assert k in RXN_LIST


@pytest.mark.unit
def test_unit_consistency(model):
    model.fs.rxns = model.fs.leach_rxns.build_reaction_block(model.fs.time)

    assert_units_consistent(model)
