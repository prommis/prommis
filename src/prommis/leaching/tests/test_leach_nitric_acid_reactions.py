#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re
import pytest

from pyomo.environ import Block, ConcreteModel, Expression, Var, units, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from prommis.leaching.leach_nitric_acid_reactions import (
    NitricAcidLeachingReactionParameterBlock,
    _default_aqueous_aliases,
    _metal_list,
)

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

EXPECTED_A = 0.37
EXPECTED_B = {
    "Al2O3": 0.183,
    "Fe2O3": 0.0440,
    "CaO": 1.160,
    "Sc2O3": 2.64e-5,
    "Y2O3": 2.23e-5,
    "La2O3": 2.86e-5,
    "Ce2O3": 3.95e-5,
    "Pr2O3": 1.19e-6,
    "Nd2O3": 2.52e-5,
    "Sm2O3": 4.68e-6,
    "Gd2O3": 3.07e-6,
    "Dy2O3": 2.29e-6,
}


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
    m.fs.solid[0].conversion_comp = Var(RXN_LIST, units=units.dimensionless)

    m.fs.solid[0].params = Block()
    m.fs.solid[0].params.dens_mass = Var(units=units.kg / units.liter)

    # Leaching reaction parameters
    m.fs.leach_rxns = NitricAcidLeachingReactionParameterBlock()

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
def test_parameter_values(model):
    for k in RXN_LIST:
        assert value(model.fs.leach_rxns.A[k]) == pytest.approx(EXPECTED_A, rel=1e-8)
        assert value(model.fs.leach_rxns.B[k]) == pytest.approx(EXPECTED_B[k], rel=1e-8)


@pytest.mark.unit
def test_stoichiometry_all_oxides(model):
    stoich = model.fs.leach_rxns.reaction_stoichiometry

    for metal in _metal_list:
        if metal == "Ca":
            oxide, n_metal, n_water, n_acid = "CaO", 1, 1, -2
        else:
            oxide, n_metal, n_water, n_acid = f"{metal}2O3", 2, 3, -6

        assert stoich[(oxide, "liquid", metal)] == n_metal
        assert stoich[(oxide, "liquid", "H2O")] == n_water
        assert stoich[(oxide, "liquid", "H")] == n_acid
        assert stoich[(oxide, "solid", oxide)] == -1


@pytest.mark.unit
def test_build_reaction_block(model):
    model.fs.rxns = model.fs.leach_rxns.build_reaction_block(model.fs.time)

    assert len(model.fs.rxns) == 1

    assert isinstance(model.fs.rxns[0].reaction_rate, Expression)
    assert len(model.fs.rxns[0].reaction_rate) == 12
    for k in model.fs.rxns[0].reaction_rate:
        assert k in RXN_LIST


@pytest.mark.unit
def test_unit_consistency(model):
    model.fs.rxns = model.fs.leach_rxns.build_reaction_block(model.fs.time)

    assert_units_consistent(model)


@pytest.mark.unit
def test_aqueous_aliases_default_config():
    """Verify default values and type of aqueous_aliases config argument."""
    m = ConcreteModel()
    m.leach_rxns = NitricAcidLeachingReactionParameterBlock()

    cfg = m.leach_rxns.config.aqueous_aliases

    # Must be a dict and match the module default
    assert isinstance(cfg, dict)
    assert cfg == _default_aqueous_aliases

    # Check that H and H2O defaults are present
    assert cfg["H"] == "H"
    assert cfg["H2O"] == "H2O"

    # Check that all metals in _metal_list default to identity mapping
    for metal in _metal_list:
        assert metal in cfg
        assert cfg[metal] == metal


@pytest.mark.unit
def test_aqueous_aliases_default_stoichiometry():
    """Verify stoichiometry keys use default names when no aliases provided."""
    m = ConcreteModel()
    m.leach_rxns = NitricAcidLeachingReactionParameterBlock()

    stoich = m.leach_rxns.reaction_stoichiometry

    # Verify H and H2O mapping in liquid phase
    assert stoich[("CaO", "liquid", "H")] == -2
    assert stoich[("CaO", "liquid", "H2O")] == 1
    assert stoich[("Al2O3", "liquid", "H")] == -6
    assert stoich[("Al2O3", "liquid", "H2O")] == 3

    # Verify metals mapped to default names
    assert stoich[("CaO", "liquid", "Ca")] == 1
    for metal in _metal_list:
        oxide = "CaO" if metal == "Ca" else f"{metal}2O3"
        assert (oxide, "liquid", metal) in stoich


@pytest.mark.unit
def test_aqueous_aliases_custom_config():
    """Verify custom aliases (ionic names, as used by other property packages)."""
    custom_aliases = {"H": "H_+", "H2O": "H2O"}
    for metal in _metal_list:
        charge = "2+" if metal == "Ca" else "3+"
        custom_aliases[metal] = f"{metal}_{charge}"

    m = ConcreteModel()
    m.leach_rxns = NitricAcidLeachingReactionParameterBlock(
        aqueous_aliases=custom_aliases
    )

    stoich = m.leach_rxns.reaction_stoichiometry

    # Verify H was remapped to H_+ and H2O remains H2O
    assert ("CaO", "liquid", "H_+") in stoich
    assert stoich[("CaO", "liquid", "H_+")] == -2
    assert ("CaO", "liquid", "H") not in stoich
    assert stoich[("CaO", "liquid", "H2O")] == 1

    # Verify Ca mapped to Ca_2+
    assert stoich[("CaO", "liquid", "Ca_2+")] == 1
    assert ("CaO", "liquid", "Ca") not in stoich

    # Verify trivalent REEs and impurities mapped to their respective ionic charges
    assert stoich[("Al2O3", "liquid", "Al_3+")] == 2
    assert stoich[("Fe2O3", "liquid", "Fe_3+")] == 2
    assert stoich[("La2O3", "liquid", "La_3+")] == 2
    assert stoich[("Sc2O3", "liquid", "Sc_3+")] == 2

    # Solid phase remains unchanged
    assert stoich[("CaO", "solid", "CaO")] == -1
    assert stoich[("Al2O3", "solid", "Al2O3")] == -1


@pytest.mark.unit
def test_aqueous_aliases_invalid_type():
    """Verify config validation fails when aqueous_aliases is not a dict."""
    m = ConcreteModel()
    with pytest.raises(
        ValueError, match=re.escape("invalid value for configuration 'aqueous_aliases'")
    ):
        m.leach_rxns = NitricAcidLeachingReactionParameterBlock(
            aqueous_aliases=["H", "H2O"]
        )


@pytest.mark.unit
def test_aqueous_aliases_missing_required_key():
    """Verify KeyError is raised if a required element key is omitted."""
    incomplete_aliases = {
        "H": "H_+",
        "H2O": "H2O",
        # missing Ca and metal list entries
    }

    m = ConcreteModel()
    with pytest.raises(KeyError, match=re.escape("La")):
        m.leach_rxns = NitricAcidLeachingReactionParameterBlock(
            aqueous_aliases=incomplete_aliases
        )
