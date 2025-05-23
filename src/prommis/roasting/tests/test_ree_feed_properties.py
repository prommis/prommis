#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Param, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.roasting.ree_feed_properties import ReeFeedParameters


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.feed_prop = ReeFeedParameters()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.feed_prop.phase_list) == 1
    for k in model.fs.feed_prop.phase_list:
        assert k == "solid"

    for k in model.fs.feed_prop.component_list:
        assert k in [
            "C",
            "H",
            "O",
            "N",
            "S",
            "H2O",
            "Kaolinite",
            "Al2O3",
            "SiO2",
            "CaCO3",
            "CaO",
            "FeS2",
            "Fe2O3",
            "Ree2X",
            "Sc2X",
            "Y2X",
            "La2X",
            "Ce2X",
            "Pr2X",
            "Nd2X",
            "Sm2X",
            "Gd2X",
            "Dy2X",
            "Sc2O3",
            "Y2O3",
            "La2O3",
            "Ce2O3",
            "Pr2O3",
            "Nd2O3",
            "Sm2O3",
            "Gd2O3",
            "Dy2O3",
        ]
        assert k in model.fs.feed_prop.mw_comp
        assert k in model.fs.feed_prop.enth0_comp
        assert k in model.fs.feed_prop.cp0_comp
        assert k in model.fs.feed_prop.cp1_comp

    assert isinstance(model.fs.feed_prop.dens_mass, Param)
    assert value(model.fs.feed_prop.dens_mass) == pytest.approx(2000, rel=1e-8)


@pytest.mark.unit
def test_build_state(model):
    model.fs.state = model.fs.feed_prop.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_mass, Var)
    assert isinstance(model.fs.state[0].temperature, Var)
    assert isinstance(model.fs.state[0].mass_frac_comp, Var)
    assert isinstance(model.fs.state[0].flow_mol_comp, Var)
    assert isinstance(model.fs.state[0].enth_mol_comp, Var)
    assert isinstance(model.fs.state[0].enth_mol, Var)
    assert isinstance(model.fs.state[0].enth_mass, Var)

    assert isinstance(model.fs.state[0].flow_mol_comp_constraint, Constraint)
    assert isinstance(model.fs.state[0].sum_mass_frac, Constraint)
    assert isinstance(model.fs.state[0].enth_mol_comp_constraint, Constraint)
    assert isinstance(model.fs.state[0].enth_mol_constraint, Constraint)
    assert isinstance(model.fs.state[0].enth_mass_constraint, Constraint)


@pytest.mark.unit
def test_fix_state(model):
    model.fs.state = model.fs.feed_prop.build_state_block(model.fs.time)

    model.fs.state[0].flow_mass.set_value(10)
    model.fs.state[0].temperature.set_value(900)
    model.fs.state[0].mass_frac_comp[:].set_value(0.02)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_mass.fixed
    assert model.fs.state[0].temperature.fixed
    for j in model.fs.feed_prop.component_list:
        assert model.fs.state[0].mass_frac_comp[j].fixed
        assert model.fs.state[0].flow_mol_comp_constraint[j].active
        assert not model.fs.state[0].sum_mass_frac.active
        assert model.fs.state[0].enth_mol_comp_constraint[j].active
    assert model.fs.state[0].enth_mol_constraint.active
    assert model.fs.state[0].enth_mass_constraint.active
