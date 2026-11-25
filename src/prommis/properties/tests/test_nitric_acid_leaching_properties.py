#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, Param, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.properties.nitric_acid_leaching_properties import (
    NitricAcidLeachingParameters,
)


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = NitricAcidLeachingParameters()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.leach_soln.phase_list) == 1
    for k in model.fs.leach_soln.phase_list:
        assert k == "liquid"

    expected_components = [
        "H2O",
        "H",
        "NO3",
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
    assert len(model.fs.leach_soln.component_list) == len(expected_components)

    for k in model.fs.leach_soln.component_list:
        assert k in expected_components
        assert k in model.fs.leach_soln.mw

    assert isinstance(model.fs.leach_soln.dens_mass, Param)
    assert value(model.fs.leach_soln.dens_mass) == pytest.approx(1, rel=1e-8)

    assert isinstance(model.fs.leach_soln.cp_mol, Param)
    assert value(model.fs.leach_soln.cp_mol) == pytest.approx(75.327, rel=1e-8)

    assert isinstance(model.fs.leach_soln.temperature_ref, Param)
    assert value(model.fs.leach_soln.temperature_ref) == pytest.approx(298.15, rel=1e-8)


@pytest.mark.unit
def test_build_state(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_vol, Var)
    assert isinstance(model.fs.state[0].conc_mass_comp, Var)
    assert isinstance(model.fs.state[0].conc_mol_comp, Var)
    assert isinstance(model.fs.state[0].temperature, Var)
    assert isinstance(model.fs.state[0].pressure, Var)

    assert isinstance(model.fs.state[0].molar_concentration_constraint, Constraint)
    assert isinstance(model.fs.state[0].h2o_concentration, Constraint)

    assert isinstance(model.fs.state[0].dens_mass, Param)


@pytest.mark.unit
def test_build_state_defined_state(model):
    model.fs.state = model.fs.leach_soln.build_state_block(
        model.fs.time, defined_state=True
    )

    assert isinstance(model.fs.state[0].molar_concentration_constraint, Constraint)


@pytest.mark.unit
def test_state_vars(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    state_vars = model.fs.state[0].define_state_vars()

    assert set(state_vars.keys()) == {
        "flow_vol",
        "conc_mass_comp",
        "temperature",
        "pressure",
    }
    assert state_vars["flow_vol"] is model.fs.state[0].flow_vol
    assert state_vars["conc_mass_comp"] is model.fs.state[0].conc_mass_comp
    assert state_vars["temperature"] is model.fs.state[0].temperature
    assert state_vars["pressure"] is model.fs.state[0].pressure


@pytest.mark.unit
def test_pH_property(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    assert not model.fs.state[0].is_property_constructed("pH_phase")

    # Accessing pH_phase triggers construction of the variable and constraint
    assert isinstance(model.fs.state[0].pH_phase, Var)
    assert isinstance(model.fs.state[0].pH_constraint, Constraint)
    assert model.fs.state[0].is_property_constructed("pH_phase")


@pytest.mark.unit
def test_fix_state(model):
    model.fs.state = model.fs.leach_soln.build_state_block(model.fs.time)

    model.fs.state[0].flow_vol.set_value(10)
    model.fs.state[0].conc_mass_comp[:].set_value(1)
    model.fs.state[0].conc_mol_comp[:].set_value(1)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_vol.fixed
    assert model.fs.state[0].temperature.fixed
    assert model.fs.state[0].pressure.fixed
    for j in model.fs.leach_soln.component_list:
        if j == "H2O":
            assert not model.fs.state[0].conc_mass_comp[j].fixed
        else:
            assert model.fs.state[0].conc_mass_comp[j].fixed
        assert not model.fs.state[0].conc_mol_comp[j].fixed

        assert model.fs.state[0].molar_concentration_constraint[j].active

    assert model.fs.state[0].h2o_concentration.active
