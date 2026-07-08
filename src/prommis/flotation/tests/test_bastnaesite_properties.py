#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the bastnaesite dry-solids property package."""

from pyomo.environ import (
    ConcreteModel,
    Expression,
    Param,
    Var,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core.scaling import get_scaling_factor
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.flotation.bastnaesite_properties import (
    BASTNAESITE_BOND_WORK_INDEX,
    BastnaesiteParameters,
    BastnaesitePropertiesScaler,
    COMPONENTS,
    REFERENCE_SCALING_FACTORS,
)


@pytest.fixture()
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters()
    m.fs.state = m.fs.properties.build_state_block(m.fs.time, defined_state=True)
    return m


@pytest.fixture()
def particle_size_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BastnaesiteParameters(has_particle_size_distribution=True)
    m.fs.state = m.fs.properties.build_state_block(m.fs.time, defined_state=True)
    return m


@pytest.mark.unit
@pytest.mark.build
def test_build(model):
    assert list(model.fs.properties.component_list) == list(COMPONENTS)
    assert list(model.fs.properties.phase_list) == ["solid"]
    assert isinstance(model.fs.properties.rho_mass_comp, Param)
    assert not hasattr(model.fs.properties, "bond_work_index")
    assert isinstance(model.fs.state[0].flow_mass_comp, Var)
    assert isinstance(model.fs.state[0].flow_mass, Expression)
    assert isinstance(model.fs.state[0].flow_vol_phase, Expression)
    assert not model.fs.state[0].is_property_constructed("mass_frac_comp")
    assert isinstance(model.fs.state[0].mass_frac_comp, Expression)
    assert model.fs.state[0].is_property_constructed("mass_frac_comp")
    assert not hasattr(model.fs.state[0], "mass_frac_comp_eq")
    assert not hasattr(model.fs.state[0], "particle_size_median")
    assert not hasattr(model.fs.state[0], "particle_size_width")
    assert str(model.fs.state[0].flow_mass_comp["REO"].get_units()) == str(
        units.kg / units.hour
    )


@pytest.mark.unit
def test_material_flow_basis_and_terms(model):
    state = model.fs.state[0]
    assert state.get_material_flow_basis() == MaterialFlowBasis.mass
    assert state.get_material_flow_terms("solid", "REO") is state.flow_mass_comp["REO"]
    assert set(state.define_state_vars()) == {"flow_mass_comp"}


@pytest.mark.unit
def test_derived_mass_quantities(model):
    state = model.fs.state[0]
    component_flows = {
        "REO": 10.0,
        "CaO": 2.0,
        "BaO": 1.0,
        "SrO": 1.0,
        "inert_gangue": 6.0,
    }
    total_flow = sum(component_flows.values())

    for component, flow in component_flows.items():
        state.flow_mass_comp[component].set_value(flow)

    assert value(state.flow_mass) == pytest.approx(total_flow)
    assert value(state.flow_vol_phase["solid"]) == pytest.approx(
        sum(
            flow / value(model.fs.properties.rho_mass_comp[component])
            for component, flow in component_flows.items()
        )
    )
    for component, flow in component_flows.items():
        assert value(state.mass_frac_comp[component]) == pytest.approx(
            flow / total_flow
        )


@pytest.mark.unit
def test_rho_mass_comp_rejects_non_positive(model):
    with pytest.raises(ValueError):
        model.fs.properties.rho_mass_comp["REO"].set_value(0)
    with pytest.raises(ValueError):
        model.fs.properties.rho_mass_comp["REO"].set_value(-1)


@pytest.mark.unit
def test_state_fixing_dof_and_units(model):
    for component in COMPONENTS:
        model.fs.state[0].flow_mass_comp[component].fix(1)
    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model)


@pytest.mark.unit
def test_default_scaler(model):
    assert model.fs.state[0].default_scaler is BastnaesitePropertiesScaler
    scaler = model.fs.state[0].default_scaler()
    scaler.variable_scaling_routine(model.fs.state[0])
    for component in COMPONENTS:
        assert get_scaling_factor(
            model.fs.state[0].flow_mass_comp[component]
        ) == pytest.approx(REFERENCE_SCALING_FACTORS[component])


@pytest.mark.unit
@pytest.mark.build
def test_particle_size_distribution_config(particle_size_model):
    state = particle_size_model.fs.state[0]

    assert isinstance(particle_size_model.fs.properties.bond_work_index, Param)
    assert value(particle_size_model.fs.properties.bond_work_index) == pytest.approx(
        BASTNAESITE_BOND_WORK_INDEX
    )
    assert isinstance(state.particle_size_median, Var)
    assert isinstance(state.particle_size_width, Var)
    assert str(state.particle_size_median.get_units()) == str(units.um)
    assert set(state.define_state_vars()) == {
        "flow_mass_comp",
        "particle_size_median",
        "particle_size_width",
    }


@pytest.mark.unit
def test_particle_size_distribution_dof_and_units(particle_size_model):
    for component in COMPONENTS:
        particle_size_model.fs.state[0].flow_mass_comp[component].fix(1)
    particle_size_model.fs.state[0].particle_size_median.fix(80000)
    particle_size_model.fs.state[0].particle_size_width.fix(1.5)

    assert degrees_of_freedom(particle_size_model) == 0
    assert_units_consistent(particle_size_model)
