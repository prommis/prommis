#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the multi-component diafiltration property model for streams.
"""

from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError

import pytest

from prommis.nanofiltration.multi_component_diafiltration_stream_properties import (
    MultiComponentDiafiltrationStreamParameter,
)


# Test single-salt model
@pytest.fixture
def model_single_salt():
    cation_list = ["Li"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.stream_properties = MultiComponentDiafiltrationStreamParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_single_salt(model_single_salt):
    assert len(model_single_salt.fs.stream_properties.phase_list) == 1
    for k in model_single_salt.fs.stream_properties.phase_list:
        assert k == "liquid"

    for j in model_single_salt.fs.stream_properties.component_list:
        assert j in [
            "Li",
            "Cl",
        ]


@pytest.mark.unit
def test_build(model_single_salt):
    assert len(model_single_salt.fs.stream_properties.config) == 3

    model_single_salt.fs.state = (
        model_single_salt.fs.stream_properties.build_state_block(
            model_single_salt.fs.time
        )
    )

    assert len(model_single_salt.fs.state) == 1

    assert isinstance(model_single_salt.fs.state[0].flow_vol, Var)
    assert isinstance(model_single_salt.fs.state[0].conc_mol_comp, Var)

    model_single_salt.fs.state[0].flow_vol.set_value(10)
    for j in model_single_salt.fs.stream_properties.component_list:
        model_single_salt.fs.state[0].conc_mol_comp[j].set_value(1)

    model_single_salt.fs.state.fix_initialization_states()

    assert model_single_salt.fs.state[0].flow_vol.fixed
    for j in model_single_salt.fs.stream_properties.component_list:
        assert model_single_salt.fs.state[0].conc_mol_comp[j].fixed


# Test two-salt model
@pytest.fixture
def model_two_salt():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.stream_properties = MultiComponentDiafiltrationStreamParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_two_salt(model_two_salt):
    assert len(model_two_salt.fs.stream_properties.phase_list) == 1
    for k in model_two_salt.fs.stream_properties.phase_list:
        assert k == "liquid"

    for j in model_two_salt.fs.stream_properties.component_list:
        assert j in [
            "Li",
            "Co",
            "Cl",
        ]


@pytest.mark.unit
def test_build_two_salt(model_two_salt):
    test_build(model_two_salt)


# Test three-salt model
@pytest.fixture
def model_three_salt():
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.stream_properties = MultiComponentDiafiltrationStreamParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_three_salt(model_three_salt):
    assert len(model_three_salt.fs.stream_properties.phase_list) == 1
    for k in model_three_salt.fs.stream_properties.phase_list:
        assert k == "liquid"

    for j in model_three_salt.fs.stream_properties.component_list:
        assert j in [
            "Li",
            "Co",
            "Al",
            "Cl",
        ]


@pytest.mark.unit
def test_build_three_salt(model_three_salt):
    test_build(model_three_salt)


# Test common anion exception
@pytest.mark.component
def test_common_anion_exception():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl", "sulfate"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    with pytest.raises(
        ConfigurationError,
        match="The multi-component diafiltration unit model only supports systems with a common anion",
    ):
        m.fs.stream_properties = MultiComponentDiafiltrationStreamParameter(
            cation_list=cation_list,
            anion_list=anion_list,
        )
