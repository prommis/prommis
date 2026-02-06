#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the multi-component diafiltration property model for solutes.
"""

from pyomo.environ import ConcreteModel, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)


################################################################################
# Test functions for single-salt model
@pytest.fixture
def sample_model_single_salt():
    num_salts = 1
    salt_system = "lithium_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build(sample_model_single_salt):
    assert len(sample_model_single_salt.fs.solute_properties.config) == 3

    sample_model_single_salt.fs.state = (
        sample_model_single_salt.fs.solute_properties.build_state_block(
            sample_model_single_salt.fs.time
        )
    )

    assert len(sample_model_single_salt.fs.state) == 1

    assert isinstance(sample_model_single_salt.fs.state[0].flow_vol, Var)
    assert isinstance(sample_model_single_salt.fs.state[0].conc_mol_comp, Var)

    sample_model_single_salt.fs.state[0].flow_vol.set_value(10)
    for j in sample_model_single_salt.fs.solute_properties.component_list:
        sample_model_single_salt.fs.state[0].conc_mol_comp[j].set_value(1)

    sample_model_single_salt.fs.state.fix_initialization_states()

    assert sample_model_single_salt.fs.state[0].flow_vol.fixed
    for j in sample_model_single_salt.fs.solute_properties.component_list:
        assert sample_model_single_salt.fs.state[0].conc_mol_comp[j].fixed


@pytest.mark.unit
def test_parameters_single_salt(sample_model_single_salt):
    assert len(sample_model_single_salt.fs.solute_properties.phase_list) == 1
    for k in sample_model_single_salt.fs.solute_properties.phase_list:
        assert k == "liquid"

    for j in sample_model_single_salt.fs.solute_properties.component_list:
        assert j in [
            "cation_1",
            "anion",
        ]
        assert j in sample_model_single_salt.fs.solute_properties.charge
        assert j in sample_model_single_salt.fs.solute_properties.diffusion_coefficient
        assert j in sample_model_single_salt.fs.solute_properties.sigma
        assert (
            j
            in sample_model_single_salt.fs.solute_properties.partition_coefficient_retentate
        )
        assert (
            j
            in sample_model_single_salt.fs.solute_properties.partition_coefficient_permeate
        )
        assert j in sample_model_single_salt.fs.solute_properties.num_solutes


################################################################################
# Test single-salt model: lithium chloride
@pytest.fixture
def model_lithium():
    num_salts = 1
    salt_system = "lithium_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_lithium(model_lithium):
    test_build(model_lithium)


@pytest.mark.unit
def test_parameters_single_salt_lithium(model_lithium):
    test_parameters_single_salt(model_lithium)

    # check lithium values
    assert value(model_lithium.fs.solute_properties.charge["cation_1"]) == 1
    assert (
        value(model_lithium.fs.solute_properties.diffusion_coefficient["cation_1"])
        == 3.71
    )
    assert value(model_lithium.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert value(model_lithium.fs.solute_properties.num_solutes["cation_1"]) == 1

    # check chloride values (single salt with lithium)
    assert value(model_lithium.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(model_lithium.fs.solute_properties.diffusion_coefficient["anion"]) == 7.31
    )
    assert value(model_lithium.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_retentate["anion"]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_permeate["anion"]
        )
        == 0.01
    )
    assert value(model_lithium.fs.solute_properties.num_solutes["anion"]) == 1


################################################################################
# Test single-salt model: cobalt chloride
@pytest.fixture
def model_cobalt():
    num_salts = 1
    salt_system = "cobalt_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_cobalt(model_cobalt):
    test_build(model_cobalt)


@pytest.mark.unit
def test_parameters_single_salt_cobalt(model_cobalt):
    test_parameters_single_salt(model_cobalt)

    # check cobalt values
    assert value(model_cobalt.fs.solute_properties.charge["cation_1"]) == 2
    assert (
        value(model_cobalt.fs.solute_properties.diffusion_coefficient["cation_1"])
        == 2.64
    )
    assert value(model_cobalt.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_permeate["cation_1"]
        )
        == 0.04
    )
    assert value(model_cobalt.fs.solute_properties.num_solutes["cation_1"]) == 1

    # check chloride values (single salt with cobalt)
    assert value(model_cobalt.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(model_cobalt.fs.solute_properties.diffusion_coefficient["anion"]) == 7.31
    )
    assert value(model_cobalt.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_retentate["anion"]
        )
        == 0.01
    )
    assert (
        value(model_cobalt.fs.solute_properties.partition_coefficient_permeate["anion"])
        == 0.01
    )
    assert value(model_cobalt.fs.solute_properties.num_solutes["anion"]) == 2


################################################################################
# Test single-salt model: aluminum chloride
@pytest.fixture
def model_aluminum():
    num_salts = 1
    salt_system = "aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_aluminum(model_aluminum):
    test_build(model_aluminum)


@pytest.mark.unit
def test_parameters_single_salt_aluminum(model_aluminum):
    test_parameters_single_salt(model_aluminum)

    # check aluminum values
    assert value(model_aluminum.fs.solute_properties.charge["cation_1"]) == 3
    assert (
        value(model_aluminum.fs.solute_properties.diffusion_coefficient["cation_1"])
        == 2.01
    )
    assert value(model_aluminum.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.004
    )
    assert value(model_aluminum.fs.solute_properties.num_solutes["cation_1"]) == 1

    # check chloride values (single salt with aluminum)
    assert value(model_aluminum.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(model_aluminum.fs.solute_properties.diffusion_coefficient["anion"])
        == 7.31
    )
    assert value(model_aluminum.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_retentate["anion"]
        )
        == 0.01
    )
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_permeate["anion"]
        )
        == 0.01
    )
    assert value(model_aluminum.fs.solute_properties.num_solutes["anion"]) == 3


################################################################################
# Test functions for two-salt model
@pytest.fixture
def sample_model_two_salt():
    num_salts = 2
    salt_system = "lithium_cobalt_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_parameters_two_salt(sample_model_two_salt):
    assert len(sample_model_two_salt.fs.solute_properties.phase_list) == 1
    for k in sample_model_two_salt.fs.solute_properties.phase_list:
        assert k == "liquid"

    for j in sample_model_two_salt.fs.solute_properties.component_list:
        assert j in [
            "cation_1",
            "cation_2",
            "anion",
        ]
        assert j in sample_model_two_salt.fs.solute_properties.charge
        assert j in sample_model_two_salt.fs.solute_properties.diffusion_coefficient
        assert j in sample_model_two_salt.fs.solute_properties.sigma
        assert (
            j
            in sample_model_two_salt.fs.solute_properties.partition_coefficient_retentate
        )
        assert (
            j
            in sample_model_two_salt.fs.solute_properties.partition_coefficient_permeate
        )
        assert j in sample_model_two_salt.fs.solute_properties.num_solutes


################################################################################
# Test two-salt model: lithium chloride + cobalt chloride
@pytest.fixture
def model_lithium_cobalt():
    num_salts = 2
    salt_system = "lithium_cobalt_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_lithium_cobalt(model_lithium_cobalt):
    test_build(model_lithium_cobalt)


@pytest.mark.unit
def test_parameters_two_salt_lithium_cobalt(model_lithium_cobalt):
    test_parameters_two_salt(model_lithium_cobalt)

    # check lithium values
    assert value(model_lithium_cobalt.fs.solute_properties.charge["cation_1"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["cation_1"]
        )
        == 3.71
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["cation_1"]) == 1

    # check cobalt values
    assert value(model_lithium_cobalt.fs.solute_properties.charge["cation_2"]) == 2
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["cation_2"]
        )
        == 2.64
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["cation_2"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "cation_2"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "cation_2"
            ]
        )
        == 0.04
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["cation_2"]) == 1

    # check chloride values (two salt with lithium and cobalt)
    assert value(model_lithium_cobalt.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["anion"])
        == 7.31
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "anion"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "anion"
            ]
        )
        == 0.01
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["anion"]) == 3


################################################################################


# Test two-salt model: lithium chloride + aluminum chloride
@pytest.fixture
def model_lithium_aluminum():
    num_salts = 2
    salt_system = "lithium_aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_lithium_aluminum(model_lithium_aluminum):
    test_build(model_lithium_aluminum)


@pytest.mark.unit
def test_parameters_two_salt_lithium_aluminum(model_lithium_aluminum):
    test_parameters_two_salt(model_lithium_aluminum)

    # check lithium values
    assert value(model_lithium_aluminum.fs.solute_properties.charge["cation_1"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient[
                "cation_1"
            ]
        )
        == 3.71
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(model_lithium_aluminum.fs.solute_properties.num_solutes["cation_1"]) == 1
    )

    # check aluminum values
    assert value(model_lithium_aluminum.fs.solute_properties.charge["cation_2"]) == 3
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient[
                "cation_2"
            ]
        )
        == 2.01
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["cation_2"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_2"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_2"
            ]
        )
        == 0.004
    )
    assert (
        value(model_lithium_aluminum.fs.solute_properties.num_solutes["cation_2"]) == 1
    )

    # check chloride values (two salt with lithium and aluminum)
    assert value(model_lithium_aluminum.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient["anion"]
        )
        == 7.31
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "anion"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "anion"
            ]
        )
        == 0.01
    )
    assert value(model_lithium_aluminum.fs.solute_properties.num_solutes["anion"]) == 4


################################################################################


# Test two-salt model: cobalt chloride + aluminum chloride
@pytest.fixture
def model_cobalt_aluminum():
    num_salts = 2
    salt_system = "cobalt_aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_cobalt_aluminum(model_cobalt_aluminum):
    test_build(model_cobalt_aluminum)


@pytest.mark.unit
def test_parameters_two_salt_cobalt_aluminum(model_cobalt_aluminum):
    test_parameters_two_salt(model_cobalt_aluminum)

    # check cobalt values
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["cation_1"]) == 2
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["cation_1"]
        )
        == 2.64
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["cation_1"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.04
    )
    assert (
        value(model_cobalt_aluminum.fs.solute_properties.num_solutes["cation_1"]) == 1
    )

    # check aluminum values
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["cation_2"]) == 3
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["cation_2"]
        )
        == 2.01
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["cation_2"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_2"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_2"
            ]
        )
        == 0.004
    )
    assert (
        value(model_cobalt_aluminum.fs.solute_properties.num_solutes["cation_2"]) == 1
    )

    # check chloride values (two salt with cobalt and aluminum)
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["anion"]) == -1
    assert (
        value(model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["anion"])
        == 7.31
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "anion"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "anion"
            ]
        )
        == 0.01
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.num_solutes["anion"]) == 5


################################################################################
# Test functions for three-salt model
@pytest.fixture
def sample_model_three_salt():
    num_salts = 3
    salt_system = "lithium_cobalt_aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_parameters_three_salt(sample_model_three_salt):
    assert len(sample_model_three_salt.fs.solute_properties.phase_list) == 1
    for k in sample_model_three_salt.fs.solute_properties.phase_list:
        assert k == "liquid"

    for j in sample_model_three_salt.fs.solute_properties.component_list:
        assert j in [
            "cation_1",
            "cation_2",
            "cation_3",
            "anion",
        ]
        assert j in sample_model_three_salt.fs.solute_properties.charge
        assert j in sample_model_three_salt.fs.solute_properties.diffusion_coefficient
        assert j in sample_model_three_salt.fs.solute_properties.sigma
        assert (
            j
            in sample_model_three_salt.fs.solute_properties.partition_coefficient_retentate
        )
        assert (
            j
            in sample_model_three_salt.fs.solute_properties.partition_coefficient_permeate
        )
        assert j in sample_model_three_salt.fs.solute_properties.num_solutes


################################################################################
# Test three-salt model: lithium chloride + cobalt chloride + aluminum chloride
@pytest.fixture
def model_lithium_cobalt_aluminum():
    num_salts = 3
    salt_system = "lithium_cobalt_aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    return m


@pytest.mark.unit
def test_build_lithium_cobalt_aluminum(model_lithium_cobalt_aluminum):
    test_build(model_lithium_cobalt_aluminum)


@pytest.mark.unit
def test_parameters_three_salt_lithium_cobalt_aluminum(model_lithium_cobalt_aluminum):
    test_parameters_three_salt(model_lithium_cobalt_aluminum)

    # check lithium values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["cation_1"])
        == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "cation_1"
            ]
        )
        == 3.71
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["cation_1"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_1"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["cation_1"]
        )
        == 1
    )

    # check cobalt values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["cation_2"])
        == 2
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "cation_2"
            ]
        )
        == 2.64
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["cation_2"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_2"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_2"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["cation_2"]
        )
        == 1
    )

    # check aluminum values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["cation_3"])
        == 3
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "cation_3"
            ]
        )
        == 2.01
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["cation_3"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cation_3"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cation_3"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["cation_3"]
        )
        == 1
    )

    # check chloride values (two salt with lithium, cobalt, and aluminum)
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["anion"]) == -1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "anion"
            ]
        )
        == 7.31
    )
    assert value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["anion"]) == 1
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "anion"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "anion"
            ]
        )
        == 0.01
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["anion"])
        == 6
    )


################################################################################
