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
from idaes.core.util.exceptions import ConfigurationError

import pytest

from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)


################################################################################
# Test functions for single-salt model
@pytest.fixture
def sample_model_single_salt():
    cation_list = ["lithium"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
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
    cation_list = ["lithium"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_lithium(model_lithium):
    test_build(model_lithium)


@pytest.mark.unit
def test_parameters_single_salt_lithium(model_lithium):
    test_parameters_single_salt(model_lithium)
    for j in model_lithium.fs.solute_properties.component_list:
        assert j in [
            "lithium",
            "chloride",
        ]

    # check lithium values
    assert value(model_lithium.fs.solute_properties.charge["lithium"]) == 1
    assert (
        value(model_lithium.fs.solute_properties.diffusion_coefficient["lithium"])
        == 3.71
    )
    assert value(model_lithium.fs.solute_properties.sigma["lithium"]) == 1
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_retentate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_permeate["lithium"]
        )
        == 0.4
    )
    assert value(model_lithium.fs.solute_properties.num_solutes["lithium"]) == 1

    # check chloride values (single salt with lithium)
    assert value(model_lithium.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(model_lithium.fs.solute_properties.diffusion_coefficient["chloride"])
        == 7.31
    )
    assert value(model_lithium.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert value(model_lithium.fs.solute_properties.num_solutes["chloride"]) == 1


################################################################################
# Test single-salt model: cobalt chloride
@pytest.fixture
def model_cobalt():
    cation_list = ["cobalt"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_cobalt(model_cobalt):
    test_build(model_cobalt)


@pytest.mark.unit
def test_parameters_single_salt_cobalt(model_cobalt):
    test_parameters_single_salt(model_cobalt)
    for j in model_cobalt.fs.solute_properties.component_list:
        assert j in [
            "cobalt",
            "chloride",
        ]

    # check cobalt values
    assert value(model_cobalt.fs.solute_properties.charge["cobalt"]) == 2
    assert (
        value(model_cobalt.fs.solute_properties.diffusion_coefficient["cobalt"]) == 2.64
    )
    assert value(model_cobalt.fs.solute_properties.sigma["cobalt"]) == 1
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_retentate["cobalt"]
        )
        == 0.04
    )
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_permeate["cobalt"]
        )
        == 0.04
    )
    assert value(model_cobalt.fs.solute_properties.num_solutes["cobalt"]) == 1

    # check chloride values (single salt with cobalt)
    assert value(model_cobalt.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(model_cobalt.fs.solute_properties.diffusion_coefficient["chloride"])
        == 7.31
    )
    assert value(model_cobalt.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_cobalt.fs.solute_properties.partition_coefficient_permeate["chloride"]
        )
        == 0.01
    )
    assert value(model_cobalt.fs.solute_properties.num_solutes["chloride"]) == 2


################################################################################
# Test single-salt model: aluminum chloride
@pytest.fixture
def model_aluminum():
    cation_list = ["aluminum"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_aluminum(model_aluminum):
    test_build(model_aluminum)


@pytest.mark.unit
def test_parameters_single_salt_aluminum(model_aluminum):
    test_parameters_single_salt(model_aluminum)
    for j in model_aluminum.fs.solute_properties.component_list:
        assert j in [
            "aluminum",
            "chloride",
        ]

    # check aluminum values
    assert value(model_aluminum.fs.solute_properties.charge["aluminum"]) == 3
    assert (
        value(model_aluminum.fs.solute_properties.diffusion_coefficient["aluminum"])
        == 2.01
    )
    assert value(model_aluminum.fs.solute_properties.sigma["aluminum"]) == 1
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert value(model_aluminum.fs.solute_properties.num_solutes["aluminum"]) == 1

    # check chloride values (single salt with aluminum)
    assert value(model_aluminum.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(model_aluminum.fs.solute_properties.diffusion_coefficient["chloride"])
        == 7.31
    )
    assert value(model_aluminum.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert value(model_aluminum.fs.solute_properties.num_solutes["chloride"]) == 3


################################################################################
# Test functions for two-salt model
@pytest.fixture
def sample_model_two_salt():
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_two_salt(sample_model_two_salt):
    assert len(sample_model_two_salt.fs.solute_properties.phase_list) == 1
    for k in sample_model_two_salt.fs.solute_properties.phase_list:
        assert k == "liquid"

    for j in sample_model_two_salt.fs.solute_properties.component_list:
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
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_lithium_cobalt(model_lithium_cobalt):
    test_build(model_lithium_cobalt)


@pytest.mark.unit
def test_parameters_two_salt_lithium_cobalt(model_lithium_cobalt):
    test_parameters_two_salt(model_lithium_cobalt)
    for j in model_lithium_cobalt.fs.solute_properties.component_list:
        assert j in [
            "lithium",
            "cobalt",
            "chloride",
        ]

    # check lithium values
    assert value(model_lithium_cobalt.fs.solute_properties.charge["lithium"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["lithium"]
        )
        == 3.71
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["lithium"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["lithium"]) == 1

    # check cobalt values
    assert value(model_lithium_cobalt.fs.solute_properties.charge["cobalt"]) == 2
    assert (
        value(model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["cobalt"])
        == 2.64
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["cobalt"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["cobalt"]) == 1

    # check chloride values (two salt with lithium and cobalt)
    assert value(model_lithium_cobalt.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.diffusion_coefficient["chloride"]
        )
        == 7.31
    )
    assert value(model_lithium_cobalt.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_cobalt.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert value(model_lithium_cobalt.fs.solute_properties.num_solutes["chloride"]) == 3


################################################################################


# Test two-salt model: lithium chloride + aluminum chloride
@pytest.fixture
def model_lithium_aluminum():
    cation_list = ["lithium", "aluminum"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_lithium_aluminum(model_lithium_aluminum):
    test_build(model_lithium_aluminum)


@pytest.mark.unit
def test_parameters_two_salt_lithium_aluminum(model_lithium_aluminum):
    test_parameters_two_salt(model_lithium_aluminum)
    for j in model_lithium_aluminum.fs.solute_properties.component_list:
        assert j in [
            "lithium",
            "aluminum",
            "chloride",
        ]

    # check lithium values
    assert value(model_lithium_aluminum.fs.solute_properties.charge["lithium"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient["lithium"]
        )
        == 3.71
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["lithium"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(model_lithium_aluminum.fs.solute_properties.num_solutes["lithium"]) == 1
    )

    # check aluminum values
    assert value(model_lithium_aluminum.fs.solute_properties.charge["aluminum"]) == 3
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient[
                "aluminum"
            ]
        )
        == 2.01
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["aluminum"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(model_lithium_aluminum.fs.solute_properties.num_solutes["aluminum"]) == 1
    )

    # check chloride values (two salt with lithium and aluminum)
    assert value(model_lithium_aluminum.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.diffusion_coefficient[
                "chloride"
            ]
        )
        == 7.31
    )
    assert value(model_lithium_aluminum.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(model_lithium_aluminum.fs.solute_properties.num_solutes["chloride"]) == 4
    )


################################################################################


# Test two-salt model: cobalt chloride + aluminum chloride
@pytest.fixture
def model_cobalt_aluminum():
    cation_list = ["cobalt", "aluminum"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_cobalt_aluminum(model_cobalt_aluminum):
    test_build(model_cobalt_aluminum)


@pytest.mark.unit
def test_parameters_two_salt_cobalt_aluminum(model_cobalt_aluminum):
    test_parameters_two_salt(model_cobalt_aluminum)
    for j in model_cobalt_aluminum.fs.solute_properties.component_list:
        assert j in [
            "cobalt",
            "aluminum",
            "chloride",
        ]

    # check cobalt values
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["cobalt"]) == 2
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["cobalt"]
        )
        == 2.64
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["cobalt"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.num_solutes["cobalt"]) == 1

    # check aluminum values
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["aluminum"]) == 3
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["aluminum"]
        )
        == 2.01
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["aluminum"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(model_cobalt_aluminum.fs.solute_properties.num_solutes["aluminum"]) == 1
    )

    # check chloride values (two salt with cobalt and aluminum)
    assert value(model_cobalt_aluminum.fs.solute_properties.charge["chloride"]) == -1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.diffusion_coefficient["chloride"]
        )
        == 7.31
    )
    assert value(model_cobalt_aluminum.fs.solute_properties.sigma["chloride"]) == 1
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(model_cobalt_aluminum.fs.solute_properties.num_solutes["chloride"]) == 5
    )


################################################################################
# Test functions for three-salt model
@pytest.fixture
def sample_model_three_salt():
    cation_list = ["lithium", "cobalt", "aluminum"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_three_salt(sample_model_three_salt):
    assert len(sample_model_three_salt.fs.solute_properties.phase_list) == 1
    for k in sample_model_three_salt.fs.solute_properties.phase_list:
        assert k == "liquid"

    for j in sample_model_three_salt.fs.solute_properties.component_list:
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
    cation_list = ["lithium", "cobalt", "aluminum"]
    anion_list = ["chloride"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_lithium_cobalt_aluminum(model_lithium_cobalt_aluminum):
    test_build(model_lithium_cobalt_aluminum)
    for j in model_lithium_cobalt_aluminum.fs.solute_properties.component_list:
        assert j in [
            "lithium",
            "cobalt",
            "aluminum",
            "chloride",
        ]


@pytest.mark.unit
def test_parameters_three_salt_lithium_cobalt_aluminum(model_lithium_cobalt_aluminum):
    test_parameters_three_salt(model_lithium_cobalt_aluminum)

    # check lithium values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["lithium"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "lithium"
            ]
        )
        == 3.71
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["lithium"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "lithium"
            ]
        )
        == 0.4
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["lithium"])
        == 1
    )

    # check cobalt values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["cobalt"]) == 2
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "cobalt"
            ]
        )
        == 2.64
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["cobalt"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "cobalt"
            ]
        )
        == 0.04
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["cobalt"])
        == 1
    )

    # check aluminum values
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["aluminum"])
        == 3
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "aluminum"
            ]
        )
        == 2.01
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["aluminum"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "aluminum"
            ]
        )
        == 0.004
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["aluminum"]
        )
        == 1
    )

    # check chloride values (two salt with lithium, cobalt, and aluminum)
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.charge["chloride"])
        == -1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.diffusion_coefficient[
                "chloride"
            ]
        )
        == 7.31
    )
    assert (
        value(model_lithium_cobalt_aluminum.fs.solute_properties.sigma["chloride"]) == 1
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_retentate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.partition_coefficient_permeate[
                "chloride"
            ]
        )
        == 0.01
    )
    assert (
        value(
            model_lithium_cobalt_aluminum.fs.solute_properties.num_solutes["chloride"]
        )
        == 6
    )


################################################################################
# Test common anion exception


@pytest.mark.component
def test_common_anion_exception():
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride", "sulfate"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    with pytest.raises(
        ConfigurationError,
        match="The multi-component diafiltration unit model only supports systems with a common anion",
    ):
        m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
            cation_list=cation_list,
            anion_list=anion_list,
        )


################################################################################
