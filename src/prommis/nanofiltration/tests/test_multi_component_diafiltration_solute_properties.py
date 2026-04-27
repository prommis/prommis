#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
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
    cation_list = ["Li"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build(sample_model_single_salt):
    prop_package = sample_model_single_salt.fs.solute_properties

    assert len(prop_package.config) == 3

    sample_model_single_salt.fs.state = prop_package.build_state_block(
        sample_model_single_salt.fs.time
    )

    assert len(sample_model_single_salt.fs.state) == 1

    assert isinstance(sample_model_single_salt.fs.state[0].flow_vol, Var)
    assert isinstance(sample_model_single_salt.fs.state[0].conc_mol_comp, Var)

    sample_model_single_salt.fs.state[0].flow_vol.set_value(10)
    for j in prop_package.component_list:
        sample_model_single_salt.fs.state[0].conc_mol_comp[j].set_value(1)

    sample_model_single_salt.fs.state.fix_initialization_states()

    assert sample_model_single_salt.fs.state[0].flow_vol.fixed
    for j in prop_package.component_list:
        assert sample_model_single_salt.fs.state[0].conc_mol_comp[j].fixed


@pytest.mark.unit
def test_parameters_single_salt(sample_model_single_salt):
    prop_package = sample_model_single_salt.fs.solute_properties

    assert len(prop_package.phase_list) == 1
    for k in prop_package.phase_list:
        assert k == "liquid"

    for j in prop_package.component_list:
        assert j in prop_package.charge
        assert j in prop_package.boundary_layer_diffusion_coefficient
        assert j in prop_package.membrane_diffusion_coefficient
        assert j in prop_package.sigma
        assert j in prop_package.partition_coefficient_retentate
        assert j in prop_package.partition_coefficient_permeate
        assert j in prop_package.num_solutes


################################################################################
# Test single-salt model: LiCl
@pytest.fixture
def model_Li():
    cation_list = ["Li"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Li(model_Li):
    test_build(model_Li)


@pytest.mark.unit
def test_parameters_single_salt_Li(model_Li):
    test_parameters_single_salt(model_Li)

    prop_package = model_Li.fs.solute_properties

    for j in prop_package.component_list:
        assert j in ["Li", "Cl"]

    # check Li values
    assert value(prop_package.charge["Li"]) == 1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.membrane_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.sigma["Li"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Li"]) == 0.4
    assert value(prop_package.partition_coefficient_permeate["Li"]) == 0.4
    assert value(prop_package.num_solutes["Li"]) == 1

    # check Cl values (single salt with Li)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 1


################################################################################
# Test single-salt model: CoCl2
@pytest.fixture
def model_Co():
    cation_list = ["Co"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Co(model_Co):
    test_build(model_Co)


@pytest.mark.unit
def test_parameters_single_salt_Co(model_Co):
    test_parameters_single_salt(model_Co)

    prop_package = model_Co.fs.solute_properties

    for j in prop_package.component_list:
        assert j in ["Co", "Cl"]

    # check Co values
    assert value(prop_package.charge["Co"]) == 2
    assert value(prop_package.boundary_layer_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.membrane_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.sigma["Co"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Co"]) == 0.04
    assert value(prop_package.partition_coefficient_permeate["Co"]) == 0.04
    assert value(prop_package.num_solutes["Co"]) == 1

    # check Cl values (single salt with Co)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 2


################################################################################
# Test single-salt model: AlCl3
@pytest.fixture
def model_Al():
    cation_list = ["Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Al(model_Al):
    test_build(model_Al)


@pytest.mark.unit
def test_parameters_single_salt_Al(model_Al):
    test_parameters_single_salt(model_Al)

    prop_package = model_Al.fs.solute_properties

    for j in prop_package.component_list:
        assert j in ["Al", "Cl"]

    # check Al values
    assert value(prop_package.charge["Al"]) == 3
    assert value(prop_package.boundary_layer_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.membrane_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.sigma["Al"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Al"]) == 0.004
    assert value(prop_package.partition_coefficient_permeate["Al"]) == 0.004
    assert value(prop_package.num_solutes["Al"]) == 1

    # check Cl values (single salt with Al)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 3


################################################################################
# Test functions for two-salt model
@pytest.fixture
def sample_model_two_salt():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_two_salt(sample_model_two_salt):
    prop_package = sample_model_two_salt.fs.solute_properties

    assert len(prop_package.phase_list) == 1
    for k in prop_package.phase_list:
        assert k == "liquid"

    for j in prop_package.component_list:
        assert j in prop_package.charge
        assert j in prop_package.boundary_layer_diffusion_coefficient
        assert j in prop_package.membrane_diffusion_coefficient
        assert j in prop_package.sigma
        assert j in prop_package.partition_coefficient_retentate
        assert j in prop_package.partition_coefficient_permeate
        assert j in prop_package.num_solutes


################################################################################
# Test two-salt model: LiCl + CoCl2
@pytest.fixture
def model_Li_Co():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Li_Co(model_Li_Co):
    test_build(model_Li_Co)


@pytest.mark.unit
def test_parameters_two_salt_Li_Co(model_Li_Co):
    test_parameters_two_salt(model_Li_Co)

    prop_package = model_Li_Co.fs.solute_properties
    for j in prop_package.component_list:
        assert j in ["Li", "Co", "Cl"]

    # check Li values
    assert value(prop_package.charge["Li"]) == 1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.membrane_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.sigma["Li"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Li"]) == 0.4
    assert value(prop_package.partition_coefficient_permeate["Li"]) == 0.4
    assert value(prop_package.num_solutes["Li"]) == 1

    # check Co values
    assert value(prop_package.charge["Co"]) == 2
    assert value(prop_package.boundary_layer_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.membrane_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.sigma["Co"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Co"]) == 0.04
    assert value(prop_package.partition_coefficient_permeate["Co"]) == 0.04
    assert value(prop_package.num_solutes["Co"]) == 1

    # check Cl values (two salt with Li and Co)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 3


################################################################################


# Test two-salt model: LiCl + AlCl3
@pytest.fixture
def model_Li_Al():
    cation_list = ["Li", "Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Li_Al(model_Li_Al):
    test_build(model_Li_Al)


@pytest.mark.unit
def test_parameters_two_salt_Li_Al(model_Li_Al):
    test_parameters_two_salt(model_Li_Al)

    prop_package = model_Li_Al.fs.solute_properties

    for j in prop_package.component_list:
        assert j in ["Li", "Al", "Cl"]

    # check Li values
    assert value(prop_package.charge["Li"]) == 1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.membrane_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.sigma["Li"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Li"]) == 0.4
    assert value(prop_package.partition_coefficient_permeate["Li"]) == 0.4
    assert value(prop_package.num_solutes["Li"]) == 1

    # check Al values
    assert value(prop_package.charge["Al"]) == 3
    assert value(prop_package.boundary_layer_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.membrane_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.sigma["Al"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Al"]) == 0.004
    assert value(prop_package.partition_coefficient_permeate["Al"]) == 0.004
    assert value(prop_package.num_solutes["Al"]) == 1

    # check Cl values (two salt with Li and Al)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 4


################################################################################


# Test two-salt model: CoCl2 + AlCl3
@pytest.fixture
def model_Co_Al():
    cation_list = ["Co", "Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Co_Al(model_Co_Al):
    test_build(model_Co_Al)


@pytest.mark.unit
def test_parameters_two_salt_Co_Al(model_Co_Al):
    test_parameters_two_salt(model_Co_Al)

    prop_package = model_Co_Al.fs.solute_properties

    for j in prop_package.component_list:
        assert j in ["Co", "Al", "Cl"]

    # check Co values
    assert value(prop_package.charge["Co"]) == 2
    assert value(prop_package.boundary_layer_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.membrane_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.sigma["Co"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Co"]) == 0.04
    assert value(prop_package.partition_coefficient_permeate["Co"]) == 0.04
    assert value(prop_package.num_solutes["Co"]) == 1

    # check Al values
    assert value(prop_package.charge["Al"]) == 3
    assert value(prop_package.boundary_layer_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.membrane_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.sigma["Al"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Al"]) == 0.004
    assert value(prop_package.partition_coefficient_permeate["Al"]) == 0.004
    assert value(prop_package.num_solutes["Al"]) == 1

    # check Cl values (two salt with Co and Al)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 5


################################################################################
# Test functions for three-salt model
@pytest.fixture
def sample_model_three_salt():
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_parameters_three_salt(sample_model_three_salt):
    prop_package = sample_model_three_salt.fs.solute_properties

    assert len(prop_package.phase_list) == 1
    for k in prop_package.phase_list:
        assert k == "liquid"

    for j in prop_package.component_list:
        assert j in prop_package.charge
        assert j in prop_package.boundary_layer_diffusion_coefficient
        assert j in prop_package.membrane_diffusion_coefficient
        assert j in prop_package.sigma
        assert j in prop_package.partition_coefficient_retentate
        assert j in prop_package.partition_coefficient_permeate
        assert j in prop_package.num_solutes


################################################################################
# Test three-salt model: LiCl + CoCl2 + AlCl3
@pytest.fixture
def model_Li_Co_Al():
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.solute_properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    return m


@pytest.mark.unit
def test_build_Li_Co_Al(model_Li_Co_Al):
    test_build(model_Li_Co_Al)
    for j in model_Li_Co_Al.fs.solute_properties.component_list:
        assert j in ["Li", "Co", "Al", "Cl"]


@pytest.mark.unit
def test_parameters_three_salt_Li_Co_Al(model_Li_Co_Al):
    test_parameters_three_salt(model_Li_Co_Al)

    prop_package = model_Li_Co_Al.fs.solute_properties

    # check Li values
    assert value(prop_package.charge["Li"]) == 1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.membrane_diffusion_coefficient["Li"]) == 3.71
    assert value(prop_package.sigma["Li"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Li"]) == 0.4
    assert value(prop_package.partition_coefficient_permeate["Li"]) == 0.4
    assert value(prop_package.num_solutes["Li"]) == 1

    # check Co values
    assert value(prop_package.charge["Co"]) == 2
    assert value(prop_package.boundary_layer_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.membrane_diffusion_coefficient["Co"]) == 2.64
    assert value(prop_package.sigma["Co"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Co"]) == 0.04
    assert value(prop_package.partition_coefficient_permeate["Co"]) == 0.04
    assert value(prop_package.num_solutes["Co"]) == 1

    # check Al values
    assert value(prop_package.charge["Al"]) == 3
    assert value(prop_package.boundary_layer_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.membrane_diffusion_coefficient["Al"]) == 2.01
    assert value(prop_package.sigma["Al"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Al"]) == 0.004
    assert value(prop_package.partition_coefficient_permeate["Al"]) == 0.004
    assert value(prop_package.num_solutes["Al"]) == 1

    # check Cl values (two salt with Li, Co, and Al)
    assert value(prop_package.charge["Cl"]) == -1
    assert value(prop_package.boundary_layer_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.membrane_diffusion_coefficient["Cl"]) == 7.31
    assert value(prop_package.sigma["Cl"]) == 1
    assert value(prop_package.partition_coefficient_retentate["Cl"]) == 0.01
    assert value(prop_package.partition_coefficient_permeate["Cl"]) == 0.01
    assert value(prop_package.num_solutes["Cl"]) == 6


################################################################################
# Test common anion exception


@pytest.mark.component
def test_common_anion_exception():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl", "SO4"]

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
