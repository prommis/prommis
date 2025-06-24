#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import pytest

from prommis.superstructure.check_superstructure_inputs import (
    check_supe_formulation_params,
)


### Test superstructure formulation parameters
def create_mock_model():
    class MockFeedParams:
        class tracked_comps:
            @staticmethod
            def data():
                return {"Nd", "Dy", "Fe"}

    class MockModel:
        feed_params = MockFeedParams()

    return MockModel()


def create_mock_model_incomplete_tracked_comps():
    class MockFeedParams:
        class tracked_comps:
            @staticmethod
            def data():
                return {"Nd", "Dy"}  # Missing "Fe" for some tests

    class MockModel:
        feed_params = MockFeedParams()

    return MockModel()


VALID_NUM_STAGES = 5
VALID_OPTIONS_IN_STAGE = {1: 2, 2: 4, 3: 6, 4: 4, 5: 5}
VALID_OPTION_OUTLETS = {
    (1, 1): [1, 2, 3, 4],
    (1, 2): [1, 2, 3, 4],
    (2, 1): [1, 2, 3, 6],
    (2, 2): [1, 2, 3, 6],
    (2, 3): [1, 2, 3, 6],
    (2, 4): [4, 5],
    (3, 1): [1],
    (3, 2): [1],
    (3, 3): [2, 3],
    (3, 4): [2, 3],
    (3, 5): [2, 3],
    (3, 6): [4],
    (4, 1): [1],
    (4, 2): [2],
    (4, 3): [3],
    (4, 4): [4, 5],
}
VALID_OPTION_EFFICIENCIES = {
    (1, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (1, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
    (2, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    (3, 1): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
    (3, 2): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
    (3, 3): {"Nd": 0.925, "Dy": 0.98, "Fe": 0},
    (3, 4): {"Nd": 1, "Dy": 1, "Fe": 0},
    (3, 5): {"Nd": 1, "Dy": 1, "Fe": 0},
    (3, 6): {"Nd": 1, "Dy": 1, "Fe": 0.403},
    (4, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 2): {"Nd": 1, "Dy": 0.899, "Fe": 0},
    (4, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
    (4, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
    (5, 1): {"Nd": 1, "Dy": 1, "Fe": 0},
    (5, 2): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 3): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 4): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
    (5, 5): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
}


def test_valid_inputs():
    """Test no errors with valid inputs."""
    m = create_mock_model()
    check_supe_formulation_params(
        m,
        VALID_NUM_STAGES,
        VALID_OPTIONS_IN_STAGE,
        VALID_OPTION_OUTLETS,
        VALID_OPTION_EFFICIENCIES,
    )


def test_num_stages_not_int():
    """Test num_stages not of type int."""
    m = create_mock_model()
    with pytest.raises(TypeError, match="num_stages is not of type int"):
        check_supe_formulation_params(
            m,
            "5",
            VALID_OPTIONS_IN_STAGE,
            VALID_OPTION_OUTLETS,
            VALID_OPTION_EFFICIENCIES,
        )


def test_options_in_stage_not_dict():
    """Test options_in_stage not of type dict."""
    m = create_mock_model()
    with pytest.raises(TypeError, match="options_in_stage is not of type dict"):
        check_supe_formulation_params(
            m, VALID_NUM_STAGES, [], VALID_OPTION_OUTLETS, VALID_OPTION_EFFICIENCIES
        )


def test_option_outlets_not_dict():
    """Test option_outlets not of type dict."""
    m = create_mock_model()
    with pytest.raises(TypeError, match="option_outlets is not of type dict"):
        check_supe_formulation_params(
            m, VALID_NUM_STAGES, VALID_OPTIONS_IN_STAGE, [], VALID_OPTION_EFFICIENCIES
        )


def test_option_efficiencies_not_dict():
    """Test option_efficiencies not of type dict."""
    m = create_mock_model()
    with pytest.raises(TypeError, match="option_efficiencies is not of type dict"):
        check_supe_formulation_params(
            m, VALID_NUM_STAGES, VALID_OPTIONS_IN_STAGE, VALID_OPTION_OUTLETS, []
        )


def test_num_stages_less_than_2():
    """Test num_stages less than 2."""
    m = create_mock_model()
    with pytest.raises(ValueError, match="There must be at least 2 processing stages"):
        check_supe_formulation_params(m, 1, {1: 1}, {(1, 1): [1]}, {(1, 1): {"Nd": 1}})


def test_options_in_stage_keys_mismatch():
    """Test options_in_stage keys don't match num_stages."""
    m = create_mock_model()
    with pytest.raises(ValueError, match="Stages must start at 1 and count up"):
        check_supe_formulation_params(m, VALID_NUM_STAGES, {2: 1, 3: 2}, {}, {})


def test_options_in_stage_value_not_positive():
    """Test options_in_stage value not positive."""
    m = create_mock_model()
    with pytest.raises(ValueError, match="number of options in stage"):
        check_supe_formulation_params(
            m, 2, {1: 0, 2: 1}, {(1, 1): [1]}, {(1, 1): {"Nd": 1}}
        )


def test_option_outlets_empty_list():
    """Test option_outlets value is empty list."""
    m = create_mock_model()
    with pytest.raises(ValueError, match="is an empty list"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 1}, {(1, 1): []}, {(1, 1): {"Nd": 1}}
        )


def test_option_outlets_not_positive_int():
    """Test option_outlets contains non-positive int."""
    m = create_mock_model()
    with pytest.raises(TypeError, match="All elements must be positive integers"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 1}, {(1, 1): [0]}, {(1, 1): {"Nd": 1}}
        )


def test_option_outlets_undefined():
    """Test option_outlets not defined for all options."""
    m = create_mock_model()
    bad_outlets = {(1, 1): [1]}
    with pytest.raises(ValueError, match="don't have defined outlets"):
        check_supe_formulation_params(
            m, 2, {1: 2, 2: 1}, bad_outlets, {(1, 1): {"Nd": 1}, (1, 2): {"Nd": 1}}
        )


def test_option_outlets_undefined_opts():
    """Test option_outlets contains undefined options."""
    m = create_mock_model()
    bad_outlets = {(1, 1): [1], (1, 3): [1]}
    with pytest.raises(ValueError, match="don't exist in the superstructure"):
        check_supe_formulation_params(
            m, 2, {1: 2, 2: 1}, bad_outlets, {(1, 1): {"Nd": 1}}
        )


def test_option_outlets_infeasible_outlet():
    """Test option_outlets contains infeasible outlet."""
    m = create_mock_model()
    bad_outlets = {(1, 1): [2]}
    with pytest.raises(ValueError, match="not in the superstructure"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 1}, bad_outlets, {(1, 1): {"Nd": 1}}
        )


def test_disconnected_options():
    """Test options not connected from previous stage."""
    m = create_mock_model()
    bad_outlets = {(1, 1): [2]}
    with pytest.raises(ValueError, match="not connected from the previous stage"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 2}, bad_outlets, {(1, 1): {"Nd": 1}}
        )


def test_option_efficiencies_not_specified():
    """Test option_efficiencies not specified for all options."""
    m = create_mock_model()
    bad_efficiencies = {(1, 1): {"Nd": 1}}
    with pytest.raises(ValueError, match="Efficiences not specified"):
        check_supe_formulation_params(
            m, 2, {1: 2, 2: 1}, {(1, 1): [1], (1, 2): [1]}, bad_efficiencies
        )


def test_option_efficiencies_infeasible_opts():
    """Test option_efficiencies contains infeasible options."""
    m = create_mock_model()
    bad_efficiencies = {(1, 1): {"Nd": 1}, (1, 3): {"Nd": 1}}
    with pytest.raises(ValueError, match="do not exist in the superstructure"):
        check_supe_formulation_params(
            m, 2, {1: 2, 2: 1}, {(1, 1): [1], (1, 2): [1]}, bad_efficiencies
        )


def test_option_efficiencies_efficiency_not_in_0_1():
    """Test option_efficiencies efficiency not in [0, 1]."""
    m = create_mock_model()
    with pytest.raises(ValueError, match="not in [0, 1]"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 1}, {(1, 1): [1]}, {(1, 1): {"Nd": 1.1}}
        )


def test_option_efficiencies_missing_effs():
    """Test option_efficiencies missing efficiencies for tracked components."""
    m = create_mock_model()
    bad_efficiencies = {(1, 1): {"Nd": 1}, (1, 2): {"Nd": 1}}
    with pytest.raises(ValueError, match="missing components"):
        check_supe_formulation_params(
            m, 2, {1: 2, 2: 1}, {(1, 1): [1], (1, 2): [1]}, bad_efficiencies
        )


def test_option_efficiencies_undefined_component():
    """Test option_efficiencies contains undefined component."""
    m = create_mock_model_incomplete_tracked_comps()  # Only tracks Nd, Dy
    bad_efficiencies = {(1, 1): {"Nd": 1, "Dy": 1, "Fe": 1}}
    with pytest.raises(ValueError, match="is not defined as a tracked component"):
        check_supe_formulation_params(
            m, 2, {1: 1, 2: 1}, {(1, 1): [1]}, bad_efficiencies
        )
