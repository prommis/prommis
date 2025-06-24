#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import warnings

from pyomo.environ import SolverFactory

from idaes.core.solvers import get_solver

import pytest

from prommis.superstructure.check_superstructure_inputs import check_feed_params

# Helper function to create a mock model with operational years
def create_mock_model(operational_years):
    class MockPlantParams:
        class operational_range:
            @staticmethod
            def data():
                return operational_years

    class MockModel:
        plant_lifetime_params = MockPlantParams()

    return MockModel()


# Test data for valid inputs
VALID_FEED = {
    2025: 290273,
    2026: 274648,
    2027: 286512,
    2028: 487819,
    2029: 592637,
    2030: 571054,
    2031: 498472,
    2032: 506565,
    2033: 566355,
    2034: 669094,
    2035: 719057,
    2036: 762656,
    2037: 1434637,
    2038: 1697805,
}
VALID_TRACKED = ["Nd", "Dy", "Fe"]
VALID_PROD_MASS = {"Nd": 0.618, "Dy": 0.309, "Fe": 2.073}


# Positive test - no errors/warnings
def test_valid_inputs():
    """Test no errors/warnings with valid inputs"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # Treat warnings as errors
        check_feed_params(m, VALID_FEED, 0.1, VALID_TRACKED, VALID_PROD_MASS)


# Type error tests
def test_empty_available_feed():
    """Test empty available_feed dict"""
    m = create_mock_model(set())
    with pytest.raises(TypeError, match="available_feed dict is empty"):
        check_feed_params(m, {}, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_non_int_key_in_feed():
    """Test non-int key in available_feed"""
    m = create_mock_model({2025})
    with pytest.raises(TypeError, match="key 2025.0 in available_feed"):
        check_feed_params(m, {2025.0: 100}, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_non_numeric_value_in_feed():
    """Test non-numeric value in available_feed"""
    m = create_mock_model({2025})
    with pytest.raises(TypeError, match="value abc in available_feed"):
        check_feed_params(m, {2025: "abc"}, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_non_numeric_collection_rate():
    """Test non-numeric collection_rate"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="collection_rate is not of type"):
        check_feed_params(m, VALID_FEED, "0.1", VALID_TRACKED, VALID_PROD_MASS)


def test_non_list_tracked_comps():
    """Test tracked_comps not being a list"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="tracked_comps is not of type list"):
        check_feed_params(m, VALID_FEED, 0.1, {"Nd", "Dy"}, VALID_PROD_MASS)


def test_non_string_in_tracked_comps():
    """Test non-string in tracked_comps"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="Value 123 in tracked_comps"):
        check_feed_params(m, VALID_FEED, 0.1, ["Nd", 123], VALID_PROD_MASS)


def test_non_dict_prod_comp_mass():
    """Test prod_comp_mass not being a dict"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="prod_comp_mass is not of type dict"):
        check_feed_params(m, VALID_FEED, 0.1, VALID_TRACKED, [("Nd", 0.618)])


def test_non_string_key_in_prod_mass():
    """Test non-string key in prod_comp_mass"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="key 123 in prod_comp_mass"):
        check_feed_params(m, VALID_FEED, 0.1, VALID_TRACKED, {123: 0.618})


def test_non_numeric_value_in_prod_mass():
    """Test non-numeric value in prod_comp_mass"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(TypeError, match="value abc in prod_comp_mass"):
        check_feed_params(m, VALID_FEED, 0.1, VALID_TRACKED, {"Nd": "abc"})


# Value error tests
def test_feed_years_mismatch():
    """Test mismatch between feed years and operational years"""
    m = create_mock_model({2025, 2026})  # Different from VALID_FEED keys
    with pytest.raises(ValueError, match="Years of available_feed do not match"):
        check_feed_params(m, VALID_FEED, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_negative_feed_value():
    """Test negative value in available_feed"""
    m = create_mock_model({2025})
    with pytest.raises(ValueError, match="available_feed contains negative values"):
        check_feed_params(m, {2025: -100}, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_all_zero_feed():
    """Test all-zero available_feed"""
    m = create_mock_model({2025})
    with pytest.raises(ValueError, match="All values in available_feed are zero"):
        check_feed_params(m, {2025: 0, 2026: 0}, 0.1, VALID_TRACKED, VALID_PROD_MASS)


def test_non_positive_collection_rate():
    """Test non-positive collection_rate"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(ValueError, match="Collection rate must be a positive value"):
        check_feed_params(m, VALID_FEED, 0, VALID_TRACKED, VALID_PROD_MASS)


def test_empty_tracked_comps():
    """Test empty tracked_comps list"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(ValueError, match="tracked_comps list is empty"):
        check_feed_params(m, VALID_FEED, 0.1, [], VALID_PROD_MASS)


def test_mismatched_keys():
    """Test mismatch between tracked_comps and prod_comp_mass keys"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(ValueError, match="prod_comp_mass keys don't match"):
        check_feed_params(m, VALID_FEED, 0.1, ["Nd", "Dy"], {"Nd": 0.618})


def test_negative_prod_mass():
    """Test negative value in prod_comp_mass"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.raises(ValueError, match="prod_comp_mass contains negative values"):
        check_feed_params(
            m, VALID_FEED, 0.1, VALID_TRACKED, {"Nd": -0.618, "Dy": 0.309, "Fe": 2.073}
        )


# Warning test
def test_zero_prod_mass_warning():
    """Test warning for zero value in prod_comp_mass"""
    m = create_mock_model(set(VALID_FEED.keys()))
    with pytest.warns(UserWarning, match="prod_comp_mass contains zero values"):
        check_feed_params(
            m, VALID_FEED, 0.1, VALID_TRACKED, {"Nd": 0, "Dy": 0.309, "Fe": 2.073}
        )