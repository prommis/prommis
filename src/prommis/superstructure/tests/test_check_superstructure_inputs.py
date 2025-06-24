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

from prommis.superstructure.check_superstructure_inputs import (  # check_byproduct_valorization_params,; check_discretized_costing_params,; check_environmental_impact_params,; check_feed_params,; check_operating_params,; check_supe_formulation_params,
    check_feed_params,
    check_plant_lifetime_params,
    check_supe_formulation_params,
)

solver_available = SolverFactory("gurobi").available()
if solver_available:
    solver = get_solver(solver="gurobi")
else:
    solver = None


@pytest.fixture(scope="module")
def get_params():
    return {
        ### Plant Lifetime Parameters
        "plant_start": 2024,
        "plant_lifetime": 15,
        ### Feed parameters
        "available_feed": {
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
        },
        "collection_rate": 0.1,
        "tracked_comps": ["Nd", "Dy", "Fe"],
        "prod_comp_mass": {
            "Nd": 0.206 * 3,
            "Dy": 0.103 * 3,
            "Fe": 0.691 * 3,
        },
        ### Superstructure formulation parameters
        "num_stages": 5,
        "options_in_stage": {
            1: 2,
            2: 4,
            3: 6,
            4: 4,
            5: 5,
        },
        "option_outlets": {
            # level 1
            (1, 1): [1, 2, 3, 4],
            (1, 2): [1, 2, 3, 4],
            # level 2
            (2, 1): [1, 2, 3, 6],
            (2, 2): [1, 2, 3, 6],
            (2, 3): [1, 2, 3, 6],
            (2, 4): [4, 5],
            # level 3
            (3, 1): [1],
            (3, 2): [1],
            (3, 3): [2, 3],
            (3, 4): [2, 3],
            (3, 5): [2, 3],
            (3, 6): [4],
            # level 4
            (4, 1): [1],
            (4, 2): [2],
            (4, 3): [3],
            (4, 4): [4, 5],
        },
        "option_efficiencies": {
            # Level 1 yields
            (1, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
            (1, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
            # level 2 yields
            (2, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
            (2, 2): {"Nd": 1, "Dy": 1, "Fe": 1},
            (2, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
            (2, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
            # level 3 yields
            (3, 1): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
            (3, 2): {"Nd": 0.985, "Dy": 0.985, "Fe": 0},
            (3, 3): {"Nd": 0.925, "Dy": 0.98, "Fe": 0},
            (3, 4): {"Nd": 1, "Dy": 1, "Fe": 0},
            (3, 5): {"Nd": 1, "Dy": 1, "Fe": 0},
            (3, 6): {"Nd": 1, "Dy": 1, "Fe": 0.403},
            # level 4 yields
            (4, 1): {"Nd": 1, "Dy": 1, "Fe": 1},
            (4, 2): {"Nd": 1, "Dy": 0.899, "Fe": 0},
            (4, 3): {"Nd": 1, "Dy": 1, "Fe": 1},
            (4, 4): {"Nd": 1, "Dy": 1, "Fe": 1},
            # level 5 yields
            (5, 1): {
                "Nd": 1,
                "Dy": 1,
                "Fe": 0,
            },
            (5, 2): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
            (5, 3): {"Nd": 0.98, "Dy": 0.98, "Fe": 0},
            (5, 4): {
                "Nd": 0.98,
                "Dy": 0.98,
                "Fe": 0,
            },
            (5, 5): {
                "Nd": 0.98,
                "Dy": 0.98,
                "Fe": 0,
            },
        },
        ### Operating parameters
        "profit": {
            (5, 1): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 2): {"Nd": 69.888, "Dy": 263.81, "Fe": 0},
            (5, 3): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 4): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
            (5, 5): {"Nd": 45.4272, "Dy": 171.4765, "Fe": 0},
        },
        "opt_var_oc_params": {
            # level 2
            (2, 1): {"a": 0.0053, "b": 7929.7},
            (2, 2): {"a": 0.0015, "b": 2233.16},
            (2, 3): {"a": 0.0034, "b": 0},
            (2, 4): {"a": 0.0117, "b": 0},
            # level 3
            (3, 1): {"a": 15.594, "b": 4e6},
            (3, 2): {"a": 35.58463, "b": 4e6},
            (3, 3): {"a": 1.8359, "b": 0},
            (3, 4): {"a": 3.7414, "b": 2378.6},
            (3, 5): {"a": 10.35427, "b": 2378.6},
            (3, 6): {"a": 1.58, "b": 0},
            # level 4
            (4, 1): {"a": 0, "b": 0},
            (4, 2): {"a": 111.09, "b": 254606},
            (4, 3): {"a": 0, "b": 0},
            (4, 4): {"a": 0, "b": 0},
            # level 5
            (5, 1): {"a": 0.4997, "b": 89832},
            (5, 2): {"a": 9.8127, "b": 964921},
            (5, 3): {"a": 9.8127, "b": 964921},
            (5, 4): {"a": 2.17, "b": 0},
            (5, 5): {"a": 6.7063559004, "b": 0},
        },
        "operators_per_discrete_unit": {
            (1, 1): 1,
            (1, 2): 0,
        },
        "yearly_cost_per_unit": {
            (1, 1): 1,
            (1, 2): 0,
        },
        "capital_cost_per_unit": {
            (1, 1): 1,
            (1, 2): 0,
        },
        "processing_rate": {
            (1, 1): 7868,
            (1, 2): 52453,
        },
        "num_operators": {
            (2, 1): 0.65,
            (2, 2): 0.65,
            (2, 3): 0.65,
            (2, 4): 0.65,
            (3, 1): 1.6,
            (3, 2): 1.6,
            (3, 3): 1.3,
            (3, 4): 0.45,
            (3, 5): 0.45,
            (3, 6): 1.15,
            (4, 1): 0,
            (4, 2): 1.3,
            (4, 3): 0,
            (4, 4): 0,
            (5, 1): 1.05,
            (5, 2): 0.75,
            (5, 3): 0.75,
            (5, 4): 1.15,
            (5, 5): 1.15,
        },
        "labor_rate": 8000 * 38.20,
        ### Discretized costing Parameters
        "discretized_purchased_equipment_cost": {
            (2, 1): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    10130.08515,
                    31353.21173,
                    48788.84678,
                    60305.81927,
                    77063.4884,
                    117214.7546,
                    151018.0699,
                    195698.5419,
                ],
            },
            (2, 2): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    11702.08515,
                    39023.21173,
                    62134.84678,
                    77539.81927,
                    100113.4884,
                    154792.7546,
                    201326.0699,
                    263374.5419,
                ],
            },
            (2, 3): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    67799.0,
                    272739.0,
                    406588.0,
                    488109.0,
                    599667.0,
                    842665.0,
                    1028482.0,
                    1255541.0,
                ],
            },
            (2, 4): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    121352.0,
                    490387.0,
                    732185.0,
                    879652.0,
                    1081651.0,
                    1522296.0,
                    1859740.0,
                    2272553.0,
                ],
            },
            (3, 1): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    343228.652,
                    482425.4684,
                    618182.0594,
                    743750.2902,
                    844443.0443,
                    978479.5225,
                    1183834.522,
                    1440660.587,
                ],
            },
            (3, 2): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    643228.652,
                    782425.4684,
                    918182.0594,
                    1043750.2902,
                    1144443.0443,
                    1278479.5225,
                    1483834.522,
                    1740660.587,
                ],
            },
            (3, 3): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    423074.7216,
                    3042779.121,
                    5348359.01,
                    6921261.68,
                    9251002.61,
                    14933803.33,
                    19762044.75,
                    26151302.79,
                ],
            },
            (3, 4): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    226790.0,
                    446435.0,
                    713714.0,
                    1270105.0,
                    1541353.0,
                    2920751.0,
                    3652064.0,
                    5323087.0,
                ],
            },
            (3, 5): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    476790.0,
                    696435.0,
                    963714.0,
                    1520105.0,
                    1791353.0,
                    3170751.0,
                    3902064.0,
                    5573087.0,
                ],
            },
            (3, 6): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    169491.0,
                    940300.0,
                    1534578.0,
                    1919653.0,
                    2469724.0,
                    3743401.0,
                    4774426.0,
                    6089420.0,
                ],
            },
            (4, 1): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            },
            (4, 2): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    450073.0,
                    1623424.0,
                    2349169.0,
                    2778687.0,
                    3349410.0,
                    4585349.0,
                    5503177.0,
                    6590434.0,
                ],
            },
            (4, 3): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            },
            (4, 4): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            },
            (5, 1): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    197813.1853,
                    324151.478,
                    468747.756,
                    547368.2774,
                    655614.4213,
                    800184.1752,
                    974415.7068,
                    1114534.98,
                ],
            },
            (5, 2): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    222906.0,
                    354009.0,
                    490597.0,
                    562047.0,
                    679397.0,
                    912244.0,
                    1097498.0,
                    1297052.0,
                ],
            },
            (5, 3): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    222906.0,
                    354009.0,
                    490597.0,
                    562047.0,
                    679397.0,
                    912244.0,
                    1097498.0,
                    1297052.0,
                ],
            },
            (5, 4): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    154685.0,
                    858157.0,
                    1400520.0,
                    1751956.0,
                    2253973.0,
                    3416384.0,
                    4357340.0,
                    5557458.0,
                ],
            },
            (5, 5): {
                "Flowrates": [
                    0.0,
                    36480.0,
                    634240.0,
                    1434800.0,
                    2083760.0,
                    3171200.0,
                    6342400.0,
                    9513600.0,
                    14270400.0,
                ],
                "Costs": [
                    0.0,
                    404685.0,
                    1108157.0,
                    1650520.0,
                    2001956.0,
                    2503973.0,
                    3666384.0,
                    4607340.0,
                    5807458.0,
                ],
            },
        },
        ### Environmental impact parameters
        "options_environmental_impacts": {
            (1, 1): 0,
            (1, 2): 1000,
            (2, 1): 0,
            (2, 2): 1000,
            (2, 3): 600,
            (2, 4): 800,
            (3, 1): 600,
            (3, 2): 0,
            (3, 3): 600,
            (3, 4): 800,
            (3, 5): 800,
            (3, 6): 1000,
            (4, 1): 0,
            (4, 2): 800,
            (4, 3): 600,
            (4, 4): 1000,
            (5, 1): 0,
            (5, 2): 800,
            (5, 3): 600,
            (5, 4): 800,
            (5, 5): 1000,
        },
        "epsilon": 1,
        ### Byproduct valorization parameters
        "byproduct_values": {
            "Jarosite": -0.17,
            "Iron oxide": 10,
            "Residue": -0.17,
            "Iron hydroxide": -0.17,
        },
        "byproduct_opt_conversions": {
            (3, 1): {"Jarosite": 0.75},
            (3, 2): {"Iron oxide": 1},
            (3, 3): {"Residue": 0.25},
            (3, 4): {"Iron hydroxide": 0.5},
            (3, 5): {"Iron oxide": 1},
            (3, 6): {"Iron oxide": 1},
            (5, 4): {"Iron hydroxide": 0.5},
            (5, 5): {"Iron oxide": 1},
        },
    }

### Test plant lifetime parameters
def test_plant_start_type():
    """Test that TypeError is raised if plant_start is not an int."""
    with pytest.raises(TypeError, match="plant_start is not of type int."):
        check_plant_lifetime_params(plant_start="2020", plant_lifetime=10)


def test_plant_lifetime_type():
    """Test that TypeError is raised if plant_lifetime is not an int."""
    with pytest.raises(TypeError, match="plant_lifetime is not of type int."):
        check_plant_lifetime_params(plant_start=2020, plant_lifetime="10")


def test_plant_lifetime_min():
    """Test that ValueError is raised if plant_lifetime is less than 3."""
    with pytest.raises(
        ValueError, match="Plant lifetime must be a minimum of three years."
    ):
        check_plant_lifetime_params(plant_start=2020, plant_lifetime=2)


def test_valid_inputs():
    """Test that no error is raised for valid inputs."""
    check_plant_lifetime_params(plant_start=2020, plant_lifetime=5)


### Test feed parameters
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
