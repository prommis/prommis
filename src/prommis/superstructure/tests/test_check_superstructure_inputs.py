#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    SolverFactory,
)

from idaes.core.solvers import get_solver

import pytest

from prommis.superstructure.check_superstructure_inputs import (  # check_byproduct_valorization_params,; check_discretized_costing_params,; check_environmental_impact_params,; check_feed_params,; check_operating_params,; check_supe_formulation_params,
    check_plant_lifetime_params,
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


# class TestSuperstructureParams(object):
#     @pytest.fixture(scope="class")
# def get_superstructure_params(self, get_params):
# common_params = get_params
# plant_start=common_params["plant_start"]
# plant_lifetime=common_params["plant_lifetime"]
# available_feed=common_params["available_feed"]
# collection_rate=common_params["collection_rate"]
# tracked_comps=common_params["tracked_comps"]
# prod_comp_mass=common_params["prod_comp_mass"]
# num_stages=common_params["num_stages"]
# options_in_stage=common_params["options_in_stage"]
# option_outlets=common_params["option_outlets"]
# option_efficiencies=common_params["option_efficiencies"]
# profit=common_params["profit"]
# opt_var_oc_params=common_params["opt_var_oc_params"]
# operators_per_discrete_unit=common_params["operators_per_discrete_unit"]
# yearly_cost_per_unit=common_params["yearly_cost_per_unit"]
# capital_cost_per_unit=common_params["capital_cost_per_unit"]
# processing_rate=common_params["processing_rate"]
# num_operators=common_params["num_operators"]
# labor_rate=common_params["labor_rate"]
# discretized_purchased_equipment_cost=common_params["discretized_purchased_equipment_cost"]
# options_environmental_impacts=common_params["options_environmental_impacts"]
# epsilon=common_params["epsilon"]
# byproduct_values=byproduct_values["byproduct_values"]


@pytest.mark.unit
def test_plant_lifetime_params(get_params):
    # params = get_params
    # plant_start=params["plant_start"]
    # plant_lifetime=params["plant_lifetime"]

    with pytest.raises(TypeError, match="plant_start is not of type int."):
        check_plant_lifetime_params(plant_start="2020", plant_lifetime=15)
