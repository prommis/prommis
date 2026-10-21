#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# WaterTAP License Agreement

# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from prommis.ion_exchange.ix_parameter_estimation import (
    IXParmest,
    run_parmest_workflow,
    main,
)

import pytest

pytest.importorskip("watertap", reason="WaterTAP dependency not available")



@pytest.fixture(scope="module")
def parmest_obj():

    regenerant = "single_use"
    num_traps = 30
    c_trap_min = 1e-3

    data_name = "breakthrough_literature_data"
    curve_file = f"data/{data_name}.csv"

    n_val = 1.00001
    kt_val = 55
    curve = 1

    # Call run_parmest_workflow exactly as called in
    # `ix_parameter_estimation`.
    c = run_parmest_workflow(
        curve_file=curve_file,
        data_name=data_name,
        curve=curve,
        regenerant=regenerant,
        num_traps=num_traps,
        c_trap_min=c_trap_min,
        n_val=n_val,
        kt_val=kt_val,
    )
    return c


def test_full_parmest_run_and_data():

    regenerant = "single_use"
    num_traps = 30
    c_trap_min = 1e-3

    data_name = "breakthrough_literature_data"
    curve_file = f"data/{data_name}.csv"

    parmest_data = {
        "c0_ref": {},
        "loading_rate": {},
        "mass_transfer_coeff": {},
        "freundlich_n": {},
        "bv_50": {},
    }

    kt_val = 55
    n_val = 1.00001
    curve = [1]

    for cnum in curve:
        c = run_parmest_workflow(
            curve_file=curve_file,
            data_name=data_name,
            curve=cnum,
            regenerant=regenerant,
            num_traps=num_traps,
            c_trap_min=c_trap_min,
            n_val=n_val,
            kt_val=kt_val,
        )

        # Check that c is returned and is the expected type.
        assert c is not None
        assert hasattr(c, "parmest_theta")
        assert hasattr(c, "target_component")
        assert hasattr(c, "curve")

        # Extract variables and assert that they are not None.
        variables = vars(c)
        compound = None
        loading_rate = None
        c0 = None
        for key, val in variables.items():
            if key == "component_data":
                compound = val["name"]
            elif key == "loading_rate":
                loading_rate = val
            elif key == "c0":
                c0 = val
        assert compound is not None
        assert loading_rate is not None
        assert c0 is not None

        # Check parmest_theta structure.
        parmest_variables = vars(c.parmest_theta)
        assert "_mgr" in parmest_variables
        mgr = parmest_variables["_mgr"]
        items = mgr.items
        values = mgr.blocks[0].values
        assert len(items) == len(values)
        parmest_values = dict(zip(items, values))

        # Check that expected keys are in parmest_values.
        for key in [
            f"fs.ix.mass_transfer_coeff[{c.target_component}]",
            f"fs.ix.freundlich_n[{c.target_component}]",
            f"fs.ix.bv_50[{c.target_component}]",
        ]:
            assert key in parmest_values


def test_main(parmest_obj):

    # Test main() function.
    parmest_data = main(
        regenerant="single_use",
        num_traps=30,
        c_trap_min=1e-3,
        data_name="breakthrough_literature_data",
        kt_val=55,
        n_val=1.00001,
        curve_list=[1],
    )

    # Check that we return a dictionary with expected keys.
    for key in [
        "c0_ref",
        "loading_rate",
        "mass_transfer_coeff",
        "freundlich_n",
        "bv_50",
    ]:
        assert key in parmest_data


def test_ixparmest_methods(parmest_obj):
    # Test all major methods of IXParmest.

    m = parmest_obj.build_model()
    assert m is not None

    parmest_obj.add_bounds(m)
    parmest_obj.add_operating_conditions(m)
    parmest_obj.scale_it(m)
    parmest_obj.deactivate_trapezoids(m)

    # Test _get_ix_blk.
    ix_blk = parmest_obj._get_ix_blk(m)
    assert ix_blk is not None

    parmest_obj.add_initialization(m)
    parmest_obj.run_all_things()
    parmest_obj.get_bvs_cbs_breakthrough()

    # Test build_parmest.
    m_parmest = parmest_obj.build_parmest()
    assert m_parmest is not None

    # Test build_parmest_experiment_list.
    parmest_obj.build_parmest_experiment_list()
    assert hasattr(parmest_obj, "experiment_list")

    # Test run_parmest.
    solver_options = {
        "print_level": 5,
        "halt_on_ampl_error": "yes",
        "nlp_scaling_method": "gradient-based",
        "max_iter": 1000,
        "tol": 1e-6,
        "bound_push": 1e-8,
        "constr_viol_tol": 1e-6,
        "acceptable_constr_viol_tol": 1e-6,
    }
    parmest_obj.run_parmest(solver_options)

    # Test that parmest_obj is populated.
    assert hasattr(parmest_obj, "parmest_theta")
    assert hasattr(parmest_obj, "thetas")
    theta_dict = parmest_obj.parmest_theta.to_dict()
    assert isinstance(theta_dict, dict)


def test_missing_initial_guess_raises():

    kwargs = {
        "fix_initial_guess": {},
        "theta_initial_guess": {},
        "thetas": ["freundlich_n", "mass_transfer_coeff", "bv_50"],
        "set_bounds": {},
        "calc_from_constr": {},
        "component_data": {"charge": 1, "mw": 1, "molar_vol": 1, "diffusivity": 1},
        "resin_data": {
            "S950": {
                "type": "cation",
                "functional_group": "sulfonic",
                "density_kgm3": {"value": 1},
                "diameter": {"value": 1},
                "porosity": {"value": 0.8},
            }
        },
        "regenerant": "single_use",
        "data_name": "breakthrough_literature_data",
        "curve": 1,
        "num_traps": 30,
        "c_trap_min": 1e-3,
        "diff_calculation": "HaydukLaudie",
        "c0_min_thresh": 0,
    }
    with pytest.raises(ValueError, match="need to provide initial guess"):
        IXParmest(**kwargs)


def test_diff_calculation_branches(parmest_obj):
    # Test HaydukLaudie branch
    parmest_obj.diff_calculation = "HaydukLaudie"
    parmest_obj.get_curve_conditions(num_traps=1, c_trap_min=1e-3)
    assert "diffus_calculation" in parmest_obj.ion_props

    # Test non-HaydukLaudie branch
    parmest_obj.diff_calculation = "Other"
    parmest_obj.get_curve_conditions(num_traps=1, c_trap_min=1e-3)
    assert "diffusivity_data" in parmest_obj.ion_props
