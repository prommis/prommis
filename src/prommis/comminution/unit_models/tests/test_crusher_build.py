#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Build, export and configuration-guard tests for CrusherSolidPSD."""

import math

import pytest

from pyomo.environ import units
from pyomo.util.check_units import assert_units_equivalent

from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom

from prommis.comminution.unit_models.tests.crusher_test_support import (
    CLASSIFICATION_RECYCLE_REF,
    SINGLE_STRESS_REF,
    crusher,
    fix_feed,
    one_mineral_flowsheet,
)

COARSE_ONE_MINERAL = {"OreA": [0.0, 0.0, 0.0, 10.0]}


@pytest.mark.unit
def test_distribution_function_d80_out_of_mesh_rejected_at_build():
    for bad_d80 in (20e-3, 1.0e-3):  # above x_N, below the bin-0 floor 1.8e-3
        m = one_mineral_flowsheet()
        with pytest.raises(ConfigurationError, match="attainable mesh band"):
            crusher(
                m,
                psd_method="distribution_function",
                shape_exponent=1.5,
                target_product_p80=bad_d80,
            )
    m = one_mineral_flowsheet()
    with pytest.raises(ConfigurationError, match="d_max must be <= the top mesh edge"):
        crusher(
            m,
            psd_method="distribution_function",
            distribution_shape="gates_gaudin_schuhmann",
            shape_exponent=1.5,
            target_product_p80=15e-3,  # d_max = 17.4e-3 > 16e-3
        )
    # A valid target gives each characteristic size a top-edge validity bound:
    # Rosin-Rammler uses 0.95 passing; GGS uses the top mesh edge.
    m = one_mineral_flowsheet()
    cr = crusher(
        m,
        psd_method="distribution_function",
        shape_exponent=1.5,
        target_product_p80=5e-3,
    )
    expected_ub = 16e-3 / ((-math.log(0.05)) ** (1.0 / 1.5))
    assert cr.d_63[0].ub == pytest.approx(expected_ub, rel=1e-12, abs=0)
    m2 = one_mineral_flowsheet()
    cr2 = crusher(
        m2,
        psd_method="distribution_function",
        distribution_shape="gates_gaudin_schuhmann",
        shape_exponent=1.5,
        target_product_p80=5e-3,
    )
    assert cr2.d_max[0].ub == pytest.approx(16e-3, rel=1e-12, abs=0)


OUTLET_P80_SCOPES = {
    "tabular",
    "single_stress",
    "repeated_stress",
    "classification_recycle",
}
# Each case lists its configuration and degrees of freedom with the feed fixed.
# The GGS power path has one free degree of freedom until power is fixed.
PATH_SHAPE_SCOPES = {
    "tabular": (dict(psd_method="tabular", tabular_psd=[0.25, 0.25, 0.25, 0.25]), 0),
    "rr_size_spec": (
        dict(
            psd_method="distribution_function",
            distribution_shape="rosin_rammler",
            shape_exponent=1.5,
            target_product_p80=5e-3,
        ),
        0,
    ),
    "ggs_power_spec": (
        dict(
            psd_method="distribution_function",
            distribution_shape="gates_gaudin_schuhmann",
            shape_exponent=1.5,
        ),
        1,
    ),
    "single_stress": (dict(SINGLE_STRESS_REF), 0),
    "repeated_stress": (dict(SINGLE_STRESS_REF, n_stress_events=3), 0),
    "classification_recycle": (dict(CLASSIFICATION_RECYCLE_REF), 0),
}


@pytest.mark.build
@pytest.mark.component
@pytest.mark.parametrize("scope", sorted(PATH_SHAPE_SCOPES))
def test_every_path_shape_scope_is_square_and_consistent(scope):
    """Check degrees of freedom, units, and structure for each product path."""
    cfg, expected_dof = PATH_SHAPE_SCOPES[scope]
    m = one_mineral_flowsheet()
    cr = crusher(m, bond_work_index=12.0, **cfg)
    fix_feed(cr.properties_in[0], COARSE_ONE_MINERAL)
    assert cr.power[0].lb == 0
    assert_units_equivalent(cr.power[0], units.kW)
    assert degrees_of_freedom(m) == expected_dof
    if expected_dof:
        cr.power[0].fix(3.0)
    # Bound checks cannot prove that outlet P80 is positive for tabular, stress,
    # and recycle paths, so ignore their potential evaluation errors.
    DiagnosticsToolbox(m).assert_no_structural_warnings(
        ignore_evaluation_errors=scope in OUTLET_P80_SCOPES
    )
