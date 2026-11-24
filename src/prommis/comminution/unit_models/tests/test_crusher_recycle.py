#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for Whiten classification recycle, including loads and product flows."""

import pytest

from pyomo.environ import value

from idaes.core.initialization import InitializationStatus
from idaes.core.util.exceptions import ConfigurationError

from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.tests.crusher_test_support import (
    CLASSIFICATION_RECYCLE_REF,
    EDGES_X,
    FEED_X,
    crusher,
    fix_feed,
    flowsheet,
    one_mineral_flowsheet,
    x_individual,
)


def coarsest_recycled_unit(breakage):
    """Build a three-interval Whiten crusher with full coarsest-bin recycle."""
    m = one_mineral_flowsheet(edges=[1e-3, 2e-3, 4e-3, 8e-3])
    cr = crusher(
        m,
        psd_method="selection_breakage",
        recycle="classification_recycle",
        selection_function="whiten",
        selection_params={"K1": 3e-3, "K2": 4e-3, "K3": 2.3},
        breakage_function="user",
        breakage_params=breakage,
    )
    return m, cr


@pytest.mark.unit
def test_recycle_rejects_singular_denominator():
    # The coarsest interval is fully recycled. If it retains all its mass,
    # its recycle denominator 1 - b_22 C_2 is zero.
    with pytest.raises(
        ConfigurationError,
        match=r"Whiten recycle for 'OreA': 1 - b_jj C_j = \S+ at interval 2 "
        r"must be at least",
    ):
        coarsest_recycled_unit([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    m, cr = coarsest_recycled_unit([[1.0, 0.0, 0.5], [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]])
    assert [value(cr.classification[k]) for k in range(3)] == [0.0, 0.0, 1.0]


@pytest.mark.component
@pytest.mark.solver
def test_classification_recycle_per_mineral_loads_solved():
    m = flowsheet(edges=EDGES_X)
    cr = crusher(m, bond_work_index=14.0, **x_individual("classification_recycle"))
    fix_feed(cr.properties_in[0], FEED_X)
    native.CrusherSolidPSDInitializer().initialize(cr)
    for s, loads, products in (
        ("OreA", [40 / 13, 64 / 13, 128 / 13], [40 / 13, 32 / 13, 32 / 13]),
        ("OreB", [72 / 35, 96 / 35, 128 / 7], [72 / 35, 48 / 35, 32 / 7]),
    ):
        for k in range(3):
            assert value(cr.recycle_load[0, s, k]) == pytest.approx(loads[k], rel=1e-8)
            assert value(
                cr.properties_out[0].flow_mass_sized_comp_size[s, k]
            ) == pytest.approx(products[k], rel=1e-8)


@pytest.mark.component
@pytest.mark.solver
def test_zero_mineral_row_classification_recycle_matches_one_species_control():
    # OreB has a different density but no feed.
    subject_m = flowsheet()
    subject = crusher(subject_m, bond_work_index=14.0, **CLASSIFICATION_RECYCLE_REF)
    rows = {"OreA": [0.5, 1.0, 2.0, 6.5], "OreB": [0.0, 0.0, 0.0, 0.0]}
    fix_feed(subject.properties_in[0], rows, liquid=1.0)
    control_m = flowsheet(species=["OreA"], density={"OreA": 2800.0})
    control = crusher(control_m, bond_work_index=14.0, **CLASSIFICATION_RECYCLE_REF)
    fix_feed(control.properties_in[0], {"OreA": rows["OreA"]}, liquid=1.0)
    for unit in (subject, control):
        init = native.CrusherSolidPSDInitializer()
        init.initialize(unit)
        assert init.summary[unit]["status"] == InitializationStatus.Ok
    for k in range(4):
        assert (
            abs(value(subject.properties_out[0].flow_mass_sized_comp_size["OreB", k]))
            <= 1e-12
        )
        assert abs(value(subject.recycle_load[0, "OreB", k])) <= 1e-12
    assert value(subject.product_p80[0]) > 0
    assert value(subject.product_p80[0]) == pytest.approx(
        value(control.product_p80[0]), rel=1e-8
    )
    assert value(subject.power[0]) > 0
    assert value(subject.power[0]) == pytest.approx(value(control.power[0]), rel=1e-8)
