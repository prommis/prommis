"""Tests for the mineral-by-size WHIMS separator.

The input and magnetics-product arrays are transcribed from the Met Dynamics
WHIMS example [3]. The fitted M50 is specific to the current reconstructed
force expression; it is not a value reported by that source.

References:
    [1] King, R. P. "8 - Magnetic Separation." In Modeling and Simulation of
        Mineral Processing Systems. Butterworth-Heinemann, 2001.
        https://doi.org/10.1016/B978-0-08-051184-9.50012-2
    [2] Dobby, G. and J. A. Finch. "Capture of Mineral Particles in a High
        Gradient Magnetic Field." Powder Technology 17, no. 1 (1977): 73-82.
        https://doi.org/10.1016/0032-5910(77)85044-4
    [3] Met Dynamics. "WHIMS (Dobby and Finch)."
        https://wiki.metdynamics.com.au/view/WHIMS_(Dobby_and_Finch)
"""

from pyomo.environ import ConcreteModel, value

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import get_scaling_factor

import pytest
from prommis.psd.properties.mineral_size_psd import MineralSizePSDParameterBlock

from prommis.magnetic_separator.whims_separator import (
    WHIMSSeparator,
    WHIMSSeparatorScaler,
)

MINERALS = ("Ore1", "Ore2", "Ore3", "Ore4", "Ore5")
FEED_TPH = (
    (2.02, 1.55, 2.33, 6.98, 2.64),
    (0.79, 0.61, 0.92, 2.75, 1.04),
    (2.37, 1.82, 2.73, 8.19, 3.09),
    (3.37, 2.59, 3.89, 11.66, 4.40),
    (3.52, 2.71, 4.07, 12.20, 4.61),
    (0.94, 0.72, 1.08, 3.24, 1.22),
)
MAGNETICS_PRODUCT_TPH = (
    (1.57, 0.40, 0.01, 0.03, 0.01),
    (0.79, 0.31, 0.08, 0.03, 0.01),
    (2.37, 1.11, 0.53, 0.16, 0.06),
    (3.37, 1.84, 1.14, 1.04, 0.13),
    (3.52, 2.15, 1.53, 2.11, 0.18),
    (0.94, 0.64, 0.51, 0.86, 0.06),
)
SIZE_EDGES_MM = (0.0, 0.011, 0.015, 0.021, 0.028, 0.035, 0.050)
ORE_DENSITIES_T_M3 = (4.60, 4.60, 4.50, 4.50, 4.30)
ORE_SUSCEPTIBILITIES_M3_KG = (5.9e-2, 3.7e-3, 4.1e-4, 1.4e-4, 0.0)
PHYSICAL_RECOVERY = (0.01, 0.01, 0.02, 0.03, 0.04, 0.05)
# Effective fit for this example using the reconstructed magnetic-force term
# with G=1; it is not a vendor-calibrated parameter.
FITTED_M50 = 5.086829e-10


@pytest.mark.unit
def test_whims_scaling_routines():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.properties = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 4e-3],
        component_list=["Ore1", "Ore2"],
    )
    model.fs.unit = WHIMSSeparator(property_package=model.fs.properties)
    unit = model.fs.unit
    scaler = WHIMSSeparatorScaler()

    scaler.variable_scaling_routine(unit)

    feed_flow = unit.feed_state[0].flow_mass_size_comp[0, "Ore1"]
    feed_factor = get_scaling_factor(feed_flow)
    assert feed_factor is not None and feed_factor > 0
    for size in model.fs.properties.size_interval_set:
        for mineral in model.fs.properties.solid_component_set:
            assert get_scaling_factor(unit.recovery[0, size, mineral]) == pytest.approx(
                scaler.RECOVERY_SCALING_FACTOR
            )
            assert get_scaling_factor(
                unit.recovery_M[0, size, mineral]
            ) == pytest.approx(scaler.RECOVERY_SCALING_FACTOR)

    scaler.constraint_scaling_routine(unit)

    for size in model.fs.properties.size_interval_set:
        for mineral in model.fs.properties.solid_component_set:
            assert get_scaling_factor(
                unit.mags_solid_balance[0, size, mineral]
            ) == pytest.approx(
                get_scaling_factor(
                    unit.feed_state[0].flow_mass_size_comp[size, mineral]
                )
            )
            assert get_scaling_factor(
                unit.nonmags_solid_balance[0, size, mineral]
            ) == pytest.approx(
                get_scaling_factor(
                    unit.feed_state[0].flow_mass_size_comp[size, mineral]
                )
            )
            assert get_scaling_factor(
                unit.recovery_M_eqn[0, size, mineral]
            ) == pytest.approx(1.0 / scaler.RECOVERY_SCALING_FACTOR)
            assert get_scaling_factor(
                unit.recovery_eqn[0, size, mineral]
            ) == pytest.approx(1.0 / scaler.RECOVERY_SCALING_FACTOR)


def test_whims_separator_with_dobby_finch_example_data():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.properties = MineralSizePSDParameterBlock(
        size_edges=[edge * 1e-3 for edge in SIZE_EDGES_MM],
        bottom_size=0.0055e-3,
        component_list=list(MINERALS),
    )
    model.fs.unit = WHIMSSeparator(property_package=model.fs.properties)
    unit = model.fs.unit

    for parameter, setting in {
        "H": 0.210,
        "Hs": 0.210,
        "u": 0.156,
        "Lm": 0.084,
        "M50": FITTED_M50,
        "B": 0.413,
        "beta": 1.044,
    }.items():
        getattr(unit, parameter).set_value(setting)

    for mineral, density, susceptibility in zip(
        MINERALS, ORE_DENSITIES_T_M3, ORE_SUSCEPTIBILITIES_M3_KG
    ):
        for size in model.fs.properties.size_interval_set:
            # The example gives one density per ore, so reuse it for each size.
            unit.rho[size, mineral].set_value(density * 1000)
            unit.chi[size, mineral].set_value(susceptibility)

    for size, recovery in enumerate(PHYSICAL_RECOVERY):
        unit.Rp[size].set_value(recovery)
        for mineral, flow_tph in zip(MINERALS, FEED_TPH[size]):
            unit.feed_state[0].flow_mass_size_comp[size, mineral].fix(
                flow_tph * 1000 / 3600
            )

    unit.default_initializer().initialize(unit)

    for mineral, expected_tph in zip(MINERALS, (13.01, 10.00, 15.02, 45.02, 17.00)):
        assert value(unit.feed_state[0].flow_mass_comp[mineral]) == pytest.approx(
            expected_tph * 1000 / 3600
        )

    for size in model.fs.properties.size_interval_set:
        for mineral, expected_magnetics_tph in zip(
            MINERALS, MAGNETICS_PRODUCT_TPH[size]
        ):
            feed = value(unit.feed_state[0].flow_mass_size_comp[size, mineral])
            magnetics = value(unit.mags_state[0].flow_mass_size_comp[size, mineral])
            non_magnetics = value(
                unit.nonmags_state[0].flow_mass_size_comp[size, mineral]
            )
            recovery = value(unit.recovery[0, size, mineral])
            assert recovery == pytest.approx(magnetics / feed if feed else 0.0)
            assert magnetics + non_magnetics == pytest.approx(feed)
            assert magnetics == pytest.approx(
                expected_magnetics_tph * 1000 / 3600,
                abs=0.27 * 1000 / 3600,
            )

    total_feed_tph = sum(sum(row) for row in FEED_TPH)
    assert value(unit.feed_state[0].flow_mass) == pytest.approx(
        total_feed_tph * 1000 / 3600
    )
    assert value(unit.mags_state[0].flow_mass) + value(
        unit.nonmags_state[0].flow_mass
    ) == pytest.approx(value(unit.feed_state[0].flow_mass))
