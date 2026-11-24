#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Shared builders, test data, and reference arithmetic for crusher tests."""

from pyomo.environ import ConcreteModel, units, value

from idaes.core import FlowsheetBlock

from prommis.comminution.core.size_mesh import geometric_series
from prommis.comminution.properties.solid_psd_properties import (
    SolidPSDParameterBlock,
)
from prommis.comminution.unit_models import crusher as native
from prommis.comminution.unit_models.crusher import CrusherSolidPSD

# Reference mesh: four intervals from 1 to 16 mm.
EDGES_REF = geometric_series(16e-3, 1e-3, 2)
# Test mesh with geometric-mean interval sizes of 1, 3, and 4 mm.
EDGES_X = [0.4e-3, 2.5e-3, 3.6e-3, (40.0 / 9.0) * 1e-3]
SPECIES_TWO = ["OreA", "OreB"]
DENSITY_TWO = {"OreA": 2800.0, "OreB": 4200.0}


def flowsheet(
    edges=None, species=None, density=None, unsized=None, vapors=None, bottom_size=None
):
    """A flowsheet with one solid PSD property package."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    kwargs = dict(
        size_edges=list(EDGES_REF if edges is None else edges),
        sized_solid_component_list=list(SPECIES_TWO if species is None else species),
        solid_density=dict(DENSITY_TWO if density is None else density),
    )
    if unsized:
        kwargs["unsized_solid_component_list"] = list(unsized)
    if vapors:
        kwargs["vapor_component_list"] = list(vapors)
    if bottom_size is not None:
        kwargs["bottom_size"] = bottom_size
    m.fs.pp = SolidPSDParameterBlock(**kwargs)
    return m


def one_mineral_flowsheet(edges=None):
    """A flowsheet whose package carries the one component OreA."""
    return flowsheet(edges=edges, species=["OreA"], density={"OreA": 2800.0})


def two_time_flowsheet(edges=None, species=None, density=None):
    """A flowsheet on the two-point time set [0.0, 1.0]."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0.0, 1.0], time_units=units.s)
    m.fs.pp = SolidPSDParameterBlock(
        size_edges=list(EDGES_REF if edges is None else edges),
        sized_solid_component_list=list(species or ["OreA"]),
        solid_density=dict(density or {"OreA": 2800.0}),
    )
    return m


def crusher(m, **cfg):
    """Build ``m.fs.cr`` on the flowsheet's package with ``cfg``."""
    m.fs.cr = CrusherSolidPSD(property_package=m.fs.pp, **cfg)
    return m.fs.cr


def fix_feed(
    feed, rows, liquid=0.0, unsized=None, temperature=298.15, pressure=101325.0
):
    """Fix a feed from ``rows = {component: [kg/s per interval]}``."""
    for s, row in rows.items():
        for k, v in enumerate(row):
            feed.flow_mass_sized_comp_size[s, k].fix(v)
    for j in feed.params.liquid_list:
        feed.flow_mass_liquid_comp[j].fix(liquid)
    for u, v in (unsized or {}).items():
        feed.flow_mass_unsized_comp[u].fix(v)
    feed.temperature.fix(temperature)
    feed.pressure.fix(pressure)


WHITEN_REF = {"K1": 1.0e-3, "K2": 16.0e-3, "K3": 2.3}
LUCKIE_AUSTIN_REF = {"phi": 0.4, "gamma": 0.9, "beta": 4.0}
SINGLE_STRESS_REF = dict(
    psd_method="selection_breakage",
    recycle="single_pass",
    n_stress_events=1,
    selection_function="whiten",
    selection_params=WHITEN_REF,
    breakage_function="luckie_austin",
    breakage_params=LUCKIE_AUSTIN_REF,
)
CLASSIFICATION_RECYCLE_REF = dict(SINGLE_STRESS_REF, recycle="classification_recycle")

RR = native.DistributionShape.rosin_rammler
GGS = native.DistributionShape.gates_gaudin_schuhmann

# one configuration per product path, on the reference mesh
PATH_CONFIGS = {
    "tabular": dict(psd_method="tabular", tabular_psd=[0.25, 0.5, 0.25, 0.0]),
    "distribution_size_spec": dict(
        psd_method="distribution_function", shape_exponent=1.5, target_product_p80=5e-3
    ),
    "distribution_power_spec": dict(
        psd_method="distribution_function", shape_exponent=1.5
    ),
    "single_stress": dict(SINGLE_STRESS_REF),
    "repeated_stress": dict(SINGLE_STRESS_REF, n_stress_events=3),
    "classification_recycle": dict(CLASSIFICATION_RECYCLE_REF),
}

# Synthetic drop-weight data for pendulum power tests.
PENDULUM_APPEARANCE = {
    "10": {"t75": 1.0, "t50": 1.5, "t25": 3.0, "t4": 20.0, "t2": 40.0},
    "20": {"t75": 2.0, "t50": 3.0, "t25": 6.0, "t4": 35.0, "t2": 60.0},
}
PENDULUM_ECS_A = {
    "sizes": [2e-3, 8e-3],
    "t10_rows": [10, 20, 30],
    "ecs": [[0.4, 0.3], [0.8, 0.6], [1.1, 0.9]],
}
PENDULUM_REF = dict(
    CLASSIFICATION_RECYCLE_REF,
    breakage_function="t10_appearance",
    breakage_params={
        "appearance_table": PENDULUM_APPEARANCE,
        "t10": 15.0,
        "ecs_table": PENDULUM_ECS_A,
    },
    power_law="pendulum",
    pendulum_power_params={"power_factor": 1.3, "no_load_power": 20.0},
)

# Three-interval mesh for hand calculations.
HC_EDGES = [1e-3, 2e-3, 4e-3, 8e-3]

# reference feed on EDGES_REF: all mass in the top interval
FEED_REF = {"OreA": [0.0, 0.0, 0.0, 6.0], "OreB": [0.0, 0.0, 0.0, 4.0]}
# Each component enters with 8 kg/s in the coarsest interval.
FEED_X = {"OreA": [0.0, 0.0, 8.0], "OreB": [0.0, 0.0, 8.0]}

# User selection and breakage data for the EDGES_X test mesh.
S_A = (0.0, 0.5, 1.0)
S_B = (0.0, 0.25, 0.5)
B_A = ((1.0, 0.5, 0.25), (0.0, 0.5, 0.5), (0.0, 0.0, 0.25))
B_B = ((1.0, 0.25, 0.125), (0.0, 0.75, 0.125), (0.0, 0.0, 0.75))
WHITEN_X = {"K1": 1.0e-3, "K2": 5.0e-3, "K3": 1.0}  # C = [0, 0.5, 0.75] on X


def x_common(path):
    """Return common inputs for ``path`` on the EDGES_X test mesh."""
    base = dict(
        psd_method="selection_breakage",
        selection_function="user",
        breakage_function="user",
    )
    if path == "single_stress":
        return dict(
            base,
            recycle="single_pass",
            n_stress_events=1,
            selection_params=list(S_A),
            breakage_params=[list(r) for r in B_A],
        )
    if path == "repeated_stress":
        return dict(
            base,
            recycle="single_pass",
            n_stress_events=2,
            selection_params=list(S_A),
            breakage_params=[list(r) for r in B_A],
        )
    return dict(
        psd_method="selection_breakage",
        recycle="classification_recycle",
        selection_function="whiten",
        selection_params=WHITEN_X,
        breakage_function="user",
        breakage_params=[list(r) for r in B_A],
    )


def x_individual(path):
    """Return EDGES_X inputs with breakage data for each component.

    Single-pass paths also use selection data for each component.
    """
    cfg = x_common(path)
    if path in ("single_stress", "repeated_stress"):
        cfg.pop("selection_params")
        cfg["selection_params_by_comp"] = {"OreA": list(S_A), "OreB": list(S_B)}
    cfg.pop("breakage_params")
    cfg["breakage_params_by_comp"] = {
        "OreA": [list(r) for r in B_A],
        "OreB": [list(r) for r in B_B],
    }
    return cfg


# A zero first edge needs a positive bottom size for characteristic-size calculations.
ZERO_EDGES = [0.0, 1.0e-3, 2.0e-3, 4.0e-3]
ZERO_BOTTOM = 0.5e-3
ZERO_WHITEN = {"K1": 0.5e-3, "K2": 4.0e-3, "K3": 2.3}

ZERO_WHITEN_BY_COMP = {
    "OreA": ZERO_WHITEN,
    "OreB": {"K1": 0.5e-3, "K2": 3.0e-3, "K3": 2.0},
}
ZERO_VOGEL_BY_COMP = {
    "OreA": {"q": 0.5, "dprime": 2.0e-3},
    "OreB": {"q": 0.7, "dprime": 1.5e-3},
}


def feed_percentile(cr, t=0):
    """The live feed 80%-passing size in metres the property package reports."""
    return value(cr.properties_in[t].percentile_size[0.8])


# Expected stress and recycle products.
def stress_reference(selection, breakage, feed, events):
    """Return the product row after ``events`` stress events.

    Each event balances the finest bin against its incoming total.
    """
    n = len(feed)
    a = [
        [
            float(k == j) * (1.0 - selection[j]) + breakage[k][j] * selection[j]
            for j in range(n)
        ]
        for k in range(n)
    ]
    row = list(feed)
    for _ in range(events):
        product = [sum(a[k][j] * row[j] for j in range(n)) for k in range(n)]
        product[0] = sum(row) - sum(product[1:])
        row = product
    return row


def recycle_reference(breakage, classification, feed):
    """Return recycle loads and product rows.

    The finest product bin balances the total product against the feed.
    """
    n = len(feed)
    load = [0.0] * n
    for k in reversed(range(n)):
        coarser = sum(
            breakage[k][j] * classification[j] * load[j] for j in range(k + 1, n)
        )
        load[k] = (feed[k] + coarser) / (1.0 - breakage[k][k] * classification[k])
    product = [(1.0 - classification[k]) * load[k] for k in range(n)]
    product[0] = sum(feed) - sum(product[1:])
    return load, product


# Shared four-interval selection and breakage data for two-component tests.
COMMON_S = [0.0, 0.3, 0.7, 1.0]
COMMON_B = [
    [1.0, 0.5, 0.25, 0.125],
    [0.0, 0.5, 0.5, 0.125],
    [0.0, 0.0, 0.25, 0.25],
    [0.0, 0.0, 0.0, 0.5],
]
