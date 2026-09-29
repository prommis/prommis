#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Stream-propagation and pass-through validation for the bulk PSD package.

Reusable "fix feed -> propagate -> compare" pattern that the Phase 1 stages copy.
A pass-through harness (no shipped pass-through unit) connects two state blocks
with an Arc; the downstream block is ``defined_state=True`` (mirroring a unit
inlet) so the Arc equalities make a square system without a redundant consistency
constraint.
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.network import Arc, Port

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import propagate_state

from prommis.psd.properties.bulk_psd import BulkPSDParameterBlock

# Pass-through / degenerate feeds place many PSD variables exactly at their 0
# lower bound; with default IPOPT bound_push these get perturbed off the bound by
# ~1e-6, swamping a tight post-solve comparison.  Tightening bound_push/bound_frac
# lets them converge to the bound, so the pass-through agreement is exact.
solver = get_solver()
solver.options["bound_push"] = 1e-10
solver.options["bound_frac"] = 1e-10


def _passthrough_model(edges_m, components, bottom_size=None):
    """Two Arc-connected state blocks (both defined_state=True) sharing a mesh."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = BulkPSDParameterBlock(
        size_edges=edges_m, component_list=components, bottom_size=bottom_size
    )
    m.fs.feed = m.fs.props.build_state_block(m.fs.time, defined_state=True)
    m.fs.prod = m.fs.props.build_state_block(m.fs.time, defined_state=True)
    t = m.fs.time.first()
    # harness documents the pattern: connected blocks must share the mesh
    m.fs.feed[t].params.assert_same_mesh(m.fs.prod[t].params)

    m.fs.feed_port = Port()
    m.fs.feed_port.add(m.fs.feed[t].flow_mass_comp, "flow_mass_comp")
    m.fs.feed_port.add(m.fs.feed[t].flow_mass_size, "flow_mass_size")
    m.fs.prod_port = Port()
    m.fs.prod_port.add(m.fs.prod[t].flow_mass_comp, "flow_mass_comp")
    m.fs.prod_port.add(m.fs.prod[t].flow_mass_size, "flow_mass_size")
    m.fs.stream = Arc(source=m.fs.feed_port, destination=m.fs.prod_port)
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


# -----------------------------------------------------------------------------
@pytest.mark.component
@pytest.mark.solver
def test_arc_propagation_handcalc():
    edges = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]
    m = _passthrough_model(edges, ["Ore1", "Ore2", "Ore3"])
    t = m.fs.time.first()
    feed = m.fs.feed[t]
    prod = m.fs.prod[t]
    for j, val in zip(["Ore1", "Ore2", "Ore3"], [5.0, 3.0, 2.0]):
        feed.flow_mass_comp[j].fix(val)
    for k, val in zip([0, 1, 2, 3], [1.0, 2.0, 3.0, 4.0]):
        feed.flow_mass_size[k].fix(val)

    propagate_state(arc=m.fs.stream)
    results = solver.solve(m)
    assert_optimal_termination(results)
    for j in ["Ore1", "Ore2", "Ore3"]:
        assert value(prod.flow_mass_comp[j]) == pytest.approx(
            value(feed.flow_mass_comp[j]), rel=1e-8
        )
    for k in [0, 1, 2, 3]:
        assert value(prod.flow_mass_size[k]) == pytest.approx(
            value(feed.flow_mass_size[k]), rel=1e-8
        )


@pytest.mark.component
@pytest.mark.solver
def test_degenerate_propagation():
    # all mass in one interval propagates unchanged; P80 equals the augmented-curve value
    edges = [1e-3, 2e-3, 4e-3, 8e-3, 16e-3]
    m = _passthrough_model(edges, ["Ore1"])
    t = m.fs.time.first()
    feed = m.fs.feed[t]
    prod = m.fs.prod[t]
    feed.flow_mass_comp["Ore1"].fix(10.0)
    for k, val in zip([0, 1, 2, 3], [0.0, 10.0, 0.0, 0.0]):  # all in bin 1
        feed.flow_mass_size[k].fix(val)
    propagate_state(arc=m.fs.stream)
    results = solver.solve(m)
    assert_optimal_termination(results)
    for k in [0, 1, 2, 3]:
        assert value(prod.flow_mass_size[k]) == pytest.approx(
            value(feed.flow_mass_size[k]), abs=1e-8
        )
    # step distribution: cum = [0,1,1,1]; P80 in [x_1, x_2] = 2 + 0.8*(4-2) = 3.6 mm
    assert value(units.convert(prod.P80, units.mm)) == pytest.approx(3.6, rel=2e-3)


@pytest.mark.unit
def test_assert_same_mesh_negative():
    m = ConcreteModel()
    m.p1 = BulkPSDParameterBlock(size_edges=[1e-3, 2e-3, 4e-3], component_list=["A"])
    # same interval count, different top edge
    m.p2 = BulkPSDParameterBlock(size_edges=[1e-3, 2e-3, 5e-3], component_list=["A"])
    with pytest.raises(ConfigurationError):
        m.p1.assert_same_mesh(m.p2)
    # identical mesh passes
    m.p3 = BulkPSDParameterBlock(size_edges=[1e-3, 2e-3, 4e-3], component_list=["A"])
    m.p1.assert_same_mesh(m.p3)
