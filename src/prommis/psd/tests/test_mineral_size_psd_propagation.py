"""Stream-propagation checks for the mineral-by-size PSD package."""

from pyomo.environ import ConcreteModel, TransformationFactory, value
from pyomo.network import Arc, Port

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import propagate_state

import pytest

from prommis.psd.properties.mineral_size_psd import MineralSizePSDParameterBlock


def _passthrough_model(edges_m, components, bottom_size=None):
    """Create two same-mesh state blocks connected by an Arc."""
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.props = MineralSizePSDParameterBlock(
        size_edges=edges_m,
        component_list=components,
        bottom_size=bottom_size,
    )
    model.fs.feed = model.fs.props.build_state_block(model.fs.time, defined_state=True)
    model.fs.product = model.fs.props.build_state_block(
        model.fs.time, defined_state=True
    )
    time = model.fs.time.first()
    model.fs.feed[time].params.assert_same_mesh(model.fs.product[time].params)

    for name, state in (
        ("feed", model.fs.feed[time]),
        ("product", model.fs.product[time]),
    ):
        port = Port()
        port.add(state.flow_mass_size_comp, "flow_mass_size_comp")
        port.add(state.flow_mass_liquid, "flow_mass_liquid")
        port.add(state.temperature, "temperature")
        port.add(state.pressure, "pressure")
        setattr(model.fs, f"{name}_port", port)

    model.fs.stream = Arc(source=model.fs.feed_port, destination=model.fs.product_port)
    TransformationFactory("network.expand_arcs").apply_to(model)
    return model


@pytest.mark.unit
def test_joint_flow_arc_propagation_and_derived_marginals():
    model = _passthrough_model([1e-3, 2e-3, 4e-3, 8e-3], ["Ore1", "Ore2"])
    time = model.fs.time.first()
    feed = model.fs.feed[time]
    product = model.fs.product[time]
    flows = {
        (0, "Ore1"): 1.0,
        (1, "Ore1"): 2.0,
        (2, "Ore1"): 3.0,
        (0, "Ore2"): 3.0,
        (1, "Ore2"): 2.0,
        (2, "Ore2"): 1.0,
    }
    for index, flow in flows.items():
        feed.flow_mass_size_comp[index].fix(flow)
    feed.flow_mass_liquid["H2O"].fix(7.0)
    feed.temperature.fix(320.0)
    feed.pressure.fix(120000.0)

    propagate_state(arc=model.fs.stream)

    for index, flow in flows.items():
        assert value(product.flow_mass_size_comp[index]) == pytest.approx(flow)
    for mineral in ("Ore1", "Ore2"):
        assert value(product.flow_mass_comp[mineral]) == pytest.approx(6.0)
    for size in range(3):
        assert value(product.flow_mass_size[size]) == pytest.approx(4.0)
    assert value(product.flow_mass_liquid["H2O"]) == pytest.approx(7.0)
    assert value(product.temperature) == pytest.approx(320.0)
    assert value(product.pressure) == pytest.approx(120000.0)
    assert value(product.cum_passing_mineral[0, "Ore1"]) == pytest.approx(1 / 6)
    assert value(product.cum_passing_mineral[1, "Ore2"]) == pytest.approx(5 / 6)


@pytest.mark.unit
def test_degenerate_joint_distribution_propagates_unchanged():
    model = _passthrough_model([1e-3, 2e-3, 4e-3, 8e-3], ["Ore1"])
    time = model.fs.time.first()
    feed = model.fs.feed[time]
    product = model.fs.product[time]
    for size, flow in enumerate((0.0, 10.0, 0.0)):
        feed.flow_mass_size_comp[size, "Ore1"].fix(flow)

    propagate_state(arc=model.fs.stream)

    for size, flow in enumerate((0.0, 10.0, 0.0)):
        assert value(product.flow_mass_size_comp[size, "Ore1"]) == pytest.approx(
            flow, abs=1e-8
        )
    assert value(product.flow_mass_comp["Ore1"]) == pytest.approx(10.0)
    assert value(product.flow_mass_size[1]) == pytest.approx(10.0)
    assert value(product.cum_passing_mineral[0, "Ore1"]) == pytest.approx(0.0)
    assert value(product.cum_passing_mineral[1, "Ore1"]) == pytest.approx(1.0)


@pytest.mark.unit
def test_assert_same_mesh_rejects_mismatched_mesh():
    model = ConcreteModel()
    model.params = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 4e-3], component_list=["Ore1"]
    )
    model.other_params = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 5e-3], component_list=["Ore1"]
    )

    with pytest.raises(ConfigurationError):
        model.params.assert_same_mesh(model.other_params)

    model.same_params = MineralSizePSDParameterBlock(
        size_edges=[1e-3, 2e-3, 4e-3], component_list=["Ore1"]
    )
    model.params.assert_same_mesh(model.same_params)
