#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Test the size-reducer bases and helpers."""

import pytest

from pyomo.common.config import ConfigDict, ConfigValue

from pyomo.environ import ConcreteModel
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.scaling.util import get_scaling_factor, set_scaling_factor

from prommis.comminution.properties.solid_psd_properties import (
    SolidPSDParameterBlock,
)
from prommis.comminution.unit_models import _size_reducer_common as sr
from prommis.comminution.unit_models import _solid_psd_common as spc

EDGES = [1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3]
SIZED = ["OreA", "OreB"]
UNSIZED = ["Inert"]
DENSITY = {"OreA": 3000.0, "OreB": 5000.0, "Inert": 2000.0}


@declare_process_block_class("_StubUnit")
class _StubUnitData(UnitModelBlockData):
    """Identity unit exercising shared state-block and pass-through helpers."""

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "property_package",
        ConfigValue(default=useDefault, domain=is_physical_parameter_block),
    )
    CONFIG.declare("property_package_args", ConfigDict(implicit=True))

    def build(self):
        super().build()
        spc.require_solid_psd_package(self.config, "_StubUnit")
        sr.add_inlet_outlet_state_blocks(self)
        sr.add_passthrough_constraints(self)

        pp = self.config.property_package

        @self.Constraint(
            self.flowsheet().time,
            pp.sized_solid_list,
            pp.size_interval_set,
            doc="Outlet sized-solid flow equals inlet flow.",
        )
        def product_rows_eqn(b, t, s, k):
            return (
                b.properties_out[t].flow_mass_sized_comp_size[s, k]
                == b.properties_in[t].flow_mass_sized_comp_size[s, k]
            )


def _model(vapors=()):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    kwargs = dict(
        size_edges=EDGES,
        sized_solid_component_list=list(SIZED),
        unsized_solid_component_list=list(UNSIZED),
        solid_density=dict(DENSITY),
    )
    if vapors:
        kwargs["vapor_component_list"] = list(vapors)
    m.fs.params = SolidPSDParameterBlock(**kwargs)
    m.fs.unit = _StubUnit(property_package=m.fs.params)
    return m


def _fix_feed(feed):
    for s in SIZED:
        for k in range(3):
            feed.flow_mass_sized_comp_size[s, k].fix(0.5 + 0.1 * k)
    feed.flow_mass_unsized_comp["Inert"].fix(0.2)
    feed.flow_mass_liquid_comp["H2O"].fix(1.5)
    feed.temperature.fix(298.15)
    feed.pressure.fix(101325.0)


@pytest.mark.unit
def test_scaler_declared_factor_propagates_to_constraint_factor():
    m = _model(vapors=("Air",))
    u = m.fs.unit
    feed = u.properties_in[0]
    _fix_feed(feed)  # sized 3.6 kg/s, Inert 0.2, H2O 1.5
    feed.flow_mass_vapor_comp["Air"].fix(0.7)
    assert degrees_of_freedom(m) == 0
    assert_units_consistent(m)
    scaler = sr.SizeReducerScalerBase()
    # Constraint scaling needs variable factors, even when the feed is fixed.
    with pytest.raises(ConfigurationError, match="is missing or non-positive"):
        scaler.constraint_scaling_routine(u)
    set_scaling_factor(feed.flow_mass_liquid_comp["H2O"], 0.5)  # 2.0 kg/s nominal
    scaler.variable_scaling_routine(u)
    scaler.constraint_scaling_routine(u)
    # The water row uses the assigned 2.0 kg/s nominal, not its 1.5 kg/s flow.
    for con, expected in (
        (u.liquid_passthrough_eqn[0, "H2O"], 0.5),
        (u.unsized_passthrough_eqn[0, "Inert"], 1.0 / 0.2),
        (u.vapor_passthrough_eqn[0, "Air"], 1.0 / 0.7),
    ):
        assert get_scaling_factor(con) == pytest.approx(expected, rel=1e-12), con.name
    assert get_scaling_factor(u.temperature_eqn[0]) == pytest.approx(
        1.0 / 300.0, rel=1e-12, abs=0
    )
    assert get_scaling_factor(u.pressure_eqn[0]) == pytest.approx(
        1.0e-5, rel=1e-12, abs=0
    )
