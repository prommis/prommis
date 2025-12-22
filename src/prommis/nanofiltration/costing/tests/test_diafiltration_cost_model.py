#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests the diafiltration cost model
"""

from pyomo.environ import ConcreteModel, SolverFactory, Var, units

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import pytest

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCosting,
    DiafiltrationCostingData,
)


@pytest.mark.component
def test_diafiltration_costing():
    # create a model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = DiafiltrationCosting()

    # create unit models to cost
    m.fs.membrane = UnitModelBlock()
    m.fs.membrane.length = Var(
        initialize=2400,
        doc="Total membrane wound length",
        units=units.m,
    )
    m.fs.membrane.length.fix()
    m.fs.membrane.width = Var(
        initialize=1.5,
        doc="Membrane width (module/tube length)",
        units=units.m,
    )
    m.fs.membrane.width.fix()
    m.fs.membrane.water_flux = Var(
        initialize=0.1,
        doc="Water flux through membrane",
        units=units.m**3 / units.m**2 / units.h,
    )
    m.fs.membrane.water_flux.fix()

    m.fs.cascade = UnitModelBlock()
    m.fs.cascade.feed_volume_flow = Var(
        initialize=100,
        doc="Volumetric flow rate of feed",
        units=units.m**3 / units.h,
    )
    m.fs.cascade.feed_volume_flow.fix()
    m.fs.cascade.permeate_volume_flow = Var(
        initialize=110,
        doc="Volumetric flow rate of permeate",
        units=units.m**3 / units.h,
    )
    m.fs.cascade.permeate_volume_flow.fix()

    m.fs.pump = UnitModelBlock()
    m.fs.pump.atmospheric_pressure = Var(
        initialize=101.325,
        doc="Atmospheric pressure in kilo Pascal",
        units=units.kPa,
    )
    m.fs.pump.atmospheric_pressure.fix()
    m.fs.pump.operating_pressure = Var(
        initialize=145,
        doc="Membrane operating pressure",
        units=units.psi,
    )
    m.fs.pump.operating_pressure.fix()

    m.fs.precipitator = UnitModelBlock()
    m.fs.precipitator.V = Var(
        initialize=20,
        doc="Precipitator volume",
        units=units.m**3,
    )
    m.fs.precipitator.V.fix()

    m.fs.membrane.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.membrane.length,
            "membrane_width": m.fs.membrane.width,
        },
    )
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop_utility,
        costing_method_arguments={
            "water_flux": m.fs.membrane.water_flux,
            "vol_flow_feed": m.fs.cascade.feed_volume_flow,
            "vol_flow_perm": m.fs.cascade.permeate_volume_flow,
        },
    )
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.fs.pump.atmospheric_pressure,
            "outlet_pressure": m.fs.pump.operating_pressure,
            "inlet_vol_flow": m.fs.cascade.feed_volume_flow,
        },
    )
    m.fs.precipitator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_precipitator,
        costing_method_arguments={"precip_volume": m.fs.precipitator.V},
    )

    # set bound to numerically zero to avoid potential log(0)
    m.fs.precipitator.costing.weight.setlb(10 ** (-4))

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    m.fs.costing.cost_process()

    # unfix the precipitator diameter such that both L and D are manipulated to satisfy V
    m.fs.precipitator.costing.precipitator_diameter.unfix()

    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    dt.assert_no_numerical_warnings()


@pytest.mark.component
def test_simple_costing():
    # create a model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = DiafiltrationCosting()

    # create unit models to cost
    m.fs.membrane = UnitModelBlock()
    m.fs.membrane.length = Var(
        initialize=2400,
        doc="Total membrane wound length",
        units=units.m,
    )
    m.fs.membrane.length.fix()
    m.fs.membrane.width = Var(
        initialize=1.5,
        doc="Membrane width (module/tube length)",
        units=units.m,
    )
    m.fs.membrane.width.fix()
    m.fs.membrane.water_flux = Var(
        initialize=0.1,
        doc="Water flux through membrane",
        units=units.m**3 / units.m**2 / units.h,
    )
    m.fs.membrane.water_flux.fix()

    m.fs.cascade = UnitModelBlock()
    m.fs.cascade.feed_volume_flow = Var(
        initialize=100,
        doc="Volumetric flow rate of feed",
        units=units.m**3 / units.h,
    )
    m.fs.cascade.feed_volume_flow.fix()
    m.fs.cascade.permeate_volume_flow = Var(
        initialize=110,
        doc="Volumetric flow rate of permeate",
        units=units.m**3 / units.h,
    )
    m.fs.cascade.permeate_volume_flow.fix()

    m.fs.pump = UnitModelBlock()
    m.fs.pump.atmospheric_pressure = Var(
        initialize=101.325,
        doc="Atmospheric pressure in kilo Pascal",
        units=units.kPa,
    )
    m.fs.pump.atmospheric_pressure.fix()
    m.fs.pump.operating_pressure = Var(
        initialize=145,
        doc="Membrane operating pressure",
        units=units.psi,
    )
    m.fs.pump.operating_pressure.fix()

    m.fs.precipitator = UnitModelBlock()
    m.fs.precipitator.V = Var(
        initialize=50,
        doc="Precipitator volume",
        units=units.m**3,
    )
    m.fs.precipitator.V.fix()

    m.fs.membrane.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membranes,
        costing_method_arguments={
            "membrane_length": m.fs.membrane.length,
            "membrane_width": m.fs.membrane.width,
        },
    )
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop_utility,
        costing_method_arguments={
            "water_flux": m.fs.membrane.water_flux,
            "vol_flow_feed": m.fs.cascade.feed_volume_flow,
            "vol_flow_perm": m.fs.cascade.permeate_volume_flow,
        },
    )
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.fs.pump.atmospheric_pressure,
            "outlet_pressure": m.fs.pump.operating_pressure,
            "inlet_vol_flow": m.fs.cascade.feed_volume_flow,
            "simple_costing": True,
        },
    )
    m.fs.precipitator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_precipitator,
        costing_method_arguments={
            "precip_volume": m.fs.precipitator.V,
            "simple_costing": True,
        },
    )

    # fix pump installation power (design decision) for simulation
    m.fs.pump.costing.pump_installation_power_simple.fix()

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()

    m.fs.costing.cost_process()
    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    dt.assert_no_numerical_warnings()
