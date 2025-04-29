#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests the diafiltration cost model
"""

from pyomo.environ import ConcreteModel, Param, SolverFactory, Var, units

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

    # define global parameters
    m.w = Param(
        initialize=1.5,
        mutable=True,
        doc="Membrane width (module/tube length)",
        units=units.m,
    )
    m.Jw = Param(
        initialize=0.1,
        mutable=True,
        doc="Water flux",
        units=units.m**3 / units.m**2 / units.h,
    )
    m.atmospheric_pressure = Param(
        initialize=101325,
        doc="Atmospheric pressure in Pascal",
        units=units.Pa,
    )
    m.operating_pressure = Param(
        initialize=145,
        mutable=True,
        doc="Membrane operating pressure",
        units=units.psi,
    )
    m.Q_feed = Param(
        initialize=100,
        mutable=True,
        doc="Cascade feed flow",
        units=units.m**3 / units.h,
    )
    m.Q_perm = Param(
        initialize=110,
        mutable=True,
        doc="Cascade permeate flow",
        units=units.m**3 / units.h,
    )

    # create unit models to cost
    m.fs.membrane = UnitModelBlock()
    m.fs.membrane.length = Var(
        initialize=2400,
        doc="Total membrane wound length",
        units=units.m,
    )
    m.fs.membrane.length.fix()
    m.fs.cascade = UnitModelBlock()
    m.fs.pump = UnitModelBlock()
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
            "membrane_width": m.w,
        },
    )
    m.fs.cascade.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop,
        costing_method_arguments={
            "water_flux": m.Jw,
            "vol_flow_feed": m.Q_feed,
            "vol_flow_perm": m.Q_perm,
        },
    )
    m.fs.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=DiafiltrationCostingData.cost_pump,
        costing_method_arguments={
            "inlet_pressure": m.atmospheric_pressure,
            "outlet_pressure": m.operating_pressure,
            "inlet_vol_flow": m.Q_feed,
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
