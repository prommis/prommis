#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for REE costing.

"""

import os
from contextlib import nullcontext as does_not_raise

import pyomo.environ as pyo
from pyomo.common.config import ConfigDict
from pyomo.common.dependencies import attempt_import
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.units_container import UnitsError
from pyomo.environ import assert_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    badly_scaled_var_generator,
    calculate_scaling_factors,
    unscaled_variables_generator,
)

import pytest

from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCostingData,
)
from prommis.uky.costing.custom_costing_example import CustomCostingData
from prommis.uky.costing.ree_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
    custom_REE_plant_currency_units,
)

_, watertap_costing_available = attempt_import("watertap.costing")
if watertap_costing_available:
    import watertap.property_models.NaCl_prop_pack as props
    from watertap.core.solvers import get_solver as get_watertap_solver
    from watertap.core.util.initialization import check_dof
    from watertap.core.wt_database import Database
    from watertap.core.zero_order_properties import WaterParameterBlock
    from watertap.property_models.multicomp_aq_sol_prop_pack import (
        ActivityCoefficientModel,
        DensityCalculation,
        MCASParameterBlock,
    )
    from watertap.unit_models.ion_exchange_0D import IonExchange0D
    from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
    from watertap.unit_models.reverse_osmosis_1D import (
        ConcentrationPolarizationType,
        MassTransferCoefficient,
        PressureChangeType,
        ReverseOsmosis1D,
    )
    from watertap.unit_models.zero_order import NanofiltrationZO

_log = idaeslog.getLogger(__name__)


@pytest.mark.unit
def test_register_REE_currency_units_twice(caplog):
    # check that units exist - they are registered when QGESSCosting imports
    assert hasattr(pyunits, "USD_UKy_2019")
    assert hasattr(pyunits, "USD_2022")
    assert hasattr(pyunits, "USD_2025")

    # register units again
    custom_REE_plant_currency_units()
    msg = (
        "Custom REE plant currency units (USD_2022, USD_2025, USD_UKy_2019) "
        "already appear in Pyomo unit registry. Assuming repeated call of "
        "custom_power_plant_currency_units."
    )
    for record in caplog.records:
        assert msg in record.message


def base_model():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting(
        # arguments for NPV
        discount_percentage=10,  # percent
        plant_lifetime=20,  # years
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[10, 60, 30],
        capital_escalation_percentage=3.6,
        capital_loan_interest_percentage=6,
        capital_loan_repayment_period=10,
        debt_percentage_of_CAPEX=50,
        operating_inflation_percentage=3,
        revenue_inflation_percentage=3,
    )
    CE_index_year = "UKy_2019"

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # accounts 1.x are crushing and screening
    # accounts 2.x are dry grinding
    # accounts 3.x are roasting
    # accounts 4.x are leaching
    # accounts 5.x are rougher solvent extraction
    # accounts 6.x are cleaner solvent extraction
    # accounts 7.x are solvent extraction wash & saponification
    # accounts 9.x are REE precipitation
    # accounts 11.x are water treatment

    # 1.1 is Front End Loader (2 cuyd)
    # this is a constant-cost unit, where n_equip is the scaling parameter
    CS_front_end_loader_2yd3_accounts = ["1.1"]
    m.fs.CS_front_end_loader_2yd3 = UnitModelBlock()
    m.fs.CS_front_end_loader_2yd3.n_equip = pyo.Var(
        initialize=1, units=pyunits.dimensionless
    )
    m.fs.CS_front_end_loader_2yd3.n_equip.fix()

    m.fs.CS_front_end_loader_2yd3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_front_end_loader_2yd3_accounts,
            "scaled_param": m.fs.CS_front_end_loader_2yd3.n_equip,  # 1 loader
            "source": 1,
            # no. units is the scaling parameter for constant-cost units,
            #     so use n_equip below to specify the number of loaders
            "n_equip": 5,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 1.5 is CS Roll Crusher
    CS_roll_crusher_accounts = ["1.5"]
    m.fs.CS_roll_crusher = UnitModelBlock()
    m.fs.CS_roll_crusher.power = pyo.Var(initialize=430, units=pyunits.hp)
    m.fs.CS_roll_crusher.power.fix()
    m.fs.CS_roll_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_roll_crusher_accounts,
            "scaled_param": m.fs.CS_roll_crusher.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 1.6 is CS Vibrating Screens
    CS_vibrating_screen_accounts = ["1.6"]
    m.fs.CS_vibrating_screens = UnitModelBlock()
    m.fs.CS_vibrating_screens.area = pyo.Var(initialize=124, units=pyunits.ft**2)
    m.fs.CS_vibrating_screens.area.fix()
    m.fs.CS_vibrating_screens.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_vibrating_screen_accounts,
            "scaled_param": m.fs.CS_vibrating_screens.area,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 1.7 is CS Conveyors
    CS_conveyors_accounts = ["1.7"]
    m.fs.CS_conveyors = UnitModelBlock()
    m.fs.CS_conveyors.throughput = pyo.Var(
        initialize=575, units=pyunits.ton / pyunits.hr
    )
    m.fs.CS_conveyors.throughput.fix()
    m.fs.CS_conveyors.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_conveyors_accounts,
            "scaled_param": m.fs.CS_conveyors.throughput,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 2.1 is DG Vibrating Screens
    DG_vibrating_screen_accounts = ["2.1"]
    m.fs.DG_vibrating_screens = UnitModelBlock()
    m.fs.DG_vibrating_screens.area = pyo.Var(initialize=332, units=pyunits.ft**2)
    m.fs.DG_vibrating_screens.area.fix()
    m.fs.DG_vibrating_screens.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": DG_vibrating_screen_accounts,
            "scaled_param": m.fs.DG_vibrating_screens.area,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 2.2 is DG Storage Bins
    DG_storage_bins_accounts = ["2.2"]
    m.fs.DG_storage_bins = UnitModelBlock()
    m.fs.DG_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.ton)
    m.fs.DG_storage_bins.capacity.fix()
    m.fs.DG_storage_bins.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": DG_storage_bins_accounts,
            "scaled_param": m.fs.DG_storage_bins.capacity,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 2.3 is DG Air Swept Ball Mill
    DG_air_swept_ball_mill_accounts = ["2.3"]
    m.fs.DG_air_swept_ball_mill = UnitModelBlock()
    m.fs.DG_air_swept_ball_mill.power = pyo.Var(initialize=5609, units=pyunits.hp)
    m.fs.DG_air_swept_ball_mill.power.fix()
    m.fs.DG_air_swept_ball_mill.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": DG_air_swept_ball_mill_accounts,
            "scaled_param": m.fs.DG_air_swept_ball_mill.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 2.4 is DG Bucket Elevator
    DG_bucket_elevator_accounts = ["2.4"]
    m.fs.DG_bucket_elevator = UnitModelBlock()
    m.fs.DG_bucket_elevator.n_equip = pyo.Var(initialize=1, units=pyunits.dimensionless)
    m.fs.DG_bucket_elevator.n_equip.fix()
    m.fs.DG_bucket_elevator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": DG_bucket_elevator_accounts,
            "scaled_param": m.fs.DG_bucket_elevator.n_equip,  # 1 elevator
            "source": 1,
            # no. units is the scaling parameter for constant-cost units,
            #     so use n_equip below to specify the number of elevators
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 2.5 is DG Elevator Motor
    DG_elevator_motor_accounts = ["2.5"]
    m.fs.DG_elevator_motor = UnitModelBlock()
    m.fs.DG_elevator_motor.power = pyo.Var(initialize=58.0, units=pyunits.hp)
    m.fs.DG_elevator_motor.power.fix()
    m.fs.DG_elevator_motor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": DG_elevator_motor_accounts,
            "scaled_param": m.fs.DG_elevator_motor.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.1 is R Storage Bins
    R_storage_bins_accounts = ["3.1"]
    m.fs.R_storage_bins = UnitModelBlock()
    m.fs.R_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.ton)
    m.fs.R_storage_bins.capacity.fix()
    m.fs.R_storage_bins.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_storage_bins_accounts,
            "scaled_param": m.fs.R_storage_bins.capacity,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.2 is R Conveyors
    R_conveyors_accounts = ["3.2"]
    m.fs.R_conveyors = UnitModelBlock()
    m.fs.R_conveyors.throughput = pyo.Var(
        initialize=575, units=pyunits.ton / pyunits.hr
    )
    m.fs.R_conveyors.throughput.fix()
    m.fs.R_conveyors.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_conveyors_accounts,
            "scaled_param": m.fs.R_conveyors.throughput,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.3 is R Roaster
    R_roaster_accounts = ["3.3"]
    m.fs.R_roaster = UnitModelBlock()
    m.fs.R_roaster.duty = pyo.Var(initialize=737, units=pyunits.MBTU / pyunits.hr)
    m.fs.R_roaster.duty.fix()
    m.fs.R_roaster.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_roaster_accounts,
            "scaled_param": m.fs.R_roaster.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.4 is R Gas Scrubber
    R_gas_scrubber_accounts = ["3.4"]
    m.fs.R_gas_scrubber = UnitModelBlock()
    m.fs.R_gas_scrubber.gas_rate = pyo.Var(
        initialize=11500, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.R_gas_scrubber.gas_rate.fix()
    m.fs.R_gas_scrubber.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_gas_scrubber_accounts,
            "scaled_param": m.fs.R_gas_scrubber.gas_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.5 is R Spray Chamber Quencher (7-60 kcfm)
    R_spray_chamber_quencher_accounts = ["3.5"]
    m.fs.R_spray_chamber_quencher = UnitModelBlock()
    m.fs.R_spray_chamber_quencher.gas_rate = pyo.Var(
        initialize=11500, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.R_spray_chamber_quencher.gas_rate.fix()
    m.fs.R_spray_chamber_quencher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_spray_chamber_quencher_accounts,
            "scaled_param": m.fs.R_spray_chamber_quencher.gas_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 3.7 is R Chiller
    R_chiller_accounts = ["3.7"]
    m.fs.R_chiller = UnitModelBlock()
    m.fs.R_chiller.duty = pyo.Var(initialize=131, units=pyunits.MBTU / pyunits.hr)
    m.fs.R_chiller.duty.fix()
    m.fs.R_chiller.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": R_chiller_accounts,
            "scaled_param": m.fs.R_chiller.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.2 is L PE Tanks
    L_pe_tanks_accounts = ["4.2"]
    m.fs.L_pe_tanks = UnitModelBlock()
    m.fs.L_pe_tanks.capacity = pyo.Var(initialize=164805, units=pyunits.gal)
    m.fs.L_pe_tanks.capacity.fix()
    m.fs.L_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pe_tanks_accounts,
            "scaled_param": m.fs.L_pe_tanks.capacity,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.3 is L Tank Mixer
    L_tank_mixer_accounts = ["4.3"]
    m.fs.L_tank_mixers = UnitModelBlock()
    m.fs.L_tank_mixers.power = pyo.Var(initialize=474, units=pyunits.hp)
    m.fs.L_tank_mixers.power.fix()
    m.fs.L_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_tank_mixer_accounts,
            "scaled_param": m.fs.L_tank_mixers.power,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.4 is L Process Pump
    L_pump_accounts = ["4.4"]
    m.fs.L_pump = UnitModelBlock()
    m.fs.L_pump.feed_rate = pyo.Var(initialize=10987, units=pyunits.gal / pyunits.min)
    m.fs.L_pump.feed_rate.fix()
    m.fs.L_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_pump_accounts,
            "scaled_param": m.fs.L_pump.feed_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.5 is L Thickener
    L_thickener_accounts = ["4.5"]
    m.fs.L_thickener = UnitModelBlock()
    m.fs.L_thickener.area = pyo.Var(initialize=22590, units=pyunits.ft**2)
    m.fs.L_thickener.area.fix()
    m.fs.L_thickener.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_thickener_accounts,
            "scaled_param": m.fs.L_thickener.area,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.6 is L Solid Waste Filter Press
    L_filter_press_accounts = ["4.6"]
    m.fs.L_filter_press = UnitModelBlock()
    m.fs.L_filter_press.volume = pyo.Var(initialize=3600, units=pyunits.ft**3)
    m.fs.L_filter_press.volume.fix()
    m.fs.L_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_filter_press_accounts,
            "scaled_param": m.fs.L_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 4.8 is L Solution Heater
    L_solution_heater_accounts = ["4.8"]
    m.fs.L_solution_heater = UnitModelBlock()
    m.fs.L_solution_heater.duty = pyo.Var(
        initialize=2.4, units=pyunits.MBTU / pyunits.hr
    )
    m.fs.L_solution_heater.duty.fix()
    m.fs.L_solution_heater.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": L_solution_heater_accounts,
            "scaled_param": m.fs.L_solution_heater.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.1 is RSX PE Tanks
    RSX_pe_tanks_accounts = ["5.1"]
    m.fs.RSX_pe_tanks = UnitModelBlock()
    m.fs.RSX_pe_tanks.capacity = pyo.Var(initialize=35136, units=pyunits.gal)
    m.fs.RSX_pe_tanks.capacity.fix()
    m.fs.RSX_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pe_tanks_accounts,
            "scaled_param": m.fs.RSX_pe_tanks.capacity,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.2 is RSX Tank Mixer
    RSX_tank_mixer_accounts = ["5.2"]
    m.fs.RSX_tank_mixers = UnitModelBlock()
    m.fs.RSX_tank_mixers.power = pyo.Var(initialize=20, units=pyunits.hp)
    m.fs.RSX_tank_mixers.power.fix()
    m.fs.RSX_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_tank_mixer_accounts,
            "scaled_param": m.fs.RSX_tank_mixers.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.3 is RSX Process Pump
    RSX_pump_accounts = ["5.3"]
    m.fs.RSX_pump = UnitModelBlock()
    m.fs.RSX_pump.feed_rate = pyo.Var(initialize=7027, units=pyunits.gal / pyunits.min)
    m.fs.RSX_pump.feed_rate.fix()
    m.fs.RSX_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pump_accounts,
            "scaled_param": m.fs.RSX_pump.feed_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 5.4 is RSX Mixer Settler
    RSX_mixer_settler_accounts = ["5.4"]
    m.fs.RSX_mixer_settler = UnitModelBlock()
    m.fs.RSX_mixer_settler.volume = pyo.Var(initialize=61107, units=pyunits.gal)
    m.fs.RSX_mixer_settler.volume.fix()
    m.fs.RSX_mixer_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": RSX_mixer_settler_accounts,
            "scaled_param": m.fs.RSX_mixer_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.1 is CSX PE Tanks
    CSX_pe_tanks_accounts = ["6.1"]
    m.fs.CSX_pe_tanks = UnitModelBlock()
    m.fs.CSX_pe_tanks.capacity = pyo.Var(initialize=1405, units=pyunits.gal)
    m.fs.CSX_pe_tanks.capacity.fix()
    m.fs.CSX_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pe_tanks_accounts,
            "scaled_param": m.fs.CSX_pe_tanks.capacity,
            "source": 1,
            "n_equip": 5,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.2 is CSX Tank Mixer
    CSX_tank_mixer_accounts = ["6.2"]
    m.fs.CSX_tank_mixers = UnitModelBlock()
    m.fs.CSX_tank_mixers.power = pyo.Var(initialize=0.8, units=pyunits.hp)
    m.fs.CSX_tank_mixers.power.fix()
    m.fs.CSX_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_tank_mixer_accounts,
            "scaled_param": m.fs.CSX_tank_mixers.power,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.3 is CSX Process Pump
    CSX_pump_accounts = ["6.3"]
    m.fs.CSX_pump = UnitModelBlock()
    m.fs.CSX_pump.feed_rate = pyo.Var(initialize=281, units=pyunits.gal / pyunits.min)
    m.fs.CSX_pump.feed_rate.fix()
    m.fs.CSX_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pump_accounts,
            "scaled_param": m.fs.CSX_pump.feed_rate,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 6.4 is CSX Mixer Settler
    CSX_mixer_settler_accounts = ["6.4"]
    m.fs.CSX_mixer_settler = UnitModelBlock()
    m.fs.CSX_mixer_settler.volume = pyo.Var(initialize=2444, units=pyunits.gal)
    m.fs.CSX_mixer_settler.volume.fix()
    m.fs.CSX_mixer_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CSX_mixer_settler_accounts,
            "scaled_param": m.fs.CSX_mixer_settler.volume,
            "source": 1,
            "n_equip": 6,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 7.1 is SX Wash PE Tanks
    SX_wash_pe_tanks_accounts = ["7.1"]
    m.fs.SX_wash_pe_tanks = UnitModelBlock()
    m.fs.SX_wash_pe_tanks.capacity = pyo.Var(initialize=3514, units=pyunits.gal)
    m.fs.SX_wash_pe_tanks.capacity.fix()
    m.fs.SX_wash_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_pe_tanks_accounts,
            "scaled_param": m.fs.SX_wash_pe_tanks.capacity,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 7.2 is SX Wash Tank Mixer
    SX_wash_tank_mixer_accounts = ["7.2"]
    m.fs.SX_wash_tank_mixers = UnitModelBlock()
    m.fs.SX_wash_tank_mixers.power = pyo.Var(initialize=2, units=pyunits.hp)
    m.fs.SX_wash_tank_mixers.power.fix()
    m.fs.SX_wash_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_tank_mixer_accounts,
            "scaled_param": m.fs.SX_wash_tank_mixers.power,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 7.3 is SX Wash Process Pump
    SX_wash_pump_accounts = ["7.3"]
    m.fs.SX_wash_pump = UnitModelBlock()
    m.fs.SX_wash_pump.feed_rate = pyo.Var(
        initialize=703, units=pyunits.gal / pyunits.min
    )
    m.fs.SX_wash_pump.feed_rate.fix()
    m.fs.SX_wash_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_pump_accounts,
            "scaled_param": m.fs.SX_wash_pump.feed_rate,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 7.4 is SX Wash Mixer Settler
    SX_wash_mixer_settler_accounts = ["7.4"]
    m.fs.SX_wash_mixer_settler = UnitModelBlock()
    m.fs.SX_wash_mixer_settler.volume = pyo.Var(initialize=18332, units=pyunits.gal)
    m.fs.SX_wash_mixer_settler.volume.fix()
    m.fs.SX_wash_mixer_settler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_mixer_settler_accounts,
            "scaled_param": m.fs.SX_wash_mixer_settler.volume,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 7.5 is SX Wash Filter Press
    SX_wash_filter_press_accounts = ["7.5"]
    m.fs.SX_wash_filter_press = UnitModelBlock()
    m.fs.SX_wash_filter_press.volume = pyo.Var(initialize=0.26, units=pyunits.ft**3)
    m.fs.SX_wash_filter_press.volume.fix()
    m.fs.SX_wash_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_filter_press_accounts,
            "scaled_param": m.fs.SX_wash_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.2 is REE Precipitation PE Tanks
    reep_pe_tanks_accounts = ["9.2"]
    m.fs.reep_pe_tanks = UnitModelBlock()
    m.fs.reep_pe_tanks.capacity = pyo.Var(initialize=1504, units=pyunits.gal)
    m.fs.reep_pe_tanks.capacity.fix()
    m.fs.reep_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pe_tanks_accounts,
            "scaled_param": m.fs.reep_pe_tanks.capacity,
            "source": 1,
            "n_equip": 4,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.3 is REE Precipitation Tank Mixer
    reep_tank_mixer_accounts = ["9.3"]
    m.fs.reep_tank_mixers = UnitModelBlock()
    m.fs.reep_tank_mixers.power = pyo.Var(initialize=0.61, units=pyunits.hp)
    m.fs.reep_tank_mixers.power.fix()
    m.fs.reep_tank_mixers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_tank_mixer_accounts,
            "scaled_param": m.fs.reep_tank_mixers.power,
            "source": 1,
            "n_equip": 3,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.4 is REE Precipitation Process Pump
    reep_pump_accounts = ["9.4"]
    m.fs.reep_pump = UnitModelBlock()
    m.fs.reep_pump.feed_rate = pyo.Var(initialize=70, units=pyunits.gal / pyunits.min)
    m.fs.reep_pump.feed_rate.fix()
    m.fs.reep_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_pump_accounts,
            "scaled_param": m.fs.reep_pump.feed_rate,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.5 is REE Precipitation Filter Press
    reep_filter_press_accounts = ["9.5"]
    m.fs.reep_filter_press = UnitModelBlock()
    m.fs.reep_filter_press.volume = pyo.Var(initialize=0.405, units=pyunits.ft**3)
    m.fs.reep_filter_press.volume.fix()
    m.fs.reep_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_filter_press_accounts,
            "scaled_param": m.fs.reep_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 9.8 is REE Precipitation Roaster
    reep_roaster_accounts = ["9.8"]
    m.fs.reep_roaster = UnitModelBlock()
    m.fs.reep_roaster.duty = pyo.Var(initialize=0.35, units=pyunits.MBTU / pyunits.hr)
    m.fs.reep_roaster.duty.fix()
    m.fs.reep_roaster.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": reep_roaster_accounts,
            "scaled_param": m.fs.reep_roaster.duty,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 11.1 is Water Treatment PE Tanks
    WT_pe_tanks_accounts = ["11.1"]
    m.fs.WT_pe_tanks = UnitModelBlock()
    m.fs.WT_pe_tanks.capacity = pyo.Var(initialize=453131, units=pyunits.gal)
    m.fs.WT_pe_tanks.capacity.fix()
    m.fs.WT_pe_tanks.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": WT_pe_tanks_accounts,
            "scaled_param": m.fs.WT_pe_tanks.capacity,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 11.2 is Water Treatment Process Pump
    WT_pump_accounts = ["11.2"]
    m.fs.WT_pump = UnitModelBlock()
    m.fs.WT_pump.feed_rate = pyo.Var(initialize=78805, units=pyunits.gal / pyunits.min)
    m.fs.WT_pump.feed_rate.fix()
    m.fs.WT_pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": WT_pump_accounts,
            "scaled_param": m.fs.WT_pump.feed_rate,
            "source": 1,
            "n_equip": 2,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 11.3 is Water Treatment Filter Press
    WT_filter_press_accounts = ["11.3"]
    m.fs.WT_filter_press = UnitModelBlock()
    m.fs.WT_filter_press.volume = pyo.Var(initialize=469, units=pyunits.ft**3)
    m.fs.WT_filter_press.volume.fix()
    m.fs.WT_filter_press.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": WT_filter_press_accounts,
            "scaled_param": m.fs.WT_filter_press.volume,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    # 11.4 is Water Treatment Conveyors
    WT_conveyors_accounts = ["11.4"]
    m.fs.WT_conveyors = UnitModelBlock()
    m.fs.WT_conveyors.throughput = pyo.Var(
        initialize=569, units=pyunits.ton / pyunits.hr
    )
    m.fs.WT_conveyors.throughput.fix()
    m.fs.WT_conveyors.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": WT_conveyors_accounts,
            "scaled_param": m.fs.WT_conveyors.throughput,
            "source": 1,
            "n_equip": 1,
            "scale_down_parallel_equip": False,
            "CE_index_year": CE_index_year,
        },
    )

    return m


class TestREECosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = base_model()
        return model

    @pytest.mark.unit
    def test_base_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.unit
    def test_REE_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included
        CE_index_year = "UKy_2019"

        # add plant-level cost constraints

        model.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)

        hours_per_shift = 8
        shifts_per_day = 3
        operating_days_per_year = 336

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
            mutable=False,
            units=pyunits.hours / pyunits.year,
        )

        model.fs.recovery_rate_per_year = pyo.Var(
            initialize=39.3
            * pyunits.kg
            / pyunits.hr
            * 0.8025  # TREO (total rare earth oxide), 80.25% REE in REO
            * model.fs.annual_operating_hours,
            units=pyunits.kg / pyunits.yr,
        )

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * getattr(pyunits, "MUSD_" + CE_index_year)
            / pyunits.ton
            * pyunits.convert(model.fs.feed_input, to_units=pyunits.ton / pyunits.hr)
            * hours_per_shift
            * pyunits.hr
            * shifts_per_day
            * pyunits.day**-1
            * operating_days_per_year
            * pyunits.day
        )

        # dummy reagent with cost of 1 USD/kg for each section
        reagent_costs = (
            (  # all USD/year
                302962  # Crushing and Screening
                + 0  # Dry Grinding
                + 5767543  # Roasting
                + 199053595  # Leaching
                + 152303329  # Rougher Solvent Extraction
                + 43702016  # Cleaner Solvent Extraction
                + 7207168  # Solvent Extraction Wash and Saponification
                + 1233763  # Rare Earth Element Precipiation
                + 18684816  # Water Treatment
            )
            * pyunits.kg
            / pyunits.year
        )

        model.fs.reagents = pyo.Var(
            model.fs.time,
            initialize=reagent_costs / (model.fs.annual_operating_hours),
            units=pyunits.kg / pyunits.hr,
        )

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)

        resources = [
            "dummy",
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
            "power",
        ]

        rates = [
            model.fs.reagents,
            model.fs.solid_waste,
            model.fs.precipitate,
            model.fs.dust_and_volatiles,
            model.fs.power,
        ]

        # define product flowrates

        pure_product_output_rates = {
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.4 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.5 * pyunits.kg / pyunits.hr,
        }

        mixed_product_output_rates = {
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            "Y2O3": 0.05418 * pyunits.kg / pyunits.hr,
            "La2O3": 0.13770 * pyunits.kg / pyunits.hr,
            "CeO2": 0.37383 * pyunits.kg / pyunits.hr,
            "Pr6O11": 0.03941 * pyunits.kg / pyunits.hr,
            "Nd2O3": 0.17289 * pyunits.kg / pyunits.hr,
            "Sm2O3": 0.02358 * pyunits.kg / pyunits.hr,
            "Eu2O3": 0.00199 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tb4O7": 0.00801 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Ho2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Er2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tm2O3": 0.00130 * pyunits.kg / pyunits.hr,
            "Yb2O3": 0.00373 * pyunits.kg / pyunits.hr,
            "Lu2O3": 0.00105 * pyunits.kg / pyunits.hr,
        }

        model.fs.costing.build_process_costs(
            # arguments related to installation costs
            piping_materials_and_labor_percentage=20,
            electrical_materials_and_labor_percentage=20,
            instrumentation_percentage=8,
            plants_services_percentage=10,
            process_buildings_percentage=40,
            auxiliary_buildings_percentage=15,
            site_improvements_percentage=10,
            equipment_installation_percentage=17,
            field_expenses_percentage=12,
            project_management_and_construction_percentage=30,
            process_contingency_percentage=15,
            # argument related to Fixed OM costs
            labor_types=[
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
            ],
            labor_rate=[24.98, 19.08, 30.39, 22.73, 21.97, 45.85],  # USD/hr
            labor_burden=25,  # % fringe benefits
            operators_per_shift=[4, 9, 2, 2, 2, 3],
            hours_per_shift=hours_per_shift,
            shifts_per_day=shifts_per_day,
            operating_days_per_year=operating_days_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            mixed_product_sale_price_realization_factor=0.65,  # 65% price realization for mixed products
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=resources,
            rates=rates,
            prices={
                "dummy": 1 * getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            },
            fixed_OM=True,
            variable_OM=True,
            feed_input=model.fs.feed_input,
            efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
            chemicals=["dummy"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
            ],
            recovery_rate_per_year=model.fs.recovery_rate_per_year,
            CE_index_year=CE_index_year,
            transport_cost_per_ton_product=0,
            # arguments related to NPV calculation
            calculate_NPV=True,
        )

        # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
        model.fs.costing.other_plant_costs.unfix()
        model.fs.costing.other_plant_costs_eq = pyo.Constraint(
            expr=(
                model.fs.costing.other_plant_costs
                == pyunits.convert(
                    1218073 * pyunits.USD_2016  # Rougher Solvent Extraction
                    + 48723 * pyunits.USD_2016  # Cleaner Solvent Extraction
                    + 182711
                    * pyunits.USD_2016,  # Solvent Extraction Wash and Saponification
                    to_units=getattr(pyunits, "MUSD_" + CE_index_year),
                )
            )
        )

        # fix costing vars that shouldn't change
        model.fs.feed_input.fix()
        model.fs.feed_grade.fix()
        model.fs.recovery_rate_per_year.fix()
        model.fs.reagents.fix()
        model.fs.solid_waste.fix()
        model.fs.precipitate.fix()
        model.fs.dust_and_volatiles.fix()
        model.fs.power.fix()

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_full_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_initialize(self, model):
        # add initialize
        QGESSCostingData.costing_initialization(model.fs.costing)
        QGESSCostingData.initialize_fixed_OM_costs(model.fs.costing)
        QGESSCostingData.initialize_variable_OM_costs(model.fs.costing)

    @pytest.mark.component
    def test_solve(self, model):
        # try solving
        solver = get_solver()
        results = solver.solve(model, tee=True)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solved_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_results(self, model):
        # check some overall cost results

        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            133.23, rel=1e-4
        )
        assert value(model.fs.costing.total_BEC) == pytest.approx(44.308, rel=1e-4)
        assert value(model.fs.costing.total_installation_cost) == pytest.approx(
            87.287, rel=1e-4
        )
        assert value(model.fs.costing.other_plant_costs) == pytest.approx(
            1.6309, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            10.916, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            525.71, rel=1e-4
        )
        assert value(model.fs.costing.land_cost) == pytest.approx(
            1.2247, rel=1e-4
        )  # Expression, not Var
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(
            27.654, rel=1e-4
        )
        assert value(model.fs.costing.pv_capital_cost) == pytest.approx(
            -112.78144, rel=1e-4
        )
        assert value(model.fs.costing.pv_loan_interest) == pytest.approx(
            -11.001142, rel=1e-4
        )
        assert value(model.fs.costing.pv_operating_cost) == pytest.approx(
            -4614.5826, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(237.25943, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(-4501.10576, rel=1e-4)

    @pytest.mark.unit
    def test_report(self, model):
        # test report methods
        QGESSCostingData.report(model.fs.costing, export=True)
        assert os.path.exists(os.path.join(os.getcwd(), "costing_report.csv"))
        # cleanup
        os.remove(os.path.join(os.getcwd(), "costing_report.csv"))
        assert not os.path.exists(os.path.join(os.getcwd(), "costing_report.csv"))

        model.fs.costing.variable_operating_costs.display()  # results will be in t = 0
        QGESSCostingData.display_total_plant_costs(model.fs.costing)
        QGESSCostingData.display_bare_erected_costs(model.fs.costing)
        QGESSCostingData.display_flowsheet_cost(model.fs.costing)

    @pytest.mark.unit
    def test_costing_bounding_build_diagnostics(self, model):
        # test costing bounding method
        CE_index_year = "UKy_2019"
        QGESSCostingData.calculate_REE_costing_bounds(
            b=model.fs.costing,
            capacity=model.fs.feed_input
            * model.fs.annual_operating_hours
            * 20
            * pyunits.year,
            grade=model.fs.feed_grade,
            CE_index_year=CE_index_year,
        )

        dt = DiagnosticsToolbox(model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_costing_bounding_solve(self, model):

        # solve new variables and constraints
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing_bounding_solve_diagnostics(self, model):

        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_costing_bounding_results(self, model):

        expected_costing_lower_bound = {
            "Beneficiation": 0.12109,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 5.3203,
            "Chemical Extraction": 0.062054,
            "Chemical Extraction, Enrichment and Separation": 0.0000,
            "Enrichment and Separation": 0.18702,
            "Mining": 0.24696,
            "Total Capital": 0.29132,
            "Total Operating": 5.4747,
        }

        expected_costing_upper_bound = {
            "Beneficiation": 1.2209,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 15.235,
            "Chemical Extraction": 1.2373,
            "Chemical Extraction, Enrichment and Separation": 55.9608,
            "Enrichment and Separation": 4.2676,
            "Mining": 2.0840,
            "Total Capital": 1.0861,
            "Total Operating": 12.648,
        }

        for key in model.fs.costing.costing_lower_bound.keys():

            assert value(model.fs.costing.costing_lower_bound[key]) == pytest.approx(
                expected_costing_lower_bound[key], abs=1e-8, rel=1e-4
            )

        for key in model.fs.costing.costing_upper_bound.keys():
            assert value(model.fs.costing.costing_upper_bound[key]) == pytest.approx(
                expected_costing_upper_bound[key], rel=1e-4
            )

    @pytest.mark.component
    def test_costing_bounding_rerun_norecalculate(self, model):

        # call again, should give same results
        CE_index_year = "UKy_2019"

        QGESSCostingData.calculate_REE_costing_bounds(
            b=model.fs.costing,
            capacity=model.fs.feed_input
            * model.fs.annual_operating_hours
            * 20
            * pyunits.year,
            grade=model.fs.feed_grade,
            CE_index_year=CE_index_year,
            recalculate=False,
        )

        expected_costing_lower_bound = {
            "Beneficiation": 0.12109,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 5.3203,
            "Chemical Extraction": 0.062054,
            "Chemical Extraction, Enrichment and Separation": 0.0000,
            "Enrichment and Separation": 0.18702,
            "Mining": 0.24696,
            "Total Capital": 0.29132,
            "Total Operating": 5.4747,
        }

        expected_costing_upper_bound = {
            "Beneficiation": 1.2209,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 15.235,
            "Chemical Extraction": 1.2373,
            "Chemical Extraction, Enrichment and Separation": 55.9608,
            "Enrichment and Separation": 4.2676,
            "Mining": 2.0840,
            "Total Capital": 1.0861,
            "Total Operating": 12.648,
        }

        for key in model.fs.costing.costing_lower_bound.keys():

            assert value(model.fs.costing.costing_lower_bound[key]) == pytest.approx(
                expected_costing_lower_bound[key], abs=1e-8, rel=1e-4
            )

        for key in model.fs.costing.costing_upper_bound.keys():
            assert value(model.fs.costing.costing_upper_bound[key]) == pytest.approx(
                expected_costing_upper_bound[key], rel=1e-4
            )

    @pytest.mark.component
    def test_costing_bounding_rerun_recalculate(self, model):

        # call again, should give same results
        CE_index_year = "UKy_2019"

        QGESSCostingData.calculate_REE_costing_bounds(
            b=model.fs.costing,
            capacity=model.fs.feed_input
            * model.fs.annual_operating_hours
            * 20
            * pyunits.year,
            grade=model.fs.feed_grade,
            CE_index_year=CE_index_year,
            recalculate=True,
        )

        expected_costing_lower_bound = {
            "Beneficiation": 0.12109,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 5.3203,
            "Chemical Extraction": 0.062054,
            "Chemical Extraction, Enrichment and Separation": 0.0000,
            "Enrichment and Separation": 0.18702,
            "Mining": 0.24696,
            "Total Capital": 0.29132,
            "Total Operating": 5.4747,
        }

        expected_costing_upper_bound = {
            "Beneficiation": 1.2209,
            "Beneficiation, Chemical Extraction, Enrichment and Separation": 15.235,
            "Chemical Extraction": 1.2373,
            "Chemical Extraction, Enrichment and Separation": 55.9608,
            "Enrichment and Separation": 4.2676,
            "Mining": 2.0840,
            "Total Capital": 1.0861,
            "Total Operating": 12.648,
        }

        for key in model.fs.costing.costing_lower_bound.keys():

            assert value(model.fs.costing.costing_lower_bound[key]) == pytest.approx(
                expected_costing_lower_bound[key], abs=1e-8, rel=1e-4
            )

        for key in model.fs.costing.costing_upper_bound.keys():
            assert value(model.fs.costing.costing_upper_bound[key]) == pytest.approx(
                expected_costing_upper_bound[key], rel=1e-4
            )


class TestWaterTAPCosting(object):
    @pytest.fixture(scope="class")
    def solver(self):
        pytest.importorskip("watertap", reason="WaterTAP dependency not available")
        return get_watertap_solver()

    @pytest.fixture(scope="class")
    def model(self, solver):
        model = base_model()
        model.fs_membrane = FlowsheetBlock(dynamic=False)
        # Nanofiltration

        model.fs_membrane.nf_properties = MCASParameterBlock(
            solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
            diffusivity_data={
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "SO4_2-"): 1.06e-09,
                ("Liq", "Mg_2+"): 7.06e-10,
                ("Liq", "Na_+"): 1.33e-09,
                ("Liq", "Cl_-"): 2.03e-09,
            },
            mw_data={
                "H2O": 0.018,
                "Ca_2+": 0.04,
                "Mg_2+": 0.024,
                "SO4_2-": 0.096,
                "Na_+": 0.023,
                "Cl_-": 0.035,
            },
            stokes_radius_data={
                "Ca_2+": 3.09e-10,
                "Mg_2+": 3.47e-10,
                "SO4_2-": 2.3e-10,
                "Cl_-": 1.21e-10,
                "Na_+": 1.84e-10,
            },
            charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
            activity_coefficient_model=ActivityCoefficientModel.davies,
            density_calculation=DensityCalculation.constant,
        )

        model.fs_membrane.nfunit = NanofiltrationDSPMDE0D(
            property_package=model.fs_membrane.nf_properties
        )
        mass_flow_in = 1 * pyunits.kg / pyunits.s
        feed_mass_frac = {
            "Ca_2+": 382e-6,
            "Mg_2+": 1394e-6,
            "SO4_2-": 2136e-6,
            "Cl_-": 20101.6e-6,
            "Na_+": 11122e-6,
        }

        # Fix mole flow rates of each ion and water
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / model.fs_membrane.nfunit.feed_side.properties_in[0].mw_comp[ion]
            )
            model.fs_membrane.nfunit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(
                mol_comp_flow
            )

        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / model.fs_membrane.nfunit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        model.fs_membrane.nfunit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
            H2O_mol_comp_flow
        )

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        model.fs_membrane.nfunit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        model.fs_membrane.nfunit.inlet.temperature[0].fix(298.15)
        model.fs_membrane.nfunit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        model.fs_membrane.nfunit.radius_pore.fix(0.5e-9)
        model.fs_membrane.nfunit.membrane_thickness_effective.fix(1.33e-6)
        model.fs_membrane.nfunit.membrane_charge_density.fix(-27)
        model.fs_membrane.nfunit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        model.fs_membrane.nfunit.mixed_permeate[0].pressure.fix(101325)

        model.fs_membrane.nfunit.spacer_porosity.fix(0.85)
        model.fs_membrane.nfunit.channel_height.fix(5e-4)
        model.fs_membrane.nfunit.velocity[0, 0].fix(0.25)
        model.fs_membrane.nfunit.area.fix(50)
        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        model.fs_membrane.nfunit.spacer_mixing_efficiency.fix()
        model.fs_membrane.nfunit.spacer_mixing_length.fix()

        check_dof(model.fs_membrane, fail_flag=True)

        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Ca_2+")
        )
        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "SO4_2-")
        )
        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Mg_2+")
        )
        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        model.fs_membrane.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(model.fs_membrane)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(model.fs_membrane.nfunit))
        assert len(unscaled_var_list) == 0

        model.fs_membrane.nfunit.initialize(optarg=solver.options)

        badly_scaled_var_lst = list(
            badly_scaled_var_generator(model.fs_membrane.nfunit, small=1e-5, zero=1e-12)
        )
        for var, val in badly_scaled_var_lst:
            print(var.name, val)
        assert len(badly_scaled_var_lst) == 0

        results = solver.solve(model.fs_membrane, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Reverse Osmosis

        model.fs_membrane.ro_properties = props.NaClParameterBlock()

        model.fs_membrane.rounit = ReverseOsmosis1D(
            property_package=model.fs_membrane.ro_properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
            has_full_reporting=True,
        )

        # fully specify system
        feed_flow_mass = 1000 / 3600
        feed_mass_frac_NaCl = 0.034283
        feed_pressure = 70e5

        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 1e5
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        model.fs_membrane.rounit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        model.fs_membrane.rounit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        model.fs_membrane.rounit.inlet.pressure[0].fix(feed_pressure)
        model.fs_membrane.rounit.inlet.temperature[0].fix(feed_temperature)
        model.fs_membrane.rounit.A_comp.fix(A)
        model.fs_membrane.rounit.B_comp.fix(B)
        model.fs_membrane.rounit.permeate.pressure[0].fix(pressure_atmospheric)
        model.fs_membrane.rounit.feed_side.N_Re[0, 0].fix(400)
        model.fs_membrane.rounit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
        model.fs_membrane.rounit.feed_side.spacer_porosity.fix(0.97)
        model.fs_membrane.rounit.feed_side.channel_height.fix(0.001)

        check_dof(model.fs_membrane, fail_flag=True)

        model.fs_membrane.ro_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        model.fs_membrane.ro_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(model.fs_membrane)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(model.fs_membrane.rounit))
        assert len(unscaled_var_list) == 0

        model.fs_membrane.rounit.initialize(optarg=solver.options)

        badly_scaled_var_lst = list(
            badly_scaled_var_generator(model.fs_membrane.rounit)
        )
        assert badly_scaled_var_lst == []

        results = solver.solve(model.fs_membrane, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Ion Exchange

        target_ion = "Ca_2+"
        ion_props = {
            "solute_list": [target_ion],
            "diffusivity_data": {("Liq", target_ion): 9.2e-10},
            "mw_data": {"H2O": 0.018, target_ion: 0.04},
            "charge": {target_ion: 2},
        }

        model.fs_membrane.ro_properties = MCASParameterBlock(**ion_props)

        ix_config = {
            "property_package": model.fs_membrane.ro_properties,
            "target_ion": target_ion,
        }
        model.fs_membrane.ixunit = IonExchange0D(**ix_config)
        model.fs_membrane.ixunit.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 0.1,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        model.fs_membrane.ixunit.process_flow.properties_in[0].flow_mass_phase_comp[...]
        model.fs_membrane.ixunit.process_flow.properties_out[0].flow_mass_phase_comp[
            ...
        ]
        model.fs_membrane.ixunit.regeneration_stream[0].flow_mass_phase_comp[...]

        model.fs_membrane.ixunit.service_flow_rate.fix(15)
        model.fs_membrane.ixunit.langmuir[target_ion].fix(0.9)
        model.fs_membrane.ixunit.resin_max_capacity.fix(3)
        model.fs_membrane.ixunit.bed_depth.fix(1.7)
        model.fs_membrane.ixunit.dimensionless_time.fix()
        model.fs_membrane.ixunit.number_columns.fix(8)
        model.fs_membrane.ixunit.resin_diam.fix()
        model.fs_membrane.ixunit.resin_bulk_dens.fix()
        model.fs_membrane.ixunit.bed_porosity.fix()

        check_dof(model.fs_membrane, fail_flag=True)

        model.fs_membrane.ro_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        model.fs_membrane.ro_properties.set_default_scaling(
            "flow_mol_phase_comp", 10, index=("Liq", "Ca_2+")
        )

        calculate_scaling_factors(model.fs_membrane)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(model.fs_membrane.ixunit))
        assert len(unscaled_var_list) == 0

        model.fs_membrane.ixunit.initialize(optarg=solver.options)

        badly_scaled_var_lst = list(
            badly_scaled_var_generator(model.fs_membrane.ixunit)
        )
        assert badly_scaled_var_lst == []

        results = solver.solve(model.fs_membrane, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Nanofiltration Zero Order

        model.fs_membrane.db = Database()
        model.fs_membrane.nfzo_params = WaterParameterBlock(solute_list=["tds", "dye"])

        model.fs_membrane.nfzounit = NanofiltrationZO(
            property_package=model.fs_membrane.nfzo_params,
            database=model.fs_membrane.db,
            process_subtype="rHGO_dye_rejection",
        )

        model.fs_membrane.nfzounit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        model.fs_membrane.nfzounit.inlet.flow_mass_comp[0, "tds"].fix(1)
        model.fs_membrane.nfzounit.inlet.flow_mass_comp[0, "dye"].fix(2)

        model.fs_membrane.nfzounit.load_parameters_from_database(
            use_default_removal=True
        )

        results = solver.solve(model.fs_membrane, tee=True)

        # # Check for optimal solution
        assert_optimal_termination(results)

        return model

    @pytest.mark.component
    def test_REE_watertap_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included
        CE_index_year = "UKy_2019"

        # add plant-level cost constraints

        model.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)

        hours_per_shift = 8
        shifts_per_day = 3
        operating_days_per_year = 336

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
            mutable=False,
            units=pyunits.hours / pyunits.year,
        )

        model.fs.recovery_rate_per_year = pyo.Var(
            initialize=39.3
            * pyunits.kg
            / pyunits.hr
            * 0.8025  # TREO (total rare earth oxide), 80.25% REE in REO
            * model.fs.annual_operating_hours,
            units=pyunits.kg / pyunits.yr,
        )

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * getattr(pyunits, "MUSD_" + CE_index_year)
            / pyunits.ton
            * pyunits.convert(model.fs.feed_input, to_units=pyunits.ton / pyunits.hr)
            * hours_per_shift
            * pyunits.hr
            * shifts_per_day
            * pyunits.day**-1
            * operating_days_per_year
            * pyunits.day
        )

        # dummy reagent with cost of 1 USD/kg for each section
        reagent_costs = (
            (  # all USD/year
                302962  # Crushing and Screening
                + 0  # Dry Grinding
                + 5767543  # Roasting
                + 199053595  # Leaching
                + 152303329  # Rougher Solvent Extraction
                + 43702016  # Cleaner Solvent Extraction
                + 7207168  # Solvent Extraction Wash and Saponification
                + 1233763  # Rare Earth Element Precipiation
                + 18684816  # Water Treatment
            )
            * pyunits.kg
            / pyunits.year
        )

        model.fs.reagents = pyo.Var(
            model.fs.time,
            initialize=reagent_costs / (model.fs.annual_operating_hours),
            units=pyunits.kg / pyunits.hr,
        )

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)

        resources = [
            "dummy",
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
            "power",
        ]

        rates = [
            model.fs.reagents,
            model.fs.solid_waste,
            model.fs.precipitate,
            model.fs.dust_and_volatiles,
            model.fs.power,
        ]

        # define product flowrates

        pure_product_output_rates = {
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.4 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.5 * pyunits.kg / pyunits.hr,
        }

        mixed_product_output_rates = {
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            "Y2O3": 0.05418 * pyunits.kg / pyunits.hr,
            "La2O3": 0.13770 * pyunits.kg / pyunits.hr,
            "CeO2": 0.37383 * pyunits.kg / pyunits.hr,
            "Pr6O11": 0.03941 * pyunits.kg / pyunits.hr,
            "Nd2O3": 0.17289 * pyunits.kg / pyunits.hr,
            "Sm2O3": 0.02358 * pyunits.kg / pyunits.hr,
            "Eu2O3": 0.00199 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tb4O7": 0.00801 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Ho2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Er2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tm2O3": 0.00130 * pyunits.kg / pyunits.hr,
            "Yb2O3": 0.00373 * pyunits.kg / pyunits.hr,
            "Lu2O3": 0.00105 * pyunits.kg / pyunits.hr,
        }

        model.fs.costing.build_process_costs(
            # arguments related to installation costs
            piping_materials_and_labor_percentage=20,
            electrical_materials_and_labor_percentage=20,
            instrumentation_percentage=8,
            plants_services_percentage=10,
            process_buildings_percentage=40,
            auxiliary_buildings_percentage=15,
            site_improvements_percentage=10,
            equipment_installation_percentage=17,
            field_expenses_percentage=12,
            project_management_and_construction_percentage=30,
            process_contingency_percentage=15,
            # argument related to Fixed OM costs
            labor_types=[
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
            ],
            labor_rate=[24.98, 19.08, 30.39, 22.73, 21.97, 45.85],  # USD/hr
            labor_burden=25,  # % fringe benefits
            operators_per_shift=[4, 9, 2, 2, 2, 3],
            hours_per_shift=hours_per_shift,
            shifts_per_day=shifts_per_day,
            operating_days_per_year=operating_days_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            mixed_product_sale_price_realization_factor=0.65,  # 65% price realization for mixed products
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=resources,
            rates=rates,
            prices={
                "dummy": 1 * getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            },
            fixed_OM=True,
            variable_OM=True,
            feed_input=model.fs.feed_input,
            efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
            chemicals=["dummy"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
            ],
            recovery_rate_per_year=model.fs.recovery_rate_per_year,
            CE_index_year=CE_index_year,
            watertap_blocks=[
                model.fs_membrane.nfunit,
                model.fs_membrane.rounit,
                model.fs_membrane.ixunit,
                model.fs_membrane.nfzounit,
            ],
        )

        # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
        model.fs.costing.other_plant_costs.unfix()
        model.fs.costing.other_plant_costs_eq = pyo.Constraint(
            expr=(
                model.fs.costing.other_plant_costs
                == pyunits.convert(
                    1218073 * pyunits.USD_2016  # Rougher Solvent Extraction
                    + 48723 * pyunits.USD_2016  # Cleaner Solvent Extraction
                    + 182711
                    * pyunits.USD_2016,  # Solvent Extraction Wash and Saponification
                    to_units=getattr(pyunits, "MUSD_" + CE_index_year),
                )
            )
        )

        # fix costing vars that shouldn't change
        model.fs.feed_input.fix()
        model.fs.feed_grade.fix()
        model.fs.recovery_rate_per_year.fix()
        model.fs.reagents.fix()
        model.fs.solid_waste.fix()
        model.fs.precipitate.fix()
        model.fs.dust_and_volatiles.fix()
        model.fs.power.fix()

    @pytest.mark.component
    def test_REE_watertap_costing_initialize(self, model, solver):

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

        QGESSCostingData.costing_initialization(model.fs.costing)
        QGESSCostingData.initialize_fixed_OM_costs(model.fs.costing)
        QGESSCostingData.initialize_variable_OM_costs(model.fs.costing)

    @pytest.mark.component
    def test_REE_watertap_costing_solve(self, model, solver):

        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_REE_watertap_costing_results_totalCAPEX(self, model):

        assert value(model.fs.costing.total_BEC) == pytest.approx(50.686, rel=1e-4)

    @pytest.mark.component
    def test_REE_watertap_costing_results_equipmentCAPEX(self, model):

        CE_index_year = "UKy_2019"

        CE_index_units = getattr(
            pyunits, "MUSD_" + CE_index_year
        )  # millions of USD, for base year

        assert value(
            pyunits.convert(
                model.fs_membrane.nfunit.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0015159, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.rounit.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0016148, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.ixunit.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(4.0354, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.nfzounit.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(2.3391, rel=1e-4)

        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs_membrane.nfunit.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs_membrane.rounit.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs_membrane.ixunit.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs_membrane.nfzounit.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(44.308, rel=1e-4)

    @pytest.mark.component
    def test_REE_watertap_costing_results_fixedOPEX(self, model):

        CE_index_year = "UKy_2019"

        CE_index_units = getattr(
            pyunits, "MUSD_" + CE_index_year
        )  # millions of USD, for base year

        assert value(
            pyunits.convert(
                model.fs_membrane.nfunit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00015159, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.rounit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00016148, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.ixunit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.037284, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.nfzounit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.467810, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs_membrane.nfunit.costing.fixed_operating_cost
                + model.fs_membrane.rounit.costing.fixed_operating_cost
                + model.fs_membrane.ixunit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
            + pyunits.convert(
                model.fs_membrane.nfzounit.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.505407, rel=1e-4)

        assert value(model.fs.costing.watertap_fixed_costs) == pytest.approx(
            0.505407, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            11.98986, rel=1e-4
        )

    # TODO commented as no WaterTAP models currently use this, may change in the future
    # @pytest.mark.component
    # def test_REE_watertap_costing_variableOPEX(self, model):

    #     assert value(model.fs.costing.watertap_variable_costs) == pytest.approx(
    #         0, abs=1e-4
    #     )

    #     assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
    #         533.10082, rel=1e-4
    #     )


class TestCustomCosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = base_model()

        # create model
        model.fs.custom_vessel = UnitModelBlock()
        model.fs.custom_vessel.volume = pyo.Var(initialize=5000, units=pyunits.m**3)
        model.fs.custom_vessel.volume.fix()
        model.fs.custom_vessel.water_injection_rate = pyo.Var(
            initialize=0.1, units=pyunits.m**3 / pyunits.s
        )
        model.fs.custom_vessel.water_injection_rate.fix()

        # add costing
        model.fs.custom_vessel.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=CustomCostingData.cost_custom_vessel,
            costing_method_arguments={
                "volume_per_unit": model.fs.custom_vessel.volume,
                "material": "carbonsteel",
                "water_injection_rate_per_unit": model.fs.custom_vessel.water_injection_rate,
                "number_of_units": 3,
            },
        )

        return model

    @pytest.mark.unit
    def test_model(self, model):

        # confirm that base units match the QGESS costing block
        assert (
            model.fs.custom_vessel.costing.costing_package.base_currency
            == pyunits.USD_2021
        )
        assert (
            model.fs.custom_vessel.costing.costing_package.base_period == pyunits.year
        )

        assert isinstance(model.fs.custom_vessel.costing.number_of_units, pyo.Param)
        assert isinstance(model.fs.custom_vessel.costing.capital_cost_per_unit, pyo.Var)
        assert isinstance(model.fs.custom_vessel.costing.capital_cost, pyo.Var)
        assert isinstance(model.fs.custom_vessel.costing.fixed_operating_cost, pyo.Var)
        assert isinstance(
            model.fs.custom_vessel.costing.variable_operating_cost_per_unit, pyo.Var
        )
        assert isinstance(
            model.fs.custom_vessel.costing.variable_operating_cost, pyo.Var
        )

        assert isinstance(
            model.fs.custom_vessel.costing.capital_cost_per_unit_eq, pyo.Constraint
        )
        assert isinstance(
            model.fs.custom_vessel.costing.capital_cost_constraint, pyo.Constraint
        )
        assert isinstance(
            model.fs.custom_vessel.costing.fixed_operating_cost_constraint,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.custom_vessel.costing.variable_operating_cost_per_unit_eq,
            pyo.Constraint,
        )
        assert isinstance(
            model.fs.custom_vessel.costing.variable_operating_cost_constraint,
            pyo.Constraint,
        )

    @pytest.mark.component
    def test_REE_custom_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included
        CE_index_year = "UKy_2019"

        CE_index_units = getattr(
            pyunits, "MUSD_" + CE_index_year
        )  # millions of USD, for base year

        # add plant-level cost constraints

        model.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)

        hours_per_shift = 8
        shifts_per_day = 3
        operating_days_per_year = 336

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
            mutable=False,
            units=pyunits.hours / pyunits.year,
        )

        model.fs.recovery_rate_per_year = pyo.Var(
            initialize=39.3
            * pyunits.kg
            / pyunits.hr
            * 0.8025  # TREO (total rare earth oxide), 80.25% REE in REO
            * model.fs.annual_operating_hours,
            units=pyunits.kg / pyunits.yr,
        )

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * getattr(pyunits, "MUSD_" + CE_index_year)
            / pyunits.ton
            * pyunits.convert(model.fs.feed_input, to_units=pyunits.ton / pyunits.hr)
            * hours_per_shift
            * pyunits.hr
            * shifts_per_day
            * pyunits.day**-1
            * operating_days_per_year
            * pyunits.day
        )

        # dummy reagent with cost of 1 USD/kg for each section
        reagent_costs = (
            (  # all USD/year
                302962  # Crushing and Screening
                + 0  # Dry Grinding
                + 5767543  # Roasting
                + 199053595  # Leaching
                + 152303329  # Rougher Solvent Extraction
                + 43702016  # Cleaner Solvent Extraction
                + 7207168  # Solvent Extraction Wash and Saponification
                + 1233763  # Rare Earth Element Precipiation
                + 18684816  # Water Treatment
            )
            * pyunits.kg
            / pyunits.year
        )

        model.fs.reagents = pyo.Var(
            model.fs.time,
            initialize=reagent_costs / (model.fs.annual_operating_hours),
            units=pyunits.kg / pyunits.hr,
        )

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)

        resources = [
            "dummy",
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
            "power",
        ]

        rates = [
            model.fs.reagents,
            model.fs.solid_waste,
            model.fs.precipitate,
            model.fs.dust_and_volatiles,
            model.fs.power,
        ]

        # define product flowrates

        pure_product_output_rates = {
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.4 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.5 * pyunits.kg / pyunits.hr,
        }

        mixed_product_output_rates = {
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            "Y2O3": 0.05418 * pyunits.kg / pyunits.hr,
            "La2O3": 0.13770 * pyunits.kg / pyunits.hr,
            "CeO2": 0.37383 * pyunits.kg / pyunits.hr,
            "Pr6O11": 0.03941 * pyunits.kg / pyunits.hr,
            "Nd2O3": 0.17289 * pyunits.kg / pyunits.hr,
            "Sm2O3": 0.02358 * pyunits.kg / pyunits.hr,
            "Eu2O3": 0.00199 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tb4O7": 0.00801 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Ho2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Er2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tm2O3": 0.00130 * pyunits.kg / pyunits.hr,
            "Yb2O3": 0.00373 * pyunits.kg / pyunits.hr,
            "Lu2O3": 0.00105 * pyunits.kg / pyunits.hr,
        }

        model.fs.costing.build_process_costs(
            # arguments related to installation costs
            piping_materials_and_labor_percentage=20,
            electrical_materials_and_labor_percentage=20,
            instrumentation_percentage=8,
            plants_services_percentage=10,
            process_buildings_percentage=40,
            auxiliary_buildings_percentage=15,
            site_improvements_percentage=10,
            equipment_installation_percentage=17,
            field_expenses_percentage=12,
            project_management_and_construction_percentage=30,
            process_contingency_percentage=15,
            # argument related to Fixed OM costs
            labor_types=[
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
            ],
            labor_rate=[24.98, 19.08, 30.39, 22.73, 21.97, 45.85],  # USD/hr
            labor_burden=25,  # % fringe benefits
            operators_per_shift=[4, 9, 2, 2, 2, 3],
            hours_per_shift=hours_per_shift,
            shifts_per_day=shifts_per_day,
            operating_days_per_year=operating_days_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            mixed_product_sale_price_realization_factor=0.65,  # 65% price realization for mixed products
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=resources,
            rates=rates,
            prices={
                "dummy": 1 * getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            },
            fixed_OM=True,
            variable_OM=True,
            feed_input=model.fs.feed_input,
            efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
            chemicals=["dummy"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
            ],
            recovery_rate_per_year=model.fs.recovery_rate_per_year,
            CE_index_year=CE_index_year,
        )

        # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
        model.fs.costing.other_plant_costs.unfix()
        model.fs.costing.other_plant_costs_rule = pyo.Constraint(
            expr=(
                model.fs.costing.other_plant_costs
                == pyunits.convert(
                    1218073 * pyunits.USD_2016  # Rougher Solvent Extraction
                    + 48723 * pyunits.USD_2016  # Cleaner Solvent Extraction
                    + 182711
                    * pyunits.USD_2016,  # Solvent Extraction Wash and Saponification
                    to_units=getattr(pyunits, "MUSD_" + CE_index_year),
                )
            )
        )

        # fix costing vars that shouldn't change
        model.fs.feed_input.fix()
        model.fs.feed_grade.fix()
        model.fs.recovery_rate_per_year.fix()
        model.fs.reagents.fix()
        model.fs.solid_waste.fix()
        model.fs.precipitate.fix()
        model.fs.dust_and_volatiles.fix()
        model.fs.power.fix()

        # check model structural diagnostics
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

        QGESSCostingData.costing_initialization(model.fs.costing)
        QGESSCostingData.initialize_fixed_OM_costs(model.fs.costing)
        QGESSCostingData.initialize_variable_OM_costs(model.fs.costing)

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

        # check model numerical diagnostics
        dt.assert_no_numerical_warnings()

        assert value(model.fs.costing.total_BEC) == pytest.approx(44.377, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0686081, rel=1e-4)
        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs.custom_vessel.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(44.308, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00343040, rel=1e-4)

        assert value(model.fs.costing.custom_fixed_costs) == pytest.approx(
            0.00343040, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            10.92575, rel=1e-4
        )

        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.variable_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(4.13750, rel=1e-4)

        assert value(model.fs.costing.custom_variable_costs) == pytest.approx(
            4.13750, rel=1e-4
        )

        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            530.67815, rel=1e-4
        )


class TestDiafiltrationCosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = base_model()

        # create dummy blocks to store the UnitModelCostingBlocks
        m.fs.stage1 = UnitModelBlock()
        m.fs.stage2 = UnitModelBlock()
        m.fs.stage3 = UnitModelBlock()
        m.fs.cascade = UnitModelBlock()  # to cost the pressure drop
        m.fs.feed_pump = UnitModelBlock()  # to cost feed pump
        m.fs.diafiltrate_pump = UnitModelBlock()  # to cost diafiltrate pump

        m.fs.stage1.length = pyo.Var(initialize=10, units=pyunits.m)
        m.fs.stage2.length = pyo.Var(initialize=10, units=pyunits.m)
        m.fs.stage3.length = pyo.Var(initialize=10, units=pyunits.m)
        m.fs.w = pyo.Var(initialize=10, units=pyunits.m)
        m.fs.Jw = pyo.Var(initialize=10, units=pyunits.m / pyunits.h)
        m.fs.stage3.retentate_flow_vol = pyo.Var(
            initialize=10, units=pyunits.m**3 / pyunits.h
        )
        m.fs.stage3.permeate_flow_vol = pyo.Var(
            initialize=10, units=pyunits.m**3 / pyunits.h
        )
        m.fs.P_atm = pyo.Var(initialize=101325, units=pyunits.Pa)
        m.fs.P_op = pyo.Var(initialize=201325, units=pyunits.Pa)

        m.fs.stage1.length.fix()
        m.fs.stage2.length.fix()
        m.fs.stage3.length.fix()
        m.fs.w.fix()
        m.fs.Jw.fix()
        m.fs.stage3.retentate_flow_vol.fix()
        m.fs.stage3.permeate_flow_vol.fix()
        m.fs.P_atm.fix()
        m.fs.P_op.fix()

        # m.fs.costing = DiafiltrationCosting()
        m.fs.stage1.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_membranes,
            costing_method_arguments={
                "membrane_length": m.fs.stage1.length,
                "membrane_width": m.fs.w,
            },
        )
        m.fs.stage2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_membranes,
            costing_method_arguments={
                "membrane_length": m.fs.stage2.length,
                "membrane_width": m.fs.w,
            },
        )
        m.fs.stage3.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_membranes,
            costing_method_arguments={
                "membrane_length": m.fs.stage3.length,
                "membrane_width": m.fs.w,
            },
        )
        m.fs.cascade.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop,
            costing_method_arguments={
                "water_flux": m.fs.Jw,
                "vol_flow_feed": m.fs.stage3.retentate_flow_vol,  # cascade feed
                "vol_flow_perm": m.fs.stage3.permeate_flow_vol,  # cascade permeate
            },
        )
        m.fs.feed_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_pump,
            costing_method_arguments={
                "inlet_pressure": m.fs.P_atm
                + pyunits.convert(
                    m.fs.cascade.costing.pressure_drop, to_units=pyunits.Pa
                ),
                "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
                * pyunits.psi,  # this should make m.fs.feed_pump.costing.fixed_operating_cost ~0
                "inlet_vol_flow": m.fs.stage3.retentate_flow_vol,  # feed
            },
        )
        m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_pump,
            costing_method_arguments={
                "inlet_pressure": m.fs.P_atm,
                "outlet_pressure": m.fs.P_op,
                "inlet_vol_flow": m.fs.stage3.retentate_flow_vol,  # diafiltrate
            },
        )

        return m

    @pytest.mark.unit
    def test_model(self, model):

        # confirm that base units match the QGESS costing block
        for blk in [
            model.fs.stage1,
            model.fs.stage2,
            model.fs.stage3,
            model.fs.cascade,
            model.fs.feed_pump,
            model.fs.diafiltrate_pump,
        ]:
            assert blk.costing.costing_package.base_currency == pyunits.USD_2021
            assert blk.costing.costing_package.base_period == pyunits.year

        assert isinstance(model.fs.stage1.costing.capital_cost, pyo.Var)
        assert isinstance(model.fs.stage1.costing.fixed_operating_cost, pyo.Var)
        assert not hasattr(model.fs.stage1.costing, "variable_operating_cost")
        assert isinstance(
            model.fs.stage1.costing.capital_cost_constraint, pyo.Constraint
        )
        assert isinstance(
            model.fs.stage1.costing.fixed_operating_cost_constraint, pyo.Constraint
        )
        assert not hasattr(
            model.fs.stage1.costing, "variable_operating_cost_constraint"
        )

        assert isinstance(model.fs.stage2.costing.capital_cost, pyo.Var)
        assert isinstance(model.fs.stage2.costing.fixed_operating_cost, pyo.Var)
        assert not hasattr(model.fs.stage2.costing, "variable_operating_cost")
        assert isinstance(
            model.fs.stage2.costing.capital_cost_constraint, pyo.Constraint
        )
        assert isinstance(
            model.fs.stage2.costing.fixed_operating_cost_constraint, pyo.Constraint
        )
        assert not hasattr(
            model.fs.stage2.costing, "variable_operating_cost_constraint"
        )

        assert isinstance(model.fs.stage3.costing.capital_cost, pyo.Var)
        assert isinstance(model.fs.stage3.costing.fixed_operating_cost, pyo.Var)
        assert not hasattr(model.fs.stage3.costing, "variable_operating_cost")
        assert isinstance(
            model.fs.stage3.costing.capital_cost_constraint, pyo.Constraint
        )
        assert isinstance(
            model.fs.stage3.costing.fixed_operating_cost_constraint, pyo.Constraint
        )
        assert not hasattr(
            model.fs.stage3.costing, "variable_operating_cost_constraint"
        )

        assert not hasattr(model.fs.cascade.costing, "capital_cost")
        assert not hasattr(model.fs.cascade.costing, "fixed_operating_cost")
        assert isinstance(model.fs.cascade.costing.variable_operating_cost, pyo.Var)
        assert not hasattr(model.fs.cascade.costing, "capital_cost_constraint")
        assert not hasattr(model.fs.cascade.costing, "fixed_operating_cost_constraint")
        assert isinstance(
            model.fs.cascade.costing.variable_operating_cost_constraint, pyo.Constraint
        )

        assert isinstance(model.fs.feed_pump.costing.capital_cost, pyo.Var)
        assert not hasattr(model.fs.feed_pump.costing, "fixed_operating_cost")
        assert isinstance(model.fs.feed_pump.costing.variable_operating_cost, pyo.Var)
        assert isinstance(
            model.fs.feed_pump.costing.capital_cost_constraint, pyo.Constraint
        )
        assert not hasattr(
            model.fs.feed_pump.costing, "fixed_operating_cost_constraint"
        )
        assert isinstance(
            model.fs.feed_pump.costing.variable_operating_cost_constraint,
            pyo.Constraint,
        )

        assert isinstance(model.fs.diafiltrate_pump.costing.capital_cost, pyo.Var)
        assert not hasattr(model.fs.diafiltrate_pump.costing, "fixed_operating_cost")
        assert isinstance(
            model.fs.diafiltrate_pump.costing.variable_operating_cost, pyo.Var
        )
        assert isinstance(
            model.fs.diafiltrate_pump.costing.capital_cost_constraint, pyo.Constraint
        )
        assert not hasattr(
            model.fs.diafiltrate_pump.costing, "fixed_operating_cost_constraint"
        )
        assert isinstance(
            model.fs.diafiltrate_pump.costing.variable_operating_cost_constraint,
            pyo.Constraint,
        )

    @pytest.mark.component
    def test_REE_custom_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included
        CE_index_year = "UKy_2019"

        CE_index_units = getattr(
            pyunits, "MUSD_" + CE_index_year
        )  # millions of USD, for base year

        # add plant-level cost constraints

        model.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)

        hours_per_shift = 8
        shifts_per_day = 3
        operating_days_per_year = 336

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
            mutable=False,
            units=pyunits.hours / pyunits.year,
        )

        model.fs.recovery_rate_per_year = pyo.Var(
            initialize=39.3
            * pyunits.kg
            / pyunits.hr
            * 0.8025  # TREO (total rare earth oxide), 80.25% REE in REO
            * model.fs.annual_operating_hours,
            units=pyunits.kg / pyunits.yr,
        )

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * getattr(pyunits, "MUSD_" + CE_index_year)
            / pyunits.ton
            * pyunits.convert(model.fs.feed_input, to_units=pyunits.ton / pyunits.hr)
            * hours_per_shift
            * pyunits.hr
            * shifts_per_day
            * pyunits.day**-1
            * operating_days_per_year
            * pyunits.day
        )

        # dummy reagent with cost of 1 USD/kg for each section
        reagent_costs = (
            (  # all USD/year
                302962  # Crushing and Screening
                + 0  # Dry Grinding
                + 5767543  # Roasting
                + 199053595  # Leaching
                + 152303329  # Rougher Solvent Extraction
                + 43702016  # Cleaner Solvent Extraction
                + 7207168  # Solvent Extraction Wash and Saponification
                + 1233763  # Rare Earth Element Precipiation
                + 18684816  # Water Treatment
            )
            * pyunits.kg
            / pyunits.year
        )

        model.fs.reagents = pyo.Var(
            model.fs.time,
            initialize=reagent_costs / (model.fs.annual_operating_hours),
            units=pyunits.kg / pyunits.hr,
        )

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)

        resources = [
            "dummy",
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
            "power",
        ]

        rates = [
            model.fs.reagents,
            model.fs.solid_waste,
            model.fs.precipitate,
            model.fs.dust_and_volatiles,
            model.fs.power,
        ]

        # define product flowrates

        pure_product_output_rates = {
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.4 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.5 * pyunits.kg / pyunits.hr,
        }

        mixed_product_output_rates = {
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            "Y2O3": 0.05418 * pyunits.kg / pyunits.hr,
            "La2O3": 0.13770 * pyunits.kg / pyunits.hr,
            "CeO2": 0.37383 * pyunits.kg / pyunits.hr,
            "Pr6O11": 0.03941 * pyunits.kg / pyunits.hr,
            "Nd2O3": 0.17289 * pyunits.kg / pyunits.hr,
            "Sm2O3": 0.02358 * pyunits.kg / pyunits.hr,
            "Eu2O3": 0.00199 * pyunits.kg / pyunits.hr,
            "Gd2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tb4O7": 0.00801 * pyunits.kg / pyunits.hr,
            "Dy2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Ho2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Er2O3": 0.00000 * pyunits.kg / pyunits.hr,
            "Tm2O3": 0.00130 * pyunits.kg / pyunits.hr,
            "Yb2O3": 0.00373 * pyunits.kg / pyunits.hr,
            "Lu2O3": 0.00105 * pyunits.kg / pyunits.hr,
        }

        model.fs.costing.build_process_costs(
            # arguments related to installation costs
            piping_materials_and_labor_percentage=20,
            electrical_materials_and_labor_percentage=20,
            instrumentation_percentage=8,
            plants_services_percentage=10,
            process_buildings_percentage=40,
            auxiliary_buildings_percentage=15,
            site_improvements_percentage=10,
            equipment_installation_percentage=17,
            field_expenses_percentage=12,
            project_management_and_construction_percentage=30,
            process_contingency_percentage=15,
            # argument related to Fixed OM costs
            labor_types=[
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
            ],
            labor_rate=[24.98, 19.08, 30.39, 22.73, 21.97, 45.85],  # USD/hr
            labor_burden=25,  # % fringe benefits
            operators_per_shift=[4, 9, 2, 2, 2, 3],
            hours_per_shift=hours_per_shift,
            shifts_per_day=shifts_per_day,
            operating_days_per_year=operating_days_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            mixed_product_sale_price_realization_factor=0.65,  # 65% price realization for mixed products
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=resources,
            rates=rates,
            prices={
                "dummy": 1 * getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
            },
            fixed_OM=True,
            variable_OM=True,
            feed_input=model.fs.feed_input,
            efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
            chemicals=["dummy"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
            ],
            recovery_rate_per_year=model.fs.recovery_rate_per_year,
            CE_index_year=CE_index_year,
        )

        # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
        model.fs.costing.other_plant_costs.unfix()
        model.fs.costing.other_plant_costs_rule = pyo.Constraint(
            expr=(
                model.fs.costing.other_plant_costs
                == pyunits.convert(
                    1218073 * pyunits.USD_2016  # Rougher Solvent Extraction
                    + 48723 * pyunits.USD_2016  # Cleaner Solvent Extraction
                    + 182711
                    * pyunits.USD_2016,  # Solvent Extraction Wash and Saponification
                    to_units=getattr(pyunits, "MUSD_" + CE_index_year),
                )
            )
        )

        # fix costing vars that shouldn't change
        model.fs.feed_input.fix()
        model.fs.feed_grade.fix()
        model.fs.recovery_rate_per_year.fix()
        model.fs.reagents.fix()
        model.fs.solid_waste.fix()
        model.fs.precipitate.fix()
        model.fs.dust_and_volatiles.fix()
        model.fs.power.fix()

        # check model structural diagnostics
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

        QGESSCostingData.costing_initialization(model.fs.costing)
        QGESSCostingData.initialize_fixed_OM_costs(model.fs.costing)
        QGESSCostingData.initialize_variable_OM_costs(model.fs.costing)

        solver = get_solver()
        results = solver.solve(model, tee=False)
        assert_optimal_termination(results)

        # check model numerical diagnostics
        dt.assert_no_numerical_warnings()

        assert value(model.fs.costing.total_BEC) == pytest.approx(44.684, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage1.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage2.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage3.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.feed_pump.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.34788, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.diafiltrate_pump.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(0.014780, rel=1e-4)
        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs.stage1.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs.stage2.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs.stage3.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs.feed_pump.costing.capital_cost, to_units=CE_index_units
            )
            - pyunits.convert(
                model.fs.diafiltrate_pump.costing.capital_cost, to_units=CE_index_units
            )
        ) == pytest.approx(44.308, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.stage1.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage2.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage3.costing.fixed_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)

        assert value(model.fs.costing.custom_fixed_costs) == pytest.approx(
            0.0025826, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            10.95225, rel=1e-4
        )

        assert value(
            pyunits.convert(
                model.fs.cascade.costing.variable_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(0.985221, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.feed_pump.costing.variable_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(8.09845e-17, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.diafiltrate_pump.costing.variable_operating_cost,
                to_units=CE_index_units / pyunits.year,
            )
        ) == pytest.approx(2.36473e-10, rel=1e-4)

        assert value(model.fs.costing.custom_variable_costs) == pytest.approx(
            2.36473e-10, rel=1e-4
        )

        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            525.7185, rel=1e-4
        )


class TestHDDRecyclingCosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()

        # add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = QGESSCosting()
        CE_index_year = "2019"

        # Source 1, 1.1 is Front End Loader (2 cuyd)
        # this is a constant-cost unit, where n_equip is the scaling parameter
        CS_front_end_loader_2yd3_accounts = ["1.1"]
        m.fs.CS_front_end_loader_2yd3 = UnitModelBlock()
        m.fs.CS_front_end_loader_2yd3.n_equip = pyo.Var(
            initialize=1, units=pyunits.dimensionless
        )
        m.fs.CS_front_end_loader_2yd3.n_equip.fix()

        m.fs.CS_front_end_loader_2yd3.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_front_end_loader_2yd3_accounts,
                "scaled_param": m.fs.CS_front_end_loader_2yd3.n_equip,  # 1 loader
                "source": 1,
                # no. units is the scaling parameter for constant-cost units,
                #     so use n_equip below to specify the number of loaders
                "n_equip": 5,
                "scale_down_parallel_equip": False,
                "CE_index_year": CE_index_year,
            },
        )

        # Source 2, 2.1 is HDD shredder
        # this is a constant-cost unit, where n_equip is the scaling parameter
        HDD_Recycling_shredder_accounts = ["2.1"]
        m.fs.HDD_Recycling_shredder = UnitModelBlock()
        m.fs.HDD_Recycling_shredder.n_equip = pyo.Var(
            initialize=1, units=pyunits.dimensionless
        )
        m.fs.HDD_Recycling_shredder.n_equip.fix()

        m.fs.HDD_Recycling_shredder.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": HDD_Recycling_shredder_accounts,
                "scaled_param": m.fs.HDD_Recycling_shredder.n_equip,  # 1 shredder
                "source": 2,
                # no. units is the scaling parameter for constant-cost units,
                # so use n_equip below to specify the number of loaders
                "n_equip": 1,
                "scale_down_parallel_equip": False,
                "CE_index_year": CE_index_year,
            },
        )

        # Source 2, 2.2 is Hydrogen Decrepitation
        HDD_Recycling_HD_accounts = ["2.2"]
        m.fs.HDD_Recycling_HD = UnitModelBlock()
        m.fs.HDD_Recycling_HD.duty = pyo.Var(
            initialize=10, units=pyunits.MBTU / pyunits.hr
        )
        m.fs.HDD_Recycling_HD.duty.fix()
        m.fs.HDD_Recycling_HD.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": HDD_Recycling_HD_accounts,
                "scaled_param": m.fs.HDD_Recycling_HD.duty,
                "source": 2,
                "n_equip": 1,
                "scale_down_parallel_equip": False,
                "CE_index_year": CE_index_year,
            },
        )

        return m

    @pytest.mark.unit
    def test_HDD_Recycling_costing_build_diagnostics(self, model):
        CE_index_year = "2019"
        model.fs.costing.build_process_costs(
            CE_index_year=CE_index_year,
            # defaults to fixed_OM=True, so explicitly set to False
            # defaults to variable_OM=False, so let that use the default
            fixed_OM=False,
        )

        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_HDD_Recycling_costing_initialize(self, model):
        QGESSCostingData.costing_initialization(model.fs.costing)

    @pytest.mark.component
    def test_HDD_Recycling_costing_solve(self, model):
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_HDD_Recycling_costing_solve_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_HDD_Recycling_costing_results(self, model):
        CS_front_end_loader_2yd3_accounts = ["1.1"]
        HDD_Recycling_shredder_accounts = ["2.1"]
        HDD_Recycling_HD_accounts = ["2.2"]

        assert value(
            model.fs.CS_front_end_loader_2yd3.costing.bare_erected_cost[
                CS_front_end_loader_2yd3_accounts
            ]
        ) == pytest.approx(0.82652, rel=1e-4)
        assert value(
            model.fs.HDD_Recycling_shredder.costing.bare_erected_cost[
                HDD_Recycling_shredder_accounts
            ]
        ) == pytest.approx(0.05000, rel=1e-4)
        assert value(
            model.fs.HDD_Recycling_HD.costing.bare_erected_cost[
                HDD_Recycling_HD_accounts
            ]
        ) == pytest.approx(0.41035, rel=1e-4)

        assert value(model.fs.costing.total_plant_cost) == pytest.approx(
            3.8220, rel=1e-4
        )


class TestNPVCostingBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = QGESSCosting(
            discount_percentage=10,  # percent
            plant_lifetime=20,  # years
            has_capital_expenditure_period=True,
        )

        # 1.3 is CS Jaw Crusher
        CS_jaw_crusher_accounts = ["1.3"]
        m.fs.CS_jaw_crusher = UnitModelBlock()
        m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
        m.fs.CS_jaw_crusher.power.fix()
        m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "source": 1,
            },
        )

        m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        m.fs.feed_input.fix()

        m.fs.water = pyo.Var(
            m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr
        )
        m.fs.water.fix()

        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
            calculate_NPV=True,
        )

        return m

    @pytest.mark.unit
    def test_NPV_costingblock_build(self, model):
        # check that some objects are built as expected
        assert isinstance(model.fs.costing.pv_capital_cost, pyo.Var)
        assert isinstance(model.fs.costing.loan_debt, pyo.Var)
        assert isinstance(model.fs.costing.pv_loan_interest, pyo.Var)
        assert isinstance(model.fs.costing.pv_operating_cost, pyo.Var)
        assert isinstance(model.fs.costing.pv_revenue, pyo.Var)
        assert isinstance(model.fs.costing.npv, pyo.Var)
        assert not hasattr(model.fs.costing, "pv_taxes")
        assert isinstance(model.fs.costing.discount_percentage, pyo.Param)
        assert isinstance(model.fs.costing.plant_lifetime, pyo.Param)
        assert isinstance(model.fs.costing.capital_expenditure_percentages, pyo.Param)
        assert isinstance(model.fs.costing.capital_escalation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.capital_loan_interest_percentage, pyo.Param)
        assert isinstance(model.fs.costing.capital_loan_repayment_period, pyo.Param)
        assert isinstance(model.fs.costing.debt_percentage_of_CAPEX, pyo.Param)
        assert isinstance(model.fs.costing.operating_inflation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.revenue_inflation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.pv_capital_cost_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.loan_debt_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_loan_interest_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_operating_cost_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_revenue_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.npv_constraint, pyo.Constraint)
        assert not hasattr(model.fs.costing, "pv_taxes_constraint")

    @pytest.mark.unit
    def test_NPV_costingblock_build_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_NPV_costingblock_initialize(self, model):
        QGESSCostingData.costing_initialization(model.fs.costing)

    @pytest.mark.component
    def test_NPV_costingblock_solve(self, model):
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_NPV_costingblock_solve_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_NPV_costingblock_results(self, model):
        # check some NPV results
        assert value(model.fs.costing.pv_capital_cost) == pytest.approx(
            -6.3162037, rel=1e-4
        )
        assert value(model.fs.costing.pv_loan_interest) == pytest.approx(
            -0.61610712, rel=1e-4
        )
        assert value(model.fs.costing.pv_operating_cost) == pytest.approx(
            -57.28850, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(239.63402, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(175.4132, rel=1e-4)


class TestNPVFixedInputs(object):
    @pytest.fixture(scope="class")
    def model(self):
        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = QGESSCosting(
            discount_percentage=10,
            plant_lifetime=20,
            # use CAPEX, OPEX, REVENUE from CostingBlock test to verify results are the same
            total_capital_cost=7.461172417869669,
            annual_operating_cost=6.880158261340908,
            annual_revenue=64.38220104959998,
            has_capital_expenditure_period=True,
            cost_year="2021",
        )

        m.fs.costing.build_process_costs(
            fixed_OM=False,
            calculate_NPV=True,
        )

        return m

    @pytest.mark.unit
    def test_NPV_fixedinputs_build(self, model):
        # check that some objects are built as expected
        assert isinstance(model.fs.costing.pv_capital_cost, pyo.Var)
        assert isinstance(model.fs.costing.loan_debt, pyo.Var)
        assert isinstance(model.fs.costing.pv_loan_interest, pyo.Var)
        assert isinstance(model.fs.costing.pv_operating_cost, pyo.Var)
        assert isinstance(model.fs.costing.pv_revenue, pyo.Var)
        assert isinstance(model.fs.costing.npv, pyo.Var)
        assert not hasattr(model.fs.costing, "pv_taxes")
        assert isinstance(model.fs.costing.discount_percentage, pyo.Param)
        assert isinstance(model.fs.costing.plant_lifetime, pyo.Param)
        assert isinstance(model.fs.costing.capital_expenditure_percentages, pyo.Param)
        assert isinstance(model.fs.costing.capital_escalation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.capital_loan_interest_percentage, pyo.Param)
        assert isinstance(model.fs.costing.capital_loan_repayment_period, pyo.Param)
        assert isinstance(model.fs.costing.debt_percentage_of_CAPEX, pyo.Param)
        assert isinstance(model.fs.costing.operating_inflation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.revenue_inflation_percentage, pyo.Param)
        assert isinstance(model.fs.costing.pv_capital_cost_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.loan_debt_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_loan_interest_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_operating_cost_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.pv_revenue_constraint, pyo.Constraint)
        assert isinstance(model.fs.costing.npv_constraint, pyo.Constraint)
        assert not hasattr(model.fs.costing, "pv_taxes_constraint")

    @pytest.mark.unit
    def test_NPV_fixedinputs_build_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_NPV_fixedinputs_solve(self, model):
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_NPV_fixedinputs_solve_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_NPV_fixedinputs_results(self, model):
        # check some NPV results
        assert value(model.fs.costing.pv_capital_cost) == pytest.approx(
            -6.3162037, rel=1e-4
        )
        assert value(model.fs.costing.pv_loan_interest) == pytest.approx(
            -0.61610712, rel=1e-4
        )
        assert value(model.fs.costing.pv_operating_cost) == pytest.approx(
            -59.02935, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(552.37672, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(486.41505, rel=1e-4)


@pytest.mark.component
def test_REE_costing_CE_index_year():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_2021 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2021.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2021.power.fix()
    m.fs.CS_jaw_crusher_2021.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2021.power,
            "source": 1,
            "CE_index_year": "2021",
        },
    )

    m.fs.CS_jaw_crusher_UKy_2019 = UnitModelBlock()
    m.fs.CS_jaw_crusher_UKy_2019.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_UKy_2019.power.fix()
    m.fs.CS_jaw_crusher_UKy_2019.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_UKy_2019.power,
            "source": 1,
            "CE_index_year": "UKy_2019",
        },
    )

    m.fs.CS_jaw_crusher_2022 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2022.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2022.power.fix()
    m.fs.CS_jaw_crusher_2022.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2022.power,
            "source": 1,
            "CE_index_year": "2022",
        },
    )

    m.fs.CS_jaw_crusher_2025 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2025.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2025.power.fix()
    m.fs.CS_jaw_crusher_2025.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2025.power,
            "source": 1,
            "CE_index_year": "2025",
        },
    )

    m.fs.CS_jaw_crusher_CE500 = UnitModelBlock()
    m.fs.CS_jaw_crusher_CE500.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_CE500.power.fix()
    m.fs.CS_jaw_crusher_CE500.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_CE500.power,
            "source": 1,
            "CE_index_year": "CE500",
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that total plant cost ratios match expected currency conversions
    for account in CS_jaw_crusher_accounts:
        # USD_2021 = CE 708.0
        assert (
            pytest.approx(
                value(
                    m.fs.CS_jaw_crusher_2021.costing.bare_erected_cost[account]
                    / m.fs.CS_jaw_crusher_CE500.costing.bare_erected_cost[account]
                ),
                abs=1e-1,
            )
            == 708.0 / 500
        )

        # USD_UKy_2019 = CE 609.495
        assert (
            pytest.approx(
                value(
                    m.fs.CS_jaw_crusher_UKy_2019.costing.bare_erected_cost[account]
                    / m.fs.CS_jaw_crusher_CE500.costing.bare_erected_cost[account]
                ),
                abs=1e-1,
            )
            == 609.495 / 500
        )

        # USD_2022 = CE 816.0
        assert (
            pytest.approx(
                value(
                    m.fs.CS_jaw_crusher_2022.costing.bare_erected_cost[account]
                    / m.fs.CS_jaw_crusher_CE500.costing.bare_erected_cost[account]
                ),
                abs=1e-1,
            )
            == 816.0 / 500
        )

        # USD_2025 = CE 815.59
        assert (
            pytest.approx(
                value(
                    m.fs.CS_jaw_crusher_2025.costing.bare_erected_cost[account]
                    / m.fs.CS_jaw_crusher_CE500.costing.bare_erected_cost[account]
                ),
                abs=1e-1,
            )
            == 815.59 / 500
        )


@pytest.mark.unit
def test_REE_costing_disallowedunitmodelcostunits():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        AttributeError,
        match="CE_index_year notayear is not a valid currency base option. "
        "Valid CE index options include CE500, CE394 and years from 1990 to 2020.",
    ):
        m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "source": 1,
                "CE_index_year": "notayear",
            },
        )


@pytest.mark.unit
def test_REE_costing_nonexistentcostaccount():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["notanaccount"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        KeyError,
        match="Account notanaccount could not be found in the dictionary for source 1",
    ):
        m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "source": 1,
            },
        )


@pytest.mark.component
def test_REE_costing_multipleaccountssameparameter():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher, 1.5 is Roll Crusher
    CS_crusher_accounts = ["1.3", "1.5"]
    m.fs.CS_crusher = UnitModelBlock()
    m.fs.CS_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_crusher.power.fix()

    m.fs.CS_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_crusher_accounts,
            "scaled_param": m.fs.CS_crusher.power,
            "source": 1,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # 1.5 is CS Roll Crusher
    CS_roll_crusher_accounts = ["1.5"]
    m.fs.CS_roll_crusher = UnitModelBlock()
    m.fs.CS_roll_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_roll_crusher.power.fix()
    m.fs.CS_roll_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_roll_crusher_accounts,
            "scaled_param": m.fs.CS_roll_crusher.power,
            "source": 1,
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    assert value(m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3"]) == pytest.approx(
        2.5122, rel=1e-4
    )
    assert value(m.fs.CS_crusher.costing.bare_erected_cost["1.3"]) == pytest.approx(
        2.5122, rel=1e-4
    )
    assert value(
        m.fs.CS_roll_crusher.costing.bare_erected_cost["1.5"]
    ) == pytest.approx(0.32769, rel=1e-4)
    assert value(m.fs.CS_crusher.costing.bare_erected_cost["1.5"]) == pytest.approx(
        0.32769, rel=1e-4
    )


@pytest.mark.unit
def test_REE_costing_multipleaccountsdifferentparameters():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher, 1.6 is Vibrating Screen
    CS_crusher_accounts = ["1.3", "1.6"]
    m.fs.CS_crusher = UnitModelBlock()
    m.fs.CS_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_crusher.power.fix()

    with pytest.raises(
        ValueError,
        match="fs.CS_crusher.costing cost accounts selected do not use the same process parameter",
    ):
        m.fs.CS_crusher.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_crusher_accounts,
                "scaled_param": m.fs.CS_crusher.power,
                "source": 1,
            },
        )


@pytest.mark.component
def test_REE_costing_additionalcostingparams_newaccount():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "1": {
            "1.3new": {
                "Account Name": "UKy Crushing and Screening - Jaw Crusher",
                "BEC": 1922101.0,
                "BEC_units": "$2016",
                "Exponent": 1.25,
                "Process Parameter": "Power Draw (hp)",
                "RP Value": 589.0,
                "Units": "hp",
            },
        },
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3new"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
            "additional_costing_params": additional_costing_params,
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # adding a check just to make sure everything works as expected
    assert value(
        m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3new"]
    ) == pytest.approx(2.5122, rel=1e-4)


@pytest.mark.component
def test_REE_costing_additionalcostingparams_overwrite():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "1": {
            "1.3": {
                "Account Name": "UKy Crushing and Screening - Jaw Crusher",
                "BEC": 1922101.0,
                "BEC_units": "$2016",
                "Exponent": 1.25,
                "Process Parameter": "Power Draw (hp)",
                "RP Value": 689.0,  # changed this value from 589.0
                "Units": "hp",
            },
        },
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
            "additional_costing_params": additional_costing_params,
            "use_additional_costing_params": True,
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # adding a check just to make sure 1.3 was overwritten before it was used
    assert value(m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3"]) == pytest.approx(
        2.0650, rel=1e-4
    )


@pytest.mark.unit
def test_REE_costing_additionalcostingparams_nooverwrite():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "1": {
            "1.3": {
                "Account Name": "UKy Crushing and Screening - Jaw Crusher",
                "BEC": 1922101.0,
                "BEC_units": "$2016",
                "Exponent": 1.25,
                "Process Parameter": "Power Draw (hp)",
                "RP Value": 689.0,  # changed this value from 589.0
                "Units": "hp",
            },
        },
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        ValueError,
        match="Data already exists for Account 1.3 using source 1. Please "
        "confirm that the custom account dictionary is correct, or add the "
        "new parameters as a new account. To use the custom account dictionary "
        "for all conflicts, please pass the argument use_additional_costing_params "
        "as True.",
    ):
        m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_REE_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "source": 1,
                "additional_costing_params": additional_costing_params,
                "use_additional_costing_params": False,
            },
        )


@pytest.mark.component
def test_REE_costing_scaledownparallelequip():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_1 = UnitModelBlock()
    m.fs.CS_jaw_crusher_1.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_1.power.fix()
    m.fs.CS_jaw_crusher_1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_1.power,
            "source": 1,
            "scale_down_parallel_equip": True,
            "n_equip": 1,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_2 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2.power.fix()
    m.fs.CS_jaw_crusher_2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2.power,
            "source": 1,
            "scale_down_parallel_equip": False,
            "n_equip": 1,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_3 = UnitModelBlock()
    m.fs.CS_jaw_crusher_3.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_3.power.fix()
    m.fs.CS_jaw_crusher_3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_3.power,
            "source": 1,
            "scale_down_parallel_equip": True,
            "n_equip": 2,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_4 = UnitModelBlock()
    m.fs.CS_jaw_crusher_4.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_4.power.fix()
    m.fs.CS_jaw_crusher_4.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_4.power,
            "source": 1,
            "scale_down_parallel_equip": False,
            "n_equip": 2,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_5 = UnitModelBlock()
    m.fs.CS_jaw_crusher_5.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_5.power.fix()
    m.fs.CS_jaw_crusher_5.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_5.power,
            "source": 1,
            "scale_down_parallel_equip": True,
            "n_equip": 5,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_6 = UnitModelBlock()
    m.fs.CS_jaw_crusher_6.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_6.power.fix()
    m.fs.CS_jaw_crusher_6.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_6.power,
            "source": 1,
            "scale_down_parallel_equip": False,
            "n_equip": 5,
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # base case
    assert value(
        m.fs.CS_jaw_crusher_1.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(1.0000, rel=1e-4)

    # only one unit, parallel doesn't change result
    assert value(
        m.fs.CS_jaw_crusher_2.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(1.0000, rel=1e-4)

    # same capacity over two units should be slightly less expensive than one unit
    assert value(
        m.fs.CS_jaw_crusher_3.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(0.84095, rel=1e-4)

    # two units with base case capacity should be double the cost
    assert value(
        m.fs.CS_jaw_crusher_4.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(2.0000, rel=1e-4)

    # same capacity over five units should be much less expensive than one unit
    assert value(
        m.fs.CS_jaw_crusher_5.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(0.66878, rel=1e-4)

    # five units with base case capacity should be five times the cost
    assert value(
        m.fs.CS_jaw_crusher_6.costing.bare_erected_cost["1.3"]
    ) == pytest.approx(5.0000, rel=1e-4)


@pytest.mark.unit
def test_REE_costing_disallowedbuildprocesscostunits():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        AttributeError,
        match="CE_index_year notayear is not a valid currency base option. "
        "Valid CE index options include CE500, CE394 and years from 1990 to 2020.",
    ):
        # defaults to fixed_OM=True, so explicitly set to False
        # defaults to variable_OM=False, so let that use the default
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            CE_index_year="notayear",
        )


@pytest.mark.component
def test_REE_costing_usersetTPC_noOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        total_purchase_cost=1,
        fixed_OM=False,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")  # built using passed TPC
    assert not hasattr(m.fs.costing, "total_BEC_eq")  # shouldn't be built
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert value(m.fs.costing.total_BEC) == pytest.approx(1.0000, rel=1e-4)
    assert value(m.fs.costing.total_plant_cost) == pytest.approx(2.9700, rel=1e-4)


@pytest.mark.component
def test_REE_costing_usersetTPC_withOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        total_purchase_cost=1,
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")  # built using passed TPC
    assert not hasattr(m.fs.costing, "total_BEC_eq")  # shouldn't be built
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert value(m.fs.costing.total_BEC) == pytest.approx(1.0000, rel=1e-4)
    assert value(m.fs.costing.total_plant_cost) == pytest.approx(2.9700, rel=1e-4)


@pytest.mark.component
def test_REE_costing_useLangfactor():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        Lang_factor=2.97,
        fixed_OM=False,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_BEC_eq")
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")

    # if piping materials and labor do not exist, all other factors do not exist as well
    assert not hasattr(m.fs.costing, "piping_materials_and_labor_percentage")
    # instead of percentages, Lang factor should be built
    assert hasattr(m.fs.costing, "Lang_factor")

    assert value(m.fs.costing.total_BEC) == pytest.approx(2.5122, rel=1e-4)
    assert value(m.fs.costing.total_plant_cost) == pytest.approx(7.4612, rel=1e-4)


# optional cost arguments - land, additional chemicals, additional waste
cost_obj_dict = {
    "Expression_withunits": pyo.Expression(expr=1e-3 * pyunits.MUSD_2021),
    "Expression_nounits": pyo.Expression(expr=1e-3),
    "Nonexpression_withunits": pyo.Var(initialize=1e-3, units=pyunits.MUSD_2021),
    "Nonexpression_nounits": pyo.Var(initialize=1e-3),
}


@pytest.mark.parametrize("cost_obj", cost_obj_dict.keys())
@pytest.mark.parametrize(
    "argument", ["land_cost", "additional_chemicals_cost", "additional_waste_cost"]
)
@pytest.mark.component
def test_REE_costing_optionalexpressionarguments(argument, cost_obj):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # define optional argument - land_cost, additional_chemicals_cost, or additional_waste_cost
    setattr(m.fs, argument, cost_obj_dict[cost_obj])
    if isinstance(getattr(m.fs, argument), pyo.Var):
        getattr(m.fs, argument).fix()  # pyo.Var objects must be fixed

    # required for variable cost calculations
    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()
    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        # pass the test argument - land_cost, additional_chemicals_cost, or waste_cost
        **{argument: getattr(m.fs, argument)},
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert hasattr(m.fs.costing, argument)
    assert value(getattr(m.fs.costing, argument)) == pytest.approx(1e-3, rel=1e-4)
    assert pyunits.get_units(getattr(m.fs.costing, argument)) == pyunits.MUSD_2021

    # clean up for subsequent parameterization runs
    delattr(m.fs, argument)


@pytest.mark.component
def test_REE_costing_fixedOM_defaults():
    # use as many default values as possible, only pass required arguments
    # will build fixed O&M costs but not variable O&M costs

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # fixed_OM defaults to True, use default
    m.fs.costing.build_process_costs(
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "annual_operating_labor_cost")
    assert hasattr(m.fs.costing, "annual_technical_labor_cost")
    assert hasattr(m.fs.costing, "annual_labor_cost")
    assert hasattr(m.fs.costing, "maintenance_and_material_cost")
    assert hasattr(m.fs.costing, "quality_assurance_and_control_cost")
    assert hasattr(m.fs.costing, "sales_patenting_and_research_cost")
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost")
    assert hasattr(m.fs.costing, "property_taxes_and_insurance_cost")
    assert hasattr(m.fs.costing, "other_fixed_costs")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_sales_revenue")
    assert not hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert not hasattr(m.fs.costing, "cost_of_recovery")
    assert not hasattr(m.fs.costing, "total_variable_OM_cost")
    assert not hasattr(m.fs.costing, "plant_overhead_cost")

    assert value(m.fs.costing.annual_operating_labor_cost) == pytest.approx(
        3.0730, rel=1e-4
    )
    assert value(m.fs.costing.annual_technical_labor_cost) == pytest.approx(
        1.1801, rel=1e-4
    )
    assert value(m.fs.costing.annual_labor_cost) == pytest.approx(4.2531, rel=1e-4)
    assert value(m.fs.costing.maintenance_and_material_cost) == pytest.approx(
        0.14922, rel=1e-4
    )
    assert value(m.fs.costing.quality_assurance_and_control_cost) == pytest.approx(
        0.30730, rel=1e-4
    )
    assert value(m.fs.costing.sales_patenting_and_research_cost) == pytest.approx(
        0.13965, rel=1e-4
    )
    assert value(m.fs.costing.admin_and_support_labor_cost) == pytest.approx(
        0.61460, rel=1e-4
    )
    assert value(m.fs.costing.property_taxes_and_insurance_cost) == pytest.approx(
        0.074612, rel=1e-4
    )
    assert value(m.fs.costing.other_fixed_costs) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.total_fixed_OM_cost) == pytest.approx(5.5384, rel=1e-4)
    assert value(m.fs.costing.total_sales_revenue) == pytest.approx(27.931, rel=1e-4)


@pytest.mark.unit
def test_REE_costing_fixedOM_twiceonsamemodel():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    # fixed_OM defaults to True, use default
    m.fs.costing.build_process_costs(
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
    )

    with pytest.raises(
        RuntimeError,
        match="Costing for the block fs.costing already exists. Please ensure "
        "that the costing build method is not called twice on the same "
        "model.",
    ):
        # call costing a second time
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_pureproductnotpassed():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_process_costs()


@pytest.mark.unit
def test_REE_costing_fixedOM_pureproductnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates=[
                1.9 * pyunits.kg / pyunits.hr,
            ],
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_mixedproductnotpassed():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_mixedproductnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates=[
                0.00143 * pyunits.kg / pyunits.hr,
            ],
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_salesprice():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.costing.build_process_costs(
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        sale_prices={"newprod": 1e-6 * pyunits.MUSD_2021 / pyunits.kg},
    )


@pytest.mark.unit
def test_REE_costing_fixedOM_salespricenotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        TypeError,
        match="Dictionary of custom sale_prices must be a dict object.",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            sale_prices=[1e-6 * pyunits.MUSD_2021 / pyunits.kg],
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_pureproductnoprice():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        AttributeError,
        match="A pure product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Al', 'Sb', 'As', 'Ba', 'Be', 'Bi', 'Ce', 'Cs', 'Cr', 'Co', 'Dy', 'Er', 'Eu', 'CaF2', 'Gd', 'Ga', 'Ge', 'C', 'Ha', 'Ho', 'In', 'Ir', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Nd', 'Ni', 'Nb', 'Pd', 'Pt', 'Pr', 'Rh', 'Rb', 'Ru', 'Sm', 'Sc', 'Ta', 'Te', 'Tb', 'Tm', 'Sn', 'Ti', 'W', 'V', 'Yb', 'Y', 'Zn', 'Zr', 'CeO2', 'Dy2O3', 'Eu2O3', 'La2O3', 'Nd2O3', 'Sc2O3', 'Ta2O5', 'Tb4O7', 'TiO2', 'WO3', 'Y2O3', 'Er2O3', 'Ho2O3', 'Gd2O3', 'Lu2O3', 'Pr6O11', 'Sm2O3', 'Tm2O3', 'Yb2O3'\\]",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "newprod": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_mixedproductnoprice():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        ValueError,
        match="Value newtype for labor_type is not allowed. "
        "Allowed labor types for operating labor include skilled,"
        "unskilled, supervisor and maintenance. Allowed labor types "
        "for direct labor include technician and engineer.",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            labor_types=[
                "skilled",
                "unskilled",
                "supervisor",
                "maintenance",
                "technician",
                "engineer",
                "newtype",
            ],
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_disallowedlabortype():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    with pytest.raises(
        AttributeError,
        match="A mixed product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Al', 'Sb', 'As', 'Ba', 'Be', 'Bi', 'Ce', 'Cs', 'Cr', 'Co', 'Dy', 'Er', 'Eu', 'CaF2', 'Gd', 'Ga', 'Ge', 'C', 'Ha', 'Ho', 'In', 'Ir', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Nd', 'Ni', 'Nb', 'Pd', 'Pt', 'Pr', 'Rh', 'Rb', 'Ru', 'Sm', 'Sc', 'Ta', 'Te', 'Tb', 'Tm', 'Sn', 'Ti', 'W', 'V', 'Yb', 'Y', 'Zn', 'Zr', 'CeO2', 'Dy2O3', 'Eu2O3', 'La2O3', 'Nd2O3', 'Sc2O3', 'Ta2O5', 'Tb4O7', 'TiO2', 'WO3', 'Y2O3', 'Er2O3', 'Ho2O3', 'Gd2O3', 'Lu2O3', 'Pr6O11', 'Sm2O3', 'Tm2O3', 'Yb2O3'\\]",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "newprod": 0.00143 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.component
def test_REE_costing_variableOM_defaults():
    # use as many default values as possible, only pass required arguments

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "feed_input_rate")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert hasattr(m.fs.costing, "plant_overhead_cost")
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert hasattr(m.fs.costing, "land_cost")
    assert hasattr(m.fs.costing, "additional_chemicals_cost")
    assert hasattr(m.fs.costing, "additional_waste_cost")
    assert not hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.feed_input_rate)) == "ton/h"
    assert value(m.fs.costing.feed_input_rate) == pytest.approx(500.00, rel=1e-4)
    assert value(m.fs.costing.total_fixed_OM_cost) == pytest.approx(5.5384, rel=1e-4)
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.1388, rel=1e-4
    )
    assert value(m.fs.costing.plant_overhead_cost[0]) == pytest.approx(1.1077, rel=1e-4)
    assert value(m.fs.costing.other_variable_costs[0]) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.land_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.additional_chemicals_cost) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.additional_waste_cost) == pytest.approx(0.0000, abs=1e-4)


@pytest.mark.component
def test_REE_costing_variableOM_steadystateflowsheet():
    # use as many default values as possible, only pass required arguments

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "feed_input_rate")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert hasattr(m.fs.costing, "plant_overhead_cost")
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert hasattr(m.fs.costing, "land_cost")
    assert hasattr(m.fs.costing, "additional_chemicals_cost")
    assert hasattr(m.fs.costing, "additional_waste_cost")
    assert not hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.feed_input_rate)) == "ton/h"
    assert value(m.fs.costing.feed_input_rate) == pytest.approx(500.00, rel=1e-4)
    assert value(m.fs.costing.total_fixed_OM_cost) == pytest.approx(5.5384, rel=1e-4)
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.1388, rel=1e-4
    )
    assert value(m.fs.costing.plant_overhead_cost[0]) == pytest.approx(1.1077, rel=1e-4)
    assert value(m.fs.costing.other_variable_costs[0]) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.land_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.additional_chemicals_cost) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.additional_waste_cost) == pytest.approx(0.0000, abs=1e-4)


@pytest.mark.unit
def test_REE_costing_variableOM_nofixedOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        AttributeError,
        match="_ScalarQGESSCosting' object has no attribute 'hours_per_shift'",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_nofeedinput():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check some cost results
    assert not hasattr(m.fs.costing, "feed_input_rate")


@pytest.mark.component
def test_REE_costing_variableOM_feedinputnounits():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        UnitsError,
        match="The argument feed_input was passed as a dimensionless quantity "
        "with no units. Please ensure that the feed rate is passed in units "
        "of mass / time.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_resourcesnotalist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        TypeError,
        match="resources argument must be a list",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources={
                "water",
            },
            rates=[
                m.fs.water,
            ],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_ratesnotalist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        TypeError,
        match="rates argument must be a list",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates={
                m.fs.water,
            },
        )


@pytest.mark.component
def test_REE_costing_variableOM_customprices():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        prices={"water": 1.90e-3 * 1e-6 * pyunits.MUSD_2021 / pyunits.gallon},
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check some cost results
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.12301, rel=1e-4
    )


@pytest.mark.unit
def test_REE_costing_variableOM_custompricesnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        TypeError,
        match="prices argument must be a dictionary",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
            prices=[1e-3 * 1e-6 * pyunits.MUSD_2021 / pyunits.gallon],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_resourcesratesdifflengths():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        AttributeError,
        match="resources and rates must be lists of the same length",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=["water", "water"],
            rates=[
                m.fs.water,
            ],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_resourcenotinpricelist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.H2O = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.H2O.fix()

    with pytest.raises(
        AttributeError,
        match="A resource was included that does not contain a price. "
        "Prices exist for the following resources: \\['power', 'water', "
        "'diesel', 'bioleaching_solution', 'H2SO4', 'natural_gas', "
        "'polymer', 'NAOH', 'CACO3', 'coal_calcite', 'HCL', 'oxalic_acid',"
        " 'ascorbic_acid', 'kerosene', 'D2EHPA', 'NA2S', "
        "'nonhazardous_solid_waste', 'nonhazardous_precipitate_waste', "
        "'dust_and_volatiles'\\]",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "H2O",
            ],
            rates=[
                m.fs.H2O,
            ],
        )


# recovery rate
recovery_rate_units_dict = {
    "base_case": pyunits.kg / pyunits.year,
    "no_units": pyunits.dimensionless,
    "not_mass_units": pyunits.mol / pyunits.year,
    "not_per_year_units": pyunits.kg / pyunits.h,
}


recovery_rate_value_dict = {
    "base_case": 254324.44799999997,
    "no_units": 254324.44799999997,
    "not_mass_units": 254324.44799999997,
    "not_per_year_units": value(
        pyunits.convert(
            254324.44799999997 * pyunits.kg / pyunits.year,
            to_units=pyunits.kg / pyunits.h,
        )
    ),
}


@pytest.mark.parametrize(
    "recovery_rate_units, expectation",
    [
        ("base_case", does_not_raise()),
        (
            "no_units",
            pytest.raises(
                UnitsError,
                match="The argument recovery_rate_per_year was passed as a "
                "dimensionless quantity with no units. Please ensure that the "
                "feed rate is passed in units of mass / time.",
            ),
        ),
        (
            "not_mass_units",
            pytest.raises(
                UnitsError,
                match="The argument recovery_rate_per_year was passed with units of "
                "mol/a which cannot be converted to units of mass per year. Please "
                "ensure that recovery_rate_per_year is passed with rate units "
                "of mass per year \\(mass/a\\).",
            ),
        ),
        ("not_per_year_units", does_not_raise()),
    ],
)
@pytest.mark.component
def test_REE_costing_recovery(recovery_rate_units, expectation):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=recovery_rate_value_dict[recovery_rate_units],
        units=recovery_rate_units_dict[recovery_rate_units],
    )
    m.fs.recovery_rate_per_year.fix()

    with expectation:
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
            recovery_rate_per_year=m.fs.recovery_rate_per_year,
        )

        dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
        dt.assert_no_structural_warnings()

        QGESSCostingData.costing_initialization(m.fs.costing)
        QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
        QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
        solver = get_solver()
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)
        dt.assert_no_numerical_warnings()

        # check that some objects are built as expected
        assert hasattr(m.fs.costing, "recovery_rate_per_year")
        assert hasattr(m.fs.costing, "additional_cost_of_recovery")
        assert hasattr(m.fs.costing, "cost_of_recovery")

        # check some cost results
        assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
        assert value(m.fs.costing.recovery_rate_per_year) == pytest.approx(
            254324, rel=1e-4
        )
        assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
        assert value(m.fs.costing.cost_of_recovery) == pytest.approx(29.6178, rel=1e-4)
        assert value(m.fs.costing.additional_cost_of_recovery) == pytest.approx(
            0.0000, abs=1e-4
        )


@pytest.mark.component
def test_REE_costing_recovery_passedinmethodcall():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    # for convenience
    m.fs.annual_operating_hours = pyo.Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=False,
        units=pyunits.hours / pyunits.year,
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        recovery_rate_per_year=39.3
        * 0.8025
        * value(m.fs.annual_operating_hours)
        * pyunits.kg
        / pyunits.year,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert value(m.fs.costing.recovery_rate_per_year) == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert value(m.fs.costing.cost_of_recovery) == pytest.approx(29.6178, rel=1e-4)
    assert value(m.fs.costing.additional_cost_of_recovery) == pytest.approx(
        0.0000, abs=1e-4
    )


# transport cost
transport_cost_obj_dict = {
    "Expression_withunits": pyo.Expression(expr=10 * pyunits.USD_2021 / pyunits.ton),
    "Expression_nounits": pyo.Expression(expr=10),
    "Param_withunits": pyo.Param(
        initialize=10, units=pyunits.USD_2021 / pyunits.ton, mutable=False
    ),
    "Param_nounits": pyo.Param(initialize=10, mutable=False),
    "Var_withunits": pyo.Var(initialize=10, units=pyunits.USD_2021 / pyunits.ton),
    "Var_nounits": pyo.Var(initialize=10),
}


@pytest.mark.parametrize("transport_cost_obj", transport_cost_obj_dict.keys())
@pytest.mark.component
def test_REE_costing_recovery_transportcost(transport_cost_obj):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=39.3
        * 0.8025
        * 8
        * 3
        * 336,  # TREO (total rare earth oxide), 80.25% REE in REO
        units=pyunits.kg / pyunits.year,
    )
    m.fs.recovery_rate_per_year.fix()

    m.fs.transport_cost_per_ton_product = transport_cost_obj_dict[transport_cost_obj]
    if isinstance(m.fs.transport_cost_per_ton_product, pyo.Var):
        m.fs.transport_cost_per_ton_product.fix()  # pyo.Var must be fixed

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        recovery_rate_per_year=m.fs.recovery_rate_per_year,
        transport_cost_per_ton_product=m.fs.transport_cost_per_ton_product,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert value(m.fs.costing.recovery_rate_per_year) == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert value(m.fs.costing.cost_of_recovery) == pytest.approx(29.6178, rel=1e-4)
    assert value(m.fs.costing.additional_cost_of_recovery) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.unit
def test_REE_costing_recovery_Nonewithtransportcost():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        AttributeError,
        match="If transport_cost_per_ton_product is not None, "
        "recovery_rate_per_year cannot be None.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=True,
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            variable_OM=True,
            feed_input=m.fs.feed_input,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
            transport_cost_per_ton_product=10,
        )


@pytest.mark.unit
def test_REE_costing_config_defaults():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    assert isinstance(m.fs.costing.config, ConfigDict)

    assert m.fs.costing.config.discount_percentage is None
    assert m.fs.costing.config.plant_lifetime is None

    assert m.fs.costing.config.total_capital_cost is None
    assert m.fs.costing.config.annual_operating_cost is None
    assert m.fs.costing.config.annual_revenue is None
    assert m.fs.costing.config.cost_year is None

    assert m.fs.costing.config.has_capital_expenditure_period is False
    assert m.fs.costing.config.capital_expenditure_percentages is None

    assert m.fs.costing.config.capital_escalation_percentage == 3.6
    assert m.fs.costing.config.capital_loan_interest_percentage == 6
    assert m.fs.costing.config.capital_loan_repayment_period == 10
    assert m.fs.costing.config.debt_percentage_of_CAPEX == 50
    assert m.fs.costing.config.debt_expression is None
    assert m.fs.costing.config.operating_inflation_percentage == 3
    assert m.fs.costing.config.revenue_inflation_percentage == 3


@pytest.mark.unit
def test_REE_costing_config_kwargs():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.debt_expression = pyo.Var(initialize=1, units=pyunits.m)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[10, 60, 30],
        capital_escalation_percentage=2,
        capital_loan_interest_percentage=2,
        capital_loan_repayment_period=2,
        debt_percentage_of_CAPEX=2,
        debt_expression=m.fs.debt_expression,
        operating_inflation_percentage=2,
        revenue_inflation_percentage=2,
    )

    assert isinstance(m.fs.costing.config, ConfigDict)
    assert m.fs.costing.config.discount_percentage == 10
    assert m.fs.costing.config.plant_lifetime == 20
    assert m.fs.costing.config.total_capital_cost == 100
    assert m.fs.costing.config.annual_operating_cost == 100
    assert m.fs.costing.config.annual_revenue == 100
    assert m.fs.costing.config.cost_year == "2021"
    assert m.fs.costing.config.has_capital_expenditure_period is True
    assert m.fs.costing.config.capital_expenditure_percentages == [10, 60, 30]
    assert m.fs.costing.config.capital_escalation_percentage == 2
    assert m.fs.costing.config.capital_loan_interest_percentage == 2
    assert m.fs.costing.config.capital_loan_repayment_period == 2
    assert m.fs.costing.config.debt_percentage_of_CAPEX == 2
    assert isinstance(m.fs.costing.config.debt_expression, pyo.Var)
    assert m.fs.costing.config.operating_inflation_percentage == 2
    assert m.fs.costing.config.revenue_inflation_percentage == 2


@pytest.mark.unit
def test_REE_costing_NPV_costing_block_success():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        calculate_NPV=True,
    )


@pytest.mark.unit
def test_REE_costing_verify_NPV_costing_block_failure():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    with pytest.raises(
        AttributeError,
        match="Expected FlowsheetCostingBlockData object "
        "with attributes total_BEC, total_installation_cost, "
        "total_fixed_OM_cost, total_variable_OM_cost, "
        "other_plant_costs, land_cost, and total_sales_revenue. "
        "Please confirm that b is a FlowsheetCostingBlockData object "
        "and that all expected attributes exist.",
    ):
        QGESSCostingData.verify_calculate_from_costing_block(m.fs.costing)


# npv cost objects
npv_cost_obj_dict = {
    "Expression_withunits": pyo.Expression(expr=100 * pyunits.MUSD_2021),
    "Expression_nounits": pyo.Expression(expr=100),
    "Param_withunits": pyo.Param(
        initialize=100, units=pyunits.MUSD_2021, mutable=False
    ),
    "Param_nounits": pyo.Param(initialize=100, mutable=False),
    "Var_withunits": pyo.Var(initialize=100, units=pyunits.MUSD_2021),
    "Var_nounits": pyo.Var(initialize=100),
    "Scalar": 100,
}


@pytest.mark.parametrize("npv_cost_obj", npv_cost_obj_dict.keys())
@pytest.mark.unit
def test_REE_costing_verify_NPV_fixed_inputs_success(npv_cost_obj):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.obj = npv_cost_obj_dict[npv_cost_obj]

    if npv_cost_obj == "Scalar":
        # set directly
        m.fs.total_capital_cost = m.fs.obj
        m.fs.annual_operating_cost = m.fs.obj
        m.fs.annual_revenue = m.fs.obj
    else:
        # need to make References to original object
        m.fs.total_capital_cost = pyo.Reference(m.fs.obj)
        m.fs.annual_operating_cost = pyo.Reference(m.fs.obj)
        m.fs.annual_revenue = pyo.Reference(m.fs.obj)

    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=m.fs.total_capital_cost,
        annual_operating_cost=m.fs.annual_operating_cost,
        annual_revenue=m.fs.annual_revenue,
        cost_year="2021",
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )


@pytest.mark.unit
def test_REE_costing_verify_NPV_fixed_inputs_failure():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        total_capital_cost=100,
    )

    with pytest.raises(
        AttributeError,
        match="If capital, fixed O&M, or variable O&M costs are not calculated, "
        "then inputs for total_capital_cost, annual_operating_cost, and annual_revenue "
        "must be passed, and cost_year must be passed as a string, e.g. '2021'.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            calculate_NPV=True,
        )


@pytest.mark.unit
def test_REE_costing_verify_NPV_no_OM_useful_exception():
    # if users try to calculate NPV when fixed or variable OM don't exist,
    # suggest that they set fixed inputs
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    with pytest.raises(
        AttributeError,
        match="If capital, fixed O&M, or variable O&M costs are not calculated, "
        "then inputs for total_capital_cost, annual_operating_cost, and annual_revenue "
        "must be passed, and cost_year must be passed as a string, e.g. '2021'."
        "Alternatively, set fixed_OM and variable_OM to True to calculate O&M results.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            calculate_NPV=True,
        )


@pytest.mark.parametrize(
    "block_passed_as_None",
    [
        "total_capital_cost",
        "annual_operating_cost",
        "annual_revenue",
        "cost_year",
    ],
)
@pytest.mark.unit
def test_REE_costing_verify_NPV_fixed_input_failure_input_not_set(block_passed_as_None):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    for i in ["total_capital_cost", "annual_operating_cost", "annual_revenue"]:
        if i == block_passed_as_None:
            setattr(m.fs, i, None)
        else:
            setattr(m.fs, i, 100)

    if "cost_year" == block_passed_as_None:
        setattr(m.fs, "cost_year", None)
    else:
        setattr(m.fs, "cost_year", "2021")

    m.fs.costing = QGESSCosting(
        total_capital_cost=m.fs.total_capital_cost,
        annual_operating_cost=m.fs.annual_operating_cost,
        annual_revenue=m.fs.annual_revenue,
        cost_year=m.fs.cost_year,
    )

    with pytest.raises(
        AttributeError,
        match="If capital, fixed O&M, or variable O&M costs are not calculated, "
        "then inputs for total_capital_cost, annual_operating_cost, and annual_revenue "
        "must be passed, and cost_year must be passed as a string, e.g. '2021'.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            calculate_NPV=True,
        )


@pytest.mark.unit
def test_REE_costing_assert_has_config_argument_success():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(discount_percentage=10)

    QGESSCostingData.assert_config_argument_set(m.fs.costing, "discount_percentage")


@pytest.mark.unit
def test_REE_costing_assert_has_config_argument_failure_missing_argument():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    with pytest.raises(
        AttributeError,
        match="Required argument discount_percentage not set",
    ):
        QGESSCostingData.assert_config_argument_set(m.fs.costing, "discount_percentage")


@pytest.mark.unit
def test_REE_costing_assert_has_config_argument_failure_argument_not_set():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
    )

    with pytest.raises(
        AttributeError,
        match="Required argument discount_percentage not set",
    ):

        m.fs.costing.build_process_costs(
            fixed_OM=False,
            calculate_NPV=True,
        )


@pytest.mark.unit
def test_REE_costing_has_capital_expenditure_period_percentagesset():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[
            100,
        ],
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    assert m.fs.costing.config.capital_expenditure_percentages == [
        100,
    ]
    assert isinstance(m.fs.costing.capital_expenditure_percentages, pyo.Param)
    assert m.fs.costing.capital_expenditure_percentages.index_set() == [
        0,
    ]


@pytest.mark.component
def test_REE_costing_has_capital_expenditure_period_percentagesset_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[
            100,
        ],
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    def series_present_worth_factor(r, g, N):
        """
        Returns expression for series present worth factor.
        """
        return (1 - ((1 + g) ** (N)) * ((1 + r) ** (-N))) / (r - g)

    assert value(m.fs.costing.pv_capital_cost) == pytest.approx(
        value(
            -pyunits.convert(
                sum(
                    pyunits.convert(
                        m.fs.costing.config.capital_expenditure_percentages[idx]
                        * pyunits.percent,
                        to_units=pyunits.dimensionless,
                    )
                    * m.fs.costing.CAPEX
                    * (  # P/A_year(i) - P/A_year(i-1))
                        series_present_worth_factor(
                            pyunits.convert(
                                m.fs.costing.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            pyunits.convert(
                                m.fs.costing.capital_escalation_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            idx + 1,
                        )
                        - series_present_worth_factor(
                            pyunits.convert(
                                m.fs.costing.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            pyunits.convert(
                                m.fs.costing.capital_escalation_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            idx,
                        )
                    )
                    for idx in range(
                        len(m.fs.costing.config.capital_expenditure_percentages)
                    )
                ),
                to_units=m.fs.costing.cost_units,
            )
        ),
        rel=1e-4,
    )


@pytest.mark.unit
def test_REE_costing_has_capital_expenditure_period_percentagesnotset():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=True,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    assert m.fs.costing.config.capital_expenditure_percentages == [10, 60, 30]
    assert isinstance(m.fs.costing.capital_expenditure_percentages, pyo.Param)
    assert m.fs.costing.capital_expenditure_percentages.index_set() == [0, 1, 2]


@pytest.mark.component
def test_REE_costing_has_capital_expenditure_period_percentagesnotset_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=True,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    def series_present_worth_factor(r, g, N):
        """
        Returns expression for series present worth factor.
        """
        return (1 - ((1 + g) ** (N)) * ((1 + r) ** (-N))) / (r - g)

    assert value(m.fs.costing.pv_capital_cost) == pytest.approx(
        value(
            -pyunits.convert(
                sum(
                    pyunits.convert(
                        m.fs.costing.config.capital_expenditure_percentages[idx]
                        * pyunits.percent,
                        to_units=pyunits.dimensionless,
                    )
                    * m.fs.costing.CAPEX
                    * (  # P/A_year(i) - P/A_year(i-1))
                        series_present_worth_factor(
                            pyunits.convert(
                                m.fs.costing.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            pyunits.convert(
                                m.fs.costing.capital_escalation_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            idx + 1,
                        )
                        - series_present_worth_factor(
                            pyunits.convert(
                                m.fs.costing.discount_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            pyunits.convert(
                                m.fs.costing.capital_escalation_percentage,
                                to_units=pyunits.dimensionless,
                            ),
                            idx,
                        )
                    )
                    for idx in range(
                        len(m.fs.costing.config.capital_expenditure_percentages)
                    )
                ),
                to_units=m.fs.costing.cost_units,
            )
        ),
        rel=1e-4,
    )


@pytest.mark.unit
def test_REE_costing_not_has_capital_expenditure_period():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=False,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    assert m.fs.costing.config.capital_expenditure_percentages == []
    assert not hasattr(m.fs, "capital_expenditure_percentages")


@pytest.mark.component
def test_REE_costing_not_has_capital_expenditure_period_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=100,
        annual_operating_cost=100,
        annual_revenue=100,
        cost_year="2021",
        has_capital_expenditure_period=False,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    def series_present_worth_factor(r, g, N):
        """
        Returns expression for series present worth factor.
        """
        return (1 - ((1 + g) ** (N)) * ((1 + r) ** (-N))) / (r - g)

    assert value(m.fs.costing.pv_capital_cost) == pytest.approx(-100, rel=1e-4)


@pytest.mark.unit
def test_REE_costing_verify_percentages_list():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.capital_expenditure_percentages = [10, 60, 30]

    QGESSCostingData.verify_percentages_list(m.fs, "capital_expenditure_percentages")


@pytest.mark.unit
def test_REE_costing_verify_percentages_list_notalist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.capital_expenditure_percentages = 10

    with pytest.raises(
        TypeError,
        match="10 is not a list. "
        "Argument capital_expenditure_percentages must be passed as a list.",
    ):
        QGESSCostingData.verify_percentages_list(
            m.fs, "capital_expenditure_percentages"
        )


@pytest.mark.unit
def test_REE_costing_verify_percentages_list_zerolength():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.capital_expenditure_percentages = []

    with pytest.raises(
        AttributeError,
        match="Argument capital_expenditure_percentages has a length of "
        "zero. List must have a nonzero length.",
    ):
        QGESSCostingData.verify_percentages_list(
            m.fs, "capital_expenditure_percentages"
        )


@pytest.mark.unit
def test_REE_costing_verify_percentages_list_doesnotsumto100():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.capital_expenditure_percentages = [10, 60, 20]

    with pytest.raises(
        AttributeError,
        match="Argument capital_expenditure_percentages has a sum of 90. "
        "List must sum to 100 percent.",
    ):
        QGESSCostingData.verify_percentages_list(
            m.fs, "capital_expenditure_percentages"
        )


# assert cost objects
assert_cost_obj_dict = {
    "Param": pyo.Param(initialize=100, units=pyunits.MUSD_2021, mutable=False),
    "Var": pyo.Var(initialize=100, units=pyunits.MUSD_2021),
    "Expression": pyo.Expression(expr=100 * pyunits.MUSD_2021),
    "ScalarExpression": ScalarExpression(expr=100 * pyunits.MUSD_2021),
    "Scalar": 100,
    "None": None,
}


@pytest.mark.parametrize(
    "assert_cost_obj, expectation",
    [
        ("Param", does_not_raise()),
        ("Var", does_not_raise()),
        ("Expression", does_not_raise()),
        ("ScalarExpression", does_not_raise()),
        (
            "Scalar",
            pytest.raises(
                TypeError,
                match="Argument 100 of type <class 'int'> is not a supported object type. "
                "Ensure x is a Pyomo Param, Var, Expression, or ScalarExpression.",
            ),
        ),
        (
            "None",
            does_not_raise(),  # method assumes argument was not passed, since None is the default
        ),
    ],
)
@pytest.mark.unit
def test_REE_costing_assert_Pyomo_object(assert_cost_obj, expectation):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.x = assert_cost_obj_dict[assert_cost_obj]

    with expectation:
        QGESSCostingData.assert_Pyomo_object(m.fs, "x")


@pytest.mark.unit
def test_REE_costing_debt_expression():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.total_capital_cost = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.annual_operating_cost = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.annual_revenue = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.cost_year = "2021"

    # suppose debt is 50% of CAPEX (default), plus 5% of OPEX
    m.fs.debt_formula = pyo.Expression(
        expr=50 / 100 * m.fs.total_capital_cost + 5 / 100 * m.fs.annual_operating_cost
    )

    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=m.fs.total_capital_cost,
        annual_operating_cost=m.fs.annual_operating_cost,
        annual_revenue=m.fs.annual_revenue,
        cost_year=m.fs.cost_year,
        debt_expression=m.fs.debt_formula,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    assert isinstance(m.fs.costing.loan_debt[None], ScalarExpression)
    assert not hasattr(m.fs.costing, "loan_debt_constraint")
    assert m.fs.costing.loan_debt[None].expr == m.fs.debt_formula.expr


@pytest.mark.component
def test_REE_costing_debt_expression_solve():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.total_capital_cost = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.total_capital_cost.fix()
    m.fs.annual_operating_cost = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.annual_operating_cost.fix()
    m.fs.annual_revenue = pyo.Var(initialize=100, units=pyunits.MUSD_2021)
    m.fs.annual_revenue.fix()
    m.fs.cost_year = "2021"

    # suppose debt is 50% of CAPEX (default), plus 5% of OPEX
    m.fs.debt_formula = pyo.Expression(
        expr=50 / 100 * m.fs.total_capital_cost + 5 / 100 * m.fs.annual_operating_cost
    )

    m.fs.costing = QGESSCosting(
        discount_percentage=10,
        plant_lifetime=20,
        total_capital_cost=m.fs.total_capital_cost,
        annual_operating_cost=m.fs.annual_operating_cost,
        annual_revenue=m.fs.annual_revenue,
        cost_year=m.fs.cost_year,
        debt_expression=m.fs.debt_formula,
    )

    m.fs.costing.build_process_costs(
        fixed_OM=False,
        calculate_NPV=True,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    assert value(m.fs.costing.loan_debt[None]) == pytest.approx(
        value(m.fs.debt_formula), rel=1e-4
    )


@pytest.mark.component
def test_REE_costing_economy_of_numbers():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.costing = QGESSCosting()

    m.fs.base_cost = pyo.Var(initialize=41.2, units=pyunits.MUSD_2007)
    m.fs.base_cost.fix()

    QGESSCostingData.economy_of_numbers(
        m.fs.costing,
        cum_num_units=5,
        cost_FOAK=m.fs.base_cost,
        CE_index_year="2007",
        learning_rate=0.05,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    assert value(m.fs.costing.cost_NOAK) == pytest.approx(36.574, rel=1e-4)


@pytest.mark.component
def test_REE_costing_consider_taxes():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting(
        # arguments for NPV
        discount_percentage=10,  # percent
        plant_lifetime=20,  # years
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[10, 60, 30],
        capital_escalation_percentage=3.6,
        capital_loan_interest_percentage=6,
        capital_loan_repayment_period=10,
        debt_percentage_of_CAPEX=50,
        operating_inflation_percentage=3,
        revenue_inflation_percentage=3,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_REE_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "source": 1,
        },
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_process_costs(
        fixed_OM=True,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        variable_OM=True,
        feed_input=m.fs.feed_input,
        resources=[
            "water",
        ],
        rates=[
            m.fs.water,
        ],
        consider_taxes=True,
        calculate_NPV=True,
    )

    dt = DiagnosticsToolbox(model=m, variable_bounds_violation_tolerance=1e-4)
    dt.assert_no_structural_warnings()

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "feed_input_rate")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert hasattr(m.fs.costing, "plant_overhead_cost")
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert hasattr(m.fs.costing, "land_cost")
    assert hasattr(m.fs.costing, "additional_chemicals_cost")
    assert hasattr(m.fs.costing, "additional_waste_cost")
    assert not hasattr(m.fs.costing, "cost_of_recovery")

    assert hasattr(m.fs.costing, "income_tax_percentage")
    assert hasattr(m.fs.costing, "mineral_depletion_percentage")
    assert hasattr(m.fs.costing, "production_incentive_percentage")
    assert hasattr(m.fs.costing, "min_net_tax_owed")
    assert hasattr(m.fs.costing, "net_tax_owed")
    assert hasattr(m.fs.costing, "income_tax")
    assert hasattr(m.fs.costing, "additional_tax_credit")
    assert hasattr(m.fs.costing, "additional_tax_owed")
    assert hasattr(m.fs.costing, "pv_taxes")
    assert hasattr(m.fs.costing, "npv")

    assert isinstance(m.fs.costing.mineral_depletion_charge, pyo.Expression)
    assert isinstance(m.fs.costing.production_incentive_charge, pyo.Expression)
    assert isinstance(m.fs.costing.income_tax_eq, pyo.Constraint)
    assert isinstance(m.fs.costing.net_tax_owed_eq, pyo.Constraint)
    assert isinstance(m.fs.costing.pv_taxes_constraint, pyo.Constraint)
    assert isinstance(m.fs.costing.npv_constraint, pyo.Constraint)

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.feed_input_rate)) == "ton/h"
    assert value(m.fs.costing.feed_input_rate) == pytest.approx(500.00, rel=1e-4)
    assert value(m.fs.costing.total_fixed_OM_cost) == pytest.approx(5.5384, rel=1e-4)
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.1388, rel=1e-4
    )
    assert value(m.fs.costing.plant_overhead_cost[0]) == pytest.approx(1.1077, rel=1e-4)
    assert value(m.fs.costing.other_variable_costs[0]) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.land_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.additional_chemicals_cost) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.additional_waste_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.income_tax) == pytest.approx(5.303479, abs=1e-4)
    assert value(m.fs.costing.net_tax_owed) == pytest.approx(2.709606, abs=1e-4)
    assert value(m.fs.costing.pv_taxes) == pytest.approx(-17.33163, abs=1e-4)
    assert value(m.fs.costing.npv) == pytest.approx(158.08158, abs=1e-4)
