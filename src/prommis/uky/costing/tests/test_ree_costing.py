#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for REE costing.

"""

import pyomo.environ as pyo
from pyomo.environ import check_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.uky.costing.ree_plant_capcost import QGESSCosting, QGESSCostingData


# fixture so other tests don't need to explicitly re-build unit blocks
@pytest.fixture(scope="module")
def m():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()
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


@pytest.mark.component
def test_REE_costing(m):
    # full smoke test with all components, O&M costs, and extra costs included
    CE_index_year = "UKy_2019"

    # add plant-level cost constraints

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)

    hours_per_shift = 8
    shifts_per_day = 3
    operating_days_per_year = 336

    # for convenience
    m.fs.annual_operating_hours = pyo.Param(
        initialize=hours_per_shift * shifts_per_day * operating_days_per_year,
        mutable=False,
        units=pyunits.hours / pyunits.a,
    )

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=39.3
        * 0.8025
        * pyo.value(
            m.fs.annual_operating_hours
        ),  # TREO (total rare earth oxide), 80.25% REE in REO
        units=pyunits.kg / pyunits.year,
    )

    # the land cost is the lease cost, or refining cost of REO produced
    m.fs.land_cost = pyo.Expression(
        expr=0.303736
        * 1e-6
        * getattr(pyunits, "MUSD_" + CE_index_year)
        / pyunits.ton
        * pyunits.convert(m.fs.feed_input, to_units=pyunits.ton / pyunits.hr)
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
        / pyunits.a
    )

    m.fs.reagents = pyo.Var(
        m.fs.time,
        initialize=reagent_costs / (m.fs.annual_operating_hours),
        units=pyunits.kg / pyunits.hr,
    )

    m.fs.solid_waste = pyo.Var(
        m.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
    )  # non-hazardous solid waste
    m.fs.precipitate = pyo.Var(
        m.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
    )  # non-hazardous precipitate
    m.fs.dust_and_volatiles = pyo.Var(
        m.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
    )  # dust and volatiles
    m.fs.power = pyo.Var(m.fs.time, initialize=14716, units=pyunits.hp)

    resources = [
        "dummy",
        "nonhazardous_solid_waste",
        "nonhazardous_precipitate_waste",
        "dust_and_volatiles",
        "power",
    ]

    rates = [
        m.fs.reagents,
        m.fs.solid_waste,
        m.fs.precipitate,
        m.fs.dust_and_volatiles,
        m.fs.power,
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

    m.fs.costing.build_process_costs(
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
        nameplate_capacity=500,  # short (US) ton/hr
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
        land_cost=m.fs.land_cost,
        resources=resources,
        rates=rates,
        prices={
            "dummy": 1 * getattr(pyunits, "USD_" + CE_index_year) / pyunits.kg,
        },
        fixed_OM=True,
        variable_OM=True,
        feed_input=m.fs.feed_input,
        efficiency=0.80,  # power usage efficiency, or fixed motor/distribution efficiency
        chemicals=["dummy"],
        waste=[
            "nonhazardous_solid_waste",
            "nonhazardous_precipitate_waste",
            "dust_and_volatiles",
        ],
        recovery_rate_per_year=m.fs.recovery_rate_per_year,
        CE_index_year=CE_index_year,
    )

    # define reagent fill costs as an other plant cost so framework adds this to TPC calculation
    m.fs.costing.other_plant_costs.unfix()
    m.fs.costing.other_plant_costs_rule = pyo.Constraint(
        expr=(
            m.fs.costing.other_plant_costs
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
    m.fs.feed_input.fix()
    m.fs.feed_grade.fix()
    m.fs.recovery_rate_per_year.fix()
    m.fs.reagents.fix()
    m.fs.solid_waste.fix()
    m.fs.precipitate.fix()
    m.fs.dust_and_volatiles.fix()
    m.fs.power.fix()

    # check that the model is set up properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0


@pytest.mark.component
def test_initialize(m):
    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)


@pytest.mark.component
def test_solve(m):
    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)


@pytest.mark.component
def test_results(m):
    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_UKy_2019
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    assert hasattr(m.fs.costing, "land_cost")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert hasattr(m.fs.costing, "plant_overhead_cost")
    # if piping materials and labor exists, all other factors exist as well
    assert hasattr(m.fs.costing, "piping_materials_and_labor_percentage")
    # using percentages, Lang factor should not be built or used
    assert not hasattr(m.fs.costing, "Lang_factor")

    # check some overall cost results
    assert m.fs.costing.total_plant_cost.value == pytest.approx(133.23, rel=1e-4)
    assert m.fs.costing.total_BEC.value == pytest.approx(44.308, rel=1e-4)
    assert m.fs.costing.total_installation_cost.value == pytest.approx(87.287, rel=1e-4)
    assert m.fs.costing.other_plant_costs.value == pytest.approx(1.6309, rel=1e-4)
    assert m.fs.costing.total_fixed_OM_cost.value == pytest.approx(11.105, rel=1e-4)
    assert m.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        532.90, rel=1e-4
    )
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(1.2247, rel=1e-4)
    assert m.fs.costing.total_sales_revenue.value == pytest.approx(65.333, rel=1e-4)
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(2199.1, rel=1e-4)


@pytest.mark.component
def test_units_consistency(m):
    # check unit consistency
    assert_units_consistent(m)


@pytest.mark.component
def test_report(m):
    # test report methods
    QGESSCostingData.report(m.fs.costing)
    m.fs.costing.variable_operating_costs.display()  # results will be in t = 0
    print()
    QGESSCostingData.display_bare_erected_costs(m.fs.costing)
    QGESSCostingData.display_flowsheet_cost(m.fs.costing)


@pytest.mark.component
def test_costing_bounding(m):
    # test costing bounding method
    CE_index_year = "UKy_2019"
    QGESSCostingData.calculate_REE_costing_bounds(
        b=m.fs.costing,
        capacity=m.fs.feed_input * m.fs.annual_operating_hours * 20 * pyunits.a,
        grade=m.fs.feed_grade,
        CE_index_year=CE_index_year,
        recalculate=True,
    )

    m.fs.costing.costing_lower_bound.pprint()
    m.fs.costing.costing_upper_bound.pprint()

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

    for key in m.fs.costing.costing_lower_bound.keys():
        assert m.fs.costing.costing_lower_bound[key].value == pytest.approx(
            expected_costing_lower_bound[key], rel=1e-4
        )

    for key in m.fs.costing.costing_upper_bound.keys():
        assert m.fs.costing.costing_upper_bound[key].value == pytest.approx(
            expected_costing_upper_bound[key], rel=1e-4
        )


@pytest.mark.component
def test_REE_costing_noOM_usedefaults():
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
        fixed_OM=False,
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")

    assert hasattr(m.fs.costing, "land_cost")
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(0, abs=1e-4)
    assert pyunits.get_units(m.fs.costing.land_cost) == pyunits.MUSD_2021

    assert not hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert not hasattr(m.fs.costing, "cost_of_recovery")
    assert not hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert not hasattr(m.fs.costing, "total_variable_OM_cost")
    assert not hasattr(m.fs.costing, "plant_overhead_cost")

    # if piping materials and labor exists, all other factors exist as well
    assert hasattr(m.fs.costing, "piping_materials_and_labor_percentage")
    # using percentages, Lang factor should not be built or used
    assert not hasattr(m.fs.costing, "Lang_factor")

    assert m.fs.costing.total_plant_cost.value == pytest.approx(7.4612, rel=1e-4)


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

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    # check that total plant cost ratios match expected currency conversions
    for account in CS_jaw_crusher_accounts:
        # USD_2021 = CE 708.0
        assert (
            pytest.approx(
                pyo.value(
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
                pyo.value(
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
                pyo.value(
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
                pyo.value(
                    m.fs.CS_jaw_crusher_2025.costing.bare_erected_cost[account]
                    / m.fs.CS_jaw_crusher_CE500.costing.bare_erected_cost[account]
                ),
                abs=1e-1,
            )
            == 815.59 / 500
        )


@pytest.mark.component
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


@pytest.mark.component
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

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_units_consistent(m)

    assert m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3"].value == pytest.approx(
        2.5122, rel=1e-4
    )
    assert m.fs.CS_crusher.costing.bare_erected_cost["1.3"].value == pytest.approx(
        2.5122, rel=1e-4
    )
    assert m.fs.CS_roll_crusher.costing.bare_erected_cost["1.5"].value == pytest.approx(
        0.32769, rel=1e-4
    )
    assert m.fs.CS_crusher.costing.bare_erected_cost["1.5"].value == pytest.approx(
        0.32769, rel=1e-4
    )


@pytest.mark.component
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

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_units_consistent(m)

    # adding a check just to make sure everything works as expected
    assert m.fs.CS_jaw_crusher.costing.bare_erected_cost[
        "1.3new"
    ].value == pytest.approx(2.5122, rel=1e-4)


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

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_units_consistent(m)

    # adding a check just to make sure 1.3 was overwritten before it was used
    assert m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3"].value == pytest.approx(
        2.0650, rel=1e-4
    )


@pytest.mark.component
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

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    # base case
    assert m.fs.CS_jaw_crusher_1.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(1.0000, rel=1e-4)

    # only one unit, parallel doesn't change result
    assert m.fs.CS_jaw_crusher_2.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(1.0000, rel=1e-4)

    # same capacity over two units should be slightly less expensive than one unit
    assert m.fs.CS_jaw_crusher_3.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(0.84095, rel=1e-4)

    # two units with base case capacity should be double the cost
    assert m.fs.CS_jaw_crusher_4.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(2.0000, rel=1e-4)

    # same capacity over five units should be much less expensive than one unit
    assert m.fs.CS_jaw_crusher_5.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(0.66878, rel=1e-4)

    # five units with base case capacity should be five times the cost
    assert m.fs.CS_jaw_crusher_6.costing.bare_erected_cost[
        "1.3"
    ].value == pytest.approx(5.0000, rel=1e-4)


@pytest.mark.component
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
def test_REE_costing_usersetTPC():
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

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert m.fs.costing.total_plant_cost.get_units() == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")  # built using passed TPC
    assert not hasattr(m.fs.costing, "total_BEC_eq")  # shouldn't be built
    assert hasattr(m.fs.costing, "total_installation_cost")
    assert hasattr(m.fs.costing, "total_plant_cost")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert m.fs.costing.total_BEC.value == pytest.approx(1.0000, rel=1e-4)
    assert m.fs.costing.total_plant_cost.value == pytest.approx(2.9700, rel=1e-4)


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

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

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

    assert m.fs.costing.total_BEC.value == pytest.approx(2.5122, rel=1e-4)
    assert m.fs.costing.total_plant_cost.value == pytest.approx(7.4612, rel=1e-4)


@pytest.mark.component
def test_REE_costing_landcostExpression_withunits():
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

    # define land_cost
    m.fs.land_cost = pyo.Expression(expr=1e-3 * pyunits.MUSD_2021)

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        land_cost=m.fs.land_cost,
        fixed_OM=False,
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert hasattr(m.fs.costing, "land_cost")
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(1e-3, rel=1e-4)
    assert pyunits.get_units(m.fs.costing.land_cost) == pyunits.MUSD_2021


@pytest.mark.component
def test_REE_costing_landcostExpression_nounits():
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

    # define land_cost
    m.fs.land_cost = pyo.Expression(expr=1e-3)

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        land_cost=m.fs.land_cost,
        fixed_OM=False,
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert hasattr(m.fs.costing, "land_cost")
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(1e-3, rel=1e-4)
    assert pyunits.get_units(m.fs.costing.land_cost) == pyunits.MUSD_2021


@pytest.mark.component
def test_REE_costing_landcostnonExpression_withunits():
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

    # define land_cost
    m.fs.land_cost = pyo.Var(initialize=1e-3, units=pyunits.MUSD_2021)

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        land_cost=m.fs.land_cost,
        fixed_OM=False,
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert hasattr(m.fs.costing, "land_cost")
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(1e-3, rel=1e-4)
    assert pyunits.get_units(m.fs.costing.land_cost) == pyunits.MUSD_2021


@pytest.mark.component
def test_REE_costing_landcostnonExpression_nounits():
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

    # define land_cost
    m.fs.land_cost = pyo.Var(initialize=1e-3)

    # defaults to fixed_OM=True, so explicitly set to False
    # defaults to variable_OM=False, so let that use the default
    m.fs.costing.build_process_costs(
        land_cost=m.fs.land_cost,
        fixed_OM=False,
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that the cost units are as expected
    assert hasattr(m.fs.costing, "land_cost")
    assert pyo.value(m.fs.costing.land_cost) == pytest.approx(1e-3, rel=1e-4)
    assert pyunits.get_units(m.fs.costing.land_cost) == pyunits.MUSD_2021


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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

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

    assert m.fs.costing.annual_operating_labor_cost.value == pytest.approx(
        3.0730, rel=1e-4
    )
    assert m.fs.costing.annual_technical_labor_cost.value == pytest.approx(
        1.1801, rel=1e-4
    )
    assert m.fs.costing.annual_labor_cost.value == pytest.approx(4.2531, rel=1e-4)
    assert m.fs.costing.maintenance_and_material_cost.value == pytest.approx(
        0.14922, rel=1e-4
    )
    assert m.fs.costing.quality_assurance_and_control_cost.value == pytest.approx(
        0.30730, rel=1e-4
    )
    assert m.fs.costing.sales_patenting_and_research_cost.value == pytest.approx(
        0.32191, rel=1e-4
    )
    assert m.fs.costing.admin_and_support_labor_cost.value == pytest.approx(
        0.61460, rel=1e-4
    )
    assert m.fs.costing.property_taxes_and_insurance_cost.value == pytest.approx(
        0.074612, rel=1e-4
    )
    assert m.fs.costing.other_fixed_costs.value == pytest.approx(0.0000, abs=1e-4)
    assert m.fs.costing.total_fixed_OM_cost.value == pytest.approx(5.7207, rel=1e-4)
    assert m.fs.costing.total_sales_revenue.value == pytest.approx(64.382, rel=1e-4)


@pytest.mark.component
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
        match="Cannot add component 'labor_rate_index' \\(type \\<class "
        "'pyomo.core.base.set.OrderedScalarSet'\\>\\) to block "
        "'fs.costing': a component by that name \\(type \\<class "
        "'pyomo.core.base.set.OrderedScalarSet'\\>\\) is already "
        "defined.",
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


@pytest.mark.component
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


@pytest.mark.component
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


@pytest.mark.component
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


@pytest.mark.component
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


@pytest.mark.component
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


@pytest.mark.component
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


@pytest.mark.component
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
        Exception,
        match="A pure product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Sc', 'Y', 'La', "
        "'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', "
        "'Yb', 'Lu', 'Sc2O3', 'Y2O3', 'La2O3', 'CeO2', 'Pr6O11', 'Nd2O3', "
        "'Sm2O3', 'Eu2O3', 'Gd2O3', 'Tb4O7', 'Dy2O3', 'Ho2O3', 'Er2O3', "
        "'Tm2O3', 'Yb2O3', 'Lu2O3'\\]",
    ):
        m.fs.costing.build_process_costs(
            pure_product_output_rates={
                "newprod": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.component
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


@pytest.mark.component
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
        Exception,
        match="A mixed product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Sc', 'Y', 'La', "
        "'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', "
        "'Yb', 'Lu', 'Sc2O3', 'Y2O3', 'La2O3', 'CeO2', 'Pr6O11', 'Nd2O3', "
        "'Sm2O3', 'Eu2O3', 'Gd2O3', 'Tb4O7', 'Dy2O3', 'Ho2O3', 'Er2O3', "
        "'Tm2O3', 'Yb2O3', 'Lu2O3'\\]",
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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "feed_input_rate")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert hasattr(m.fs.costing, "plant_overhead_cost")
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert not hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.feed_input_rate)) == "ton/h"
    assert pyo.value(m.fs.costing.feed_input_rate) == pytest.approx(500.00, rel=1e-4)
    assert m.fs.costing.total_fixed_OM_cost.value == pytest.approx(5.7207, rel=1e-4)
    assert m.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        1.1595, rel=1e-4
    )
    assert m.fs.costing.plant_overhead_cost[0].value == pytest.approx(1.1441, rel=1e-4)
    assert m.fs.costing.other_variable_costs[0].value == pytest.approx(0.0000, abs=1e-4)


@pytest.mark.component
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


@pytest.mark.component
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

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        AttributeError,
        match="No feed_input rate variable passed to main costing block.",
    ):
        m.fs.costing.build_process_costs(
            fixed_OM=False,
            variable_OM=True,
            resources=[
                "water",
            ],
            rates=[
                m.fs.water,
            ],
        )


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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.feed_input_rate)) == "ton/h"
    assert pyo.value(m.fs.costing.feed_input_rate) == pytest.approx(500.00, rel=1e-4)


@pytest.mark.component
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


@pytest.mark.component
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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check some cost results
    assert m.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        1.1595, rel=1e-4
    )


@pytest.mark.component
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


@pytest.mark.component
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
        Exception,
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


@pytest.mark.component
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
        Exception,
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


@pytest.mark.component
def test_REE_costing_recovery_basecase():
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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )


@pytest.mark.component
def test_REE_costing_recovery_nounits():
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
        units=pyunits.hours / pyunits.a,
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=39.3
        * 0.8025
        * pyo.value(
            m.fs.annual_operating_hours
        ),  # TREO (total rare earth oxide), 80.25% REE in REO,
        # no units passed - should be parsed as kg/a (multipliers give numerical rate per year)
    )
    m.fs.recovery_rate_per_year.fix()

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )


@pytest.mark.component
def test_REE_costing_recovery_notmassunits():
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
        units=pyunits.hours / pyunits.a,
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=39.3
        * 0.8025
        * pyo.value(
            m.fs.annual_operating_hours
        ),  # TREO (total rare earth oxide), 80.25% REE in REO,
        units=pyunits.mol / pyunits.year,
    )
    m.fs.recovery_rate_per_year.fix()

    with pytest.raises(
        Exception,
        match="The argument recovery_rate_per_year was passed with units of "
        "mol/a which cannot be converted to units of mass per year. Please "
        "ensure that recovery_rate_per_year is passed with rate units "
        "of mass per year \\(mass/a\\) or dimensionless.",
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
            recovery_rate_per_year=m.fs.recovery_rate_per_year,
        )


@pytest.mark.component
def test_REE_costing_recovery_notannualbasisunits():
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
        units=pyunits.hours / pyunits.a,
    )

    m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feed_input.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.recovery_rate_per_year = pyo.Var(
        initialize=39.3
        * 0.8025
        * pyo.value(
            m.fs.annual_operating_hours
        ),  # TREO (total rare earth oxide), 80.25% REE in REO,
        units=pyunits.kg / pyunits.h,
    )
    m.fs.recovery_rate_per_year.fix()

    with pytest.raises(
        Exception,
        match="The argument recovery_rate_per_year was passed with units of "
        "kg/h and must be on an anuual basis. Please "
        "ensure that recovery_rate_per_year is passed with rate units "
        "of mass per year \\(mass/a\\) or dimensionless.",
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
            recovery_rate_per_year=m.fs.recovery_rate_per_year,
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
        units=pyunits.hours / pyunits.a,
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
        recovery_rate_per_year=39.3 * 0.8025 * pyo.value(m.fs.annual_operating_hours),
    )

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )


@pytest.mark.component
def test_REE_costing_recovery_transportcostExpression():
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

    m.fs.transport_cost_per_ton_product = pyo.Expression(
        expr=10 * pyunits.USD_2021 / pyunits.ton
    )

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
def test_REE_costing_recovery_transportcostExpressionnounits():
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

    m.fs.transport_cost_per_ton_product = pyo.Expression(expr=10)

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
def test_REE_costing_recovery_transportcostParam():
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

    m.fs.transport_cost_per_ton_product = pyo.Param(
        initialize=10, units=pyunits.USD_2021 / pyunits.ton, mutable=False
    )

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
def test_REE_costing_recovery_transportcostParamnounits():
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

    m.fs.transport_cost_per_ton_product = pyo.Param(initialize=10, mutable=False)

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
def test_REE_costing_recovery_transportcostVar():
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

    m.fs.transport_cost_per_ton_product = pyo.Var(
        initialize=10, units=pyunits.USD_2021 / pyunits.ton
    )
    m.fs.transport_cost_per_ton_product.fix()

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
def test_REE_costing_recovery_transportcostVarnounits():
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

    m.fs.transport_cost_per_ton_product = pyo.Var(initialize=10)
    m.fs.transport_cost_per_ton_product.fix()

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

    QGESSCostingData.costing_initialization(m.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert check_optimal_termination(results)
    assert_units_consistent(m)

    # check that some objects builts as expected
    assert hasattr(m.fs.costing, "recovery_rate_per_year")
    assert hasattr(m.fs.costing, "additional_cost_of_recovery")
    assert hasattr(m.fs.costing, "cost_of_recovery")
    assert hasattr(m.fs.costing, "transport_cost")

    # check some cost results
    assert str(pyunits.get_units(m.fs.costing.recovery_rate_per_year)) == "kg/a"
    assert m.fs.costing.recovery_rate_per_year.value == pytest.approx(254324, rel=1e-4)
    assert str(pyunits.get_units(m.fs.costing.cost_of_recovery)) == "USD_2021/kg"
    assert pyo.value(m.fs.costing.cost_of_recovery) == pytest.approx(30.416, rel=1e-4)
    assert m.fs.costing.additional_cost_of_recovery.value == pytest.approx(
        0.0000, abs=1e-4
    )
    assert pyo.value(m.fs.costing.transport_cost) == pytest.approx(0.0028034, rel=1e-4)


@pytest.mark.component
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
        Exception,
        match="If a transport_cost_per_ton_product is not None, "
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
