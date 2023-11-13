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

__author__ = "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, and M. Zamarripa)"
__version__ = "1.0.0"

import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import check_optimal_termination

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from ree_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)

# Create a Concrete Model as the top level object
m = pyo.ConcreteModel()

# Add a flowsheet object to the model
m.fs = FlowsheetBlock(dynamic=False)
m.fs.costing = QGESSCosting()
CE_index_year = "2016"

# check that the model solved properly and has 0 degrees of freedom
assert degrees_of_freedom(m) == 0

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
m.fs.CS_front_end_loader_2yd3.n_equip = pyo.Var(initialize=1,
                                                units=pyunits.dimensionless)
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
m.fs.CS_conveyors.throughput = pyo.Var(initialize=575,
                                       units=pyunits.tonne/pyunits.hr)
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
m.fs.DG_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.tonne)
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
m.fs.R_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.tonne)
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
m.fs.R_conveyors.throughput = pyo.Var(initialize=575,
                                      units=pyunits.tonne/pyunits.hr)
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
m.fs.R_roaster.duty = pyo.Var(initialize=737,
                              units=pyunits.MBTU/pyunits.hr)
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
m.fs.R_gas_scrubber.gas_rate = pyo.Var(initialize=11500,
                                    units=pyunits.ft**3/pyunits.min)
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
m.fs.R_spray_chamber_quencher.gas_rate = pyo.Var(initialize=11500,
                                                 units=pyunits.ft**3/pyunits.min)
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
m.fs.R_chiller.duty = pyo.Var(initialize=131,
                              units=pyunits.MBTU/pyunits.hr)
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
m.fs.L_pe_tanks.capacity = pyo.Var(initialize=164805,
                                units=pyunits.gal)
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
m.fs.L_tank_mixers.power = pyo.Var(initialize=474,
                                    units=pyunits.hp)
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
m.fs.L_pump.feed_rate = pyo.Var(initialize=10987,
                                    units=pyunits.gal/pyunits.min)
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
m.fs.L_thickener.area = pyo.Var(initialize=22590,
                                units=pyunits.ft**2)
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
m.fs.L_filter_press.volume = pyo.Var(initialize=3600,
                                units=pyunits.ft**3)
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
m.fs.L_solution_heater.duty = pyo.Var(initialize=2.4,
                                      units=pyunits.MBTU/pyunits.hr)
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
m.fs.RSX_pe_tanks.capacity = pyo.Var(initialize=35136,
                                     units=pyunits.gal)
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
m.fs.RSX_tank_mixers.power = pyo.Var(initialize=20,
                                     units=pyunits.hp)
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
m.fs.RSX_pump.feed_rate = pyo.Var(initialize=7027,
                                  units=pyunits.gal/pyunits.min)
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
m.fs.RSX_mixer_settler.volume = pyo.Var(initialize=61107,
                                        units=pyunits.gal)
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
m.fs.CSX_pe_tanks.capacity = pyo.Var(initialize=1405,
                                     units=pyunits.gal)
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
m.fs.CSX_tank_mixers.power = pyo.Var(initialize=0.8,
                                     units=pyunits.hp)
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
m.fs.CSX_pump.feed_rate = pyo.Var(initialize=281,
                                  units=pyunits.gal/pyunits.min)
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
m.fs.CSX_mixer_settler.volume = pyo.Var(initialize=2444,
                                        units=pyunits.gal)
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
m.fs.SX_wash_pe_tanks.capacity = pyo.Var(initialize=3514,
                                         units=pyunits.gal)
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
m.fs.SX_wash_tank_mixers.power = pyo.Var(initialize=2,
                                         units=pyunits.hp)
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
m.fs.SX_wash_pump.feed_rate = pyo.Var(initialize=703,
                                      units=pyunits.gal/pyunits.min)
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
m.fs.SX_wash_mixer_settler.volume = pyo.Var(initialize=18332,
                                            units=pyunits.gal)
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
m.fs.SX_wash_filter_press.volume = pyo.Var(initialize=0.26,
                                           units=pyunits.ft**3)
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
m.fs.reep_pe_tanks.capacity = pyo.Var(initialize=1405,
                                      units=pyunits.gal)
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
m.fs.reep_tank_mixers.power = pyo.Var(initialize=0.61,
                                      units=pyunits.hp)
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
m.fs.reep_pump.feed_rate = pyo.Var(initialize=70,
                                   units=pyunits.gal/pyunits.min)
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
m.fs.reep_filter_press.volume = pyo.Var(initialize=0.405,
                                        units=pyunits.ft**3)
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
m.fs.reep_roaster.duty = pyo.Var(initialize=0.35,
                              units=pyunits.MBTU/pyunits.hr)
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
m.fs.WT_pe_tanks.capacity = pyo.Var(initialize=453131,
                                    units=pyunits.gal)
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
m.fs.WT_pump.feed_rate = pyo.Var(initialize=78805,
                                 units=pyunits.gal/pyunits.min)
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
m.fs.WT_filter_press.volume = pyo.Var(initialize=469,
                                        units=pyunits.ft**3)
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
m.fs.WT_conveyors.throughput = pyo.Var(initialize=569,
                                      units=pyunits.tonne/pyunits.hr)
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

# add plant-level cost constraints
m.fs.feed_input = pyo.Var(initialize=500, units=pyunits.tonne/pyunits.hr)
m.fs.feed_grade = pyo.Var(initialize=329, units=pyunits.ppm)
tonnes_REE_capture = pyo.value(pyunits.convert(m.fs.feed_input, to_units=pyunits.tonne/pyunits.hr)
                           * 8 * 3 * 336 * pyunits.hr  # 8-hr shifts, 3 shifts/day, 336 days/year
                           * pyunits.convert(m.fs.feed_grade, to_units=pyunits.dimensionless)
                           * 0.10  # 10% total REE recovery
                           )
m.fs.land_cost_eq = pyo.Expression(expr=0.30e-6 * m.fs.feed_input
                               * 8 * 3 * 336 * pyunits.hr)  # 0.30 USD/tonne per year, in MUSD

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
        nameplate_capacity=500,  # tonne/hr
        capacity_factor=0.92,  # % of maximum capacity that is utilized
        labor_types = ["skilled", "unskilled", "supervisor", "maintenance", "technician", "engineer"],
        labor_rate=[27.90, 23.26, 30.29, 24.06, 23.43, 46.82],  # USD/hr
        labor_burden=25,  # % fringe benefits
        operators_per_shift=[2, 5, 2, 3, 1, 2],
        hours_per_shift=8,
        shifts_per_day=3,
        operating_days_per_year=336,
        # arguments related to total owners costs
        land_cost=m.fs.land_cost_eq,
        fixed_OM=True,
        tonne_REE_capture=tonnes_REE_capture,
        CE_index_year=CE_index_year,
    )

# add initialize
QGESSCostingData.costing_initialization(m.fs.costing)

# try solving
solver = get_solver()
results = solver.solve(m, tee=True)
assert check_optimal_termination(results)

# check unit consistency
assert_units_consistent(m)

# report results
QGESSCostingData.report(m.fs.costing)
QGESSCostingData.display_bare_erected_costs(m.fs.costing)
QGESSCostingData.display_flowsheet_cost(m.fs.costing)

# test costing bounding method
QGESSCostingData.calculate_REE_costing_bounds(
    b=m.fs.costing,
    capacity=m.fs.feed_input * 8064 * pyunits.hr/pyunits.a * 20 * pyunits.a,
    grade=m.fs.feed_grade,
    CE_index_year=CE_index_year,
    recalculate=True,
    )
