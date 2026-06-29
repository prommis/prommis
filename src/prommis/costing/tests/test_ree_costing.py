#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for REE costing.

"""

import os

import pyomo.environ as pyo
from pyomo.common.dependencies import attempt_import
from pyomo.environ import assert_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.environ import value

import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    get_jacobian,
    jacobian_cond,
)

import pytest

from prommis.costing.custom_costing_example import CustomCostingData
from prommis.costing.ree_costing import (
    REECosting,
    REECostingData,
    REEUnitModelCostingBlock,
)
from prommis.costing.ree_costing_dictionaries import (
    load_REE_costing_dictionary,
    register_ree_currency_units,
)
from prommis.nanofiltration.costing.diafiltration_cost_model import (
    DiafiltrationCostingData,
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
def test_no_config_set():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # QGESSCosting that REECosting inherits from requires tech to be set explicitly
    # but REECosting modifies to set the default to 10 (UKy)
    assert m.fs.costing.config.tech == 10


@pytest.mark.component
def test_REEUnitModelCostingBlock_defaults():
    REE_costing_params = load_REE_costing_dictionary()  # UKy

    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()
    m.fs.unit = UnitModelBlock()
    m.fs.unit.scaled_var = pyo.Var(initialize=1)
    m.fs.unit.scaled_var.fix()

    m.fs.unit.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": [
                "1.1",
            ],
            "scaled_param": m.fs.unit.scaled_var,
        },
    )

    assert m.fs.unit.costing.costing_method_arguments["cost_accounts"] == ["1.1"]
    assert (
        m.fs.unit.costing.costing_method_arguments["scaled_param"]
        == m.fs.unit.scaled_var
    )

    # defaults set by REEUnitModelCostingBlock() wrapper on UnitModelCostingBlock
    assert m.fs.unit.costing.costing_method_arguments["tech"] == 10
    assert m.fs.unit.costing.costing_method_arguments["ccs"] == "A"
    assert m.fs.unit.costing.costing_method_arguments["additional_costing_params"] == [
        REE_costing_params,
    ]
    assert m.fs.unit.costing.costing_method_arguments["CEPCI_year"] == "2021"

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    # indirect test that it worked as expected - default year is 2021
    assert value(m.fs.unit.costing.bare_erected_cost["1.1"]) == pytest.approx(
        value(
            pyunits.convert(
                147400.0 * pyunits.USD_2016,  # reference value
                to_units=pyunits.MUSD_2021,  # default CEPCI plant units
            )
        )
    )
    assert (
        pyunits.get_units(m.fs.unit.costing.bare_erected_cost["1.1"])
        == pyunits.MUSD_2021
    )


@pytest.mark.unit
def test_register_REE_currency_units_twice(caplog):
    # check that units exist - they are registered when REECosting imports
    assert hasattr(pyunits, "USD_UKy_2019")
    assert hasattr(pyunits, "USD_2025")

    # register units again
    register_ree_currency_units()
    msg = (
        "Custom REE plant currency units (USD_2025, USD_UKy_2019) "
        "already appear in Pyomo unit registry. Assuming repeated call of "
        "register_ree_currency_units."
    )
    for record in caplog.records:
        assert msg in record.message


def base_model():

    CEPCI_year = "UKy_2019"

    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REECosting(
        Lang_factor=2.97,
        has_fixed_OM=True,
        has_variable_OM=True,
        has_net_present_value=True,
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[10, 60, 30],
        CEPCI_year=CEPCI_year,
    )

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

    m.fs.CS_front_end_loader_2yd3.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_front_end_loader_2yd3_accounts,
            "scaled_param": m.fs.CS_front_end_loader_2yd3.n_equip,  # 1 loader
            # no. units is the scaling parameter for constant-cost units,
            #     so use n_equip below to specify the number of loaders
            "n_equip": 5,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 1.5 is CS Roll Crusher
    CS_roll_crusher_accounts = ["1.5"]
    m.fs.CS_roll_crusher = UnitModelBlock()
    m.fs.CS_roll_crusher.power = pyo.Var(initialize=430, units=pyunits.hp)
    m.fs.CS_roll_crusher.power.fix()
    m.fs.CS_roll_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_roll_crusher_accounts,
            "scaled_param": m.fs.CS_roll_crusher.power,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 1.6 is CS Vibrating Screens
    CS_vibrating_screen_accounts = ["1.6"]
    m.fs.CS_vibrating_screens = UnitModelBlock()
    m.fs.CS_vibrating_screens.area = pyo.Var(initialize=124, units=pyunits.ft**2)
    m.fs.CS_vibrating_screens.area.fix()
    m.fs.CS_vibrating_screens.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_vibrating_screen_accounts,
            "scaled_param": m.fs.CS_vibrating_screens.area,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 1.7 is CS Conveyors
    CS_conveyors_accounts = ["1.7"]
    m.fs.CS_conveyors = UnitModelBlock()
    m.fs.CS_conveyors.throughput = pyo.Var(
        initialize=575, units=pyunits.ton / pyunits.hr
    )
    m.fs.CS_conveyors.throughput.fix()
    m.fs.CS_conveyors.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_conveyors_accounts,
            "scaled_param": m.fs.CS_conveyors.throughput,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 2.1 is DG Vibrating Screens
    DG_vibrating_screen_accounts = ["2.1"]
    m.fs.DG_vibrating_screens = UnitModelBlock()
    m.fs.DG_vibrating_screens.area = pyo.Var(initialize=332, units=pyunits.ft**2)
    m.fs.DG_vibrating_screens.area.fix()
    m.fs.DG_vibrating_screens.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": DG_vibrating_screen_accounts,
            "scaled_param": m.fs.DG_vibrating_screens.area,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 2.2 is DG Storage Bins
    DG_storage_bins_accounts = ["2.2"]
    m.fs.DG_storage_bins = UnitModelBlock()
    m.fs.DG_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.ton)
    m.fs.DG_storage_bins.capacity.fix()
    m.fs.DG_storage_bins.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": DG_storage_bins_accounts,
            "scaled_param": m.fs.DG_storage_bins.capacity,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 2.3 is DG Air Swept Ball Mill
    DG_air_swept_ball_mill_accounts = ["2.3"]
    m.fs.DG_air_swept_ball_mill = UnitModelBlock()
    m.fs.DG_air_swept_ball_mill.power = pyo.Var(initialize=5609, units=pyunits.hp)
    m.fs.DG_air_swept_ball_mill.power.fix()
    m.fs.DG_air_swept_ball_mill.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": DG_air_swept_ball_mill_accounts,
            "scaled_param": m.fs.DG_air_swept_ball_mill.power,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 2.4 is DG Bucket Elevator
    DG_bucket_elevator_accounts = ["2.4"]
    m.fs.DG_bucket_elevator = UnitModelBlock()
    m.fs.DG_bucket_elevator.n_equip = pyo.Var(initialize=1, units=pyunits.dimensionless)
    m.fs.DG_bucket_elevator.n_equip.fix()
    m.fs.DG_bucket_elevator.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": DG_bucket_elevator_accounts,
            "scaled_param": m.fs.DG_bucket_elevator.n_equip,  # 1 elevator
            "CEPCI_year": CEPCI_year,
        },
    )

    # 2.5 is DG Elevator Motor
    DG_elevator_motor_accounts = ["2.5"]
    m.fs.DG_elevator_motor = UnitModelBlock()
    m.fs.DG_elevator_motor.power = pyo.Var(initialize=58.0, units=pyunits.hp)
    m.fs.DG_elevator_motor.power.fix()
    m.fs.DG_elevator_motor.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": DG_elevator_motor_accounts,
            "scaled_param": m.fs.DG_elevator_motor.power,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.1 is R Storage Bins
    R_storage_bins_accounts = ["3.1"]
    m.fs.R_storage_bins = UnitModelBlock()
    m.fs.R_storage_bins.capacity = pyo.Var(initialize=100, units=pyunits.ton)
    m.fs.R_storage_bins.capacity.fix()
    m.fs.R_storage_bins.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_storage_bins_accounts,
            "scaled_param": m.fs.R_storage_bins.capacity,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.2 is R Conveyors
    R_conveyors_accounts = ["3.2"]
    m.fs.R_conveyors = UnitModelBlock()
    m.fs.R_conveyors.throughput = pyo.Var(
        initialize=575, units=pyunits.ton / pyunits.hr
    )
    m.fs.R_conveyors.throughput.fix()
    m.fs.R_conveyors.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_conveyors_accounts,
            "scaled_param": m.fs.R_conveyors.throughput,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.3 is R Roaster
    R_roaster_accounts = ["3.3"]
    m.fs.R_roaster = UnitModelBlock()
    m.fs.R_roaster.duty = pyo.Var(initialize=737, units=pyunits.MBTU / pyunits.hr)
    m.fs.R_roaster.duty.fix()
    m.fs.R_roaster.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_roaster_accounts,
            "scaled_param": m.fs.R_roaster.duty,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.4 is R Gas Scrubber
    R_gas_scrubber_accounts = ["3.4"]
    m.fs.R_gas_scrubber = UnitModelBlock()
    m.fs.R_gas_scrubber.gas_rate = pyo.Var(
        initialize=11500, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.R_gas_scrubber.gas_rate.fix()
    m.fs.R_gas_scrubber.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_gas_scrubber_accounts,
            "scaled_param": m.fs.R_gas_scrubber.gas_rate,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.5 is R Spray Chamber Quencher (7-60 kcfm)
    R_spray_chamber_quencher_accounts = ["3.5"]
    m.fs.R_spray_chamber_quencher = UnitModelBlock()
    m.fs.R_spray_chamber_quencher.gas_rate = pyo.Var(
        initialize=11500, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.R_spray_chamber_quencher.gas_rate.fix()
    m.fs.R_spray_chamber_quencher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_spray_chamber_quencher_accounts,
            "scaled_param": m.fs.R_spray_chamber_quencher.gas_rate,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 3.7 is R Chiller
    R_chiller_accounts = ["3.7"]
    m.fs.R_chiller = UnitModelBlock()
    m.fs.R_chiller.duty = pyo.Var(initialize=131, units=pyunits.MBTU / pyunits.hr)
    m.fs.R_chiller.duty.fix()
    m.fs.R_chiller.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": R_chiller_accounts,
            "scaled_param": m.fs.R_chiller.duty,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.2 is L PE Tanks
    L_pe_tanks_accounts = ["4.2"]
    m.fs.L_pe_tanks = UnitModelBlock()
    m.fs.L_pe_tanks.capacity = pyo.Var(initialize=164805, units=pyunits.gal)
    m.fs.L_pe_tanks.capacity.fix()
    m.fs.L_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_pe_tanks_accounts,
            "scaled_param": m.fs.L_pe_tanks.capacity,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.3 is L Tank Mixer
    L_tank_mixer_accounts = ["4.3"]
    m.fs.L_tank_mixers = UnitModelBlock()
    m.fs.L_tank_mixers.power = pyo.Var(initialize=474, units=pyunits.hp)
    m.fs.L_tank_mixers.power.fix()
    m.fs.L_tank_mixers.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_tank_mixer_accounts,
            "scaled_param": m.fs.L_tank_mixers.power,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.4 is L Process Pump
    L_pump_accounts = ["4.4"]
    m.fs.L_pump = UnitModelBlock()
    m.fs.L_pump.feed_rate = pyo.Var(initialize=10987, units=pyunits.gal / pyunits.min)
    m.fs.L_pump.feed_rate.fix()
    m.fs.L_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_pump_accounts,
            "scaled_param": m.fs.L_pump.feed_rate,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.5 is L Thickener
    L_thickener_accounts = ["4.5"]
    m.fs.L_thickener = UnitModelBlock()
    m.fs.L_thickener.area = pyo.Var(initialize=22590, units=pyunits.ft**2)
    m.fs.L_thickener.area.fix()
    m.fs.L_thickener.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_thickener_accounts,
            "scaled_param": m.fs.L_thickener.area,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.6 is L Solid Waste Filter Press
    L_filter_press_accounts = ["4.6"]
    m.fs.L_filter_press = UnitModelBlock()
    m.fs.L_filter_press.volume = pyo.Var(initialize=3600, units=pyunits.ft**3)
    m.fs.L_filter_press.volume.fix()
    m.fs.L_filter_press.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_filter_press_accounts,
            "scaled_param": m.fs.L_filter_press.volume,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 4.8 is L Solution Heater
    L_solution_heater_accounts = ["4.8"]
    m.fs.L_solution_heater = UnitModelBlock()
    m.fs.L_solution_heater.duty = pyo.Var(
        initialize=2.4, units=pyunits.MBTU / pyunits.hr
    )
    m.fs.L_solution_heater.duty.fix()
    m.fs.L_solution_heater.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": L_solution_heater_accounts,
            "scaled_param": m.fs.L_solution_heater.duty,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 5.1 is RSX PE Tanks
    RSX_pe_tanks_accounts = ["5.1"]
    m.fs.RSX_pe_tanks = UnitModelBlock()
    m.fs.RSX_pe_tanks.capacity = pyo.Var(initialize=35136, units=pyunits.gal)
    m.fs.RSX_pe_tanks.capacity.fix()
    m.fs.RSX_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pe_tanks_accounts,
            "scaled_param": m.fs.RSX_pe_tanks.capacity,
            "n_equip": 6,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 5.2 is RSX Tank Mixer
    RSX_tank_mixer_accounts = ["5.2"]
    m.fs.RSX_tank_mixers = UnitModelBlock()
    m.fs.RSX_tank_mixers.power = pyo.Var(initialize=20, units=pyunits.hp)
    m.fs.RSX_tank_mixers.power.fix()
    m.fs.RSX_tank_mixers.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": RSX_tank_mixer_accounts,
            "scaled_param": m.fs.RSX_tank_mixers.power,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 5.3 is RSX Process Pump
    RSX_pump_accounts = ["5.3"]
    m.fs.RSX_pump = UnitModelBlock()
    m.fs.RSX_pump.feed_rate = pyo.Var(initialize=7027, units=pyunits.gal / pyunits.min)
    m.fs.RSX_pump.feed_rate.fix()
    m.fs.RSX_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": RSX_pump_accounts,
            "scaled_param": m.fs.RSX_pump.feed_rate,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 5.4 is RSX Mixer Settler
    RSX_mixer_settler_accounts = ["5.4"]
    m.fs.RSX_mixer_settler = UnitModelBlock()
    m.fs.RSX_mixer_settler.volume = pyo.Var(initialize=61107, units=pyunits.gal)
    m.fs.RSX_mixer_settler.volume.fix()
    m.fs.RSX_mixer_settler.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": RSX_mixer_settler_accounts,
            "scaled_param": m.fs.RSX_mixer_settler.volume,
            "n_equip": 6,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 6.1 is CSX PE Tanks
    CSX_pe_tanks_accounts = ["6.1"]
    m.fs.CSX_pe_tanks = UnitModelBlock()
    m.fs.CSX_pe_tanks.capacity = pyo.Var(initialize=1405, units=pyunits.gal)
    m.fs.CSX_pe_tanks.capacity.fix()
    m.fs.CSX_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pe_tanks_accounts,
            "scaled_param": m.fs.CSX_pe_tanks.capacity,
            "n_equip": 5,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 6.2 is CSX Tank Mixer
    CSX_tank_mixer_accounts = ["6.2"]
    m.fs.CSX_tank_mixers = UnitModelBlock()
    m.fs.CSX_tank_mixers.power = pyo.Var(initialize=0.8, units=pyunits.hp)
    m.fs.CSX_tank_mixers.power.fix()
    m.fs.CSX_tank_mixers.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CSX_tank_mixer_accounts,
            "scaled_param": m.fs.CSX_tank_mixers.power,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 6.3 is CSX Process Pump
    CSX_pump_accounts = ["6.3"]
    m.fs.CSX_pump = UnitModelBlock()
    m.fs.CSX_pump.feed_rate = pyo.Var(initialize=281, units=pyunits.gal / pyunits.min)
    m.fs.CSX_pump.feed_rate.fix()
    m.fs.CSX_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CSX_pump_accounts,
            "scaled_param": m.fs.CSX_pump.feed_rate,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 6.4 is CSX Mixer Settler
    CSX_mixer_settler_accounts = ["6.4"]
    m.fs.CSX_mixer_settler = UnitModelBlock()
    m.fs.CSX_mixer_settler.volume = pyo.Var(initialize=2444, units=pyunits.gal)
    m.fs.CSX_mixer_settler.volume.fix()
    m.fs.CSX_mixer_settler.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CSX_mixer_settler_accounts,
            "scaled_param": m.fs.CSX_mixer_settler.volume,
            "n_equip": 6,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 7.1 is SX Wash PE Tanks
    SX_wash_pe_tanks_accounts = ["7.1"]
    m.fs.SX_wash_pe_tanks = UnitModelBlock()
    m.fs.SX_wash_pe_tanks.capacity = pyo.Var(initialize=3514, units=pyunits.gal)
    m.fs.SX_wash_pe_tanks.capacity.fix()
    m.fs.SX_wash_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_pe_tanks_accounts,
            "scaled_param": m.fs.SX_wash_pe_tanks.capacity,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 7.2 is SX Wash Tank Mixer
    SX_wash_tank_mixer_accounts = ["7.2"]
    m.fs.SX_wash_tank_mixers = UnitModelBlock()
    m.fs.SX_wash_tank_mixers.power = pyo.Var(initialize=2, units=pyunits.hp)
    m.fs.SX_wash_tank_mixers.power.fix()
    m.fs.SX_wash_tank_mixers.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_tank_mixer_accounts,
            "scaled_param": m.fs.SX_wash_tank_mixers.power,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 7.3 is SX Wash Process Pump
    SX_wash_pump_accounts = ["7.3"]
    m.fs.SX_wash_pump = UnitModelBlock()
    m.fs.SX_wash_pump.feed_rate = pyo.Var(
        initialize=703, units=pyunits.gal / pyunits.min
    )
    m.fs.SX_wash_pump.feed_rate.fix()
    m.fs.SX_wash_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_pump_accounts,
            "scaled_param": m.fs.SX_wash_pump.feed_rate,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 7.4 is SX Wash Mixer Settler
    SX_wash_mixer_settler_accounts = ["7.4"]
    m.fs.SX_wash_mixer_settler = UnitModelBlock()
    m.fs.SX_wash_mixer_settler.volume = pyo.Var(initialize=18332, units=pyunits.gal)
    m.fs.SX_wash_mixer_settler.volume.fix()
    m.fs.SX_wash_mixer_settler.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_mixer_settler_accounts,
            "scaled_param": m.fs.SX_wash_mixer_settler.volume,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 7.5 is SX Wash Filter Press
    SX_wash_filter_press_accounts = ["7.5"]
    m.fs.SX_wash_filter_press = UnitModelBlock()
    m.fs.SX_wash_filter_press.volume = pyo.Var(initialize=0.26, units=pyunits.ft**3)
    m.fs.SX_wash_filter_press.volume.fix()
    m.fs.SX_wash_filter_press.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": SX_wash_filter_press_accounts,
            "scaled_param": m.fs.SX_wash_filter_press.volume,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 9.2 is REE Precipitation PE Tanks
    reep_pe_tanks_accounts = ["9.2"]
    m.fs.reep_pe_tanks = UnitModelBlock()
    m.fs.reep_pe_tanks.capacity = pyo.Var(initialize=1504, units=pyunits.gal)
    m.fs.reep_pe_tanks.capacity.fix()
    m.fs.reep_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": reep_pe_tanks_accounts,
            "scaled_param": m.fs.reep_pe_tanks.capacity,
            "n_equip": 4,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 9.3 is REE Precipitation Tank Mixer
    reep_tank_mixer_accounts = ["9.3"]
    m.fs.reep_tank_mixers = UnitModelBlock()
    m.fs.reep_tank_mixers.power = pyo.Var(initialize=0.61, units=pyunits.hp)
    m.fs.reep_tank_mixers.power.fix()
    m.fs.reep_tank_mixers.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": reep_tank_mixer_accounts,
            "scaled_param": m.fs.reep_tank_mixers.power,
            "n_equip": 3,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 9.4 is REE Precipitation Process Pump
    reep_pump_accounts = ["9.4"]
    m.fs.reep_pump = UnitModelBlock()
    m.fs.reep_pump.feed_rate = pyo.Var(initialize=70, units=pyunits.gal / pyunits.min)
    m.fs.reep_pump.feed_rate.fix()
    m.fs.reep_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": reep_pump_accounts,
            "scaled_param": m.fs.reep_pump.feed_rate,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 9.5 is REE Precipitation Filter Press
    reep_filter_press_accounts = ["9.5"]
    m.fs.reep_filter_press = UnitModelBlock()
    m.fs.reep_filter_press.volume = pyo.Var(initialize=0.405, units=pyunits.ft**3)
    m.fs.reep_filter_press.volume.fix()
    m.fs.reep_filter_press.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": reep_filter_press_accounts,
            "scaled_param": m.fs.reep_filter_press.volume,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 9.8 is REE Precipitation Roaster
    reep_roaster_accounts = ["9.8"]
    m.fs.reep_roaster = UnitModelBlock()
    m.fs.reep_roaster.duty = pyo.Var(initialize=0.35, units=pyunits.MBTU / pyunits.hr)
    m.fs.reep_roaster.duty.fix()
    m.fs.reep_roaster.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": reep_roaster_accounts,
            "scaled_param": m.fs.reep_roaster.duty,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 11.1 is Water Treatment PE Tanks
    WT_pe_tanks_accounts = ["11.1"]
    m.fs.WT_pe_tanks = UnitModelBlock()
    m.fs.WT_pe_tanks.capacity = pyo.Var(initialize=453131, units=pyunits.gal)
    m.fs.WT_pe_tanks.capacity.fix()
    m.fs.WT_pe_tanks.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": WT_pe_tanks_accounts,
            "scaled_param": m.fs.WT_pe_tanks.capacity,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 11.2 is Water Treatment Process Pump
    WT_pump_accounts = ["11.2"]
    m.fs.WT_pump = UnitModelBlock()
    m.fs.WT_pump.feed_rate = pyo.Var(initialize=78805, units=pyunits.gal / pyunits.min)
    m.fs.WT_pump.feed_rate.fix()
    m.fs.WT_pump.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": WT_pump_accounts,
            "scaled_param": m.fs.WT_pump.feed_rate,
            "n_equip": 2,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 11.3 is Water Treatment Filter Press
    WT_filter_press_accounts = ["11.3"]
    m.fs.WT_filter_press = UnitModelBlock()
    m.fs.WT_filter_press.volume = pyo.Var(initialize=469, units=pyunits.ft**3)
    m.fs.WT_filter_press.volume.fix()
    m.fs.WT_filter_press.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": WT_filter_press_accounts,
            "scaled_param": m.fs.WT_filter_press.volume,
            "CEPCI_year": CEPCI_year,
        },
    )

    # 11.4 is Water Treatment Conveyors
    WT_conveyors_accounts = ["11.4"]
    m.fs.WT_conveyors = UnitModelBlock()
    m.fs.WT_conveyors.throughput = pyo.Var(
        initialize=569, units=pyunits.ton / pyunits.hr
    )
    m.fs.WT_conveyors.throughput.fix()
    m.fs.WT_conveyors.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": WT_conveyors_accounts,
            "scaled_param": m.fs.WT_conveyors.throughput,
            "CEPCI_year": CEPCI_year,
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
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.unit
    def test_REE_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included

        CEPCI_year_units = getattr(pyunits, "USD_" + model.fs.costing.config.CEPCI_year)

        # set some parameter values
        labor_types = [
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ]

        labor_rates = [24.98, 19.08, 30.39, 22.73, 21.97, 45.85]  # USD/hr
        labor_rates_dict = dict(zip(labor_types, labor_rates))

        operators_per_shift = [4, 9, 2, 2, 2, 3]
        operators_per_shift_dict = dict(zip(labor_types, operators_per_shift))

        model.fs.costing.labor_burden.set_value(
            25 * pyunits.percent
        ),  # % fringe benefits

        for l in labor_types:
            model.fs.costing.labor_rates[l].set_value(
                labor_rates_dict[l] * CEPCI_year_units / pyunits.hr
            )
            model.fs.costing.operators_per_shift[l].set_value(
                operators_per_shift_dict[l]
            )

        # add plant-level cost constraints

        model.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feedstock.fix()
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)
        model.fs.feed_grade.fix()

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=8 * 3 * 336,
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
        model.fs.recovery_rate_per_year.fix()

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * model.fs.costing.CEPCI_units
            / pyunits.ton
            * pyunits.convert(model.fs.feedstock, to_units=pyunits.ton / pyunits.hr)
            * model.fs.annual_operating_hours
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
        model.fs.reagents.fix()

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.solid_waste.fix()

        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.precipitate.fix()

        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.dust_and_volatiles.fix()

        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)
        model.fs.power.fix()

        resources = [
            "reagents",
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

        # argument related to Fixed OM costs
        model.fs.costing.build_REE_process_costs(
            feedstock_rate=model.fs.feedstock,
            production_rate=model.fs.recovery_rate_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=dict(zip(resources, rates)),
            resource_prices={
                "reagents": 1 * CEPCI_year_units / pyunits.kg,
            },
            chemicals=["reagents"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
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
                    to_units=model.fs.costing.CEPCI_units,
                )
            )
        )

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_full_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_initialize(self, model):
        # add initialize
        REECostingData.initialize(model.fs.costing)

    @pytest.mark.component
    def test_solve(self, model):
        # try solving
        solver = get_solver()
        results = solver.solve(model, tee=True)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solved_model_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_results(self, model):
        # check some overall cost results

        assert value(model.fs.costing.total_TPC) == pytest.approx(133.23, rel=1e-4)
        assert value(model.fs.costing.total_BEC) == pytest.approx(44.308, rel=1e-4)
        assert value(model.fs.costing.total_overnight_capital) == pytest.approx(
            133.23, rel=1e-4
        )
        assert value(model.fs.costing.other_plant_costs) == pytest.approx(
            1.6309, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            11.519, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            576.07, rel=1e-4
        )
        assert value(model.fs.costing.land_cost) == pytest.approx(
            1.2247, rel=1e-4
        )  # Expression, not Var
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(
            30.061, rel=1e-4
        )
        assert value(model.fs.costing.pv_capital_cost) == pytest.approx(
            -114.1620885,
            rel=1e-4,  # TODO update to -112.78144 once IDAES side is patched
        )
        assert value(model.fs.costing.pv_loan_interest) == pytest.approx(
            -11.1358148,
            rel=1e-4,  # TODO update to -11.001142 once IDAES side is patched
        )
        assert value(model.fs.costing.pv_operating_cost) == pytest.approx(
            -5051.7729, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(257.91371, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(
            -4919.1571, rel=1e-4
        )  # TODO update to -4917.6417 once IDAES side is patched

    @pytest.mark.unit
    def test_report(self, model):
        # test report methods
        REECostingData.report(model.fs.costing, export=True)
        assert os.path.exists(os.path.join(os.getcwd(), "costing_report.csv"))
        # cleanup
        os.remove(os.path.join(os.getcwd(), "costing_report.csv"))
        assert not os.path.exists(os.path.join(os.getcwd(), "costing_report.csv"))

        REECostingData.display_total_plant_costs(model.fs.costing)
        REECostingData.display_bare_erected_costs(model.fs.costing)

    @pytest.mark.unit
    def test_costing_bounding_build_diagnostics(self, model):
        # test costing bounding method
        REECostingData.calculate_REE_costing_bounds(
            b=model.fs.costing,
            capacity=model.fs.feedstock
            * model.fs.annual_operating_hours
            * 20
            * pyunits.year,
            grade=model.fs.feed_grade,
        )

        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_costing_bounding_solve(self, model):

        # solve new variables and constraints
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing_bounding_solve_diagnostics(self, model):

        dt = DiagnosticsToolbox(model=model)
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

            if key == "Chemical Extraction, Enrichment and Separation":
                assert value(
                    model.fs.costing.costing_lower_bound[key]
                ) == pytest.approx(
                    expected_costing_lower_bound[key],
                    abs=1e-8,
                )
            else:
                assert value(
                    model.fs.costing.costing_lower_bound[key]
                ) == pytest.approx(
                    expected_costing_lower_bound[key], abs=1e-8, rel=1e-4
                )

        for key in model.fs.costing.costing_upper_bound.keys():
            assert value(model.fs.costing.costing_upper_bound[key]) == pytest.approx(
                expected_costing_upper_bound[key], rel=1e-4
            )

    @pytest.mark.component
    def test_costing_bounding_rerun(self, model):

        # call again, should give same results

        REECostingData.calculate_REE_costing_bounds(
            b=model.fs.costing,
            capacity=model.fs.feedstock
            * model.fs.annual_operating_hours
            * 20
            * pyunits.year,
            grade=model.fs.feed_grade,
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

            if key == "Chemical Extraction, Enrichment and Separation":
                assert value(
                    model.fs.costing.costing_lower_bound[key]
                ) == pytest.approx(
                    expected_costing_lower_bound[key],
                    abs=1e-8,
                )
            else:
                assert value(
                    model.fs.costing.costing_lower_bound[key]
                ) == pytest.approx(
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

        # Nanofiltration

        model.fs.nf_properties = MCASParameterBlock(
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

        model.fs.nfunit = NanofiltrationDSPMDE0D(
            property_package=model.fs.nf_properties
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
                / model.fs.nfunit.feed_side.properties_in[0].mw_comp[ion]
            )
            model.fs.nfunit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)

        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / model.fs.nfunit.feed_side.properties_in[0].mw_comp["H2O"]
        )
        model.fs.nfunit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
            H2O_mol_comp_flow
        )

        # Use assert electroneutrality method from property model to ensure the ion concentrations provided
        # obey electroneutrality condition
        model.fs.nfunit.feed_side.properties_in[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl_-",
            get_property="mass_frac_phase_comp",
        )

        # Fix other inlet state variables
        model.fs.nfunit.inlet.temperature[0].fix(298.15)
        model.fs.nfunit.inlet.pressure[0].fix(4e5)

        # Fix the membrane variables that are usually fixed for the DSPM-DE model
        model.fs.nfunit.radius_pore.fix(0.5e-9)
        model.fs.nfunit.membrane_thickness_effective.fix(1.33e-6)
        model.fs.nfunit.membrane_charge_density.fix(-27)
        model.fs.nfunit.dielectric_constant_pore.fix(41.3)

        # Fix final permeate pressure to be ~atmospheric
        model.fs.nfunit.mixed_permeate[0].pressure.fix(101325)

        model.fs.nfunit.spacer_porosity.fix(0.85)
        model.fs.nfunit.channel_height.fix(5e-4)
        model.fs.nfunit.velocity[0, 0].fix(0.25)
        model.fs.nfunit.area.fix(50)
        # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
        model.fs.nfunit.spacer_mixing_efficiency.fix()
        model.fs.nfunit.spacer_mixing_length.fix()

        check_dof(model.fs, fail_flag=True)

        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Ca_2+")
        )
        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "SO4_2-")
        )
        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Mg_2+")
        )
        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
        )
        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
        )
        model.fs.nf_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
        )

        calculate_scaling_factors(model.fs)

        model.fs.nfunit.initialize(optarg=solver.options)

        results = solver.solve(model.fs, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Reverse Osmosis

        model.fs.ro_properties = props.NaClParameterBlock()

        model.fs.rounit = ReverseOsmosis1D(
            property_package=model.fs.ro_properties,
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

        model.fs.rounit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        model.fs.rounit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        model.fs.rounit.inlet.pressure[0].fix(feed_pressure)
        model.fs.rounit.inlet.temperature[0].fix(feed_temperature)
        model.fs.rounit.A_comp.fix(A)
        model.fs.rounit.B_comp.fix(B)
        model.fs.rounit.permeate.pressure[0].fix(pressure_atmospheric)
        model.fs.rounit.feed_side.N_Re[0, 0].fix(400)
        model.fs.rounit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
        model.fs.rounit.feed_side.spacer_porosity.fix(0.97)
        model.fs.rounit.feed_side.channel_height.fix(0.001)

        check_dof(model.fs, fail_flag=True)

        model.fs.ro_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        model.fs.ro_properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(model.fs)

        model.fs.rounit.initialize(optarg=solver.options)

        results = solver.solve(model.fs, tee=True)

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

        model.fs.ro_properties = MCASParameterBlock(**ion_props)

        ix_config = {
            "property_package": model.fs.ro_properties,
            "target_ion": target_ion,
        }
        model.fs.ixunit = IonExchange0D(**ix_config)
        model.fs.ixunit.process_flow.properties_in.calculate_state(
            var_args={
                ("flow_vol_phase", "Liq"): 0.5,
                ("conc_mass_phase_comp", ("Liq", target_ion)): 0.1,
                ("pressure", None): 101325,
                ("temperature", None): 298,
            },
            hold_state=True,
        )

        model.fs.ixunit.process_flow.properties_in[0].flow_mass_phase_comp[...]
        model.fs.ixunit.process_flow.properties_out[0].flow_mass_phase_comp[...]
        model.fs.ixunit.regeneration_stream[0].flow_mass_phase_comp[...]

        model.fs.ixunit.service_flow_rate.fix(15)
        model.fs.ixunit.langmuir[target_ion].fix(0.9)
        model.fs.ixunit.resin_max_capacity.fix(3)
        model.fs.ixunit.bed_depth.fix(1.7)
        model.fs.ixunit.dimensionless_time.fix()
        model.fs.ixunit.number_columns.fix(8)
        model.fs.ixunit.resin_diam.fix()
        model.fs.ixunit.resin_bulk_dens.fix()
        model.fs.ixunit.bed_porosity.fix()

        check_dof(model.fs, fail_flag=True)

        model.fs.ro_properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        model.fs.ro_properties.set_default_scaling(
            "flow_mol_phase_comp", 10, index=("Liq", "Ca_2+")
        )

        calculate_scaling_factors(model.fs)

        model.fs.ixunit.initialize(optarg=solver.options)

        results = solver.solve(model.fs, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Nanofiltration Zero Order

        model.fs.db = Database()
        model.fs.nfzo_params = WaterParameterBlock(solute_list=["tds", "dye"])

        model.fs.nfzounit = NanofiltrationZO(
            property_package=model.fs.nfzo_params,
            database=model.fs.db,
            process_subtype="rHGO_dye_rejection",
        )

        model.fs.nfzounit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        model.fs.nfzounit.inlet.flow_mass_comp[0, "tds"].fix(1)
        model.fs.nfzounit.inlet.flow_mass_comp[0, "dye"].fix(2)

        model.fs.nfzounit.load_parameters_from_database(use_default_removal=True)

        results = solver.solve(model.fs, tee=True)

        # # Check for optimal solution
        assert_optimal_termination(results)

        return model

    @pytest.mark.solver
    @pytest.mark.component
    def test_condition_number(self, model):

        # Check condition number to confirm scaling for nanofiltration
        jac, _ = get_jacobian(model.fs.nfunit, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            5.6905e14, rel=1e-3
        )

        # Check condition number to confirm scaling for RO
        jac, _ = get_jacobian(model.fs.rounit, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            4.29867e16, rel=1e-3
        )

        # Check condition number to confirm scaling for ion exchange
        jac, _ = get_jacobian(model.fs.ixunit, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            4.59436e15, rel=1e-3
        )

    @pytest.mark.component
    def test_REE_watertap_equipment_costing(self, model):

        from watertap.costing import WaterTAPCosting
        from watertap.costing.zero_order_costing import ZeroOrderCosting

        # loop through WaterTAP blocks and import costing from WaterTAP
        model.fs.watertap_costing = WaterTAPCosting()
        model.fs.zeroorder_costing = ZeroOrderCosting()

        model.fs.nfunit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.watertap_costing
        )
        model.fs.costing._registered_unit_costing.append(
            model.fs.nfunit.costing
        )  # pylint: disable=protected-access

        model.fs.rounit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.watertap_costing
        )
        model.fs.costing._registered_unit_costing.append(
            model.fs.rounit.costing
        )  # pylint: disable=protected-access

        model.fs.ixunit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.watertap_costing
        )
        model.fs.costing._registered_unit_costing.append(
            model.fs.ixunit.costing
        )  # pylint: disable=protected-access

        model.fs.nfzounit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.zeroorder_costing
        )
        model.fs.costing._registered_unit_costing.append(
            model.fs.nfzounit.costing
        )  # pylint: disable=protected-access

    @pytest.mark.component
    def test_REE_watertap_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included

        CEPCI_year_units = getattr(pyunits, "USD_" + model.fs.costing.config.CEPCI_year)

        # set some parameter values
        labor_types = [
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ]

        labor_rates = [24.98, 19.08, 30.39, 22.73, 21.97, 45.85]  # USD/hr
        labor_rates_dict = dict(zip(labor_types, labor_rates))

        operators_per_shift = [4, 9, 2, 2, 2, 3]
        operators_per_shift_dict = dict(zip(labor_types, operators_per_shift))

        model.fs.costing.labor_burden.set_value(
            25 * pyunits.percent
        ),  # % fringe benefits

        for l in labor_types:
            model.fs.costing.labor_rates[l].set_value(
                labor_rates_dict[l] * CEPCI_year_units / pyunits.hr
            )
            model.fs.costing.operators_per_shift[l].set_value(
                operators_per_shift_dict[l]
            )

        # add plant-level cost constraints

        model.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feedstock.fix()
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)
        model.fs.feed_grade.fix()

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=8 * 3 * 336,
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
        model.fs.recovery_rate_per_year.fix()

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * model.fs.costing.CEPCI_units
            / pyunits.ton
            * pyunits.convert(model.fs.feedstock, to_units=pyunits.ton / pyunits.hr)
            * model.fs.annual_operating_hours
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
        model.fs.reagents.fix()

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.solid_waste.fix()

        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.precipitate.fix()

        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.dust_and_volatiles.fix()

        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)
        model.fs.power.fix()

        resources = [
            "reagents",
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

        # argument related to Fixed OM costs
        model.fs.costing.build_REE_process_costs(
            feedstock_rate=model.fs.feedstock,
            production_rate=model.fs.recovery_rate_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=dict(zip(resources, rates)),
            resource_prices={
                "reagents": 1 * CEPCI_year_units / pyunits.kg,
            },
            chemicals=["reagents"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
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
                    to_units=model.fs.costing.CEPCI_units,
                )
            )
        )

    @pytest.mark.component
    def test_REE_watertap_costing_initialize(self, model, solver):

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

        REECostingData.initialize(model.fs.costing)

    @pytest.mark.component
    def test_REE_watertap_costing_solve(self, model, solver):

        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_REE_watertap_costing_results_totalcapex(self, model):

        assert value(model.fs.costing.total_BEC) == pytest.approx(56.5333, rel=1e-4)

    @pytest.mark.component
    def test_REE_watertap_costing_results_equipmentcapex(self, model):

        assert value(
            pyunits.convert(
                model.fs.nfunit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0015159, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.rounit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0016148, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.ixunit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(4.0354, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.nfzounit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(8.1867, rel=1e-4)

        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs.nfunit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.rounit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.ixunit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.nfzounit.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(44.308, rel=1e-4)

    @pytest.mark.component
    def test_REE_watertap_costing_results_fixedopex(self, model):

        assert value(
            pyunits.convert(
                model.fs.nfunit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00015159, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.rounit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00016148, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.ixunit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.037284, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.nfzounit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(1.63734, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.nfunit.costing.fixed_operating_cost
                + model.fs.rounit.costing.fixed_operating_cost
                + model.fs.ixunit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
            + pyunits.convert(
                model.fs.nfzounit.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(1.6749, rel=1e-4)

        assert value(model.fs.costing.custom_fixed_costs) == pytest.approx(
            1.67493, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            14.2828, rel=1e-4
        )

    @pytest.mark.component
    def test_REE_watertap_costing_variableopex(self, model):

        assert value(model.fs.costing.custom_variable_costs) == pytest.approx(
            0, abs=1e-4
        )

        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            576.61838, rel=1e-4
        )


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
            == pyunits.USD_UKy_2019
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

        CEPCI_year_units = getattr(pyunits, "USD_" + model.fs.costing.config.CEPCI_year)

        # set some parameter values
        labor_types = [
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ]

        labor_rates = [24.98, 19.08, 30.39, 22.73, 21.97, 45.85]  # USD/hr
        labor_rates_dict = dict(zip(labor_types, labor_rates))

        operators_per_shift = [4, 9, 2, 2, 2, 3]
        operators_per_shift_dict = dict(zip(labor_types, operators_per_shift))

        model.fs.costing.labor_burden.set_value(
            25 * pyunits.percent
        ),  # % fringe benefits

        for l in labor_types:
            model.fs.costing.labor_rates[l].set_value(
                labor_rates_dict[l] * CEPCI_year_units / pyunits.hr
            )
            model.fs.costing.operators_per_shift[l].set_value(
                operators_per_shift_dict[l]
            )

        # add plant-level cost constraints

        model.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feedstock.fix()
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)
        model.fs.feed_grade.fix()

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=8 * 3 * 336,
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
        model.fs.recovery_rate_per_year.fix()

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * model.fs.costing.CEPCI_units
            / pyunits.ton
            * pyunits.convert(model.fs.feedstock, to_units=pyunits.ton / pyunits.hr)
            * model.fs.annual_operating_hours
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
        model.fs.reagents.fix()

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.solid_waste.fix()

        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.precipitate.fix()

        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.dust_and_volatiles.fix()

        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)
        model.fs.power.fix()

        resources = [
            "reagents",
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

        # argument related to Fixed OM costs
        model.fs.costing.build_REE_process_costs(
            feedstock_rate=model.fs.feedstock,
            production_rate=model.fs.recovery_rate_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=dict(zip(resources, rates)),
            resource_prices={
                "reagents": 1 * CEPCI_year_units / pyunits.kg,
            },
            chemicals=["reagents"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
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
                    to_units=model.fs.costing.CEPCI_units,
                )
            )
        )

    @pytest.mark.component
    def test_REE_custom_costing_initialize(self, model):

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

        # check model structural diagnostics
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

        REECostingData.initialize(model.fs.costing)

    @pytest.mark.component
    def test_REE_custom_costing_solve(self, model):

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

        # check model numerical diagnostics
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_REE_custom_costing_results(self, model):

        assert value(model.fs.costing.total_BEC) == pytest.approx(44.377, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0686081, rel=1e-4)
        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs.custom_vessel.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(44.308, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00343040, rel=1e-4)

        assert value(model.fs.costing.custom_fixed_costs) == pytest.approx(
            0.00343040, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            11.52811, rel=1e-4
        )

        assert value(
            pyunits.convert(
                model.fs.custom_vessel.costing.variable_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(4.13750, rel=1e-4)

        assert value(model.fs.costing.custom_variable_costs) == pytest.approx(
            4.13750, rel=1e-4
        )

        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            581.03245, rel=1e-4
        )


class TestDiafiltrationCosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = base_model()

        # create dummy blocks to store the REEUnitModelCostingBlocks
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
            costing_method=DiafiltrationCostingData.cost_membrane_pressure_drop_utility,
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
                "inlet_pressure": pyunits.convert(
                    m.fs.P_atm, to_units=pyunits.kPa
                ),  # units of kPa
                "outlet_pressure": 1e-5  # assume numerically 0 since SEC accounts for feed pump OPEX
                * pyunits.psi,  # this should make m.fs.feed_pump.costing.variable_operating_cost ~0
                "inlet_vol_flow": m.fs.stage3.retentate_flow_vol,  # feed
            },
        )
        m.fs.diafiltrate_pump.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=DiafiltrationCostingData.cost_pump,
            costing_method_arguments={
                "inlet_pressure": pyunits.convert(
                    m.fs.P_atm, to_units=pyunits.kPa
                ),  # units of kPa
                "outlet_pressure": pyunits.convert(
                    m.fs.P_op, to_units=pyunits.psi
                ),  # units of psi
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
            assert blk.costing.costing_package.base_currency == pyunits.USD_UKy_2019
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
    def test_REE_diafiltration_costing(self, model):
        # full smoke test with all components, O&M costs, and extra costs included

        CEPCI_year_units = getattr(pyunits, "USD_" + model.fs.costing.config.CEPCI_year)

        # set some parameter values
        labor_types = [
            "skilled",
            "unskilled",
            "supervisor",
            "maintenance",
            "technician",
            "engineer",
        ]

        labor_rates = [24.98, 19.08, 30.39, 22.73, 21.97, 45.85]  # USD/hr
        labor_rates_dict = dict(zip(labor_types, labor_rates))

        operators_per_shift = [4, 9, 2, 2, 2, 3]
        operators_per_shift_dict = dict(zip(labor_types, operators_per_shift))

        model.fs.costing.labor_burden.set_value(
            25 * pyunits.percent
        ),  # % fringe benefits

        for l in labor_types:
            model.fs.costing.labor_rates[l].set_value(
                labor_rates_dict[l] * CEPCI_year_units / pyunits.hr
            )
            model.fs.costing.operators_per_shift[l].set_value(
                operators_per_shift_dict[l]
            )

        # add plant-level cost constraints

        model.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        model.fs.feedstock.fix()
        model.fs.feed_grade = pyo.Var(initialize=356.64, units=pyunits.ppm)
        model.fs.feed_grade.fix()

        # for convenience
        model.fs.annual_operating_hours = pyo.Param(
            initialize=8 * 3 * 336,
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
        model.fs.recovery_rate_per_year.fix()

        # the land cost is the lease cost, or refining cost of REO produced
        model.fs.land_cost = pyo.Expression(
            expr=0.303736
            * 1e-6
            * model.fs.costing.CEPCI_units
            / pyunits.ton
            * pyunits.convert(model.fs.feedstock, to_units=pyunits.ton / pyunits.hr)
            * model.fs.annual_operating_hours
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
        model.fs.reagents.fix()

        model.fs.solid_waste = pyo.Var(
            model.fs.time, initialize=11136 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous solid waste
        model.fs.solid_waste.fix()

        model.fs.precipitate = pyo.Var(
            model.fs.time, initialize=732 / 24, units=pyunits.ton / pyunits.hr
        )  # non-hazardous precipitate
        model.fs.precipitate.fix()

        model.fs.dust_and_volatiles = pyo.Var(
            model.fs.time, initialize=120 / 24, units=pyunits.ton / pyunits.hr
        )  # dust and volatiles
        model.fs.dust_and_volatiles.fix()

        model.fs.power = pyo.Var(model.fs.time, initialize=14716, units=pyunits.hp)
        model.fs.power.fix()

        resources = [
            "reagents",
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

        # argument related to Fixed OM costs
        model.fs.costing.build_REE_process_costs(
            feedstock_rate=model.fs.feedstock,
            production_rate=model.fs.recovery_rate_per_year,
            pure_product_output_rates=pure_product_output_rates,
            mixed_product_output_rates=mixed_product_output_rates,
            # arguments related to total owners costs
            land_cost=model.fs.land_cost,
            resources=dict(zip(resources, rates)),
            resource_prices={
                "reagents": 1 * CEPCI_year_units / pyunits.kg,
            },
            chemicals=["reagents"],
            waste=[
                "nonhazardous_solid_waste",
                "nonhazardous_precipitate_waste",
                "dust_and_volatiles",
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
                    to_units=model.fs.costing.CEPCI_units,
                )
            )
        )

        # check that the model is set up properly and has 0 degrees of freedom
        assert degrees_of_freedom(model) == 0

        # check model structural diagnostics
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

        REECostingData.initialize(model.fs.costing)

        solver = get_solver()
        results = solver.solve(model, tee=False)
        assert_optimal_termination(results)

        # check model numerical diagnostics
        dt.assert_no_numerical_warnings()

        assert value(model.fs.costing.total_BEC) == pytest.approx(44.351, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage1.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage2.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage3.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.0043043, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.feed_pump.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.014780, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.diafiltrate_pump.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(0.014780, rel=1e-4)
        assert value(
            model.fs.costing.total_BEC
            - pyunits.convert(
                model.fs.stage1.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.stage2.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.stage3.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.feed_pump.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
            - pyunits.convert(
                model.fs.diafiltrate_pump.costing.capital_cost,
                to_units=model.fs.costing.CEPCI_units,
            )
        ) == pytest.approx(44.308, rel=1e-4)

        assert value(
            pyunits.convert(
                model.fs.stage1.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage2.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.stage3.costing.fixed_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.00086087, rel=1e-4)

        assert value(model.fs.costing.custom_fixed_costs) == pytest.approx(
            0.0025826, rel=1e-4
        )

        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            11.52493, rel=1e-4
        )

        assert value(
            pyunits.convert(
                model.fs.cascade.costing.variable_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(0.985221, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.feed_pump.costing.variable_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(2.91544e-10, rel=1e-4)
        assert value(
            pyunits.convert(
                model.fs.diafiltrate_pump.costing.variable_operating_cost,
                to_units=model.fs.costing.CEPCI_units / pyunits.year,
            )
        ) == pytest.approx(8.51301e-4, rel=1e-4)

        assert value(model.fs.costing.custom_variable_costs) == pytest.approx(
            0.986072, rel=1e-4
        )

        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            577.2501, rel=1e-4
        )


class TestHDDRecyclingCosting(object):
    @pytest.fixture(scope="class")
    def model(self):
        CEPCI_year = "2019"

        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()

        # add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = REECosting(
            Lang_factor=2.97,
            has_fixed_OM=True,
            CEPCI_year=CEPCI_year,
        )

        # 1.1 is Front End Loader (2 cuyd)
        # this is a constant-cost unit, where n_equip is the scaling parameter
        CS_front_end_loader_2yd3_accounts = ["1.1"]
        m.fs.CS_front_end_loader_2yd3 = UnitModelBlock()
        m.fs.CS_front_end_loader_2yd3.n_equip = pyo.Var(
            initialize=1, units=pyunits.dimensionless
        )
        m.fs.CS_front_end_loader_2yd3.n_equip.fix()

        m.fs.CS_front_end_loader_2yd3.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_front_end_loader_2yd3_accounts,
                "scaled_param": m.fs.CS_front_end_loader_2yd3.n_equip,  # 1 loader
                "n_equip": 5,
                "CEPCI_year": CEPCI_year,
            },
        )

        # 13.1 is HDD shredder
        # this is a constant-cost unit, where n_equip is the scaling parameter
        HDD_Recycling_shredder_accounts = ["13.1"]
        m.fs.HDD_Recycling_shredder = UnitModelBlock()
        m.fs.HDD_Recycling_shredder.n_equip = pyo.Var(
            initialize=1, units=pyunits.dimensionless
        )
        m.fs.HDD_Recycling_shredder.n_equip.fix()

        m.fs.HDD_Recycling_shredder.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": HDD_Recycling_shredder_accounts,
                "scaled_param": m.fs.HDD_Recycling_shredder.n_equip,  # 1 shredder
                "CEPCI_year": CEPCI_year,
            },
        )

        # 13.2 is Hydrogen Decrepitation
        HDD_Recycling_HD_accounts = ["13.2"]
        m.fs.HDD_Recycling_HD = UnitModelBlock()
        m.fs.HDD_Recycling_HD.duty = pyo.Var(
            initialize=10, units=pyunits.MBTU / pyunits.hr
        )
        m.fs.HDD_Recycling_HD.duty.fix()
        m.fs.HDD_Recycling_HD.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": HDD_Recycling_HD_accounts,
                "scaled_param": m.fs.HDD_Recycling_HD.duty,
                "CEPCI_year": CEPCI_year,
            },
        )

        return m

    @pytest.mark.unit
    def test_HDD_Recycling_costing_build_diagnostics(self, model):
        model.fs.costing.build_REE_process_costs()

        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_HDD_Recycling_costing_initialize(self, model):
        REECostingData.initialize(model.fs.costing)

    @pytest.mark.component
    def test_HDD_Recycling_costing_solve(self, model):

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

        # check model numerical diagnostics
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_HDD_Recycling_costing_results(self, model):
        CS_front_end_loader_2yd3_accounts = ["1.1"]
        HDD_Recycling_shredder_accounts = ["13.1"]
        HDD_Recycling_HD_accounts = ["13.2"]

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

        assert value(model.fs.costing.total_TPC) == pytest.approx(3.8220, rel=1e-4)


class TestNPVCostingBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        CEPCI_year = "2021"

        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = REECosting(
            Lang_factor=2.97,
            has_fixed_OM=True,
            has_variable_OM=True,
            has_net_present_value=True,
            has_capital_expenditure_period=True,
            capital_expenditure_percentages=[10, 60, 30],
            CEPCI_year=CEPCI_year,
        )

        # 1.3 is CS Jaw Crusher
        CS_jaw_crusher_accounts = ["1.3"]
        m.fs.CS_jaw_crusher = UnitModelBlock()
        m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
        m.fs.CS_jaw_crusher.power.fix()
        m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "CEPCI_year": CEPCI_year,
            },
        )

        m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
        m.fs.feedstock.fix()

        m.fs.water = pyo.Var(
            m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr
        )
        m.fs.water.fix()

        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            feedstock_rate=m.fs.feedstock,
            resources={
                "water": m.fs.water,
            },
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
        assert isinstance(model.fs.costing.debt_percentage_of_capex, pyo.Param)
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
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_NPV_costingblock_initialize(self, model):
        REECostingData.initialize(model.fs.costing)

    @pytest.mark.component
    def test_NPV_costingblock_solve(self, model):
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_NPV_costingblock_solve_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_NPV_costingblock_results(self, model):
        # values that are hard-coded in the following test class
        assert value(model.fs.costing.total_BEC) == pytest.approx(2.512179, rel=1e-4)
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            6.001094, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            1.23406, rel=1e-4
        )
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(
            30.36196, rel=1e-4
        )

        # check some NPV results
        assert value(model.fs.costing.pv_capital_cost) == pytest.approx(
            -6.3162037, rel=1e-4
        )
        assert value(model.fs.costing.pv_loan_interest) == pytest.approx(
            -0.61610712, rel=1e-4
        )
        assert value(model.fs.costing.pv_operating_cost) == pytest.approx(
            -62.07505, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(260.49501, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(191.4877, rel=1e-4)


class TestNPVFixedInputs(object):
    @pytest.fixture(scope="class")
    def model(self):
        # Create a concrete model as the top level object
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
        m.fs.costing = REECosting(
            Lang_factor=2.97,
            has_fixed_OM=True,
            has_net_present_value=True,
            has_capital_expenditure_period=True,
            capital_expenditure_percentages=[10, 60, 30],
            CEPCI_year="2021",
        )

        # values from previous test class
        m.fs.costing.build_REE_process_costs(
            total_purchase_cost=2.512179 * m.fs.costing.CEPCI_units,
            annual_fixed_operating_cost=(6.001094 + 1.23406)
            * m.fs.costing.CEPCI_units
            / pyunits.year,
            annual_revenue=30.36196 * m.fs.costing.CEPCI_units / pyunits.year,
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
        assert isinstance(model.fs.costing.debt_percentage_of_capex, pyo.Param)
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
        dt = DiagnosticsToolbox(model=model)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_NPV_fixedinputs_solve(self, model):
        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_NPV_fixedinputs_solve_diagnostics(self, model):
        dt = DiagnosticsToolbox(model=model)
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
            -62.07505, rel=1e-4
        )
        assert value(model.fs.costing.pv_revenue) == pytest.approx(260.49501, rel=1e-4)
        assert value(model.fs.costing.npv) == pytest.approx(191.4877, rel=1e-4)


@pytest.mark.component
def test_REE_costing_CEPCI_year():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_2021 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2021.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2021.power.fix()
    m.fs.CS_jaw_crusher_2021.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2021.power,
            "CEPCI_year": "2021",
        },
    )

    m.fs.CS_jaw_crusher_UKy_2019 = UnitModelBlock()
    m.fs.CS_jaw_crusher_UKy_2019.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_UKy_2019.power.fix()
    m.fs.CS_jaw_crusher_UKy_2019.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_UKy_2019.power,
            "CEPCI_year": "UKy_2019",
        },
    )

    m.fs.CS_jaw_crusher_2022 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2022.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2022.power.fix()
    m.fs.CS_jaw_crusher_2022.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2022.power,
            "CEPCI_year": "2022",
        },
    )

    m.fs.CS_jaw_crusher_2025 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2025.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2025.power.fix()
    m.fs.CS_jaw_crusher_2025.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2025.power,
            "CEPCI_year": "2025",
        },
    )

    m.fs.CS_jaw_crusher_CE500 = UnitModelBlock()
    m.fs.CS_jaw_crusher_CE500.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher_CE500.power.fix()
    m.fs.CS_jaw_crusher_CE500.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_CE500.power,
            "CEPCI_year": "CE500",
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
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

    with pytest.raises(
        AttributeError,
        match="CEPCI_year notayear is not a valid currency base option. "
        "Valid CEPCI options include CE500, CE394, years from 1990 to 2023, or "
        "user-defined units such as 2019_Sep and UKy_2019.",
    ):
        m.fs.costing = REECosting(CEPCI_year="notayear")


@pytest.mark.unit
def test_REE_costing_nonexistentcostaccount():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["notanaccount"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        KeyError,
        # TODO fails as expected but on wrong part of the account loading, look into this later
        # match="Account notanaccount could not be found in the dictionary for source 1",
    ):
        m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
            },
        )


@pytest.mark.component
def test_REE_costing_multipleaccountssameparameter():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher, 1.5 is Roll Crusher
    CS_crusher_accounts = ["1.3", "1.5"]
    m.fs.CS_crusher = UnitModelBlock()
    m.fs.CS_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_crusher.power.fix()

    m.fs.CS_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_crusher_accounts,
            "scaled_param": m.fs.CS_crusher.power,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    # 1.5 is CS Roll Crusher
    CS_roll_crusher_accounts = ["1.5"]
    m.fs.CS_roll_crusher = UnitModelBlock()
    m.fs.CS_roll_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_roll_crusher.power.fix()
    m.fs.CS_roll_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_roll_crusher_accounts,
            "scaled_param": m.fs.CS_roll_crusher.power,
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
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
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher, 1.6 is Vibrating Screen
    CS_crusher_accounts = ["1.3", "1.6"]
    m.fs.CS_crusher = UnitModelBlock()
    m.fs.CS_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_crusher.power.fix()

    with pytest.raises(
        ValueError,
        match="fs.CS_crusher.costing cost accounts selected do not use the same process parameter",
    ):
        m.fs.CS_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_crusher_accounts,
                "scaled_param": m.fs.CS_crusher.power,
            },
        )


@pytest.mark.component
def test_REE_costing_additionalcostingparams_newaccount():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "10": {
            "A": {
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
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3new"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "additional_costing_params": [
                additional_costing_params,
            ],
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # adding a check just to make sure everything works as expected
    assert value(
        m.fs.CS_jaw_crusher.costing.bare_erected_cost["1.3new"]
    ) == pytest.approx(2.5122, rel=1e-4)


@pytest.mark.component
def test_REE_costing_additionalcostingparams_newaccount_notlist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "10": {
            "A": {
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
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3new"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        TypeError,
        match="additional_costing_params must be a list of dicts, not a single dict, "
        "e.g. \\[{'1': data}, {'2': data},\\] or \\[{'1': data},\\] and not {'1': data}.",
    ):
        m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "additional_costing_params": additional_costing_params,
            },
        )


@pytest.mark.component
def test_REE_costing_additionalcostingparams_newaccount_notlistofdicts():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "10": {
            "A": {
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
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3new"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        TypeError,
        match="additional_costing_params must be a list of dicts, not a single dict, "
        "e.g. \\[{'1': data}, {'2': data},\\] or \\[{'1': data},\\] and not {'1': data}.",
    ):
        m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "additional_costing_params": [additional_costing_params, 5],
            },
        )


@pytest.mark.component
def test_REE_costing_additionalcostingparams_overwrite():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "10": {
            "A": {
                "1.3": {
                    "Account Name": "UKy Crushing and Screening - Jaw Crusher",
                    "BEC": 1922101.0,
                    "BEC_units": "$2016",
                    "Exponent": 1.25,
                    "Process Parameter": "Power Draw (hp)",
                    "RP Value": 689.0,  # changed this value from 589.0
                    "Units": "hp",
                },
            }
        },
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
            "additional_costing_params": [
                additional_costing_params,
            ],
            "use_additional_costing_params": True,
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
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
    m.fs.costing = REECosting()

    # create new account to test that doesn't exist
    additional_costing_params = {
        "10": {
            "A": {
                "1.3": {
                    "Account Name": "UKy Crushing and Screening - Jaw Crusher",
                    "BEC": 1922101.0,
                    "BEC_units": "$2016",
                    "Exponent": 1.25,
                    "Process Parameter": "Power Draw (hp)",
                    "RP Value": 689.0,  # changed this value from 589.0
                    "Units": "hp",
                },
            }
        },
    }

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()

    with pytest.raises(
        ValueError,
        match="Data already exists for Account 1.3 using technology 10 with CCS A. Please "
        "confirm that the custom account dictionary is correct, or add the "
        "new parameters as a new account. To use the custom account dictionary "
        "for all conflicts, please pass the argument use_additional_costing_params "
        "as True.",
    ):
        m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=REECostingData.get_equipment_costing,
            costing_method_arguments={
                "cost_accounts": CS_jaw_crusher_accounts,
                "scaled_param": m.fs.CS_jaw_crusher.power,
                "additional_costing_params": [
                    additional_costing_params,
                ],
                "use_additional_costing_params": False,
            },
        )


@pytest.mark.component
def test_REE_costing_scaledownparallelequip():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_1 = UnitModelBlock()
    m.fs.CS_jaw_crusher_1.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_1.power.fix()
    m.fs.CS_jaw_crusher_1.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_1.power,
            "scale_down_parallel_equip": True,
            "n_equip": 1,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_2 = UnitModelBlock()
    m.fs.CS_jaw_crusher_2.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_2.power.fix()
    m.fs.CS_jaw_crusher_2.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_2.power,
            "scale_down_parallel_equip": False,
            "n_equip": 1,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_3 = UnitModelBlock()
    m.fs.CS_jaw_crusher_3.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_3.power.fix()
    m.fs.CS_jaw_crusher_3.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_3.power,
            "scale_down_parallel_equip": True,
            "n_equip": 2,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_4 = UnitModelBlock()
    m.fs.CS_jaw_crusher_4.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_4.power.fix()
    m.fs.CS_jaw_crusher_4.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_4.power,
            "scale_down_parallel_equip": False,
            "n_equip": 2,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_5 = UnitModelBlock()
    m.fs.CS_jaw_crusher_5.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_5.power.fix()
    m.fs.CS_jaw_crusher_5.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_5.power,
            "scale_down_parallel_equip": True,
            "n_equip": 5,
        },
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher_6 = UnitModelBlock()
    m.fs.CS_jaw_crusher_6.power = pyo.Var(initialize=281.9, units=pyunits.hp)
    m.fs.CS_jaw_crusher_6.power.fix()
    m.fs.CS_jaw_crusher_6.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher_6.power,
            "scale_down_parallel_equip": False,
            "n_equip": 5,
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
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


@pytest.mark.component
def test_REE_costing_usersetTPC_noOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.costing.build_REE_process_costs(
        total_purchase_cost=1 * m.fs.costing.CEPCI_units,
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert pyunits.get_units(m.fs.costing.total_TPC) == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_BEC_eq")
    assert hasattr(m.fs.costing, "total_TPC")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert value(m.fs.costing.total_BEC) == pytest.approx(1.0000, rel=1e-4)
    assert value(m.fs.costing.total_TPC) == pytest.approx(2.9700, rel=1e-4)


@pytest.mark.component
def test_REE_costing_usersetTPC_withOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_REE_process_costs(
        total_purchase_cost=1 * m.fs.costing.CEPCI_units,
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        feedstock_rate=m.fs.feedstock,
        resources={"water": m.fs.water},
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert pyunits.get_units(m.fs.costing.total_TPC) == pyunits.MUSD_2021
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_BEC_eq")
    assert hasattr(m.fs.costing, "total_TPC")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert value(m.fs.costing.total_BEC) == pytest.approx(1.0000, rel=1e-4)
    assert value(m.fs.costing.total_TPC) == pytest.approx(2.9700, rel=1e-4)


@pytest.mark.component
def test_REE_costing_fixedOM_productratesnotpassed():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.costing.build_REE_process_costs()

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    assert value(m.fs.costing.total_sales_revenue) == pytest.approx(0, abs=1e-6)


@pytest.mark.unit
def test_REE_costing_fixedOM_pureproductnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates=[
                1.9 * pyunits.kg / pyunits.hr,
            ],
        )


@pytest.mark.unit
def test_REE_costing_fixedOM_mixedproductnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    with pytest.raises(
        TypeError,
        match="product_output_rates argument must be a dict",
    ):
        m.fs.costing.build_REE_process_costs(
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
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.costing.build_REE_process_costs(
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
    m.fs.costing = REECosting()

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    with pytest.raises(
        TypeError,
        match="Dictionary of custom sale_prices must be a dict object.",
    ):
        m.fs.costing.build_REE_process_costs(
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
    m.fs.costing = REECosting(
        has_fixed_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    with pytest.raises(
        AttributeError,
        match="A pure product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Al', 'Sb', 'As', 'Ba', 'Be', 'Bi', 'Ce', 'Cs', 'Cr', 'Co', 'Dy', 'Er', 'Eu', 'CaF2', 'Gd', 'Ga', 'Ge', 'C', 'Ha', 'Ho', 'In', 'Ir', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Nd', 'Ni', 'Nb', 'Pd', 'Pt', 'Pr', 'Rh', 'Rb', 'Ru', 'Sm', 'Sc', 'Ta', 'Te', 'Tb', 'Tm', 'Sn', 'Ti', 'W', 'V', 'Yb', 'Y', 'Zn', 'Zr', 'CeO2', 'Dy2O3', 'Eu2O3', 'La2O3', 'Nd2O3', 'Sc2O3', 'Ta2O5', 'Tb4O7', 'TiO2', 'WO3', 'Y2O3', 'Er2O3', 'Ho2O3', 'Gd2O3', 'Lu2O3', 'Pr6O11', 'Sm2O3', 'Tm2O3', 'Yb2O3'\\]",
    ):
        m.fs.costing.build_REE_process_costs(
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
    m.fs.costing = REECosting(
        has_fixed_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    with pytest.raises(
        AttributeError,
        match="A mixed product was included that does not contain a sale price. "
        "Sale prices exist for the following products: \\['Al', 'Sb', 'As', 'Ba', 'Be', 'Bi', 'Ce', 'Cs', 'Cr', 'Co', 'Dy', 'Er', 'Eu', 'CaF2', 'Gd', 'Ga', 'Ge', 'C', 'Ha', 'Ho', 'In', 'Ir', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Nd', 'Ni', 'Nb', 'Pd', 'Pt', 'Pr', 'Rh', 'Rb', 'Ru', 'Sm', 'Sc', 'Ta', 'Te', 'Tb', 'Tm', 'Sn', 'Ti', 'W', 'V', 'Yb', 'Y', 'Zn', 'Zr', 'CeO2', 'Dy2O3', 'Eu2O3', 'La2O3', 'Nd2O3', 'Sc2O3', 'Ta2O5', 'Tb4O7', 'TiO2', 'WO3', 'Y2O3', 'Er2O3', 'Ho2O3', 'Gd2O3', 'Lu2O3', 'Pr6O11', 'Sm2O3', 'Tm2O3', 'Yb2O3'\\]",
    ):
        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "newprod": 0.00143 * pyunits.kg / pyunits.hr,
            },
        )


@pytest.mark.unit
def test_REE_costing_variableOM_nofixedOM():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_REE_process_costs(
        feedstock_rate=m.fs.feedstock,
        resources={
            "water": m.fs.water,
        },
    )


@pytest.mark.unit
def test_REE_costing_variableOM_resourcesnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        TypeError,
        match="resources argument must be a dictionary",
    ):
        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            feedstock_rate=m.fs.feedstock,
            resources=[
                m.fs.water,
            ],
        )


# TODO why did this value change?
@pytest.mark.component
def test_REE_costing_variableOM_customprices():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_REE_process_costs(
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        feedstock_rate=m.fs.feedstock,
        resources={
            "water": m.fs.water,
        },
        resource_prices={"water": 1.90e-3 * 1e-6 * pyunits.MUSD_2021 / pyunits.gallon},
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()
    m.fs.costing.variable_operating_costs.display()
    c = m.fs.costing
    t = 0
    resources = {
        "water": m.fs.water,
    }
    print(value(sum(c.variable_operating_costs[t, r] for r in resources)))
    print(
        value(
            c.plant_overhead_cost[t]
            if hasattr(c, "plant_overhead_cost")
            else 0 * c.CEPCI_units / pyunits.year
        )
    )
    print(
        value(
            c.land_cost
            if c.land_cost_reoccurrence == "annual"
            else 0 * c.CEPCI_units / pyunits.year
        )
    )
    print(
        value(
            c.additional_chemicals_cost
            if c.additional_chemicals_cost_reoccurrence == "annual"
            else 0 * c.CEPCI_units / pyunits.year
        )
    )
    print(
        value(
            c.additional_waste_cost
            if c.additional_waste_cost_reoccurrence == "annual"
            else 0 * c.CEPCI_units / pyunits.year
        )
    )
    print(
        value(
            c.maintenance_material_cost
            if hasattr(c, "maintenance_material_cost")
            else 0 * c.CEPCI_units / pyunits.year
        )
    )
    print(value(c.other_variable_costs[t]))
    print(value(c.custom_variable_costs))

    # check some cost results
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.23406, rel=1e-4
    )


@pytest.mark.unit
def test_REE_costing_variableOM_custompricesnotadict():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    with pytest.raises(
        TypeError,
        match="Dictionary of resource prices must be a dict object",
    ):
        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            feedstock_rate=m.fs.feedstock,
            resources={
                "water": m.fs.water,
            },
            resource_prices=[1e-3 * 1e-6 * pyunits.MUSD_2021 / pyunits.gallon],
        )


@pytest.mark.unit
def test_REE_costing_variableOM_resourcenotinpricelist():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        has_fixed_OM=True,
        has_variable_OM=True,
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.H2O = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.H2O.fix()

    with pytest.raises(
        AttributeError,
        match="A resource was included that does not contain a price. Prices "
        "exist for the following resources: \\['natural_gas', 'coal', 'water', "
        "'water_treatment_chemicals', 'ammonia', 'SCR_catalyst', 'triethylene_glycol', "
        "'SCR_catalyst_waste', 'triethylene_glycol_waste', 'amine_purification_unit_waste', "
        "'thermal_reclaimer_unit_waste', 'PSA_adsorbent', 'amine_entrainer', "
        "'electricity', 'reactive_membrane_replacement', 'PSA_adsorbent_waste_disposal', "
        "'nonharzardous_waste_disposal', 'cooling_water', 'chilled_water', "
        "'lp_heating_steam', 'mp_heating_steam', 'hp_heating_steam', 'sulfolane_entrainer', "
        "'ZnO_sulfur_guard_catalyst', 'prereformer_catalyst', 'reformer_catalyst', "
        "'ZnO_sulfur_guard_catalyst_waste_disposal', 'prereformer_catalyst_waste_disposal', "
        "'reformer_catalyst_waste_disposal', 'power', 'diesel', 'bioleaching_solution', "
        "'H2SO4', 'polymer', 'NAOH', 'CACO3', 'coal_calcite', 'HCl', 'oxalic_acid', "
        "'ascorbic_acid', 'kerosene', 'D2EHPA', 'NA2S', 'nonhazardous_solid_waste', "
        "'nonhazardous_precipitate_waste', 'dust_and_volatiles'\\]",
    ):
        m.fs.costing.build_REE_process_costs(
            pure_product_output_rates={
                "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
            },
            mixed_product_output_rates={
                "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
            },
            feedstock_rate=m.fs.feedstock,
            resources={
                "H2O": m.fs.H2O,
            },
        )


@pytest.mark.component
def test_REE_costing_economy_of_numbers():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        Lang_factor=1,
        has_economy_of_numbers=True,
        CEPCI_year="2007",
    )

    m.fs.costing.cum_num_units.set_value(5)
    m.fs.costing.learning_rate.set_value(0.05)

    m.fs.costing.build_REE_process_costs(
        total_purchase_cost=41.2 * m.fs.costing.CEPCI_units,
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check that the cost units are as expected
    assert pyunits.get_units(m.fs.costing.total_TPC) == pyunits.MUSD_2007
    # check that some objects are built as expected
    assert hasattr(m.fs.costing, "total_BEC")
    assert hasattr(m.fs.costing, "total_BEC_eq")
    assert hasattr(m.fs.costing, "total_TPC")
    assert hasattr(m.fs.costing, "total_overnight_capital")
    # check some results
    assert value(m.fs.costing.total_BEC) == pytest.approx(41.2, rel=1e-4)
    assert value(m.fs.costing.total_TPC) == pytest.approx(36.574, rel=1e-4)


@pytest.mark.component
def test_REE_costing_consider_taxes():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = REECosting(
        # arguments for NPV
        Lang_factor=2.97,
        has_fixed_OM=True,
        has_variable_OM=True,
        has_taxes_and_credits=True,
        has_net_present_value=True,
        has_capital_expenditure_period=True,
        capital_expenditure_percentages=[10, 60, 30],
    )

    # 1.3 is CS Jaw Crusher
    CS_jaw_crusher_accounts = ["1.3"]
    m.fs.CS_jaw_crusher = UnitModelBlock()
    m.fs.CS_jaw_crusher.power = pyo.Var(initialize=589, units=pyunits.hp)
    m.fs.CS_jaw_crusher.power.fix()
    m.fs.CS_jaw_crusher.costing = REEUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=REECostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": CS_jaw_crusher_accounts,
            "scaled_param": m.fs.CS_jaw_crusher.power,
        },
    )

    m.fs.feedstock = pyo.Var(initialize=500, units=pyunits.ton / pyunits.hr)
    m.fs.feedstock.fix()

    m.fs.water = pyo.Var(m.fs.time, initialize=1000, units=pyunits.gallon / pyunits.hr)
    m.fs.water.fix()

    m.fs.costing.build_REE_process_costs(
        pure_product_output_rates={
            "Sc2O3": 1.9 * pyunits.kg / pyunits.hr,
        },
        mixed_product_output_rates={
            "Sc2O3": 0.00143 * pyunits.kg / pyunits.hr,
        },
        feedstock_rate=m.fs.feedstock,
        resources={
            "water": m.fs.water,
        },
    )

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    REECostingData.initialize(m.fs.costing)

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    # check some cost results
    assert value(m.fs.costing.total_fixed_OM_cost) == pytest.approx(6.00109, rel=1e-4)
    assert value(m.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
        1.234056, rel=1e-4
    )
    assert value(m.fs.costing.plant_overhead_cost[0]) == pytest.approx(
        1.20022, rel=1e-4
    )
    assert value(m.fs.costing.other_variable_costs[0]) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.land_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.additional_chemicals_cost) == pytest.approx(
        0.0000, abs=1e-4
    )
    assert value(m.fs.costing.additional_waste_cost) == pytest.approx(0.0000, abs=1e-4)
    assert value(m.fs.costing.income_tax) == pytest.approx(5.790603, abs=1e-4)
    assert value(m.fs.costing.net_tax_owed) == pytest.approx(3.066235, abs=1e-4)
    assert value(m.fs.costing.pv_taxes) == pytest.approx(-19.612760, abs=1e-4)
    assert value(m.fs.costing.npv) == pytest.approx(171.874896, abs=1e-4)
