#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re

from pyomo.environ import (
    Objective,
    SolverFactory,
    Suffix,
    TransformationFactory,
    Var,
    assert_optimal_termination,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.nanofiltration.membrane_cascade_flowsheet.diafiltration_flowsheet_model import (
    DiafiltrationModel,
)

# -----------------------------------------------------------------------------
# Test settings

# Get default solver for testing
solver = get_solver()

# solutes
solutes = ["Li", "Co"]

# yields
yields = {
    "permeate": {"Li": 0.81, "Co": 0.01},
    "retentate": {"Li": 0.20, "Co": 0.89},
}

# membrane performance
flux = 0.1
sieving_coefficient = {"Li": 1.3, "Co": 0.5}

# inlet conditions
feed = {
    "solvent": 100,  # m^3/hr of water
    "Li": 1.7 * 100,  # kg/hr
    "Co": 17 * 100,  # kg/hr
}
diaf = {
    "solvent": 30,  # m^3/hr of water
    "Li": 0.1 * 30,  # kg/hr
    "Co": 0.2 * 30,  # kg/hr
}

# precipitators
use_precipitators = True

# model sizes
sizes = [(1, 3), (2, 5), (3, 10)]

# superstructure configuration
mixing = ["stage", "tube"]


@pytest.fixture(scope="module")
def flowsheet():
    """Construct a flowsheet for a membrane cascade with diafiltration and downstream precipitators."""
    flowsheet_setup = DiafiltrationModel(
        NS=1,
        NT=3,
        solutes=solutes,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        feed=feed,
        diafiltrate=diaf,
        precipitate=use_precipitators,
        precipitate_yield=yields,
    )

    return flowsheet_setup


# -----------------------------------------------------------------------------
class TestFlowsheet(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, flowsheet):
        # check all units are constructed across model sizes
        # and superstructure configuration
        for (
            stages,
            tubes,
        ) in sizes:
            for mix in mixing:
                flowsheet.ns = stages
                flowsheet.nt = tubes
                m = flowsheet.build_flowsheet(mixing=mix)

                # sets
                assert len(m.fs.stages) == stages
                assert len(m.fs.tubes) == tubes

                # check each unit model is indexed correctly
                # feed
                assert len(m.fs.split_feed) == 1
                # diafiltrate
                assert len(m.fs.split_diafiltrate) == 1
                # stages
                assert len(m.fs.stage) == stages
                # NT tubes per stage
                for i in m.fs.stages:
                    assert len(m.fs.stage[i].elements) == tubes
                # 2 splitters (retentate/permeate) per stage
                assert len(m.fs.split_retentate) == stages
                assert len(m.fs.split_permeate) == stages
                if mix == "tube":
                    # mixers into all tubes in all stages
                    assert len(m.fs.inlet_mixers) == stages * tubes
                    # recycle splitters into every previous tube inlet
                    # for every stage after the first
                    if stages != 1:
                        assert len(m.fs.recycle_splitters) == stages - 1
                if mix == "stage":
                    # mixers only into each stage
                    assert len(m.fs.inlet_mixers) == stages
                    # splitters for mixed inlet for every tube
                    assert len(m.fs.splitters) == stages
                    # no recycle split needed for stage mixing
                # mix products to send to precipitator
                assert len(m.fs.mix_product_retentate) == 1
                assert len(m.fs.mix_product_permeate) == 1
                # precipitators for each product stream
                assert len(m.fs.precipitator) == 2
                # mix precipitator recycle stream
                assert len(m.fs.mix_precipitate_recycle) == 1
                # split recycle into diafiltrate and waste
                assert len(m.fs.split_precipitate_recycle) == 1

        # recovery lower bounds
        assert hasattr(m, "recovery_li")  # Li recovery LB parameter
        assert hasattr(m, "recovery_co")  # Co recovery LB parameter
        assert hasattr(m, "li_lb")  # Li recovery LB constraint
        assert hasattr(m, "co_lb")  # Co recovery LB constraint
        assert hasattr(m, "prec_li_lb")  # LB constraint for precipitators
        assert hasattr(m, "prec_co_lb")  # LB constraint for precipitators

        # objective
        assert (
            len([obj for obj in m.component_data_objects(Objective, active=True)]) == 1
        )
        # contains alternative objectives for flowsheet with/without precipitators
        assert hasattr(m, "co_obj")
        assert hasattr(m, "li_obj")
        assert hasattr(m, "prec_co_obj")
        assert hasattr(m, "prec_li_obj")

    @pytest.mark.component
    def test_units(self, flowsheet):
        flowsheet.ns = 1
        flowsheet.nt = 3
        m = flowsheet.build_flowsheet(mixing="stage")
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof(self, flowsheet):
        # check dof across model sizes and superstructure configuration
        for (
            stages,
            tubes,
        ) in sizes:
            for mix in mixing:
                flowsheet.ns = stages
                flowsheet.nt = tubes
                m = flowsheet.build_flowsheet(mixing=mix)
                flowsheet.unfix_dof(m, mixing=mix, precipitate=use_precipitators)

                # assume fixed diafiltrate flow, precipitator volumes
                m.fs.split_diafiltrate.mixed_state[0].flow_vol.fix(30)
                m.fs.precipitator["retentate"].volume.fix(500)
                m.fs.precipitator["permeate"].volume.fix(500)

                if mix == "tube":
                    # DOF = 3*NS*NT + NS - NT - 2
                    assert (
                        degrees_of_freedom(m) == 3 * stages * tubes + stages - tubes - 2
                    )
                if mix == "stage":
                    # DOF = NS*(NT + 3) - 3
                    assert degrees_of_freedom(m) == stages * (tubes + 3) - 3

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, flowsheet):
        # initialization starts from the feed inlet and propagates state
        # across the entire model, initializing each unit.
        # fixed point iteration is used to converge recycles
        #     - this leads to multiple cycles of propagate state across the flowsheet
        #       until recycles are converged
        flowsheet.ns = 2
        flowsheet.nt = 5
        mix = "stage"
        m = flowsheet.build_flowsheet(mixing=mix)
        flowsheet.initialize(m, mixing=mix, precipitate=use_precipitators, info=False)

        # check converged recycles
        # check internal membrane cascade recycles
        for i in m.fs.stages:
            if i != flowsheet.ns:
                assert value(
                    m.fs.inlet_mixers[i].recycle_state[0].flow_vol
                ) == pytest.approx(
                    value(m.fs.split_retentate[i + 1].recycle_state[0].flow_vol),
                    rel=1e-6,
                    abs=1e-6,
                )
                for sol in solutes:
                    assert value(
                        m.fs.inlet_mixers[i].recycle_state[0].flow_mass_solute[sol]
                    ) == pytest.approx(
                        value(
                            m.fs.split_retentate[i + 1]
                            .recycle_state[0]
                            .flow_mass_solute[sol]
                        ),
                        rel=1e-6,
                        abs=1e-6,
                    )
        # check diafiltrate recycle
        assert value(m.fs.split_diafiltrate.mixed_state[0].flow_vol) == pytest.approx(
            value(m.fs.split_precipitate_recycle.recycle_state[0].flow_vol),
            rel=1e-6,
            abs=1e-6,
        )
        for sol in solutes:
            assert value(
                m.fs.split_diafiltrate.mixed_state[0].flow_mass_solute[sol]
            ) == pytest.approx(
                value(
                    m.fs.split_precipitate_recycle.recycle_state[0].flow_mass_solute[
                        sol
                    ]
                ),
                rel=1e-6,
                abs=1e-6,
            )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, flowsheet):
        flowsheet.ns = 2
        flowsheet.nt = 5
        mix = "stage"
        m = flowsheet.build_flowsheet(mixing=mix)
        flowsheet.initialize(m, mixing=mix, precipitate=use_precipitators, info=False)
        flowsheet.unfix_dof(m, mixing=mix, precipitate=use_precipitators)

        # assume fixed diafiltrate flow, precipitator volumes
        m.fs.split_diafiltrate.mixed_state[0].flow_vol.fix(30)
        m.fs.precipitator["retentate"].volume.fix(500)
        m.fs.precipitator["permeate"].volume.fix(500)

        # set lower bound
        m.recovery_li = 0.8

        solver = get_solver()
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_numerical_issues(self, flowsheet):
        flowsheet.ns = 2
        flowsheet.nt = 5
        mix = "stage"
        m = flowsheet.build_flowsheet(mixing=mix)
        flowsheet.unfix_dof(m, mixing=mix, precipitate=use_precipitators)

        # assume fixed diafiltrate flow, precipitator volumes
        m.fs.split_diafiltrate.mixed_state[0].flow_vol.fix(30)
        m.fs.precipitator["retentate"].volume.fix(500)
        m.fs.precipitator["permeate"].volume.fix(500)

        # set lower bound
        m.recovery_li = 0.8

        solver = get_solver()
        solver.solve(m)

        # expecting one warning: variables at or outside bounds
        # this is due to some unused flows being fixed to their lower bound of 0
        dt = DiagnosticsToolbox(model=m)
        with pytest.raises(
            AssertionError, match=re.escape("Numerical issues found (1).")
        ):
            dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, flowsheet):
        # solve model
        flowsheet.ns = 2
        flowsheet.nt = 5
        mix = "stage"
        m = flowsheet.build_flowsheet(mixing=mix)
        flowsheet.initialize(m, mixing=mix, precipitate=use_precipitators, info=False)
        flowsheet.unfix_dof(m, mixing=mix, precipitate=use_precipitators)

        # assume fixed diafiltrate flow, precipitator volumes
        m.fs.split_diafiltrate.mixed_state[0].flow_vol.fix(30)
        m.fs.precipitator["retentate"].volume.fix(500)
        m.fs.precipitator["permeate"].volume.fix(500)

        # set lower bound
        m.recovery_li = 0.8

        solver = get_solver()
        solver.solve(m)

        # precipitator outlets
        assert pytest.approx(0, abs=1e-6) == value(
            m.fs.precipitator["retentate"].solid.flow_vol[0]
        )
        assert pytest.approx(893.679, abs=1e-3) == value(
            m.fs.precipitator["retentate"].solid.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(2.320, abs=1e-3) == value(
            m.fs.precipitator["retentate"].solid.flow_mass_solute[0, "Li"]
        )
        assert pytest.approx(0, abs=1e-6) == value(
            m.fs.precipitator["permeate"].solid.flow_vol[0]
        )
        assert pytest.approx(9.349, abs=1e-3) == value(
            m.fs.precipitator["permeate"].solid.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(135.999, abs=1e-3) == value(
            m.fs.precipitator["permeate"].solid.flow_mass_solute[0, "Li"]
        )

        # solvent exit outlet
        assert pytest.approx(100, abs=1e-6) == value(
            m.fs.split_precipitate_recycle.waste.flow_vol[0]
        )
        assert pytest.approx(796.971, abs=1e-3) == value(
            m.fs.split_precipitate_recycle.waste.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(31.679, abs=1e-3) == value(
            m.fs.split_precipitate_recycle.waste.flow_mass_solute[0, "Li"]
        )

        # system recoveries
        assert pytest.approx(0.525, abs=1e-3) == value(value(m.prec_perc_co))
        assert pytest.approx(0.799, abs=1e-3) == value(value(m.prec_perc_li))

        # membrane length
        assert pytest.approx(1105.388, abs=1e-3) == value(value(m.fs.stage[1].length))

        # objective
        assert pytest.approx(893.679, abs=1e-3) == value(value(m.prec_co_obj))

        # flows all greater than 0
        for i in m.component_data_objects(Var):
            if "flow_vol" in i.name or "flow_mass_solute" in i.name:
                assert value(i) >= 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, flowsheet):
        # solve model
        flowsheet.ns = 2
        flowsheet.nt = 5
        mix = "stage"
        m = flowsheet.build_flowsheet(mixing=mix)
        flowsheet.initialize(m, mixing=mix, precipitate=use_precipitators, info=False)
        flowsheet.unfix_dof(m, mixing=mix, precipitate=use_precipitators)

        # assume fixed diafiltrate flow, precipitator volumes
        m.fs.split_diafiltrate.mixed_state[0].flow_vol.fix(30)
        m.fs.precipitator["retentate"].volume.fix(500)
        m.fs.precipitator["permeate"].volume.fix(500)

        # set lower bound
        m.recovery_li = 0.8

        solver = get_solver()
        solver.solve(m)

        # this model only uses mass balances
        # solvent mass balance
        assert value(
            m.fs.split_feed.mixed_state[0].flow_vol * m.fs.properties.dens_H2O
        ) == pytest.approx(
            value(
                (
                    m.fs.precipitator["permeate"].solid_state[0].flow_vol
                    + m.fs.precipitator["retentate"].solid_state[0].flow_vol
                    + m.fs.split_precipitate_recycle.waste_state[0].flow_vol
                )
                * m.fs.properties.dens_H2O
            ),
            rel=1e-6,
            abs=1e-6,
        )

        # solute mass balances
        for sol in solutes:
            assert value(
                m.fs.split_feed.mixed_state[0].flow_mass_solute[sol]
            ) == pytest.approx(
                value(
                    m.fs.precipitator["permeate"].solid_state[0].flow_mass_solute[sol]
                    + m.fs.precipitator["retentate"]
                    .solid_state[0]
                    .flow_mass_solute[sol]
                    + m.fs.split_precipitate_recycle.waste_state[0].flow_mass_solute[
                        sol
                    ]
                ),
                rel=1e-6,
                abs=1e-6,
            )


# -----------------------------------------------------------------------------
@pytest.fixture(scope="module")
def flowsheet_with_costing():
    """Construct a flowsheet for a membrane cascade with diafiltration and downstream precipitators."""
    flowsheet_with_costing_setup = DiafiltrationModel(
        NS=3,
        NT=10,
        solutes=solutes,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        feed=feed,
        diafiltrate=diaf,
        precipitate=use_precipitators,
        precipitate_yield=yields,
    )

    return flowsheet_with_costing_setup


class TestFlowsheetCosting(object):
    @pytest.mark.component
    def test_simple_costing(self, flowsheet_with_costing):
        # build a 3x10 stage mixing cascade with simple costing
        mix_style = "stage"
        m_simple = flowsheet_with_costing.build_flowsheet(mixing=mix_style)
        flowsheet_with_costing.initialize(
            m_simple, mixing=mix_style, precipitate=use_precipitators
        )
        flowsheet_with_costing.unfix_dof(
            m_simple, mixing=mix_style, precipitate=use_precipitators
        )
        m_simple.fs.split_diafiltrate.inlet.flow_vol.setub(200)
        flowsheet_with_costing.add_costing(
            m_simple,
            NS=3,
            flux=flux,
            feed=feed,
            diaf=diaf,
            precipitate=use_precipitators,
            atmospheric_pressure=101.325,  # ambient pressure, kPa
            operating_pressure=145,  # nanofiltration operating pressure, psi
            simple_costing=True,
        )
        flowsheet_with_costing.add_costing_objectives(m_simple)
        flowsheet_with_costing.add_costing_scaling(m_simple, NS=3, simple_costing=True)

        assert_units_consistent(m_simple)

        # set recovery lower bounds
        m_simple.recovery_li = 0.8
        m_simple.recovery_co = 0.8

        # check scaling factors
        assert (
            m_simple.scaling_factor[m_simple.fs.costing.aggregate_capital_cost] == 1e-5
        )
        assert m_simple.scaling_factor[
            m_simple.fs.costing.aggregate_fixed_operating_cost
        ] == (1e-4)
        assert (
            m_simple.scaling_factor[
                m_simple.fs.costing.aggregate_variable_operating_cost
            ]
            == 1e-5
        )
        assert m_simple.scaling_factor[m_simple.fs.costing.total_capital_cost] == 1e-5
        assert m_simple.scaling_factor[m_simple.fs.costing.total_operating_cost] == 1e-5
        assert (
            m_simple.scaling_factor[
                m_simple.fs.costing.maintenance_labor_chemical_operating_cost
            ]
            == 1e-4
        )
        assert (
            m_simple.scaling_factor[m_simple.fs.costing.total_annualized_cost] == 1e-5
        )

        for n in range(1, 3 + 1):
            assert (
                m_simple.scaling_factor[m_simple.fs.stage[n].costing.capital_cost]
                == 1e-5
            )
            assert (
                m_simple.scaling_factor[
                    m_simple.fs.stage[n].costing.fixed_operating_cost
                ]
                == 1e-4
            )
            assert (
                m_simple.scaling_factor[m_simple.fs.stage[n].costing.membrane_area]
                == 1e-3
            )

        assert (
            m_simple.scaling_factor[m_simple.fs.feed_pump.costing.capital_cost] == 1e-4
        )
        assert m_simple.scaling_factor[
            m_simple.fs.diafiltrate_pump.costing.capital_cost
        ] == (1e-3)
        assert (
            m_simple.scaling_factor[
                m_simple.fs.diafiltrate_pump.costing.variable_operating_cost
            ]
            == 1e-3
        )

        assert (
            m_simple.scaling_factor[
                m_simple.fs.feed_pump.costing.pump_power_factor_simple
            ]
            == 1e-2
        )
        assert (
            m_simple.scaling_factor[
                m_simple.fs.feed_pump.costing.variable_operating_cost
            ]
            == 1e-4
        )
        assert (
            m_simple.scaling_factor[
                m_simple.fs.diafiltrate_pump.costing.pump_power_factor_simple
            ]
            == 1e-2
        )

        for prod in ["retentate", "permeate"]:
            assert (
                m_simple.scaling_factor[
                    m_simple.fs.precipitator[prod].costing.capital_cost
                ]
                == 1e-5
            )

        # solve scaled model
        scaling = TransformationFactory("core.scale_model")
        solver = SolverFactory("ipopt")
        scaled_model = scaling.create_using(m_simple, rename=False)
        result = solver.solve(scaled_model, tee=True)
        assert_optimal_termination(result)
        scaling.propagate_solution(scaled_model, m_simple)

        dt = DiagnosticsToolbox(m_simple)
        # some flows are at their bounds of zero
        with pytest.raises(
            AssertionError, match=re.escape("Numerical issues found (1).")
        ):
            dt.assert_no_numerical_warnings()

        # check key results
        test_dict = {
            "lithium_recovery": [value(m_simple.prec_perc_li), 0.8],
            "cobalt_recovery": [value(m_simple.prec_perc_co), 0.8],
            "stage_1_area": [value(m_simple.fs.stage[1].length), 863.52],
            "stage_2_area": [value(m_simple.fs.stage[2].length), 863.52],
            "stage_3_area": [value(m_simple.fs.stage[3].length), 863.52],
            "retentate_precipitator_volume": [
                value(m_simple.fs.precipitator["retentate"].volume),
                123.26,
            ],
            "permeate_precipitator_volume": [
                value(m_simple.fs.precipitator["permeate"].volume),
                161.76,
            ],
            "stage_1_capex": [
                value(m_simple.fs.stage[1].costing.capital_cost),
                43176.19,
            ],
            "stage_2_capex": [
                value(m_simple.fs.stage[2].costing.capital_cost),
                43176.19,
            ],
            "stage_3_capex": [
                value(m_simple.fs.stage[3].costing.capital_cost),
                43176.19,
            ],
            "stage_1_opex": [
                value(m_simple.fs.stage[1].costing.fixed_operating_cost),
                8635.24,
            ],
            "stage_2_opex": [
                value(m_simple.fs.stage[2].costing.fixed_operating_cost),
                8635.24,
            ],
            "stage_3_opex": [
                value(m_simple.fs.stage[3].costing.fixed_operating_cost),
                8635.24,
            ],
            "feed_pump_capex": [
                value(m_simple.fs.feed_pump.costing.capital_cost),
                20507.67,
            ],
            "diafiltrate_pump_capex": [
                value(m_simple.fs.diafiltrate_pump.costing.capital_cost),
                6152.30,
            ],
            "feed_pump_opex": [
                value(m_simple.fs.feed_pump.costing.variable_operating_cost),
                28150.33,
            ],
            "diafiltrate_pump_opex": [
                value(m_simple.fs.diafiltrate_pump.costing.variable_operating_cost),
                8445.10,
            ],
            "retentate_precipitator_capex": [
                value(m_simple.fs.precipitator["retentate"].costing.capital_cost),
                122936.46,
            ],
            "permeate_precipitator_capex": [
                value(m_simple.fs.precipitator["permeate"].costing.capital_cost),
                157927.27,
            ],
            "mainentance_and_labor": [
                value(m_simple.fs.costing.maintenance_labor_chemical_operating_cost),
                26223.14,
            ],
            "total_capex": [
                value(m_simple.fs.costing.total_capital_cost),
                874104.54,
            ],
            "total_opex": [
                value(m_simple.fs.costing.total_operating_cost),
                88724.28,
            ],
            "total_annualized_cost": [
                value(m_simple.fs.costing.total_annualized_cost),
                176134.73,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_default_costing(self, flowsheet_with_costing):
        # build a 3x10 stage mixing cascade with default costing
        mix_style = "stage"
        m_default = flowsheet_with_costing.build_flowsheet(mixing=mix_style)
        flowsheet_with_costing.initialize(
            m_default, mixing=mix_style, precipitate=use_precipitators
        )
        flowsheet_with_costing.unfix_dof(
            m_default, mixing=mix_style, precipitate=use_precipitators
        )
        m_default.fs.split_diafiltrate.inlet.flow_vol.setub(200)
        flowsheet_with_costing.add_costing(
            m_default,
            NS=3,
            flux=flux,
            feed=feed,
            diaf=diaf,
            precipitate=use_precipitators,
            atmospheric_pressure=101.325,  # ambient pressure, kPa
            operating_pressure=145,  # nanofiltration operating pressure, psi
            simple_costing=False,
        )
        flowsheet_with_costing.add_costing_objectives(m_default)
        flowsheet_with_costing.add_costing_scaling(
            m_default, NS=3, simple_costing=False
        )

        assert_units_consistent(m_default)

        # set recovery lower bounds
        m_default.recovery_li = 0.8
        m_default.recovery_co = 0.8

        # check scaling
        assert (
            m_default.scaling_factor[m_default.fs.costing.aggregate_capital_cost]
            == 1e-5
        )
        assert (
            m_default.scaling_factor[
                m_default.fs.costing.aggregate_fixed_operating_cost
            ]
            == 1e-4
        )
        assert (
            m_default.scaling_factor[
                m_default.fs.costing.aggregate_variable_operating_cost
            ]
            == 1e-5
        )
        assert m_default.scaling_factor[m_default.fs.costing.total_capital_cost] == 1e-5
        assert (
            m_default.scaling_factor[m_default.fs.costing.total_operating_cost] == 1e-5
        )
        assert (
            m_default.scaling_factor[
                m_default.fs.costing.maintenance_labor_chemical_operating_cost
            ]
            == 1e-4
        )
        assert (
            m_default.scaling_factor[m_default.fs.costing.total_annualized_cost] == 1e-5
        )

        for n in range(1, 3 + 1):
            assert (
                m_default.scaling_factor[m_default.fs.stage[n].costing.capital_cost]
                == 1e-5
            )
            assert (
                m_default.scaling_factor[
                    m_default.fs.stage[n].costing.fixed_operating_cost
                ]
                == 1e-4
            )
            assert (
                m_default.scaling_factor[m_default.fs.stage[n].costing.membrane_area]
                == 1e-3
            )

        assert (
            m_default.scaling_factor[m_default.fs.feed_pump.costing.capital_cost]
            == 1e-4
        )
        assert m_default.scaling_factor[
            m_default.fs.diafiltrate_pump.costing.capital_cost
        ] == (1e-3)
        assert (
            m_default.scaling_factor[
                m_default.fs.diafiltrate_pump.costing.variable_operating_cost
            ]
            == 1e-3
        )

        assert (
            m_default.scaling_factor[
                m_default.fs.cascade.costing.variable_operating_cost
            ]
            == 1e-5
        )
        assert (
            m_default.scaling_factor[m_default.fs.cascade.costing.pressure_drop] == 1e-2
        )
        assert (
            m_default.scaling_factor[
                m_default.fs.feed_pump.costing.variable_operating_cost
            ]
            == 1e-4
        )
        assert m_default.scaling_factor[m_default.fs.feed_pump.costing.pump_head] == 1e6
        assert (
            m_default.scaling_factor[m_default.fs.feed_pump.costing.pump_power] == 1e2
        )
        assert (
            m_default.scaling_factor[m_default.fs.diafiltrate_pump.costing.pump_head]
            == 1e-2
        )
        assert m_default.scaling_factor[
            m_default.fs.diafiltrate_pump.costing.pump_power
        ] == (1e-5)

        for prod in ["retentate", "permeate"]:
            assert (
                m_default.scaling_factor[
                    m_default.fs.precipitator[prod].costing.capital_cost
                ]
                == 1e-5
            )
            assert (
                m_default.scaling_factor[
                    m_default.fs.precipitator[prod].costing.base_cost_per_unit
                ]
                == 1e-4
            )

        # solve scaled model
        scaling = TransformationFactory("core.scale_model")
        solver = SolverFactory("ipopt")
        scaled_model = scaling.create_using(m_default, rename=False)
        result = solver.solve(scaled_model, tee=True)
        assert_optimal_termination(result)
        scaling.propagate_solution(scaled_model, m_default)

        dt = DiagnosticsToolbox(m_default)
        # some flows are at their bounds of zero
        with pytest.raises(
            AssertionError, match=re.escape("Numerical issues found (1).")
        ):
            dt.assert_no_numerical_warnings()

        # check key results
        test_dict = {
            "lithium_recovery": [value(m_default.prec_perc_li), 0.8],
            "cobalt_recovery": [value(m_default.prec_perc_co), 0.8],
            "stage_1_area": [value(m_default.fs.stage[1].length), 1106.92],
            "stage_2_area": [value(m_default.fs.stage[2].length), 1106.92],
            "stage_3_area": [value(m_default.fs.stage[3].length), 1106.92],
            "retentate_precipitator_volume": [
                value(m_default.fs.precipitator["retentate"].volume),
                88.733,
            ],
            "permeate_precipitator_volume": [
                value(m_default.fs.precipitator["permeate"].volume),
                112.44,
            ],
            "stage_1_capex": [
                value(m_default.fs.stage[1].costing.capital_cost),
                55346.07,
            ],
            "stage_2_capex": [
                value(m_default.fs.stage[2].costing.capital_cost),
                55346.07,
            ],
            "stage_3_capex": [
                value(m_default.fs.stage[3].costing.capital_cost),
                55346.07,
            ],
            "stage_1_opex": [
                value(m_default.fs.stage[1].costing.fixed_operating_cost),
                11069.21,
            ],
            "stage_2_opex": [
                value(m_default.fs.stage[2].costing.fixed_operating_cost),
                11069.21,
            ],
            "stage_3_opex": [
                value(m_default.fs.stage[3].costing.fixed_operating_cost),
                11069.21,
            ],
            "casacde_opex": [
                value(m_default.fs.cascade.costing.variable_operating_cost),
                114445.00,
            ],
            "feed_pump_capex": [
                value(m_default.fs.feed_pump.costing.capital_cost),
                42144.68,
            ],
            "diafiltrate_pump_capex": [
                value(m_default.fs.diafiltrate_pump.costing.capital_cost),
                26352.39,
            ],
            "feed_pump_opex": [
                value(m_default.fs.feed_pump.costing.variable_operating_cost),
                0.0033866,
            ],
            "diafiltrate_pump_opex": [
                value(m_default.fs.diafiltrate_pump.costing.variable_operating_cost),
                14731.84,
            ],
            "retentate_precipitator_capex": [
                value(m_default.fs.precipitator["retentate"].costing.capital_cost),
                209492.90,
            ],
            "permeate_precipitator_capex": [
                value(m_default.fs.precipitator["permeate"].costing.capital_cost),
                251679.20,
            ],
            "mainentance_and_labor": [
                value(m_default.fs.costing.maintenance_labor_chemical_operating_cost),
                53253.30,
            ],
            "total_capex": [
                value(m_default.fs.costing.total_capital_cost),
                1775109.96,
            ],
            "total_opex": [
                value(m_default.fs.costing.total_operating_cost),
                215637.78,
            ],
            "total_annualized_cost": [
                value(m_default.fs.costing.total_annualized_cost),
                393148.78,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)
