#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import re

from pyomo.environ import Objective, Var, assert_optimal_termination, value
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
