#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

import pytest

from prommis.nanofiltration.membrane_cascade_flowsheet.membrane import Membrane
from prommis.nanofiltration.membrane_cascade_flowsheet.solute_property import (
    SoluteParameters,
)

# -----------------------------------------------------------------------------
# Test settings

# Get default solver for testing
solver = get_solver()

# solutes
solutes = ["Li", "Co"]

# membrane performance
flux = 0.1
sieving_coefficient = {"Li": 1.3, "Co": 0.5}


@pytest.fixture(scope="module")
def membrane():
    """Construct a flowsheet with a membrane unit."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteParameters(solutes=solutes)

    m.fs.unit = Membrane(
        number_of_finite_elements=1,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        streams={
            "retentate": {
                "property_package": m.fs.properties,
                "side_streams": [1],
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            "permeate": {
                "property_package": m.fs.properties,
                "has_feed": False,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
        },
    )

    # inlet flow
    m.fs.unit.retentate_inlet.flow_vol[0].fix(100)
    for sol in solutes:
        m.fs.unit.retentate_inlet.flow_mass_solute[0, sol].fix(100)

    # side stream flow
    m.fs.unit.retentate_side_stream_1.flow_vol[0].fix(100)
    for sol in solutes:
        m.fs.unit.retentate_side_stream_1.flow_mass_solute[0, sol].fix(100)

    # membrane length [m^2]
    m.fs.unit.length.fix(100)

    return m


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config(membrane):
    # extends MSContactor with additional flux and sieving args
    assert len(membrane.fs.unit.config) == 9

    # MSContactor model, physical separation (no reactions)
    assert not membrane.fs.unit.config.dynamic
    assert not membrane.fs.unit.config.has_holdup
    assert not membrane.fs.unit.config.interacting_streams
    assert not membrane.fs.unit.config.heterogeneous_reactions
    for stream in membrane.fs.unit.config.streams:
        assert (
            membrane.fs.unit.config.streams[stream]["property_package"]
            is membrane.fs.properties
        )
        assert not membrane.fs.unit.config.streams[stream]["reaction_package"]
        assert (
            membrane.fs.unit.config.streams[stream]["flow_direction"]
            is FlowDirection.forward
        )
        assert not membrane.fs.unit.config.streams[stream]["has_rate_reactions"]
        assert not membrane.fs.unit.config.streams[stream]["has_equilibrium_reactions"]
        assert not membrane.fs.unit.config.streams[stream]["has_energy_balance"]
        assert not membrane.fs.unit.config.streams[stream]["has_heat_transfer"]
        assert not membrane.fs.unit.config.streams[stream]["has_heat_of_reaction"]
        assert not membrane.fs.unit.config.streams[stream]["has_pressure_balance"]
        assert not membrane.fs.unit.config.streams[stream]["has_pressure_change"]


# -----------------------------------------------------------------------------
class TestMembrane(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, membrane):
        # inlet/outlet flows
        assert hasattr(membrane.fs.unit, "retentate_inlet")
        assert len(membrane.fs.unit.retentate_inlet.vars) == 2
        assert hasattr(membrane.fs.unit.retentate_inlet, "flow_vol")
        assert hasattr(membrane.fs.unit.retentate_inlet, "flow_mass_solute")

        assert hasattr(membrane.fs.unit, "retentate_side_stream_1")
        assert len(membrane.fs.unit.retentate_side_stream_1.vars) == 2
        assert hasattr(membrane.fs.unit.retentate_side_stream_1, "flow_vol")
        assert hasattr(membrane.fs.unit.retentate_side_stream_1, "flow_mass_solute")

        assert hasattr(membrane.fs.unit, "retentate_outlet")
        assert len(membrane.fs.unit.retentate_outlet.vars) == 2
        assert hasattr(membrane.fs.unit.retentate_outlet, "flow_vol")
        assert hasattr(membrane.fs.unit.retentate_outlet, "flow_mass_solute")

        assert hasattr(membrane.fs.unit, "permeate_outlet")
        assert len(membrane.fs.unit.permeate_outlet.vars) == 2
        assert hasattr(membrane.fs.unit.permeate_outlet, "flow_vol")
        assert hasattr(membrane.fs.unit.permeate_outlet, "flow_mass_solute")

        # auxiliary variables for exponential reformulation
        assert hasattr(membrane.fs.unit, "LN_F_in")
        assert hasattr(membrane.fs.unit, "LN_F_out")
        assert hasattr(membrane.fs.unit, "LN_M_in")
        assert hasattr(membrane.fs.unit, "LN_M_out")
        assert hasattr(membrane.fs.unit, "LN_F_in_exp")
        assert hasattr(membrane.fs.unit, "LN_F_out_exp")
        assert hasattr(membrane.fs.unit, "LN_M_in_exp")
        assert hasattr(membrane.fs.unit, "LN_M_out_exp")

        # performance equations: linearized solute sieving, flux
        assert hasattr(membrane.fs.unit, "solute_rule_lin")
        assert hasattr(membrane.fs.unit, "solvent_rule")

        # material balances
        assert hasattr(membrane.fs.unit, "permeate_material_balance")
        assert hasattr(membrane.fs.unit, "retentate_material_balance")

        # membrane stage upper bounds
        assert hasattr(membrane.fs.unit, "conc_limits")
        assert hasattr(membrane.fs.unit, "flow_limits")

        assert number_variables(membrane.fs.unit) == 26
        assert number_total_constraints(membrane.fs.unit) == 21
        assert number_unused_variables(membrane.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, membrane):
        assert_units_consistent(membrane.fs.unit)

    @pytest.mark.component
    def test_dof(self, membrane):
        dt = DiagnosticsToolbox(model=membrane)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, membrane):
        initializer = BlockTriangularizationInitializer()
        initializer.initialize(membrane.fs.unit)
        assert (
            initializer.summary[membrane.fs.unit]["status"] == InitializationStatus.Ok
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, membrane):
        solver = get_solver()
        results = solver.solve(membrane)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_numerical_issues(self, membrane):
        dt = DiagnosticsToolbox(model=membrane)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, membrane):
        # retentate outlet
        membrane.fs.pprint()
        assert pytest.approx(190, abs=1e-6) == value(
            membrane.fs.unit.retentate_outlet.flow_vol[0]
        )
        assert pytest.approx(194.935, abs=1e-3) == value(
            membrane.fs.unit.retentate_outlet.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(187.098, abs=1e-3) == value(
            membrane.fs.unit.retentate_outlet.flow_mass_solute[0, "Li"]
        )
        # permeate outlet
        assert pytest.approx(10, abs=1e-6) == value(
            membrane.fs.unit.permeate_outlet.flow_vol[0]
        )
        assert pytest.approx(5.064, abs=1e-3) == value(
            membrane.fs.unit.permeate_outlet.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(12.901, abs=1e-3) == value(
            membrane.fs.unit.permeate_outlet.flow_mass_solute[0, "Li"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, membrane):
        # this model only uses mass balances
        # solvent mass balance
        assert value(
            (
                membrane.fs.unit.retentate_inlet.flow_vol[0]
                + membrane.fs.unit.retentate_side_stream_1.flow_vol[0]
            )
            * membrane.fs.properties.dens_H2O
        ) == pytest.approx(
            value(
                (
                    membrane.fs.unit.retentate_outlet.flow_vol[0]
                    + membrane.fs.unit.permeate_outlet.flow_vol[0]
                )
                * membrane.fs.properties.dens_H2O
            ),
            rel=1e-6,
            abs=1e-6,
        )

        # solute mass balances
        for sol in solutes:
            assert value(
                membrane.fs.unit.retentate_inlet.flow_mass_solute[0, sol]
                + membrane.fs.unit.retentate_side_stream_1.flow_mass_solute[0, sol]
            ) == pytest.approx(
                value(
                    membrane.fs.unit.retentate_outlet.flow_mass_solute[0, sol]
                    + membrane.fs.unit.permeate_outlet.flow_mass_solute[0, sol]
                ),
                rel=1e-6,
                abs=1e-6,
            )
