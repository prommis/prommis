#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, MaterialBalanceType, MomentumBalanceType
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
from idaes.models.unit_models import EnergySplittingType

import pytest

from prommis.nanofiltration.membrane_cascade_flowsheet.precipitator import Precipitator
from prommis.nanofiltration.membrane_cascade_flowsheet.solute_property import (
    SoluteParameters,
)

# -----------------------------------------------------------------------------
# Test settings

# Get default solver for testing
solver = get_solver()

# solutes
solutes = ["Li", "Co"]

# yields
yields = {"Li": 0.81, "Co": 0.05}


@pytest.fixture(scope="module")
def precip():
    """Construct a flowsheet with a precipitator unit."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteParameters(solutes=solutes)

    m.fs.unit = Precipitator(
        outlet_list=["solid", "recycle"],
        yields=yields,
        property_package=m.fs.properties,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    # inlet flow
    m.fs.unit.inlet.flow_vol[0].fix(100)
    for sol in solutes:
        m.fs.unit.inlet.flow_mass_solute[0, sol].fix(100)

    # precipitator volume [m^3]
    m.fs.unit.volume.fix(100)

    return m


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config(precip):
    # same config as idaes.models.unit_models.separator with added "yields" arg
    assert len(precip.fs.unit.config) == 16

    # yield-based model with only material balances
    assert not precip.fs.unit.config.dynamic
    assert not precip.fs.unit.config.has_holdup
    assert not precip.fs.unit.config.has_phase_equilibrium
    assert precip.fs.unit.config.property_package is precip.fs.properties
    assert (
        precip.fs.unit.config.material_balance_type
        is MaterialBalanceType.componentTotal
    )
    assert precip.fs.unit.config.momentum_balance_type is MomentumBalanceType.none
    assert precip.fs.unit.config.energy_split_basis is EnergySplittingType.none


# -----------------------------------------------------------------------------
class TestPrecip(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, precip):
        # inlet/outlet flows
        assert hasattr(precip.fs.unit, "inlet")
        assert len(precip.fs.unit.inlet.vars) == 2
        assert hasattr(precip.fs.unit.inlet, "flow_vol")
        assert hasattr(precip.fs.unit.inlet, "flow_mass_solute")

        assert hasattr(precip.fs.unit, "solid")
        assert len(precip.fs.unit.solid.vars) == 2
        assert hasattr(precip.fs.unit.solid, "flow_vol")
        assert hasattr(precip.fs.unit.solid, "flow_mass_solute")

        assert hasattr(precip.fs.unit, "recycle")
        assert len(precip.fs.unit.recycle.vars) == 2
        assert hasattr(precip.fs.unit.recycle, "flow_vol")
        assert hasattr(precip.fs.unit.recycle, "flow_mass_solute")

        # yield constraints and residence time
        assert hasattr(precip.fs.unit, "outlet_yield_eqn")
        assert hasattr(precip.fs.unit, "prec_res_time")
        assert hasattr(precip.fs.unit, "prec_vol")
        assert hasattr(precip.fs.unit, "yields_eqn")

        assert number_variables(precip.fs.unit) == 19
        assert number_total_constraints(precip.fs.unit) == 19
        assert number_unused_variables(precip.fs.unit) == 2

    @pytest.mark.component
    def test_units(self, precip):
        assert_units_consistent(precip.fs.unit)

    @pytest.mark.component
    def test_dof(self, precip):
        dt = DiagnosticsToolbox(model=precip)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, precip):
        initializer = BlockTriangularizationInitializer()
        initializer.initialize(precip.fs.unit)
        assert initializer.summary[precip.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, precip):
        solver = get_solver()
        results = solver.solve(precip)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_numerical_issues(self, precip):
        dt = DiagnosticsToolbox(model=precip)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, precip):
        # note that the maximum yield is not reached as the volume is not extremely large
        # aqueous outlet
        assert pytest.approx(100, abs=1e-6) == value(precip.fs.unit.recycle.flow_vol[0])
        assert pytest.approx(95.05, abs=1e-3) == value(
            precip.fs.unit.recycle.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(19.814, abs=1e-3) == value(
            precip.fs.unit.recycle.flow_mass_solute[0, "Li"]
        )
        # solid outlet
        assert pytest.approx(0, abs=1e-6) == value(precip.fs.unit.solid.flow_vol[0])
        assert pytest.approx(4.949, abs=1e-3) == value(
            precip.fs.unit.solid.flow_mass_solute[0, "Co"]
        )
        assert pytest.approx(80.185, abs=1e-3) == value(
            precip.fs.unit.solid.flow_mass_solute[0, "Li"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, precip):
        # this model only uses mass balances
        # solvent mass balance
        assert value(
            precip.fs.unit.inlet.flow_vol[0] * precip.fs.properties.dens_H2O
        ) == pytest.approx(
            value(precip.fs.unit.recycle.flow_vol[0] * precip.fs.properties.dens_H2O),
            rel=1e-6,
            abs=1e-6,
        )

        # solute mass balances
        for sol in solutes:
            assert value(
                precip.fs.unit.inlet.flow_mass_solute[0, sol]
            ) == pytest.approx(
                value(
                    precip.fs.unit.recycle.flow_mass_solute[0, sol]
                    + precip.fs.unit.solid.flow_mass_solute[0, sol]
                ),
                rel=1e-6,
                abs=1e-6,
            )
