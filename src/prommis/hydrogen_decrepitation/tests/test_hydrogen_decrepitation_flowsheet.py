#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    SolverFactory,
    TransformationFactory,
    Var,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core.util.math import smooth_max
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)

from idaes.models.unit_models import Feed

import pytest

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_furnace import (
    REPMHydrogenDecrepitationFurnace,
)

from prommis.hydrogen_decrepitation.hydrogen_decrepitation_flowsheet import main


@pytest.fixture(scope="module")
def model():
    m = main()
    
    return m


@pytest.mark.unit
def test_build(model):
    assert hasattr(model.fs, "shredder")
    assert isinstance(
        model.fs.shredder, Feed
    )
    assert hasattr(model.fs, "hydrogen_decrepitation_furnace")
    assert isinstance(
        model.fs.hydrogen_decrepitation_furnace, REPMHydrogenDecrepitationFurnace
    )
    assert hasattr(model.fs, "plant_basis_year")
    assert isinstance(
        model.fs.plant_basis_year, Param
    )
    assert hasattr(model.fs, "flow_2_5inch_HDDs")
    assert isinstance(
        model.fs.flow_2_5inch_HDDs, Var
    )
    assert hasattr(model.fs, "flow_3_5inch_HDDs")
    assert isinstance(
        model.fs.flow_3_5inch_HDDs, Var
    )
    assert hasattr(model.fs.shredder, "HDD_to_REPM_conversion_constraint")
    assert isinstance(
        model.fs.shredder.HDD_to_REPM_conversion_constraint, Constraint
    )
    assert hasattr(model.fs.hydrogen_decrepitation_furnace, "flow_mol_gas_constraint")
    assert isinstance(
        model.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint, Constraint
    )

    assert number_variables(model.fs) == 68
    assert number_total_constraints(model.fs) == 43
    assert_units_consistent(model.fs)


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve(model):
    initializer = BlockTriangularizationInitializer()

    initializer = BlockTriangularizationInitializer()
    initializer.initialize(model.fs.shredder)
    propagate_state(model.fs.shredded_REPM)
    
    model.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.deactivate()  # flow mol will be fixed by initializer
    model.fs.hydrogen_decrepitation_furnace.solid_in[
        0
    ].sum_mass_frac.deactivate()  # mass frac will be fixed by initializer
    
    initializer.initialize(model.fs.hydrogen_decrepitation_furnace)
    
    model.fs.hydrogen_decrepitation_furnace.flow_mol_gas_constraint.activate()
    model.fs.hydrogen_decrepitation_furnace.solid_in[0].sum_mass_frac.activate()
    
    assert (
        initializer.summary[model.fs.hydrogen_decrepitation_furnace]["status"]
        == InitializationStatus.Ok
    )

    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=True)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    model.fs.hydrogen_decrepitation_furnace.report()
    assert value(model.fs.hydrogen_decrepitation_furnace.solid_out[0].flow_mass) == pytest.approx(0.0057367, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd"]) == pytest.approx(0.0100001, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.solid_out[0].mass_frac_comp["Nd2Fe14B"]) == pytest.approx(0.99, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].flow_mol) == pytest.approx(0.0056509, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].mole_frac_comp["H2"]) == pytest.approx(1.0000, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].temperature) == pytest.approx(443.15, rel=1e-5)
    assert value(model.fs.hydrogen_decrepitation_furnace.gas_out[0].pressure) == pytest.approx(1.01325e5, rel=1e-5)
