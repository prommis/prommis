from pyomo.environ import (
    SolverFactory,
    assert_optimal_termination,
    value,
    ConcreteModel,
    Var,
    units,
    TransformationFactory,
)
from pyomo.dae.flatten import flatten_dae_components

from idaes.core.util import DiagnosticsToolbox
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core import FlowDirection, FlowsheetBlock
from idaes.core.util import from_json

import pytest

import numpy as np

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    time_duration = 24

    m.fs = FlowsheetBlock(
        dynamic=True, time_set=[0, time_duration], time_units=units.hour
    )

    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()

    number_of_stages = 3

    m.fs.solex = SolventExtraction(
        number_of_finite_elements=number_of_stages,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
    )

    """
    Discretization of the time domain, and specification of the partition coefficients,
    volume, volume fractions, and the initial conditions of state variables for the components
    for all the stages.

    """

    m.discretizer = TransformationFactory("dae.collocation")
    m.discretizer.apply_to(m, nfe=3, ncp=2, wrt=m.fs.time, scheme="LAGRANGE-RADAU")

    """
    Initialization of the model, which gives a good starting point.

    """

    from_json(m, fname="solvent_extraction.json")

    def copy_first_steady_state(m):
        # Function that propagates initial steady state guess to future time points
        # regular_vars
        regular_vars, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
        # Copy initial conditions forward
        for var in time_vars:
            for t in m.fs.time:
                if t == m.fs.time.first():
                    continue
                else:
                    var[t].value = var[m.fs.time.first()].value

    copy_first_steady_state(m)

    """
    Specifications of the partition coefficients, volume and volume fractions for all
    the stages.

    """

    m.fs.solex.mscontactor.volume[:].fix(0.4)

    m.fs.solex.mscontactor.volume_frac_stream[:, :, "organic"].fix(0.4)

    number_of_stages = 3
    stage_number = np.arange(1, number_of_stages + 1)

    for s in stage_number:
        if s == 1:
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 5.2 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 3 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 24.7 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 99.1 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = 99.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 32.4 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = 58.2 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 58.2 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = 87.6 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = 99.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = 69.8 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = 96.6 / 100
        else:
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Al"] = 4.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ca"] = 12.3 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Fe"] = 6.4 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sc"] = 16.7 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Y"] = 99.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "La"] = 23.2 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Ce"] = 24.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Pr"] = 15.1 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Nd"] = 99.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Sm"] = 99.9 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Gd"] = 7.6 / 100
            m.fs.solex.partition_coefficient[s, "aqueous", "organic", "Dy"] = 5 / 100

    """
    Fixation of the inlet conditions and the initial state values for all the components.

    """

    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H2O"].fix(1e6)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["H"].fix(1.755)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["SO4"].fix(3999.818)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["HSO4"].fix(693.459)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Al"].fix(422.375)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ca"].fix(109.542)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Cl"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Fe"].fix(688.266)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sc"].fix(0.032)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Y"].fix(0.124)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["La"].fix(0.986)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Ce"].fix(2.277)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Pr"].fix(0.303)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Nd"].fix(0.946)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Sm"].fix(0.097)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Gd"].fix(0.2584)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp["Dy"].fix(0.047)

    m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al"].fix(1.267e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca"].fix(2.684e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe"].fix(2.873e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc"].fix(1.734)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y"].fix(2.179e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La"].fix(0.000105)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce"].fix(0.00031)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr"].fix(3.711e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd"].fix(0.000165)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm"].fix(1.701e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd"].fix(3.357e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy"].fix(8.008e-6)

    m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Cl"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(1e-7)

    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["SO4"].fix(3999.818)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["HSO4"].fix(693.459)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Al"].fix(422.375)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ca"].fix(109.542)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Cl"].fix(1e-7)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Fe"].fix(688.266)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sc"].fix(0.032)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Y"].fix(0.124)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["La"].fix(0.986)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Ce"].fix(2.277)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Pr"].fix(0.303)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Nd"].fix(0.946)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Sm"].fix(0.097)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Gd"].fix(0.2584)
    m.fs.solex.mscontactor.aqueous[0, :].conc_mass_comp["Dy"].fix(0.047)

    m.fs.solex.mscontactor.aqueous_inherent_reaction_extent[0, :, "Ka2"].fix(0)
    m.fs.solex.mscontactor.aqueous[0, :].flow_vol.fix(62.01)

    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Al"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ca"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Fe"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sc"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Y"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["La"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Ce"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Pr"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Nd"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Sm"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Gd"].fix(1e-7)
    m.fs.solex.mscontactor.organic[0, :].conc_mass_comp["Dy"].fix(1e-7)

    m.fs.solex.mscontactor.organic[0, :].flow_vol.fix(62.01)

    return m


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings(ignore_unit_consistency=True)


@pytest.mark.component
@pytest.mark.solver
def test_solve(model):
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=False)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    time_duration = 24
    number_of_stages = 3
    aqueous_outlet = {
        "H2O": 1000000,
        "H": 1.75564,
        "SO4": 3999.8761,
        "HSO4": 693.4001,
        "Al": 363.321,
        "Ca": 81.9851,
        "Cl": 9.9999e-8,
        "Ce": 0.57801,
        "Dy": 3.4488e-03,
        "Fe": 461.7292,
        "Gd": 0.07428,
        "La": 0.40433,
        "Nd": 9.7359e-04,
        "Pr": 0.09789,
        "Sc": 1.3159e-03,
        "Sm": 9.9817e-05,
        "Y": 1.27602e-4,
    }

    organic_outlet = {
        "Al": 60.42,
        "Ca": 27.923,
        "Ce": 1.7771,
        "Dy": 0.045844,
        "Fe": 235.631,
        "Gd": 0.19347,
        "La": 0.6021,
        "Nd": 0.9821,
        "Pr": 0.21488,
        "Sc": 1.76797,
        "Sm": 0.10081,
        "Y": 0.12888,
    }

    for k, v in model.fs.solex.mscontactor.organic[
        time_duration, 1
    ].conc_mass_comp.items():
        assert value(v) == pytest.approx(organic_outlet[k], rel=1e-4)

    for k, v in model.fs.solex.mscontactor.aqueous[
        time_duration, number_of_stages
    ].conc_mass_comp.items():
        assert value(v) == pytest.approx(aqueous_outlet[k], rel=1e-4)
