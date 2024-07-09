#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
    value,
    Var,
)

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import MSContactor
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)


import pytest

from prommis.leaching.leach_train import LeachingTrain
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactions()

    m.fs.leach = LeachingTrain(
        number_of_tanks=1,
        liquid_phase={
            "property_package": m.fs.leach_soln,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        solid_phase={
            "property_package": m.fs.coal,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        reaction_package=m.fs.leach_rxns,
    )

    # Liquid feed state
    m.fs.leach.liquid_inlet.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach.liquid_inlet.conc_mass_comp.fix(1e-10 * units.mg / units.L)

    m.fs.leach.liquid_inlet.conc_mass_comp[0, "H"].fix(
        2 * 0.05 * 1e3 * units.mg / units.L
    )
    m.fs.leach.liquid_inlet.conc_mass_comp[0, "HSO4"].fix(1e-8 * units.mg / units.L)
    m.fs.leach.liquid_inlet.conc_mass_comp[0, "SO4"].fix(
        0.05 * 96e3 * units.mg / units.L
    )

    # Solid feed state
    m.fs.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sc2O3"].fix(
        2.77966e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Y2O3"].fix(
        3.28653e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "La2O3"].fix(
        6.77769e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Ce2O3"].fix(
        0.000156161 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Pr2O3"].fix(
        1.71438e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Nd2O3"].fix(
        6.76618e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sm2O3"].fix(
        1.47926e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Gd2O3"].fix(
        1.0405e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Dy2O3"].fix(
        7.54827e-06 * units.kg / units.kg
    )

    m.fs.leach.volume.fix(100 * units.gallon)

    return m


@pytest.mark.unit
def test_build(model):
    assert hasattr(model.fs, "leach")

    assert isinstance(model.fs.leach.mscontactor, MSContactor)
    assert isinstance(model.fs.leach.volume, Var)
    assert len(model.fs.leach.volume) == 1
    assert isinstance(
        model.fs.leach.mscontactor.heterogeneous_reaction_extent_constraint, Constraint
    )
    assert (
        len(model.fs.leach.mscontactor.heterogeneous_reaction_extent_constraint) == 12
    )

    assert number_variables(model.fs.leach) == 197
    assert number_total_constraints(model.fs.leach) == 164
    assert number_unused_variables(model.fs.leach) == 0


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve(model):
    # Cannot separate initialization, as it fails to converge in first pass
    model.scaling_factor = Suffix(direction=Suffix.EXPORT)

    for j in model.fs.coal.component_list:
        if j not in ["Al2O3", "Fe2O3", "CaO", "inerts"]:
            model.scaling_factor[
                model.fs.leach.mscontactor.solid[0.0, 1].mass_frac_comp[j]
            ] = 1e5
            model.scaling_factor[
                model.fs.leach.mscontactor.solid_inlet_state[0.0].mass_frac_comp[j]
            ] = 1e5
            model.scaling_factor[
                model.fs.leach.mscontactor.heterogeneous_reactions[
                    0.0, 1
                ].reaction_rate[j]
            ] = 1e5
            model.scaling_factor[
                model.fs.leach.mscontactor.solid[0.0, 1].conversion_eq[j]
            ] = 1e3
            model.scaling_factor[
                model.fs.leach.mscontactor.solid_inlet_state[0.0].conversion_eq[j]
            ] = 1e3
            model.scaling_factor[
                model.fs.leach.mscontactor.heterogeneous_reactions[
                    0.0, 1
                ].reaction_rate_eq[j]
            ] = 1e5

    # Create a scaled version of the model to solve
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model, rename=False)

    initializer = model.fs.leach.default_initializer()
    try:
        initializer.initialize(scaled_model.fs.leach)
    except:
        pass

    # Solve scaled model
    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=False)

    # Propagate results back to unscaled model
    scaling.propagate_solution(scaled_model, model)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    conversion = {
        "Al2O3": 0.03003765546930494,
        "CaO": 0.5145052459223426,
        "Ce2O3": 0.4207435931354908,
        "Dy2O3": 0.22365403062817607,
        "Fe2O3": 0.16338507097541013,
        "Gd2O3": 0.6106358012228825,
        "La2O3": 0.38342494088341084,
        "Nd2O3": 0.45795409973702533,
        "Pr2O3": 0.4842509098265882,
        "Sc2O3": 0.05486783341762388,
        "Sm2O3": 0.23377407833771635,
        "Y2O3": 0.15092538416884166,
        "inerts": 0.0,
    }
    recovery = {
        "inerts": 3.3306690738754696e-14,
        "Sc2O3": 5.486783341761903,
        "Y2O3": 15.092538416882373,
        "La2O3": 38.342494088335066,
        "Ce2O3": 42.074359313539645,
        "Pr2O3": 48.42509098259142,
        "Nd2O3": 45.79540997369025,
        "Sm2O3": 23.377407833768427,
        "Gd2O3": 61.06358012227954,
        "Dy2O3": 22.365403062815247,
        "Al2O3": 3.0037655469309033,
        "CaO": 51.45052459223044,
        "Fe2O3": 16.33850709753919,
    }

    for k, v in model.fs.leach.mscontactor.solid[0, 1].conversion.items():
        f_in = model.fs.leach.solid_inlet.flow_mass[0]
        f_out = model.fs.leach.solid_outlet.flow_mass[0]
        x_in = model.fs.leach.solid_inlet.mass_frac_comp[0, k]
        x_out = model.fs.leach.solid_outlet.mass_frac_comp[0, k]

        r = value(1 - f_out * x_out / (f_in * x_in)) * 100

        assert value(v) == pytest.approx(conversion[k], rel=1e-5, abs=1e-6)
        assert r == pytest.approx(recovery[k], rel=1e-5, abs=1e-6)
