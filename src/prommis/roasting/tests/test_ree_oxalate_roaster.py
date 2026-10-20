#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.core.expr import identify_variables
from pyomo.environ import (
    ComponentMap,
    ConcreteModel,
    Constraint,
    SolverFactory,
    Var,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.scaling.util import jacobian_cond
from idaes.core.util.testing import assert_solution_equivalent
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

import pytest

from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.properties import HClStrippingParameterBlock, REEOxalateParameterBlock
from prommis.properties.mixed_acid_properties import (
    MixedAcidParameterBlock,
    get_aliases,
)
from prommis.roasting.ree_oxalate_roaster import (
    REEOxalateRoaster,
    REEOxalateRoasterScaler,
)

""""
Reference:

[1] Udara S. P. R. Arachchige, Muhammad Mohsin, Morten C. Melaaen,
Optimization of post combustion carbon capture process solvent selection,
International Journal of Energy and Environment. Volume 3, Issue 6, 861-870. 
https://www.academia.edu/38249324/Optimization_of_post_combustion_carbon_capture_process_solvent_selection

"""


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )

    m.fs.prop_solid = REEOxalateParameterBlock()
    m.fs.prop_liquid = HClStrippingParameterBlock()
    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate_solid=m.fs.prop_solid,
        property_package_precipitate_liquid=m.fs.prop_liquid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        metal_list=[
            "Al",
            "Fe",
            "Ca",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Gd",
            "Dy",
        ],
    )

    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # inlet flue gas mole flow rate
    fgas = 0.00781
    # inlet flue gas composition, typical flue gas by burning CH4 with air with stoichiometric ratio of 2.3
    gas_comp = {
        "O2": 0.1118,
        "H2O": 0.1005,
        "CO2": 0.0431,
        "N2": 0.7446,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

    # fix outlet product temperature
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature, needs to be higher than 298.15 which is the lower bound
    # defined by the solid precipitate property package
    m.fs.roaster.solid_in[0].temperature.fix(299.15)
    m.fs.roaster.solid_in[0].flow_mol_comp.fix(6.1e-5)
    m.fs.roaster.liquid_in[0].flow_vol.fix(6.75e-4 * 0.018 * 3600)  # in L/hr
    m.fs.roaster.liquid_in[0].conc_mass_comp.fix(1e-5)
    m.fs.roaster.liquid_in[0].conc_mass_comp["H2O"].fix(1e6)  # mg/L
    m.fs.roaster.frac_comp_recovery.fix(0.95)

    generic_prop_scaler = m.fs.roaster.gas_in.default_scaler()
    generic_prop_scaler.default_scaling_factors["flow_mol_phase"] = 1 / fgas

    HCl_prop_scaler = m.fs.roaster.liquid_in.default_scaler()
    HCl_prop_scaler.default_scaling_factors["flow_vol"] = 1 / (6.75e-4 * 0.018 * 3600)

    submodel_scalers = ComponentMap()
    submodel_scalers[m.fs.roaster.gas_in] = generic_prop_scaler
    submodel_scalers[m.fs.roaster.gas_out] = generic_prop_scaler

    submodel_scalers[m.fs.roaster.liquid_in] = HCl_prop_scaler

    roaster_scaler = m.fs.roaster.default_scaler()

    roaster_scaler.scale_model(m.fs.roaster, submodel_scalers=submodel_scalers)

    return m


@pytest.mark.unit
def test_build(model):
    assert hasattr(model.fs, "roaster")
    assert isinstance(model.fs.roaster, REEOxalateRoaster)
    assert len(model.fs.roaster.config) == 11
    assert not model.fs.roaster.config.dynamic
    assert not model.fs.roaster.config.has_holdup
    assert model.fs.roaster.config.has_heat_transfer
    assert model.fs.roaster.config.has_pressure_change
    assert model.fs.roaster.config.property_package_gas is model.fs.prop_gas
    assert (
        model.fs.roaster.config.property_package_precipitate_solid
        is model.fs.prop_solid
    )
    assert (
        model.fs.roaster.config.property_package_precipitate_liquid
        is model.fs.prop_liquid
    )
    assert len(model.fs.prop_gas.component_list) == 4
    assert len(model.fs.roaster.metal_list) == 12
    assert isinstance(model.fs.roaster.heat_duty, Var)
    assert isinstance(model.fs.roaster.deltaP, Var)
    assert isinstance(model.fs.roaster.flow_mol_outlet_eqn, Constraint)
    assert len(model.fs.roaster.flow_mol_outlet_eqn) == 4
    assert number_variables(model.fs.roaster) == 161
    assert number_total_constraints(model.fs.roaster) == 111
    assert number_unused_variables(model.fs.roaster) == 1
    assert_units_consistent(model.fs.roaster)

    assert model.fs.roaster.default_scaler == REEOxalateRoasterScaler


@pytest.mark.unit
def test_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve(model):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(model.fs.roaster)
    assert initializer.summary[model.fs.roaster]["status"] == InitializationStatus.Ok
    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model, tee=False)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()

    assert jacobian_cond(model, scaled=False) == pytest.approx(1.6981e14, rel=1e-3)
    assert jacobian_cond(model, scaled=True) == pytest.approx(2103.0, rel=1e-3)


@pytest.mark.component
@pytest.mark.solver
def test_solution(model):
    model = model
    tol = 1e-4

    expected_results = {}

    expected_results["roaster.gas_outlet.flow_mol"] = {0: (0.008487, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "H2O"): (0.172195, tol, None),
        (0, "O2"): (0.102842, tol, None),
        (0, "CO2"): (0.039795, tol, None),
    }
    expected_results["roaster.heat_duty"] = {0: (-82.248, tol, None)}
    expected_results["roaster.flow_mass_product"] = {0: (4.9676e-8, tol, None)}
    expected_results["roaster.mass_frac_comp_product"] = {
        (0, "Al"): (0.03303983600788328, tol, None),
        (0, "Fe"): (0.051745591869350606, tol, None),
        (0, "Ca"): (0.018171407536352827, tol, None),
        (0, "Sc"): (0.04468856468660738, tol, None),
        (0, "Y"): (0.0731720199792481, tol, None),
        (0, "La"): (0.1055789983285816, tol, None),
        (0, "Ce"): (0.10636318446974416, tol, None),
        (0, "Pr"): (0.10687517376851974, tol, None),
        (0, "Nd"): (0.10903330587601676, tol, None),
        (0, "Sm"): (0.11299960272222746, tol, None),
        (0, "Gd"): (0.11746492711281438, tol, None),
        (0, "Dy"): (0.12086738764265201, tol, None),
    }

    assert_solution_equivalent(model.fs, expected_results)


@pytest.fixture(scope="module")
def model_coal_fired():
    """
    Uses the flue gas composition from a typical coal-fired plant
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_gas.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.prop_gas.set_default_scaling("pressure", 1e-5)
    m.fs.prop_gas.set_default_scaling("temperature", 1e-2)
    m.fs.prop_gas.set_default_scaling("flow_mol", 1e1)
    m.fs.prop_gas.set_default_scaling("flow_mol_phase", 1e1)
    m.fs.prop_gas.set_default_scaling("_energy_density_term", 1e-4)
    m.fs.prop_gas.set_default_scaling("phase_frac", 1)

    _mf_scale = {
        "O2": 5,
        "CO2": 10,
        "H2O": 5,
        "N2": 1,
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )

    m.fs.prop_solid = REEOxalateParameterBlock()
    m.fs.prop_liquid = HClStrippingParameterBlock()
    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate_solid=m.fs.prop_solid,
        property_package_precipitate_liquid=m.fs.prop_liquid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        metal_list=[
            "Al",
            "Fe",
            "Ca",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Gd",
            "Dy",
        ],
    )

    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # inlet flue gas mole flow rate
    fgas = 0.00781
    # inlet flue gas composition from coal-fired plant
    gas_comp = {
        "O2": 0.05,
        "H2O": 0.08,
        "CO2": 0.14,
        "N2": 0.73,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

    # fix outlet product temperature
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature, needs to be higher than 298.15 which is the lower bound
    # defined by the solid precipitate property package
    m.fs.roaster.solid_in[0].temperature.fix(299.15)
    m.fs.roaster.solid_in[0].flow_mol_comp.fix(6.1e-5)
    m.fs.roaster.liquid_in[0].flow_vol.fix(6.75e-4 * 0.018 * 3600)  # in L/hr
    m.fs.roaster.liquid_in[0].conc_mass_comp.fix(1e-5)
    m.fs.roaster.liquid_in[0].conc_mass_comp["H2O"].fix(1e6)  # mg/L
    m.fs.roaster.frac_comp_recovery.fix(0.95)

    return m


@pytest.mark.unit
def test_build_coal_fired(model_coal_fired):
    assert hasattr(model_coal_fired.fs, "roaster")
    assert isinstance(model_coal_fired.fs.roaster, REEOxalateRoaster)
    assert len(model_coal_fired.fs.roaster.config) == 11
    assert not model_coal_fired.fs.roaster.config.dynamic
    assert not model_coal_fired.fs.roaster.config.has_holdup
    assert model_coal_fired.fs.roaster.config.has_heat_transfer
    assert model_coal_fired.fs.roaster.config.has_pressure_change
    assert (
        model_coal_fired.fs.roaster.config.property_package_gas
        is model_coal_fired.fs.prop_gas
    )
    assert (
        model_coal_fired.fs.roaster.config.property_package_precipitate_solid
        is model_coal_fired.fs.prop_solid
    )
    assert (
        model_coal_fired.fs.roaster.config.property_package_precipitate_liquid
        is model_coal_fired.fs.prop_liquid
    )
    assert len(model_coal_fired.fs.prop_gas.component_list) == 4
    assert len(model_coal_fired.fs.roaster.metal_list) == 12
    assert isinstance(model_coal_fired.fs.roaster.heat_duty, Var)
    assert isinstance(model_coal_fired.fs.roaster.deltaP, Var)
    assert isinstance(model_coal_fired.fs.roaster.flow_mol_outlet_eqn, Constraint)
    assert len(model_coal_fired.fs.roaster.flow_mol_outlet_eqn) == 4
    assert number_variables(model_coal_fired.fs.roaster) == 161
    assert number_total_constraints(model_coal_fired.fs.roaster) == 111
    assert number_unused_variables(model_coal_fired.fs.roaster) == 1
    assert_units_consistent(model_coal_fired.fs.roaster)


@pytest.mark.unit
def test_structural_issues_coal_fired(model_coal_fired):
    dt = DiagnosticsToolbox(model_coal_fired)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve_coal_fired(model_coal_fired):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(model_coal_fired.fs.roaster)
    assert (
        initializer.summary[model_coal_fired.fs.roaster]["status"]
        == InitializationStatus.Ok
    )
    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model_coal_fired, tee=False)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues_coal_fired(model_coal_fired):
    dt = DiagnosticsToolbox(model_coal_fired)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution_coal_fired(model_coal_fired):
    model = model_coal_fired
    tol = 1e-4

    expected_results = {}

    expected_results["roaster.gas_outlet.flow_mol"] = {0: (0.008487, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "H2O"): (0.153331, tol, None),
        (0, "O2"): (0.0459751, tol, None),
        (0, "CO2"): (0.128961, tol, None),
    }
    expected_results["roaster.heat_duty"] = {0: (-88.75058, tol, None)}
    expected_results["roaster.flow_mass_product"] = {0: (4.9676e-8, tol, None)}
    expected_results["roaster.mass_frac_comp_product"] = {
        (0, "Al"): (0.03303983600788328, tol, None),
        (0, "Fe"): (0.051745591869350606, tol, None),
        (0, "Ca"): (0.018171407536352827, tol, None),
        (0, "Sc"): (0.04468856468660738, tol, None),
        (0, "Y"): (0.0731720199792481, tol, None),
        (0, "La"): (0.1055789983285816, tol, None),
        (0, "Ce"): (0.10636318446974416, tol, None),
        (0, "Pr"): (0.10687517376851974, tol, None),
        (0, "Nd"): (0.10903330587601676, tol, None),
        (0, "Sm"): (0.11299960272222746, tol, None),
        (0, "Gd"): (0.11746492711281438, tol, None),
        (0, "Dy"): (0.12086738764265201, tol, None),
    }

    assert_solution_equivalent(model.fs, expected_results)


@pytest.fixture(scope="module")
def model_gas_fired():
    """
    Uses the flue gas composition from a typical gas-fired plant
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    gas_species = {"O2", "H2O", "CO2", "N2"}
    m.fs.prop_gas = GenericParameterBlock(
        **get_prop(gas_species, ["Vap"], EosType.IDEAL),
        doc="gas property",
    )
    m.fs.prop_gas.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.prop_gas.set_default_scaling("pressure", 1e-5)
    m.fs.prop_gas.set_default_scaling("temperature", 1e-2)
    m.fs.prop_gas.set_default_scaling("flow_mol", 1e1)
    m.fs.prop_gas.set_default_scaling("flow_mol_phase", 1e1)
    m.fs.prop_gas.set_default_scaling("_energy_density_term", 1e-4)
    m.fs.prop_gas.set_default_scaling("phase_frac", 1)

    _mf_scale = {
        "O2": 5,
        "CO2": 10,
        "H2O": 5,
        "N2": 1,
    }
    for comp, s in _mf_scale.items():
        m.fs.prop_gas.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.prop_gas.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.prop_gas.set_default_scaling(
            "flow_mol_phase_comp", s * 1e1, index=("Vap", comp)
        )

    m.fs.prop_solid = REEOxalateParameterBlock()
    m.fs.prop_liquid = HClStrippingParameterBlock()
    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate_solid=m.fs.prop_solid,
        property_package_precipitate_liquid=m.fs.prop_liquid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        metal_list=[
            "Al",
            "Fe",
            "Ca",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Gd",
            "Dy",
        ],
    )

    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # inlet flue gas mole flow rate
    fgas = 0.00781
    # inlet flue gas composition from gas-fired plant
    gas_comp = {
        "O2": 0.12,
        "H2O": 0.08,
        "CO2": 0.04,
        "N2": 0.76,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

    # fix outlet product temperature
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature, needs to be higher than 298.15 which is the lower bound
    # defined by the solid precipitate property package
    m.fs.roaster.solid_in[0].temperature.fix(299.15)
    m.fs.roaster.solid_in[0].flow_mol_comp.fix(6.1e-5)
    m.fs.roaster.liquid_in[0].flow_vol.fix(6.75e-4 * 0.018 * 3600)  # in L/hr
    m.fs.roaster.liquid_in[0].conc_mass_comp.fix(1e-5)
    m.fs.roaster.liquid_in[0].conc_mass_comp["H2O"].fix(1e6)  # mg/L
    m.fs.roaster.frac_comp_recovery.fix(0.95)

    return m


@pytest.mark.unit
def test_build_gas_fired(model_gas_fired):
    assert hasattr(model_gas_fired.fs, "roaster")
    assert isinstance(model_gas_fired.fs.roaster, REEOxalateRoaster)
    assert len(model_gas_fired.fs.roaster.config) == 11
    assert not model_gas_fired.fs.roaster.config.dynamic
    assert not model_gas_fired.fs.roaster.config.has_holdup
    assert model_gas_fired.fs.roaster.config.has_heat_transfer
    assert model_gas_fired.fs.roaster.config.has_pressure_change
    assert (
        model_gas_fired.fs.roaster.config.property_package_gas
        is model_gas_fired.fs.prop_gas
    )
    assert (
        model_gas_fired.fs.roaster.config.property_package_precipitate_solid
        is model_gas_fired.fs.prop_solid
    )
    assert (
        model_gas_fired.fs.roaster.config.property_package_precipitate_liquid
        is model_gas_fired.fs.prop_liquid
    )
    assert len(model_gas_fired.fs.prop_gas.component_list) == 4
    assert len(model_gas_fired.fs.roaster.metal_list) == 12
    assert isinstance(model_gas_fired.fs.roaster.heat_duty, Var)
    assert isinstance(model_gas_fired.fs.roaster.deltaP, Var)
    assert isinstance(model_gas_fired.fs.roaster.flow_mol_outlet_eqn, Constraint)
    assert len(model_gas_fired.fs.roaster.flow_mol_outlet_eqn) == 4
    assert number_variables(model_gas_fired.fs.roaster) == 161
    assert number_total_constraints(model_gas_fired.fs.roaster) == 111
    assert number_unused_variables(model_gas_fired.fs.roaster) == 1
    assert_units_consistent(model_gas_fired.fs.roaster)


@pytest.mark.unit
def test_structural_issues_gas_fired(model_gas_fired):
    dt = DiagnosticsToolbox(model_gas_fired)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_initialize_and_solve_gas_fired(model_gas_fired):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(model_gas_fired.fs.roaster)
    assert (
        initializer.summary[model_gas_fired.fs.roaster]["status"]
        == InitializationStatus.Ok
    )
    # Solve model
    solver = SolverFactory("ipopt")
    results = solver.solve(model_gas_fired, tee=False)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues_gas_fired(model_gas_fired):
    dt = DiagnosticsToolbox(model_gas_fired)
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution_gas_fired(model_gas_fired):
    model = model_gas_fired
    tol = 1e-4

    expected_results = {}

    expected_results["roaster.gas_outlet.flow_mol"] = {0: (0.008487, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "H2O"): (0.153331, tol, None),
        (0, "O2"): (0.110388, tol, None),
        (0, "CO2"): (0.036943, tol, None),
    }
    expected_results["roaster.heat_duty"] = {0: (-81.3816, tol, None)}
    expected_results["roaster.flow_mass_product"] = {0: (4.9676e-8, tol, None)}
    expected_results["roaster.mass_frac_comp_product"] = {
        (0, "Al"): (0.03303983600788328, tol, None),
        (0, "Fe"): (0.051745591869350606, tol, None),
        (0, "Ca"): (0.018171407536352827, tol, None),
        (0, "Sc"): (0.04468856468660738, tol, None),
        (0, "Y"): (0.0731720199792481, tol, None),
        (0, "La"): (0.1055789983285816, tol, None),
        (0, "Ce"): (0.10636318446974416, tol, None),
        (0, "Pr"): (0.10687517376851974, tol, None),
        (0, "Nd"): (0.10903330587601676, tol, None),
        (0, "Sm"): (0.11299960272222746, tol, None),
        (0, "Gd"): (0.11746492711281438, tol, None),
        (0, "Dy"): (0.12086738764265201, tol, None),
    }

    assert_solution_equivalent(model.fs, expected_results)


@pytest.mark.parametrize(
    "use_mixed_acid",
    [False, True],
    scope="class",
)
class Test_REE_Oxalate_Roaster_model:

    @pytest.fixture(scope="class")
    @classmethod
    def model_mixed_acid(cls, use_mixed_acid):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        gas_species = {"O2", "H2O", "CO2", "N2"}
        m.fs.prop_gas = GenericParameterBlock(
            **get_prop(gas_species, ["Vap"], EosType.IDEAL),
            doc="gas property",
        )

        if use_mixed_acid:
            m.fs.prop_solid = PrecipitateParameters(
                aqueous_aliases=get_aliases(include_sulfates=True)
            )
            m.fs.prop_liquid = MixedAcidParameterBlock(include_sulfates=True)
            # The precipitate liquid properties use a molecular weight of
            # water with fewer significant figures than the mixed acid
            # properties, which then leads to a discrepancy in the model
            # solution. Use the version with fewer significant figures here
            # to make sure that the models are structurally equivalent.
            m.fs.prop_liquid.mw["H2O"] = 18e-3
        else:
            m.fs.prop_solid = PrecipitateParameters()
            m.fs.prop_liquid = AqueousParameter()

        m.fs.roaster = REEOxalateRoaster(
            property_package_gas=m.fs.prop_gas,
            property_package_precipitate_solid=m.fs.prop_solid,
            property_package_precipitate_liquid=m.fs.prop_liquid,
            has_holdup=False,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.roaster.deltaP.fix(0)
        m.fs.roaster.gas_inlet.temperature.fix(1330)
        m.fs.roaster.gas_inlet.pressure.fix(101325)
        # inlet flue gas mole flow rate
        fgas = 0.00781
        # inlet flue gas composition, typical flue gas by burning CH4 with air with stoichiometric ratio of 2.3
        gas_comp = {
            "O2": 0.1118,
            "H2O": 0.1005,
            "CO2": 0.0431,
            "N2": 0.7446,
        }
        for i, v in gas_comp.items():
            m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
        m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

        # fix outlet product temperature
        m.fs.roaster.gas_outlet.temperature.fix(873.15)

        m.fs.roaster.solid_in[0].temperature.fix(299.15)
        m.fs.roaster.solid_in[0].flow_mol_comp.fix(6.1e-5)
        m.fs.roaster.liquid_in[0].flow_vol.fix(6.75e-4 * 0.018 * 3600)  # in L/hr
        m.fs.roaster.liquid_in[0].conc_mass_comp.fix(1e-5)
        m.fs.roaster.liquid_in[0].conc_mass_comp["H2O"].fix(1e6)  # mg/L
        m.fs.roaster.frac_comp_recovery.fix(0.95)

        return m

    @pytest.mark.unit
    def test_build(self, model_mixed_acid, use_mixed_acid):
        m = model_mixed_acid

        assert isinstance(m.fs.roaster, REEOxalateRoaster)
        assert len(m.fs.roaster.config) == 11
        assert len(m.fs.roaster.metal_list) == 12
        # The mixed acid property package has more variables that are not used
        # by the roaster
        if use_mixed_acid:
            assert number_variables(m.fs.roaster) == 169
            assert number_unused_variables(m.fs.roaster) == 3
        else:
            assert number_variables(m.fs.roaster) == 167
            assert number_unused_variables(m.fs.roaster) == 1
        assert number_total_constraints(m.fs.roaster) == 115
        assert_units_consistent(m.fs.roaster)

        # The oxalate feed of each metal must come from the matching solid
        solid_props = m.fs.roaster.config.property_package_precipitate_solid
        reversed_react = dict(map(reversed, solid_props.react.items()))
        for i in m.fs.roaster.metal_list:
            solid = reversed_react[i]

            feed_vars = {
                v.name
                for v in identify_variables(
                    m.fs.roaster.flow_mol_comp_feed_eqn[0, i].body
                )
            }

            assert m.fs.roaster.solid_in[0].flow_mol_comp[solid].name in feed_vars
            assert m.fs.roaster.flow_mol_comp_feed[0, i].name in feed_vars

    @pytest.mark.unit
    def test_structural_issues(self, model_mixed_acid):
        dt = DiagnosticsToolbox(model_mixed_acid)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    @pytest.mark.solver
    def test_initialize_and_solve(self, model_mixed_acid):
        m = model_mixed_acid

        initializer = BlockTriangularizationInitializer()
        initializer.initialize(m.fs.roaster)
        assert initializer.summary[m.fs.roaster]["status"] == InitializationStatus.Ok
        # Solve model
        solver = SolverFactory("ipopt")
        results = solver.solve(m, tee=False)
        assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.solver
    def test_solution(self, model_mixed_acid):
        # The names of the species do not change the results, so these are
        # the same as the results using the other property packages.
        m = model_mixed_acid

        flow_mol_out_gas = value(m.fs.roaster.gas_out[0].flow_mol)
        assert flow_mol_out_gas == pytest.approx(0.008487, rel=1e-5, abs=1e-6)
        mole_frac_h2o = value(m.fs.roaster.gas_out[0].mole_frac_comp["H2O"])
        assert mole_frac_h2o == pytest.approx(0.172195, rel=1e-5, abs=1e-6)
        mole_frac_o2 = value(m.fs.roaster.gas_out[0].mole_frac_comp["O2"])
        assert mole_frac_o2 == pytest.approx(0.102842, rel=1e-5, abs=1e-6)
        mole_frac_co2 = value(m.fs.roaster.gas_out[0].mole_frac_comp["CO2"])
        assert mole_frac_co2 == pytest.approx(0.039795, rel=1e-5, abs=1e-6)
        heat_duty = value(m.fs.roaster.heat_duty[0])
        assert heat_duty == pytest.approx(-82.248, rel=1e-5, abs=1e-6)
        flow_mass_product = value(m.fs.roaster.flow_mass_product[0])
        assert flow_mass_product == pytest.approx(4.9676e-08, rel=1e-5, abs=1e-9)
        mass_frac_comp_product = {
            "Al": 0.03303983600788328,
            "Fe": 0.051745591869350606,
            "Ca": 0.018171407536352827,
            "Sc": 0.04468856468660738,
            "Y": 0.0731720199792481,
            "La": 0.1055789983285816,
            "Ce": 0.10636318446974416,
            "Pr": 0.10687517376851974,
            "Nd": 0.10903330587601676,
            "Sm": 0.11299960272222746,
            "Gd": 0.11746492711281438,
            "Dy": 0.12086738764265201,
        }

        for i in m.fs.roaster.metal_list:
            base_metal = m.fs.roaster.base_metal[i]
            assert value(m.fs.roaster.mass_frac_comp_product[0, i]) == pytest.approx(
                mass_frac_comp_product[base_metal], rel=1e-5, abs=1e-6
            )
