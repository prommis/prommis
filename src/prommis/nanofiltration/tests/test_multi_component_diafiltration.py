#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the multi-component diafiltration unit model.

Author: Molly Dougher
"""

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    Set,
    SolverFactory,
    TransformationFactory,
    Var,
    assert_optimal_termination,
    value,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)
from prommis.nanofiltration.multi_component_diafiltration import (
    MultiComponentDiafiltration,
)
from prommis.util import assert_solution_equivalent

# TODO: test positive and neutral membrane cases
# currently, the property package only supports a negative fixed charge


################################################################################
# Test functions for single-salt model
@pytest.fixture(scope="module")
def sample_single_salt_model():
    cation_list = ["Li"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Cl": 245},
        "diafiltrate": {"Li": 14, "Cl": 14},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.fixture(scope="module")
def sample_single_salt_model_no_boundary_layer():
    cation_list = ["Li"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Cl": 245},
        "diafiltrate": {"Li": 14, "Cl": 14},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.mark.unit
def test_config(sample_single_salt_model):
    membrane = sample_single_salt_model.fs.unit

    assert len(membrane.config) == 10

    assert not membrane.config.dynamic
    assert not membrane.config.has_holdup

    assert membrane.config.property_package is sample_single_salt_model.fs.properties

    assert len(membrane.config.anion_list) == 1
    assert membrane.config.NFE_module_length == 10
    assert membrane.config.NFE_boundary_layer_thickness == 5
    assert membrane.config.NFE_membrane_thickness == 5


@pytest.mark.unit
def test_config_boundary_layer(sample_single_salt_model):
    assert sample_single_salt_model.fs.unit.config.include_boundary_layer == True


@pytest.mark.unit
def test_config_no_boundary_layer(sample_single_salt_model_no_boundary_layer):
    assert (
        sample_single_salt_model_no_boundary_layer.fs.unit.config.include_boundary_layer
        == False
    )


@pytest.mark.unit
def test_config_single_salt(sample_single_salt_model):
    assert len(sample_single_salt_model.fs.unit.config.cation_list) == 1


@pytest.mark.build
@pytest.mark.unit
def test_build(sample_single_salt_model):
    membrane = sample_single_salt_model.fs.unit
    # parameters
    assert isinstance(membrane.numerical_zero_tolerance, Param)
    assert value(membrane.numerical_zero_tolerance) == 1e-10

    assert isinstance(membrane.total_membrane_thickness, Param)
    assert value(membrane.total_membrane_thickness) == 1e-7

    assert isinstance(membrane.membrane_fixed_charge, Param)
    assert value(membrane.membrane_fixed_charge) == -44

    assert isinstance(membrane.membrane_permeability, Param)
    assert value(membrane.membrane_permeability) == 0.01

    assert isinstance(membrane.temperature, Param)
    assert value(membrane.temperature) == 298

    # sets
    assert isinstance(membrane.dimensionless_module_length, ContinuousSet)
    assert len(membrane.dimensionless_module_length) == 11

    assert isinstance(
        membrane.dimensionless_membrane_thickness,
        ContinuousSet,
    )
    assert len(membrane.dimensionless_membrane_thickness) == 6

    assert isinstance(membrane.time, Set)
    assert len(membrane.time) == 1

    # variables
    assert isinstance(membrane.total_module_length, Var)
    assert len(membrane.total_module_length) == 1

    assert isinstance(membrane.total_membrane_length, Var)
    assert len(membrane.total_membrane_length) == 1

    assert isinstance(membrane.applied_pressure, Var)
    assert len(membrane.applied_pressure) == 1

    assert isinstance(membrane.feed_flow_volume, Var)
    assert len(membrane.feed_flow_volume) == 1

    assert isinstance(membrane.diafiltrate_flow_volume, Var)
    assert len(membrane.diafiltrate_flow_volume) == 1

    assert isinstance(membrane.membrane_D_tilde, Var)
    assert len(membrane.membrane_D_tilde) == 66

    assert isinstance(membrane.volume_flux_water, Var)
    assert len(membrane.volume_flux_water) == 11

    assert isinstance(membrane.retentate_flow_volume, Var)
    assert len(membrane.retentate_flow_volume) == 11

    assert isinstance(membrane.permeate_flow_volume, Var)
    assert len(membrane.permeate_flow_volume) == 11

    assert isinstance(membrane.osmotic_pressure, Var)
    assert len(membrane.osmotic_pressure) == 11

    assert isinstance(membrane.d_retentate_flow_volume_dx, DerivativeVar)
    assert len(membrane.d_retentate_flow_volume_dx) == 11

    # constraints
    assert isinstance(membrane.overall_mol_balance, Constraint)
    assert len(membrane.overall_mol_balance) == 10

    assert isinstance(membrane.overall_bulk_flux_equation, Constraint)
    assert len(membrane.overall_bulk_flux_equation) == 10

    assert isinstance(membrane.lumped_water_flux, Constraint)
    assert len(membrane.lumped_water_flux) == 10

    assert isinstance(membrane.anion_flux_membrane, Constraint)
    assert len(membrane.anion_flux_membrane) == 10

    assert isinstance(membrane.osmotic_pressure_calculation, Constraint)
    assert len(membrane.osmotic_pressure_calculation) == 10

    assert isinstance(membrane.electroneutrality_retentate, Constraint)
    assert len(membrane.electroneutrality_retentate) == 11

    assert isinstance(membrane.electroneutrality_permeate, Constraint)
    assert len(membrane.electroneutrality_permeate) == 10

    assert isinstance(membrane.electroneutrality_membrane, Constraint)
    assert len(membrane.electroneutrality_membrane) == 60

    assert isinstance(membrane.membrane_D_tilde_calculation, Constraint)
    assert len(membrane.membrane_D_tilde_calculation) == 60

    assert isinstance(
        membrane.retentate_flow_volume_boundary_condition,
        Constraint,
    )
    assert len(membrane.retentate_flow_volume_boundary_condition) == 1

    assert isinstance(
        membrane.permeate_flow_volume_boundary_condition,
        Constraint,
    )
    assert len(membrane.permeate_flow_volume_boundary_condition) == 1

    assert isinstance(
        membrane.d_retentate_flow_volume_dx_boundary_condition,
        Constraint,
    )
    assert len(membrane.d_retentate_flow_volume_dx_boundary_condition) == 1

    assert isinstance(
        membrane.volume_flux_water_boundary_condition,
        Constraint,
    )
    assert len(membrane.volume_flux_water_boundary_condition) == 1

    for t in membrane.time:
        for x in membrane.dimensionless_module_length:
            assert membrane.d_retentate_conc_mol_comp_dx[
                t, x, membrane.config.anion_list[0]
            ].fixed
            if x != 0:
                assert not membrane.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, membrane.config.anion_list[0]
                ].active

            for z in membrane.dimensionless_membrane_thickness:
                assert membrane.d_membrane_conc_mol_comp_dz[
                    t, x, z, membrane.config.anion_list[0]
                ].fixed
                if z != 0:
                    assert not membrane.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, membrane.config.anion_list[0]
                    ].active

    # scaling factors
    assert membrane.scaling_factor[membrane.volume_flux_water] == 1e2
    assert membrane.scaling_factor[membrane.membrane_D_tilde] == 1e-1
    assert (
        membrane.scaling_factor[membrane.membrane_cross_diffusion_coefficient_bilinear]
        == 1e-2
    )
    assert (
        membrane.scaling_factor[membrane.membrane_convection_coefficient_bilinear]
        == 1e-1
    )
    assert membrane.scaling_factor[membrane.membrane_cross_diffusion_coefficient] == 1e1
    assert membrane.scaling_factor[membrane.membrane_convection_coefficient] == 1e1

    # ports
    assert isinstance(membrane.feed_inlet, Port)
    assert len(membrane.feed_inlet.flow_vol) == 1

    assert isinstance(membrane.diafiltrate_inlet, Port)
    assert len(membrane.diafiltrate_inlet.flow_vol) == 1

    assert isinstance(membrane.retentate_outlet, Port)
    assert len(membrane.retentate_outlet.flow_vol) == 1

    assert isinstance(membrane.permeate_outlet, Port)
    assert len(membrane.permeate_outlet.flow_vol) == 1


@pytest.mark.build
@pytest.mark.unit
def test_build_boundary_layer(sample_single_salt_model):
    membrane = sample_single_salt_model.fs.unit
    # parameters
    assert isinstance(membrane.total_boundary_layer_thickness, Param)
    assert value(membrane.total_boundary_layer_thickness) == 2e-5

    # sets
    assert isinstance(
        membrane.dimensionless_boundary_layer_thickness,
        ContinuousSet,
    )
    assert len(membrane.dimensionless_boundary_layer_thickness) == 6

    # variables
    assert isinstance(membrane.boundary_layer_D_tilde, Var)
    assert len(membrane.boundary_layer_D_tilde) == 66

    # constraints
    assert isinstance(membrane.electroneutrality_boundary_layer, Constraint)
    assert len(membrane.electroneutrality_boundary_layer) == 60

    assert isinstance(membrane.boundary_layer_D_tilde_calculation, Constraint)
    assert len(membrane.boundary_layer_D_tilde_calculation) == 60

    for t in membrane.time:
        for x in membrane.dimensionless_module_length:
            for z in membrane.dimensionless_boundary_layer_thickness:
                assert membrane.d_boundary_layer_conc_mol_comp_dz[
                    t, x, z, membrane.config.anion_list[0]
                ].fixed
                if z != 0:
                    assert not membrane.d_boundary_layer_conc_mol_comp_dz_disc_eq[
                        t, x, z, membrane.config.anion_list[0]
                    ].active

    # scaling factors
    assert membrane.scaling_factor[membrane.boundary_layer_D_tilde] == 1e-3
    assert (
        membrane.scaling_factor[
            membrane.boundary_layer_cross_diffusion_coefficient_bilinear
        ]
        == 1e-4
    )
    assert (
        membrane.scaling_factor[membrane.boundary_layer_cross_diffusion_coefficient]
        == 1e1
    )


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt(sample_single_salt_model):
    membrane = sample_single_salt_model.fs.unit
    # sets
    assert isinstance(membrane.solutes, Set)
    assert len(membrane.solutes) == 2

    assert isinstance(membrane.cations, Set)
    assert len(membrane.cations) == 1

    # variables
    assert isinstance(membrane.feed_conc_mol_comp, Var)
    assert len(membrane.feed_conc_mol_comp) == 2

    assert isinstance(membrane.diafiltrate_conc_mol_comp, Var)
    assert len(membrane.diafiltrate_conc_mol_comp) == 2

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 66

    assert isinstance(
        membrane.membrane_convection_coefficient_bilinear,
        Var,
    )
    assert len(membrane.membrane_convection_coefficient_bilinear) == 66

    assert isinstance(membrane.membrane_cross_diffusion_coefficient, Var)
    assert len(membrane.membrane_cross_diffusion_coefficient) == 66

    assert isinstance(membrane.membrane_convection_coefficient, Var)
    assert len(membrane.membrane_convection_coefficient) == 66

    assert isinstance(membrane.molar_ion_flux, Var)
    assert len(membrane.molar_ion_flux) == 22

    assert isinstance(membrane.retentate_conc_mol_comp, Var)
    assert len(membrane.retentate_conc_mol_comp) == 22

    assert isinstance(membrane.permeate_conc_mol_comp, Var)
    assert len(membrane.permeate_conc_mol_comp) == 22

    assert isinstance(membrane.membrane_conc_mol_comp, Var)
    assert len(membrane.membrane_conc_mol_comp) == 132

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx) == 22

    assert isinstance(
        membrane.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_membrane_conc_mol_comp_dz) == 132

    # constraints
    assert isinstance(membrane.cation_mol_balance, Constraint)
    assert len(membrane.cation_mol_balance) == 10

    assert isinstance(membrane.cation_bulk_flux_equation, Constraint)
    assert len(membrane.cation_bulk_flux_equation) == 10

    assert isinstance(
        membrane.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_membrane_permeate_interface) == 10

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation) == 60

    assert isinstance(
        membrane.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 60

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 60

    assert isinstance(
        membrane.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_calculation) == 60

    assert isinstance(membrane.cation_flux_membrane, Constraint)
    assert len(membrane.cation_flux_membrane) == 60

    assert isinstance(
        membrane.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.retentate_conc_mol_comp_boundary_condition) == 1

    assert isinstance(
        membrane.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 6

    assert isinstance(
        membrane.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 2

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx_boundary_condition) == 1

    assert isinstance(
        membrane.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(membrane.molar_ion_flux_boundary_condition) == 2

    # ports
    assert len(membrane.feed_inlet.conc_mol_comp) == 2
    assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 2
    assert len(membrane.retentate_outlet.conc_mol_comp) == 2
    assert len(membrane.permeate_outlet.conc_mol_comp) == 2


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt_boundary_layer(sample_single_salt_model):
    membrane = sample_single_salt_model.fs.unit
    # variables
    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 66

    assert isinstance(membrane.boundary_layer_cross_diffusion_coefficient, Var)
    assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 66

    assert isinstance(membrane.boundary_layer_conc_mol_comp, Var)
    assert len(membrane.boundary_layer_conc_mol_comp) == 132

    assert isinstance(
        membrane.d_boundary_layer_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 132

    # constraints
    assert isinstance(membrane.retentate_boundary_layer_interface, Constraint)
    assert len(membrane.retentate_boundary_layer_interface) == 10

    assert isinstance(
        membrane.cation_equilibrium_boundary_layer_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_boundary_layer_membrane_interface) == 10

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation)
        == 60
    )

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 60

    assert isinstance(membrane.cation_flux_boundary_layer, Constraint)
    assert len(membrane.cation_flux_boundary_layer) == 50

    assert isinstance(
        membrane.boundary_layer_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 6


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt_no_boundary_layer(
    sample_single_salt_model_no_boundary_layer,
):
    membrane = sample_single_salt_model_no_boundary_layer.fs.unit
    assert isinstance(
        membrane.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_retentate_membrane_interface) == 10


@pytest.mark.component
def test_diagnostics(sample_single_salt_model):
    dt = DiagnosticsToolbox(sample_single_salt_model.fs.unit)
    dt.assert_no_structural_warnings()


@pytest.mark.solver
@pytest.mark.component
def test_solve(sample_single_salt_model):
    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(sample_single_salt_model, rename=False)

    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=True)

    scaling.propagate_solution(scaled_model, sample_single_salt_model)

    assert_optimal_termination(results)


@pytest.mark.component
def test_numerical_issues(sample_single_salt_model):
    dt = DiagnosticsToolbox(sample_single_salt_model.fs.unit)
    dt.assert_no_numerical_warnings()


################################################################################
# Test single-salt model: LiCl


@pytest.fixture(scope="module")
def diafiltration_single_salt_Li():
    """
    Build single-salt diafiltration unit model for LiCl
    with boundary layer.
    """
    cation_list = ["Li"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Cl": 245},
        "diafiltrate": {"Li": 14, "Cl": 14},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li(diafiltration_single_salt_Li):
    test_config(diafiltration_single_salt_Li)
    test_config_boundary_layer(diafiltration_single_salt_Li)
    test_config_single_salt(diafiltration_single_salt_Li)


class TestDiafiltrationSingleSaltLithium(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li(self, diafiltration_single_salt_Li):
        test_build(diafiltration_single_salt_Li)
        test_build_boundary_layer(diafiltration_single_salt_Li)
        test_build_single_salt(diafiltration_single_salt_Li)
        test_build_single_salt_boundary_layer(diafiltration_single_salt_Li)

    @pytest.mark.component
    def test_diagnostics_Li(self, diafiltration_single_salt_Li):
        test_diagnostics(diafiltration_single_salt_Li)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li(self, diafiltration_single_salt_Li):
        test_solve(diafiltration_single_salt_Li)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (8.3858, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (194.01, 1e-4, None),
                (0, 1, "Cl"): (194.01, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (7.8636, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (190.62, 1e-4, None),
                (0, 1, "Cl"): (190.62, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_single_salt_Li.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Li(self, diafiltration_single_salt_Li):
        test_numerical_issues(diafiltration_single_salt_Li)


@pytest.fixture(scope="module")
def diafiltration_single_salt_Li_no_boundary_layer():
    """
    Build single-salt diafiltration unit model for LiCl
    without boundary layer.
    """
    cation_list = ["Li"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Cl": 245},
        "diafiltrate": {"Li": 14, "Cl": 14},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_no_boundary_layer(
    diafiltration_single_salt_Li_no_boundary_layer,
):
    test_config(diafiltration_single_salt_Li_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_single_salt_Li_no_boundary_layer)
    test_config_single_salt(diafiltration_single_salt_Li_no_boundary_layer)


class TestDiafiltrationSingleSaltLithiumNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_no_boundary_layer(
        self, diafiltration_single_salt_Li_no_boundary_layer
    ):
        test_build(diafiltration_single_salt_Li_no_boundary_layer)
        test_build_single_salt(diafiltration_single_salt_Li_no_boundary_layer)
        test_build_single_salt_no_boundary_layer(
            diafiltration_single_salt_Li_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Li_no_boundary_layer(
        self, diafiltration_single_salt_Li_no_boundary_layer
    ):
        test_diagnostics(diafiltration_single_salt_Li_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_no_boundary_layer(
        self, diafiltration_single_salt_Li_no_boundary_layer
    ):
        test_solve(diafiltration_single_salt_Li_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (8.3856, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (194.52, 1e-4, None),
                (0, 1, "Cl"): (194.52, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (7.8637, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (190.38, 1e-4, None),
                (0, 1, "Cl"): (190.38, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_single_salt_Li_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Li_no_boundary_layer(
        self, diafiltration_single_salt_Li_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_single_salt_Li_no_boundary_layer)


################################################################################
# Test single-salt model: CoCl2


@pytest.fixture(scope="module")
def diafiltration_single_salt_Co():
    """
    Build single-salt diafiltration unit model for CoCl2
    with boundary layer.
    """
    cation_list = ["Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Co": 288, "Cl": 576},
        "diafiltrate": {"Co": 3, "Cl": 6},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Co(diafiltration_single_salt_Co):
    test_config(diafiltration_single_salt_Co)
    test_config_boundary_layer(diafiltration_single_salt_Co)
    test_config_single_salt(diafiltration_single_salt_Co)


class TestDiafiltrationSingleSaltCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co(self, diafiltration_single_salt_Co):
        test_build(diafiltration_single_salt_Co)
        test_build_boundary_layer(diafiltration_single_salt_Co)
        test_build_single_salt(diafiltration_single_salt_Co)
        test_build_single_salt_boundary_layer(diafiltration_single_salt_Co)

    @pytest.mark.component
    def test_diagnostics_Co(self, diafiltration_single_salt_Co):
        test_diagnostics(diafiltration_single_salt_Co)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co(self, diafiltration_single_salt_Co):
        test_solve(diafiltration_single_salt_Co)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (10.475, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Co"): (226.82, 1e-4, None),
                (0, 1, "Cl"): (453.64, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (5.7640, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Co"): (216.55, 1e-4, None),
                (0, 1, "Cl"): (433.11, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_single_salt_Co.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Co(self, diafiltration_single_salt_Co):
        test_numerical_issues(diafiltration_single_salt_Co)


@pytest.fixture(scope="module")
def diafiltration_single_salt_Co_no_boundary_layer():
    """
    Build single-salt diafiltration unit model for CoCl2
    without boundary layer.
    """
    cation_list = ["Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Co": 288, "Cl": 576},
        "diafiltrate": {"Co": 3, "Cl": 6},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Co_no_boundary_layer(
    diafiltration_single_salt_Co_no_boundary_layer,
):
    test_config(diafiltration_single_salt_Co_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_single_salt_Co_no_boundary_layer)
    test_config_single_salt(diafiltration_single_salt_Co_no_boundary_layer)


class TestDiafiltrationSingleSaltCobaltNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co_no_boundary_layer(
        self, diafiltration_single_salt_Co_no_boundary_layer
    ):
        test_build(diafiltration_single_salt_Co_no_boundary_layer)
        test_build_single_salt(diafiltration_single_salt_Co_no_boundary_layer)
        test_build_single_salt_no_boundary_layer(
            diafiltration_single_salt_Co_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Co_no_boundary_layer(
        self, diafiltration_single_salt_Co_no_boundary_layer
    ):
        test_diagnostics(diafiltration_single_salt_Co_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co_no_boundary_layer(
        self, diafiltration_single_salt_Co_no_boundary_layer
    ):
        test_solve(diafiltration_single_salt_Co_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (10.468, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Co"): (227.59, 1e-4, None),
                (0, 1, "Cl"): (455.17, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (5.7687, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Co"): (215.62, 1e-4, None),
                (0, 1, "Cl"): (431.24, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_single_salt_Co_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Co_no_boundary_layer(
        self, diafiltration_single_salt_Co_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_single_salt_Co_no_boundary_layer)


################################################################################
# Test single-salt model: AlCl3


@pytest.fixture(scope="module")
def diafiltration_single_salt_Al():
    """
    Build single-salt diafiltration unit model for AlCl3
    with boundary layer.
    """
    cation_list = ["Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Al": 20, "Cl": 60},
        "diafiltrate": {"Al": 3, "Cl": 9},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Al(diafiltration_single_salt_Al):
    test_config(diafiltration_single_salt_Al)
    test_config_boundary_layer(diafiltration_single_salt_Al)
    test_config_single_salt(diafiltration_single_salt_Al)


class TestDiafiltrationSingleSaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Al(self, diafiltration_single_salt_Al):
        test_build(diafiltration_single_salt_Al)
        test_build_boundary_layer(diafiltration_single_salt_Al)
        test_build_single_salt(diafiltration_single_salt_Al)
        test_build_single_salt_boundary_layer(diafiltration_single_salt_Al)

    @pytest.mark.component
    def test_diagnostics_Al(self, diafiltration_single_salt_Al):
        test_diagnostics(diafiltration_single_salt_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Al(self, diafiltration_single_salt_Al):
        test_solve(diafiltration_single_salt_Al)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (9.5153, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Al"): (17.722, 1e-4, None),
                (0, 1, "Cl"): (53.166, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (6.6918, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Al"): (14.649, 1e-4, None),
                (0, 1, "Cl"): (43.947, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_single_salt_Al.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Al(self, diafiltration_single_salt_Al):
        test_numerical_issues(diafiltration_single_salt_Al)


@pytest.fixture(scope="module")
def diafiltration_single_salt_Al_no_boundary_layer():
    """
    Build single-salt diafiltration unit model for AlCl3
    without boundary layer.
    """
    cation_list = ["Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Al": 20, "Cl": 60},
        "diafiltrate": {"Al": 3, "Cl": 9},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure of single salts
    m.fs.unit.applied_pressure.fix(5)

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Al(diafiltration_single_salt_Al_no_boundary_layer):
    test_config(diafiltration_single_salt_Al_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_single_salt_Al_no_boundary_layer)
    test_config_single_salt(diafiltration_single_salt_Al_no_boundary_layer)


class TestDiafiltrationSingleSaltAluminumNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Al_no_boundary_layer(
        self, diafiltration_single_salt_Al_no_boundary_layer
    ):
        test_build(diafiltration_single_salt_Al_no_boundary_layer)
        test_build_single_salt(diafiltration_single_salt_Al_no_boundary_layer)
        test_build_single_salt_no_boundary_layer(
            diafiltration_single_salt_Al_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Al_no_boundary_layer(
        self, diafiltration_single_salt_Al_no_boundary_layer
    ):
        test_diagnostics(diafiltration_single_salt_Al_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Al_no_boundary_layer(
        self, diafiltration_single_salt_Al_no_boundary_layer
    ):
        test_solve(diafiltration_single_salt_Al_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (9.4910, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Al"): (18.044, 1e-4, None),
                (0, 1, "Cl"): (54.131, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (6.7079, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Al"): (14.372, 1e-4, None),
                (0, 1, "Cl"): (43.115, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_single_salt_Al_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Al_no_boundary_layer(
        self, diafiltration_single_salt_Al_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_single_salt_Al_no_boundary_layer)


################################################################################
# Test functions for two-salt model
@pytest.fixture(scope="module")
def sample_two_salt_model():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Cl": 821},
        "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.fixture(scope="module")
def sample_two_salt_model_no_boundary_layer():
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Cl": 821},
        "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.mark.unit
def test_config_two_salt(sample_two_salt_model):
    assert len(sample_two_salt_model.fs.unit.config.cation_list) == 2


@pytest.mark.build
@pytest.mark.unit
def test_build_two_salt(sample_two_salt_model):
    membrane = sample_two_salt_model.fs.unit
    # sets
    assert isinstance(membrane.solutes, Set)
    assert len(membrane.solutes) == 3

    assert isinstance(membrane.cations, Set)
    assert len(membrane.cations) == 2

    # variables
    assert isinstance(membrane.feed_conc_mol_comp, Var)
    assert len(membrane.feed_conc_mol_comp) == 3

    assert isinstance(membrane.diafiltrate_conc_mol_comp, Var)
    assert len(membrane.diafiltrate_conc_mol_comp) == 3

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 264

    assert isinstance(membrane.membrane_convection_coefficient_bilinear, Var)
    assert len(membrane.membrane_convection_coefficient_bilinear) == 132

    assert isinstance(membrane.membrane_cross_diffusion_coefficient, Var)
    assert len(membrane.membrane_cross_diffusion_coefficient) == 264

    assert isinstance(membrane.membrane_convection_coefficient, Var)
    assert len(membrane.membrane_convection_coefficient) == 132

    assert isinstance(membrane.molar_ion_flux, Var)
    assert len(membrane.molar_ion_flux) == 33

    assert isinstance(membrane.retentate_conc_mol_comp, Var)
    assert len(membrane.retentate_conc_mol_comp) == 33

    assert isinstance(membrane.permeate_conc_mol_comp, Var)
    assert len(membrane.permeate_conc_mol_comp) == 33

    assert isinstance(membrane.membrane_conc_mol_comp, Var)
    assert len(membrane.membrane_conc_mol_comp) == 198

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx) == 33

    assert isinstance(
        membrane.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_membrane_conc_mol_comp_dz) == 198

    # constraints
    assert isinstance(membrane.cation_mol_balance, Constraint)
    assert len(membrane.cation_mol_balance) == 20

    assert isinstance(membrane.cation_bulk_flux_equation, Constraint)
    assert len(membrane.cation_bulk_flux_equation) == 20

    assert isinstance(
        membrane.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_membrane_permeate_interface) == 20

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation) == 240
    )

    assert isinstance(
        membrane.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 120

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 240

    assert isinstance(
        membrane.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_calculation) == 120

    assert isinstance(membrane.cation_flux_membrane, Constraint)
    assert len(membrane.cation_flux_membrane) == 120

    assert isinstance(
        membrane.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.retentate_conc_mol_comp_boundary_condition) == 2

    assert isinstance(
        membrane.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 12

    assert isinstance(
        membrane.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 3

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx_boundary_condition) == 2

    assert isinstance(
        membrane.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(membrane.molar_ion_flux_boundary_condition) == 3

    # ports
    assert len(membrane.feed_inlet.conc_mol_comp) == 3
    assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 3
    assert len(membrane.retentate_outlet.conc_mol_comp) == 3
    assert len(membrane.permeate_outlet.conc_mol_comp) == 3


@pytest.mark.build
@pytest.mark.unit
def test_build_two_salt_boundary_layer(sample_two_salt_model):
    membrane = sample_two_salt_model.fs.unit
    # variables
    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 264

    assert isinstance(membrane.boundary_layer_cross_diffusion_coefficient, Var)
    assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 264

    assert isinstance(membrane.boundary_layer_conc_mol_comp, Var)
    assert len(membrane.boundary_layer_conc_mol_comp) == 198

    assert isinstance(
        membrane.d_boundary_layer_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 198

    # constraints
    assert isinstance(membrane.retentate_boundary_layer_interface, Constraint)
    assert len(membrane.retentate_boundary_layer_interface) == 20

    assert isinstance(
        membrane.cation_equilibrium_boundary_layer_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_boundary_layer_membrane_interface) == 20

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation)
        == 240
    )

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 240

    assert isinstance(membrane.cation_flux_boundary_layer, Constraint)
    assert len(membrane.cation_flux_boundary_layer) == 100

    assert isinstance(
        membrane.boundary_layer_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 12


@pytest.mark.build
@pytest.mark.unit
def test_build_two_salt_no_boundary_layer(sample_two_salt_model_no_boundary_layer):
    membrane = sample_two_salt_model_no_boundary_layer.fs.unit
    assert isinstance(
        membrane.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_retentate_membrane_interface) == 20


################################################################################
# Test a two-salt model: LiCl + CoCl2


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Co():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2
    with boundary layer.
    """
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Cl": 821},
        "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Co(diafiltration_two_salt_Li_Co):
    test_config(diafiltration_two_salt_Li_Co)
    test_config_boundary_layer(diafiltration_two_salt_Li_Co)
    test_config_two_salt(diafiltration_two_salt_Li_Co)


class TestDiafiltrationTwoSaltLithiumCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_build(diafiltration_two_salt_Li_Co)
        test_build_boundary_layer(diafiltration_two_salt_Li_Co)
        test_build_two_salt(diafiltration_two_salt_Li_Co)
        test_build_two_salt_boundary_layer(diafiltration_two_salt_Li_Co)

    @pytest.mark.component
    def test_diagnostics_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_diagnostics(diafiltration_two_salt_Li_Co)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_solve(diafiltration_two_salt_Li_Co)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (6.0854, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (190.89, 1e-4, None),
                (0, 1, "Co"): (239.83, 1e-4, None),
                (0, 1, "Cl"): (670.55, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (10.035, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (191.70, 1e-4, None),
                (0, 1, "Co"): (222.48, 1e-4, None),
                (0, 1, "Cl"): (636.67, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_two_salt_Li_Co.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_numerical_issues(diafiltration_two_salt_Li_Co)


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Co_no_boundary_layer():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2
    without boundary layer.
    """
    cation_list = ["Li", "Co"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Cl": 821},
        "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Co_no_boundary_layer(
    diafiltration_two_salt_Li_Co_no_boundary_layer,
):
    test_config(diafiltration_two_salt_Li_Co_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_two_salt_Li_Co_no_boundary_layer)
    test_config_two_salt(diafiltration_two_salt_Li_Co_no_boundary_layer)


class TestDiafiltrationTwoSaltLithiumCobaltNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co(self, diafiltration_two_salt_Li_Co_no_boundary_layer):
        test_build(diafiltration_two_salt_Li_Co_no_boundary_layer)
        test_build_two_salt(diafiltration_two_salt_Li_Co_no_boundary_layer)
        test_build_two_salt_no_boundary_layer(
            diafiltration_two_salt_Li_Co_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Li_Co_no_boundary_layer(
        self, diafiltration_two_salt_Li_Co_no_boundary_layer
    ):
        test_diagnostics(diafiltration_two_salt_Li_Co_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co_no_boundary_layer(
        self, diafiltration_two_salt_Li_Co_no_boundary_layer
    ):
        test_solve(diafiltration_two_salt_Li_Co_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (6.0465, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (188.88, 1e-4, None),
                (0, 1, "Co"): (246.67, 1e-4, None),
                (0, 1, "Cl"): (682.22, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (10.033, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (191.55, 1e-4, None),
                (0, 1, "Co"): (222.76, 1e-4, None),
                (0, 1, "Cl"): (637.06, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_two_salt_Li_Co_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Li_Co_no_boundary_layer(
        self, diafiltration_two_salt_Li_Co_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_two_salt_Li_Co_no_boundary_layer)


################################################################################
# Test a two-salt model: LiCl + AlCl3


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Al():
    """
    Build two-salt diafiltration unit model for LiCl + AlCl3
    with boundary layer.
    """
    cation_list = ["Li", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Al": 20, "Cl": 305},
        "diafiltrate": {"Li": 14, "Al": 3, "Cl": 23},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Al(diafiltration_two_salt_Li_Al):
    test_config(diafiltration_two_salt_Li_Al)
    test_config_boundary_layer(diafiltration_two_salt_Li_Al)
    test_config_two_salt(diafiltration_two_salt_Li_Al)


class TestDiafiltrationTwoSaltLithiumAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_build(diafiltration_two_salt_Li_Al)
        test_build_boundary_layer(diafiltration_two_salt_Li_Al)
        test_build_two_salt(diafiltration_two_salt_Li_Al)
        test_build_two_salt_boundary_layer(diafiltration_two_salt_Li_Al)

    @pytest.mark.component
    def test_diagnostics_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_diagnostics(diafiltration_two_salt_Li_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_solve(diafiltration_two_salt_Li_Al)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (6.2797, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (181.35, 1e-4, None),
                (0, 1, "Al"): (28.631, 1e-4, None),
                (0, 1, "Cl"): (267.24, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (9.0427, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (194.19, 1e-4, None),
                (0, 1, "Al"): (13.722, 1e-4, None),
                (0, 1, "Cl"): (235.36, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_two_salt_Li_Al.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_numerical_issues(diafiltration_two_salt_Li_Al)


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Al_no_boundary_layer():
    """
    Build two-salt diafiltration unit model for LiCl + AlCl3
    without boundary layer.
    """
    cation_list = ["Li", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Al": 20, "Cl": 305},
        "diafiltrate": {"Li": 14, "Al": 3, "Cl": 23},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Al_no_boundary_layer(
    diafiltration_two_salt_Li_Al_no_boundary_layer,
):
    test_config(diafiltration_two_salt_Li_Al_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_two_salt_Li_Al_no_boundary_layer)
    test_config_two_salt(diafiltration_two_salt_Li_Al_no_boundary_layer)


class TestDiafiltrationTwoSaltLithiumAluminumNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Al_no_boundary_layer(
        self, diafiltration_two_salt_Li_Al_no_boundary_layer
    ):
        test_build(diafiltration_two_salt_Li_Al_no_boundary_layer)
        test_build_two_salt(diafiltration_two_salt_Li_Al_no_boundary_layer)
        test_build_two_salt_no_boundary_layer(
            diafiltration_two_salt_Li_Al_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Li_Al_no_boundary_layer(
        self, diafiltration_two_salt_Li_Al_no_boundary_layer
    ):
        test_diagnostics(diafiltration_two_salt_Li_Al_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Al_no_boundary_layer(
        self, diafiltration_two_salt_Li_Al_no_boundary_layer
    ):
        test_solve(diafiltration_two_salt_Li_Al_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (5.6461, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (171.18, 1e-4, None),
                (0, 1, "Al"): (36.668, 1e-4, None),
                (0, 1, "Cl"): (281.18, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (9.0248, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (194.44, 1e-4, None),
                (0, 1, "Al"): (13.760, 1e-4, None),
                (0, 1, "Cl"): (235.72, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_two_salt_Li_Al_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Li_Al(
        self, diafiltration_two_salt_Li_Al_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_two_salt_Li_Al_no_boundary_layer)


################################################################################
# Test a two-salt model: CoCl2 + AlCl3


@pytest.fixture(scope="module")
def diafiltration_two_salt_Co_Al():
    """
    Build two-salt diafiltration unit model for CoCl2 + AlCl3
    with boundary layer.
    """
    cation_list = ["Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Co": 288, "Al": 20, "Cl": 636},
        "diafiltrate": {"Co": 3, "Al": 3, "Cl": 15},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Co_Al(diafiltration_two_salt_Co_Al):
    test_config(diafiltration_two_salt_Co_Al)
    test_config_boundary_layer(diafiltration_two_salt_Co_Al)
    test_config_two_salt(diafiltration_two_salt_Co_Al)


class TestDiafiltrationTwoSaltCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_build(diafiltration_two_salt_Co_Al)
        test_build_two_salt(diafiltration_two_salt_Co_Al)
        test_build_two_salt_boundary_layer(diafiltration_two_salt_Co_Al)

    @pytest.mark.component
    def test_diagnostics_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_diagnostics(diafiltration_two_salt_Co_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_solve(diafiltration_two_salt_Co_Al)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (8.3116, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Co"): (230.55, 1e-4, None),
                (0, 1, "Al"): (17.615, 1e-4, None),
                (0, 1, "Cl"): (513.94, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (7.8653, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Co"): (218.56, 1e-4, None),
                (0, 1, "Al"): (15.305, 1e-4, None),
                (0, 1, "Cl"): (483.04, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_two_salt_Co_Al.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_numerical_issues(diafiltration_two_salt_Co_Al)


@pytest.fixture(scope="module")
def diafiltration_two_salt_Co_Al_no_boundary_layer():
    """
    Build two-salt diafiltration unit model for CoCl2 + AlCl3
    without boundary layer.
    """
    cation_list = ["Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Co": 288, "Al": 20, "Cl": 636},
        "diafiltrate": {"Co": 3, "Al": 3, "Cl": 15},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Co_Al_no_boundary_layer(
    diafiltration_two_salt_Co_Al_no_boundary_layer,
):
    test_config(diafiltration_two_salt_Co_Al_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_two_salt_Co_Al_no_boundary_layer)
    test_config_two_salt(diafiltration_two_salt_Co_Al_no_boundary_layer)


class TestDiafiltrationTwoSaltCobaltAluminumNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co_Al_no_boundary_layer(
        self, diafiltration_two_salt_Co_Al_no_boundary_layer
    ):
        test_build(diafiltration_two_salt_Co_Al_no_boundary_layer)
        test_build_two_salt(diafiltration_two_salt_Co_Al_no_boundary_layer)
        test_build_two_salt_no_boundary_layer(
            diafiltration_two_salt_Co_Al_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Co_Al_no_boundary_layer(
        self, diafiltration_two_salt_Co_Al_no_boundary_layer
    ):
        test_diagnostics(diafiltration_two_salt_Co_Al_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co_Al_no_boundary_layer(
        self, diafiltration_two_salt_Co_Al_no_boundary_layer
    ):
        test_solve(diafiltration_two_salt_Co_Al_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (8.2716, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Co"): (232.23, 1e-4, None),
                (0, 1, "Al"): (18.276, 1e-4, None),
                (0, 1, "Cl"): (519.28, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (7.8847, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Co"): (218.01, 1e-4, None),
                (0, 1, "Al"): (14.954, 1e-4, None),
                (0, 1, "Cl"): (480.88, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_two_salt_Co_Al_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Co_Al_no_boundary_layer(
        self, diafiltration_two_salt_Co_Al_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_two_salt_Co_Al_no_boundary_layer)


################################################################################
# Test functions for three-salt model


@pytest.fixture(scope="module")
def sample_three_salt_model():
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
        "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.fixture(scope="module")
def sample_three_salt_model_no_boundary_layer():
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
        "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    return m


@pytest.mark.unit
def test_config_three_salt(sample_three_salt_model):
    assert len(sample_three_salt_model.fs.unit.config.cation_list) == 3


@pytest.mark.build
@pytest.mark.unit
def test_build_three_salt(sample_three_salt_model):
    membrane = sample_three_salt_model.fs.unit
    assert isinstance(membrane.solutes, Set)
    assert len(membrane.solutes) == 4

    assert isinstance(membrane.cations, Set)
    assert len(membrane.cations) == 3

    # variables
    assert isinstance(membrane.feed_conc_mol_comp, Var)
    assert len(membrane.feed_conc_mol_comp) == 4

    assert isinstance(membrane.diafiltrate_conc_mol_comp, Var)
    assert len(membrane.diafiltrate_conc_mol_comp) == 4

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 594

    assert isinstance(membrane.membrane_convection_coefficient_bilinear, Var)
    assert len(membrane.membrane_convection_coefficient_bilinear) == 198

    assert isinstance(membrane.membrane_cross_diffusion_coefficient, Var)
    assert len(membrane.membrane_cross_diffusion_coefficient) == 594

    assert isinstance(membrane.membrane_convection_coefficient, Var)
    assert len(membrane.membrane_convection_coefficient) == 198

    assert isinstance(membrane.molar_ion_flux, Var)
    assert len(membrane.molar_ion_flux) == 44

    assert isinstance(membrane.retentate_conc_mol_comp, Var)
    assert len(membrane.retentate_conc_mol_comp) == 44

    assert isinstance(membrane.permeate_conc_mol_comp, Var)
    assert len(membrane.permeate_conc_mol_comp) == 44

    assert isinstance(membrane.membrane_conc_mol_comp, Var)
    assert len(membrane.membrane_conc_mol_comp) == 264

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx) == 44

    assert isinstance(
        membrane.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_membrane_conc_mol_comp_dz) == 264

    # constraints
    assert isinstance(membrane.cation_mol_balance, Constraint)
    assert len(membrane.cation_mol_balance) == 30

    assert isinstance(membrane.cation_bulk_flux_equation, Constraint)
    assert len(membrane.cation_bulk_flux_equation) == 30

    assert isinstance(
        membrane.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_membrane_permeate_interface) == 30

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation) == 540
    )

    assert isinstance(
        membrane.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 180

    assert isinstance(
        membrane.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 540

    assert isinstance(
        membrane.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.membrane_convection_coefficient_calculation) == 180

    assert isinstance(membrane.cation_flux_membrane, Constraint)
    assert len(membrane.cation_flux_membrane) == 180

    assert isinstance(
        membrane.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.retentate_conc_mol_comp_boundary_condition) == 3

    assert isinstance(
        membrane.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 18

    assert isinstance(
        membrane.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 4

    assert isinstance(
        membrane.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert len(membrane.d_retentate_conc_mol_comp_dx_boundary_condition) == 3

    assert isinstance(
        membrane.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(membrane.molar_ion_flux_boundary_condition) == 4

    # ports
    assert len(membrane.feed_inlet.conc_mol_comp) == 4
    assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 4
    assert len(membrane.retentate_outlet.conc_mol_comp) == 4
    assert len(membrane.permeate_outlet.conc_mol_comp) == 4


@pytest.mark.build
@pytest.mark.unit
def test_build_three_salt_boundary_layer(sample_three_salt_model):
    membrane = sample_three_salt_model.fs.unit
    # variables
    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 594

    assert isinstance(membrane.boundary_layer_cross_diffusion_coefficient, Var)
    assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 594

    assert isinstance(membrane.boundary_layer_conc_mol_comp, Var)
    assert len(membrane.boundary_layer_conc_mol_comp) == 264

    assert isinstance(
        membrane.d_boundary_layer_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 264

    # constraints
    assert isinstance(membrane.retentate_boundary_layer_interface, Constraint)
    assert len(membrane.retentate_boundary_layer_interface) == 30

    assert isinstance(
        membrane.cation_equilibrium_boundary_layer_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_boundary_layer_membrane_interface) == 30

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation)
        == 540
    )

    assert isinstance(
        membrane.boundary_layer_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 540

    assert isinstance(membrane.cation_flux_boundary_layer, Constraint)
    assert len(membrane.cation_flux_boundary_layer) == 150

    assert isinstance(
        membrane.boundary_layer_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 18


@pytest.mark.build
@pytest.mark.unit
def test_build_three_salt_no_boundary_layer(sample_three_salt_model_no_boundary_layer):
    membrane = sample_three_salt_model_no_boundary_layer.fs.unit
    assert isinstance(
        membrane.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert len(membrane.cation_equilibrium_retentate_membrane_interface) == 30


################################################################################
# Test a three-salt model: LiCl + CoCl2 + AlCl3


@pytest.fixture(scope="module")
def diafiltration_three_salt_Li_Co_Al():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2
    + AlCl3 with boundary layer.
    """
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
        "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=True,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 11

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Co_Al(diafiltration_three_salt_Li_Co_Al):
    test_config(diafiltration_three_salt_Li_Co_Al)
    test_config_boundary_layer(diafiltration_three_salt_Li_Co_Al)
    test_config_three_salt(diafiltration_three_salt_Li_Co_Al)


class TestDiafiltrationThreeSaltLithiumCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_build(diafiltration_three_salt_Li_Co_Al)
        test_build_boundary_layer(diafiltration_three_salt_Li_Co_Al)
        test_build_three_salt(diafiltration_three_salt_Li_Co_Al)
        test_build_three_salt_boundary_layer(diafiltration_three_salt_Li_Co_Al)

    @pytest.mark.component
    def test_diagnostics_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_diagnostics(diafiltration_three_salt_Li_Co_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_solve(diafiltration_three_salt_Li_Co_Al)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (10.528, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (189.94, 1e-4, None),
                (0, 1, "Co"): (224.29, 1e-4, None),
                (0, 1, "Al"): (20.430, 1e-4, None),
                (0, 1, "Cl"): (699.80, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (5.5393, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (194.13, 1e-4, None),
                (0, 1, "Co"): (220.51, 1e-4, None),
                (0, 1, "Al"): (9.7103, 1e-4, None),
                (0, 1, "Cl"): (664.29, 1e-4, None),
            },
        }

        assert_solution_equivalent(diafiltration_three_salt_Li_Co_Al.fs.unit, test_dict)

    @pytest.mark.component
    def test_numerical_issues_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_numerical_issues(diafiltration_three_salt_Li_Co_Al)


@pytest.fixture(scope="module")
def diafiltration_three_salt_Li_Co_Al_no_boundary_layer():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2
    + AlCl3 without boundary layer.
    """
    cation_list = ["Li", "Co", "Al"]
    anion_list = ["Cl"]
    inlet_flow_volume = {"feed": 12.5, "diafiltrate": 3.75}
    inlet_concentration = {
        "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
        "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        cation_list=cation_list,
        anion_list=anion_list,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        cation_list=cation_list,
        anion_list=anion_list,
        include_boundary_layer=False,
        NFE_module_length=10,
        NFE_boundary_layer_thickness=5,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 11

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
    m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

    for t in m.fs.unit.time:
        for j in m.fs.unit.solutes:
            m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
            m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                inlet_concentration["diafiltrate"][j]
            )

    initialized_model = m.fs.unit.default_initializer()
    initialized_model.initialize(m.fs.unit)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_Li_Co_Al_no_boundary_layer(
    diafiltration_three_salt_Li_Co_Al_no_boundary_layer,
):
    test_config(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)
    test_config_no_boundary_layer(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)
    test_config_three_salt(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)


class TestDiafiltrationThreeSaltLithiumCobaltAluminumNoBoundaryLayer(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co_Al_no_boundary_layer(
        self, diafiltration_three_salt_Li_Co_Al_no_boundary_layer
    ):
        test_build(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)
        test_build_three_salt(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)
        test_build_three_salt_no_boundary_layer(
            diafiltration_three_salt_Li_Co_Al_no_boundary_layer
        )

    @pytest.mark.component
    def test_diagnostics_Li_Co_Al_no_boundary_layer(
        self, diafiltration_three_salt_Li_Co_Al_no_boundary_layer
    ):
        test_diagnostics(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co_Al_no_boundary_layer(
        self, diafiltration_three_salt_Li_Co_Al_no_boundary_layer
    ):
        test_solve(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)

        test_dict = {
            "retentate_flow_volume": {(0, 1): (10.295, 1e-4, None)},
            "retentate_conc_mol_comp": {
                (0, 1, "Li"): (189.13, 1e-4, None),
                (0, 1, "Co"): (224.68, 1e-4, None),
                (0, 1, "Al"): (21.673, 1e-4, None),
                (0, 1, "Cl"): (703.51, 1e-4, None),
            },
            "permeate_flow_volume": {(0, 1): (5.7276, 1e-4, None)},
            "permeate_conc_mol_comp": {
                (0, 1, "Li"): (194.91, 1e-4, None),
                (0, 1, "Co"): (220.89, 1e-4, None),
                (0, 1, "Al"): (8.3159, 1e-4, None),
                (0, 1, "Cl"): (661.63, 1e-4, None),
            },
        }

        assert_solution_equivalent(
            diafiltration_three_salt_Li_Co_Al_no_boundary_layer.fs.unit, test_dict
        )

    @pytest.mark.component
    def test_numerical_issues_Li_Co_Al_no_boundary_layer(
        self, diafiltration_three_salt_Li_Co_Al_no_boundary_layer
    ):
        test_numerical_issues(diafiltration_three_salt_Li_Co_Al_no_boundary_layer)


################################################################################
