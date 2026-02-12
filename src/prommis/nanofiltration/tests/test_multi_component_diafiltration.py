#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
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

# TODO: test positive and neutral membrane cases
# currently, the property package only supports a negative fixed charge


################################################################################
# Test functions for single-salt model
@pytest.fixture(scope="module")
def sample_single_salt_model():
    cation_list = ["lithium"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "chloride": 245},
            "diafiltrate": {"lithium": 14, "chloride": 14},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    return m


@pytest.mark.unit
def test_config_single_salt(sample_single_salt_model):
    assert len(sample_single_salt_model.fs.unit.config) == 10

    assert not sample_single_salt_model.fs.unit.config.dynamic
    assert not sample_single_salt_model.fs.unit.config.has_holdup

    assert (
        sample_single_salt_model.fs.unit.config.property_package
        is sample_single_salt_model.fs.properties
    )
    assert len(sample_single_salt_model.fs.unit.config.cation_list) == 1
    assert len(sample_single_salt_model.fs.unit.config.anion_list) == 1
    assert sample_single_salt_model.fs.unit.config.NFE_module_length == 10
    assert sample_single_salt_model.fs.unit.config.NFE_membrane_thickness == 5


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt(sample_single_salt_model):
    # parameters
    assert isinstance(sample_single_salt_model.fs.unit.total_membrane_thickness, Param)
    assert value(sample_single_salt_model.fs.unit.total_membrane_thickness) == 1e-7

    assert isinstance(sample_single_salt_model.fs.unit.membrane_fixed_charge, Param)
    assert value(sample_single_salt_model.fs.unit.membrane_fixed_charge) == -44

    assert isinstance(sample_single_salt_model.fs.unit.membrane_permeability, Param)
    assert value(sample_single_salt_model.fs.unit.membrane_permeability) == 0.01

    assert isinstance(sample_single_salt_model.fs.unit.temperature, Param)
    assert value(sample_single_salt_model.fs.unit.temperature) == 298

    # sets
    assert isinstance(
        sample_single_salt_model.fs.unit.dimensionless_module_length, ContinuousSet
    )
    assert len(sample_single_salt_model.fs.unit.dimensionless_module_length) == 11

    assert isinstance(
        sample_single_salt_model.fs.unit.dimensionless_membrane_thickness,
        ContinuousSet,
    )
    assert len(sample_single_salt_model.fs.unit.dimensionless_membrane_thickness) == 6

    assert isinstance(sample_single_salt_model.fs.unit.time, Set)
    assert len(sample_single_salt_model.fs.unit.time) == 1

    assert isinstance(sample_single_salt_model.fs.unit.solutes, Set)
    assert len(sample_single_salt_model.fs.unit.solutes) == 2

    assert isinstance(sample_single_salt_model.fs.unit.cations, Set)
    assert len(sample_single_salt_model.fs.unit.cations) == 1

    # variables
    assert isinstance(sample_single_salt_model.fs.unit.total_module_length, Var)
    assert len(sample_single_salt_model.fs.unit.total_module_length) == 1

    assert isinstance(sample_single_salt_model.fs.unit.total_membrane_length, Var)
    assert len(sample_single_salt_model.fs.unit.total_membrane_length) == 1

    assert isinstance(sample_single_salt_model.fs.unit.applied_pressure, Var)
    assert len(sample_single_salt_model.fs.unit.applied_pressure) == 1

    assert isinstance(sample_single_salt_model.fs.unit.feed_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.feed_flow_volume) == 1

    assert isinstance(sample_single_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.feed_conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_flow_volume) == 1

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.membrane_D_tilde, Var)
    assert len(sample_single_salt_model.fs.unit.membrane_D_tilde) == 66

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear
        )
        == 66
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_convection_coefficient_bilinear,
        Var,
    )
    assert (
        len(sample_single_salt_model.fs.unit.membrane_convection_coefficient_bilinear)
        == 66
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient, Var
    )
    assert (
        len(sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient) == 66
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_convection_coefficient, Var
    )
    assert len(sample_single_salt_model.fs.unit.membrane_convection_coefficient) == 66

    assert isinstance(sample_single_salt_model.fs.unit.volume_flux_water, Var)
    assert len(sample_single_salt_model.fs.unit.volume_flux_water) == 11

    assert isinstance(sample_single_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_single_salt_model.fs.unit.molar_ion_flux) == 22

    assert isinstance(sample_single_salt_model.fs.unit.retentate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.retentate_flow_volume) == 11

    assert isinstance(sample_single_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.retentate_conc_mol_comp) == 22

    assert isinstance(sample_single_salt_model.fs.unit.permeate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.permeate_flow_volume) == 11

    assert isinstance(sample_single_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.permeate_conc_mol_comp) == 22

    assert isinstance(sample_single_salt_model.fs.unit.osmotic_pressure, Var)
    assert len(sample_single_salt_model.fs.unit.osmotic_pressure) == 11

    assert isinstance(sample_single_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.membrane_conc_mol_comp) == 132

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 22

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
    )
    assert len(sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx) == 11

    assert isinstance(
        sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 132

    # constraints
    assert isinstance(sample_single_salt_model.fs.unit.overall_mol_balance, Constraint)
    assert len(sample_single_salt_model.fs.unit.overall_mol_balance) == 10

    assert isinstance(sample_single_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_single_salt_model.fs.unit.cation_mol_balance) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.overall_bulk_flux_equation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.overall_bulk_flux_equation) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.cation_bulk_flux_equation) == 10

    assert isinstance(sample_single_salt_model.fs.unit.lumped_water_flux, Constraint)
    assert len(sample_single_salt_model.fs.unit.lumped_water_flux) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.chloride_flux_membrane, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.chloride_flux_membrane) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.osmotic_pressure_calculation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.osmotic_pressure_calculation) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.electroneutrality_retentate, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.electroneutrality_retentate) == 11

    assert isinstance(
        sample_single_salt_model.fs.unit.electroneutrality_permeate, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.electroneutrality_permeate) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface
        )
        == 10
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface
        )
        == 10
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_D_tilde_calculation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.membrane_D_tilde_calculation) == 60

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation
        )
        == 60
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation
        )
        == 60
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation
        )
        == 60
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.membrane_convection_coefficient_calculation
        )
        == 60
    )

    assert isinstance(sample_single_salt_model.fs.unit.cation_flux_membrane, Constraint)
    assert len(sample_single_salt_model.fs.unit.cation_flux_membrane) == 60

    assert isinstance(
        sample_single_salt_model.fs.unit.electroneutrality_membrane, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.electroneutrality_membrane) == 60

    assert isinstance(
        sample_single_salt_model.fs.unit.retentate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.retentate_flow_volume_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition)
        == 6
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.permeate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.permeate_flow_volume_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 2
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition
        )
        == 1
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(
            sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition
        )
        == 1
    )
    assert isinstance(
        sample_single_salt_model.fs.unit.volume_flux_water_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.volume_flux_water_boundary_condition) == 1
    )

    assert isinstance(
        sample_single_salt_model.fs.unit.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(sample_single_salt_model.fs.unit.molar_ion_flux_boundary_condition) == 2

    for t in sample_single_salt_model.fs.unit.time:
        for x in sample_single_salt_model.fs.unit.dimensionless_module_length:
            assert sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx[
                t, x, "chloride"
            ].fixed
            if x != 0:
                assert not sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, "chloride"
                ].active

            for z in sample_single_salt_model.fs.unit.dimensionless_membrane_thickness:
                assert sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz[
                    t, x, z, "chloride"
                ].fixed
                if z != 0:
                    assert not sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, "chloride"
                    ].active

    # scaling factors
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.volume_flux_water
        ]
        == 1e2
    )
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.membrane_D_tilde
        ]
        == 1e-1
    )
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear
        ]
        == 1e-2
    )
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.membrane_convection_coefficient_bilinear
        ]
        == 1e-1
    )
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.membrane_cross_diffusion_coefficient
        ]
        == 1e1
    )
    assert (
        sample_single_salt_model.fs.unit.scaling_factor[
            sample_single_salt_model.fs.unit.membrane_convection_coefficient
        ]
        == 1e1
    )

    # ports
    assert isinstance(sample_single_salt_model.fs.unit.feed_inlet, Port)
    assert len(sample_single_salt_model.fs.unit.feed_inlet.flow_vol) == 1
    assert len(sample_single_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_inlet, Port)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_inlet.flow_vol) == 1
    assert len(sample_single_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.retentate_outlet, Port)
    assert len(sample_single_salt_model.fs.unit.retentate_outlet.flow_vol) == 1
    assert len(sample_single_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.permeate_outlet, Port)
    assert len(sample_single_salt_model.fs.unit.permeate_outlet.flow_vol) == 1
    assert len(sample_single_salt_model.fs.unit.permeate_outlet.conc_mol_comp) == 2


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
# Test single-salt model: lithium chloride


@pytest.fixture(scope="module")
def diafiltration_single_salt_lithium():
    """
    Build single-salt diafiltration unit model for lithium chloride.
    """
    cation_list = ["lithium"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "chloride": 245},
            "diafiltrate": {"lithium": 14, "chloride": 14},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure
    m.fs.unit.applied_pressure.fix(5)
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_lithium(diafiltration_single_salt_lithium):
    test_config_single_salt(diafiltration_single_salt_lithium)


class TestDiafiltrationSingleSaltLithium(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_lithium(self, diafiltration_single_salt_lithium):
        assert isinstance(
            diafiltration_single_salt_lithium.fs.unit.numerical_zero_tolerance, Param
        )
        assert (
            value(diafiltration_single_salt_lithium.fs.unit.numerical_zero_tolerance)
            == 1e-10
        )

        test_build_single_salt(diafiltration_single_salt_lithium)

    @pytest.mark.component
    def test_diagnostics_lithium(self, diafiltration_single_salt_lithium):
        test_diagnostics(diafiltration_single_salt_lithium)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_lithium(self, diafiltration_single_salt_lithium):
        test_solve(diafiltration_single_salt_lithium)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                8.3856,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.retentate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                194.52,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                194.52,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_flow_volume[0, 1]
                ),
                7.8637,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                190.38,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                190.38,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_lithium(self, diafiltration_single_salt_lithium):
        test_numerical_issues(diafiltration_single_salt_lithium)


################################################################################
# Test single-salt model: cobalt chloride


@pytest.fixture(scope="module")
def diafiltration_single_salt_cobalt():
    """
    Build single-salt diafiltration unit model for cobalt chloride.
    """
    cation_list = ["cobalt"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"cobalt": 288, "chloride": 576},
            "diafiltrate": {"cobalt": 3, "chloride": 6},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure
    m.fs.unit.applied_pressure.fix(5)
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_cobalt(diafiltration_single_salt_cobalt):
    test_config_single_salt(diafiltration_single_salt_cobalt)


class TestDiafiltrationSingleSaltCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_cobalt(self, diafiltration_single_salt_cobalt):
        assert isinstance(
            diafiltration_single_salt_cobalt.fs.unit.numerical_zero_tolerance, Param
        )
        assert (
            value(diafiltration_single_salt_cobalt.fs.unit.numerical_zero_tolerance)
            == 1e-10
        )

        test_build_single_salt(diafiltration_single_salt_cobalt)

    @pytest.mark.component
    def test_diagnostics_cobalt(self, diafiltration_single_salt_cobalt):
        test_diagnostics(diafiltration_single_salt_cobalt)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_cobalt(self, diafiltration_single_salt_cobalt):
        test_solve(diafiltration_single_salt_cobalt)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.retentate_flow_volume[0, 1]
                ),
                10.468,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                227.59,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                455.17,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_flow_volume[0, 1]
                ),
                5.7687,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                215.62,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                431.24,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_cobalt(self, diafiltration_single_salt_cobalt):
        test_numerical_issues(diafiltration_single_salt_cobalt)


################################################################################
# Test single-salt model: aluminum chloride


@pytest.fixture(scope="module")
def diafiltration_single_salt_aluminum():
    """
    Build single-salt diafiltration unit model for aluminum chloride.
    """
    cation_list = ["aluminum"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"aluminum": 20, "chloride": 60},
            "diafiltrate": {"aluminum": 3, "chloride": 9},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 7

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    # reduce pressure to accommodate lower osmotic pressure
    m.fs.unit.applied_pressure.fix(5)
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_aluminum(diafiltration_single_salt_aluminum):
    test_config_single_salt(diafiltration_single_salt_aluminum)


class TestDiafiltrationSingleSaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_aluminum(self, diafiltration_single_salt_aluminum):
        assert isinstance(
            diafiltration_single_salt_aluminum.fs.unit.numerical_zero_tolerance, Param
        )
        assert (
            value(diafiltration_single_salt_aluminum.fs.unit.numerical_zero_tolerance)
            == 1e-10
        )

        test_build_single_salt(diafiltration_single_salt_aluminum)

    @pytest.mark.component
    def test_diagnostics_aluminum(self, diafiltration_single_salt_aluminum):
        test_diagnostics(diafiltration_single_salt_aluminum)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_aluminum(self, diafiltration_single_salt_aluminum):
        test_solve(diafiltration_single_salt_aluminum)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                9.4910,
            ],
            "aluminum_retentate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                18.044,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                54.131,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                6.7079,
            ],
            "aluminum_permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                14.372,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                43.115,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_aluminum(self, diafiltration_single_salt_aluminum):
        test_numerical_issues(diafiltration_single_salt_aluminum)


################################################################################
# Test functions for two-salt model
@pytest.fixture(scope="module")
def sample_two_salt_model():
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "cobalt": 288, "chloride": 822},
            "diafiltrate": {"lithium": 14, "cobalt": 3, "chloride": 21},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    return m


@pytest.mark.unit
def test_config_two_salt(sample_two_salt_model):
    assert len(sample_two_salt_model.fs.unit.config) == 10

    assert not sample_two_salt_model.fs.unit.config.dynamic
    assert not sample_two_salt_model.fs.unit.config.has_holdup

    assert (
        sample_two_salt_model.fs.unit.config.property_package
        is sample_two_salt_model.fs.properties
    )
    assert len(sample_two_salt_model.fs.unit.config.cation_list) == 2
    assert len(sample_two_salt_model.fs.unit.config.anion_list) == 1
    assert sample_two_salt_model.fs.unit.config.NFE_module_length == 10
    assert sample_two_salt_model.fs.unit.config.NFE_membrane_thickness == 5


@pytest.mark.build
@pytest.mark.unit
def test_build_two_salt(sample_two_salt_model):
    # parameters
    assert isinstance(sample_two_salt_model.fs.unit.total_membrane_thickness, Param)
    assert value(sample_two_salt_model.fs.unit.total_membrane_thickness) == 1e-7

    assert isinstance(sample_two_salt_model.fs.unit.membrane_fixed_charge, Param)
    assert value(sample_two_salt_model.fs.unit.membrane_fixed_charge) == -44

    assert isinstance(sample_two_salt_model.fs.unit.membrane_permeability, Param)
    assert value(sample_two_salt_model.fs.unit.membrane_permeability) == 0.01

    assert isinstance(sample_two_salt_model.fs.unit.temperature, Param)
    assert value(sample_two_salt_model.fs.unit.temperature) == 298

    # sets
    assert isinstance(
        sample_two_salt_model.fs.unit.dimensionless_module_length, ContinuousSet
    )
    assert len(sample_two_salt_model.fs.unit.dimensionless_module_length) == 11

    assert isinstance(
        sample_two_salt_model.fs.unit.dimensionless_membrane_thickness,
        ContinuousSet,
    )
    assert len(sample_two_salt_model.fs.unit.dimensionless_membrane_thickness) == 6

    assert isinstance(sample_two_salt_model.fs.unit.time, Set)
    assert len(sample_two_salt_model.fs.unit.time) == 1

    assert isinstance(sample_two_salt_model.fs.unit.solutes, Set)
    assert len(sample_two_salt_model.fs.unit.solutes) == 3

    assert isinstance(sample_two_salt_model.fs.unit.cations, Set)
    assert len(sample_two_salt_model.fs.unit.cations) == 2

    # variables
    assert isinstance(sample_two_salt_model.fs.unit.total_module_length, Var)
    assert len(sample_two_salt_model.fs.unit.total_module_length) == 1

    assert isinstance(sample_two_salt_model.fs.unit.total_membrane_length, Var)
    assert len(sample_two_salt_model.fs.unit.total_membrane_length) == 1

    assert isinstance(sample_two_salt_model.fs.unit.applied_pressure, Var)
    assert len(sample_two_salt_model.fs.unit.applied_pressure) == 1

    assert isinstance(sample_two_salt_model.fs.unit.feed_flow_volume, Var)
    assert len(sample_two_salt_model.fs.unit.feed_flow_volume) == 1

    assert isinstance(sample_two_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.feed_conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.diafiltrate_flow_volume, Var)
    assert len(sample_two_salt_model.fs.unit.diafiltrate_flow_volume) == 1

    assert isinstance(sample_two_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.membrane_D_tilde, Var)
    assert len(sample_two_salt_model.fs.unit.membrane_D_tilde) == 66

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert (
        len(sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear)
        == 264
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_convection_coefficient_bilinear, Var
    )
    assert (
        len(sample_two_salt_model.fs.unit.membrane_convection_coefficient_bilinear)
        == 132
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient, Var
    )
    assert (
        len(sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient) == 264
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_convection_coefficient, Var
    )
    assert len(sample_two_salt_model.fs.unit.membrane_convection_coefficient) == 132

    assert isinstance(sample_two_salt_model.fs.unit.volume_flux_water, Var)
    assert len(sample_two_salt_model.fs.unit.volume_flux_water) == 11

    assert isinstance(sample_two_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_two_salt_model.fs.unit.molar_ion_flux) == 33

    assert isinstance(sample_two_salt_model.fs.unit.retentate_flow_volume, Var)
    assert len(sample_two_salt_model.fs.unit.retentate_flow_volume) == 11

    assert isinstance(sample_two_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.retentate_conc_mol_comp) == 33

    assert isinstance(sample_two_salt_model.fs.unit.permeate_flow_volume, Var)
    assert len(sample_two_salt_model.fs.unit.permeate_flow_volume) == 11

    assert isinstance(sample_two_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.permeate_conc_mol_comp) == 33

    assert isinstance(sample_two_salt_model.fs.unit.osmotic_pressure, Var)
    assert len(sample_two_salt_model.fs.unit.osmotic_pressure) == 11

    assert isinstance(sample_two_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.membrane_conc_mol_comp) == 198

    assert isinstance(
        sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 33

    assert isinstance(
        sample_two_salt_model.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
    )
    assert len(sample_two_salt_model.fs.unit.d_retentate_flow_volume_dx) == 11

    assert isinstance(
        sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 198

    # constraints
    assert isinstance(sample_two_salt_model.fs.unit.overall_mol_balance, Constraint)
    assert len(sample_two_salt_model.fs.unit.overall_mol_balance) == 10

    assert isinstance(sample_two_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_two_salt_model.fs.unit.cation_mol_balance) == 20

    assert isinstance(
        sample_two_salt_model.fs.unit.overall_bulk_flux_equation, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.overall_bulk_flux_equation) == 10

    assert isinstance(
        sample_two_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.cation_bulk_flux_equation) == 20

    assert isinstance(sample_two_salt_model.fs.unit.lumped_water_flux, Constraint)
    assert len(sample_two_salt_model.fs.unit.lumped_water_flux) == 10

    assert isinstance(sample_two_salt_model.fs.unit.chloride_flux_membrane, Constraint)
    assert len(sample_two_salt_model.fs.unit.chloride_flux_membrane) == 10

    assert isinstance(
        sample_two_salt_model.fs.unit.osmotic_pressure_calculation, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.osmotic_pressure_calculation) == 10

    assert isinstance(
        sample_two_salt_model.fs.unit.electroneutrality_retentate, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.electroneutrality_retentate) == 11

    assert isinstance(
        sample_two_salt_model.fs.unit.electroneutrality_permeate, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.electroneutrality_permeate) == 10

    assert isinstance(
        sample_two_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface
        )
        == 20
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface
        )
        == 20
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_D_tilde_calculation, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.membrane_D_tilde_calculation) == 60

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation
        )
        == 240
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation
        )
        == 120
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation
        )
        == 240
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.membrane_convection_coefficient_calculation)
        == 120
    )

    assert isinstance(sample_two_salt_model.fs.unit.cation_flux_membrane, Constraint)
    assert len(sample_two_salt_model.fs.unit.cation_flux_membrane) == 120

    assert isinstance(
        sample_two_salt_model.fs.unit.electroneutrality_membrane, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.electroneutrality_membrane) == 60

    assert isinstance(
        sample_two_salt_model.fs.unit.retentate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.retentate_flow_volume_boundary_condition) == 1
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition)
        == 2
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition)
        == 12
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.permeate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.permeate_flow_volume_boundary_condition) == 1
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 3
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(
            sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition
        )
        == 2
    )
    assert isinstance(
        sample_two_salt_model.fs.unit.volume_flux_water_boundary_condition,
        Constraint,
    )
    assert len(sample_two_salt_model.fs.unit.volume_flux_water_boundary_condition) == 1

    assert isinstance(
        sample_two_salt_model.fs.unit.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(sample_two_salt_model.fs.unit.molar_ion_flux_boundary_condition) == 3

    for t in sample_two_salt_model.fs.unit.time:
        for x in sample_two_salt_model.fs.unit.dimensionless_module_length:
            assert sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx[
                t, x, "chloride"
            ].fixed
            if x != 0:
                assert not sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, "chloride"
                ].active

            for z in sample_two_salt_model.fs.unit.dimensionless_membrane_thickness:
                assert sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz[
                    t, x, z, "chloride"
                ].fixed
                if z != 0:
                    assert not sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, "chloride"
                    ].active

    # scaling factors
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.volume_flux_water
        ]
        == 1e2
    )
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.membrane_D_tilde
        ]
        == 1e-1
    )
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear
        ]
        == 1e-2
    )
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.membrane_convection_coefficient_bilinear
        ]
        == 1e-1
    )
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.membrane_cross_diffusion_coefficient
        ]
        == 1e1
    )
    assert (
        sample_two_salt_model.fs.unit.scaling_factor[
            sample_two_salt_model.fs.unit.membrane_convection_coefficient
        ]
        == 1e1
    )

    # ports
    assert isinstance(sample_two_salt_model.fs.unit.feed_inlet, Port)
    assert len(sample_two_salt_model.fs.unit.feed_inlet.flow_vol) == 1
    assert len(sample_two_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.diafiltrate_inlet, Port)
    assert len(sample_two_salt_model.fs.unit.diafiltrate_inlet.flow_vol) == 1
    assert len(sample_two_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.retentate_outlet, Port)
    assert len(sample_two_salt_model.fs.unit.retentate_outlet.flow_vol) == 1
    assert len(sample_two_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.permeate_outlet, Port)
    assert len(sample_two_salt_model.fs.unit.permeate_outlet.flow_vol) == 1
    assert len(sample_two_salt_model.fs.unit.permeate_outlet.conc_mol_comp) == 3


################################################################################
# Test a two-salt model: lithium chloride + cobalt chloride


@pytest.fixture(scope="module")
def diafiltration_two_salt_lithium_cobalt():
    """
    Build two-salt diafiltration unit model for lithium chloride + cobalt chloride.
    """
    cation_list = ["lithium", "cobalt"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "cobalt": 288, "chloride": 822},
            "diafiltrate": {"lithium": 14, "cobalt": 3, "chloride": 21},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_lithium_cobalt(diafiltration_two_salt_lithium_cobalt):
    test_config_two_salt(diafiltration_two_salt_lithium_cobalt)


class TestDiafiltrationTwoSaltLithiumCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_lithium_cobalt(self, diafiltration_two_salt_lithium_cobalt):
        assert isinstance(
            diafiltration_two_salt_lithium_cobalt.fs.unit.numerical_zero_tolerance,
            Param,
        )
        assert (
            value(
                diafiltration_two_salt_lithium_cobalt.fs.unit.numerical_zero_tolerance
            )
            == 1e-10
        )

        test_build_two_salt(diafiltration_two_salt_lithium_cobalt)

    @pytest.mark.component
    def test_diagnostics_lithium_cobalt(self, diafiltration_two_salt_lithium_cobalt):
        test_diagnostics(diafiltration_two_salt_lithium_cobalt)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_lithium_cobalt(self, diafiltration_two_salt_lithium_cobalt):
        test_solve(diafiltration_two_salt_lithium_cobalt)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                6.0465,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                188.88,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                246.67,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                682.22,
            ],
            "permeate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                10.033,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                191.55,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                222.76,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                637.06,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_lithium_cobalt(
        self, diafiltration_two_salt_lithium_cobalt
    ):
        test_numerical_issues(diafiltration_two_salt_lithium_cobalt)


################################################################################
# Test a two-salt model: lithium chloride + aluminum chloride


@pytest.fixture(scope="module")
def diafiltration_two_salt_lithium_aluminum():
    """
    Build two-salt diafiltration unit model for lithium chloride + aluminum chloride.
    """
    cation_list = ["lithium", "aluminum"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "aluminum": 20, "chloride": 305},
            "diafiltrate": {"lithium": 14, "aluminum": 3, "chloride": 23},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    # reduce numerical tolerance
    m.fs.unit.numerical_zero_tolerance.set_value(1e-8)

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_lithium_aluminum(diafiltration_two_salt_lithium_aluminum):
    test_config_two_salt(diafiltration_two_salt_lithium_aluminum)


class TestDiafiltrationTwoSaltLithiumAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_lithium_aluminum(self, diafiltration_two_salt_lithium_aluminum):
        assert isinstance(
            diafiltration_two_salt_lithium_aluminum.fs.unit.numerical_zero_tolerance,
            Param,
        )
        assert (
            value(
                diafiltration_two_salt_lithium_aluminum.fs.unit.numerical_zero_tolerance
            )
            == 1e-8
        )

        test_build_two_salt(diafiltration_two_salt_lithium_aluminum)

    @pytest.mark.component
    def test_diagnostics_lithium_aluminum(
        self, diafiltration_two_salt_lithium_aluminum
    ):
        test_diagnostics(diafiltration_two_salt_lithium_aluminum)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_lithium_aluminum(self, diafiltration_two_salt_lithium_aluminum):
        test_solve(diafiltration_two_salt_lithium_aluminum)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                5.6461,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                171.18,
            ],
            "aluminum_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                36.668,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                281.18,
            ],
            "permeate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                9.0248,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                194.44,
            ],
            "aluminum_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                13.760,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_two_salt_lithium_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                235.72,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_lithium_aluminum(
        self, diafiltration_two_salt_lithium_aluminum
    ):
        test_numerical_issues(diafiltration_two_salt_lithium_aluminum)


################################################################################
# Test a two-salt model: cobalt chloride + aluminum chloride


@pytest.fixture(scope="module")
def diafiltration_two_salt_cobalt_aluminum():
    """
    Build two-salt diafiltration unit model for cobalt chloride + aluminum chloride.
    """
    cation_list = ["cobalt", "aluminum"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"cobalt": 288, "aluminum": 20, "chloride": 636},
            "diafiltrate": {"cobalt": 3, "aluminum": 3, "chloride": 15},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_cobalt_aluminum(diafiltration_two_salt_cobalt_aluminum):
    test_config_two_salt(diafiltration_two_salt_cobalt_aluminum)


class TestDiafiltrationTwoSaltCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_cobalt_aluminum(self, diafiltration_two_salt_cobalt_aluminum):
        assert isinstance(
            diafiltration_two_salt_cobalt_aluminum.fs.unit.numerical_zero_tolerance,
            Param,
        )
        assert (
            value(
                diafiltration_two_salt_cobalt_aluminum.fs.unit.numerical_zero_tolerance
            )
            == 1e-10
        )

        test_build_two_salt(diafiltration_two_salt_cobalt_aluminum)

    @pytest.mark.component
    def test_diagnostics_cobalt_aluminum(self, diafiltration_two_salt_cobalt_aluminum):
        test_diagnostics(diafiltration_two_salt_cobalt_aluminum)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_cobalt_aluminum(self, diafiltration_two_salt_cobalt_aluminum):
        test_solve(diafiltration_two_salt_cobalt_aluminum)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                8.2716,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                232.23,
            ],
            "aluminum_retentate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                18.276,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                519.28,
            ],
            "permeate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                7.8847,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                218.01,
            ],
            "aluminum_permeate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                14.954,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_two_salt_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                480.88,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_cobalt_aluminum(
        self, diafiltration_two_salt_cobalt_aluminum
    ):
        test_numerical_issues(diafiltration_two_salt_cobalt_aluminum)


################################################################################
# Test functions for three-salt model
@pytest.fixture(scope="module")
def sample_three_salt_model():
    cation_list = ["lithium", "cobalt", "aluminum"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "cobalt": 288, "aluminum": 20, "chloride": 881},
            "diafiltrate": {"lithium": 14, "cobalt": 3, "aluminum": 3, "chloride": 29},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    return m


@pytest.mark.unit
def test_config_three_salt(sample_three_salt_model):
    assert len(sample_three_salt_model.fs.unit.config) == 10

    assert not sample_three_salt_model.fs.unit.config.dynamic
    assert not sample_three_salt_model.fs.unit.config.has_holdup

    assert (
        sample_three_salt_model.fs.unit.config.property_package
        is sample_three_salt_model.fs.properties
    )
    assert len(sample_three_salt_model.fs.unit.config.cation_list) == 3
    assert len(sample_three_salt_model.fs.unit.config.anion_list) == 1
    assert sample_three_salt_model.fs.unit.config.NFE_module_length == 10
    assert sample_three_salt_model.fs.unit.config.NFE_membrane_thickness == 5


@pytest.mark.build
@pytest.mark.unit
def test_build_three_salt(sample_three_salt_model):
    # parameters
    assert isinstance(sample_three_salt_model.fs.unit.total_membrane_thickness, Param)
    assert value(sample_three_salt_model.fs.unit.total_membrane_thickness) == 1e-7

    assert isinstance(sample_three_salt_model.fs.unit.membrane_fixed_charge, Param)
    assert value(sample_three_salt_model.fs.unit.membrane_fixed_charge) == -44

    assert isinstance(sample_three_salt_model.fs.unit.membrane_permeability, Param)
    assert value(sample_three_salt_model.fs.unit.membrane_permeability) == 0.01

    assert isinstance(sample_three_salt_model.fs.unit.temperature, Param)
    assert value(sample_three_salt_model.fs.unit.temperature) == 298

    # sets
    assert isinstance(
        sample_three_salt_model.fs.unit.dimensionless_module_length, ContinuousSet
    )
    assert len(sample_three_salt_model.fs.unit.dimensionless_module_length) == 11

    assert isinstance(
        sample_three_salt_model.fs.unit.dimensionless_membrane_thickness,
        ContinuousSet,
    )
    assert len(sample_three_salt_model.fs.unit.dimensionless_membrane_thickness) == 6

    assert isinstance(sample_three_salt_model.fs.unit.time, Set)
    assert len(sample_three_salt_model.fs.unit.time) == 1

    assert isinstance(sample_three_salt_model.fs.unit.solutes, Set)
    assert len(sample_three_salt_model.fs.unit.solutes) == 4

    assert isinstance(sample_three_salt_model.fs.unit.cations, Set)
    assert len(sample_three_salt_model.fs.unit.cations) == 3

    # variables
    assert isinstance(sample_three_salt_model.fs.unit.total_module_length, Var)
    assert len(sample_three_salt_model.fs.unit.total_module_length) == 1

    assert isinstance(sample_three_salt_model.fs.unit.total_membrane_length, Var)
    assert len(sample_three_salt_model.fs.unit.total_membrane_length) == 1

    assert isinstance(sample_three_salt_model.fs.unit.applied_pressure, Var)
    assert len(sample_three_salt_model.fs.unit.applied_pressure) == 1

    assert isinstance(sample_three_salt_model.fs.unit.feed_flow_volume, Var)
    assert len(sample_three_salt_model.fs.unit.feed_flow_volume) == 1

    assert isinstance(sample_three_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.feed_conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.diafiltrate_flow_volume, Var)
    assert len(sample_three_salt_model.fs.unit.diafiltrate_flow_volume) == 1

    assert isinstance(sample_three_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.membrane_D_tilde, Var)
    assert len(sample_three_salt_model.fs.unit.membrane_D_tilde) == 66

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear,
        Var,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear
        )
        == 594
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_convection_coefficient_bilinear, Var
    )
    assert (
        len(sample_three_salt_model.fs.unit.membrane_convection_coefficient_bilinear)
        == 198
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient, Var
    )
    assert (
        len(sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient) == 594
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_convection_coefficient, Var
    )
    assert len(sample_three_salt_model.fs.unit.membrane_convection_coefficient) == 198

    assert isinstance(sample_three_salt_model.fs.unit.volume_flux_water, Var)
    assert len(sample_three_salt_model.fs.unit.volume_flux_water) == 11

    assert isinstance(sample_three_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_three_salt_model.fs.unit.molar_ion_flux) == 44

    assert isinstance(sample_three_salt_model.fs.unit.retentate_flow_volume, Var)
    assert len(sample_three_salt_model.fs.unit.retentate_flow_volume) == 11

    assert isinstance(sample_three_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.retentate_conc_mol_comp) == 44

    assert isinstance(sample_three_salt_model.fs.unit.permeate_flow_volume, Var)
    assert len(sample_three_salt_model.fs.unit.permeate_flow_volume) == 11

    assert isinstance(sample_three_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.permeate_conc_mol_comp) == 44

    assert isinstance(sample_three_salt_model.fs.unit.osmotic_pressure, Var)
    assert len(sample_three_salt_model.fs.unit.osmotic_pressure) == 11

    assert isinstance(sample_three_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.membrane_conc_mol_comp) == 264

    assert isinstance(
        sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 44

    assert isinstance(
        sample_three_salt_model.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
    )
    assert len(sample_three_salt_model.fs.unit.d_retentate_flow_volume_dx) == 11

    assert isinstance(
        sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 264

    # constraints
    assert isinstance(sample_three_salt_model.fs.unit.overall_mol_balance, Constraint)
    assert len(sample_three_salt_model.fs.unit.overall_mol_balance) == 10

    assert isinstance(sample_three_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_three_salt_model.fs.unit.cation_mol_balance) == 30

    assert isinstance(
        sample_three_salt_model.fs.unit.overall_bulk_flux_equation, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.overall_bulk_flux_equation) == 10

    assert isinstance(
        sample_three_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.cation_bulk_flux_equation) == 30

    assert isinstance(sample_three_salt_model.fs.unit.lumped_water_flux, Constraint)
    assert len(sample_three_salt_model.fs.unit.lumped_water_flux) == 10

    assert isinstance(
        sample_three_salt_model.fs.unit.chloride_flux_membrane, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.chloride_flux_membrane) == 10

    assert isinstance(
        sample_three_salt_model.fs.unit.osmotic_pressure_calculation, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.osmotic_pressure_calculation) == 10

    assert isinstance(
        sample_three_salt_model.fs.unit.electroneutrality_retentate, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.electroneutrality_retentate) == 11

    assert isinstance(
        sample_three_salt_model.fs.unit.electroneutrality_permeate, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.electroneutrality_permeate) == 10

    assert isinstance(
        sample_three_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface
        )
        == 30
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface
        )
        == 30
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_D_tilde_calculation, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.membrane_D_tilde_calculation) == 60

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation
        )
        == 540
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.membrane_convection_coefficient_bilinear_calculation
        )
        == 180
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_calculation
        )
        == 540
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_convection_coefficient_calculation,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.membrane_convection_coefficient_calculation)
        == 180
    )

    assert isinstance(sample_three_salt_model.fs.unit.cation_flux_membrane, Constraint)
    assert len(sample_three_salt_model.fs.unit.cation_flux_membrane) == 180

    assert isinstance(
        sample_three_salt_model.fs.unit.electroneutrality_membrane, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.electroneutrality_membrane) == 60

    assert isinstance(
        sample_three_salt_model.fs.unit.retentate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.retentate_flow_volume_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.retentate_conc_mol_comp_boundary_condition)
        == 3
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.membrane_conc_mol_comp_boundary_condition)
        == 18
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.permeate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.permeate_flow_volume_boundary_condition)
        == 1
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 4
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.d_retentate_flow_volume_dx_boundary_condition
        )
        == 1
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition,
        Constraint,
    )
    assert (
        len(
            sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition
        )
        == 3
    )
    assert isinstance(
        sample_three_salt_model.fs.unit.volume_flux_water_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.volume_flux_water_boundary_condition) == 1
    )

    assert isinstance(
        sample_three_salt_model.fs.unit.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(sample_three_salt_model.fs.unit.molar_ion_flux_boundary_condition) == 4

    for t in sample_three_salt_model.fs.unit.time:
        for x in sample_three_salt_model.fs.unit.dimensionless_module_length:
            assert sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx[
                t, x, "chloride"
            ].fixed
            if x != 0:
                assert not sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, "chloride"
                ].active

            for z in sample_three_salt_model.fs.unit.dimensionless_membrane_thickness:
                assert sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz[
                    t, x, z, "chloride"
                ].fixed
                if z != 0:
                    assert not sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, "chloride"
                    ].active

    # scaling factors
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.volume_flux_water
        ]
        == 1e2
    )
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.membrane_D_tilde
        ]
        == 1e-1
    )
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient_bilinear
        ]
        == 1e-2
    )
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.membrane_convection_coefficient_bilinear
        ]
        == 1e-1
    )
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.membrane_cross_diffusion_coefficient
        ]
        == 1e1
    )
    assert (
        sample_three_salt_model.fs.unit.scaling_factor[
            sample_three_salt_model.fs.unit.membrane_convection_coefficient
        ]
        == 1e1
    )
    for t in sample_three_salt_model.fs.unit.time:
        for x in sample_three_salt_model.fs.unit.dimensionless_module_length:
            if x != 0:
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[2]
                        ]
                    ]
                    == 1e-5
                )
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[1]
                        ]
                    ]
                    == 1e-5
                )
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[0]
                        ]
                    ]
                    == 1e-3
                )
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[2]
                        ]
                    ]
                    == 1e-5
                )
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[1]
                        ]
                    ]
                    == 1e-5
                )
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface[
                            t, x, sample_three_salt_model.fs.unit.config.cation_list[0]
                        ]
                    ]
                    == 1e-3
                )

    # ports
    assert isinstance(sample_three_salt_model.fs.unit.feed_inlet, Port)
    assert len(sample_three_salt_model.fs.unit.feed_inlet.flow_vol) == 1
    assert len(sample_three_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.diafiltrate_inlet, Port)
    assert len(sample_three_salt_model.fs.unit.diafiltrate_inlet.flow_vol) == 1
    assert len(sample_three_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.retentate_outlet, Port)
    assert len(sample_three_salt_model.fs.unit.retentate_outlet.flow_vol) == 1
    assert len(sample_three_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.permeate_outlet, Port)
    assert len(sample_three_salt_model.fs.unit.permeate_outlet.flow_vol) == 1
    assert len(sample_three_salt_model.fs.unit.permeate_outlet.conc_mol_comp) == 4


################################################################################
# Test a three-salt model: lithium chloride + cobalt chloride + aluminum chloride


@pytest.fixture(scope="module")
def diafiltration_three_salt_lithium_cobalt_aluminum():
    """
    Build two-salt diafiltration unit model for lithium chloride + cobalt chloride + aluminum chloride.
    """
    cation_list = ["lithium", "cobalt", "aluminum"]
    anion_list = ["chloride"]

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
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"lithium": 245, "cobalt": 288, "aluminum": 20, "chloride": 881},
            "diafiltrate": {"lithium": 14, "cobalt": 3, "aluminum": 3, "chloride": 29},
        },
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -44

    assert degrees_of_freedom(m.fs.unit) == 11

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.feed_flow_volume.fix()
    m.fs.unit.feed_conc_mol_comp.fix()
    m.fs.unit.diafiltrate_flow_volume.fix()
    m.fs.unit.diafiltrate_conc_mol_comp.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config_lithium_cobalt_aluminum(
    diafiltration_three_salt_lithium_cobalt_aluminum,
):
    test_config_three_salt(diafiltration_three_salt_lithium_cobalt_aluminum)


class TestDiafiltrationThreeSaltLithiumCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_lithium_cobalt_aluminum(
        self, diafiltration_three_salt_lithium_cobalt_aluminum
    ):
        assert isinstance(
            diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.numerical_zero_tolerance,
            Param,
        )
        assert (
            value(
                diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.numerical_zero_tolerance
            )
            == 1e-10
        )

        test_build_three_salt(diafiltration_three_salt_lithium_cobalt_aluminum)

    @pytest.mark.component
    def test_diagnostics_lithium_cobalt_aluminum(
        self, diafiltration_three_salt_lithium_cobalt_aluminum
    ):
        test_diagnostics(diafiltration_three_salt_lithium_cobalt_aluminum)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_lithium_cobalt_aluminum(
        self, diafiltration_three_salt_lithium_cobalt_aluminum
    ):
        test_solve(diafiltration_three_salt_lithium_cobalt_aluminum)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                10.2946,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                189.13,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                224.68,
            ],
            "aluminum_retentate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                21.673,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                703.51,
            ],
            "permeate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                5.7276,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "lithium"
                    ]
                ),
                194.91,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cobalt"
                    ]
                ),
                220.89,
            ],
            "aluminum_permeate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "aluminum"
                    ]
                ),
                8.3159,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_three_salt_lithium_cobalt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "chloride"
                    ]
                ),
                661.63,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_lithium_cobalt_aluminum(
        self, diafiltration_three_salt_lithium_cobalt_aluminum
    ):
        test_numerical_issues(diafiltration_three_salt_lithium_cobalt_aluminum)


################################################################################
