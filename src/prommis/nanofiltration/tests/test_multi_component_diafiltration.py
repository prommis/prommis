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
    num_salts = 1
    salt_system = "lithium_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )
    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        num_salts=num_salts,
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
    assert len(sample_single_salt_model.fs.unit.config) == 7

    assert not sample_single_salt_model.fs.unit.config.dynamic
    assert not sample_single_salt_model.fs.unit.config.has_holdup

    assert (
        sample_single_salt_model.fs.unit.config.property_package
        is sample_single_salt_model.fs.properties
    )
    assert sample_single_salt_model.fs.unit.config.num_salts == 1
    assert sample_single_salt_model.fs.unit.config.NFE_module_length == 10
    assert sample_single_salt_model.fs.unit.config.NFE_membrane_thickness == 5


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt(sample_single_salt_model):
    # parameters
    assert isinstance(sample_single_salt_model.fs.unit.numerical_zero_tolerance, Param)
    assert value(sample_single_salt_model.fs.unit.numerical_zero_tolerance) == 1e-10

    assert isinstance(sample_single_salt_model.fs.unit.total_membrane_thickness, Param)
    assert value(sample_single_salt_model.fs.unit.total_membrane_thickness) == 1e-7

    assert isinstance(sample_single_salt_model.fs.unit.membrane_fixed_charge, Param)
    assert value(sample_single_salt_model.fs.unit.membrane_fixed_charge) == -140

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
                t, x, "anion"
            ].fixed
            if x != 0:
                assert not sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, "anion"
                ].active

            for z in sample_single_salt_model.fs.unit.dimensionless_membrane_thickness:
                assert sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz[
                    t, x, z, "anion"
                ].fixed
                if z != 0:
                    assert not sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, "anion"
                    ].active

    # scaling factors
    for t in sample_single_salt_model.fs.unit.time:
        for x in sample_single_salt_model.fs.unit.dimensionless_module_length:
            if x != 0:
                assert (
                    sample_single_salt_model.fs.unit.scaling_factor[
                        sample_single_salt_model.fs.unit.cation_equilibrium_retentate_membrane_interface[
                            t, x, "cation_1"
                        ]
                    ]
                    == 1e-3
                )
                assert (
                    sample_single_salt_model.fs.unit.scaling_factor[
                        sample_single_salt_model.fs.unit.cation_equilibrium_membrane_permeate_interface[
                            t, x, "cation_1"
                        ]
                    ]
                    == 1e-3
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
    num_salts = 1
    salt_system = "lithium_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        num_salts=num_salts,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -140

    assert degrees_of_freedom(m.fs.unit) == 7

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
def test_config_lithium(diafiltration_single_salt_lithium):
    test_config_single_salt(diafiltration_single_salt_lithium)


class TestDiafiltrationSingleSaltLithium(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_lithium(self, diafiltration_single_salt_lithium):
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
                1.5194,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                251.17,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.retentate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                251.17,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_flow_volume[0, 1]
                ),
                14.708,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                230.35,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_lithium.fs.unit.permeate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                230.35,
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
    num_salts = 1
    salt_system = "cobalt_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        num_salts=num_salts,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -140

    assert degrees_of_freedom(m.fs.unit) == 7

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
def test_config_cobalt(diafiltration_single_salt_cobalt):
    test_config_single_salt(diafiltration_single_salt_cobalt)


class TestDiafiltrationSingleSaltCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_cobalt(self, diafiltration_single_salt_cobalt):
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
                6.0411,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                223.67,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                447.34,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_flow_volume[0, 1]
                ),
                10.020,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                192.27,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_cobalt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                384.54,
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
    num_salts = 1
    salt_system = "aluminum_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        num_salts=num_salts,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -140

    assert degrees_of_freedom(m.fs.unit) == 7

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
def test_config_aluminum(diafiltration_single_salt_aluminum):
    test_config_single_salt(diafiltration_single_salt_aluminum)


class TestDiafiltrationSingleSaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_aluminum(self, diafiltration_single_salt_aluminum):
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
                10.303,
            ],
            "aluminum_retentate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                203.70,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.retentate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                611.09,
            ],
            "permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_flow_volume[
                        0, 1
                    ]
                ),
                5.8608,
            ],
            "aluminum_permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                177.76,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_single_salt_aluminum.fs.unit.permeate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                533.28,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_aluminum(self, diafiltration_single_salt_aluminum):
        test_numerical_issues(diafiltration_single_salt_aluminum)


################################################################################
# Test a two-salt model: lithium chloride + cobalt chloride


@pytest.fixture(scope="module")
def diafiltration_two_salt():
    """
    Build two-salt diafiltration unit model.
    """
    num_salts = 2
    salt_system = "lithium_cobalt_chloride"

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
        num_salts=num_salts,
        salt_system=salt_system,
    )

    m.fs.unit = MultiComponentDiafiltration(
        property_package=m.fs.properties,
        num_salts=num_salts,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -140

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
def test_config(diafiltration_two_salt):
    assert len(diafiltration_two_salt.fs.unit.config) == 7

    assert not diafiltration_two_salt.fs.unit.config.dynamic
    assert not diafiltration_two_salt.fs.unit.config.has_holdup

    assert (
        diafiltration_two_salt.fs.unit.config.property_package
        is diafiltration_two_salt.fs.properties
    )
    assert diafiltration_two_salt.fs.unit.config.num_salts == 2
    assert diafiltration_two_salt.fs.unit.config.NFE_module_length == 10
    assert diafiltration_two_salt.fs.unit.config.NFE_membrane_thickness == 5


class TestDiafiltrationTwoSalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, diafiltration_two_salt):
        # parameters
        assert isinstance(
            diafiltration_two_salt.fs.unit.numerical_zero_tolerance, Param
        )
        assert value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance) == 1e-10

        assert isinstance(
            diafiltration_two_salt.fs.unit.total_membrane_thickness, Param
        )
        assert value(diafiltration_two_salt.fs.unit.total_membrane_thickness) == 1e-7

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_fixed_charge, Param)
        assert value(diafiltration_two_salt.fs.unit.membrane_fixed_charge) == -140

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_permeability, Param)
        assert value(diafiltration_two_salt.fs.unit.membrane_permeability) == 0.01

        assert isinstance(diafiltration_two_salt.fs.unit.temperature, Param)
        assert value(diafiltration_two_salt.fs.unit.temperature) == 298

        # sets
        assert isinstance(
            diafiltration_two_salt.fs.unit.dimensionless_module_length, ContinuousSet
        )
        assert len(diafiltration_two_salt.fs.unit.dimensionless_module_length) == 11

        assert isinstance(
            diafiltration_two_salt.fs.unit.dimensionless_membrane_thickness,
            ContinuousSet,
        )
        assert len(diafiltration_two_salt.fs.unit.dimensionless_membrane_thickness) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.time, Set)
        assert len(diafiltration_two_salt.fs.unit.time) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.solutes, Set)
        assert len(diafiltration_two_salt.fs.unit.solutes) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.cations, Set)
        assert len(diafiltration_two_salt.fs.unit.cations) == 2

        # variables
        assert isinstance(diafiltration_two_salt.fs.unit.total_module_length, Var)
        assert len(diafiltration_two_salt.fs.unit.total_module_length) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.total_membrane_length, Var)
        assert len(diafiltration_two_salt.fs.unit.total_membrane_length) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.applied_pressure, Var)
        assert len(diafiltration_two_salt.fs.unit.applied_pressure) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.feed_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.feed_flow_volume) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.feed_conc_mol_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.feed_conc_mol_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_conc_mol_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_conc_mol_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_D_tilde, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_D_tilde) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_bilinear,
            Var,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_bilinear
            )
            == 264
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_convection_coefficient_bilinear, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_convection_coefficient_bilinear)
            == 132
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient)
            == 264
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_convection_coefficient, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_convection_coefficient) == 132
        )

        assert isinstance(diafiltration_two_salt.fs.unit.volume_flux_water, Var)
        assert len(diafiltration_two_salt.fs.unit.volume_flux_water) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.molar_ion_flux, Var)
        assert len(diafiltration_two_salt.fs.unit.molar_ion_flux) == 33

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.retentate_flow_volume) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_conc_mol_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.retentate_conc_mol_comp) == 33

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.permeate_flow_volume) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_conc_mol_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.permeate_conc_mol_comp) == 33

        assert isinstance(diafiltration_two_salt.fs.unit.osmotic_pressure, Var)
        assert len(diafiltration_two_salt.fs.unit.osmotic_pressure) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_conc_mol_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mol_comp) == 198

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx,
            DerivativeVar,
        )
        assert len(diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx) == 33

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
        )
        assert len(diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx) == 11

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mol_comp_dz,
            DerivativeVar,
        )
        assert len(diafiltration_two_salt.fs.unit.d_membrane_conc_mol_comp_dz) == 198

        # constraints
        assert isinstance(
            diafiltration_two_salt.fs.unit.overall_mol_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.overall_mol_balance) == 10

        assert isinstance(diafiltration_two_salt.fs.unit.cation_mol_balance, Constraint)
        assert len(diafiltration_two_salt.fs.unit.cation_mol_balance) == 20

        assert isinstance(
            diafiltration_two_salt.fs.unit.overall_bulk_flux_equation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.overall_bulk_flux_equation) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.cation_bulk_flux_equation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.cation_bulk_flux_equation) == 20

        assert isinstance(diafiltration_two_salt.fs.unit.lumped_water_flux, Constraint)
        assert len(diafiltration_two_salt.fs.unit.lumped_water_flux) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.chloride_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.chloride_flux_membrane) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.osmotic_pressure_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.osmotic_pressure_calculation) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_retentate, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_retentate) == 11

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_permeate, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_permeate) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.cation_equilibrium_retentate_membrane_interface,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.cation_equilibrium_retentate_membrane_interface
            )
            == 20
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.cation_equilibrium_membrane_permeate_interface,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.cation_equilibrium_membrane_permeate_interface
            )
            == 20
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_D_tilde_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.membrane_D_tilde_calculation) == 60

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_bilinear_calculation
            )
            == 240
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_convection_coefficient_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_convection_coefficient_bilinear_calculation
            )
            == 120
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_cross_diffusion_coefficient_calculation
            )
            == 240
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_convection_coefficient_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_convection_coefficient_calculation
            )
            == 120
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.cation_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.cation_flux_membrane) == 120

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_membrane) == 60

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_flow_volume_boundary_condition,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_flow_volume_boundary_condition)
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.retentate_conc_mol_comp_boundary_condition
            )
            == 2
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.membrane_conc_mol_comp_boundary_condition
            )
            == 12
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.permeate_flow_volume_boundary_condition,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.permeate_flow_volume_boundary_condition)
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.permeate_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.permeate_conc_mol_comp_boundary_condition
            )
            == 3
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx_boundary_condition,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx_boundary_condition
            )
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx_boundary_condition
            )
            == 2
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.volume_flux_water_boundary_condition,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.volume_flux_water_boundary_condition)
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.molar_ion_flux_boundary_condition,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.molar_ion_flux_boundary_condition) == 3
        )

        for t in diafiltration_two_salt.fs.unit.time:
            for x in diafiltration_two_salt.fs.unit.dimensionless_module_length:
                assert diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[
                    t, x, "anion"
                ].fixed
                if x != 0:
                    assert not diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                        t, x, "anion"
                    ].active

                for (
                    z
                ) in diafiltration_two_salt.fs.unit.dimensionless_membrane_thickness:
                    assert diafiltration_two_salt.fs.unit.d_membrane_conc_mol_comp_dz[
                        t, x, z, "anion"
                    ].fixed
                    if z != 0:
                        assert not diafiltration_two_salt.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                            t, x, z, "anion"
                        ].active

        # scaling factors
        for t in diafiltration_two_salt.fs.unit.time:
            for x in diafiltration_two_salt.fs.unit.dimensionless_module_length:
                if x != 0:
                    assert (
                        diafiltration_two_salt.fs.unit.scaling_factor[
                            diafiltration_two_salt.fs.unit.cation_equilibrium_retentate_membrane_interface[
                                t, x, "cation_2"
                            ]
                        ]
                        == 1e-5
                    )
                    assert (
                        diafiltration_two_salt.fs.unit.scaling_factor[
                            diafiltration_two_salt.fs.unit.cation_equilibrium_retentate_membrane_interface[
                                t, x, "cation_1"
                            ]
                        ]
                        == 1e-3
                    )
                    assert (
                        diafiltration_two_salt.fs.unit.scaling_factor[
                            diafiltration_two_salt.fs.unit.cation_equilibrium_membrane_permeate_interface[
                                t, x, "cation_2"
                            ]
                        ]
                        == 1e-5
                    )
                    assert (
                        diafiltration_two_salt.fs.unit.scaling_factor[
                            diafiltration_two_salt.fs.unit.cation_equilibrium_membrane_permeate_interface[
                                t, x, "cation_1"
                            ]
                        ]
                        == 1e-3
                    )

        # ports
        assert isinstance(diafiltration_two_salt.fs.unit.feed_inlet, Port)
        assert len(diafiltration_two_salt.fs.unit.feed_inlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.feed_inlet.conc_mol_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_inlet, Port)
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_inlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_inlet.conc_mol_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_outlet, Port)
        assert len(diafiltration_two_salt.fs.unit.retentate_outlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.retentate_outlet.conc_mol_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_outlet, Port)
        assert len(diafiltration_two_salt.fs.unit.permeate_outlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.permeate_outlet.conc_mol_comp) == 3

    @pytest.mark.component
    def test_diagnostics(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve(self, diafiltration_two_salt):
        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(diafiltration_two_salt, rename=False)

        solver = SolverFactory("ipopt")
        results = solver.solve(scaled_model, tee=True)

        scaling.propagate_solution(scaled_model, diafiltration_two_salt)

        assert_optimal_termination(results)

        test_dict = {
            "retentate_final": [
                value(diafiltration_two_salt.fs.unit.retentate_flow_volume[0, 1]),
                6.6914,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                197.90,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "cation_2"
                    ]
                ),
                241.11,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[
                        0, 1, "anion"
                    ]
                ),
                680.12,
            ],
            "permeate_final": [
                value(diafiltration_two_salt.fs.unit.permeate_flow_volume[0, 1]),
                9.4615,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cation_1"
                    ]
                ),
                191.35,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[
                        0, 1, "cation_2"
                    ]
                ),
                220.46,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 1, "anion"]
                ),
                632.27,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        dt.assert_no_numerical_warnings()


################################################################################
