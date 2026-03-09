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
        NFE_module_length=10,
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
    assert len(sample_single_salt_model.fs.unit.config) == 8

    assert not sample_single_salt_model.fs.unit.config.dynamic
    assert not sample_single_salt_model.fs.unit.config.has_holdup

    assert (
        sample_single_salt_model.fs.unit.config.property_package
        is sample_single_salt_model.fs.properties
    )

    assert len(sample_single_salt_model.fs.unit.config.anion_list) == 1
    assert sample_single_salt_model.fs.unit.config.NFE_module_length == 10
    assert sample_single_salt_model.fs.unit.config.NFE_membrane_thickness == 5


@pytest.mark.unit
def test_config_single_salt(sample_single_salt_model):
    assert len(sample_single_salt_model.fs.unit.config.cation_list) == 1


@pytest.mark.build
@pytest.mark.unit
def test_build(sample_single_salt_model):
    # parameters
    assert isinstance(sample_single_salt_model.fs.unit.numerical_zero_tolerance, Param)
    assert value(sample_single_salt_model.fs.unit.numerical_zero_tolerance) == 1e-10

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

    # variables
    assert isinstance(sample_single_salt_model.fs.unit.total_module_length, Var)
    assert len(sample_single_salt_model.fs.unit.total_module_length) == 1

    assert isinstance(sample_single_salt_model.fs.unit.total_membrane_length, Var)
    assert len(sample_single_salt_model.fs.unit.total_membrane_length) == 1

    assert isinstance(sample_single_salt_model.fs.unit.applied_pressure, Var)
    assert len(sample_single_salt_model.fs.unit.applied_pressure) == 1

    assert isinstance(sample_single_salt_model.fs.unit.feed_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.feed_flow_volume) == 1

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_flow_volume) == 1

    assert isinstance(sample_single_salt_model.fs.unit.membrane_D_tilde, Var)
    assert len(sample_single_salt_model.fs.unit.membrane_D_tilde) == 66

    assert isinstance(sample_single_salt_model.fs.unit.volume_flux_water, Var)
    assert len(sample_single_salt_model.fs.unit.volume_flux_water) == 11

    assert isinstance(sample_single_salt_model.fs.unit.retentate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.retentate_flow_volume) == 11

    assert isinstance(sample_single_salt_model.fs.unit.permeate_flow_volume, Var)
    assert len(sample_single_salt_model.fs.unit.permeate_flow_volume) == 11

    assert isinstance(sample_single_salt_model.fs.unit.osmotic_pressure, Var)
    assert len(sample_single_salt_model.fs.unit.osmotic_pressure) == 11

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
    )
    assert len(sample_single_salt_model.fs.unit.d_retentate_flow_volume_dx) == 11

    # constraints
    assert isinstance(sample_single_salt_model.fs.unit.overall_mol_balance, Constraint)
    assert len(sample_single_salt_model.fs.unit.overall_mol_balance) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.overall_bulk_flux_equation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.overall_bulk_flux_equation) == 10

    assert isinstance(sample_single_salt_model.fs.unit.lumped_water_flux, Constraint)
    assert len(sample_single_salt_model.fs.unit.lumped_water_flux) == 10

    assert isinstance(sample_single_salt_model.fs.unit.anion_flux_membrane, Constraint)
    assert len(sample_single_salt_model.fs.unit.anion_flux_membrane) == 10

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
        sample_single_salt_model.fs.unit.membrane_D_tilde_calculation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.membrane_D_tilde_calculation) == 60

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
        sample_single_salt_model.fs.unit.permeate_flow_volume_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.permeate_flow_volume_boundary_condition)
        == 1
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
        sample_single_salt_model.fs.unit.volume_flux_water_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.volume_flux_water_boundary_condition) == 1
    )

    for t in sample_single_salt_model.fs.unit.time:
        for x in sample_single_salt_model.fs.unit.dimensionless_module_length:
            assert sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx[
                t, x, sample_single_salt_model.fs.unit.config.anion_list[0]
            ].fixed
            if x != 0:
                assert not sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    t, x, sample_single_salt_model.fs.unit.config.anion_list[0]
                ].active

            for z in sample_single_salt_model.fs.unit.dimensionless_membrane_thickness:
                assert sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz[
                    t, x, z, sample_single_salt_model.fs.unit.config.anion_list[0]
                ].fixed
                if z != 0:
                    assert not sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz_disc_eq[
                        t, x, z, sample_single_salt_model.fs.unit.config.anion_list[0]
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

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_inlet, Port)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_inlet.flow_vol) == 1

    assert isinstance(sample_single_salt_model.fs.unit.retentate_outlet, Port)
    assert len(sample_single_salt_model.fs.unit.retentate_outlet.flow_vol) == 1

    assert isinstance(sample_single_salt_model.fs.unit.permeate_outlet, Port)
    assert len(sample_single_salt_model.fs.unit.permeate_outlet.flow_vol) == 1


@pytest.mark.build
@pytest.mark.unit
def test_build_single_salt(sample_single_salt_model):
    # sets
    assert isinstance(sample_single_salt_model.fs.unit.solutes, Set)
    assert len(sample_single_salt_model.fs.unit.solutes) == 2

    assert isinstance(sample_single_salt_model.fs.unit.cations, Set)
    assert len(sample_single_salt_model.fs.unit.cations) == 1

    # variables
    assert isinstance(sample_single_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.feed_conc_mol_comp) == 2

    assert isinstance(sample_single_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 2

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

    assert isinstance(sample_single_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_single_salt_model.fs.unit.molar_ion_flux) == 22

    assert isinstance(sample_single_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.retentate_conc_mol_comp) == 22

    assert isinstance(sample_single_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.permeate_conc_mol_comp) == 22

    assert isinstance(sample_single_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_single_salt_model.fs.unit.membrane_conc_mol_comp) == 132

    assert isinstance(
        sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_single_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 22

    assert isinstance(
        sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_single_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 132

    # constraints
    assert isinstance(sample_single_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_single_salt_model.fs.unit.cation_mol_balance) == 10

    assert isinstance(
        sample_single_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_single_salt_model.fs.unit.cation_bulk_flux_equation) == 10

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
        sample_single_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_single_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 2
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
        sample_single_salt_model.fs.unit.molar_ion_flux_boundary_condition,
        Constraint,
    )
    assert len(sample_single_salt_model.fs.unit.molar_ion_flux_boundary_condition) == 2

    # ports
    assert len(sample_single_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 2
    assert len(sample_single_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 2
    assert len(sample_single_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 2
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
# Test single-salt model: LiCl


@pytest.fixture(scope="module")
def diafiltration_single_salt_Li():
    """
    Build single-salt diafiltration unit model for LiCl.
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
def test_config_Li(diafiltration_single_salt_Li):
    test_config(diafiltration_single_salt_Li)
    test_config_single_salt(diafiltration_single_salt_Li)


class TestDiafiltrationSingleSaltLithium(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li(self, diafiltration_single_salt_Li):
        test_build(diafiltration_single_salt_Li)
        test_build_single_salt(diafiltration_single_salt_Li)

    @pytest.mark.component
    def test_diagnostics_Li(self, diafiltration_single_salt_Li):
        test_diagnostics(diafiltration_single_salt_Li)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li(self, diafiltration_single_salt_Li):
        test_solve(diafiltration_single_salt_Li)

        test_dict = {
            "retentate_final": [
                value(diafiltration_single_salt_Li.fs.unit.retentate_flow_volume[0, 1]),
                8.3856,
            ],
            "Li_retentate_final": [
                value(
                    diafiltration_single_salt_Li.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                194.52,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_single_salt_Li.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                194.52,
            ],
            "permeate_final": [
                value(diafiltration_single_salt_Li.fs.unit.permeate_flow_volume[0, 1]),
                7.8637,
            ],
            "Li_permeate_final": [
                value(
                    diafiltration_single_salt_Li.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                190.38,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_single_salt_Li.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                190.38,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Li(self, diafiltration_single_salt_Li):
        test_numerical_issues(diafiltration_single_salt_Li)


################################################################################
# Test single-salt model: CoCl2


@pytest.fixture(scope="module")
def diafiltration_single_salt_Co():
    """
    Build single-salt diafiltration unit model for CoCl2.
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
def test_config_Co(diafiltration_single_salt_Co):
    test_config(diafiltration_single_salt_Co)
    test_config_single_salt(diafiltration_single_salt_Co)


class TestDiafiltrationSingleSaltCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co(self, diafiltration_single_salt_Co):
        test_build(diafiltration_single_salt_Co)
        test_build_single_salt(diafiltration_single_salt_Co)

    @pytest.mark.component
    def test_diagnostics_Co(self, diafiltration_single_salt_Co):
        test_diagnostics(diafiltration_single_salt_Co)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co(self, diafiltration_single_salt_Co):
        test_solve(diafiltration_single_salt_Co)

        test_dict = {
            "retentate_final": [
                value(diafiltration_single_salt_Co.fs.unit.retentate_flow_volume[0, 1]),
                10.468,
            ],
            "Co_retentate_final": [
                value(
                    diafiltration_single_salt_Co.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                227.59,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_single_salt_Co.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                455.17,
            ],
            "permeate_final": [
                value(diafiltration_single_salt_Co.fs.unit.permeate_flow_volume[0, 1]),
                5.7687,
            ],
            "Co_permeate_final": [
                value(
                    diafiltration_single_salt_Co.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                215.62,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_single_salt_Co.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                431.24,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Co(self, diafiltration_single_salt_Co):
        test_numerical_issues(diafiltration_single_salt_Co)


################################################################################
# Test single-salt model: AlCl3


@pytest.fixture(scope="module")
def diafiltration_single_salt_Al():
    """
    Build single-salt diafiltration unit model for AlCl3.
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
def test_config_Al(diafiltration_single_salt_Al):
    test_config(diafiltration_single_salt_Al)
    test_config_single_salt(diafiltration_single_salt_Al)


class TestDiafiltrationSingleSaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Al(self, diafiltration_single_salt_Al):
        test_build(diafiltration_single_salt_Al)
        test_build_single_salt(diafiltration_single_salt_Al)

    @pytest.mark.component
    def test_diagnostics_Al(self, diafiltration_single_salt_Al):
        test_diagnostics(diafiltration_single_salt_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Al(self, diafiltration_single_salt_Al):
        test_solve(diafiltration_single_salt_Al)

        test_dict = {
            "retentate_final": [
                value(diafiltration_single_salt_Al.fs.unit.retentate_flow_volume[0, 1]),
                9.4910,
            ],
            "Al_retentate_final": [
                value(
                    diafiltration_single_salt_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                18.044,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_single_salt_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                54.131,
            ],
            "permeate_final": [
                value(diafiltration_single_salt_Al.fs.unit.permeate_flow_volume[0, 1]),
                6.7079,
            ],
            "Al_permeate_final": [
                value(
                    diafiltration_single_salt_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                14.372,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_single_salt_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                43.115,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Al(self, diafiltration_single_salt_Al):
        test_numerical_issues(diafiltration_single_salt_Al)


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
        NFE_module_length=10,
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
    # sets
    assert isinstance(sample_two_salt_model.fs.unit.solutes, Set)
    assert len(sample_two_salt_model.fs.unit.solutes) == 3

    assert isinstance(sample_two_salt_model.fs.unit.cations, Set)
    assert len(sample_two_salt_model.fs.unit.cations) == 2

    # variables
    assert isinstance(sample_two_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.feed_conc_mol_comp) == 3

    assert isinstance(sample_two_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 3

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

    assert isinstance(sample_two_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_two_salt_model.fs.unit.molar_ion_flux) == 33

    assert isinstance(sample_two_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.retentate_conc_mol_comp) == 33

    assert isinstance(sample_two_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.permeate_conc_mol_comp) == 33

    assert isinstance(sample_two_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_two_salt_model.fs.unit.membrane_conc_mol_comp) == 198

    assert isinstance(
        sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_two_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 33

    assert isinstance(
        sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_two_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 198

    # constraints
    assert isinstance(sample_two_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_two_salt_model.fs.unit.cation_mol_balance) == 20

    assert isinstance(
        sample_two_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_two_salt_model.fs.unit.cation_bulk_flux_equation) == 20

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
        sample_two_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_two_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 3
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

    # scaling factors
    for t in sample_two_salt_model.fs.unit.time:
        for x in sample_two_salt_model.fs.unit.dimensionless_module_length:
            if x != 0:
                assert (
                    sample_two_salt_model.fs.unit.scaling_factor[
                        sample_two_salt_model.fs.unit.lumped_water_flux[t, x]
                    ]
                    == 1e3
                )

    # ports
    assert len(sample_two_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 3
    assert len(sample_two_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 3
    assert len(sample_two_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 3
    assert len(sample_two_salt_model.fs.unit.permeate_outlet.conc_mol_comp) == 3


################################################################################
# Test a two-salt model: LiCl + CoCl2


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Co():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2.
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
        NFE_module_length=10,
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
    test_config_two_salt(diafiltration_two_salt_Li_Co)


class TestDiafiltrationTwoSaltLithiumCobalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_build(diafiltration_two_salt_Li_Co)
        test_build_two_salt(diafiltration_two_salt_Li_Co)

    @pytest.mark.component
    def test_diagnostics_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_diagnostics(diafiltration_two_salt_Li_Co)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_solve(diafiltration_two_salt_Li_Co)

        test_dict = {
            "retentate_final": [
                value(diafiltration_two_salt_Li_Co.fs.unit.retentate_flow_volume[0, 1]),
                6.0465,
            ],
            "Li_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                188.88,
            ],
            "Co_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                246.67,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                682.22,
            ],
            "permeate_final": [
                value(diafiltration_two_salt_Li_Co.fs.unit.permeate_flow_volume[0, 1]),
                10.033,
            ],
            "Li_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                191.55,
            ],
            "Co_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                222.76,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Co.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                637.06,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Li_Co(self, diafiltration_two_salt_Li_Co):
        test_numerical_issues(diafiltration_two_salt_Li_Co)


################################################################################
# Test a two-salt model: LiCl + AlCl3


@pytest.fixture(scope="module")
def diafiltration_two_salt_Li_Al():
    """
    Build two-salt diafiltration unit model for LiCl + AlCl3.
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
        NFE_module_length=10,
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
    test_config_two_salt(diafiltration_two_salt_Li_Al)


class TestDiafiltrationTwoSaltLithiumAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_build(diafiltration_two_salt_Li_Al)
        test_build_two_salt(diafiltration_two_salt_Li_Al)

    @pytest.mark.component
    def test_diagnostics_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_diagnostics(diafiltration_two_salt_Li_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_solve(diafiltration_two_salt_Li_Al)

        test_dict = {
            "retentate_final": [
                value(diafiltration_two_salt_Li_Al.fs.unit.retentate_flow_volume[0, 1]),
                5.6461,
            ],
            "Li_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                171.18,
            ],
            "Al_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                36.668,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                281.18,
            ],
            "permeate_final": [
                value(diafiltration_two_salt_Li_Al.fs.unit.permeate_flow_volume[0, 1]),
                9.0248,
            ],
            "Li_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                194.44,
            ],
            "Al_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                13.760,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_two_salt_Li_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                235.72,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Li_Al(self, diafiltration_two_salt_Li_Al):
        test_numerical_issues(diafiltration_two_salt_Li_Al)


################################################################################
# Test a two-salt model: CoCl2 + AlCl3


@pytest.fixture(scope="module")
def diafiltration_two_salt_Co_Al():
    """
    Build two-salt diafiltration unit model for CoCl2 + AlCl3.
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
        NFE_module_length=10,
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
    test_config_two_salt(diafiltration_two_salt_Co_Al)


class TestDiafiltrationTwoSaltCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_build(diafiltration_two_salt_Co_Al)
        test_build_two_salt(diafiltration_two_salt_Co_Al)

    @pytest.mark.component
    def test_diagnostics_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_diagnostics(diafiltration_two_salt_Co_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_solve(diafiltration_two_salt_Co_Al)

        test_dict = {
            "retentate_final": [
                value(diafiltration_two_salt_Co_Al.fs.unit.retentate_flow_volume[0, 1]),
                8.2716,
            ],
            "Co_retentate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                232.23,
            ],
            "Al_retentate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                18.276,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                519.28,
            ],
            "permeate_final": [
                value(diafiltration_two_salt_Co_Al.fs.unit.permeate_flow_volume[0, 1]),
                7.8847,
            ],
            "Co_permeate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                218.01,
            ],
            "Al_permeate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                14.954,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_two_salt_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                480.88,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Co_Al(self, diafiltration_two_salt_Co_Al):
        test_numerical_issues(diafiltration_two_salt_Co_Al)


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
        NFE_module_length=10,
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
    # sets
    assert isinstance(sample_three_salt_model.fs.unit.solutes, Set)
    assert len(sample_three_salt_model.fs.unit.solutes) == 4

    assert isinstance(sample_three_salt_model.fs.unit.cations, Set)
    assert len(sample_three_salt_model.fs.unit.cations) == 3

    # variables
    assert isinstance(sample_three_salt_model.fs.unit.feed_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.feed_conc_mol_comp) == 4

    assert isinstance(sample_three_salt_model.fs.unit.diafiltrate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.diafiltrate_conc_mol_comp) == 4

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

    assert isinstance(sample_three_salt_model.fs.unit.molar_ion_flux, Var)
    assert len(sample_three_salt_model.fs.unit.molar_ion_flux) == 44

    assert isinstance(sample_three_salt_model.fs.unit.retentate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.retentate_conc_mol_comp) == 44

    assert isinstance(sample_three_salt_model.fs.unit.permeate_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.permeate_conc_mol_comp) == 44

    assert isinstance(sample_three_salt_model.fs.unit.membrane_conc_mol_comp, Var)
    assert len(sample_three_salt_model.fs.unit.membrane_conc_mol_comp) == 264

    assert isinstance(
        sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx,
        DerivativeVar,
    )
    assert len(sample_three_salt_model.fs.unit.d_retentate_conc_mol_comp_dx) == 44

    assert isinstance(
        sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz,
        DerivativeVar,
    )
    assert len(sample_three_salt_model.fs.unit.d_membrane_conc_mol_comp_dz) == 264

    # constraints
    assert isinstance(sample_three_salt_model.fs.unit.cation_mol_balance, Constraint)
    assert len(sample_three_salt_model.fs.unit.cation_mol_balance) == 30

    assert isinstance(
        sample_three_salt_model.fs.unit.cation_bulk_flux_equation, Constraint
    )
    assert len(sample_three_salt_model.fs.unit.cation_bulk_flux_equation) == 30

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
        sample_three_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition,
        Constraint,
    )
    assert (
        len(sample_three_salt_model.fs.unit.permeate_conc_mol_comp_boundary_condition)
        == 4
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

    # scaling factors
    for t in sample_three_salt_model.fs.unit.time:
        for x in sample_three_salt_model.fs.unit.dimensionless_module_length:
            if x != 0:
                assert (
                    sample_three_salt_model.fs.unit.scaling_factor[
                        sample_three_salt_model.fs.unit.lumped_water_flux[t, x]
                    ]
                    == 1e3
                )

    # ports
    assert len(sample_three_salt_model.fs.unit.feed_inlet.conc_mol_comp) == 4
    assert len(sample_three_salt_model.fs.unit.diafiltrate_inlet.conc_mol_comp) == 4
    assert len(sample_three_salt_model.fs.unit.retentate_outlet.conc_mol_comp) == 4
    assert len(sample_three_salt_model.fs.unit.permeate_outlet.conc_mol_comp) == 4


################################################################################
# Test a three-salt model: LiCl + CoCl2 + AlCl3


@pytest.fixture(scope="module")
def diafiltration_three_salt_Li_Co_Al():
    """
    Build two-salt diafiltration unit model for LiCl + CoCl2 + AlCl3.
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
        NFE_module_length=10,
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
def test_config_Li_Co_Al(
    diafiltration_three_salt_Li_Co_Al,
):
    test_config(diafiltration_three_salt_Li_Co_Al)
    test_config_three_salt(diafiltration_three_salt_Li_Co_Al)


class TestDiafiltrationThreeSaltLithiumCobaltAluminum(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_build(diafiltration_three_salt_Li_Co_Al)
        test_build_three_salt(diafiltration_three_salt_Li_Co_Al)

    @pytest.mark.component
    def test_diagnostics_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_diagnostics(diafiltration_three_salt_Li_Co_Al)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_solve(diafiltration_three_salt_Li_Co_Al)

        test_dict = {
            "retentate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.retentate_flow_volume[
                        0, 1
                    ]
                ),
                10.2946,
            ],
            "Li_retentate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                189.13,
            ],
            "Co_retentate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                224.68,
            ],
            "Al_retentate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                21.673,
            ],
            "Cl_retentate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.retentate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                703.51,
            ],
            "permeate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.permeate_flow_volume[0, 1]
                ),
                5.7276,
            ],
            "Li_permeate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Li"
                    ]
                ),
                194.91,
            ],
            "Co_permeate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Co"
                    ]
                ),
                220.89,
            ],
            "Al_permeate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Al"
                    ]
                ),
                8.3159,
            ],
            "Cl_permeate_final": [
                value(
                    diafiltration_three_salt_Li_Co_Al.fs.unit.permeate_conc_mol_comp[
                        0, 1, "Cl"
                    ]
                ),
                661.63,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-4) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues_Li_Co_Al(self, diafiltration_three_salt_Li_Co_Al):
        test_numerical_issues(diafiltration_three_salt_Li_Co_Al)


################################################################################
