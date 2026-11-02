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
    Expression,
    Param,
    Set,
    SolverFactory,
    TransformationFactory,
    Var,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import assert_solution_equivalent

import pytest

from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)
from prommis.nanofiltration.multi_component_diafiltration import (
    MultiComponentDiafiltration,
)

# TODO: test positive and neutral membrane cases (low priority)


################################################################################
# Define fixture for membrane model
@pytest.fixture(scope="module")
def diafiltration_model():
    """
    Build single-salt diafiltration unit model for LiCl.
    """

    def _diafiltration_model(
        cation_list,
        anion_list,
        inlet_flow_volume,
        inlet_concentration,
        non_Donnan_partition_dict,
        boundary_layer,
    ):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MultiComponentDiafiltrationSoluteParameter(
            cation_list=cation_list,
            anion_list=anion_list,
            non_Donnan_partition_dict=non_Donnan_partition_dict,
        )

        Dm_Cl = units.convert(
            m.fs.properties.membrane_diffusion_coefficient[anion_list[0]],
            to_units=units.um**2 / units.s,
        )
        Dm_over_l_value = 60  # um/s
        l_um = value(Dm_Cl) / Dm_over_l_value  # um
        l_m = l_um / 1e6  # m

        m.fs.unit = MultiComponentDiafiltration(
            property_package=m.fs.properties,
            cation_list=cation_list,
            anion_list=anion_list,
            include_boundary_layer=boundary_layer,
            total_membrane_thickness=l_m,
            NFE_module_length=10,
            NFE_boundary_layer_thickness=10,
            NFE_membrane_thickness=10,
        )

        assert degrees_of_freedom(m.fs.unit) == 5 + 2 * len(m.fs.unit.cations)

        m.fs.unit.total_module_length.fix()
        m.fs.unit.total_membrane_length.fix()

        m.fs.unit.feed_flow_volume.fix(inlet_flow_volume["feed"])
        m.fs.unit.diafiltrate_flow_volume.fix(inlet_flow_volume["diafiltrate"])

        for t in m.fs.unit.time:
            for j in m.fs.unit.solutes:
                m.fs.unit.feed_conc_mol_comp[t, j].fix(inlet_concentration["feed"][j])
                m.fs.unit.diafiltrate_conc_mol_comp[t, j].fix(
                    inlet_concentration["diafiltrate"][j]
                )

        feed_ionic_strength = value(m.fs.unit.feed_ionic_strength[0])

        if feed_ionic_strength < 201:
            m.fs.unit.applied_pressure.fix(5)
        elif feed_ionic_strength >= 201:
            m.fs.unit.applied_pressure.fix(15)

        initialized_model = m.fs.unit.default_initializer(
            multiplier_H_feed=1.2,  # increase can help solver performance
            multiplier_H_perm=1,
        )
        initialized_model.initialize(m.fs.unit)

        assert degrees_of_freedom(m.fs.unit) == 0

        return m

    return _diafiltration_model


################################################################################
# Build general test functions
@pytest.mark.unit
def test_config(diafiltration_model):
    def _test_config(m):
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties
        assert len(m.fs.unit.config.anion_list) == 1
        assert m.fs.unit.config.NFE_module_length == 10
        assert m.fs.unit.config.NFE_boundary_layer_thickness == 10
        assert m.fs.unit.config.NFE_membrane_thickness == 10

    model = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 200},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_config(model)
    assert model.fs.unit.config.include_boundary_layer == True
    assert len(model.fs.unit.config.cation_list) == 1

    model_no_boundary_layer = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 200},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_config(model_no_boundary_layer)
    assert model_no_boundary_layer.fs.unit.config.include_boundary_layer == False

    model_two_salt = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Cl": 600},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_config(model_two_salt)
    assert len(model_two_salt.fs.unit.config.cation_list) == 2

    model_three_salt = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Al": 20, "Cl": 660},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_config(model_three_salt)
    assert len(model_three_salt.fs.unit.config.cation_list) == 3


@pytest.mark.build
@pytest.mark.unit
def test_build(diafiltration_model):
    def _test_build(membrane):
        # parameters
        assert isinstance(membrane.numerical_zero_tolerance, Param)
        assert value(membrane.numerical_zero_tolerance) == 1e-10

        assert isinstance(membrane.total_membrane_thickness, Param)
        assert (
            value(membrane.total_membrane_thickness)
            == membrane.config.total_membrane_thickness
        )

        assert isinstance(membrane.membrane_fixed_charge, Param)
        assert value(membrane.membrane_fixed_charge) == -44

        assert isinstance(membrane.membrane_permeability, Param)
        assert value(membrane.membrane_permeability) == 11

        assert isinstance(membrane.temperature, Param)
        assert value(membrane.temperature) == 298

        # sets
        assert isinstance(membrane.dimensionless_module_length, ContinuousSet)
        assert len(membrane.dimensionless_module_length) == 11

        assert isinstance(
            membrane.dimensionless_membrane_thickness,
            ContinuousSet,
        )
        assert len(membrane.dimensionless_membrane_thickness) == 11

        assert isinstance(membrane.time, Set)
        assert len(membrane.time) == 1

        # dependent on number of solutes
        assert isinstance(membrane.solutes, Set)
        assert isinstance(membrane.cations, Set)

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
        assert len(membrane.membrane_D_tilde) == 121

        assert isinstance(membrane.volume_flux_water, Var)
        assert len(membrane.volume_flux_water) == 11

        assert isinstance(membrane.retentate_flow_volume, Var)
        assert len(membrane.retentate_flow_volume) == 11

        assert isinstance(membrane.osmotic_pressure, Var)
        assert len(membrane.osmotic_pressure) == 11

        assert isinstance(membrane.Donnan_potential_feed_side, Var)
        assert len(membrane.Donnan_potential_feed_side) == 11

        assert isinstance(membrane.Donnan_potential_permeate_side, Var)
        assert len(membrane.Donnan_potential_permeate_side) == 11

        assert isinstance(membrane.d_retentate_flow_volume_dx, DerivativeVar)
        assert len(membrane.d_retentate_flow_volume_dx) == 11

        # dependent on number of solutes
        assert isinstance(membrane.feed_conc_mol_comp, Var)
        assert isinstance(membrane.diafiltrate_conc_mol_comp, Var)
        assert isinstance(membrane.membrane_cross_diffusion_coefficient_bilinear, Var)
        assert isinstance(membrane.membrane_convection_coefficient_bilinear, Var)
        assert isinstance(membrane.membrane_cross_diffusion_coefficient, Var)
        assert isinstance(membrane.membrane_convection_coefficient, Var)
        assert isinstance(membrane.membrane_conc_mol_comp, Var)
        assert isinstance(membrane.molar_ion_flux, Var)
        assert isinstance(membrane.retentate_conc_mol_comp, Var)
        assert isinstance(membrane.permeate_local_conc_mol_comp, Var)
        assert isinstance(membrane.partitioning_term_bilinear_feed, Var)
        assert isinstance(membrane.partitioning_term_bilinear_permeate, Var)
        assert isinstance(membrane.d_retentate_conc_mol_comp_dx, DerivativeVar)
        assert isinstance(membrane.d_membrane_conc_mol_comp_dz, DerivativeVar)

        # constraints
        assert isinstance(membrane.overall_mol_balance, Constraint)
        assert len(membrane.overall_mol_balance) == 10

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
        assert len(membrane.electroneutrality_membrane) == 110

        assert isinstance(membrane.membrane_D_tilde_calculation, Constraint)
        assert len(membrane.membrane_D_tilde_calculation) == 110

        assert isinstance(
            membrane.feed_flow_volume_boundary_condition,
            Constraint,
        )
        assert len(membrane.feed_flow_volume_boundary_condition) == 1

        # dependent on number of solutes
        assert isinstance(membrane.cation_mol_balance, Constraint)
        assert isinstance(membrane.cation_bulk_flux_equation, Constraint)
        assert isinstance(
            membrane.membrane_cross_diffusion_coefficient_bilinear_calculation,
            Constraint,
        )
        assert isinstance(
            membrane.membrane_convection_coefficient_bilinear_calculation,
            Constraint,
        )
        assert isinstance(
            membrane.membrane_cross_diffusion_coefficient_calculation,
            Constraint,
        )
        assert isinstance(
            membrane.membrane_convection_coefficient_calculation,
            Constraint,
        )
        assert isinstance(membrane.cation_flux_membrane, Constraint)
        assert isinstance(
            membrane.partitioning_term_bilinear_feed_constraint, Constraint
        )
        assert isinstance(
            membrane.partitioning_term_bilinear_permeate_constraint, Constraint
        )
        assert isinstance(membrane.membrane_permeate_interface, Constraint)
        assert isinstance(
            membrane.feed_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert isinstance(
            membrane.membrane_conc_mol_comp_boundary_condition,
            Constraint,
        )

        for t in membrane.time:
            for x in membrane.dimensionless_module_length:
                assert membrane.d_retentate_conc_mol_comp_dx[
                    t, x, membrane.config.anion_list[0]
                ].fixed
                if x != 0:
                    assert not membrane.d_retentate_conc_mol_comp_dx_disc_eq[
                        t, x, membrane.config.anion_list[0]
                    ].active

        # scaling factors
        assert membrane.scaling_factor[membrane.volume_flux_water] == 1e2
        assert membrane.scaling_factor[membrane.lumped_water_flux] == 1e3
        assert membrane.scaling_factor[membrane.membrane_D_tilde] == 1e1
        assert (
            membrane.scaling_factor[
                membrane.membrane_cross_diffusion_coefficient_bilinear
            ]
            == 1e3
        )
        assert (
            membrane.scaling_factor[
                membrane.membrane_cross_diffusion_coefficient_bilinear_calculation
            ]
            == 1e3
        )
        assert (
            membrane.scaling_factor[membrane.membrane_convection_coefficient_bilinear]
            == 1e2
        )
        assert (
            membrane.scaling_factor[membrane.membrane_cross_diffusion_coefficient]
            == 1e5
        )
        assert (
            membrane.scaling_factor[
                membrane.membrane_cross_diffusion_coefficient_calculation
            ]
            == 1e5
        )
        assert membrane.scaling_factor[membrane.membrane_convection_coefficient] == 1e3

        # ports
        assert isinstance(membrane.feed_inlet, Port)
        assert len(membrane.feed_inlet.flow_vol) == 1

        assert isinstance(membrane.diafiltrate_inlet, Port)
        assert len(membrane.diafiltrate_inlet.flow_vol) == 1

        assert isinstance(membrane.retentate_outlet, Port)
        assert len(membrane.retentate_outlet.flow_vol) == 1

        assert isinstance(membrane.permeate_outlet, Port)
        assert len(membrane.permeate_outlet.flow_vol) == 1

        # expressions
        assert isinstance(membrane.permeate_outlet_flow, Expression)
        assert len(membrane.permeate_outlet_flow) == 1

        assert isinstance(membrane.feed_ionic_strength, Expression)
        assert len(membrane.feed_ionic_strength) == 1

        assert isinstance(membrane.total_feed_flow_volume, Expression)
        assert len(membrane.total_feed_flow_volume) == 1

        assert isinstance(membrane.membrane_electric_potential_gradient, Expression)
        assert len(membrane.membrane_electric_potential_gradient) == 110

        # dependent on number of solutes
        assert isinstance(membrane.permeate_outlet_conc_mol_comp, Expression)
        assert isinstance(membrane.total_feed_conc_mol_comp, Expression)
        assert isinstance(membrane.local_partition_coefficient_feed_side, Expression)
        assert isinstance(
            membrane.local_partition_coefficient_permeate_side, Expression
        )
        assert isinstance(membrane.local_observed_rejection_percent, Expression)
        assert isinstance(membrane.overall_observed_rejection_percent, Expression)
        assert isinstance(membrane.local_observed_sieving_coefficient, Expression)
        assert isinstance(membrane.overall_observed_sieving_coefficient, Expression)
        assert isinstance(membrane.membrane_convective_flux, Expression)
        assert isinstance(membrane.membrane_diffusive_flux, Expression)
        assert isinstance(membrane.membrane_electromigrative_flux, Expression)

    def _test_build_boundary_layer(membrane):
        # parameters
        assert isinstance(membrane.total_boundary_layer_thickness, Param)
        assert value(membrane.total_boundary_layer_thickness) == 2e-5

        # sets
        assert isinstance(
            membrane.dimensionless_boundary_layer_thickness,
            ContinuousSet,
        )
        assert len(membrane.dimensionless_boundary_layer_thickness) == 11

        # variables
        assert isinstance(membrane.boundary_layer_D_tilde, Var)
        assert len(membrane.boundary_layer_D_tilde) == 121

        # dependent on number of solutes
        assert isinstance(
            membrane.boundary_layer_cross_diffusion_coefficient_bilinear,
            Var,
        )
        assert isinstance(membrane.boundary_layer_cross_diffusion_coefficient, Var)
        assert isinstance(membrane.boundary_layer_conc_mol_comp, Var)
        assert isinstance(
            membrane.d_boundary_layer_conc_mol_comp_dz,
            DerivativeVar,
        )

        # constraints
        assert isinstance(membrane.electroneutrality_boundary_layer, Constraint)
        assert len(membrane.electroneutrality_boundary_layer) == 110

        assert isinstance(membrane.boundary_layer_D_tilde_calculation, Constraint)
        assert len(membrane.boundary_layer_D_tilde_calculation) == 110

        # dependent on number of solutes
        assert isinstance(membrane.retentate_boundary_layer_interface, Constraint)
        assert isinstance(membrane.boundary_layer_membrane_interface, Constraint)
        assert isinstance(
            membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation,
            Constraint,
        )
        assert isinstance(
            membrane.boundary_layer_cross_diffusion_coefficient_calculation,
            Constraint,
        )
        assert isinstance(membrane.cation_flux_boundary_layer, Constraint)
        assert isinstance(
            membrane.boundary_layer_conc_mol_comp_boundary_condition,
            Constraint,
        )

        # scaling factors
        assert membrane.scaling_factor[membrane.boundary_layer_D_tilde] == 1e-2
        assert (
            membrane.scaling_factor[
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear
            ]
            == 1e-3
        )
        assert (
            membrane.scaling_factor[
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            ]
            == 1e-2
        )
        assert (
            membrane.scaling_factor[
                membrane.boundary_layer_cross_diffusion_coefficient_calculation
            ]
            == 1e-2
        )

        # expressions
        assert isinstance(
            membrane.boundary_layer_electric_potential_gradient, Expression
        )
        assert len(membrane.boundary_layer_electric_potential_gradient) == 110

        # dependent on number of solutes
        assert isinstance(membrane.local_actual_rejection_percent, Expression)
        assert isinstance(membrane.local_actual_sieving_coefficient, Expression)
        assert isinstance(membrane.boundary_layer_convective_flux, Expression)
        assert isinstance(membrane.boundary_layer_diffusive_flux, Expression)
        assert isinstance(membrane.boundary_layer_electromigrative_flux, Expression)

    def _test_build_single_salt(membrane):
        # sets
        assert len(membrane.solutes) == 2
        assert len(membrane.cations) == 1

        # variables
        assert len(membrane.feed_conc_mol_comp) == 2
        assert len(membrane.diafiltrate_conc_mol_comp) == 2
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 121
        assert len(membrane.membrane_convection_coefficient_bilinear) == 121
        assert len(membrane.membrane_cross_diffusion_coefficient) == 121
        assert len(membrane.membrane_convection_coefficient) == 121
        assert len(membrane.molar_ion_flux) == 22
        assert len(membrane.retentate_conc_mol_comp) == 22
        assert len(membrane.permeate_local_conc_mol_comp) == 22
        assert len(membrane.partitioning_term_bilinear_feed) == 22
        assert len(membrane.partitioning_term_bilinear_permeate) == 22
        assert len(membrane.membrane_conc_mol_comp) == 242
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 22
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 242

        # constraints
        assert len(membrane.cation_mol_balance) == 10
        assert len(membrane.cation_bulk_flux_equation) == 10
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 20
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 20
        assert len(membrane.membrane_permeate_interface) == 20
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 110
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 110
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 110
        assert len(membrane.membrane_convection_coefficient_calculation) == 110
        assert len(membrane.cation_flux_membrane) == 110
        assert len(membrane.feed_conc_mol_comp_boundary_condition) == 1
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 22

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 2
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 2
        assert len(membrane.retentate_outlet.conc_mol_comp) == 2
        assert len(membrane.permeate_outlet.conc_mol_comp) == 2

        # expressions
        assert len(membrane.permeate_outlet_conc_mol_comp) == 2
        assert len(membrane.total_feed_conc_mol_comp) == 2
        assert len(membrane.local_partition_coefficient_feed_side) == 22
        assert len(membrane.local_partition_coefficient_permeate_side) == 22
        assert len(membrane.local_observed_rejection_percent) == 22
        assert len(membrane.overall_observed_rejection_percent) == 2
        assert len(membrane.local_observed_sieving_coefficient) == 22
        assert len(membrane.overall_observed_sieving_coefficient) == 2
        assert len(membrane.membrane_convective_flux) == 220
        assert len(membrane.membrane_diffusive_flux) == 220
        assert len(membrane.membrane_electromigrative_flux) == 220

    def _test_build_single_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 121
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 121
        assert len(membrane.boundary_layer_conc_mol_comp) == 242
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 242

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 10
        assert len(membrane.boundary_layer_membrane_interface) == 20
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 110
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 110
        )
        assert len(membrane.cation_flux_boundary_layer) == 100
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 22

        # expressions
        assert len(membrane.local_actual_rejection_percent) == 22
        assert len(membrane.local_actual_sieving_coefficient) == 22
        assert len(membrane.boundary_layer_convective_flux) == 220
        assert len(membrane.boundary_layer_diffusive_flux) == 220
        assert len(membrane.boundary_layer_electromigrative_flux) == 220

    def _test_build_single_salt_no_boundary_layer(membrane):
        assert isinstance(membrane.retentate_membrane_interface, Constraint)
        assert len(membrane.retentate_membrane_interface) == 20

    def _test_build_two_salt(membrane):
        # sets
        assert len(membrane.solutes) == 3
        assert len(membrane.cations) == 2

        # variables
        assert len(membrane.feed_conc_mol_comp) == 3
        assert len(membrane.diafiltrate_conc_mol_comp) == 3
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 484
        assert len(membrane.membrane_convection_coefficient_bilinear) == 242
        assert len(membrane.membrane_cross_diffusion_coefficient) == 484
        assert len(membrane.membrane_convection_coefficient) == 242
        assert len(membrane.molar_ion_flux) == 33
        assert len(membrane.retentate_conc_mol_comp) == 33
        assert len(membrane.permeate_local_conc_mol_comp) == 33
        assert len(membrane.partitioning_term_bilinear_feed) == 33
        assert len(membrane.partitioning_term_bilinear_permeate) == 33
        assert len(membrane.membrane_conc_mol_comp) == 363
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 33
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 363

        # constraints
        assert len(membrane.cation_mol_balance) == 20
        assert len(membrane.cation_bulk_flux_equation) == 20
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 30
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 30
        assert len(membrane.membrane_permeate_interface) == 30
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 440
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 220
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 440
        assert len(membrane.membrane_convection_coefficient_calculation) == 220
        assert len(membrane.cation_flux_membrane) == 220
        assert len(membrane.feed_conc_mol_comp_boundary_condition) == 2
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 33

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 3
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 3
        assert len(membrane.retentate_outlet.conc_mol_comp) == 3
        assert len(membrane.permeate_outlet.conc_mol_comp) == 3

        # expressions
        assert len(membrane.permeate_outlet_conc_mol_comp) == 3
        assert len(membrane.total_feed_conc_mol_comp) == 3
        assert len(membrane.local_partition_coefficient_feed_side) == 33
        assert len(membrane.local_partition_coefficient_permeate_side) == 33
        assert len(membrane.local_observed_rejection_percent) == 33
        assert len(membrane.overall_observed_rejection_percent) == 3
        assert len(membrane.local_observed_sieving_coefficient) == 33
        assert len(membrane.overall_observed_sieving_coefficient) == 3
        assert len(membrane.membrane_convective_flux) == 330
        assert len(membrane.membrane_diffusive_flux) == 330
        assert len(membrane.membrane_electromigrative_flux) == 330

    def _test_build_two_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 484
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 484
        assert len(membrane.boundary_layer_conc_mol_comp) == 363
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 363

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 20
        assert len(membrane.boundary_layer_membrane_interface) == 30
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 440
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 440
        )
        assert len(membrane.cation_flux_boundary_layer) == 200
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 33

        # expressions
        assert len(membrane.local_actual_rejection_percent) == 33
        assert len(membrane.local_actual_sieving_coefficient) == 33
        assert len(membrane.boundary_layer_convective_flux) == 330
        assert len(membrane.boundary_layer_diffusive_flux) == 330
        assert len(membrane.boundary_layer_electromigrative_flux) == 330

    def _test_build_two_salt_no_boundary_layer(membrane):
        assert isinstance(membrane.retentate_membrane_interface, Constraint)
        assert len(membrane.retentate_membrane_interface) == 30

    def _test_build_three_salt(membrane):
        # sets
        assert len(membrane.solutes) == 4
        assert len(membrane.cations) == 3

        # variables
        assert len(membrane.feed_conc_mol_comp) == 4
        assert len(membrane.diafiltrate_conc_mol_comp) == 4
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 1089
        assert len(membrane.membrane_convection_coefficient_bilinear) == 363
        assert len(membrane.membrane_cross_diffusion_coefficient) == 1089
        assert len(membrane.membrane_convection_coefficient) == 363
        assert len(membrane.molar_ion_flux) == 44
        assert len(membrane.retentate_conc_mol_comp) == 44
        assert len(membrane.permeate_local_conc_mol_comp) == 44
        assert len(membrane.partitioning_term_bilinear_feed) == 44
        assert len(membrane.partitioning_term_bilinear_permeate) == 44
        assert len(membrane.membrane_conc_mol_comp) == 484
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 44
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 484

        # constraints
        assert len(membrane.cation_mol_balance) == 30
        assert len(membrane.cation_bulk_flux_equation) == 30
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 40
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 40
        assert len(membrane.membrane_permeate_interface) == 40
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 990
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 330
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 990
        assert len(membrane.membrane_convection_coefficient_calculation) == 330
        assert len(membrane.cation_flux_membrane) == 330
        assert len(membrane.feed_conc_mol_comp_boundary_condition) == 3
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 44

        # scaling factors
        assert membrane.scaling_factor[membrane.molar_ion_flux] == 1e2

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 4
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 4
        assert len(membrane.retentate_outlet.conc_mol_comp) == 4
        assert len(membrane.permeate_outlet.conc_mol_comp) == 4

        # expressions
        assert len(membrane.permeate_outlet_conc_mol_comp) == 4
        assert len(membrane.total_feed_conc_mol_comp) == 4
        assert len(membrane.local_partition_coefficient_feed_side) == 44
        assert len(membrane.local_partition_coefficient_permeate_side) == 44
        assert len(membrane.local_observed_rejection_percent) == 44
        assert len(membrane.overall_observed_rejection_percent) == 4
        assert len(membrane.local_observed_sieving_coefficient) == 44
        assert len(membrane.overall_observed_sieving_coefficient) == 4
        assert len(membrane.membrane_convective_flux) == 440
        assert len(membrane.membrane_diffusive_flux) == 440
        assert len(membrane.membrane_electromigrative_flux) == 440

    def _test_build_three_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 1089
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 1089
        assert len(membrane.boundary_layer_conc_mol_comp) == 484
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 484

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 30
        assert len(membrane.boundary_layer_membrane_interface) == 40
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 990
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 990
        )
        assert len(membrane.cation_flux_boundary_layer) == 300
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 44

        # expressions
        assert len(membrane.local_actual_rejection_percent) == 44
        assert len(membrane.local_actual_sieving_coefficient) == 44
        assert len(membrane.boundary_layer_convective_flux) == 440
        assert len(membrane.boundary_layer_diffusive_flux) == 440
        assert len(membrane.boundary_layer_electromigrative_flux) == 440

    def _test_build_three_salt_no_boundary_layer(membrane):
        assert isinstance(membrane.retentate_membrane_interface, Constraint)
        assert len(membrane.retentate_membrane_interface) == 40

    model_single_salt = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 200},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_build(model_single_salt.fs.unit)
    _test_build_boundary_layer(model_single_salt.fs.unit)
    _test_build_single_salt(model_single_salt.fs.unit)
    _test_build_single_salt_boundary_layer(model_single_salt.fs.unit)

    model_single_salt_no_boundary_layer = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 2400},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_build(model_single_salt_no_boundary_layer.fs.unit)
    _test_build_single_salt(model_single_salt_no_boundary_layer.fs.unit)
    _test_build_single_salt_no_boundary_layer(
        model_single_salt_no_boundary_layer.fs.unit
    )

    model_two_salt = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Cl": 600},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_build(model_two_salt.fs.unit)
    _test_build_boundary_layer(model_two_salt.fs.unit)
    _test_build_two_salt(model_two_salt.fs.unit)
    _test_build_two_salt_boundary_layer(model_two_salt.fs.unit)

    model_two_salt_no_boundary_layer = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Cl": 600},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_build(model_two_salt_no_boundary_layer.fs.unit)
    _test_build_two_salt(model_two_salt_no_boundary_layer.fs.unit)
    _test_build_two_salt_no_boundary_layer(model_two_salt_no_boundary_layer.fs.unit)

    model_three_salt = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Al": 20, "Cl": 660},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_build(model_three_salt.fs.unit)
    _test_build_boundary_layer(model_three_salt.fs.unit)
    _test_build_three_salt(model_three_salt.fs.unit)
    _test_build_three_salt_boundary_layer(model_three_salt.fs.unit)

    model_three_salt = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Al": 20, "Cl": 660},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_build(model_three_salt.fs.unit)
    _test_build_three_salt(model_three_salt.fs.unit)
    _test_build_three_salt_no_boundary_layer(model_three_salt.fs.unit)


################################################################################
# Test model solves
@pytest.mark.solver
@pytest.mark.component
def test_solve(diafiltration_model):
    def _test_diagnostics(membrane):
        dt = DiagnosticsToolbox(membrane)
        dt.assert_no_structural_warnings()

    def _solve_model(model):
        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(model, rename=False)
        solver = SolverFactory("ipopt")
        results = solver.solve(scaled_model, tee=True)
        assert_optimal_termination(results)
        scaling.propagate_solution(scaled_model, model)

    def _set_water_flux_target(model):
        model.fs.unit.applied_pressure.unfix()

        def _water_flux_constraint(m):
            return (
                sum(
                    m.fs.unit.volume_flux_water[0, x]
                    for x in m.fs.unit.dimensionless_module_length
                    if x != 0
                )
                / (len(m.fs.unit.dimensionless_module_length) - 1)
                == 0.02
            )

        model.water_flux_constraint = Constraint(rule=_water_flux_constraint)

    def _test_solve(model):
        _solve_model(model)
        _set_water_flux_target(model)
        _solve_model(model)

    def _test_numerical_issues(membrane):
        dt = DiagnosticsToolbox(membrane)
        dt.assert_no_numerical_warnings()

    model_LiCl = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 200},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiCl.fs.unit)
    _test_solve(model_LiCl)
    _test_numerical_issues(model_LiCl.fs.unit)
    LiCl_test_dict = {
        "applied_pressure": {(0): (4.8035, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (169.74, 1e-4, None),
            (0, 1, "Cl"): (169.74, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (108.25, 1e-4, None),
            (0, "Cl"): (108.25, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiCl.fs.unit, LiCl_test_dict)

    model_LiCl_no_boundary_layer = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Cl": 200},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_diagnostics(model_LiCl_no_boundary_layer.fs.unit)
    _test_solve(model_LiCl_no_boundary_layer)
    _test_numerical_issues(model_LiCl_no_boundary_layer.fs.unit)
    LiCl_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (4.7434, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (170.54, 1e-4, None),
            (0, 1, "Cl"): (170.54, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (105.17, 1e-4, None),
            (0, "Cl"): (105.17, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_LiCl_no_boundary_layer.fs.unit, LiCl_no_boundary_layer_test_dict
    )

    model_CoCl2 = diafiltration_model(
        cation_list=["Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Co": 200, "Cl": 400},
            "diafiltrate": {"Co": 3, "Cl": 6},
        },
        non_Donnan_partition_dict={"Co": 0.4, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_CoCl2.fs.unit)
    _test_solve(model_CoCl2)
    _test_numerical_issues(model_CoCl2.fs.unit)
    CoCl2_test_dict = {
        "applied_pressure": {(0): (8.4877, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Co"): (173.29, 1e-4, None),
            (0, 1, "Cl"): (346.59, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Co"): (82.248, 1e-4, None),
            (0, "Cl"): (164.50, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_CoCl2.fs.unit, CoCl2_test_dict)

    model_CoCl2_no_boundary_layer = diafiltration_model(
        cation_list=["Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Co": 200, "Cl": 400},
            "diafiltrate": {"Co": 3, "Cl": 6},
        },
        non_Donnan_partition_dict={"Co": 0.4, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_diagnostics(model_CoCl2_no_boundary_layer.fs.unit)
    _test_solve(model_CoCl2_no_boundary_layer)
    _test_numerical_issues(model_CoCl2_no_boundary_layer.fs.unit)
    CoCl2_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (8.2255, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Co"): (174.20, 1e-4, None),
            (0, 1, "Cl"): (348.40, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Co"): (78.761, 1e-4, None),
            (0, "Cl"): (157.52, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_CoCl2_no_boundary_layer.fs.unit, CoCl2_no_boundary_layer_test_dict
    )

    model_AlCl3 = diafiltration_model(
        cation_list=["Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Al": 20, "Cl": 60},
            "diafiltrate": {"Al": 3, "Cl": 9},
        },
        non_Donnan_partition_dict={"Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_AlCl3.fs.unit)
    _test_solve(model_AlCl3)
    _test_numerical_issues(model_AlCl3.fs.unit)
    AlCl3_test_dict = {
        "applied_pressure": {(0): (3.7290, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Al"): (20.092, 1e-4, None),
            (0, 1, "Cl"): (60.277, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Al"): (0.59824, 1e-4, None),
            (0, "Cl"): (1.7947, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_AlCl3.fs.unit, AlCl3_test_dict)

    model_AlCl3_no_boundary_layer = diafiltration_model(
        cation_list=["Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Al": 20, "Cl": 60},
            "diafiltrate": {"Al": 3, "Cl": 9},
        },
        non_Donnan_partition_dict={"Al": 0.0001, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_diagnostics(model_AlCl3_no_boundary_layer.fs.unit)
    _test_solve(model_AlCl3_no_boundary_layer)
    _test_numerical_issues(model_AlCl3_no_boundary_layer.fs.unit)
    AlCl3_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (3.5704, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Al"): (20.109, 1e-4, None),
            (0, 1, "Cl"): (60.328, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Al"): (0.53274, 1e-4, None),
            (0, "Cl"): (1.5982, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_AlCl3_no_boundary_layer.fs.unit, AlCl3_no_boundary_layer_test_dict
    )

    model_LiCoCl3 = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Cl": 600},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiCoCl3.fs.unit)
    _test_solve(model_LiCoCl3)
    _test_numerical_issues(model_LiCoCl3.fs.unit)
    LiCoCl3_test_dict = {
        "applied_pressure": {(0): (11.356, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (162.79, 1e-4, None),
            (0, 1, "Co"): (177.57, 1e-4, None),
            (0, 1, "Cl"): (517.93, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (135.05, 1e-4, None),
            (0, "Co"): (65.776, 1e-4, None),
            (0, "Cl"): (266.61, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiCoCl3.fs.unit, LiCoCl3_test_dict)

    model_LiCoCl3_no_boundary_layer = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Cl": 600},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_diagnostics(model_LiCoCl3_no_boundary_layer.fs.unit)
    _test_solve(model_LiCoCl3_no_boundary_layer)
    _test_numerical_issues(model_LiCoCl3_no_boundary_layer.fs.unit)
    LiCoCl3_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (10.960, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (163.00, 1e-4, None),
            (0, 1, "Co"): (178.63, 1e-4, None),
            (0, 1, "Cl"): (520.27, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (134.23, 1e-4, None),
            (0, "Co"): (61.667, 1e-4, None),
            (0, "Cl"): (257.57, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_LiCoCl3_no_boundary_layer.fs.unit, LiCoCl3_no_boundary_layer_test_dict
    )

    # TODO: For completeness, debug numerical issues in no boundary layer two-salt solutions with Al
    # +3 ions are almost 100% rejected, likely causing near-zero instabilities
    # This is low priority

    model_LiAlCl4 = diafiltration_model(
        cation_list=["Li", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Al": 20, "Cl": 260},
            "diafiltrate": {"Li": 14, "Al": 3, "Cl": 23},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiAlCl4.fs.unit)
    _test_solve(model_LiAlCl4)
    _test_numerical_issues(model_LiAlCl4.fs.unit)
    LiAlCl4_test_dict = {
        "applied_pressure": {(0): (5.7847, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (165.48, 1e-4, None),
            (0, 1, "Al"): (20.248, 1e-4, None),
            (0, 1, "Cl"): (226.23, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (124.67, 1e-4, None),
            (0, "Al"): (0.0011345, 1e-4, None),
            (0, "Cl"): (124.68, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiAlCl4.fs.unit, LiAlCl4_test_dict)

    model_CoAlCl5 = diafiltration_model(
        cation_list=["Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Co": 200, "Al": 20, "Cl": 460},
            "diafiltrate": {"Co": 3, "Al": 3, "Cl": 15},
        },
        non_Donnan_partition_dict={"Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_CoAlCl5.fs.unit)
    _test_solve(model_CoAlCl5)
    _test_numerical_issues(model_CoAlCl5.fs.unit)
    CoAlCl5_test_dict = {
        "applied_pressure": {(0): (9.6966, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Co"): (171.11, 1e-4, None),
            (0, 1, "Al"): (20.247, 1e-4, None),
            (0, 1, "Cl"): (402.96, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Co"): (90.669, 1e-4, None),
            (0, "Al"): (0.0031789, 1e-4, None),
            (0, "Cl"): (181.35, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_CoAlCl5.fs.unit, CoAlCl5_test_dict)

    model_LiCoAlCl6 = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Al": 20, "Cl": 660},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiCoAlCl6.fs.unit)
    _test_solve(model_LiCoAlCl6)
    _test_numerical_issues(model_LiCoAlCl6.fs.unit)
    LiCoAlCl6_test_dict = {
        "applied_pressure": {(0): (12.547, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (161.36, 1e-4, None),
            (0, 1, "Co"): (176.28, 1e-4, None),
            (0, 1, "Al"): (20.247, 1e-4, None),
            (0, 1, "Cl"): (574.66, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (140.55, 1e-4, None),
            (0, "Co"): (70.752, 1e-4, None),
            (0, "Al"): (0.0011973, 1e-4, None),
            (0, "Cl"): (282.06, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiCoAlCl6.fs.unit, LiCoAlCl6_test_dict)

    model_LiCoAlCl6_no_boundary_layer = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 200, "Co": 200, "Al": 20, "Cl": 660},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.75, "Co": 0.4, "Al": 0.0001, "Cl": 0.05},
        boundary_layer=False,
    )
    _test_diagnostics(model_LiCoAlCl6_no_boundary_layer.fs.unit)
    _test_solve(model_LiCoAlCl6_no_boundary_layer)
    _test_numerical_issues(model_LiCoAlCl6_no_boundary_layer.fs.unit)
    LiCoAlCl6_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (12.052, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (161.51, 1e-4, None),
            (0, 1, "Co"): (177.42, 1e-4, None),
            (0, 1, "Al"): (20.247, 1e-4, None),
            (0, 1, "Cl"): (577.09, 1e-4, None),
        },
        "permeate_outlet_flow": {(0): (3.2800, 1e-4, None)},
        "permeate_outlet_conc_mol_comp": {
            (0, "Li"): (140.00, 1e-4, None),
            (0, "Co"): (66.345, 1e-4, None),
            (0, "Al"): (0.0010161, 1e-4, None),
            (0, "Cl"): (272.70, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_LiCoAlCl6_no_boundary_layer.fs.unit, LiCoAlCl6_no_boundary_layer_test_dict
    )
