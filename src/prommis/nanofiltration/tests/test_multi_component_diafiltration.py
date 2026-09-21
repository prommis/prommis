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
    units,
    value,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util.diagnostics_tools.diagnostics_toolbox import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import assert_solution_equivalent

import pytest

from prommis.nanofiltration.multi_component_diafiltration_solute_properties import (
    MultiComponentDiafiltrationSoluteParameter,
)
from prommis.nanofiltration.multi_component_diafiltration import (
    MultiComponentDiafiltration,
)

# TODO: test positive and neutral membrane cases


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
        Dm_over_l_value = 40  # um/s
        l_um = value(Dm_Cl) / Dm_over_l_value  # um
        l_m = l_um / 1e6  # m

        m.fs.unit = MultiComponentDiafiltration(
            property_package=m.fs.properties,
            cation_list=cation_list,
            anion_list=anion_list,
            include_boundary_layer=boundary_layer,
            total_membrane_thickness=l_m,
            NFE_module_length=10,
            NFE_boundary_layer_thickness=5,
            NFE_membrane_thickness=5,
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

        if feed_ionic_strength < 199:
            m.fs.unit.applied_pressure.fix(5)
        elif (feed_ionic_strength >= 199) and (feed_ionic_strength < 299):
            m.fs.unit.applied_pressure.fix(15)
        elif feed_ionic_strength >= 299:
            m.fs.unit.applied_pressure.fix(20)

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
        assert m.fs.unit.config.NFE_boundary_layer_thickness == 5
        assert m.fs.unit.config.NFE_membrane_thickness == 5

    model = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
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
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
        boundary_layer=False,
    )
    _test_config(model_no_boundary_layer)
    assert model_no_boundary_layer.fs.unit.config.include_boundary_layer == False

    model_two_salt = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Co": 288, "Cl": 821},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Cl": 0.1},
        boundary_layer=True,
    )
    _test_config(model_two_salt)
    assert len(model_two_salt.fs.unit.config.cation_list) == 2

    model_three_salt = diafiltration_model(
        cation_list=["Li", "Co", "Al"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Al": 0.0005, "Cl": 0.1},
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
        assert len(membrane.dimensionless_membrane_thickness) == 6

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
        assert len(membrane.membrane_D_tilde) == 66

        assert isinstance(membrane.volume_flux_water, Var)
        assert len(membrane.volume_flux_water) == 11

        assert isinstance(membrane.retentate_flow_volume, Var)
        assert len(membrane.retentate_flow_volume) == 11

        assert isinstance(membrane.permeate_flow_volume, Var)
        assert len(membrane.permeate_flow_volume) == 11

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
        assert isinstance(membrane.permeate_conc_mol_comp, Var)
        assert isinstance(membrane.partitioning_term_bilinear_feed, Var)
        assert isinstance(membrane.partitioning_term_bilinear_permeate, Var)
        assert isinstance(membrane.d_retentate_conc_mol_comp_dx, DerivativeVar)
        assert isinstance(membrane.d_membrane_conc_mol_comp_dz, DerivativeVar)

        # constraints
        assert isinstance(membrane.differential_overall_mass_balance, Constraint)
        assert len(membrane.differential_overall_mass_balance) == 10

        assert isinstance(membrane.overall_mass_balance, Constraint)
        assert len(membrane.overall_mass_balance) == 11

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
            membrane.permeate_flow_volume_boundary_condition,
            Constraint,
        )
        assert len(membrane.permeate_flow_volume_boundary_condition) == 1

        assert isinstance(
            membrane.volume_flux_water_boundary_condition,
            Constraint,
        )
        assert len(membrane.volume_flux_water_boundary_condition) == 1

        # dependent on number of solutes
        assert isinstance(membrane.differential_cation_mol_balance, Constraint)
        assert isinstance(membrane.cation_mol_balance, Constraint)
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
            membrane.membrane_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert isinstance(
            membrane.permeate_conc_mol_comp_boundary_condition,
            Constraint,
        )
        assert isinstance(
            membrane.molar_ion_flux_boundary_condition,
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

    def _test_build_boundary_layer(membrane):
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
        assert len(membrane.electroneutrality_boundary_layer) == 60

        assert isinstance(membrane.boundary_layer_D_tilde_calculation, Constraint)
        assert len(membrane.boundary_layer_D_tilde_calculation) == 60

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

    def _test_build_single_salt(membrane):
        # sets
        assert len(membrane.solutes) == 2
        assert len(membrane.cations) == 1

        # variables
        assert len(membrane.feed_conc_mol_comp) == 2
        assert len(membrane.diafiltrate_conc_mol_comp) == 2
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 66
        assert len(membrane.membrane_convection_coefficient_bilinear) == 66
        assert len(membrane.membrane_cross_diffusion_coefficient) == 66
        assert len(membrane.membrane_convection_coefficient) == 66
        assert len(membrane.molar_ion_flux) == 22
        assert len(membrane.retentate_conc_mol_comp) == 22
        assert len(membrane.permeate_conc_mol_comp) == 22
        assert len(membrane.partitioning_term_bilinear_feed) == 22
        assert len(membrane.partitioning_term_bilinear_permeate) == 22
        assert len(membrane.membrane_conc_mol_comp) == 132
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 22
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 132

        # constraints
        assert len(membrane.differential_cation_mol_balance) == 10
        assert len(membrane.cation_mol_balance) == 11
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 20
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 20
        assert len(membrane.membrane_permeate_interface) == 20
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 60
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 60
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 60
        assert len(membrane.membrane_convection_coefficient_calculation) == 60
        assert len(membrane.cation_flux_membrane) == 60
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 12
        assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 2
        assert len(membrane.molar_ion_flux_boundary_condition) == 2

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 2
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 2
        assert len(membrane.retentate_outlet.conc_mol_comp) == 2
        assert len(membrane.permeate_outlet.conc_mol_comp) == 2

    def _test_build_single_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 66
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 66
        assert len(membrane.boundary_layer_conc_mol_comp) == 132
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 132

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 10
        assert len(membrane.boundary_layer_membrane_interface) == 20
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 60
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 60
        )
        assert len(membrane.cation_flux_boundary_layer) == 50
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 12

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
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 264
        assert len(membrane.membrane_convection_coefficient_bilinear) == 132
        assert len(membrane.membrane_cross_diffusion_coefficient) == 264
        assert len(membrane.membrane_convection_coefficient) == 132
        assert len(membrane.molar_ion_flux) == 33
        assert len(membrane.retentate_conc_mol_comp) == 33
        assert len(membrane.permeate_conc_mol_comp) == 33
        assert len(membrane.partitioning_term_bilinear_feed) == 33
        assert len(membrane.partitioning_term_bilinear_permeate) == 33
        assert len(membrane.membrane_conc_mol_comp) == 198
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 33
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 198

        # constraints
        assert len(membrane.differential_cation_mol_balance) == 20
        assert len(membrane.cation_mol_balance) == 22
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 30
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 30
        assert len(membrane.membrane_permeate_interface) == 30
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 240
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 120
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 240
        assert len(membrane.membrane_convection_coefficient_calculation) == 120
        assert len(membrane.cation_flux_membrane) == 120
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 18
        assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 3
        assert len(membrane.molar_ion_flux_boundary_condition) == 3

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 3
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 3
        assert len(membrane.retentate_outlet.conc_mol_comp) == 3
        assert len(membrane.permeate_outlet.conc_mol_comp) == 3

    def _test_build_two_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 264
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 264
        assert len(membrane.boundary_layer_conc_mol_comp) == 198
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 198

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 20
        assert len(membrane.boundary_layer_membrane_interface) == 30
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 240
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 240
        )
        assert len(membrane.cation_flux_boundary_layer) == 100
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 18

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
        assert len(membrane.membrane_cross_diffusion_coefficient_bilinear) == 594
        assert len(membrane.membrane_convection_coefficient_bilinear) == 198
        assert len(membrane.membrane_cross_diffusion_coefficient) == 594
        assert len(membrane.membrane_convection_coefficient) == 198
        assert len(membrane.molar_ion_flux) == 44
        assert len(membrane.retentate_conc_mol_comp) == 44
        assert len(membrane.permeate_conc_mol_comp) == 44
        assert len(membrane.partitioning_term_bilinear_feed) == 44
        assert len(membrane.partitioning_term_bilinear_permeate) == 44
        assert len(membrane.membrane_conc_mol_comp) == 264
        assert len(membrane.d_retentate_conc_mol_comp_dx) == 44
        assert len(membrane.d_membrane_conc_mol_comp_dz) == 264

        # constraints
        assert len(membrane.differential_cation_mol_balance) == 30
        assert len(membrane.cation_mol_balance) == 33
        assert len(membrane.partitioning_term_bilinear_feed_constraint) == 40
        assert len(membrane.partitioning_term_bilinear_permeate_constraint) == 40
        assert len(membrane.membrane_permeate_interface) == 40
        assert (
            len(membrane.membrane_cross_diffusion_coefficient_bilinear_calculation)
            == 540
        )
        assert len(membrane.membrane_convection_coefficient_bilinear_calculation) == 180
        assert len(membrane.membrane_cross_diffusion_coefficient_calculation) == 540
        assert len(membrane.membrane_convection_coefficient_calculation) == 180
        assert len(membrane.cation_flux_membrane) == 180
        assert len(membrane.membrane_conc_mol_comp_boundary_condition) == 24
        assert len(membrane.permeate_conc_mol_comp_boundary_condition) == 4
        assert len(membrane.molar_ion_flux_boundary_condition) == 4

        # ports
        assert len(membrane.feed_inlet.conc_mol_comp) == 4
        assert len(membrane.diafiltrate_inlet.conc_mol_comp) == 4
        assert len(membrane.retentate_outlet.conc_mol_comp) == 4
        assert len(membrane.permeate_outlet.conc_mol_comp) == 4

    def _test_build_three_salt_boundary_layer(membrane):
        # variables
        assert len(membrane.boundary_layer_cross_diffusion_coefficient_bilinear) == 594
        assert len(membrane.boundary_layer_cross_diffusion_coefficient) == 594
        assert len(membrane.boundary_layer_conc_mol_comp) == 264
        assert len(membrane.d_boundary_layer_conc_mol_comp_dz) == 264

        # constraints
        assert len(membrane.retentate_boundary_layer_interface) == 30
        assert len(membrane.boundary_layer_membrane_interface) == 40
        assert (
            len(
                membrane.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            == 540
        )
        assert (
            len(membrane.boundary_layer_cross_diffusion_coefficient_calculation) == 540
        )
        assert len(membrane.cation_flux_boundary_layer) == 150
        assert len(membrane.boundary_layer_conc_mol_comp_boundary_condition) == 24

    def _test_build_three_salt_no_boundary_layer(membrane):
        assert isinstance(membrane.retentate_membrane_interface, Constraint)
        assert len(membrane.retentate_membrane_interface) == 40

    model_single_salt = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
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
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
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
            "feed": {"Li": 245, "Co": 288, "Cl": 821},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Cl": 0.1},
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
            "feed": {"Li": 245, "Co": 288, "Cl": 821},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Cl": 0.1},
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
            "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Al": 0.0005, "Cl": 0.1},
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
            "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
            "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Al": 0.0005, "Cl": 0.1},
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
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiCl.fs.unit)
    _test_solve(model_LiCl)
    _test_numerical_issues(model_LiCl.fs.unit)
    LiCl_test_dict = {
        "applied_pressure": {(0): (5.2828, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (205.21, 1e-4, None),
            (0, 1, "Cl"): (205.21, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Li"): (138.22, 1e-4, None),
            (0, 1, "Cl"): (138.22, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiCl.fs.unit, LiCl_test_dict)

    model_LiCl_no_boundary_layer = diafiltration_model(
        cation_list=["Li"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Cl": 245},
            "diafiltrate": {"Li": 14, "Cl": 14},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Cl": 0.1},
        boundary_layer=False,
    )
    _test_diagnostics(model_LiCl_no_boundary_layer.fs.unit)
    _test_solve(model_LiCl_no_boundary_layer)
    _test_numerical_issues(model_LiCl_no_boundary_layer.fs.unit)
    LiCl_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (5.2255, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (205.96, 1e-4, None),
            (0, 1, "Cl"): (205.96, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Li"): (135.27, 1e-4, None),
            (0, 1, "Cl"): (135.27, 1e-4, None),
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
            "feed": {"Co": 288, "Cl": 576},
            "diafiltrate": {"Co": 3, "Cl": 6},
        },
        non_Donnan_partition_dict={"Co": 0.3, "Cl": 0.1},
        boundary_layer=True,
    )
    _test_diagnostics(model_CoCl2.fs.unit)
    _test_solve(model_CoCl2)
    _test_numerical_issues(model_CoCl2.fs.unit)
    CoCl2_test_dict = {
        "applied_pressure": {(0): (18.145, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Co"): (247.51, 1e-4, None),
            (0, 1, "Cl"): (495.02, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Co"): (122.28, 1e-4, None),
            (0, 1, "Cl"): (244.56, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_CoCl2.fs.unit, CoCl2_test_dict)

    model_CoCl2_no_boundary_layer = diafiltration_model(
        cation_list=["Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Co": 288, "Cl": 576},
            "diafiltrate": {"Co": 3, "Cl": 6},
        },
        non_Donnan_partition_dict={"Co": 0.3, "Cl": 0.1},
        boundary_layer=False,
    )
    _test_diagnostics(model_CoCl2_no_boundary_layer.fs.unit)
    _test_solve(model_CoCl2_no_boundary_layer)
    _test_numerical_issues(model_CoCl2_no_boundary_layer.fs.unit)
    CoCl2_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (17.668, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Co"): (248.50, 1e-4, None),
            (0, 1, "Cl"): (497.01, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Co"): (118.34, 1e-4, None),
            (0, 1, "Cl"): (236.67, 1e-4, None),
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
        non_Donnan_partition_dict={"Al": 0.0005, "Cl": 0.1},
        boundary_layer=True,
    )
    _test_diagnostics(model_AlCl3.fs.unit)
    _test_solve(model_AlCl3)
    _test_numerical_issues(model_AlCl3.fs.unit)
    AlCl3_test_dict = {
        "applied_pressure": {(0): (6.5109, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Al"): (19.889, 1e-4, None),
            (0, 1, "Cl"): (59.668, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Al"): (1.0023, 1e-4, None),
            (0, 1, "Cl"): (3.0068, 1e-4, None),
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
        non_Donnan_partition_dict={"Al": 0.0005, "Cl": 0.1},
        boundary_layer=False,
    )
    _test_diagnostics(model_AlCl3_no_boundary_layer.fs.unit)
    _test_solve(model_AlCl3_no_boundary_layer)
    _test_numerical_issues(model_AlCl3_no_boundary_layer.fs.unit)
    AlCl3_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (6.1631, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Al"): (19.926, 1e-4, None),
            (0, 1, "Cl"): (59.778, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Al"): (0.85716, 1e-4, None),
            (0, 1, "Cl"): (2.5715, 1e-4, None),
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
            "feed": {"Li": 245, "Co": 288, "Cl": 821},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Cl": 0.1},
        boundary_layer=True,
    )
    _test_diagnostics(model_LiCoCl3.fs.unit)
    _test_solve(model_LiCoCl3)
    _test_numerical_issues(model_LiCoCl3.fs.unit)
    LiCoCl3_test_dict = {
        "applied_pressure": {(0): (33.533, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (195.37, 1e-4, None),
            (0, 1, "Co"): (255.02, 1e-4, None),
            (0, 1, "Cl"): (705.41, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Li"): (177.15, 1e-4, None),
            (0, 1, "Co"): (92.576, 1e-4, None),
            (0, 1, "Cl"): (362.30, 1e-4, None),
        },
    }
    assert_solution_equivalent(model_LiCoCl3.fs.unit, LiCoCl3_test_dict)

    model_LiCoCl3_no_boundary_layer = diafiltration_model(
        cation_list=["Li", "Co"],
        anion_list=["Cl"],
        inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
        inlet_concentration={
            "feed": {"Li": 245, "Co": 288, "Cl": 821},
            "diafiltrate": {"Li": 14, "Co": 3, "Cl": 20},
        },
        non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Cl": 0.1},
        boundary_layer=False,
    )
    _test_diagnostics(model_LiCoCl3_no_boundary_layer.fs.unit)
    _test_solve(model_LiCoCl3_no_boundary_layer)
    _test_numerical_issues(model_LiCoCl3_no_boundary_layer.fs.unit)
    LiCoCl3_no_boundary_layer_test_dict = {
        "applied_pressure": {(0): (32.532, 1e-4, None)},
        "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
        "retentate_conc_mol_comp": {
            (0, 1, "Li"): (194.96, 1e-4, None),
            (0, 1, "Co"): (256.47, 1e-4, None),
            (0, 1, "Cl"): (707.91, 1e-4, None),
        },
        "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
        "permeate_conc_mol_comp": {
            (0, 1, "Li"): (178.77, 1e-4, None),
            (0, 1, "Co"): (86.825, 1e-4, None),
            (0, 1, "Cl"): (352.42, 1e-4, None),
        },
    }
    assert_solution_equivalent(
        model_LiCoCl3_no_boundary_layer.fs.unit, LiCoCl3_no_boundary_layer_test_dict
    )

    # TODO: debug numerical issues in multi-salt solutions with Al (+3 ions almost 100% rejected, causes near-zero instabilities)
    # model_LiAlCl4 = diafiltration_model(
    #     cation_list=["Li", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Li": 245, "Al": 20, "Cl": 305},
    #         "diafiltrate": {"Li": 14, "Al": 3, "Cl": 23},
    #     },
    #     non_Donnan_partition_dict={"Li": 0.7, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=True,
    # )
    # _test_diagnostics(model_LiAlCl4.fs.unit)
    # _test_solve(model_LiAlCl4)
    # _test_numerical_issues(model_LiAlCl4.fs.unit)
    # LiAlCl4_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(model_LiAlCl4.fs.unit, LiAlCl4_test_dict)

    # model_LiAlCl4_no_boundary_layer = diafiltration_model(
    #     cation_list=["Li", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Li": 245, "Al": 20, "Cl": 305},
    #         "diafiltrate": {"Li": 14, "Al": 3, "Cl": 23},
    #     },
    #     non_Donnan_partition_dict={"Li": 0.7, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=False,
    # )
    # _test_diagnostics(model_LiAlCl4_no_boundary_layer.fs.unit)
    # _test_solve(model_LiAlCl4_no_boundary_layer)
    # _test_numerical_issues(model_LiAlCl4_no_boundary_layer.fs.unit)
    # LiAlCl4_no_boundary_layer_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(
    #     model_LiAlCl4_no_boundary_layer.fs.unit, LiAlCl4_no_boundary_layer_test_dict
    # )

    # model_CoAlCl5 = diafiltration_model(
    #     cation_list=["Co", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Co": 288, "Al": 20, "Cl": 636},
    #         "diafiltrate": {"Co": 3, "Al": 3, "Cl": 15},
    #     },
    #     non_Donnan_partition_dict={"Co": 0.3, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=True,
    # )
    # _test_diagnostics(model_CoAlCl5.fs.unit)
    # _test_solve(model_CoAlCl5)
    # _test_numerical_issues(model_CoAlCl5.fs.unit)
    # CoAlCl5_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(model_CoAlCl5.fs.unit, CoAlCl5_test_dict)

    # model_CoAlCl5_no_boundary_layer = diafiltration_model(
    #     cation_list=["Co", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Co": 288, "Al": 20, "Cl": 636},
    #         "diafiltrate": {"Co": 3, "Al": 3, "Cl": 15},
    #     },
    #     non_Donnan_partition_dict={"Co": 0.3, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=False,
    # )
    # _test_diagnostics(model_CoAlCl5_no_boundary_layer.fs.unit)
    # _test_solve(model_CoAlCl5_no_boundary_layer)
    # _test_numerical_issues(model_CoAlCl5_no_boundary_layer.fs.unit)
    # CoAlCl5_no_boundary_layer_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(
    #     model_CoAlCl5_no_boundary_layer.fs.unit, CoAlCl5_no_boundary_layer_test_dict
    # )

    # model_LiCoAlCl6 = diafiltration_model(
    #     cation_list=["Li", "Co", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
    #         "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    #     },
    #     non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=True,
    # )
    # _test_diagnostics(model_LiCoAlCl6.fs.unit)
    # _test_solve(model_LiCoAlCl6)
    # _test_numerical_issues(model_LiCoAlCl6.fs.unit)
    # LiCoAlCl6_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(model_LiCoAlCl6.fs.unit, LiCoAlCl6_test_dict)

    # model_LiCoAlCl6_no_boundary_layer = diafiltration_model(
    #     cation_list=["Li", "Co", "Al"],
    #     anion_list=["Cl"],
    #     inlet_flow_volume={"feed": 12.5, "diafiltrate": 3.75},
    #     inlet_concentration={
    #         "feed": {"Li": 245, "Co": 288, "Al": 20, "Cl": 881},
    #         "diafiltrate": {"Li": 14, "Co": 3, "Al": 3, "Cl": 29},
    #     },
    #     non_Donnan_partition_dict={"Li": 0.7, "Co": 0.3, "Al": 0.0005, "Cl": 0.1},
    #     boundary_layer=False,
    # )
    # _test_diagnostics(model_LiCoAlCl6_no_boundary_layer.fs.unit)
    # _test_solve(model_LiCoAlCl6_no_boundary_layer)
    # _test_numerical_issues(model_LiCoAlCl6_no_boundary_layer.fs.unit)
    # LiCoAlCl6_no_boundary_layer_test_dict = {
    #     "applied_pressure": {(0): (, 1e-4, None)},
    #     "retentate_flow_volume": {(0, 1): (12.970, 1e-4, None)},
    #     "retentate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    #     "permeate_flow_volume": {(0, 1): (3.2800, 1e-4, None)},
    #     "permeate_conc_mol_comp": {
    #         (0, 1, "Li"): (, 1e-4, None),
    #         (0, 1, "Co"): (, 1e-4, None),
    #         (0, 1, "Al"): (, 1e-4, None),
    #         (0, 1, "Cl"): (, 1e-4, None),
    #     },
    # }
    # assert_solution_equivalent(
    #     model_LiCoAlCl6_no_boundary_layer.fs.unit, LiCoAlCl6_no_boundary_layer_test_dict
    # )
