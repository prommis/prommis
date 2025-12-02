#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the two-salt diafiltration unit model.

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

from prommis.nanofiltration.diafiltration_solute_properties import (
    SoluteParameter,
)
from prommis.nanofiltration.diafiltration_two_salt import TwoSaltDiafiltration


@pytest.mark.unit
def test_zero_chi_implementation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteParameter()

    m.fs.unit = TwoSaltDiafiltration(
        property_package=m.fs.properties,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()
    m.fs.unit.membrane_fixed_charge.set_value(0)

    dt = DiagnosticsToolbox(m.fs.unit)
    dt.assert_no_structural_warnings()

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)
    solver = SolverFactory("ipopt")
    results = solver.solve(scaled_model, tee=True)
    scaling.propagate_solution(scaled_model, m)

    assert_optimal_termination(results)
    dt.assert_no_numerical_warnings()

    test_dict_zero_chi = {
        "retentate_final": [
            value(m.fs.unit.retentate_flow_volume[0, 1]),
            4.2918725888627725,
        ],
        "lithium_retentate_final": [
            value(m.fs.unit.retentate_conc_mol_comp[0, 1, "Li"]),
            192.59324298587342,
        ],
        "cobalt_retentate_final": [
            value(m.fs.unit.retentate_conc_mol_comp[0, 1, "Co"]),
            223.44760568527482,
        ],
        "chloride_retentate_final": [
            value(m.fs.unit.retentate_conc_mol_comp[0, 1, "Cl"]),
            639.4884543564231,
        ],
        "permeate_final": [
            value(m.fs.unit.permeate_flow_volume[0, 1]),
            11.957147777311302,
        ],
        "lithium_permeate_final": [
            value(m.fs.unit.permeate_conc_mol_comp[0, 1, "Li"]),
            191.96627751421204,
        ],
        "cobalt_permeate_final": [
            value(m.fs.unit.permeate_conc_mol_comp[0, 1, "Co"]),
            222.60051761723855,
        ],
        "chloride_permeate_final": [
            value(m.fs.unit.permeate_conc_mol_comp[0, 1, "Cl"]),
            637.1673127486891,
        ],
    }

    for model_result, test_val in test_dict_zero_chi.values():
        assert pytest.approx(test_val, rel=1e-5) == value(model_result)


@pytest.fixture(scope="module")
def diafiltration_two_salt():
    """
    Build a flowsheet with the two-salt diafiltration unit model.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteParameter()

    m.fs.unit = TwoSaltDiafiltration(
        property_package=m.fs.properties,
        NFE_module_length=10,
        NFE_membrane_thickness=5,
    )

    assert value(m.fs.unit.membrane_fixed_charge) == -140

    assert degrees_of_freedom(m.fs.unit) == 3

    m.fs.unit.total_module_length.fix()
    m.fs.unit.total_membrane_length.fix()
    m.fs.unit.applied_pressure.fix()

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config(diafiltration_two_salt):
    assert len(diafiltration_two_salt.fs.unit.config) == 6

    assert not diafiltration_two_salt.fs.unit.config.dynamic
    assert not diafiltration_two_salt.fs.unit.config.has_holdup

    assert (
        diafiltration_two_salt.fs.unit.config.property_package
        is diafiltration_two_salt.fs.properties
    )
    assert diafiltration_two_salt.fs.unit.config.NFE_module_length is 10
    assert diafiltration_two_salt.fs.unit.config.NFE_membrane_thickness is 5


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

        assert isinstance(diafiltration_two_salt.fs.unit.D_tilde, Var)
        assert len(diafiltration_two_salt.fs.unit.D_tilde) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_bilinear,
            Var,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_bilinear
            )
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_bilinear,
            Var,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_bilinear
            )
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_bilinear,
            Var,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_bilinear
            )
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_bilinear,
            Var,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_bilinear
            )
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_lithium_bilinear, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.convection_coefficient_lithium_bilinear)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_bilinear, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_bilinear)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt, Var
        )
        assert (
            len(diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt)
            == 66
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_lithium, Var
        )
        assert len(diafiltration_two_salt.fs.unit.convection_coefficient_lithium) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_cobalt, Var
        )
        assert len(diafiltration_two_salt.fs.unit.convection_coefficient_cobalt) == 66

        assert isinstance(diafiltration_two_salt.fs.unit.volume_flux_water, Var)
        assert len(diafiltration_two_salt.fs.unit.volume_flux_water) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.mol_flux_lithium, Var)
        assert len(diafiltration_two_salt.fs.unit.mol_flux_lithium) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.mol_flux_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.mol_flux_cobalt) == 11

        assert isinstance(diafiltration_two_salt.fs.unit.mol_flux_chloride, Var)
        assert len(diafiltration_two_salt.fs.unit.mol_flux_chloride) == 11

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

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_conc_mol_lithium, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mol_lithium) == 66

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_conc_mol_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mol_cobalt) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mol_chloride, Var
        )
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mol_chloride) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_lithium, Var
        )
        assert len(diafiltration_two_salt.fs.unit.convection_coefficient_lithium) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_cobalt, Var
        )
        assert len(diafiltration_two_salt.fs.unit.convection_coefficient_cobalt) == 66

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
            diafiltration_two_salt.fs.unit.d_membrane_conc_mol_lithium_dz,
            DerivativeVar,
        )
        assert len(diafiltration_two_salt.fs.unit.d_membrane_conc_mol_lithium_dz) == 66

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mol_cobalt_dz, DerivativeVar
        )
        assert len(diafiltration_two_salt.fs.unit.d_membrane_conc_mol_cobalt_dz) == 66

        # constraints
        assert isinstance(
            diafiltration_two_salt.fs.unit.overall_mol_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.overall_mol_balance) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_mol_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.lithium_mol_balance) == 10

        assert isinstance(diafiltration_two_salt.fs.unit.cobalt_mol_balance, Constraint)
        assert len(diafiltration_two_salt.fs.unit.cobalt_mol_balance) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.bulk_flux_equation_overall, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.bulk_flux_equation_overall) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.bulk_flux_equation_lithium, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.bulk_flux_equation_lithium) == 10

        assert isinstance(
            diafiltration_two_salt.fs.unit.bulk_flux_equation_cobalt, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.bulk_flux_equation_cobalt) == 10

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
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_permeate) == 11

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_membrane_interface_lithium)
            == 10
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_cobalt,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_membrane_interface_cobalt)
            == 10
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_permeate_interface_lithium)
            == 10
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_cobalt,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_permeate_interface_cobalt) == 10
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.D_tilde_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.D_tilde_calculation) == 60

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_lithium_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.convection_coefficient_lithium_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_bilinear_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_bilinear_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_lithium_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_lithium_cobalt_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_lithium_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.diffusion_coefficient_cobalt_cobalt_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_lithium_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.convection_coefficient_lithium_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_calculation,
            Constraint,
        )
        assert (
            len(
                diafiltration_two_salt.fs.unit.convection_coefficient_cobalt_calculation
            )
            == 60
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.lithium_flux_membrane) == 60

        assert isinstance(
            diafiltration_two_salt.fs.unit.cobalt_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.cobalt_flux_membrane) == 60

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_membrane) == 60

        for x in diafiltration_two_salt.fs.unit.dimensionless_module_length:
            assert diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[
                0, x, "Cl"
            ].fixed
            if x != 0:
                assert not diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx_disc_eq[
                    0, x, "Cl"
                ].active

        assert diafiltration_two_salt.fs.unit.retentate_flow_volume[0, 0].fixed
        assert value(diafiltration_two_salt.fs.unit.retentate_flow_volume[0, 0]) == (
            value(diafiltration_two_salt.fs.unit.feed_flow_volume[0])
            + value(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume[0])
        )

        assert diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 0, "Li"].fixed
        assert value(
            diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 0, "Li"]
        ) == (
            (
                value(diafiltration_two_salt.fs.unit.feed_flow_volume[0])
                * value(diafiltration_two_salt.fs.unit.feed_conc_mol_comp[0, "Li"])
                + value(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume[0])
                * value(
                    diafiltration_two_salt.fs.unit.diafiltrate_conc_mol_comp[0, "Li"]
                )
            )
            / (
                value(diafiltration_two_salt.fs.unit.feed_flow_volume[0])
                + value(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume[0])
            )
        )

        assert diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 0, "Co"].fixed
        assert value(
            diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 0, "Co"]
        ) == (
            (
                value(diafiltration_two_salt.fs.unit.feed_flow_volume[0])
                * value(diafiltration_two_salt.fs.unit.feed_conc_mol_comp[0, "Co"])
                + value(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume[0])
                * value(
                    diafiltration_two_salt.fs.unit.diafiltrate_conc_mol_comp[0, "Co"]
                )
            )
            / (
                value(diafiltration_two_salt.fs.unit.feed_flow_volume[0])
                + value(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume[0])
            )
        )

        assert diafiltration_two_salt.fs.unit.permeate_flow_volume[0, 0].fixed
        assert value(
            diafiltration_two_salt.fs.unit.permeate_flow_volume[0, 0]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 0, "Li"].fixed
        assert value(
            diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 0, "Li"]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 0, "Co"].fixed
        assert value(
            diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 0, "Co"]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        for z in diafiltration_two_salt.fs.unit.dimensionless_membrane_thickness:
            assert diafiltration_two_salt.fs.unit.membrane_conc_mol_lithium[0, z].fixed
            assert value(
                diafiltration_two_salt.fs.unit.membrane_conc_mol_lithium[0, z]
            ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

            assert diafiltration_two_salt.fs.unit.membrane_conc_mol_cobalt[0, z].fixed
            assert value(
                diafiltration_two_salt.fs.unit.membrane_conc_mol_cobalt[0, z]
            ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

            assert diafiltration_two_salt.fs.unit.membrane_conc_mol_chloride[0, z].fixed
            assert value(
                diafiltration_two_salt.fs.unit.membrane_conc_mol_chloride[0, z]
            ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[
            0, 0, "Li"
        ].fixed
        assert value(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[0, 0, "Li"]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[
            0, 0, "Co"
        ].fixed
        assert value(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mol_comp_dx[0, 0, "Co"]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx[0, 0].fixed
        assert value(
            diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx[0, 0]
        ) == value(diafiltration_two_salt.fs.unit.numerical_zero_tolerance)

        assert diafiltration_two_salt.fs.unit.volume_flux_water[0].fixed
        assert value(diafiltration_two_salt.fs.unit.volume_flux_water[0]) == value(
            diafiltration_two_salt.fs.unit.numerical_zero_tolerance
        )

        assert diafiltration_two_salt.fs.unit.mol_flux_lithium[0].fixed
        assert value(diafiltration_two_salt.fs.unit.mol_flux_lithium[0]) == value(
            diafiltration_two_salt.fs.unit.numerical_zero_tolerance
        )

        assert diafiltration_two_salt.fs.unit.mol_flux_cobalt[0].fixed
        assert value(diafiltration_two_salt.fs.unit.mol_flux_cobalt[0]) == value(
            diafiltration_two_salt.fs.unit.numerical_zero_tolerance
        )

        assert diafiltration_two_salt.fs.unit.mol_flux_chloride[0].fixed
        assert value(diafiltration_two_salt.fs.unit.mol_flux_chloride[0]) == value(
            diafiltration_two_salt.fs.unit.numerical_zero_tolerance
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
                4.247958521740219,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 1, "Li"]
                ),
                192.6044474143687,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 1, "Co"]
                ),
                223.2315850704493,
            ],
            "chloride_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mol_comp[0, 1, "Cl"]
                ),
                639.0676175552672,
            ],
            "permeate_final": [
                value(diafiltration_two_salt.fs.unit.permeate_flow_volume[0, 1]),
                12.001247151436175,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 1, "Li"]
                ),
                191.975023527502,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 1, "Co"]
                ),
                222.54094050056483,
            ],
            "chloride_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mol_comp[0, 1, "Cl"]
                ),
                637.0569045286317,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-5) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        dt.assert_no_numerical_warnings()
