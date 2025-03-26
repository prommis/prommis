#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the two-salt diafiltration unit model.
"""

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Param,
    SolverFactory,
    Suffix,
    TransformationFactory,
    Var,
    assert_optimal_termination,
)
from pyomo.dae import ContinuousSet, DerivativeVar

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

import pytest

from prommis.nanofiltration.diafiltration_solute_properties import SoluteParameter
from prommis.nanofiltration.diafiltration_two_salt import TwoSaltDiafiltration


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
        NFEx=5,
        NFEz=5,
    )

    return m


@pytest.mark.unit
def test_config(diafiltration_two_salt):
    assert len(diafiltration_two_salt.fs.unit.config) == 6

    assert not diafiltration_two_salt.fs.unit.config.dynamic
    assert not diafiltration_two_salt.fs.unit.config.has_holdup

    # TODO: add additional relevant assertions


class TestDiafiltrationTwoSalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, diafiltration_two_salt):
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_thickness, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_width, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_length, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.applied_pressure, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_permeability, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.feed_flow_volume, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.feed_conc_mass_lithium, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.feed_conc_mass_cobalt, Param)
        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume, Param)
        assert isinstance(
            diafiltration_two_salt.fs.unit.diafiltrate_conc_mass_lithium, Param
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.diafiltrate_conc_mass_cobalt, Param
        )
        assert isinstance(diafiltration_two_salt.fs.unit.temperature, Param)

        assert isinstance(diafiltration_two_salt.fs.unit.x_bar, ContinuousSet)
        assert isinstance(diafiltration_two_salt.fs.unit.z_bar, ContinuousSet)

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_flow_volume, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_permeate_flow_volume, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_lithium,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_cobalt,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_lithium,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_cobalt, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_chlorine,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_permeate_conc_mass_lithium,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_permeate_conc_mass_cobalt, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_lithium_dx,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_cobalt_dx,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit._initial_d_retentate_flow_volume_dx,
            Constraint,
        )

        assert isinstance(diafiltration_two_salt.fs.unit.volume_flux_water, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_lithium, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_cobalt, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_chlorine, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.retentate_flow_volume, Var)
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_conc_mass_lithium, Var
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_conc_mass_cobalt, Var
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_conc_mass_chlorine, Var
        )
        assert isinstance(diafiltration_two_salt.fs.unit.permeate_flow_volume, Var)
        assert isinstance(
            diafiltration_two_salt.fs.unit.permeate_conc_mass_lithium, Var
        )
        assert isinstance(diafiltration_two_salt.fs.unit.permeate_conc_mass_cobalt, Var)
        assert isinstance(
            diafiltration_two_salt.fs.unit.permeate_conc_mass_chlorine, Var
        )
        assert isinstance(diafiltration_two_salt.fs.unit.osmotic_pressure, Var)
        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mass_lithium, Var
        )
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_conc_mass_cobalt, Var)
        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mass_chlorine, Var
        )
        assert isinstance(diafiltration_two_salt.fs.unit.D_lithium_lithium, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.D_lithium_cobalt, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.D_cobalt_lithium, Var)
        assert isinstance(diafiltration_two_salt.fs.unit.D_cobalt_cobalt, Var)

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mass_lithium_dx,
            DerivativeVar,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mass_cobalt_dx,
            DerivativeVar,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mass_lithium_dz,
            DerivativeVar,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mass_cobalt_dz, DerivativeVar
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.overall_mass_balance, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_mass_balance, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.cobalt_mass_balance, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.general_mass_balance_lithium, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.general_mass_balance_cobalt, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_overall, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_lithium, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_cobalt, Constraint
        )
        assert isinstance(diafiltration_two_salt.fs.unit.lumped_water_flux, Constraint)
        assert isinstance(
            diafiltration_two_salt.fs.unit.chlorine_flux_membrane, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.osmotic_pressure_calculation, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_retentate, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_lithium,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_cobalt,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_chlorine,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_lithium,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_cobalt,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_chlorine,
            Constraint,
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.D_lithium_lithium_calculation, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.D_lithium_cobalt_calculation, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.D_cobalt_lithium_calculation, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.D_cobalt_cobalt_calculation, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_flux_membrane, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.cobalt_flux_membrane, Constraint
        )
        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_membrane, Constraint
        )

        # TODO: add additional relevant assertions

    @pytest.mark.component
    def test_diagnostics(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve(self, diafiltration_two_salt):
        diafiltration_two_salt.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add scaling factors for poorly scaled variables
        for x in diafiltration_two_salt.fs.unit.x_bar:
            diafiltration_two_salt.scaling_factor[
                diafiltration_two_salt.fs.unit.retentate_flow_volume[x]
            ] = 1e-2
            diafiltration_two_salt.scaling_factor[
                diafiltration_two_salt.fs.unit.retentate_conc_mass_cobalt[x]
            ] = 1e-1
            diafiltration_two_salt.scaling_factor[
                diafiltration_two_salt.fs.unit.retentate_conc_mass_chlorine[x]
            ] = 1e-1
            for z in diafiltration_two_salt.fs.unit.z_bar:
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_lithium_lithium[x, z]
                ] = 1e8
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_lithium_cobalt[x, z]
                ] = 1e8
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_cobalt_lithium[x, z]
                ] = 1e8
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_cobalt_cobalt[x, z]
                ] = 1e8

                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.volume_flux_water[x]
                ] = 1e2
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.mass_flux_lithium[x]
                ] = 1e2
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.mass_flux_cobalt[x]
                ] = 1e2
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.mass_flux_chlorine[x]
                ] = 1e2

        # Add scaling factors for poorly scaled constraints
        for x in diafiltration_two_salt.fs.unit.x_bar:
            for z in diafiltration_two_salt.fs.unit.z_bar:
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_lithium_lithium_calculation[x, z]
                ] = 1e12
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_lithium_cobalt_calculation[x, z]
                ] = 1e12
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_cobalt_lithium_calculation[x, z]
                ] = 1e12
                diafiltration_two_salt.scaling_factor[
                    diafiltration_two_salt.fs.unit.D_cobalt_cobalt_calculation[x, z]
                ] = 1e12

                if z != 0:
                    diafiltration_two_salt.scaling_factor[
                        diafiltration_two_salt.fs.unit.lithium_flux_membrane[x, z]
                    ] = 1e8
                    diafiltration_two_salt.scaling_factor[
                        diafiltration_two_salt.fs.unit.cobalt_flux_membrane[x, z]
                    ] = 1e8

        scaling = TransformationFactory("core.scale_model")
        scaled_model = scaling.create_using(diafiltration_two_salt, rename=False)

        solver = SolverFactory("ipopt")
        results = solver.solve(scaled_model, tee=True)

        scaling.propagate_solution(scaled_model, diafiltration_two_salt)

        assert_optimal_termination(results)

    @pytest.mark.component
    def test_numerical_issues(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        dt.report_numerical_issues()  # some vars are purposely at their lower bound (0)

    # TODO: add test for the solution
