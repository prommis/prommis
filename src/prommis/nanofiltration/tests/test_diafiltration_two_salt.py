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

    assert degrees_of_freedom(m.fs.unit) == 9

    m.fs.unit.membrane_width.fix(1)
    m.fs.unit.membrane_length.fix(100)
    m.fs.unit.applied_pressure.fix(10)

    m.fs.unit.feed_flow_volume.fix(100)
    m.fs.unit.feed_conc_mass_comp[0, "Li"].fix(1.7)
    m.fs.unit.feed_conc_mass_comp[0, "Co"].fix(17)
    m.fs.unit.feed_conc_mass_comp[0, "Cl"].fix(10.7)

    m.fs.unit.diafiltrate_flow_volume.fix(30)
    m.fs.unit.diafiltrate_conc_mass_comp[0, "Li"].fix(0.1)
    m.fs.unit.diafiltrate_conc_mass_comp[0, "Co"].fix(0.2)
    m.fs.unit.diafiltrate_conc_mass_comp[0, "Cl"].fix(0.2)

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
    assert diafiltration_two_salt.fs.unit.config.NFEx is 5
    assert diafiltration_two_salt.fs.unit.config.NFEz is 5


class TestDiafiltrationTwoSalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, diafiltration_two_salt):
        # parameters
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_thickness, Param)
        assert value(diafiltration_two_salt.fs.unit.membrane_thickness) == 1e-7

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_permeability, Param)
        assert value(diafiltration_two_salt.fs.unit.membrane_permeability) == 0.03

        assert isinstance(diafiltration_two_salt.fs.unit.temperature, Param)
        assert value(diafiltration_two_salt.fs.unit.temperature) == 298

        # sets
        assert isinstance(diafiltration_two_salt.fs.unit.x_bar, ContinuousSet)
        assert len(diafiltration_two_salt.fs.unit.x_bar) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.z_bar, ContinuousSet)
        assert len(diafiltration_two_salt.fs.unit.z_bar) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.time, Set)
        assert len(diafiltration_two_salt.fs.unit.time) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.solutes, Set)
        assert len(diafiltration_two_salt.fs.unit.solutes) == 3

        # variables
        assert isinstance(diafiltration_two_salt.fs.unit.membrane_width, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_width) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_length, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_length) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.applied_pressure, Var)
        assert len(diafiltration_two_salt.fs.unit.applied_pressure) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.feed_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.feed_flow_volume) == 1

        assert isinstance(diafiltration_two_salt.fs.unit.feed_conc_mass_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.feed_conc_mass_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_flow_volume) == 1

        assert isinstance(
            diafiltration_two_salt.fs.unit.diafiltrate_conc_mass_comp, Var
        )
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_conc_mass_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.volume_flux_water, Var)
        assert len(diafiltration_two_salt.fs.unit.volume_flux_water) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_lithium, Var)
        assert len(diafiltration_two_salt.fs.unit.mass_flux_lithium) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.mass_flux_cobalt) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.mass_flux_chlorine, Var)
        assert len(diafiltration_two_salt.fs.unit.mass_flux_chlorine) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.retentate_flow_volume) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_conc_mass_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.retentate_conc_mass_comp) == 18

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_flow_volume, Var)
        assert len(diafiltration_two_salt.fs.unit.permeate_flow_volume) == 6

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_conc_mass_comp, Var)
        assert len(diafiltration_two_salt.fs.unit.permeate_conc_mass_comp) == 18

        assert isinstance(diafiltration_two_salt.fs.unit.osmotic_pressure, Var)
        assert len(diafiltration_two_salt.fs.unit.osmotic_pressure) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mass_lithium, Var
        )
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mass_lithium) == 36

        assert isinstance(diafiltration_two_salt.fs.unit.membrane_conc_mass_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mass_cobalt) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_conc_mass_chlorine, Var
        )
        assert len(diafiltration_two_salt.fs.unit.membrane_conc_mass_chlorine) == 36

        assert isinstance(diafiltration_two_salt.fs.unit.D_lithium_lithium, Var)
        assert len(diafiltration_two_salt.fs.unit.D_lithium_lithium) == 36

        assert isinstance(diafiltration_two_salt.fs.unit.D_lithium_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.D_lithium_cobalt) == 36

        assert isinstance(diafiltration_two_salt.fs.unit.D_cobalt_lithium, Var)
        assert len(diafiltration_two_salt.fs.unit.D_cobalt_lithium) == 36

        assert isinstance(diafiltration_two_salt.fs.unit.D_cobalt_cobalt, Var)
        assert len(diafiltration_two_salt.fs.unit.D_cobalt_cobalt) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_conc_mass_comp_dx,
            DerivativeVar,
        )
        assert len(diafiltration_two_salt.fs.unit.d_retentate_conc_mass_comp_dx) == 18

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx, DerivativeVar
        )
        assert len(diafiltration_two_salt.fs.unit.d_retentate_flow_volume_dx) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mass_lithium_dz,
            DerivativeVar,
        )
        assert len(diafiltration_two_salt.fs.unit.d_membrane_conc_mass_lithium_dz) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.d_membrane_conc_mass_cobalt_dz, DerivativeVar
        )
        assert len(diafiltration_two_salt.fs.unit.d_membrane_conc_mass_cobalt_dz) == 36

        # constraints
        assert isinstance(
            diafiltration_two_salt.fs.unit.overall_mass_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.overall_mass_balance) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_mass_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.lithium_mass_balance) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.cobalt_mass_balance, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.cobalt_mass_balance) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_overall, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.geometric_flux_equation_overall) == 5

        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_lithium, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.geometric_flux_equation_lithium) == 5

        assert isinstance(
            diafiltration_two_salt.fs.unit.geometric_flux_equation_cobalt, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.geometric_flux_equation_cobalt) == 5

        assert isinstance(diafiltration_two_salt.fs.unit.lumped_water_flux, Constraint)
        assert len(diafiltration_two_salt.fs.unit.lumped_water_flux) == 5

        assert isinstance(
            diafiltration_two_salt.fs.unit.chlorine_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.chlorine_flux_membrane) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.osmotic_pressure_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.osmotic_pressure_calculation) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_retentate, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_retentate) == 6

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_membrane_interface_lithium)
            == 5
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_cobalt,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_membrane_interface_cobalt) == 5
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.retentate_membrane_interface_chlorine,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.retentate_membrane_interface_chlorine)
            == 5
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_permeate_interface_lithium) == 6
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_cobalt,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_permeate_interface_cobalt) == 6
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.membrane_permeate_interface_chlorine,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.membrane_permeate_interface_chlorine)
            == 6
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.D_lithium_lithium_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.D_lithium_lithium_calculation) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.D_lithium_cobalt_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.D_lithium_cobalt_calculation) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.D_cobalt_lithium_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.D_cobalt_lithium_calculation) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.D_cobalt_cobalt_calculation, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.D_cobalt_cobalt_calculation) == 36

        assert isinstance(
            diafiltration_two_salt.fs.unit.lithium_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.lithium_flux_membrane) == 30

        assert isinstance(
            diafiltration_two_salt.fs.unit.cobalt_flux_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.cobalt_flux_membrane) == 30

        assert isinstance(
            diafiltration_two_salt.fs.unit.electroneutrality_membrane, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.electroneutrality_membrane) == 30

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_flow_volume, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.initial_retentate_flow_volume) == 1

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_permeate_flow_volume, Constraint
        )
        assert len(diafiltration_two_salt.fs.unit.initial_permeate_flow_volume) == 1

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_lithium) == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_cobalt,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_retentate_conc_mass_cobalt) == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_lithium,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_membrane_interface_lithium) == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_cobalt, Constraint
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_membrane_interface_cobalt) == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_membrane_interface_chlorine,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_membrane_interface_chlorine) == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_lithium_dx,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_lithium_dx)
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_cobalt_dx,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit.initial_d_retentate_conc_mass_cobalt_dx)
            == 1
        )

        assert isinstance(
            diafiltration_two_salt.fs.unit._initial_d_retentate_flow_volume_dx,
            Constraint,
        )
        assert (
            len(diafiltration_two_salt.fs.unit._initial_d_retentate_flow_volume_dx) == 1
        )

        # ports
        assert isinstance(diafiltration_two_salt.fs.unit.feed_inlet, Port)
        assert len(diafiltration_two_salt.fs.unit.feed_inlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.feed_inlet.conc_mass_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.diafiltrate_inlet, Port)
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_inlet.flow_vol) == 1
        assert len(diafiltration_two_salt.fs.unit.diafiltrate_inlet.conc_mass_comp) == 3

        assert isinstance(diafiltration_two_salt.fs.unit.retentate_outlet, Port)
        assert len(diafiltration_two_salt.fs.unit.retentate_outlet.flow_vol) == 6
        assert len(diafiltration_two_salt.fs.unit.retentate_outlet.conc_mass_comp) == 18

        assert isinstance(diafiltration_two_salt.fs.unit.permeate_outlet, Port)
        assert len(diafiltration_two_salt.fs.unit.permeate_outlet.flow_vol) == 6
        assert len(diafiltration_two_salt.fs.unit.permeate_outlet.conc_mass_comp) == 18

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
                99.99999997000003,
            ],
            "lithium_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mass_comp[0, "Li", 1]
                ),
                1.330769230688393,
            ],
            "cobalt_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mass_comp[0, "Co", 1]
                ),
                13.12307692306405,
            ],
            "chlorine_retentate_final": [
                value(
                    diafiltration_two_salt.fs.unit.retentate_conc_mass_comp[0, "Cl", 1]
                ),
                22.58534905091955,
            ],
            "permeate_final": [
                value(diafiltration_two_salt.fs.unit.permeate_flow_volume[0, 1]),
                30.000000029999928,
            ],
            "lithium_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mass_comp[0, "Li", 1]
                ),
                1.330769231005046,
            ],
            "cobalt_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mass_comp[0, "Co", 1]
                ),
                13.123076923114663,
            ],
            "chlorine_permeate_final": [
                value(
                    diafiltration_two_salt.fs.unit.permeate_conc_mass_comp[0, "Cl", 1]
                ),
                22.585349052131857,
            ],
        }

        for model_result, test_val in test_dict.values():
            assert pytest.approx(test_val, rel=1e-5) == value(model_result)

    @pytest.mark.component
    def test_numerical_issues(self, diafiltration_two_salt):
        dt = DiagnosticsToolbox(diafiltration_two_salt.fs.unit)
        # TODO: resolve numerical warnings
        # some variables are hitting their lower bound of 0 (expected)
        # some residuals are large (unexpected)
        # dt.assert_no_numerical_warnings()
        dt.report_numerical_issues()
