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
    Param,
    SolverFactory,
    Suffix,
    TransformationFactory,
    assert_optimal_termination,
    units,
)
from pyomo.util.check_units import assert_units_consistent

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

    # Define global parameters
    m.membrane_thickness = Param(
        initialize=1e-7,
        units=units.m,
        doc="Thickness of membrane (z-direction)",
    )
    m.membrane_width = Param(
        initialize=1,
        units=units.m,
        doc="Width of the membrane (x-direction)",
    )
    m.membrane_length = Param(
        initialize=100,
        units=units.m,
        doc="Length of the membrane, wound radially",
    )
    m.dP = Param(
        initialize=10,
        units=units.bar,
        doc="Pressure applied to membrane",
    )
    m.Lp = Param(
        initialize=0.01,
        units=units.m / units.h / units.bar,
        doc="Hydraulic permeability coefficient",
    )
    m.feed_flow_volume = Param(
        initialize=100,
        units=units.m**3 / units.h,
        doc="Volumetric flow rate of the feed",
    )
    m.feed_conc_mass_lithium = Param(
        initialize=1.7,
        units=units.kg / units.m**3,
        doc="Mass concentration of lithium in the feed",
    )
    m.feed_conc_mass_cobalt = Param(
        initialize=17,
        units=units.kg / units.m**3,
        doc="Mass concentration of cobalt in the feed",
    )
    m.diafiltrate_flow_volume = Param(
        initialize=30,
        units=units.m**3 / units.h,
        doc="Volumetric flow rate of the diafiltrate",
    )
    m.diafiltrate_conc_mass_lithium = Param(
        initialize=0.1,
        units=units.kg / units.m**3,
        doc="Mass concentration of lithium in the diafiltrate",
    )
    m.diafiltrate_conc_mass_cobalt = Param(
        initialize=0.2,
        units=units.kg / units.m**3,
        doc="Mass concentration of cobalt in the diafiltrate",
    )

    m.fs.unit = TwoSaltDiafiltration(
        property_package=m.fs.properties,
        membrane_length=m.membrane_length,
        membrane_width=m.membrane_width,
        membrane_thickness=m.membrane_thickness,
        membrane_permeability=m.Lp,
        applied_pressure=m.dP,
        feed_flow_volume=m.feed_flow_volume,
        feed_conc_mass_lithium=m.feed_conc_mass_lithium,
        feed_conc_mass_cobalt=m.feed_conc_mass_cobalt,
        diafiltrate_flow_volume=m.diafiltrate_flow_volume,
        diafiltrate_conc_mass_lithium=m.diafiltrate_conc_mass_lithium,
        diafiltrate_conc_mass_cobalt=m.diafiltrate_conc_mass_cobalt,
        NFEx=5,
        NFEz=5,
    )

    return m


@pytest.mark.unit
def test_config(diafiltration_two_salt):
    assert len(diafiltration_two_salt.fs.unit.config) == 17

    assert not diafiltration_two_salt.fs.unit.config.dynamic
    assert not diafiltration_two_salt.fs.unit.config.has_holdup

    # TODO: add additional relevant assertions


class TestDiafiltrationTwoSalt(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, diafiltration_two_salt):
        assert hasattr(diafiltration_two_salt.fs.unit, "volume_flux_water")
        # TODO: add additional relevant assertions

    @pytest.mark.component
    def test_diagnostics(self, diafiltration_two_salt):
        assert_units_consistent(diafiltration_two_salt.fs.unit)

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
