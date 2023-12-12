#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for UKy flowsheet.

"""

import pytest
from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models.feed import Feed
from idaes.models.unit_models.mixer import Mixer
from idaes.models.unit_models.mscontactor import MSContactor
from idaes.models.unit_models.product import Product
from idaes.models.unit_models.separator import Separator
from idaes.models.unit_models.solid_liquid import SLSeparator
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from prommis.Solvent_Extraction.REEOgdistribution import REESolExOgParameters
from prommis.Solvent_Extraction.SolventExtraction import SolventExtraction
from prommis.UKy_flowsheet.UKy_flowsheet import (
    build,
    initialize_system,
    set_operating_conditions,
    set_scaling,
    solve,
)

# @pytest.fixture(scope="module")
# def model():
#     m = build()
#
#     return


class TestUKyFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m = build()
        set_operating_conditions(m)
        return m

    @pytest.mark.known_issue(6)
    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.known_issue(6)
    @pytest.mark.unit
    def test_build_flowsheet(self, model):
        assert isinstance(model.fs, FlowsheetBlock)

        # Leaching section property packages and unit models
        assert isinstance(model.fs.leach_soln, LeachSolutionParameters)
        assert isinstance(model.fs.coal, CoalRefuseParameters)
        assert isinstance(model.fs.leach_rxns, CoalRefuseLeachingReactions)

        assert isinstance(model.fs.leach, MSContactor)
        assert isinstance(model.fs.sl_sep1, SLSeparator)
        assert isinstance(model.fs.leach_mixer, Mixer)
        assert isinstance(model.fs.leach_liquid_feed, Feed)
        assert isinstance(model.fs.leach_solid_feed, Feed)
        assert isinstance(model.fs.leach_filter_cake, Product)
        assert isinstance(model.fs.leach_filter_cake_liquid, Product)

        # Solvent extraction section property packages and unit models
        assert isinstance(model.fs.prop_a, REESolExAqParameters)
        assert isinstance(model.fs.prop_o, REESolExOgParameters)

        assert isinstance(model.fs.solex_rougher, SolventExtraction)
        assert isinstance(model.fs.solex_cleaner, SolventExtraction)
        assert isinstance(model.fs.sep1, Separator)
        assert isinstance(model.fs.sx_mixer, Mixer)
        assert isinstance(model.fs.recycle1_purge, Product)

        # Precipitation property packages and unit models
        assert isinstance(model.fs.properties_aq, AqueousParameter)
        assert isinstance(model.fs.properties_solid, PrecipitateParameters)

        assert isinstance(model.fs.precipitator, Precipitator)
        assert isinstance(model.fs.sl_sep2, SLSeparator)
        assert isinstance(model.fs.sep2, Separator)
        assert isinstance(model.fs.recycle2_purge, Product)

        # Roasting property packages and unit models
        assert isinstance(model.fs.prop_gas, GenericParameterBlock)
        assert isinstance(model.fs.prop_solid, PrecipitateParameters)

        assert isinstance(model.fs.roaster, REEOxalateRoaster)

        # Flowsheet connections
        assert isinstance(model.fs.sol_feed, Arc)
        assert isinstance(model.fs.liq_feed, Arc)
        assert isinstance(model.fs.feed_mixture, Arc)
        assert isinstance(model.fs.s01, Arc)
        assert isinstance(model.fs.s02, Arc)
        assert isinstance(model.fs.sep1_solid, Arc)
        assert isinstance(model.fs.sep1_retained_liquid, Arc)
        assert isinstance(model.fs.sep1_liquid, Arc)
        assert isinstance(model.fs.recycle1, Arc)
        assert isinstance(model.fs.purge1, Arc)
        assert isinstance(model.fs.recycle_feed, Arc)
        assert isinstance(model.fs.s03, Arc)
        assert isinstance(model.fs.s04, Arc)
        assert isinstance(model.fs.s05, Arc)
        assert isinstance(model.fs.s06, Arc)
        assert isinstance(model.fs.s07, Arc)
        assert isinstance(model.fs.s08, Arc)
        assert isinstance(model.fs.sep2_solid, Arc)
        assert isinstance(model.fs.sep2_recovered_liquid, Arc)
        assert isinstance(model.fs.purge2, Arc)
        assert isinstance(model.fs.recycle2, Arc)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.known_issue(6)
    @pytest.mark.unit
    def test_set_dof(self, model):
        set_scaling(model)
        set_operating_conditions(model)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.known_issue(6)
    @pytest.mark.unit
    def test_initialize_flowsheet(self, model):
        initialize_system(model)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.integration
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.known_issue(6)
    @pytest.mark.unit
    def test_solve_flowsheet(self, model):
        solve(model)

        assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
            22.304218, 1e-4
        )
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Al2O3"
        ].value == pytest.approx(0.23304405, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "CaO"
        ].value == pytest.approx(0.0024095, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Ce2O3"
        ].value == pytest.approx(0.0001138787, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Dy2O3"
        ].value == pytest.approx(6.7167941e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Fe2O3"
        ].value == pytest.approx(0.0573268, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Gd2O3"
        ].value == pytest.approx(4.7948411e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "La2O3"
        ].value == pytest.approx(5.2969083e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Nd2O3"
        ].value == pytest.approx(4.6071439e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Pr2O3"
        ].value == pytest.approx(1.1902727e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sc2O3"
        ].value == pytest.approx(2.7237136e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sm2O3"
        ].value == pytest.approx(1.2870297e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Y2O3"
        ].value == pytest.approx(3.0385806e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "inerts"
        ].value == pytest.approx(0.7069127, 1e-4)

        assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
            606.6516829, 1e-4
        )
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(355.0736527, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(57.4481733, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(1.40988819, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.03071039, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(500.3169574, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(0.1845452, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(4.6227523, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(1508.1269700, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3303.975195, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.5048725, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.716497, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(0.173730, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.024634, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.0688498, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.087815, 1e-4)

        assert model.fs.solex_rougher.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.01, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(253.273371, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(42.0733366, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(9.6544128, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.21031067, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(211.353471, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.2637967, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(3.4047256, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(4.9067072, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1.1896853, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(321.508700, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.4714962, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.60137392, 1e-4)

        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(424.656178, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(4.6227523, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3303.975195, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(1508.12697, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(318.089658, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(51.3044558, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.000112096, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(3.07104e-11, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(469.4542750, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(5.0639192e-7, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.007700816, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(3.6684669e-7, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(7.4486762e-6, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(2.463428e-11, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(6.884975e-11, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(8.781502e-11, 1e-4)

        assert model.fs.solex_cleaner.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.010, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(297.115596, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(49.734607, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(12.0853316, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.228820, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(232.147087, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.367060, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(4.879478, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(5.552049, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1.381098, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(564.89968, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.512326, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.7166909, 1e-4)

        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(105.5845946, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(221.457043, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(37.5738603, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.01135198e-2, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(1.08708e-11, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(185.75940857, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.66415e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(1.3415618e-2, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(1.9405347e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(4.820097e-6, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(1.4294391e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(2.397942e-11, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(6.77259e-11, 1e-4)

        assert model.fs.precipitator.cv_aqueous.properties_out[
            0
        ].flow_vol.value == pytest.approx(167.5946, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Al"
        ].value == pytest.approx(247.20565109, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ca"
        ].value == pytest.approx(42.0733366, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ce"
        ].value == pytest.approx(1.42779599, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Dy"
        ].value == pytest.approx(1.08708e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Fe"
        ].value == pytest.approx(197.9715322, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Gd"
        ].value == pytest.approx(6.06468988e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "La"
        ].value == pytest.approx(87.9540105e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Nd"
        ].value == pytest.approx(37.9010675e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Pr"
        ].value == pytest.approx(11.2422085e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sc"
        ].value == pytest.approx(142.94391044, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sm"
        ].value == pytest.approx(2.39794e-2, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Y"
        ].value == pytest.approx(6.772586e-2, 1e-4)

        assert model.fs.precipitator.precipitate_outlet.temperature[
            0
        ].value == pytest.approx(348.15, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Al2(C2O4)3(s)"
        ].value == pytest.approx(6.972413e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Ce2(C2O4)3(s)"
        ].value == pytest.approx(1.820393e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Dy2(C2O4)3(s)"
        ].value == pytest.approx(0.0380531e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Fe2(C2O4)3(s)"
        ].value == pytest.approx(7.429618e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Gd2(C2O4)3(s)"
        ].value == pytest.approx(0.23722518e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "La2(C2O4)3(s)"
        ].value == pytest.approx(0.5636469e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Nd2(C2O4)3(s)"
        ].value == pytest.approx(0.973234e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Pr2(C2O4)3(s)"
        ].value == pytest.approx(0.2370397e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sc2(C2O4)3(s)"
        ].value == pytest.approx(123.178955e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sm2(C2O4)3(s)"
        ].value == pytest.approx(0.0922802e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Y2(C2O4)3(s)"
        ].value == pytest.approx(0.1861060e-3, 1e-4)

        assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
            9.021911e-3, 1e-4
        )
        assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
            873.15, 1e-4
        )
        assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(
            101325, 1e-4
        )
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "CO2"
        ].value == pytest.approx(0.0634927, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "H2O"
        ].value == pytest.approx(0.2016929, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "N2"
        ].value == pytest.approx(0.644578, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "O2"
        ].value == pytest.approx(0.090236, 1e-4)
