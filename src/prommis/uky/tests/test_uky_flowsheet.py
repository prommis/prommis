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

from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

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

import pytest

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.uky.uky_flowsheet import (
    main,
    build,
    initialize_system,
    set_operating_conditions,
    set_scaling,
    solve,
)


class TestUKyFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m = build()
        set_operating_conditions(m)
        return m

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

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
        assert isinstance(model.fs.prop_o, REESolExOgParameters)

        assert isinstance(model.fs.solex_rougher_load, SolventExtraction)
        assert isinstance(model.fs.solex_rougher_scrub, SolventExtraction)
        assert isinstance(model.fs.solex_rougher_strip, SolventExtraction)
        assert isinstance(model.fs.solex_cleaner_load, SolventExtraction)
        assert isinstance(model.fs.solex_cleaner_strip, SolventExtraction)
        assert isinstance(model.fs.rougher_sep, Separator)
        assert isinstance(model.fs.load_sep, Separator)
        assert isinstance(model.fs.scrub_sep, Separator)
        assert isinstance(model.fs.cleaner_sep, Separator)
        assert isinstance(model.fs.rougher_mixer, Mixer)
        assert isinstance(model.fs.cleaner_mixer, Mixer)
        assert isinstance(model.fs.leach_sx_mixer, Mixer)
        assert isinstance(model.fs.sc_circuit_purge, Product)
        assert isinstance(model.fs.cleaner_purge, Product)
        assert isinstance(model.fs.cleaner_org_make_up, Feed)
        assert isinstance(model.fs.rougher_org_make_up, Feed)
        assert isinstance(model.fs.acid_feed1, Feed)
        assert isinstance(model.fs.acid_feed2, Feed)
        assert isinstance(model.fs.acid_feed3, Feed)

        # Precipitation property packages and unit models
        assert isinstance(model.fs.properties_aq, AqueousParameter)
        assert isinstance(model.fs.properties_solid, PrecipitateParameters)

        assert isinstance(model.fs.precipitator, Precipitator)
        assert isinstance(model.fs.sl_sep2, SLSeparator)
        assert isinstance(model.fs.precip_sx_mixer, Mixer)
        assert isinstance(model.fs.precip_purge, Product)

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
        assert isinstance(model.fs.mixed_aq_feed, Arc)
        assert isinstance(model.fs.org_feed, Arc)
        assert isinstance(model.fs.mixed_org_feed, Arc)
        assert isinstance(model.fs.s03, Arc)
        assert isinstance(model.fs.load_recycle, Arc)
        assert isinstance(model.fs.s04, Arc)
        assert isinstance(model.fs.s05, Arc)
        assert isinstance(model.fs.s06, Arc)
        assert isinstance(model.fs.scrub_recycle, Arc)
        assert isinstance(model.fs.s07, Arc)
        assert isinstance(model.fs.s08, Arc)
        assert isinstance(model.fs.s09, Arc)
        assert isinstance(model.fs.s10, Arc)
        assert isinstance(model.fs.s11, Arc)
        assert isinstance(model.fs.s12, Arc)
        assert isinstance(model.fs.s13, Arc)
        assert isinstance(model.fs.org_feed2, Arc)
        assert isinstance(model.fs.s14, Arc)
        assert isinstance(model.fs.s15, Arc)
        assert isinstance(model.fs.s16, Arc)
        assert isinstance(model.fs.s17, Arc)
        assert isinstance(model.fs.s18, Arc)
        assert isinstance(model.fs.s19, Arc)
        assert isinstance(model.fs.s20, Arc)
        assert isinstance(model.fs.s21, Arc)
        assert isinstance(model.fs.s22, Arc)
        assert isinstance(model.fs.s23, Arc)
        assert isinstance(model.fs.sep2_solid, Arc)
        assert isinstance(model.fs.sep2_recovered_liquid, Arc)
        assert isinstance(model.fs.s24, Arc)
        assert isinstance(model.fs.s25, Arc)

        assert degrees_of_freedom(model) == 0

    @pytest.fixture(scope="class")
    def scaled_model(self, model):
        scaled_model = set_scaling(model)

        return scaled_model

    @pytest.mark.unit
    def test_initialize_flowsheet(self, scaled_model):
        initialize_system(scaled_model)

        assert degrees_of_freedom(scaled_model) == 0

    @pytest.mark.integration
    def test_unit_consistency(self, scaled_model):
        assert_units_consistent(scaled_model)

    @pytest.mark.unit
    def test_solve_flowsheet(self, scaled_model):
        solve(scaled_model)

        assert scaled_model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
            22.234694, 1e-4
        )
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Al2O3"
        ].value == pytest.approx(0.23350023, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "CaO"
        ].value == pytest.approx(0.0017923, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Ce2O3"
        ].value == pytest.approx(9.590437e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Dy2O3"
        ].value == pytest.approx(6.2170773e-6, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Fe2O3"
        ].value == pytest.approx(0.0553200, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Gd2O3"
        ].value == pytest.approx(3.6507596e-6, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "La2O3"
        ].value == pytest.approx(4.50720462e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Nd2O3"
        ].value == pytest.approx(3.8214143e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Pr2O3"
        ].value == pytest.approx(9.4563816e-6, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sc2O3"
        ].value == pytest.approx(2.6963267e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sm2O3"
        ].value == pytest.approx(1.1943966e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Y2O3"
        ].value == pytest.approx(2.91101636e-5, 1e-4)
        assert scaled_model.fs.leach.solid_outlet.mass_frac_comp[
            0, "inerts"
        ].value == pytest.approx(0.7091231, 1e-4)

        assert scaled_model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
            620.9470213, 1e-4
        )
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(422.37519528, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(109.542284, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(2.277073662, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.04719091, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(688.2668940, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(0.2584, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(1.7558201, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(693.459103, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3999.81895, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.9865420, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.94624, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(0.303450, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.032623, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.0971730, 1e-4)
        assert scaled_model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.124519, 1e-4)

        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.01, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(1.267028e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(2.6847018e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.0003046844, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(8.008099e-6, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(2.8731444e-6, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(3.3573517e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.00010512, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.00016566, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(3.7119487e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(1.7348667, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(1.7011849e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(2.1799337e-5, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(41.44, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3.940233e-8, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(1.6122889e-7, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(349.184092, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(739.885332, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(8396.893147, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(220.697679, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(79.18188, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(925.263005, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(2897.036916, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(4565.36775, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1022.9874915, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(36.682984, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(468.834838, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(600.77471, 1e-4)

        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.010, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(7.05737918e-9, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(1.53922935e-8, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(1.9068033e-6, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(4.35634416e-8, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(9.211715e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.8157884e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(7.402615e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(9.371941e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(2.154441e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.4212117, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(9.242277e-8, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(1.2989422e-7, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(0.389002, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(0.8484226, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(105.10298449, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(2.4012161, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(0.0507687, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(10.008624, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(40.803207, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(51.65813, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(11.8752795, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.1346149, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(5.094342, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(7.15976825, 1e-4)

        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[
            0
        ].flow_vol.value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Al"
        ].value == pytest.approx(38.5501172, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ca"
        ].value == pytest.approx(84.84226098, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ce"
        ].value == pytest.approx(3355.93828, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Dy"
        ].value == pytest.approx(30.83161, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Fe"
        ].value == pytest.approx(0.049530, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Gd"
        ].value == pytest.approx(1.20003, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "La"
        ].value == pytest.approx(19.7854751, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Nd"
        ].value == pytest.approx(9.530926, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Pr"
        ].value == pytest.approx(2.61256, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sc"
        ].value == pytest.approx(0.0920631, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sm"
        ].value == pytest.approx(0.644434, 1e-4)
        assert scaled_model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Y"
        ].value == pytest.approx(1.8286048, 1e-4)

        assert scaled_model.fs.precipitator.precipitate_outlet.temperature[
            0
        ].value == pytest.approx(348.15, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Al2(C2O4)3(s)"
        ].value == pytest.approx(5.838925e-7, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Ce2(C2O4)3(s)"
        ].value == pytest.approx(0.00229770, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Dy2(C2O4)3(s)"
        ].value == pytest.approx(5.7957229e-5, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Fe2(C2O4)3(s)"
        ].value == pytest.approx(9.981925e-8, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Gd2(C2O4)3(s)"
        ].value == pytest.approx(0.0002520741, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "La2(C2O4)3(s)"
        ].value == pytest.approx(0.000680896, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Nd2(C2O4)3(s)"
        ].value == pytest.approx(0.00131427, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Pr2(C2O4)3(s)"
        ].value == pytest.approx(0.0002958, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sc2(C2O4)3(s)"
        ].value == pytest.approx(4.260289e-6, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sm2(C2O4)3(s)"
        ].value == pytest.approx(0.000133178, 1e-4)
        assert scaled_model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Y2(C2O4)3(s)"
        ].value == pytest.approx(0.00026984, 1e-4)

        assert scaled_model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
            0.00850637, 1e-4
        )
        assert scaled_model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
            873.15, 1e-4
        )
        assert scaled_model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(
            101325, 1e-4
        )
        assert scaled_model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "CO2"
        ].value == pytest.approx(0.04061137, 1e-4)
        assert scaled_model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "H2O"
        ].value == pytest.approx(0.173358, 1e-4)
        assert scaled_model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "N2"
        ].value == pytest.approx(0.683643, 1e-4)
        assert scaled_model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "O2"
        ].value == pytest.approx(0.102388, 1e-4)

        assert scaled_model.fs.leach_mixer.outlet.flow_vol[0].value == pytest.approx(
            620.76162, 1e-4
        )
        assert scaled_model.fs.rougher_org_make_up.outlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert scaled_model.fs.solex_rougher_load.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(440.422915, 1e-4)
        assert scaled_model.fs.solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.acid_feed1.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.acid_feed2.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.acid_feed3.outlet.flow_vol[0].value == pytest.approx(9, 1e-4)
        assert scaled_model.fs.rougher_sep.inlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert scaled_model.fs.load_sep.inlet.flow_vol[0].value == pytest.approx(
            440.422915, 1e-4
        )
        assert scaled_model.fs.scrub_sep.inlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert scaled_model.fs.rougher_mixer.outlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert scaled_model.fs.sc_circuit_purge.inlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert scaled_model.fs.solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(5.76, 1e-4)
        assert scaled_model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(9.0, 1e-4)
        assert scaled_model.fs.cleaner_org_make_up.outlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert scaled_model.fs.cleaner_mixer.outlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert scaled_model.fs.cleaner_sep.inlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert scaled_model.fs.leach_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
            440.422915, 1e-4
        )
        assert scaled_model.fs.cleaner_purge.inlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert scaled_model.fs.sl_sep1.recovered_liquid_outlet.flow_vol[
            0
        ].value == pytest.approx(434.662915, 1e-4)
        assert scaled_model.fs.sl_sep2.recovered_liquid_outlet.flow_vol[
            0
        ].value == pytest.approx(6.3, 1e-4)
        assert scaled_model.fs.precip_sep.inlet.flow_vol[0].value == pytest.approx(6.3, 1e-4)
        assert scaled_model.fs.precip_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
            5.76, 1e-4
        )
        assert scaled_model.fs.precip_purge.inlet.flow_vol[0].value == pytest.approx(
            0.63, 1e-4
        )

    def test_main(self):
        main()
