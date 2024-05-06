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

from pyomo.environ import value
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models.feed import Feed
from idaes.models.unit_models.mixer import Mixer
from idaes.models.unit_models.product import Product
from idaes.models.unit_models.separator import Separator
from idaes.models.unit_models.solid_liquid import SLSeparator

import pytest

from prommis.leaching.leach_train import LeachingTrain
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
    add_costing,
    display_costing,
)


class TestUKyFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main()
        m.results = res
        return m

    @pytest.mark.component
    def test_structural_issues(self, system_frame):
        model = system_frame
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.unit
    def test_build_flowsheet(self, system_frame):
        model = system_frame
        assert isinstance(model.fs, FlowsheetBlock)

        # Leaching section property packages and unit models
        assert isinstance(model.fs.leach_soln, LeachSolutionParameters)
        assert isinstance(model.fs.coal, CoalRefuseParameters)
        assert isinstance(model.fs.leach_rxns, CoalRefuseLeachingReactions)

        assert isinstance(model.fs.leach, LeachingTrain)
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

    @pytest.mark.unit
    def test_solution(self, system_frame):
        model = system_frame
        assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
            22.234694, 1e-4
        )
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Al2O3"
        ].value == pytest.approx(0.23350023, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "CaO"
        ].value == pytest.approx(0.0017923, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Ce2O3"
        ].value == pytest.approx(9.590437e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Dy2O3"
        ].value == pytest.approx(6.2170773e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Fe2O3"
        ].value == pytest.approx(0.0553200, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Gd2O3"
        ].value == pytest.approx(3.6507596e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "La2O3"
        ].value == pytest.approx(4.50720462e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Nd2O3"
        ].value == pytest.approx(3.8214143e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Pr2O3"
        ].value == pytest.approx(9.4563816e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sc2O3"
        ].value == pytest.approx(2.6963267e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sm2O3"
        ].value == pytest.approx(1.1943966e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Y2O3"
        ].value == pytest.approx(2.91101636e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "inerts"
        ].value == pytest.approx(0.7091231, 1e-4)

        assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
            620.9470213, 1e-4
        )
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(422.37519528, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(109.542284, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(2.277073662, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.04719091, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(688.2668940, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(0.2584, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(1.7558201, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(693.459103, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3999.81895, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.9865420, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.94624, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(0.303450, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.032623, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.0971730, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.124519, 1e-4)

        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.01, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(1.267028e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(2.6847018e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.0003046844, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(8.008099e-6, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(2.8731444e-6, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(3.3573517e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.00010512, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.00016566, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(3.7119487e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(1.7348667, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(1.7011849e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(2.1799337e-5, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(41.44, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(3.940233e-8, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(1.6122889e-7, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(349.184092, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(739.885332, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(8396.893147, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(220.697679, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(79.18188, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(925.263005, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(2897.036916, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(4565.36775, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1022.9874915, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(36.682984, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(468.834838, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(600.77471, 1e-4)

        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.010, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(7.05737918e-9, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(1.53922935e-8, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(1.9068033e-6, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(4.35634416e-8, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(9.211715e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.8157884e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(7.402615e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(9.371941e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(2.154441e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.4212117, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(9.242277e-8, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(1.2989422e-7, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(9, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(0.389002, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(0.8484226, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(105.10298449, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(2.4012161, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(0.0507687, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(10.008624, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(40.803207, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(51.65813, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(11.8752795, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.1346149, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(5.094342, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(7.15976825, 1e-4)

        assert model.fs.precipitator.cv_aqueous.properties_out[
            0
        ].flow_vol.value == pytest.approx(9, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Al"
        ].value == pytest.approx(0.3855012, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ca"
        ].value == pytest.approx(0.84842261, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ce"
        ].value == pytest.approx(33.55938, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Dy"
        ].value == pytest.approx(0.30831, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Fe"
        ].value == pytest.approx(0.049530, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Gd"
        ].value == pytest.approx(1.20003, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "La"
        ].value == pytest.approx(19.78548, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Nd"
        ].value == pytest.approx(9.530926, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Pr"
        ].value == pytest.approx(2.61256, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sc"
        ].value == pytest.approx(0.0920631, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sm"
        ].value == pytest.approx(0.644434, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Y"
        ].value == pytest.approx(1.8286048, 1e-4)

        assert model.fs.precipitator.precipitate_outlet.temperature[
            0
        ].value == pytest.approx(348.15, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Al2(C2O4)3(s)"
        ].value == pytest.approx(5.838925e-7, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Ce2(C2O4)3(s)"
        ].value == pytest.approx(0.00229770, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Dy2(C2O4)3(s)"
        ].value == pytest.approx(5.7957229e-5, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Fe2(C2O4)3(s)"
        ].value == pytest.approx(9.981925e-8, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Gd2(C2O4)3(s)"
        ].value == pytest.approx(0.0002520741, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "La2(C2O4)3(s)"
        ].value == pytest.approx(6.80896e-4, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Nd2(C2O4)3(s)"
        ].value == pytest.approx(1.31427e-3, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Pr2(C2O4)3(s)"
        ].value == pytest.approx(2.958e-4, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sc2(C2O4)3(s)"
        ].value == pytest.approx(4.260289e-6, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sm2(C2O4)3(s)"
        ].value == pytest.approx(1.33178e-4, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Y2(C2O4)3(s)"
        ].value == pytest.approx(2.6984e-4, 1e-4)

        assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
            0.0085064, 1e-4
        )
        assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
            873.15, 1e-4
        )
        assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(
            101325, 1e-4
        )
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "CO2"
        ].value == pytest.approx(0.04061137, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "H2O"
        ].value == pytest.approx(0.173358, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "N2"
        ].value == pytest.approx(0.683643, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "O2"
        ].value == pytest.approx(0.102388, 1e-4)

        assert model.fs.leach_mixer.outlet.flow_vol[0].value == pytest.approx(
            620.76162, 1e-4
        )
        assert model.fs.rougher_org_make_up.outlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert model.fs.solex_rougher_load.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(440.422915, 1e-4)
        assert model.fs.solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(0.09, 1e-4)
        assert model.fs.acid_feed1.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert model.fs.acid_feed2.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert model.fs.acid_feed3.outlet.flow_vol[0].value == pytest.approx(9, 1e-4)
        assert model.fs.rougher_sep.inlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert model.fs.load_sep.inlet.flow_vol[0].value == pytest.approx(
            440.422915, 1e-4
        )
        assert model.fs.scrub_sep.inlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
        assert model.fs.rougher_mixer.outlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert model.fs.sc_circuit_purge.inlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert model.fs.solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(5.76, 1e-4)
        assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(9, 1e-4)
        assert model.fs.cleaner_org_make_up.outlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert model.fs.cleaner_mixer.outlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert model.fs.cleaner_sep.inlet.flow_vol[0].value == pytest.approx(
            62.01, 1e-4
        )
        assert model.fs.leach_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
            440.422915, 1e-4
        )
        assert model.fs.cleaner_purge.inlet.flow_vol[0].value == pytest.approx(
            6.201, 1e-4
        )
        assert model.fs.sl_sep1.recovered_liquid_outlet.flow_vol[
            0
        ].value == pytest.approx(434.662915, 1e-4)
        assert model.fs.sl_sep2.recovered_liquid_outlet.flow_vol[
            0
        ].value == pytest.approx(6.3, 1e-4)
        assert model.fs.precip_sep.inlet.flow_vol[0].value == pytest.approx(6.3, 1e-4)
        assert model.fs.precip_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
            5.76, 1e-4
        )
        assert model.fs.precip_purge.inlet.flow_vol[0].value == pytest.approx(
            0.63, 1e-4
        )

    @pytest.mark.unit
    def test_costing(self, system_frame):
        model = add_costing(system_frame)
        display_costing(model)

        assert value(model.fs.costing.total_plant_cost) == pytest.approx(15.576, rel=1e-4)
        assert value(model.fs.costing.total_BEC) == pytest.approx(5.244, rel=1e-4)
        assert value(model.fs.costing.total_installation_cost) == pytest.approx(
            10.33, rel=1e-4
        )
        assert value(model.fs.costing.other_plant_costs) == pytest.approx(
            0.0016309, rel=1e-4
        )
        assert value(model.fs.costing.total_fixed_OM_cost) == pytest.approx(
            7.2484, rel=1e-4
        )
        assert value(model.fs.costing.total_variable_OM_cost[0]) == pytest.approx(
            10.7422, rel=1e-4
        )
        assert value(model.fs.costing.land_cost) == pytest.approx(
            6.1234e-5, rel=1e-4
        )  # Expression, not Var
        assert value(model.fs.costing.total_sales_revenue) == pytest.approx(
            0.00018136, rel=1e-4
        )
        assert value(model.fs.L_pe_tanks.costing.bare_erected_cost["4.2"]) == pytest.approx(
            0.0003690, rel=1e-4
        )
        assert value(model.fs.L_tank_mixers.costing.bare_erected_cost["4.3"]) == pytest.approx(
            0.08625, rel=1e-4
        )
        assert value(model.fs.L_pump.costing.bare_erected_cost["4.4"]) == pytest.approx(
            0.01446, rel=1e-4
        )
        assert value(model.fs.L_thickener.costing.bare_erected_cost["4.5"]) == pytest.approx(
            0.0243540, rel=1e-4
        )
        assert value(model.fs.L_filter_press.costing.bare_erected_cost["4.6"]) == pytest.approx(
            0.08919, rel=1e-4
        )
        assert value(model.fs.L_solution_heater.costing.bare_erected_cost["4.8"]) == pytest.approx(
            0.00750107, rel=1e-4
        )
        assert value(model.fs.RSX_pe_tanks.costing.bare_erected_cost["5.1"]) == pytest.approx(
            0.000273433, rel=1e-4
        )
        assert value(model.fs.RSX_tank_mixers.costing.bare_erected_cost["5.2"]) == pytest.approx(
            0.03558, rel=1e-4
        )
        assert value(model.fs.RSX_pump.costing.bare_erected_cost["5.3"]) == pytest.approx(
            0.00623416, rel=1e-4
        )
        assert value(model.fs.RSX_mixer_settler.costing.bare_erected_cost["5.4"]) == pytest.approx(
            0.3876, rel=1e-4
        )
        assert value(model.fs.CSX_pe_tanks.costing.bare_erected_cost["6.1"]) == pytest.approx(
            9.54821e-5, rel=1e-4
        )
        assert value(model.fs.CSX_tank_mixers.costing.bare_erected_cost["6.2"]) == pytest.approx(
            0.0057947, rel=1e-4
        )
        assert value(model.fs.CSX_pump.costing.bare_erected_cost["6.3"]) == pytest.approx(
            0.00357727, rel=1e-4
        )
        assert value(model.fs.CSX_mixer_settler.costing.bare_erected_cost["6.4"]) == pytest.approx(
            0.25796, rel=1e-4
        )
        assert value(model.fs.reep_pe_tanks.costing.bare_erected_cost["9.2"]) == pytest.approx(
            1.453344e-5, rel=1e-4
        )
        assert value(model.fs.reep_tank_mixers.costing.bare_erected_cost["9.3"]) == pytest.approx(
            0.0090271, rel=1e-4
        )
        assert value(model.fs.reep_pump.costing.bare_erected_cost["9.4"]) == pytest.approx(
            0.0014137, rel=1e-4
        )
        assert value(model.fs.reep_filter_press.costing.bare_erected_cost["9.5"]) == pytest.approx(
            0.00356786, rel=1e-4
        )
        assert value(model.fs.reep_roaster.costing.bare_erected_cost["9.8"]) == pytest.approx(
            0.0877, rel=1e-4
        )
        assert value(model.fs.R_storage_bins.costing.bare_erected_cost["3.1"]) == pytest.approx(
            0.04151, rel=1e-4
        )
        assert value(model.fs.R_conveyors.costing.bare_erected_cost["3.2"]) == pytest.approx(
            1.7062e-6, rel=1e-4
        )
        assert value(model.fs.R_roaster.costing.bare_erected_cost["3.3"]) == pytest.approx(
            3.4557, rel=1e-4
        )
        assert value(model.fs.R_gas_scrubber.costing.bare_erected_cost["3.4"]) == pytest.approx(
            7.4055e-5, rel=1e-4
        )
        assert value(model.fs.R_spray_chamber_quencher.costing.bare_erected_cost["3.5"]) == pytest.approx(
            0.21234, rel=1e-4
        )
        assert value(model.fs.R_chiller.costing.bare_erected_cost["3.7"]) == pytest.approx(
            0.51352, rel=1e-4
        )
