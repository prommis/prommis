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
from pyomo.environ import value
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
from prommis.solvent_extraction.ree_aq_distribution import REESolExAqParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.uky.uky_flowsheet import (
    build,
    initialize_system,
    set_operating_conditions,
    set_scaling,
    solve,
    add_costing,
    display_costing,
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

    @pytest.mark.unit
    def test_set_dof(self, model):
        set_scaling(model)
        set_operating_conditions(model)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_initialize_flowsheet(self, model):
        initialize_system(model)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.integration
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_solve_flowsheet(self, model):
        solve(model)

        assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
            22.237084, 1e-4
        )
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Al2O3"
        ].value == pytest.approx(0.23359591, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "CaO"
        ].value == pytest.approx(0.0017714, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Ce2O3"
        ].value == pytest.approx(9.557201e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Dy2O3"
        ].value == pytest.approx(6.2014028e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Fe2O3"
        ].value == pytest.approx(0.0553200, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Gd2O3"
        ].value == pytest.approx(3.6701135e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "La2O3"
        ].value == pytest.approx(4.48588583e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Nd2O3"
        ].value == pytest.approx(3.8120246e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Pr2O3"
        ].value == pytest.approx(9.4020634e-6, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sc2O3"
        ].value == pytest.approx(2.6956793e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Sm2O3"
        ].value == pytest.approx(1.1922270e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "Y2O3"
        ].value == pytest.approx(2.90770689e-5, 1e-4)
        assert model.fs.leach.solid_outlet.mass_frac_comp[
            0, "inerts"
        ].value == pytest.approx(0.7090469, 1e-4)

        assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
            606.7136886, 1e-4
        )
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(361.79296237, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(96.0984875, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(1.993387898, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.04781406, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(636.9422878, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(0.2207516, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(1.6990424, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(690.045112, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(4113.133015, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.7657892, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(0.970650, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(0.253155, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(0.033295, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.1000351, 1e-4)
        assert model.fs.leach.liquid_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.128225, 1e-4)

        assert model.fs.solex_rougher.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.01, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(258.092619, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(70.3868727, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(13.651406, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.32747326, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(269.096862, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.5118992, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(5.1648061, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(6.6478682, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1.7337573, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(321.5680337, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.6851294, 1e-4)
        assert model.fs.solex_rougher.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(0.87820016, 1e-4)

        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(424.699582, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H2O"
        ].value == pytest.approx(1000000.0, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "H"
        ].value == pytest.approx(1.6990424, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "SO4"
        ].value == pytest.approx(4113.133015, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "HSO4"
        ].value == pytest.approx(690.0451119, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(324.109094, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(85.8213642, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.000158488, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(4.78215e-11, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(597.6516996, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(6.05742461e-7, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.011680577, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(4.96972640e-7, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(1.0854028e-5, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(3.330248e-11, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(1.0004260e-10, 1e-4)
        assert model.fs.solex_rougher.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(1.2823268e-10, 1e-4)

        assert model.fs.solex_cleaner.mscontactor.organic_outlet.flow_vol[
            0
        ].value == pytest.approx(62.010, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(302.769067, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(83.203847, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(17.0887422, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(0.35629, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(295.571453, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.63543, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(7.401935, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(7.522212, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(2.012708, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(565.00393, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(0.744459, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.organic_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(1.0466002, 1e-4)

        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.flow_vol[
            0
        ].value == pytest.approx(105.5845946, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Al"
        ].value == pytest.approx(225.670894, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ca"
        ].value == pytest.approx(62.8594434, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Ce"
        ].value == pytest.approx(0.0001605177, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Dy"
        ].value == pytest.approx(1.692685e-11, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Fe"
        ].value == pytest.approx(236.51030445, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Gd"
        ].value == pytest.approx(1.99085e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "La"
        ].value == pytest.approx(0.020350852, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Nd"
        ].value == pytest.approx(2.6291397e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Pr"
        ].value == pytest.approx(7.024444e-6, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sc"
        ].value == pytest.approx(1.4297029e-7, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Sm"
        ].value == pytest.approx(3.484441e-11, 1e-4)
        assert model.fs.solex_cleaner.mscontactor.aqueous_outlet.conc_mass_comp[
            0, "Y"
        ].value == pytest.approx(9.89016e-11, 1e-4)

        assert model.fs.precipitator.cv_aqueous.properties_out[
            0
        ].flow_vol.value == pytest.approx(167.5946, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Al"
        ].value == pytest.approx(251.90944330, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ca"
        ].value == pytest.approx(70.3868727, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Ce"
        ].value == pytest.approx(2.01891339, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Dy"
        ].value == pytest.approx(0.0169268, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Fe"
        ].value == pytest.approx(252.0588739, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Gd"
        ].value == pytest.approx(0.0725528, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "La"
        ].value == pytest.approx(1.3342203, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Nd"
        ].value == pytest.approx(0.5135038, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Pr"
        ].value == pytest.approx(0.16383543, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sc"
        ].value == pytest.approx(142.97029033, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Sm"
        ].value == pytest.approx(0.0348444, 1e-4)
        assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
            "Y"
        ].value == pytest.approx(0.0989016, 1e-4)

        assert model.fs.precipitator.precipitate_outlet.temperature[
            0
        ].value == pytest.approx(348.15, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Al2(C2O4)3(s)"
        ].value == pytest.approx(0.007105083, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Ce2(C2O4)3(s)"
        ].value == pytest.approx(0.0025740, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Dy2(C2O4)3(s)"
        ].value == pytest.approx(5.9252256e-5, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Fe2(C2O4)3(s)"
        ].value == pytest.approx(0.0094594, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Gd2(C2O4)3(s)"
        ].value == pytest.approx(0.000283796, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "La2(C2O4)3(s)"
        ].value == pytest.approx(0.000855025, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Nd2(C2O4)3(s)"
        ].value == pytest.approx(0.001318589, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Pr2(C2O4)3(s)"
        ].value == pytest.approx(0.00034544, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sc2(C2O4)3(s)"
        ].value == pytest.approx(0.123201687, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Sm2(C2O4)3(s)"
        ].value == pytest.approx(0.00013409, 1e-4)
        assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
            0, "Y2(C2O4)3(s)"
        ].value == pytest.approx(0.00027177, 1e-4)

        assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
            0.00903269, 1e-4
        )
        assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
            873.15, 1e-4
        )
        assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(
            101325, 1e-4
        )
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "CO2"
        ].value == pytest.approx(0.0641327, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "H2O"
        ].value == pytest.approx(0.2021090, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "N2"
        ].value == pytest.approx(0.643809, 1e-4)
        assert model.fs.roaster.gas_outlet.mole_frac_comp[
            0, "O2"
        ].value == pytest.approx(0.089950, 1e-4)

    @pytest.mark.unit
    def test_costing(self, model):
        m = add_costing(model)

        assert m.fs.costing.total_plant_cost.value == pytest.approx(16.2919, rel=1e-4)
        assert m.fs.costing.total_BEC.value == pytest.approx(5.485, rel=1e-4)
        assert m.fs.costing.total_installation_cost.value == pytest.approx(
            10.805, rel=1e-4
        )
        assert m.fs.costing.other_plant_costs.value == pytest.approx(
            0.0016309, rel=1e-4
        )
        assert m.fs.costing.total_fixed_OM_cost.value == pytest.approx(7.2716, rel=1e-4)
        assert m.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
            10.7472, rel=1e-4
        )
        assert value(m.fs.costing.land_cost) == pytest.approx(
            6.1234e-5, rel=1e-4
        )  # Expression, not Var
        assert m.fs.costing.total_sales_revenue.value == pytest.approx(0.35555, rel=1e-4)

    @pytest.mark.component
    def test_report(self, model):
        m = add_costing(model)
        display_costing(m)
