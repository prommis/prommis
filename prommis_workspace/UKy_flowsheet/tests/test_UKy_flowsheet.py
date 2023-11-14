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

from workspace.prommis_workspace.UKy_flowsheet.UKy_flowsheet_without_recycle import (
    build,
    set_scaling,
    set_operating_conditions,
    initialize_system,
    solve,
)

from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.mscontactor import (
    MSContactor,
)

from idaes.models.unit_models.feed import (
    Feed,
)
from idaes.models.unit_models.solid_liquid import SLSeparator
from idaes.models.unit_models.product import (
    Product,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from workspace.leaching.leach_solution_properties import (
    LeachSolutionParameters,
)
from workspace.leaching.leach_solids_properties import (
    CoalRefuseParameters,
)
from workspace.leaching.leach_reactions import (
    CoalRefuseLeachingReactions,
)

from workspace.prommis_workspace.Solvent_Extraction.REESXmodel import REESX
from workspace.prommis_workspace.Solvent_Extraction.REEAqdistribution import REESolExAqParameters
from workspace.prommis_workspace.Solvent_Extraction.REEOgdistribution import REESolExOgParameters

from workspace.prommis_workspace.precipitate.precipitator import Precipitator
from workspace.prommis_workspace.precipitate.precip_prop import (
    AqueousStateParameterBlock,
    PrecipitateStateParameterBlock,
)

from workspace.prommis_workspace.roasting.ree_oxalate_roaster import REEOxalateRoaster

from workspace.prommis_workspace.UKy_flowsheet.Translators.translator_leaching_SX import (
    Translator_leaching_SX,
)
from workspace.prommis_workspace.UKy_flowsheet.Translators.translator_SX_precipitator import (
    Translator_SX_precipitator,
)

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)


@pytest.fixture(scope="module")
def model():
    m = build()

    return m


@pytest.mark.unit
def test_build_flowsheet(model):
    assert isinstance(model.fs, FlowsheetBlock)

    # Leaching property packages and unit models
    assert isinstance(model.fs.leach_soln, LeachSolutionParameters)
    assert isinstance(model.fs.coal, CoalRefuseParameters)
    assert isinstance(model.fs.leach_rxns, CoalRefuseLeachingReactions)

    assert isinstance(model.fs.leach, MSContactor)
    assert isinstance(model.fs.separator1, SLSeparator)
    assert isinstance(model.fs.leach_liquid_feed, Feed)
    assert isinstance(model.fs.leach_solid_feed, Feed)
    assert isinstance(model.fs.leach_filter_cake, Product)
    assert isinstance(model.fs.leach_filter_cake_liquid, Product)

    # Solvent extraction property packages and unit models
    assert isinstance(model.fs.prop_a, REESolExAqParameters)
    assert isinstance(model.fs.prop_o, REESolExOgParameters)

    assert isinstance(model.fs.leach_to_SX, Translator_leaching_SX)
    assert isinstance(model.fs.solex, REESX)
    assert isinstance(model.fs.sx_leach_acid, Product)

    # Precipitation property packages and unit models
    assert isinstance(model.fs.properties_aq, AqueousStateParameterBlock)
    assert isinstance(model.fs.properties_solid, PrecipitateStateParameterBlock)

    assert isinstance(model.fs.SX_to_precipitator, Translator_SX_precipitator)
    assert isinstance(model.fs.precipitator, Precipitator)
    assert isinstance(model.fs.separator2, SLSeparator)
    assert isinstance(model.fs.precipitate_feed, Feed)

    # Roasting property packages and unit models
    assert isinstance(model.fs.prop_gas, GenericParameterBlock)
    assert isinstance(model.fs.prop_solid, PrecipitateStateParameterBlock)

    assert isinstance(model.fs.roaster, REEOxalateRoaster)


    # Flowsheet connections
    assert isinstance(model.fs.s01, Arc)
    assert isinstance(model.fs.s02, Arc)
    assert isinstance(model.fs.sep1_solid, Arc)
    assert isinstance(model.fs.sep1_retained_liquid, Arc)
    assert isinstance(model.fs.sep1_liquid, Arc)
    assert isinstance(model.fs.s03, Arc)
    assert isinstance(model.fs.s04, Arc)
    assert isinstance(model.fs.s05, Arc)
    assert isinstance(model.fs.s06, Arc)
    assert isinstance(model.fs.s07, Arc)
    assert isinstance(model.fs.s08, Arc)
    assert isinstance(model.fs.s09, Arc)
    assert isinstance(model.fs.sep2_solid, Arc)

    assert degrees_of_freedom(model) == 63


@pytest.mark.unit
def test_set_dof(model):
    set_scaling(model)
    set_operating_conditions(model)

    assert degrees_of_freedom(model) == 0


@pytest.mark.unit
def test_initialize_flowsheet(model):
    initialize_system(model)

    assert degrees_of_freedom(model) == 0

@pytest.mark.integration
def test_unit_consistency(model):
    assert_units_consistent(model)


@pytest.mark.unit
def test_solve_flowsheet(model):
    solve(model)

    assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(22.24119, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Al2O3"].value == pytest.approx(0.234416, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "CaO"].value == pytest.approx(0.00163869, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Ce2O3"].value == pytest.approx(9.224196e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Dy2O3"].value == pytest.approx(5.975687e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Fe2O3"].value == pytest.approx(0.054770, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Gd2O3"].value == pytest.approx(4.131267e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "La2O3"].value == pytest.approx(4.261405e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Nd2O3"].value == pytest.approx(3.739941e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Pr2O3"].value == pytest.approx(9.016348e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Sc2O3"].value == pytest.approx(2.678979e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Sm2O3"].value == pytest.approx(1.155810e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(2.845565e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "inerts"].value == pytest.approx(0.708916, 1e-4)

    assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(224.47892, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(380.672853, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(122.97209, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(5.667589, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(0.148617, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(741.241302, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(0.556940, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(2.274231, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "HSO4"].value == pytest.approx(881.193491, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(3924.065118, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(2.238809, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(2.684072, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(0.716708, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(0.100454, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(0.301299, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(0.394626, 1e-4)

    assert model.fs.separator1.solid_outlet.flow_mass[0].value == pytest.approx(22.24119, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Al2O3"].value == pytest.approx(0.234416, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "CaO"].value == pytest.approx(0.00163869, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Ce2O3"].value == pytest.approx(9.224196e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Dy2O3"].value == pytest.approx(5.975687e-6, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Fe2O3"].value == pytest.approx(0.054770, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Gd2O3"].value == pytest.approx(4.131267e-6, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "La2O3"].value == pytest.approx(4.261405e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Nd2O3"].value == pytest.approx(3.739941e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Pr2O3"].value == pytest.approx(9.016348e-6, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Sc2O3"].value == pytest.approx(2.678979e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Sm2O3"].value == pytest.approx(1.155810e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(2.845565e-5, 1e-4)
    assert model.fs.separator1.solid_outlet.mass_frac_comp[0, "inerts"].value == pytest.approx(0.708916, 1e-4)

    assert model.fs.separator1.retained_liquid_outlet.flow_vol[0].value == pytest.approx(67.343677, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(380.672853, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(122.972089, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(5.667589, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(0.148617, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(741.241302, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(0.556940, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(2.274231, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "HSO4"].value == pytest.approx(881.193491, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(2.238809, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(2.684072, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(0.716708, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(3924.065118, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(0.100454, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(0.301299, 1e-4)
    assert model.fs.separator1.retained_liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(0.394626, 1e-4)

    assert model.fs.separator1.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(157.135245, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(380.672853, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(122.972089, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(5.667589, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(0.148617, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(741.241302, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(0.556940, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(2.274231, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "HSO4"].value == pytest.approx(881.193491, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(2.238809, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(2.684072, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(0.716708, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(3924.065118, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(0.100454, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(0.301299, 1e-4)
    assert model.fs.separator1.recovered_liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(0.394626, 1e-4)

    assert model.fs.leach_to_SX.properties_out[0].flow_vol.value == pytest.approx(157.135245, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Al"].value == pytest.approx(380.672853, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Ca"].value == pytest.approx(122.972089, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Ce"].value == pytest.approx(5.667589, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Dy"].value == pytest.approx(0.148617, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Fe"].value == pytest.approx(741.241302, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Gd"].value == pytest.approx(0.556940, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["La"].value == pytest.approx(2.238809, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Nd"].value == pytest.approx(2.684072, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Pr"].value == pytest.approx(0.716708, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Sc"].value == pytest.approx(0.100454, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Sm"].value == pytest.approx(0.301299, 1e-4)
    assert model.fs.leach_to_SX.properties_out[0].conc_mass_comp["Y"].value == pytest.approx(0.394626, 1e-4)

    assert model.fs.solex.Orgacid_outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(100.475263, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(33.325252, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(14.360703, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(0.376601, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(115.867046, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(1.411298, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(5.586677, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(6.801518, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(1.816083, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(321.594554, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(0.763501, 1e-4)
    assert model.fs.solex.Orgacid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(0.999996, 1e-4)

    assert model.fs.solex.Acidsoln_outlet.flow_vol[0].value == pytest.approx(157.135245, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(341.022481, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(109.821005, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(0.000450613, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(1.486174e-10, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(695.516898, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(1.528242e-6, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "La"].value == pytest.approx(0.034149, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(1.374245e-6, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(3.072887e-5, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(1.004542e-10, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(3.012988e-10, 1e-4)
    assert model.fs.solex.Acidsoln_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(3.946264e-10, 1e-4)

    assert model.fs.SX_to_precipitator.outlet.flow_mass[0].value == pytest.approx(0.042592, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(C2O4)2^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(C2O4)3^3-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(C2O4)^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(HC2O4)^2+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)(C2O4)"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)(HC2O4)2^2-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)2(C2O4)^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)2^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)3"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)4^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al(OH)^2+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al2(OH)2^4+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al3(OH)4^5+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Al^3+"].value == pytest.approx(-1.887652, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "C2O4^2-"].value == pytest.approx(-4, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ca(C2O4)"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ca(OH)^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ca^2+"].value == pytest.approx(-2.551581, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(C2O4)2^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(C2O4)3^3-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(C2O4)^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(OH)2^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(OH)3"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(OH)4^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce(OH)^2+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Ce^3+"].value == pytest.approx(-8.482045, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(C2O4)2^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(C2O4)3^3-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(C2O4)^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(OH)2^+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(OH)3"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(OH)4^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe(OH)^2+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe2(OH)2^4+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe3(OH)4^5+"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "Fe^3+"].value == pytest.approx(-1.894037, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "H2C2O4"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "HC2O4^-"].value == pytest.approx(-20, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "H^+"].value == pytest.approx(-2.5, 1e-4)
    assert model.fs.SX_to_precipitator.outlet.log10_molality_comp[0, "OH^-"].value == pytest.approx(-11.5, 1e-4)

    assert model.fs.precipitator.aqueous_outlet.flow_mass[0].value == pytest.approx(0.042592, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(C2O4)2^-"].value == pytest.approx(-7.296814, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(C2O4)3^3-"].value == pytest.approx(-11.412392, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(C2O4)^+"].value == pytest.approx(-4.106161, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(HC2O4)^2+"].value == pytest.approx(-5.455648, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)(C2O4)"].value == pytest.approx(-7.914906, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)(HC2O4)2^2-"].value == pytest.approx(-14.943300, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)2(C2O4)^-"].value == pytest.approx(-11.984883, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)2^+"].value == pytest.approx(-10.284230, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)3"].value == pytest.approx(-15.332975, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)4^-"].value == pytest.approx(-20.022952, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al(OH)^2+"].value == pytest.approx(-6.076717, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al2(OH)2^4+"].value == pytest.approx(-8.730122, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al3(OH)4^5+"].value == pytest.approx(-13.811041, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Al^3+"].value == pytest.approx(-1.890434, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "C2O4^2-"].value == pytest.approx(-8.371352, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ca(C2O4)"].value == pytest.approx(-8.769771, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ca(OH)^+"].value == pytest.approx(-14.169095, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ca^2+"].value == pytest.approx(-2.551582, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(C2O4)2^-"].value == pytest.approx(-16.818587, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(C2O4)3^3-"].value == pytest.approx(-23.784165, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(C2O4)^+"].value == pytest.approx(-11.907935, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(OH)2^+"].value == pytest.approx(-22.886004, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(OH)3"].value == pytest.approx(-31.234749, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(OH)4^-"].value == pytest.approx(-41.614725, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce(OH)^2+"].value == pytest.approx(-16.008490, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Ce^3+"].value == pytest.approx(-8.482208, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(C2O4)2^-"].value == pytest.approx(-7.353076, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(C2O4)3^3-"].value == pytest.approx(-10.768653, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(C2O4)^+"].value == pytest.approx(-4.782423, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(OH)2^+"].value == pytest.approx(-7.840492, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(OH)3"].value == pytest.approx(-15.739237, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(OH)4^-"].value == pytest.approx(-21.819213, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe(OH)^2+"].value == pytest.approx(-5.192978, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe2(OH)2^4+"].value == pytest.approx(-8.122645, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe3(OH)4^5+"].value == pytest.approx(-12.499826, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "Fe^3+"].value == pytest.approx(-3.986696, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "H2C2O4"].value == pytest.approx(-6.354168, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "HC2O4^-"].value == pytest.approx(-5.987144, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "H^+"].value == pytest.approx(-1.386492, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.log10_molality_comp[0, "OH^-"].value == pytest.approx(-12.418213, 1e-4)

    assert model.fs.precipitator.precipitate_outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Al(OH)3(s)"].value == pytest.approx(1.635349e-19, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Ca(C2O4)*3H2O(s)"].value == pytest.approx(2.930532e-19, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Ca(C2O4)*H2O(s)"].value == pytest.approx(3.339473e-19, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Ca(OH)2(s)"].value == pytest.approx(3.261348e-20, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Ce(OH)3(s)"].value == pytest.approx(2.272350e-20, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Ce2(C2O4)3(s)"].value == pytest.approx(8.749139e-20, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mass[0, "Fe2O3(s)"].value == pytest.approx(0.000268984, 1e-4)

    assert model.fs.separator2.solid_outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Al(OH)3(s)"].value == pytest.approx(1.635349e-19, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Ca(C2O4)*3H2O(s)"].value == pytest.approx(2.930532e-19, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Ca(C2O4)*H2O(s)"].value == pytest.approx(3.339473e-19, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Ca(OH)2(s)"].value == pytest.approx(3.261348e-20, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Ce(OH)3(s)"].value == pytest.approx(2.272350e-20, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Ce2(C2O4)3(s)"].value == pytest.approx(8.749139e-20, 1e-4)
    assert model.fs.separator2.solid_outlet.flow_mass[0, "Fe2O3(s)"].value == pytest.approx(0.000268984, 1e-4)

    assert model.fs.separator2.retained_liquid_outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.flow_mass[0].value == pytest.approx(1, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(C2O4)2^-"].value == pytest.approx(-2.189044, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(C2O4)3^3-"].value == pytest.approx(-3.423718, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(C2O4)^+"].value == pytest.approx(-1.231848, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(HC2O4)^2+"].value == pytest.approx(-1.636694, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)(C2O4)"].value == pytest.approx(-2.374472, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)(HC2O4)2^2-"].value == pytest.approx(-4.482990, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)2(C2O4)^-"].value == pytest.approx(-3.595465, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)2^+"].value == pytest.approx(-3.085269, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)3"].value == pytest.approx(-4.599893, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)4^-"].value == pytest.approx(-6.006885, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al(OH)^2+"].value == pytest.approx(-1.823015, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al2(OH)2^4+"].value == pytest.approx(-2.619037, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al3(OH)4^5+"].value == pytest.approx(-4.143312, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Al^3+"].value == pytest.approx(-0.567130, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "C2O4^2-"].value == pytest.approx(-2.511406, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ca(C2O4)"].value == pytest.approx(-2.63093, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ca(OH)^+"].value == pytest.approx(-4.250729, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ca^2+"].value == pytest.approx(-0.765475, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)2^-"].value == pytest.approx(-5.045576, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)3^3-"].value == pytest.approx(-7.135250, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)^+"].value == pytest.approx(-3.572380, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(OH)2^+"].value == pytest.approx(-6.865801, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(OH)3"].value == pytest.approx(-9.370425, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(OH)4^-"].value == pytest.approx(-12.484417, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce(OH)^2+"].value == pytest.approx(-4.802547, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Ce^3+"].value == pytest.approx(-2.544662, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)2^-"].value == pytest.approx(-2.205923, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)3^3-"].value == pytest.approx(-3.230596, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)^+"].value == pytest.approx(-1.434727, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(OH)2^+"].value == pytest.approx(-2.352148, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(OH)3"].value == pytest.approx(-4.721771, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(OH)4^-"].value == pytest.approx(-6.54576467, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe(OH)^2+"].value == pytest.approx(-1.557894, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe2(OH)2^4+"].value == pytest.approx(-2.436794, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe3(OH)4^5+"].value == pytest.approx(-3.749948, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "Fe^3+"].value == pytest.approx(-1.196009, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "H2C2O4"].value == pytest.approx(-1.906250, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "HC2O4^-"].value == pytest.approx(-1.796143, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "H^+"].value == pytest.approx(-0.415948, 1e-4)
    assert model.fs.separator2.retained_liquid_outlet.log10_molality_comp[0, "OH^-"].value == pytest.approx(-3.725464, 1e-4)

    assert model.fs.separator2.recovered_liquid_outlet.temperature[0].value == pytest.approx(348.15, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.flow_mass[0].value == pytest.approx(1, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(C2O4)2^-"].value == pytest.approx(-5.107770, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(C2O4)3^3-"].value == pytest.approx(-7.988674, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(C2O4)^+"].value == pytest.approx(-2.874313, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(HC2O4)^2+"].value == pytest.approx(-3.818953, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)(C2O4)"].value == pytest.approx(-5.540434, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)(HC2O4)2^2-"].value == pytest.approx(-10.460310, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)2(C2O4)^-"].value == pytest.approx(-8.389418, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)2^+"].value == pytest.approx(-7.198961, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)3"].value == pytest.approx(-10.733083, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)4^-"].value == pytest.approx(-14.016066, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al(OH)^2+"].value == pytest.approx(-4.253702, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al2(OH)2^4+"].value == pytest.approx(-6.111085, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al3(OH)4^5+"].value == pytest.approx(-9.667729, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Al^3+"].value == pytest.approx(-1.323304, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "C2O4^2-"].value == pytest.approx(-5.859946, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ca(C2O4)"].value == pytest.approx(-6.138840, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ca(OH)^+"].value == pytest.approx(-9.918367, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ca^2+"].value == pytest.approx(-1.786107, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)2^-"].value == pytest.approx(-11.773011, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)3^3-"].value == pytest.approx(-16.648916, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(C2O4)^+"].value == pytest.approx(-8.335554, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(OH)2^+"].value == pytest.approx(-16.020202, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(OH)3"].value == pytest.approx(-21.864324, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(OH)4^-"].value == pytest.approx(-29.130307, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce(OH)^2+"].value == pytest.approx(-11.205943, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Ce^3+"].value == pytest.approx(-5.937545, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)2^-"].value == pytest.approx(-5.147153, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)3^3-"].value == pytest.approx(-7.538057, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(C2O4)^+"].value == pytest.approx(-3.347696, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(OH)2^+"].value == pytest.approx(-5.488344, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(OH)3"].value == pytest.approx(-11.017466, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(OH)4^-"].value == pytest.approx(-15.273449, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe(OH)^2+"].value == pytest.approx(-3.635085, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe2(OH)2^4+"].value == pytest.approx(-5.685851, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe3(OH)4^5+"].value == pytest.approx(-8.749878, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "Fe^3+"].value == pytest.approx(-2.790687, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "H2C2O4"].value == pytest.approx(-4.447918, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "HC2O4^-"].value == pytest.approx(-4.191001, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "H^+"].value == pytest.approx(-0.970544, 1e-4)
    assert model.fs.separator2.recovered_liquid_outlet.log10_molality_comp[0, "OH^-"].value == pytest.approx(-8.692749, 1e-4)

    assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(0.008485, 1e-4)
    assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(873.15, 1e-4)
    assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "CO2"].value == pytest.approx(0.0396713, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "H2O"].value == pytest.approx(0.172057, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "N2"].value == pytest.approx(0.685365, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "O2"].value == pytest.approx(0.102906, 1e-4)
