#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for UKy flowsheet.

"""

from pyomo.environ import (
    TransformationFactory,
    assert_optimal_termination,
    units,
    value,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models.feed import Feed
from idaes.models.unit_models.mixer import Mixer
from idaes.models.unit_models.product import Product
from idaes.models.unit_models.separator import Separator
from idaes.models.unit_models.solid_liquid import SLSeparator

import pytest

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactions
from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.leaching.leach_train import LeachingTrain
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitate_reactions import OxalatePrecipitationReactions
from prommis.precipitate.precipitator import OxalatePrecipitator
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.uky.costing.ree_plant_capcost import QGESSCostingData
from prommis.uky.uky_flowsheet import (
    add_costing,
    build,
    display_costing,
    display_results,
    fix_organic_recycle,
    initialize_system,
    set_operating_conditions,
    set_scaling,
    solve_system,
)


@pytest.fixture(scope="module")
def system_frame():
    m = build()
    set_operating_conditions(m)

    return m


@pytest.mark.component
def test_structural_issues(system_frame):
    model = system_frame
    dt = DiagnosticsToolbox(model)
    dt.report_structural_issues()
    dt.display_overconstrained_set()

    dt.assert_no_structural_warnings()


@pytest.mark.unit
def test_build_flowsheet(system_frame):
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
    assert isinstance(model.fs.precip_rxns, OxalatePrecipitationReactions)

    assert isinstance(model.fs.precipitator, OxalatePrecipitator)
    assert isinstance(model.fs.sl_sep2, SLSeparator)
    assert isinstance(model.fs.precip_sx_mixer, Mixer)
    assert isinstance(model.fs.precip_purge, Product)

    # Roasting property packages and unit models
    assert isinstance(model.fs.prop_gas, GenericParameterBlock)
    assert isinstance(model.fs.prop_solid, PrecipitateParameters)

    assert isinstance(model.fs.roaster, REEOxalateRoaster)

    # Flowsheet connections
    assert isinstance(model.fs.leaching_sol_feed, Arc)
    assert isinstance(model.fs.leaching_liq_feed, Arc)
    assert isinstance(model.fs.leaching_feed_mixture, Arc)
    assert isinstance(model.fs.leaching_solid_outlet, Arc)
    assert isinstance(model.fs.leaching_liquid_outlet, Arc)
    assert isinstance(model.fs.sl_sep1_solid_outlet, Arc)
    assert isinstance(model.fs.sl_sep1_retained_liquid_outlet, Arc)
    assert isinstance(model.fs.sl_sep1_liquid_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_load_aq_feed, Arc)
    assert isinstance(model.fs.sx_rougher_org_feed, Arc)
    assert isinstance(model.fs.sx_rougher_mixed_org_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_load_aq_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_load_aq_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_load_org_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_scrub_acid_feed, Arc)
    assert isinstance(model.fs.sx_rougher_scrub_aq_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_scrub_aq_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_scrub_org_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_strip_acid_feed, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_purge, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_strip_aq_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_load_aq_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_org_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_mixed_org_recycle, Arc)
    assert isinstance(model.fs.sx_cleaner_load_aq_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_acid_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_load_org_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_purge, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_recycle, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_aq_outlet, Arc)
    assert isinstance(model.fs.precip_solid_outlet, Arc)
    assert isinstance(model.fs.precip_aq_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_solid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_liquid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_aq_purge, Arc)
    assert isinstance(model.fs.sl_sep2_aq_recycle, Arc)


@pytest.mark.component
@pytest.mark.solver
def test_solve(system_frame):
    model = system_frame

    set_scaling(model)

    scaling = TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(model, rename=False)

    initialize_system(scaled_model)

    solve_system(scaled_model)

    fix_organic_recycle(scaled_model)

    results = solve_system(scaled_model)

    scaling.propagate_solution(scaled_model, model)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_solution(system_frame):
    model = system_frame

    assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
        22.2454484, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Al2O3"
    ].value == pytest.approx(0.2326854366, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "CaO"].value == pytest.approx(
        0.00214479, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Ce2O3"
    ].value == pytest.approx(1.04406e-4, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Dy2O3"
    ].value == pytest.approx(6.508519e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Fe2O3"
    ].value == pytest.approx(0.056103923, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Gd2O3"
    ].value == pytest.approx(3.7753056e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "La2O3"
    ].value == pytest.approx(4.931877e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Nd2O3"
    ].value == pytest.approx(4.1483573e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Pr2O3"
    ].value == pytest.approx(1.065413119e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sc2O3"
    ].value == pytest.approx(2.7121932e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sm2O3"
    ].value == pytest.approx(1.243956e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(
        2.9827698e-5, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "inerts"
    ].value == pytest.approx(0.70878, 1e-4)

    assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
        926.337724, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(
        307.2412917, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(
        57.04876439, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(
        3.0360766, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Cl"].value == pytest.approx(
        146.5275185, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(
        0.0350188, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(
        424.4005869, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(
        0.36763938, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(
        2.7608537, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[
        0, "H2C2O4"
    ].value == pytest.approx(527.4485164, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(
        1000000.0, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(674.38469, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(
        2473.79559, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(
        1.09416279, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(
        1.52977089, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(
        0.3783432, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(
        0.0210985, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(
        0.147079, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(
        0.0887758, 1e-4
    )

    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.01, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al_o"
    ].value == pytest.approx(0.0483158, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca_o"
    ].value == pytest.approx(0.0189648, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce_o"
    ].value == pytest.approx(0.0053627, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy_o"
    ].value == pytest.approx(1.005819, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe_o"
    ].value == pytest.approx(1.91696046, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd_o"
    ].value == pytest.approx(0.130396, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La_o"
    ].value == pytest.approx(0.003697, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd_o"
    ].value == pytest.approx(0.00319497, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr_o"
    ].value == pytest.approx(0.001004536, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc_o"
    ].value == pytest.approx(1.8641447, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm_o"
    ].value == pytest.approx(0.00442077, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y_o"
    ].value == pytest.approx(4.086477, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(39.537047, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2C2O4"
    ].value == pytest.approx(4.293345e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(5.698696e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(5.879417e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(3.3969304, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.9239557, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.2401165, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(0.9811165, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(25.2456133, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(1.4885115, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.174653, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.129265, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(0.032062, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.00294795, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(0.11509334, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(1.829172, 1e-4)

    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.010, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al_o"
    ].value == pytest.approx(0.00132359, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca_o"
    ].value == pytest.approx(0.00076523, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce_o"
    ].value == pytest.approx(0.000319098, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy_o"
    ].value == pytest.approx(0.43903026, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe_o"
    ].value == pytest.approx(0.24553399, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd_o"
    ].value == pytest.approx(0.03692634, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La_o"
    ].value == pytest.approx(0.00025373, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd_o"
    ].value == pytest.approx(0.00021787289, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr_o"
    ].value == pytest.approx(9.67546e-5, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc_o"
    ].value == pytest.approx(0.0042601826, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm_o"
    ].value == pytest.approx(0.000595803, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y_o"
    ].value == pytest.approx(1.2759757, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(0.27885617, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.09313549, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.0327623, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(0.4677576, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(4.778415, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(0.63568078, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(41.10369, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2C2O4"
    ].value == pytest.approx(6.02203e-6, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(2.832189e-5, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(2.734739e-5, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.0285543, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.01859269, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(0.005211387, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(6.779843e-6, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(0.030301444, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(0.61573168, 1e-4)

    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(0.0626532, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.02095606, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(8.710016e-5, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(323.676, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(0.02853685, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(1.0623022, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(0.00731239, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(9.25731, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "H2C2O4"
    ].value == pytest.approx(6199.59893, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(7.41622e-6, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(7.259714e-6, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.001157344, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(5.917e-5, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(3.42338e-5, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(4.76e-6, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(0.000174336, 1e-4)
    assert model.fs.precipitator.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(0.02911317, 1e-4)

    assert model.fs.precipitator.precipitate_outlet.temperature[
        0
    ].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Al2(C2O4)3(s)"
    ].value == pytest.approx(6.741e-8, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ca(C2O4)(s)"
    ].value == pytest.approx(6.8603e-10, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ce2(C2O4)3(s)"
    ].value == pytest.approx(1.0413874e-6, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Dy2(C2O4)3(s)"
    ].value == pytest.approx(9.4421e-6, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Fe2(C2O4)3(s)"
    ].value == pytest.approx(4.6e-6, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Gd2(C2O4)3(s)"
    ].value == pytest.approx(1.72622e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "La2(C2O4)3(s)"
    ].value == pytest.approx(7.593e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Nd2(C2O4)3(s)"
    ].value == pytest.approx(5.7345e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Pr2(C2O4)3(s)"
    ].value == pytest.approx(1.63043e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sc2(C2O4)3(s)"
    ].value == pytest.approx(1.28184e-9, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sm2(C2O4)3(s)"
    ].value == pytest.approx(8.8516e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Y2(C2O4)3(s)"
    ].value == pytest.approx(2.46175e-5, 1e-4)

    assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
        6.95386e-2, 1e-4
    )
    assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
        873.15, 1e-4
    )
    assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "CO2"].value == pytest.approx(
        0.0048422, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "H2O"].value == pytest.approx(
        0.898975, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "N2"].value == pytest.approx(
        0.0836273, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "O2"].value == pytest.approx(
        0.0125561, 1e-4
    )


@pytest.mark.component
@pytest.mark.solver
def test_REE_conservation(system_frame):
    model = system_frame

    REE_mass_frac = {
        "Y2O3": 88.906 * 2 / (88.906 * 2 + 16 * 3),
        "La2O3": 138.91 * 2 / (138.91 * 2 + 16 * 3),
        "Ce2O3": 140.12 * 2 / (140.12 * 2 + 16 * 3),
        "Pr2O3": 140.91 * 2 / (140.91 * 2 + 16 * 3),
        "Nd2O3": 144.24 * 2 / (144.24 * 2 + 16 * 3),
        "Sm2O3": 150.36 * 2 / (150.36 * 2 + 16 * 3),
        "Gd2O3": 157.25 * 2 / (157.25 * 2 + 16 * 3),
        "Dy2O3": 162.5 * 2 / (162.5 * 2 + 16 * 3),
    }

    molar_mass = {
        "Y2O3": (88.906 * 2 + 16 * 3) * units.g / units.mol,
        "La2O3": (138.91 * 2 + 16 * 3) * units.g / units.mol,
        "Ce2O3": (140.12 * 2 + 16 * 3) * units.g / units.mol,
        "Pr2O3": (140.91 * 2 + 16 * 3) * units.g / units.mol,
        "Nd2O3": (144.24 * 2 + 16 * 3) * units.g / units.mol,
        "Sm2O3": (150.36 * 2 + 16 * 3) * units.g / units.mol,
        "Gd2O3": (157.25 * 2 + 16 * 3) * units.g / units.mol,
        "Dy2O3": (162.5 * 2 + 16 * 3) * units.g / units.mol,
    }

    component_list = [
        ("Y2O3", "Y", "Y_o"),
        ("La2O3", "La", "La_o"),
        ("Ce2O3", "Ce", "Ce_o"),
        ("Pr2O3", "Pr", "Pr_o"),
        ("Nd2O3", "Nd", "Nd_o"),
        ("Sm2O3", "Sm", "Sm_o"),
        ("Gd2O3", "Gd", "Gd_o"),
        ("Dy2O3", "Dy", "Dy_o"),
    ]

    for REO, REE, Org in component_list:
        solid_feed = value(
            units.convert(
                model.fs.leach_solid_feed.flow_mass[0]
                * model.fs.leach_solid_feed.mass_frac_comp[0, REO]
                * REE_mass_frac[REO],
                to_units=units.kg / units.hr,
            )
        )

        liquid_feed = value(
            units.convert(
                model.fs.leach_liquid_feed.conc_mass_comp[0, REE]
                * model.fs.leach_liquid_feed.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        acid_feeds = value(
            units.convert(
                model.fs.acid_feed1.conc_mass_comp[0, REE]
                * model.fs.acid_feed1.flow_vol[0]
                + model.fs.acid_feed2.conc_mass_comp[0, REE]
                * model.fs.acid_feed2.flow_vol[0]
                + model.fs.acid_feed3.conc_mass_comp[0, REE]
                * model.fs.acid_feed3.flow_vol[0]
                + model.fs.oxalic_acid_feed.conc_mass_comp[0, REE]
                * model.fs.oxalic_acid_feed.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        filter_cake = value(
            units.convert(
                model.fs.sl_sep1.solid_outlet.flow_mass[0]
                * model.fs.sl_sep1.solid_outlet.mass_frac_comp[0, REO]
                * REE_mass_frac[REO],
                to_units=units.kg / units.hr,
            )
        )

        filter_cake_liquid = value(
            units.convert(
                model.fs.sl_sep1.retained_liquid_outlet.conc_mass_comp[0, REE]
                * model.fs.sl_sep1.retained_liquid_outlet.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        load_purge = value(
            units.convert(
                model.fs.load_sep.purge.conc_mass_comp[0, REE]
                * model.fs.load_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        scrub_purge = value(
            units.convert(
                model.fs.scrub_sep.purge.conc_mass_comp[0, REE]
                * model.fs.scrub_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        precip_purge = value(
            units.convert(
                model.fs.precip_sep.purge.conc_mass_comp[0, REE]
                * model.fs.precip_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        roaster_retained_liquid = value(
            units.convert(
                model.fs.sl_sep2.retained_liquid_outlet.conc_mass_comp[0, REE]
                * model.fs.sl_sep2.retained_liquid_outlet.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        roaster_product = value(
            units.convert(
                model.fs.roaster.flow_mol_comp_product[0, REE]
                * molar_mass[REO]
                * REE_mass_frac[REO],
                to_units=units.kg / units.hr,
            )
        )

        roaster_dust = value(
            units.convert(
                model.fs.roaster.flow_mol_comp_dust[0, REE]
                * molar_mass[REO]
                * REE_mass_frac[REO],
                to_units=units.kg / units.hr,
            )
        )

        rougher_purge = value(
            units.convert(
                model.fs.rougher_sep.purge.conc_mass_comp[0, Org]
                * model.fs.rougher_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        cleaner_purge = value(
            units.convert(
                model.fs.cleaner_sep.purge.conc_mass_comp[0, Org]
                * model.fs.cleaner_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        assert solid_feed + liquid_feed + acid_feeds == pytest.approx(
            filter_cake
            + filter_cake_liquid
            + load_purge
            + rougher_purge
            + scrub_purge
            + cleaner_purge
            + precip_purge
            + roaster_retained_liquid
            + roaster_product
            + roaster_dust,
            rel=1e-6,
            abs=1e-6,
        )


@pytest.mark.component
@pytest.mark.solver
def test_costing(system_frame):
    model = system_frame
    add_costing(model)

    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_costing_diagnostics(system_frame):
    model = system_frame
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_costing_initialize(system_frame):
    model = system_frame
    QGESSCostingData.costing_initialization(model.fs.costing)
    QGESSCostingData.initialize_fixed_OM_costs(model.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(model.fs.costing)


@pytest.mark.component
@pytest.mark.solver
def test_costing_solve(system_frame):
    model = system_frame
    results = solve_system(model)
    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_costing_solution(system_frame):
    model = system_frame

    assert model.fs.costing.total_plant_cost.value == pytest.approx(1.4426, rel=1e-4)
    assert model.fs.costing.total_BEC.value == pytest.approx(0.48573, rel=1e-4)
    assert model.fs.costing.total_installation_cost.value == pytest.approx(
        0.95688, rel=1e-4
    )
    assert model.fs.costing.other_plant_costs.value == pytest.approx(
        1.04575e-5, rel=1e-4
    )
    assert model.fs.costing.total_fixed_OM_cost.value == pytest.approx(6.8244, rel=1e-4)
    assert model.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        1.36498, rel=1e-4
    )
    assert value(model.fs.costing.land_cost) == pytest.approx(6.1234e-5, rel=1e-4)
    assert model.fs.costing.total_sales_revenue.value == pytest.approx(
        2.9486e-5, rel=1e-4
    )


@pytest.mark.unit
def test_display(system_frame):
    model = system_frame
    display_results(model)
    display_costing(model)
