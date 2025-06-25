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
from prommis.precipitate.precipitator import Precipitator
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

    assert isinstance(model.fs.precipitator, Precipitator)
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
        22.238837, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Al2O3"
    ].value == pytest.approx(0.232736, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "CaO"].value == pytest.approx(
        0.0020873, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Ce2O3"
    ].value == pytest.approx(1.027533e-4, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Dy2O3"
    ].value == pytest.approx(6.459123e-6, abs=1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Fe2O3"
    ].value == pytest.approx(0.05590266, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Gd2O3"
    ].value == pytest.approx(3.71286e-6, abs=1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "La2O3"
    ].value == pytest.approx(4.8569e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Nd2O3"
    ].value == pytest.approx(4.07967e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Pr2O3"
    ].value == pytest.approx(1.043202e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sc2O3"
    ].value == pytest.approx(2.70925995e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sm2O3"
    ].value == pytest.approx(1.234895e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(
        2.969969e-5, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "inerts"
    ].value == pytest.approx(0.708991, 1e-4)

    assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
        865.2653, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(
        329.3326, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(
        63.95648, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(
        3.349581, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Cl"].value == pytest.approx(
        156.86977, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(
        0.03894, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(
        464.83, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(
        0.394755, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(
        2.58245, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(
        1000000.0, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(684.7322, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(
        2685.26996, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(
        1.216377, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(
        1.67917, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(
        0.418309, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(
        0.023101, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(
        0.162865, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(
        0.097377, 1e-4
    )

    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.01, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al_o"
    ].value == pytest.approx(0.053654, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca_o"
    ].value == pytest.approx(0.022027, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce_o"
    ].value == pytest.approx(0.0063643, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy_o"
    ].value == pytest.approx(1.19399, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe_o"
    ].value == pytest.approx(2.17512, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd_o"
    ].value == pytest.approx(0.15928, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La_o"
    ].value == pytest.approx(0.0045008, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd_o"
    ].value == pytest.approx(0.0037576, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr_o"
    ].value == pytest.approx(0.0011516, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc_o"
    ].value == pytest.approx(1.9319, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm_o"
    ].value == pytest.approx(0.0054847, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y_o"
    ].value == pytest.approx(4.6218, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9.0, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(39.28586, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(2.41459e-7, abs=1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(2.70217e-7, abs=1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(3.77224, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(1.07311, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.28459, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(1.15087, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(28.64549, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(1.80357, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.212375, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.15231, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(0.037087, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.00305517, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(0.14207, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(2.04473, 1e-4)

    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.010, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al_o"
    ].value == pytest.approx(0.0041203, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca_o"
    ].value == pytest.approx(0.0024165, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce_o"
    ].value == pytest.approx(0.0006731, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy_o"
    ].value == pytest.approx(0.373165, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe_o"
    ].value == pytest.approx(0.834055, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd_o"
    ].value == pytest.approx(0.0455115, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La_o"
    ].value == pytest.approx(0.00046228, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd_o"
    ].value == pytest.approx(0.0004254, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr_o"
    ].value == pytest.approx(0.00016404, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc_o"
    ].value == pytest.approx(0.00441393, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm_o"
    ].value == pytest.approx(0.00089498, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y_o"
    ].value == pytest.approx(1.15377, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(0.869364, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.29638, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.07427, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(0.38478, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(16.23179, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(0.76294, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(40.4119, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(1.22141e-6, abs=1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(1.13409e-6, abs=1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.05575, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.040607, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(0.011155, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(6.993204e-6, abs=1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(0.04601, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(0.53946, 1e-4)

    assert model.fs.precipitator.cv_aqueous.properties_out[
        0
    ].flow_vol.value == pytest.approx(9, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Al"
    ].value == pytest.approx(0.86154, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ca"
    ].value == pytest.approx(0.23563, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ce"
    ].value == pytest.approx(0.023715, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Dy"
    ].value == pytest.approx(0.0494056, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Fe"
    ].value == pytest.approx(15.835737, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Gd"
    ].value == pytest.approx(0.091476, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "H"
    ].value == pytest.approx(40.41, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "HSO4"
    ].value == pytest.approx(1.2214e-6, abs=1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "SO4"
    ].value == pytest.approx(1.134087e-6, abs=1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "La"
    ].value == pytest.approx(0.027033, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Nd"
    ].value == pytest.approx(0.007492, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Pr"
    ].value == pytest.approx(0.002454, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sc"
    ].value == pytest.approx(4.78265e-6, abs=1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sm"
    ].value == pytest.approx(0.0058206, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Y"
    ].value == pytest.approx(0.137778, 1e-4)

    assert model.fs.precipitator.precipitate_outlet.temperature[
        0
    ].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Al2(C2O4)3(s)"
    ].value == pytest.approx(1.3049e-6, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ca(C2O4)(s)"
    ].value == pytest.approx(1.36438e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ce2(C2O4)3(s)"
    ].value == pytest.approx(1.6235e-6, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Dy2(C2O4)3(s)"
    ].value == pytest.approx(9.2872e-6, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Fe2(C2O4)3(s)"
    ].value == pytest.approx(3.1914e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Gd2(C2O4)3(s)"
    ].value == pytest.approx(1.92151e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "La2(C2O4)3(s)"
    ].value == pytest.approx(9.30326e-7, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Nd2(C2O4)3(s)"
    ].value == pytest.approx(1.0331e-6, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Pr2(C2O4)3(s)"
    ].value == pytest.approx(2.77874e-7, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sc2(C2O4)3(s)"
    ].value == pytest.approx(2.28046e-10, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sm2(C2O4)3(s)"
    ].value == pytest.approx(1.20288e-6, abs=1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Y2(C2O4)3(s)"
    ].value == pytest.approx(2.03315e-5, 1e-4)

    assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
        0.021699, 1e-4
    )
    assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
        873.15, 1e-4
    )
    assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "CO2"].value == pytest.approx(
        0.0155197, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "H2O"].value == pytest.approx(
        0.676245, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "N2"].value == pytest.approx(
        0.267998, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "O2"].value == pytest.approx(
        0.040237, 1e-4
    )

    assert model.fs.leach_mixer.outlet.flow_vol[0].value == pytest.approx(865.078, 1e-4)
    assert model.fs.rougher_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_rougher_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(621.9757, 1e-4)
    assert model.fs.solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(90, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.acid_feed1.outlet.flow_vol[0].value == pytest.approx(90, 1e-4)
    assert model.fs.acid_feed2.outlet.flow_vol[0].value == pytest.approx(9, 1e-4)
    assert model.fs.acid_feed3.outlet.flow_vol[0].value == pytest.approx(9, 1e-4)
    assert model.fs.rougher_sep.inlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.load_sep.inlet.flow_vol[0].value == pytest.approx(621.9757, 1e-4)
    assert model.fs.scrub_sep.inlet.flow_vol[0].value == pytest.approx(90, 1e-4)
    assert model.fs.rougher_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.sc_circuit_purge.inlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(16.29, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.cleaner_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.cleaner_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.cleaner_sep.inlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.leach_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        621.97569, 1e-4
    )
    assert model.fs.cleaner_purge.inlet.flow_vol[0].value == pytest.approx(6.201, 1e-4)
    assert model.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        605.68569, 1e-4
    )
    assert model.fs.sl_sep2.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        8.1, 1e-4
    )
    assert model.fs.precip_sep.inlet.flow_vol[0].value == pytest.approx(8.1, 1e-4)
    assert model.fs.precip_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        16.29, 1e-4
    )
    assert model.fs.precip_purge.inlet.flow_vol[0].value == pytest.approx(0.81, 1e-4)


@pytest.mark.component
@pytest.mark.solver
def test_conservation(system_frame):
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

    for REO, REE, REE_o in component_list:
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
                model.fs.rougher_sep.purge.conc_mass_comp[0, REE_o]
                * model.fs.rougher_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        cleaner_purge = value(
            units.convert(
                model.fs.cleaner_sep.purge.conc_mass_comp[0, REE_o]
                * model.fs.cleaner_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        assert solid_feed + liquid_feed == pytest.approx(
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
            rel=1e-7,
            abs=1e-7,
        )


@pytest.mark.component
@pytest.mark.solver
def test_costing(system_frame):
    model = system_frame
    add_costing(model)


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

    assert model.fs.costing.total_plant_cost.value == pytest.approx(1.18071, rel=1e-4)
    assert model.fs.costing.total_BEC.value == pytest.approx(0.39754, rel=1e-4)
    assert model.fs.costing.total_installation_cost.value == pytest.approx(
        0.78316, rel=1e-4
    )
    assert model.fs.costing.other_plant_costs.value == pytest.approx(
        7.5987e-06, rel=1e-4
    )
    assert model.fs.costing.total_fixed_OM_cost.value == pytest.approx(6.8166, rel=1e-4)
    assert model.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        1.3634, rel=1e-4
    )
    assert value(model.fs.costing.land_cost) == pytest.approx(6.1234e-5, rel=1e-4)
    assert model.fs.costing.total_sales_revenue.value == pytest.approx(
        2.69422e-5, abs=1e-4
    )


@pytest.mark.component
@pytest.mark.solver
def test_costing_solution_diagnostics(system_frame):

    model = system_frame
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()


@pytest.mark.unit
def test_display(system_frame):
    model = system_frame
    display_results(model)
    display_costing(model)
