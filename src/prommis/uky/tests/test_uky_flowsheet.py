#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for UKy flowsheet.

"""

from pyomo.network import Arc
from pyomo.environ import (
    assert_optimal_termination,
    value,
    units,
    TransformationFactory,
)

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
    build,
    set_partition_coefficients,
    set_operating_conditions,
    set_scaling,
    initialize_system,
    solve_system,
    fix_organic_recycle,
    display_results,
    add_costing,
    display_costing,
)


@pytest.fixture(scope="module")
def system_frame():
    m = build()
    set_partition_coefficients(m)
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

    scaling = TransformationFactory("core.scale_model")
    scaling.propagate_solution(scaled_model, model)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_numerical_issues(system_frame):
    dt = DiagnosticsToolbox(system_frame)
    dt.report_numerical_issues()
    dt.display_variables_at_or_outside_bounds()
    dt.assert_no_numerical_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_solution(system_frame):
    model = system_frame

    assert model.fs.leach.solid_outlet.flow_mass[0].value == pytest.approx(
        22.234694, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Al2O3"
    ].value == pytest.approx(0.23347629, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "CaO"].value == pytest.approx(
        0.0017975, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Ce2O3"
    ].value == pytest.approx(9.598532e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Dy2O3"
    ].value == pytest.approx(6.2209025e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Fe2O3"
    ].value == pytest.approx(0.0553200, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Gd2O3"
    ].value == pytest.approx(3.6459069e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "La2O3"
    ].value == pytest.approx(4.51241048e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Nd2O3"
    ].value == pytest.approx(3.8236866e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Pr2O3"
    ].value == pytest.approx(9.4696496e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sc2O3"
    ].value == pytest.approx(2.6963267e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sm2O3"
    ].value == pytest.approx(1.1949254e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(
        2.91182383e-5, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "inerts"
    ].value == pytest.approx(0.7091231, 1e-4)

    assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
        624.4944177, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(
        421.52139874, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(
        108.574479, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(
        2.261334816, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(
        0.04680713, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(
        684.508471, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(
        0.2571, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(
        1.7700757, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(
        1000000.0, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(694.2926543, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(
        3972.37476, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(
        0.977964, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(
        0.94020, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(
        0.30123, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(
        0.032418, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(
        0.096468, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(
        0.123607, 1e-4
    )

    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.01, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(1.271694e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(2.67625214e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.0003043025, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(7.988351e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(2.8737799e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(3.3597390e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.00010451, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.00016553, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(3.7059156e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(1.7338035, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(1.6985006e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(2.1763374e-5, 1e-4)
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
    ].value == pytest.approx(3.940239e-9, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(1.6122889e-8, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(350.470145, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(737.556650, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(8386.369603, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(220.153435, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(79.19939, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(925.920992, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(2880.493852, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(4562.15409, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(1021.3248275, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(36.660503, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(468.095076, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(599.78359, 1e-4)

    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.010, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(7.21706637e-9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(1.56342660e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(1.5333975e-6, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(4.92297979e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(9.29630099e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(1.741709e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(5.865032e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(8.3523683e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(1.88444340e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.4329824, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(8.6417247e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(1.273777e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(0.3978046, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.861760, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(84.520856, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(2.713546, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(0.0512348, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(9.6003014, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(32.328056, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(46.03824, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(10.387050, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.1383766, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(4.763318, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(7.021059, 1e-4)

    assert model.fs.precipitator.cv_aqueous.properties_out[
        0
    ].flow_vol.value == pytest.approx(9, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Al"
    ].value == pytest.approx(0.3971172, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ca"
    ].value == pytest.approx(0.8617606, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ce"
    ].value == pytest.approx(0.8382345, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Dy"
    ].value == pytest.approx(0.6464795, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Fe"
    ].value == pytest.approx(0.0487210, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Gd"
    ].value == pytest.approx(0.430704, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "La"
    ].value == pytest.approx(5.08062, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Nd"
    ].value == pytest.approx(0.526176, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Pr"
    ].value == pytest.approx(0.22001, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sc"
    ].value == pytest.approx(0.088505, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sm"
    ].value == pytest.approx(0.103999, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Y"
    ].value == pytest.approx(1.291948, 1e-4)

    assert model.fs.precipitator.precipitate_outlet.temperature[
        0
    ].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Al2(C2O4)3(s)"
    ].value == pytest.approx(1.1464429e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ce2(C2O4)3(s)"
    ].value == pytest.approx(0.0026875, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Dy2(C2O4)3(s)"
    ].value == pytest.approx(5.724184e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Fe2(C2O4)3(s)"
    ].value == pytest.approx(2.0255644e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Gd2(C2O4)3(s)"
    ].value == pytest.approx(0.0002624, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "La2(C2O4)3(s)"
    ].value == pytest.approx(8.82714e-4, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Nd2(C2O4)3(s)"
    ].value == pytest.approx(1.41986e-3, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Pr2(C2O4)3(s)"
    ].value == pytest.approx(3.24694e-4, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sc2(C2O4)3(s)"
    ].value == pytest.approx(4.993087e-6, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sm2(C2O4)3(s)"
    ].value == pytest.approx(1.39444e-4, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Y2(C2O4)3(s)"
    ].value == pytest.approx(2.89983e-4, 1e-4)

    assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
        0.0245011, 1e-4
    )
    assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
        873.15, 1e-4
    )
    assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "CO2"].value == pytest.approx(
        0.01415145, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "H2O"].value == pytest.approx(
        0.712964, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "N2"].value == pytest.approx(
        0.237349, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "O2"].value == pytest.approx(
        0.035534, 1e-4
    )

    assert model.fs.leach_mixer.outlet.flow_vol[0].value == pytest.approx(
        624.30868, 1e-4
    )
    assert model.fs.rougher_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_rougher_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(444.364092, 1e-4)
    assert model.fs.solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(0.09, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(0.09, 1e-4)
    assert model.fs.acid_feed1.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
    assert model.fs.acid_feed2.outlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
    assert model.fs.acid_feed3.outlet.flow_vol[0].value == pytest.approx(9, 1e-4)
    assert model.fs.rougher_sep.inlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.load_sep.inlet.flow_vol[0].value == pytest.approx(444.364092, 1e-4)
    assert model.fs.scrub_sep.inlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
    assert model.fs.rougher_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.sc_circuit_purge.inlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(7.218, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.cleaner_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.cleaner_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.cleaner_sep.inlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.leach_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        444.364092, 1e-4
    )
    assert model.fs.cleaner_purge.inlet.flow_vol[0].value == pytest.approx(6.201, 1e-4)
    assert model.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        437.146092, 1e-4
    )
    assert model.fs.sl_sep2.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        7.9199, 1e-4
    )
    assert model.fs.precip_sep.inlet.flow_vol[0].value == pytest.approx(7.9199, 1e-4)
    assert model.fs.precip_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        7.218, 1e-4
    )
    assert model.fs.precip_purge.inlet.flow_vol[0].value == pytest.approx(0.79199, 1e-4)


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
        ("Y2O3", "Y"),
        ("La2O3", "La"),
        ("Ce2O3", "Ce"),
        ("Pr2O3", "Pr"),
        ("Nd2O3", "Nd"),
        ("Sm2O3", "Sm"),
        ("Gd2O3", "Gd"),
        ("Dy2O3", "Dy"),
    ]

    for REO, REE in component_list:
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
                model.fs.rougher_sep.purge.conc_mass_comp[0, REE]
                * model.fs.rougher_sep.purge.flow_vol[0],
                to_units=units.kg / units.hr,
            )
        )

        cleaner_purge = value(
            units.convert(
                model.fs.cleaner_sep.purge.conc_mass_comp[0, REE]
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

    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


@pytest.mark.component
@pytest.mark.solver
def test_costing_solution(system_frame):
    model = system_frame

    assert model.fs.costing.total_plant_cost.value == pytest.approx(16.0115, rel=1e-4)
    assert model.fs.costing.total_BEC.value == pytest.approx(5.3905, rel=1e-4)
    assert model.fs.costing.total_installation_cost.value == pytest.approx(
        10.619, rel=1e-4
    )
    assert model.fs.costing.other_plant_costs.value == pytest.approx(
        0.0016309, rel=1e-4
    )
    assert model.fs.costing.total_fixed_OM_cost.value == pytest.approx(7.2615, rel=1e-4)
    assert model.fs.costing.total_variable_OM_cost[0].value == pytest.approx(
        1.4573, rel=1e-4
    )
    assert value(model.fs.costing.land_cost) == pytest.approx(6.1234e-5, rel=1e-4)
    assert model.fs.costing.total_sales_revenue.value == pytest.approx(
        0.00019575, rel=1e-4
    )


@pytest.mark.unit
def test_display(system_frame):
    model = system_frame
    display_results(model)
    display_costing(model)
