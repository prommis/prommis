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
from prommis.precipitate.precipitate_reactions import OxalatePrecipitationReactions
from prommis.precipitate.precipitator import OxalatePrecipitator
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
        0.0017995, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Ce2O3"
    ].value == pytest.approx(9.601668e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Dy2O3"
    ].value == pytest.approx(6.22238e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Fe2O3"
    ].value == pytest.approx(0.0553200, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Gd2O3"
    ].value == pytest.approx(3.64403e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "La2O3"
    ].value == pytest.approx(4.51443e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Nd2O3"
    ].value == pytest.approx(3.82457e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Pr2O3"
    ].value == pytest.approx(9.474788e-6, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sc2O3"
    ].value == pytest.approx(2.69655e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "Sm2O3"
    ].value == pytest.approx(1.19513e-5, 1e-4)
    assert model.fs.leach.solid_outlet.mass_frac_comp[0, "Y2O3"].value == pytest.approx(
        2.9121365e-5, 1e-4
    )
    assert model.fs.leach.solid_outlet.mass_frac_comp[
        0, "inerts"
    ].value == pytest.approx(0.70915, 1e-4)

    assert model.fs.leach.liquid_outlet.flow_vol[0].value == pytest.approx(
        625.87396, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Al"].value == pytest.approx(
        421.19129, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ca"].value == pytest.approx(
        108.196656, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Ce"].value == pytest.approx(
        2.255278, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Cl"].value == pytest.approx(
        43.65097, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Dy"].value == pytest.approx(
        0.046659, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Fe"].value == pytest.approx(
        683.0582, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Gd"].value == pytest.approx(
        0.2566, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H"].value == pytest.approx(
        1.7756, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "H2O"].value == pytest.approx(
        1000000.0, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(694.61357, 1e-4)
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "SO4"].value == pytest.approx(
        3961.7867, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "La"].value == pytest.approx(
        0.9767, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Nd"].value == pytest.approx(
        0.93787, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Pr"].value == pytest.approx(
        0.30038, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sc"].value == pytest.approx(
        0.03234, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Sm"].value == pytest.approx(
        0.0962, 1e-4
    )
    assert model.fs.leach.liquid_outlet.conc_mass_comp[0, "Y"].value == pytest.approx(
        0.12326, 1e-4
    )

    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.01, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(1.2735e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(2.6728e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(0.00030416, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(7.9807e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(2.874e-6, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(3.3607e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(0.00010515, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(0.00016549, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(3.7036e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(1.7334, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(1.6975e-5, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(2.1749e-5, 1e-4)
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
    ].value == pytest.approx(3.940238e-9, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(1.6122889e-8, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(350.9697, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(736.601, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(8382.571, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(219.943, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(79.2062, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(926.1761, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(2897.92095, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(4560.9097, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(1020.6912, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(36.6518, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(467.8087, 1e-4)
    assert model.fs.solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(599.4000, 1e-4)

    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.flow_vol[
        0
    ].value == pytest.approx(62.010, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(7.27414e-9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(1.541183e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(2.0916e-6, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(4.4823e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(9.34607e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(1.8722e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(8.7493e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(9.8243e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(2.2807e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.44604, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(9.5164e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(1.3912e-7, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Al"
    ].value == pytest.approx(0.40095, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ca"
    ].value == pytest.approx(0.8495, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Ce"
    ].value == pytest.approx(115.2912, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Dy"
    ].value == pytest.approx(2.47066, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Fe"
    ].value == pytest.approx(0.051509, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Gd"
    ].value == pytest.approx(10.3196, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H"
    ].value == pytest.approx(41.44, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "HSO4"
    ].value == pytest.approx(1.61229e-8, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "SO4"
    ].value == pytest.approx(3.9402e-9, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "La"
    ].value == pytest.approx(48.2264, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Nd"
    ].value == pytest.approx(54.1513, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Pr"
    ].value == pytest.approx(12.57097, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sc"
    ].value == pytest.approx(0.14255, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Sm"
    ].value == pytest.approx(5.24542, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp[
        0, "Y"
    ].value == pytest.approx(7.66856, 1e-4)

    assert model.fs.precipitator.cv_aqueous.properties_out[
        0
    ].flow_vol.value == pytest.approx(9, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Al"
    ].value == pytest.approx(0.39734, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ca"
    ].value == pytest.approx(0.67535, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Ce"
    ].value == pytest.approx(36.8124, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Cl"
    ].value == pytest.approx(1438.56, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Dy"
    ].value == pytest.approx(0.317233, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Fe"
    ].value == pytest.approx(0.05025, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Gd"
    ].value == pytest.approx(1.2373, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "H"
    ].value == pytest.approx(41.44, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "H2O"
    ].value == pytest.approx(1000000.0, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "HSO4"
    ].value == pytest.approx(1.6122e-8, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "SO4"
    ].value == pytest.approx(3.9402e-9, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "La"
    ].value == pytest.approx(23.38496, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Nd"
    ].value == pytest.approx(9.9909, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Pr"
    ].value == pytest.approx(2.7656, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sc"
    ].value == pytest.approx(0.09749, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Sm"
    ].value == pytest.approx(0.663545, 1e-4)
    assert model.fs.precipitator.cv_aqueous.properties_out[0].conc_mass_comp[
        "Y"
    ].value == pytest.approx(1.95855, 1e-4)

    assert model.fs.precipitator.precipitate_outlet.temperature[
        0
    ].value == pytest.approx(348.15, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Al2(C2O4)3(s)"
    ].value == pytest.approx(6.018e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ca(C2O4)(s)"
    ].value == pytest.approx(3.9106e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Ce2(C2O4)3(s)"
    ].value == pytest.approx(0.0025204, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Dy2(C2O4)3(s)"
    ].value == pytest.approx(5.9633e-5, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Fe2(C2O4)3(s)"
    ].value == pytest.approx(1.01275e-7, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Gd2(C2O4)3(s)"
    ].value == pytest.approx(0.00025991, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "La2(C2O4)3(s)"
    ].value == pytest.approx(0.0008048, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Nd2(C2O4)3(s)"
    ].value == pytest.approx(1.3777e-3, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Pr2(C2O4)3(s)"
    ].value == pytest.approx(3.13144e-4, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sc2(C2O4)3(s)"
    ].value == pytest.approx(4.5114e-6, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Sm2(C2O4)3(s)"
    ].value == pytest.approx(1.3713e-4, 1e-4)
    assert model.fs.precipitator.precipitate_outlet.flow_mol_comp[
        0, "Y2(C2O4)3(s)"
    ].value == pytest.approx(2.890e-4, 1e-4)

    assert model.fs.roaster.gas_outlet.flow_mol[0].value == pytest.approx(
        0.014778, 1e-4
    )
    assert model.fs.roaster.gas_outlet.temperature[0].value == pytest.approx(
        873.15, 1e-4
    )
    assert model.fs.roaster.gas_outlet.pressure[0].value == pytest.approx(101325, 1e-4)
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "CO2"].value == pytest.approx(
        0.02343, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "H2O"].value == pytest.approx(
        0.524126, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "N2"].value == pytest.approx(
        0.39352, 1e-4
    )
    assert model.fs.roaster.gas_outlet.mole_frac_comp[0, "O2"].value == pytest.approx(
        0.05892, 1e-4
    )

    assert model.fs.leach_mixer.outlet.flow_vol[0].value == pytest.approx(
        625.688095, 1e-4
    )
    assert model.fs.rougher_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_rougher_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(445.8968, 1e-4)
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
    assert model.fs.load_sep.inlet.flow_vol[0].value == pytest.approx(445.89677, 1e-4)
    assert model.fs.scrub_sep.inlet.flow_vol[0].value == pytest.approx(0.09, 1e-4)
    assert model.fs.rougher_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.sc_circuit_purge.inlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(7.785, 1e-4)
    assert model.fs.solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol[
        0
    ].value == pytest.approx(9, 1e-4)
    assert model.fs.cleaner_org_make_up.outlet.flow_vol[0].value == pytest.approx(
        6.201, 1e-4
    )
    assert model.fs.cleaner_mixer.outlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.cleaner_sep.inlet.flow_vol[0].value == pytest.approx(62.01, 1e-4)
    assert model.fs.leach_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        445.89677, 1e-4
    )
    assert model.fs.cleaner_purge.inlet.flow_vol[0].value == pytest.approx(6.201, 1e-4)
    assert model.fs.sl_sep1.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        438.11177, 1e-4
    )
    assert model.fs.sl_sep2.recovered_liquid_outlet.flow_vol[0].value == pytest.approx(
        8.55, 1e-4
    )
    assert model.fs.precip_sep.inlet.flow_vol[0].value == pytest.approx(8.55, 1e-4)
    assert model.fs.precip_sx_mixer.outlet.flow_vol[0].value == pytest.approx(
        7.785, 1e-4
    )
    assert model.fs.precip_purge.inlet.flow_vol[0].value == pytest.approx(0.855, 1e-4)


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
        0.00019047, rel=1e-4
    )


@pytest.mark.unit
def test_display(system_frame):
    model = system_frame
    display_results(model)
    display_costing(model)
