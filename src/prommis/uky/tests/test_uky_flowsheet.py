#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for UKy flowsheet.

"""

from pyomo.environ import assert_optimal_termination, units, value
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

from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.leaching.leach_train import LeachingTrain
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.properties.hcl_stripping_properties import HClStrippingParameterBlock
from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.uky.uky_flowsheet import (
    add_costing,
    add_result_expressions,
    build,
    data_reconcilliation,
    display_costing,
    display_results,
    fix_organic_recycle,
    initialize_costing,
    initialize_system,
    optimize_model,
    set_operating_conditions,
    set_scaling,
    solve_system,
)
from prommis.util import assert_solution_equivalent


@pytest.fixture(scope="module")
def system_frame():
    m = build()
    set_operating_conditions(m)
    add_result_expressions(m)

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
    assert isinstance(model.fs.leach_soln, SulfuricAcidLeachingParameters)
    assert isinstance(model.fs.coal, CoalRefuseParameters)
    assert isinstance(model.fs.leach_rxns, CoalRefuseLeachingReactionParameterBlock)

    assert isinstance(model.fs.leach, LeachingTrain)
    assert isinstance(model.fs.sl_sep1, SLSeparator)
    assert isinstance(model.fs.leach_mixer, Mixer)
    assert isinstance(model.fs.leach_liquid_feed, Feed)
    assert isinstance(model.fs.leach_solid_feed, Feed)
    assert isinstance(model.fs.leach_filter_cake, Product)
    assert isinstance(model.fs.leach_filter_cake_liquid, Product)

    # Solvent extraction section property packages and unit models
    assert isinstance(model.fs.HCl_stripping_params, HClStrippingParameterBlock)
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
    assert isinstance(model.fs.rougher_organic_purge, Product)
    assert isinstance(model.fs.cleaner_organic_purge, Product)
    assert isinstance(model.fs.cleaner_org_make_up, Feed)
    assert isinstance(model.fs.rougher_org_make_up, Feed)
    assert isinstance(model.fs.acid_feed1, Feed)
    assert isinstance(model.fs.acid_feed2, Feed)
    assert isinstance(model.fs.acid_feed3, Feed)

    # Precipitation property packages and unit models

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
    assert isinstance(model.fs.sx_rougher_scrub_aq_translator, Arc)
    assert isinstance(model.fs.translator_scrub_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_scrub_org_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_strip_acid_feed, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_outlet, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_purge, Arc)
    assert isinstance(model.fs.sx_rougher_strip_org_recycle, Arc)
    assert isinstance(model.fs.sx_rougher_strip_aq_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_load_aq_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_org_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_mixed_org_recycle, Arc)
    assert isinstance(model.fs.sx_cleaner_load_aq_outlet_translator, Arc)
    assert isinstance(model.fs.sx_cleaner_load_translator_leach_sx_mixer, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_acid_feed, Arc)
    assert isinstance(model.fs.sx_cleaner_load_org_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_outlet, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_purge, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_org_recycle, Arc)
    assert isinstance(model.fs.sx_cleaner_strip_aq_precip, Arc)
    assert isinstance(model.fs.precip_solid_outlet, Arc)
    assert isinstance(model.fs.precip_aq_sl_sep2, Arc)
    assert isinstance(model.fs.sl_sep2_solid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_liquid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_aq_purge, Arc)
    assert isinstance(model.fs.sl_sep2_aq_recycle, Arc)


@pytest.mark.component
@pytest.mark.solver
def test_solve(system_frame):
    model = system_frame

    set_scaling(model)

    initialize_system(model)

    solve_system(model, tee=True)

    fix_organic_recycle(model)

    results = solve_system(model)

    assert_optimal_termination(results)


@pytest.mark.component
@pytest.mark.solver
def test_solution(system_frame):
    model = system_frame
    tol = 1e-4

    expected_results = {}

    # --- Data Definition ---

    # Leach Unit Results
    expected_results["leach.solid_outlet.flow_mass"] = {0: (2.22213e01, tol, None)}
    expected_results["leach.solid_outlet.mass_frac_comp"] = {
        (0, "Al2O3"): (0.233339, tol, None),
        (0, "CaO"): (0.00176827, tol, None),
        (0, "Ce2O3"): (9.47052e-5, tol, None),
        (0, "Dy2O3"): (6.459123e-6, None, tol),
        (0, "Fe2O3"): (0.0550767, tol, None),
        (0, "Gd2O3"): (3.71286e-6, None, tol),
        (0, "La2O3"): (4.46256e-5, tol, None),
        (0, "Nd2O3"): (3.76359e-5, tol, None),
        (0, "Pr2O3"): (9.31775e-06, tol, None),
        (0, "Sc2O3"): (2.69415e-5, tol, None),
        (0, "Sm2O3"): (1.18787e-5, tol, None),
        (0, "Y2O3"): (2.90205e-5, tol, None),
        (0, "inerts"): (0.709552, tol, None),
    }
    expected_results["leach.liquid_outlet.flow_vol"] = {0: (607.327, tol, None)}
    expected_results["leach.liquid_outlet.conc_mass_comp"] = {
        (0, "Al"): (446.962, tol, None),
        (0, "Ca"): (113.572, tol, None),
        (0, "Ce"): (5.39646, tol, None),
        (0, "Cl"): (124.023, tol, None),
        (0, "Dy"): (0.105014, tol, None),
        (0, "Fe"): (719.756, tol, None),
        (0, "Gd"): (0.475374, tol, None),
        (0, "H"): (1.75316, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "HSO4"): (708.093, tol, None),
        (0, "SO4"): (4.09043e03, tol, None),
        (0, "La"): (2.04625, tol, None),
        (0, "Nd"): (2.6378, tol, None),
        (0, "Pr"): (0.686617, tol, None),
        (0, "Sc"): (0.036107, tol, None),
        (0, "Y"): (0.19884, tol, None),
    }

    # Solex Rougher Strip Results
    expected_results["solex_rougher_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (62.01, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (2.21854e01, tol, None),
        (0, "Ce_o"): (1.00019e00, tol, None),
        (0, "Dy_o"): (9.46796e-03, tol, None),
        (0, "Fe_o"): (1.05113e02, tol, None),
        (0, "Gd_o"): (4.20017e-01, tol, None),
        (0, "La_o"): (3.75435e-01, tol, None),
        (0, "Nd_o"): (4.34448e-01, tol, None),
        (0, "Pr_o"): (1.16491e-01, tol, None),
        (0, "Sc_o"): (2.25521e00, tol, None),
        (0, "Sm_o"): (8.91130e-02, tol, None),
        (0, "Y_o"): (1.16378e-01, tol, None),
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (1.30000e-01, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "H2O"): (1000000.0, tol, None),
        (0, "H"): (4.84108e03, tol, None),
        (0, "Al"): (4.34299e02, tol, None),
        (0, "Ca"): (1.10980e02, tol, None),
        (0, "Ce"): (2.71989e02, tol, None),
        (0, "Cl"): (1.79700e05, tol, None),
        (0, "Dy"): (9.39532e02, tol, None),
        (0, "Fe"): (7.12470e02, tol, None),
        (0, "Gd"): (1.23220e03, tol, None),
        (0, "La"): (1.06115e02, tol, None),
        (0, "Nd"): (9.56514e01, tol, None),
        (0, "Pr"): (6.10344e00, tol, None),
        (0, "Sc"): (3.56558e-03, tol, None),
        (0, "Sm"): (4.51837e01, tol, None),
        (0, "Y"): (3.63620e03, tol, None),
    }

    # Solex Cleaner Strip Results
    expected_results["solex_cleaner_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (6.20100e01, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (6.32280e00, tol, None),
        (0, "Ca_o"): (1.70227e00, tol, None),
        (0, "Ce_o"): (8.77313e-01, tol, None),
        (0, "Dy_o"): (1.98951e-03, tol, None),
        (0, "Fe_o"): (1.29152e01, tol, None),
        (0, "Gd_o"): (5.74821e-01, tol, None),
        (0, "La_o"): (3.52850e-01, tol, None),
        (0, "Nd_o"): (3.52635e-01, tol, None),
        (0, "Pr_o"): (5.63599e-02, tol, None),
        (0, "Sc_o"): (7.47573e-05, tol, None),
        (0, "Sm_o"): (9.54875e-02, tol, None),
        (0, "Y_o"): (5.77284e-02, tol, None),
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (0.03, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "Al"): (1.20168e02, tol, None),
        (0, "Ca"): (2.60620e01, tol, None),
        (0, "Ce"): (1.98303e02, tol, None),
        (0, "Cl"): (1.79700e05, tol, None),
        (0, "Dy"): (2.35258e02, tol, None),
        (0, "Fe"): (8.66114e01, tol, None),
        (0, "Gd"): (9.88843e02, tol, None),
        (0, "H"): (5.00580e03, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "La"): (7.90549e01, tol, None),
        (0, "Nd"): (6.60411e01, tol, None),
        (0, "Pr"): (2.77363e00, tol, None),
        (0, "Sc"): (6.993204e-6, None, tol),
        (0, "Sm"): (3.44977e01, tol, None),
        (0, "Y"): (1.56414e03, tol, None),
    }

    # Precipitator Results
    expected_results["precipitator.cv_aqueous.properties_out[0].flow_vol"] = {
        None: (0.03, tol, None)
    }
    expected_results["precipitator.cv_aqueous.properties_out[0].conc_mass_comp"] = {
        "Al": (1.19087e02, tol, None),
        "Ca": (2.07193e01, tol, None),
        "Ce": (6.33183e01, tol, None),
        "Cl": (1.79700e05, tol, None),
        "Dy": (3.02072e01, tol, None),
        "Fe": (8.44981e01, tol, None),
        "Gd": (1.18562e02, tol, None),
        "H": (5.00580e03, tol, None),
        "H2O": (1000000.0, tol, None),
        "La": (3.83337e01, tol, None),
        "Nd": (1.21846e01, tol, None),
        "Pr": (6.10199e-01, tol, None),
        "Sc": (4.78265e-6, None, tol),
        "Sm": (4.36397e00, tol, None),
        "Y": (3.99481e02, tol, None),
    }
    expected_results["precipitator.precipitate_outlet.temperature"] = {
        0: (348.15, tol, None)
    }
    expected_results["precipitator.precipitate_outlet.flow_mol_comp"] = {
        (0, "Al2(C2O4)3(s)"): (1.3049e-6, None, tol),
        (0, "Ca(C2O4)(s)"): (3.99923e-06, tol, None),
        (0, "Ce2(C2O4)3(s)"): (1.6235e-6, None, tol),
        (0, "Dy2(C2O4)3(s)"): (9.2872e-6, None, tol),
        (0, "Fe2(C2O4)3(s)"): (5.67638e-07, tol, None),
        (0, "Gd2(C2O4)3(s)"): (8.30157e-05, tol, None),
        (0, "La2(C2O4)3(s)"): (9.30326e-7, None, tol),
        (0, "Nd2(C2O4)3(s)"): (1.0331e-6, None, tol),
        (0, "Pr2(C2O4)3(s)"): (2.77874e-7, None, tol),
        (0, "Sc2(C2O4)3(s)"): (1.24687e-14, tol, None),
        (0, "Sm2(C2O4)3(s)"): (1.20288e-6, None, tol),
        (0, "Y2(C2O4)3(s)"): (1.96500e-04, tol, None),
    }

    # Roaster Results
    expected_results["roaster.gas_outlet.flow_mol"] = {0: (7.85761e-03, tol, None)}
    expected_results["roaster.gas_outlet.temperature"] = {0: (873.15, tol, None)}
    expected_results["roaster.gas_outlet.pressure"] = {0: (101325, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "CO2"): (4.29085e-02, tol, None),
        (0, "H2O"): (1.05898e-01, tol, None),
        (0, "N2"): (7.40088e-01, tol, None),
        (0, "O2"): (1.11105e-01, tol, None),
    }

    # Volumetric flow rates
    expected_results["leach_mixer.outlet.flow_vol"] = {0: (6.07136e02, tol, None)}
    expected_results["rougher_org_make_up.outlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["solex_rougher_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (4.25283e02, tol, None)
    }
    expected_results["solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (9.00000e-02, tol, None)
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (1.30000e-01, tol, None)
    }
    expected_results["acid_feed1.outlet.flow_vol"] = {0: (9.00000e-02, tol, None)}
    expected_results["acid_feed2.outlet.flow_vol"] = {0: (1.30000e-01, tol, None)}
    expected_results["acid_feed3.outlet.flow_vol"] = {0: (3.00000e-02, tol, None)}
    expected_results["rougher_sep.inlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["load_sep.inlet.flow_vol"] = {0: (4.25283e02, tol, None)}
    expected_results["scrub_sep.inlet.flow_vol"] = {0: (9.00000e-02, tol, None)}
    expected_results["rougher_mixer.outlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["rougher_organic_purge.inlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (1.54300e-01, tol, None)
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (3.00000e-02, tol, None)
    }
    expected_results["cleaner_org_make_up.outlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["cleaner_mixer.outlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["cleaner_sep.inlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["leach_sx_mixer.outlet.flow_vol"] = {0: (4.25283e02, tol, None)}
    expected_results["cleaner_organic_purge.inlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["sl_sep1.recovered_liquid_outlet.flow_vol"] = {
        0: (4.25129e02, tol, None)
    }
    expected_results["sl_sep2.recovered_liquid_outlet.flow_vol"] = {
        0: (2.70000e-02, tol, None)
    }
    expected_results["precip_sep.inlet.flow_vol"] = {0: (2.70000e-02, tol, None)}
    expected_results["precip_sx_mixer.outlet.flow_vol"] = {0: (1.54300e-01, tol, None)}
    expected_results["precip_purge.inlet.flow_vol"] = {0: (2.70000e-03, tol, None)}

    assert_solution_equivalent(model.fs, expected_results)


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
            rel=1e-4,
        )


@pytest.mark.component
@pytest.mark.solver
def test_costing(system_frame):
    model = system_frame
    add_costing(model)
    initialize_costing(model)


@pytest.mark.component
@pytest.mark.solver
def test_costing_diagnostics(system_frame):
    model = system_frame
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()


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
    tol = 1e-4

    expected_results = {
        "costing.total_plant_cost": {None: (8.32020e-01, tol, None)},
        "costing.total_BEC": {None: (2.80140e-01, tol, None)},
        "costing.total_installation_cost": {None: (5.51875e-01, tol, None)},
        "costing.other_plant_costs": {None: (4.78618e-06, tol, None)},
        "costing.total_fixed_OM_cost": {None: (6.80612e00, tol, None)},
        "costing.total_variable_OM_cost": {0: (1.36131e00, tol, None)},
        "costing.total_sales_revenue": {None: (2.69422e-5, None, tol)},
        "costing.land_cost": {None: (6.1234e-5, tol, None)},
    }
    assert_solution_equivalent(model.fs, expected_results)


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


# Smoke tests to make sure data reconcilliation and optimization solve
@pytest.mark.integration
@pytest.mark.solver
def test_data_reconcilliation(system_frame):
    data_reconcilliation(system_frame)


@pytest.mark.integration
@pytest.mark.solver
def test_optimize_model(system_frame):
    optimize_model(system_frame)
