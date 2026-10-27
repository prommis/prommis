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
from idaes.core.util.testing import assert_solution_equivalent
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
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.properties.mixed_acid_properties import (
    MixedAcidParameterBlock,
)
from prommis.uky.uky_flowsheet import (
    add_costing,
    add_result_expressions,
    build,
    data_reconciliation,
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
    assert isinstance(model.fs.leach_soln, MixedAcidParameterBlock)
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
    assert isinstance(model.fs.sx_soln, MixedAcidParameterBlock)
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
    assert isinstance(model.fs.sx_cleaner_strip_aq_translator, Arc)
    assert isinstance(model.fs.precip_aq_inlet, Arc)
    assert isinstance(model.fs.precip_solid_outlet, Arc)
    assert isinstance(model.fs.precip_aq_sl_sep2, Arc)
    assert isinstance(model.fs.sl_sep2_solid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_liquid_outlet, Arc)
    assert isinstance(model.fs.sl_sep2_retained_liquid_roaster, Arc)
    assert isinstance(model.fs.sl_sep2_aq_purge, Arc)
    assert isinstance(model.fs.precip_sx_mixer_translator, Arc)
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
    expected_results["leach.solid_outlet.flow_mass"] = {0: (9.98692, tol, None)}
    expected_results["leach.solid_outlet.mass_frac_comp"] = {
        (0, "Al2O3"): (2.00334e-01, tol, None),
        (0, "CaO"): (8.51935e-4, tol, None),
        (0, "Ce2O3"): (3.09813e-5, tol, None),
        (0, "Dy2O3"): (6.459123e-6, None, tol),
        (0, "Fe2O3"): (3.29703e-2, tol, None),
        (0, "Gd2O3"): (3.71286e-6, None, tol),
        (0, "La2O3"): (2.08603e-5, tol, None),
        (0, "Nd2O3"): (8.74798e-6, tol, None),
        (0, "Pr2O3"): (2.85399e-06, tol, None),
        (0, "Sc2O3"): (2.4983e-5, tol, None),
        (0, "Sm2O3"): (7.3565e-6, tol, None),
        (0, "Y2O3"): (2.23177e-5, tol, None),
        (0, "inerts"): (0.765721, tol, None),
    }
    expected_results["leach.liquid_outlet.flow_vol"] = {0: (286.891, tol, None)}
    expected_results["leach.liquid_outlet.conc_mass_comp"] = {
        (0, "Al_3+"): (2984.03, tol, None),
        (0, "Ca_2+"): (184.801, tol, None),
        (0, "Ce_3+"): (5.54825, tol, None),
        (0, "Cl_-"): (1900.08, tol, None),
        (0, "Dy_3+"): (0.125034, tol, None),
        (0, "Fe_3+"): (2410.41, tol, None),
        (0, "Gd_3+"): (0.39232, tol, None),
        (0, "H_+"): (4.76635, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "HSO4_-"): (7726.15, tol, None),
        (0, "SO4_2-"): (1.65466e4, tol, None),
        (0, "La_3+"): (4.07358, tol, None),
        (0, "Nd_3+"): (3.41841, tol, None),
        (0, "Pr_3+"): (1.23246, tol, None),
        (0, "Sc_3+"): (0.129667, tol, None),
        (0, "Y_3+"): (0.38354, tol, None),
    }

    # Solex Rougher Strip Results
    expected_results["solex_rougher_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (128.9, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (8.79553e1, tol, None),
        (0, "Ce_o"): (2.01508, tol, None),
        (0, "Dy_o"): (8.49862e-02, tol, None),
        (0, "Fe_o"): (2.90217e02, tol, None),
        (0, "Gd_o"): (1.61545e-01, tol, None),
        (0, "La_o"): (1.17401e-01, tol, None),
        (0, "Nd_o"): (6.0895e-01, tol, None),
        (0, "Pr_o"): (7.69972e-02, tol, None),
        (0, "Sc_o"): (1.97423e00, tol, None),
        (0, "Sm_o"): (1.00326e-02, tol, None),
        (0, "Y_o"): (3.23259, tol, None),
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (3.37500e00, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "H2O"): (1000000.0, tol, None),
        (0, "H_+"): (5.82914e02, tol, None),
        (0, "Al_3+"): (2.47317e3, tol, None),
        (0, "Ca_2+"): (1.60318e2, tol, None),
        (0, "Ce_3+"): (5.00699e02, tol, None),
        (0, "Cl_-"): (3.54530e04, tol, None),
        (0, "Dy_3+"): (1.7572e01, tol, None),
        (0, "Fe_3+"): (2.27948e03, tol, None),
        (0, "Gd_3+"): (4.24106e01, tol, None),
        (0, "La_3+"): (2.10862e01, tol, None),
        (0, "Nd_3+"): (1.49683e02, tol, None),
        (0, "Pr_3+"): (3.39757e00, tol, None),
        (0, "Sc_3+"): (3.12145e-03, tol, None),
        (0, "Sm_3+"): (1.80673e00, tol, None),
        (0, "Y_3+"): (2.21579e01, tol, None),
    }

    # Solex Cleaner Strip Results
    expected_results["solex_cleaner_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (6.03300e02, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (4.60449e01, tol, None),
        (0, "Ca_o"): (3.24556e00, tol, None),
        (0, "Ce_o"): (8.92687e-01, tol, None),
        (0, "Dy_o"): (1.18607e-02, tol, None),
        (0, "Fe_o"): (7.4743e01, tol, None),
        (0, "Gd_o"): (6.0861e-02, tol, None),
        (0, "La_o"): (7.97732e-02, tol, None),
        (0, "Nd_o"): (2.87135e-01, tol, None),
        (0, "Pr_o"): (3.51405e-02, tol, None),
        # (0, "Sc_o"): (1.746e-04, tol, None),  # This value fluctuates on Linux
        (0, "Sm_o"): (6.39721e-03, tol, None),
        (0, "Y_o"): (2.17061e-01, tol, None),
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (3.51700e00, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "Al_3+"): (9.72786e02, tol, None),
        (0, "Ca_2+"): (5.40818e01, tol, None),
        (0, "Ce_3+"): (2.30808e02, tol, None),
        (0, "Cl_-"): (3.54530e04, tol, None),
        (0, "Dy_3+"): (6.77404, tol, None),
        (0, "Fe_3+"): (5.19847e02, tol, None),
        (0, "Gd_3+"): (1.79659e01, tol, None),
        (0, "H_+"): (8.60741e02, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "La_3+"): (1.03905e01, tol, None),
        (0, "Nd_3+"): (6.60005e01, tol, None),
        (0, "Pr_3+"): (1.21321e00, tol, None),
        (0, "Sc_3+"): (6.993204e-6, None, tol),
        (0, "Sm_3+"): (7.71446e-1, tol, None),
        (0, "Y_3+"): (7.63437e00, tol, None),
    }

    # Precipitator Results
    expected_results["precipitator.cv_aqueous.properties_out[0].flow_vol"] = {
        None: (3.51700e00, tol, None)
    }
    expected_results["precipitator.cv_aqueous.properties_out[0].conc_mass_comp"] = {
        "Al_3+": (9.64031e02, tol, None),
        "Ca_2+": (4.2995e01, tol, None),
        "Ce_3+": (7.3697e01, tol, None),
        "Cl_-": (3.54530e04, tol, None),
        "Dy_3+": (8.69787e-1, tol, None),
        "Fe_3+": (5.07162e02, tol, None),
        "Gd_3+": (2.15411e00, tol, None),
        "H_+": (8.60741e02, tol, None),
        "H2O": (1000000.0, tol, None),
        "La_3+": (5.03834, tol, None),
        "Nd_3+": (1.21771e01, tol, None),
        "Pr_3+": (2.66907e-01, tol, None),
        "Sc_3+": (4.78265e-6, None, tol),
        "Sm_3+": (9.75879e-02, tol, None),
        "Y_3+": (1.94982e00, tol, None),
    }
    expected_results["precipitator.precipitate_outlet.temperature"] = {
        0: (348.15, tol, None)
    }
    expected_results["precipitator.precipitate_outlet.flow_mol_comp"] = {
        (0, "Al2(C2O4)3(s)"): (0.00053, None, tol),
        (0, "Ca(C2O4)(s)"): (9.72906e-04, tol, None),
        (0, "Ce2(C2O4)3(s)"): (0.00197, None, tol),
        (0, "Dy2(C2O4)3(s)"): (9.2872e-6, None, tol),
        (0, "Fe2(C2O4)3(s)"): (3.99414e-04, tol, None),
        (0, "Gd2(C2O4)3(s)"): (1.7682e-04, tol, None),
        (0, "La2(C2O4)3(s)"): (0.000153, None, tol),
        (0, "Nd2(C2O4)3(s)"): (0.00066, None, tol),
        (0, "Pr2(C2O4)3(s)"): (2.77874e-7, None, tol),
        (0, "Sc2(C2O4)3(s)"): (3.11363e-18, None, 1e-8),
        (0, "Sm2(C2O4)3(s)"): (1.20288e-6, None, tol),
        (0, "Y2(C2O4)3(s)"): (1.12438e-04, tol, None),
    }

    # Roaster Results
    expected_results["roaster.gas_outlet.flow_mol"] = {0: (1.32476e-02, tol, None)}
    expected_results["roaster.gas_outlet.temperature"] = {0: (873.15, tol, None)}
    expected_results["roaster.gas_outlet.pressure"] = {0: (101325, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "CO2"): (2.59581e-02, tol, None),
        (0, "H2O"): (4.69296e-01, tol, None),
        (0, "N2"): (4.38973e-01, tol, None),
        (0, "O2"): (6.57735e-02, tol, None),
    }

    # Volumetric flow rates
    expected_results["leach_mixer.outlet.flow_vol"] = {0: (2.86433e02, tol, None)}
    expected_results["rougher_org_make_up.outlet.flow_vol"] = {
        0: (1.28900e01, tol, None)
    }
    expected_results["solex_rougher_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (2.07048e02, tol, None)
    }
    expected_results["solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (1.00000e-01, tol, None)
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (3.37500e00, tol, None)
    }
    expected_results["acid_feed1.outlet.flow_vol"] = {0: (1.00000e-01, tol, None)}
    expected_results["acid_feed2.outlet.flow_vol"] = {0: (3.37500e00, tol, None)}
    expected_results["acid_feed3.outlet.flow_vol"] = {0: (3.51700e00, tol, None)}
    expected_results["rougher_sep.inlet.flow_vol"] = {0: (1.28900e02, tol, None)}
    expected_results["load_sep.inlet.flow_vol"] = {0: (2.07048e02, tol, None)}
    expected_results["scrub_sep.inlet.flow_vol"] = {0: (1.00000e-01, tol, None)}
    expected_results["rougher_mixer.outlet.flow_vol"] = {0: (1.28900e02, tol, None)}
    expected_results["rougher_organic_purge.inlet.flow_vol"] = {
        0: (1.28900e01, tol, None)
    }
    expected_results["solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (6.22377e00, tol, None)
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (3.51700e00, tol, None)
    }
    expected_results["cleaner_org_make_up.outlet.flow_vol"] = {
        0: (6.03300e01, tol, None)
    }
    expected_results["cleaner_mixer.outlet.flow_vol"] = {0: (6.03300e02, tol, None)}
    expected_results["cleaner_sep.inlet.flow_vol"] = {0: (6.03300e02, tol, None)}
    expected_results["leach_sx_mixer.outlet.flow_vol"] = {0: (2.07048e02, tol, None)}
    expected_results["cleaner_organic_purge.inlet.flow_vol"] = {
        0: (6.03300e01, tol, None)
    }
    expected_results["sl_sep1.recovered_liquid_outlet.flow_vol"] = {
        0: (2.00824e02, tol, None)
    }
    expected_results["sl_sep2.recovered_liquid_outlet.flow_vol"] = {
        0: (3.16530e00, tol, None)
    }
    expected_results["precip_sep.inlet.flow_vol"] = {0: (3.16530e00, tol, None)}
    expected_results["precip_sx_mixer.outlet.flow_vol"] = {0: (6.22377e00, tol, None)}
    expected_results["precip_purge.inlet.flow_vol"] = {0: (3.16530e-01, tol, None)}

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
        ("Y2O3", "Y_3+", "Y_o"),
        ("La2O3", "La_3+", "La_o"),
        ("Ce2O3", "Ce_3+", "Ce_o"),
        ("Pr2O3", "Pr_3+", "Pr_o"),
        ("Nd2O3", "Nd_3+", "Nd_o"),
        ("Sm2O3", "Sm_3+", "Sm_o"),
        ("Gd2O3", "Gd_3+", "Gd_o"),
        ("Dy2O3", "Dy_3+", "Dy_o"),
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
        "costing.total_plant_cost": {None: (7.27569e-01, tol, None)},
        "costing.total_BEC": {None: (2.44971e-01, tol, None)},
        "costing.total_installation_cost": {None: (4.82594e-01, tol, None)},
        "costing.other_plant_costs": {None: (3.59344e-06, tol, None)},
        "costing.total_fixed_OM_cost": {None: (6.80364e00, tol, None)},
        "costing.total_variable_OM_cost": {0: (1.36139e00, tol, None)},
        "costing.total_sales_revenue": {None: (3.3e-4, None, tol)},
        "costing.land_cost": {None: (2.96991e-5, tol, None)},
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


# Smoke tests to make sure data reconciliation and optimization solve
@pytest.mark.integration
@pytest.mark.solver
def test_data_reconciliation(system_frame):
    data_reconciliation(system_frame)


@pytest.mark.integration
@pytest.mark.solver
def test_optimize_model(system_frame):
    optimize_model(system_frame)
