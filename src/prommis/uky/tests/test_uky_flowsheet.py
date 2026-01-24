#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests for UKy flowsheet.

"""
import pytest

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

from prommis.util import assert_solution_equivalent
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
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
    tol = 1e-4

    expected_results = {}

    # --- Data Definition ---

    # Leach Unit Results
    expected_results["leach.solid_outlet.flow_mass"] = {0: (22.238837, tol, None)}
    expected_results["leach.solid_outlet.mass_frac_comp"] = {
        (0, "Al2O3"): (0.232736, tol, None),
        (0, "CaO"): (0.0020873, tol, None),
        (0, "Ce2O3"): (1.027533e-4, tol, None),
        (0, "Dy2O3"): (6.459123e-6, None, tol),
        (0, "Fe2O3"): (0.05590266, tol, None),
        (0, "Gd2O3"): (3.71286e-6, None, tol),
        (0, "La2O3"): (4.8569e-5, tol, None),
        (0, "Nd2O3"): (4.07967e-5, tol, None),
        (0, "Pr2O3"): (1.043202e-5, tol, None),
        (0, "Sc2O3"): (2.70925995e-5, tol, None),
        (0, "Sm2O3"): (1.234895e-5, tol, None),
        (0, "Y2O3"): (2.969969e-5, tol, None),
        (0, "inerts"): (0.708991, tol, None),
    }
    expected_results["leach.liquid_outlet.flow_vol"] = {0: (865.2653, tol, None)}
    expected_results["leach.liquid_outlet.conc_mass_comp"] = {
        (0, "Al"): (329.3326, tol, None),
        (0, "Ca"): (63.95648, tol, None),
        (0, "Ce"): (3.349581, tol, None),
        (0, "Cl"): (156.86977, tol, None),
        (0, "Dy"): (0.03894, tol, None),
        (0, "Fe"): (464.83, tol, None),
        (0, "Gd"): (0.394755, tol, None),
        (0, "H"): (2.58245, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "HSO4"): (684.7322, tol, None),
        (0, "SO4"): (2685.26996, tol, None),
        (0, "La"): (1.216377, tol, None),
        (0, "Nd"): (1.67917, tol, None),
        (0, "Pr"): (0.418309, tol, None),
        (0, "Sc"): (0.023101, tol, None),
        (0, "Y"): (0.097377, tol, None),
    }

    # Solex Rougher Strip Results
    expected_results["solex_rougher_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (62.01, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (0.053654, tol, None),
        (0, "Ce_o"): (0.0063643, tol, None),
        (0, "Dy_o"): (1.19399, tol, None),
        (0, "Fe_o"): (2.17512, tol, None),
        (0, "Gd_o"): (0.15928, tol, None),
        (0, "La_o"): (0.0045008, tol, None),
        (0, "Nd_o"): (0.0037576, tol, None),
        (0, "Pr_o"): (0.0011516, tol, None),
        (0, "Sc_o"): (1.9319, tol, None),
        (0, "Sm_o"): (0.0054847, tol, None),
        (0, "Y_o"): (4.6218, tol, None),
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (9.0, tol, None)
    }
    expected_results[
        "solex_rougher_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "H2O"): (1000000.0, tol, None),
        (0, "H"): (39.28586, tol, None),
        (0, "SO4"): (2.41459e-7, None, tol),
        (0, "HSO4"): (2.70217e-7, None, tol),
        (0, "Al"): (3.77224, tol, None),
        (0, "Ca"): (1.07311, tol, None),
        (0, "Ce"): (0.28459, tol, None),
        (0, "Cl"): (1438.56, tol, None),
        (0, "Dy"): (1.15087, tol, None),
        (0, "Fe"): (28.64549, tol, None),
        (0, "Gd"): (1.80357, tol, None),
        (0, "La"): (0.212375, tol, None),
        (0, "Nd"): (0.15231, tol, None),
        (0, "Pr"): (0.037087, tol, None),
        (0, "Sc"): (0.00305517, tol, None),
        (0, "Sm"): (0.14207, tol, None),
        (0, "Y"): (2.04473, tol, None),
    }

    # Solex Cleaner Strip Results
    expected_results["solex_cleaner_strip.mscontactor.organic_outlet.flow_vol"] = {
        0: (62.010, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.organic_outlet.conc_mass_comp"
    ] = {
        (0, "Al_o"): (0.0041203, tol, None),
        (0, "Ca_o"): (0.0024165, tol, None),
        (0, "Ce_o"): (0.0006731, tol, None),
        (0, "Dy_o"): (0.373165, tol, None),
        (0, "Fe_o"): (0.834055, tol, None),
        (0, "Gd_o"): (0.0455115, tol, None),
        (0, "La_o"): (0.00046228, tol, None),
        (0, "Nd_o"): (0.0004254, tol, None),
        (0, "Pr_o"): (0.00016404, tol, None),
        (0, "Sc_o"): (0.00441393, tol, None),
        (0, "Sm_o"): (0.00089498, tol, None),
        (0, "Y_o"): (1.15377, tol, None),
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (9, tol, None)
    }
    expected_results[
        "solex_cleaner_strip.mscontactor.aqueous_outlet.conc_mass_comp"
    ] = {
        (0, "Al"): (0.869364, tol, None),
        (0, "Ca"): (0.29638, tol, None),
        (0, "Ce"): (0.07427, tol, None),
        (0, "Cl"): (1438.56, tol, None),
        (0, "Dy"): (0.38478, tol, None),
        (0, "Fe"): (16.23179, tol, None),
        (0, "Gd"): (0.76294, tol, None),
        (0, "H"): (40.4119, tol, None),
        (0, "H2O"): (1000000.0, tol, None),
        (0, "HSO4"): (1.22141e-6, None, tol),
        (0, "SO4"): (1.13409e-6, None, tol),
        (0, "La"): (0.05575, tol, None),
        (0, "Nd"): (0.040607, tol, None),
        (0, "Pr"): (0.011155, tol, None),
        (0, "Sc"): (6.993204e-6, None, tol),
        (0, "Sm"): (0.04601, tol, None),
        (0, "Y"): (0.53946, tol, None),
    }

    # Precipitator Results
    expected_results["precipitator.cv_aqueous.properties_out[0].flow_vol"] = {
        None: (9, tol, None)
    }
    expected_results["precipitator.cv_aqueous.properties_out[0].conc_mass_comp"] = {
        "Al": (0.86154, tol, None),
        "Ca": (0.23563, tol, None),
        "Ce": (0.023715, tol, None),
        "Cl": (1438.56, tol, None),
        "Dy": (0.0494056, tol, None),
        "Fe": (15.835737, tol, None),
        "Gd": (0.091476, tol, None),
        "H": (40.41, tol, None),
        "H2O": (1000000.0, tol, None),
        "HSO4": (1.2214e-6, None, tol),
        "SO4": (1.134087e-6, None, tol),
        "La": (0.027033, tol, None),
        "Nd": (0.007492, tol, None),
        "Pr": (0.002454, tol, None),
        "Sc": (4.78265e-6, None, tol),
        "Sm": (0.0058206, tol, None),
        "Y": (0.137778, tol, None),
    }
    expected_results["precipitator.precipitate_outlet.temperature"] = {
        0: (348.15, tol, None)
    }
    expected_results["precipitator.precipitate_outlet.flow_mol_comp"] = {
        (0, "Al2(C2O4)3(s)"): (1.3049e-6, None, tol),
        (0, "Ca(C2O4)(s)"): (1.36438e-5, tol, None),
        (0, "Ce2(C2O4)3(s)"): (1.6235e-6, None, tol),
        (0, "Dy2(C2O4)3(s)"): (9.2872e-6, None, tol),
        (0, "Fe2(C2O4)3(s)"): (3.1914e-5, tol, None),
        (0, "Gd2(C2O4)3(s)"): (1.92151e-5, tol, None),
        (0, "La2(C2O4)3(s)"): (9.30326e-7, None, tol),
        (0, "Nd2(C2O4)3(s)"): (1.0331e-6, None, tol),
        (0, "Pr2(C2O4)3(s)"): (2.77874e-7, None, tol),
        (0, "Sc2(C2O4)3(s)"): (2.28046e-10, tol, None),
        (0, "Sm2(C2O4)3(s)"): (1.20288e-6, None, tol),
        (0, "Y2(C2O4)3(s)"): (2.03315e-5, tol, None),
    }

    # Roaster Results
    expected_results["roaster.gas_outlet.flow_mol"] = {0: (0.021699, tol, None)}
    expected_results["roaster.gas_outlet.temperature"] = {0: (873.15, tol, None)}
    expected_results["roaster.gas_outlet.pressure"] = {0: (101325, tol, None)}
    expected_results["roaster.gas_outlet.mole_frac_comp"] = {
        (0, "CO2"): (0.0155197, tol, None),
        (0, "H2O"): (0.676245, tol, None),
        (0, "N2"): (0.267998, tol, None),
        (0, "O2"): (0.040237, tol, None),
    }

    # Volumetric flow rates
    expected_results["leach_mixer.outlet.flow_vol"] = {0: (865.078, tol, None)}
    expected_results["rougher_org_make_up.outlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["solex_rougher_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (621.9757, tol, None)
    }
    expected_results["solex_rougher_scrub.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (90, tol, None)
    }
    expected_results["solex_rougher_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (9, tol, None)
    }
    expected_results["acid_feed1.outlet.flow_vol"] = {0: (90, tol, None)}
    expected_results["acid_feed2.outlet.flow_vol"] = {0: (9, tol, None)}
    expected_results["acid_feed3.outlet.flow_vol"] = {0: (9, tol, None)}
    expected_results["rougher_sep.inlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["load_sep.inlet.flow_vol"] = {0: (621.9757, tol, None)}
    expected_results["scrub_sep.inlet.flow_vol"] = {0: (90, tol, None)}
    expected_results["rougher_mixer.outlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["rougher_organic_purge.inlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["solex_cleaner_load.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (16.29, tol, None)
    }
    expected_results["solex_cleaner_strip.mscontactor.aqueous_outlet.flow_vol"] = {
        0: (9, tol, None)
    }
    expected_results["cleaner_org_make_up.outlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["cleaner_mixer.outlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["cleaner_sep.inlet.flow_vol"] = {0: (62.01, tol, None)}
    expected_results["leach_sx_mixer.outlet.flow_vol"] = {0: (621.97569, tol, None)}
    expected_results["cleaner_organic_purge.inlet.flow_vol"] = {0: (6.201, tol, None)}
    expected_results["sl_sep1.recovered_liquid_outlet.flow_vol"] = {
        0: (605.68569, tol, None)
    }
    expected_results["sl_sep2.recovered_liquid_outlet.flow_vol"] = {0: (8.1, tol, None)}
    expected_results["precip_sep.inlet.flow_vol"] = {0: (8.1, tol, None)}
    expected_results["precip_sx_mixer.outlet.flow_vol"] = {0: (16.29, tol, None)}
    expected_results["precip_purge.inlet.flow_vol"] = {0: (0.81, tol, None)}

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
    tol = 1e-4

    expected_results = {
        "costing.total_plant_cost": {None: (1.18071, tol, None)},
        "costing.total_BEC": {None: (0.39754, tol, None)},
        "costing.total_installation_cost": {None: (0.78316, tol, None)},
        "costing.other_plant_costs": {None: (7.5987e-06, tol, None)},
        "costing.total_fixed_OM_cost": {None: (6.8166, tol, None)},
        "costing.total_variable_OM_cost": {0: (1.3634, tol, None)},
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
