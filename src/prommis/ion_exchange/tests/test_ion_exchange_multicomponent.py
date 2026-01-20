#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This test is a modified version from
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/tests/test_ion_exchange_clark.py

# WaterTAP License Agreement

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import os
import pytest
import pandas as pd
import pyomo.environ as pyo

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.exceptions import ConfigurationError

pytest.importorskip("watertap", reason="WaterTAP dependency not available")

# pylint: disable=wrong-import-position
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from prommis.ion_exchange.ion_exchange_multicomponent import (
    IonExchangeMultiComp,
)
from prommis.ion_exchange.ix_freundlich_multicomponent_example import (
    main,
    build_model,
    add_data,
    build_clark,
    set_bounds,
    set_operating_conditions,
    set_scaling,
    add_costing,
    run_optimization,
)

"""
modified by: Soraya Rawlings
"""


@pytest.fixture(scope="module")
def m():

    path = os.path.dirname(os.path.realpath(__file__))
    resin_file = os.path.join(path, "..", "data", "resin_data.json")
    comp_prop_file = os.path.join(path, "..", "data", "properties_data.json")
    parmest_file = os.path.join(path, "..", "data", "parmest_data.json")
    curve_file = os.path.join(path, "..", "data", "breakthrough_literature_data.csv")
    curve_data = pd.read_csv(curve_file)

    solver = get_solver()

    resin = "S950"
    target_component = "La"
    num_traps = 30
    c_trap_min = 1e-3
    regenerant = "single_use"
    hazardous_waste = False

    # Add sets for solvent and ion species
    m = build_model()

    add_data(
        m,
        target_component=target_component,
        resin=resin,
        curve_data=curve_data,
        resin_file=resin_file,
        comp_prop_file=comp_prop_file,
        parmest_file=parmest_file,
    )

    parmest_data = m.fs.parmest_data
    build_clark(
        m,
        resin=resin,
        regenerant=regenerant,
        target_component=target_component,
        num_traps=num_traps,
        c_trap_min=c_trap_min,
        resin_file=resin_file,
        hazardous_waste=hazardous_waste,
    )

    set_bounds(m)

    set_operating_conditions(
        m, parmest_data=parmest_data, target_component=target_component, resin=resin
    )

    set_scaling(m)

    return m


def build_clark_with_costing(m, regenerant, target_component):

    add_costing(m, regenerant=regenerant, target_component=target_component)

    # Set values for variables needed during the test for staticmethod
    # in the costing model file.
    m.fs.costing.aggregate_capital_cost.set_value(1e3)
    m.fs.costing.aggregate_fixed_operating_cost.set_value(1e3)
    m.fs.costing.aggregate_variable_operating_cost.set_value(1e3)
    m.fs.costing.aggregate_flow_costs["electricity"].set_value(1e-3)
    if regenerant != "single_use":
        m.fs.costing.aggregate_flow_costs["NaCl"].set_value(1e3)
    m.fs.costing.total_capital_cost.set_value(1e3)
    m.fs.costing.total_operating_cost.set_value(1e3)
    m.fs.costing.initialize_build(m.fs.costing)

    return m


@pytest.mark.unit
def test_config_error_in_ix_type():

    # Set up the model with parameters that will trigger the
    # ConfigurationError
    with pytest.raises(
        ConfigurationError,
        match="The current ion exchange model is limited to cation exchange methods and alternative techniques are not addressed at this time.",
    ):

        path = os.path.dirname(os.path.realpath(__file__))
        resin_file = os.path.join(path, "..", "data", "resin_data.json")
        resin = "S950"
        target_component = "Cl"
        regenerant = "single_use"
        list_solvent = ["H2O"]
        list_reactive_ions = ["Cl"]
        hazardous_waste = False
        num_traps = 30
        c_trap_min = 1e-3

        # Add sets for solvent and ion species
        m = build_model()

        m.fs.set_solvent = pyo.Set(initialize=list_solvent)
        m.fs.set_reactive_ions = pyo.Set(initialize=list_reactive_ions)
        m.fs.set_all = pyo.Set(initialize=list_solvent + list_reactive_ions)

        ion_props = {
            "solute_list": [],
            "diffusivity_data": {},
            "molar_volume_data": {},
            "mw_data": {},  # in kg/mol
            "charge": {},
        }

        ion_props["solute_list"] = list_reactive_ions
        ion_props["diffusivity_data"] = {}
        ion_props["molar_volume_data"] = {
            ("Liq", "Cl"): 0.0006818,
        }
        ion_props["mw_data"] = {
            "H2O": 0.018,
            "Cl": 0.035453,
        }
        ion_props["charge"] = {"Cl": -1}
        ion_props["diffus_calculation"] = "HaydukLaudie"
        m.fs.properties = MCASParameterBlock(**ion_props)

        ix_config = {
            "property_package": m.fs.properties,
            "regenerant": regenerant,
            "target_component": target_component,
            "reactive_ions": list_reactive_ions,
            "number_traps": num_traps,
            "c_trap_min": c_trap_min,
            "resin_data_path": resin_file,
            "resin": resin,
            "hazardous_waste": hazardous_waste,
        }

        m.fs.unit = ix = IonExchangeMultiComp(**ix_config)


@pytest.mark.unit
def test_structural_issues(m):

    dt = DiagnosticsToolbox(m)
    dt.assert_no_structural_warnings()


@pytest.mark.component
def test_resin_specific_data(m):

    ix = m.fs.unit_ix
    resin = "S950"

    # Test that data coming from .json file is correct
    test_json_data = {
        "S950": {
            "bed_expansion": {"param_a": 0, "param_b": 0.0596, "param_c": -4e-5},
            "pressure_drop": {"param_a": 0, "param_b": 0.1668, "param_c": 2e-4},
        }
    }

    assert test_json_data[resin]["bed_expansion"]["param_a"] == pyo.value(
        ix.bed_expansion_frac_param_A
    )
    assert test_json_data[resin]["bed_expansion"]["param_b"] == pyo.value(
        ix.bed_expansion_frac_param_B
    )
    assert test_json_data[resin]["bed_expansion"]["param_c"] == pyo.value(
        ix.bed_expansion_frac_param_C
    )

    assert test_json_data[resin]["pressure_drop"]["param_a"] == pyo.value(
        ix.pressure_drop_param_A
    )
    assert test_json_data[resin]["pressure_drop"]["param_b"] == pyo.value(
        ix.pressure_drop_param_B
    )
    assert test_json_data[resin]["pressure_drop"]["param_c"] == pyo.value(
        ix.pressure_drop_param_C
    )


@pytest.mark.component
def test_initialization(m):

    # Scale model
    init_scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = init_scaling.create_using(m, rename=False)

    check_dof(m, fail_flag=True)

    # Solve scaled model with zero degrees of freedom
    solver = get_solver()
    scaled_results = solver.solve(scaled_model)
    pyo.assert_optimal_termination(scaled_results)

    # Propagate the solution back to the original model
    results = init_scaling.propagate_solution(scaled_model, m)

    # Test expected results for ix variables. NOTE: bed_expansion_frac
    # and pressure_drop are the ones calculated with resin-specific
    # data from .json file
    ix_vars_results = {
        "bed_depth": 1,
        "column_height": 2.2969999999999997,
        "bed_diameter": 0.015,
        "number_columns": 1.0,
        "loading_rate": 9.450285573574311e-05,
        "bed_volume_total": 0.00017671458676442585,
        "target_breakthrough_time": 1693392.3126900392,
        "service_flow_rate": 0.3402102806486752,
        "N_Re": 0.07087714180180733,
        "N_Pe_particle": 0.01403501814152805,
        "N_Pe_bed": 18.7133575220374,
    }

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit_ix, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)


@pytest.mark.component
def test_optimization_single_use(m):

    regenerant = "single_use"
    target_component = "La"
    solver = get_solver()

    # Scale model
    init_scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = init_scaling.create_using(m, rename=False)

    check_dof(m, fail_flag=True)

    # Solve scaled model with zero degrees of freedom
    scaled_results = solver.solve(scaled_model)
    pyo.assert_optimal_termination(scaled_results)

    # Propagate the solution back to the original model
    init_scaling.propagate_solution(scaled_model, m)

    build_clark_with_costing(
        m, regenerant=regenerant, target_component=target_component
    )

    check_dof(m, fail_flag=True)

    # Solve scaled model
    results = solver.solve(m)
    pyo.assert_optimal_termination(results)

    run_optimization(m, target_component=target_component)

    results_opt = solver.solve(m)
    pyo.assert_optimal_termination(results_opt)

    # Test expected results for ix variables. NOTE: bed_expansion_frac
    # and pressure_drop are the ones calculated with resin-specific
    # data from .json file
    ix_vars_results = {
        "resin_diam": 0.00075,
        "resin_density": 1126.61,
        "bed_volume": 0.00017671456829028136,
        "bed_volume_total": 0.00017671456829028136,
        "bed_depth": 2.249999764670581,
        "bed_porosity": 0.80,
        "column_height": 3.918249694777743,
        "bed_diameter": 0.010000000000242901,
        "number_columns": 1.0,
        "number_columns_redundant": 1.0,
        "target_breakthrough_time": 2291425.7780593033,
        "ebct": 10581.689703155758,
        "loading_rate": 0.00021263142539509238,
        "service_flow_rate": 0.34021031621503484,
        "N_Re": 0.1594735690463193,
        "N_Pe_particle": 0.02071383855568819,
        "N_Pe_bed": 62.141509167630446,
    }

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit_ix, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)

    # Expected results for system costs
    sys_cost_results = {
        "utilization_factor": 1,
        "TPEC": 4.121212121212121,
        "TIC": 2.0,
        "total_capital_cost": 2398.5398832175597,
        "total_operating_cost": 88.9915519881193,
        "electricity_cost": 0.07,
        "aggregate_capital_cost": 2398.5398832175597,
        "aggregate_fixed_operating_cost": 17.03535548053501,
        "aggregate_variable_operating_cost": 0.0,
        "aggregate_flow_electricity": 4.1445910422021544e-08,
        "aggregate_flow_costs": {
            "electricity": 2.5432039553160866e-05,
        },
        "maintenance_labor_chemical_operating_cost": 71.95619649459807,
    }

    for v, r in sys_cost_results.items():
        mv = getattr(m.fs.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)

    # Expected results for ion exchange costs considering "single_use"
    # for regenerant use
    ix_cost_results = {
        "capital_cost": 2398.539883153269,
        "fixed_operating_cost": 17.03535548053501,
        "capital_cost_vessel": 598.398018274306,
        "capital_cost_resin": 1.2369525140112791,
        "capital_cost_backwash_tank": 47.95379556887687,
        "operating_cost_hazardous": 0,
        "total_pumping_power": 4.1445910422021544e-08,
        "flow_vol_resin": 0.0024337195286942945,
        "single_use_resin_replacement_cost": 17.03535548053501,
        "cost_factor": 2.0,
    }

    for v, r in ix_cost_results.items():
        mv = getattr(m.fs.unit_ix.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)


@pytest.mark.component
def test_scaling(m):

    m.fs.unit_ix.calculate_scaling_factors()


@pytest.mark.unit
def test_get_stream_table_contents(m):

    stable = m.fs.unit_ix._get_stream_table_contents()

    expected = pd.DataFrame.from_dict(
        {
            "Units": {
                "flow_mol_phase_comp ('Liq', 'H2O')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'La')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Dy')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Ho')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Er')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Yb')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Sm')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "temperature": getattr(pyo.units.pint_registry, "kelvin"),
                "pressure": getattr(pyo.units.pint_registry, "pascal"),
            },
            "Feed Inlet": {
                "flow_mol_phase_comp ('Liq', 'H2O')": 0.00092778,
                "flow_mol_phase_comp ('Liq', 'La')": 1.2023e-10,
                "flow_mol_phase_comp ('Liq', 'Dy')": 3.2886e-11,
                "flow_mol_phase_comp ('Liq', 'Ho')": 5.2653e-12,
                "flow_mol_phase_comp ('Liq', 'Er')": 1.1981e-11,
                "flow_mol_phase_comp ('Liq', 'Yb')": 7.7203e-12,
                "flow_mol_phase_comp ('Liq', 'Sm')": 4.5537e-11,
                "temperature": 298.15,
                "pressure": 1.0132e5,
            },
            "Liquid Outlet": {
                "flow_mol_phase_comp ('Liq', 'H2O')": 0.00092778,
                "flow_mol_phase_comp ('Liq', 'La')": 8.5727e-11,
                "flow_mol_phase_comp ('Liq', 'Dy')": 1.3930e-11,
                "flow_mol_phase_comp ('Liq', 'Ho')": 2.1721e-12,
                "flow_mol_phase_comp ('Liq', 'Er')": 5.0967e-12,
                "flow_mol_phase_comp ('Liq', 'Yb')": 3.2839e-12,
                "flow_mol_phase_comp ('Liq', 'Sm')": 1.9374e-11,
                "temperature": 298.15,
                "pressure": 1.0132e5,
            },
        }
    )

    pd.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)


@pytest.mark.unit
def test_get_performance_contents(m):

    ix = m.fs.unit_ix

    perf_dict = ix._get_performance_contents()

    assert perf_dict == {
        "vars": {
            "Max Breakthrough Time": ix.target_breakthrough_time,
            "EBCT": ix.ebct,
            "Number Columns": ix.number_columns,
            "Bed Volume Total": ix.bed_volume_total,
            "Bed Depth": ix.bed_depth,
            "Bed Porosity": ix.bed_porosity,
            "Service Flow Rate [BV/hr]": ix.service_flow_rate,
            "Bed Velocity": ix.loading_rate,
            "Resin Particle Diameter": ix.resin_diam,
            "Resin Bulk Density": ix.resin_density,
            "Reynolds Number": ix.N_Re,
            "Peclet Number (bed)": ix.N_Pe_bed,
            "Peclet Number (particle)": ix.N_Pe_particle,
        }
    }


@pytest.fixture(scope="module")
def m_nacl():
    """
    Uses the NaCl regenerant IX model
    """

    path = os.path.dirname(os.path.realpath(__file__))
    resin_file = os.path.join(path, "..", "data", "resin_data.json")
    comp_prop_file = os.path.join(path, "..", "data", "properties_data.json")
    parmest_file = os.path.join(path, "..", "data", "parmest_data.json")
    curve_file = os.path.join(path, "..", "data", "breakthrough_literature_data.csv")
    curve_data = pd.read_csv(curve_file)

    solver = get_solver()

    resin = "S950"
    target_component = "La"
    num_traps = 30
    c_trap_min = 1e-3
    regenerant = "NaCl"
    hazardous_waste = True

    # Add sets for solvent and ion species
    m = build_model()

    add_data(
        m,
        target_component=target_component,
        resin=resin,
        curve_data=curve_data,
        resin_file=resin_file,
        comp_prop_file=comp_prop_file,
        parmest_file=parmest_file,
    )

    parmest_data = m.fs.parmest_data
    build_clark(
        m,
        resin=resin,
        regenerant=regenerant,
        target_component=target_component,
        num_traps=num_traps,
        c_trap_min=c_trap_min,
        resin_file=resin_file,
        hazardous_waste=hazardous_waste,
    )

    set_bounds(m)

    set_operating_conditions(
        m, parmest_data=parmest_data, target_component=target_component, resin=resin
    )

    set_scaling(m)

    return m


@pytest.mark.component
def test_optimization_nacl(m_nacl):

    regenerant = "NaCl"
    target_component = "La"
    solver = get_solver()

    # Scale model
    init_scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = init_scaling.create_using(m_nacl, rename=False)

    check_dof(scaled_model, fail_flag=True)

    # Solve scaled model with zero degrees of freedom
    scaled_results = solver.solve(scaled_model)
    pyo.assert_optimal_termination(scaled_results)

    # Propagate the solution back to the original model
    init_scaling.propagate_solution(scaled_model, m_nacl)

    m = build_clark_with_costing(
        m_nacl, regenerant=regenerant, target_component=target_component
    )

    check_dof(m, fail_flag=True)

    # Solve scaled model
    results = solver.solve(m)
    pyo.assert_optimal_termination(results)

    run_optimization(m, target_component=target_component)

    results_opt = solver.solve(m)
    pyo.assert_optimal_termination(results_opt)

    # Test expected results for ix variables. NOTE: bed_expansion_frac
    # and pressure_drop are the ones calculated with resin-specific
    # data from .json file
    ix_vars_results = {
        "resin_diam": 0.00075,
        "resin_density": 1126.61,
        "bed_volume": 0.00017671456829883997,
        "bed_volume_total": 0.00017671456829883997,
        "bed_depth": 2.2499997646005814,
        "bed_porosity": 0.8,
        "column_height": 3.918249694686954,
        "bed_diameter": 0.010000000000640615,
        "number_columns": 1.0,
        "number_columns_redundant": 1.0,
        "target_breakthrough_time": 2291425.7741685435,
        "ebct": 10581.689703668251,
        "loading_rate": 0.00021263142537817904,
        "service_flow_rate": 0.34021031619855785,
        "N_Re": 0.15947356903363427,
        "N_Pe_particle": 0.020713838554897324,
        "N_Pe_bed": 62.14150916332458,
    }

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit_ix, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)

    for v, r in ix_vars_results.items():
        mv = getattr(m.fs.unit_ix, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-6) == pyo.value(mv)

    # Expected results for system costs
    sys_cost_results = {
        "utilization_factor": 1,
        "TPEC": 4.121212121212121,
        "TIC": 2.0,
        "total_capital_cost": 2777.5679769146095,
        "total_operating_cost": 3932.8064452419,
        "electricity_cost": 0.07,
        "NaCl_cost": 0.09,
        "aggregate_capital_cost": 2777.5679769146095,
        "aggregate_fixed_operating_cost": 3849.3271494414466,
        "aggregate_variable_operating_cost": 0.0,
        "aggregate_flow_electricity": 4.189304217061552e-08,
        "aggregate_flow_costs": {
            "NaCl": 0.1522307890750307,
            "electricity": 2.5706408536733103e-05,
        },
        "maintenance_labor_chemical_operating_cost": 83.32703930743828,
    }

    for v, r in sys_cost_results.items():
        mv = getattr(m.fs.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)

    # Expected results for ion exchange costs considering "single_use"
    # for regenerant use
    ix_cost_results = {
        "capital_cost": 2777.5679769146095,
        "fixed_operating_cost": 3849.3271494414466,
        "capital_cost_vessel": 598.398018289634,
        "capital_cost_resin": 1.2369525140701472,
        "capital_cost_backwash_tank": 183.3206369055396,
        "operating_cost_hazardous": 3849.2034541900457,
        "total_pumping_power": 4.189304217061552e-08,
        "capital_cost_regen_tank": 6.193409944356747,
        "flow_mass_regen_soln": 1.4243565041828834,
        "cost_factor": 2.0,
    }

    for v, r in ix_cost_results.items():
        mv = getattr(m.fs.unit_ix.costing, v)
        if mv.is_indexed():
            for i, s in r.items():
                assert pytest.approx(s, rel=1e-3) == pyo.value(mv[i])
        else:
            assert pytest.approx(r, rel=1e-3) == pyo.value(mv)


@pytest.mark.unit
def test_get_stream_table_contents_nacl(m_nacl):

    stable = m_nacl.fs.unit_ix._get_stream_table_contents()

    expected = pd.DataFrame.from_dict(
        {
            "Units": {
                "flow_mol_phase_comp ('Liq', 'H2O')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'La')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Dy')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Ho')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Er')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Yb')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "flow_mol_phase_comp ('Liq', 'Sm')": getattr(
                    pyo.units.pint_registry, "mole / second"
                ),
                "temperature": getattr(pyo.units.pint_registry, "kelvin"),
                "pressure": getattr(pyo.units.pint_registry, "pascal"),
            },
            "Feed Inlet": {
                "flow_mol_phase_comp ('Liq', 'H2O')": 0.00092778,
                "flow_mol_phase_comp ('Liq', 'La')": 1.2023e-10,
                "flow_mol_phase_comp ('Liq', 'Dy')": 3.2886e-11,
                "flow_mol_phase_comp ('Liq', 'Ho')": 5.2653e-12,
                "flow_mol_phase_comp ('Liq', 'Er')": 1.1981e-11,
                "flow_mol_phase_comp ('Liq', 'Yb')": 7.7203e-12,
                "flow_mol_phase_comp ('Liq', 'Sm')": 4.5537e-11,
                "temperature": 298.15,
                "pressure": 1.0132e5,
            },
            "Liquid Outlet": {
                "flow_mol_phase_comp ('Liq', 'H2O')": 0.00092778,
                "flow_mol_phase_comp ('Liq', 'La')": 8.5727e-11,
                "flow_mol_phase_comp ('Liq', 'Dy')": 1.3930e-11,
                "flow_mol_phase_comp ('Liq', 'Ho')": 2.1721e-12,
                "flow_mol_phase_comp ('Liq', 'Er')": 5.0967e-12,
                "flow_mol_phase_comp ('Liq', 'Yb')": 3.2839e-12,
                "flow_mol_phase_comp ('Liq', 'Sm')": 1.9374e-11,
                "temperature": 298.15,
                "pressure": 1.0132e5,
            },
            "Regen Outlet": {
                "flow_mol_phase_comp ('Liq', 'H2O')": 0.0000,
                "flow_mol_phase_comp ('Liq', 'La')": 3.4499e-11,
                "flow_mol_phase_comp ('Liq', 'Dy')": 1.8956e-11,
                "flow_mol_phase_comp ('Liq', 'Ho')": 3.0931e-12,
                "flow_mol_phase_comp ('Liq', 'Er')": 6.8847e-12,
                "flow_mol_phase_comp ('Liq', 'Yb')": 4.4365e-12,
                "flow_mol_phase_comp ('Liq', 'Sm')": 2.6163e-11,
                "temperature": 298.15,
                "pressure": 1.0132e5,
            },
        }
    )

    pd.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)


@pytest.mark.component
def test_main_in_ix_example():
    """Tests the execution of main function in example script."""

    m = main()
