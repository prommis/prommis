####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Tests for Mountain Pass kinetic flotation parameters."""

import json

import pytest

from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation import (
    mountain_pass_flotation_kinetic_closed_form_flowsheet as kinetic,
)
from prommis.flotation.bastnaesite_properties import COMPONENTS


def _forward_recovery(k_cf_per_min, tau_min, number_of_cells):
    return 1 - (1 + k_cf_per_min * tau_min / number_of_cells) ** (-number_of_cells)


@pytest.mark.unit
def test_json_file_loads_from_source_tree():
    assert kinetic.KINETIC_PARAMETERS_FILE.is_file()
    with kinetic.KINETIC_PARAMETERS_FILE.open(encoding="utf-8") as file:
        data = json.load(file)

    assert data == kinetic.load_kinetic_parameters()


@pytest.mark.unit
def test_json_schema():
    data = kinetic.load_kinetic_parameters()

    assert set(data) == {
        "schema_version",
        "source",
        "property_package",
        "bank_constants",
        "basis",
        "assumptions",
        "units",
        "banks",
    }
    assert data["schema_version"] == 1
    assert "rho_mass_comp_kg_per_m3" in data["property_package"]
    assert "rho_water_kg_per_m3" not in data["property_package"]
    assert "rho_water_kg_per_m3" in data["bank_constants"]
    assert set(data["property_package"]["rho_mass_comp_kg_per_m3"]) == set(COMPONENTS)
    assert set(data["banks"]) == set(recovery_flowsheet.BANK_NAMES)
    for bank_data in data["banks"].values():
        assert {
            "number_of_cells",
            "cell_volume",
            "cell_volume_basis",
            "air_holdup",
            "pulp_solids_mass_fraction",
            "tau_required_min",
            "Q_slurry_m3_per_min",
            "k_cf_per_min",
            "R_ref",
        } <= set(bank_data)
        assert set(bank_data["k_cf_per_min"]) == set(COMPONENTS)
        assert set(bank_data["R_ref"]) == set(COMPONENTS)


@pytest.mark.unit
def test_basis_values_are_consistent():
    basis = kinetic.load_kinetic_parameters()["basis"]

    assert (
        basis["json_fresh_feed_total_kg_per_h"] * 24 / basis["shift_duration_h"]
    ) == pytest.approx(basis["plant_throughput_t_per_day"], rel=0.02)
    assert (
        basis["scale_jsonbasis_to_plant_h"]
        * recovery_flowsheet.stream_total(recovery_flowsheet.FEED_INPUTS["fresh_feed"])
    ) == pytest.approx(basis["plant_dry_feed_kg_per_h"], abs=1)


@pytest.mark.unit
def test_cell_volume_back_calculation():
    data = kinetic.load_kinetic_parameters()

    for bank_name, bank_data in data["banks"].items():
        effective_volume = (
            bank_data["tau_required_min"] / 60 * bank_data["Q_slurry_m3_per_min"] * 60
        )
        expected_cell_volume = effective_volume / (
            bank_data["number_of_cells"] * (1 - bank_data["air_holdup"])
        )
        assert bank_data["cell_volume"] == pytest.approx(expected_cell_volume, abs=1e-3)
        if bank_name == "rougher":
            assert bank_data["cell_volume_basis"] == "paper_anchored"
        else:
            assert bank_data["cell_volume_basis"] == "inferred_from_assumed_tau"


@pytest.mark.unit
def test_reference_recoveries_match_closed_form_equation():
    data = kinetic.load_kinetic_parameters()

    for bank_data in data["banks"].values():
        for component in COMPONENTS:
            assert _forward_recovery(
                bank_data["k_cf_per_min"][component],
                bank_data["tau_required_min"],
                bank_data["number_of_cells"],
            ) == pytest.approx(bank_data["R_ref"][component], abs=1e-3)


@pytest.mark.unit
def test_tau_rougher_in_expected_range():
    data = kinetic.load_kinetic_parameters()

    assert 2 <= data["banks"]["rougher"]["tau_required_min"] <= 15
