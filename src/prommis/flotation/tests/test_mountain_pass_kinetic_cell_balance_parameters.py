####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
####################################################################################################
"""Tests for Mountain Pass kinetic cell-balance parameters."""

import json

import pytest

from prommis.flotation import (
    mountain_pass_flotation_fixed_recovery_flowsheet as recovery_flowsheet,
)
from prommis.flotation import (
    mountain_pass_flotation_kinetic_cell_balance_flowsheet as cell_balance,
)
from prommis.flotation.bastnaesite_properties import COMPONENTS
from prommis.flotation.initializer import (
    audit_cascade_solvability,
    cell_cascade_forward,
)


@pytest.mark.unit
def test_json_file_loads_from_source_tree():
    assert cell_balance.KINETIC_CELL_BALANCE_PARAMETERS_FILE.is_file()
    with cell_balance.KINETIC_CELL_BALANCE_PARAMETERS_FILE.open(
        encoding="utf-8"
    ) as file:
        data = json.load(file)

    assert data == cell_balance.load_kinetic_cell_balance_parameters()


@pytest.mark.unit
def test_json_schema():
    data = cell_balance.load_kinetic_cell_balance_parameters()

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
    assert set(data["property_package"]["rho_mass_comp_kg_per_m3"]) == set(COMPONENTS)
    assert set(data["banks"]) == set(recovery_flowsheet.BANK_NAMES)
    for bank_data in data["banks"].values():
        assert {
            "number_of_cells",
            "cell_volume",
            "cell_total_solid_holdup",
            "F_in_reference_kg_per_h",
            "k_cb_per_hour",
            "R_ref",
            "calibration",
        } <= set(bank_data)
        assert set(bank_data["F_in_reference_kg_per_h"]) == set(COMPONENTS)
        assert set(bank_data["k_cb_per_hour"]) == set(COMPONENTS)
        assert set(bank_data["R_ref"]) == set(COMPONENTS)


@pytest.mark.unit
def test_reference_recovery_round_trip():
    data = cell_balance.load_kinetic_cell_balance_parameters()

    for bank_name, bank_data in data["banks"].items():
        M_total_by_cell = [
            bank_data["cell_total_solid_holdup"]
            for _ in range(bank_data["number_of_cells"])
        ]
        _, _, F_float = cell_cascade_forward(
            bank_data["F_in_reference_kg_per_h"],
            M_total_by_cell,
            bank_data["k_cb_per_hour"],
            context=f"bank={bank_name}",
        )
        for component in COMPONENTS:
            recovery = (
                sum(cell[component] for cell in F_float)
                / bank_data["F_in_reference_kg_per_h"][component]
            )
            assert recovery == pytest.approx(bank_data["R_ref"][component], abs=1e-9)


@pytest.mark.unit
def test_cascade_audit_walks_complete():
    data = cell_balance.load_kinetic_cell_balance_parameters()

    for bank_name, bank_data in data["banks"].items():
        M_total_by_cell = [
            bank_data["cell_total_solid_holdup"]
            for _ in range(bank_data["number_of_cells"])
        ]
        audit = audit_cascade_solvability(
            bank_data["F_in_reference_kg_per_h"],
            M_total_by_cell,
            bank_data["k_cb_per_hour"],
            bank_name=bank_name,
        )
        assert audit["walk_completed"]
        assert audit["first_failure_mode"] is None


@pytest.mark.unit
def test_legacy_generator_keys_are_not_committed():
    data = cell_balance.load_kinetic_cell_balance_parameters()

    legacy_keys = {
        "R_inf",
        "warm_start",
        "feasibility_margin_min",
        "calibration_bound_proximity",
    }
    for bank_data in data["banks"].values():
        assert legacy_keys.isdisjoint(bank_data)
