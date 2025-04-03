#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Tests the diafiltration flowsheet
"""

from pyomo.environ import value

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.models.unit_models import MSContactor

import pytest

from prommis.nanofiltration.diafiltration import (
    add_costing,
    add_product_constraints,
    build_model,
    main,
)


@pytest.mark.component
def test_main():
    """
    Tests the execution of the main function in diafiltration.py
    """
    m = main()
    dt = DiagnosticsToolbox(m)
    dt.assert_no_numerical_warnings()

    assert isinstance(
        m.fs.stage3, MSContactor
    )  # check that stage3 exists and is an MSContactor
    # Retentate side checks
    assert hasattr(
        m.fs.stage3, "retentate_inlet_state"
    )  # check that there is a retentate feed
    assert hasattr(
        m.fs.stage3, "retentate_side_stream_state"
    )  # check that a side stream exists
    for k in m.fs.stage3.retentate_side_stream_state:
        assert k == (0, 10)  # check that the side stream only exists at element 10
    assert not hasattr(
        m.fs.stage3, "retentate_energy_balance"
    )  # check that there are no energy balances
    assert not hasattr(
        m.fs.stage3, "retentate_pressure_balance"
    )  # check that there are no pressure balances
    # Permeate side checks
    assert not hasattr(
        m.fs.stage3, "permeate_inlet_state"
    )  # check that there is no permeate feed
    assert not hasattr(
        m.fs.stage3, "permeate_side_stream_state"
    )  # check that there are no side streams on permeate side
    assert not hasattr(
        m.fs.stage3, "permeate_energy_balance"
    )  # check that there are no energy balances
    assert not hasattr(
        m.fs.stage3, "permeate_pressure_balance"
    )  # check that there are no pressure balances

    test_dict = {
        "lithium_recovery": [value(m.Li_recovery) * 100, 94.99999900000033],
        "lithium_purity": [value(m.Li_purity) * 100, 20.8820384057719],
        "cobalt_recovery": [value(m.Co_recovery) * 100, 63.49999900000085],
        "cobalt_purity": [value(m.Co_purity) * 100, 99.20784627170993],
        "membrane_area": [
            (value(m.membrane_length) * value(m.w)),
            3502.558324664463,
        ],
        "total_annualized_cost": [
            value(m.fs.costing.total_annualized_cost),
            251240.71142955605,
        ],
    }

    for model_result, test_val in test_dict.values():
        assert pytest.approx(test_val, rel=1e-5) == value(model_result)


@pytest.mark.component
def test_Li_purity_constraint_exception():
    m = build_model()
    add_costing(m)

    with pytest.raises(
        ValueError, match="A lithium product purity bound was not provided"
    ):
        add_product_constraints(
            m,
            Li_recovery_bound=0.95,
            Co_recovery_bound=0.635,
            recovery=True,
            Co_purity_bound=0.5,
            purity=True,
        )


@pytest.mark.component
def test_Co_purity_constraint_exception():
    m = build_model()
    add_costing(m)

    with pytest.raises(
        ValueError, match="A cobalt product purity bound was not provided"
    ):
        add_product_constraints(
            m,
            Li_recovery_bound=0.95,
            Co_recovery_bound=0.635,
            recovery=True,
            Li_purity_bound=0.5,
            purity=True,
        )
