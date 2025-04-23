#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the two-salt diafiltration property model for feed streams.
"""

from pyomo.environ import ConcreteModel, Set, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.nanofiltration.diafiltration_solute_product_properties import (
    SoluteProductParameter,
)


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.product_properties = SoluteProductParameter()

    return m


@pytest.mark.unit
def test_parameters(model):
    assert len(model.fs.product_properties.phase_list) == 1
    for k in model.fs.product_properties.phase_list:
        assert k == "liquid"

    for j in model.fs.product_properties.component_list:
        assert j in [
            "Li",
            "Co",
            "Cl",
        ]

    assert model.fs.product_properties.nfe == 5
    assert isinstance(model.fs.product_properties.x, Set)


@pytest.mark.unit
def test_build(model):
    model.fs.state = model.fs.product_properties.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_vol, Var)
    assert isinstance(model.fs.state[0].conc_mass_comp, Var)

    for x in model.fs.product_properties.x:
        model.fs.state[0].flow_vol[x].set_value(10)
        for j in model.fs.product_properties.component_list:
            model.fs.state[0].conc_mass_comp[j, x].set_value(1)

    model.fs.state.fix_initialization_states()

    for x in model.fs.product_properties.x:
        assert model.fs.state[0].flow_vol[x].fixed
        for j in model.fs.product_properties.component_list:
            assert model.fs.state[0].conc_mass_comp[j, x].fixed
