#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the two-salt diafiltration property model for feed streams.
"""

from pyomo.environ import ConcreteModel, Param, Var, value

from idaes.core import FlowsheetBlock

import pytest

from prommis.nanofiltration.diafiltration_solute_properties import (
    SoluteParameter,
)


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.product_properties = SoluteParameter()

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
        assert j in model.fs.product_properties.charge
        assert j in model.fs.product_properties.molar_mass
        assert j in model.fs.product_properties.sigma

    assert isinstance(model.fs.product_properties.num_solutes, Param)
    assert value(model.fs.product_properties.num_solutes) == 5


@pytest.mark.unit
def test_build(model):
    model.fs.state = model.fs.product_properties.build_state_block(model.fs.time)

    assert len(model.fs.state) == 1

    assert isinstance(model.fs.state[0].flow_vol, Var)
    assert isinstance(model.fs.state[0].conc_mass_comp, Var)

    model.fs.state[0].flow_vol.set_value(10)
    for j in model.fs.product_properties.component_list:
        model.fs.state[0].conc_mass_comp[j].set_value(1)

    model.fs.state.fix_initialization_states()

    assert model.fs.state[0].flow_vol.fixed
    for j in model.fs.product_properties.component_list:
        assert model.fs.state[0].conc_mass_comp[j].fixed
