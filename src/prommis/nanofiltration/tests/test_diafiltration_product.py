#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Diagnostic tests for the custom product block used with the two-salt diafiltration unit model.

Author: Molly Dougher
"""

from pyomo.environ import (
    ConcreteModel,
)
from pyomo.network import Port

from idaes.core import FlowsheetBlock
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom

import pytest

from prommis.nanofiltration.diafiltration_product import DiafiltrationProduct
from prommis.nanofiltration.diafiltration_solute_product_properties import (
    SoluteProductParameter,
)


@pytest.fixture(scope="module")
def diafiltration_product():
    """
    Build a flowsheet with the custom diafiltration product unit model.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SoluteProductParameter()

    m.fs.unit = DiafiltrationProduct(property_package=m.fs.properties)

    assert degrees_of_freedom(m.fs.unit) == 0

    return m


@pytest.mark.unit
def test_config(diafiltration_product):
    # product.py has 4 config args
    assert len(diafiltration_product.fs.unit.config) == 4

    assert not diafiltration_product.fs.unit.config.dynamic
    assert not diafiltration_product.fs.unit.config.has_holdup

    assert (
        diafiltration_product.fs.unit.config.property_package
        is diafiltration_product.fs.properties
    )


class TestDiafiltrationProduct(object):
    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, diafiltration_product):
        assert isinstance(diafiltration_product.fs.unit.inlet, Port)
        assert hasattr(diafiltration_product.fs.unit.inlet, "flow_vol")
        assert hasattr(diafiltration_product.fs.unit.inlet, "conc_mass_lithium")
        assert hasattr(diafiltration_product.fs.unit.inlet, "conc_mass_cobalt")
        assert hasattr(diafiltration_product.fs.unit.inlet, "conc_mass_chlorine")

    @pytest.mark.component
    def test_diagnostics(self, diafiltration_product):
        dt = DiagnosticsToolbox(diafiltration_product.fs.unit)
        dt.assert_no_structural_warnings()
