import pytest

import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.scaling.util import get_jacobian
from idaes.core.util import DiagnosticsToolbox

from prommis.properties import (
    HClStrippingParameterBlock,
    SulfuricAcidLeachingParameters,
    TranslatorHClLeach
)
from prommis.properties.hcl_stripping_properties import HClStrippingStateBlock, HClStrippingPropertiesScaler
from prommis.properties.sulfuric_acid_leaching_properties import SulfuricAcidLeachingStateBlock, SulfuricAcidLeachingPropertiesScaler
from prommis.properties.translator_hcl_leach import TranslatorHClLeachScaler

@pytest.fixture
def model():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.hcl = HClStrippingParameterBlock()
    m.fs.h2so4 = SulfuricAcidLeachingParameters()
    m.fs.unit = TranslatorHClLeach(
        inlet_property_package=m.fs.hcl,
        outlet_property_package=m.fs.h2so4
    )

    # Taken from tear guess in UKy flowsheet
    port_values = {
        "flow_vol": {0: 16.70},
        "conc_mass_comp": {
            (0, "Al"): 2.42,
            (0, "Ca"): 0.68,
            (0, "Ce"): 0.16,
            (0, "Cl"): 1438.56,
            (0, "Dy"): 0.64,
            (0, "Fe"): 22.67,
            (0, "Gd"): 1.01,
            (0, "H"): 39.81,
            (0, "H2O"): 1000000,
            (0, "La"): 0.13,
            (0, "Nd"): 8.52e-2,
            (0, "Pr"): 2.10e-2,
            (0, "Sc"): 1.65e-3,
            (0, "Sm"): 7.88e-2,
            (0, "Y"): 1.17,
        },
    }
    for varname, var_dict in port_values.items():
        var_obj = getattr(m.fs.unit.inlet, varname)
        for idx, val in var_dict.items():
            var_obj[idx].fix(val)

    hcl_scaler = HClStrippingPropertiesScaler()
    hcl_scaler.default_scaling_factors["flow_vol"] = 1 / 16.70

    h2so4_scaler = SulfuricAcidLeachingPropertiesScaler()
    h2so4_scaler.default_scaling_factors["flow_vol"] = 1 / 16.70

    submodel_scalers = pyo.ComponentMap()
    submodel_scalers[m.fs.unit.properties_in] = hcl_scaler
    submodel_scalers[m.fs.unit.properties_out] = h2so4_scaler

    scaler_obj = m.fs.unit.default_scaler()
    scaler_obj.scale_model(m.fs.unit, submodel_scalers=submodel_scalers)

    return m
    

@pytest.mark.unit
def test_build(model):
    unit = model.fs.unit

    assert isinstance(unit.properties_in, HClStrippingStateBlock)
    assert len(unit.properties_in) == 1
    
    assert isinstance(unit.properties_out, SulfuricAcidLeachingStateBlock)
    assert len(unit.properties_out) == 1
    
    assert unit.default_scaler is TranslatorHClLeachScaler
    assert unit.default_initializer is BlockTriangularizationInitializer

    assert isinstance(unit.HCl_components, pyo.Set)
    assert len(unit.HCl_components) == 15

    assert isinstance(unit.sulfate_components, pyo.Set)
    assert len(unit.sulfate_components) == 2

    assert isinstance(unit.eps_conc_mass, pyo.Param)
    assert len(unit.eps_conc_mass) == 1
    assert pyo.value(unit.eps_conc_mass) == 1e-15
    assert_units_equivalent(unit.eps_conc_mass, pyo.units.mg / pyo.units.L)

    assert isinstance(unit.flow_vol_eqn, pyo.Constraint)
    assert len(unit.flow_vol_eqn) == 1

    assert isinstance(unit.conc_mass_comp_hcl_eqn, pyo.Constraint)
    assert len(unit.conc_mass_comp_hcl_eqn) == 15
    
    assert isinstance(unit.conc_mass_sulfates_eqn, pyo.Constraint)
    assert len(unit.conc_mass_sulfates_eqn) == 2

@pytest.mark.unit
def test_no_structural_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_structural_warnings()

@pytest.mark.unit
@pytest.mark.solver
def test_initialize_model(model):
    init_obj = model.default_initializer()
    init_obj.initialize(model)

@pytest.mark.unit
@pytest.mark.solver
def test_no_numerical_issues(model):
    dt = DiagnosticsToolbox(model)
    dt.assert_no_numerical_warnings()

    assert get_jacobian(model, scaled=False) == pytest.approx(1e12, rel=1e-3)
    assert get_jacobian(model, scaled=True) == pytest.approx(1e14, rel=1e-3)
