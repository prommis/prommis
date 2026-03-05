import pytest

from pyomo.environ import ConcreteModel, units

from idaes.core import (
    EnergyBalanceType,
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json


from idaes.core.solvers import get_solver

from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction import SettlerTank
from prommis.solvent_extraction.settler_tank import SettlerTankInitializer, SettlerTankScaler
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.organic_parameters = REESolExOgParameters()
    m.fs.aqueous_parameters = SulfuricAcidLeachingParameters()

    m.fs.unit = SettlerTank(
        light_phase_alias="organic",
        heavy_phase_alias="aqueous",
        light_phase_config={
            "property_package": m.fs.leach_soln,
            "energy_balance_type": EnergyBalanceType.isothermal,
            "has_pressure_balance": False,
        },
        heavy_phase_config={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.backward,
            "energy_balance_type": EnergyBalanceType.isothermal,
            "has_pressure_balance": False,
        },
        has_holdup=True,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=4,
    )

    # Specify aqueous inlet
    m.fs.unit.aqueous_inlet.flow_vol.fix(50)
    m.fs.unit.aqueous_inlet.temperature.fix(300)
    m.fs.unit.aqueous_inlet.pressure.fix(1e5)

    # Fix most concentrations to effectively zero
    # and then change a handfulto be nonzero
    m.fs.unit.aqueous_inlet.conc_mass_comp.fix(1e-20)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H"].fix(10.75)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(100)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e4)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Y"].fix(8.89)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(1752.34)
    
    # Specify organic inlet
    m.fs.unit.organic_inlet.flow_vol.fix(50)
    m.fs.unit.organic_inlet.temperature.fix(300)
    m.fs.unit.organic_inlet.pressure.fix(1e5)

    # Fix most concentrations to zero and then change a handful
    m.fs.unit.organic_inlet.conc_mass_comp.fix(1e-20)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(0.05 * 975.8e3)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(20)
    m.fs.unit.organic_inlet.conc_mass_comp[0, "Y_o"].fix(10)

    m.fs.unit.length.fix(1)
    m.fs.unit.settler_width(0.5)
    m.fs.unit.light_phase_weir_height(0.5)

@pytest.mark.unit
def test_model_construction(model):
    pass