"""
Demonstration flowsheet for EvaporationPond unit model.

Authors: Andrew Lee
"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
)

from idaes.core import FlowsheetBlock

from prommis.evaporation_pond.evaporation_pond import EvaporationPond
from prommis.evaporation_pond.brine_properties import BrineParameters
from prommis.evaporation_pond.brine_precipitation_reactions import (
    BrineReactionParameters,
)


def build_model():
    """
    Method to build a single stage evaporation pond model for testing
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.brine_props = BrineParameters()
    m.fs.brine_rxns = BrineReactionParameters(property_package=m.fs.brine_props)

    m.fs.pond_1 = EvaporationPond(
        property_package=m.fs.brine_props,
        reaction_package=m.fs.brine_rxns,
    )

    return m


def set_inputs(m):
    """
    Inlet and design conditions.
    """
    # TODO: Get real input values
    m.fs.pond_1.inlet.flow_vol.fix(240 * units.L / units.s)
    m.fs.pond_1.inlet.conc_mass_comp[0, "Li"].fix(650 * units.mg / units.L)
    m.fs.pond_1.inlet.conc_mass_comp[0, "Na"].fix(650 * units.mg / units.L)
    m.fs.pond_1.inlet.conc_mass_comp[0, "Cl"].fix(1300 * units.mg / units.L)
    m.fs.pond_1.inlet.conc_mass_comp[0, "H2O"].fix(1 * units.kg / units.L)

    m.fs.pond_1.surface_area.fix(50000 * units.m**2)
    m.fs.pond_1.evaporation_rate.fix(4.75 * units.mm / units.day)


def set_scaling(m):
    """
    Apply scaling factors to improve solver performance.
    """
    pass


# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    # Call build model function
    m = build_model()
    set_inputs(m)
    set_scaling(m)

    # Solve scaled model
    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    # Display some results
    m.fs.pond_1.report()
