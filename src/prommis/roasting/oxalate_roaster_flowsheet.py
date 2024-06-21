#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

"""
flowsheet to test the roster unit model
authors: J. Ma
"""
# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES standard unit model
import idaes.logger as idaeslog

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    EosType,
    get_prop,
)

from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.roasting.ree_oxalate_roaster import REEOxalateRoaster

_log = idaeslog.getModelLogger(__name__)


def main(m=None):
    """Create concrete model, make the flowsheet object, fix some variables,
    and solve the problem."""

    if m is None:
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        # Add property packages to flowsheet library
        gas_species = {"O2", "H2O", "CO2", "N2"}
        m.fs.prop_gas = GenericParameterBlock(
            **get_prop(gas_species, ["Vap"], EosType.IDEAL),
            doc="gas property",
        )
        key_components = {
            "H^+",
            "Ce^3+",
            "Al^3+",
            "Fe^3+",
            # "Fe^2+",
            # "Ca^2+",
            # "Mg^2+",
            "C2O4^2-",
            # "NO3^-",
            # "SO4^2-",
            # "Cl^-",
        }
        m.fs.prop_solid = PrecipitateParameters(
            key_components=key_components,
        )

    create_model(m)
    set_inputs(m)
    initialize_system(m)
    solver = get_solver(options={"max_iter": 50})
    dof = degrees_of_freedom(m)
    print("dof=", dof)
    result = solver.solve(m, tee=True)
    print("Gas feed mole flow =", pyo.value(m.fs.roaster.gas_in[0].flow_mol), "mol/s")
    print(
        "Solid feed Ce mole flow =",
        m.fs.roaster.solid_in[0].flow_mol_comp["Ce2(C2O4)3(s)"].value,
        "mol/s",
    )
    print("heat_duty=", m.fs.roaster.heat_duty[0].value)
    print("mass fraction of metal oxide in solid product:")
    for x in m.fs.roaster.metal_list:
        print(x, pyo.value(m.fs.roaster.mass_frac_comp_product[0, x]))
    return m


def create_model(m):
    """Create unit models"""
    m.fs.roaster = REEOxalateRoaster(
        property_package_gas=m.fs.prop_gas,
        property_package_precipitate=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        # metal_list = ["Sc","Y","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Tm","Yb","Lu"], default is ["Ce"] only
        metal_list=["Ce"],
    )


def set_inputs(m):
    """fix variables for geometry and design data"""
    m.fs.roaster.deltaP.fix(0)
    m.fs.roaster.gas_inlet.temperature.fix(1330)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    # inlet flue gas mole flow rate
    fgas = 0.00781
    # inlet flue gas composition, typical flue gas by burning CH4 with air with stoichiometric ratio 0f 2.3
    gas_comp = {
        "O2": 0.1118,
        "H2O": 0.1005,
        "CO2": 0.0431,
        "N2": 0.7446,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    m.fs.roaster.gas_inlet.flow_mol.fix(fgas)

    # fix outlet product temperature
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature
    m.fs.roaster.solid_in[0].temperature.fix(298.15)
    m.fs.roaster.solid_in[0].flow_mol_comp["Ce2(C2O4)3(s)"].fix(6.1e-5)
    m.fs.roaster.flow_mol_moist_feed.fix(6.75e-4)
    # total solid mass flow rate including surface moisture

    """
    m.fs.roaster.mass_frac_feed_dry[0,'Sc'].fix(0.001648997)
    m.fs.roaster.mass_frac_feed_dry[0,'Y'].fix(0.0619823)
    m.fs.roaster.mass_frac_feed_dry[0,'La'].fix(0.160501197)
    m.fs.roaster.mass_frac_feed_dry[0,'Ce'].fix(0.415971651)
    m.fs.roaster.mass_frac_feed_dry[0,'Pr'].fix(0.072397859)
    m.fs.roaster.mass_frac_feed_dry[0,'Nd'].fix(0.200619089)
    m.fs.roaster.mass_frac_feed_dry[0,'Sm'].fix(0.036538614)
    m.fs.roaster.mass_frac_feed_dry[0,'Eu'].fix(0.002551071)
    m.fs.roaster.mass_frac_feed_dry[0,'Gd'].fix(0.01729078)
    m.fs.roaster.mass_frac_feed_dry[0,'Tb'].fix(0.00745766)
    m.fs.roaster.mass_frac_feed_dry[0,'Dy'].fix(0.01002913)
    m.fs.roaster.mass_frac_feed_dry[0,'Tm'].fix(0.004159939)
    m.fs.roaster.mass_frac_feed_dry[0,'Yb'].fix(0.005879861)
    m.fs.roaster.mass_frac_feed_dry[0,'Lu'].fix(0.002971853)
    """
    m.fs.roaster.frac_comp_recovery.fix(0.95)


def initialize_system(m):
    initializer = BlockTriangularizationInitializer()
    initializer.initialize(m.fs.roaster)
    assert initializer.summary[m.fs.roaster]["status"] == InitializationStatus.Ok


if __name__ == "__main__":
    """
    Main function to to run simulation
    To run steady-state model, call main_steady()
    to run dynamic model, call main_dyn()
    """
    m = main()
