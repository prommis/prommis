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

import idaes.core.util.scaling as iscale

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

from prommis.leaching.leach_solids_properties import CoalRefuseParameters
from prommis.roasting.ree_feed_roaster import REEFeedRoaster

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
        m.fs.prop_solid = CoalRefuseParameters(
            doc="solid property",
        )

    create_model(m)
    set_inputs(m)
    iscale.calculate_scaling_factors(m)
    initialize_system(m)
    solver = get_solver(options={"max_iter": 50})
    dof = degrees_of_freedom(m)
    print("dof=", dof)
    result = solver.solve(m, tee=True)
    print("heat_duty=", m.fs.roaster.heat_duty[0].value)
    print(
        "Total mass flow of dry organic free impurity in feed stream=",
        pyo.value(m.fs.roaster.flow_mass_impurity_feed[0]),
    )
    print(
        "Un mass fraction in dry organic free feed",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_feed[0, "Un"]),
    )
    print(
        "Al mass fraction in dry organic free feed",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_feed[0, "Al"]),
    )
    print(
        "Ca mass fraction in dry organic free feed",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_feed[0, "Ca"]),
    )
    print(
        "Fe mass fraction in dry organic free feed",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_feed[0, "Fe"]),
    )
    print(
        "Si mass fraction in dry organic free feed",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_feed[0, "Si"]),
    )
    print(
        "Total mass flow of recovered solid product=",
        pyo.value(m.fs.roaster.flow_mass_product_recovered[0]),
    )
    print(
        "Un mass fraction in recovered solid product",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_product[0, "Un"]),
    )
    print(
        "Al mass fraction in recovered solid product",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_product[0, "Al"]),
    )
    print(
        "Ca mass fraction in recovered solid product",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_product[0, "Ca"]),
    )
    print(
        "Fe mass fraction in recovered solid product",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_product[0, "Fe"]),
    )
    print(
        "Si mass fraction in recovered solid product",
        pyo.value(m.fs.roaster.mass_frac_comp_impurity_ele_product[0, "Si"]),
    )
    print("insoluble REE in the solid recovered product in ppm of total mass:")
    for i in m.fs.roaster.ree_list:
        print(i, pyo.value(m.fs.roaster.ppm_comp_ree_ins_product[0, i]))
    print("dissovable REE in the solid recovered product in ppm of total mass:")
    for i in m.fs.roaster.ree_list:
        print(i, pyo.value(m.fs.roaster.ppm_comp_ree_dis_product[0, i]))
    print(
        "mass flow rate of solid outlet=",
        pyo.value(m.fs.roaster.solid_out[0].flow_mass),
    )
    print("mass fractions of species in solid outlet:")
    for i in m.fs.prop_solid.component_list:
        print(i, pyo.value(m.fs.roaster.solid_out[0].mass_frac_comp[i]))
    return m


def create_model(m):
    """Create unit models"""
    m.fs.roaster = REEFeedRoaster(
        gas_property_package=m.fs.prop_gas,
        solid_property_package=m.fs.prop_solid,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
        ree_list=[
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
        ],
    )


def set_inputs(m):
    """fix variables for geometry and design data"""
    # assuming no pressure drop
    m.fs.roaster.deltaP.fix(0)
    # inlet flue gas composition, typical flue gas by burning CH4 with air to get excess O2 at 15.5%
    m.fs.roaster.gas_inlet.temperature.fix(895)
    m.fs.roaster.gas_inlet.pressure.fix(101325)
    gas_comp = {
        "O2": 0.155,
        "H2O": 0.0617,
        "CO2": 0.0235,
        "N2": 0.7598,
    }
    for i, v in gas_comp.items():
        m.fs.roaster.gas_inlet.mole_frac_comp[0, i].fix(v)
    # inlet flue gas mole flow rate
    m.fs.roaster.gas_inlet.flow_mol.fix(80)

    # fix outlet temperature at 600 C
    m.fs.roaster.gas_outlet.temperature.fix(873.15)

    # solid feed temperature
    m.fs.roaster.temp_feed.fix(298.15)
    # limestone conversion
    m.fs.roaster.xconv_caco3.fix(0.5)

    # solid feed rate at 1030 kg/hr
    m.fs.roaster.flow_mass_feed.fix(1030 / 3600)
    # surface moisture in solid feed
    m.fs.roaster.mass_frac_moist_feed.fix(0.029126)
    # combustible organic material
    # Note that the heat of combustion will cause very large heat release
    m.fs.roaster.mass_frac_organic_feed.fix(0.294854)

    # impurity minerals modeled plus some small amount of unknown mineral that is modeled as Un2O3 here
    # note that the unknown in UK's excel sheet also include O, C, and S elements in minerals
    # assume all Al2O3 is in Kaolinite form
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Un2O3"].fix(0.00793)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "CaCO3"].fix(0.041217)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "SiO2"].fix(0.143882)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Al2O3"].fix(0.0)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Kaolinite"].fix(0.611613)
    m.fs.roaster.mass_frac_comp_impurity_feed[0, "Pyrite"].fix(0.195359)

    # organic composition assumed as a typical coal
    m.fs.roaster.mass_frac_comp_organic_feed[0, "C"].fix(0.804)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "H"].fix(0.055)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "O"].fix(0.111)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "N"].fix(0.018)
    m.fs.roaster.mass_frac_comp_organic_feed[0, "S"].fix(0.012)

    # Insoluble REE based on UK's excel file
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Sc"].fix(10.5836902)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Y"].fix(18.25358553)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "La"].fix(31.88774016)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Ce"].fix(67.77169805)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Pr"].fix(12.74306914)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Nd"].fix(30.06183493)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Sm"].fix(6.93187974)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Eu"].fix(0.902019051)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Gd"].fix(5.559717425)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Tb"].fix(0.492010392)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Dy"].fix(3.394871702)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Ho"].fix(0.984020783)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Er"].fix(3.044997646)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Tm"].fix(0.69428133)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Yb"].fix(2.20857998)
    m.fs.roaster.ppm_comp_ree_ins_feed[0, "Lu"].fix(0.579478906)

    # Dissovable REE based on UK's excel file
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Sc"].fix(2.896677798)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Y"].fix(4.995871471)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "La"].fix(8.727438841)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Ce"].fix(18.54861295)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Pr"].fix(3.487683857)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Nd"].fix(8.227702072)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Sm"].fix(1.89720426)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Eu"].fix(0.246875949)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Gd"].fix(1.521653575)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Tb"].fix(0.134659608)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Dy"].fix(0.929151298)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Ho"].fix(0.269319217)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Er"].fix(0.833393354)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Tm"].fix(0.19001967)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Yb"].fix(0.60447202)
    m.fs.roaster.ppm_comp_ree_dis_feed[0, "Lu"].fix(0.158599094)

    # conversion of insoluble REE to dissovable REE
    m.fs.roaster.xconv_comp_ins[0, "Sc"].fix(0.850340479)
    m.fs.roaster.xconv_comp_ins[0, "Y"].fix(0.737455969)
    m.fs.roaster.xconv_comp_ins[0, "La"].fix(0.784164924)
    m.fs.roaster.xconv_comp_ins[0, "Ce"].fix(0.706035503)
    m.fs.roaster.xconv_comp_ins[0, "Pr"].fix(0.478582677)
    m.fs.roaster.xconv_comp_ins[0, "Nd"].fix(0.674653216)
    m.fs.roaster.xconv_comp_ins[0, "Sm"].fix(0.678933333)
    m.fs.roaster.xconv_comp_ins[0, "Eu"].fix(0.1)
    m.fs.roaster.xconv_comp_ins[0, "Gd"].fix(0.023783784)
    m.fs.roaster.xconv_comp_ins[0, "Tb"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Dy"].fix(0.712053571)
    m.fs.roaster.xconv_comp_ins[0, "Ho"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Er"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Tm"].fix(1)
    m.fs.roaster.xconv_comp_ins[0, "Yb"].fix(0.766633166)
    m.fs.roaster.xconv_comp_ins[0, "Lu"].fix(1)

    # recovery fraction of all REE
    m.fs.roaster.frac_comp_ree_recovery.fix(0.99)

    # recovery fraction of impurity minerals
    m.fs.roaster.frac_impurity_recovery.fix(0.99)


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
