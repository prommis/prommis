#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import pyomo.environ as pyo
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import unscaled_variables_generator

import pytest

from prommis.precipitate.CMR_liquid_properties import AqueousParameter
from prommis.precipitate.CMR_solid_properties import PrecipitateParameters
from prommis.precipitate.precipitator import Precipitator

def build_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_aq = AqueousParameter()
    m.fs.properties_solid = PrecipitateParameters()

    m.fs.unit = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )

    m.fs.unit.aqueous_inlet.flow_vol[0].fix(10000)

    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(11.3742)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Co"].fix(105.732)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(17454.9915)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Mg"].fix(93.6369)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Mn"].fix(40.2903)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Na"].fix(6.3279)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ni"].fix(15.9399)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Zn"].fix(668.3544)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.801)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(2.6433)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(521.0505)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(134.4078)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(1762.2)
    m.fs.unit.aqueous_inlet.conc_mass_comp[0, "OH"].fix(69.12e-2)
    # m.fs.unit.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1000)

    m.fs.unit.cv_precipitate[0].temperature.fix(298.15)

    assert_units_consistent(m.fs.unit)
    assert degrees_of_freedom(m) == 0
    initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
    initializer.initialize(m.fs.unit)
    solver = get_solver()
    results = solver.solve(m)
    assert_optimal_termination(results)

    return m
def display_performance_metrics(m):
    print("---- System Performance Metrics ----")
    print("---- Feed Metrics ----")
    f_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.flow_vol[0],
        to_units=pyo.units.m ** 3 / pyo.units.hr,
    )
    print(f"Influent flow: " f"{pyo.value(f_in):.3g}" f"{pyo.units.get_units(f_in)}")
    Ca_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ca"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ca2+ feed mass concentration: "
        f"{pyo.value(Ca_in):.3g}"
        f"{pyo.units.get_units(Ca_in)}"
    )
    Co_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Co"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Co2+ feed mass concentration: "
        f"{pyo.value(Co_in):.3g}"
        f"{pyo.units.get_units(Co_in)}"
    )
    Fe_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Fe"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Fe3+ feed mass concentration: "
        f"{pyo.value(Fe_in):.3g}"
        f"{pyo.units.get_units(Fe_in)}"
    )
    Mg_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Mg"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mg2+ feed mass concentration: "
        f"{pyo.value(Mg_in):.3g}"
        f"{pyo.units.get_units(Mg_in)}"
    )
    Mn_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Mn"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mn2+ feed mass concentration: "
        f"{pyo.value(Mn_in):.3g}"
        f"{pyo.units.get_units(Mn_in)}"
    )
    Na_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Na"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Na+ feed mass concentration: "
        f"{pyo.value(Na_in):.3g}"
        f"{pyo.units.get_units(Na_in)}"
    )
    Ni_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Ni"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ni2+ feed mass concentration: "
        f"{pyo.value(Ni_in):.3g}"
        f"{pyo.units.get_units(Ni_in)}"
    )
    Zn_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Zn"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Zn2+ feed mass concentration: "
        f"{pyo.value(Zn_in):.3g}"
        f"{pyo.units.get_units(Zn_in)}"
    )
    Dy_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Dy"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Dy3+ feed mass concentration: "
        f"{pyo.value(Dy_in):.3g}"
        f"{pyo.units.get_units(Dy_in)}"
    )
    Gd_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Gd"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Gd3+ feed mass concentration: "
        f"{pyo.value(Gd_in):.3g}"
        f"{pyo.units.get_units(Gd_in)}"
    )
    Nd_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Nd"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Nd3+ feed mass concentration: "
        f"{pyo.value(Nd_in):.3g}"
        f"{pyo.units.get_units(Nd_in)}"
    )
    Pr_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "Pr"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Pr3+ feed mass concentration: "
        f"{pyo.value(Pr_in):.3g}"
        f"{pyo.units.get_units(Pr_in)}"
    )
    SO4_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "SO4"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"SO4 2- feed mass concentration: "
        f"{pyo.value(SO4_in):.3g}"
        f"{pyo.units.get_units(SO4_in)}"
    )
    OH_in = pyo.units.convert(
        m.fs.unit.aqueous_inlet.conc_mass_comp[0, "OH"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"OH- feed mass concentration: "
        f"{pyo.value(OH_in):.3g}"
        f"{pyo.units.get_units(OH_in)}"
    )

    print("\n---- Liquid Outlet Metrics ----")
    f_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.flow_vol[0],
        to_units=pyo.units.m ** 3 / pyo.units.hr,
    )
    print(f"Effluent flow: " f"{pyo.value(f_out):.3g}" f"{pyo.units.get_units(f_out)}")
    Ca_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ca"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ca2+ feed mass concentration: "
        f"{pyo.value(Ca_out):.3g}"
        f"{pyo.units.get_units(Ca_out)}"
    )
    Co_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Co"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Co2+ effluent mass concentration: "
        f"{pyo.value(Co_out):.3g}"
        f"{pyo.units.get_units(Co_out)}"
    )
    Fe_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Fe"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Fe3+ effluent mass concentration: "
        f"{pyo.value(Fe_out):.3g}"
        f"{pyo.units.get_units(Fe_out)}"
    )
    Mg_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Mg"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mg2+ effluent mass concentration: "
        f"{pyo.value(Mg_out):.3g}"
        f"{pyo.units.get_units(Mg_out)}"
    )
    Mn_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Mn"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Mn2+ effluent mass concentration: "
        f"{pyo.value(Mn_out):.3g}"
        f"{pyo.units.get_units(Mn_out)}"
    )
    Na_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Na"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Na+ effluent mass concentration: "
        f"{pyo.value(Na_out):.3g}"
        f"{pyo.units.get_units(Na_out)}"
    )
    Ni_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Ni"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Ni2+ effluent mass concentration: "
        f"{pyo.value(Ni_out):.3g}"
        f"{pyo.units.get_units(Ni_out)}"
    )
    Zn_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Zn"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Zn2+ effluent mass concentration: "
        f"{pyo.value(Zn_out):.3g}"
        f"{pyo.units.get_units(Zn_out)}"
    )
    Dy_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Dy"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Dy3+ effluent mass concentration: "
        f"{pyo.value(Dy_out):.3g}"
        f"{pyo.units.get_units(Dy_out)}"
    )
    Gd_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Gd"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Gd3+ effluent mass concentration: "
        f"{pyo.value(Gd_out):.3g}"
        f"{pyo.units.get_units(Gd_out)}"
    )
    Nd_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Nd"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Nd3+ effluent mass concentration: "
        f"{pyo.value(Nd_out):.3g}"
        f"{pyo.units.get_units(Nd_out)}"
    )
    Pr_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "Pr"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"Pr3+ effluent mass concentration: "
        f"{pyo.value(Pr_out):.3g}"
        f"{pyo.units.get_units(Pr_out)}"
    )
    SO4_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "SO4"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"SO4 2- effluent mass concentration: "
        f"{pyo.value(SO4_out):.3g}"
        f"{pyo.units.get_units(SO4_out)}"
    )
    OH_out = pyo.units.convert(
        m.fs.unit.aqueous_outlet.conc_mass_comp[0, "OH"],
        to_units=pyo.units.mg / pyo.units.L,
    )
    print(
        f"OH- effluent mass concentration: "
        f"{pyo.value(OH_out):.3g}"
        f"{pyo.units.get_units(OH_out)}"
    )

    print("\n---- Precipitate Metrics ----")
    CoFe2O4_out = pyo.units.convert(
        m.fs.unit.precipitate_outlet.flow_mol_comp[0, "CoFe2O4(s)"],
        to_units=pyo.units.mol / pyo.units.hr,
    )
    print(
        f"CoFe2O4 mol flow rate: "
        f"{pyo.value(CoFe2O4_out):.3g}"
        f"{pyo.units.get_units(CoFe2O4_out)}"
    )
    Fe2O3_out = pyo.units.convert(
        m.fs.unit.precipitate_outlet.flow_mol_comp[0, "Fe2O3(s)"],
        to_units=pyo.units.mol / pyo.units.hr,
    )
    print(
        f"Fe2O3 mol flow rate: "
        f"{pyo.value(Fe2O3_out):.3g}"
        f"{pyo.units.get_units(Fe2O3_out)}"
    )
    Zn4_OH6_SO4_out = pyo.units.convert(
        m.fs.unit.precipitate_outlet.flow_mol_comp[0, "Zn4(OH)6SO4(s)"],
        to_units=pyo.units.mol / pyo.units.hr,
    )
    print(
        f"Zn4(OH)6SO4 mol flow rate: "
        f"{pyo.value(Zn4_OH6_SO4_out):.3g}"
        f"{pyo.units.get_units(Zn4_OH6_SO4_out)}"
    )

if __name__ == "__main__":
    # Call build model function
    m = build_flowsheet()
    display_performance_metrics(m)
