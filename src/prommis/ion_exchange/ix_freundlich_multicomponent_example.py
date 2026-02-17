#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This model is derived from
# https://github.com/watertap-org/watertap/blob/main/watertap/flowsheets/ion_exchange/ion_exchange_demo.py

# WaterTAP License Agreement

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

r"""This model is an example on how to use the ion exchange multicomponent
model (IXMC) for the removal of REEs.

"""

# Import Python libraries
import os
import json
import logging
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

# Import Pyomo components
import pyomo.environ as pyo

# Import IDAES models and libraries
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.scaling.scaling_base import ScalerBase
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.exceptions import ConfigurationError

# Import WaterTAP models and libraries
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

# Import PrOMMiS ion exchange models
from prommis.ion_exchange.ion_exchange_multicomponent import (
    IonExchangeMultiComp,
)
from prommis.ion_exchange.costing.ion_exchange_cost_model import (
    IXCosting,
    IXCostingData,
)


logging.basicConfig(level=logging.INFO)
logging.getLogger("pyomo.repn.plugins.nl_writer").setLevel(logging.ERROR)


""" This model solves an example for the separation of multiple rare
earth elements (REEs) using an Ion Exchange (IX) model


REFERENCES: 

[1] M. Hermassi, M. Granados, C. Valderrama, N. Skoglund, C. Ayora,
J.L. Cortina, Impact of functional group types in ion exchange resins
on rare earth element recovery from treated acid mine waters, Journal
of Cleaner Production, Volume 379, Part 2, 2022, 134742,
https://doi.org/10.1016/j.jclepro.2022.134742.

[2] WebPlotDigitizer, https://automeris.io/docs/

[3] Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle,
R., Oppenheimer, J., & Jacangelo, J. G. (2023).  Fundamental kinetic
constants for breakthrough of per- and polyfluoroalkyl substances at
varying empty bed contact times: Theoretical analysis and pilot scale
demonstration.  Chemical Engineering Journal,
464. doi:10.1016/j.cej.2023.142587

[4] Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012).
Chapter 16: Ion Exchange. MWH's Water Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

[5] Ana T. Lima and Lisbeth M. Ottosen 2022 J. Electrochem. Soc. 169 033501

"""

__author__ = "Soraya_Rawlings"


def main():
    """Method to build an ion exchange unit model using data from the
    literature.

    """

    # Define path for the data files for resins, properties of
    # components, and calculated parameters from Parmest model
    # solution. The pressure drop and bed expansion parameters in the
    # resin_data file where obtained using ref[2].
    # curr_dir = dirname(abspath(__file__))
    # data_path = abspath(join(curr_dir, "data"))
    # print(data_ptah)

    path = os.path.dirname(os.path.realpath(__file__))
    resin_file = os.path.join(path, "data", "resin_data.json")
    comp_prop_file = os.path.join(path, "data", "properties_data.json")
    parmest_file = os.path.join(path, "data", "parmest_data.json")

    # Read original breakthrough data for multiple components. The
    # data is from ref[1], Figure 4, using ref[2].
    curve_file = os.path.join(path, "data", "breakthrough_literature_data.csv")
    curve_data = pd.read_csv(curve_file)

    # Define solver to use from WaterTAP
    solver = get_solver()

    # Add relevant data for IX model: resin name, target component,
    # number of trapezoidal points, and the minimum concentration for
    # the period 0 in the trapezoidal rule equations.
    resin = "S950"
    target_component = "La"
    num_traps = 30
    c_trap_min = 1e-3
    regenerant = "single_use"
    save_plot = False
    hazardous_waste = False

    # Build the model by calling the Pyomo ConcreteModel() component
    m = build_model()

    # Add all relevant data
    add_data(
        m,
        resin=resin,
        curve_data=curve_data,
        resin_file=resin_file,
        comp_prop_file=comp_prop_file,
        parmest_file=parmest_file,
    )

    # Build IX model
    parmest_data = m.fs.parmest_data
    build_clark(
        m,
        resin=resin,
        regenerant=regenerant,
        target_component=target_component,
        num_traps=num_traps,
        c_trap_min=c_trap_min,
        resin_file=resin_file,
        hazardous_waste=hazardous_waste,
    )

    # Add lower and upper bounds to relevant variables in the IX model
    set_bounds(m)

    # Add the operating conditions for the column
    set_operating_conditions(
        m, parmest_data=parmest_data, target_component=target_component, resin=resin
    )

    # Call IDAES diagnostics tool box to make sure there are no
    # structural issues in the model
    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    print()
    print("=========== Start initialization")
    initialize_system(m, solver=solver)

    # Add scaling factors to the IX model variables and constraints
    set_scaling(m)

    # Scale model
    scaling = pyo.TransformationFactory("core.scale_model")
    scaled_model = scaling.create_using(m, rename=False)

    # Solve scaled initialization model
    init_scaled_results = model_solve(scaled_model, solver=solver)
    pyo.assert_optimal_termination(init_scaled_results)

    # Propagate the solution back to the original model
    init_results = scaling.propagate_solution(scaled_model, m)

    print("=========== End initialization")

    # Add relevant costing metrics for the IX model
    add_costing(m)

    print()
    print("=========== Re-run initialization with costing")

    # Solve and re-initialize the model after adding costing variables
    initialize_system(m, solver=solver)

    dt = DiagnosticsToolbox(model=m)
    dt.assert_no_structural_warnings()

    init_costing_results = model_solve(m, solver=solver)
    pyo.assert_optimal_termination(init_costing_results)

    print("=========== End initialization(s)")
    print()
    print()
    print("=========== Start optimization run")
    print()
    # Solve an optimization model
    run_optimization(m, target_component=target_component)

    results = model_solve(m)
    pyo.assert_optimal_termination(results)

    print()
    print()
    print("****** Optimization results")
    m.fs.unit_ix.report()

    # Plot the breakthrough profiles with optimization solution
    output_plot = "breakthrough_curves_optimization"
    plot_traps(
        m,
        results=results,
        curve_data=curve_data,
        output_plot=output_plot,
        save_plot=save_plot,
        show_plot=False,
    )


def add_data(
    m,
    resin=None,
    curve_data=None,
    resin_file=None,
    comp_prop_file=None,
    parmest_file=None,
):
    """
    Add component properties, resin, and parameter estimation data

    """

    # Declare all components in the feed stream as lists
    m.list_solvent = ["H2O"]
    m.list_reactive_ions = ["La", "Dy", "Ho", "Er", "Yb", "Sm"]
    m.list_all = m.list_solvent + m.list_reactive_ions

    # Declare lists and dictionaries to be used during scaling
    m.fs.calc_from_constr_dict = {
        "service_flow_rate": "eq_service_flow_rate",
        "bed_volume_total": "eq_bed_design",
        "column_height": "eq_column_height",
    }
    m.fs.scale_from_value = [
        "bed_volume_total",
        "bed_diameter",
        "column_height",
        "service_flow_rate",
        "ebct",
        "N_Re",
        "N_Sh",
        "N_Pe_bed",
        "N_Pe_particle",
    ]

    # Read data for resins, properties of components, and calculated
    # parameters from Parmest model solution using .json files.
    with open(resin_file, "r") as file:
        m.fs.resin_data = json.load(file)
    with open(comp_prop_file, "r") as file:
        m.fs.props_data = json.load(file)
    with open(parmest_file, "r") as file:
        m.fs.parmest_data = json.load(file)

    # Save the data of breakthrough curves in a new empty dictionary
    # "data_init" and assign the relevant values to new elements in
    # the dictionary.
    m.fs.data_init = {}
    m.fs.data_init["conc_mass"] = curve_data.set_index("compound")[
        "c0"
    ].to_dict()  # in mg/L
    m.fs.data_init["c_norm"] = curve_data.set_index("compound")[
        "c_norm"
    ].to_dict()  # dimensionless
    m.fs.data_init["flow_in"] = curve_data["flow_in"][0]  # in m3/s
    m.fs.data_init["bed_diam"] = curve_data["bed_diam"][0]  # in m
    m.fs.data_init["bed_depth"] = curve_data["bed_depth"][0]  # in m
    m.fs.data_init["loading_rate"] = curve_data["vel_bed"][0]  # in m/s
    m.fs.data_init["density_solution"] = 1e3  # in kg/m3, assumed

    return m


def build_model():
    """
    Method to build a flowsheet
    """

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    return m


def build_clark(
    m,
    resin=None,
    regenerant=None,
    target_component=None,
    num_traps=None,
    c_trap_min=None,
    resin_file=None,
    hazardous_waste=None,
):

    # Add sets for solvent and ion species
    m.fs.set_solvent = pyo.Set(initialize=m.list_solvent)
    m.fs.set_reactive_ions = pyo.Set(initialize=m.list_reactive_ions)
    m.fs.set_all = pyo.Set(initialize=m.list_all)

    # Declare properties for ions needed in the IX unit model
    ion_props = {
        "solute_list": [],
        "diffusivity_data": {},
        "molar_volume_data": {},
        "mw_data": {"H2O": m.fs.props_data["mw"]["H2O"]},  # in kg/mol
        "charge": {},
    }
    for ion in m.fs.set_reactive_ions:
        ion_props["solute_list"].append(ion)
        # Diffusivity data from ref[5]
        if ion in m.fs.props_data["diffusivity"]:
            ion_props["diffusivity_data"][("Liq", ion)] = m.fs.props_data[
                "diffusivity"
            ][ion]
        else:
            ion_props["molar_volume_data"][("Liq", ion)] = m.fs.props_data["molar_vol"][
                ion
            ]
        ion_props["mw_data"][ion] = m.fs.props_data["mw"][ion]
        ion_props["charge"][ion] = m.fs.props_data["charge"][ion]
    # NOTE: Use Hayduk Laudie to calculate diffusivity
    ion_props["diffus_calculation"] = "HaydukLaudie"

    m.fs.ix_properties = MCASParameterBlock(**ion_props)

    ix_config = {
        "property_package": m.fs.ix_properties,
        "regenerant": regenerant,
        "hazardous_waste": hazardous_waste,
        "target_component": target_component,
        "reactive_ions": m.list_reactive_ions,
        "number_traps": num_traps,
        "c_trap_min": c_trap_min,
        "resin_data_path": resin_file,
        "resin": resin,
    }

    ix = m.fs.unit_ix = IonExchangeMultiComp(**ix_config)

    # Touch relevant variables
    ix.process_flow.properties_in[0].flow_vol_phase["Liq"]
    ix.process_flow.properties_in[0].pressure
    ix.process_flow.properties_in[0].temperature
    for i in m.fs.set_all:
        ix.process_flow.properties_in[0].flow_mass_phase_comp["Liq", i]
        ix.process_flow.properties_out[0].flow_mass_phase_comp["Liq", i]
        ix.process_flow.properties_in[0].conc_mass_phase_comp["Liq", i]
        ix.process_flow.properties_out[0].conc_mass_phase_comp["Liq", i]

    @m.fs.unit_ix.Expression(
        m.fs.set_reactive_ions, doc="Percentage of recovery for all the REEs"
    )
    def recovery_comp(b, s):
        conc_in = b.process_flow.properties_in[0].conc_mass_phase_comp["Liq", s]
        conc_out = b.process_flow.properties_out[0].conc_mass_phase_comp["Liq", s]
        return ((conc_in - conc_out) / conc_in) * 100


def set_bounds(m):
    """This method adds bounds to relevant variables"""

    ix = m.fs.unit_ix

    ix.bed_diameter.setlb(0.01)
    ix.bed_diameter.setub(100)
    ix.bed_depth.setlb(0.01)
    ix.bed_depth.setub(100)
    ix.service_flow_rate.setlb(1e-10)
    ix.process_flow.properties_in[0.0].visc_k_phase["Liq"].setlb(1e-16)
    ix.process_flow.properties_in[0].flow_vol_phase["Liq"].setlb(1e-16)
    for c in m.fs.set_reactive_ions:
        ix.process_flow.properties_in[0.0].diffus_phase_comp["Liq", c].setlb(1e-16)
        ix.bv_50[c].setlb(1e-3)
        ix.loading_rate.setlb(1e-16)
        ix.process_flow.properties_in[0].flow_mass_phase_comp["Liq", c].setlb(1e-16)
        ix.process_flow.properties_out[0].flow_mass_phase_comp["Liq", c].setlb(1e-16)
        for i in range(1, ix.num_traps + 1):
            ix.c_traps[c, i].setlb(1e-16)


def set_operating_conditions(m, parmest_data=None, target_component=None, resin=None):
    """Set design parameters and operating conditions for the ion exchange
    column.

    """

    ix = m.fs.unit_ix
    pf = ix.process_flow

    # Add water density
    rho_m3 = 1000 * (pyo.units.kg / pyo.units.m**3)

    # Calculate the mass of water based on the given data in ref[1]
    L_solution = pyo.units.convert(50 * pyo.units.cm**3, to_units=pyo.units.L)
    dens = m.fs.data_init["density_solution"] * (pyo.units.kg / pyo.units.m**3)
    kg_solution = dens * pyo.units.convert(L_solution, to_units=pyo.units.m**3)

    # Add volumetric flow, initial concentration for all REEs, and
    # calculate missing parameters needed in the IX model. All values
    # are obtained from ref[1].
    flow_vol = m.fs.data_init["flow_in"] * (pyo.units.m**3 / pyo.units.seconds)

    init_mass_comp = {}
    for c in m.fs.set_reactive_ions:

        # Calculate total and individual mass for REEs
        init_mass_comp[c] = L_solution * pyo.units.convert(
            m.fs.data_init["conc_mass"][c] * (pyo.units.mg / pyo.units.L),
            to_units=pyo.units.kg / pyo.units.L,
        )
    sum_init_mass_comp = sum(init_mass_comp[c] for c in m.list_reactive_ions)  # in kg/L

    # Save the initial concentration for water in the data_init
    # dictionary
    m.fs.data_init["conc_mass"]["H2O"] = pyo.value(
        pyo.units.convert(
            (kg_solution - sum_init_mass_comp) / L_solution,
            to_units=pyo.units.mg / pyo.units.L,
        )
    )

    # Convert initial concentrations from data to kg/m3
    conc_mass_init_conv = {}
    for c in m.fs.set_all:
        conc = m.fs.data_init["conc_mass"][c] * (pyo.units.mg / pyo.units.L)
        conc_mass_init_conv[c] = pyo.units.convert(
            conc, to_units=pyo.units.kg / pyo.units.m**3
        )

    mw = {}
    flow_mol = {}
    for i in m.fs.set_all:
        mw[i] = m.fs.props_data["mw"][i] * (pyo.units.kg / pyo.units.mol)
        if i == "H2O":
            flow_mol[i] = pyo.units.convert(
                (flow_vol * rho_m3) / mw[i], to_units=pyo.units.mol / pyo.units.seconds
            )
        else:
            flow_mol[i] = pyo.units.convert(
                (conc_mass_init_conv[i] * flow_vol) / mw[i],  # when conc_mass in kg/m3
                to_units=pyo.units.mol / pyo.units.seconds,
            )

    # Add inlet pressure, temperature, and flow to IX
    pf.properties_in[0].pressure.fix(101325)
    pf.properties_in[0].temperature.fix(298.15)
    for i in m.fs.set_all:
        pf.properties_in[0].flow_mol_phase_comp["Liq", i].fix(flow_mol[i])

    # Add resin parameters from .json file
    ix.resin_diam.fix(m.fs.resin_data[resin]["diameter"]["value"])
    ix.resin_density.fix(m.fs.resin_data[resin]["density_kgm3"]["value"])

    # Add column design parameters. bed_depth and bed_diameter are
    # taken from ref[1].
    ix.bed_depth.fix(m.fs.data_init["bed_depth"])
    ix.bed_diameter.fix(m.fs.data_init["bed_diam"])
    ix.bed_porosity.fix(0.8)
    ix.number_columns.fix(1)
    ix.number_columns_redundant.fix(1)

    loading_rate0 = parmest_data["loading_rate"][target_component]
    ix.loading_rate.set_value(loading_rate0)

    ebct0 = (
        (np.pi * (m.fs.data_init["bed_diam"] * pyo.units.m / 2) ** 2)
        * m.fs.data_init["bed_depth"]
        * pyo.units.m
        / flow_vol
    )  # in seconds
    ix.ebct.set_value(ebct0)

    # Equilibrium parameters
    for c in m.fs.set_reactive_ions:

        # Fix Freundlich parameters using the data from parameter
        # estimation
        ix.freundlich_n[c].fix(parmest_data["freundlich_n"][c])
        ix.mass_transfer_coeff[c].fix(parmest_data["mass_transfer_coeff"][c])
        ix.bv_50[c].fix(parmest_data["bv_50"][c])
        ix.c_norm[c].fix(0.99)

    return m


def get_comp_list(blk, comp=pyo.Var):
    cs = []
    split_name = blk.name + "."
    skip_list = ["ref", "process_flow", "regeneration"]
    for c in blk.component_objects(comp):
        if any(s in c.name for s in skip_list):
            continue
        cs.append(c.name.split(split_name)[1])
    return cs


def set_scaling(m):
    """Set scaling factors for all variables in flowsheet"""

    ix = m.fs.unit_ix

    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

    sb = ScalerBase()

    for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
        if "temperature" in var.name:
            sb.set_variable_scaling_factor(var, 1e-1, overwrite=True)
        if "pressure" in var.name:
            sb.set_variable_scaling_factor(var, 1e-5)
        if "flow_vol" in var.name:
            sb.set_variable_scaling_factor(var, 1e8)
        if "flow_mol_phase_comp" in var.name:
            if "H2O" in var.name:
                sb.set_variable_scaling_factor(var, 1e4)
            else:
                sb.set_variable_scaling_factor(var, 1e8)

    for v in get_comp_list(ix):
        ixv = getattr(ix, v)
        if ixv.is_indexed():
            idx = [*ixv.index_set()]
            for i in idx:
                if ixv[i].is_fixed():
                    if ixv[i]() == 0:
                        continue
                    sb.set_variable_scaling_factor(ixv[i], 1 / ixv[i]())
        else:
            if ixv.is_fixed():
                if ixv() == 0:
                    continue
                sb.set_variable_scaling_factor(ixv, 1 / ixv())

    for v in m.fs.scale_from_value:
        ixv = getattr(ix, v)
        if ixv.is_indexed():
            idx = [*ixv.index_set()]
            for i in idx:
                if ixv[i]() == 0:
                    continue
                sb.set_variable_scaling_factor(ixv[i], 1 / ixv[i]())
        else:
            if ixv() == 0:
                continue
            sb.set_variable_scaling_factor(ixv, 1 / ixv())

    return m


def initialize_system(
    m,
    solver=None,
    **kwargs,
):

    ix = m.fs.unit_ix

    ix.initialize()

    # Check and raise an error if the degrees of freedom are not 0
    if degrees_of_freedom(m) != 0:
        raise ConfigurationError(
            "The degrees of freedom after building the model are not 0. "
            "You have {} degrees of freedom. "
            "Please check your inputs to ensure a square problem "
            "before initializing the model.".format(degrees_of_freedom(m))
        )


def add_costing(m):

    flow_out = m.fs.unit_ix.process_flow.properties_out[0].flow_vol_phase["Liq"]

    m.fs.costing = IXCosting()
    m.fs.costing.base_currency = pyo.units.USD_2021
    m.fs.costing.base_period = pyo.units.year

    # Add costs related only to IX unit
    m.fs.unit_ix.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=IXCostingData.cost_ion_exchange,
    )

    # Calculate costs of entire process
    m.fs.costing.cost_process()
    m.fs.costing.utilization_factor.fix(1)
    m.fs.costing.aggregate_fixed_operating_cost()
    m.fs.costing.aggregate_variable_operating_cost()
    m.fs.costing.aggregate_capital_cost()
    m.fs.costing.total_capital_cost()
    m.fs.costing.total_operating_cost()

    # Touch relevant cost variable
    m.fs.costing.total_annualized_cost

    # Add costs for REEs. References are: [a]
    # https://www.metal.com/Rare-Earth-Metals/ and [b]
    # https://www.metal.com/price/Rare%20Earth/Rare-Earth-Oxides.
    market_prices = {
        # From ref[a]
        "La": 2.642,
        "Dy": 247.07,
        "Ho": 63.610,
        # From [b]
        "Er": 39.641,
        "Yb": 12.291,
        "Sm": 8.973,
    }
    m.fs.market_price = pyo.Param(
        m.fs.set_reactive_ions,
        initialize=market_prices,
        units=pyo.units.USD_2021 / pyo.units.kg,
        doc="Market price for REEs",
    )
    m.fs.operational_daily_hours = pyo.Param(
        initialize=8, units=pyo.units.hours / pyo.units.day, doc="IX operational hours"
    )
    m.fs.operational_yearly_days = pyo.Param(
        initialize=365, units=pyo.units.day / pyo.units.year, doc="IX operational hours"
    )
    m.fs.ix_lifetime = pyo.Param(
        initialize=15, units=pyo.units.years, doc="IX lifetime"
    )

    @m.fs.unit_ix.Expression()
    def expected_annual_profit_ree(b):
        return (
            sum(
                m.fs.market_price[s]
                * pyo.units.convert(
                    (
                        b.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", s]
                        - b.process_flow.properties_out[0.0].flow_mass_phase_comp[
                            "Liq", s
                        ]
                    ),
                    to_units=pyo.units.kg / pyo.units.hour,
                )
                for s in m.fs.set_reactive_ions
            )
            * m.fs.operational_daily_hours
            * m.fs.operational_yearly_days
        )

    m.fs.unit_ix.target_breakthrough_time.setlb(1e-3)
    m.fs.unit_ix.target_breakthrough_time.setub(1e10)


def run_optimization(m, target_component=None):
    """This method unfixes variables and add constraints to solve an
    optimization problem

    """

    ix = m.fs.unit_ix

    for c in m.fs.set_reactive_ions:
        ix.c_norm[c].unfix()

    ix.bed_depth.unfix()
    ix.bed_diameter.unfix()

    # For this example, we optimize the model to have an effluent with
    # a very small concentration of the multiple REEs (c_norm closer
    # to 1 or final concentration closer to initial concentration)
    @m.fs.unit_ix.Constraint(m.fs.set_reactive_ions)
    def components_specifications(b, c):
        if c == target_component:
            return b.c_norm[c] == 0.999
        else:
            return b.c_norm[c] >= 0.9

    @m.fs.unit_ix.costing.Expression()
    def annualized_capital_cost(b):
        return b.capital_cost / m.fs.ix_lifetime

    m.fs.obj = pyo.Objective(
        expr=(
            ix.costing.annualized_capital_cost
            + ix.costing.fixed_operating_cost
            - ix.expected_annual_profit_ree
        )
    )


def plot_traps(
    m, results=None, curve_data=None, output_plot=None, save_plot=None, show_plot=None
):

    distinct_colors = [
        "#1f77b4",  # Blue
        "#ff7f0e",  # Orange
        "#2ca02c",  # Green
        "#d62728",  # Red
        "#9467bd",  # Purple
        "#8c564b",  # Brown
    ]

    # Define a list of distinct markers
    distinct_markers = [
        "o",  # Circle
        "s",  # Square
        "D",  # Diamond
        "^",  # Triangle
        "v",  # Inverted Triangle
        "<",  # Left Triangle
    ]

    ix = m.fs.unit_ix
    color = ["b", "k"]
    num_traps = ix.num_traps

    plt.figure(figsize=(14, 12))
    plt.minorticks_on()
    plt.grid(linestyle=":", which="both", color="gray", alpha=0.50)

    for idx, c in enumerate(m.fs.set_reactive_ions):

        comp = curve_data["compound"]
        x_orig_values = curve_data[comp == c]["bv"].values
        y_orig_values = curve_data[comp == c]["c_norm"].values

        plt.plot(
            x_orig_values,
            y_orig_values,
            marker=distinct_markers[idx % len(distinct_markers)],
            ms=8,
            linestyle="",
            linewidth=0.5,
            color=distinct_colors[idx % len(distinct_colors)],
            alpha=0.8,
            label=f"Literature {c}",
        )

    for idx, c in enumerate(m.fs.set_reactive_ions):

        bvs_values = []
        c_breakthru_values = []
        for i in range(1, num_traps + 1):
            bvs_values.append(pyo.value(ix.bv_traps[c, i]))
            c_breakthru_values.append(pyo.value(ix.c_traps[c, i]))

        plt.plot(
            bvs_values,
            c_breakthru_values,
            linestyle="-",
            linewidth=1,
            color=distinct_colors[idx % len(distinct_colors)],
            alpha=0.6,
            label=f"Calculated {c}",
        )

    plt.xlabel("Bed Volume (BV)", fontsize=16)
    plt.ylabel("C/C0", fontsize=16)
    plt.xticks(np.arange(0, 241, 20), fontsize=14)
    plt.yticks(np.arange(0, 1.1, 0.2), fontsize=14)
    plt.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.80, 0.75),
        ncol=1,
        fontsize=12,
    )
    if save_plot:
        plt.savefig("breakthru_trapezoidal_all_components.png", bbox_inches="tight")
    if show_plot:
        plt.show()
    plt.close()


def model_solve(m, solver=None):

    solver = get_solver()
    results = solver.solve(
        m,
        tee=True,
        symbolic_solver_labels=True,
        options={
            "halt_on_ampl_error": "yes",
            "max_iter": 1000,
        },
    )

    return results


if __name__ == "__main__":

    m = main()
