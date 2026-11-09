#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import (
    ComponentMap,
    ConcreteModel,
    units,
)

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver


from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)
from prommis.properties.mixed_acid_properties import (
    get_aliases,
    MixedAcidParameterBlock,
)
from prommis.solvent_extraction.ree_og_distribution import (
    REESolExOgParameters,
    ree_list,
)
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)


def build_model(dosage, number_of_stages, has_holdup, use_mixed_acid=False):
    """
    Method to build a steady state model for solvent extraction
    Args:
        dosage: Percentage dosage of extractant to the system.
        number_of_stages: Number of stages in the model.
        has_holdup: Boolean flag about whether or not to create terms
            associated with material holdup and hydrostatic pressure
        use_mixed_acid: Boolean flag to use mixed acid properties instead
            of the old sulfuric acid properties.
    Returns:
        m: ConcreteModel object with the solvent extraction system.
    """

    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_o = REESolExOgParameters()
    if use_mixed_acid:
        m.fs.leach_soln = MixedAcidParameterBlock(include_sulfates=True)
        m.fs.reaxn = SolventExtractionReactions(
            aqueous_aliases=get_aliases(include_sulfates=True)
        )
    else:
        m.fs.leach_soln = SulfuricAcidLeachingParameters()
        m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.solex = SolventExtraction(
        number_of_finite_elements=number_of_stages,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": not has_holdup,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": not has_holdup,
        },
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=has_holdup,
        create_hydrostatic_pressure_terms=has_holdup,
    )

    return m


def set_inputs(m, dosage, has_holdup, use_mixed_acid=False):
    """
    Set inlet conditions to the solvent extraction model and fixing the parameters
    of the model.
    Args:
        m: ConcreteModel object with the solvent extraction system.
        dosage: Percentage dosage of extractant to the system.
        has_holdup: Boolean flag about whether or not to fix the terms
            associated with material holdup and hydrostatic pressure
        use_mixed_acid: Boolean flag to use the component names of the
            mixed acid properties instead of the old sulfuric acid properties.
    Returns:
        None

    """

    if has_holdup:
        m.fs.solex.mscontactor.volume[:].fix(0.4 * units.m**3)
        m.fs.solex.area_cross_stage[:] = 1
        m.fs.solex.elevation[:] = 0

    aqueous_inlet_comp = {
        "H2O": 1e6,
        "H": 10.75,
        "SO4": 100,
        "HSO4": 1e4,
        "Al": 422.375,
        "Ca": 109.542,
        "Cl": 1e-7,
        "Fe": 688.266,
        "Sc": 0.032,
        "Y": 0.124,
        "La": 0.986,
        "Ce": 2.277,
        "Pr": 0.303,
        "Nd": 0.946,
        "Sm": 0.097,
        "Gd": 0.2584,
        "Dy": 0.047,
    }
    if use_mixed_acid:
        for j1, j2 in get_aliases(include_sulfates=True).items():
            m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp[j2].fix(
                aqueous_inlet_comp[j1]
            )
    else:
        for idx, val in aqueous_inlet_comp.items():
            m.fs.solex.mscontactor.aqueous_inlet_state[:].conc_mass_comp[idx].fix(val)

    m.fs.solex.mscontactor.aqueous_inlet_state[:].flow_vol.fix(
        62.01 * units.L / units.h
    )
    m.fs.solex.mscontactor.aqueous_inlet_state[:].pressure.fix(101300)
    m.fs.solex.mscontactor.aqueous[:, :].temperature.fix(305.15 * units.K)
    m.fs.solex.mscontactor.aqueous_inlet_state[:].temperature.fix(305.15 * units.K)

    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Kerosene"].fix(820e3)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Al_o"].fix(1.267e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ca_o"].fix(2.684e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Fe_o"].fix(2.873e-6)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sc_o"].fix(1.734)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Y_o"].fix(2.179e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["La_o"].fix(0.000105)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Ce_o"].fix(0.00031)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Pr_o"].fix(3.711e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Nd_o"].fix(0.000165)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Sm_o"].fix(1.701e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Gd_o"].fix(3.357e-5)
    m.fs.solex.mscontactor.organic_inlet_state[:].conc_mass_comp["Dy_o"].fix(8.008e-6)

    m.fs.solex.mscontactor.organic_inlet_state[:].flow_vol.fix(62.01)
    m.fs.solex.mscontactor.organic_inlet_state[:].pressure.fix(101300)
    m.fs.solex.mscontactor.organic[:, :].temperature.fix(305.15 * units.K)
    m.fs.solex.mscontactor.organic_inlet_state[:].temperature.fix(305.15 * units.K)


def scale_model(m):
    """
    Apply scaling factors to improve solver performance.
    """

    aqueous_scaler = m.fs.solex.mscontactor.aqueous.default_scaler()
    aqueous_scaler.default_scaling_factors["flow_vol"] = 1 / 62.01

    organic_scaler = m.fs.solex.mscontactor.organic.default_scaler()
    organic_scaler.default_scaling_factors["flow_vol"] = 1 / 62.01
    for ree in ree_list:
        if ree == "Sc_o":
            organic_scaler.default_scaling_factors[f"conc_mass_comp[{ree}]"] = 1
        else:
            organic_scaler.default_scaling_factors[f"conc_mass_comp[{ree}]"] = 100

    submodel_scalers = ComponentMap()
    submodel_scalers[m.fs.solex.mscontactor.aqueous_inlet_state] = aqueous_scaler
    submodel_scalers[m.fs.solex.mscontactor.aqueous] = aqueous_scaler
    submodel_scalers[m.fs.solex.mscontactor.organic_inlet_state] = organic_scaler
    submodel_scalers[m.fs.solex.mscontactor.organic] = organic_scaler

    scaler_obj = m.fs.solex.default_scaler(
        max_variable_scaling_factor=1e12,
        max_constraint_scaling_factor=1e12,
        max_expression_scaling_hint=1e12,
        min_variable_scaling_factor=1e-12,
        min_constraint_scaling_factor=1e-12,
        min_expression_scaling_hint=1e-12,
    )
    scaler_obj.scale_model(m.fs.solex, submodel_scalers=submodel_scalers)


def model_buildup_and_set_inputs(
    dosage, number_of_stages, has_holdup, use_mixed_acid=False
):
    """
    A function to build up the solvent extraction model and set inlet streams
    to the model.
    Args:
        dosage: Percentage dosage of extractant to the system.
        number_of_stages: Number of stages in the model.
        has_holdup: Boolean flag about whether or not to create terms
            associated with material holdup and hydrostatic pressure
        use_mixed_acid: Boolean flag to use mixed acid properties instead
            of the old sulfuric acid properties.
    Returns:
        m: ConcreteModel object with the solvent extraction system.
    """
    m = build_model(
        dosage, number_of_stages, has_holdup=has_holdup, use_mixed_acid=use_mixed_acid
    )
    set_inputs(m, dosage, has_holdup=has_holdup, use_mixed_acid=use_mixed_acid)
    scale_model(m)

    return m


def initialize_steady_model(m):
    """
    A function to initialize the solvent extraction model with the default initializer
    after setting input conditions.
    Args:
        m: ConcreteModel object with the solvent extraction system.
    Returns:
        None
    """
    initializer = m.fs.solex.default_initializer()
    initializer.initialize(m.fs.solex)


def solve_model(m):
    """
    A function to solve the initialized solvent extraction model.
    Args:
        m: ConcreteModel object with the solvent extraction system.
    Returns:
        None
    """
    solver = get_solver("ipopt_v2")
    results = solver.solve(m, tee=True)
    return results


def main(dosage, number_of_stages, has_holdup, used_mixed_acid=False):
    """
    The main function used to build a solvent extraction model, set inlets to the model,
    initialize the model, solve the model and export the results to a json file.
    Args:
        dosage: Percentage dosage of extractant to the system.
        number_of_stages: Number of stages in the model.
        has_holdup: Boolean flag about whether or not to create
            material holdup and hydrostatic pressure terms
    Returns:
        m: ConcreteModel object with the solvent extraction system.
    """
    m = model_buildup_and_set_inputs(
        dosage, number_of_stages, has_holdup, used_mixed_acid
    )
    initialize_steady_model(m)
    results = solve_model(m)

    return m, results


dosage = 5
number_of_stages = 3

if __name__ == "__main__":
    m, results = main(dosage, number_of_stages, has_holdup=False, used_mixed_acid=True)
