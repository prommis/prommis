#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Demonstration flowsheet for LeachTrain unit model using
parameters and data for West Kentucky No. 13 coal refuse.

Authors: Andrew Lee
"""

from pyomo.environ import (
    ComponentMap,
    ConcreteModel,
    SolverFactory,
    TransformationFactory,
    units,
)

from idaes.core import FlowsheetBlock
from idaes.core.util import to_json
from idaes.core.util.scaling import set_scaling_factor

from prommis.leaching.leach_train import LeachingTrain, LeachingTrainInitializer
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)


def build_model(has_holdup=True, n_tanks=1):
    """
    Method to build a single stage leaching system using data for
    West Kentucky No. 13 coal refuse.
    Args:
        has_holdup: Boolean flag to determine whether to create
            material holdup terms
        n_tanks: Number of tanks in leaching train
    Returns:
        m: ConcreteModel containing the leaching flowsheet
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.leach_soln = SulfuricAcidLeachingParameters()
    m.fs.coal = CoalRefuseParameters()
    m.fs.leach_rxns = CoalRefuseLeachingReactionParameterBlock()

    m.fs.leach = LeachingTrain(
        number_of_tanks=n_tanks,
        liquid_phase={
            "property_package": m.fs.leach_soln,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        solid_phase={
            "property_package": m.fs.coal,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        reaction_package=m.fs.leach_rxns,
        has_holdup=has_holdup,
    )

    return m


def set_inputs(m):
    """
    Set inlet conditions to leach reactor based on one case study from
    University of Kentucky pilot plant study.
    """
    # Liquid feed state
    m.fs.leach.liquid_inlet.flow_vol.fix(224.3 * units.L / units.hour)
    m.fs.leach.liquid_inlet.conc_mass_comp.fix(1e-10 * units.mg / units.L)
    m.fs.leach.liquid_inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)
    m.fs.leach.liquid_inlet.conc_mass_comp[0, "H"].fix(
        2 * 0.05 * 1e3 * units.mg / units.L
    )
    m.fs.leach.liquid_inlet.conc_mass_comp[0, "HSO4"].fix(1e-8 * units.mg / units.L)
    m.fs.leach.liquid_inlet.conc_mass_comp[0, "SO4"].fix(
        0.05 * 96e3 * units.mg / units.L
    )

    # Solid feed state
    m.fs.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "inerts"].fix(0.6952 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.237 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.0642 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "CaO"].fix(3.31e-3 * units.kg / units.kg)
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sc2O3"].fix(
        2.77966e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Y2O3"].fix(
        3.28653e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "La2O3"].fix(
        6.77769e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Ce2O3"].fix(
        0.000156161 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Pr2O3"].fix(
        1.71438e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Nd2O3"].fix(
        6.76618e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Sm2O3"].fix(
        1.47926e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Gd2O3"].fix(
        1.0405e-05 * units.kg / units.kg
    )
    m.fs.leach.solid_inlet.mass_frac_comp[0, "Dy2O3"].fix(
        7.54827e-06 * units.kg / units.kg
    )

    m.fs.leach.volume.fix(100 * units.gallon)

    if m.fs.leach.config.has_holdup:
        m.fs.leach.liquid_solid_residence_time_ratio.fix(1 / 32)


def set_scaling(m):
    """
    Apply scaling factors to improve solver performance.
    """

    solid_scaler = m.fs.leach.mscontactor.solid.default_scaler()
    solid_scaler.default_scaling_factors["flow_mass"] = 1 / 22.68

    liquid_scaler = m.fs.leach.mscontactor.liquid.default_scaler()
    liquid_scaler.default_scaling_factors["flow_vol"] = 1 / 224.3
    liquid_scaler.default_scaling_factors["conc_mass_comp[Ce]"] = 1 / 5
    liquid_scaler.default_scaling_factors["conc_mass_comp[Nd]"] = 1 / 2
    liquid_scaler.default_scaling_factors["conc_mass_comp[La]"] = 1
    liquid_scaler.default_scaling_factors["conc_mass_comp[SO4]"] = 1e-3

    submodel_scalers = ComponentMap()
    submodel_scalers[m.fs.leach.mscontactor.liquid_inlet_state] = liquid_scaler
    submodel_scalers[m.fs.leach.mscontactor.liquid] = liquid_scaler
    submodel_scalers[m.fs.leach.mscontactor.solid_inlet_state] = solid_scaler
    submodel_scalers[m.fs.leach.mscontactor.solid] = solid_scaler

    scaler_obj = m.fs.leach.default_scaler()
    scaler_obj.default_scaling_factors["liquid_phase_fraction"] = 1
    scaler_obj.default_scaling_factors["solid_phase_fraction"] = 1
    scaler_obj.scale_model(m.fs.leach, submodel_scalers=submodel_scalers)


# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    # Call build model function
    m = build_model(has_holdup=True, n_tanks=1)
    set_inputs(m)
    set_scaling(m)

    # Initialize model
    # This is likely to fail to converge, but gives a good enough starting point
    initializer = LeachingTrainInitializer()
    initializer.initialize(m.fs.leach)

    # Solve scaled model
    solver = SolverFactory("ipopt_v2")
    solver.solve(m, tee=True)

    # Store steady state values in a json file
    to_json(m, fname="leaching.json", human_read=True)

    # Display some results
    m.fs.leach.liquid_outlet.display()
    m.fs.leach.solid_outlet.display()

    m.fs.leach.report()

    # Call build model function
    m2 = build_model(has_holdup=True, n_tanks=2)
    set_inputs(m2)
    m2.fs.leach.volume[:].set_value(50 * units.gallon)
    set_scaling(m2)

    # Initialize model
    # This is likely to fail to converge, but gives a good enough starting point
    initializer = LeachingTrainInitializer()
    initializer.initialize(m2.fs.leach)

    # Solve scaled model
    solver = SolverFactory("ipopt_v2")
    solver.solve(m2, tee=True)

    # Store steady state values in a json file
    to_json(m2, fname="leaching2.json", human_read=True)

    # Display some results
    m2.fs.leach.liquid_outlet.display()
    m2.fs.leach.solid_outlet.display()

    m2.fs.leach.report()
