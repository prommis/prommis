#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Demonstration flowsheet for LeachTrain unit model using
parameters and data for West Kentucky No. 13 coal refuse.

Authors: Andrew Lee, Arkoprabho Dasgupta, Akintomiwa Ojo, Douglas Allan
"""
import matplotlib.pyplot as plt
from pyomo.common.config import Bool, ConfigDict, ConfigValue, In, PositiveInt
from pyomo.environ import (
    Block,
    ComponentMap,
    ConcreteModel,
    TransformationFactory,
    units,
    Var,
    Param,
)
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import FlowsheetBlockData, declare_process_block_class
from idaes.core.initialization import ModularInitializerBase
from idaes.core.util import to_json, from_json

from prommis.leaching.leach_train import LeachingTrain, LeachingTrainInitializer
from prommis.leaching.leach_reactions import CoalRefuseLeachingReactionParameterBlock
from prommis.properties.coal_refuse_properties import CoalRefuseParameters
from prommis.properties.sulfuric_acid_leaching_properties import (
    SulfuricAcidLeachingParameters,
)


def copy_first_steady_state(m):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    """

    _, time_vars = flatten_dae_components(m, m.fs.time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in m.fs.time:
            if t == m.fs.time.first():
                continue
            else:
                var[t].value = var[m.fs.time.first()].value


class CocurrentSlurryLeachingFlowsheetInitializer(ModularInitializerBase):
    """
    This is an initializer for the CocurrentSlurryLeachingFlowsheet.
    This method calls the
    """

    CONFIG = ModularInitializerBase.CONFIG()

    CONFIG.declare(
        "ssc_solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments for solver calls by ssc_solver",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "scc_solver method.",
        ),
    )

    def initialization_routine(
        self,
        model: Block,
    ):
        """
        Initialization routine for MSContactor Blocks.

        Args:
            model: model to be initialized

        Returns:
            None
        """
        # Initialize MSContactor
        leach_init = model.leach.default_initializer(**self.config)
        leach_init.initialize(model.leach)

        solver = self._get_solver()
        results = solver.solve(model)

        return results


@declare_process_block_class("CocurrentSlurryLeachingFlowsheet")
class CocurrentSlurryLeachingFlowsheetData(FlowsheetBlockData):
    """
    Flowsheet containing leaching train for critical mineral recovery
    from West Kentucky No. 13 coal refuse.
    """

    default_initializer = CocurrentSlurryLeachingFlowsheetInitializer

    CONFIG = FlowsheetBlockData.CONFIG()

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Determines whether to create holdup terms in leach train.",
        ),
    )
    CONFIG.declare(
        "number_of_tanks",
        ConfigValue(
            default=1,
            domain=PositiveInt,
            description="Number of tanks in leaching train.",
        ),
    )
    CONFIG.declare(
        "dynamic_transformation_method",
        ConfigValue(
            default="finite_difference", domain=In(["finite_difference", "collocation"])
        ),
    )
    CONFIG.declare(
        "dynamic_transformation_scheme",
        ConfigValue(
            default="BACKWARD",
        ),
    )

    def build(self):
        super().build()
        self._add_parameter_blocks()
        self._add_units()
        if self.config.dynamic:
            TransformationFactory(
                f"dae.{self.config.dynamic_transformation_method}"
            ).apply_to(
                self,
                nfe=len(self.time) - 1,
                wrt=self.time,
                scheme=self.config.dynamic_transformation_scheme,
            )
        self._set_unit_geometry()
        self._set_inputs()

    def _add_parameter_blocks(self):
        self.leach_soln = SulfuricAcidLeachingParameters()
        self.coal = CoalRefuseParameters()
        self.leach_rxns = CoalRefuseLeachingReactionParameterBlock()

    def _add_units(self):
        self.leach = LeachingTrain(
            number_of_tanks=self.config.number_of_tanks,
            liquid_phase={
                "property_package": self.leach_soln,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            solid_phase={
                "property_package": self.coal,
                "has_energy_balance": False,
                "has_pressure_balance": False,
            },
            reaction_package=self.leach_rxns,
            has_holdup=self.config.has_holdup,
        )

    def _set_unit_geometry(self):
        self.leach.volume.fix(100 * units.gallon)

        if self.leach.config.has_holdup:
            self.leach.liquid_solid_residence_time_ratio.fix(1 / 32)

    def _set_inputs(self):
        """
        Set inlet conditions to leach reactor based on one case study from
        University of Kentucky pilot plant study.
        """
        # Liquid feed state
        self.leach.liquid_inlet.flow_vol.fix(224.3 * units.L / units.hour)
        self.leach.liquid_inlet.conc_mass_comp.fix(1e-10 * units.mg / units.L)
        self.leach.liquid_inlet.conc_mass_comp[:, "H2O"].fix(1e6 * units.mg / units.L)
        self.leach.liquid_inlet.conc_mass_comp[:, "H"].fix(
            2 * 0.05 * 1e3 * units.mg / units.L
        )
        self.leach.liquid_inlet.conc_mass_comp[:, "HSO4"].fix(1e-8 * units.mg / units.L)
        self.leach.liquid_inlet.conc_mass_comp[:, "SO4"].fix(
            0.05 * 96e3 * units.mg / units.L
        )

        # Solid feed state
        self.leach.solid_inlet.flow_mass.fix(22.68 * units.kg / units.hour)
        self.leach.solid_inlet.mass_frac_comp[:, "inerts"].fix(
            0.6952 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Al2O3"].fix(
            0.237 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Fe2O3"].fix(
            0.0642 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "CaO"].fix(
            3.31e-3 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Sc2O3"].fix(
            2.77966e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Y2O3"].fix(
            3.28653e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "La2O3"].fix(
            6.77769e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Ce2O3"].fix(
            0.000156161 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Pr2O3"].fix(
            1.71438e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Nd2O3"].fix(
            6.76618e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Sm2O3"].fix(
            1.47926e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Gd2O3"].fix(
            1.0405e-05 * units.kg / units.kg
        )
        self.leach.solid_inlet.mass_frac_comp[:, "Dy2O3"].fix(
            7.54827e-06 * units.kg / units.kg
        )

    def fix_initial_conditions(self, t=None):
        if t is None:
            t = self.time.first()

        self.leach.mscontactor.liquid[t, :].flow_vol.fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["H"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["HSO4"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Cl"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Sc"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Y"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["La"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Ce"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Pr"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Nd"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Sm"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Gd"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Dy"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Al"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Ca"].fix()
        self.leach.mscontactor.liquid[t, :].conc_mass_comp["Fe"].fix()

        self.leach.mscontactor.solid[t, :].mass_frac_comp["inerts"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Y2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["La2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Ce2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Pr2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Nd2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Sm2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Gd2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Dy2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Al2O3"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["CaO"].fix()
        self.leach.mscontactor.solid[t, :].mass_frac_comp["Fe2O3"].fix()
        self.leach.mscontactor.solid[t, :].flow_mass.fix()

        self.leach.mscontactor.liquid_inherent_reaction_extent[t, :, "Ka2"].fix()

    def scale_model(self):
        """
        Apply scaling factors to improve solver performance.
        """

        solid_scaler = self.leach.mscontactor.solid.default_scaler()
        solid_scaler.default_scaling_factors["flow_mass"] = 1 / 22.68

        liquid_scaler = self.leach.mscontactor.liquid.default_scaler()
        liquid_scaler.default_scaling_factors["flow_vol"] = 1 / 224.3
        liquid_scaler.default_scaling_factors["conc_mass_comp[Ce]"] = 1 / 5
        liquid_scaler.default_scaling_factors["conc_mass_comp[Nd]"] = 1 / 2
        liquid_scaler.default_scaling_factors["conc_mass_comp[La]"] = 1
        liquid_scaler.default_scaling_factors["conc_mass_comp[SO4]"] = 1e-3

        submodel_scalers = ComponentMap()
        submodel_scalers[self.leach.mscontactor.liquid_inlet_state] = liquid_scaler
        submodel_scalers[self.leach.mscontactor.liquid] = liquid_scaler
        submodel_scalers[self.leach.mscontactor.solid_inlet_state] = solid_scaler
        submodel_scalers[self.leach.mscontactor.solid] = solid_scaler

        scaler_obj = self.leach.default_scaler()
        scaler_obj.default_scaling_factors["liquid_phase_fraction"] = 1
        scaler_obj.default_scaling_factors["solid_phase_fraction"] = 1
        scaler_obj.scale_model(self.leach, submodel_scalers=submodel_scalers)


def plot_results(m, perturb_time):
    # Plotting the results for REE oxides
    REE_set = m.fs.coal.component_list - ["inerts", "Al2O3", "Fe2O3", "CaO"]
    for e in REE_set:
        plt.plot(m.fs.time, m.fs.leach.recovery[:, e]())
    plt.legend(REE_set)
    plt.xlabel("Time (h)")
    plt.ylabel("Recovery %")
    plt.axvline(
        x=perturb_time,
        color="r",
        linestyle="--",
        label=f"Perturbation at t={perturb_time}h",
    )
    plt.title("REE oxide recovery variation wrt time, with perturbation at t=12 ")
    plt.figure()
    plt.plot(m.fs.time, m.fs.leach.mscontactor.solid[:, :].flow_mass())
    plt.title("Solid mass flow rate variation wrt time")
    plt.xlabel("Time (h)")
    plt.ylabel("Solid mass flow rate (kg/h)")


def create_one_tank_json():
    m = ConcreteModel()
    m.fs = CocurrentSlurryLeachingFlowsheet(
        has_holdup=True,
        number_of_tanks=1,
    )
    m.fs.scale_model()
    init_obj = m.fs.default_initializer()
    init_obj.initialize(m.fs)
    to_json(m, fname="leaching_one_tank.json", human_read=True)
    return m


def create_two_tanks_json():
    m = ConcreteModel()
    m.fs = CocurrentSlurryLeachingFlowsheet(
        has_holdup=True,
        number_of_tanks=2,
    )
    m.fs.leach.volume[:].fix(50 * units.gallon)
    m.fs.scale_model()
    init_obj = m.fs.default_initializer()
    init_obj.initialize(m.fs)
    to_json(m, fname="leaching_two_tanks.json", human_read=True)
    return m


# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    m = create_one_tank_json()
    # m = create_two_tanks_json()

    m.fs.leach.liquid_outlet.display()
    m.fs.leach.solid_outlet.display()

    m.fs.leach.report()
