#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Oxalate Precipitator Unit Model
===================================

Author: Alejandro Garciadiego, Bo-Xun Wang

The Precipitator Unit Model represents an Equilibrium reactor unit model with the equilibrium equation derived from Minteq data.

Configuration Arguments
-----------------------

The precipitator unit model needs an aqueous property package which includes stoichiometric values for solids being
created in the precipitator and the parameters used in the equilibrium equation.

Model Structure
---------------

The Precitator unit model has hard coded stream names (``aqueous`` and ``precipitate`` respectively). The Precipitator 
model also has one inlet and two outlets named ``aqueous_inlet``, ``aqueous_outlet`` and ``precipitate_outlet`` respectively.

Additional Constraints
----------------------

The Precipitator unit adds one additional constraint to define the conversion.

.. math:: Conversion_{c} = \exp(-\frac{\epsilon}{Oxalic Acid Dosage}^{n_{DA}})

where :math:`Conversion` is the conversion of component c, :math:`\epsilon` and :math:`n_{DA}` are the parameters 
estimated based on Minteq data, :math:`Oxalic Acid Dosage` is the amount of oxalic acid added into the precipitator.

"""

from pyomo.common.config import Bool, ConfigDict, ConfigValue

from pyomo.environ import (
    Var,
    Block,
    log,
    NonNegativeReals,
    units as pyunits,
)

from pyomo.network import Port

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale

from idaes.models.unit_models.mscontactor import MSContactor
from idaes.core.initialization import ModularInitializerBase


# -----------------------------------------------------------------------------
# Precipitator unit model
class OxalatePrecipitatorInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer  for the Oxalate Precipitator unit model.

    This routine calls the initializer for the internal MSContactor model.

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

    def initialize_main_model(
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
        model.mscontactor.heterogeneous_reaction_extent.fix()

        msc_init = model.mscontactor.default_initializer(
            ssc_solver_options=self.config.ssc_solver_options,
            calculate_variable_options=self.config.calculate_variable_options,
        )

        try:
            msc_init.initialize(model.mscontactor)
        except:
            pass

        model.mscontactor.heterogeneous_.unfix()

        solver = self._get_solver()
        results = solver.solve(model, tee=True)

        return results


StreamCONFIG = ConfigDict()
StreamCONFIG.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for given stream",
        doc="""Property parameter object used to define property calculations for given stream,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
    ),
)
StreamCONFIG.declare(
    "property_package_args",
    ConfigDict(
        implicit=True,
        description="Dict of arguments to use for constructing property package",
        doc="""A ConfigDict with arguments to be passed to property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
StreamCONFIG.declare(
    "has_energy_balance",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include energy balance for stream. Default=True.",
    ),
)
StreamCONFIG.declare(
    "has_pressure_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include pressure balance for stream. Default=True.",
    ),
)


@declare_process_block_class("OxalatePrecipitator")
class OxalatePrecipitatorData(UnitModelBlockData):
    """
    Oxalate Precipitator Unit Model Class
    """

    # Set default initializer
    default_initializer = OxalatePrecipitatorInitializer

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "liquid_phase",
        StreamCONFIG(
            description="Liquid phase properties",
        ),
    )
    CONFIG.declare(
        "solid_phase",
        StreamCONFIG(
            description="Solid phase properties",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            # TODO: Add a domain validator for this
            description="Heterogeneous reaction package for precipitaction.",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for heterogeneous reaction package for precipitaction.",
        ),
    )
    CONFIG.declare(
        "number_of_2s",
        ConfigValue(
            default=1, domain=int, description="Number of tanks in precipitaction"
        ),
    )

    def build(self):
        """
        Build method for OxalatePrecipitator unit model.
        """
        super().build()

        self.mscontactor = MSContactor(
            number_of_finite_elements=self.config.number_of_tanks,
            streams={
                "liquid": {
                    "property_package": self.config.liquid_phase.property_package,
                    "property_package_args": self.config.liquid_phase.property_package_args,
                    "has_energy_balance": self.config.liquid_phase.has_energy_balance,
                    "has_pressure_balance": self.config.liquid_phase.has_pressure_balance,
                },
                "solid": {
                    "property_package": self.config.solid_phase.property_package,
                    "property_package_args": self.config.solid_phase.property_package_args,
                    "has_energy_balance": self.config.solid_phase.has_energy_balance,
                    "has_pressure_balance": self.config.solid_phase.has_pressure_balance,
                },
            },
            heterogeneous_reactions=self.config.reaction_package,
            heterogeneous_reactions_args=self.config.reaction_package_args,
        )

        # Get units of measurement from MSContactor
        flow_basis = self.mscontactor.flow_basis
        uom = self.mscontactor.uom

        self.hydraulic_retention_time = Var(
            self.flowsheet().time,
            initialize=2,
            domain=NonNegativeReals,
            units=pyunits.h,
            doc="Hydraulic retention time",
        )

        self.volume = Var(
            self.flowsheet().time,
            initialize=1800,
            domain=NonNegativeReals,
            units=pyunits.l,
            doc="Volume of precipitator",
        )

        self.conversion = Var(
            self.config.reaction_package.reaction_idx,
            initialize=0.5,
            units=pyunits.dimensionless,
            bounds=(1e-20, 0.999999),
        )

        # Create unit level Ports
        self.aqueous_inlet = Port(extends=self.mscontactor.liquid_inlet)
        self.aqueous_outlet = Port(extends=self.mscontactor.liquid_outlet)
        self.precipitate_outlet = Port(extends=self.mscontactor.solid_outlet)

        @self.Constraint(self.flowsheet().time, doc="Hydraulic retention time equation")
        def eq_hydraulic_retention(blk, t):
            return (
                self.hydraulic_retention_time[t] * self.aqueous_inlet.flow_vol[t]
                == self.volume[t]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.reaction_package.reaction_idx,
            doc="Reaction extent constraint",
        )
        def heterogeneous_reaction_extent_constraint(blk, t, s, r):
            return blk.mscontactor.heterogeneous_reaction_extent[t, s, r] == (
                blk.mscontactor.heterogeneous_reactions[t, s].reaction_rate[r]
                - (
                    self.conversion[r]
                    * blk.mscontactor.liquid_inlet_state[t].flow_mol_comp[
                        blk.mscontactor.config.streams.solid.property_package.react[r]
                    ]
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.reaction_package.reaction_idx,
            doc="conversion constraint",
        )
        def conversion_constraint(blk, t, s, r):
            if r == "Ca(C2O4)(s)":
                return self.conversion[r] == 1e-6
            else:
                return log(self.conversion[r]) == (
                    -(
                        (blk.config.reaction_package.E_D[r])
                        ** blk.config.reaction_package.N_D[r]
                    )
                ) / (
                    (
                        (
                            (blk.aqueous_inlet.conc_mass_comp[0, "H2C2O4"])
                            / (1000 * pyunits.mg / pyunits.l)
                        )
                        ** blk.config.reaction_package.N_D[r]
                    )
                )

        @self.Constraint(self.flowsheet().time, doc="temperature equation")
        def temp_constraint(blk, t):
            return (
                blk.mscontactor.solid_inlet_state[t].temperature
                == blk.mscontactor.solid_outlet.temperature[t]
            )

        @self.Constraint(self.flowsheet().time, doc="liquid temperature equation")
        def liq_temp_constraint(blk, t):
            return blk.aqueous_inlet.temperature[t] == blk.aqueous_outlet.temperature[t]

        @self.Constraint(self.flowsheet().time, doc="pressure equation")
        def press_constraint(blk, t):
            return blk.aqueous_inlet.pressure[t] == blk.aqueous_outlet.pressure[t]

        @self.Constraint(
            self.flowsheet().time,
            self.config.solid_phase.property_package.component_list,
            doc="Initial solids",
        )
        def init_solid_constraint(blk, t, r):
            return (
                blk.mscontactor.solid_inlet_state[t].flow_mol_comp[r]
                == 1e-9 * pyunits.mole / pyunits.hour
            )

        iscale.set_scaling_factor(self.hydraulic_retention_time, 1e0)
        iscale.set_scaling_factor(self.conversion, 1e1)
        iscale.set_scaling_factor(self.mscontactor.heterogeneous_reaction_extent, 1e3)
        iscale.set_scaling_factor(
            self.mscontactor.solid_heterogeneous_reactions_generation, 1e3
        )
        iscale.set_scaling_factor(
            self.mscontactor.liquid_heterogeneous_reactions_generation, 1e3
        )
        iscale.set_scaling_factor(self.volume, 1e-3)

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        param_dict = {}
        var_dict["Unit Volume"] = self.volume[time_point]
        var_dict["Hydraulic Retention Time"] = self.hydraulic_retention_time[time_point]
        var_dict["Unit Height"] = self.height
        var_dict["Unit Diameter"] = self.diameter
        expr_dict["Surface Area"] = self.surface_area
        return {"vars": var_dict, "params": param_dict, "exprs": expr_dict}
