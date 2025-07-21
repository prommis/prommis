#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Solvent Extraction Model

========================

Author: Arkoprabho Dasgupta

The Solvent Extraction unit model is used to perform the solvent extraction unit operation.
It represents a series of tanks, referred to as stages, through which the aqueous and organic
phases are passed, and the desired components are extracted subsequently.

Configuration Arguments
-----------------------

The user must specify the following configurations in a solvent extraction model to be able to
use it.

The user must specify the aqueous feed input in the ``aqueous_stream`` configuration, with a
configuration that describes the aqueous feed's properties.

The user must specify the organic feed input in the ``organic_stream`` configuration, with a
configuration that describes the organic feed's properties.

The number of stages in the solvent extraction process has to be specified by the user through
the ``number_of_finite_elements`` configuration. It takes an integer value.


Stream configurations
---------------------

Each of the feed streams has to have a dictionary that specifies the property packages and other
details as mentioned below.

The ``property_package`` configuration is the property package that describes the state conditions
and properties of a particular stream.

The ``property_package_args`` configuration is any specific set of arguments that has to be passed
to the property block for the unit operation.

The user can specify the direction of the flow of the stream through the stages through the
configuration ``flow_direction``. This is a configuration, that uses FlowDirection Enum, which
can have two possible values.

Degrees of freedom
------------------

When the solvent extraction model is operated in steady state, the number of degrees of freedom of
the model is equal to the sum of the number of distribution coefficients of the total components
involved in the mass transfer operation and the volumes and volume fractions, for all the stages.

If the model is operated in dynamic state, the number of degrees of freedom is equal to the sum
of the distribution coefficient of all components involved in the mass transfer operation, values
of the state block variables of all the components of the system at the start of the operation, the
volumes and the volume fractions, for all the stages.

Model structure
---------------

The core model consists of a MSContactor model, with stream names hard coded as 'aqueous' and
'organic', and the stream dictionaries and number of finite elements are the same as those provided
by the user.

This model uses the heterogeneous reaction term defined in the MSContactor to calculate the amount of
material transferred between the phases for each of the rare earth elements. The distribution coefficients
used for this quantification are defined in the reaction package, and the constraint pertaining to the
distribution coefficient is defined in the solvent extraction model.

The pressure buildup in each of the stages has been defined in the model. For defining the pressure, we
need the volume of the phases, so the configuration ``has_holdup`` has to be set to True to obtain the
pressure of the phases.

"""

from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Param, Block, units
from pyomo.network import Port

from idaes.core import (
    FlowDirection,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.initialization import ModularInitializerBase

from idaes.models.unit_models.mscontactor import MSContactor


class SolventExtractionInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer  for the Solvent Extraction unit model.

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

        model.mscontactor.heterogeneous_reaction_extent.fix(1e-8)

        model.mscontactor.volume.fix()
        model.mscontactor.volume_frac_stream[:, :, "aqueous"].fix()

        # Initialize MSContactor
        msc_init = model.mscontactor.default_initializer(
            ssc_solver_options=self.config.ssc_solver_options,
            calculate_variable_options=self.config.calculate_variable_options,
        )
        msc_init.initialize(model.mscontactor)

        model.mscontactor.heterogeneous_reaction_extent.unfix()

        solver = self._get_solver()
        results = solver.solve(model)

        return results


Stream_Config = ConfigDict()

Stream_Config.declare(
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

Stream_Config.declare(
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

Stream_Config.declare(
    "flow_direction",
    ConfigValue(
        default=FlowDirection.forward,
        domain=In(FlowDirection),
        doc="Direction of flow for stream",
        description="FlowDirection Enum indicating direction of "
        "flow for given stream. Default=FlowDirection.forward.",
    ),
)

Stream_Config.declare(
    "has_energy_balance",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include energy balance for stream. Default=False.",
    ),
)

Stream_Config.declare(
    "has_pressure_balance",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include pressure balance for stream. Default=False.",
    ),
)


@declare_process_block_class("SolventExtraction")
class SolventExtractionData(UnitModelBlockData):

    default_initializer = SolventExtractionInitializer

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "aqueous_stream",
        Stream_Config(
            description="Aqueous stream properties",
        ),
    )

    CONFIG.declare(
        "organic_stream",
        Stream_Config(
            description="Organic stream properties",
        ),
    )

    CONFIG.declare(
        "number_of_finite_elements",
        ConfigValue(domain=int, description="Number of finite elements to use"),
    )

    CONFIG.declare(
        "heterogeneous_reaction_package",
        ConfigValue(
            description="Heterogeneous reaction package for solvent extraction.",
        ),
    )
    CONFIG.declare(
        "heterogeneous_reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for heterogeneous reaction package for solvent extractiong.",
        ),
    )

    def build(self):
        super().build()

        streams_dict = {
            "aqueous": self.config.aqueous_stream,
            "organic": self.config.organic_stream,
        }
        self.mscontactor = MSContactor(
            streams=streams_dict,
            number_of_finite_elements=self.config.number_of_finite_elements,
            heterogeneous_reactions=self.config.heterogeneous_reaction_package,
            heterogeneous_reactions_args=self.config.heterogeneous_reaction_package_args,
            has_holdup=self.config.has_holdup,
        )

        def distribution_ratio_rule(b, t, s, e):
            return (
                b.mscontactor.organic[t, s].conc_mol_comp[f"{e}_o"]
                == b.mscontactor.heterogeneous_reactions[t, s].distribution_coefficient[
                    e
                ]
                * b.mscontactor.aqueous[t, s].conc_mol_comp[e]
            )

        self.distribution_extent_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.heterogeneous_reaction_package.element_list,
            rule=distribution_ratio_rule,
        )

        self.area_cross_stage = Param(
            self.mscontactor.elements,
            units=units.m**2,
            doc="Cross sectional area stage",
            initialize=1,
            mutable=True,
        )

        self.elevation = Param(
            self.mscontactor.elements,
            units=units.m,
            doc="Height of settler tank base above outflow valve level",
            initialize=1,
            mutable=True,
        )

        def organic_pressure_calculation(b, t, s):

            g = Constants.acceleration_gravity
            P_atm = 101325 * units.Pa

            rho_og = b.config.organic_stream["property_package"].dens_mass

            P_org = units.convert(
                (
                    rho_og
                    * g
                    * b.mscontactor.volume[s]
                    * b.mscontactor.volume_frac_stream[t, s, "organic"]
                    / b.area_cross_stage[s]
                ),
                to_units=units.Pa,
            )

            return b.mscontactor.organic[t, s].pressure == P_org + P_atm

        self.organic_pressure_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            rule=organic_pressure_calculation,
        )

        def aqueous_pressure_calculation(b, t, s):
            g = Constants.acceleration_gravity

            rho_aq = b.config.aqueous_stream["property_package"].dens_mass
            P_aq = units.convert(
                (
                    rho_aq
                    * g
                    * (
                        b.mscontactor.volume[s]
                        * b.mscontactor.volume_frac_stream[t, s, "aqueous"]
                        / b.area_cross_stage[s]
                        + b.elevation[s]
                    )
                ),
                to_units=units.Pa,
            )

            return (
                b.mscontactor.aqueous[t, s].pressure
                == P_aq + b.mscontactor.organic[t, s].pressure
            )

        self.aqueous_pressure_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            rule=aqueous_pressure_calculation,
        )

        self.aqueous_inlet = Port(extends=self.mscontactor.aqueous_inlet)
        self.aqueous_outlet = Port(extends=self.mscontactor.aqueous_outlet)
        self.organic_inlet = Port(extends=self.mscontactor.organic_inlet)
        self.organic_outlet = Port(extends=self.mscontactor.organic_outlet)
