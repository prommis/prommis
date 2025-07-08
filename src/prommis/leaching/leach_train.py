#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Leach Train
===========

Author: Andrew Lee

The Leach Train unit model represents a series of leaching tanks in which a solid and liquid stream are contacted and
undergo heterogeneous chemical reactions.

Configuration Arguments
-----------------------

When creating an instance of a Leach Train model, the user must specify the number of tanks in train using the
``number_of_tanks`` configuration argument. This value cannot be changed later.

A Leach Train model also requires a special "heterogeneous reaction package" to define the reactions which occur in the
system and the constraints defining the rate of reaction or equilibrium condition associated with each. These custom
reaction packages are similar to normal IDAES reaction packages with the following differences:

* the reaction package is not coupled to any specific thermophysical property package (i.e., the liquid or solid phase
  properties). Instead, it is expected that the reaction package will link to both the solid and liquid phase properties
  using ``self.parent_block().liquid`` and ``self.parent_block().solid`` as required.

Degrees of Freedom
------------------

A Leach Train module has a number of degrees of freedom equal to the number of tanks in the train. These are the volume
of each tank, and is specified using the ``volume`` attribute on the unit model (indexed by tank number).

Model Structure
---------------

The core Leach Train unit model consists of a single ``MSContactor`` model (named ``mscontactor``) with hard coded
stream names (``liquid`` and ``solid`` respectively). The Leach Train model also has two inlets and two outlets named
``liquid_inlet``, ``solid_inlet``, ``liquid_outlet`` and ``solid_outlet`` respectively.

Additional Constraints
----------------------

Leach Train units add one additional constraint to define the extent of reaction for rate reactions in the system
(indexed by number of tanks and list of rate reactions).

.. math:: X_{t,e,r} = V_{t,e} \times r_{t,e,r}

where :math:`X_{t,e,r}` is the extent of reaction for reaction :math:`r` in tank :math:`e` at time :math:`t`,
:math:`V_{t,e}` is the volume of the reacting material in tank :math:`e` at time :math:`t` (allows for varying reactor
volume with time) and :math:`r_{t,e,r}` is the volumetric rate of reaction for reaction :math:`r` in tank :math:`e` at
time :math:`t` (from the reaction package).

Variables
---------

Leach Train units add the following additional Variables beyond those created by the MSContactor Block.

================ ====== ================================================================================================
Variable         Name   Notes
================ ====== ================================================================================================
:math:`V_{t,e}`  volume
================ ====== ================================================================================================

"""

from pyomo.common.config import Bool, ConfigDict, ConfigValue
from pyomo.environ import Block, Constraint, Var, units
from pyomo.network import Port

from idaes.core import (
    MaterialFlowBasis,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.unit_models.mscontactor import MSContactor


class LeachingTrainInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer  for the Leaching Train unit model.

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
        msc_init = model.mscontactor.default_initializer(
            ssc_solver_options=self.config.ssc_solver_options,
            calculate_variable_options=self.config.calculate_variable_options,
        )
        msc_init.initialize(model.mscontactor)

        solver = self._get_solver()
        results = solver.solve(model)

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


@declare_process_block_class("LeachingTrain")
class LeachingTrainData(UnitModelBlockData):
    """
    Leaching Train Unit Model Class
    """

    # Set default initializer
    default_initializer = LeachingTrainInitializer

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
            description="Heterogeneous reaction package for leaching.",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for heterogeneous reaction package for leaching.",
        ),
    )
    CONFIG.declare(
        "number_of_tanks",
        ConfigValue(
            default=1, domain=int, description="Number of tanks in leaching train"
        ),
    )

    def build(self):
        """
        Build method for LeachingTrain unit model.
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
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
        )

        # Get units of measurement from MSContactor
        flow_basis = self.mscontactor.flow_basis
        uom = self.mscontactor.uom

        # Reactor volume
        self.volume = Var(
            self.flowsheet().time,
            self.mscontactor.elements,
            initialize=1,
            units=uom.VOLUME,
            doc="Volume of each tank.",
        )

        # Note that this is being added to the MSContactor block
        def rule_heterogeneous_reaction_extent(b, t, s, r):
            volume = b.parent_block().volume

            if flow_basis == MaterialFlowBasis.mass:
                m_units = uom.MASS
                x_units = m_units / uom.TIME
            elif flow_basis == MaterialFlowBasis.molar:
                m_units = uom.AMOUNT
                x_units = m_units / uom.TIME
            else:
                # Undefined
                x_units = None

            return units.convert(
                b.heterogeneous_reaction_extent[t, s, r], to_units=x_units
            ) == units.convert(
                b.heterogeneous_reactions[t, s].reaction_rate[r] * volume[t, s],
                to_units=x_units,
            )

        self.mscontactor.heterogeneous_reaction_extent_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.reaction_package.reaction_idx,
            rule=rule_heterogeneous_reaction_extent,
        )

        # Create unit level Ports
        self.liquid_inlet = Port(extends=self.mscontactor.liquid_inlet)
        self.liquid_outlet = Port(extends=self.mscontactor.liquid_outlet)
        self.solid_inlet = Port(extends=self.mscontactor.solid_inlet)
        self.solid_outlet = Port(extends=self.mscontactor.solid_outlet)

        # Create Expression for element recovery
        @self.Expression(
            self.flowsheet().time,
            self.mscontactor.solid.component_list,
            doc="Percent recovery of solid species",
        )
        def recovery(b, t, j):
            f_in = b.solid_inlet.flow_mass[t]
            f_out = b.solid_outlet.flow_mass[t]
            x_in = b.solid_inlet.mass_frac_comp[t, j]
            x_out = b.solid_outlet.mass_frac_comp[t, j]

            return (1 - f_out * x_out / (f_in * x_in)) * 100

    def _get_stream_table_contents(self, time_point=0):
        return self.mscontactor._get_stream_table_contents(time_point)

    def _get_performance_contents(self, time_point=0):
        exprs = {}
        for j in self.mscontactor.solid.component_list:
            exprs[f"Recovery % {j}"] = self.recovery[time_point, j]

        return {"exprs": exprs}
