#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
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

The material transfer can happen from either of the phases to the other. To specify the direction 
of the transfer, the ``aqueous_to_organic`` configuration is to be used by the user. This is a boolean 
configuration. The default value is True, which means the material transfer is happening from the 
aqueous phase to the organic phase, like in the loading operation. For scrubbing and stripping, the 
reverse happens, so the value of the configuration will be False.

Stream configurations
---------------------

Each of the feed streams has to have a dictionary that specifies the property packages and other 
details as mentioned below.

The ``property_package`` configuration is the property package that describes the state conditions 
and properties of a particular stream.

The ``property_package_args`` configuration is any specific set of arguments that has to be passed to 
the property block for the unit operation.

The user can specify the direction of the flow of the stream through the stages through the 
configuration ``flow_direction``. This is a configuration, that uses FlowDirection Enum, which
can have two possible values.

Degrees of freedom 
------------------

When the solvent extraction model is operated in steady state, the number of degrees of freedom of
the model is equal to the number of partition coefficients of the total components involved in the
mass transfer operation, for all the stages.

If the model is operated in dynamic state, the number of degrees of freedom is equal to the sum
of the partition coefficient of all components involved in the mass transfer operation, 
values of the state block variables of all the components of the system at the start of the
operation, the volumes and the volume fractions, for all the stages.

Model structure
---------------

The core model consists of a MSContactor model, with stream names hard coded as 'aqueous' and 
'organic', and the stream dictionaries and number of finite elements are the same as those provided
by the user.

This model defines the material transfer term defined in the MSContactor and expresses it as a
function of the parameter of partition coefficient defined by the user.
"""

from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Param, Block
from pyomo.network import Port

from idaes.core import (
    FlowDirection,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.initialization import ModularInitializerBase
from idaes.models.unit_models.mscontactor import MSContactor


class SolventExtractionInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer  for the Solvent Extraction unit model.

    This routine calls the initializer for the internal MSContactor model.

    """

    CONFIG = ModularInitializerBase.CONFIG()

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
        """
        This model adds an additional constraint of the material transfer term.
        This constraint is present outside the main MSContactor model, so this term is 
        fixed to give a square model to the MScontactor Initializer.

        """

        model.mscontactor.material_transfer_term.fix(1e-8)

        msc_init = self.get_submodel_initializer(model.mscontactor)
        msc_init.initialize(model.mscontactor)

        model.mscontactor.material_transfer_term.unfix()

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
        "aqueous_to_organic",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Direction of the transfer between two phases",
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
        )

        self.partition_coefficient = Param(
            self.mscontactor.elements,
            self.mscontactor.stream_component_interactions,
            initialize=1,
            mutable=True,
            doc="The fraction of component that goes from aqueous to organic phase",
        )

        self.aqueous_inlet = Port(extends=self.mscontactor.aqueous_inlet)
        self.aqueous_outlet = Port(extends=self.mscontactor.aqueous_outlet)
        self.organic_inlet = Port(extends=self.mscontactor.organic_inlet)
        self.organic_outlet = Port(extends=self.mscontactor.organic_outlet)

        def mass_transfer_term(b, t, s, k, l, m):
            if self.config.aqueous_to_organic:
                stream_state = b.mscontactor.aqueous
                in_state = b.mscontactor.aqueous_inlet_state
                stream_name = self.config.aqueous_stream
                sign = -1
            else:
                stream_state = b.mscontactor.organic
                in_state = b.mscontactor.organic_inlet_state
                stream_name = self.config.organic_stream
                sign = 1

            if stream_name.flow_direction == FlowDirection.forward:
                if s == b.mscontactor.elements.first():
                    state = in_state[t]
                else:
                    state = stream_state[t, b.mscontactor.elements.prev(s)]
            else:
                if s == b.mscontactor.elements.last():
                    state = in_state[t]
                else:
                    state = stream_state[t, b.mscontactor.elements.next(s)]

            return (
                b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                == sign
                * state.get_material_flow_terms(stream_state.phase_list, m)
                * b.partition_coefficient[s, (k, l, m)]
            )

        self.mass_transfer_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.mscontactor.stream_component_interactions,
            rule=mass_transfer_term,
        )
