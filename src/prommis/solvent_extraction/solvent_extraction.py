from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Param
from pyomo.network import Port

from idaes.core import (
    FlowDirection,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.unit_models.mscontactor import MSContactor

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

        def param_init(b, s, k, l, m):
            b.partition_coefficient[s, (k, l, m)] = 1

        self.partition_coefficient = Param(
            self.mscontactor.elements,
            self.mscontactor.stream_component_interactions,
            initialize=param_init,
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
