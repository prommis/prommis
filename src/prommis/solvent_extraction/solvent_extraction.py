from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Param, Var, units
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
        doc="Bool indicating whether to include energy balance for stream. Default=True.",
    ),
)

Stream_Config.declare(
    "has_pressure_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include pressure balance for stream. Default=True.",
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
            mutable=True,
            doc="The fraction of component that goes from aqueous to organic phase",
        )

        self.aqueous_inlet = Port(extends=self.mscontactor.aqueous_inlet)
        self.aqueous_outlet = Port(extends=self.mscontactor.aqueous_outlet)
        self.organic_inlet = Port(extends=self.mscontactor.organic_inlet)
        self.organic_outlet = Port(extends=self.mscontactor.organic_outlet)

        def mass_transfer_term(b, t, s, k, l, m):
            if self.config.aqueous_to_organic == True:
                aqueous = b.mscontactor.aqueous
                if self.config.aqueous_stream.flow_direction == FlowDirection.forward:
                    if s == b.mscontactor.elements.first():
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == -b.mscontactor.aqueous_inlet_state[
                                t
                            ].get_material_flow_terms(aqueous.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                    else:
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == -b.mscontactor.aqueous[
                                t, b.mscontactor.elements.prev(s)
                            ].get_material_flow_terms(aqueous.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                else:
                    if s == b.mscontactor.elements.last():
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == -b.mscontactor.aqueous_inlet_state[
                                t
                            ].get_material_flow_terms(aqueous.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                    else:
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == -b.mscontactor.aqueous[
                                t, b.mscontactor.elements.next(s)
                            ].get_material_flow_terms(aqueous.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
            else:
                organic = b.mscontactor.organic
                if self.config.organic_stream.flow_direction == FlowDirection.forward:
                    if s == b.mscontactor.elements.first():
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == b.mscontactor.organic_inlet_state[
                                t
                            ].get_material_flow_terms(organic.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                    else:
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == b.mscontactor.organic[
                                t, b.mscontactor.elements.prev(s)
                            ].get_material_flow_terms(organic.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                else:
                    if s == b.mscontactor.elements.last():
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == b.mscontactor.organic_inlet_state[
                                t
                            ].get_material_flow_terms(organic.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )
                    else:
                        return (
                            b.mscontactor.material_transfer_term[t, s, (k, l, m)]
                            == b.mscontactor.organic[
                                t, b.mscontactor.elements.next(s)
                            ].get_material_flow_terms(organic.phase_list, m)
                            * b.partition_coefficient[s, (k, l, m)]
                        )

        self.mass_transfer_constraint = Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.mscontactor.stream_component_interactions,
            rule=mass_transfer_term,
        )

        self.clean_sx_pe_tank_cap = Var(
            initialize=1405, units=units.gal, doc="Volume of polyethylene tanks"
        )
        self.num_sx_pe_tank = Var(initialize=5, doc="Number of polyethylene tanks")

        self.clean_sx_tank_mixer_pow = Var(
            initialize=1, units=units.hp, doc="Power of tank mixers"
        )
        self.num_sx_tank_mixer = Var(initialize=2, doc="Number of mixers")

        self.clean_sx_process_pump_feed = Var(
            initialize=281, units=units.gal / units.min, doc="Process pump feed rates"
        )
        self.num_sx_process_pump = Var(initialize=3, doc="Number of process pumps")

        self.clean_sx_mix_set_cap = Var(
            initialize=2444, units=units.gal, doc="Volume of mixer settler tanks"
        )
        self.num_sx_mix_set = Var(initialize=6, doc="Number of mixer settler tanks")