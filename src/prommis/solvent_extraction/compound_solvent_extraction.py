from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import (
    Constraint,
    Param,
    Block,
    units,
    RangeSet,
    TransformationFactory,
)
from pyomo.network import Port, Arc

from idaes.core import (
    FlowDirection,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
    ControlVolume1DBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.initialization import ModularInitializerBase

from prommis.solvent_extraction.solvent_extraction import SolventExtraction


# class CompoundSolventExtractionInitializer(ModularInitializerBase):

#     def initialize_main_model(
#         self,
#         model: Block,
#     ):
#         for i in model.elements:
#             model.solvent_extraction_mixer[i]


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


@declare_process_block_class("CompoundSolventExtraction")
class CompoundSolventExtractionData(UnitModelBlockData):

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
            description="Arguments for heterogeneous reaction package for solvent extraction.",
        ),
    )

    CONFIG.declare(
        "settler_finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
            domain (default=20)""",
        ),
    )

    CONFIG.declare(
        "settler_collocation_points",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=3)""",
        ),
    )

    CONFIG.declare(
        "settler_transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See
        Pyomo documentation for supported transformations.""",
        ),
    )

    CONFIG.declare(
        "settler_transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transformating domain. See
        Pyomo documentation for supported schemes.""",
        ),
    )

    def build(self):
        super().build()

        self.elements = RangeSet(
            1,
            self.config.number_of_finite_elements,
            doc="Set of finite elements in cascade (1 to number of elements)",
        )

        self.solvent_extraction_mixer = Block(self.elements)
        self.solvent_extraction_aqueous_settler = Block(self.elements)
        self.solvent_extraction_organic_settler = Block(self.elements)

        # Declare the mixer and settler tanks
        for i in self.elements:
            self.solvent_extraction_mixer[i].unit = SolventExtraction(
                number_of_finite_elements=1,
                aqueous_stream=self.config.aqueous_stream,
                organic_stream=self.config.organic_stream,
                heterogeneous_reaction_package=self.config.heterogeneous_reaction_package,
                heterogeneous_reaction_package_args=self.config.heterogeneous_reaction_package_args,
                has_holdup=self.config.has_holdup,
                dynamic=self.config.dynamic,
            )

        # Declare aqueous settler tanks
        for i in self.elements:
            self.solvent_extraction_aqueous_settler[i].unit = ControlVolume1DBlock(
                dynamic=self.config.dynamic,
                has_holdup=self.config.has_holdup,
                property_package=self.config.aqueous_stream.property_package,
                property_package_args=self.config.aqueous_stream.property_package_args,
                transformation_method=self.config.settler_transformation_method,
                transformation_scheme=self.config.settler_transformation_scheme,
                finite_elements=self.config.settler_finite_elements,
                collocation_points=self.config.settler_collocation_points,
            )
            self.solvent_extraction_aqueous_settler[i].unit.add_geometry(
                flow_direction=FlowDirection.forward,
            )
            self.solvent_extraction_aqueous_settler[i].unit.add_state_blocks(
                information_flow=FlowDirection.forward,
                has_phase_equilibrium=False,
            )
            self.solvent_extraction_aqueous_settler[i].unit.add_material_balances(
                balance_type=MaterialBalanceType.componentTotal,
                has_phase_equilibrium=False,
                has_mass_transfer=False,
            )
            self.solvent_extraction_aqueous_settler[i].unit.add_energy_balances(
                balance_type=EnergyBalanceType.none,
            )
            self.solvent_extraction_aqueous_settler[i].unit.add_momentum_balances(
                balance_type=MomentumBalanceType.none,
            )
            self.solvent_extraction_aqueous_settler[i].unit.apply_transformation()
            self.add_inlet_port(
                name=f"aqueous_settler_{i}_inlet",
                block=self.solvent_extraction_aqueous_settler[i].unit,
            )
            self.add_outlet_port(
                name=f"aqueous_settler_{i}_outlet",
                block=self.solvent_extraction_aqueous_settler[i].unit,
            )

        # Declare organic settler tanks
        for i in self.elements:

            self.solvent_extraction_organic_settler[i].unit = ControlVolume1DBlock(
                dynamic=self.config.dynamic,
                has_holdup=self.config.has_holdup,
                property_package=self.config.organic_stream.property_package,
                property_package_args=self.config.organic_stream.property_package_args,
                transformation_method=self.config.settler_transformation_method,
                transformation_scheme=self.config.settler_transformation_scheme,
                finite_elements=self.config.settler_finite_elements,
                collocation_points=self.config.settler_collocation_points,
            )
            self.solvent_extraction_organic_settler[i].unit.add_geometry(
                flow_direction=FlowDirection.forward,
            )
            self.solvent_extraction_organic_settler[i].unit.add_state_blocks(
                information_flow=FlowDirection.forward,
                has_phase_equilibrium=False,
            )
            self.solvent_extraction_organic_settler[i].unit.add_material_balances(
                balance_type=MaterialBalanceType.componentTotal,
                has_phase_equilibrium=False,
                has_mass_transfer=False,
            )
            self.solvent_extraction_organic_settler[i].unit.add_energy_balances(
                balance_type=EnergyBalanceType.none,
            )
            self.solvent_extraction_organic_settler[i].unit.add_momentum_balances(
                balance_type=MomentumBalanceType.none,
            )
            self.solvent_extraction_organic_settler[i].unit.apply_transformation()
            self.add_inlet_port(
                name=f"organic_settler_{i}_inlet",
                block=self.solvent_extraction_organic_settler[i].unit,
            )
            self.add_outlet_port(
                name=f"organic_settler_{i}_outlet",
                block=self.solvent_extraction_organic_settler[i].unit,
            )

        # Declare the mixer settler arcs
        for i in self.elements:
            for j in ["aqueous", "organic"]:
                setattr(
                    self,
                    f"solvent_extraction_mixer_{j}_settler_arc_{i}",
                    Arc(
                        source=getattr(
                            self.solvent_extraction_mixer[i].unit, f"{j}_outlet"
                        ),
                        destination=getattr(self, f"{j}_settler_{i}_inlet"),
                    ),
                )

        # Declare directional arcs between mixer-settler tanks
        for i in self.elements:
            for j in ["aqueous", "organic"]:
                if (
                    getattr(self.config, f"{j}_stream").flow_direction
                    == FlowDirection.forward
                ):
                    if i != self.elements.first():
                        previous_stage = i - 1
                        current_stage = i
                    else:
                        continue
                else:
                    if i != self.elements.last():
                        previous_stage = i + 1
                        current_stage = i
                    else:
                        continue
                setattr(
                    self,
                    f"solvent_extraction_{j}_settler_{previous_stage}_to_mixer_{current_stage}_arc",
                    Arc(
                        source=getattr(
                            self,
                            f"{j}_settler_{previous_stage}_outlet",
                        ),
                        destination=getattr(
                            self.solvent_extraction_mixer[current_stage].unit,
                            f"{j}_inlet",
                        ),
                    ),
                )

        TransformationFactory("network.expand_arcs").apply_to(self)

        # define ports
        for j in ["aqueous", "organic"]:
            if (
                getattr(self.config, f"{j}_stream").flow_direction
                == FlowDirection.forward
            ):
                inlet_stage = self.elements.first()
                outlet_stage = self.elements.last()
            else:
                inlet_stage = self.elements.last()
                outlet_stage = self.elements.first()
            setattr(
                self,
                f"{j}_inlet",
                Port(
                    extends=getattr(
                        self.solvent_extraction_mixer[inlet_stage].unit,
                        f"{j}_inlet",
                    ),
                ),
            )
            setattr(
                self,
                f"{j}_outlet",
                Port(
                    extends=getattr(
                        self,
                        f"{j}_settler_{outlet_stage}_outlet",
                    ),
                ),
            )
