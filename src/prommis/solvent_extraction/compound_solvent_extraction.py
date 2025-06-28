from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import (
    Block,
    RangeSet,
    TransformationFactory,
    value,
    Var,
)
from pyomo.network import Port, Arc
from pyomo.dae.flatten import flatten_dae_components
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
from idaes.core.initialization import (
    ModularInitializerBase,
)
from idaes.core.util.initialization import propagate_state

from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
)


class CompoundSolventExtractionInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer  for the compound solvent extraction unit model.

    This routine calls the initializer for the internal SolventExtraction model.

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

        # Deconstructing the model and initializing the mixer and settler tanks separately

        for e in model.elements:
            for j in ["aqueous", "organic"]:

                getattr(model, f"mixer_{j}_settler_arc_{e}_expanded").deactivate()
                if (
                    getattr(model.config, f"{j}_stream").flow_direction
                    == FlowDirection.backward
                ):
                    if e != model.elements.last():
                        getattr(
                            model,
                            f"{j}_settler_{model.elements.next(e)}_to_mixer_{e}_arc_expanded",
                        ).deactivate()
                        for k, v in getattr(model, f"{j}_inlet").vars.items():
                            for i in v:
                                for k1, v1 in getattr(
                                    model.mixer[e].unit, f"{j}_inlet"
                                ).vars.items():
                                    if k1 == k:
                                        v1[i].fix(value(v[i]))
                else:
                    if e != model.elements.first():
                        getattr(
                            model,
                            f"{j}_settler_{model.elements.prev(e)}_to_mixer_{e}_arc_expanded",
                        ).deactivate()
                        setattr(
                            model,
                            f"{j}_{model.elements.prev(e)}_bypass_{e}_arc",
                            Arc(
                                source=model.mixer[e - 1].unit.aqueous_outlet,
                                destination=model.mixer[e].unit.aqueous_inlet,
                            ),
                        ),

            model.aqueous_settler[e].unit.deactivate()
            model.organic_settler[e].unit.deactivate()

        TransformationFactory("network.expand_arcs").apply_to(model)

        sx_initializer = model.mixer[model.elements.first()].unit.default_initializer(
            ssc_solver_options=self.config.ssc_solver_options,
            calculate_variable_options=self.config.calculate_variable_options,
        )
        for e in model.elements:
            sx_initializer.initialize(model.mixer[e].unit)
            for j in ["aqueous", "organic"]:
                if (
                    getattr(model.config, f"{j}_stream").flow_direction
                    == FlowDirection.forward
                ):
                    if e != model.elements.last():
                        propagate_state(
                            arc=getattr(
                                model, f"{j}_{e}_bypass_{model.elements.next(e)}_arc"
                            )
                        )

        solver = self._get_solver()
        results_mx = solver.solve(model)

        for e in model.elements:

            model.aqueous_settler[e].unit.activate()
            model.organic_settler[e].unit.activate()

            for j in ["aqueous", "organic"]:
                for k, v in getattr(model, f"{j}_inlet").vars.items():
                    for i in v:
                        for k1, v1 in getattr(
                            model, f"{j}_settler_{e}_in"
                        ).vars.items():
                            if k1 == k:
                                v1[i].fix(value(v[i]))

                target_model = getattr(model, f"{j}_settler")[e].unit
                target_x = getattr(model, f"{j}_settler")[e].unit.length_domain

                regular_vars, length_vars = flatten_dae_components(
                    target_model,
                    target_x,
                    Var,
                    active=True,
                )
                for var in length_vars:
                    for x in target_x:
                        if x == target_x.first():
                            continue
                        else:
                            var[x].value = var[target_x.first()].value

                result_sx = solver.solve(target_model)

        # Reconstructing the model and solving

        for e in model.elements:
            for j in ["aqueous", "organic"]:

                getattr(model, f"{j}_settler_{e}_in").unfix()
                getattr(model, f"mixer_{j}_settler_arc_{e}_expanded").activate()

                if (
                    getattr(model.config, f"{j}_stream").flow_direction
                    == FlowDirection.backward
                ):
                    if e != model.elements.last():
                        getattr(model.mixer[e].unit, f"{j}_inlet").unfix()
                        getattr(
                            model,
                            f"{j}_settler_{model.elements.next(e)}_to_mixer_{e}_arc_expanded",
                        ).activate()
                else:
                    if e != model.elements.first():
                        getattr(
                            model,
                            f"{j}_settler_{model.elements.prev(e)}_to_mixer_{e}_arc_expanded",
                        ).activate()
                        model.del_component(
                            f"{j}_{model.elements.prev(e)}_bypass_{e}_arc_expanded"
                        )
                        model.del_component(
                            f"{j}_{model.elements.prev(e)}_bypass_{e}_arc"
                        )

        TransformationFactory("network.expand_arcs").apply_to(model)

        final_results = solver.solve(model)

        return final_results


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

    default_initializer = CompoundSolventExtractionInitializer

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

        self.mixer = Block(self.elements)
        self.aqueous_settler = Block(self.elements)
        self.organic_settler = Block(self.elements)

        for i in self.elements:

            # Declare the mixer tank
            self.mixer[i].unit = SolventExtraction(
                number_of_finite_elements=1,
                aqueous_stream=self.config.aqueous_stream,
                organic_stream=self.config.organic_stream,
                heterogeneous_reaction_package=self.config.heterogeneous_reaction_package,
                heterogeneous_reaction_package_args=self.config.heterogeneous_reaction_package_args,
                has_holdup=self.config.has_holdup,
                dynamic=self.config.dynamic,
            )

            # Declare aqueous settler tank
            self.aqueous_settler[i].unit = ControlVolume1DBlock(
                dynamic=self.config.dynamic,
                has_holdup=self.config.has_holdup,
                property_package=self.config.aqueous_stream.property_package,
                property_package_args=self.config.aqueous_stream.property_package_args,
                transformation_method=self.config.settler_transformation_method,
                transformation_scheme=self.config.settler_transformation_scheme,
                finite_elements=self.config.settler_finite_elements,
                collocation_points=self.config.settler_collocation_points,
            )
            self.aqueous_settler[i].unit.add_geometry(
                flow_direction=FlowDirection.forward,
            )
            self.aqueous_settler[i].unit.add_state_blocks(
                information_flow=FlowDirection.forward,
                has_phase_equilibrium=False,
            )
            self.aqueous_settler[i].unit.add_material_balances(
                balance_type=MaterialBalanceType.componentTotal,
                has_phase_equilibrium=False,
                has_mass_transfer=False,
            )
            self.aqueous_settler[i].unit.add_energy_balances(
                balance_type=EnergyBalanceType.none,
            )
            self.aqueous_settler[i].unit.add_momentum_balances(
                balance_type=MomentumBalanceType.none,
            )
            self.aqueous_settler[i].unit.apply_transformation()
            self.add_inlet_port(
                name=f"aqueous_settler_{i}_in",
                block=self.aqueous_settler[i].unit,
            )
            self.add_outlet_port(
                name=f"aqueous_settler_{i}_out",
                block=self.aqueous_settler[i].unit,
            )

            # Declare organic settler tanks
            self.organic_settler[i].unit = ControlVolume1DBlock(
                dynamic=self.config.dynamic,
                has_holdup=self.config.has_holdup,
                property_package=self.config.organic_stream.property_package,
                property_package_args=self.config.organic_stream.property_package_args,
                transformation_method=self.config.settler_transformation_method,
                transformation_scheme=self.config.settler_transformation_scheme,
                finite_elements=self.config.settler_finite_elements,
                collocation_points=self.config.settler_collocation_points,
            )
            self.organic_settler[i].unit.add_geometry(
                flow_direction=FlowDirection.forward,
            )
            self.organic_settler[i].unit.add_state_blocks(
                information_flow=FlowDirection.forward,
                has_phase_equilibrium=False,
            )
            self.organic_settler[i].unit.add_material_balances(
                balance_type=MaterialBalanceType.componentTotal,
                has_phase_equilibrium=False,
                has_mass_transfer=False,
            )
            self.organic_settler[i].unit.add_energy_balances(
                balance_type=EnergyBalanceType.none,
            )
            self.organic_settler[i].unit.add_momentum_balances(
                balance_type=MomentumBalanceType.none,
            )
            self.organic_settler[i].unit.apply_transformation()
            self.add_inlet_port(
                name=f"organic_settler_{i}_in",
                block=self.organic_settler[i].unit,
            )
            self.add_outlet_port(
                name=f"organic_settler_{i}_out",
                block=self.organic_settler[i].unit,
            )

        for i in self.elements:
            for j in ["aqueous", "organic"]:

                setattr(
                    self,
                    f"mixer_{j}_settler_arc_{i}",
                    Arc(
                        source=getattr(self.mixer[i].unit, f"{j}_outlet"),
                        destination=getattr(self, f"{j}_settler_{i}_in"),
                    ),
                )

                # Declare directional arcs between mixer-settler tanks
                if (
                    getattr(self.config, f"{j}_stream").flow_direction
                    == FlowDirection.forward
                ):
                    if i != self.elements.first():
                        setattr(
                            self,
                            f"{j}_settler_{self.elements.prev(i)}_to_mixer_{i}_arc",
                            Arc(
                                source=getattr(
                                    self,
                                    f"{j}_settler_{self.elements.prev(i)}_out",
                                ),
                                destination=getattr(
                                    self.mixer[i].unit,
                                    f"{j}_inlet",
                                ),
                            ),
                        )

                else:
                    if i != self.elements.last():
                        setattr(
                            self,
                            f"{j}_settler_{self.elements.next(i)}_to_mixer_{i}_arc",
                            Arc(
                                source=getattr(
                                    self,
                                    f"{j}_settler_{self.elements.next(i)}_out",
                                ),
                                destination=getattr(
                                    self.mixer[i].unit,
                                    f"{j}_inlet",
                                ),
                            ),
                        )

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
                        self.mixer[inlet_stage].unit,
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
                        f"{j}_settler_{outlet_stage}_out",
                    ),
                ),
            )

        TransformationFactory("network.expand_arcs").apply_to(self)
