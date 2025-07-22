#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

r"""
Mixer Settler Extraction Model

========================

Author: Arkoprabho Dasgupta

The Mixer Settler Extraction unit is an extension of the Solvent Extraction unit model from
a pure mixer tank-based model to a mixer-settler based model. This model gives outlet concentration
and flow values identical to the Solvent Extraction model, but it takes time delays into consideration
which is important for dynamics.
The mixer models are the Solvent Extraction model with 1 tank, and the settler tanks are modeled
as an independent process flow reactor without reaction (an empty pipe) for each phase, using 1D Control Volume blocks.

Configuration Arguments
-----------------------

The user must specify the following configurations in a solvent extraction model to be able to
use it.

The user must specify the aqueous feed input in the ``aqueous_stream`` configuration, with a
dictionary containing the aqueous feed's properties.

The user must specify the organic feed input in the ``organic_stream`` configuration, with a
dictionary containing the organic feed's properties.

The number of stages in the solvent extraction process has to be specified by the user through
the ``number_of_stages`` configuration. It takes an integer value.

The user must give a heterogeneous reaction package in the ``heterogeneous_reaction_package``
argument for the extraction reaction between the two phases, and any additional arguments can
be given in the ``heterogeneous_reaction_package_args`` argument.

The user needs to specify the discretization method of the settler tank length in the argument
``settler_transformation_method``. Once the discretization method is specified, the user needs
to specify the corresponding discretization scheme in the argument ``settler_transformation_scheme``.

After the discretization, the user needs to specify the number of finite elements in the argument
``settler_finite_elements``. If the discretization method is collocation, then the user has to
additionally specify the number of collocation points in the ``settler_collocation_points``, else
the user does not need to specify anything in this argument.

Stream Configurations
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

The stream has two more arguments, ``has_energy_balance`` and ``has_pressure_balance`` in lieu with
the base MSContactor model. However, our model does not consider energy balance and pressure balance
yet, so the default arguments will be ``has_energy_balance = False`` and ``has_pressure_balance = False``.

Degrees of Freedom
------------------

When the model is operated in steady state, the number of degrees of freedom of each unit is the sum
of the volume of the mixer tank, and the length and area of both the settler tanks. The total degrees
of freedom is the degrees of freedom of one unit multiplied by the number of stages.

If the model is operated in dynamic state, the number of degrees of freedom of each unit is equal to
the sum of all the variables in the mixer and the two settler tanks whose values have to be fixed at
time=0, the volume of the mixer tank, and the area and length of the two settler tanks. The total degrees
of freedom is the degrees of freedom of each unit multiplied by the number of tanks.

Model Structure
---------------

The core model is a combination of the Solvent Extraction models and 1D Control Volume models.
Each unit consists of a mixer tank model, which is a Solvent Extraction unit model, accompanied by
two 1D Control Volume unit models for the aqueous and organic phases, respectively. This entire unit
is repeated as many times as the number of stages in the unit model.
The material transfer between the two phases happens only in the mixer tanks, ie. the Solvent Extraction
models, and they are quantified by the heterogeneous reaction package.
Each of the mixer and settler tanks are connected extensively by Arcs, named after the sequence and
the stage number of the tanks.

Additional Constraints
----------------------

In addition to the MSContactor model, the model declares the following constraints.

1. temperature_constraint = This constraint equates the inlet and outlet temperatures of all the
settler tanks. These constraints are added directly to the 1D Control Volume blocks.

2. pressure_constraint = This constraint equates the inlet and outlet pressures of all the
settler tanks. These constraints are added directly to the 1D Control Volume blocks.

"""

from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import Block, RangeSet, TransformationFactory, value, Var, Constraint
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

__author__ = "Arkoprabho Dasgupta"


class MixerSettlerExtractionInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer for the mixer settler extraction unit model.

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
        """
        This initializer model can be decomposed into the following steps.
        Step 1:
            Deactivation of the arcs connecting the mixer tank to their corresponding settler
            tanks, and the arcs connecting the settler tank to the next mixer tank.
        Step 2:
            If a stream is flowing backward, then the values of the state variables of that
            corresponding inlet stream to each of the mixer tank are fixed to the values of
            the state variables of that inlet stream to the overall model.
            If a stream is flowing forward, then a bypass arc connecting the corresponding
            outlet stream of the mixer to the inlet stream of the next mixer tank is created.
            The settler tanks are all deactivated, and the newly created arcs are transformed.
        Step 3:
            The default initializer model for SolventExtraction model is called and the mixer
            tanks are all initialized sequentially. While initializing the mixer tanks, if a
            stream is forward flowing, and the tank is not the first tank, then the values of
            the state variables of the stream outlet of the previous mixer is propagated along
            the bypass arcs to the inlet stream of the mixer, then the mixer is initialized.
            After the individual mixer tanks are initialized, the overall model is solved.
        Step 4:
            The settler tanks are all activated. The values of the state variables of the feed
            streams to all the settlers are fixed to the values of state variables of the corresponding
            inlet streams to the overall model, and the values at the inlet of the settler tanks
            are propagated throughout the settler tank length. Then the overall model is solved.
        Step 5:
            The state variables of the inlet streams to the settler tanks and the backward flowing
            inlet stream to the mixer tanks are all unfixed. All the bypass arcs are deleted and
            the mixer-settler arcs are all reactivated. Then the overall model is solved. The
            resulting solution is the initialized model.

        """

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
                        src_port = getattr(model, f"{j}_inlet")
                        dest_port = getattr(model.mixer[e].unit, f"{j}_inlet")
                        for var_name, src_var_obj in src_port.vars.items():
                            dest_var_obj = dest_port.vars[var_name]
                            for idx in src_var_obj:
                                dest_var_obj[idx].fix(value(src_var_obj[idx]))

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
                                source=model.mixer[
                                    model.elements.prev(e)
                                ].unit.aqueous_outlet,
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
                src_port = getattr(model, f"{j}_inlet")
                dest_port = getattr(model, f"{j}_settler_{e}_in")
                for var_name, src_var_obj in src_port.vars.items():
                    dest_var_obj = dest_port.vars[var_name]
                    for idx in src_var_obj:
                        dest_var_obj[idx].fix(value(src_var_obj[idx]))

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


@declare_process_block_class("MixerSettlerExtraction")
class MixerSettlerExtractionData(UnitModelBlockData):

    default_initializer = MixerSettlerExtractionInitializer

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
        "number_of_stages",
        ConfigValue(domain=int, description="Number of stages in the model"),
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

    def build(self):
        super().build()

        # Declare indexed blocks

        self.elements = RangeSet(
            1,
            self.config.number_of_stages,
            doc="Set of number of stages in cascade (1 to number of elements)",
        )

        # Declare indexed mixer tanks and settler tanks

        self.mixer = Block(self.elements)
        self.aqueous_settler = Block(self.elements)
        self.organic_settler = Block(self.elements)

        # Define the individual mixer and settler models

        for i in self.elements:

            # Declare the mixer tank
            self.mixer[i].unit = SolventExtraction(
                dynamic=self.config.dynamic,
                number_of_finite_elements=1,
                aqueous_stream=self.config.aqueous_stream,
                organic_stream=self.config.organic_stream,
                heterogeneous_reaction_package=self.config.heterogeneous_reaction_package,
                heterogeneous_reaction_package_args=self.config.heterogeneous_reaction_package_args,
                has_holdup=self.config.has_holdup,
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
            if self.config.aqueous_stream.has_energy_balance == True:
                aqueous_energy_balance_type = EnergyBalanceType.enthalpyTotal
            else:
                aqueous_energy_balance_type = EnergyBalanceType.none
            self.aqueous_settler[i].unit.add_energy_balances(
                balance_type=aqueous_energy_balance_type,
            )
            if self.config.aqueous_stream.has_pressure_balance == True:
                aqueous_pressure_balance_type = MomentumBalanceType.pressureTotal
            else:
                aqueous_pressure_balance_type = MomentumBalanceType.none
            self.aqueous_settler[i].unit.add_momentum_balances(
                balance_type=aqueous_pressure_balance_type,
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
            if self.config.organic_stream.has_energy_balance == True:
                organic_energy_balance_type = EnergyBalanceType.enthalpyTotal
            else:
                organic_energy_balance_type = EnergyBalanceType.none
            self.organic_settler[i].unit.add_energy_balances(
                balance_type=organic_energy_balance_type,
            )
            if self.config.organic_stream.has_pressure_balance == True:
                organic_pressure_balance_type = MomentumBalanceType.pressureTotal
            else:
                organic_pressure_balance_type = MomentumBalanceType.none
            self.organic_settler[i].unit.add_momentum_balances(
                balance_type=organic_pressure_balance_type,
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

        # if there is no energy balance, define function for temperature outlet
        if (
            self.config.aqueous_stream.has_energy_balance == False
            or self.config.organic_stream.has_energy_balance == False
        ):

            def settler_temperature_outlet(b, t):
                x0 = b.length_domain.first()
                xf = b.length_domain.last()
                return (
                    b.properties[t, x0].temperature == b.properties[t, xf].temperature
                )

        # if there is no pressure balance, define function for pressure outlet
        if (
            self.config.aqueous_stream.has_pressure_balance == False
            or self.config.organic_stream.has_pressure_balance == False
        ):

            def settler_pressure_outlet(b, t):
                x0 = b.length_domain.first()
                xf = b.length_domain.last()
                return b.properties[t, x0].pressure == b.properties[t, xf].pressure

        for i in self.elements:

            for j in ["aqueous", "organic"]:

                # Declare temperature constraint for the settlers if it does not have energy balance
                if getattr(self.config, f"{j}_stream").has_energy_balance == False:
                    model = getattr(self, f"{j}_settler")[i].unit
                    setattr(
                        model,
                        f"temperature_constraint",
                        Constraint(
                            self.flowsheet().time, rule=settler_temperature_outlet
                        ),
                    )

                # Declare pressure constraint for the settlers if it does not have pressure balance
                if getattr(self.config, f"{j}_stream").has_pressure_balance == False:
                    model = getattr(self, f"{j}_settler")[i].unit
                    setattr(
                        model,
                        f"pressure_constraint",
                        Constraint(self.flowsheet().time, rule=settler_pressure_outlet),
                    )

                # Declare arcs connecting mixers to their corresponding settler tanks
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
