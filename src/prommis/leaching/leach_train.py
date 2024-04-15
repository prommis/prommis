"""
Unit model for a leaching train.

Train is modeled as a series of well mixed tank reactors.

Authors: Andrew Lee
"""

from pyomo.environ import (
    Block,
    Constraint,
    Var,
    units,
)
from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.network import Port

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.unit_models.mscontactor import MSContactor
from idaes.core.initialization import ModularInitializerBase


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
        return msc_init.initialize(model.mscontactor)


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
            description="Solid stream properties",
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
