r"""
Solar Evaporation Pond
======================

Author: Andrew Lee

"""

from pyomo.environ import (
    Reference,
    Var,
    units,
)
from pyomo.common.config import ConfigDict, ConfigValue

from idaes.core import (
    ControlVolume0DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.config import is_reaction_parameter_block
from idaes.core.util.exceptions import ConfigurationError


@declare_process_block_class("EvaporationPond")
class EvaporationPondData(UnitModelBlockData):
    """
    Leaching Train Unit Model Class
    """

    # TODO: Set default initializer

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
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
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            domain=is_reaction_parameter_block,
            description="Reaction package defining precipitation equilibria.",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for reaction package.",
        ),
    )
    CONFIG.declare(
        "solvent_id",
        ConfigValue(
            default="H2O",
            domain=str,
            description="Name of solvent component (default='H2O').",
        ),
    )

    def build(self):
        """
        Build method for EvaporationPond unit model.
        """
        super().build()

        # Verify solvent ID
        if self.config.solvent_id not in self.config.property_package.component_list:
            raise ConfigurationError(
                f"{self.name} - {self.config.solvent_id} was set as the solvent_id for "
                f"an EvaporationPond model, however this is not a valid component name "
                f"in the property package provided."
            )

        # Verify single phase
        if len(self.config.property_package.phase_list) > 1:
            raise ConfigurationError(
                f"{self.name} - EvaporationPond model only supports a single phase. However, "
                f"the property package provided includes "
                f"{len(self.config.property_package.phase_list)} phases."
            )
        else:
            phase_name = self.config.property_package.phase_list.at(1)

        # Build control volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            reaction_package=self.config.reaction_package,
            reaction_package_args=self.config.reaction_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)
        self.control_volume.add_reaction_blocks(has_equilibrium=True)
        self.control_volume.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_equilibrium_reactions=True,
            has_mass_transfer=True,
        )

        # TODO: Energy and momentum balances
        # TODO: Rate based reactions?

        # Construct Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Get units of measurement
        flow_basis = self.control_volume.properties_in[
            self.flowsheet().time.first()
        ].get_material_flow_basis()
        uom = self.config.property_package.get_metadata().derived_units

        if flow_basis is MaterialFlowBasis.molar:
            mb_units = uom.FLOW_MOLE
        elif flow_basis is MaterialFlowBasis.mass:
            mb_units = uom.FLOW_MASS
        else:
            mb_units = None

        # Build unit model level variables
        self.surface_area = Var(
            self.flowsheet().time,
            units=uom.AREA,
            doc="Surface area of evaporation pond",
        )

        self.evaporation_rate = Var(
            self.flowsheet().time,
            units=units.mm / units.day,
            doc="Evaporation rate of water",
        )

        self.water_loss_rate = Reference(
            self.control_volume.mass_transfer_term[
                :, phase_name, self.config.solvent_id
            ]
        )

        # Add unit level constraints
        @self.Constraint(
            self.flowsheet().time,
            self.control_volume.properties_out.phase_component_set,
            doc="Evaporation constraint",
        )
        def evaporation_constraint(b, t, p, j):
            if j == self.config.solvent_id:
                rhs = -b.evaporation_rate[t] * b.surface_area[t]
                if flow_basis is MaterialFlowBasis.molar:
                    rhs = (
                        rhs * 1000 * units.kg / units.m**3 / (18 * units.g / units.mol)
                    )
                elif flow_basis is MaterialFlowBasis.mass:
                    rhs = rhs * 1000 * units.kg / units.m**3
                # No action if flow basis other

                return b.control_volume.mass_transfer_term[t, p, j] == units.convert(
                    rhs,
                    to_units=mb_units,
                )
            else:
                # TODO: User defined loss terms for other species?
                return b.control_volume.mass_transfer_term[t, p, j] == 0

        if hasattr(self.control_volume, "volume"):
            # Add constraints relating area to volume
            self.average_pond_depth = Var(
                self.flowsheet().time,
                units=uom.LENGTH,
                doc="Average depth of evaporation pond",
            )

            @self.Constraint(self.flowsheet().time, doc="Volume constraint")
            def volume_constraint(b, t):
                return (
                    b.control_volume.volume[t]
                    == b.surface_area[t] * b.average_pond_depth[t]
                )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Surface Area"] = self.surface_area[time_point]
        if hasattr(self.control_volume, "volume"):
            var_dict["Average Depth"] = self.average_pond_depth[time_point]
            var_dict["Volume"] = self.control_volume.volume[time_point]
        var_dict["Evaporation Rate"] = self.evaporation_rate[time_point]
        var_dict["Water Loss Rate"] = self.water_loss_rate[time_point]
        for p, j in self.control_volume.properties_out.phase_component_set:
            var_dict[f"Precipitation Rate {j}"] = (
                self.control_volume.equilibrium_reaction_generation[time_point, p, j]
            )

        return {"vars": var_dict}
