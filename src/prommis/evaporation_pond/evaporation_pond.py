r"""
Solar Evaporation Pond
======================

Author: Andrew Lee

"""

from pyomo.environ import (
    Block,
    log,
    Param,
    units,
    Var,
)
from pyomo.common.config import ConfigDict, ConfigValue

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_max
from idaes.core.initialization.initializer_base import ModularInitializerBase
from idaes.core.util import to_json, from_json, StoreSpec
import idaes.logger as idaeslog

__author__ = "Andrew Lee"


class EvaporationPondInitializer(ModularInitializerBase):
    """
    Initializer object for EvaporationPond unit models.

    This Initializer initializes the inlet StateBlock, maps the solution of
    this onto the outlet StateBlock, and then solves the full model.

    """

    CONFIG = ModularInitializerBase.CONFIG()

    def initialization_routine(
        self,
        model: Block,
        plugin_initializer_args: dict = None,
    ):
        """
        Common initialization routine for EvaporationPond unit models.

        Args:
            model: Pyomo Block to be initialized
            plugin_initializer_args: dict-of-dicts containing arguments to be passed to plug-in Initializers.
                Keys should be submodel components.

        Returns:
            Pyomo solver results object
        """
        # The default initialization_routine is sufficient
        return super().initialization_routine(
            model=model,
            plugin_initializer_args=plugin_initializer_args,
        )

    def initialize_main_model(
        self,
        model: Block,
    ):
        """
        Initialization routine for EvaporationPond unit models.

        Args:
            model: current model being initialized

        Returns:
            Pyomo solver results object from solve of main model

        """
        # Get logger
        _log = self.get_logger(model)

        # Initialize inlet properties - inlet state should already be fixed
        prop_init = self.get_submodel_initializer(model.properties_in)

        if prop_init is not None:
            prop_init.initialize(
                model=model.properties_in,
                output_level=self.get_output_level(),
            )
            _log.info_high("Inlet properties initialization complete.")

        # Map solution from inlet properties to outlet properties
        state = to_json(
            model.properties_in,
            wts=StoreSpec().value(),
            return_dict=True,
        )
        from_json(
            model.properties_out,
            sd=state,
            wts=StoreSpec().value(only_not_fixed=True),
        )

        # Solve main model
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = self._get_solver().solve(model, tee=slc.tee)

        _log.info_high(
            f"Evaporation Pond initialization {idaeslog.condition(results)}."
        )

        return results


@declare_process_block_class("EvaporationPond")
class EvaporationPondData(UnitModelBlockData):
    """
    Evaporation Pond Unit Model Class
    """

    default_initializer = EvaporationPondInitializer

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
            description="Reaction package to use",
            doc="""Reaction parameter object used to define precipitation reactions,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **ReactionParameterObject** - a ReactionParameterBlock object.}""",
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
        phase_name = self.config.property_package.phase_list.at(1)

        # Build property blocks
        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )
        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=False,
            **self.config.property_package_args,
        )

        # Construct Ports
        self.add_port("inlet", self.properties_in)
        self.add_port("outlet", self.properties_out)

        # Get units of measurement
        flow_basis = self.properties_in[
            self.flowsheet().time.first()
        ].get_material_flow_basis()
        uom = self.config.property_package.get_metadata().derived_units

        if flow_basis is MaterialFlowBasis.molar:
            mb_units = uom.FLOW_MOLE
        elif flow_basis is MaterialFlowBasis.mass:
            mb_units = uom.FLOW_MASS
        else:
            mb_units = None

        # Get indexing sets
        time = self.flowsheet().time
        comp_set = self.config.property_package.component_list
        pc_set = self.config.property_package.get_phase_component_set()

        # Build unit model level variables
        self.surface_area = Var(
            time,
            initialize=1000,
            units=uom.AREA,
            doc="Surface area of evaporation pond",
        )
        self.average_pond_depth = Var(
            time,
            initialize=1,
            units=uom.LENGTH,
            doc="Average depth of evaporation pond",
        )
        self.volume = Var(
            time, initialize=1000, units=uom.VOLUME, doc="Volume of liquid in pond"
        )

        self.evaporation_rate = Var(
            time,
            initialize=0,
            units=units.mm / units.day,
            doc="Evaporation rate of water",
        )
        self.water_loss_rate = Var(
            time,
            initialize=0,
            units=mb_units,
            doc="Flowrate of water lost to evaporation",
        )

        self.precipitation_rate = Var(
            time,
            comp_set,
            initialize=0,
            units=mb_units,
            doc="Rate of precipitation of species",
        )

        # Add unit level constraints
        @self.Constraint(
            time,
            doc="Evaporation constraint",
        )
        def evaporation_constraint(b, t):
            rhs = -b.evaporation_rate[t] * b.surface_area[t]
            if flow_basis is MaterialFlowBasis.molar:
                rhs = rhs * 1000 * units.kg / units.m**3 / (18 * units.g / units.mol)
            elif flow_basis is MaterialFlowBasis.mass:
                rhs = rhs * 1000 * units.kg / units.m**3
            # No action if flow basis other

            return b.water_loss_rate[t] == units.convert(
                rhs,
                to_units=mb_units,
            )

        # Add constraint relating area to volume
        @self.Constraint(time, doc="Volume constraint")
        def volume_constraint(b, t):
            return b.volume[t] == b.surface_area[t] * b.average_pond_depth[t]

        # Material balances
        # TODO: For dynamics, need to include tracking of solids
        @self.Constraint(time, pc_set, doc="Component balances")
        def component_balances(b, t, p, j):
            rhs = (
                b.properties_in[t].get_material_flow_terms(p, j)
                - b.properties_out[t].get_material_flow_terms(p, j)
                + b.precipitation_rate[t, j]
            )

            if j == self.config.solvent_id:
                rhs += b.water_loss_rate[t]

            # TODO: Only support steady-state for now
            lhs = 0

            return lhs == rhs

        # Equilibrium reactions and precipitation rate
        # Reaction parameters
        rxn_idx = self.config.reaction_package.equilibrium_reaction_idx
        rxn_stoic = self.config.reaction_package.equilibrium_reaction_stoichiometry

        self.reaction_extent = Var(
            time,
            rxn_idx,
            initialize=0,
            units=mb_units,
            doc="Extent of equilibrium reaction",
        )

        self.eps = Param(
            mutable=True, initialize=1e-4, doc="Smoothing parameter for smooth maximum"
        )

        self.s_norm = Param(
            rxn_idx,
            mutable=True,
            initialize=1e6,  # TODO: Value should match Ksp
            units=mb_units,  # Units are based on mb_units... how to match value?
            doc="Normalizing factor for solid precipitation term",
        )

        self.s_scale = Param(
            rxn_idx,
            mutable=True,
            initialize=1,
            doc="Scaling factor for solid precipitation term w.r.t saturated status Q = Ksp - f(C)",
        )

        @self.Constraint(time, pc_set, doc="Stoichiometry constraints")
        def stoichiometry_constraints(b, t, p, j):
            # TODO: Handle mass-mole conversions as necessary
            return b.precipitation_rate[t, j] == sum(
                rxn_stoic[r, p, j] * b.reaction_extent[t, r] for r in rxn_idx
            )

        @self.Constraint(time, rxn_idx, doc="Equilibrium constraint")
        def equilibrium_constraint(b, t, r):
            # TODO: Time indexing/varying Ksp?
            Ksp = getattr(b.config.reaction_package, "solubility_constant_" + r)
            k_units = Ksp.get_units()
            c_units = b.properties_out[t].conc_mass_comp.get_units()

            # Equilibrium expression
            # TODO: Different bases
            e = sum(
                -rxn_stoic[r, phase_name, j]
                * log(b.properties_out[t].conc_mass_comp[j] / c_units)
                for j in comp_set
                if rxn_stoic[r, phase_name, j] != 0.0
            )

            # Complementarity formulation to support conditional precipitation
            k = log(Ksp / k_units)

            # TODO: For steady-state, it is sufficient to use extent of precipitation reactions
            # TODO: For dynamics, need to track amount of solids in basin and use that instead
            # This is because solids can redissolve in dynamic mode, thus we need to check
            # for the presence of solids, not just the occurance of precipitation
            s = (
                b.s_scale[r]
                * b.reaction_extent[t, r]
                / (b.reaction_extent[t, r] + b.s_norm[r])
            )
            Q = k - e

            return Q - smooth_max(0, Q - s, b.eps) == 0

        # TODO: Energy and momentum balances
        # TODO: Rate based reactions?

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Surface Area"] = self.surface_area[time_point]
        var_dict["Average Depth"] = self.average_pond_depth[time_point]
        var_dict["Volume"] = self.volume[time_point]
        var_dict["Evaporation Rate"] = self.evaporation_rate[time_point]
        var_dict["Water Loss Rate"] = self.water_loss_rate[time_point]
        for j in self.properties_out.component_list:
            var_dict[f"Precipitation Rate {j}"] = self.precipitation_rate[time_point, j]

        return {"vars": var_dict}
