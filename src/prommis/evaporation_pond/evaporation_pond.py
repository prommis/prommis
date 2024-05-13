r"""
Solar Evaporation Pond
======================

Author: Andrew Lee

"""

from pyomo.environ import (
    Reference,
    Var,
    units,
    Param,
    Set,
    log,
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
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_max


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
        # Reaction Index
        self.equilibrium_reaction_idx = Set(initialize=["P1", "P2"])

        # Reaction Stoichiometry
        self.equilibrium_reaction_stoichiometry = {
            ("P1", "liquid", "Li"): 0,
            ("P1", "liquid", "Na"): -1,
            ("P1", "liquid", "Cl"): -1,
            ("P1", "liquid", "H2O"): 0,
            ("P2", "liquid", "Li"): -1,
            ("P2", "liquid", "Na"): 0,
            ("P2", "liquid", "Cl"): -1,
            ("P2", "liquid", "H2O"): 0,
        }

        # TODO: Time indexing?
        self.solubility_constant = Param(
            self.equilibrium_reaction_idx,
            default={"P1": 1.5e6, "P2": 3e6},
            units=units.mg**2 / units.L**2,
            mutable=True,
        )

        self.reaction_extent = Var(
            time,
            self.equilibrium_reaction_idx,
            initialize=0,
            units=mb_units,
            doc="Extent of equilibrium reaction",
        )

        self.eps = Param(
            mutable=True, initialize=1e-4, doc="Smoothing parameter for smooth maximum"
        )

        self.s_norm = Param(
            mutable=True,
            initialize=1e6,  # TODO: Value should match Ksp
            units=mb_units,  # Units are based on mb_units... how to match value?
            doc="Normalizing factor for solid precipitation term",
        )

        self.s_scale = Param(
            mutable=True,
            initialize=1,
            doc="Scaling factor for solid precipitation term w.r.t saturated status Q = Ksp - f(C)",
        )

        @self.Constraint(time, pc_set, doc="Stoichiometry constraints")
        def stoichiometry_constraints(b, t, p, j):
            # TODO: Handle mass-mole conversions as necessary
            return b.precipitation_rate[t, j] == sum(
                b.equilibrium_reaction_stoichiometry[r, p, j]
                * b.reaction_extent[t, r]
                for r in b.equilibrium_reaction_idx
            )

        @self.Constraint(time, self.equilibrium_reaction_idx, doc="Equilibrium constraint")
        def equilibrium_constraint(b, t, r):
            k_units = b.solubility_constant.get_units()
            c_units = b.properties_out[t].conc_mass_comp.get_units()
            # TODO: different Ks to allow for different units
            # TODO: phase name
            e = sum(
                -b.equilibrium_reaction_stoichiometry[r, "liquid", j]
                * log(b.properties_out[t].conc_mass_comp[j]/c_units)
                for j in comp_set
                if b.equilibrium_reaction_stoichiometry[r, "liquid", j] != 0.0
            )

            # Complementarity formulation to support conditional precipitation
            k = log(b.solubility_constant[r] / k_units)
            s = b.s_scale * b.reaction_extent[t, r] / (b.reaction_extent[t, r] + b.s_norm)
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
