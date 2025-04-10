#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Solar Evaporation Pond
======================

Author: Andrew Lee

The Evaporation Pond model is a general purpose model for representing evaporation ponds
(both solar and enhanced). Currently, the model only supports steady-state operation
(although it could be easily extended to dynamics). The model also only supports
material balances with equilibrium reactions, as energy and momentum balances are less
meaningful given the open-air nature of these operations and the long-duration time
constants involved.

Configuration Arguments
-----------------------

When creating an instance of an Evaporation Pond model, the user may specify the name of the solvent phase
used by the associated property package, which is used to determine which component requires an evaporation
term in the material balances. The default name is 'H2O', and only a single solvent is supported.

Evaporation Pond models also require a reaction package to define the precipitation reactions which occur in the
system. However, unlike normal reaction packages, the equilibrium constraints are written by the unit model as
it is necessary to include a conditional check for sub-saturated reactions. Thus, the reaction package should not
create any equilibrium constraints, and need only meet the following conditions:

Reaction Parameter Block:

* Define the set of equilibrium reaction names (``equilibrium_reaction_idx``), and
* define the stoichiometry coefficients for the equilibrium reactions (``equilibrium_reaction_stoichiometry``).

Reaction Block:

* Define the solubility products for each reaction. Each reaction requires a separate Var or Param
  named ``solubility_product_X`` where ``X`` is the identifier for each reaction in ``equilibrium_reaction_idx``.
  Separate components are required as the units of measurement for each reaction may be different.
* All reactions are assumed to use a standard power law form.

Degrees of Freedom
------------------

Evaporation Pond models have three degrees of freedom in addition to the feed stream state. The most common
variables to specify are:

* the surface area of the evaporation pond, ``surface_area``,
* the average depth of liquid in the pond, ``average_pond_depth``, and
* the rate of evaporation in units of length/time (or volume evaporated  per unit area per time), ``evaporation_rate``.

Model Structure
---------------

The Evaporation Pond unit model does not use Control Volumes, and writes custom material balances and
reaction constraints. The Evaporation Pond model has one inlet and one outlet ``Port`` (named
``inlet`` and ``outlet`` respectively).

Variables
---------

Evaporation Pond models have the following Variables.

================= =================== =====================
Variable          Name                Indexing Set(s)
================= =================== =====================
:math:`A_{t}`     surface_area        time
:math:`D_{t}`     average_pond_depth  time
:math:`V_{t}`     volume              time
:math:`e_{t}`     evaporation_rate    time
:math:`W_{t}`     water_loss_rate     time
:math:`P_{t, j}`  precipitation_rate  time, component list
:math:`X_{t, r}`  reaction_extent     time, reaction list
================= =================== =====================

Params
------

Evaporation Pond models have the following Params, which are used in the complementarity constraint
for precipitation.

================== ======== ============== =================================================================
Parameter          Name     Indexing Set   Notes
================== ======== ============== =================================================================
:math:`\epsilon`   eps                     Smoothing parameter for smooth_max
:math:`norm_{r}`   s_norm   reaction list  Normalizing factor, should match magnitude of solubility product
:math:`scale_{r}`  s_scale  reaction list  Scaling factor for reaction product, Q = Ksp - f(C)
================== ======== ============== =================================================================

Constraints
-----------

Evaporation Pond models have the following constraints.

Component balances:

For the solvent:

.. math:: 0 = F_{in,t,j} - F_{out,t,j} + P_{t_j} + W_{t}

For other species:

.. math:: 0 = F_{in,t,j} - F_{out,t,j} + P_{t_j}

where :math:`F_{in,t,j}` and :math:`F_{out,t,j}` are the flow rate of component :math:`j` in to and out of
the pond at time :math:`t`.

Stoichiometry Constraint:

.. math:: P_{t, j} = \sum_r{n_{r,j} \times X_{t,r}}

where :math:`n_{r,j}` is the stoichiometric coefficient for component :math:`j` in reaction :math:`r`.

Equilibrium Constraint:

.. math:: Q_{t, r} - \max{0, Q_{t,r}-\bar{X}_{t,r}} = 0

where :math:`Q_{t,r} = \ln{K_{t,r}} - \sum_j{-n_{r,j} \times \ln{C_{t,j}}}` and :math:`K_{t,r}` is the
solubility product for reaction :math:`r` at time :math:`t`, :math:`C_{t,j}` is the concentration
of species :math:`j` at time :math:`t`, and :math:`\bar{X}_{t,r} = scale_{r} \times \frac{X_{t,r}}{X_{t,r} + norm_{r}}`
The math:`max` operator is approximated using a smooth maximum,

Evaporation Constraint:

.. math:: W_{t} = e_{t} \times A_{t}

Volume Constraint:

.. math:: V_{t} = D_{t} \times A_{t}

"""

from pyomo.common.config import ConfigDict, ConfigValue, In
from pyomo.environ import Block, Param, Var, log, units

import idaes.logger as idaeslog
from idaes.core import (
    MaterialFlowBasis,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.initialization.initializer_base import ModularInitializerBase
from idaes.core.util import StoreSpec, from_json, to_json
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
from idaes.core.util.math import smooth_max

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

    CONFIG = ConfigDict()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. Equilibrium Reactors do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Equilibrium reactors do not have defined volume, thus
    this must be False.""",
        ),
    )
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
            description="Dict of arguments to use for constructing property blocks",
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
        "reaction_package_args",
        ConfigDict(
            implicit=True,
            description="Dict of arguments to use for constructing reaction blocks",
            doc="""A ConfigDict with arguments to be passed to reaction block(s)
        and used when constructing these,
        **default** - None.
        **Valid values:** {
        see reaction package for documentation.}""",
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
                f"EvaporationPond model, however this is not a valid component name "
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
        self.reactions = self.config.reaction_package.build_reaction_block(
            self.flowsheet().time,
            state_block=self.properties_out,
            **self.config.reaction_package_args,
        )

        # Construct Ports
        self.add_port("inlet", self.properties_in)
        self.add_port("outlet", self.properties_out)

        # Get indexing sets
        time = self.flowsheet().time
        comp_set = self.config.property_package.component_list
        pc_set = self.config.property_package.get_phase_component_set()
        rxn_basis = self.reactions[time.first()].get_reaction_rate_basis()

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
            raise PropertyPackageError(
                f"{self.name} - Property package uses a flow basis of 'other'. "
                "EvaporationPond model only supports mass or molar bases."
            )

        if rxn_basis is MaterialFlowBasis.molar:
            rxn_units = uom.FLOW_MOLE
        elif rxn_basis is MaterialFlowBasis.mass:
            rxn_units = uom.FLOW_MASS
        else:
            raise PropertyPackageError(
                f"{self.name} - Reaction package uses a reaction basis other than mass or "
                "mole. EvaporationPond model only supports mass or molar bases."
            )

        # Build unit model level variables
        self.surface_area = Var(
            time,
            initialize=1000,
            bounds=(1, None),
            units=uom.AREA,
            doc="Surface area of evaporation pond",
        )
        self.average_pond_depth = Var(
            time,
            initialize=1,
            bounds=(1e-6, None),
            units=uom.LENGTH,
            doc="Average depth of evaporation pond",
        )
        self.volume = Var(
            time,
            initialize=1000,
            bounds=(1e-6, None),
            units=uom.VOLUME,
            doc="Volume of liquid in pond",
        )

        self.evaporation_rate = Var(
            time,
            initialize=0,
            bounds=(0, None),  # TODO: Need to relax this for dynamics
            units=units.mm / units.day,
            doc="Evaporation rate of water",
        )
        self.water_loss_rate = Var(
            time,
            initialize=0,
            bounds=(None, 0),  # TODO: Need to relax this for dynamics
            units=mb_units,
            doc="Flowrate of water lost to evaporation",
        )

        self.precipitation_rate = Var(
            time,
            comp_set,
            initialize=0,
            bounds=(None, 0),  # TODO: Need to relax this for dynamics
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
            bounds=(0, None),  # TODO: Need to relax this for dynamics
            units=rxn_units,
            doc="Extent of equilibrium reaction",
        )

        self.eps = Param(
            mutable=True,
            initialize=1e-4,
            units=units.dimensionless,
            doc="Smoothing parameter for smooth maximum",
        )

        self.s_norm = Param(
            rxn_idx,
            mutable=True,
            initialize=1,
            units=rxn_units,
            doc="Normalizing factor for solid precipitation term. Should match magnitude of Ksp",
        )

        self.s_scale = Param(
            rxn_idx,
            mutable=True,
            initialize=1,
            units=units.dimensionless,
            doc="Scaling factor for solid precipitation term w.r.t saturated status Q = Ksp - f(C)",
        )

        @self.Constraint(time, pc_set, doc="Stoichiometry constraints")
        def stoichiometry_constraint(b, t, p, j):
            # Catch inerts and avoid sum expression
            if all(rxn_stoic[r, p, j] == 0 for r in rxn_idx):
                return b.precipitation_rate[t, j] == 0 * mb_units

            exp = sum(
                rxn_stoic[r, p, j] * b.reaction_extent[t, r]
                for r in rxn_idx
                if rxn_stoic[r, p, j] != 0
            )
            # Convert units if required
            if (
                flow_basis is MaterialFlowBasis.mass
                and rxn_basis is MaterialFlowBasis.molar
            ):
                exp = units.convert(exp * b.properties_out[t].mw[j], to_units=mb_units)
            elif (
                flow_basis is MaterialFlowBasis.molar
                and rxn_basis is MaterialFlowBasis.mass
            ):
                exp = units.convert(exp / b.properties_out[t].mw[j], to_units=mb_units)
            return b.precipitation_rate[t, j] == exp

        @self.Constraint(time, rxn_idx, doc="Equilibrium constraint")
        def equilibrium_constraint(b, t, r):
            Ksp = getattr(b.reactions[t], "solubility_product_" + r)
            k_units = Ksp.get_units()

            if rxn_basis == MaterialFlowBasis.molar:
                conc = b.properties_out[t].conc_mole_comp
            else:
                conc = b.properties_out[t].conc_mass_comp

            # Equilibrium expression
            # Assume equilibrium reactions will always be defined on a molar basis
            e = sum(
                -rxn_stoic[r, phase_name, j] * log(conc[j] / units.get_units(conc[j]))
                for j in comp_set
                if rxn_stoic[r, phase_name, j] != 0.0
            )

            # Complementarity formulation to support conditional precipitation
            k = log(Ksp / k_units)

            # TODO: For steady-state, it is sufficient to use extent of precipitation reactions
            # TODO: For dynamics, need to track amount of solids in basin and use that instead
            # This is because solids can redissolve in dynamic mode, thus we need to check
            # for the presence of solids, not just the occurrence of precipitation
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
        # TODO: For dynamics, show inventory of solids as well

        return {"vars": var_dict}
