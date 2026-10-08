#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Per-component slurry particle-size-distribution property package.

This package represents a slurry through component mass flows, with sized
solids resolved over a common particle-size mesh. Unit models balance these
flows; the package derives PSD, volume, and density properties from the state
and configured component densities.

In the tables, ``s`` is a sized solid, ``u`` an unsized solid, ``j`` a component
of the indicated phase, ``v`` a vapor component, ``p`` a phase, and ``q`` a
passing fraction. Size interval ``k`` follows the ascending mesh, with bin 0
finest. Unsized solids have no PSD and default to an empty list. Liquids
default to ``["H2O"]``; vapor is optional and defaults to off.

State variables
---------------

=================================== ====== ==================================================
Name                                Units  Meaning and basis
=================================== ====== ==================================================
``flow_mass_sized_comp_size[s, k]`` kg/s   State: retained mass flow of sized solid ``s`` in
                                           interval ``k``.
``flow_mass_unsized_comp[u]``       kg/s   State: mass flow of unsized solid ``u``, when
                                           configured.
``flow_mass_liquid_comp[j]``        kg/s   State: mass flow of liquid component ``j``.
``flow_mass_vapor_comp[v]``         kg/s   State: mass flow of vapor component ``v``, when
                                           configured.
``temperature``                     K      State temperature.
``pressure``                        Pa     State pressure.
=================================== ====== ==================================================

Parameters and derived properties
---------------------------------

=================================== ====== ==================================================
Name                                Units  Meaning and basis
=================================== ====== ==================================================
``size_edges[k]``                   m      Parameter: interval boundaries; ``N+1`` edges for
                                           ``N`` intervals.
``size_char[k]``                    m      Parameter: geometric-mean interval size; uses
                                           ``bottom_size`` for a zero finest edge.
``dens_mass_solid_comp[j]``         kg/m^3 Parameter: configured density of sized or unsized
                                           solid ``j``.
``dens_mass_liquid_comp[j]``        kg/m^3 Parameter: configured density of liquid component
                                           ``j``.
``flow_mass_phase[p]``              kg/s   Total phase flow; ``Sol`` includes sized and
                                           unsized solids, ``Liq`` all liquids, and optional
                                           ``Vap`` all vapors.
``flow_mass_sized_comp[s]``         kg/s   Flow of sized solid ``s``, summed over its
                                           intervals.
``flow_mass_solid_comp[j]``         kg/s   Flow of solid ``j``, whether sized or unsized.
``flow_mass_size[k]``               kg/s   Retained flow in interval ``k``, summed over sized
                                           solids.
``flow_mass_sized``                 kg/s   Total sized-solids flow; excludes unsized solids.
``mass_frac_solid_comp[j]``         –      Composition of solid ``j`` on a total-solids mass
                                           basis.
``mass_frac_size[k]``               –      Retained fraction in interval ``k`` on a
                                           sized-solids mass basis.
``mass_frac_size_comp[s, k]``       –      Fraction of sized solid ``s`` retained in interval
                                           ``k``; its PSD.
``mass_frac_comp_size[s, k]``       –      Component composition within interval ``k`` on a
                                           sized-solids basis.
``cum_passing_size[k]``             –      Sized-solids mass fraction passing the upper edge
                                           of interval ``k``.
``cum_passing_comp_size[s, k]``     –      Mass fraction of sized solid ``s`` passing that
                                           upper edge.
``percentile_size[q]``              m      Smooth size at passing fraction ``q`` for all
                                           sized solids.
``percentile_size_comp[s, q]``      m      Smooth size at passing fraction ``q`` for sized
                                           solid ``s``.
``flow_vol_solid``                  m^3/s  Sum of each sized and unsized solid's mass flow
                                           divided by its density.
``flow_vol_liquid``                 m^3/s  Sum of each liquid component's mass flow divided
                                           by its density.
``flow_vol_slurry``                 m^3/s  Solid plus liquid volume flow; excludes vapor.
``vol_frac_solid_slurry``           –      All-solids volume fraction of the solid-liquid
                                           slurry.
``dens_mass_solid``                 kg/m^3 Total solids mass flow divided by total solids
                                           volume flow.
``dens_mass_solid_mass_weighted``   kg/m^3 Sum of each solid's mass flow times density
                                           divided by total solids mass flow; differs from
                                           ``dens_mass_solid``.
``dens_mass_slurry``                kg/m^3 Solid plus liquid mass flow divided by slurry
                                           volume flow; excludes vapor.
=================================== ====== ==================================================

Mass fractions are calculated from component flows, without a separate
constraint forcing them to sum to one. Unit models supply the full
component-by-size flow balances. Retained mass fractions, per-component PSD
properties, and ``dens_mass_solid_mass_weighted`` are built when first accessed.

Small fixed terms in mass and volume flow denominators keep ratios finite at
zero flow, but can bias them near zero. Percentiles have no physical meaning
when the selected PSD contains no sized solids. ``validate_feed`` accepts
water-only feeds but requires finite, nonnegative flows, positive finite
temperature and pressure, and positive finite total flow across all phases.
For reporting, ``percentile_size_report`` sets ``size_m`` only when the selected
PSD contains sized solids and the smooth percentile agrees with linear bin
interpolation within ``rel_tol``.

Solid and liquid volumes are added using constant component densities; vapor
and volume changes on mixing are excluded from slurry volume. Solid densities
must describe the solid material, excluding interparticle voids and pores
occupied by liquid counted in ``flow_mass_liquid_comp``. A density that includes
those pore volumes would count that liquid volume twice. Porosity is compatible
with this model when solid and liquid volumes do not overlap. The package does
not model pore uptake, swelling, dissolution, or changes in liquid density
with dissolved composition. Those effects need a property model that
represents them and, when mass moves between phases, a phase-transfer model.
"""

__author__ = "Daison Yancy Caballero"

import math

from pyomo.common.config import ConfigValue
from pyomo.environ import NonNegativeReals, Param, RangeSet, Set, Var, units, value

from idaes.core import (
    Component,
    LiquidPhase,
    MaterialBalanceType,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    SolidPhase,
    StateBlock,
    StateBlockData,
    VaporPhase,
    declare_process_block_class,
)
from idaes.core.initialization import InitializerBase
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.math import smooth_min

from prommis.comminution.core import size_mesh
from prommis.comminution.core.config_utils import _cfg_float, _is_real_number
from prommis.comminution.core.psd_math import (
    cumulative_passing,
    size_at_passing as size_at_passing_numeric,
    size_at_passing_smooth,
)
from prommis.comminution.core.smooth_functions import EPS_SMOOTH
from prommis.comminution.core.unit_utils import (
    MASS_FLOW_EPS,
    VOLUME_FLOW_EPS,
    FactorSource,
    declare_factor_source,
    fixed_value_or_none,
)

# Density of water at 298.15 K and 101325 Pa (kg/m^3), IAPWS-95 (Wagner, W., and
# Pruss, A. (2002). J. Phys. Chem. Ref. Data, 31(2), 387-535).
_WATER_DENSITY_25C = 997.048

# Reject names that would collide with phase or parameter-block attributes.
_RESERVED_NAMES = (
    "component_list",
    "Sol",
    "Liq",
    "Vap",
    "sized_solid_list",
    "unsized_solid_list",
    "all_solid_list",
    "liquid_list",
    "vapor_list",
    "size_interval_set",
    "size_edge_index",
    "size_edges",
    "size_char",
    "dens_mass_solid_comp",
    "dens_mass_liquid_comp",
    "assert_same_mesh",
    "size_edges_m",
    "size_char_m",
    "size_at_passing_from_flows",
    "percentile_targets",
    "_edge_values",
    "_size_char_values",
)

_FLOW_VAR_NAMES = (
    "flow_mass_sized_comp_size",
    "flow_mass_unsized_comp",
    "flow_mass_liquid_comp",
    "flow_mass_vapor_comp",
)


class SolidPSDScaler(CustomScalerBase):
    """Scale slurry PSD state variables from current values or fixed inputs.

    Current-values mode uses current flow magnitudes with a stream-relative
    floor. Scaling factors depends on model initialization.
    Input-based mode reads only fixed flows; unfixed flows receive a factor of
    1.0.
    Both modes use fixed temperature and pressure nominals.
    """

    CONFIG = declare_factor_source(CustomScalerBase.CONFIG())

    INPUT_BASED_DEFAULT_FLOW_FACTOR = 1.0

    # Floor on each flow scale: EPS_REL times the sum of the absolute stream flows
    # (fixed values only in input-based mode; a zero total counts as 1.0).
    EPS_REL = 1e-8
    TEMPERATURE_NOMINAL = 300.0  # K
    PRESSURE_NOMINAL = 1.0e5  # Pa

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """Scale state variables using the configured factor source."""
        if self.config.factor_source == FactorSource.input_based:
            self._flow_factors(
                model,
                fixed_value_or_none,
                overwrite,
                unset_factor=self.INPUT_BASED_DEFAULT_FLOW_FACTOR,
            )
        else:
            self._flow_factors(
                model, lambda var: value(var, exception=False), overwrite
            )
        self.set_variable_scaling_factor(
            model.temperature,
            self.temperature_nominal_factor(model),
            overwrite=overwrite,
        )
        self.set_variable_scaling_factor(
            model.pressure,
            self.pressure_nominal_factor(model),
            overwrite=overwrite,
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """No constraint scaling is needed for this state block."""
        pass

    # Scaling helpers

    @staticmethod
    def temperature_nominal_factor(blk):
        """1/nominal for a temperature Var, in the Var's own declared units."""
        return 1.0 / value(
            units.convert(
                SolidPSDScaler.TEMPERATURE_NOMINAL * units.K,
                to_units=units.get_units(blk.temperature),
            )
        )

    @staticmethod
    def pressure_nominal_factor(blk):
        """1/nominal for a pressure Var, in the Var's own declared units."""
        return 1.0 / value(
            units.convert(
                SolidPSDScaler.PRESSURE_NOMINAL * units.Pa,
                to_units=units.get_units(blk.pressure),
            )
        )

    @staticmethod
    def iter_flow_var_data(blk):
        """Yield every flow VarData of a state-block member (all phases)."""
        for name in _FLOW_VAR_NAMES:
            var = getattr(blk, name, None)
            if var is not None:
                yield from var.values()

    def _flow_factors(self, model, read, overwrite, unset_factor=None):
        """Set flow scaling factors from the values returned by ``read``.

        The factor is ``1 / max(abs(value), floor)``, where ``floor`` is
        ``EPS_REL`` times the sum of available absolute flow values. Use
        ``EPS_REL`` as the floor when that sum is zero. A second lower bound on
        the denominator keeps the scaling factor finite, even for extremely
        small flows. If ``read`` returns ``None``, use
        ``unset_factor`` when provided; otherwise use ``1 / floor``.
        """
        flows = list(self.iter_flow_var_data(model))
        values = [read(var) for var in flows]
        total = sum(abs(v) for v in values if v is not None)
        # 1e10: IDAES's default maximum variable scaling factor.
        maximum_factor = min(self.config.max_variable_scaling_factor, 1e10)
        floor = max(
            self.EPS_REL * (total if total > 0.0 else 1.0),
            1.0 / maximum_factor,
        )
        for var, v in zip(flows, values):
            if v is None and unset_factor is not None:
                sf = unset_factor
            else:
                sf = 1.0 / max(abs(v or 0.0), floor)
            self.set_variable_scaling_factor(var, sf, overwrite=overwrite)


class SolidPSDInitializer(InitializerBase):
    """Initializer for the slurry PSD state block.

    No solve is needed. Fully specified blocks are validated while the base
    initializer manages state fixing and restoration.
    """

    def initialization_routine(self, model):
        """Validate defined feed states without solving."""
        for idx in model:
            sbd = model[idx]
            if sbd.config.defined_state:
                sbd.validate_feed()
        return None


@declare_process_block_class("SolidPSDParameterBlock")
class SolidPSDParameterData(PhysicalParameterBlock):
    """Parameter block for the per-component slurry PSD package.

    Configuration:

        - ``size_edges``: ascending list of ``N+1`` size edges in meters (required).
        - ``bottom_size``: size in meters that replaces a zero finest edge for
          characteristic sizes; required when ``size_edges[0] == 0`` and then
          must satisfy ``0 < bottom_size < size_edges[1]``; rejected when
          ``size_edges[0] > 0``.
        - ``sized_solid_component_list``: non-empty list of sized solid components (required).
        - ``unsized_solid_component_list``: optional unsized solid components (default
          empty); disjoint from the sized list.
        - ``solid_density``: REQUIRED dict, one positive finite value per solid
          component (sized and unsized), kg/m^3.
        - ``liquid_component_list``: default ``["H2O"]``.
        - ``liquid_density``: per-liquid dict of densities, kg/m^3; defaults to
          ``{"H2O": 997.048}``. Keys must match the liquid list exactly;
          other liquid lists require a complete dict.
        - ``vapor_component_list``: default empty; when empty no vapor state exists.
        - ``percentile_targets``: passing fractions for the indexed percentile
          properties; defaults to ``(0.8,)``. The 80% passing size is always
          included.

    Use one parameter block instance across connected units to keep components,
    densities, and bottom-size definitions consistent. ``assert_same_mesh``
    compares cached size edges of separate instances within tolerance; it does
    not check those other properties or enforce shared-instance identity.
    The mesh Params are ``mutable=True`` only because Pyomo
    requires it for unit-bearing Params; post-build changes are unsupported.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "size_edges",
        ConfigValue(
            default=None,
            description="Strictly increasing finite size-mesh edges in meters "
            "(at least two).",
        ),
    )
    CONFIG.declare(
        "bottom_size",
        ConfigValue(
            default=None,
            description="Size in meters that replaces a zero first mesh edge "
            "when computing characteristic sizes; the mesh edge remains zero. "
            "Required when size_edges[0] == 0, and then must satisfy "
            "0 < bottom_size < size_edges[1]; rejected when size_edges[0] > 0.",
        ),
    )
    CONFIG.declare(
        "sized_solid_component_list",
        ConfigValue(
            default=None,
            description="Non-empty list of sized solid component names.",
        ),
    )
    CONFIG.declare(
        "unsized_solid_component_list",
        ConfigValue(
            default=None,
            description="Unsized solid component names (default: none).",
        ),
    )
    CONFIG.declare(
        "solid_density",
        ConfigValue(
            default=None,
            description="Required dict of positive densities (kg/m^3) for "
            "every sized and unsized solid component.",
        ),
    )
    CONFIG.declare(
        "liquid_component_list",
        ConfigValue(
            default=None,
            description="Liquid component names (default ['H2O']).",
        ),
    )
    CONFIG.declare(
        "liquid_density",
        ConfigValue(
            default={"H2O": _WATER_DENSITY_25C},
            description="Dict of positive liquid densities (kg/m^3); defaults "
            f"to {{'H2O': {_WATER_DENSITY_25C}}}. Keys must match the liquid "
            "component list exactly.",
        ),
    )
    CONFIG.declare(
        "vapor_component_list",
        ConfigValue(
            default=None,
            description="Vapor component names (default: none; no vapor phase "
            "or mass-flow variable).",
        ),
    )
    CONFIG.declare(
        "percentile_targets",
        ConfigValue(
            default=(0.8,),
            description="Non-empty passing fractions in (0, 1] for indexed "
            "size properties; 0.8 is always included (default: (0.8,)).",
        ),
    )

    def build(self):
        """Build the mesh, components, phases, and density parameters."""
        super().build()

        # Mesh
        if self.config.size_edges is None:
            raise ConfigurationError(
                "SolidPSDParameterBlock requires a 'size_edges' list."
            )
        bottom = self.config.bottom_size
        try:
            edges = list(size_mesh.validate_mesh(self.config.size_edges, bottom))
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(str(exc)) from exc
        if bottom is not None:
            bottom = float(bottom)
        n_int = len(edges) - 1

        targets = _validate_percentile_targets(self.config.percentile_targets)
        self.percentile_targets = Set(initialize=targets, ordered=True)

        # Component lists
        sized = _validated_name_tuple(
            self.config.sized_solid_component_list,
            "sized_solid_component_list",
            allow_empty=False,
        )
        unsized = _validated_name_tuple(
            self.config.unsized_solid_component_list,
            "unsized_solid_component_list",
            allow_empty=True,
        )
        liquids = self.config.liquid_component_list
        liquids = (
            ("H2O",)
            if liquids is None
            else _validated_name_tuple(
                liquids, "liquid_component_list", allow_empty=False
            )
        )
        vapors = _validated_name_tuple(
            self.config.vapor_component_list,
            "vapor_component_list",
            allow_empty=True,
        )
        overlap = set(sized).intersection(unsized)
        if overlap:
            raise ConfigurationError(
                f"component names {sorted(overlap)!r} appear in both sized and "
                "unsized solid lists."
            )
        union = dict.fromkeys(sized + unsized + liquids + vapors)

        # Densities
        dens_solid = _validated_density_dict(
            self.config.solid_density, sized + unsized, "solid_density"
        )
        dens_liquid = _validated_density_dict(
            self.config.liquid_density, liquids, "liquid_density"
        )

        # Phases and components
        # Restrict the IDAES phase-component set to declared pairs.
        self.Sol = SolidPhase(component_list=list(sized + unsized))
        self.Liq = LiquidPhase(component_list=list(liquids))
        if vapors:
            self.Vap = VaporPhase(component_list=list(vapors))
        for c in union:
            if hasattr(self, c):
                raise ConfigurationError(
                    f"component name {c!r} collides with an existing block "
                    "attribute."
                )
            setattr(self, c, Component())

        self.sized_solid_list = sized
        self.unsized_solid_list = unsized
        self.all_solid_list = sized + unsized
        self.liquid_list = liquids
        self.vapor_list = vapors

        # Mesh parameters
        self.size_interval_set = RangeSet(0, n_int - 1)
        self.size_edge_index = RangeSet(0, n_int)
        self.size_edges = Param(
            self.size_edge_index,
            initialize={i: edges[i] for i in range(n_int + 1)},
            units=units.m,
            mutable=True,
            doc="Strictly increasing size-mesh edges (m).",
        )
        dchars = size_mesh.characteristic_sizes(edges, bottom)
        self.size_char = Param(
            self.size_interval_set,
            initialize={k: dchars[k] for k in range(n_int)},
            units=units.m,
            mutable=True,
            doc="Geometric-mean size per interval (m); uses bottom_size for a zero "
            "lower edge.",
        )
        # Density parameters
        self.dens_mass_solid_comp = Param(
            self.all_solid_list,
            initialize=dens_solid,
            units=units.kg / units.m**3,
            mutable=True,
            doc="Configured density of each solid component (kg/m^3).",
        )
        self.dens_mass_liquid_comp = Param(
            self.liquid_list,
            initialize=dens_liquid,
            units=units.kg / units.m**3,
            mutable=True,
            doc="Configured density of each liquid component (kg/m^3).",
        )

        self._edge_values = tuple(edges)
        self._size_char_values = tuple(dchars)

        self._state_block_class = SolidPSDStateBlock

    def assert_same_mesh(self, other):
        """Compare cached size edges within the mesh tolerance.

        ``bottom_size``, components, and densities are not compared. Raise
        ``ConfigurationError`` if the edges differ.
        """
        if not hasattr(other, "_edge_values"):
            raise ConfigurationError(
                "assert_same_mesh requires another PSD parameter block "
                f"(got {type(other).__name__})."
            )
        if not size_mesh.meshes_equal(self._edge_values, other._edge_values):
            raise ConfigurationError(
                "size-mesh mismatch between parameter blocks:\n"
                f"  this : {self._edge_values}\n"
                f"  other: {other._edge_values}"
            )

    @property
    def size_edges_m(self):
        """Mesh edges in meters as floats, as validated at build."""
        return self._edge_values

    @property
    def size_char_m(self):
        """Characteristic interval sizes in meters as floats, as built."""
        return self._size_char_values

    def size_at_passing_from_flows(self, size_bin_flows, target=0.8, eps=EPS_SMOOTH):
        """Compute a smooth passing size from proposed size-bin mass flows (m).

        ``size_bin_flows`` maps each label to ``N`` kg/s flows, one per size
        interval; flows are summed across labels. The calculation uses the same
        smoothing and fixed mass-flow floor as ``size_at_passing``. The target
        is smoothly limited to top-edge passing. An empty stream yields a finite
        value with no physical meaning.
        """
        n_int = len(self._edge_values) - 1
        size = [sum(row[k] for row in size_bin_flows.values()) for k in range(n_int)]
        total = sum(size)
        flow_mass_eps = value(MASS_FLOW_EPS)
        cum = [
            sum(size[i] for i in range(k + 1)) / (total + flow_mass_eps)
            for k in range(n_int)
        ]
        target = smooth_min(target, cum[-1], eps=eps)
        return size_at_passing_smooth(self._edge_values, cum, target, eps=eps)

    @classmethod
    def define_metadata(cls, obj):
        """Register state properties and their units.

        Every property the state block offers is registered, whether or not a
        unit model reads it. ``flow_mass_solid_comp`` and ``mass_frac_solid_comp`` are
        custom because they cover solids only.
        """
        obj.add_properties(
            {
                "temperature": {"method": None},
                "pressure": {"method": None},
                "flow_mass_phase": {"method": None},
            }
        )
        flow_mass = units.kg / units.s
        flow_vol = units.m**3 / units.s
        density = units.kg / units.m**3
        obj.define_custom_properties(
            {
                # Properties constructed in build()
                "flow_mass_sized_comp_size": {"method": None, "units": flow_mass},
                "flow_mass_unsized_comp": {"method": None, "units": flow_mass},
                "flow_mass_liquid_comp": {"method": None, "units": flow_mass},
                "flow_mass_vapor_comp": {"method": None, "units": flow_mass},
                "flow_mass_size": {"method": None, "units": flow_mass},
                "flow_mass_sized": {"method": None, "units": flow_mass},
                "flow_mass_sized_comp": {"method": None, "units": flow_mass},
                "flow_mass_solid_comp": {"method": None, "units": flow_mass},
                "cum_passing_size": {"method": None, "units": units.dimensionless},
                "percentile_size": {"method": None, "units": units.m},
                "flow_vol_solid": {"method": None, "units": flow_vol},
                "flow_vol_liquid": {"method": None, "units": flow_vol},
                "flow_vol_slurry": {"method": None, "units": flow_vol},
                "vol_frac_solid_slurry": {"method": None, "units": units.dimensionless},
                "dens_mass_solid": {"method": None, "units": density},
                "dens_mass_slurry": {"method": None, "units": density},
                # Properties constructed on first access
                "mass_frac_solid_comp": {
                    "method": "_mass_frac_solid_comp",
                    "units": units.dimensionless,
                },
                "mass_frac_size": {
                    "method": "_mass_frac_size",
                    "units": units.dimensionless,
                },
                "mass_frac_size_comp": {
                    "method": "_mass_frac_size_comp",
                    "units": units.dimensionless,
                },
                "mass_frac_comp_size": {
                    "method": "_mass_frac_comp_size",
                    "units": units.dimensionless,
                },
                "cum_passing_comp_size": {
                    "method": "_cum_passing_comp_size",
                    "units": units.dimensionless,
                },
                "percentile_size_comp": {
                    "method": "_percentile_size_comp",
                    "units": units.m,
                },
                "dens_mass_solid_mass_weighted": {
                    "method": "_dens_mass_solid_mass_weighted",
                    "units": density,
                },
            }
        )
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _SolidPSDStateBlock(StateBlock):
    """Indexed state block for the slurry PSD package."""

    default_initializer = SolidPSDInitializer
    default_scaler = SolidPSDScaler

    def fix_initialization_states(self):
        """Fix all state vars (no consistency constraint exists to deactivate)."""
        fix_state_vars(self)


@declare_process_block_class("SolidPSDStateBlock", block_class=_SolidPSDStateBlock)
class SolidPSDStateBlockData(StateBlockData):
    """State block holding the per-component slurry PSD state."""

    def build(self):
        """Build state variables and main PSD and slurry properties."""
        super().build()
        p = self.params
        sized = p.sized_solid_list
        unsized = p.unsized_solid_list
        sset = p.size_interval_set

        # State variables
        self.flow_mass_sized_comp_size = Var(
            sized,
            sset,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=1e-3,
            units=units.kg / units.s,
            doc="Retained mass flow per sized solid component and size interval (kg/s).",
        )
        if unsized:
            self.flow_mass_unsized_comp = Var(
                unsized,
                domain=NonNegativeReals,
                bounds=(0, None),
                initialize=1e-3,
                units=units.kg / units.s,
                doc="Mass flow of an unsized solid component (kg/s).",
            )
        self.flow_mass_liquid_comp = Var(
            p.liquid_list,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=1e-3,
            units=units.kg / units.s,
            doc="Liquid component mass flow (kg/s).",
        )
        if p.vapor_list:
            self.flow_mass_vapor_comp = Var(
                p.vapor_list,
                domain=NonNegativeReals,
                bounds=(0, None),
                initialize=1e-3,
                units=units.kg / units.s,
                doc="Vapor component mass flow (kg/s).",
            )
        self.temperature = Var(
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=298.15,
            units=units.K,
            doc="State temperature (K).",
        )
        self.pressure = Var(
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=101325.0,
            units=units.Pa,
            doc="State pressure (Pa).",
        )

        # Per-component PSD expressions
        @self.Expression(
            sized, doc="Total mass flow of a sized solid component (kg/s)."
        )
        def flow_mass_sized_comp(b, s):
            return sum(b.flow_mass_sized_comp_size[s, k] for k in sset)

        # Aggregate properties
        @self.Expression(p.phase_list, doc="Total mass flow of each phase (kg/s).")
        def flow_mass_phase(b, phase):
            if phase == "Sol":
                total = sum(b.flow_mass_sized_comp[s] for s in sized)
                if unsized:
                    total = total + sum(b.flow_mass_unsized_comp[u] for u in unsized)
                return total
            if phase == "Liq":
                return sum(b.flow_mass_liquid_comp[j] for j in p.liquid_list)
            return sum(b.flow_mass_vapor_comp[v] for v in p.vapor_list)

        @self.Expression(p.all_solid_list, doc="Per-solid-component mass flow (kg/s).")
        def flow_mass_solid_comp(b, j):
            if j in sized:
                return b.flow_mass_sized_comp[j]
            return b.flow_mass_unsized_comp[j]

        @self.Expression(
            sset, doc="Retained mass flow per interval, sized solids (kg/s)."
        )
        def flow_mass_size(b, k):
            return sum(b.flow_mass_sized_comp_size[s, k] for s in sized)

        @self.Expression(doc="Total sized-solids flow (kg/s).")
        def flow_mass_sized(b):
            return sum(b.flow_mass_size[k] for k in sset)

        @self.Expression(
            sset,
            doc="Cumulative passing fraction of sized solids at each interval "
            "upper edge.",
        )
        def cum_passing_size(b, k):
            return sum(b.flow_mass_size[i] for i in sset if i <= k) / (
                b.flow_mass_sized + MASS_FLOW_EPS
            )

        @self.Expression(
            p.percentile_targets,
            doc="Smooth passing size of sized solids at each configured fraction (m).",
        )
        def percentile_size(b, target):
            return b.size_at_passing(target)

        # Slurry volume and density expressions
        @self.Expression(doc="Volumetric flow of all solids (m^3/s).")
        def flow_vol_solid(b):
            return sum(
                b.flow_mass_solid_comp[j] / b.params.dens_mass_solid_comp[j]
                for j in b.params.all_solid_list
            )

        @self.Expression(doc="Volumetric flow of liquid (m^3/s).")
        def flow_vol_liquid(b):
            return sum(
                b.flow_mass_liquid_comp[j] / b.params.dens_mass_liquid_comp[j]
                for j in b.params.liquid_list
            )

        @self.Expression(doc="Volumetric slurry flow, solid + liquid (m^3/s).")
        def flow_vol_slurry(b):
            return b.flow_vol_solid + b.flow_vol_liquid

        @self.Expression(doc="Solids volume fraction of the slurry (Cv).")
        def vol_frac_solid_slurry(b):
            return b.flow_vol_solid / (b.flow_vol_slurry + VOLUME_FLOW_EPS)

        @self.Expression(
            doc="Solid-mixture density from total mass and volume flows (kg/m^3)."
        )
        def dens_mass_solid(b):
            return b.flow_mass_phase["Sol"] / (b.flow_vol_solid + VOLUME_FLOW_EPS)

        @self.Expression(
            doc="Slurry density over the solid + liquid phases (kg/m^3); "
            "vapor excluded."
        )
        def dens_mass_slurry(b):
            return (b.flow_mass_phase["Sol"] + b.flow_mass_phase["Liq"]) / (
                b.flow_vol_slurry + VOLUME_FLOW_EPS
            )

    # Build-on-demand properties
    def _mass_frac_size_comp(self):
        """Construct retained size fractions for each sized solid component."""
        p = self.params

        @self.Expression(
            p.sized_solid_list,
            p.size_interval_set,
            doc="Fraction of each sized solid component's mass flow retained in an interval.",
        )
        def mass_frac_size_comp(b, s, k):
            return b.flow_mass_sized_comp_size[s, k] / (
                b.flow_mass_sized_comp[s] + MASS_FLOW_EPS
            )

    def _mass_frac_comp_size(self):
        """Construct component mass fractions within each size interval."""
        p = self.params

        @self.Expression(
            p.sized_solid_list,
            p.size_interval_set,
            doc="Fraction of an interval's sized-solids mass flow in each component.",
        )
        def mass_frac_comp_size(b, s, k):
            return b.flow_mass_sized_comp_size[s, k] / (
                b.flow_mass_size[k] + MASS_FLOW_EPS
            )

    def _cum_passing_comp_size(self):
        """Construct cumulative passing fractions for each sized solid component."""
        p = self.params
        sset = p.size_interval_set

        @self.Expression(
            p.sized_solid_list,
            sset,
            doc="Cumulative passing fraction for each sized solid component at the "
            "interval upper edge.",
        )
        def cum_passing_comp_size(b, s, k):
            return sum(b.flow_mass_sized_comp_size[s, i] for i in sset if i <= k) / (
                b.flow_mass_sized_comp[s] + MASS_FLOW_EPS
            )

    def _mass_frac_solid_comp(self):
        """Construct composition fractions on a total-solids basis."""

        @self.Expression(
            self.params.all_solid_list,
            doc="Fraction of total solid mass flow in each solid component.",
        )
        def mass_frac_solid_comp(b, j):
            return b.flow_mass_solid_comp[j] / (
                b.flow_mass_phase["Sol"] + MASS_FLOW_EPS
            )

    def _mass_frac_size(self):
        """Construct retained size fractions on a sized-solids basis."""

        @self.Expression(
            self.params.size_interval_set,
            doc="Fraction of sized-solids mass flow retained in an interval.",
        )
        def mass_frac_size(b, k):
            return b.flow_mass_size[k] / (b.flow_mass_sized + MASS_FLOW_EPS)

    def _percentile_size_comp(self):
        """Construct passing sizes for each sized component and target."""

        @self.Expression(
            self.params.sized_solid_list,
            self.params.percentile_targets,
            doc="Smooth passing size by sized component and configured fraction (m).",
        )
        def percentile_size_comp(b, comp, target):
            return b.size_at_passing(target, comp=comp)

    def _dens_mass_solid_mass_weighted(self):
        """Construct the mass-weighted mean solids density."""

        @self.Expression(doc="Mass-weighted mean solids density (kg/m^3).")
        def dens_mass_solid_mass_weighted(b):
            return sum(
                b.flow_mass_solid_comp[j] * b.params.dens_mass_solid_comp[j]
                for j in b.params.all_solid_list
            ) / (b.flow_mass_phase["Sol"] + MASS_FLOW_EPS)

    # Percentile helpers
    def size_at_passing(self, target, comp=None, eps=EPS_SMOOTH):
        """Return a smooth estimate of the size at ``target`` passing (m).

        ``target`` is the fraction of sized-solid mass passing, in (0, 1];
        0.8 requests the 80% passing size. With ``comp=None``, use the bulk
        sized PSD;
        otherwise use the named sized-solid component. The first component
        call builds ``cum_passing_comp_size``. ``eps`` controls smoothing,
        not the size error. The target is smoothly limited to top-edge passing.
        For an empty sized stream, the result is finite but has no physical
        meaning.
        """
        p = self.params
        edges = [p.size_edges[i] for i in p.size_edge_index]
        if comp is None:
            cum = [self.cum_passing_size[k] for k in p.size_interval_set]
        else:
            cum = [self.cum_passing_comp_size[comp, k] for k in p.size_interval_set]
        target = smooth_min(target, cum[-1], eps=eps)
        return size_at_passing_smooth(edges, cum, target, eps=eps)

    def percentile_size_report(self, target, comp=None, *, rel_tol):
        """Report whether a smooth passing size agrees with linear bin interpolation.

        ``comp=None`` uses all sized solids; otherwise select a sized solid.
        ``target`` is in (0, 1] and need not be configured.
        ``rel_tol`` is the allowed difference relative to the interpolated size.

        Return a dict with ``status`` (``reliable``, ``unreliable``, or
        ``no_particles``), ``smooth_size_m``, ``reference_size_m``, and
        ``relative_error``. ``size_m`` is the smooth size only when reliable;
        the reference and error are None for an empty sized PSD. Reliability
        measures numerical agreement, not physical accuracy.
        """
        target = _cfg_float(target, "target", positive=True, hi=1.0)
        tolerance = _cfg_float(rel_tol, "rel_tol", lo=0.0)
        p = self.params
        if comp is not None and comp not in p.sized_solid_list:
            raise KeyError(f"{comp!r} is not a sized solid component.")
        components = p.sized_solid_list if comp is None else (comp,)
        flows = []
        for k in p.size_interval_set:
            bin_flows = [
                value(
                    units.convert(
                        self.flow_mass_sized_comp_size[s, k], units.kg / units.s
                    ),
                    exception=False,
                )
                for s in components
            ]
            if any(v is None or not math.isfinite(v) or v < 0.0 for v in bin_flows):
                raise ConfigurationError(
                    f"{self.name}: selected sized flows must be set, finite, "
                    "and nonnegative."
                )
            flows.append(sum(bin_flows))
        total = sum(flows)
        if not math.isfinite(total):
            raise ConfigurationError(
                f"{self.name}: selected sized flow total is not finite."
            )
        smooth_size = value(
            units.convert(self.size_at_passing(target, comp=comp), units.m)
        )
        report = {
            "status": "no_particles",
            "size_m": None,
            "smooth_size_m": smooth_size,
            "reference_size_m": None,
            "relative_error": None,
        }
        if total == 0.0:
            return report
        edges = [
            value(units.convert(p.size_edges[i], units.m)) for i in p.size_edge_index
        ]
        reference = size_at_passing_numeric(edges, cumulative_passing(flows), target)
        error = abs(smooth_size - reference) / reference
        reliable = math.isfinite(smooth_size) and error <= tolerance
        report.update(
            status="reliable" if reliable else "unreliable",
            size_m=smooth_size if reliable else None,
            reference_size_m=reference,
            relative_error=error,
        )
        return report

    # Feed validation
    def validate_feed(self):
        """Validate current feed values and return True; fixedness is not checked.

        Raise ``ConfigurationError`` for unset/non-finite or negative mass flows
        and unset/non-finite or non-positive temperature/pressure. The total
        flow across all phases must be positive and finite. Water-only feeds
        are valid.
        """
        flows = [
            value(var, exception=False)
            for var in SolidPSDScaler.iter_flow_var_data(self)
        ]
        for v in flows:
            if v is None or not math.isfinite(v):
                raise ConfigurationError(
                    f"{self.name}: feed contains an unset or non-finite flow value."
                )
            if v < 0.0:
                raise ConfigurationError(
                    f"{self.name}: feed contains a negative flow value."
                )
        for name in ("temperature", "pressure"):
            v = value(getattr(self, name), exception=False)
            if v is None or not math.isfinite(v) or v <= 0.0:
                raise ConfigurationError(
                    f"{self.name}: {name} must be a positive finite value."
                )
        total = sum(flows)
        if not math.isfinite(total) or total <= 0.0:
            raise ConfigurationError(
                f"{self.name}: total feed flow across all phases must be "
                "positive and finite."
            )
        return True

    # Material-flow interface
    def get_material_flow_terms(self, p, j):
        """Return the mass-flow term for a declared phase and component."""
        if p == "Sol":
            return self.flow_mass_solid_comp[j]
        if p == "Liq":
            return self.flow_mass_liquid_comp[j]
        if p == "Vap" and hasattr(self, "flow_mass_vapor_comp"):
            return self.flow_mass_vapor_comp[j]
        raise KeyError(f"phase {p!r} is not declared on this state block")

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def default_material_balance_type(self):
        """Balance each component over its phases (no size resolution)."""
        return MaterialBalanceType.componentTotal

    def define_state_vars(self):
        """Return the independent flow, temperature, and pressure variables."""
        sv = {"flow_mass_sized_comp_size": self.flow_mass_sized_comp_size}
        if hasattr(self, "flow_mass_unsized_comp"):
            sv["flow_mass_unsized_comp"] = self.flow_mass_unsized_comp
        sv["flow_mass_liquid_comp"] = self.flow_mass_liquid_comp
        if hasattr(self, "flow_mass_vapor_comp"):
            sv["flow_mass_vapor_comp"] = self.flow_mass_vapor_comp
        sv["temperature"] = self.temperature
        sv["pressure"] = self.pressure
        return sv


def _validate_percentile_targets(cfg_targets):
    """Return validated passing fractions, appending 0.8 if absent."""
    if not isinstance(cfg_targets, (list, tuple)) or not cfg_targets:
        raise ConfigurationError(
            "percentile_targets must be a non-empty list or tuple."
        )
    targets = tuple(
        _cfg_float(target, f"percentile_targets[{i}]", positive=True, hi=1.0)
        for i, target in enumerate(cfg_targets)
    )
    if len(set(targets)) != len(targets):
        raise ConfigurationError("percentile_targets must not contain duplicates.")
    if 0.8 not in targets:
        targets += (0.8,)
    return targets


def _validated_name_tuple(raw_names, option_name, *, allow_empty):
    """Validate a component-name list (identifier-safe, unique, concrete)."""
    if raw_names is None:
        if allow_empty:
            return ()
        raise ConfigurationError(f"{option_name} is required and must be non-empty.")
    if isinstance(raw_names, (str, bytes)) or not isinstance(raw_names, (list, tuple)):
        raise ConfigurationError(
            f"{option_name} must be a list or tuple of component names "
            f"(got {type(raw_names).__name__})."
        )
    if not raw_names and not allow_empty:
        raise ConfigurationError(f"{option_name} must be non-empty.")
    seen = set()
    for c in raw_names:
        if not isinstance(c, str) or not c.isidentifier():
            raise ConfigurationError(
                f"component name {c!r} in {option_name} must be a valid Python "
                "identifier; rename source labels before configuration."
            )
        if c in _RESERVED_NAMES:
            raise ConfigurationError(
                f"component name {c!r} in {option_name} is reserved (collides with a "
                "phase name or an internal parameter-block attribute)."
            )
        if c in seen:
            raise ConfigurationError(
                f"duplicate component name {c!r} in {option_name}."
            )
        seen.add(c)
    return tuple(raw_names)


def _validated_density_dict(raw_densities, keys, option_name):
    """Validate exact component keys and positive, finite density values."""
    if raw_densities is None or not isinstance(raw_densities, dict):
        raise ConfigurationError(
            f"{option_name} is required as a dict of component densities in kg/m^3."
        )
    missing = [
        component_name for component_name in keys if component_name not in raw_densities
    ]
    extra = [
        component_name for component_name in raw_densities if component_name not in keys
    ]
    if missing or extra:
        raise ConfigurationError(
            f"{option_name} keys must match the declared components exactly "
            f"(missing: {missing}, unexpected: {extra})."
        )
    validated_densities = {}
    for component_name in keys:
        density = raw_densities[component_name]
        if not _is_real_number(density):
            raise ConfigurationError(
                f"{option_name}[{component_name!r}] must be an int or float "
                "(not bool)."
            )
        density = float(density)
        if not math.isfinite(density) or density <= 0.0:
            raise ConfigurationError(
                f"{option_name}[{component_name!r}] must be a positive finite "
                "density (kg/m^3)."
            )
        validated_densities[component_name] = density
    return validated_densities
