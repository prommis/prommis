#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Mineral-by-size particle-size-distribution (PSD) property package.

One joint (size, mineral) distribution per stream, plus a lumped liquid
carrier flow. The state variables include:

* ``flow_mass_size_comp[k, j]`` -- retained mass flow of mineral ``j`` within
  size interval ``k`` (kg/s).
* ``flow_mass_liquid[j]`` -- mass flow of liquid component ``j`` (kg/s).
* ``flow_mass_vapor[j]`` -- mass flow of vapor component ``j`` (kg/s), when
    vapor components are configured.

Unlike :mod:`bulk_psd.bulk_psd`, where composition (``flow_mass_comp``) and PSD
(``flow_mass_size``) are independent vectors reconciled by a consistency
constraint, this package tracks the joint distribution directly.  The bulk
marginals are **derived**, not stated:

* ``flow_mass_comp[j]     = sum_k flow_mass_size_comp[k, j]`` for solids only
* ``flow_mass_size[k]     = sum_j flow_mass_size_comp[k, j]``
* ``grade_by_size[k, j]   = flow_mass_size_comp[k, j] / flow_mass_size[k]``

so marginal consistency is automatic by construction -- there is no
``flow_consistency_eqn`` on this package, on any block, regardless of
``defined_state``.  This is the intended difference from bulk PSD: bulk PSD
represents a unit (e.g. a crusher) that reshapes the size vector while leaving
composition untouched exactly; mineral-by-size represents streams where
composition genuinely varies by size (liberation, differential settling,
size-dependent magnetic/gravity/flotation response) and the two axes cannot be
treated independently.

This is a solid + liquid package by default: ``flow_mass_size_comp`` ranges
over solid minerals only, while liquid components are tracked in a separate
component-indexed state. Vapor components can optionally be configured and
are tracked in their own component-indexed state. PSD and grade diagnostics
remain solids-only, matching the mineral-by-size representation.
``temperature`` and ``pressure`` are carried as ordinary state variables, but this package adds no energy or momentum balance of its own: there is no enthalpy, heat
capacity, or pressure-drop correlation here, so a unit model built on this
package that wants isothermal/isobaric behavior must impose that itself, e.g.
via a ``ControlVolume0DBlockData`` built with
``energy_balance_type=EnergyBalanceType.none`` /
``has_pressure_change=False`` and its own
``properties_out[t].temperature == properties_in[t].temperature`` /
``properties_out[t].pressure == properties_in[t].pressure`` constraints.
This package only supplies ``temperature``/``pressure`` as free state vars
for such a constraint to reference; it does not add that constraint itself.
Size-mesh handling
(:mod:`psd.core.size_mesh`), percentile math
(:mod:`psd.core.psd_math`), and scaler read idioms
(:mod:`psd.core.unit_utils`) are reused directly from ``psd.core``,
which is documented as representation-independent.

Derived mass fractions, grades, and percentiles use a unit-bearing denominator
floor and are diagnostic-only at zero flow. ``get_material_flow_terms``
exposes the derived component-total flow, so unit models needing size
resolution (e.g. a WHIMS separator) must read ``flow_mass_size_comp``
directly rather than going through the generic material-flow-terms hook.
"""

__author__ = "Carolina Tristan"

import math

from pyomo.common.config import ConfigValue
from pyomo.environ import (
    NonNegativeReals,
    Param,
    RangeSet,
    Reals,
    Set,
    Var,
    units,
    value,
)

from idaes.core import (
    Component,
    LiquidPhase,
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

from prommis.psd.core import size_mesh
from prommis.psd.core.psd_math import size_at_passing_smooth
from prommis.psd.core.smooth_functions import EPS_SMOOTH
from prommis.psd.core.unit_utils import (
    FactorSource,
    declare_factor_source,
    fixed_value_or_none,
)

# Absolute numeric threshold (kg/s) used by validate_feed to detect a
# degenerate zero-flow feed.
_EPS_NUMERIC = 1e-12

# Reject names that would be overwritten by parameter-block attributes created
# after component declaration, mirroring bulk_psd's reserved-name gate plus
# this package's own derived-Expression names.
_RESERVED_NAMES = (
    "Sol",
    "Liq",
    "H2O",
    "liquid_component_set",
    "vapor_component_set",
    "solid",
    "component_list",
    "solid_component_set",
    "size_interval_set",
    "size_edge_index",
    "size_edges",
    "d_char",
    "flow_eps",
    "assert_same_mesh",
    "_edge_values",
    "_bottom_size",
    "flow_mass_size_comp",
    "flow_mass_comp",
    "flow_mass_size",
    "flow_mass",
    "mass_frac_comp",
    "mass_frac_size",
    "grade_by_size",
    "cum_passing",
    "cum_passing_mineral",
    "P80",
    "P50",
    "P80_mineral",
    "temperature",
    "pressure",
)


class MineralSizePSDScaler(CustomScalerBase):
    """Per-(size, mineral)-cell scaler for the mineral-by-size PSD state.

    The single joint flow Var uses inverse-magnitude scaling with a
    stream-relative floor, which avoids extreme factors on empty/near-empty (size, mineral) cells -- these are common in liberation data, where most of a coarse bin's mass sits in one or two host minerals and everything else is exactly or near zero.

    There is no consistency constraint in this package (marginals are derived
    Expressions, not Vars reconciled by a Constraint), so
    ``constraint_scaling_routine`` has nothing package-specific to do; it is
    a no-op provided for CustomScalerBase interface completeness and so a
    calling unit scaler can invoke it unconditionally without a hasattr
    check.
    """

    CONFIG = declare_factor_source(CustomScalerBase.CONFIG())

    EPS_REL = 1e-8

    #: Documented coarse default for an unfixed flow Var scaled standalone in
    #: baseline mode (a 1 kg/s-order nominal in the Var's declared units).
    BASELINE_DEFAULT_FLOW_FACTOR = 1.0

    #: Fixed factors for the pass-through state vars (same constants as
    #: ``SulfuricAcidLeachingPropertiesScaler``), independent of factor_source
    #: since these aren't magnitude-scaled off the joint flow state.
    TEMPERATURE_SCALING_FACTOR = 1 / 300
    PRESSURE_SCALING_FACTOR = 1e-5
    BASELINE_DEFAULT_LIQUID_FACTOR = 1.0
    BASELINE_DEFAULT_VAPOR_FACTOR = 1.0

    @staticmethod
    def _magnitude(var):
        v = value(var, exception=False)
        return abs(v) if v is not None else 0.0

    def _stream_total(self, model, magnitude=None):
        magnitude = magnitude or self._magnitude
        return (
            sum(
                magnitude(model.flow_mass_size_comp[k, j])
                for (k, j) in model.flow_mass_size_comp
            )
            or 1.0
        )

    @staticmethod
    def _fixed_magnitude(var):
        v = fixed_value_or_none(var)
        return abs(v) if v is not None else 0.0

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if self.config.factor_source == FactorSource.baseline:
            floor = self.EPS_REL * self._stream_total(model, self._fixed_magnitude)
            for var in model.flow_mass_size_comp.values():
                v = fixed_value_or_none(var)
                sf = (
                    1.0 / max(abs(v), floor)
                    if v is not None
                    else self.BASELINE_DEFAULT_FLOW_FACTOR
                )
                self.set_variable_scaling_factor(var, sf, overwrite=overwrite)
        else:
            floor = self.EPS_REL * self._stream_total(model)
            for var in model.flow_mass_size_comp.values():
                mag = max(self._magnitude(var), floor)
                self.set_variable_scaling_factor(var, 1.0 / mag, overwrite=overwrite)

        flow_floor = self.EPS_REL * self._stream_total(model)
        for flow_var, default_factor in (
            (model.flow_mass_liquid, self.BASELINE_DEFAULT_LIQUID_FACTOR),
            (
                model.flow_mass_vapor if hasattr(model, "flow_mass_vapor") else (),
                self.BASELINE_DEFAULT_VAPOR_FACTOR,
            ),
        ):
            for var in getattr(flow_var, "values", lambda: ())():
                fixed_flow = fixed_value_or_none(var)
                if self.config.factor_source == FactorSource.baseline:
                    factor = (
                        1.0 / max(abs(fixed_flow), flow_floor)
                        if fixed_flow is not None
                        else default_factor
                    )
                else:
                    factor = 1.0 / max(self._magnitude(var), flow_floor)
                self.set_variable_scaling_factor(var, factor, overwrite=overwrite)

        self.set_variable_scaling_factor(
            model.temperature, self.TEMPERATURE_SCALING_FACTOR, overwrite=overwrite
        )
        self.set_variable_scaling_factor(
            model.pressure, self.PRESSURE_SCALING_FACTOR, overwrite=overwrite
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict | None = None
    ):
        # No package-declared Constraints to scale (see class docstring).
        # A downstream unit model (e.g. WHIMSSeparator) that adds its own
        # recovery_eqn / mags_solid_balance Constraints owns scaling those
        # itself, using `required_scaling_factor(self, ...)` against the
        # factors this routine sets on flow_mass_size_comp.
        pass


class MineralSizePSDInitializer(InitializerBase):
    """Initializer for the mineral-by-size PSD state block.

    No solve is needed. Because marginals are derived (no consistency
    constraint to satisfy), validation reduces to checking the fixed
    (or set) joint-flow values are finite and non-negative with a nonzero
    total -- simpler than ``BulkPSDInitializer``, which must also reconcile
    two independently-fixed vectors.
    """

    CONFIG = InitializerBase.CONFIG()

    def initialization_routine(self, model):
        for idx in model:
            sbd = model[idx]
            if sbd.config.defined_state:
                sbd.validate_feed()
        return None


@declare_process_block_class("MineralSizePSDParameterBlock")
class MineralSizePSDParameterData(PhysicalParameterBlock):
    """Parameter block for the mineral-by-size PSD package.

    Configuration:
        size_edges: ascending list of ``N+1`` size edges in meters (required).
        bottom_size: positive replacement for a zero finest edge; required iff
            ``size_edges[0] == 0``.
        component_list: non-empty list of unique, identifier-safe mineral
            component names (required).
        liquid_component_list: liquid species names; defaults to ``["H2O"]``.
        vapor_component_list: optional vapor species names; an omitted or
            empty list disables the vapor phase. Species may occur in both
            fluid phases.

    Connected blocks must share **the same parameter block instance** -- mesh
    identity is by object; ``assert_same_mesh`` is provided for translator/
    diagnostic use, exactly as in ``BulkPSDParameterBlock``.

    The mesh Params (``size_edges``, ``d_char``) are fixed at build time, for
    the same reasons documented on ``BulkPSDParameterBlock``: mutation after
    building state/unit blocks is unsupported and will silently desynchronize
    cached unit-model data from the live state-block PSD.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "size_edges",
        ConfigValue(
            default=None,
            description="Ascending list of N+1 size-mesh edges in meters.",
        ),
    )
    CONFIG.declare(
        "bottom_size",
        ConfigValue(
            default=None,
            description="Positive bottom edge for characteristic sizes when "
            "the finest edge is zero (meters).",
        ),
    )
    CONFIG.declare(
        "component_list",
        ConfigValue(
            default=None,
            description="Non-empty list of identifier-safe mineral/solid "
            "component names.",
        ),
    )
    CONFIG.declare(
        "liquid_component_list",
        ConfigValue(
            default=None,
            description="Component names in the liquid phase; defaults " "to ['H2O'].",
        ),
    )
    CONFIG.declare(
        "vapor_component_list",
        ConfigValue(
            default=None,
            description="Optional component names in the vapor phase. "
            "An empty or omitted list disables the vapor phase.",
        ),
    )

    def build(self):
        super().build()

        edges, bottom, n_int = size_mesh.validate_property_mesh(
            self.config,
            "MineralSizePSDParameterBlock requires a 'size_edges' list.",
        )

        comps = self.config.component_list
        if isinstance(comps, (str, bytes)) or not isinstance(comps, (list, tuple)):
            raise ConfigurationError(
                "component_list must be a list or tuple of component names (got "
                f"{type(comps).__name__})."
            )
        if not comps:
            raise ConfigurationError(
                "MineralSizePSDParameterBlock requires a non-empty "
                "'component_list' (an empty list would leave "
                "flow_mass_size_comp with no second index)."
            )
        seen = set()
        for c in comps:
            if not isinstance(c, str) or not c.isidentifier():
                raise ConfigurationError(
                    f"component name {c!r} is not a valid Python identifier; map "
                    "raw artifact labels (e.g. 'Monazite(s)') to canonical IDs."
                )
            if c in _RESERVED_NAMES:
                raise ConfigurationError(
                    f"component name {c!r} is reserved (collides with the solid "
                    "phase name or an internal parameter-block attribute)."
                )
            if c in seen:
                raise ConfigurationError(f"duplicate component name {c!r}")
            seen.add(c)

        liquid_comps = self.config.liquid_component_list
        if liquid_comps is None:
            liquid_comps = ["H2O"]
        vapor_comps = self.config.vapor_component_list
        if vapor_comps is None:
            vapor_comps = []
        for label, phase_comps in (
            ("liquid_component_list", liquid_comps),
            ("vapor_component_list", vapor_comps),
        ):
            if isinstance(phase_comps, (str, bytes)) or not isinstance(
                phase_comps, (list, tuple)
            ):
                raise ConfigurationError(
                    f"{label} must be a list or tuple of component names."
                )
            for c in phase_comps:
                if not isinstance(c, str) or not c.isidentifier():
                    raise ConfigurationError(
                        f"component name {c!r} in {label} is not a valid "
                        "Python identifier."
                    )
                if c in _RESERVED_NAMES and c != "H2O":
                    raise ConfigurationError(
                        f"component name {c!r} in {label} is reserved."
                    )
            if len(set(phase_comps)) != len(phase_comps):
                raise ConfigurationError(f"{label} contains duplicate names.")
        if not liquid_comps:
            raise ConfigurationError("liquid_component_list cannot be empty.")

        # Each component belongs to exactly one phase; only solids are
        # represented in the size-resolved PSD state.
        self.Sol = SolidPhase()
        self.Liq = LiquidPhase()
        if vapor_comps:
            self.Vap = VaporPhase()
        self.solid_component_set = Set(initialize=comps, ordered=True)
        self.liquid_component_set = Set(initialize=liquid_comps, ordered=True)
        self.vapor_component_set = Set(initialize=vapor_comps, ordered=True)
        all_components = dict.fromkeys(
            list(comps) + list(liquid_comps) + list(vapor_comps)
        )
        for c in all_components:
            setattr(self, c, Component())
        self.phase_component_set = Set(
            initialize=(
                [("Sol", c) for c in comps]
                + [("Liq", c) for c in liquid_comps]
                + [("Vap", c) for c in vapor_comps]
            ),
            doc="Valid phase-component pairs.",
        )

        self.size_interval_set = RangeSet(0, n_int - 1)
        self.size_edge_index = RangeSet(0, n_int)
        self.size_edges = Param(
            self.size_edge_index,
            initialize={i: edges[i] for i in range(n_int + 1)},
            units=units.m,
            mutable=True,
            doc="Size-mesh edges (ascending, meters).",
        )
        dchars = size_mesh.characteristic_sizes(edges, bottom)
        self.d_char = Param(
            self.size_interval_set,
            initialize={k: dchars[k] for k in range(n_int)},
            units=units.m,
            mutable=True,
            doc="Geometric-mean characteristic size of each interval (meters).",
        )
        self.flow_eps = Param(
            initialize=1e-12,
            units=units.kg / units.s,
            mutable=True,
            doc="Additive unit-bearing floor on derived-Expression denominators.",
        )

        self._edge_values = tuple(edges)
        self._bottom_size = bottom

        self._state_block_class = MineralSizePSDStateBlock

    def assert_same_mesh(self, other):
        """Raise unless ``other`` carries an identical size mesh.

        ``other`` must be a MineralSizePSDParameterBlock or a MineralSizePSDStateBlock with a ``params`` attribute pointing to one.

        Raises ``ConfigurationError`` if the meshes differ in length or any edge value, rather than silently returning a boolean.  This is the same behavior as ``BulkPSDParameterBlock.assert_same_mesh`` and is intended to be used in translator/diagnostic code, not in a unit model's normal execution path.
        """
        size_mesh.assert_same_mesh(self, other, "mineral-by-size PSD parameter block")

    @classmethod
    def define_metadata(cls, obj):
        # These flow states do not use the standard (phase, component) shape.
        obj.define_custom_properties(
            {
                "flow_mass_size_comp": {"method": None},
                "flow_mass_liquid": {"method": None},
                "flow_mass_vapor": {"method": None},
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


class _MineralSizePSDStateBlock(StateBlock):
    default_initializer = MineralSizePSDInitializer
    default_scaler = MineralSizePSDScaler

    def fix_initialization_states(self):
        """Fix the joint state var.

        No consistency constraint exists in this package (see module
        docstring), so unlike ``_BulkPSDStateBlock.fix_initialization_states``
        there is nothing to deactivate: fixing every ``flow_mass_size_comp``
        entry alone drives every derived Expression (marginals, grades,
        percentiles) to a fully determined value with no risk of an
        over-determined system.
        """
        fix_state_vars(self)


@declare_process_block_class(
    "MineralSizePSDStateBlock", block_class=_MineralSizePSDStateBlock
)
class MineralSizePSDStateBlockData(StateBlockData):
    """State block holding the joint (size, mineral) mass-flow distribution."""

    default_scaler = MineralSizePSDScaler

    def build(self):
        super().build()
        comps = self.params.solid_component_set
        sset = self.params.size_interval_set

        self.flow_mass_size_comp = Var(
            sset,
            comps,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=1e-3,
            units=units.kg / units.s,
            doc="Retained mass flow of each mineral within each size "
            "interval (kg/s) -- the joint mineral-by-size state.",
        )
        self.flow_mass_liquid = Var(
            self.params.liquid_component_set,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=10.0,
            units=units.kg / units.s,
            doc="Mass flow of each liquid component [kg/s].",
        )
        if len(self.params.vapor_component_set):
            self.flow_mass_vapor = Var(
                self.params.vapor_component_set,
                domain=NonNegativeReals,
                bounds=(0, None),
                initialize=0.0,
                units=units.kg / units.s,
                doc="Mass flow of each vapor component [kg/s].",
            )

        self.temperature = Var(
            domain=Reals,
            initialize=298.15,
            bounds=(298.1, None),
            doc="State temperature (pass-through: no energy balance is "
            "declared by this package) [K]",
            units=units.K,
        )
        self.pressure = Var(
            domain=Reals,
            initialize=101325.0,
            bounds=(1e3, 1e6),
            doc="State pressure (pass-through: no momentum balance is "
            "declared by this package) [Pa]",
            units=units.Pa,
        )

        # ------------------------------------------------------------------
        # Derived marginals and grades. All Expressions: consistency with
        # the joint state is automatic, not enforced by a Constraint.
        @self.Expression(comps, doc="Bulk component (mineral) mass flow.")
        def flow_mass_comp(b, j):
            return sum(b.flow_mass_size_comp[k, j] for k in sset)

        @self.Expression(sset, doc="Bulk retained mass flow per size interval.")
        def flow_mass_size(b, k):
            return sum(b.flow_mass_size_comp[k, j] for j in comps)

        @self.Expression(doc="Total solids flow (kg/s).")
        def flow_mass(b):
            return sum(b.flow_mass_size_comp[k, j] for k in sset for j in comps)

        @self.Expression(comps, doc="Bulk component (mineral) mass fraction.")
        def mass_frac_comp(b, j):
            return b.flow_mass_comp[j] / (b.flow_mass + b.params.flow_eps)

        @self.Expression(sset, doc="Bulk retained mass fraction per interval.")
        def mass_frac_size(b, k):
            return b.flow_mass_size[k] / (b.flow_mass + b.params.flow_eps)

        @self.Expression(
            sset,
            comps,
            doc="Grade-by-size: fraction of size interval k's mass that is "
            "mineral j (each interval's grades sum to 1).",
        )
        def grade_by_size(b, k, j):
            return b.flow_mass_size_comp[k, j] / (
                b.flow_mass_size[k] + b.params.flow_eps
            )

        @self.Expression(sset, doc="Cumulative passing at each interval upper edge.")
        def cum_passing(b, k):
            return sum(b.flow_mass_size[i] for i in sset if i <= k) / (
                b.flow_mass + b.params.flow_eps
            )

        @self.Expression(
            sset,
            comps,
            doc="Cumulative passing by interval and mineral.",
        )
        def cum_passing_mineral(b, k, j):
            return sum(b.flow_mass_size_comp[i, j] for i in sset if i <= k) / (
                b.flow_mass_comp[j] + b.params.flow_eps
            )

        @self.Expression(doc="80% passing size (smooth interpolation, meters).")
        def P80(b):
            return b.Pxx(0.8)

        @self.Expression(doc="50% passing size (smooth interpolation, meters).")
        def P50(b):
            return b.Pxx(0.5)

        @self.Expression(
            comps,
            doc="80% passing size for each mineral (meters).",
        )
        def P80_mineral(b, j):
            return b.mineral_Pxx(j, 0.8)

        @self.Expression(
            comps,
            doc="50% passing size for each mineral (meters).",
        )
        def P50_mineral(b, j):
            return b.mineral_Pxx(j, 0.5)

    def Pxx(self, target, eps=EPS_SMOOTH):
        """Smooth linear-in-size percentile of the block's bulk PSD (meters).

        Uses the bulk-rolled-up ``cum_passing`` (summed over minerals), i.e.
        this is the size-only Pxx of the stream, not a per-mineral Pxx.  For
        a per-mineral percentile, build the same
        ``size_at_passing_smooth(edges, cum_j, target)`` call against a
        mineral-specific cumulative curve constructed from
        ``flow_mass_size_comp[:, j] / flow_mass_comp[j]``.
        """
        edges = [self.params.size_edges[i] for i in self.params.size_edge_index]
        cum = [self.cum_passing[k] for k in self.params.size_interval_set]
        return size_at_passing_smooth(edges, cum, target, eps=eps)

    def mineral_Pxx(self, j, target, eps=EPS_SMOOTH):
        """Smooth linear-in-size percentile for mineral ``j`` (meters)."""
        if j not in self.params.solid_component_set:
            raise KeyError(f"Unknown solid mineral {j!r}.")
        edges = [self.params.size_edges[i] for i in self.params.size_edge_index]
        cum = [self.cum_passing_mineral[k, j] for k in self.params.size_interval_set]
        return size_at_passing_smooth(edges, cum, target, eps=eps)

    def validate_feed(self):
        """Validate a fixed feed specification.

        Returns ``True`` for finite, non-negative, nonzero-total joint flow
        values. Raises ``ConfigurationError`` otherwise. There is no
        cross-vector reconciliation to check (unlike
        ``BulkPSDStateBlockData.validate_feed``) because marginals are
        derived from the single joint state, not independently specified.
        """
        vals = [
            value(self.flow_mass_size_comp[k, j], exception=False)
            for k in self.params.size_interval_set
            for j in self.params.solid_component_set
        ]

        for v in vals:
            if v is None or not math.isfinite(v):
                raise ConfigurationError(
                    f"{self.name}: feed contains a non-finite flow value."
                )
        for v in vals:
            if v < 0.0:
                raise ConfigurationError(
                    f"{self.name}: feed contains a negative flow value."
                )

        if sum(vals) <= _EPS_NUMERIC:
            raise ConfigurationError(
                f"{self.name}: degenerate feed -- total solids flow is zero."
            )
        return True

    def get_material_flow_terms(self, p, j):
        phase_name = p if isinstance(p, str) else p.local_name
        if phase_name == "Sol" and j in self.params.solid_component_set:
            return self.flow_mass_comp[j]
        if phase_name == "Liq" and j in self.params.liquid_component_set:
            return self.flow_mass_liquid[j]
        if (
            phase_name == "Vap"
            and hasattr(self, "flow_mass_vapor")
            and j in self.params.vapor_component_set
        ):
            return self.flow_mass_vapor[j]
        raise KeyError(f"Unrecognised (phase, component) pair ({p}, {j})")

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def default_material_balance_type(self):
        from idaes.core import MaterialBalanceType

        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        from idaes.core import EnergyBalanceType

        return EnergyBalanceType.none

    def define_state_vars(self):
        return {
            "flow_mass_size_comp": self.flow_mass_size_comp,
            "flow_mass_liquid": self.flow_mass_liquid,
            "temperature": self.temperature,
            "pressure": self.pressure,
            **(
                {"flow_mass_vapor": self.flow_mass_vapor}
                if hasattr(self, "flow_mass_vapor")
                else {}
            ),
        }
