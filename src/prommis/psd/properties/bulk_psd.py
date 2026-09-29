#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Bulk particle-size-distribution (PSD) property package.

One bulk PSD per stream plus a separate bulk composition vector.  The state is:

* ``flow_mass_comp[j]`` -- mass flow of each solid component ``j`` (kg/s),
* ``flow_mass_size[k]`` -- retained mass flow in each size interval ``k`` (kg/s),

with a consistency constraint ``sum_j flow_mass_comp == sum_k flow_mass_size``
built only on outlet/product blocks (``defined_state=False``).  Composition and
PSD are independent vectors that happen to share a total; a unit that transforms
the size vector leaves composition untouched (exact for a crusher).

This is a **solids-only** package: ``flow_mass_comp`` ranges over solid
components only, with no liquid phase and no temperature/pressure/energy state.

Derived mass fractions and percentiles use a unit-bearing denominator floor and
are diagnostic-only at zero flow. ``get_material_flow_terms`` exposes component
flows, so unit models must add explicit equations for the independent size
vector.
"""

__author__ = "Daison Yancy Caballero"

import math

from pyomo.common.config import ConfigValue
from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Param,
    RangeSet,
    Var,
    units,
    value,
)

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.initialization import InitializerBase
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import fix_state_vars

from prommis.psd.core import size_mesh
from prommis.psd.core.psd_math import size_at_passing_smooth
from prommis.psd.core.unit_utils import (
    FactorSource,
    declare_factor_source,
    fixed_value_or_none,
    required_scaling_factor,
)
from prommis.psd.core.smooth_functions import EPS_SMOOTH

# Absolute numeric threshold (kg/s) used by validate_feed to detect a degenerate
# zero-flow feed and to floor the consistency reference total.
_EPS_NUMERIC = 1e-12

# Reject names that would be overwritten by parameter-block attributes created
# after component declaration.
_RESERVED_NAMES = (
    "solid",
    "component_list",
    "size_interval_set",
    "size_edge_index",
    "size_edges",
    "d_char",
    "flow_eps",
    "assert_same_mesh",
    "_edge_values",
    "_bottom_size",
)


class BulkPSDScaler(CustomScalerBase):
    """Per-size-class scaler for the bulk PSD state.

    Flow variables use inverse-magnitude scaling with a stream-relative floor,
    which avoids extreme factors on empty bins. The consistency constraint uses
    the inverse stream total.

    Under ``factor_source="baseline"`` the flow reads are gated by
    ``var.fixed``: a fixed (feed-spec) Var keeps the identical inverse-magnitude
    arithmetic, while an unfixed Var gets the documented coarse default 1.0 --
    the owning unit scaler propagates real factors onto outlet blocks instead.

    ``constraint_scaling_routine`` derives its factors from the declared
    variable scaling factors and raises ``ConfigurationError`` when a required
    factor is missing -- run ``variable_scaling_routine`` (or set factors
    manually) first.
    """

    CONFIG = declare_factor_source(CustomScalerBase.CONFIG())

    EPS_REL = 1e-8

    #: Documented coarse default for an unfixed flow Var scaled standalone in
    #: baseline mode (a 1 kg/s-order nominal in the Var's declared units).
    BASELINE_DEFAULT_FLOW_FACTOR = 1.0

    @staticmethod
    def _magnitude(var):
        # tolerate uninitialized values (the MineralSizePSDScaler read idiom)
        v = value(var, exception=False)
        return abs(v) if v is not None else 0.0

    def _stream_total(self, model, magnitude=None):
        magnitude = magnitude or self._magnitude
        comp_total = sum(
            magnitude(model.flow_mass_comp[j]) for j in model.flow_mass_comp
        )
        size_total = sum(
            magnitude(model.flow_mass_size[k]) for k in model.flow_mass_size
        )
        ref = max(comp_total, size_total)
        return ref if ref > 0.0 else 1.0

    @staticmethod
    def _fixed_magnitude(var):
        v = fixed_value_or_none(var)
        return abs(v) if v is not None else 0.0

    def _flow_vars(self, model):
        yield from model.flow_mass_comp.values()
        yield from model.flow_mass_size.values()

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if self.config.factor_source == FactorSource.baseline:
            # same arithmetic as seeded, but only fixed (feed-spec) values are
            # data; unfixed Vars take the documented coarse default
            floor = self.EPS_REL * self._stream_total(model, self._fixed_magnitude)
            for var in self._flow_vars(model):
                v = fixed_value_or_none(var)
                sf = (
                    1.0 / max(abs(v), floor)
                    if v is not None
                    else self.BASELINE_DEFAULT_FLOW_FACTOR
                )
                self.set_variable_scaling_factor(var, sf, overwrite=overwrite)
        else:
            floor = self.EPS_REL * self._stream_total(model)
            for var in self._flow_vars(model):
                mag = max(self._magnitude(var), floor)
                self.set_variable_scaling_factor(var, 1.0 / mag, overwrite=overwrite)

    def _side_nominal(self, var):
        """Summed nominal of one consistency-row side from declared factors."""
        return sum(1.0 / required_scaling_factor(self, v) for v in var.values())

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if hasattr(model, "flow_consistency_eqn"):
            # max-of-the-two-sides basis (the _stream_total rule), derived from
            # the declared flow-Var factors set one call earlier.
            nominal = max(
                self._side_nominal(model.flow_mass_comp),
                self._side_nominal(model.flow_mass_size),
            )
            self.set_constraint_scaling_factor(
                model.flow_consistency_eqn,
                1.0 / nominal,
                overwrite=overwrite,
            )


class BulkPSDInitializer(InitializerBase):
    """Initializer for the bulk PSD state block.

    No solve is needed. Fully specified blocks are validated while the base
    initializer manages state fixing and restoration.
    """

    CONFIG = InitializerBase.CONFIG()

    def initialization_routine(self, model):
        # ``model`` is the time-indexed StateBlock; validate every defined-state
        # member.
        for idx in model:
            sbd = model[idx]
            if sbd.config.defined_state:
                sbd.validate_feed()
        return None


@declare_process_block_class("BulkPSDParameterBlock")
class BulkPSDParameterData(PhysicalParameterBlock):
    """Parameter block for the bulk PSD package.

    Configuration:
        size_edges: ascending list of ``N+1`` size edges in meters (required).
        bottom_size: positive replacement for a zero finest edge; required iff
            ``size_edges[0] == 0``.
        component_list: non-empty list of unique, identifier-safe solid component
            names (required).

    Connected blocks must share **the same parameter block instance** -- mesh
    identity is by object; ``assert_same_mesh`` is provided for translator/
    diagnostic use.

    The mesh Params (``size_edges``, ``d_char``) are fixed at build time.  They
    are declared ``mutable=True`` only because Pyomo requires unit-bearing Params
    to be mutable, not as a signal that mutation is supported.  Mutating them
    after building state/unit blocks is **unsupported**: it silently desynchronizes
    cached unit-model data (selection/breakage matrices, distribution anchor
    bounds) from the live state-block P80.
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
            description="Non-empty list of identifier-safe solid component names.",
        ),
    )

    def build(self):
        super().build()

        edges, bottom, n_int = size_mesh.validate_property_mesh(
            self.config,
            "BulkPSDParameterBlock requires a 'size_edges' list.",
        )

        comps = self.config.component_list
        # require a concrete, reusable, ordered sequence: a bare string would be
        # iterated into single-character components and a one-shot iterator would be
        # exhausted by validation, leaving zero components created downstream.
        if isinstance(comps, (str, bytes)) or not isinstance(comps, (list, tuple)):
            raise ConfigurationError(
                "component_list must be a list or tuple of component names (got "
                f"{type(comps).__name__})."
            )
        if not comps:
            raise ConfigurationError(
                "BulkPSDParameterBlock requires a non-empty 'component_list' "
                "(an empty list would force the whole PSD to zero through the "
                "consistency constraint)."
            )
        seen = set()
        for c in comps:
            if not isinstance(c, str) or not c.isidentifier():
                raise ConfigurationError(
                    f"component name {c!r} is not a valid Python identifier; map "
                    "raw artifact labels (e.g. 'Ore1(s)') to canonical IDs."
                )
            if c in _RESERVED_NAMES:
                raise ConfigurationError(
                    f"component name {c!r} is reserved (collides with the solid "
                    "phase name or an internal parameter-block attribute)."
                )
            if c in seen:
                raise ConfigurationError(f"duplicate component name {c!r}")
            seen.add(c)

        # Single solid phase; one Component per requested name (populates
        # component_list in declaration order).
        self.solid = Phase()
        for c in comps:
            setattr(self, c, Component())

        # Ordered index sets and pinned mesh Params (read by the unit models and
        # the function-library builders).
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

        # Stash plain-float mesh data for the numeric helpers / assert_same_mesh.
        self._edge_values = tuple(edges)
        self._bottom_size = bottom

        self._state_block_class = BulkPSDStateBlock

    def assert_same_mesh(self, other):
        """Raise unless ``other`` carries an identical size mesh.

        ``other`` must be a bulk-PSD parameter block (it is matched by its mesh,
        an object attribute), a wrong-type argument raises a clear
        ConfigurationError rather than an incidental ``AttributeError``.
        """
        size_mesh.assert_same_mesh(self, other, "bulk-PSD parameter block")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mass_comp": {"method": None},
            }
        )
        obj.define_custom_properties(
            {
                "flow_mass_size": {"method": None},
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


class _BulkPSDStateBlock(StateBlock):
    default_initializer = BulkPSDInitializer
    default_scaler = BulkPSDScaler

    def fix_initialization_states(self):
        """Fix all state vars; deactivate the consistency constraint while fixed.

        On a ``defined_state=False`` block the consistency constraint would leave
        DOF = -1 once all state vars are fixed, so it is deactivated for the
        duration of initialization (restored by the Initializer's state-restore).
        """
        fix_state_vars(self)
        for sbd in self.values():
            if not sbd.config.defined_state and hasattr(sbd, "flow_consistency_eqn"):
                sbd.flow_consistency_eqn.deactivate()


@declare_process_block_class("BulkPSDStateBlock", block_class=_BulkPSDStateBlock)
class BulkPSDStateBlockData(StateBlockData):
    """State block holding solid composition and retained PSD mass-flow vectors."""

    default_scaler = BulkPSDScaler

    def build(self):
        super().build()
        comps = self.params.component_list
        sset = self.params.size_interval_set

        self.flow_mass_comp = Var(
            comps,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=1e-3,
            units=units.kg / units.s,
            doc="Component mass flow (kg/s).",
        )
        self.flow_mass_size = Var(
            sset,
            domain=NonNegativeReals,
            bounds=(0, None),
            initialize=1e-3,
            units=units.kg / units.s,
            doc="Retained mass flow per size interval (kg/s).",
        )

        if not self.config.defined_state:
            self.flow_consistency_eqn = Constraint(
                expr=sum(self.flow_mass_comp[j] for j in comps)
                == sum(self.flow_mass_size[k] for k in sset),
                doc="Composition total equals PSD total.",
            )

        # ------------------------------------------------------------------
        # Derived expressions.
        @self.Expression(doc="Total solids flow from composition (kg/s).")
        def flow_mass(b):
            return sum(b.flow_mass_comp[j] for j in comps)

        @self.Expression(doc="Total solids flow from PSD (kg/s).")
        def flow_mass_sized(b):
            return sum(b.flow_mass_size[k] for k in sset)

        @self.Expression(comps, doc="Component mass fraction.")
        def mass_frac_comp(b, j):
            return b.flow_mass_comp[j] / (b.flow_mass + b.params.flow_eps)

        @self.Expression(sset, doc="Retained mass fraction per interval.")
        def mass_frac_size(b, k):
            return b.flow_mass_size[k] / (b.flow_mass_sized + b.params.flow_eps)

        @self.Expression(sset, doc="Cumulative passing at each interval upper edge.")
        def cum_passing(b, k):
            return sum(b.flow_mass_size[i] for i in sset if i <= k) / (
                b.flow_mass_sized + b.params.flow_eps
            )

        @self.Expression(doc="80% passing size (smooth interpolation, meters).")
        def P80(b):
            return b.Pxx(0.8)

        @self.Expression(doc="50% passing size (smooth interpolation, meters).")
        def P50(b):
            return b.Pxx(0.5)

    def Pxx(self, target, eps=EPS_SMOOTH):
        """Smooth linear-in-size percentile of the block's PSD (meters).

        Diagnostic-only on degenerate states (the smooth form returns a finite
        but physically meaningless value when the PSD is near-zero).
        """
        edges = [self.params.size_edges[i] for i in self.params.size_edge_index]
        cum = [self.cum_passing[k] for k in self.params.size_interval_set]
        return size_at_passing_smooth(edges, cum, target, eps=eps)

    def validate_feed(self, tol_rel=1e-6):
        """Validate a fixed feed specification.

        Returns ``True`` for finite, non-negative, nonzero composition and PSD
        totals that agree within ``tol_rel``.  Raises ``ConfigurationError`` for
        non-finite values, negative flows, zero total solids flow, or a
        composition/PSD total mismatch.

        The default tolerance admits totals reconstructed from rounded fraction
        data while rejecting material composition/PSD mismatches.
        """
        # exception=False: uninitialized Vars come back None and hit the check
        # below, instead of value() raising a bare ValueError.
        comp_vals = [
            value(self.flow_mass_comp[j], exception=False) for j in self.flow_mass_comp
        ]
        size_vals = [
            value(self.flow_mass_size[k], exception=False) for k in self.flow_mass_size
        ]
        all_vals = comp_vals + size_vals

        for v in all_vals:
            if v is None or not math.isfinite(v):
                raise ConfigurationError(
                    f"{self.name}: feed contains a non-finite flow value."
                )
        for v in all_vals:
            if v < 0.0:
                raise ConfigurationError(
                    f"{self.name}: feed contains a negative flow value."
                )

        sum_comp = sum(comp_vals)
        sum_size = sum(size_vals)
        if sum_comp <= _EPS_NUMERIC or sum_size <= _EPS_NUMERIC:
            raise ConfigurationError(
                f"{self.name}: degenerate feed -- total solids flow is zero."
            )

        residual = abs(sum_comp - sum_size)
        reference = max(sum_comp, _EPS_NUMERIC)
        if residual > tol_rel * reference:
            raise ConfigurationError(
                f"{self.name}: composition total ({sum_comp:.8g}) and PSD total "
                f"({sum_size:.8g}) are inconsistent (relative residual "
                f"{residual / reference:.3e} > {tol_rel:.1e}). A PSD fraction "
                "vector that does not sum to 1 manifests as this residual once "
                "scaled to flows."
            )
        return True

    def get_material_flow_terms(self, p, j):
        return self.flow_mass_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "flow_mass_comp": self.flow_mass_comp,
            "flow_mass_size": self.flow_mass_size,
        }
