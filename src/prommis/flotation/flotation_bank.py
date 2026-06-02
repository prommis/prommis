#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Dry-solids Flotation Bank
-------------------------

The Flotation Bank unit model represents a bank of flotation cells as a
component-wise dry-solids separation between concentrate and tails streams. The
unit supports three recovery bases:

* ``fixed``: component recoveries are user-fixed or otherwise constrained by the
  surrounding flowsheet, and the unit applies those recoveries as direct stream
  split fractions.
* ``kinetic_closed_form``: component recoveries are calculated from closed-form
  first-order tanks-in-series kinetics using fitted apparent rate constants and
  bank residence time.
* ``kinetic_cell_balance``: component recoveries emerge from an explicit
  fixed-inventory CSTR cascade. Each cell carries component solid holdup, and
  component flotation flow is first order in that holdup.

The model is intended for dry-solids beneficiation flowsheets. It conserves
solid component mass but does not include water as a stream component. Kinetic
modes infer bank-local apparent water or slurry quantities from dry-solids flow,
``pulp_solids_mass_fraction``, and fixed density assumptions only for residence
time and inventory calculations. These inferred quantities are not conserved
stream states.

Configuration Arguments
-----------------------

The following configuration arguments are used by the Flotation Bank model.

``property_package``
    Physical property package used to construct inlet, concentrate, and tails
    state blocks. The package must provide component dry-solid mass flow terms.

``property_package_args``
    Optional arguments passed to each state block created by the property
    package.

``recovery_basis``
    Selects how the local component recovery is determined. Use ``"fixed"`` for
    a recovery-split bank with free recovery variables. Use
    ``"kinetic_closed_form"`` for a bank where recovery is constrained by
    apparent first-order kinetics (rate constant ``k_cf``). Use
    ``"kinetic_cell_balance"`` for a bank with explicit per-cell solid holdup,
    pulp outlet flow, and flotation flow balances (rate constant ``k_cb``).

``number_of_cells``
    Number of perfectly mixed cells in series. This is used in
    ``recovery_basis="kinetic_closed_form"`` and ``"kinetic_cell_balance"``. In
    fixed mode it is accepted but does not add kinetic variables, expressions,
    or constraints.

Degrees of Freedom
------------------

For each time point and component, fixed mode has one degree of freedom in the
local recovery variable, ``recovery[t, j]``. Flowsheets typically fix these
recoveries from plant data or add process constraints that determine them.

Closed-form kinetic mode (``"kinetic_closed_form"``) is a prediction form.
Users fix ``k_cf[t, j]`` and leave ``recovery[t, j]`` free. They also fix
``cell_volume``, ``air_holdup``, and ``pulp_solids_mass_fraction`` for each
bank.

Kinetic-cell-balance mode (``"kinetic_cell_balance"``) is a prediction form.
Users fix ``cell_volume``, ``air_holdup``, ``pulp_solids_mass_fraction``,
``rho_solid``, ``rho_water``, and ``k_cb[t, j]``; the cell balances then
determine cell holdups, outlet flows, concentrate and tails flows, and
``recovery[t, j]``. This mode is steady-state only and rejects ``dynamic=True``.

The bank's ``default_initializer`` is
``FlotationBankKineticCellBalanceInitializer``, a dispatcher that selects
between two interchangeable strategies via its ``strategy`` config option:

* ``strategy="analytical"`` (the default) runs
  ``FlotationBankKineticCellBalanceAnalyticalInitializer``, which uses
  ``cell_cascade_forward`` to compute a constraint-satisfying cascade in pure
  Python and seed it directly (no solver call).
* ``strategy="staged"`` runs
  ``FlotationBankKineticCellBalanceStagedInitializer``, which applies the
  standard IDAES staged fix / deactivate / solve / unfix / reactivate pattern,
  solving the cell cascade through the Pyomo equations themselves.

Select the staged strategy with
``unit.default_initializer(strategy="staged").initialize(unit)``, or import and
use either strategy class directly.

Dynamic-flowsheet compatibility
-------------------------------

The ``fixed`` and ``kinetic`` modes are algebraic in time — every constraint
is evaluated independently at each time point and the unit carries no
accumulation state. They can be embedded in dynamic flowsheets without
restriction; the state blocks propagate the flowsheet's dynamic context. Only
the ``kinetic_cell_balance`` mode rejects ``dynamic=True``, because its
``cell_total_solid_holdup`` and ``cell_solid_holdup`` Vars are named like
inventory accumulators but are algebraically closed by the geometric inventory
expression — a future dynamic variant of that mode will add the corresponding
``dM/dt`` terms before lifting the restriction.

Model Structure
---------------

The Flotation Bank creates three state blocks from the configured property
package:

* ``properties_in``: inlet dry-solids state.
* ``properties_concentrate``: concentrate dry-solids state.
* ``properties_tails``: tails dry-solids state.

The unit exposes corresponding ``inlet``, ``concentrate``, and ``tails`` ports.
All recovery bases expose the same ports, so downstream flowsheets can switch
formulations without changing stream connections. The fixed and closed-form
kinetic modes use direct stream-split constraints. The kinetic-cell-balance mode
instead uses cell material balances and outlet equations; it does not construct
``concentrate_split_eq`` or ``tails_split_eq``.

Additional Variables and Parameters
-----------------------------------

The following variables and parameters are added beyond the property package.

================================= ======================== ======================
Name                              Symbol                   Notes
================================= ======================== ======================
``recovery[t, j]``                :math:`R_{t,j}`          Local recovery of
                                                           component ``j`` to
                                                           concentrate.
``cell_volume``                   :math:`V_{cell}`         Volume per flotation
                                                           cell; kinetic modes.
``air_holdup``                    :math:`\epsilon_{air}`   Volumetric air hold-up
                                                           fraction; kinetic
                                                           modes.
``pulp_solids_mass_fraction``     :math:`w_s`              Inlet pulp solids mass
                                                           fraction used to infer
                                                           bank-local slurry
                                                           quantities; kinetic
                                                           modes.
``rho_solid``                     :math:`\rho_s`           Bank-local dry-solid
                                                           mixture density;
                                                           kinetic-cell-balance
                                                           mode.
``rho_water``                     :math:`\rho_w`           Water density for the
                                                           bank-local slurry
                                                           estimate; kinetic
                                                           modes.
``k_cf[t, j]``                    :math:`k^{cf}_{t,j}`     Closed-form flotation
                                                           rate constant;
                                                           ``kinetic_closed_form``
                                                           mode only.
``k_cb[t, j]``                    :math:`k^{cb}_{t,j}`     Cell-balance flotation
                                                           rate constant;
                                                           ``kinetic_cell_balance``
                                                           mode only.
``R_inf[t, j]``                   :math:`R_{\infty,t,j}`   Ultimate recoverable
                                                           fraction, fixed at
                                                           1.0 by default;
                                                           closed-form kinetic
                                                           mode only.
``cell_total_solid_holdup``       :math:`M_{t,i}`          Total dry-solid
                                                           inventory in cell
                                                           ``i``; kinetic-cell-
                                                           balance mode.
``cell_solid_holdup``             :math:`M_{t,i,j}`        Component dry-solid
                                                           inventory in cell
                                                           ``i``; kinetic-cell-
                                                           balance mode.
``cell_pulp_out_flow``            :math:`P_{t,i,j}`        Component dry-solid
                                                           pulp outlet flow from
                                                           cell ``i``; kinetic-
                                                           cell-balance mode.
``cell_float_flow``               :math:`G_{t,i,j}`        Component dry-solid
                                                           flotation flow from
                                                           cell ``i``; kinetic-
                                                           cell-balance mode.
================================= ======================== ======================

Rate constant note
~~~~~~~~~~~~~~~~~~

``k_cf`` and ``k_cb`` both have units of ``1/hour`` and both drive first-order
flotation via ``F_float = k * M``, but they are **not interchangeable**:

* ``k_cf`` is calibrated against the closed-form tanks-in-series recovery where
  cell inventory is proportional to outlet flow (``M = tau * F_out``).
* ``k_cb`` is calibrated against the fixed-inventory cell cascade where
  inventory is set by cell geometry (``M_total = V * (1 - eps) * rho * w_s``).

The same plant recovery requires different numerical values of ``k`` in each
mode because the underlying inventory model is different. Calibrate each
independently using the appropriate preprocessing script. Do not convert or
reuse ``k_cf`` values for ``kinetic_cell_balance`` or vice versa.

Additional Constraints
----------------------

Fixed and closed-form kinetic modes include the dry-solids split constraints:

.. math::

    F^{conc}_{t,j} = R_{t,j} F^{in}_{t,j}

.. math::

    F^{tails}_{t,j} = F^{in}_{t,j} - F^{conc}_{t,j}

where :math:`F^{in}_{t,j}`, :math:`F^{conc}_{t,j}`, and
:math:`F^{tails}_{t,j}` are component dry-solid mass flows.

Kinetic mode also includes a tanks-in-series first-order model recovery constraint:

.. math::

    R_{t,j} = R_{\infty,t,j}
    \left(1 - \left(1 + \frac{k_{t,j}\tau_t}{N}\right)^{-N}\right)

where :math:`N` is ``number_of_cells`` and :math:`\tau_t` is the active slurry
residence time,  :math:`R_{\infty,t,j}` is the ultimate or maximum recoverable
fraction, :math:`k_{t,j}` is the flotation rate constant.

Kinetic-cell-balance mode replaces the direct split equations with a cell
cascade. For each cell ``i`` and component ``j``:

.. math::

    F^{feed}_{t,i,j} = P_{t,i,j} + G_{t,i,j}

.. math::

    G_{t,i,j} = k_{t,j} M_{t,i,j}

The total cell inventory is fixed by the geometric closure:

.. math::

    M_{t,i} = V_{cell}(1 - \epsilon_{air})\rho_{slurry}w_s

and component holdups sum to ``M[t, i]``. A well-mixed algebraic closure is
written for all but one component; the omitted component follows from the total
holdup and material balances. Concentrate flow is the sum of all cell flotation
flows, tails flow is the last cell pulp outlet, and recovery is still reported
through ``recovery[t, j]``.

The forward cascade helper used by the initializer validates finite and numeric
inputs, applies ``FEED_FLOOR_KG_PER_H`` to component feeds at each cell entry,
treats only exactly ``k == 0`` as a zero-rate carrier, and raises
``CascadeInfeasibleError`` with a ``failure_mode`` when no positive cell root
exists.

Expressions and Performance Variables
-------------------------------------

The unit reports solid mass pull in all modes. Closed-form kinetic mode also
reports residence time, inlet solid volumetric flow, inferred inlet water
volumetric flow, apparent inlet slurry flow, and component rate constants.
Kinetic-cell-balance mode reports apparent bank-inlet water flow, apparent
slurry flow, apparent cell and bank residence times, the first-cell inventory
residence time, and component rate constants.

In closed-form kinetic mode, residence time is calculated as:

.. math::

    \tau_t =
    \frac{V_{cell} N (1 - \epsilon_{air})}
         {Q^{solid}_{t,in} + Q^{water}_{t,in}}

The apparent inlet water flow is inferred from dry solids as:

.. math::

    \dot{m}^{water}_{t,in} =
    \dot{m}^{solid}_{t,in}\frac{1 - w_s}{w_s}

This water estimate is bank-local. It should not be interpreted as a conserved
water balance through recycle mixers or as a replacement for a slurry property
package.
"""

import math

from pyomo.common.config import ConfigBlock, ConfigValue, In, PositiveInt
from pyomo.environ import (
    Expression,
    Param,
    PositiveReals,
    RangeSet,
    Var,
    units,
    value,
)

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe

from prommis.flotation.initializer import (
    FlotationBankKineticCellBalanceInitializer,
    FlotationBankKineticClosedFormInitializer,
    _cell_feed_nominal,
    _safe_value,
)

__author__ = "Daison Yancy Caballero"


class FlotationBankScaler(CustomScalerBase):
    """Scaler for flotation-bank recovery modes."""

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for state_name in (
            "properties_in",
            "properties_concentrate",
            "properties_tails",
        ):
            self.call_submodel_scaler_method(
                submodel=getattr(model, state_name),
                submodel_scalers=submodel_scalers,
                method="variable_scaling_routine",
                overwrite=overwrite,
            )
        for var in model.recovery.values():
            self.set_variable_scaling_factor(var, 1.0, overwrite=overwrite)
        if model.config.recovery_basis == "kinetic_closed_form":
            k_floor = 0.01
            for var in model.k_cf.values():
                nominal = _safe_value(var, default=1.0)
                self.set_variable_scaling_factor(
                    var, 1 / max(abs(nominal), k_floor), overwrite=overwrite
                )
            cell_volume = _safe_value(model.cell_volume, default=1.0)
            self.set_variable_scaling_factor(
                model.cell_volume, 1 / max(abs(cell_volume), 1e-6), overwrite=overwrite
            )
            self.set_variable_scaling_factor(model.air_holdup, 1.0, overwrite=overwrite)
            self.set_variable_scaling_factor(
                model.pulp_solids_mass_fraction, 1.0, overwrite=overwrite
            )
        if model.config.recovery_basis == "kinetic_cell_balance":
            cell_volume = _safe_value(model.cell_volume, default=1.0)
            self.set_variable_scaling_factor(
                model.cell_volume,
                1 / max(abs(cell_volume), 1e-6),
                overwrite=overwrite,
            )
            self.set_variable_scaling_factor(model.air_holdup, 1.0, overwrite=overwrite)
            self.set_variable_scaling_factor(
                model.pulp_solids_mass_fraction, 1.0, overwrite=overwrite
            )
            for var in model.k_cb.values():
                nominal = _safe_value(var, default=0.01)
                self.set_variable_scaling_factor(
                    var, 1 / max(abs(nominal), 0.01), overwrite=overwrite
                )
            for (time, cell), var in model.cell_total_solid_holdup.items():
                nominal = _safe_value(var, default=1.0)
                self.set_variable_scaling_factor(
                    var, 1 / max(abs(nominal), 1.0), overwrite=overwrite
                )
            for var in model.cell_solid_holdup.values():
                nominal = _safe_value(var, default=1e-3)
                self.set_variable_scaling_factor(
                    var, 1 / max(abs(nominal), 1e-3), overwrite=overwrite
                )
            for var_obj in (model.cell_pulp_out_flow, model.cell_float_flow):
                for (time, cell, component), var in var_obj.items():
                    feed_nominal = _cell_feed_nominal(model, time, cell, component)
                    nominal = _safe_value(var, default=feed_nominal)
                    self.set_variable_scaling_factor(
                        var, 1 / max(abs(nominal), 1e-12), overwrite=overwrite
                    )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        if model.config.recovery_basis in ("fixed", "kinetic_closed_form"):
            for (time, component), condata in model.concentrate_split_eq.items():
                self.scale_constraint_by_component(
                    condata,
                    model.properties_concentrate[time].flow_mass_comp[component],
                    overwrite=overwrite,
                )
            for (time, component), condata in model.tails_split_eq.items():
                self.scale_constraint_by_component(
                    condata,
                    model.properties_tails[time].flow_mass_comp[component],
                    overwrite=overwrite,
                )
        if model.config.recovery_basis == "kinetic_closed_form":
            for (time, component), condata in model.kinetic_recovery_eq.items():
                self.scale_constraint_by_component(
                    condata, model.recovery[time, component], overwrite=overwrite
                )
        if model.config.recovery_basis == "kinetic_cell_balance":
            component_list = list(model.config.property_package.component_list)
            for (time, cell), condata in model.geometric_holdup_eq.items():
                nominal = _safe_value(
                    model.cell_total_solid_holdup[time, cell], default=1.0
                )
                self.set_constraint_scaling_factor(
                    condata, 1 / max(abs(nominal), 1.0), overwrite=overwrite
                )
            for (time, cell), condata in model.total_component_holdup_eq.items():
                nominal = _safe_value(
                    model.cell_total_solid_holdup[time, cell], default=1.0
                )
                self.set_constraint_scaling_factor(
                    condata, 1 / max(abs(nominal), 1.0), overwrite=overwrite
                )
            for (time, cell, component), condata in model.well_mixed_eq.items():
                holdup = _safe_value(
                    model.cell_total_solid_holdup[time, cell], default=1.0
                )
                outlet_flow = sum(
                    _safe_value(
                        model.cell_pulp_out_flow[time, cell, local_component],
                        default=1e-3,
                    )
                    for local_component in component_list
                )
                self.set_constraint_scaling_factor(
                    condata,
                    1 / max(abs(holdup * outlet_flow), 1e-3),
                    overwrite=overwrite,
                )
            for (
                time,
                cell,
                component,
            ), condata in model.cell_material_balance_eq.items():
                feed_nominal = _cell_feed_nominal(model, time, cell, component)
                self.set_constraint_scaling_factor(
                    condata, 1 / max(abs(feed_nominal), 1e-12), overwrite=overwrite
                )
            for (time, cell, component), condata in model.flotation_removal_eq.items():
                k_nominal = _safe_value(model.k_cb[time, component], default=0.01)
                holdup = _safe_value(
                    model.cell_solid_holdup[time, cell, component], default=1e-3
                )
                self.set_constraint_scaling_factor(
                    condata,
                    1 / max(abs(k_nominal * holdup), 1e-3),
                    overwrite=overwrite,
                )
            for constraint_obj in (
                model.concentrate_outlet_eq,
                model.tails_outlet_eq,
                model.recovery_eq,
            ):
                for (time, component), condata in constraint_obj.items():
                    feed_nominal = _safe_value(
                        model.properties_in[time].flow_mass_comp[component],
                        default=1e-3,
                    )
                    self.set_constraint_scaling_factor(
                        condata,
                        1 / max(abs(feed_nominal), 1e-3),
                        overwrite=overwrite,
                    )


@declare_process_block_class("FlotationBank")
class FlotationBankData(UnitModelBlockData):
    """Dry-solids flotation bank with fixed, kinetic, and cell-balance modes.

    All modes expose inlet, concentrate, and tails dry-solids ports plus local
    component recovery variables. Fixed mode leaves recoveries free for the
    flowsheet to fix or constrain. Closed-form kinetic mode constrains
    recoveries with first-order tanks-in-series kinetics. Kinetic-cell-balance
    mode replaces the direct split equations with a steady-state fixed-inventory
    CSTR cascade that determines cell holdups, outlet flows, and recoveries.
    """

    default_scaler = FlotationBankScaler
    default_initializer = BlockTriangularizationInitializer

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for the bank",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use when constructing property blocks",
        ),
    )
    CONFIG.declare(
        "recovery_basis",
        ConfigValue(
            default="fixed",
            domain=In(["fixed", "kinetic_closed_form", "kinetic_cell_balance"]),
            description="How recovery is determined",
        ),
    )
    CONFIG.declare(
        "number_of_cells",
        ConfigValue(
            default=1,
            domain=PositiveInt,
            description="Number of perfectly mixed cells in series",
        ),
    )

    def build(self):
        super().build()

        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )
        self.properties_concentrate = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=False,
            **self.config.property_package_args,
        )
        self.properties_tails = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=False,
            **self.config.property_package_args,
        )

        self.add_port("inlet", self.properties_in)
        self.add_port("concentrate", self.properties_concentrate)
        self.add_port("tails", self.properties_tails)

        component_list = self.config.property_package.component_list
        self.recovery = Var(
            self.flowsheet().time,
            component_list,
            bounds=(0, 1),
            units=units.dimensionless,
            initialize=0.5,
            doc="Local component recovery to concentrate",
        )

        if self.config.recovery_basis == "kinetic_closed_form":
            number_of_cells = self.config.number_of_cells

            self.cell_volume = Var(
                initialize=1.7,
                bounds=(1e-6, None),
                units=units.m**3,
                doc="Volume per flotation cell",
            )
            self.air_holdup = Var(
                initialize=0.10,
                bounds=(0, 0.5),
                units=units.dimensionless,
                doc="Volumetric air hold-up fraction in the cells",
            )
            self.pulp_solids_mass_fraction = Var(
                initialize=0.325,
                bounds=(1e-3, 1 - 1e-3),
                units=units.dimensionless,
                doc="Solids mass fraction of the inlet pulp",
            )
            self.rho_water = Param(
                initialize=1000.0,
                mutable=True,
                domain=PositiveReals,
                units=units.kg / units.m**3,
                doc="Water mass density used for apparent inlet slurry flow",
            )
            self.k_cf = Var(
                self.flowsheet().time,
                component_list,
                initialize=1.0,
                bounds=(0, None),
                units=1 / units.hour,
                doc="Closed-form apparent first-order flotation rate constant.",
            )
            self.R_inf = Var(
                self.flowsheet().time,
                component_list,
                initialize=1.0,
                # The tiny upper tolerance avoids diagnostics warning on the
                # default fixed-at-one value while preserving the intended
                # ultimate-recovery bound for estimation workflows.
                bounds=(0, 1 + 1e-8),
                units=units.dimensionless,
                doc="Ultimate recoverable fraction fixed at 1.0 by default",
            )
            for time in self.flowsheet().time:
                for component in component_list:
                    # Fix R_inf for current implementation. IF rich recovery-time
                    # data are available, it can be estimated.
                    self.R_inf[time, component].fix(1.0)
            self.effective_volume = Expression(
                expr=self.cell_volume * number_of_cells * (1 - self.air_holdup),
                doc="Active slurry volume across the bank",
            )

            # Water is not a conserved stream component in this dry-solids
            # property package. Kinetic mode infers a bank-local apparent
            # inlet water flow from the dry-solids flow and this bank's
            # pulp-solids fraction only to size residence time. Mixer water
            # imbalances therefore indicate unmodeled dilution, dewatering, or
            # water adjustment between banks, not a dry-solids mass-balance
            # failure.
            # TODO: Replace this bank-local water estimate with a slurry-aware
            # property package when water conservation is needed. That future
            # package should carry water flow or pulp solids fraction as a
            # stream state, expose slurry volumetric flow on the state block,
            # and let mixers balance water through the recycle circuit.
            @self.Expression(self.flowsheet().time)
            def flow_mass_water(b, t):
                mass_solid = sum(
                    b.properties_in[t].flow_mass_comp[j] for j in component_list
                )
                return (
                    mass_solid
                    * (1 - b.pulp_solids_mass_fraction)
                    / b.pulp_solids_mass_fraction
                )

            @self.Expression(self.flowsheet().time)
            def flow_vol_water(b, t):
                return b.flow_mass_water[t] / b.rho_water

            @self.Expression(self.flowsheet().time)
            def flow_vol_slurry(b, t):
                return b.properties_in[t].flow_vol_solid + b.flow_vol_water[t]

            @self.Expression(self.flowsheet().time)
            def tau(b, t):
                return b.effective_volume / b.flow_vol_slurry[t]

            @self.Constraint(self.flowsheet().time, component_list)
            def kinetic_recovery_eq(b, t, j):
                return b.recovery[t, j] == b.R_inf[t, j] * (
                    1
                    - (1 + b.k_cf[t, j] * b.tau[t] / number_of_cells)
                    ** (-number_of_cells)
                )

            self.default_initializer = FlotationBankKineticClosedFormInitializer

        elif self.config.recovery_basis == "kinetic_cell_balance":
            self._validate_kinetic_cell_balance_config()
            self._build_kinetic_cell_balance(component_list)
            self.default_initializer = FlotationBankKineticCellBalanceInitializer

        if self.config.recovery_basis in ("fixed", "kinetic_closed_form"):

            @self.Constraint(self.flowsheet().time, component_list)
            def concentrate_split_eq(b, t, j):
                return (
                    b.properties_concentrate[t].flow_mass_comp[j]
                    == b.recovery[t, j] * b.properties_in[t].flow_mass_comp[j]
                )

            @self.Constraint(self.flowsheet().time, component_list)
            def tails_split_eq(b, t, j):
                return (
                    b.properties_tails[t].flow_mass_comp[j]
                    == b.properties_in[t].flow_mass_comp[j]
                    - b.properties_concentrate[t].flow_mass_comp[j]
                )

        @self.Expression(self.flowsheet().time)
        def solid_mass_pull(b, t):
            return sum(
                b.properties_concentrate[t].flow_mass_comp[j] for j in component_list
            ) / sum(b.properties_in[t].flow_mass_comp[j] for j in component_list)

    def _validate_kinetic_cell_balance_config(self):
        if self.config.dynamic is True:
            raise ConfigurationError(
                "recovery_basis='kinetic_cell_balance' uses geometric "
                "(algebraic) holdup and does not integrate cell_solid_holdup "
                "over time. dynamic=True is rejected to avoid producing a "
                "model that names a Var 'holdup' but behaves steady-state. "
                "Use recovery_basis='fixed' or 'kinetic_closed_form' (both "
                "algebraic in time) inside dynamic flowsheets, or wait for a "
                "future dynamic kinetic-cell-balance mode with proper "
                "dM/dt terms."
            )

    def _build_kinetic_cell_balance(self, component_list):
        component_names = list(component_list)
        self.cells = RangeSet(1, self.config.number_of_cells)
        self._zero_feed_recovery_fixes = set()
        self._well_mixed_omitted_component = tuple(component_names)[-1]

        self.cell_volume = Var(
            initialize=1.7,
            bounds=(1e-6, None),
            units=units.m**3,
            doc="Volume per flotation cell",
        )
        self.air_holdup = Var(
            initialize=0.10,
            bounds=(0, 0.5),
            units=units.dimensionless,
            doc="Volumetric air hold-up fraction in the cells",
        )
        self.pulp_solids_mass_fraction = Var(
            initialize=0.325,
            bounds=(1e-3, 1 - 1e-3),
            units=units.dimensionless,
            doc="Solids mass fraction of the inlet pulp",
        )
        self.rho_solid = Param(
            initialize=3000.0,
            mutable=True,
            domain=PositiveReals,
            units=units.kg / units.m**3,
            doc="Reference mixture density of dry solids in the bank",
        )
        self.rho_water = Param(
            initialize=1000.0,
            mutable=True,
            domain=PositiveReals,
            units=units.kg / units.m**3,
            doc="Water density used for bank-local slurry estimates",
        )
        self.k_cb = Var(
            self.flowsheet().time,
            component_names,
            initialize=0.01,
            bounds=(0, None),
            units=1 / units.hour,
            doc="Cell-balance first-order flotation rate constant.",
        )
        self.cell_total_solid_holdup = Var(
            self.flowsheet().time,
            self.cells,
            initialize=1.0,
            bounds=(0, None),
            units=units.kg,
            doc="Total dry-solid holdup in each cell",
        )
        self.cell_solid_holdup = Var(
            self.flowsheet().time,
            self.cells,
            component_names,
            initialize=1e-3,
            bounds=(0, None),
            units=units.kg,
            doc="Component dry-solid holdup in each cell",
        )
        self.cell_pulp_out_flow = Var(
            self.flowsheet().time,
            self.cells,
            component_names,
            initialize=1e-3,
            bounds=(0, None),
            units=units.kg / units.hour,
            doc="Component dry-solid pulp outlet flow from each cell",
        )
        self.cell_float_flow = Var(
            self.flowsheet().time,
            self.cells,
            component_names,
            initialize=1e-3,
            bounds=(0, None),
            units=units.kg / units.hour,
            doc="Component dry-solid flotation flow from each cell",
        )

        @self.Expression(self.flowsheet().time, self.cells, component_names)
        def cell_feed(b, t, i, j):
            if i == b.cells.first():
                return b.properties_in[t].flow_mass_comp[j]
            return b.cell_pulp_out_flow[t, b.cells.prev(i), j]

        @self.Expression(self.flowsheet().time, self.cells)
        def cell_total_pulp_out_flow(b, t, i):
            return sum(b.cell_pulp_out_flow[t, i, j] for j in component_names)

        @self.Expression()
        def rho_slurry(b):
            w_s = b.pulp_solids_mass_fraction
            return 1 / (w_s / b.rho_solid + (1 - w_s) / b.rho_water)

        @self.Expression(self.flowsheet().time)
        def flow_mass_water_bank_inlet(b, t):
            mass_solid = sum(
                b.properties_in[t].flow_mass_comp[j] for j in component_names
            )
            return (
                mass_solid
                * (1 - b.pulp_solids_mass_fraction)
                / b.pulp_solids_mass_fraction
            )

        @self.Expression(self.flowsheet().time)
        def Q_slurry_in(b, t):
            mass_solid = sum(
                b.properties_in[t].flow_mass_comp[j] for j in component_names
            )
            return (mass_solid + b.flow_mass_water_bank_inlet[t]) / b.rho_slurry

        @self.Expression(self.flowsheet().time)
        def tau_apparent_cell(b, t):
            return b.cell_volume * (1 - b.air_holdup) / b.Q_slurry_in[t]

        @self.Expression(self.flowsheet().time)
        def tau_apparent_bank(b, t):
            return b.config.number_of_cells * b.tau_apparent_cell[t]

        @self.Expression(self.flowsheet().time, self.cells)
        def tau_inventory_cell(b, t, i):
            return b.cell_total_solid_holdup[t, i] / sum(
                b.cell_feed[t, i, j] for j in component_names
            )

        @self.Constraint(self.flowsheet().time, self.cells)
        def geometric_holdup_eq(b, t, i):
            return (
                b.cell_total_solid_holdup[t, i]
                == b.cell_volume
                * (1 - b.air_holdup)
                * b.rho_slurry
                * b.pulp_solids_mass_fraction
            )

        @self.Constraint(self.flowsheet().time, self.cells)
        def total_component_holdup_eq(b, t, i):
            return (
                sum(b.cell_solid_holdup[t, i, j] for j in component_names)
                == b.cell_total_solid_holdup[t, i]
            )

        @self.Constraint(self.flowsheet().time, self.cells, component_names)
        def cell_material_balance_eq(b, t, i, j):
            return (
                b.cell_feed[t, i, j]
                == b.cell_pulp_out_flow[t, i, j] + b.cell_float_flow[t, i, j]
            )

        @self.Constraint(self.flowsheet().time, self.cells, component_names)
        def flotation_removal_eq(b, t, i, j):
            return (
                b.cell_float_flow[t, i, j]
                == b.k_cb[t, j] * b.cell_solid_holdup[t, i, j]
            )

        well_mixed_components = [
            j for j in component_names if j != self._well_mixed_omitted_component
        ]

        @self.Constraint(self.flowsheet().time, self.cells, well_mixed_components)
        def well_mixed_eq(b, t, i, j):
            return (
                b.cell_solid_holdup[t, i, j] * b.cell_total_pulp_out_flow[t, i]
                == b.cell_total_solid_holdup[t, i] * b.cell_pulp_out_flow[t, i, j]
            )

        @self.Constraint(self.flowsheet().time, component_names)
        def concentrate_outlet_eq(b, t, j):
            return b.properties_concentrate[t].flow_mass_comp[j] == sum(
                b.cell_float_flow[t, i, j] for i in b.cells
            )

        @self.Constraint(self.flowsheet().time, component_names)
        def tails_outlet_eq(b, t, j):
            return (
                b.properties_tails[t].flow_mass_comp[j]
                == b.cell_pulp_out_flow[t, b.cells.last(), j]
            )

        @self.Constraint(self.flowsheet().time, component_names)
        def recovery_eq(b, t, j):
            return (
                b.properties_concentrate[t].flow_mass_comp[j]
                == b.recovery[t, j] * b.properties_in[t].flow_mass_comp[j]
            )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Inlet": self.properties_in,
                "Concentrate": self.properties_concentrate,
                "Tails": self.properties_tails,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        def numeric(obj):
            val = value(obj, exception=False)
            return math.nan if val is None else val

        variables = {
            "Mass Pull": numeric(self.solid_mass_pull[time_point]),
        }
        if self.config.recovery_basis == "kinetic_closed_form":
            variables.update(
                {
                    "Residence Time": numeric(self.tau[time_point]),
                    "Inlet Solid Volumetric Flow": numeric(
                        self.properties_in[time_point].flow_vol_solid
                    ),
                    "Inlet Water Volumetric Flow": numeric(
                        self.flow_vol_water[time_point]
                    ),
                    "Apparent Inlet Slurry Flow": numeric(
                        self.flow_vol_slurry[time_point]
                    ),
                }
            )
            variables.update(
                {
                    f"k_cf[{component}]": numeric(self.k_cf[time_point, component])
                    for component in self.config.property_package.component_list
                }
            )
        if self.config.recovery_basis == "kinetic_cell_balance":
            variables.update(
                {
                    "Inlet Water Mass Flow": numeric(
                        self.flow_mass_water_bank_inlet[time_point]
                    ),
                    "Apparent Inlet Slurry Flow": numeric(self.Q_slurry_in[time_point]),
                    "Apparent Cell Residence Time": numeric(
                        self.tau_apparent_cell[time_point]
                    ),
                    "Apparent Bank Residence Time": numeric(
                        self.tau_apparent_bank[time_point]
                    ),
                    "Cell 1 Inventory Residence Time": numeric(
                        self.tau_inventory_cell[time_point, self.cells.first()]
                    ),
                }
            )
            variables.update(
                {
                    f"k_cb[{component}]": numeric(self.k_cb[time_point, component])
                    for component in self.config.property_package.component_list
                }
            )
        return {"vars": variables}
