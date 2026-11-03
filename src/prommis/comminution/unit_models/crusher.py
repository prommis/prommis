#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Particle Size Distribution Crusher
----------------------------------

``CrusherSolidPSD`` models size reduction in a stream whose solids are tracked by
component and particle size. Each sized solid component has its own feed particle size
distribution (PSD). The unit moves each component's mass among size intervals, usually
toward finer sizes, while conserving its total mass. It also calculates the power
required for crushing.

The unit uses a ``SolidPSDParameterBlock``. A sized solid component is any solid in the
package's ``sized_solid_list``, whether or not it is a mineral. Unsized solids,
liquid, vapor, temperature, and pressure pass through unchanged. The crusher is
steady state (set ``dynamic=False`` and ``has_holdup=False``).

Product Methods
---------------

Set ``psd_method`` to match the data you have:

* ``"tabular"``: use measured or specified product size fractions.
* ``"distribution_function"``: use a product P80 or measured power with an assumed
  distribution curve.
* ``"selection_breakage"``: use data that describes how the feed breaks. This is the
  default ``psd_method``.

P80 is the particle size below which 80% of the stream's sized solid mass lies.

The full feed size distribution affects the product only in ``selection_breakage``. The
other methods prescribe a product shape. The feed affects power in every method. A
power-specified distribution also uses the feed P80 to find the product P80.

Method: ``tabular``
-------------------

Give the product fractions with ``tabular_psd``. The crusher multiplies each component's
total inlet flow by these fractions. The same row applies to every sized solid component. Use
``tabular_psd_by_comp`` for a different row per component. Give one form, not both.

Each row needs ``N`` nonnegative fractions, finest interval first. Its sum must be
within ``1e-5`` of one. The unit normalizes it to one as needed. A per-component
dictionary must name every sized solid component. The initializer checks that each component's
cumulative product is no coarser than its feed. This method does not use the feed
fractions to shape the product.

Method: ``distribution_function``
---------------------------------

This method puts an analytic curve on the product. Use
``distribution_shape="rosin_rammler"`` for Rosin-Rammler [8], the default, or
``"gates_gaudin_schuhmann"`` for Gates-Gaudin-Schuhmann [9]. Set a positive
``shape_exponent`` to control the spread. The same product fractions apply to every
sized solid component. The mesh needs at least two size intervals.

There are two ways to place the curve:

* Set ``target_product_p80`` to the analytic product P80 you want, in meters. The
  crusher calculates the power.
* Omit ``target_product_p80`` and fix ``power[t]`` at every time point. The crusher
  solves for the curve's size parameter and product P80.

The target P80 must be between the point 80 percent through the finest interval and the
top mesh edge. The curve must also fit the mesh, which caps the target at
:math:`x_N (\ln 5/\ln 20)^{1/n}` for Rosin-Rammler and :math:`x_N \cdot 0.8^{1/n}` for
Gates-Gaudin-Schuhmann, where :math:`x_N` is the top mesh edge and :math:`n` is
``shape_exponent``. A fixed power must imply a P80 in the same range. In power mode,
the size parameter is ``d_63[t]`` for Rosin-Rammler or ``d_max[t]`` for
Gates-Gaudin-Schuhmann.

Method: ``selection_breakage``
------------------------------

This method solves a mass-based population balance for each sized solid component [2]. It is the
default product method. For each parent size interval, the selection fraction :math:`S`
specifies how much mass breaks, from zero to one. The breakage fractions :math:`B` send
that mass to the same or finer intervals and sum to one for each parent interval. The
product therefore depends on each component's full feed size distribution.

**There are two ways to run it:**

* ``recycle="classification_recycle"``, the default, sends classified oversize back to
  the crushing zone. It uses the Whiten method described below.
* ``recycle="single_pass"`` applies the selection and breakage step ``n_stress_events``
  times. Use an integer from 1 to 1000. Each extra event adds intermediate ``psd_stage``
  flows and makes the model larger.

**Selection functions** (``selection_function``):

Put the function inputs in ``selection_params``. Sizes are in meters. The unit checks
for missing and unknown keys.

* ``"whiten"``, the default: selection rises from zero below ``K1`` to one
  above ``K2`` [1, 7]. Give ``K1`` and ``K2``; ``K3`` is optional and defaults
  to 2.3. This is the only function allowed with recycle.
* ``"king"``: a ramp between two multiples of the closed side setting [2]. Give
  ``CSS``, ``alpha1``, ``alpha2``, and ``n``.
* ``"austin"``: a size-based power law [3]. Give ``S1``, ``d1``, and ``a``.
* ``"vogel_peukert"``: selection from impact energy [5]. Give ``f_mat``,
  ``xw_min``, and ``v``.
* ``"user"``: give a list of ``N`` values from zero to one, in mesh order.

**Breakage functions** (``breakage_function``):

Put the function inputs in ``breakage_params``.

* ``"luckie_austin"``, the default, or ``"reid_stewart"``: a two-term
  cumulative curve [3, 4]. Give ``phi``, ``gamma``, and ``beta``. The two
  functions use different size edges for the second term.
* ``"vogel"``: a power law with a smooth cut near fragment size ``dprime`` [5].
  Give ``q`` and ``dprime``.
* ``"t10_appearance"`` [6, 7]: use only with classification recycle. Give
  ``t10`` and an ``appearance_table``. Pendulum power also needs an ``ecs_table``.
* ``"user"``: give an :math:`N` by :math:`N` matrix :math:`B_{i,j}`, with daughter
  interval :math:`i` and parent interval :math:`j`. Entries must be nonnegative and
  zero when :math:`i > j`. Each column must sum to one. The unit normalizes each
  column.

**Different values for each component:**

Use ``selection_params_by_comp`` or ``breakage_params_by_comp`` when components need
different values. Each dictionary must name every sized solid component. For each function, give
either its common input or its ``_by_comp`` input, not both.

In single-pass mode, the selection parameters may differ by component. ``CSS`` for King
and ``v`` for Vogel-Peukert are machine settings, so each component must have the same
value. Classification recycle requires common ``selection_params``, but allows
per-component breakage.

Classification Recycle (Whiten Method)
--------------------------------------

The Whiten circuit [1] uses one classification fraction :math:`C_k` for each size
interval. The fraction sent back through the crusher is common to all components. The unit
solves the internal load ``recycle_load[t, s, k]`` and sends the remaining mass to the
product. It accepts only ``selection_function="whiten"`` with common
``selection_params``. Keep ``n_stress_events=1``.

Power
-----

``power_law`` selects ``"bond"``, ``"rittinger"``, ``"kick"``, or ``"pendulum"``. Bond
is the default. The first three laws use ``bond_work_index`` in kWh/t, from 1 to 150.
Its default is 14.0. They calculate power from the total sized feed flow and the feed
and product P80 values [10-12]. The Rittinger and Kick coefficients are chosen so
each law's slope in size approximately matches Bond's at its boundary size, 50 µm for
Rittinger and 50 mm for Kick [14].

Use ``"rittinger"`` for product P80 up to 0.05 mm, ``"bond"`` from 0.05 to 50 mm,
and ``"kick"`` from 50 mm upward [13]. These are approximate applicability ranges;
they do not constrain the model's solution.

Pendulum power uses the specific comminution energy of the classified load [7]. It
requires selection and breakage with classification recycle,
``breakage_function="t10_appearance"``, and an ``ecs_table`` in the breakage data. Set
``pendulum_power_params`` to a dictionary with a nonnegative ``power_factor`` and
nonnegative ``no_load_power`` in kW. Do not set ``bond_work_index`` for this law.

Equipment and Stage
-------------------

Both ``crusher_equipment`` and ``crusher_stage`` are optional. Use them to supply
defaults and to check whether the product P80 is typical for the chosen stage.

Set ``crusher_equipment`` to ``"jaw"``, ``"gyratory"``, ``"cone"``, ``"roll1"``,
``"roll2"``, ``"short_head_cone"``, or ``"hammer_mill"``. Equipment does two separate
things:

* Stage inference, for every ``psd_method``: it sets the matching ``crusher_stage``.
* Function presets, for ``selection_breakage`` only: it fills in ``selection_function``
  and ``breakage_function`` if you did not set them. Functions you choose yourself take
  precedence. Equipment supplies function names only; you must still give the inputs
  for the chosen functions. The presets are engineering defaults based on the
  equipment and function pairings documented for a commercial crusher simulator; they
  are not fitted to data.

``"gyratory"``, ``"roll1"``, and ``"roll2"`` have no built-in function choices. The unit
warns and keeps the general function defaults unless you select the functions yourself.
``"cone"``, ``"short_head_cone"``, and ``"hammer_mill"`` choose non-Whiten selection
functions. Use ``recycle="single_pass"`` with those defaults.

You can set ``crusher_stage`` on its own to ``"primary"``, ``"secondary"``, or
``"tertiary"``. The expected product P80 ranges are at least 10 cm for primary, 2 to
under 10 cm for secondary, and 0.5 to under 2 cm for tertiary [15]. After solving,
``check_applicability()`` warns if the product P80 falls outside the range. The stage
adds no constraint. If you set both options, the stage must match the equipment.

Degrees of Freedom
------------------

With the inlet state fixed, the tabular method, a distribution with
``target_product_p80``, and the selection and breakage methods determine the product
size distribution. Their power follows from the power equation. For a distribution
without ``target_product_p80``, fix ``power[t]`` at every time point. The model then
solves for the distribution's size parameter.

``bond_work_index`` is fixed at build time for Bond, Rittinger, and Kick. At one time
point, a study can unfix it and fix measured power to solve for the work index. Pendulum
power fixes its power factor and no-load power at build time. To fit the power factor,
unfix it and fix measured power.

To change the crushing data, build a new unit. The initializer also rejects a sized feed
below 1e-6 kg/s. ``CrusherSolidPSDScaler`` supplies scaling
factors from the configured data and fixed feed when ``factor_source="input_based"``.

Model Structure
---------------

The unit creates the inlet and outlet state blocks, ``properties_in`` and
``properties_out``, from the property package. The ``inlet`` and ``outlet`` ports expose
these states to a flowsheet. The crusher has no control volume. It adds product size
equations for each sized solid component and constraints that copy the pass-through state
variables to the outlet. It calculates power from the crushing duty; it has no energy or
momentum balance.

Additional Variables and Parameters
-----------------------------------

The crusher adds the following Pyomo components alongside its inlet and outlet states. It
creates method-specific components only when their method is active. ``bond_work_index``
is always present, although pendulum power does not use it. The crushing-data
``Param`` objects (``tabular_psd``, ``selection``, ``breakage``, ``classification`` and
``ecs``) are mutable internally; changing them after build is unsupported.

============================= ============= ======================================
Name                          Units         Description
============================= ============= ======================================
``power[t]``                  kW            Crusher power at each time point.
``bond_work_index``           kWh/t         Work index; fixed at build.
``tabular_psd[s, k]``         dimensionless Product row for tabular mode.
``selection[s, k]``           dimensionless Selection by component and size, single
                                            pass only. Recycle uses its common
                                            row as ``classification``.
``breakage[s, k, i]``         dimensionless Daughter size ``k`` from parent size
                                            ``i``.
``classification[k]``         dimensionless Fraction returned to the crusher in
                                            Whiten recycle.
``psd_stage[t, r, s, k]``     kg/s          Intermediate rows between stress
                                            events when ``n_stress_events > 1``.
``recycle_load[t, s, k]``     kg/s          Internal load for classification
                                            recycle.
``d_63[t]`` or ``d_max[t]``   m             Size parameter for a distribution
                                            method.
``target_product_p80``        m             Target analytic P80 for a
                                            size-specified distribution; may be
                                            changed after build.
``feed_p80[t]``               m             Smoothed feed P80 used by the
                                            size-specified distribution's
                                            no-coarsening rule.
``ecs[s, k]``                 kWh/t         Specific comminution energy for
                                            pendulum power.
``pendulum_power_factor``     dimensionless Multiplier fixed at build.
``no_load_power``             kW            Idle power fixed at build.
============================= ============= ======================================

Additional Constraints
----------------------

Every product method conserves each sized solid component separately. The method-specific
equations determine the outlet flows above the finest interval. The flow in interval 0
is then the inlet component total minus the outlet flows in all other intervals:

.. math::

    m^{out}_{t,s,0}
    = m^{in}_{t,s,total} - \sum_{k=1}^{N-1} m^{out}_{t,s,k}

Here :math:`m^{in}_{t,s,total}` is the inlet flow of component :math:`s`, summed over all
size intervals at time :math:`t`. This closure also collects any material finer than the
mesh in interval 0.

For ``tabular``, each outlet flow above interval 0 equals the component's inlet total
multiplied by its specified product fraction. For ``distribution_function``, the unit
calculates the fractions from the selected cumulative curve, normalizes them over the
mesh, and applies the same fractions to every component. The distribution's size parameter
is set by ``target_product_p80`` or determined by fixed ``power[t]`` through the power
equation.

For ``selection_breakage`` with ``recycle="single_pass"``, the flow from parent interval
:math:`i` is split into an unbroken fraction :math:`1-S_{s,i}` and a broken fraction
:math:`S_{s,i}`. The breakage fraction :math:`B_{s,k,i}` sends the broken part to
daughter interval :math:`k`. One stress event has the transfer coefficient

.. math::

    A_{s,k,i}
    = \delta_{k,i}(1-S_{s,i}) + B_{s,k,i}S_{s,i}

where :math:`\delta_{k,i}` is one when :math:`k=i` and zero otherwise. Only the same or
finer intervals can receive mass. With more than one stress event, the unit applies this
transfer to each intermediate ``psd_stage`` row before calculating the product row.

Classification recycle solves the internal load before calculating the product. The
common fraction :math:`C_k` returns material in interval :math:`k` to the crusher, while
the remaining fraction leaves as product. For each component, the recycle and product
equations are

.. math::

    (I-B_s C)X_{t,s} = F_{t,s},
    \qquad Y_{t,s} = (I-C)X_{t,s}

Here :math:`F_{t,s}` is the inlet size-flow vector for component :math:`s` at time
:math:`t`, :math:`X_{t,s}` is the internal load, and :math:`Y_{t,s}` is the product
size-flow vector. :math:`C` is the diagonal matrix of classification fractions.

The tabular method also constrains the cumulative product flow of each component to be at
least its cumulative inlet flow at every interior size boundary. This prevents a coarser
product. A distribution specified by ``target_product_p80`` instead checks the
aggregate sized stream: it accepts a product P80 no greater than the feed P80, or
product cumulative passing fractions no lower than the feed's at every interior
boundary.

For Bond, Rittinger, and Kick, the power equation multiplies total sized inlet flow by
the specific energy calculated from the feed and product P80 values. With flow in kg/s
and energy in kWh/t, this is

.. math::

    P_t = 3.6\,m^{in}_{t,sized} E_t

The unit floors negative specific energy at zero, except when fixed power determines a
distribution curve. Pendulum power uses the energy of the classified internal load
instead:

.. math::

    P_t = A\left(3.6\sum_{s,k} E^{cs}_{s,k} C_k X_{t,s,k}\right) + P_0

:math:`A` is ``pendulum_power_factor``, :math:`P_0` is ``no_load_power``, and
:math:`E^{cs}_{s,k}` comes from the configured ``ecs_table``.

Expressions and Reported Values
-------------------------------

The outlet state contains the product flow of each component in each size interval.
``product_p80[t]`` gives the product size at 80 percent passing. It is analytic for a
distribution function and comes from the outlet state's smoothed PSD for the other
methods. The distribution method also exposes ``shape_cdf[t, i]`` at the mesh edges and
``product_mass_frac[t, k]`` for the normalized product fractions. Bond, Rittinger, and
Kick use the inlet state's smoothed P80 for the total sized feed. A smoothed P80
approximates the size at 80 percent passing.

For a distribution function, ``product_p80[t]`` comes from the full size
curve. ``target_product_p80`` sets this value when supplied. Power and
``check_applicability()`` use it too.

The outlet stream uses size fractions from that curve, but leaves out any
material above the mesh's largest size. The model scales the remaining
fractions to sum to one, then calculates ``properties_out[t].percentile_size[0.8]``
from them. That is why the two P80 values can differ.

To reduce the difference, extend the mesh to larger sizes and use smaller
intervals near P80.

``report()`` shows ``power[t]``, ``product_p80[t]`` and an inlet and outlet stream
table. Pendulum mode also reports ``pendulum_power[t]``, the calculated power before the
power factor and no-load power. After solving, call ``check_applicability()`` to check
the product P80 against the selected crusher stage's range. It returns any warnings
and, by default, logs them. It does not change the solution.

References
----------

[1] Whiten, W. J. (1972). The simulation of crushing plants with models
developed using multiple spline regression. Proc. 10th APCOM Symposium,
Johannesburg, SAIMM, 317-323. Classification recycle and the Whiten
function.

[2] King, R. P. (2001). Modeling and Simulation of Mineral Processing
Systems. Butterworth-Heinemann, Oxford. King selection and population
balance.

[3] Austin, L. G., Klimpel, R. R., and Luckie, P. T. (1984). Process
Engineering of Size Reduction: Ball Milling. SME-AIME, New York. Austin
selection and Luckie-Austin breakage.

[4] Reid, K. J. (1965). A solution to the batch grinding equation. Chemical
Engineering Science, 20(11), 953-963. Batch-grinding population balance.

[5] Vogel, L., and Peukert, W. (2003). Breakage behaviour of different
materials: construction of a mastercurve for the breakage probability.
Powder Technology, 129(1-3), 101-110. Vogel-Peukert selection and Vogel
breakage.

[6] Narayanan, S. S., and Whiten, W. J. (1988). Determination of comminution
characteristics from single-particle breakage tests and its application
to ball-mill scale-up. Trans. IMM Sect. C, 97, C115-C124. The t10
appearance method.

[7] Napier-Munn, T. J., Morrell, S., Morrison, R. D., and Kojovic, T. (1996).
Mineral Comminution Circuits: Their Operation and Optimisation. JKMRC,
University of Queensland. Variable-exponent Whiten function, appearance
functions, and pendulum power.

[8] Rosin, P., and Rammler, E. (1933). The laws governing the fineness of
powdered coal. Journal of the Institute of Fuel, 7, 29-36.
Rosin-Rammler product distribution.

[9] Schuhmann, R. (1940). Principles of comminution, I: size distribution
and surface calculations. AIME Technical Publication No. 1189.
Gates-Gaudin-Schuhmann product distribution.

[10] von Rittinger, P. R. (1867). Lehrbuch der Aufbereitungskunde.
Ernst & Korn, Berlin. Rittinger energy law.

[11] Kick, F. (1885). Das Gesetz der proportionalen Widerstände und seine
Anwendung. Arthur Felix, Leipzig. Kick energy law.

[12] Bond, F. C. (1952). The third theory of comminution. Trans. AIME, 193,
484-494. Bond energy law and work index.

[13] Hukki, R. T. (1961). Proposal for a Solomonic settlement between the
theories of von Rittinger, Kick and Bond. Trans. AIME, 220, 403-408.
Applicable size ranges for the three energy laws.

[14] Zogg, M. (1993). Einführung in die Mechanische Verfahrenstechnik,
3rd ed. B. G. Teubner, Stuttgart, sec. 2.1.2.2. Coefficients that
match the Rittinger and Kick laws to Bond at their size boundaries.

[15] Upadhyay, R. K. (2025). Mining, Mineral Beneficiation, and Environment.
In Geology and Mineral Resources, 799-858. Springer Nature Singapore.
"""

__author__ = "Daison Yancy Caballero"

__all__ = [
    "CrusherSolidPSD",
    "CrusherSolidPSDData",
    "CrusherSolidPSDInitializer",
    "CrusherSolidPSDScaler",
]

import math

from pyomo.common.config import ConfigDict, ConfigValue, In
from pyomo.environ import (
    Expr_if,
    NonNegativeReals,
    Param,
    PositiveReals,
    UnitInterval,
    Var,
    units,
    value,
)

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_max

from prommis.comminution.core.psd_math import P80_PASSING, fractions_from_passing
from prommis.comminution.core.smooth_functions import EPS_SMOOTH
from prommis.comminution.core.unit_utils import MASS_FLOW_EPS
from prommis.comminution.functions import distributions, power_laws
from prommis.comminution.functions.power_laws import specific_energy
from prommis.comminution.unit_models import _size_reducer_common as sr
from prommis.comminution.unit_models import _solid_psd_common as spc
from prommis.comminution.unit_models._crusher.shared import (
    PSDMethod,
    DistributionShape,
    SelectionFunction,
    BreakageFunction,
    RecycleMode,
    CrusherPowerLaw,
    CrusherEquipment,
    CrusherStage,
    ProductPath,
    PreparedCrusherData,
    WORK_INDEX_BOUNDS,
    N_STRESS_EVENTS_MAX,
    MIN_SIZED_FEED_KG_S,
    _DISTRIBUTION_PATHS,
    _MIN_SHAPE_SIZE,
    _TABULAR_SUM_TOL,
    _shape_size_upper_bound,
)
from prommis.comminution.unit_models._crusher.config_processing import (
    validate_and_copy_configuration,
    prepare_crushing_data,
)
from prommis.comminution.unit_models._crusher.initializer import (
    CrusherSolidPSDInitializer,
)
from prommis.comminution.unit_models._crusher.scaler import CrusherSolidPSDScaler

_log = idaeslog.getLogger(__name__)

STAGE_VARIABLE_WARNING_THRESHOLD = 1_000_000
# Allow cumulative product passing to fall slightly below feed passing due to rounding.
_NO_COARSENING_TOL = 1e-12


# Product-P80 applicability ranges in cm: lower inclusive, upper exclusive
# (Upadhyay 2025). They support the warning-only check_applicability check.
CRUSHER_STAGE_DATA = {
    CrusherStage.primary: {"p80_lower_cm": 10.0, "p80_upper_cm": None},
    CrusherStage.secondary: {"p80_lower_cm": 2.0, "p80_upper_cm": 10.0},
    CrusherStage.tertiary: {"p80_lower_cm": 0.5, "p80_upper_cm": 2.0},
}


@declare_process_block_class("CrusherSolidPSD")
class CrusherSolidPSDData(UnitModelBlockData):
    """Steady-state crusher on the solid PSD property package."""

    CONFIG = ConfigDict()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc=(
                "This unit is quasi-steady; dynamic=True is not supported. "
                "**default** - False."
            ),
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Holdup construction flag - must be False",
            doc=(
                "This unit has no holdup terms; has_holdup=True is not supported. "
                "**default** - False."
            ),
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description=(
                "SolidPSD property parameter block. **default** - useDefault."
            ),
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigDict(
            implicit=True,
            description="Arguments for the property block. **default** - {}.",
        ),
    )
    CONFIG.declare(
        "psd_method",
        ConfigValue(
            default=PSDMethod.selection_breakage,
            domain=In(PSDMethod),
            description=(
                "Product PSD method: 'selection_breakage', 'tabular', or "
                "'distribution_function'. **default** - 'selection_breakage'."
            ),
        ),
    )
    CONFIG.declare(
        "selection_function",
        ConfigValue(
            default=SelectionFunction.whiten,
            domain=In(SelectionFunction),
            description=(
                "Selection function for 'selection_breakage': 'whiten', "
                "'austin', 'vogel_peukert', 'king', or 'user'. Classification "
                "recycle requires 'whiten'. **default** - 'whiten'."
            ),
        ),
    )
    CONFIG.declare(
        "breakage_function",
        ConfigValue(
            default=BreakageFunction.luckie_austin,
            domain=In(BreakageFunction),
            description=(
                "Breakage function for 'selection_breakage': 'luckie_austin', "
                "'reid_stewart', 'vogel', 't10_appearance', or 'user'. "
                "'t10_appearance' requires classification recycle. "
                "**default** - 'luckie_austin'."
            ),
        ),
    )
    CONFIG.declare(
        "recycle",
        ConfigValue(
            default=RecycleMode.classification_recycle,
            domain=In(RecycleMode),
            description=(
                "Selection-breakage circuit: 'single_pass' or "
                "'classification_recycle' (Whiten closed circuit; recrushes "
                "classified oversize). **default** - 'classification_recycle'."
            ),
        ),
    )
    CONFIG.declare(
        "power_law",
        ConfigValue(
            default=CrusherPowerLaw.bond,
            domain=In(CrusherPowerLaw),
            description=(
                "Power law: 'bond', 'rittinger', 'kick', or 'pendulum'. The first "
                "three use bond_work_index. 'pendulum' requires Whiten recycle, "
                "t10_appearance breakage, pendulum_power_params, and an "
                "'ecs_table' in the breakage parameters; do not set "
                "bond_work_index. **default** - 'bond'."
            ),
        ),
    )
    CONFIG.declare(
        "pendulum_power_params",
        ConfigValue(
            default=None,
            description=(
                "power_law='pendulum' only: {'power_factor': A (-), "
                "'no_load_power': P_no_load (kW)}, both finite and nonnegative. "
                "**default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "distribution_shape",
        ConfigValue(
            default=DistributionShape.rosin_rammler,
            domain=In(DistributionShape),
            description=(
                "Product PSD shape for 'distribution_function': 'rosin_rammler' "
                "or 'gates_gaudin_schuhmann'. **default** - 'rosin_rammler'."
            ),
        ),
    )
    CONFIG.declare(
        "crusher_equipment",
        ConfigValue(
            default=None,
            domain=In([None] + list(CrusherEquipment)),
            description=(
                "Optional equipment preset; sets the stage and supplies "
                "selection and breakage function defaults when available. "
                "**default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "crusher_stage",
        ConfigValue(
            default=None,
            domain=In([None] + list(CrusherStage)),
            description=(
                "Optional crushing stage for equipment consistency and "
                "product-P80 applicability checks. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "tabular_psd",
        ConfigValue(
            default=None,
            description=(
                "Shared tabular product fractions: N nonnegative values, finest "
                f"first, with sum within {_TABULAR_SUM_TOL:g} of one; normalized "
                "at build. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "tabular_psd_by_comp",
        ConfigValue(
            default=None,
            description=(
                "Tabular product fractions indexed by every sized solid component; "
                "each row has N values. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "target_product_p80",
        ConfigValue(
            default=None,
            description=(
                "Target product P80 in meters for 'distribution_function'; "
                "omit and fix power[t] to infer P80. After build, the "
                "target_product_p80 Param may be changed; the initializer rechecks it. "
                "**default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "shape_exponent",
        ConfigValue(
            default=None,
            description=(
                "Exponent controlling product distribution spread, at least "
                "1e-3; required for 'distribution_function'. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "n_stress_events",
        ConfigValue(
            default=1,
            description=(
                "Number of single-pass stress events; integer from 1 to "
                f"{N_STRESS_EVENTS_MAX}, excluding booleans. **default** - 1."
            ),
        ),
    )
    CONFIG.declare(
        "selection_params",
        ConfigValue(
            default=None,
            description=(
                "Shared selection-function inputs: a parameter dictionary, or N "
                "finest-first values for 'user'. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "selection_params_by_comp",
        ConfigValue(
            default=None,
            description=(
                "Selection inputs indexed by every sized solid component "
                "for single-pass selection-breakage. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "breakage_params",
        ConfigValue(
            default=None,
            description=(
                "Shared breakage-function inputs: a parameter dictionary, or an "
                "N-by-N matrix for 'user' (daughter rows, parent columns; "
                "finest first). **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "breakage_params_by_comp",
        ConfigValue(
            default=None,
            description=(
                "Breakage inputs indexed by every sized solid component; values are "
                "parameter dictionaries or user matrices. **default** - None."
            ),
        ),
    )
    CONFIG.declare(
        "bond_work_index",
        ConfigValue(
            default=14.0,
            description=(
                "Bond work index for Bond, Rittinger, and Kick power laws "
                f"(kWh/t); finite and within [{WORK_INDEX_BOUNDS[0]:g}, "
                f"{WORK_INDEX_BOUNDS[1]:g}]. Do not set explicitly for pendulum. "
                "**default** - 14.0."
            ),
        ),
    )

    default_initializer = CrusherSolidPSDInitializer
    default_scaler = CrusherSolidPSDScaler

    def build(self):
        """Validate configuration and build the crusher's states and equations."""
        super().build()
        self._get_property_package()
        pp = spc.require_solid_psd_package(self.config, "CrusherSolidPSD")
        self._resolve_configuration()
        time = self.flowsheet().time
        n_int = len(pp.size_interval_set)
        if self.product_path == ProductPath.repeated_stress:
            positions, expression_terms = self._stage_resource_estimate(
                self.n_stress_events, n_int, len(pp.sized_solid_list), len(time)
            )
            if positions > STAGE_VARIABLE_WARNING_THRESHOLD:
                _log.warning(
                    "CrusherSolidPSD repeated-stress build declares %s stage "
                    "variable entries, above the warning threshold of %s, and about "
                    "%s expression terms; construction may require substantial "
                    "time and memory.",
                    f"{positions:,}",
                    f"{STAGE_VARIABLE_WARNING_THRESHOLD:,}",
                    f"{expression_terms:,}",
                )
        sr.add_inlet_outlet_state_blocks(self)
        self.bond_work_index = Var(
            initialize=float(self.config.bond_work_index),
            bounds=WORK_INDEX_BOUNDS,
            units=units.kWh / units.metric_ton,
            doc="Bond work index (kWh/t); fixed at build.",
        )
        self.bond_work_index.fix()
        self.power = Var(
            time,
            units=units.kW,
            bounds=(0, None),
            initialize=1.0,
            doc="Crusher power (kW).",
        )
        # The crushing data live in mutable Params; the initializer seeds from them.
        prepared = self._prepare_crushing_data()
        sset = pp.size_interval_set
        sized = pp.sized_solid_list
        if prepared.tabular_fractions is not None:
            self.tabular_psd = Param(
                sized,
                sset,
                mutable=True,
                domain=UnitInterval,
                units=units.dimensionless,
                initialize={
                    (s, k): prepared.tabular_fractions[s][k]
                    for s in sized
                    for k in sset
                },
                doc=(
                    "Normalized tabulated product mass fraction by component and "
                    "size interval."
                ),
            )
        # Classification recycle uses the selection row only as classification.
        if prepared.selection is not None and prepared.classification is None:
            self.selection = Param(
                sized,
                sset,
                mutable=True,
                domain=UnitInterval,
                units=units.dimensionless,
                initialize={
                    (s, k): prepared.selection[s][k] for s in sized for k in sset
                },
                doc="Fraction selected for breakage by component and size interval.",
            )
        if prepared.breakage is not None:
            self.breakage = Param(
                sized,
                sset,
                sset,
                mutable=True,
                domain=UnitInterval,
                units=units.dimensionless,
                initialize={
                    (s, i, j): prepared.breakage[s][i][j]
                    for s in sized
                    for i in sset
                    for j in sset
                },
                doc=(
                    "Fraction of broken parent interval j entering daughter "
                    "interval i, by component."
                ),
            )
        if prepared.classification is not None:
            self.classification = Param(
                sset,
                mutable=True,
                domain=UnitInterval,
                units=units.dimensionless,
                initialize={k: prepared.classification[k] for k in sset},
                doc=(
                    "Fraction sent for breakage by size interval, shared across "
                    "components."
                ),
            )
        if prepared.ecs is not None:
            self.ecs = Param(
                sized,
                sset,
                mutable=True,
                domain=NonNegativeReals,
                units=units.kWh / units.metric_ton,
                initialize={(s, k): prepared.ecs[s][k] for s in sized for k in sset},
                doc="Specific comminution energy per component and interval (kWh/t).",
            )

        sr.add_passthrough_constraints(self)
        interior = list(sset)[1:]
        boundary = list(sset)[: n_int - 1]

        if self.product_path == ProductPath.tabular:
            self._build_tabular(time, sized, sset, interior, boundary)
        elif self.product_path in _DISTRIBUTION_PATHS:
            self._build_distribution(time, sized, sset, interior, boundary)
        elif self.product_path == ProductPath.single_stress:
            self._build_single_stress(time, sized, sset, interior)
        elif self.product_path == ProductPath.repeated_stress:
            self._build_repeated_stress(time, sized, sset, interior)
        elif self.product_path == ProductPath.classification_recycle:
            self._build_classification_recycle(time, sized, sset, interior)

        if self.product_path not in _DISTRIBUTION_PATHS:

            @self.Expression(
                time, doc="Aggregate sized-solids outlet P80 (package size units)."
            )
            def product_p80(b, t):
                return b.properties_out[t].percentile_size[P80_PASSING]

        self._build_power(time)

    def initialize(self, *args, **kwargs):
        """Reject the inherited initializer before it changes model state.

        The inherited method expects a ``control_volume``, which this crusher
        does not have. This legacy method could deactivate an attached block
        and then fail, leaving that block inactive.
        """
        raise NotImplementedError(
            f"{self.name}: CrusherSolidPSD does not support the legacy "
            "initialize() method. Use CrusherSolidPSDInitializer instead."
        )

    # Configuration and input validation

    def _resolve_configuration(self):
        """Check the options and settle the product path.

        Validation runs first because it writes equipment defaults into
        ``config``.
        """
        self.n_stress_events = validate_and_copy_configuration(self.config)
        self.product_path = self._select_product_path(self.config)

    def _select_product_path(self, config):
        """Return the product path selected by the PSD method and its options."""
        if config.psd_method == PSDMethod.tabular:
            return ProductPath.tabular
        if config.psd_method == PSDMethod.distribution_function:
            if config.target_product_p80 is not None:
                return ProductPath.distribution_size_spec
            return ProductPath.distribution_power_spec
        if config.recycle == RecycleMode.classification_recycle:
            return ProductPath.classification_recycle
        if self.n_stress_events == 1:
            return ProductPath.single_stress
        return ProductPath.repeated_stress

    # Crushing-data preparation and size calculations

    def _prepare_crushing_data(self):
        """Prepare validated crushing inputs for this unit's selected product path."""
        return prepare_crushing_data(self.config, self.product_path, self.name)

    def _crushing_data_from_params(self):
        """Return the crushing data the equations use, read from the Params.

        The Params hold the data built from the configuration; post-build edits
        are unsupported. The Param domains reject out-of-range values on
        assignment and on ``from_json``.
        """
        pp = self.config.property_package
        sized = list(pp.sized_solid_list)
        sset = list(pp.size_interval_set)

        def entry(param, key):
            number = value(param[key])
            if not math.isfinite(number):
                raise ConfigurationError(
                    f"{self.name}: {param.local_name}[{key}] = {number!r} must be "
                    "finite."
                )
            return number

        def rows(name):
            param = self.component(name)
            if param is None:
                return None
            return {s: tuple(entry(param, (s, k)) for k in sset) for s in sized}

        breakage = None
        if self.component("breakage") is not None:
            breakage = {
                s: tuple(
                    tuple(entry(self.breakage, (s, i, j)) for j in sset) for i in sset
                )
                for s in sized
            }
        classification = None
        if self.component("classification") is not None:
            classification = tuple(entry(self.classification, k) for k in sset)
        return PreparedCrusherData(
            rows("tabular_psd"),
            rows("selection"),
            breakage,
            classification,
            rows("ecs"),
        )

    def _shape_size_var(self):
        """Return the characteristic-size Var (``d_63`` or ``d_max``), or ``None``."""
        if self.product_path not in _DISTRIBUTION_PATHS:
            return None
        return getattr(
            self, distributions.SHAPE_SIZE_NAME[self.config.distribution_shape]
        )

    # Product and power equations

    def _build_tabular(self, time, sized, sset, interior, boundary):
        """Build tabular product equations and per-component no-coarsening checks."""

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Product mass flow from tabulated fractions; the finest interval "
                "enforces component mass balance."
            ),
        )
        def product_psd_eqn(b, t, s, k):
            if k == 0:
                return b._finest_product_balance(t, s, interior)
            return (
                b.properties_out[t].flow_mass_sized_comp_size[s, k]
                == b.tabular_psd[s, k] * b.properties_in[t].flow_mass_sized_comp[s]
            )

        @self.Constraint(
            time,
            sized,
            boundary,
            doc="Per-component no-coarsening at every interior boundary.",
        )
        def no_coarsening_comp_eqn(b, t, s, k):
            return (
                sum(
                    b.properties_out[t].flow_mass_sized_comp_size[s, j]
                    for j in range(k + 1)
                )
                - sum(
                    b.properties_in[t].flow_mass_sized_comp_size[s, j]
                    for j in range(k + 1)
                )
                >= 0
            )

    def _build_distribution(self, time, sized, sset, interior, boundary):
        """Build distribution product rows for a target P80 or specified power."""
        cfg = self.config
        pp = cfg.property_package
        shape = cfg.distribution_shape
        n = float(cfg.shape_exponent)
        n_int = len(sset)
        x_top = pp.size_edges_m[n_int]
        shape_size_ub = _shape_size_upper_bound(shape, x_top, n)
        name = distributions.SHAPE_SIZE_NAME[shape]
        setattr(
            self,
            name,
            Var(
                time,
                units=units.m,
                bounds=(_MIN_SHAPE_SIZE, shape_size_ub),
                initialize=max(_MIN_SHAPE_SIZE, 0.5 * shape_size_ub),
                doc=f"{distributions.SHAPES[shape].shape_size_doc} (m).",
            ),
        )
        shape_size_var = getattr(self, name)

        @self.Expression(
            time,
            pp.size_edge_index,
            doc=(
                "Analytic cumulative passing fraction at each mesh edge before "
                "normalization."
            ),
        )
        def shape_cdf(b, t, i):
            edge = units.convert(pp.size_edges[i], to_units=units.m)
            return distributions.shape_cdf(shape, edge, shape_size_var[t], n)

        @self.Expression(
            time,
            sset,
            doc=(
                "Product mass fraction by size interval, normalized at the top "
                "mesh edge."
            ),
        )
        def product_mass_frac(b, t, k):
            passing = [b.shape_cdf[t, i] for i in pp.size_edge_index]
            return fractions_from_passing(passing)[k]

        p80_factor = distributions.p80_ratio(shape, n)

        @self.Expression(
            time,
            doc=(
                "P80 of the analytic product distribution before mesh-top "
                "normalization (m)."
            ),
        )
        def product_p80(b, t):
            return shape_size_var[t] * p80_factor

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Product mass flow from distribution fractions; the finest "
                "interval enforces component mass balance."
            ),
        )
        def product_psd_eqn(b, t, s, k):
            if k == 0:
                return b._finest_product_balance(t, s, interior)
            return (
                b.properties_out[t].flow_mass_sized_comp_size[s, k]
                == b.product_mass_frac[t, k]
                * b.properties_in[t].flow_mass_sized_comp[s]
            )

        if self.product_path == ProductPath.distribution_size_spec:
            self.target_product_p80 = Param(
                initialize=float(cfg.target_product_p80),
                units=units.m,
                mutable=True,
                domain=PositiveReals,
                doc=(
                    "Target analytic product P80 (m); may be changed after build, "
                    "and the initializer rechecks it."
                ),
            )

            @self.Constraint(time, doc="Match analytic product P80 to its target.")
            def p80_shape_size_eqn(b, t):
                return b.product_p80[t] == b.target_product_p80

            top_edge = units.convert(pp.size_edges[n_int], to_units=units.m)
            # Every no-coarsening row compares against the feed P80, so hold it
            # in one variable rather than repeat the percentile expression.
            self.feed_p80 = Var(
                time,
                units=units.m,
                bounds=(0.0, x_top),
                initialize=0.5 * x_top,
                doc="Smoothed P80 of the sized feed (m).",
            )

            @self.Constraint(time, doc="Set feed_p80 from the inlet state.")
            def feed_p80_eqn(b, t):
                return b.feed_p80[t] == units.convert(
                    b.properties_in[t].percentile_size[P80_PASSING], to_units=units.m
                )

            # With positive sized feed, the power equation prevents analytic
            # product P80 from exceeding feed P80 on the power-spec path.

            @self.Constraint(
                time,
                boundary,
                doc=(
                    "Require analytic product P80 <= feed P80, or aggregate "
                    "cumulative outlet mass flow >= inlet, at each interior "
                    "boundary."
                ),
            )
            def no_coarsening_eqn(b, t, k):
                # Both margins are dimensionless, so the scaler uses a factor of 1.
                p80_margin = (b.feed_p80[t] - b.product_p80[t]) / top_edge
                cumulative_gain = sum(
                    b.properties_out[t].flow_mass_sized_comp_size[s, j]
                    - b.properties_in[t].flow_mass_sized_comp_size[s, j]
                    for s in sized
                    for j in range(k + 1)
                ) / (
                    b.properties_in[t].flow_mass_sized
                    + units.convert(
                        MASS_FLOW_EPS,
                        to_units=units.get_units(b.properties_in[t].flow_mass_sized),
                    )
                )
                return smooth_max(p80_margin, cumulative_gain, EPS_SMOOTH) >= 0

    def _build_single_stress(self, time, sized, sset, interior):
        """Build product equations for one selection and breakage event."""
        order = list(sset)

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Product mass flow after one stress event; the finest interval "
                "enforces component mass balance."
            ),
        )
        def product_psd_eqn(b, t, s, k):
            if k == 0:
                return b._finest_product_balance(t, s, interior)
            return b.properties_out[t].flow_mass_sized_comp_size[s, k] == sum(
                b._calculate_crushing_coefficient(s, k, i)
                * b.properties_in[t].flow_mass_sized_comp_size[s, i]
                for i in order
                if i >= k
            )

    def _build_repeated_stress(self, time, sized, sset, interior):
        """Build intermediate stress-event flows and the final product equations."""
        order = list(sset)
        n = self.n_stress_events
        stages = range(1, n)
        self.psd_stage = Var(
            time,
            stages,
            sized,
            sset,
            units=units.kg / units.s,
            bounds=(0, None),
            initialize=1e-3,
            doc=(
                "Mass flow after each intermediate stress event, by component and "
                "size interval (kg/s)."
            ),
        )

        def _incoming(b, t, r, s, i):
            if r == 1:
                return units.convert(
                    b.properties_in[t].flow_mass_sized_comp_size[s, i],
                    to_units=units.kg / units.s,
                )
            return b.psd_stage[t, r - 1, s, i]

        @self.Constraint(
            time,
            stages,
            sized,
            sset,
            doc=(
                "Mass flow after each intermediate stress event; the finest "
                "interval enforces component mass balance."
            ),
        )
        def psd_stage_eqn(b, t, r, s, k):
            if k == 0:
                return b.psd_stage[t, r, s, 0] == sum(
                    _incoming(b, t, r, s, i) for i in order
                ) - sum(b.psd_stage[t, r, s, j] for j in interior)
            return b.psd_stage[t, r, s, k] == sum(
                b._calculate_crushing_coefficient(s, k, i) * _incoming(b, t, r, s, i)
                for i in order
                if i >= k
            )

        last = n - 1

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Product mass flow after the final stress event; the finest "
                "interval enforces component mass balance."
            ),
        )
        def product_psd_eqn(b, t, s, k):
            if k == 0:
                return b._finest_product_balance(t, s, interior)
            return units.convert(
                b.properties_out[t].flow_mass_sized_comp_size[s, k],
                to_units=units.kg / units.s,
            ) == sum(
                b._calculate_crushing_coefficient(s, k, i) * b.psd_stage[t, last, s, i]
                for i in order
                if i >= k
            )

    @staticmethod
    def _stage_resource_estimate(n_stress_events, n_intervals, n_species, n_time):
        """Return repeated-stress stage positions and estimated expression terms.

        The estimate counts terms in the intermediate stage and final product
        equations for each component at each time point.
        """
        n, N = n_stress_events, n_intervals
        positions = (n - 1) * N * n_species * n_time
        expression_terms = (
            n_species * n_time * (n * N * (N - 1) // 2 + (n - 1) * (2 * N - 1) + N)
        )
        return positions, expression_terms

    def _build_classification_recycle(self, time, sized, sset, interior):
        """Build Whiten recycle-load and product equations for each component."""
        order = list(sset)
        self.recycle_load = Var(
            time,
            sized,
            sset,
            units=units.kg / units.s,
            bounds=(0, None),
            initialize=1e-3,
            doc=(
                "Internal classifier feed X, including fresh feed and crushed "
                "recycle, by component and size interval (kg/s)."
            ),
        )

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Whiten closed-circuit balance for the internal classifier feed "
                "by component and size interval."
            ),
        )
        def recycle_load_eqn(b, t, s, k):
            # Row k: X_k - sum_{j>=k} B_kj C_j X_j == F_k (diagonal inside the sum).
            return b.recycle_load[t, s, k] - sum(
                b.breakage[s, k, j] * b.classification[j] * b.recycle_load[t, s, j]
                for j in order
                if j >= k
            ) == units.convert(
                b.properties_in[t].flow_mass_sized_comp_size[s, k],
                to_units=units.kg / units.s,
            )

        @self.Constraint(
            time,
            sized,
            sset,
            doc=(
                "Whiten closed-circuit product flow by component and size interval; "
                "the finest interval enforces component mass balance."
            ),
        )
        def product_psd_eqn(b, t, s, k):
            if k == 0:
                return b._finest_product_balance(t, s, interior)
            return (
                units.convert(
                    b.properties_out[t].flow_mass_sized_comp_size[s, k],
                    to_units=units.kg / units.s,
                )
                == (1 - b.classification[k]) * b.recycle_load[t, s, k]
            )

    def _finest_product_balance(self, t, s, interior):
        """Constrain the finest outlet flow to conserve each component's mass."""
        return self.properties_out[t].flow_mass_sized_comp_size[
            s, 0
        ] == self.properties_in[t].flow_mass_sized_comp[s] - sum(
            self.properties_out[t].flow_mass_sized_comp_size[s, k] for k in interior
        )

    def _calculate_crushing_coefficient(self, s, k, i):
        """Return the symbolic one-event coefficient :math:`A_{s,k,i}` of the
        module docstring, from parent interval ``i`` to daughter interval ``k``
        of component ``s``.

        Symbolic twin of ``_numeric_crushing_coefficients``; keep the two aligned.
        """
        broken = self.breakage[s, k, i] * self.selection[s, i]
        if k == i:
            return 1 - self.selection[s, i] + broken
        return broken

    def _build_power(self, time):
        """Build the selected crusher power equation."""
        if self.config.power_law == CrusherPowerLaw.pendulum:
            self._build_pendulum_power(time)
            return
        # With positive sized feed, do not floor energy on the power-spec path.
        # Otherwise, zero power could accept a product P80 above the feed P80.
        floored = self.product_path != ProductPath.distribution_power_spec
        # Read units from the mesh edges, the sized flows and the work index,
        # then use constant factors. units.convert would scan the large
        # percentile expressions for units.
        size_units = self.config.property_package.size_edges.get_units()
        inlet_state = self.properties_in[time.first()]
        flow_units = inlet_state.flow_mass_sized_comp_size.get_units()
        work_index_units = units.get_units(self.bond_work_index)
        zero_specific_energy = 0.0 * work_index_units
        to_kw = units.convert_value(
            1.0, from_units=flow_units * work_index_units, to_units=units.kW
        )

        def microns(size, from_units):
            factor = units.convert_value(1.0, from_units=from_units, to_units=units.um)
            return size * (factor / from_units)

        @self.Constraint(
            time,
            doc="Crusher power from sized feed mass flow and specific energy (kW).",
        )
        def power_eqn(b, t):
            inlet = b.properties_in[t]
            if b.product_path in _DISTRIBUTION_PATHS:
                # The analytic product P80 has the shape size's units, which
                # need not match the package size units of the percentiles.
                product_units = units.get_units(b.product_p80[t])
            else:
                product_units = size_units
            e_raw = specific_energy(
                b.config.power_law,
                b.bond_work_index,
                microns(b.product_p80[t], product_units),
                microns(inlet.percentile_size[P80_PASSING], size_units),
                sizes_are_microns=True,
            )
            e_spec = (
                Expr_if(
                    IF=e_raw >= zero_specific_energy,
                    THEN=e_raw,
                    ELSE=zero_specific_energy,
                )
                if floored
                else e_raw
            )
            return b.power[t] == (
                inlet.flow_mass_sized / flow_units * e_spec / work_index_units
            ) * (to_kw * units.kW)

    def _build_pendulum_power(self, time):
        """Build pendulum power variables, expression, and constraint."""
        params = self.config.pendulum_power_params
        pp = self.config.property_package
        self.pendulum_power_factor = Var(
            initialize=float(params["power_factor"]),
            bounds=(0, None),
            units=units.dimensionless,
            doc="Dimensionless pendulum-to-crusher power factor A; fixed at build.",
        )
        self.pendulum_power_factor.fix()
        self.no_load_power = Var(
            initialize=float(params["no_load_power"]),
            bounds=(0, None),
            units=units.kW,
            doc="No-load (idle) power (kW); fixed at build.",
        )
        self.no_load_power.fix()

        @self.Expression(
            time, doc="Pendulum power, sum of Ecs * C * X converted to kW."
        )
        def pendulum_power(b, t):
            return power_laws.pendulum_power(
                (b.ecs[s, k], b.classification[k], b.recycle_load[t, s, k])
                for s in pp.sized_solid_list
                for k in pp.size_interval_set
            )

        @self.Constraint(
            time,
            doc="Crusher power P = A*P_p + P_0 from scaled pendulum power plus "
            "no-load power (kW) (Napier-Munn et al., 1996, ch. 11).",
        )
        def power_eqn(b, t):
            return (
                b.power[t]
                == b.pendulum_power_factor * b.pendulum_power[t] + b.no_load_power
            )

    # Required feed and product checks

    def validate_sized_feed_flow(self):
        """Require total sized feed to be at least ``MIN_SIZED_FEED_KG_S`` (1e-6 kg/s).

        Call after ``validate_feed`` has checked that all feed flow values are
        set and finite.
        """
        for t in self.flowsheet().time:
            m_sized = value(
                units.convert(
                    self.properties_in[t].flow_mass_sized, to_units=units.kg / units.s
                )
            )
            if m_sized < MIN_SIZED_FEED_KG_S:
                raise ConfigurationError(
                    f"{self.name}: sized feed mass flow at t = {t} ({m_sized!r} "
                    f"kg/s) is below the minimum {MIN_SIZED_FEED_KG_S:.6g} kg/s."
                )

    def _check_no_coarsening(self, product_fractions, feed_rows, n_int, at_time=None):
        """Require cumulative product fractions to be no lower than the feed's.

        Inputs map components to full rows; ``None`` denotes an aggregate row.
        Compare positive-feed rows within ``_NO_COARSENING_TOL`` and skip
        zero-feed rows. Supplied data are unchanged. The tolerance applies to
        fractions, so the allowed mass-flow difference grows with the feed rate.
        """
        remedy = (
            "change the table or the feed; the table is never altered."
            if self.product_path == ProductPath.tabular
            else "review the prescribed distribution or the feed; supplied data are not altered."
        )
        when = "" if at_time is None else f" at t = {at_time}"
        for mineral, feed in feed_rows.items():
            total = math.fsum(feed)
            if total == 0:
                continue
            fractions = product_fractions[mineral]
            normalized_feed = [entry / total for entry in feed]
            for k in range(n_int - 1):
                cum_product = math.fsum(fractions[: k + 1])
                cum_feed = math.fsum(normalized_feed[: k + 1])
                difference = cum_product - cum_feed
                if difference < -_NO_COARSENING_TOL:
                    context = (
                        "aggregate" if mineral is None else f"component {mineral!r}"
                    )
                    prefix = f"{self.name}: no-coarsening{when} for {context}"
                    raise ConfigurationError(
                        f"{prefix}: cumulative product fraction {cum_product!r} "
                        f"is below cumulative feed fraction {cum_feed!r} at boundary "
                        f"{k} by {-difference:.3e}, exceeding the tolerance "
                        f"{_NO_COARSENING_TOL:.3e}; {remedy}"
                    )

    # Optional checks

    def check_applicability(self, emit_warning=True):
        """Return messages when product P80 is unavailable or outside the stage range.

        The lower limit is inclusive and the upper limit exclusive. Log messages
        only when ``emit_warning`` is true; return an empty list when no stage is
        set or no issue is found.
        """
        stage = self.config.crusher_stage
        if stage is None:
            return []
        data = CRUSHER_STAGE_DATA[stage]
        lower = data["p80_lower_cm"]
        upper = data["p80_upper_cm"]
        messages = []
        for t in self.flowsheet().time:
            p80_cm = value(
                units.convert(self.product_p80[t], units.cm), exception=False
            )
            if p80_cm is None:
                messages.append(
                    f"{self.name}: product P80 at t = {t} has no value; "
                    "the crusher-stage range "
                    "was not checked."
                )
                continue
            if (lower is not None and p80_cm < lower) or (
                upper is not None and p80_cm >= upper
            ):
                rng = (
                    f"[{lower}, {upper}) cm" if upper is not None else f">= {lower} cm"
                )
                messages.append(
                    f"{self.name}: at t = {t}, crusher_stage='{stage}' has product "
                    f"P80 = {p80_cm:.4g} cm outside its applicable range {rng}."
                )
        return spc.report_applicability_warnings(
            _log, "CrusherSolidPSD", messages, emit_warning
        )

    # Reporting

    def _get_performance_contents(self, time_point=0):
        """Return crusher power, product P80, and any pendulum power for reporting."""
        contents = {
            "vars": {"Power": self.power[time_point]},
            "exprs": {"Product P80": self.product_p80[time_point]},
        }
        if hasattr(self, "pendulum_power"):
            contents["exprs"]["Pendulum power"] = self.pendulum_power[time_point]
        return contents
