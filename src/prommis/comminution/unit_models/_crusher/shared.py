#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Shared crusher types, numerical limits, and helper functions."""

__author__ = "Daison Yancy Caballero"

import math
from collections import namedtuple

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum

from prommis.comminution.functions import distributions
from prommis.comminution.functions.power_laws import invert_specific_energy

WORK_INDEX_BOUNDS = (1.0, 150.0)  # kWh/t
N_STRESS_EVENTS_MAX = 1000
MIN_SIZED_FEED_KG_S = 1e-6
_TABULAR_SUM_TOL = 1e-5

# For valid Whiten inputs, this floor bounds internal classifier flow
# relative to fresh sized feed. It does not bound solver residuals.
_RECYCLE_FLOOR_PER_INTERVAL = 1e-5
# Minimum characteristic size (m): d_63 for Rosin-Rammler,
# d_max for Gates-Gaudin-Schuhmann.
_MIN_SHAPE_SIZE = 1e-9


class PSDMethod(StrEnum):
    """Product-PSD methods."""

    distribution_function = "distribution_function"
    tabular = "tabular"
    selection_breakage = "selection_breakage"


class DistributionShape(StrEnum):
    """Product-PSD distribution shapes."""

    rosin_rammler = "rosin_rammler"
    gates_gaudin_schuhmann = "gates_gaudin_schuhmann"


class SelectionFunction(StrEnum):
    """Crusher selection functions, including Whiten classification."""

    whiten = "whiten"
    austin = "austin"
    vogel_peukert = "vogel_peukert"
    king = "king"
    user = "user"


class BreakageFunction(StrEnum):
    """Breakage (daughter-distribution) functions."""

    luckie_austin = "luckie_austin"
    reid_stewart = "reid_stewart"
    vogel = "vogel"
    t10_appearance = "t10_appearance"
    user = "user"


class RecycleMode(StrEnum):
    """Crusher circuit options."""

    single_pass = "single_pass"
    classification_recycle = "classification_recycle"


class CrusherPowerLaw(StrEnum):
    """Crusher power models."""

    bond = "bond"
    rittinger = "rittinger"
    kick = "kick"
    pendulum = "pendulum"


class CrusherEquipment(StrEnum):
    """Equipment presets."""

    jaw = "jaw"
    gyratory = "gyratory"
    cone = "cone"
    roll1 = "roll1"
    roll2 = "roll2"
    short_head_cone = "short_head_cone"
    hammer_mill = "hammer_mill"


class CrusherStage(StrEnum):
    """Crushing stages."""

    primary = "primary"
    secondary = "secondary"
    tertiary = "tertiary"


class ProductPath(StrEnum):
    """Calculation routes selected from the product method and its options."""

    tabular = "tabular"
    distribution_size_spec = "distribution_size_spec"
    distribution_power_spec = "distribution_power_spec"
    single_stress = "single_stress"
    repeated_stress = "repeated_stress"
    classification_recycle = "classification_recycle"


_DISTRIBUTION_PATHS = (
    ProductPath.distribution_size_spec,
    ProductPath.distribution_power_spec,
)


# Crushing data by component s, interval k (finest first); fields unused by the
# product path are None. tabular_fractions[s][k] and selection[s][k] are
# fractions; breakage[s][i][j] is the fraction of broken parent j entering
# daughter i; classification[k] is shared by all components; ecs[s][k] is in kWh/t.
PreparedCrusherData = namedtuple(
    "PreparedCrusherData",
    ["tabular_fractions", "selection", "breakage", "classification", "ecs"],
)


def _shape_size_upper_bound(shape, x_top, n):
    """Return the mesh-based upper limit for a distribution's characteristic size.

    Raise ``ConfigurationError`` if the limit cannot be calculated, is
    nonfinite, or is below ``_MIN_SHAPE_SIZE`` (1e-9 m).
    """
    context = (
        f"{shape} exponent {n!r}, upper limit for characteristic size "
        f"at mesh top {x_top!r} m"
    )
    try:
        upper = distributions.shape_size_upper_bound(shape, x_top, n)
    except (ArithmeticError, ValueError) as exc:
        raise ConfigurationError(f"{context} could not be calculated: {exc}") from exc
    if not math.isfinite(upper) or upper < _MIN_SHAPE_SIZE:
        raise ConfigurationError(
            f"{context} must be finite and at least {_MIN_SHAPE_SIZE} m "
            f"(got {upper!r})."
        )
    return upper


# Explain each shape's mesh-based upper limit in configuration errors.
_SHAPE_SIZE_LIMIT_RULE = {
    DistributionShape.rosin_rammler: (
        "a larger d_63 gives Rosin-Rammler Q(x_N) below "
        f"{distributions.RR_VALIDITY_THRESHOLD}"
    ),
    DistributionShape.gates_gaudin_schuhmann: (
        "Gates-Gaudin-Schuhmann d_max must be <= the top mesh edge"
    ),
}


def _p80_from_power(power_law, p_kw, m_sized, work_index, f80_live):
    """Return the analytic product P80 (m) implied by specified power.

    ``p_kw`` is in kW, ``m_sized`` is sized feed mass flow in kg/s,
    ``work_index`` is in kWh/t, and ``f80_live`` is feed P80 in m.
    Zero power returns feed P80. Otherwise, with positive sized feed,
    convert power to specific energy and invert the selected power law.
    """
    if p_kw == 0.0:
        return f80_live
    return invert_specific_energy(
        power_law, p_kw / (m_sized * 3.6), work_index, f80_live
    )


def _close_finest(row, total=1.0):
    """Return ``row`` with its finest entry set to ``total`` minus the coarser
    entries."""
    return [math.fsum([total, *(-q for q in row[1:])]), *row[1:]]


def _numeric_crushing_coefficients(selection, breakage, n_int):
    """Return numeric one-event transfer coefficients, indexed
    ``[daughter][parent]``.

    This function does not adjust the finest row to make each column sum to one.
    Keep it aligned with the symbolic ``_calculate_crushing_coefficient``.
    """
    return [
        [
            (1 if i == j else 0) * (1.0 - selection[j]) + breakage[i][j] * selection[j]
            for j in range(n_int)
        ]
        for i in range(n_int)
    ]


def validate_recycle_denominators(breakage, classification, minerals, *, unit_name):
    """Require each recycle denominator ``1 - B[j][j] C[j]`` to be at least
    the number of size intervals times ``_RECYCLE_FLOOR_PER_INTERVAL``.

    Raise ``ConfigurationError`` if a denominator is too close to zero for
    the recycle calculation. Callers must provide finite, validated inputs.
    """
    n = len(classification)
    floor = _RECYCLE_FLOOR_PER_INTERVAL * n
    for j in range(n):
        delta = 1.0 - breakage[j][j] * classification[j]
        if delta < floor:
            shared = ", ".join(repr(s) for s in minerals)
            raise ConfigurationError(
                f"{unit_name}: Whiten recycle for {shared}: 1 - b_jj C_j = "
                f"{delta!r} at interval {j} must be at least {floor!r}."
            )
