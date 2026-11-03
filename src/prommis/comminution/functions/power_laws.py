#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Specific-energy laws and power correlations.

Bond, Rittinger, and Kick return mass-specific energy. They are cases of one
differential energy-size law with different size exponents (Hoeffl, 1986).
Work indices use ``W*h/kg`` (numerically equal to kWh/metric ton), and size
terms use the dimensionless magnitudes of sizes in microns. Unit models
multiply the result by solids mass flow to obtain power.

The Rittinger and Kick coefficients approximately match Bond's local size slope
at the crossover sizes ``RITTINGER_BOUNDARY_UM`` and ``KICK_BOUNDARY_UM``
(Zogg, 1993). ``APPLICABLE_RANGE_M`` holds the size ranges used for warnings:
Rittinger below the first crossover, Bond between the two, and Kick above the
second.

The pendulum crusher helpers calculate ``P_p = sum(Ecs * C * X)`` and
``P = A * P_p + P_0``. Here ``X`` is incoming classifier mass flow, and
``C * X`` is the flow selected for breakage. With ``Ecs`` in kWh/t and
``X`` in kg/s, multiplying the sum by 3.6 gives kW.
:func:`pendulum_power` relies on Pyomo unit
conversion; :func:`numeric_pendulum_power` takes ``Ecs`` in kWh/t and ``X``
in kg/s. These helpers return power rather than mass-specific energy.

References:

[1] von Rittinger, P. R. (1867). Lehrbuch der Aufbereitungskunde.
    Ernst & Korn, Berlin.

[2] Kick, F. (1885). Das Gesetz der proportionalen Widerstaende und seine
    Anwendung. Arthur Felix, Leipzig.

[3] Bond, F. C. (1952). The third theory of comminution. Trans. AIME,
    193, 484-494.

[4] Hukki, R. T. (1961). Proposal for a Solomonic settlement between the
    theories of von Rittinger, Kick and Bond. Trans. AIME, 220, 403-408.

[5] Hoeffl, K. (1986). Zerkleinerungs- und Klassiermaschinen.
    Springer, Berlin.

[6] Napier-Munn, T. J., Morrell, S., Morrison, R. D., & Kojovic, T. (1996).
    Mineral Comminution Circuits: Their Operation and Optimisation. JKMRC,
    Brisbane. Chapter 11: the pendulum crusher power model.

[7] Zogg, M. (1993). Einfuehrung in die Mechanische Verfahrenstechnik
    (3rd ed.). B. G. Teubner, Stuttgart.
"""

__author__ = "Daison Yancy Caballero"

import math

from pyomo.environ import log10
from pyomo.environ import units as pyunits
from pyomo.environ import value

# Bond-law crossover sizes are micron numbers used by the energy formulas.
RITTINGER_BOUNDARY_UM = 50.0
KICK_BOUNDARY_UM = 50000.0
# Zogg's factor (eq. 2.14).
_ZOGG_KICK_FACTOR = 1.151
# Suggested size ranges in meters; used for warnings, not hard limits (Hukki, 1961).
APPLICABLE_RANGE_M = {
    "rittinger": (0.0, RITTINGER_BOUNDARY_UM / 1e6),
    "bond": (RITTINGER_BOUNDARY_UM / 1e6, KICK_BOUNDARY_UM / 1e6),
    "kick": (KICK_BOUNDARY_UM / 1e6, None),
}


def _magnitude(quantity, unit):
    """Dimensionless magnitude of ``quantity`` expressed in ``unit``."""
    return pyunits.convert(quantity, to_units=unit) / unit


def _rittinger_coeff(bond_work_index):
    """Return ``C_R = 0.5*(10*WI)*sqrt(RITTINGER_BOUNDARY_UM)`` (micron numbers).

    This is Zogg's relation (2.15); it matches Bond's local size slope at the
    lower crossover. ``bond_work_index`` may be numeric or Pyomo object.
    """
    return 0.5 * (10 * bond_work_index) * RITTINGER_BOUNDARY_UM ** (0.5)


def _kick_coeff(bond_work_index):
    """Return ``C_K = _ZOGG_KICK_FACTOR*(10*WI)/sqrt(KICK_BOUNDARY_UM)``.

    Sizes are micron numbers. This is Zogg's relation (2.14); it approximately
    matches Bond's local size slope at the upper crossover.
    """
    return _ZOGG_KICK_FACTOR * (10 * bond_work_index) * KICK_BOUNDARY_UM ** (-0.5)


def specific_energy(power_law, bond_work_index, p80, f80, *, sizes_are_microns=False):
    """Return mass-specific energy from the selected Bond, Rittinger, or Kick law.

    ``bond_work_index`` has specific-energy units (``W*h/kg`` = kWh/tonne).
    By default, ``p80`` and ``f80`` are lengths with Pyomo units. With
    ``sizes_are_microns=True``, they are dimensionless micron magnitudes.
    The crusher needs the latter to avoid walking its large percentile
    expressions during unit conversion.

    The laws use sizes in microns:

    - Bond (1952): ``10*WI*(1/sqrt(P80_um) - 1/sqrt(F80_um))``.
    - Rittinger (1867): ``C_R*(1/P80_um - 1/F80_um)`` for fine sizes;
      :func:`_rittinger_coeff` supplies ``C_R``.
    - Kick (1885): ``C_K*log10(F80_um/P80_um)`` for coarse sizes;
      :func:`_kick_coeff` supplies ``C_K``.

    Return an expression in the units of ``bond_work_index``.
    """
    if sizes_are_microns:
        p80_um, f80_um = p80, f80
    else:
        p80_um = _magnitude(p80, pyunits.um)
        f80_um = _magnitude(f80, pyunits.um)
    if power_law == "bond":
        return 10 * bond_work_index * (1 / p80_um**0.5 - 1 / f80_um**0.5)
    if power_law == "rittinger":
        return _rittinger_coeff(bond_work_index) * (1 / p80_um - 1 / f80_um)
    if power_law == "kick":
        return _kick_coeff(bond_work_index) * log10(f80_um / p80_um)
    raise NotImplementedError(f"power_law {power_law!r} is not implemented")


def numeric_specific_energy(power_law, bond_work_index, p80, f80):
    """Evaluate :func:`specific_energy` as a float.

    ``bond_work_index`` is a number in ``W*h/kg``; ``p80``/``f80`` are lengths in meters
    (converted to the micron-number basis here, matching the Pyomo forms).
    Returns ``e_spec`` as a float in ``W*h/kg``.
    """
    return value(
        specific_energy(
            power_law,
            bond_work_index * (pyunits.W * pyunits.hr / pyunits.kg),
            p80 * 1e6 * pyunits.um,
            f80 * 1e6 * pyunits.um,
        )
    )


def invert_specific_energy(power_law, e_spec, bond_work_index, f80):
    """Closed-form product P80 (meters) from a specific energy ``e_spec``.

    Inverse of :func:`numeric_specific_energy`: ``e_spec`` in ``W*h/kg``,
    ``bond_work_index`` in ``W*h/kg``, ``f80`` in meters.  Shares the Rittinger/Kick
    coefficients with the forward forms.

    ``e_spec`` must be finite and nonnegative. The inverse assumes positive,
    finite ``bond_work_index`` and ``f80``; under these conditions, the returned
    P80 is no larger than ``f80``.

    Raises:
        ValueError: if ``e_spec`` is negative or not finite.
        NotImplementedError: for an unknown ``power_law``.
    """
    if not (math.isfinite(e_spec) and e_spec >= 0.0):
        raise ValueError(
            f"invert_specific_energy: e_spec must be finite and >= 0 (got {e_spec!r})."
        )
    f80_um = f80 * 1e6
    if power_law == "bond":
        inv_sqrt = e_spec / (10.0 * bond_work_index) + 1.0 / math.sqrt(f80_um)
        return (1.0 / inv_sqrt**2) * 1e-6
    if power_law == "rittinger":
        inv = e_spec / _rittinger_coeff(bond_work_index) + 1.0 / f80_um
        return (1.0 / inv) * 1e-6
    if power_law == "kick":
        return f80_um * 10.0 ** (-e_spec / _kick_coeff(bond_work_index)) * 1e-6
    raise NotImplementedError(f"power_law {power_law!r} inversion is not implemented")


def pendulum_power(terms):
    """Pendulum power ``P_p = sum Ecs * C * X`` in kW.

    Napier-Munn et al. (1996, ch. 11).

    ``terms`` yields ``(ecs, classification, load)`` triples, one per component
    and size interval. ``load`` is the incoming classifier mass flow ``X``;
    ``classification * load`` is the flow selected for breakage.
    """
    return pyunits.convert(
        sum(ecs * classification * load for ecs, classification, load in terms),
        to_units=pyunits.kW,
    )


def numeric_pendulum_power(terms):
    """Evaluate pendulum power as a float in kW.

    ``ecs`` is in kWh/metric ton and incoming ``load`` is in kg/s.
    """
    return value(
        pendulum_power(
            (
                ecs * (pyunits.kWh / pyunits.metric_ton),
                classification,
                load * pyunits.kg / pyunits.s,
            )
            for ecs, classification, load in terms
        )
    )
