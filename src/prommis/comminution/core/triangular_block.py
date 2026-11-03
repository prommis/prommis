#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Triangular population-balance equations and numeric back-substitution.

Each row is the balance equation for one size interval. Callers supply bins
from finest to coarsest. The seeders solve the explicit-diagonal form,
``d_i*x_i - sum_{j>i} c_ij*x_j = rhs_i``, or the recycle form,
``x_i - sum_{j>=i} c_ij*x_j = rhs_i``, from coarsest to finest.
"""

__author__ = "Daison Yancy Caballero"

import math

# Minimum diagonal magnitude for back substitution and row scaling.
DIAGONAL_TOL = 1e-8


def back_substitution(index_list, diagonal_value, coupling_value, rhs_value):
    """Solve diagonal-explicit rows by float back substitution, coarsest first.

    Return ``{i: x_i}``; raise ``ValueError`` for a non-finite or near-zero diagonal.
    """
    order = list(index_list)
    out = {}
    for i in reversed(order):
        d = float(diagonal_value(i))
        if not math.isfinite(d):
            raise ValueError(
                f"triangular system has a non-finite diagonal at interval "
                f"{i}: diagonal = {d}."
            )
        if abs(d) < DIAGONAL_TOL:
            raise ValueError(
                f"triangular system has a near-zero diagonal at interval {i}: "
                f"|diagonal| = {abs(d):.3e} < {DIAGONAL_TOL:.3e} "
                f"(diagonal = {d:.3e})."
            )
        coupled = sum(float(coupling_value(i, j)) * out[j] for j in order if j > i)
        out[i] = (float(rhs_value(i)) + coupled) / d
    return out


def back_substitution_recycle(index_list, coupling_value, rhs_value):
    """Solve recycle rows by float back substitution (diagonal ``1 - c_ii``)."""
    return back_substitution(
        index_list,
        lambda i: 1.0 - float(coupling_value(i, i)),
        coupling_value,
        rhs_value,
    )
