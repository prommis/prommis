#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for triangular population-balance rows and back substitution."""

import pytest


from prommis.comminution.core.triangular_block import (
    back_substitution,
    back_substitution_recycle,
)

# Shared bins and right-hand side for the recycle and explicit-diagonal cases.
_IDX = [0, 1, 2]
_RHS = {0: 1.0, 1: 2.0, 2: 3.0}


@pytest.mark.unit
def test_back_substitution_recycle_matches_hand_solution():
    # 3-bin recycle coupling c_ij = b_ij * C_j with the diagonal c_ii included;
    # each breakage column sums to one
    idx = _IDX
    brk = {(0, 0): 1.0, (0, 1): 0.3, (1, 1): 0.7, (0, 2): 0.2, (1, 2): 0.3, (2, 2): 0.5}
    cls = {0: 0.0, 1: 0.4, 2: 0.9}
    feed = _RHS
    coup = lambda i, j: brk.get((i, j), 0.0) * cls[j]

    # Hand solution of the recycle balance, including 1 - b_ii C_i:
    # x2 = 3 / 0.55, x1 = (2 + 0.27 x2) / 0.72 and x0 = 1 + 0.12 x1 + 0.18 x2.
    seed = back_substitution_recycle(idx, coup, lambda k: feed[k])
    hand = {0: 169.0 / 66.0, 1: 955.0 / 198.0, 2: 60.0 / 11.0}
    for i in idx:
        assert seed[i] == pytest.approx(hand[i], rel=1e-14)


_DIAG = {0: 0.9, 1: 0.8, 2: 1.0}
_COUP = {(0, 1): 0.2, (0, 2): 0.1, (1, 2): 0.3}


def _diag(i):
    return _DIAG[i]


def _coup(i, j):
    return _COUP.get((i, j), 0.0)


def _rhs(i):
    return _RHS[i]


@pytest.mark.unit
def test_back_substitution_hand_values():
    x = back_substitution(_IDX, _diag, _coup, _rhs)
    assert x[2] == pytest.approx(3.0, rel=1e-14)
    assert x[1] == pytest.approx(3.625, rel=1e-14)
    assert x[0] == pytest.approx(2.25, rel=1e-14)

    # Both signs are checked against the absolute diagonal floor of 1e-8.
    for small in (5e-9, -5e-9, 0.0):
        with pytest.raises(ValueError, match="near-zero diagonal at interval 1"):
            back_substitution(_IDX, lambda i: small if i == 1 else 1.0, _coup, _rhs)
    for small in (2e-8, -2e-8):
        x = back_substitution(_IDX, lambda i: small if i == 1 else 1.0, _coup, _rhs)
        assert x[1] == pytest.approx((2.0 + 0.3 * 3.0) / small, rel=1e-14)
    # Reject nonfinite diagonals before they propagate through the calculation.
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ValueError, match="non-finite diagonal at interval 1"):
            back_substitution(_IDX, lambda i: bad if i == 1 else 1.0, _coup, _rhs)
