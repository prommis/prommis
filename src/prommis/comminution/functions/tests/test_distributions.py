#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Tests for the assumed-shape distribution forms and their discretization."""

import pytest

from prommis.comminution.core.psd_math import cumulative_passing, size_at_passing
from prommis.comminution.core.size_mesh import geometric_series
from prommis.comminution.functions.distributions import (
    distribution_fractions,
)


@pytest.mark.unit
def test_ggs_happy_path_sum_and_p80():
    edges = list(geometric_series(top=40e-3, bottom=0.1e-3, ratio=400 ** (1.0 / 60)))
    d_max, n = 20e-3, 2.0
    fracs = distribution_fractions("gates_gaudin_schuhmann", edges, d_max, n)
    assert sum(fracs) == pytest.approx(1.0, rel=1e-12)
    cum = cumulative_passing(fracs)
    p80 = size_at_passing(edges, cum, 0.8)
    analytic = d_max * 0.8 ** (1.0 / n)
    assert p80 == pytest.approx(analytic, rel=1e-3)
