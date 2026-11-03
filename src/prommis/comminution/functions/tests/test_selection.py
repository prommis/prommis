#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Reference values for selection and classification functions."""

import pytest

from prommis.comminution.functions.selection import whiten_selection


@pytest.mark.unit
def test_whiten_anchor():
    k1, k2, k3 = 0.200, 0.250, 2.3
    assert whiten_selection(0.200, k1, k2, k3) == 0.0
    assert whiten_selection(0.250, k1, k2, k3) == 1.0
    # C(0.225) = 1 - 0.5^2.3
    assert whiten_selection(0.225, k1, k2, k3) == pytest.approx(
        1.0 - 0.5**2.3, rel=1e-12
    )
