#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Choice of Objective Function Code
=================================

Author: Chris Laliwala
"""

from enum import Enum


class ObjectiveFunctionChoice(Enum):
    """
    Enum of supported objective functions.

    Note:
        The element-flow (bilinear stream-mixing) superstructure currently implements
        only ``NET_PRESENT_VALUE``. ``COST_OF_RECOVERY`` is reserved for parity with the
        fixed-yield PrOMMiS superstructure and raises ``NotImplementedError`` if selected.
    """

    NET_PRESENT_VALUE = 1
    COST_OF_RECOVERY = 2
