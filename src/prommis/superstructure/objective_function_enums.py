#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Superstructure Code
===================

Author: Chris Laliwala
"""

from enum import Enum


class ObjectiveFunction(Enum):
    """
    Enum of supported objective functions
    """

    NET_PRESENT_VALUE = 1
    COST_OF_RECOVERY = 2
