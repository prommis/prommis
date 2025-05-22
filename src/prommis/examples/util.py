#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Utility functions for the examples
"""

import logging
from typing import Union


class DropAll(logging.Filter):
    def filter(self, record):
        return False


def silence_log(obj: Union[logging.Logger, logging.LoggerAdapter]):
    log = getattr(obj, "logger", obj)
    log.addFilter(DropAll())
