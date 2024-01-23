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
