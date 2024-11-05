# coding: utf-8
"""
Boilerplate for usable and pretty logging.

If Colorama is installed, add some color to log messages.
"""
import argparse
import logging

try:
    from colorama import Fore, Style
except ImportError:

    class Fore:
        RED, YELLOW, CYAN, GREEN, RESET = "", "", "", "", ""

    class Style:
        DIM, RESET_ALL = "", ""


__author__ = "Dan Gunter (LBNL)"


class ColorFormatter(logging.Formatter):
    """Custom formatter that prints different levels in different colors."""

    colors = {
        logging.DEBUG: "",
        logging.INFO: Fore.GREEN,
        logging.WARNING: Fore.YELLOW,
        logging.ERROR: Fore.RED,
        logging.CRITICAL: Fore.RED,
    }
    _formatters = None

    @classmethod
    def _get_formatter(cls, levelno: int) -> logging.Formatter:
        if cls._formatters is None:
            formats = {
                k: color
                + "[{levelname}] "
                + Style.DIM
                + "{asctime} ({name}) "
                + Style.RESET_ALL
                + color
                + "{message}"
                + Fore.RESET
                for k, color in cls.colors.items()
            }
            cls._formatters = {
                k: logging.Formatter(v, style="{") for k, v in formats.items()
            }
        return cls._formatters[levelno]

    def format(self, record):
        formatter = self._get_formatter(record.levelno)
        return formatter.format(record)

    @classmethod
    def set_colors(cls, colors: dict):
        """Set replacement or additional colors.

        Args:
            colors (dict): Key is a logging level and value is a Colorama code to set
                           the color, e.g., `{logging.DEBUG: Fore.CYAN}`
        """
        cls.colors.update(colors)
        cls._formatters = None  # create on next logging call


def module_logger(name: str, color: bool = True) -> logging.Logger:
    """Get the main logger for this module.

    First time it's called, sets up the logger.

    Returns:
        logging.Logger: An initialized logger with a handler writing to stderr
    """
    log = logging.getLogger(name)
    if not log.handlers:
        h = logging.StreamHandler()
        if color:
            h.setFormatter(ColorFormatter())
        else:
            fmt = "[{levelname}] {asctime} ({name}) {message}"
            h.setFormatter(logging.Formatter(fmt, style="{"))
        log.addHandler(h)
    return log


def add_log_options(parser: argparse.ArgumentParser) -> None:
    """Add logging-specific options to the argument parser

    Args:
        parser (argparse.ArgumentParser): Parser to modify
    """
    parser.add_argument(
        "--no-color", action="store_true", help="Do not add colors to the logs"
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Minimal logging")
    parser.add_argument(
        "-v",
        action="count",
        dest="vb",
        default=0,
        help="Increase verbosity (repeatable)",
    )


def process_log_options(module_name: str, args: argparse.Namespace) -> logging.Logger:
    log = module_logger(module_name, color=not args.no_color)
    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel(
            (logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG)[
                min(args.vb, 3)
            ]
        )
    return log
