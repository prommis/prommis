"""
Visualize UKy flowsheet with IDAES Flowsheet Visualizer

Usage:
    python uky_fv.py
"""

# system
import argparse
import logging
import sys

# idaes
try:
    from idaes_ui.fv import visualize
except ImportError:
    visualize = None
from idaes import logger as idaeslog

# prommis
from prommis.uky.uky_flowsheet_ui import build_flowsheet
from prommis.examples import util

_log = idaeslog.getLogger(__name__)


def main():
    if visualize is None:
        print("cannot continue: idaes_ui not installed")
        return -1

    result = 0

    p = argparse.ArgumentParser()
    p.add_argument("-q", "--quiet", action="store_true")
    args = p.parse_args()

    if args.quiet:
        util.silence_log(_log)
        util.silence_log(idaeslog.getLogger("idaes_ui.fv.fsvis"))
    else:
        _log.setLevel(logging.WARNING)

    fs = build_flowsheet().fs
    try:
        visualize(name="UKy", flowsheet=fs, loop_forever=True)
    except KeyboardInterrupt:
        pass
    except Exception as err:
        _log.critical(f"in visualize(): {err}")
        result = 1

    return result


if __name__ == "__main__":
    sys.exit(main())
