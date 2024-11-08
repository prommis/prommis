"""
Calculate and display model connectivity.

The main thing this module does is create and process a connectivity
matrix, where each cell at the intersection of an Arc and Unit.
See `USAGE` for details.

This module can be run as a script or used programmatically, using the
public functions `create_from_matrix` and `create_from_model`.
"""

import argparse
import csv
from dataclasses import dataclass, field
import enum
import importlib
from io import StringIO
import logging
from pathlib import Path
import re
import sys
from tempfile import TemporaryFile
from traceback import format_stack
from typing import TextIO
import warnings

# only external dependency -> Pyomo
try:
    import pyomo
    from pyomo.network import Arc
except ImportError as err:
    pyomo = None
    warnings.warn(f"Could not import pyomo: {err}")

# Constants
AS_STRING = "-"

# Logging
# This variable is reassigned if run as script
_log = logging.getLogger(__name__)


##############
# Classes
##############


class OutputFormats(enum.Enum):
    MARKDOWN = "markdown"
    HTML = "html"
    MERMAID = "mermaid"
    CSV = "csv"


class Mermaid:

    # URL to load Mermaid from, for the HTML output
    CDN_URL = "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs"

    def __init__(
        self,
        connectivity,
        stream_labels: bool = False,
        direction: str = "LR",
        indent="   ",
    ):
        self._conn = connectivity
        self.indent = indent
        self._stream_labels = stream_labels
        self._direction = direction

    def write(self, output_file: str | TextIO | None, output_format: str = None):
        if output_file is None:
            f = StringIO()
        elif hasattr(output_file, "write"):
            f = output_file
        else:
            f = open(output_file, "w")

        try:
            fmt = OutputFormats(output_format)
        except ValueError:
            raise ValueError(f"Bad output format: {output_format}")

        if fmt == OutputFormats.MARKDOWN:
            f.write("# Graph\n```mermaid\n")
            self._body(f)
            f.write("\n```\n")
        elif fmt == OutputFormats.MERMAID:
            self._body(f)
        elif fmt == OutputFormats.HTML:
            f.write(
                "<!doctype html>\n<html lang='en'>\n<body>\n<pre class='mermaid'>\n"
            )
            self._body(f)
            f.write("</pre>\n<script type='module'>\n")
            # XXX: May want to make this URL configurable
            url = self.CDN_URL
            f.write(f"import mermaid from '{url}';\n")
            f.write("mermaid.initialize({securityLevel: 'loose', maxEdges: 2000});\n")
            f.write("await mermaid.run();\n")
            f.write("</script></body></html>")
        else:  # !! should not get here
            raise RuntimeError(f"Output format unaccounted for: {output_format}")

        if output_file is None:
            return f.getvalue()

    def _body(self, outfile):
        i = self.indent
        outfile.write(f"flowchart {self._direction}\n")
        # Get connections first, so we know which streams to show
        connections, show_streams = self._get_connections()
        # Units
        for s in self._get_mermaid_units():
            outfile.write(f"{i}{s}\n")
        # Streams
        for abbr, s in self._get_mermaid_streams():
            if abbr in show_streams:
                outfile.write(f"{i}{s}\n")
        # Connections
        for s in connections:
            outfile.write(f"{i}{s}\n")

    def _get_mermaid_units(self):
        for name, abbr in self._conn.units.items():
            yield f"{abbr}[{name}]"

    def _get_mermaid_streams(self):
        for name, abbr in self._conn.streams.items():
            yield abbr, f"{abbr}([{name}])"

    def _get_connections(self):
        connections = []
        show_streams = set()
        # if we are going to label the streams,
        # construct a mapping to get names from abbreviations
        if self._stream_labels:
            abbr_name = {v: k for k, v in self._conn.streams.items()}
        for stream_abbr, values in self._conn.connections.items():
            if values[0] is not None and values[1] is not None:
                if self._stream_labels:
                    label = self._clean_stream_label(abbr_name[stream_abbr])
                    connections.append(f"{values[0]} --|{label}| -->{values[1]}")
                else:
                    connections.append(f"{values[0]} --> {values[1]}")
            elif values[0] is not None:
                connections.append(f"{values[0]} --> {stream_abbr}")
                show_streams.add(stream_abbr)
            elif values[1] is not None:
                connections.append(f" {stream_abbr} --> {values[1]}")
                show_streams.add(stream_abbr)
        return connections, show_streams

    @staticmethod
    def _clean_stream_label(label):
        if label.endswith("_outlet"):
            label = label[:-7]
        elif label.endswith("_feed"):
            label = label[:-5]
        label = label.replace("_", " ")
        return label


@dataclass
class Connectivity:
    units: dict = field(default_factory=dict)
    streams: dict = field(default_factory=dict)
    connections: dict = field(default_factory=dict)


class ConnectivityFromFile:
    """Build connectivity information from input data."""

    def __init__(self, input_file: str | Path | TextIO):
        """Constructor

        Args:
            input_file (str | Path | TextIO): Input CSV file
        """
        if isinstance(input_file, str) or isinstance(input_file, Path):
            datafile = open(input_file, "r")
        else:
            datafile = input_file
        reader = csv.reader(datafile)
        self._header = next(reader)
        self._rows = list(reader)
        self._c = None

    @property
    def connectivity(self):
        if self._c is None:
            units = self._build_units()
            streams = self._build_streams()
            connections = self._build_connections(units, streams)
            self._c = Connectivity(
                units=units, streams=streams, connections=connections
            )
        return self._c

    def _build_units(self):
        units = {}
        c1, c2 = 1, -1
        for s in self._header[1:]:
            abbr = "Unit_"
            # Pick abbreviations that match spreadsheet column names,
            # for easier comparison and debugging.
            # i.e. A-Z then AA-ZZ. chr(x) returns ASCII character at x,
            # and ord(x) is the reverse function. e.g., chr(ord("A") + 1) == "B"
            if c2 > -1:
                abbr += chr(ord("A") + c2)
            abbr += chr(ord("A") + c1)
            units[s] = abbr
            c1 += 1
            if c1 == 26:  # after Z, start on AA..ZZ
                c1 = 0
                c2 += 1
        return units

    def _build_streams(self):
        streams = {}
        n = 3  # pick numbers that match spreadsheet row numbers
        for row in self._rows[1:]:
            s = row[0]
            abbr = f"Stream_{n}"
            streams[s] = abbr
            n += 1
        return streams

    def _build_connections(self, units, streams):
        connections = {s: [None, None] for s in streams.values()}
        for row in self._rows[1:]:
            stream_name = row[0]
            col = 1
            for conn in row[1:]:
                if conn not in ("", "0"):
                    conn = max(0, int(conn))  # -1 -> 0, 1 -> 1
                    try:
                        unit_name = self._header[col]
                    except IndexError:
                        print(f"col={col} :: header-len={len(self._header)}")
                        raise
                    unit_abbr = units[unit_name]
                    stream_abbr = streams[stream_name]
                    connections[stream_abbr][conn] = unit_abbr
                col += 1
        return connections


class ModelConnectivity:
    """Build connectivity information from a model."""

    def __init__(self, model, flowsheet_attr: str = "fs"):
        """Constructor

        Args:
            model: Pyomo ConcreteModel instance with an attribute, ".fs"
                   that is (or acts like) an IDAES flowsheet

        Raises:
            NotImplementedError: If Pyomo isn't installed
        """
        if pyomo is None:
            raise NotImplementedError(
                "Trying to build from a Pyomo model, but Pyomo is not installed"
            )
        fa = flowsheet_attr.strip("'").strip('"').strip()
        self._fs = eval(f"model.{fa}")
        self._units = []
        self._streams = []
        self._build()

    def _build(self):
        fs = self._fs  # alias
        units_ord, units_idx = {}, 0
        streams_ord, streams_idx = {}, 0
        rows, empty = [], True
        for comp in self._arcs_sorted_by_name(fs):
            stream_name = comp.getname()
            src, dst = comp.source.parent_block(), comp.dest.parent_block()
            src_name, dst_name = src.getname(), dst.getname()

            src_i, dst_i, stream_i = -1, -1, -1

            try:
                idx = streams_ord[stream_name]
            except KeyError:
                self._streams.append(stream_name)
                idx = streams_ord[stream_name] = streams_idx
                streams_idx += 1
                if empty:
                    rows = [[]]
                    empty = False
                else:
                    rows.append([0] * len(rows[0]))
            stream_i = idx

            for unit_name, is_src in (src_name, True), (dst_name, False):
                try:
                    idx = units_ord[unit_name]
                except KeyError:
                    self._units.append(unit_name)
                    idx = units_ord[unit_name] = units_idx
                    units_idx += 1
                    for row in rows:
                        row.append(0)
                if is_src:
                    src_i = idx
                else:
                    dst_i = idx

            rows[stream_i][src_i] = -1
            rows[stream_i][dst_i] = 1

        self._rows = rows

    @staticmethod
    def _arcs_sorted_by_name(fs):
        """Try and make the output stable by looping through the Pyomo
        Arcs in alphabetical order, by their name.
        """
        arcs = fs.component_objects(Arc, descend_into=False)
        return sorted(arcs, key=lambda arc: arc.getname())

    def write(self, f: TextIO):
        """Write the CSV file."""
        header = self._units.copy()
        header.insert(0, "Units")
        f.write(",".join(header))
        f.write("\n")
        header_sep = ["" for _ in self._units]
        header_sep.insert(0, "Arcs")
        f.write(",".join(header_sep))
        f.write("\n")
        for row_idx, row in enumerate(self._rows):
            row.insert(0, self._streams[row_idx])
            f.write(",".join((str(value) for value in row)))
            f.write("\n")


############
# Utility
############


def _get_model(module_name, build_func):
    _log.info("[begin] load and build model")
    mod = importlib.import_module(module_name)
    build_function = getattr(mod, build_func)
    _log.debug(f"[begin] build model function={build_func}")
    model = build_function()
    _log.debug("[ end ] build model")
    _log.info("[ end ] load and build model")
    return model


def _real_output_file(ifile: str, to_, input_file=None):
    try:
        to_fmt = OutputFormats(to_)
    except ValueError:
        raise ValueError(f"Bad format: {to_}")
    ext = {
        OutputFormats.HTML: "html",
        OutputFormats.MARKDOWN: "md",
        OutputFormats.MERMAID: "mmd",
        OutputFormats.CSV: "csv",
    }[to_fmt]
    i = ifile.rfind(".")
    if i > 0:
        filename = ifile[:i] + "." + ext
    else:
        filename = ifile + +"." + ext
    _log.info(f"No output file specified, using '{filename}'")
    return filename


##################
# Public interface
##################


def create_from_matrix(
    ifile: str, ofile: str | None, to_format: str, mermaid_options: dict = None
) -> Connectivity | None:
    """Programmatic interface to create a graph of the model from a connectivity matrix.

    Args:
        ifile (str): Input filename
        ofile: Output file name. If this is the special value defined by `AS_STRING`, then
               The output will go to the console. If None, then no output will
               be created and the connectivity will be returned as an object.
        to_format: Output format, which should match one of the values in `OutputFormat`

    Raises:
        RuntimeError: For all errors captured during Mermaid processing
    """
    conn_file = ConnectivityFromFile(ifile)

    if ofile is None:
        return conn_file.connectivity

    try:
        mermaid_kw = mermaid_options or {}
        mermaid = Mermaid(conn_file.connectivity, **mermaid_kw)
    except Exception as err:
        err_msg = f"Could not parse connectivity information: {err}. "
        err_msg += f"Stack trace:\n{format_stack()}"
        raise RuntimeError(err_msg)

    if ofile == AS_STRING:
        print(mermaid.write(None, output_format=to_format))
    else:
        mermaid.write(ofile, output_format=to_format)


def create_from_model(
    model: object = None,
    module_name: str = None,
    ofile: str | TextIO = None,
    to_format: str = None,
    build_func: str = "build",
    flowsheet_attr: str = "fs",
    mermaid_options: dict = None,
) -> Connectivity | None:
    """Programmatic interface to create the connectivity or mermaid output from a python model.

    Arguments:
        model: If present, the model touse
        module_name: Dotted Python module name (absolute, e.g. package.subpackage.module).
                     The protocol is to call the `build()` function in the module to get
                     back a model.
        ofile: Output file name. If this is the special value defined by `AS_STRING`, then
               The output will go to the console. If None, then no output will
               be created and the connectivity will be returned as an object.
        to_format: Output format, which should match one of the values in `OutputFormat`

    Returns:
        Connectivity instance, if ofile is None

    Raises:
        RuntimeError: For all errors captured while building the model
    """
    if model is None:
        try:
            model = _get_model(module_name, build_func)
        except Exception as err:
            raise RuntimeError(f"Could not load model: {err}")
    else:
        pass  # assume it's already loaded into this variable

    model_conn = ModelConnectivity(model, flowsheet_attr=flowsheet_attr)

    output_file = ofile
    conn_file = None
    with TemporaryFile(mode="w+t") as tempfile:
        model_conn.write(tempfile)
        tempfile.flush()  # make sure all data is written
        tempfile.seek(0)  # reset to start of file for reading
        # if no output, just parse the file and return the result
        if ofile is None:
            conn_file = ConnectivityFromFile(tempfile)
            return conn_file.connectivity
        # for CSV, just copy the temporary file
        elif to_format == "csv":
            tempfile.seek(0)
            if ofile == AS_STRING:
                for line in tempfile:
                    print(line, end="")
            elif hasattr(ofile, "write"):
                for line in tempfile:
                    ofile.write(line)
            else:
                with open(output_file, "w") as ofile:
                    for line in tempfile:
                        ofile.write(line)
        # for Mermaid, etc., process the temporary file
        else:
            conn_file = ConnectivityFromFile(tempfile)
            mermaid_kw = mermaid_options or {}
            mermaid = Mermaid(conn_file.connectivity, **mermaid_kw)
            if output_file == AS_STRING:
                print(mermaid.write(None, output_format=to_format))
            else:
                mermaid.write(output_file, output_format=to_format)


###############
# Command-line
###############


def csv_main(args) -> int:
    """CLI function for creating a graph from connectivity matrix specifying units and streams.

    Args:
        args: Parsed args from ArgumentParser

    Returns:
        int: Code for sys.exit()
    """
    _log.info("[begin] create from matrix")

    if args.ofile is None:
        args.ofile = _real_output_file(args.source, args.to)
        print(f"Output in: {args.ofile}")

    mermaid_options = {"stream_labels": args.labels, "direction": args.direction}

    try:
        create_from_matrix(
            args.source, args.ofile, args.to, mermaid_options=mermaid_options
        )
    except Exception as err:
        _log.info("[ end ] create from matrix (1)")
        _log.error(f"{err}")
        return 1
    _log.info("[ end ] create from matrix")

    return 0


def module_main(args) -> int:
    """CLI function for creating connectivity/graph from a Python model.

    Args:
        args: Parsed args from ArgumentParser

    Returns:
        int: Code for sys.exit()
    """
    _log.info("[begin] create from Python model")

    if args.ofile is None:
        args.ofile = AS_STRING

    mermaid_options = {"stream_labels": args.labels, "direction": args.direction}

    try:
        create_from_model(
            module_name=args.source,
            ofile=args.ofile,
            to_format=args.to,
            flowsheet_attr=args.fs,
            build_func=args.build,
            mermaid_options=mermaid_options,
        )
    except RuntimeError as err:
        _log.info("[ end ] create from Python model (1)")
        _log.error(f"{err}")
        return 1
    _log.info("[ end ] create from Python model")


SCRIPT_NAME = "connectivity"
USAGE = f"""
This script generates connectivity information from models,
or graphs of the model structure from connectivity information,
or both together -- i.e. graphs of the model structure from the model.

It has two main 'modes': csv and module.

The *csv* mode starts from a connectivity matrix stored in a file, and 
generates MermaidJS code (optionally wrapped by Markdown or HTML) that can be
rendered as a simple diagram of the flowsheet.

The *module* mode starts from a Python module, calls its 'build()' function to
build a model, then either (a) writes out the connectivity CSV for that model, or
(b) generates MermaidJS code as in the csv mode above. The second function
is a convenience and is equivalent to generating the CSV file then running again
in 'csv' mode with that file as input, which will be shown below.

You can explicitly indicate the mode with the --type/-t argument, though the
program will try to infer it as well (anything ending in ".csv" will be assumed to
be a CSV file, for example).

Example command-lines (showing the two modes):

    # Generate the connectivity matrix in uky_conn.csv
    {SCRIPT_NAME} prommis.uky.uky_flowsheet -O uky_conn.csv --to csv

    # Generate the MermaidJS code wrapped in a HTML page that can be viewed in a
    # browser without any further installation (MermaidJS is fetched from the network)
    # The page will be called 'uky_conn.html' (since no filename was specified).
    {SCRIPT_NAME} uky_conn.csv --to html

    # Print the 'raw' MermaidJS code to the console instead of to a file
    {SCRIPT_NAME}  uky_conn.csv --to mermaid --output-file "-"

    # Print mermaid info to default file, with streams labeled
    {SCRIPT_NAME} uky_conn.csv --to mermaid --labels
    # (console)> Output in: uky_conn.mmd

For more information about MermaidJS, see http://mermaid.js.org

The connectivity matrix format is:

    Units,Unit 1 name,Unit 2 Name, ..., Unit N Name
    Arcs,,,	... ,
    Arc1 Name,	-1,	0, 0, ..., 0
    Arc2 Name, 0 , 1, 0 ,... , 0

Where each cell at the intersection of an Arc (row i) and Unit (column j)
is either:
  *  -1 meaning Arc(i) is an outlet of Unit(j), 
  *  1 meaning Arc(i) is an inlet for Unit(j),
  *  0 meaning there is no connection

"""


def _add_log_options(parser: argparse.ArgumentParser) -> None:
    """Add logging-specific options to the argument parser

    Args:
        parser (argparse.ArgumentParser): Parser to modify
    """
    parser.add_argument("-q", "--quiet", action="store_true", help="Minimal logging")
    parser.add_argument(
        "-v",
        action="count",
        dest="vb",
        default=0,
        help="Increase verbosity (repeatable)",
    )


def _process_log_options(module_name: str, args: argparse.Namespace) -> logging.Logger:
    log = logging.getLogger(module_name)
    if not log.handlers:
        h = logging.StreamHandler()
        fmt = "[{levelname}] {asctime} ({name}) {message}"
        h.setFormatter(logging.Formatter(fmt, style="{"))
        log.addHandler(h)
    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel(
            (logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG)[
                min(args.vb, 3)
            ]
        )
    return log


def main():
    p = argparse.ArgumentParser(description="Process and/or generate model connectivity information")
    p.add_argument("--usage", action="store_true", help="Print usage with examples")
    # set nargs=? so --usage works without any other argument; though
    # this will require more checks later
    p.add_argument("source", help="Build source", metavar="FILE or MODULE", nargs="?")
    p.add_argument(
        "--type",
        "-t",
        choices=("csv", "module"),
        help="Build source type: csv=CSV file, module=Python module",
        default=None,
    )
    p.add_argument(
        "-O",
        "--output-file",
        dest="ofile",
        help=f"Output file",
        default=None,
    )
    p.add_argument(
        "--to",
        help="Output format for mermaid graph (default=csv)",
        choices=("markdown", "mermaid", "html", "csv"),
        default="csv",
    )
    p.add_argument(
        "--fs",
        help="Name of flowsheet attribute on model object (default=fs)",
        default="fs",
    )
    p.add_argument(
        "--build",
        help="Name of build function in module (dafault=build)",
        default="build",
    )
    p.add_argument(
        "--labels", "-L", help="Add stream labels to diagram", action="store_true"
    )
    p.add_argument(
        "--direction",
        "-D",
        help="Direction of diagram",
        choices=("LR", "TD"),
        default="LR",
    )
    _add_log_options(p)
    args = p.parse_args()
    if args.usage:
        print(USAGE)
        return 0
    if args.source is None:
        print("File or module source is required. Try --usage for details.\n")
        p.print_help()
        return 0
    _log = _process_log_options("idaes_ui.conn.connectivity", args)
    if args.type is None:
        main_method = None
        if args.source.lower().endswith(".csv"):
            path = Path(args.source)
            if not path.exists():
                p.error(
                    f"Source looks like a CSV file, but does not exist: {args.source}"
                )
            main_method = csv_main
        elif "/" in args.source:
            path = Path(args.source)
            if path.exists():
                _log.warning(
                    "File path given, but suffix is not .csv; assuming CSV mode"
                )
                main_method = csv_main
            else:
                p.error(
                    f"Source looks like file path, but does not exist: {args.source}"
                )
        else:
            m = re.match(
                r"[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)*", args.source
            )
            if m.span() != (0, len(args.source)):
                p.error(
                    "Source looks like a module name, but is not valid: {srgs.source}"
                )
            main_method = module_main
    else:
        if args.type == "csv":
            if not Path(args.source).exists():
                p.error(f"Source file path does not exist: {args.source}")
            main_method = csv_main
        elif args.type == "module":
            main_method = module_main

    return main_method(args)


if __name__ == "__main__":
    sys.exit(main())
