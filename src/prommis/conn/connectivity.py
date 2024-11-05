"""
Calculate and display model connectivity.

The main thing this module does is create and process a connectivity
matrix encoded as follows:

    Units,Unit 1 name,Unit 2 Name, ..., Unit N Name
    Arcs,,,,,,,	...,
    Arc1 Name,	-1,	0, 0, ..., 0
    Arc2 Name, 0 , 1, 0 ,... , 0

Where each cell at the intersection of an Arc (row i) and Unit (column j)
is either:
  *  -1 meaning Arc(i) is an outlet of Unit(j), 
  *  1 meaning Arc(i) is an inlet for Unit(j),
  *  0 meaning there is no connection


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
import pprint
import sys
from tempfile import TemporaryFile
from traceback import format_stack
from typing import TextIO
import warnings

try:
    import pyomo
    from pyomo.environ import Block, value
    from pyomo.network import Arc
    from pyomo.network.port import Port
except ImportError as err:
    pyomo = None
    warnings.warn(f"Could not import pyomo: {err}")    

# For logging
from .util.yoder import add_log_options, process_log_options

# Constants
AS_STRING = "-"

# Logging
# This variable is reassigned if run as script
_log = logging.getLogger(__name__)


##############
# Classes
##############

class OutputFormats(enum.Enum):
    markdown = "markdown"
    html = "html"
    mermaid = "mermaid"
    csv = "csv"

    @classmethod
    def get_ext(cls, fmt):
        match fmt:
            case cls.html:
                result = "html"
            case cls.markdown:
                result = "md"
            case cls.mermaid:
                result = "mmd"
            case _:
                raise ValueError("Bad format")
        return result

class Mermaid:
    def __init__(self, connectivity, indent="   "):
        self._conn = connectivity
        self.indent = indent

    def write(self, output_file: str | None, output_format: str = None):
        if output_file is None:
            f = StringIO()
        else:
            f = open(output_file, "w")
        match output_format:
            case OutputFormats.markdown.value:
                f.write("# Graph\n```mermaid\n")
                self._body(f)
                f.write("\n```\n")
            case OutputFormats.mermaid.value:
                self._body(f)
            case OutputFormats.html.value:
                f.write("<!doctype html>\n<html lang='en'>\n<body>\n<pre class='mermaid'>\n")
                self._body(f)
                f.write("</pre>\n<script type='module'>\n")
                f.write("import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs';\n")
                f.write("mermaid.initialize({securityLevel: 'loose', maxEdges: 2000});\n")
                f.write("await mermaid.run();\n")
                f.write("</script></body></html>")
            case _:
                raise ValueError(f"Bad output format: {output_format}")
        if output_file is None:
            return(f.getvalue())

    def _body(self, outfile):
        i = self.indent
        outfile.write("flowchart TD\n")
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
        for stream_abbr, values in self._conn.connections.items():
            if values[0] is not None and values[1] is not None:
                connections.append(f"{values[0]} --> {values[1]}")
            elif values[0] is not None:
                connections.append(f"{values[0]} --> {stream_abbr}")
                show_streams.add(stream_abbr)
            elif values[1] is not None:
                connections.append(f" {stream_abbr} --> {values[1]}")
                show_streams.add(stream_abbr)
        return connections, show_streams


@dataclass
class Connectivity:
    units: dict = field(default_factory=dict)
    streams: dict = field(default_factory=dict)
    connections: dict = field(default_factory=dict)

class ConnectivityFromFile:
    """Build connectivity information from input data.
    """
    def __init__(self, input_file: str | Path | TextIO):
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
            self._c = Connectivity(units=units, streams=streams, connections=connections)
        return self._c

    def _build_units(self):
        units = {}
        c1, c2 = 1, -1
        for s in self._header[1:]:
            abbr = "Unit_"
            if c2 > -1:
                abbr += chr(ord("A") + c2)
            abbr += chr(ord("A") + c1)
            units[s] = abbr
            c1 += 1
            if c1 == 26:
                c1 = 0
                c2 += 1
        return units

    def _build_streams(self):
        streams = {}
        n = 3
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
    """Build connectivity information from a model.
    """
    def __init__(self, model):
        if pyomo is None:
            raise NotImplementedError("Trying to build from a Pyomo model, but Pyomo is not installed")
        self._fs = model.fs
        self._units = []
        self._streams = []
        self._build()


    def _build(self):
        fs = self._fs  # alias
        units_ord, units_idx = {}, 0
        streams_ord, streams_idx = {}, 0
        rows, empty = [], True
        for comp in fs.component_objects(Arc, descend_into=False):
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
                # if empty, there are no columns to lengthen; defer
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
        

    def write(self, f: TextIO):
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

def _get_model(module_name):
    _log.info("[begin] load and build model")
    mod = importlib.import_module(module_name)
    build_function = mod.build
    _log.debug("[begin] build model")
    model = build_function()
    _log.debug("[ end ] build model")
    _log.info("[ end ] load and build model")
    return model
    
def _get_default_filename(fname: str, ext: str) -> str:
    i = fname.rfind(".")
    if i > 0:
        filename = fname[:i] + "." + ext
    else:
        filename = fname + + "." + ext
    return filename

def _real_output_file(file_, to_, input_file=None):
        if file_ is None:
            if input_file is None:
                return None
            ext = OutputFormats.get_ext(to_)
            output_file = _get_default_filename(file_, ext)
            print(f"No output file specified, using '{output_file}'")
        else:
            output_file = file_
        return output_file

##################
# Public interface
##################

def create_from_matrix(ifile: str, ofile: str | None, to_format: str) -> Connectivity | None:
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
        mermaid = Mermaid(conn_file.connectivity)
    except Exception as err:
        err_msg = f"Could not parse connectivity information: {err}. "
        err_msg += f"Stack trace:\n{format_stack()}"
        raise RuntimeError(err_msg)
    
    if args.ofile == AS_STRING:
        print(mermaid.write(None, output_format=args.to))
    else:
        output_file = _real_output_file(args.ofile, args.to, input_file=args.input_file)
        mermaid.write(output_file, output_format=args.to)


    
def create_from_model(module_name: str, ofile: str, to_format: str) -> Connectivity | None:
    """Programmatic interface to create the connectivity or mermaid output from a python model.

    Arguments:
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
    try:
        model = _get_model(module_name)
    except Exception as err:
        raise RuntimeError(f"Could not load model: {err}")

    model_conn = ModelConnectivity(model)
    
    output_file = _real_output_file(ofile, to_format)
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
            else:
                with open(output_file, "w") as ofile:
                    for line in tempfile:
                        ofile.write(line)
        # for Mermaid, etc., process the temporary file
        else:
            conn_file = ConnectivityFromFile(tempfile)
            mermaid = Mermaid(conn_file.connectivity)
            if output_file == AS_STRING:
                print(mermaid.write(None, output_format=to_format))
            else:
                mermaid.write(output_file, output_format=to_format)

##############
# Command-line
##############

def csv_main(args) -> int:
    """CLI function for creating a graph from connectivity matrix specifying units and streams.

    Args:
        args: Parsed args from ArgumentParser

    Returns:
        int: Code for sys.exit()
    """
    _log.info("[begin] create from matrix")
    try:
        create_from_matrix(args.input_file, args.ofile, args.to)
    except Exception as err:
        _log.info("[ end ] create from matrix (1)")
        _log.error("{err}")
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
    try:
        create_from_model(module_name=args.module_name, ofile=args.ofile, to_format=args.to)
    except RuntimeError as err:
        _log.info("[ end ] create from Python model (1)")
        _log.error("{err}")
        return 1
    _log.info("[ end ] create from Python model")

if __name__ == "__main__":
    p_top = argparse.ArgumentParser()
    add_log_options(p_top)
    subp = p_top.add_subparsers()
    p = subp.add_parser("csv", help="Read from a CSV file and generate a graph")
    p.add_argument("input_file", help="Input CSV file")
    p.add_argument(
        "-O",
        "--output-file",
        dest="ofile",
        help=f"Output file",
        default=AS_STRING,
    )
    p.add_argument("--to", help="Output format for mermaid graph (default=mermaid)",
                   choices=("markdown", "mermaid", "html"), default="mermaid")
    p.set_defaults(main_method=csv_main)
    p = subp.add_parser("module", help="Read from a Python model and generate a CSV file or a graph")
    p.add_argument(
        "-O",
        "--output-file",
        dest="ofile",
        help=f"Output file",
        default=AS_STRING,
    )
    p.add_argument("module_name", help="A dotted Python module name")
    p.add_argument("--to", help="Output format for CSV or mermaid graph (default=csv)",
                   choices=("csv", "markdown", "mermaid", "html"), default="csv")
    p.set_defaults(main_method=module_main)
    args = p_top.parse_args()
    _log = process_log_options("idaes_ui.conn.connectivity", args)

    sys.exit(args.main_method(args))
