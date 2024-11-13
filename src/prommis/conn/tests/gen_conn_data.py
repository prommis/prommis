"""
Run this script to regenerate the file
`connectivity_data.py`, which will be imported
in `test_connectivity.py` to perform regression tests.
"""

from io import StringIO
from prommis.uky.uky_flowsheet import build
from prommis.conn.connectivity import (
    ModelConnectivity,
    ConnectivityFromFile,
    Mermaid,
    OutputFormats,
)

TAB = "    "


def csv_data(fs, filename):
    mc = ModelConnectivity(fs)
    csv = StringIO()
    mc.write(csv)
    lines = csv.getvalue().split("\n")
    fp = open(filename, "w")
    write_header(fp)
    write_list(fp, lines, "csv")
    fp.write("\n\n")
    csv.seek(0)
    return csv, fp


def mermaid_data(csv, fp):
    conn = ConnectivityFromFile(csv)
    mmd = Mermaid(conn.connectivity)
    mermaid = StringIO()
    mmd.write(mermaid, output_format=OutputFormats.MERMAID)
    lines = mermaid.getvalue().split("\n")
    write_list(fp, lines, "mermaid")


def write_header(fp):
    fp.write("import pytest\n\n")


def write_list(fp, lines, name):
    while lines[-1] == "":
        lines = lines[:-1]
    fp.write("@pytest.fixture\n")
    fp.write(f"def uky_{name}():\n")
    fp.write(f"{TAB}{TAB}return [\n")
    fp.write(",\n".join((f'{TAB}{TAB}{TAB}"{x}"' for x in lines)))
    fp.write(f"\n{TAB}{TAB}]\n")


if __name__ == "__main__":
    filename = "connectivity_data.py"
    fs = build()
    csv, fp = csv_data(fs, filename)
    mermaid_data(csv, fp)
    print(f"Updated file: {filename}")
