"""
Run this script to regenerate the file
`connectivity_data.py`, which will be imported
in `test_connectivity.py` to perform regression tests.
"""

from io import StringIO
from prommis.uky.uky_flowsheet import build
from prommis.conn.connectivity import (
    ModelConnectivity,
    ConnectivityBuilder,
    Mermaid,
    OutputFormats,
)

TAB = "    "


def csv_data(fs, filename):
    mc = ModelConnectivity(fs)
    data = mc.get_data()
    fp = open(filename, "w")
    write_header(fp)
    lines = (",".join(map(str, values)) for values in data)
    write_list(fp, lines, "csv")
    fp.write("\n\n")
    return data, fp


def mermaid_data(csv, fp):
    conn = ConnectivityBuilder(input_data=csv)
    mmd = Mermaid(conn.connectivity)
    mermaid = StringIO()
    mmd.write(mermaid, output_format=OutputFormats.MERMAID)
    lines = mermaid.getvalue().split("\n")
    write_list(fp, lines, "mermaid")


def write_header(fp):
    fp.write("import pytest\n\n")


def write_list(fp, lines, name):
    fp.write("@pytest.fixture\n")
    fp.write(f"def uky_{name}():\n")
    fp.write(f"{TAB}return [\n")
    fp.write(",\n".join((f'{TAB}{TAB}"{x}"' for x in lines)))
    fp.write(f"\n{TAB}]\n")


if __name__ == "__main__":
    filename = "connectivity_data.py"
    fs = build()
    csv, fp = csv_data(fs, filename)
    mermaid_data(csv, fp)
    print(f"Updated file: {filename}")
