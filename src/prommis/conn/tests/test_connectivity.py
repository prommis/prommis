"""
Tests for `connectivity` module.
"""

from io import StringIO
import re
import pytest

from prommis.conn import connectivity
from prommis.conn.connectivity import OutputFormats as OF
from prommis.uky.uky_flowsheet import build
from prommis.conn.tests import connectivity_data as cdata


def test_outputformats():
    assert OF.MARKDOWN == OF("markdown")
    assert OF.HTML == OF("html")
    assert OF.MERMAID == OF("mermaid")
    assert OF.CSV == OF("csv")


def test__output_format():
    assert connectivity._output_format(OF.CSV) == OF.CSV
    with pytest.raises(ValueError):
        connectivity._output_format("--bad--")


@pytest.fixture
def example_conn():
    #
    # Connectivity for this little loop:
    #
    # UnitA -- Stream1 --> UnitB ---+
    #  ^                            |
    #  +----- Stream 2 -----<-------+
    #
    u = {"Unit A": "U-A", "Unit B": "U-B"}
    s = {"Stream 1": "S-1", "Stream 2": "S-2"}
    c = {"S-1": ["U-A", "U-B"], "S-2": ["U-B", "U-A"]}
    conn = connectivity.Connectivity(units=u, streams=s, connections=c)
    yield conn


@pytest.mark.unit
def test_mermaid(example_conn):
    mmd = connectivity.Mermaid(example_conn)
    s = mmd.write(None, output_format=connectivity.OutputFormats.MERMAID.value)
    unit_patterns = [r"U-A..?Unit A..?", r"U-B..?Unit B..?"]
    connection_patterns = [r"U-A\s*-->\s*U-B", r"U-B\s*-->\s*U-A"]

    def find_and_remove(text, patterns):
        match_idx, match_item = -1, None
        for j, pat in enumerate(patterns):
            m = re.search(pat, line)
            if m:
                match_idx, match_item = j, pat
                break
        if match_idx >= 0:
            patterns.remove(match_item)
        return match_idx

    for i, line in enumerate(s.split("\n")):
        line = line.strip()
        if not line:
            continue
        if i == 0:
            assert line.startswith("flowchart")
        else:
            patterns = unit_patterns if i < 3 else connection_patterns
            match_idx = find_and_remove(line, patterns)
            assert match_idx >= 0
    # everything was found
    assert len(unit_patterns) == 0
    assert len(connection_patterns) == 0


@pytest.mark.unit
def test_mermaid_containers(example_conn):
    ofmt = connectivity.OutputFormats
    mmd = connectivity.Mermaid(example_conn)
    base = mmd.write(None, output_format=ofmt.MERMAID)
    for container in ofmt.HTML, ofmt.MARKDOWN:
        text = mmd.write(None, output_format=container)
        for line in text.split("\n"):
            assert line.strip() in text


@pytest.mark.unit
def test_mermaid_options(example_conn):
    kwargs_list = [{}, {"direction": "TD"}, {"direction": "TD", "stream_labels": True}]
    for kwargs in kwargs_list:
        mmd = connectivity.Mermaid(example_conn, **kwargs)
        _ = mmd.write(None, output_format=connectivity.OutputFormats.MERMAID)


@pytest.mark.unit
def test_ordering():
    # is the ordering consistent?
    model = build()
    conn1 = connectivity.create_from_model(model=model)
    conn2 = connectivity.create_from_model(model=model)

    for attribute in "units", "streams", "connections":
        keys1 = getattr(conn1, attribute).keys()
        keys2 = getattr(conn2, attribute).keys()
        assert list(keys1) == list(keys2)


# this roundabout method avoids pylint warnings
uky_csv_data, uky_mermaid_data = cdata.uky_csv, cdata.uky_mermaid


@pytest.mark.unit
def test_uky_data(uky_csv_data, uky_mermaid_data):
    model = build()
    for fmt, ref_lines in (("csv", uky_csv_data), ("mermaid", uky_mermaid_data)):
        print(f"Format={fmt}")
        sio = StringIO()
        csv = connectivity.create_from_model(model=model, to_format=fmt, ofile=sio)
        data = sio.getvalue()
        lines = data.split("\n")
        for line, ref_line in zip(lines, ref_lines):
            assert line == ref_line


# With a parametrized fixture, create a set of CSV files
# representing different connectivities.
# First element is expected exception (None if OK)
matrix_list = [
    (None, [[1, 0], [-1, 1], [0, -1]]),  # (S1) -> U1 -(S2)-> U2 -> (S3)
    (None, [[1, 0], [-1, 1], [0, -1, 1]]),  # Bad rowlen / ignored
    (None, [[1, 0], [-1, 1], [0, -1.0]]),  # Float OK
    (ValueError, [[1, 9], [-1, 1], [0, -1]]),  # Bad value
    (ValueError, [[1, "?"], [-1, 1], [0, -1]]),  # Bad value
]


@pytest.fixture(scope="function", params=matrix_list)
def connectivity_info(tmp_path, request):
    exc, m = request.param
    n_streams = len(m)
    n_units = len(m[0])
    csv_file = tmp_path / f"matrix.csv"
    with open(csv_file, "w") as f:
        unit_list = ",".join((f"Unit {n + 1}" for n in range(n_units)))
        f.write(f"Arcs,{unit_list}\n")
        for i in range(n_streams):
            f.write(f"Stream {i + 1},")
            values = ",".join((f"{v}" for v in m[i]))
            f.write(values)
            f.write("\n")
    yield (exc, csv_file)


@pytest.mark.unit
def test_connectivity_builder(connectivity_info):
    expect_exc, csv = connectivity_info
    # csv_dump = "".join(csv.open())
    # print(f"exc={expect_exc}; csv={csv_dump}")
    builder = connectivity.ConnectivityBuilder(input_file=csv)
    if expect_exc is None:
        conn = builder.connectivity
    else:
        with pytest.raises(expect_exc):
            _ = builder.connectivity
    builder = connectivity.ConnectivityBuilder(input_file=csv.open())


@pytest.mark.unit
def test_connectivity_builder_args():
    with pytest.raises(ValueError):
        connectivity.ConnectivityBuilder(input_file=None, input_data=None)


@pytest.mark.unit
@pytest.mark.parametrize("output_format", [OF.HTML, OF.MARKDOWN, OF.MERMAID])
def test_create_from_matrix(connectivity_info, output_format):
    expect_exc, csv = connectivity_info
    if expect_exc is None:
        conn = connectivity.create_from_matrix(ifile=csv)
        connectivity.create_from_matrix(
            ifile=csv, ofile=connectivity.AS_STRING, to_format=output_format
        )
    else:
        with pytest.raises(expect_exc):
            _ = connectivity.create_from_matrix(ifile=csv)


@pytest.mark.unit
@pytest.mark.parametrize(
    "args,code",
    [
        (["--usage"], 0),
        (["prommis.uky.uky_flowsheet", "-O", "{path}/uky_conn.csv", "--to", "csv"], 0),
        (["prommis.uky.uky_flowsheet", "-O", "-", "--to", "csv"], 0),
        (
            [
                "prommis.uky.uky_flowsheet",
                "-tmodule",
                "-O",
                "{path}/uky_conn.csv",
                "--to",
                "csv",
            ],
            0,
        ),
        (["invalidmodule.1.name"], 2),
        (
            ["prommis.me.this"],
            1,
        ),
        (["{path}/uky_conn.csv", "--to", "html"], 0),
        (["{path}/uky_conn.csv"], 0),
        (["{path}/uky_conn.csv", "-v"], 0),
        (["{path}/uky_conn.csv", "-vv"], 0),
        (["{path}/uky_conn.csv", "-q"], 0),
        (["{path}/uky_conn.csv", "-q", "-v"], 0),
        (["{path}/uky_conn.txt"], 0),
        (["{path}/uky_conn.csv", "-tcsv"], 0),
        (["{path}/uky_conn.csv", "--to", "mermaid", "--output-file", "-"], 0),
        (["{path}/uky_conn.csv", "--to", "mermaid", "--labels"], 0),
        (["{path}/nope.csv", "--to", "mermaid", "--labels"], 2),
        (["{path}/nope.csv", "--to", "mermaid", "--type", "csv"], 2),
        (["{path}/nope.mmd", "--to", "mermaid", "--labels"], 2),
        (["--labels"], 2),
        (["{path}/junk.csv", "--to", "html"], 1),
    ],
)
def test_main(tmp_path, args, code, uky_csv_data):
    from_model = args and "uky_flowsheet" in args[0]
    if not from_model:
        csv_file = tmp_path / "uky_conn.csv"
        with csv_file.open("w") as f:
            for line in uky_csv_data:
                f.write(line)
                f.write("\n")
        csv_txt_file = tmp_path / "uky_conn.txt"
        with csv_txt_file.open("w") as f:
            for line in uky_csv_data:
                f.write(line)
                f.write("\n")
        (tmp_path / "junk.csv").open(mode="w").write("This,is,some\njunk\n")
    for i in range(len(args)):
        if "{" in args[i]:
            args[i] = args[i].format(path=tmp_path)
    print(f"Run CLI with args: {args}")
    ret_code = connectivity.main(command_line=args)
    assert ret_code == code
