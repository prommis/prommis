"""
Tests for `connectivity` module.
"""

from io import StringIO
import re
import pytest

from prommis.conn import connectivity
from prommis.uky.uky_flowsheet import build
from prommis.conn.tests.connectivity_data import uky_csv, uky_mermaid


def test_outputformats():
    ofmt = connectivity.OutputFormats
    assert ofmt.MARKDOWN == ofmt("markdown")
    assert ofmt.HTML == ofmt("html")
    assert ofmt.MERMAID == ofmt("mermaid")
    assert ofmt.CSV == ofmt("csv")


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
    s = {"Stream 1": "S-1", "Stream2": "S-2"}
    c = {"Stream1": ["UnitA", "UnitB"], "Stream2": ["UnitB", "UnitA"]}
    conn = connectivity.Connectivity(units=u, streams=s, connections=c)
    yield conn


@pytest.mark.unit
def test_mermaid(example_conn):
    mmd = connectivity.Mermaid(example_conn)
    s = mmd.write(None, output_format=connectivity.OutputFormats.MERMAID.value)
    unit_patterns = [r"U-A..?Unit A..?", r"U-B..?Unit B..?"]
    connection_patterns = [r"UnitA\s*-->\s*UnitB", r"UnitB\s*-->\s*UnitA"]

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
def test_ordering():
    # is the ordering consistent?
    model = build()
    conn1 = connectivity.create_from_model(model=model)
    conn2 = connectivity.create_from_model(model=model)

    for attribute in "units", "streams", "connections":
        keys1 = getattr(conn1, attribute).keys()
        keys2 = getattr(conn2, attribute).keys()
        assert list(keys1) == list(keys2)


@pytest.mark.unit
def test_uky_data(uky_csv, uky_mermaid):
    model = build()
    for fmt, ref_lines in (("csv", uky_csv), ("mermaid", uky_mermaid)):
        sio = StringIO()
        csv = connectivity.create_from_model(model=model, to_format=fmt, ofile=sio)
        data = sio.getvalue()
        lines = data.split("\n")
        for line, ref_line in zip(lines, ref_lines):
            assert line == ref_line
