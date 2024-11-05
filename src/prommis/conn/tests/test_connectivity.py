"""
Tests for `connectivity` module.
"""
import re
from idaes_ui.conn import connectivity
import pytest

def test_outputformats():
    ofmt = connectivity.OutputFormats
    assert ofmt.get_ext(ofmt.markdown) == "md"
    assert ofmt.get_ext(ofmt.html) == "html"

@pytest.fixture
def example_conn():
    #
    # Connectivity for this little loop:
    #
    # UnitA -- Stream1 --> UnitB ---+
    #  ^                            |
    #  +----- Stream 2 -------------+
    #
    u = {"Unit A": "U-A", "Unit B": "U-B"}
    s = {"Stream 1": "S-1", "Stream2": "S-2"}
    c = {"Stream1": ["UnitA", "UnitB"], "Stream2": ["UnitB", "UnitA"]}
    conn = connectivity.Connectivity(units=u, streams=s, connections=c)
    yield conn

def test_mermaid(example_conn):
    mmd = connectivity.Mermaid(example_conn)
    s = mmd.write(None, output_format=connectivity.OutputFormats.mermaid.value)
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
