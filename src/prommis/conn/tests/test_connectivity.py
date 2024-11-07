"""
Tests for `connectivity` module.
"""

from io import StringIO
import re
from prommis.conn import connectivity
from prommis.uky.uky_flowsheet import build
import pytest


def test_outputformats():
    ofmt = connectivity.OutputFormats
    assert ofmt.markdown == ofmt("markdown")
    assert ofmt.html == ofmt("html")
    assert ofmt.mermaid == ofmt("mermaid")
    assert ofmt.csv == ofmt("csv")


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
def test_uky_csv(uky_csv, uky_mermaid):
    model = build()
    for fmt, ref_lines in (("csv", uky_csv), ("mermaid", uky_mermaid)):
        sio = StringIO()
        csv = connectivity.create_from_model(model=model, to_format=fmt, ofile=sio)
        data = sio.getvalue()
        lines = data.split("\n")
        for line, ref_line in zip(lines, ref_lines):
            assert line == ref_line


# Reference data


@pytest.fixture
def uky_csv():
    return [
        "Units,leach_mixer,leach,leach_liquid_feed,sl_sep1,leach_solid_feed,precipitator,sl_sep2,leach_sx_mixer,leach_filter_cake_liquid,leach_filter_cake,precip_sep,precip_purge,precip_sx_mixer,roaster,solex_cleaner_load,solex_cleaner_strip,cleaner_mixer,cleaner_org_make_up,acid_feed3,cleaner_sep,cleaner_purge,solex_rougher_load,load_sep,solex_rougher_scrub,rougher_mixer,rougher_org_make_up,acid_feed1,scrub_sep,solex_rougher_strip,acid_feed2,rougher_sep,sc_circuit_purge",
        "Arcs,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
        "leaching_feed_mixture,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "leaching_liq_feed,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "leaching_liquid_outlet,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "leaching_sol_feed,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "leaching_solid_outlet,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "precip_aq_outlet,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "precip_solid_outlet,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep1_liquid_outlet,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep1_retained_liquid_outlet,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep1_solid_outlet,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep2_aq_purge,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep2_aq_recycle,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep2_liquid_outlet,0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep2_retained_liquid_outlet,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sl_sep2_solid_outlet,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_load_aq_feed,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_load_aq_outlet,0,0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_load_org_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_mixed_org_recycle,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_org_feed,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_strip_acid_feed,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_strip_aq_outlet,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_strip_org_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_strip_org_purge,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0",
        "sx_cleaner_strip_org_recycle,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0",
        "sx_rougher_load_aq_feed,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0",
        "sx_rougher_load_aq_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0",
        "sx_rougher_load_aq_recycle,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0",
        "sx_rougher_load_org_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0",
        "sx_rougher_mixed_org_recycle,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0",
        "sx_rougher_org_feed,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0",
        "sx_rougher_scrub_acid_feed,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0",
        "sx_rougher_scrub_aq_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0",
        "sx_rougher_scrub_aq_recycle,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0",
        "sx_rougher_scrub_org_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0",
        "sx_rougher_strip_acid_feed,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0",
        "sx_rougher_strip_aq_outlet,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0",
        "sx_rougher_strip_org_outlet,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0",
        "sx_rougher_strip_org_purge,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1",
        "sx_rougher_strip_org_recycle,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1,0",
    ]


@pytest.fixture
def uky_mermaid():
    return [
        "flowchart LR",
        "   Unit_B[leach_mixer]",
        "   Unit_C[leach]",
        "   Unit_D[leach_liquid_feed]",
        "   Unit_E[sl_sep1]",
        "   Unit_F[leach_solid_feed]",
        "   Unit_G[precipitator]",
        "   Unit_H[sl_sep2]",
        "   Unit_I[leach_sx_mixer]",
        "   Unit_J[leach_filter_cake_liquid]",
        "   Unit_K[leach_filter_cake]",
        "   Unit_L[precip_sep]",
        "   Unit_M[precip_purge]",
        "   Unit_N[precip_sx_mixer]",
        "   Unit_O[roaster]",
        "   Unit_P[solex_cleaner_load]",
        "   Unit_Q[solex_cleaner_strip]",
        "   Unit_R[cleaner_mixer]",
        "   Unit_S[cleaner_org_make_up]",
        "   Unit_T[acid_feed3]",
        "   Unit_U[cleaner_sep]",
        "   Unit_V[cleaner_purge]",
        "   Unit_W[solex_rougher_load]",
        "   Unit_X[load_sep]",
        "   Unit_Y[solex_rougher_scrub]",
        "   Unit_Z[rougher_mixer]",
        "   Unit_AA[rougher_org_make_up]",
        "   Unit_AB[acid_feed1]",
        "   Unit_AC[scrub_sep]",
        "   Unit_AD[solex_rougher_strip]",
        "   Unit_AE[acid_feed2]",
        "   Unit_AF[rougher_sep]",
        "   Unit_AG[sc_circuit_purge]",
        "   Unit_B --> Unit_C",
        "   Unit_D --> Unit_B",
        "   Unit_C --> Unit_E",
        "   Unit_F --> Unit_C",
        "   Unit_C --> Unit_E",
        "   Unit_G --> Unit_H",
        "   Unit_G --> Unit_H",
        "   Unit_E --> Unit_I",
        "   Unit_E --> Unit_J",
        "   Unit_E --> Unit_K",
        "   Unit_L --> Unit_M",
        "   Unit_L --> Unit_N",
        "   Unit_H --> Unit_L",
        "   Unit_H --> Unit_O",
        "   Unit_H --> Unit_O",
        "   Unit_N --> Unit_P",
        "   Unit_P --> Unit_I",
        "   Unit_P --> Unit_Q",
        "   Unit_R --> Unit_P",
        "   Unit_S --> Unit_R",
        "   Unit_T --> Unit_Q",
        "   Unit_Q --> Unit_G",
        "   Unit_Q --> Unit_U",
        "   Unit_U --> Unit_V",
        "   Unit_U --> Unit_R",
        "   Unit_I --> Unit_W",
        "   Unit_W --> Unit_X",
        "   Unit_X --> Unit_B",
        "   Unit_W --> Unit_Y",
        "   Unit_Z --> Unit_W",
        "   Unit_AA --> Unit_Z",
        "   Unit_AB --> Unit_Y",
        "   Unit_Y --> Unit_AC",
        "   Unit_AC --> Unit_B",
        "   Unit_Y --> Unit_AD",
        "   Unit_AE --> Unit_AD",
        "   Unit_AD --> Unit_N",
        "   Unit_AD --> Unit_AF",
        "   Unit_AF --> Unit_AG",
        "   Unit_AF --> Unit_Z",
        "",
    ]
