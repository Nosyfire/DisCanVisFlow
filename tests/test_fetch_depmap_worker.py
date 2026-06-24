"""Unit tests for fetch_depmap_worker.py pure helpers (no network)."""
import importlib.util
import io
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "fetch_depmap_worker.py"
spec = importlib.util.spec_from_file_location("fetch_depmap_worker", BIN)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_pos_from_protein_change():
    assert m.pos_from_protein_change("p.G12D") == "12"
    assert m.pos_from_protein_change("p.(Gly12Asp)") == "12"
    assert m.pos_from_protein_change("p.Gly12Aspfs*5") == "12"
    assert m.pos_from_protein_change("p.*42L") == "42"
    assert m.pos_from_protein_change("") == ""
    assert m.pos_from_protein_change("silent") == ""


def test_normalise_maps_columns_and_skips_unparseable():
    raw = io.StringIO(
        "HugoSymbol,ProteinChange,ModelID,Pos,EntrezGeneID,Hotspot\n"
        "KRAS,p.G12D,ACH-000001,25245350,3845,Y\n"
        "TP53,p.R175H,ACH-000002,7675088,7157,\n"
        "NOISE,,ACH-000003,1,0,\n"          # no protein change → skipped
    )
    rows = m.normalise(raw)
    assert len(rows) == 2
    kras = rows[0]
    assert kras["HugoSymbol"] == "KRAS"
    assert kras["HGVSp_Short"] == "p.G12D"
    assert kras["Protein_position"] == "12"
    assert kras["ModelID"] == "ACH-000001"
    assert kras["Start_Position"] == "25245350"
    assert rows[1]["HugoSymbol"] == "TP53"
    assert rows[1]["Protein_position"] == "175"


def test_normalise_accepts_alt_column_names():
    raw = io.StringIO(
        "Hugo_Symbol,HGVSp_Short,DepMap_ID\n"
        "BRAF,p.V600E,ACH-000009\n"
    )
    rows = m.normalise(raw)
    assert len(rows) == 1
    assert rows[0]["HugoSymbol"] == "BRAF"
    assert rows[0]["Protein_position"] == "600"
    assert rows[0]["ModelID"] == "ACH-000009"


def test_resolve_file_url_picks_newest_release(monkeypatch):
    catalogue = (
        "release,release_date,filename,url,md5_hash\n"
        "DepMap Public 25Q4,2025-10-01,OmicsSomaticMutations.csv,http://x/old,aaa\n"
        "DepMap Public 26Q1,2026-04-01,OmicsSomaticMutations.csv,http://x/new,bbb\n"
        "DepMap Public 26Q1,2026-04-01,ScreenSequenceMap.csv,http://x/other,ccc\n"
    )

    class FakeResp:
        text = catalogue

    monkeypatch.setattr(m, "_req", lambda *a, **k: FakeResp())
    url = m.resolve_file_url(None, "ignored", "OmicsSomaticMutations.csv", "")
    assert url == "http://x/new"
    # explicit release wins
    url2 = m.resolve_file_url(None, "ignored", "OmicsSomaticMutations.csv",
                              "DepMap Public 25Q4")
    assert url2 == "http://x/old"
