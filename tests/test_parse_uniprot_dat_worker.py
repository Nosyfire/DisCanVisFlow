"""
Tests for bin/parse_uniprot_dat_worker.py — bulk UniProt feature parser.
"""
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "parse_uniprot_dat_worker.py"


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_dat(path: Path, entries: list[str]) -> Path:
    """Write a minimal .dat.gz with given entry blocks."""
    content = "\n".join(entries) + "\n"
    gz = path / "test.dat.gz"
    with gzip.open(gz, 'wt', encoding='latin-1') as f:
        f.write(content)
    return gz


def _entry(acc: str, features: str = "", dr_lines: str = "") -> str:
    """Minimal UniProt flat-file entry."""
    return (
        f"ID   TEST_{acc}              Reviewed;        100 AA.\n"
        f"AC   {acc};\n"
        f"OX   NCBI_TaxID=9606;\n"
        + features
        + dr_lines
        + "//\n"
    )


def _signal_feat(start=1, end=25) -> str:
    return (
        f"FT   SIGNAL          {start}..{end}\n"
        f"FT                   /evidence=\"ECO:0000269\"\n"
    )


def _transmem_feat(start=30, end=55, note="TM helix") -> str:
    return (
        f"FT   TRANSMEM        {start}..{end}\n"
        f"FT                   /note=\"{note}\"\n"
    )


def _binding_feat(start=70, end=70, ligand="ATP") -> str:
    return (
        f"FT   BINDING         {start}..{end}\n"
        f"FT                   /ligand=\"{ligand}\"\n"
    )


def _pfam_dr(pfam_acc="PF00001", pfam_name="7tm_1") -> str:
    return f"DR   Pfam; {pfam_acc}; {pfam_name}; 1.\n"


def _make_interpro(path: Path, rows: list[tuple]) -> Path:
    """Write minimal protein2ipr.dat.gz in the actual EBI 6-column format.
    rows = [(acc, pfam_acc, pfam_name, start, end), ...]
    Columns: UniProt_acc  InterPro_acc  InterPro_desc  sig_acc  start  end
    """
    gz = path / "protein2ipr.dat.gz"
    with gzip.open(gz, 'wt', encoding='utf-8') as f:
        for acc, pfam_acc, pfam_name, start, end in rows:
            f.write(f"{acc}\tIPR000001\t{pfam_name}\t{pfam_acc}\t{start}\t{end}\n")
    return gz


def _run(tmp_path: Path, dat_gz: Path, interpro_gz: Path | None = None,
         accessions: list[str] | None = None) -> subprocess.CompletedProcess:
    acc_file = None
    if accessions:
        acc_file = tmp_path / "accs.txt"
        acc_file.write_text("\n".join(accessions))

    cmd = [sys.executable, str(WORKER),
           "--uniprot_dat", str(dat_gz),
           "--outdir",      str(tmp_path / "out")]
    if interpro_gz:
        cmd += ["--interpro_pfam", str(interpro_gz)]
    if acc_file:
        cmd += ["--accessions", str(acc_file)]

    return subprocess.run(cmd, capture_output=True, text=True)


# ── Tests: parse_uniprot_dat ─────────────────────────────────────────────────

class TestUniprotDat:
    def test_signal_peptide_extracted(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P12345", features=_signal_feat(1, 25))])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert len(df) == 1
        row = df.iloc[0]
        assert row["Accession"] == "P12345"
        assert row["Type"] == "Signal peptide"
        assert row["Start"] == 1
        assert row["End"] == 25

    def test_transmembrane_and_binding(self, tmp_path):
        feats = _transmem_feat(30, 55, "TM1") + _binding_feat(70, 70, "Zn(2+)")
        dat = _make_dat(tmp_path, [_entry("P12346", features=feats)])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        types = set(df["Type"])
        assert "Transmembrane" in types
        assert "Binding site" in types
        tm = df[df["Type"] == "Transmembrane"].iloc[0]
        assert tm["Start"] == 30
        assert tm["End"] == 55
        bs = df[df["Type"] == "Binding site"].iloc[0]
        assert bs["Ligand"] == "Zn(2+)"

    def test_multiple_entries(self, tmp_path):
        dat = _make_dat(tmp_path, [
            _entry("P00001", features=_signal_feat(1, 20)),
            _entry("P00002", features=_transmem_feat(10, 30)),
        ])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert set(df["Accession"]) == {"P00001", "P00002"}

    def test_accession_filter(self, tmp_path):
        dat = _make_dat(tmp_path, [
            _entry("P00001", features=_signal_feat()),
            _entry("P00002", features=_transmem_feat()),
        ])
        r = _run(tmp_path, dat, accessions=["P00001"])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert list(df["Accession"]) == ["P00001"]

    def test_no_matching_features_empty_output(self, tmp_path):
        # Entry with only MOD_RES (not in ROI or BINDING types)
        feats = "FT   MOD_RES         50..50\nFT                   /note=\"Phosphoserine\"\n"
        dat = _make_dat(tmp_path, [_entry("P99999", features=feats)])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert len(df) == 0

    def test_uncertain_positions_skipped(self, tmp_path):
        # Uncertain start/end (with ?) should be skipped
        feats = "FT   SIGNAL          ?..25\n"
        dat = _make_dat(tmp_path, [_entry("P77777", features=feats)])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert len(df) == 0

    def test_lt_gt_positions_stripped(self, tmp_path):
        # <1..>25 means uncertain but present — we strip < and > and keep
        feats = "FT   TRANSMEM        <10..>30\n"
        dat = _make_dat(tmp_path, [_entry("P66666", features=feats)])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["Start"] == 10
        assert df.iloc[0]["End"] == 30

    def test_empty_dat_produces_empty_output(self, tmp_path):
        dat = _make_dat(tmp_path, [])
        r = _run(tmp_path, dat)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "uniprot_features.tsv", sep="\t")
        assert len(df) == 0


# ── Tests: parse_interpro_pfam ────────────────────────────────────────────────

class TestInterproPfam:
    def test_pfam_rows_extracted(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P12345")])
        ipr = _make_interpro(tmp_path, [("P12345", "PF00001", "7tm_1", 5, 90)])
        r = _run(tmp_path, dat, interpro_gz=ipr)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "pfam_domains.tsv", sep="\t")
        assert len(df) == 1
        row = df.iloc[0]
        assert row["Accession"] == "P12345"
        assert row["hmm_acc"] == "PF00001"
        assert row["hmm_name"] == "7tm_1"
        assert row["start"] == 5
        assert row["end"] == 90
        assert row["type"] == "Pfam"

    def test_non_pfam_analysis_skipped(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P22222")])
        gz = tmp_path / "ipr_nopfam.dat.gz"
        with gzip.open(gz, 'wt') as f:
            f.write("P22222\tmd5\t100\tCDD\tcd00001\tCDD domain\t1\t50\t1e-5\tT\t01-01-2024\tIPR000001\tTest\n")
        r = _run(tmp_path, dat, interpro_gz=gz)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "pfam_domains.tsv", sep="\t")
        assert len(df) == 0

    def test_accession_filter_interpro(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P11111"), _entry("P22222")])
        ipr = _make_interpro(tmp_path, [
            ("P11111", "PF00001", "Domain_A", 1, 50),
            ("P22222", "PF00002", "Domain_B", 1, 60),
        ])
        r = _run(tmp_path, dat, interpro_gz=ipr, accessions=["P11111"])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "pfam_domains.tsv", sep="\t")
        assert list(df["Accession"]) == ["P11111"]

    def test_multiple_domains_per_protein(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P33333")])
        ipr = _make_interpro(tmp_path, [
            ("P33333", "PF00001", "Domain_A", 1, 40),
            ("P33333", "PF00002", "Domain_B", 50, 90),
        ])
        r = _run(tmp_path, dat, interpro_gz=ipr)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "pfam_domains.tsv", sep="\t")
        assert len(df) == 2

    def test_no_interpro_file_gives_empty_pfam(self, tmp_path):
        dat = _make_dat(tmp_path, [_entry("P44444", dr_lines=_pfam_dr("PF00001", "7tm_1"))])
        r = _run(tmp_path, dat)  # no interpro_gz
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "out" / "pfam_domains.tsv", sep="\t")
        # No positions available without interpro file
        assert len(df) == 0
