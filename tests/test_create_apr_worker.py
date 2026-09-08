"""Tests for bin/create_apr_worker.py (AGGRESCAN a3v aggregation-prone regions)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_apr_worker.py"
OUT_COLS = ["Protein_ID", "start", "end", "length", "mean_a3v", "peak_a3v"]

# Hydrophobic (high a3v) and hydrophilic (very negative a3v) building blocks.
HYDRO = "VILFVILFVILF"          # 12 aa, a3v ~ +1.5
POLAR = "DEKRNQDEKRNQDEKRNQ"    # 18 aa, a3v ~ -1.3


def _run(args):
    return subprocess.run([sys.executable, str(BIN)] + args,
                          capture_output=True, text=True)


def _seq_table(path: Path, rows):
    pd.DataFrame(rows, columns=["Protein_ID", "main_isoform", "Sequence"]).to_csv(
        path, sep="\t", index=False)


def _read(tmp_path):
    return pd.read_csv(tmp_path / "aggregation_prone_regions.tsv", sep="\t")


def test_finds_hydrophobic_run(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["APR-201", "yes", POLAR + HYDRO + POLAR]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--threshold", "0.5"])
    assert res.returncode == 0, res.stderr
    out = _read(tmp_path)
    assert list(out.columns) == OUT_COLS
    assert len(out) == 1
    row = out.iloc[0]
    # the hydrophobic block occupies 1-based positions 19..30
    assert 19 <= row["start"] <= 30
    assert 19 <= row["end"] <= 30
    assert row["end"] >= row["start"]
    assert row["length"] == row["end"] - row["start"] + 1
    assert row["mean_a3v"] > 0.5
    assert row["peak_a3v"] > 0.5


def test_no_apr_in_hydrophilic_protein(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["POLAR-201", "yes", POLAR * 3]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--threshold", "0.5"])
    assert res.returncode == 0, res.stderr
    out = _read(tmp_path)
    assert list(out.columns) == OUT_COLS
    assert len(out) == 0


def test_min_len_enforced(tmp_path):
    # A 3-residue hydrophobic patch cannot produce a >=5-position run.
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["SHORT-201", "yes", POLAR + "VIL" + POLAR]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--threshold", "0.5", "--min_len", "5"])
    assert res.returncode == 0, res.stderr
    assert len(_read(tmp_path)) == 0


def test_nonstandard_residue_breaks_the_run(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["CLEAN-201", "yes", POLAR + HYDRO + POLAR],
                     ["XBROKEN-201", "yes", POLAR + "VILFVXLFVILF" + POLAR]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--threshold", "0.5"])
    assert res.returncode == 0, res.stderr
    out = _read(tmp_path)
    clean = out[out["Protein_ID"] == "CLEAN-201"]
    broken = out[out["Protein_ID"] == "XBROKEN-201"]
    assert len(clean) == 1
    # the X blanks its whole window, so the broken protein yields fewer/shorter APRs
    assert broken["length"].sum() < clean["length"].sum()


def test_only_main_isoforms_filter(tmp_path):
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["A-201", "yes", POLAR + HYDRO + POLAR],
                     ["A-202", "no", POLAR + HYDRO + POLAR]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--threshold", "0.5", "--only_main_isoforms"])
    assert res.returncode == 0, res.stderr
    out = _read(tmp_path)
    assert set(out["Protein_ID"]) == {"A-201"}


def test_falls_back_to_default_threshold_on_small_input(tmp_path):
    """Too few proteins to calibrate a quantile -> DEFAULT_A3V_THRESHOLD + a warning."""
    seq = tmp_path / "seq.tsv"
    _seq_table(seq, [["APR-201", "yes", POLAR + HYDRO + POLAR]])
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path)])
    assert res.returncode == 0, res.stderr
    assert "DEFAULT" in res.stderr
    assert (tmp_path / "aggregation_prone_regions.tsv").exists()


def test_calibrates_threshold_when_given_enough_proteins(tmp_path):
    """With >= min_proteins_for_threshold sequences the q90 is computed, not defaulted."""
    seq = tmp_path / "seq.tsv"
    rows = [[f"P-{i}", "yes", POLAR + HYDRO + POLAR] for i in range(20)]
    _seq_table(seq, rows)
    res = _run(["--seq_table", str(seq), "--outdir", str(tmp_path),
                "--min_proteins_for_threshold", "10"])
    assert res.returncode == 0, res.stderr
    assert "q90" in res.stderr
    assert "DEFAULT" not in res.stderr
