"""Tests for create_conservation_worker.py (Module 7)."""

import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_conservation_worker.py"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(args: list[str], tmpdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True,
        text=True,
        cwd=tmpdir,
    )


def _seq_table(tmpdir: Path) -> Path:
    p = tmpdir / "seq.tsv"
    p.write_text(
        "Protein_ID\tEntry_Isoform\tChromosome\n"
        "GENE1-201\tP12345\tchr1\n"
        "GENE1-202\tP12345-2\tchr1\n"
        "GENE2-201\tQ99999\tchr1\n",
        encoding="utf-8",
    )
    return p


def _conservation_table(tmpdir: Path) -> Path:
    """Minimal GOPHER conservation table with 3-residue scores."""
    p = tmpdir / "conservation_table.tsv"
    # P12345 canonical: 3 residues, 2 levels
    p.write_text(
        "uniprot_acc\tlevel\tconservation_score\n"
        "P12345\tglobal\t0.5, 0.8, 1.0\n"
        "P12345\tMammalia\t1.0, 1.0, 0.9\n",
        encoding="utf-8",
    )
    return p


def _combined_map(tmpdir: Path, strand: str = "+") -> Path:
    """Two-residue combined map for GENE1-201 on chr1."""
    p = tmpdir / "combined_map.map"
    p.write_text(
        "# header|id|x|y|GENE1-201|GENE1|100|CDS:1-9| chr1 " + strand + " 1000-1009\n"
        "0 M 1,2,3 ATG M 1000,1001,1002, ATG M\n"
        "1 E 4,5,6 GAG E 1003,1004,1005, GAG E\n"
        "2 H 7,8,9 CAC H 1006,1007,1008, CAC H\n",
        encoding="utf-8",
    )
    return p


def _phastcons_bedgraph(tmpdir: Path) -> Path:
    """Fake bigWig → BedGraph output covering chr1:1000-1009."""
    p = tmpdir / "fake.bedgraph"
    p.write_text(
        "chr1\t1000\t1003\t0.6\n"   # covers pos 1000,1001,1002 → residue 0 score 0.6
        "chr1\t1003\t1006\t0.8\n"   # residue 1
        "chr1\t1006\t1009\t0.9\n",  # residue 2
        encoding="utf-8",
    )
    return p


# ---------------------------------------------------------------------------
# Tests: GOPHER multiple-level conservation
# ---------------------------------------------------------------------------

class TestGopherConservation:
    def test_output_file_created(self, tmp_path):
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        result = _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(tmp_path / "combined_map.map"),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "conservation_multiple_level.tsv").exists()

    def test_canonical_acc_lookup(self, tmp_path):
        """P12345-2 isoform should resolve to P12345 canonical in GOPHER table."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(tmp_path / "combined_map.map"),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_multiple_level.tsv", sep="\t")
        # Both GENE1-201 (P12345) and GENE1-202 (P12345-2) should have scores
        gene1_pids = df[df["Protein_ID"].str.startswith("GENE1")]["Protein_ID"].unique()
        assert "GENE1-201" in gene1_pids
        assert "GENE1-202" in gene1_pids

    def test_missing_acc_skipped(self, tmp_path):
        """GENE2-201 (Q99999) has no entry in conservation table — row must be absent."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(tmp_path / "combined_map.map"),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_multiple_level.tsv", sep="\t")
        assert "GENE2-201" not in df["Protein_ID"].values

    def test_levels_present(self, tmp_path):
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(tmp_path / "combined_map.map"),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_multiple_level.tsv", sep="\t")
        levels = set(df[df["Protein_ID"] == "GENE1-201"]["level"])
        assert "global" in levels
        assert "Mammalia" in levels

    def test_scores_string_format(self, tmp_path):
        """conservationscores should be a comma-separated float string."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(tmp_path / "combined_map.map"),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_multiple_level.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["level"] == "global")].iloc[0]
        scores = [float(x) for x in row["conservationscores"].split(",")]
        assert scores == pytest.approx([0.5, 0.8, 1.0])


# ---------------------------------------------------------------------------
# Tests: phastCons
# ---------------------------------------------------------------------------

class TestPhastCons:
    def test_phastcons_output_created(self, tmp_path, monkeypatch):
        """With a fake BedGraph file passed via --phastcons_bedgraph, output is created."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        cmap = _combined_map(tmp_path)
        bg = _phastcons_bedgraph(tmp_path)

        result = _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(cmap),
            "--outdir", str(tmp_path),
            "--skip_gopher",
            "--phastcons_bedgraph", str(bg),
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "conservation_phastcons.tsv").exists()

    def test_phastcons_scores_mapped(self, tmp_path):
        """Per-residue phastCons score = mean of 3 nucleotide positions from BedGraph."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        cmap = _combined_map(tmp_path)
        bg = _phastcons_bedgraph(tmp_path)

        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(cmap),
            "--outdir", str(tmp_path),
            "--skip_gopher",
            "--phastcons_bedgraph", str(bg),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_phastcons.tsv", sep="\t")
        row = df[df["Protein_ID"] == "GENE1-201"].iloc[0]
        scores = [float(x) for x in row["conservationscores"].split(",")]
        assert len(scores) == 3
        assert scores[0] == pytest.approx(0.6)
        assert scores[1] == pytest.approx(0.8)
        assert scores[2] == pytest.approx(0.9)

    def test_phastcons_minus_strand(self, tmp_path):
        """Minus-strand protein: genomic coords are still resolved correctly."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        cmap = _combined_map(tmp_path, strand="-")
        bg = _phastcons_bedgraph(tmp_path)

        result = _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(cmap),
            "--outdir", str(tmp_path),
            "--skip_gopher",
            "--phastcons_bedgraph", str(bg),
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        df = pd.read_csv(tmp_path / "conservation_phastcons.tsv", sep="\t")
        assert "GENE1-201" in df["Protein_ID"].values

    def test_missing_coords_zero(self, tmp_path):
        """Residues with '-' genomic coords get score 0.0."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        p = tmp_path / "combined_map_gap.map"
        p.write_text(
            "# header|id|x|y|GENE1-201|GENE1|100|CDS:1-9| chr1 + 1000-1009\n"
            "0 M 1,2,3 ATG M -,-,-, ATG M\n"
            "1 E 4,5,6 GAG E 1003,1004,1005, GAG E\n",
            encoding="utf-8",
        )
        bg = _phastcons_bedgraph(tmp_path)
        _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(p),
            "--outdir", str(tmp_path),
            "--skip_gopher",
            "--phastcons_bedgraph", str(bg),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "conservation_phastcons.tsv", sep="\t")
        row = df[df["Protein_ID"] == "GENE1-201"].iloc[0]
        scores = [float(x) for x in row["conservationscores"].split(",")]
        assert scores[0] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_empty_combined_map(self, tmp_path):
        """Worker should not crash when combined_map.map has no data lines."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        p = tmp_path / "empty.map"
        p.write_text("", encoding="utf-8")
        result = _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(p),
            "--outdir", str(tmp_path),
            "--skip_phastcons",
        ], tmp_path)
        assert result.returncode == 0, result.stderr

    def test_skip_both_flags(self, tmp_path):
        """With --skip_gopher and --skip_phastcons, outputs are empty but created."""
        seq = _seq_table(tmp_path)
        cons = _conservation_table(tmp_path)
        p = tmp_path / "combined_map.map"
        p.write_text("", encoding="utf-8")
        result = _run([
            "--seq_table", str(seq),
            "--conservation_table", str(cons),
            "--combined_map", str(p),
            "--outdir", str(tmp_path),
            "--skip_gopher",
            "--skip_phastcons",
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "conservation_multiple_level.tsv").exists()
        assert (tmp_path / "conservation_phastcons.tsv").exists()
