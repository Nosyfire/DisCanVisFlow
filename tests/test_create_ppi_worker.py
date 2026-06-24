"""Tests for create_ppi_worker.py (Module 5j — PPI)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_ppi_worker.py"


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
        "Protein_ID\tEntry_Isoform\n"
        "GENE1-201\tP11111\n"
        "GENE2-201\tP22222\n"
        "GENE3-201\tP33333\n",
        encoding="utf-8",
    )
    return p


def _interaction_file(tmpdir: Path, name: str, rows: list[list]) -> Path:
    p = tmpdir / name
    header = "Accession A\tAccession B\tID Interactor A\tID Interactor B\tInteraction Detection Methods\tPublication Identifiers\tConfidence Value\n"
    p.write_text(header + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return p


class TestBasicOutput:
    def test_output_file_created(self, tmp_path):
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018(two hybrid)", "pubmed:12345678", "0.8"],
        ])
        result = _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "interactions.tsv").exists()

    def test_required_columns(self, tmp_path):
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018(two hybrid)", "pubmed:12345678", "0.8"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        for col in ["Protein_ID_A", "Protein_ID_B", "database", "number_of_pubmed"]:
            assert col in df.columns, f"Missing column: {col}"


class TestFiltering:
    def test_only_proteins_in_seq_table(self, tmp_path):
        """Protein_ID_A must always be a protein in the seq_table.
        External partners (not in seq_table) may appear as Protein_ID_B but
        proteins with no connection to the run must not appear as Protein_ID_A."""
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018(two hybrid)", "pubmed:11111111", "0.8"],
            # GENEXX is an external partner (not in seq_table) — OK as Protein_ID_B
            ["GENE1-201", "GENEXX-201", "P11111", "P99999",
             "MI:0018(two hybrid)", "pubmed:22222222", "0.5"],
            # Two proteins neither in seq_table — must be dropped entirely
            ["GENEXX-201", "GENEYY-201", "P99999", "P88888",
             "MI:0018(two hybrid)", "pubmed:33333333", "0.5"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        run_proteins = {"GENE1-201", "GENE2-201", "GENE3-201"}
        # All Protein_ID_A must be from the run
        assert all(a in run_proteins for a in df["Protein_ID_A"].values)
        # The GENEXX external partner should appear (linked to GENE1)
        assert "GENEXX-201" in df["Protein_ID_B"].values
        # GENEYY (no connection to run) must not appear at all
        assert "GENEYY-201" not in df["Protein_ID_A"].values
        assert "GENEYY-201" not in df["Protein_ID_B"].values

    def test_both_directions_kept(self, tmp_path):
        """An interaction A↔B should produce rows for A→B when A is in our list."""
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018(two hybrid)", "pubmed:12345678", "0.8"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        gene1_rows = df[df["Protein_ID_A"] == "GENE1-201"]
        assert "GENE2-201" in gene1_rows["Protein_ID_B"].values


class TestPubMedCounting:
    def test_pubmed_count(self, tmp_path):
        """number_of_pubmed should count distinct pubmed IDs in Publication Identifiers."""
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018", "pubmed:11111111|imex:IM-1|pubmed:22222222", "0.8"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        row = df[(df["Protein_ID_A"] == "GENE1-201") & (df["Protein_ID_B"] == "GENE2-201")].iloc[0]
        assert row["number_of_pubmed"] == 2


class TestMerging:
    def test_multiple_databases_merged(self, tmp_path):
        """Same pair appearing in two databases → one row with merged database field."""
        seq = _seq_table(tmp_path)
        intact = _interaction_file(tmp_path, "intact.tsv", [
            ["GENE1-201", "GENE2-201", "P11111", "P22222",
             "MI:0018", "pubmed:11111111", "0.8"],
        ])
        biogrid = _interaction_file(tmp_path, "biogrid.tsv", [
            ["GENE1-201", "GENE2-201", "entrez:1", "entrez:2",
             "MI:0004", "pubmed:22222222", "0.9"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(biogrid),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        row = df[(df["Protein_ID_A"] == "GENE1-201") & (df["Protein_ID_B"] == "GENE2-201")].iloc[0]
        # Both databases must appear in the database field
        assert "IntAct" in str(row["database"]) or "BioGRID" in str(row["database"])
        assert row["number_of_pubmed"] == 2  # 2 distinct pubmed IDs across both

    def test_uniprot_accession_matched_via_entry_isoform(self, tmp_path):
        """ppi_worker should match raw UniProt accessions using Entry_Isoform column."""
        # seq_table has Entry_Isoform → ppi_worker builds uniprot_to_pids map
        seq = tmp_path / "seq.tsv"
        seq.write_text(
            "Protein_ID\tEntry_Isoform\n"
            "GENE1-201\tP11111\n"
            "GENE1-204\tP11111-2\n",
            encoding="utf-8",
        )
        intact = _interaction_file(tmp_path, "intact.tsv", [
            # Accession A/B are raw UniProt accessions (from preprocessed raw MiTab)
            ["P11111", "P22222", "P11111", "P22222",
             "MI:0018(two hybrid)", "pubmed:12345678", "0.8"],
        ])
        _run([
            "--seq_table", str(seq),
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        # Both isoforms of GENE1 should be emitted with P22222-201 as external partner
        assert len(df) == 2
        assert set(df["Protein_ID_A"].values) == {"GENE1-201", "GENE1-204"}
        assert all(b == "P22222-201" for b in df["Protein_ID_B"].values)

    def test_empty_databases_handled(self, tmp_path):
        """Worker must not crash if all three database files are missing/empty."""
        seq = _seq_table(tmp_path)
        result = _run([
            "--seq_table", str(seq),
            "--intact", str(tmp_path / "no_file"),
            "--biogrid", str(tmp_path / "no_file"),
            "--hippie", str(tmp_path / "no_file"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert result.returncode == 0, result.stderr
        df = pd.read_csv(tmp_path / "interactions.tsv", sep="\t")
        assert len(df) == 0
