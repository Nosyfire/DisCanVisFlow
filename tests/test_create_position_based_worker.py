"""Tests for create_position_based_worker.py (Module 5m — PositionBasedAnnotations)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_position_based_worker.py"


def _run(args: list[str], tmpdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True, text=True, cwd=tmpdir,
    )


def _write(tmpdir: Path, name: str, content: str) -> Path:
    p = tmpdir / name
    p.write_text(content, encoding="utf-8")
    return p


def _seq_table(tmpdir: Path, rows: list[dict] | None = None) -> Path:
    if rows is None:
        rows = [
            {"Protein_ID": "GENE1-201", "Sequence": "ACDE"},   # 4 AA
            {"Protein_ID": "GENE2-201", "Sequence": "MAST"},   # 4 AA
        ]
    df = pd.DataFrame(rows)
    p = tmpdir / "seq.tsv"
    df.to_csv(p, sep="\t", index=False)
    return p


def _iupred(tmpdir: Path) -> Path:
    return _write(tmpdir, "IUPredscores.tsv",
                  "Protein_ID\tIUPredscores\n"
                  "GENE1-201\t0.1, 0.2, 0.3, 0.4\n"
                  "GENE2-201\t0.5, 0.6, 0.7, 0.8\n")


def _plddt(tmpdir: Path) -> Path:
    return _write(tmpdir, "AlphaFoldTable.tsv",
                  "Protein_ID\tPlldtscores\n"
                  "GENE1-201\t80.0, 60.0, 40.0, 20.0\n"
                  "GENE2-201\t90.0, 70.0, 50.0, 30.0\n")


def _combined_pos(tmpdir: Path) -> Path:
    lines = ["Protein_ID\tPosition\tCombinedDisorder"]
    for i in range(1, 5):
        lines.append(f"GENE1-201\t{i}\t{1 if i >= 3 else 0}")
        lines.append(f"GENE2-201\t{i}\t0")
    return _write(tmpdir, "CombinedDisorderNew_Pos.tsv", "\n".join(lines) + "\n")


def _run_basic(tmpdir: Path, extra_args: list[str] | None = None) -> subprocess.CompletedProcess:
    seq   = _seq_table(tmpdir)
    iup   = _iupred(tmpdir)
    plddt = _plddt(tmpdir)
    cdis  = _combined_pos(tmpdir)
    args  = [
        "--seq_table", str(seq),
        "--iupred_tsv", str(iup),
        "--plddt_tsv", str(plddt),
        "--combined_pos_tsv", str(cdis),
        "--outdir", str(tmpdir),
    ]
    if extra_args:
        args.extend(extra_args)
    return _run(args, tmpdir)


class TestBasicOutput:
    def test_output_files_created(self, tmp_path):
        r = _run_basic(tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "position_based_annotations.tsv").exists()
        assert (tmp_path / "rsa_scores.tsv").exists()

    def test_row_count_equals_total_sequence_length(self, tmp_path):
        r = _run_basic(tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        # 2 proteins × 4 AA each = 8 rows
        assert len(df) == 8

    def test_required_columns_present(self, tmp_path):
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        for col in ["Protein_ID", "position", "plddt", "rsa", "iupred",
                    "edisorder", "combineddisorder", "pfam"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_positions_are_one_indexed(self, tmp_path):
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        g1 = df[df["Protein_ID"] == "GENE1-201"]
        assert set(g1["position"].values) == {1, 2, 3, 4}


class TestRSAScores:
    def test_rsa_derived_from_plddt(self, tmp_path):
        """RSA = (100 - pLDDT) / 100."""
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        assert abs(row["plddt"] - 80.0) < 1e-6
        assert abs(row["rsa"]   - 0.20) < 1e-4

    def test_rsa_scores_file_columns(self, tmp_path):
        _run_basic(tmp_path)
        rsa_df = pd.read_csv(tmp_path / "rsa_scores.tsv", sep="\t")
        assert "Protein_ID" in rsa_df.columns
        assert "rsascores"  in rsa_df.columns

    def test_rsa_scores_comma_separated_string(self, tmp_path):
        _run_basic(tmp_path)
        rsa_df = pd.read_csv(tmp_path / "rsa_scores.tsv", sep="\t")
        row = rsa_df[rsa_df["Protein_ID"] == "GENE1-201"].iloc[0]
        vals = [float(x.strip()) for x in str(row["rsascores"]).split(",")]
        assert len(vals) == 4
        assert abs(vals[0] - 0.20) < 1e-4    # pos1: plddt=80 → rsa=0.20
        assert abs(vals[3] - 0.80) < 1e-4    # pos4: plddt=20 → rsa=0.80


class TestDisorderColumns:
    def test_combined_disorder_binary(self, tmp_path):
        """Positions 3,4 of GENE1-201 are in combined disorder regions."""
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        g1 = df[df["Protein_ID"] == "GENE1-201"].sort_values("position")
        assert g1.iloc[0]["combineddisorder"] == 0.0
        assert g1.iloc[2]["combineddisorder"] == 1.0
        assert g1.iloc[3]["combineddisorder"] == 1.0

    def test_edisorder_bool(self, tmp_path):
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        g1 = df[df["Protein_ID"] == "GENE1-201"].sort_values("position")
        assert not bool(g1.iloc[0]["edisorder"])
        assert bool(g1.iloc[2]["edisorder"])

    def test_iupred_values_correct(self, tmp_path):
        _run_basic(tmp_path)
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 2)].iloc[0]
        assert abs(row["iupred"] - 0.2) < 1e-4


class TestConservation:
    def test_phastcons_loaded(self, tmp_path):
        phastcons = _write(tmp_path, "phastcons.tsv",
                           "Protein_ID\tEntry_Isoform\tconservationscores\n"
                           "GENE1-201\tP11111\t0.90, 0.80, 0.70, 0.60\n")
        r = _run_basic(tmp_path, ["--phastcons_tsv", str(phastcons)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        assert abs(row["phastCons"] - 0.90) < 1e-4

    def test_gopher_conservation_levels(self, tmp_path):
        cons = _write(tmp_path, "conservation.tsv",
                      "Protein_ID\tEntry_Isoform\tlevel\tconservationscores\n"
                      "GENE1-201\tP11111\tglobal\t0.9, 0.8, 0.7, 0.6\n"
                      "GENE1-201\tP11111\tMammalia\t0.5, 0.4, 0.3, 0.2\n")
        r = _run_basic(tmp_path, ["--conservation_tsv", str(cons)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        assert abs(row["conservationGlobal"] - 0.9) < 1e-4
        assert abs(row["conservationMammal"] - 0.5) < 1e-4

    def test_missing_level_is_null(self, tmp_path):
        """Protein with no conservation data gets null values."""
        cons = _write(tmp_path, "conservation.tsv",
                      "Protein_ID\tEntry_Isoform\tlevel\tconservationscores\n"
                      "GENE1-201\tP11111\tglobal\t0.9, 0.8, 0.7, 0.6\n")
        r = _run_basic(tmp_path, ["--conservation_tsv", str(cons)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE2-201") & (df["position"] == 1)].iloc[0]
        assert pd.isna(row["conservationGlobal"])


class TestPfamAnnotation:
    def test_pfam_domain_at_position(self, tmp_path):
        pfam = _write(tmp_path, "pfam_domains.tsv",
                      "Protein_ID\tAccession\thmm_acc\thmm_name\ttype\tenvelope_start\tenvelope_end\n"
                      "GENE1-201\tP11111\tPF00001\tKinase\tDomain\t2\t3\n")
        r = _run_basic(tmp_path, ["--pfam_tsv", str(pfam)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row1 = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        row2 = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 2)].iloc[0]
        assert row1["pfam"] == "-"
        assert row2["pfam"] == "Kinase"

    def test_multiple_pfam_domains_joined(self, tmp_path):
        pfam = _write(tmp_path, "pfam_domains.tsv",
                      "Protein_ID\tAccession\thmm_acc\thmm_name\ttype\tenvelope_start\tenvelope_end\n"
                      "GENE1-201\tP11111\tPF00001\tKinase\tDomain\t1\t2\n"
                      "GENE1-201\tP11111\tPF00002\tSH2\tDomain\t1\t2\n")
        r = _run_basic(tmp_path, ["--pfam_tsv", str(pfam)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        assert "Kinase" in str(row["pfam"])
        assert "SH2" in str(row["pfam"])

    def test_repeat_type_excluded(self, tmp_path):
        """Non-Domain type Pfam entries should not appear in pfam column."""
        pfam = _write(tmp_path, "pfam_domains.tsv",
                      "Protein_ID\tAccession\thmm_acc\thmm_name\ttype\tenvelope_start\tenvelope_end\n"
                      "GENE1-201\tP11111\tPF00001\tRepeatEntry\tRepeat\t1\t4\n")
        r = _run_basic(tmp_path, ["--pfam_tsv", str(pfam)])
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE1-201") & (df["position"] == 1)].iloc[0]
        assert row["pfam"] == "-"


class TestEdgeCases:
    def test_no_file_sentinel_handled(self, tmp_path):
        """Worker must not crash when optional files are NO_FILE."""
        r = _run_basic(tmp_path, [
            "--phastcons_tsv",    str(tmp_path / "NO_FILE"),
            "--conservation_tsv", str(tmp_path / "NO_FILE"),
            "--pfam_tsv",         str(tmp_path / "NO_FILE"),
        ])
        assert r.returncode == 0, r.stderr

    def test_missing_iupred_yields_null(self, tmp_path):
        """Protein with no IUPred data gets null iupred values."""
        iupred = _write(tmp_path, "IUPredscores.tsv",
                        "Protein_ID\tIUPredscores\n"
                        "GENE1-201\t0.1, 0.2, 0.3, 0.4\n")  # GENE2-201 missing
        plddt = _plddt(tmp_path)
        cdis  = _combined_pos(tmp_path)
        seq   = _seq_table(tmp_path)
        r = _run([
            "--seq_table", str(seq),
            "--iupred_tsv", str(iupred),
            "--plddt_tsv", str(plddt),
            "--combined_pos_tsv", str(cdis),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        row = df[(df["Protein_ID"] == "GENE2-201") & (df["position"] == 1)].iloc[0]
        assert pd.isna(row["iupred"])

    def test_empty_seq_table(self, tmp_path):
        """Worker produces empty output files on empty seq_table."""
        seq = _write(tmp_path, "seq.tsv", "Protein_ID\tSequence\n")
        iup = _write(tmp_path, "IUPredscores.tsv", "Protein_ID\tIUPredscores\n")
        plt = _write(tmp_path, "AlphaFoldTable.tsv", "Protein_ID\tPlldtscores\n")
        cdis = _write(tmp_path, "CombinedDisorderNew_Pos.tsv",
                      "Protein_ID\tPosition\tCombinedDisorder\n")
        r = _run([
            "--seq_table", str(seq),
            "--iupred_tsv", str(iup),
            "--plddt_tsv", str(plt),
            "--combined_pos_tsv", str(cdis),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        assert len(df) == 0

    def test_duplicate_protein_ids_deduplicated(self, tmp_path):
        """Duplicate Protein_ID rows in seq_table produce one set of rows."""
        seq = _seq_table(tmp_path, [
            {"Protein_ID": "GENE1-201", "Sequence": "ACDE"},
            {"Protein_ID": "GENE1-201", "Sequence": "ACDE"},  # duplicate
        ])
        r = _run_basic(tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "position_based_annotations.tsv", sep="\t")
        assert len(df[df["Protein_ID"] == "GENE1-201"]) == 4  # 4 positions, not 8
