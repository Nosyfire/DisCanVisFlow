"""Tests for create_ppi_preprocess_worker.py (Module 5j-prep — PPI preprocessing)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

WORKER = Path(__file__).parent.parent / "bin" / "create_ppi_preprocess_worker.py"


def _run(args: list[str], tmpdir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(WORKER)] + args,
        capture_output=True, text=True, cwd=tmpdir,
    )


def _write(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


# ── IntAct MiTab2.7 fixtures ──────────────────────────────────────────────────
INTACT_HEADER = "\t".join([
    "#ID(s) interactor A", "ID(s) interactor B",
    "Alt. ID(s) interactor A", "Alt. ID(s) interactor B",
    "Alias(es) interactor A", "Alias(es) interactor B",
    "Interaction detection method(s)", "Publication 1st author(s)",
    "Publication Identifier(s)", "Taxon interactor A", "Taxon interactor B",
    "Interaction type(s)", "Source database(s)", "Interaction identifier(s)",
    "Confidence value(s)",
])


def _intact_row(acc_a, acc_b, taxon_a="taxid:9606(human)", taxon_b="taxid:9606(human)",
                method="psi-mi:\"MI:0018\"(two hybrid)",
                pubmed="pubmed:12345678|imex:IM-1",
                conf="intact-miscore:0.37"):
    parts = [
        f"uniprotkb:{acc_a}", f"uniprotkb:{acc_b}",
        "-", "-", "-", "-",
        method, "Smith J (2021)",
        pubmed, taxon_a, taxon_b,
        "psi-mi:\"MI:0915\"(physical association)", "psi-mi:\"MI:0469\"(IntAct)",
        "intact:EBI-12345",
        conf,
    ]
    return "\t".join(parts)


class TestIntActParsing:
    def test_output_files_created(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" + _intact_row("P11111", "P22222"))
        r = _run([
            "--intact", str(intact),
            "--biogrid", str(tmp_path / "NO_FILE"),
            "--hippie", str(tmp_path / "NO_FILE"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "Interaction_intact.tsv").exists()

    def test_columns_correct(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" + _intact_row("P11111", "P22222"))
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        for col in ["Accession A", "Accession B", "Publication Identifiers", "Confidence Value"]:
            assert col in df.columns

    def test_uniprot_stripped_of_prefix(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" + _intact_row("P11111", "P22222"))
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        assert df["Accession A"].iloc[0] == "P11111"
        assert df["Accession B"].iloc[0] == "P22222"

    def test_isoform_suffix_stripped(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" + _intact_row("P11111-2", "P22222-3"))
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        assert df["Accession A"].iloc[0] == "P11111"
        assert df["Accession B"].iloc[0] == "P22222"

    def test_non_human_filtered(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        mouse = "taxid:10090(mouse)"
        human = "taxid:9606(human)"
        content = (INTACT_HEADER + "\n"
                   + _intact_row("P11111", "P22222", taxon_a=human, taxon_b=human) + "\n"
                   + _intact_row("P33333", "P44444", taxon_a=mouse, taxon_b=human))
        _write(intact, content)
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        assert len(df) == 1
        assert df["Accession A"].iloc[0] == "P11111"

    def test_confidence_extracted(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" + _intact_row("P11111", "P22222",
                                                           conf="intact-miscore:0.75"))
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        assert float(df["Confidence Value"].iloc[0]) == pytest.approx(0.75)

    def test_pubmed_retained(self, tmp_path):
        intact = tmp_path / "intact.mitab"
        _write(intact, INTACT_HEADER + "\n" +
               _intact_row("P11111", "P22222", pubmed="pubmed:99999999|imex:IM-1"))
        _run(["--intact", str(intact), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_intact.tsv", sep="\t")
        assert "pubmed:99999999" in str(df["Publication Identifiers"].iloc[0])


# ── BioGRID MiTab fixtures ────────────────────────────────────────────────────
BIOGRID_HEADER = "\t".join([
    "#BioGRID Interaction ID", "Entrez Gene Interactor A", "Entrez Gene Interactor B",
    "BioGRID ID Interactor A", "BioGRID ID Interactor B",
    "Systematic Name(s) Interactor A", "Systematic Name(s) Interactor B",
    "Official Symbol Interactor A", "Official Symbol Interactor B",
    "Synonyms Interactor A", "Synonyms Interactor B",
    "Experimental System", "Experimental System Type",
    "Author", "Pubmed ID", "Organism Interactor A", "Organism Interactor B",
    "Throughput", "Score", "Modification", "Phenotypes", "Qualifications",
    "Tags", "Source Database",
    "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B",
])


def _biogrid_row(uniprot_a, uniprot_b, gene_a="GENEA", gene_b="GENEB",
                 org_a="9606", org_b="9606", pubmed="12345678", score="0.9"):
    parts = [
        "123456", "1", "2", "111", "222",
        "-", "-", gene_a, gene_b, "-", "-",
        "Co-purification", "physical", "Smith J (2021)",
        pubmed, org_a, org_b,
        "Low Throughput", score, "-", "-", "-", "-", "BioGRID",
        uniprot_a, uniprot_b,
    ]
    return "\t".join(parts)


class TestBioGRIDParsing:
    def test_output_created(self, tmp_path):
        bg = tmp_path / "biogrid.mitab"
        _write(bg, BIOGRID_HEADER + "\n" + _biogrid_row("P11111", "P22222"))
        r = _run([
            "--intact", str(tmp_path / "NO_FILE"),
            "--biogrid", str(bg),
            "--hippie", str(tmp_path / "NO_FILE"),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "Interaction_biogrid.tsv", sep="\t")
        assert len(df) == 1
        assert df["Accession A"].iloc[0] == "P11111"

    def test_non_human_filtered(self, tmp_path):
        bg = tmp_path / "biogrid.mitab"
        content = (BIOGRID_HEADER + "\n"
                   + _biogrid_row("P11111", "P22222", org_a="9606", org_b="9606") + "\n"
                   + _biogrid_row("P33333", "P44444", org_a="10090", org_b="9606"))
        _write(bg, content)
        _run(["--intact", str(tmp_path / "NO_FILE"), "--biogrid", str(bg),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_biogrid.tsv", sep="\t")
        assert len(df) == 1

    def test_pubmed_formatted(self, tmp_path):
        bg = tmp_path / "biogrid.mitab"
        _write(bg, BIOGRID_HEADER + "\n" + _biogrid_row("P11111", "P22222", pubmed="87654321"))
        _run(["--intact", str(tmp_path / "NO_FILE"), "--biogrid", str(bg),
              "--hippie", str(tmp_path / "NO_FILE"), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_biogrid.tsv", sep="\t")
        assert "pubmed:87654321" in str(df["Publication Identifiers"].iloc[0])


# ── HIPPIE fixture ────────────────────────────────────────────────────────────

class TestHIPPIEParsing:
    def test_output_created(self, tmp_path):
        hi = tmp_path / "hippie.txt"
        _write(hi, "GENEA\tP11111\tGENEB\tP22222\t0.63\texperiments:1,score:0.63,pmids:12345678\n")
        r = _run([
            "--intact", str(tmp_path / "NO_FILE"),
            "--biogrid", str(tmp_path / "NO_FILE"),
            "--hippie", str(hi),
            "--outdir", str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0, r.stderr
        df = pd.read_csv(tmp_path / "Interaction_hippie.tsv", sep="\t")
        assert len(df) == 1
        assert df["Accession A"].iloc[0] == "P11111"
        assert df["Accession B"].iloc[0] == "P22222"

    def test_pubmed_extracted(self, tmp_path):
        hi = tmp_path / "hippie.txt"
        _write(hi, "GENEA\tP11111\tGENEB\tP22222\t0.63\texperiments:1,pmids:99887766\n")
        _run(["--intact", str(tmp_path / "NO_FILE"), "--biogrid", str(tmp_path / "NO_FILE"),
              "--hippie", str(hi), "--outdir", str(tmp_path)], tmp_path)
        df = pd.read_csv(tmp_path / "Interaction_hippie.tsv", sep="\t")
        assert "pubmed:99887766" in str(df["Publication Identifiers"].iloc[0])


# ── Robustness ────────────────────────────────────────────────────────────────

class TestRobustness:
    def test_all_no_file(self, tmp_path):
        r = _run([
            "--intact",  str(tmp_path / "NO_FILE"),
            "--biogrid", str(tmp_path / "NO_FILE"),
            "--hippie",  str(tmp_path / "NO_FILE"),
            "--outdir",  str(tmp_path),
        ], tmp_path)
        assert r.returncode == 0
        for name in ["Interaction_intact.tsv", "Interaction_biogrid.tsv", "Interaction_hippie.tsv"]:
            df = pd.read_csv(tmp_path / name, sep="\t")
            assert len(df) == 0
