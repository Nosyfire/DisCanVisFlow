"""Tests for create_disprot_worker.py (DisProt curated disorder regions)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd

WORKER = Path(__file__).parent.parent / "bin" / "create_disprot_worker.py"

# Deterministic 60-aa sequences; region sequences are derived by slicing so the
# worker's coordinate-validation always sees a genuine match/mismatch.
SEQ1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA"
SEQ1_ALT = "MKTAYIAKQRWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"  # diverges after res 10
SEQ2 = "GSHMASMTGGQQMGRGSEFMKRISTTITTTITITTGNGAGKALEEVLSKGNITTPTQINSS"

DISPROT_HEADER = (
    "UniProt ACC\tDisProt ID\tRegion ID\tStart\tEnd\t"
    "Term namespace\tTerm ID\tTerm name\tECO Term ID\tPMID\t"
    "Region sequence\tObsolete\tDataset\n"
)


def _run(args, tmpdir):
    return subprocess.run([sys.executable, str(WORKER)] + args,
                          capture_output=True, text=True, cwd=tmpdir)


def _seq(tmpdir):
    p = tmpdir / "seq.tsv"
    p.write_text(
        "Protein_ID\tEntry_Isoform\tmain_isoform\tSequence\n"
        f"GENE1-201\tP11111\tyes\t{SEQ1}\n"
        f"GENE1-204\tP11111-2\tno\t{SEQ1_ALT}\n"
        f"GENE2-201\tP22222\tyes\t{SEQ2}\n",
        encoding="utf-8",
    )
    return p


def _region(acc, dpid, rid, start, end, seq, term_ns="Structural state",
            term_id="IDPO:0000002", term_name="disorder",
            eco="ECO:0006220", pmid="pmid:123", obsolete="false",
            dataset="Human proteins"):
    region_seq = seq[start - 1:end]
    return (f"{acc}\t{dpid}\t{rid}\t{start}\t{end}\t{term_ns}\t{term_id}\t"
            f"{term_name}\t{eco}\t{pmid}\t{region_seq}\t{obsolete}\t{dataset}\n")


def _disprot(tmpdir, rows):
    p = tmpdir / "disprot_in.tsv"
    p.write_text(DISPROT_HEADER + "".join(rows), encoding="utf-8")
    return p


def _out(tmpdir):
    return pd.read_csv(tmpdir / "disprot.tsv", sep="\t", dtype=str)


class TestBasicOutput:
    def test_file_and_columns(self, tmp_path):
        dis = _disprot(tmp_path, [_region("P11111", "DP01", "DP01r001", 10, 30, SEQ1)])
        r = _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        df = _out(tmp_path)
        for col in ["Protein_ID", "Entry_Isoform", "disprot_id", "region_id",
                    "start", "end", "term_namespace", "term_id", "term_name",
                    "eco_id", "pmid", "dataset"]:
            assert col in df.columns, f"Missing: {col}"

    def test_region_mapped_to_pid(self, tmp_path):
        dis = _disprot(tmp_path, [_region("P11111", "DP01", "DP01r001", 10, 30, SEQ1)])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        df = _out(tmp_path)
        row = df[df["Protein_ID"] == "GENE1-201"].iloc[0]
        assert int(row["start"]) == 10 and int(row["end"]) == 30
        assert row["disprot_id"] == "DP01"

    def test_term_fields_preserved(self, tmp_path):
        dis = _disprot(tmp_path, [_region(
            "P22222", "DP02", "DP02r001", 5, 25, SEQ2,
            term_ns="Disorder function", term_id="GO:0005515",
            term_name="protein binding", eco="ECO:0000269", pmid="pmid:999",
            dataset="Viral proteins")])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        row = _out(tmp_path).iloc[0]
        assert row["term_namespace"] == "Disorder function"
        assert row["term_id"] == "GO:0005515"
        assert row["term_name"] == "protein binding"
        assert row["eco_id"] == "ECO:0000269"
        assert row["pmid"] == "pmid:999"
        assert row["dataset"] == "Viral proteins"


class TestCoordinateValidation:
    def test_mismatched_isoform_rejected(self, tmp_path):
        """P11111 region 15-30 matches SEQ1 (GENE1-201) but not SEQ1_ALT
        (GENE1-204, diverges after res 10) → only the canonical isoform is kept."""
        dis = _disprot(tmp_path, [_region("P11111", "DP01", "DP01r001", 15, 30, SEQ1)])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        pids = set(_out(tmp_path)["Protein_ID"].values)
        assert "GENE1-201" in pids
        assert "GENE1-204" not in pids

    def test_matching_isoform_kept(self, tmp_path):
        """A region within the identical N-term prefix (res 1-8) matches BOTH
        isoforms and is emitted for each."""
        dis = _disprot(tmp_path, [_region("P11111", "DP01", "DP01r001", 1, 8, SEQ1)])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        pids = set(_out(tmp_path)["Protein_ID"].values)
        assert {"GENE1-201", "GENE1-204"} <= pids

    def test_out_of_range_skipped(self, tmp_path):
        """End beyond sequence length → skipped (no crash)."""
        row = ("P22222\tDP02\tDP02r001\t50\t200\tStructural state\tIDPO:0000002\t"
               "disorder\tECO:0006220\tpmid:1\t\tfalse\tHuman proteins\n")
        dis = _disprot(tmp_path, [row])
        r = _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert "GENE2-201" not in set(_out(tmp_path)["Protein_ID"].values)


class TestFiltering:
    def test_obsolete_dropped(self, tmp_path):
        dis = _disprot(tmp_path, [
            _region("P11111", "DP01", "DP01r001", 10, 30, SEQ1, obsolete="true"),
            _region("P22222", "DP02", "DP02r001", 5, 25, SEQ2, obsolete="false"),
        ])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        pids = set(_out(tmp_path)["Protein_ID"].values)
        assert "GENE1-201" not in pids
        assert "GENE2-201" in pids

    def test_protein_not_in_disprot_excluded(self, tmp_path):
        dis = _disprot(tmp_path, [_region("P22222", "DP02", "DP02r001", 5, 25, SEQ2)])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path)], tmp_path)
        pids = set(_out(tmp_path)["Protein_ID"].values)
        assert "GENE1-201" not in pids
        assert "GENE2-201" in pids

    def test_only_main_isoforms(self, tmp_path):
        """--only_main_isoforms drops GENE1-204 (main_isoform=no) as a candidate."""
        dis = _disprot(tmp_path, [_region("P11111", "DP01", "DP01r001", 1, 8, SEQ1)])
        _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
              "--outdir", str(tmp_path), "--only_main_isoforms"], tmp_path)
        pids = set(_out(tmp_path)["Protein_ID"].values)
        assert "GENE1-201" in pids
        assert "GENE1-204" not in pids


class TestEdgeCases:
    def test_no_file_returns_empty(self, tmp_path):
        r = _run(["--seq_table", str(_seq(tmp_path)),
                  "--disprot_tsv", str(tmp_path / "NO_FILE"),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert len(_out(tmp_path)) == 0

    def test_empty_disprot_graceful(self, tmp_path):
        dis = tmp_path / "disprot_in.tsv"
        dis.write_text(DISPROT_HEADER, encoding="utf-8")
        r = _run(["--seq_table", str(_seq(tmp_path)), "--disprot_tsv", str(dis),
                  "--outdir", str(tmp_path)], tmp_path)
        assert r.returncode == 0, r.stderr
        assert len(_out(tmp_path)) == 0
