"""Tests for split_seq_table.py — gene-balanced scatter chunking."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

WORKER = Path(__file__).parent.parent / "bin" / "split_seq_table.py"


def _run(loc, n, outdir, prefix="chunk_"):
    return subprocess.run(
        [sys.executable, str(WORKER), "--loc_chrom", str(loc),
         "--n_chunks", str(n), "--outdir", str(outdir), "--prefix", prefix],
        capture_output=True, text=True)


def _make(tmp, rows):
    p = tmp / "seq.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _rows(genes):
    """genes = {gene: n_isoforms}"""
    out = []
    for g, n in genes.items():
        for i in range(n):
            out.append({"Gene_Gencode": g, "Protein_ID": f"{g}-{201+i}",
                        "Entry_Isoform": f"P{g}", "Sequence": "MEEK"})
    return out


def _read_chunks(outdir, prefix="chunk_"):
    return sorted(Path(outdir).glob(f"{prefix}*.tsv"))


class TestSplit:
    def test_all_isoforms_preserved(self, tmp_path):
        loc = _make(tmp_path, _rows({"TP53": 9, "RAF1": 5, "BRAF": 3}))
        r = _run(loc, 3, tmp_path / "out")
        assert r.returncode == 0, r.stderr
        chunks = _read_chunks(tmp_path / "out")
        total = sum(len(pd.read_csv(c, sep="\t")) for c in chunks)
        assert total == 17  # 9+5+3, nothing lost or duplicated

    def test_gene_not_split_across_chunks(self, tmp_path):
        loc = _make(tmp_path, _rows({"TP53": 9, "RAF1": 5, "BRAF": 3, "KRAS": 2}))
        _run(loc, 4, tmp_path / "out")
        chunks = _read_chunks(tmp_path / "out")
        # each gene must appear in exactly one chunk
        gene_to_chunks = {}
        for c in chunks:
            for g in pd.read_csv(c, sep="\t")["Gene_Gencode"].unique():
                gene_to_chunks.setdefault(g, set()).add(c.name)
        assert all(len(v) == 1 for v in gene_to_chunks.values())

    def test_header_in_every_chunk(self, tmp_path):
        loc = _make(tmp_path, _rows({"TP53": 4, "RAF1": 4}))
        _run(loc, 2, tmp_path / "out")
        for c in _read_chunks(tmp_path / "out"):
            assert pd.read_csv(c, sep="\t").columns.tolist()[0] == "Gene_Gencode"

    def test_fewer_genes_than_chunks(self, tmp_path):
        loc = _make(tmp_path, _rows({"TP53": 2, "RAF1": 2}))
        _run(loc, 50, tmp_path / "out")
        # only 2 genes → at most 2 chunks, no empties
        chunks = _read_chunks(tmp_path / "out")
        assert 1 <= len(chunks) <= 2
        assert all(len(pd.read_csv(c, sep="\t")) > 0 for c in chunks)

    def test_single_chunk_equals_input(self, tmp_path):
        loc = _make(tmp_path, _rows({"TP53": 3, "RAF1": 2}))
        _run(loc, 1, tmp_path / "out")
        chunks = _read_chunks(tmp_path / "out")
        assert len(chunks) == 1
        assert len(pd.read_csv(chunks[0], sep="\t")) == 5

    def test_balanced_distribution(self, tmp_path):
        # 4 genes of 5 isoforms each into 2 chunks → ~10 each
        loc = _make(tmp_path, _rows({"A": 5, "B": 5, "C": 5, "D": 5}))
        _run(loc, 2, tmp_path / "out")
        sizes = [len(pd.read_csv(c, sep="\t")) for c in _read_chunks(tmp_path / "out")]
        assert max(sizes) - min(sizes) <= 5  # balanced within one gene-group
