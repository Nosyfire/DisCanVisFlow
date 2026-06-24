"""Tests for create_pem_transfer_worker.py"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

WORKER = Path(__file__).parent.parent / "bin" / "create_pem_transfer_worker.py"


def _run(loc, pem, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--loc_chrom", str(loc),
         "--pem_tsv", str(pem),
         "--outdir", str(outdir)],
        capture_output=True, text=True,
    )


def test_transfers_motif_to_second_isoform(tmp_path):
    loc = tmp_path / "loc.tsv"
    pd.DataFrame({
        "Protein_ID": ["GENE-201", "GENE-202"],
        "Gene_Gencode": ["GENE", "GENE"],
        "Sequence": ["MAEAKLLPKL", "MAEAKXXX"],
    }).to_csv(loc, sep="\t", index=False)

    pem = tmp_path / "pem.tsv"
    pd.DataFrame({
        "Protein_ID": ["GENE-201"],
        "ELM_Accession": ["ELME000001"],
        "Start": [1],
        "End": [5],
    }).to_csv(pem, sep="\t", index=False)

    out = tmp_path / "out"
    out.mkdir()
    r = _run(loc, pem, out)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out / "pem_core_motifs_mapped.tsv", sep="\t")
    assert "GENE-202" in df["Protein_ID"].values
