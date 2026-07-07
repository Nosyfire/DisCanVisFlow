"""Tests for bin/create_finches_worker.py using a mocked finches backend."""
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_finches_worker.py"
OUT_COLS = ["Protein_ID", "Position", "WT_AA", "Mut_AA",
            "WT_Epsilon", "Mut_Epsilon", "Delta_Epsilon"]


def _fake_finches(libdir: Path):
    """Write a minimal fake `finches` package: epsilon = -(# of F/Y/W)."""
    pkg = libdir / "finches"
    (pkg / "forcefields").mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "forcefields" / "__init__.py").write_text("")
    (pkg / "epsilon_calculation.py").write_text(textwrap.dedent("""
        class InteractionMatrixConstructor:
            def __init__(self, parameters=None):
                self.p = parameters
            def calculate_epsilon_value(self, a, b):
                return -float(sum(a.count(x) for x in 'FYW'))
    """))
    (pkg / "forcefields" / "mpipi.py").write_text(textwrap.dedent("""
        def Mpipi_model(version='Mpipi_GGv1'):
            return {'v': version}
    """))


def test_delta_epsilon_and_main_filter(tmp_path):
    lib = tmp_path / "lib"
    lib.mkdir()
    _fake_finches(lib)

    seq = tmp_path / "seq.tsv"
    # main isoform 'AFA' (one aromatic F) + an alt isoform that must be filtered out.
    pd.DataFrame(
        [["M-201", "yes", "AFA"],
         ["M-202", "no",  "AAAA"]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)

    res = subprocess.run(
        [sys.executable, str(BIN),
         "--loc_chrom", str(seq), "--output_dir", str(tmp_path),
         "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1"],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr

    out = pd.read_csv(tmp_path / "finches_saturation.tsv", sep="\t")
    assert list(out.columns) == OUT_COLS
    # only the main isoform survives the filter
    assert set(out["Protein_ID"]) == {"M-201"}
    # WT epsilon for 'AFA' = -1 (single F)
    assert (out["WT_Epsilon"] == -1.0).all()
    # mutating the F (pos 2) to A removes the aromatic → mut epsilon 0, Δε = +1
    f_to_a = out[(out["Position"] == 2) & (out["Mut_AA"] == "A")]
    assert len(f_to_a) == 1
    assert float(f_to_a.iloc[0]["Mut_Epsilon"]) == 0.0
    assert float(f_to_a.iloc[0]["Delta_Epsilon"]) == 1.0
