"""Tests for bin/create_finches_worker.py using a mocked finches backend."""
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_finches_worker.py"
OUT_COLS = ["Protein_ID", "Position", "WT_AA", "Mut_AA",
            "WT_Epsilon", "Mut_Epsilon", "Delta_Epsilon"]


def _fake_finches(libdir: Path):
    """Write a minimal fake `finches` package: epsilon = -(# of F/Y/W).

    Mirrors the real Mpipi forcefield's alphabet: it parameterises the standard
    20 plus U (selenocysteine), and raises KeyError on anything else (e.g. X).
    """
    pkg = libdir / "finches"
    (pkg / "forcefields").mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "forcefields" / "__init__.py").write_text("")
    (pkg / "epsilon_calculation.py").write_text(textwrap.dedent("""
        class InteractionMatrixConstructor:
            def __init__(self, parameters=None):
                self.p = parameters
            def calculate_epsilon_value(self, a, b):
                for c in a:
                    if c not in 'ACDEFGHIKLMNPQRSTVWYU':
                        raise KeyError(c)
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
         "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1",
         "--engine", "full"],
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


def test_unparameterised_residue_is_skipped_not_fatal(tmp_path):
    """X has no Mpipi parameters: skip that protein, keep going.

    Regression: the worker used to let the KeyError escape the pool, killing a
    whole proteome run because of one bad isoform.
    """
    lib = tmp_path / "lib"
    lib.mkdir()
    _fake_finches(lib)

    seq = tmp_path / "seq.tsv"
    pd.DataFrame(
        [["UNK-201", "yes", "AFXA"],   # X — unknown residue, no parameters
         ["OK-201",  "yes", "AFA"]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)

    res = subprocess.run(
        [sys.executable, str(BIN),
         "--loc_chrom", str(seq), "--output_dir", str(tmp_path),
         "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1",
         "--engine", "full"],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr

    out = pd.read_csv(tmp_path / "finches_saturation.tsv", sep="\t")
    # the clean protein is still scanned; the X one is dropped entirely
    assert set(out["Protein_ID"]) == {"OK-201"}


def test_selenoprotein_is_scored_at_standard_positions(tmp_path):
    """Mpipi DOES parameterise U, so selenoproteins must be scored, not skipped.

    Regression: a blanket 'any non-standard residue → skip the protein' guard
    wrongly dropped all 25 U-containing main isoforms. U is kept in the sequence
    context (it contributes to epsilon) but is not itself mutated, since the
    substitution alphabet is the standard 20.
    """
    lib = tmp_path / "lib"
    lib.mkdir()
    _fake_finches(lib)

    seq = tmp_path / "seq.tsv"
    pd.DataFrame([["SEL-201", "yes", "AFUA"]],   # U at position 3
                 columns=["Protein_ID", "main_isoform", "Sequence"]).to_csv(
        seq, sep="\t", index=False)

    res = subprocess.run(
        [sys.executable, str(BIN),
         "--loc_chrom", str(seq), "--output_dir", str(tmp_path),
         "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1",
         "--engine", "full"],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr

    out = pd.read_csv(tmp_path / "finches_saturation.tsv", sep="\t")
    assert set(out["Protein_ID"]) == {"SEL-201"}
    # 3 standard positions (1, 2, 4) x 19 substitutions; U at 3 is never mutated
    assert sorted(out["Position"].unique()) == [1, 2, 4]
    assert len(out) == 19 * 3
    # U still counts as context: WT epsilon reflects the single F
    assert (out["WT_Epsilon"] == -1.0).all()


def test_unusable_engine_fails_fast_instead_of_hanging(tmp_path):
    """A model that can't initialise must exit, not deadlock.

    Regression: the Pool initializer called sys.exit(1) on init failure, which
    makes multiprocessing respawn workers forever — the run hung silently. The
    model is now built in the parent first. The fake backend lacks the internals
    the incremental engine needs, so it stands in for a broken env here.
    """
    lib = tmp_path / "lib"
    lib.mkdir()
    _fake_finches(lib)

    seq = tmp_path / "seq.tsv"
    pd.DataFrame([["A-201", "yes", "AFA"]],
                 columns=["Protein_ID", "main_isoform", "Sequence"]).to_csv(
        seq, sep="\t", index=False)

    res = subprocess.run(
        [sys.executable, str(BIN),
         "--loc_chrom", str(seq), "--output_dir", str(tmp_path),
         "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1",
         "--engine", "incremental"],
        capture_output=True, text=True, timeout=120)   # times out if it hangs
    assert res.returncode != 0
    assert "FINCHES init failed" in res.stderr


def test_resume_appends_and_skips_completed(tmp_path):
    """Re-running must keep existing proteins and only add the missing ones."""
    lib = tmp_path / "lib"
    lib.mkdir()
    _fake_finches(lib)

    seq = tmp_path / "seq.tsv"
    pd.DataFrame(
        [["A-201", "yes", "AFA"],
         ["B-201", "yes", "AWA"]],
        columns=["Protein_ID", "main_isoform", "Sequence"]
    ).to_csv(seq, sep="\t", index=False)

    cmd = [sys.executable, str(BIN),
           "--loc_chrom", str(seq), "--output_dir", str(tmp_path),
           "--finches_lib", str(lib), "--only_main_isoforms", "--n_cpu", "1",
           "--engine", "full"]

    # first pass: only A-201 in the table
    only_a = tmp_path / "seq_a.tsv"
    pd.DataFrame([["A-201", "yes", "AFA"]],
                 columns=["Protein_ID", "main_isoform", "Sequence"]).to_csv(
        only_a, sep="\t", index=False)
    res = subprocess.run(
        [sys.executable, str(BIN), "--loc_chrom", str(only_a),
         "--output_dir", str(tmp_path), "--finches_lib", str(lib),
         "--only_main_isoforms", "--n_cpu", "1", "--engine", "full"],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    first = pd.read_csv(tmp_path / "finches_saturation.tsv", sep="\t")
    assert set(first["Protein_ID"]) == {"A-201"}

    # second pass with both: A-201 must be preserved, B-201 appended
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    both = pd.read_csv(tmp_path / "finches_saturation.tsv", sep="\t")
    assert set(both["Protein_ID"]) == {"A-201", "B-201"}
    # A-201's rows are the originals, not recomputed duplicates
    assert len(both[both["Protein_ID"] == "A-201"]) == len(first)


def test_blas_threading_pinned_before_numpy_import(tmp_path):
    """BLAS env vars must be set before numpy is imported, or the affinity fix is dead.

    Regression: OpenBLAS pins every forked worker to the core its thread pool
    starts on, which collapsed a 64-worker proteome run onto 2 logical CPUs
    (3.4% CPU each, 27x slower). Setting these before the numpy import stops the
    thread pool — and the pinning — from ever happening.
    """
    probe = tmp_path / "probe.py"
    probe.write_text(
        "import importlib.util, os, sys\n"
        f"spec = importlib.util.spec_from_file_location('cfw', {str(BIN)!r})\n"
        "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
        "print(os.environ.get('OPENBLAS_NUM_THREADS'), os.environ.get('OMP_NUM_THREADS'))\n"
        "print(len(os.sched_getaffinity(0)), os.cpu_count())\n"
    )
    env = {k: v for k, v in os.environ.items()
           if k not in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")}
    res = subprocess.run([sys.executable, str(probe)], capture_output=True,
                         text=True, env=env)
    assert res.returncode == 0, res.stderr
    threads, affinity = res.stdout.strip().splitlines()[:2]
    assert threads == "1 1", f"BLAS threading not pinned: {threads}"
    visible, total = affinity.split()
    assert visible == total, f"module import narrowed affinity: {visible}/{total}"
