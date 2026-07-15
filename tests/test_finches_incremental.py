"""Exactness tests for the incremental FINCHES Δε engine.

The incremental engine (bin/finches_incremental.py) must reproduce, to within
floating-point noise, the Δε that the reference path produces by calling
finches.calculate_epsilon_value(mut, mut) for every variant.

Skipped when the `finches` package is not importable.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

BIN = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN))

finches_available = importlib.util.find_spec("finches") is not None
pytestmark = pytest.mark.skipif(not finches_available, reason="finches not installed")


@pytest.fixture(scope="module")
def engine_and_model():
    from finches.epsilon_calculation import InteractionMatrixConstructor
    from finches.forcefields.mpipi import Mpipi_model
    from finches_incremental import IncrementalEpsilon

    X = InteractionMatrixConstructor(parameters=Mpipi_model(version="Mpipi_GGv1"))
    return IncrementalEpsilon(X), X


# sequences chosen to exercise: charged clusters, mixed charge, aliphatic runs,
# aromatics, termini, and a longer realistic IDR-ish stretch.
SEQS = [
    "MEHIQGAWKTISNGFGFKDAVFDGSSCISPTIVQQFGYQRRASDDGKLTD",  # realistic-ish
    "KKKKEEEERRRRDDDD",                                     # charge clusters
    "AAAALLLLVVVVIIII",                                     # aliphatic clusters
    "MGSSHHHHHHSSGLVPRGSHM",                                # tag-like
    "GSGSGSGSGSGSGSGSGSGS",                                 # low complexity, no charge/ali
    "DKDKDKDKAAAALLLLKRKR",                                 # mixed
]

AA = "ACDEFGHIKLMNPQRSTVWY"


def _ref_delta(X, seq, p, mut_aa):
    """Reference Δε via full library recompute."""
    wt_eps = float(X.calculate_epsilon_value(seq, seq))
    mut_seq = seq[:p] + mut_aa + seq[p + 1:]
    mut_eps = float(X.calculate_epsilon_value(mut_seq, mut_seq))
    return wt_eps, mut_eps, mut_eps - wt_eps


@pytest.mark.parametrize("seq", SEQS)
def test_wt_epsilon_matches_library(engine_and_model, seq):
    eng, X = engine_and_model
    ctx = eng.prepare_wt(seq)
    assert ctx is not None
    ref = float(X.calculate_epsilon_value(seq, seq))
    assert abs(ctx["eps_wt"] - ref) < 1e-9, (ctx["eps_wt"], ref)


@pytest.mark.parametrize("seq", SEQS)
def test_all_variants_match_library(engine_and_model, seq):
    """Exhaustive: every position × every substitution must match reference."""
    eng, X = engine_and_model
    ctx = eng.prepare_wt(seq)
    max_abs = 0.0
    worst = None
    for p, wt_aa in enumerate(seq):
        for mut_aa in AA:
            if mut_aa == wt_aa:
                continue
            _, mut_eps, delta = _ref_delta(X, seq, p, mut_aa)
            inc_mut, inc_delta = eng.delta_for_variant(ctx, p, mut_aa)
            d = abs(inc_mut - mut_eps)
            if d > max_abs:
                max_abs, worst = d, (p, wt_aa, mut_aa, mut_eps, inc_mut)
    assert max_abs < 1e-8, f"max abs diff {max_abs} at {worst}"


def test_non_standard_residue_returns_none(engine_and_model):
    eng, _ = engine_and_model
    assert eng.prepare_wt("ACDEFXGHIK") is None  # X is non-standard → caller falls back
