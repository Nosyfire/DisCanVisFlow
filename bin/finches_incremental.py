#!/usr/bin/env python3
"""
finches_incremental.py — exact, fast site-wise FINCHES Δε saturation engine.

Background
----------
The naive saturation scan calls ``calculate_epsilon_value(mut, mut)`` once per
variant, rebuilding the entire L×L weighted interaction matrix every time. That
is O(19·L³) per protein (19 substitutions × L positions × O(L²) matrix build) and
is what makes a full-proteome run take ~2 weeks.

A single-point substitution at position ``p`` only perturbs a small, detectable
band of the L×L matrix, because the FINCHES self-epsilon reduces to a
per-element sum:

    w[i,j] = base[i,j] · (1 − repulsive_charge_mask[i,j]·cp) · aliphatic_mask[i,j]
    ε      = (1/L) · Σ_ij h(w[i,j])           with   h(x) = x·(x≠b) − 2b

    (base            : identity-only 20×20 pairwise lookup, convert_to_custom)
    (repulsive mask  : nonzero only where both residues are charged (R,K,E,D);
                       value = |NCPR/FCR| over the ±1 windows of i and j
                       = |Δcharge| / N_charged  — independent of window length)
    (aliphatic mask  : multiplier keyed on local aliphatic-cluster groups)
    (b = null_interaction_baseline, cp = charge_prefactor)

so ε is a normalised sum of an elementwise function of ``w``. A mutation at ``p``
changes only:
  * row/col ``p`` of ``base`` (residue identity),
  * rows/cols within ±1 of ``p`` of the charge mask,
  * whichever positions change aliphatic-group (found by diffing the 1-D group
    vector — exact regardless of cluster length).

Because ``w`` is symmetric for the homotypic (self) case, Δε is obtained from the
affected rows only, in O(L·|affected|) per variant → O(20·L²) per protein instead
of O(19·L³). Every number is reproduced from FINCHES' own primitives (the base
lookup, ``get_neighbors_window_of3``, ``get_aliphatic_groups``), and the result is
validated bit-for-bit against ``calculate_epsilon_value`` by the test-suite and by
the worker's ``--validate`` gate.

This module is import-only (no argparse); the worker drives it.
"""

from __future__ import annotations

import numpy as np

AA20 = "ACDEFGHIKLMNPQRSTVWY"
_AA_SET = set(AA20)
_POS = set("RK")
_NEG = set("ED")
# aliphatic-cluster multiplier table (mirrors get_aliphatic_weighted_mask)
_ALI_MULT = {(1, 1): 1.0, (1, 2): 1.0, (1, 3): 1.0, (2, 1): 1.0, (3, 1): 1.0,
             (2, 2): 1.5, (2, 3): 1.5, (3, 2): 1.5, (3, 3): 3.0}


class IncrementalEpsilon:
    """Holds model-level constants and computes per-protein saturation Δε.

    One instance per worker process (the FINCHES model is not picklable, so it
    is built inside each worker and wrapped here).
    """

    def __init__(self, X):
        # X = finches InteractionMatrixConstructor (already initialised)
        from finches.parsing_aminoacid_sequences import get_aliphatic_groups
        from finches.sequence_tools import get_neighbors_window_of3

        self.X = X
        self.b = float(X.null_interaction_baseline)
        self.cp = float(X.charge_prefactor)
        self._ali_groups = get_aliphatic_groups
        self._win3 = get_neighbors_window_of3

        # identity-only base pairwise lookup (convert_to_custom=True), 20×20
        M = np.asarray(
            X.calculate_pairwise_heterotypic_matrix(AA20, AA20, convert_to_custom=True),
            dtype=float,
        )
        self.Marr = M
        self.code = {a: i for i, a in enumerate(AA20)}

    # -- elementwise reduction kernel -------------------------------------
    def _h(self, w: np.ndarray) -> np.ndarray:
        """h(x) = x·(x≠b) − 2b  (exact match of get_attractive_repulsive_matrices)."""
        b = self.b
        return np.where(w != b, w, 0.0) - 2.0 * b

    # -- per-position charge window counts (±1) ---------------------------
    def _charge_window_counts(self, seq: str):
        """pos[k], neg[k] = #(+), #(−) charged residues in the ±1 window of k."""
        L = len(seq)
        pos = np.zeros(L, dtype=float)
        neg = np.zeros(L, dtype=float)
        for k in range(L):
            frag = self._win3(k, seq)
            pos[k] = sum(c in _POS for c in frag)
            neg[k] = sum(c in _NEG for c in frag)
        return pos, neg

    # -- weighted-matrix row for one row index ----------------------------
    def _w_row(self, i, codes, charged, win_pos, win_neg, groups):
        """Full length-L row w[i, :] of the weighted matrix for a sequence.

        codes    : int array of AA codes (len L)
        charged  : bool array, True where residue is R/K/E/D
        win_pos/win_neg : ±1-window +/- charge counts per position
        groups   : aliphatic group per position (0/1/2/3)
        """
        L = codes.shape[0]
        base = self.Marr[codes[i], codes]                    # (L,) identity lookup

        # repulsive charge mask row: nonzero only where both i and j charged
        rmask = np.zeros(L, dtype=float)
        if charged[i]:
            both = charged                                    # (L,) bool
            tot_pos = win_pos[i] + win_pos                    # (L,)
            tot_neg = win_neg[i] + win_neg
            denom = tot_pos + tot_neg                         # N_charged of fragment
            with np.errstate(divide="ignore", invalid="ignore"):
                val = np.abs(tot_pos - tot_neg) / denom
            rmask = np.where(both & (denom > 0), val, 0.0)

        # aliphatic mask row
        amask = np.ones(L, dtype=float)
        gi = groups[i]
        if gi > 0:
            for j in range(L):
                gj = groups[j]
                if gj > 0:
                    amask[j] = _ALI_MULT[(gi, gj)]

        return base * (1.0 - rmask * self.cp) * amask

    # -- WT setup (once per protein) --------------------------------------
    def prepare_wt(self, seq: str):
        """Return a context dict with cached WT quantities, or None if the
        sequence contains a non-standard residue (caller should fall back)."""
        if any(c not in _AA_SET for c in seq):
            return None
        L = len(seq)
        codes = np.fromiter((self.code[c] for c in seq), dtype=int, count=L)
        charged = np.fromiter((c in _POS or c in _NEG for c in seq), dtype=bool, count=L)
        win_pos, win_neg = self._charge_window_counts(seq)
        groups = np.asarray(self._ali_groups(seq), dtype=int)

        # exact WT weighted matrix straight from the library
        w_wt = np.asarray(
            self.X.calculate_weighted_pairwise_matrix(seq, seq), dtype=float
        )
        h_wt = self._h(w_wt)
        S_wt = float(h_wt.sum())
        eps_wt = S_wt / L
        return {
            "seq": seq, "L": L, "codes": codes, "charged": charged,
            "win_pos": win_pos, "win_neg": win_neg, "groups": groups,
            "w_wt": w_wt, "h_wt": h_wt, "S_wt": S_wt, "eps_wt": eps_wt,
        }

    # -- Δε for one variant (p 0-based, mut_aa) ---------------------------
    def delta_for_variant(self, ctx, p: int, mut_aa: str):
        """Return (eps_mut, delta) for substituting position p (0-based) → mut_aa."""
        L = ctx["L"]
        seq = ctx["seq"]
        mut_seq = seq[:p] + mut_aa + seq[p + 1:]

        # mutant per-position quantities (cheap 1-D recomputes)
        codes = ctx["codes"].copy()
        codes[p] = self.code[mut_aa]
        charged = ctx["charged"].copy()
        charged[p] = mut_aa in _POS or mut_aa in _NEG
        win_pos, win_neg = self._charge_window_counts(mut_seq)
        groups = np.asarray(self._ali_groups(mut_seq), dtype=int)

        # affected row/col set S: base row/col p, charge ±1, aliphatic-group diffs
        g_wt = ctx["groups"]
        S = {p}
        if p - 1 >= 0:
            S.add(p - 1)
        if p + 1 < L:
            S.add(p + 1)
        S.update(np.nonzero(groups != g_wt)[0].tolist())
        S = sorted(S)

        w_wt = ctx["w_wt"]
        h_wt = ctx["h_wt"]
        delta_sum = 0.0
        # rows i∈S, all columns; w is symmetric so col-block == row-block
        rows_new = {}
        for i in S:
            w_row = self._w_row(i, codes, charged, win_pos, win_neg, groups)
            rows_new[i] = w_row
            dh = self._h(w_row) - h_wt[i, :]
            delta_sum += 2.0 * float(dh.sum())
        # subtract doubly-counted S×S block
        for i in S:
            wr = rows_new[i]
            for j in S:
                delta_sum -= float(self._h(np.array([wr[j]]))[0] - h_wt[i, j])

        eps_mut = (ctx["S_wt"] + delta_sum) / L
        return eps_mut, eps_mut - ctx["eps_wt"]
