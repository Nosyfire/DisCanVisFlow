#!/usr/bin/env python3
"""
run_gopher_worker.py — GOPHER multi-level conservation RECOMPUTE.

GOPHER (SLiMSuite) builds per-protein orthologue alignments. This worker turns
those alignments into the per-residue, per-taxonomic-level conservation table the
pipeline consumes (the same schema as the pre-computed gopher conservation_table):

    uniprot_acc <TAB> level <TAB> conservation_score    (comma-separated per residue)

Two stages:
  1. (optional) run GOPHER to produce the orthologue alignments — supply a command
     template via --gopher_cmd (the user's SLiMSuite install + proteome DB). This
     is the slow part and is delegated to GOPHER itself; we don't reimplement it.
  2. score every alignment: for each taxonomic LEVEL, restrict the alignment to the
     sequences whose species belongs to that level, then compute a per-query-residue
     conservation score (Valdar-style: 1 - normalised Shannon entropy of the column,
     down-weighted by the column gap fraction), emitted only at query non-gap columns.

Levels (default, matching the legacy table): global, Mammalia, Vertebrata,
Eukaryota, Eumetazoa, Opisthokonta, Viridiplantae. `global` uses every sequence;
the others use the --taxon_map species→levels membership.

Usage:
  run_gopher_worker.py --seq_table seq.tsv --aln_dir gopher_aln/ --out conservation_table.tsv
      [--taxon_map species_levels.tsv] [--levels global,Mammalia,...]
      [--gopher_cmd "PYTHON slimsuite/tools/gopher.py seqin={fasta} orthdb=DB ..."]
"""

import argparse
import logging
import math
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

DEFAULT_LEVELS = ["global", "Mammalia", "Vertebrata", "Eukaryota",
                  "Eumetazoa", "Opisthokonta", "Viridiplantae"]

# Minimal built-in species→levels for common orthologue species codes. For a real
# run supply a complete --taxon_map; this default keeps the worker self-contained.
DEFAULT_TAXON = {
    "HUMAN": ["Mammalia", "Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "PANTR": ["Mammalia", "Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "MOUSE": ["Mammalia", "Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "RAT":   ["Mammalia", "Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "BOVIN": ["Mammalia", "Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "CHICK": ["Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "XENLA": ["Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "DANRE": ["Vertebrata", "Eukaryota", "Eumetazoa", "Opisthokonta"],
    "DROME": ["Eukaryota", "Eumetazoa", "Opisthokonta"],
    "CAEEL": ["Eukaryota", "Eumetazoa", "Opisthokonta"],
    "YEAST": ["Eukaryota", "Opisthokonta"],
    "ARATH": ["Eukaryota", "Viridiplantae"],
}


# ── FASTA + species helpers ─────────────────────────────────────────────────
def read_fasta(path):
    """Return list of (header, sequence) preserving order."""
    seqs, hdr, buf = [], None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    seqs.append((hdr, "".join(buf)))
                hdr, buf = line[1:].strip(), []
            elif hdr is not None:
                buf.append(line.strip())
    if hdr is not None:
        seqs.append((hdr, "".join(buf)))
    return seqs


def species_of(header):
    """Extract a species code from an alignment header.

    Handles 'ACC_SPECIES', 'sp|ACC|NAME_SPECIES', and trailing '[Species]'."""
    h = header.split()[0] if header else ""
    if "|" in h:                       # sp|P04637|P53_HUMAN
        h = h.split("|")[-1]
    if "_" in h:                       # NAME_HUMAN / P04637_HUMAN
        return h.rsplit("_", 1)[-1].upper()
    return h.upper()


def load_taxon_map(path):
    """species_code <TAB> comma-separated levels  → {code: set(levels)}."""
    if not path:
        return {k: set(v) for k, v in DEFAULT_TAXON.items()}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    code_col = df.columns[0]
    lvl_col = df.columns[1] if len(df.columns) > 1 else None
    out = {}
    for _, r in df.iterrows():
        code = str(r[code_col]).strip().upper()
        lvls = {x.strip() for x in str(r[lvl_col]).split(",") if x.strip()} if lvl_col else set()
        if code:
            out[code] = lvls
    return out


# ── conservation scoring ────────────────────────────────────────────────────
def column_conservation(symbols):
    """Per-column conservation in [0,1].

    Valdar-style simplification: 1 - normalised Shannon entropy over the non-gap
    residues, multiplied by (1 - gap_fraction). All-gap / single-symbol columns
    score 0 / 1 respectively. Fully conserved residue → 1.0."""
    n = len(symbols)
    if n == 0:
        return 0.0
    non_gap = [s for s in symbols if s not in ("-", ".", "")]
    gap_frac = (n - len(non_gap)) / n
    if not non_gap:
        return 0.0
    counts = defaultdict(int)
    for s in non_gap:
        counts[s.upper()] += 1
    total = len(non_gap)
    entropy = -sum((c / total) * math.log(c / total, 20) for c in counts.values())
    # entropy is in [0,1] because log base 20 (20 aa); 0 = identical column
    conservation = (1.0 - entropy) * (1.0 - gap_frac)
    return max(0.0, min(1.0, conservation))


def score_alignment(seqs, query_idx, member_indices):
    """Per-query-residue conservation over the member subsequence set.

    seqs: list of (header, aligned_seq). query_idx: index of the query (reference)
    sequence. member_indices: indices (incl. query) of sequences in this level.
    Returns a list of floats, one per NON-GAP residue of the query."""
    if not seqs:
        return []
    query = seqs[query_idx][1]
    members = [seqs[i][1] for i in member_indices] or [query]
    out = []
    for col, q_aa in enumerate(query):
        if q_aa in ("-", ".", ""):
            continue                       # query gap → not a query residue
        col_symbols = [m[col] for m in members if col < len(m)]
        out.append(round(column_conservation(col_symbols), 6))
    return out


def members_for_level(seqs, query_idx, level, taxon):
    """Indices of sequences belonging to `level` (query always included)."""
    idxs = [query_idx]
    for i, (hdr, _) in enumerate(seqs):
        if i == query_idx:
            continue
        if level == "global":
            idxs.append(i)
        else:
            sp = species_of(hdr)
            if level in taxon.get(sp, set()):
                idxs.append(i)
    return idxs


def acc_from_filename(name):
    """'P04637.orthaln.fas' / 'P04637.fas' → 'P04637'."""
    base = Path(name).name
    for suffix in (".orthaln.fas", ".orthaln.fasta", ".fas", ".fasta", ".aln"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base.rsplit(".", 1)[0]


def pick_query_index(seqs, acc):
    """Query = the sequence whose header contains the accession, else the first."""
    for i, (hdr, _) in enumerate(seqs):
        if acc and acc in hdr:
            return i
    return 0


def run_gopher(cmd_template, fasta, aln_dir):
    """Run the user-supplied GOPHER command to populate aln_dir (optional)."""
    cmd = cmd_template.format(fasta=shlex.quote(str(fasta)),
                              aln_dir=shlex.quote(str(aln_dir)))
    log.info("running GOPHER: %s", cmd)
    subprocess.run(cmd, shell=True, check=True)


def main():
    ap = argparse.ArgumentParser(description="GOPHER multi-level conservation recompute")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--out", required=True, help="conservation_table.tsv")
    ap.add_argument("--aln_dir", help="dir of per-accession orthologue alignments")
    ap.add_argument("--taxon_map", default="", help="species→levels TSV (default built-in)")
    ap.add_argument("--levels", default=",".join(DEFAULT_LEVELS))
    ap.add_argument("--gopher_cmd", default="",
                    help="command template to GENERATE --aln_dir first "
                         "({fasta}/{aln_dir} substituted); slow, runs the real GOPHER")
    ap.add_argument("--gopher_fasta", default="",
                    help="protein FASTA passed to --gopher_cmd (default: built from seq_table)")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    levels = [x.strip() for x in args.levels.split(",") if x.strip()]
    taxon = load_taxon_map(args.taxon_map)

    aln_dir = Path(args.aln_dir) if args.aln_dir else None

    # Optional: run GOPHER to build the alignments
    if args.gopher_cmd:
        seq_df = pd.read_csv(args.seq_table, sep="\t", dtype=str).fillna("")
        fasta = Path(args.gopher_fasta) if args.gopher_fasta else out_path.parent / "gopher_query.fasta"
        if not args.gopher_fasta:
            with open(fasta, "w") as fh:
                for _, r in seq_df.iterrows():
                    acc, seq = str(r.get("Entry_Isoform", "")), str(r.get("Sequence", ""))
                    if acc and seq and seq != "nan":
                        fh.write(f">{acc}\n{seq}\n")
        aln_dir = aln_dir or (out_path.parent / "gopher_aln")
        aln_dir.mkdir(parents=True, exist_ok=True)
        run_gopher(args.gopher_cmd, fasta, aln_dir)

    rows = []
    if aln_dir and aln_dir.exists():
        aln_files = sorted([p for p in aln_dir.iterdir()
                            if p.suffix in (".fas", ".fasta", ".aln")
                            or p.name.endswith(".orthaln.fas")])
        for p in aln_files:
            acc = acc_from_filename(p.name)
            seqs = read_fasta(p)
            if not seqs:
                continue
            q = pick_query_index(seqs, acc)
            for level in levels:
                members = members_for_level(seqs, q, level, taxon)
                scores = score_alignment(seqs, q, members)
                if scores:
                    rows.append({"uniprot_acc": acc, "level": level,
                                 "conservation_score": ", ".join(str(s) for s in scores)})
    else:
        log.warning("no --aln_dir with alignments — writing empty conservation table")

    df = pd.DataFrame(rows, columns=["uniprot_acc", "level", "conservation_score"])
    df.to_csv(out_path, sep="\t", index=False)
    log.info("GOPHER recompute: %d (acc,level) rows for %d accessions → %s",
             len(df), df["uniprot_acc"].nunique() if len(df) else 0, out_path)


if __name__ == "__main__":
    main()
