#!/usr/bin/env python3
"""
create_pdb_bulk_worker.py — Module 5c (BULK): PDB structure → transcript mapping
                            from a single local SIFTS flat file (no per-protein API).

This is a drop-in *fast* replacement for create_pdb_worker.py's network path.
Instead of one PDBe graph-API call per UniProt accession
(https://www.ebi.ac.uk/pdbe/api/mappings/{accession} — ~tens of seconds each,
~10 days for the full human proteome), it does ONE local join against the SIFTS
flat file:

    https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz

columns: PDB  CHAIN  SP_PRIMARY  RES_BEG  RES_END  PDB_BEG  PDB_END  SP_BEG  SP_END
  • SP_PRIMARY  = canonical UniProt accession (P04637, never P04637-9)
  • SP_BEG/END  = UniProt residue range covered  -> unp_start / unp_end
  • CHAIN       = author chain id               -> chain_id
  • RES_BEG/END = SEQRES (label) residue numbers (used internally only)
  • PDB_BEG/END = PDB author residue numbers     (ignored; SIFTS-only path)

For each isoform row in the sequence table we strip the "-N" isoform suffix to
the canonical accession, look up its chains in SIFTS, and translate the covered
UniProt region onto every transcript (Protein_ID) sequence of the same gene
using the exact / >=min_identity window search — identical convention to
create_pdb_worker.py — so the per-isoform Protein_ID rows line up.

Schema-compatibility with create_pdb_worker.py
-----------------------------------------------
pdb_structures.tsv columns are IDENTICAL:
    Protein_ID Accession pdb_id chain_id struct_asym_id entity_id
    prot_start prot_end unp_start unp_end resolution experimental_method

SIFTS pdb_chain_uniprot.tsv does NOT carry:
    struct_asym_id  -> emitted blank   (API only)
    entity_id       -> emitted blank   (API only)
    resolution      -> blank, unless --resolu_idx (wwPDB resolu.idx) is supplied
    experimental_method -> blank, unless --resolu_idx supplied (XRAY/NMR/...)

pdb_missing.tsv (unobserved residues = structure-derived disorder) is NOT
reproducible from this flat file — computing it needs the per-PDB
polymer_coverage observed-residue ranges that only the API provides. We still
write a pdb_missing.tsv with the correct header but zero data rows so the
downstream contract (file exists, same columns) is preserved.
"""

import argparse
import gzip
import logging
import sys
import time
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

STRUCT_COLS = ["Protein_ID", "Accession", "pdb_id", "chain_id", "struct_asym_id",
               "entity_id", "prot_start", "prot_end", "unp_start", "unp_end",
               "resolution", "experimental_method"]
MISSING_COLS = ["Protein_ID", "Accession", "pdb_id", "chain_id",
                "prot_start", "prot_end", "unp_start", "unp_end", "length"]


# ---------------------------------------------------------------------------
# UniProt -> transcript region mapping (same convention as create_pdb_worker.py)
# ---------------------------------------------------------------------------

def best_window(region: str, tgt_seq: str, min_identity: float):
    """Return (start_0based, identity) of best ungapped placement, or (-1, id)."""
    L = len(region)
    if L == 0 or len(tgt_seq) < L:
        return -1, 0.0
    idx = tgt_seq.find(region)
    if idx != -1:
        return idx, 1.0
    best_idx, best_id = -1, 0.0
    for i in range(0, len(tgt_seq) - L + 1):
        m = sum(1 for a, b in zip(region, tgt_seq[i:i + L]) if a == b)
        ident = m / L
        if ident > best_id:
            best_id, best_idx = ident, i
    return (best_idx, best_id) if best_id >= min_identity else (-1, best_id)


def canonical_acc(entry_isoform: str) -> str:
    """Strip the -N isoform suffix to the canonical accession for SIFTS lookup."""
    if not entry_isoform:
        return ""
    return entry_isoform.split("-", 1)[0]


def build_gene_lookup(loc_df: pd.DataFrame):
    """
    Returns:
      gene_to_rows: gene -> list of (entry_isoform, protein_id, seq)
      iso_to_seq:   entry_isoform -> sequence (reference for that isoform)
      iso_gene:     entry_isoform -> gene
    """
    gene_col = next((c for c in ["Gene_Gencode", "Gene_Uniprot", "Gene"]
                     if c in loc_df.columns), None)
    gene_to_rows: dict[str, list] = {}
    iso_to_seq: dict[str, str] = {}
    iso_gene: dict[str, str] = {}
    for _, row in loc_df.iterrows():
        iso = str(row.get("Entry_Isoform", "")).strip()
        pid = str(row.get("Protein_ID", "")).strip()
        seq = str(row.get("Sequence", "")).strip()
        gene = str(row.get(gene_col, "")).strip() if gene_col else ""
        if not iso or iso == "nan":
            continue
        seq = "" if seq == "nan" else seq
        if seq:
            iso_to_seq[iso] = seq
            iso_gene[iso] = gene
        gene_to_rows.setdefault(gene, [])
        entry = (iso, pid, seq)
        if entry not in gene_to_rows[gene]:
            gene_to_rows[gene].append(entry)
    return gene_to_rows, iso_to_seq, iso_gene


# ---------------------------------------------------------------------------
# SIFTS load (filtered to the accessions we actually need)
# ---------------------------------------------------------------------------

def _open(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_sifts(path: str, wanted_accs: set) -> dict:
    """
    Stream the SIFTS flat file, keeping only rows whose SP_PRIMARY is wanted.
    Returns acc -> list of {pdb_id, chain_id, unp_start, unp_end}.

    Each SIFTS line is ONE UniProt segment of a (pdb, chain) and is kept as a
    separate row — discontiguous segments of the same chain must NOT be merged
    (e.g. 5hou/A covers 1-61 and 69 with a gap, which the PDBe API also reports
    as two segments). Only exact-duplicate segments are de-duplicated.
    """
    out: dict[str, set] = {}         # acc -> {(pdb, chain, sp_beg, sp_end)}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("PDB\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            pdb, chain, acc = parts[0], parts[1], parts[2]
            if acc not in wanted_accs:
                continue
            try:
                sp_beg = int(parts[7])
                sp_end = int(parts[8])
            except ValueError:
                continue
            if sp_end < sp_beg:
                sp_beg, sp_end = sp_end, sp_beg
            out.setdefault(acc, set()).add((pdb, chain, sp_beg, sp_end))
    # flatten
    flat: dict[str, list] = {}
    for acc, segs in out.items():
        rows = [{"pdb_id": pdb, "chain_id": chain,
                 "unp_start": s, "unp_end": e}
                for (pdb, chain, s, e) in segs]
        flat[acc] = rows
    return flat


def load_resolu_idx(path: str) -> dict:
    """Optional wwPDB resolu.idx -> pdb_id(lower) -> resolution string.
    resolu.idx has no experimental method, so method stays blank for that path.
    -1.00 marks non-X-ray; we keep it blank rather than emit -1."""
    res: dict[str, str] = {}
    try:
        with _open(path) as fh:
            for line in fh:
                if ";" not in line:
                    continue
                parts = [p.strip() for p in line.split(";")]
                if len(parts) < 2:
                    continue
                pid = parts[0].strip().lower()
                val = parts[1].strip()
                if not pid or len(pid) != 4:
                    continue
                try:
                    fv = float(val)
                except ValueError:
                    continue
                res[pid] = "" if fv <= 0 else val
    except OSError as exc:
        log.warning("could not read resolu.idx (%s): %s", path, exc)
    return res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Module 5c (BULK): PDB structure -> transcript mapping from local SIFTS")
    p.add_argument("--seq_table", required=True,
                   help="loc_chrom_with_names_isoforms_with_seq.tsv")
    p.add_argument("--sifts_tsv", required=True,
                   help="pdb_chain_uniprot.tsv(.gz) SIFTS flat file")
    p.add_argument("--outdir", default=".")
    p.add_argument("--resolu_idx", default=None,
                   help="optional wwPDB resolu.idx for the resolution column")
    p.add_argument("--min_identity", type=float, default=0.9)
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    loc_df = pd.read_csv(args.seq_table, sep="\t", dtype=str)
    if "Entry_Isoform" not in loc_df.columns:
        log.error("seq_table missing 'Entry_Isoform' column")
        pd.DataFrame(columns=STRUCT_COLS).to_csv(outdir / "pdb_structures.tsv", sep="\t", index=False)
        pd.DataFrame(columns=MISSING_COLS).to_csv(outdir / "pdb_missing.tsv", sep="\t", index=False)
        return

    gene_to_rows, iso_to_seq, iso_gene = build_gene_lookup(loc_df)

    # canonical accession set to filter SIFTS by
    wanted = {canonical_acc(iso) for iso in iso_to_seq}
    wanted.discard("")
    log.info("Loading SIFTS for %d canonical accessions from %s …",
             len(wanted), args.sifts_tsv)
    sifts = load_sifts(args.sifts_tsv, wanted)
    log.info("SIFTS: %d/%d accessions have PDB chains",
             len(sifts), len(wanted))

    resolu = load_resolu_idx(args.resolu_idx) if args.resolu_idx else {}

    struct_rows = []

    for iso, ref_seq in iso_to_seq.items():
        if not ref_seq:
            continue
        acc = canonical_acc(iso)
        chains = sifts.get(acc)
        if not chains:
            continue
        gene = iso_gene.get(iso, "")
        siblings = gene_to_rows.get(gene, [(iso, "", ref_seq)])

        for seg in chains:
            unp_s, unp_e = seg["unp_start"], seg["unp_end"]
            # clamp to this isoform's reference sequence
            if unp_s < 1 or unp_e < unp_s or unp_e > len(ref_seq):
                continue
            region = ref_seq[unp_s - 1: unp_e]
            pdb_id = seg["pdb_id"]
            resolution = resolu.get(pdb_id, "")

            for (sib_iso, pid, sib_seq) in siblings:
                if not sib_seq or not pid:
                    continue
                idx, ident = best_window(region, sib_seq, args.min_identity)
                if idx == -1:
                    continue
                prot_off = idx + 1 - unp_s          # prot = unp + prot_off
                struct_rows.append({
                    "Protein_ID": pid, "Accession": acc, "pdb_id": pdb_id,
                    "chain_id": seg["chain_id"], "struct_asym_id": "",
                    "entity_id": "",
                    "prot_start": unp_s + prot_off, "prot_end": unp_e + prot_off,
                    "unp_start": unp_s, "unp_end": unp_e,
                    "resolution": resolution, "experimental_method": "",
                })

    # de-duplicate identical rows (a (pdb,chain) might be hit through >1 sibling iso)
    struct_df = pd.DataFrame(struct_rows, columns=STRUCT_COLS).drop_duplicates()
    struct_df.to_csv(outdir / "pdb_structures.tsv", sep="\t", index=False)

    # pdb_missing requires per-PDB observed-residue coverage (API only) -> header only
    pd.DataFrame(columns=MISSING_COLS).to_csv(
        outdir / "pdb_missing.tsv", sep="\t", index=False)

    log.info("PDB(bulk): %d structure-chain rows (%d distinct pdb_id) in %.1fs",
             len(struct_df),
             struct_df["pdb_id"].nunique() if len(struct_df) else 0,
             time.time() - t0)


if __name__ == "__main__":
    main()
