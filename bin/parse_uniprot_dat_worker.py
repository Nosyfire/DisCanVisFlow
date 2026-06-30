#!/usr/bin/env python3
"""
parse_uniprot_dat_worker.py — Bulk-extract UniProt feature annotations and Pfam
domain mappings from two FTP flat files:

  uniprot_sprot.dat.gz  — Swiss-Prot flat file, FTP from UniProt
  protein2ipr.dat.gz    — InterPro protein→domain mappings, FTP from EBI

This replaces ~37 k per-protein REST API calls in create_annotation_worker.py
(UniProt REST + InterPro REST) with two local file scans.  Outputs are
storeDir-cached by Nextflow so parsing only runs once per dat.gz version.

Outputs
-------
  uniprot_features.tsv  cols: Accession, Type, Start, End, Note, Evidence, Ligand
  pfam_domains.tsv      cols: Accession, hmm_acc, hmm_name, start, end, type

Usage
-----
  parse_uniprot_dat_worker.py \\
      --uniprot_dat   uniprot_sprot.dat.gz \\
      --interpro_pfam protein2ipr.dat.gz   \\
      --accessions    accessions.txt        \\
      --outdir        .
"""

import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

# ── Feature type mappings: UniProt dat FT keyword → display name ──────────────
# Only these are kept in uniprot_features.tsv (ROI + binding).
_FT_TYPE_MAP: dict[str, str] = {
    "SIGNAL":   "Signal peptide",
    "TRANSMEM": "Transmembrane",
    "TOPO_DOM": "Topological domain",
    "INTRAMEM": "Intramembrane",
    "PROPEP":   "Propeptide",
    "CHAIN":    "Chain",
    "COILED":   "Coiled coil",
    "COMPBIAS": "Compositional bias",
    "REGION":   "Region",
    "SITE":     "Site",
    "TRANSIT":  "Transit peptide",
    "BINDING":  "Binding site",
    "ACT_SITE": "Active site",
}

# FT feature header pattern: "FT   TYPE   start..end" or "FT   TYPE   start"
_FEAT_RE = re.compile(r'^FT   (\S+)\s+([<>?\d]+)(?:\.\.([<>?\d]+))?')
_NOTE_RE = re.compile(r'^FT\s+/note="([^"]*)"')
_LIGAND_RE = re.compile(r'^FT\s+/ligand="([^"]*)"')
_EVID_RE = re.compile(r'^FT\s+/evidence="([^"]*)"')

_OUT_FEAT_COLS = ["Accession", "Type", "Start", "End", "Note", "Evidence", "Ligand"]
_OUT_PFAM_COLS = ["Accession", "hmm_acc", "hmm_name", "start", "end", "type"]


# ── UniProt flat file parser ──────────────────────────────────────────────────

def _strip_pos(s: str) -> str | None:
    """Strip '<', '>', '?' position modifiers; return None for uncertain positions."""
    s = s.strip().lstrip('<>')
    return None if '?' in s else s


def parse_uniprot_dat(dat_gz: Path, keep: set[str] | None) -> tuple[list, list]:
    """
    Stream uniprot_sprot.dat.gz and collect:
      • FT feature lines → feat_rows  (for ROI + binding types)
      • DR Pfam lines    → pfam_dr    (accession, pfam_acc, pfam_name; NO positions)

    Returns (feat_rows, pfam_dr) both as lists of dicts.
    keep: set of accessions to keep; None = keep all.
    """
    feat_rows: list[dict] = []
    pfam_dr: list[dict] = []

    acc: str | None = None
    cur_feat: dict | None = None
    skip_entry: bool = False

    total = 0

    open_fn = gzip.open if str(dat_gz).endswith('.gz') else open

    with open_fn(dat_gz, 'rt', encoding='latin-1') as fh:
        for raw in fh:
            tag = raw[:2]

            if tag == 'AC':
                if acc is None:
                    # First AC line of this entry → primary accession
                    acc = raw[5:].strip().rstrip(';').split(';')[0].strip()
                    skip_entry = (keep is not None and acc not in keep)
                    total += 1
                    if total % 50_000 == 0:
                        log.info("  dat: %d entries read, %d features kept …",
                                 total, len(feat_rows))

            elif tag == '//' :
                # End of entry
                acc = None
                cur_feat = None
                skip_entry = False

            elif skip_entry:
                continue

            elif tag == 'FT':
                m = _FEAT_RE.match(raw)
                if m:
                    # Feature header line
                    raw_type = m.group(1)
                    display  = _FT_TYPE_MAP.get(raw_type)
                    start_s  = _strip_pos(m.group(2) or '')
                    end_s    = _strip_pos(m.group(3) or m.group(2) or '')
                    if display and start_s and end_s:
                        cur_feat = {
                            'Accession': acc,
                            'Type':      display,
                            'Start':     int(start_s),
                            'End':       int(end_s),
                            'Note':      '',
                            'Evidence':  '',
                            'Ligand':    '',
                        }
                        feat_rows.append(cur_feat)
                    else:
                        cur_feat = None
                elif cur_feat:
                    # Qualifier continuation for current feature
                    nm = _NOTE_RE.match(raw)
                    if nm:
                        cur_feat['Note'] = nm.group(1)
                    else:
                        lm = _LIGAND_RE.match(raw)
                        if lm:
                            cur_feat['Ligand'] = lm.group(1)
                        else:
                            em = _EVID_RE.match(raw)
                            if em:
                                cur_feat['Evidence'] = em.group(1)

            elif tag == 'DR':
                # DR   Pfam; PF00001; 7tm_1; 3.
                if 'Pfam;' in raw:
                    parts = [p.strip() for p in raw[5:].strip().rstrip(';').split(';')]
                    if len(parts) >= 3 and parts[0] == 'Pfam':
                        pfam_dr.append({
                            'acc':       acc,
                            'pfam_acc':  parts[1],
                            'pfam_name': parts[2],
                        })

    log.info("dat: %d entries scanned → %d feature rows, %d Pfam DR refs",
             total, len(feat_rows), len(pfam_dr))
    return feat_rows, pfam_dr


# ── InterPro protein2ipr.dat.gz parser ───────────────────────────────────────

def parse_interpro_pfam(interpro_gz: Path, keep: set[str] | None) -> list[dict]:
    """
    Stream protein2ipr.dat.gz (tab-separated, all species).

    Actual EBI FTP format (protein2ipr.dat, 6 columns):
      0  UniProt_accession      e.g. P04049
      1  InterPro_accession     e.g. IPR004839
      2  InterPro_description   e.g. "Aminotransferase, class I..."
      3  signature_accession    e.g. PF00155 (Pfam), TIGR01821 (TIGRFAM), ...
      4  start_location         1-based integer
      5  end_location           1-based integer

    Pfam entries are identified by signature_accession starting with 'PF'.

    Returns list of {Accession, hmm_acc, hmm_name, start, end, type='Pfam'}.
    """
    rows: list[dict] = []
    total = 0

    open_fn = gzip.open if str(interpro_gz).endswith('.gz') else open

    with open_fn(interpro_gz, 'rt', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            total += 1
            if total % 5_000_000 == 0:
                log.info("  interpro: %dM lines scanned, %d Pfam rows kept …",
                         total // 1_000_000, len(rows))
            if '\t' not in raw:
                continue
            parts = raw.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            acc     = parts[0]
            sig_acc = parts[3]
            if not sig_acc.startswith('PF'):
                continue
            if keep and acc not in keep:
                continue
            try:
                start = int(parts[4])
                end   = int(parts[5])
            except ValueError:
                continue
            rows.append({
                'Accession': acc,
                'hmm_acc':   sig_acc,
                'hmm_name':  parts[2],
                'start':     start,
                'end':       end,
                'type':      'Pfam',
            })

    log.info("interpro: %d lines scanned → %d Pfam domain rows", total, len(rows))
    return rows


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--uniprot_dat',   required=True, help='uniprot_sprot.dat.gz')
    p.add_argument('--interpro_pfam', default=None,  help='protein2ipr.dat.gz (optional)')
    p.add_argument('--accessions',    default=None,
                   help='text file with one accession per line to filter (optional, default=all)')
    p.add_argument('--outdir',        required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Optional accession filter
    keep: set[str] | None = None
    if args.accessions and Path(args.accessions).exists():
        keep = {l.strip() for l in Path(args.accessions).read_text().splitlines() if l.strip()}
        log.info("Filtering to %d accessions from %s", len(keep), args.accessions)
    else:
        log.info("No accession filter — keeping all Swiss-Prot entries")

    # ── Parse UniProt dat ────────────────────────────────────────────────────
    log.info("Parsing %s …", args.uniprot_dat)
    feat_rows, pfam_dr = parse_uniprot_dat(Path(args.uniprot_dat), keep)

    feat_df = pd.DataFrame(feat_rows, columns=_OUT_FEAT_COLS) if feat_rows \
              else pd.DataFrame(columns=_OUT_FEAT_COLS)
    feat_df.to_csv(outdir / 'uniprot_features.tsv', sep='\t', index=False)
    log.info("uniprot_features.tsv: %d rows (%d unique accessions)",
             len(feat_df), feat_df['Accession'].nunique() if not feat_df.empty else 0)

    # ── Parse InterPro protein2ipr (preferred Pfam source with positions) ───
    pfam_rows: list[dict] = []
    if args.interpro_pfam and Path(args.interpro_pfam).exists() \
            and Path(args.interpro_pfam).name != 'NO_FILE':
        log.info("Parsing %s for Pfam domains …", args.interpro_pfam)
        pfam_rows = parse_interpro_pfam(Path(args.interpro_pfam), keep)
    else:
        # Fallback: use DR lines from dat file (no positions, but at least Pfam accessions).
        # This mode produces pfam_domains.tsv without start/end — ANNOTATION_MAP
        # workers that need positions should supply --interpro_pfam.
        log.info("No protein2ipr.dat.gz supplied; Pfam table will be empty. "
                 "Supply --interpro_pfam for domain positions.")

    pfam_df = pd.DataFrame(pfam_rows, columns=_OUT_PFAM_COLS) if pfam_rows \
              else pd.DataFrame(columns=_OUT_PFAM_COLS)
    pfam_df.to_csv(outdir / 'pfam_domains.tsv', sep='\t', index=False)
    log.info("pfam_domains.tsv: %d rows (%d unique accessions)",
             len(pfam_df), pfam_df['Accession'].nunique() if not pfam_df.empty else 0)


if __name__ == '__main__':
    main()
