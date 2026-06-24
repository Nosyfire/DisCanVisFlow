#!/usr/bin/env python3
"""
create_polymorphism_worker.py — Module 5g: Polymorphism annotation (legacy SNP).

Source: the legacy SNP pipeline output (common_poly.out / all_poly.out), a
UCSC BED-style table whose last column is the Protein_ID (Gencode transcript)
the SNP was assigned to. Each SNP's genomic coordinate is mapped to a protein
residue position via combined_map.map (per-nucleotide genomic↔protein index),
rather than trusting any precomputed position column.

This replaces the previous UniProt-REST-API approach.

Inputs
------
--loc_chrom      loc_chrom_with_names_isoforms_with_seq.tsv (Protein_ID column)
--combined_map   combined_map.map (Module 3) — genomic↔protein index
--snp_common     common_poly.out  (or NO_FILE)
--snp_all        all_poly.out     (or NO_FILE)
--output_dir     output directory (default: .)

Output
------
polymorphism.tsv — Protein_ID | Position | rsid | ref | alt |
                   allele_frequency | Type
                   allele_frequency = maximum finite per-population minor-allele
                   frequency reported in the .out (blank if undeterminable).
                   Type: 'Common Polymorphisms' | 'All Polymorphisms'
"""

import argparse
import logging
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

OUT_COLS = ["Protein_ID", "Position", "rsid", "ref", "alt", "allele_frequency", "Type"]
_CHUNK = 500_000
# Column index (0-based) of the per-population allele-frequency CSV in the
# legacy UCSC dbSNP .out row (col 10 in 1-based terms).
_FREQ_COL = 9


def _max_allele_freq(cell) -> str:
    """Return the maximum finite per-population frequency as a string, or ''.
    Missing values are encoded as -inf / inf / empty in the legacy .out."""
    if cell is None:
        return ""
    best = None
    for tok in str(cell).split(","):
        tok = tok.strip()
        if not tok or tok in ("-inf", "inf", "nan"):
            continue
        try:
            v = float(tok)
        except ValueError:
            continue
        if v != v or v in (float("inf"), float("-inf")):  # NaN / inf guard
            continue
        if best is None or v > best:
            best = v
    return "" if best is None else f"{best:.6g}"


def parse_map(map_path: str):
    """Parse combined_map.map.

    Returns (g2p, regions):
      g2p     : {(Protein_ID, genomic_pos_str): protein_pos_1based}
      regions : {Protein_ID: (chrom, start_int, end_int)}  from the '# … chrN ± a-b' header
    """
    g2p: dict = {}
    regions: dict = {}
    cur_pid = None
    with open(map_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.lstrip("#").split()
                if len(parts) < 3:
                    cur_pid = None
                    continue
                fields = parts[0].split("|")
                cur_pid = next((f for f in fields if re.match(r".+-\d+$", f)), fields[0])
                # header tail: '<fields> <chrom> <strand> <start>-<end>'
                m = re.search(r"(chr[\w.]+)\s+[+-]\s+(\d+)-(\d+)", line)
                if m:
                    regions[cur_pid] = (m.group(1), int(m.group(2)), int(m.group(3)))
                continue
            if cur_pid is None:
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            try:
                prot_pos = int(cols[0]) + 1                       # 1-based
                gposes = [g for g in cols[5].rstrip(",").split(",") if g and g != "-"]
            except (ValueError, IndexError):
                continue
            for g in gposes:
                g2p[(cur_pid, g)] = prot_pos
    log.info("combined_map: indexed %d (protein,nucleotide) positions, %d regions",
             len(g2p), len(regions))
    return g2p, regions


def _classify_type(freq_class: str) -> str:
    """Map the dbSNP freqClass flag column to a polymorphism Type label.
    'commonAll'/'commonSome' → Common Polymorphisms, otherwise All Polymorphisms."""
    fc = (freq_class or "").lower()
    if "common" in fc:
        return "Common Polymorphisms"
    return "All Polymorphisms"


def extract_from_bigbed(dbsnp_bb: str, ucsc_bin: str, regions: dict,
                        g2p: dict, protein_ids: set, rows: list) -> int:
    """Extract every dbSNP record overlapping each selected isoform's genomic
    region directly from the dbSnp*.bb bigBed (legacy get_snp.py approach), map
    each SNV's genomic coordinate to a protein residue via combined_map, and emit
    a row carrying rsid + ref/alt + allele_frequency for ALL selected isoforms.

    The bigBed columns follow the UCSC dbSnp155 schema (same as common_poly.out):
      0 chrom  1 chromStart(0-based)  2 chromEnd(1-based)  3 name(rsid)
      4 ref    6 alts                 9 freqs              14 freqSourceCount/flags
    """
    bb = Path(dbsnp_bb)
    if not bb.exists() or bb.stat().st_size == 0:
        log.info("dbsnp_bb: missing/empty — skipped")
        return 0
    tool = Path(ucsc_bin) / "bigBedToBed" if ucsc_bin and ucsc_bin != "NO_FILE" else Path("bigBedToBed")
    tool_s = str(tool) if Path(tool).exists() else "bigBedToBed"
    # Degrade gracefully if the UCSC binary is absent: a missing optional tool
    # should not abort the whole pipeline — emit no bigBed rows and warn.
    if not Path(tool_s).exists() and shutil.which(tool_s) is None:
        log.warning("bigBedToBed not found (PATH or --ucsc_bin) — skipping dbSnp "
                    "bigBed extraction; install ucsc-bigbedtobed for the "
                    "polymorphism track")
        return 0

    # Dedupe identical genomic regions so we only run bigBedToBed once per locus,
    # then fan the SNVs out to every isoform sharing that region.
    region_to_pids: dict = {}
    for pid in protein_ids:
        reg = regions.get(pid)
        if reg:
            region_to_pids.setdefault(reg, []).append(pid)

    n = 0
    for (chrom, start, end), pids in region_to_pids.items():
        with tempfile.NamedTemporaryFile(mode="r", suffix=".bed", delete=False) as tf:
            bed_path = tf.name
        try:
            cmd = [tool_s, dbsnp_bb, f"-chrom={chrom}",
                   f"-start={start}", f"-end={end}", bed_path]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                log.warning("bigBedToBed failed for %s:%d-%d: %s",
                            chrom, start, end, res.stderr.strip()[:200])
                continue
            with open(bed_path, encoding="utf-8") as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 7:
                        continue
                    gpos = c[2].strip()             # chromEnd = SNV base (1-based)
                    rsid = c[3]
                    ref = c[4]
                    alt = c[6].rstrip(",")
                    freq = _max_allele_freq(c[9]) if len(c) > 9 else ""
                    ptype = _classify_type(c[14]) if len(c) > 14 else "All Polymorphisms"
                    for pid in pids:
                        prot_pos = g2p.get((pid, gpos))
                        if prot_pos is None:
                            # fallback: 0-based start + 1
                            try:
                                prot_pos = g2p.get((pid, str(int(c[1]) + 1)))
                            except (ValueError, TypeError):
                                prot_pos = None
                        if prot_pos is None:
                            continue   # SNV not in this isoform's CDS (e.g. UTR/intron)
                        rows.append({
                            "Protein_ID": pid,
                            "Position": prot_pos,
                            "rsid": rsid,
                            "ref": ref,
                            "alt": alt,
                            "allele_frequency": freq,
                            "Type": ptype,
                        })
                        n += 1
        finally:
            try:
                Path(bed_path).unlink()
            except OSError:
                pass
    log.info("dbsnp_bb: %d polymorphism rows extracted across %d regions",
             n, len(region_to_pids))
    return n


def parse_out(path: str, type_label: str, protein_ids: set, g2p: dict, rows: list) -> int:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        log.info("%s: missing/empty — skipped", type_label)
        return 0
    n = 0
    for chunk in pd.read_csv(p, sep="\t", header=None, dtype=str, chunksize=_CHUNK):
        sub = chunk[chunk.iloc[:, -1].isin(protein_ids)]
        for _, r in sub.iterrows():
            pid = r.iloc[-1]
            # UCSC BED: col1 = 0-based start, col2 = 1-based end → SNV base = end
            gpos = str(r.iloc[2]).strip()
            prot_pos = g2p.get((pid, gpos))
            if prot_pos is None:
                # fallback: start+1
                try:
                    prot_pos = g2p.get((pid, str(int(float(r.iloc[1])) + 1)))
                except (ValueError, TypeError):
                    prot_pos = None
            if prot_pos is None:
                continue   # SNP not in this protein's CDS (e.g. UTR) → skip
            freq = _max_allele_freq(r.iloc[_FREQ_COL]) if len(r) > _FREQ_COL else ""
            rows.append({
                "Protein_ID": pid,
                "Position": prot_pos,
                "rsid": str(r.iloc[3]),
                "ref": str(r.iloc[4]),
                "alt": str(r.iloc[6]).rstrip(","),
                "allele_frequency": freq,
                "Type": type_label,
            })
            n += 1
    log.info("%s: %d SNPs mapped to run proteins", type_label, n)
    return n


def parse_pos_tsv(path: str, protein_ids: set, covered: set, rows: list) -> int:
    """Fold the comprehensive pre-mapped polymorphism table (polymorphism_pos.tsv:
    'AccessionPosition' = Protein_ID|position, 'Polymorphism' = type) into the
    output. These are the *all* polymorphisms; positions already enriched from the
    SNP .out (allele frequency / rsid) are skipped here to avoid duplication."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        log.info("snp_pos: missing/empty — skipped")
        return 0
    n = 0
    for chunk in pd.read_csv(p, sep="\t", dtype=str, chunksize=_CHUNK):
        acc_col = "AccessionPosition" if "AccessionPosition" in chunk.columns else chunk.columns[0]
        type_col = "Polymorphism" if "Polymorphism" in chunk.columns else chunk.columns[-1]
        for _, r in chunk.iterrows():
            key = str(r[acc_col])
            if "|" not in key:
                continue
            pid, _, pos = key.rpartition("|")
            if pid not in protein_ids:
                continue
            if (pid, str(pos)) in covered:
                continue   # already have an allele-frequency-bearing row here
            rows.append({
                "Protein_ID": pid,
                "Position": pos,
                "rsid": "",
                "ref": "",
                "alt": "",
                "allele_frequency": "",
                "Type": str(r[type_col]),
            })
            n += 1
    log.info("snp_pos: %d additional polymorphism positions for run proteins", n)
    return n


def main():
    ap = argparse.ArgumentParser(description="Module 5g: polymorphism mapping "
                                             "(all polymorphisms + allele frequency)")
    ap.add_argument("--loc_chrom", required=True)
    ap.add_argument("--combined_map", default="NO_FILE")
    ap.add_argument("--snp_common", default="NO_FILE")
    ap.add_argument("--snp_all", default="NO_FILE")
    ap.add_argument("--snp_pos_tsv", default="NO_FILE",
                    help="polymorphism_pos.tsv — comprehensive pre-mapped "
                         "Protein_ID|position polymorphism table (fallback only)")
    ap.add_argument("--dbsnp_bb", default="NO_FILE",
                    help="dbSnp*.bb bigBed (e.g. dbSnp155Common.bb) — extract "
                         "polymorphisms with rsid + allele frequency per isoform")
    ap.add_argument("--ucsc_bin", default="NO_FILE",
                    help="directory containing bigBedToBed (defaults to PATH)")
    ap.add_argument("--output_dir", default=".")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "polymorphism.tsv"

    protein_ids = set(
        pd.read_csv(args.loc_chrom, sep="\t", dtype=str, usecols=["Protein_ID"])
        ["Protein_ID"].dropna()
    )
    if not protein_ids:
        pd.DataFrame(columns=OUT_COLS).to_csv(out_path, sep="\t", index=False)
        log.warning("No run proteins — wrote empty %s", out_path)
        return

    rows: list = []
    cm = Path(args.combined_map)
    bigbed_used = False
    if cm.exists() and cm.stat().st_size > 0:
        g2p, regions = parse_map(args.combined_map)
        # 1a. Preferred: extract directly from the dbSnp bigBed per isoform so that
        #     EVERY selected isoform gets its polymorphisms with rsid + allele freq.
        if args.dbsnp_bb and args.dbsnp_bb != "NO_FILE" and Path(args.dbsnp_bb).exists():
            n_bb = extract_from_bigbed(args.dbsnp_bb, args.ucsc_bin, regions,
                                       g2p, protein_ids, rows)
            bigbed_used = n_bb > 0 or Path(args.dbsnp_bb).stat().st_size > 0
        # 1b. Fallback / supplement: pre-baked legacy .out tables (if present)
        if not bigbed_used:
            parse_out(args.snp_common, "Common Polymorphisms", protein_ids, g2p, rows)
            parse_out(args.snp_all, "All Polymorphisms", protein_ids, g2p, rows)
    else:
        log.info("No combined_map — skipping allele-frequency SNP enrichment")

    covered = {(r["Protein_ID"], str(r["Position"])) for r in rows}

    # 2. Supplement with the comprehensive pre-mapped positional table only when no
    #    bigBed was available (those rows have no allele frequency / rsid).
    if not bigbed_used:
        parse_pos_tsv(args.snp_pos_tsv, protein_ids, covered, rows)

    df = pd.DataFrame(rows, columns=OUT_COLS).drop_duplicates()
    df.to_csv(out_path, sep="\t", index=False)
    n_freq = (df["allele_frequency"].astype(str) != "").sum() if not df.empty else 0
    log.info("polymorphism.tsv: %d rows across %d proteins (%d with allele frequency) → %s",
             len(df), df["Protein_ID"].nunique() if not df.empty else 0, n_freq, out_path)


if __name__ == "__main__":
    main()
