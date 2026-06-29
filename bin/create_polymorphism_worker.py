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
import json
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    stream=sys.stderr)
log = logging.getLogger(__name__)

OUT_COLS_BASE   = ["Protein_ID", "Position", "rsid", "ref", "alt", "allele_frequency", "Type"]
OUT_COLS_GNOMAD = ["Protein_ID", "Position", "rsid", "ref", "alt", "allele_frequency",
                   "gnomad_af", "gnomad_af_popmax", "Type"]
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


def extract_from_dbsnp_maf(maf_gz: str, regions: dict, g2p: dict,
                           protein_ids: set, rows: list) -> int:
    """Extract polymorphisms from the compact dbSNP MAF TSV (from FETCH_DBSNP_VCF).

    Columns: chrom | pos | rsid | ref | alt | maf | is_common
    Loads the compact TSV into memory grouped by chromosome, then for each protein
    queries its positional slice and maps hits to protein residues via g2p.
    """
    p = Path(maf_gz)
    if not p.exists() or p.stat().st_size == 0:
        log.info("dbsnp_maf: missing/empty — skipped")
        return 0

    log.info("dbsnp_maf: loading compact MAF table from %s …", p)
    import pandas as pd
    try:
        df = pd.read_csv(p, sep="\t", dtype=str,
                         names=["chrom", "pos", "rsid", "ref", "alt", "maf", "is_common"],
                         skiprows=1)
    except Exception as exc:
        log.warning("dbsnp_maf: could not read %s — %s", p, exc)
        return 0

    # Group by chromosome for fast per-chrom slice
    chrom_groups: dict = {}
    for chrom, grp in df.groupby("chrom"):
        grp_int = grp.copy()
        try:
            grp_int["_pos_int"] = grp_int["pos"].astype(int)
        except ValueError:
            grp_int["_pos_int"] = 0
        chrom_groups[chrom] = grp_int

    # Same bounding-box strategy as extract_from_bigbed
    chrom_bounds: dict = {}
    chrom_pids:   dict = {}
    for pid in protein_ids:
        reg = regions.get(pid)
        if not reg:
            continue
        chrom, start, end = reg
        if chrom not in chrom_bounds:
            chrom_bounds[chrom] = [start, end]
            chrom_pids[chrom]   = {}
        else:
            chrom_bounds[chrom][0] = min(chrom_bounds[chrom][0], start)
            chrom_bounds[chrom][1] = max(chrom_bounds[chrom][1], end)
        chrom_pids[chrom][pid] = (start, end)

    n = 0
    for chrom, (b_start, b_end) in chrom_bounds.items():
        cdf = chrom_groups.get(chrom)
        if cdf is None or cdf.empty:
            continue
        pids_on_chrom = chrom_pids[chrom]
        window = cdf[(cdf["_pos_int"] >= b_start) & (cdf["_pos_int"] <= b_end)]
        for row in window.itertuples(index=False):
            pos_str = row.pos
            try:
                pos_int = int(pos_str)
            except (ValueError, TypeError):
                continue
            ptype = "Common Polymorphisms" if str(row.is_common) == "1" else "All Polymorphisms"
            for pid, (p_start, p_end) in pids_on_chrom.items():
                if not (p_start <= pos_int <= p_end):
                    continue
                prot_pos = g2p.get((pid, pos_str))
                if prot_pos is None:
                    continue
                rows.append({
                    "Protein_ID": pid,
                    "Position":   prot_pos,
                    "rsid":       row.rsid,
                    "ref":        row.ref,
                    "alt":        row.alt,
                    "allele_frequency": row.maf,
                    "Type":       ptype,
                })
                n += 1
    log.info("dbsnp_maf: %d polymorphism rows extracted (%d isoforms)",
             n, len(protein_ids))
    return n


def extract_from_gnomad_maf(gnomad_gz: str, regions: dict, g2p: dict,
                            protein_ids: set) -> dict:
    """Load the compact gnomAD MAF TSV (from FETCH_GNOMAD_VCF) into a position lookup.

    Returns {(Protein_ID, prot_pos): {'gnomad_af': str, 'gnomad_af_popmax': str}}
    for every position that has gnomAD data.  The caller merges this into rows built
    from the dbSNP extraction step, or creates new rows for gnomAD-only variants.
    """
    p = Path(gnomad_gz)
    if not p.exists() or p.stat().st_size == 0:
        log.info("gnomad_maf: missing/empty — skipped")
        return {}

    log.info("gnomad_maf: loading compact gnomAD MAF table from %s …", p)
    try:
        df = pd.read_csv(p, sep="\t", dtype=str,
                         names=["chrom","pos","rsid","ref","alt","af","af_popmax","is_common"],
                         skiprows=1)
    except Exception as exc:
        log.warning("gnomad_maf: could not read %s — %s", p, exc)
        return {}

    chrom_groups: dict = {}
    for chrom, grp in df.groupby("chrom"):
        grp = grp.copy()
        try:
            grp["_pos_int"] = grp["pos"].astype(int)
        except ValueError:
            grp["_pos_int"] = 0
        chrom_groups[chrom] = grp

    chrom_bounds: dict = {}
    chrom_pids:   dict = {}
    for pid in protein_ids:
        reg = regions.get(pid)
        if not reg:
            continue
        chrom, start, end = reg
        if chrom not in chrom_bounds:
            chrom_bounds[chrom] = [start, end]
            chrom_pids[chrom]   = {}
        else:
            chrom_bounds[chrom][0] = min(chrom_bounds[chrom][0], start)
            chrom_bounds[chrom][1] = max(chrom_bounds[chrom][1], end)
        chrom_pids[chrom][pid] = (start, end)

    result: dict = {}   # (pid, prot_pos) → {gnomad_af, gnomad_af_popmax, rsid, ref, alt}
    n = 0
    for chrom, (b_start, b_end) in chrom_bounds.items():
        cdf = chrom_groups.get(chrom)
        if cdf is None or cdf.empty:
            continue
        pids_on_chrom = chrom_pids[chrom]
        window = cdf[(cdf["_pos_int"] >= b_start) & (cdf["_pos_int"] <= b_end)]
        for row in window.itertuples(index=False):
            pos_str = row.pos
            try:
                pos_int = int(pos_str)
            except (ValueError, TypeError):
                continue
            for pid, (p_start, p_end) in pids_on_chrom.items():
                if not (p_start <= pos_int <= p_end):
                    continue
                prot_pos = g2p.get((pid, pos_str))
                if prot_pos is None:
                    continue
                key = (pid, prot_pos, row.ref, row.alt)
                result[key] = {
                    "gnomad_af":        row.af,
                    "gnomad_af_popmax": row.af_popmax,
                    "rsid":             row.rsid if row.rsid != "." else "",
                }
                n += 1
    log.info("gnomad_maf: %d (protein, position) gnomAD entries mapped", n)
    return result


def _api_get(url: str, retries: int = 3, delay: float = 1.0) -> dict | list | None:
    """HTTP GET with retry and exponential backoff."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 429:              # rate limited
                time.sleep(delay * (2 ** attempt))
            else:
                log.warning("HTTP %d for %s", e.code, url)
                return None
        except Exception as exc:
            log.warning("Request failed (%s): %s", url, exc)
            time.sleep(delay)
    return None


def _api_post(url: str, payload: dict, retries: int = 3, delay: float = 1.0) -> dict | None:
    """HTTP POST JSON with retry."""
    data = json.dumps(payload).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"Content-Type": "application/json",
                                                  "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(delay * (2 ** attempt))
            else:
                log.warning("HTTP %d for %s", e.code, url)
                return None
        except Exception as exc:
            log.warning("Request failed (%s): %s", url, exc)
            time.sleep(delay)
    return None


def extract_from_dbsnp_api(regions: dict, g2p: dict, protein_ids: set,
                            rows: list, max_proteins: int = 200) -> int:
    """Query the Ensembl REST API for dbSNP variants in each protein's genomic region.

    Uses /overlap/region/human/{region}?feature=variation which returns variants from
    dbSNP annotated with rsid and minor allele frequency.  Rate-limited to ~15 req/s
    without an API key — suitable for small runs only (≤ max_proteins proteins).
    """
    base = "https://rest.ensembl.org/overlap/region/human"
    if len(protein_ids) > max_proteins:
        log.warning("dbsnp_api: %d proteins exceeds limit %d — skipped. "
                    "Use --fetch_dbsnp_vcf for large runs.", len(protein_ids), max_proteins)
        return 0

    n = 0
    for pid in sorted(protein_ids):
        reg = regions.get(pid)
        if not reg:
            continue
        chrom, start, end = reg
        # Ensembl uses chromosome numbers without 'chr' prefix
        ens_chrom = chrom.replace("chr", "", 1)
        url = (f"{base}/{ens_chrom}:{start}-{end}"
               f"?feature=variation&content-type=application/json")
        variants = _api_get(url)
        if not variants:
            continue
        for v in variants:
            rsid = v.get("id", "")
            if not rsid or not rsid.startswith("rs"):
                continue
            alleles = v.get("alleles", [])
            if len(alleles) < 2:
                continue
            ref  = alleles[0]
            maf  = str(v.get("minor_allele_freq") or "")
            vtype = "Common Polymorphisms" if maf and float(maf) >= 0.01 else "All Polymorphisms"
            gpos_str = str(v.get("end", ""))
            prot_pos = g2p.get((pid, gpos_str))
            if prot_pos is None:
                continue
            for alt in alleles[1:]:
                if len(ref) == 1 and len(alt) == 1:
                    rows.append({"Protein_ID": pid, "Position": prot_pos,
                                 "rsid": rsid, "ref": ref, "alt": alt,
                                 "allele_frequency": maf, "Type": vtype})
                    n += 1
        time.sleep(0.07)   # ~14 req/s to stay under the 15/s anonymous limit
    log.info("dbsnp_api (Ensembl REST): %d SNPs mapped across %d proteins",
             n, len(protein_ids))
    return n


def extract_from_gnomad_api(regions: dict, g2p: dict, protein_ids: set,
                             max_proteins: int = 200) -> dict:
    """Query the gnomAD GraphQL API for variant allele frequencies.

    Returns the same gnomAD result dict as extract_from_gnomad_maf():
      {(pid, prot_pos, ref, alt): {'gnomad_af': str, 'gnomad_af_popmax': str, 'rsid': str}}
    Suitable for small runs only.
    """
    url = "https://gnomad.broadinstitute.org/api"
    if len(protein_ids) > max_proteins:
        log.warning("gnomad_api: %d proteins exceeds limit %d — skipped. "
                    "Use --fetch_gnomad_vcf for large runs.", len(protein_ids), max_proteins)
        return {}

    result: dict = {}
    n = 0
    for pid in sorted(protein_ids):
        reg = regions.get(pid)
        if not reg:
            continue
        chrom, start, end = reg
        ens_chrom = chrom.replace("chr", "", 1)
        # gnomAD GraphQL API v4: genome/exome each have .af; no af_popmax at variant level
        query = (
            '{ region(chrom: "%s", start: %d, stop: %d, reference_genome: GRCh38) {'
            ' variants(dataset: gnomad_r4) {'
            ' variant_id pos ref alt genome { af } exome { af } } } }'
        ) % (ens_chrom, start, end)
        resp = _api_post(url, {"query": query})
        if not resp:
            continue
        variants = (resp.get("data", {}).get("region", {}) or {}).get("variants") or []
        for v in variants:
            ref = v.get("ref", "")
            alt = v.get("alt", "")
            if len(ref) != 1 or len(alt) != 1:
                continue
            pos_str = str(v.get("pos", ""))
            prot_pos = g2p.get((pid, pos_str))
            if prot_pos is None:
                continue
            genome_af = (v.get("genome") or {}).get("af")
            exome_af  = (v.get("exome")  or {}).get("af")
            # Use whichever is available; prefer the higher of the two
            af_vals = [x for x in [genome_af, exome_af] if x is not None]
            af       = max(af_vals) if af_vals else None
            af_s     = f"{af:.6g}" if af is not None else ""
            # gnomad_af_popmax: not available from region API; store empty string
            key = (pid, prot_pos, ref, alt)
            result[key] = {"gnomad_af": af_s, "gnomad_af_popmax": "", "rsid": ""}
            n += 1
        time.sleep(0.25)   # ~4 req/s to stay within gnomAD GraphQL limits
    log.info("gnomad_api: %d (protein, position) gnomAD entries mapped", n)
    return result


def extract_from_bigbed(dbsnp_bb: str, ucsc_bin: str, regions: dict,
                        g2p: dict, protein_ids: set, rows: list) -> int:
    """Extract every dbSNP record overlapping each selected isoform's genomic
    region directly from the dbSnp*.bb bigBed, map each SNV's genomic coordinate
    to a protein residue via combined_map, and emit a row carrying
    rsid + ref/alt + allele_frequency for ALL selected isoforms.

    The bigBed columns follow the UCSC dbSnp155 schema:
      0 chrom  1 chromStart(0-based)  2 chromEnd(1-based)  3 name(rsid)
      4 ref    6 alts                 9 freqs              14 freqSourceCount/flags

    Optimisation: instead of one bigBedToBed call per isoform region, sweep the
    entire bounding box of each chromosome once (at most 24 calls for a full-
    proteome run) and filter to individual isoform positions in memory.  For a
    single gene this reduces N calls (one per isoform) to 1 call.
    """
    bb = Path(dbsnp_bb)
    if not bb.exists() or bb.stat().st_size == 0:
        log.info("dbsnp_bb: missing/empty — skipped")
        return 0
    tool = Path(ucsc_bin) / "bigBedToBed" if ucsc_bin and ucsc_bin != "NO_FILE" else Path("bigBedToBed")
    tool_s = str(tool) if Path(tool).exists() else "bigBedToBed"
    if not Path(tool_s).exists() and shutil.which(tool_s) is None:
        log.warning("bigBedToBed not found (PATH or --ucsc_bin) — skipping dbSnp "
                    "bigBed extraction; install ucsc-bigbedtobed for the "
                    "polymorphism track")
        return 0

    # Build chromosome → (bounding_start, bounding_end, {pid: (start, end)}) map.
    # One bigBedToBed call per chromosome sweeps the entire locus set at once.
    chrom_bounds: dict = {}   # chrom → [min_start, max_end]
    chrom_pids:   dict = {}   # chrom → {pid: (start, end)}
    for pid in protein_ids:
        reg = regions.get(pid)
        if not reg:
            continue
        chrom, start, end = reg
        if chrom not in chrom_bounds:
            chrom_bounds[chrom] = [start, end]
            chrom_pids[chrom]   = {}
        else:
            chrom_bounds[chrom][0] = min(chrom_bounds[chrom][0], start)
            chrom_bounds[chrom][1] = max(chrom_bounds[chrom][1], end)
        chrom_pids[chrom][pid] = (start, end)

    n = 0
    for chrom, (b_start, b_end) in chrom_bounds.items():
        pids_on_chrom = chrom_pids[chrom]
        with tempfile.NamedTemporaryFile(mode="r", suffix=".bed", delete=False) as tf:
            bed_path = tf.name
        try:
            cmd = [tool_s, dbsnp_bb, f"-chrom={chrom}",
                   f"-start={b_start}", f"-end={b_end}", bed_path]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                log.warning("bigBedToBed failed for %s:%d-%d: %s",
                            chrom, b_start, b_end, res.stderr.strip()[:200])
                continue
            with open(bed_path, encoding="utf-8") as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 7:
                        continue
                    gpos = c[2].strip()             # chromEnd = SNV base (1-based)
                    rsid = c[3]
                    ref  = c[4]
                    alt  = c[6].rstrip(",")
                    freq  = _max_allele_freq(c[9]) if len(c) > 9 else ""
                    ptype = _classify_type(c[14]) if len(c) > 14 else "All Polymorphisms"
                    try:
                        snv_pos = int(gpos)
                    except (ValueError, TypeError):
                        continue
                    for pid, (p_start, p_end) in pids_on_chrom.items():
                        # Quick range check before the dict lookup
                        if not (p_start <= snv_pos <= p_end):
                            continue
                        prot_pos = g2p.get((pid, gpos))
                        if prot_pos is None:
                            try:
                                prot_pos = g2p.get((pid, str(int(c[1]) + 1)))
                            except (ValueError, TypeError):
                                prot_pos = None
                        if prot_pos is None:
                            continue
                        rows.append({
                            "Protein_ID": pid,
                            "Position":   prot_pos,
                            "rsid":       rsid,
                            "ref":        ref,
                            "alt":        alt,
                            "allele_frequency": freq,
                            "Type":       ptype,
                        })
                        n += 1
        finally:
            try:
                Path(bed_path).unlink()
            except OSError:
                pass
    log.info("dbsnp_bb: %d polymorphism rows extracted across %d chromosomes (%d isoforms)",
             n, len(chrom_bounds), len(protein_ids))
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
    ap.add_argument("--dbsnp_maf", default="NO_FILE",
                    help="compact dbSNP MAF gzip TSV from FETCH_DBSNP_VCF "
                         "(preferred over --dbsnp_bb when provided)")
    ap.add_argument("--gnomad_maf", default="NO_FILE",
                    help="compact gnomAD MAF gzip TSV from FETCH_GNOMAD_VCF")
    ap.add_argument("--use_dbsnp_api", action="store_true",
                    help="query Ensembl REST API for dbSNP variants per protein region "
                         "(for small runs only; large runs use --dbsnp_maf)")
    ap.add_argument("--use_gnomad_api", action="store_true",
                    help="query gnomAD GraphQL API for AF per protein region "
                         "(for small runs only; large runs use --gnomad_maf)")
    ap.add_argument("--api_max_proteins", type=int, default=200,
                    help="protein count above which API queries are refused (default 200)")
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
    gnomad_data: dict = {}   # (pid, prot_pos, ref, alt) → {gnomad_af, gnomad_af_popmax, rsid}
    cm = Path(args.combined_map)
    snp_map_used = False
    if cm.exists() and cm.stat().st_size > 0:
        g2p, regions = parse_map(args.combined_map)

        # ── dbSNP source: VCF compact table > Ensembl API > bigBed > legacy .out ──
        if args.dbsnp_maf and args.dbsnp_maf != "NO_FILE" and Path(args.dbsnp_maf).exists():
            n_maf = extract_from_dbsnp_maf(args.dbsnp_maf, regions, g2p, protein_ids, rows)
            snp_map_used = n_maf > 0 or Path(args.dbsnp_maf).stat().st_size > 0
        if not snp_map_used and args.use_dbsnp_api:
            n_api = extract_from_dbsnp_api(regions, g2p, protein_ids, rows,
                                           max_proteins=args.api_max_proteins)
            snp_map_used = n_api > 0
        if not snp_map_used and args.dbsnp_bb and args.dbsnp_bb != "NO_FILE" \
                and Path(args.dbsnp_bb).exists():
            n_bb = extract_from_bigbed(args.dbsnp_bb, args.ucsc_bin, regions,
                                       g2p, protein_ids, rows)
            snp_map_used = n_bb > 0 or Path(args.dbsnp_bb).stat().st_size > 0
        if not snp_map_used:
            parse_out(args.snp_common, "Common Polymorphisms", protein_ids, g2p, rows)
            parse_out(args.snp_all, "All Polymorphisms", protein_ids, g2p, rows)

        # ── gnomAD source (independent of dbSNP, always run if configured) ──
        if args.gnomad_maf and args.gnomad_maf != "NO_FILE" and Path(args.gnomad_maf).exists():
            gnomad_data = extract_from_gnomad_maf(args.gnomad_maf, regions, g2p, protein_ids)
        elif args.use_gnomad_api:
            gnomad_data = extract_from_gnomad_api(regions, g2p, protein_ids,
                                                   max_proteins=args.api_max_proteins)
    else:
        log.info("No combined_map — skipping allele-frequency SNP enrichment")

    covered = {(r["Protein_ID"], str(r["Position"])) for r in rows}

    # 2. Supplement with the comprehensive pre-mapped positional table only when no
    #    VCF/bigBed/API was available (those rows have no allele frequency / rsid).
    if not snp_map_used:
        parse_pos_tsv(args.snp_pos_tsv, protein_ids, covered, rows)

    use_gnomad = bool(gnomad_data)
    out_cols   = OUT_COLS_GNOMAD if use_gnomad else OUT_COLS_BASE

    if use_gnomad:
        # Attach gnomAD AF to existing rows; add gnomAD-only rows for variants not in dbSNP
        new_rows = []
        for r in rows:
            key = (r["Protein_ID"], r["Position"], r.get("ref",""), r.get("alt",""))
            gd  = gnomad_data.pop(key, {})
            r["gnomad_af"]        = gd.get("gnomad_af", "")
            r["gnomad_af_popmax"] = gd.get("gnomad_af_popmax", "")
            if not r.get("rsid") and gd.get("rsid"):
                r["rsid"] = gd["rsid"]
        # Remaining gnomad_data entries are gnomAD-only (not in dbSNP sources)
        for (pid, prot_pos, ref, alt), gd in gnomad_data.items():
            af_f = 0.0
            try:
                af_f = float(gd.get("gnomad_af") or 0)
            except ValueError:
                pass
            vtype = "Common Polymorphisms" if af_f >= 0.01 else "All Polymorphisms"
            new_rows.append({
                "Protein_ID": pid, "Position": prot_pos,
                "rsid": gd.get("rsid", ""), "ref": ref, "alt": alt,
                "allele_frequency": gd.get("gnomad_af", ""),
                "gnomad_af": gd.get("gnomad_af", ""),
                "gnomad_af_popmax": gd.get("gnomad_af_popmax", ""),
                "Type": vtype,
            })
        rows.extend(new_rows)
        log.info("gnomAD-only variants added: %d", len(new_rows))

    df = pd.DataFrame(rows, columns=out_cols).drop_duplicates()
    df.to_csv(out_path, sep="\t", index=False)
    n_freq = (df["allele_frequency"].astype(str) != "").sum() if not df.empty else 0
    n_gaf  = (df["gnomad_af"].astype(str) != "").sum() if use_gnomad and not df.empty else 0
    log.info("polymorphism.tsv: %d rows across %d proteins "
             "(%d with allele_frequency%s) → %s",
             len(df), df["Protein_ID"].nunique() if not df.empty else 0, n_freq,
             f", {n_gaf} with gnomad_af" if use_gnomad else "", out_path)


if __name__ == "__main__":
    main()
