#!/usr/bin/env python3
"""
create_annotation_worker.py — Module 5: Annotation Mapping

Maps functional annotations to proteins from loc_chrom_with_names_isoforms_with_seq.tsv.

Data sources (all from legacy_data/ local files or cached downloads)
----------------------------------------------------------------------
  ELM        — elm_instances-2023.tsv (legacy, skiprows=5)
  DIBS       — dibs_parsed.tsv (space-sep, no header: acc id start end)
  MFIB       — mfib_parsed.tsv (space-sep, no header: acc id start end)
  PhasePro   — phasepro_parsed.tsv (space-sep, no header: acc pubmed start end)
  UniProt    — ROI + binding site features via REST API
  PTMdb      — local TSV: Gene_name Entry_Isoform Position Type PubmedID Motifs (tab, no header)
  PhosphoSite— local TSV: skip 3 rows, filter ORGANISM==human, col ACC_ID → Entry_Isoform
  Pfam       — InterPro REST API per accession

Isoform annotation transfer
----------------------------
  For canonical annotations, check if the annotated region is 100 % identical
  in each alternative GENCODE isoform via substring search.
  Transferred rows get homology_transfer=True.

Inputs
------
  --loc_chrom     loc_chrom_with_names_isoforms_with_seq.tsv  (Module 2)
  --elm_tsv       ELM instances TSV (legacy_data/elm/elm_instances-2023.tsv)
  --dibs_tsv      DIBS parsed TSV  (legacy_data/dibs/dibs_parsed.tsv)
  --mfib_tsv      MFIB parsed TSV  (legacy_data/mfib/mfib_parsed.tsv)
  --phasepro_tsv  PhasePro parsed TSV (legacy_data/phasepro/phasepro_parsed.tsv)
  --ptmdb_dir     dir with PTMdb type files  (legacy_data/ptm/ptmdb/)
  --ptmphs_dir    dir with PhosphoSite files (legacy_data/ptm/ptmphs/)
  --output_dir    output directory (default: .)
  --request_delay seconds between API requests (default: 0.5)
  --skip_uniprot  skip UniProt REST API
  --skip_pfam     skip InterPro/Pfam API

Outputs
-------
  elm.tsv
  dibs.tsv
  mfib.tsv
  phasepro.tsv
  uniprot_roi.tsv
  uniprot_binding.tsv
  ptm_merged.tsv
  pfam_domains.tsv
  annotation_stats.tsv
"""

import argparse
import io
import logging
import sys
import time
from pathlib import Path

import pandas as pd

try:
    import requests
    _HAS_REQUESTS = True
except ImportError:
    _HAS_REQUESTS = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/{acc}.json"
INTERPRO_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/?format=json"

ROI_TYPES     = {"Region", "Site", "Propeptide", "Signal peptide",
                 "Transit peptide", "Chain", "Coiled coil", "Compositional bias",
                 "Intramembrane", "Topological domain", "Transmembrane"}
BINDING_TYPES = {"Binding site", "Active site"}
PTM_TYPES     = {"Modified residue", "Cross-link", "Glycosylation",
                 "Lipidation", "Disulfide bond", "Natural variant"}

# ---------------------------------------------------------------------------
# HTTP helper
# ---------------------------------------------------------------------------

def _get(url: str, delay: float = 0.5, timeout: int = 60):
    if not _HAS_REQUESTS:
        return None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
            log.warning("HTTP %d  %s (attempt %d)", r.status_code, url, attempt + 1)
        except requests.RequestException as e:
            log.warning("Request error: %s (attempt %d)", e, attempt + 1)
        time.sleep(delay * (attempt + 1))
    return None


# ---------------------------------------------------------------------------
# loc_chrom loader
# ---------------------------------------------------------------------------

def load_loc_chrom(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["Entry_Isoform"].notna() & (df["Entry_Isoform"] != "")]
    log.info("loc_chrom: %d rows, %d unique accessions",
             len(df), df["Entry_Isoform"].nunique())
    return df


def _main_isoform_df(df: pd.DataFrame) -> pd.DataFrame:
    if "main_isoform" in df.columns:
        return df[df["main_isoform"].str.lower() == "yes"]
    return df


# ---------------------------------------------------------------------------
# Isoform annotation transfer — 100% identity substring search
# ---------------------------------------------------------------------------

def _apply_isoform_transfer(
    annot_df:   pd.DataFrame,
    loc_df:     pd.DataFrame,
    out_path:   Path,
    out_cols:   list[str],
) -> pd.DataFrame:
    """
    Append rows for alternative isoforms where the annotated region is
    100 % identical (substring match).  Marks homology_transfer=True.
    Returns the combined DataFrame and re-saves the TSV.
    """
    if annot_df.empty or "Sequence" not in loc_df.columns:
        return annot_df

    acc_col   = next((c for c in ["Entry_Isoform", "Accession"] if c in annot_df.columns), None)
    gene_col  = next((c for c in ["Gene", "Gene_Gencode", "Gene_Uniprot"] if c in loc_df.columns), None)
    if acc_col is None:
        return annot_df

    # Build acc→sequence dict; force all values to str so NaN becomes "nan"
    acc_to_seq: dict[str, str] = {
        k: str(v) for k, v in loc_df.set_index("Entry_Isoform")["Sequence"].to_dict().items()
    }
    gene_groups: dict[str, list[str]] = {}
    if gene_col:
        for _, row in loc_df.iterrows():
            gene_groups.setdefault(str(row.get(gene_col, "")), []).append(
                str(row.get("Entry_Isoform", "")))

    extra_rows = []
    for _, row in annot_df.iterrows():
        can_acc = str(row.get(acc_col, ""))
        can_seq = acc_to_seq.get(can_acc, "")
        if not can_seq or can_seq in ("nan", ""):
            continue
        try:
            s = int(float(row.get("Start", 0)))
            e = int(float(row.get("End",   0)))
        except (ValueError, TypeError):
            continue
        if s < 1 or e < s:
            continue
        region = can_seq[s - 1: e]
        if not region:
            continue

        if gene_col:
            gene_vals = loc_df.loc[loc_df["Entry_Isoform"] == can_acc, gene_col]
            gene_name = gene_vals.iloc[0] if not gene_vals.empty else ""
            alt_accs  = [a for a in gene_groups.get(str(gene_name), []) if a != can_acc]
        else:
            alt_accs = []

        for alt_acc in alt_accs:
            alt_seq = acc_to_seq.get(alt_acc, "")
            if not alt_seq or alt_seq in ("nan", ""):
                continue
            idx = alt_seq.find(region)
            if idx != -1:
                nr = row.to_dict()
                nr[acc_col] = alt_acc
                nr["Start"]             = idx + 1
                nr["End"]               = idx + len(region)
                nr["homology_transfer"] = True
                extra_rows.append(nr)

    if extra_rows:
        result = pd.concat([annot_df, pd.DataFrame(extra_rows)], ignore_index=True)
        real_cols = [c for c in out_cols if c in result.columns]
        result[real_cols].to_csv(out_path, sep="\t", index=False)
        log.info("Isoform transfer: +%d rows → %s", len(extra_rows), out_path.name)
        return result

    return annot_df


# ---------------------------------------------------------------------------
# ELM — from legacy TSV (skip # header lines, skiprows=5)
# ---------------------------------------------------------------------------

def map_elm(loc_df: pd.DataFrame, elm_tsv: str | None,
            outdir: Path, delay: float) -> int:
    out_cols = ["Protein_ID", "Entry_Isoform", "Accession", "ELMType",
                "ELMIdentifier", "Start", "End", "References",
                "Methods", "InstanceLogic", "Organism", "homology_transfer"]

    if not elm_tsv or not Path(elm_tsv).exists():
        log.warning("ELM TSV not provided or missing; skipping ELM.")
        pd.DataFrame(columns=out_cols).to_csv(outdir / "elm.tsv", sep="\t", index=False)
        return 0

    try:
        elm = pd.read_csv(elm_tsv, sep="\t", skiprows=5, dtype=str)
        # Strip quotes from column names if present
        elm.columns = elm.columns.str.strip('"').str.strip()
    except Exception as e:
        log.warning("ELM parse error: %s", e)
        pd.DataFrame(columns=out_cols).to_csv(outdir / "elm.tsv", sep="\t", index=False)
        return 0

    elm = elm[elm["Organism"].str.contains("Homo sapiens", na=False)]
    # Strip quotes from all string values
    for col in elm.columns:
        elm[col] = elm[col].astype(str).str.strip('"')

    # Explode multi-accession rows
    acc_col = "Primary_Acc" if "Primary_Acc" in elm.columns else "Accessions"
    if elm[acc_col].str.contains(r"\s", na=False).any():
        elm = elm.assign(**{acc_col: elm[acc_col].str.split()}).explode(acc_col)

    elm = elm.rename(columns={acc_col: "Entry_Isoform", "Accession": "ELM_Accession"})

    main_df = _main_isoform_df(loc_df)
    meta_cols = [c for c in ["Entry_Isoform", "Entry_Name", "Chromosome", "Gene",
                              "Gene_Gencode", "Gene_Uniprot", "Sequence"] if c in loc_df.columns]
    merged = pd.merge(elm, main_df[meta_cols].drop_duplicates(),
                      on="Entry_Isoform", how="inner")
    if "Entry_Name" in merged.columns:
        merged = merged.rename(columns={"Entry_Name": "Protein_ID"})
    merged["homology_transfer"] = False

    out_path = outdir / "elm.tsv"
    real_cols = [c for c in out_cols if c in merged.columns]
    merged[real_cols].to_csv(out_path, sep="\t", index=False)
    merged = _apply_isoform_transfer(merged, loc_df, out_path, real_cols)
    log.info("ELM: %d rows → elm.tsv", len(merged))
    return len(merged)


# ---------------------------------------------------------------------------
# DIBS / MFIB / PhasePro — space-sep legacy files, no header
# Format: acc  id/pubmed  start  end
# ---------------------------------------------------------------------------

def _load_binding_legacy(tsv_path: str | None, source_label: str,
                          cols: list[str]) -> pd.DataFrame | None:
    if not tsv_path or not Path(tsv_path).exists():
        log.info("%s: file not provided/missing, skipping.", source_label)
        return None
    try:
        df = pd.read_csv(tsv_path, sep=r"\s+", header=None, dtype=str,
                         names=cols, engine="python")
        df["Data"] = source_label
        return df
    except Exception as e:
        log.warning("%s read error: %s", source_label, e)
        return None


def _map_binding_single(df: pd.DataFrame | None, loc_df: pd.DataFrame,
                         source: str, outdir: Path) -> int:
    out_cols = ["Accession", "Entry_Isoform", "Name", "Start", "End",
                "Data", "homology_transfer"]
    out_path = outdir / f"{source}.tsv"
    if df is None or df.empty:
        pd.DataFrame(columns=out_cols).to_csv(out_path, sep="\t", index=False)
        return 0

    # Rename acc column
    if "Entry_Isoform" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Entry_Isoform"})
    if "Name" not in df.columns and len(df.columns) > 1:
        df = df.rename(columns={df.columns[1]: "Name"})

    main_df   = _main_isoform_df(loc_df)
    meta_cols = [c for c in ["Entry_Isoform", "Entry_Name", "Chromosome",
                              "Gene", "Gene_Gencode", "Sequence"] if c in loc_df.columns]
    merged = pd.merge(df, main_df[meta_cols].drop_duplicates(),
                      on="Entry_Isoform", how="inner")
    if "Entry_Name" in merged.columns:
        merged = merged.rename(columns={"Entry_Name": "Accession"})
    if "Accession" not in merged.columns:
        merged["Accession"] = merged["Entry_Isoform"]
    merged["homology_transfer"] = False

    real_cols = [c for c in out_cols if c in merged.columns]
    merged[real_cols].to_csv(out_path, sep="\t", index=False)
    merged = _apply_isoform_transfer(merged, loc_df, out_path, real_cols)
    log.info("%s: %d rows → %s", source, len(merged), out_path.name)
    return len(merged)


def map_dibs(loc_df, dibs_tsv, outdir):
    df = _load_binding_legacy(dibs_tsv, "dibs",
                               ["Entry_Isoform", "Name", "Start", "End"])
    return _map_binding_single(df, loc_df, "dibs", outdir)


def map_mfib(loc_df, mfib_tsv, outdir):
    df = _load_binding_legacy(mfib_tsv, "mfib",
                               ["Entry_Isoform", "Name", "Start", "End"])
    return _map_binding_single(df, loc_df, "mfib", outdir)


def map_phasepro(loc_df, phasepro_tsv, outdir):
    df = _load_binding_legacy(phasepro_tsv, "phasepro",
                               ["Entry_Isoform", "Name", "Start", "End"])
    return _map_binding_single(df, loc_df, "phasepro", outdir)


# ---------------------------------------------------------------------------
# PTM — PTMdb (local tab-sep) + PhosphoSite (local TSV)
# ---------------------------------------------------------------------------

def _ptmdb_type(path: Path, ptype: str, loc_df: pd.DataFrame) -> pd.DataFrame:
    """Load one PTMdb type file and merge with loc_chrom."""
    try:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str,
                         names=["Gene_name", "Entry_Isoform", "Position",
                                "Type", "PubmedID", "Motifs"])
    except Exception as e:
        log.warning("PTMdb %s read error: %s", ptype, e)
        return pd.DataFrame()

    main_df   = _main_isoform_df(loc_df)
    meta_cols = [c for c in ["Entry_Isoform", "Sequence"] if c in loc_df.columns]
    merged = pd.merge(df, main_df[meta_cols].drop_duplicates(),
                      on="Entry_Isoform", how="inner")

    # Quality check: motif in sequence around position
    def qcheck(row):
        try:
            pos   = int(row["Position"])
            motif = str(row["Motifs"]).strip("-_").upper()
            seq   = str(row["Sequence"])
            idx   = seq.find(motif)
            return abs(idx - pos) < len(motif) + 2
        except Exception:
            return False

    merged["qc"] = merged.apply(qcheck, axis=1)
    merged = merged[merged["qc"]].drop(columns=["qc", "Sequence"], errors="ignore")
    merged["Database"] = "PTMdb"
    log.info("PTMdb %s: %d rows after QC", ptype, len(merged))
    return merged


def _phosphosite_type(path: Path, loc_df: pd.DataFrame) -> pd.DataFrame:
    """Load one PhosphoSite dataset file and merge with loc_chrom."""
    try:
        df = pd.read_csv(path, sep="\t", skiprows=3, dtype=str)
    except Exception as e:
        log.warning("PhosphoSite read error (%s): %s", path.name, e)
        return pd.DataFrame()

    if "ORGANISM" not in df.columns:
        return pd.DataFrame()
    df = df[df["ORGANISM"].str.lower() == "human"].copy()
    if df.empty:
        return pd.DataFrame()

    df = df.rename(columns={"ACC_ID": "Entry_Isoform"})

    main_df   = _main_isoform_df(loc_df)
    meta_cols = [c for c in ["Entry_Isoform", "Sequence"] if c in loc_df.columns]
    merged = pd.merge(df, main_df[meta_cols].drop_duplicates(),
                      on="Entry_Isoform", how="inner")

    def qcheck(row):
        try:
            mod_rsd = str(row.get("MOD_RSD", ""))
            pos     = int(mod_rsd.split("-")[0][1:])
            motif   = str(row.get("SITE_+/-7_AA", "")).strip("_-").upper()
            seq     = str(row.get("Sequence", ""))
            idx     = seq.find(motif)
            return abs(idx - pos) < len(motif) + 2
        except Exception:
            return False

    merged["qc"] = merged.apply(qcheck, axis=1)
    merged = merged[merged["qc"]].drop(columns=["qc", "Sequence"], errors="ignore")

    # Extract position and type from MOD_RSD (e.g. "S12-p" → pos=12, type=Phosphorylation)
    type_map = {"p": "Phosphorylation", "ac": "Acetylation", "ub": "Ubiquitination",
                "sm": "Sumoylation", "m1": "Methylation", "m2": "Methylation",
                "m3": "Methylation"}
    if "MOD_RSD" in merged.columns and not merged.empty:
        def _parse_modrsd(v):
            try:
                parts = str(v).split("-")
                pos = parts[0][1:]          # strip residue letter, keep number
                suf = parts[1] if len(parts) > 1 else ""
                return pos, type_map.get(suf, suf)
            except Exception:
                return "", ""
        _parsed        = merged["MOD_RSD"].apply(_parse_modrsd)
        merged["Position"] = _parsed.apply(lambda x: x[0])
        merged["Type"]     = _parsed.apply(lambda x: x[1])

    merged["Database"] = "PhosphoSitePlus"
    log.info("PhosphoSite %s: %d rows after QC", path.name, len(merged))
    return merged


def map_ptm(loc_df: pd.DataFrame, ptmdb_dir: str | None,
             ptmphs_dir: str | None, outdir: Path) -> int:
    out_cols = ["Entry_Isoform", "Position", "Type", "Database", "homology_transfer"]
    out_path = outdir / "ptm_merged.tsv"

    all_frames = []

    # PTMdb
    if ptmdb_dir and Path(ptmdb_dir).exists():
        for f in sorted(Path(ptmdb_dir).glob("*.tsv")):
            fr = _ptmdb_type(f, f.stem, loc_df)
            if not fr.empty:
                all_frames.append(fr)
    else:
        log.info("PTMdb dir not provided; skipping.")

    # PhosphoSite
    if ptmphs_dir and Path(ptmphs_dir).exists():
        for f in sorted(Path(ptmphs_dir).glob("*.tsv")):
            fr = _phosphosite_type(f, loc_df)
            if not fr.empty:
                all_frames.append(fr)
    else:
        log.info("PhosphoSite dir not provided; skipping.")

    if not all_frames:
        pd.DataFrame(columns=out_cols).to_csv(out_path, sep="\t", index=False)
        return 0

    combined = pd.concat(all_frames, ignore_index=True)
    combined = combined[combined["Entry_Isoform"].notna()].copy()
    combined["homology_transfer"] = False

    # Deduplicate
    combined = combined.groupby(
        ["Entry_Isoform", "Position", "Type"], as_index=False
    ).agg(Database=("Database", lambda x: " ,".join(sorted(set(x)))))
    combined["homology_transfer"] = False

    real_cols = [c for c in out_cols if c in combined.columns]
    combined[real_cols].to_csv(out_path, sep="\t", index=False)
    # Transfer to isoforms
    _apply_isoform_transfer(combined, loc_df, out_path, real_cols)
    log.info("PTM merged: %d unique sites → ptm_merged.tsv", len(combined))
    return len(combined)


# ---------------------------------------------------------------------------
# UniProt features (ROI + binding)
# ---------------------------------------------------------------------------

def _fetch_uniprot(acc: str, delay: float) -> list[dict]:
    r = _get(UNIPROT_REST.format(acc=acc), delay)
    if r is None:
        return []
    try:
        data = r.json()
    except Exception:
        return []
    rows = []
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        loc   = feat.get("location", {})
        start = loc.get("start", {}).get("value", "")
        end   = loc.get("end",   {}).get("value", "")
        rows.append({
            "Accession": acc,
            "Type":      ftype,
            "Start":     start,
            "End":       end,
            "Note":      feat.get("description", ""),
            "Evidence":  "|".join(e.get("evidenceCode", "")
                                  for e in feat.get("evidences", [])),
            "Ligand":    feat.get("ligand", {}).get("name", ""),
        })
    return rows


def map_uniprot_features_from_file(loc_df: pd.DataFrame, features_tsv: str,
                                    outdir: Path) -> tuple[int, int]:
    """Read pre-parsed uniprot_features.tsv (from parse_uniprot_dat_worker.py).

    The file contains ALL Swiss-Prot features; we filter to our accession set
    and split into ROI vs BINDING — no REST API calls needed.
    """
    src = Path(features_tsv)
    if not src.exists() or src.stat().st_size == 0:
        log.warning("--uniprot_features_tsv file missing or empty: %s", src)
        for fname in ("uniprot_roi.tsv", "uniprot_binding.tsv"):
            pd.DataFrame(columns=["Accession", "Type", "Start", "End",
                                   "Note", "Evidence", "Ligand"]).to_csv(
                outdir / fname, sep="\t", index=False)
        return 0, 0

    feat_df = pd.read_csv(src, sep="\t", dtype=str)
    accs    = set(loc_df["Entry_Isoform"].dropna().astype(str))
    # Strip isoform suffix for lookup (P04049-2 → P04049)
    canonical = {a.split("-")[0] for a in accs}
    feat_df   = feat_df[feat_df["Accession"].isin(canonical)].copy()

    roi_df  = feat_df[feat_df["Type"].isin(ROI_TYPES)]
    bind_df = feat_df[feat_df["Type"].isin(BINDING_TYPES)]

    roi_df.to_csv( outdir / "uniprot_roi.tsv",     sep="\t", index=False)
    bind_df.to_csv(outdir / "uniprot_binding.tsv", sep="\t", index=False)

    log.info("UniProt (local): ROI=%d  Binding=%d", len(roi_df), len(bind_df))
    return len(roi_df), len(bind_df)


def map_uniprot_features(loc_df: pd.DataFrame, outdir: Path,
                          delay: float) -> tuple[int, int]:
    accessions = loc_df["Entry_Isoform"].dropna().unique().tolist()
    log.info("Fetching UniProt features for %d accessions …", len(accessions))

    roi_rows, bind_rows = [], []
    for i, acc in enumerate(accessions):
        if (i + 1) % 10 == 0:
            log.info("  UniProt: %d / %d", i + 1, len(accessions))
        for f in _fetch_uniprot(acc, delay):
            if f["Type"] in ROI_TYPES:
                roi_rows.append(f)
            elif f["Type"] in BINDING_TYPES:
                bind_rows.append(f)
        time.sleep(delay)

    base_cols = ["Accession", "Type", "Start", "End", "Note", "Evidence", "Ligand"]
    for rows, fname in [(roi_rows, "uniprot_roi.tsv"),
                        (bind_rows, "uniprot_binding.tsv")]:
        pd.DataFrame(rows, columns=base_cols).to_csv(outdir / fname, sep="\t", index=False)

    log.info("UniProt: ROI=%d  Binding=%d", len(roi_rows), len(bind_rows))
    return len(roi_rows), len(bind_rows)


# ---------------------------------------------------------------------------
# Pfam via InterPro API
# ---------------------------------------------------------------------------

def map_pfam(loc_df: pd.DataFrame, outdir: Path, delay: float) -> int:
    """Fetch InterPro/Pfam per transcript (Protein_ID); one API call per Entry_Isoform."""
    acc_to_pids: dict[str, list[str]] = {}
    for _, row in loc_df.iterrows():
        acc = str(row.get("Entry_Isoform", "")).strip()
        pid = str(row.get("Protein_ID", "")).strip()
        if not acc or acc == "nan" or not pid or pid == "nan":
            continue
        acc_to_pids.setdefault(acc, [])
        if pid not in acc_to_pids[acc]:
            acc_to_pids[acc].append(pid)

    rows = []
    log.info("Fetching Pfam domains for %d accessions (%d transcripts) …",
             len(acc_to_pids), sum(len(v) for v in acc_to_pids.values()))
    for i, (acc, pids) in enumerate(acc_to_pids.items()):
        if (i + 1) % 10 == 0:
            log.info("  Pfam: %d / %d", i + 1, len(acc_to_pids))
        r = _get(INTERPRO_URL.format(acc=acc), delay)
        time.sleep(delay)
        if r is None:
            continue
        try:
            data = r.json()
        except Exception:
            continue
        for entry in data.get("results", []):
            acc_pfam = entry.get("metadata", {}).get("accession", "")
            name_pfam = entry.get("metadata", {}).get("name", {})
            if isinstance(name_pfam, dict):
                name_pfam = name_pfam.get("name", "")
            etype = entry.get("metadata", {}).get("type", "")
            for prot in entry.get("proteins", []):
                for loc_entry in prot.get("entry_protein_locations", []):
                    for frag in loc_entry.get("fragments", []):
                        base = {
                            "Accession":      acc,
                            "hmm_acc":        acc_pfam,
                            "hmm_name":       name_pfam,
                            "type":           etype,
                            "envelope_start": frag.get("start", ""),
                            "envelope_end":   frag.get("end",   ""),
                        }
                        for pid in pids:
                            rows.append({**base, "Protein_ID": pid})

    pfam_cols = ["Protein_ID", "Accession", "hmm_acc", "hmm_name", "type",
                 "envelope_start", "envelope_end"]
    pd.DataFrame(rows, columns=pfam_cols).to_csv(
        outdir / "pfam_domains.tsv", sep="\t", index=False)
    log.info("Pfam: %d domain entries across %d transcripts",
             len(rows), len({r["Protein_ID"] for r in rows}) if rows else 0)
    return len(rows)


def map_pfam_from_file(loc_df: pd.DataFrame, pfam_tsv: str,
                        outdir: Path) -> int:
    """Read pre-parsed pfam_domains.tsv (from parse_uniprot_dat_worker.py).

    Joins accession → Protein_ID using loc_df — no InterPro REST API calls.
    """
    src = Path(pfam_tsv)
    if not src.exists() or src.stat().st_size == 0:
        log.warning("--pfam_tsv file missing or empty: %s", src)
        pd.DataFrame(columns=["Protein_ID", "Accession", "hmm_acc", "hmm_name",
                               "type", "envelope_start", "envelope_end"]).to_csv(
            outdir / "pfam_domains.tsv", sep="\t", index=False)
        return 0

    pfam_df = pd.read_csv(src, sep="\t", dtype=str)
    if pfam_df.empty:
        pfam_df.to_csv(outdir / "pfam_domains.tsv", sep="\t", index=False)
        return 0

    # Build acc (canonical) → list of Protein_IDs map from loc_df
    acc_to_pids: dict[str, list[str]] = {}
    for _, row in loc_df.iterrows():
        acc = str(row.get("Entry_Isoform", "")).strip()
        pid = str(row.get("Protein_ID", "")).strip()
        if not acc or acc in ("nan", "") or not pid or pid in ("nan", ""):
            continue
        can = acc.split("-")[0]
        acc_to_pids.setdefault(can, [])
        if pid not in acc_to_pids[can]:
            acc_to_pids[can].append(pid)

    # Join: for each Pfam row, expand to all Protein_IDs sharing that accession
    rows = []
    # interpro uses 'start'/'end'; normalise to envelope_start/end for consistency
    start_col = "start" if "start" in pfam_df.columns else "envelope_start"
    end_col   = "end"   if "end"   in pfam_df.columns else "envelope_end"
    for _, r in pfam_df.iterrows():
        can  = str(r["Accession"]).split("-")[0]
        pids = acc_to_pids.get(can, [])
        for pid in pids:
            rows.append({
                "Protein_ID":     pid,
                "Accession":      r["Accession"],
                "hmm_acc":        r.get("hmm_acc", ""),
                "hmm_name":       r.get("hmm_name", ""),
                "type":           r.get("type", "Pfam"),
                "envelope_start": r.get(start_col, ""),
                "envelope_end":   r.get(end_col,   ""),
            })

    out_df = pd.DataFrame(rows, columns=["Protein_ID", "Accession", "hmm_acc",
                                          "hmm_name", "type",
                                          "envelope_start", "envelope_end"])
    out_df.to_csv(outdir / "pfam_domains.tsv", sep="\t", index=False)
    log.info("Pfam (local): %d domain entries across %d transcripts",
             len(out_df), out_df["Protein_ID"].nunique() if not out_df.empty else 0)
    return len(out_df)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Module 5: functional annotation mapping")
    p.add_argument("--loc_chrom",    required=True)
    p.add_argument("--elm_tsv",      default=None,
                   help="ELM instances TSV (legacy_data/elm/elm_instances-2023.tsv)")
    p.add_argument("--dibs_tsv",     default=None,
                   help="DIBS parsed (legacy_data/dibs/dibs_parsed.tsv)")
    p.add_argument("--mfib_tsv",     default=None,
                   help="MFIB parsed (legacy_data/mfib/mfib_parsed.tsv)")
    p.add_argument("--phasepro_tsv", default=None,
                   help="PhasePro parsed (legacy_data/phasepro/phasepro_parsed.tsv)")
    p.add_argument("--ptmdb_dir",    default=None,
                   help="Dir with PTMdb type files (legacy_data/ptm/ptmdb/)")
    p.add_argument("--ptmphs_dir",   default=None,
                   help="Dir with PhosphoSite files (legacy_data/ptm/ptmphs/)")
    p.add_argument("--output_dir",            default=".")
    p.add_argument("--request_delay",         type=float, default=0.5)
    p.add_argument("--skip_uniprot",          action="store_true", default=False)
    p.add_argument("--skip_pfam",             action="store_true", default=False)
    # Pre-parsed bulk files (from parse_uniprot_dat_worker.py).
    # When supplied, the per-protein REST API calls are replaced by local joins.
    p.add_argument("--uniprot_features_tsv",  default=None,
                   help="uniprot_features.tsv from PARSE_UNIPROT_DAT (replaces UniProt REST)")
    p.add_argument("--pfam_tsv",              default=None,
                   help="pfam_domains.tsv from PARSE_UNIPROT_DAT (replaces InterPro REST)")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    loc_df = load_loc_chrom(args.loc_chrom)

    n_elm      = map_elm(loc_df, args.elm_tsv,         outdir, args.request_delay)
    n_dibs     = map_dibs(loc_df, args.dibs_tsv,       outdir)
    n_mfib     = map_mfib(loc_df, args.mfib_tsv,       outdir)
    n_phasepro = map_phasepro(loc_df, args.phasepro_tsv, outdir)
    n_ptm      = map_ptm(loc_df, args.ptmdb_dir, args.ptmphs_dir, outdir)

    n_roi = n_bind = 0
    if args.uniprot_features_tsv and Path(args.uniprot_features_tsv).exists():
        # Fast path: pre-parsed bulk file — no REST API calls
        n_roi, n_bind = map_uniprot_features_from_file(
            loc_df, args.uniprot_features_tsv, outdir)
    elif not args.skip_uniprot:
        # Fallback: per-protein REST API (slow for full proteome)
        n_roi, n_bind = map_uniprot_features(loc_df, outdir, args.request_delay)
    else:
        for fname in ("uniprot_roi.tsv", "uniprot_binding.tsv"):
            pd.DataFrame(columns=["Accession", "Type", "Start", "End",
                                   "Note", "Evidence", "Ligand"]).to_csv(
                outdir / fname, sep="\t", index=False)

    n_pfam = 0
    if args.pfam_tsv and Path(args.pfam_tsv).exists():
        # Fast path: pre-parsed InterPro bulk file — no REST API calls
        n_pfam = map_pfam_from_file(loc_df, args.pfam_tsv, outdir)
    elif not args.skip_pfam:
        # Fallback: per-protein InterPro REST API (slow for full proteome)
        n_pfam = map_pfam(loc_df, outdir, args.request_delay)
    else:
        pd.DataFrame(columns=["Accession", "hmm_acc", "hmm_name", "type",
                               "envelope_start", "envelope_end"]).to_csv(
            outdir / "pfam_domains.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "elm_rows":      n_elm,
        "dibs_rows":     n_dibs,
        "mfib_rows":     n_mfib,
        "phasepro_rows": n_phasepro,
        "ptm_rows":      n_ptm,
        "uniprot_roi":   n_roi,
        "uniprot_bind":  n_bind,
        "pfam_rows":     n_pfam,
    }]).to_csv(outdir / "annotation_stats.tsv", sep="\t", index=False)

    log.info("Done — ELM=%d  DIBS=%d  MFIB=%d  PhasePro=%d  PTM=%d  "
             "ROI=%d  Binding=%d  Pfam=%d",
             n_elm, n_dibs, n_mfib, n_phasepro, n_ptm, n_roi, n_bind, n_pfam)


if __name__ == "__main__":
    main()
