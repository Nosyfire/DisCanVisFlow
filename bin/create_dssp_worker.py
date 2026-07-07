#!/usr/bin/env python3
"""create_dssp_worker.py — DSSP secondary structure + true RSA from AlphaFold.

For each protein, resolves an AlphaFold mmCIF (local --cif_dir cache first, then
EBI download), runs `mkdssp`, and emits per-residue 8-state SS, collapsed 3-state
SS, and true RSA = solvent-accessible area ÷ Tien (2013) theoretical max-ASA.

This is a *true* RSA track, distinct from the pLDDT-pseudo-RSA in rsa_scores.tsv.

Output
------
dssp.tsv — columns: Protein_ID  Position  aa  ss8  ss3  rsa
Proteins with no AlphaFold model are skipped (logged).
"""
import argparse
import json
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

try:
    import requests
    _HAS_REQUESTS = True
except Exception:
    _HAS_REQUESTS = False

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

OUT_COLS = ["Protein_ID", "Position", "aa", "ss8", "ss3", "rsa"]

ALPHAFOLD_SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
ALPHAFOLD_CIF_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.cif"
_ALPHAFOLD_FALLBACK_VERSIONS = (6, 5, 4)

# DSSP 8-state code -> collapsed 3-state (H=helix, E=strand, C=coil/other).
_SS8_TO_SS3 = {"H": "H", "G": "H", "I": "H",       # helices
               "E": "E", "B": "E",                  # strands
               "T": "C", "S": "C", "P": "C",        # turns/bends/PPII -> coil
               "-": "C", " ": "C", "": "C", "C": "C"}

# Tien et al. (2013) theoretical maximum ASA (A^2), 1-letter AA.
_MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


def _get(url: str, delay: float, timeout: int = 60):
    if not _HAS_REQUESTS:
        return None
    try:
        r = requests.get(url, timeout=timeout)
        time.sleep(delay)
        return r if r.status_code == 200 else None
    except Exception:
        return None


def _resolve_cif(acc_base: str, cif_dir: Path, delay: float,
                 allow_download: bool) -> "Path | None":
    """Local cache AF-<acc>-F1*.cif first, else download to cif_dir."""
    hits = sorted(cif_dir.glob(f"AF-{acc_base}-F1*.cif"))
    if hits:
        return hits[0]
    if not allow_download:
        return None
    cif_dir.mkdir(parents=True, exist_ok=True)
    # summary API -> current cifUrl
    summ = _get(ALPHAFOLD_SUMMARY_URL.format(acc=acc_base), delay=0, timeout=15)
    urls = []
    if summ is not None:
        try:
            data = json.loads(summ.text)
            u = data[0].get("cifUrl") if isinstance(data, list) and data else None
            if u:
                urls.append(u)
        except Exception:
            pass
    urls += [ALPHAFOLD_CIF_URL.format(acc=acc_base, ver=v)
             for v in _ALPHAFOLD_FALLBACK_VERSIONS]
    for u in urls:
        r = _get(u, delay)
        if r is not None:
            dest = cif_dir / f"AF-{acc_base}-F1-model.cif"
            dest.write_text(r.text)
            return dest
    return None


def _run_dssp(cif: Path, mkdssp: str, tmp: Path) -> list:
    """Return list of (resnum:int, aa:str, ss8:str, acc:float) from mkdssp mmCIF output."""
    out = tmp / "out.cif"
    proc = subprocess.run(
        [mkdssp, "--output-format", "mmcif", "--calculate-accessibility",
         str(cif), str(out)],
        capture_output=True, text=True)
    if proc.returncode != 0 or not out.exists():
        log.warning("mkdssp failed on %s: %s", cif.name, proc.stderr.strip()[:200])
        return []
    return _parse_dssp_mmcif(out.read_text())


def _parse_dssp_mmcif(text: str) -> list:
    """Parse mkdssp 4.x mmCIF output loop `_dssp_struct_summary`.

    mkdssp 4.6.1 writes an mmCIF with a `_dssp_struct_summary` category holding
    per-residue label_seq_id, label_comp_id, secondary_structure, and
    accessibility (only populated when mkdssp is run with
    --calculate-accessibility — otherwise every row is '.'). Column order is
    read from the loop header (robust to reordering across mkdssp point
    releases).
    """
    rows: list = []
    lines = text.splitlines()
    i = 0
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    while i < len(lines):
        if lines[i].strip() == "loop_":
            # collect the header keys that follow
            keys = []
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                keys.append(lines[j].strip())
                j += 1
            if any(k.startswith("_dssp_struct_summary.") for k in keys):
                idx = {k.split(".", 1)[1]: n for n, k in enumerate(keys)}
                need = ("label_seq_id", "label_comp_id",
                        "secondary_structure", "accessibility")
                if all(n in idx for n in need):
                    k = j
                    while k < len(lines):
                        row = lines[k].strip()
                        if row.startswith("#") or row.startswith("_") \
                                or row == "loop_" or row == "":
                            break
                        toks = row.split()
                        if len(toks) >= len(keys):
                            try:
                                seqid = int(toks[idx["label_seq_id"]])
                            except ValueError:
                                k += 1
                                continue
                            comp = toks[idx["label_comp_id"]]
                            ss = toks[idx["secondary_structure"]]
                            if ss in (".", "?"):
                                ss = "-"
                            try:
                                acc = float(toks[idx["accessibility"]])
                            except ValueError:
                                acc = float("nan")
                            aa = three_to_one.get(comp, "X")
                            rows.append((seqid, aa, ss, acc))
                        k += 1
                    i = k
                    continue
            i = j
            continue
        i += 1
    return rows


def main():
    ap = argparse.ArgumentParser(description="DSSP SS + true RSA from AlphaFold mmCIF")
    ap.add_argument("--seq_table", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--cif_dir", default="references/alphafold_cif")
    ap.add_argument("--only_main_isoforms", action="store_true", default=False)
    ap.add_argument("--mkdssp", default="mkdssp")
    ap.add_argument("--delay", type=float, default=0.2)
    ap.add_argument("--no_download", action="store_true", default=False)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = outdir / "dssp.tsv"
    cif_dir = Path(args.cif_dir)

    if shutil.which(args.mkdssp) is None:
        log.warning("mkdssp not found on PATH — skipping DSSP track (empty output).")
        pd.DataFrame(columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)
        sys.exit(0)

    df = pd.read_csv(args.seq_table, sep="\t", dtype=str).dropna(subset=["Sequence"])
    if args.only_main_isoforms and "main_isoform" in df.columns:
        df = df[df["main_isoform"] == "yes"]

    # AlphaFold only has canonical accessions -> strip isoform suffix, dedup.
    df = df.copy()
    df["acc_base"] = df["Entry_Isoform"].fillna("").str.split("-").str[0]

    pd.DataFrame(columns=OUT_COLS).to_csv(out_tsv, sep="\t", index=False)

    import tempfile
    ok = skipped = 0
    buffer: list = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for _, r in df.iterrows():
            pid = r["Protein_ID"]
            acc = r["acc_base"]
            if not acc:
                skipped += 1
                continue
            cif = _resolve_cif(acc, cif_dir, args.delay, not args.no_download)
            if cif is None:
                skipped += 1
                continue
            residues = _run_dssp(cif, args.mkdssp, tmp)
            if not residues:
                skipped += 1
                continue
            for resnum, aa, ss8, acc_area in residues:
                ss8n = ss8 if ss8 not in ("-", " ", "") else "C"
                ss3 = _SS8_TO_SS3.get(ss8, "C")
                maxasa = _MAX_ASA.get(aa)
                rsa = (acc_area / maxasa) if (maxasa and acc_area == acc_area) else ""
                buffer.append({
                    "Protein_ID": pid, "Position": resnum, "aa": aa,
                    "ss8": ss8n, "ss3": ss3,
                    "rsa": round(rsa, 4) if rsa != "" else "",
                })
            ok += 1
            if len(buffer) >= 50000:
                pd.DataFrame(buffer, columns=OUT_COLS).to_csv(
                    out_tsv, mode="a", header=False, sep="\t", index=False)
                buffer = []
        if buffer:
            pd.DataFrame(buffer, columns=OUT_COLS).to_csv(
                out_tsv, mode="a", header=False, sep="\t", index=False)
    log.info("Done — %d proteins with DSSP, %d skipped (no model)", ok, skipped)


if __name__ == "__main__":
    main()
