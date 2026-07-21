#!/usr/bin/env python3
"""Build a per-variant ClinVar submission-date table.

ClinVar's VCF carries no dates, so the dates come from ``submission_summary.txt.gz``
(one row per submitted record / SCV). Those rows are aggregated per VariationID
into a single "when was this variant last touched by a submitter" summary:

  last_evaluated   max(DateLastEvaluated) across all SCVs
  first_evaluated  min(DateLastEvaluated) across all SCVs
  n_submissions    number of SCV records (including undated ones)
  last_submitter / last_clinical_significance / last_review_status
                   taken from the SCV holding `last_evaluated`

The table is keyed on the ClinVar VCF's own ``CLNHGVS`` string, which is exactly
the ``Mutation`` column emitted by create_mutation_map_worker.py — so it joins
onto the mapped mutation TSVs with no coordinate reconstruction. Every variant in
the VCF is emitted, including those with no submission rows (empty dates,
``n_submissions=0``), so the table is a complete companion to the VCF.

Joining by VariationID rather than by coordinate also makes the release skew
between the (pinned) VCF and the (always-current) submission file harmless:
VariationIDs are stable, unmatched submissions are ignored.
"""
import argparse
import gzip
import logging
import re
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger("clinvar_dates")

OUT_COLS = ["VariationID", "Mutation", "Chromosome", "Position", "Ref", "Alt", "Gene",
            "last_evaluated", "first_evaluated", "n_submissions", "last_submitter",
            "last_clinical_significance", "last_review_status"]

CLNHGVS_RE  = re.compile(r"CLNHGVS=([^;]+)")
GENEINFO_RE = re.compile(r"GENEINFO=([^;:]+)")
_EMPTY = {"", "-", "na", "n/a", "none"}


def _open(path):
    p = str(path)
    return gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, errors="replace")


def parse_date(raw: str) -> str:
    """'Dec 17, 2024' → '2024-12-17'. Unparseable/absent → ''."""
    s = (raw or "").strip()
    if s.lower() in _EMPTY:
        return ""
    for fmt in ("%b %d, %Y", "%Y-%m-%d", "%b %Y", "%Y"):
        try:
            return datetime.strptime(s, fmt).strftime("%Y-%m-%d")
        except ValueError:
            continue
    return ""


def load_submissions(path) -> dict:
    """VariationID → [last, first, n, submitter, clnsig, review].

    Built for every VariationID in the file, not just those in the VCF: holding a
    second dict of VCF variants purely to pre-filter costs more memory than the
    unmatched entries it would drop. Unmatched IDs are ignored at write time.
    """
    agg = {}
    seen_header = False
    rows = 0
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                # the final comment line IS the header; everything before is prose
                if line.startswith("#VariationID"):
                    seen_header = True
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            vid = f[0].strip()
            rows += 1
            date = parse_date(f[2])
            sig, review, submitter = f[1].strip(), f[6].strip(), f[9].strip()
            rec = agg.get(vid)
            if rec is None:
                # n starts at 1; dateless rows count but leave last/first empty
                agg[vid] = [date, date, 1, submitter, sig, review] if date \
                    else ["", "", 1, "", "", ""]
                continue
            rec[2] += 1
            if not date:
                continue
            if not rec[0] or date > rec[0]:
                rec[0], rec[3], rec[4], rec[5] = date, submitter, sig, review
            if not rec[1] or date < rec[1]:
                rec[1] = date
    if not seen_header:
        log.warning("submission_summary: no '#VariationID' header line found")
    log.info("Submissions: %d rows over %d variants", rows, len(agg))
    return agg


def write_table(vcf_path, agg, out_path):
    """Stream the VCF and join each variant to its aggregated submission dates."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    n_var = n_dated = 0
    empty = ["", "", 0, "", "", ""]
    with _open(vcf_path) as src, open(out, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for line in src:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt, info = f[0], f[1], f[2], f[3], f[4], f[7]
            if not vid or vid == ".":
                continue
            m_h = CLNHGVS_RE.search(info)
            m_g = GENEINFO_RE.search(info)
            last, first, n, submitter, sig, review = agg.get(vid, empty)
            n_var += 1
            if last:
                n_dated += 1
            fh.write("\t".join([vid, m_h.group(1) if m_h else "", chrom, pos, ref, alt,
                                m_g.group(1) if m_g else "", last, first, str(n),
                                submitter, sig, review]) + "\n")
    log.info("Written %s — %d variants, %d with a submission date", out, n_var, n_dated)


def main():
    p = argparse.ArgumentParser(description="Per-variant ClinVar submission dates")
    p.add_argument("--clinvar_vcf", required=True, help="ClinVar GRCh38 VCF (.vcf[.gz])")
    p.add_argument("--submission_summary", required=True,
                   help="ClinVar submission_summary.txt[.gz]")
    p.add_argument("--output", required=True, help="output TSV path")
    args = p.parse_args()

    write_table(args.clinvar_vcf, load_submissions(args.submission_summary), args.output)


if __name__ == "__main__":
    main()
