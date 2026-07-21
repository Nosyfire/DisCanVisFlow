"""Tests for create_clinvar_dates_worker.py (ClinVar submission dates).

The worker joins per-SCV submission dates (submission_summary.txt.gz) onto the
ClinVar VCF by VariationID, emitting one compact row per variant keyed by the
same CLNHGVS string the mutation maps already carry. Exercised as a subprocess
with tiny hand-built inputs — no network, no real ClinVar.
"""
import gzip
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd

WORKER = Path(__file__).resolve().parents[1] / "bin" / "create_clinvar_dates_worker.py"

SUBM_HEADER = (
    "#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\t"
    "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\tCollectionMethod\t"
    "OriginCounts\tSubmitter\tSCV\tSubmittedGeneSymbol\tExplanationOfInterpretation\t"
    "SomaticClinicalImpact\tOncogenicity\tContributesToAggregateClassification"
)


def _subm_row(vid, sig, date, review, submitter, scv):
    return "\t".join([str(vid), sig, date, "-", "-", "-", review, "clinical testing",
                      "-", submitter, scv, "-", "-", "-", "-", "yes"])


def _vcf(tmp_path, records):
    """records: list of (chrom, pos, variation_id, ref, alt, clnhgvs, gene)."""
    p = tmp_path / "clinvar.vcf.gz"
    lines = [
        "##fileformat=VCFv4.1",
        "##fileDate=2026-03-09",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for chrom, pos, vid, ref, alt, hgvs, gene in records:
        info = (f"ALLELEID=999;CLNHGVS={hgvs};CLNSIG=Uncertain_significance;"
                f"CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO={gene}:1")
        lines.append(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}")
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return p


def _subm(tmp_path, rows):
    p = tmp_path / "submission_summary.txt.gz"
    body = ["## a comment line", "## another", SUBM_HEADER] + rows
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(body) + "\n")
    return p


def _run(tmp_path, vcf, subm):
    out = tmp_path / "clinvar_submission_dates.tsv"
    r = subprocess.run(
        [sys.executable, str(WORKER), "--clinvar_vcf", str(vcf),
         "--submission_summary", str(subm), "--output", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    return pd.read_csv(out, sep="\t", dtype=str).fillna("")


def test_aggregates_max_min_and_count_per_variation(tmp_path):
    vcf = _vcf(tmp_path, [("7", 4781213, 2, "G", "A",
                           "NC_000007.14:g.4781213G>A", "AP5Z1")])
    subm = _subm(tmp_path, [
        _subm_row(2, "Likely pathogenic", "Dec 17, 2024",
                  "criteria provided, single submitter", "Basel Group", "SCV005909190.1"),
        _subm_row(2, "Pathogenic", "Jun 25, 2024",
                  "criteria provided, single submitter", "Athena", "SCV005622007.1"),
        _subm_row(2, "Pathogenic", "Jun 29, 2010",
                  "no assertion criteria provided", "OMIM", "SCV000020155.3"),
    ])
    df = _run(tmp_path, vcf, subm)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["last_evaluated"] == "2024-12-17"
    assert row["first_evaluated"] == "2010-06-29"
    assert row["n_submissions"] == "3"
    assert row["last_submitter"] == "Basel Group"
    assert row["last_clinical_significance"] == "Likely pathogenic"
    assert row["last_review_status"] == "criteria provided, single submitter"


def test_undated_submission_counts_but_does_not_shift_dates(tmp_path):
    vcf = _vcf(tmp_path, [("1", 100, 5, "C", "T", "NC_000001.11:g.100C>T", "GENA")])
    subm = _subm(tmp_path, [
        _subm_row(5, "Pathogenic", "Mar 12, 2024", "criteria provided", "Fulgent", "SCV1.1"),
        # dateless submission — ClinVar writes a bare '-'
        _subm_row(5, "Pathogenic", "-", "criteria provided", "Paris Brain", "SCV2.1"),
    ])
    df = _run(tmp_path, vcf, subm)
    row = df.iloc[0]
    assert row["n_submissions"] == "2"
    assert row["last_evaluated"] == "2024-03-12"
    assert row["first_evaluated"] == "2024-03-12"
    assert row["last_submitter"] == "Fulgent"


def test_join_key_is_clnhgvs_and_coordinates_are_carried(tmp_path):
    vcf = _vcf(tmp_path, [("1", 69134, 2205837, "A", "G",
                           "NC_000001.11:g.69134A>G", "OR4F5")])
    subm = _subm(tmp_path, [
        _subm_row(2205837, "Likely benign", "Feb 01, 2023", "criteria provided",
                  "Lab X", "SCV003526545.1"),
    ])
    df = _run(tmp_path, vcf, subm)
    row = df.iloc[0]
    assert row["Mutation"] == "NC_000001.11:g.69134A>G"
    assert row["VariationID"] == "2205837"
    assert (row["Chromosome"], row["Position"], row["Ref"], row["Alt"]) == \
           ("1", "69134", "A", "G")
    assert row["Gene"] == "OR4F5"


def test_variant_without_submissions_is_kept_with_empty_dates(tmp_path):
    vcf = _vcf(tmp_path, [
        ("1", 100, 11, "C", "T", "NC_000001.11:g.100C>T", "GENA"),
        ("2", 200, 22, "G", "A", "NC_000002.12:g.200G>A", "GENB"),
    ])
    subm = _subm(tmp_path, [
        _subm_row(11, "Pathogenic", "Mar 12, 2024", "criteria provided", "Lab", "SCV1.1"),
    ])
    df = _run(tmp_path, vcf, subm)
    assert len(df) == 2
    orphan = df[df["VariationID"] == "22"].iloc[0]
    assert orphan["last_evaluated"] == ""
    assert orphan["n_submissions"] == "0"


def test_submissions_for_unknown_variation_are_ignored(tmp_path):
    vcf = _vcf(tmp_path, [("1", 100, 11, "C", "T", "NC_000001.11:g.100C>T", "GENA")])
    subm = _subm(tmp_path, [
        _subm_row(11, "Pathogenic", "Mar 12, 2024", "criteria provided", "Lab", "SCV1.1"),
        _subm_row(999999, "Pathogenic", "Dec 31, 2025", "criteria provided", "Ghost", "SCV9.1"),
    ])
    df = _run(tmp_path, vcf, subm)
    assert list(df["VariationID"]) == ["11"]


def test_empty_submission_file_still_emits_all_variants(tmp_path):
    vcf = _vcf(tmp_path, [("1", 100, 11, "C", "T", "NC_000001.11:g.100C>T", "GENA")])
    subm = _subm(tmp_path, [])
    df = _run(tmp_path, vcf, subm)
    assert len(df) == 1
    assert df.iloc[0]["n_submissions"] == "0"
    assert df.iloc[0]["last_evaluated"] == ""
