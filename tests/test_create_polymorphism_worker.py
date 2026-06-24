"""Tests for bin/create_polymorphism_worker.py — legacy SNP .out → protein mapping."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin" / "create_polymorphism_worker.py"


def _write_map(path: Path):
    # GENE-201: residue 1 (M) codon at genomic 100,101,102; residue 2 (F) at 200,201,202
    lines = [
        "# ENST0001.1|ENSG1|OTT1|-|GENE-201|GENE|9|CDS:1-6| chr3 - 1-300",
        "0 M 1,2,3 ATG M 100,101,102, ATG M",
        "1 F 4,5,6 TTC F 200,201,202, TTC F",
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_loc(path: Path):
    pd.DataFrame({"Protein_ID": ["GENE-201"]}).to_csv(path, sep="\t", index=False)


def _write_out(path: Path, gstart, gend, pid, rsid="rs1", ref="T", alt="C,",
               freqs="-inf,0.0195687,-inf,0.428571,-inf,"):
    # UCSC BED-style: enough columns; last column = Protein_ID, col2 = 1-based end,
    # col10 (index 9) = per-population minor allele frequency CSV (-inf = missing)
    cols = ["chr3", str(gstart), str(gend), rsid, ref, "1", alt, "0", "31",
            freqs, "refcol", "altcol", "9", "snv", "flags", "id", "999", pid]
    path.write_text("\t".join(cols) + "\n")


def test_snp_maps_to_protein_position(tmp_path):
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    common = tmp_path / "common_poly.out"
    _write_map(mp)
    _write_loc(loc)
    # SNP at genomic base 201 (BED end) → residue 2 (F)
    _write_out(common, 200, 201, "GENE-201")

    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--combined_map", str(mp), "--snp_common", str(common),
                    "--output_dir", str(tmp_path)], check=True)

    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    assert list(df.columns) == ["Protein_ID", "Position", "rsid", "ref", "alt",
                                "allele_frequency", "Type"]
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Protein_ID"] == "GENE-201"
    assert row["Position"] == "2"          # genomic 201 → residue 2 (1-based)
    assert row["rsid"] == "rs1"
    assert row["ref"] == "T" and row["alt"] == "C"
    assert row["allele_frequency"] == "0.428571"   # max finite, -inf ignored
    assert row["Type"] == "Common Polymorphisms"


def test_utr_snp_not_in_cds_is_skipped(tmp_path):
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    common = tmp_path / "common_poly.out"
    _write_map(mp)
    _write_loc(loc)
    # genomic 999 is not part of any codon → must be dropped
    _write_out(common, 998, 999, "GENE-201")
    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--combined_map", str(mp), "--snp_common", str(common),
                    "--output_dir", str(tmp_path)], check=True)
    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    assert len(df) == 0


def test_no_combined_map_writes_empty(tmp_path):
    loc = tmp_path / "loc.tsv"
    _write_loc(loc)
    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--combined_map", str(tmp_path / "missing.map"),
                    "--output_dir", str(tmp_path)], check=True)
    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    assert len(df) == 0
    assert list(df.columns) == ["Protein_ID", "Position", "rsid", "ref", "alt",
                                "allele_frequency", "Type"]


def _write_pos_tsv(path: Path, rows):
    """rows: list of (Protein_ID|pos, type)"""
    with open(path, "w") as fh:
        fh.write("AccessionPosition\tPolymorphism\n")
        for key, typ in rows:
            fh.write(f"{key}\t{typ}\n")


def test_snp_pos_tsv_folds_in_all_polymorphisms(tmp_path):
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    common = tmp_path / "common_poly.out"
    pos = tmp_path / "polymorphism_pos.tsv"
    _write_map(mp)
    _write_loc(loc)
    # one allele-frequency-bearing common SNP at residue 2
    _write_out(common, 200, 201, "GENE-201")
    # comprehensive set: residue 2 (already covered → enriched row wins) + residue 1 (new)
    _write_pos_tsv(pos, [("GENE-201|2", "Common Polymorphisms"),
                         ("GENE-201|1", "All Polymorphisms"),
                         ("OTHER-201|5", "Polymorphism")])  # filtered out (not a run protein)

    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--combined_map", str(mp), "--snp_common", str(common),
                    "--snp_pos_tsv", str(pos), "--output_dir", str(tmp_path)], check=True)
    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    # residue 2 (with freq) + residue 1 (from pos tsv, no freq) = 2 rows
    assert len(df) == 2
    assert set(df["Position"]) == {"1", "2"}
    assert "OTHER-201" not in set(df["Protein_ID"])
    # residue 2 keeps the allele frequency from common_poly (not duplicated)
    r2 = df[df["Position"] == "2"].iloc[0]
    assert r2["allele_frequency"] == "0.428571"
    # residue 1 comes from the pos tsv → no allele frequency
    r1 = df[df["Position"] == "1"].iloc[0]
    assert pd.isna(r1["allele_frequency"]) or str(r1["allele_frequency"]) == ""
    assert r1["Type"] == "All Polymorphisms"


def test_snp_pos_tsv_without_combined_map(tmp_path):
    """polymorphism_pos.tsv alone (no genome map) still yields the all-poly set."""
    loc = tmp_path / "loc.tsv"
    pos = tmp_path / "polymorphism_pos.tsv"
    _write_loc(loc)
    _write_pos_tsv(pos, [("GENE-201|10", "Common Polymorphisms"),
                         ("GENE-201|20", "All Polymorphisms")])
    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--snp_pos_tsv", str(pos), "--output_dir", str(tmp_path)], check=True)
    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    assert set(df["Position"]) == {"10", "20"}


def _write_map_two_isoforms(path: Path):
    # Two isoforms of the same gene sharing genomic codon at base 201 (residue 2 / residue 1).
    lines = [
        "# ENST0001.1|ENSG1|OTT1|-|GENE-201|GENE|9|CDS:1-6| chr3 - 1-300",
        "0 M 1,2,3 ATG M 100,101,102, ATG M",
        "1 F 4,5,6 TTC F 200,201,202, TTC F",
        "# ENST0002.1|ENSG1|OTT1|-|GENE-202|GENE|9|CDS:1-3| chr3 - 1-300",
        "0 F 1,2,3 TTC F 200,201,202, TTC F",
    ]
    path.write_text("\n".join(lines) + "\n")


def _make_fake_bigbedtobed(bindir: Path, bed_lines):
    """Create a fake bigBedToBed that ignores args and writes fixed BED to the last arg."""
    bindir.mkdir(parents=True, exist_ok=True)
    tool = bindir / "bigBedToBed"
    body = (bed_lines).replace("\n", "\\n")
    tool.write_text(
        "#!/usr/bin/env bash\n"
        "out=\"${@: -1}\"\n"
        f"printf '{body}' > \"$out\"\n"
    )
    tool.chmod(0o755)
    return bindir


def test_dbsnp_bigbed_extracts_freq_for_all_isoforms(tmp_path):
    """With --dbsnp_bb, every selected isoform that contains the SNV codon gets a
    row carrying rsid + allele frequency (bigBed path supersedes the pos table)."""
    mp = tmp_path / "combined_map.map"
    loc = tmp_path / "loc.tsv"
    _write_map_two_isoforms(mp)
    pd.DataFrame({"Protein_ID": ["GENE-201", "GENE-202"]}).to_csv(loc, sep="\t", index=False)

    # dbSnp155 bigBed BED row: chrom start end rsid ref altCount alts ... freqs(col10) ... flags(col15)
    bed = ("chr3\t200\t201\trs99\tT\t1\tC,\t0\t31\t"
           "-inf,0.0195687,-inf,0.428571,-inf,\trefcol\taltcol\t9\tsnv\tcommonSome,\n")
    bindir = _make_fake_bigbedtobed(tmp_path / "bin", bed)
    dbsnp = tmp_path / "dbSnp155Common.bb"
    dbsnp.write_bytes(b"\x00fake-bigbed")   # non-empty so worker uses it

    subprocess.run([sys.executable, str(BIN), "--loc_chrom", str(loc),
                    "--combined_map", str(mp), "--dbsnp_bb", str(dbsnp),
                    "--ucsc_bin", str(bindir), "--output_dir", str(tmp_path)], check=True)
    df = pd.read_csv(tmp_path / "polymorphism.tsv", sep="\t", dtype=str)
    # SNV at genomic 201 maps to residue 2 of GENE-201 and residue 1 of GENE-202
    assert set(df["Protein_ID"]) == {"GENE-201", "GENE-202"}
    assert df[df["Protein_ID"] == "GENE-201"].iloc[0]["Position"] == "2"
    assert df[df["Protein_ID"] == "GENE-202"].iloc[0]["Position"] == "1"
    # allele frequency populated for ALL rows
    assert (df["allele_frequency"].astype(str) == "0.428571").all()
    assert (df["rsid"] == "rs99").all()
    assert (df["Type"] == "Common Polymorphisms").all()   # commonSome flag


def test_classify_type_and_parse_map_regions():
    sys.path.insert(0, str(BIN.parent))
    import importlib
    mod = importlib.import_module("create_polymorphism_worker")
    assert mod._classify_type("refIsMinor,commonSome,") == "Common Polymorphisms"
    assert mod._classify_type("rareSome,") == "All Polymorphisms"
    assert mod._classify_type("") == "All Polymorphisms"
