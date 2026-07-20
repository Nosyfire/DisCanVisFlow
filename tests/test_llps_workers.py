"""Tests for the LLPS-database workers (parse + region map + variant map).

Each worker is exercised as a subprocess with tiny hand-built inputs, matching
the pattern of the other bin/*.py tests. No Nextflow, no network, no real DB.
"""
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[1] / "bin"
PARSE = BIN / "parse_llps_sources.py"
REGIONS = BIN / "create_llps_regions_worker.py"
VARIANTS = BIN / "create_llps_variants_worker.py"


def _loc(tmp_path):
    """Two isoforms of one accession with different sequences + an unrelated one."""
    seq_main = "MKDEFGHIKL" + "MNPQRSTVWY" + "ACDEACDEAC"   # 30 aa, Q1 main
    seq_alt  = "XX" + seq_main + "GGGGGGGGGG"                # 42 aa, shifted, Q1-2
    seq_other = "GGGGWWWWCCCC"                               # 12 aa, P9
    p = tmp_path / "loc.tsv"
    pd.DataFrame(
        [["Q1-1",  "Q1",   "GENA", "yes", seq_main],
         ["Q1-2",  "Q1-2", "GENA", "no",  seq_alt],
         ["P9-1",  "P9",   "GENB", "yes", seq_other]],
        columns=["Protein_ID", "Entry_Isoform", "Gene", "main_isoform", "Sequence"],
    ).to_csv(p, sep="\t", index=False)
    return p


def _norm(tmp_path):
    d = tmp_path / "norm"
    d.mkdir()
    pd.DataFrame(
        [["Q1", 1, 10, "PS_region", "phasepdb"],
         ["Q1", 5, 40, "IDR", "llpsdb"],          # end 40 > main len 30 → dropped for Q1-1
         ["P9", 1, 4, "LCR", "llpsdb"]],
        columns=["acc", "start", "end", "region_type", "source"],
    ).to_csv(d / "llps_regions.tsv", sep="\t", index=False)
    pd.DataFrame(
        [["Q1", "GENA", "Nucleus", "liquid", "PS-self", "phasepdb"],
         ["P9", "GENB", "Cytoplasm", "", "", "llpsdb"]],
        columns=["acc", "gene", "organelle_mlo", "material_state", "klass", "source"],
    ).to_csv(d / "llps_proteins.tsv", sep="\t", index=False)
    pd.DataFrame(
        # M at pos 1 matches main seq; K at pos 2 matches main; Z at 3 matches nothing
        [["Q1", 1, "M", "P", "experimental", "", "phasepdb"],
         ["Q1", 2, "K", "E", "disease", "myopathy", "disphasedb"],
         ["Q1", 3, "Z", "A", "experimental", "", "phasepdb"]],
        columns=["acc", "pos", "wt_aa", "mut_aa", "variant_class", "disease_term", "source"],
    ).to_csv(d / "llps_variants.tsv", sep="\t", index=False)
    return d


def test_regions_residue_validation_and_source_tag(tmp_path):
    loc, norm = _loc(tmp_path), _norm(tmp_path)
    out = tmp_path / "out"
    res = subprocess.run(
        [sys.executable, str(REGIONS), "--regions", str(norm / "llps_regions.tsv"),
         "--proteins", str(norm / "llps_proteins.tsv"), "--loc_chrom", str(loc),
         "--outdir", str(out)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    reg = pd.read_csv(out / "llps_regions.tsv", sep="\t")
    # PS_region 1-10 fits both Q1 isoforms; IDR 5-40 overflows main (len 30) but
    # fits the alt isoform (len 32) → present only for Q1-2; LCR 1-4 fits P9.
    q1_main = reg[(reg.Protein_ID == "Q1-1")]
    assert set(q1_main.region_type) == {"PS_region"}          # IDR dropped for short main
    assert "IDR" in set(reg[reg.Protein_ID == "Q1-2"].region_type)
    assert set(reg[reg.Protein_ID == "P9-1"].region_type) == {"LCR"}
    assert set(reg.source) >= {"phasepdb", "llpsdb"}
    assert (reg.mapping_type == "direct").all()


def test_proteins_map_all_isoforms_of_accession(tmp_path):
    loc, norm = _loc(tmp_path), _norm(tmp_path)
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(REGIONS), "--proteins", str(norm / "llps_proteins.tsv"),
         "--loc_chrom", str(loc), "--outdir", str(out)], check=True,
        capture_output=True, text=True)
    prot = pd.read_csv(out / "llps_proteins.tsv", sep="\t")
    # Q1 protein info attaches to BOTH Q1 isoforms
    assert set(prot[prot.Gene == "GENA"].Protein_ID) == {"Q1-1", "Q1-2"}
    assert prot[prot.Protein_ID == "Q1-1"].iloc[0].organelle_mlo == "Nucleus"


def test_variants_residue_validation_and_class(tmp_path):
    loc, norm = _loc(tmp_path), _norm(tmp_path)
    out = tmp_path / "out"
    res = subprocess.run(
        [sys.executable, str(VARIANTS), "--variants", str(norm / "llps_variants.tsv"),
         "--loc_chrom", str(loc), "--outdir", str(out)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    var = pd.read_csv(out / "llps_variants.tsv", sep="\t")
    # pos1 M→P and pos2 K→E match the main isoform; pos3 Z→A matches nothing → dropped
    main = var[var.Protein_ID == "Q1-1"]
    assert set(zip(main.Protein_position, main.WT_AA)) == {(1, "M"), (2, "K")}
    assert "Z" not in set(var.WT_AA)
    # variant_class + disease_term carried through and source tagged
    dis = var[(var.WT_AA == "K")].iloc[0]
    assert dis.variant_class == "disease" and dis.disease_term == "myopathy"
    assert dis.source_db == "disphasedb"


def test_empty_inputs_produce_headers_and_exit_zero(tmp_path):
    loc = _loc(tmp_path)
    out = tmp_path / "out"
    r1 = subprocess.run([sys.executable, str(REGIONS), "--loc_chrom", str(loc),
                         "--outdir", str(out)], capture_output=True, text=True)
    r2 = subprocess.run([sys.executable, str(VARIANTS), "--loc_chrom", str(loc),
                         "--outdir", str(out)], capture_output=True, text=True)
    assert r1.returncode == 0 and r2.returncode == 0, r1.stderr + r2.stderr
    assert list(pd.read_csv(out / "llps_regions.tsv", sep="\t").columns)[:2] == ["Protein_ID", "Gene"]
    assert "Protein_position" in pd.read_csv(out / "llps_variants.tsv", sep="\t").columns


def test_parse_hgvs_and_ranges(tmp_path):
    """The parser turns a DisPhaseDB-shaped CSV into normalised variants."""
    disp = tmp_path / "disphasedb.csv"
    disp.write_text(textwrap.dedent("""\
        uniprot,variant,disease
        Q1,p.M1P,cancer
        Q1,p.(K2E),myopathy
        Q1,notavariant,x
    """))
    out = tmp_path / "norm"
    res = subprocess.run(
        [sys.executable, str(PARSE), "--disphasedb", str(disp), "--outdir", str(out)],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    v = pd.read_csv(out / "llps_variants.tsv", sep="\t")
    assert set(zip(v.wt_aa, v.pos, v.mut_aa)) == {("M", 1, "P"), ("K", 2, "E")}
    assert (v.variant_class == "disease").all()
    assert set(v.disease_term) == {"cancer", "myopathy"}
