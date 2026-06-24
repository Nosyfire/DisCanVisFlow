"""Tests for create_cancer_driver_worker.py — combined membership + role enrichment."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

WORKER = Path(__file__).parent.parent / "bin" / "create_cancer_driver_worker.py"


def _make_seq(tmp, proteins):
    p = tmp / "seq.tsv"
    pd.DataFrame({"Protein_ID": proteins}).to_csv(p, sep="\t", index=False)
    return p


def _make_combined(tmp, rows):
    p = tmp / "cancer_driver_src.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _make_census_roles(tmp, rows):
    p = tmp / "census_roles.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _make_compendium_roles(tmp, rows):
    p = tmp / "compendium_roles.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _run(seq, combined, census_roles, compendium_roles, outdir):
    return subprocess.run(
        [sys.executable, str(WORKER),
         "--seq_table", str(seq),
         "--cancer_driver", str(combined),
         "--census_roles", str(census_roles),
         "--compendium_roles", str(compendium_roles),
         "--outdir", str(outdir)],
        capture_output=True, text=True,
    )


class TestRoleEnrichment:
    def test_census_carries_role(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201", "RAF1-262", "TP53-201"])
        combined = _make_combined(tmp_path, [
            {"Protein_ID": "RAF1-201", "Cancer Driver": "Census, Compendium"},
            {"Protein_ID": "RAF1-262", "Cancer Driver": "Census"},
            {"Protein_ID": "TP53-201", "Cancer Driver": "Compendium"},
            {"Protein_ID": "BRAF-201", "Cancer Driver": "Census, Compendium"},
        ])
        census_roles = _make_census_roles(tmp_path, [
            {"Gene": "RAF1", "Tier": "1", "Role in Cancer": "oncogene, fusion",
             "Tumour Types(Somatic)": "astrocytoma", "Tumour Types(Germline)": ""},
        ])
        comp_roles = _make_compendium_roles(tmp_path, [
            {"Gene": "RAF1", "ROLE": "Act", "CANCER_TYPE": "BLCA, CM"},
            {"Gene": "TP53", "ROLE": "LoF", "CANCER_TYPE": "BRCA"},
        ])
        r = _run(seq, combined, census_roles, comp_roles, tmp_path / "out")
        assert r.returncode == 0, r.stderr

        census = pd.read_csv(tmp_path / "out" / "census_driver.tsv", sep="\t", dtype=str)
        # BRAF filtered out (not in run)
        assert set(census["Protein_ID"]) == {"RAF1-201", "RAF1-262"}
        raf201 = census[census["Protein_ID"] == "RAF1-201"].iloc[0]
        assert raf201["Role in Cancer"] == "oncogene, fusion"
        assert raf201["Tier"] == "1"
        assert raf201["Gene"] == "RAF1"

    def test_compendium_carries_role_and_cancer_type(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201", "TP53-201"])
        combined = _make_combined(tmp_path, [
            {"Protein_ID": "RAF1-201", "Cancer Driver": "Census, Compendium"},
            {"Protein_ID": "TP53-201", "Cancer Driver": "Compendium"},
        ])
        census_roles = _make_census_roles(tmp_path, [
            {"Gene": "RAF1", "Tier": "1", "Role in Cancer": "oncogene",
             "Tumour Types(Somatic)": "", "Tumour Types(Germline)": ""},
        ])
        comp_roles = _make_compendium_roles(tmp_path, [
            {"Gene": "RAF1", "ROLE": "Act", "CANCER_TYPE": "BLCA, CM"},
            {"Gene": "TP53", "ROLE": "LoF", "CANCER_TYPE": "BRCA"},
        ])
        _run(seq, combined, census_roles, comp_roles, tmp_path / "out")
        comp = pd.read_csv(tmp_path / "out" / "compendium_driver.tsv", sep="\t", dtype=str)
        assert set(comp["Protein_ID"]) == {"RAF1-201", "TP53-201"}
        tp53 = comp[comp["Protein_ID"] == "TP53-201"].iloc[0]
        assert tp53["ROLE"] == "LoF"
        assert tp53["CANCER_TYPE"] == "BRCA"

    def test_combined_has_role_columns(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        combined = _make_combined(tmp_path, [
            {"Protein_ID": "RAF1-201", "Cancer Driver": "Census, Compendium"},
        ])
        census_roles = _make_census_roles(tmp_path, [
            {"Gene": "RAF1", "Tier": "1", "Role in Cancer": "oncogene, fusion",
             "Tumour Types(Somatic)": "", "Tumour Types(Germline)": ""},
        ])
        comp_roles = _make_compendium_roles(tmp_path, [
            {"Gene": "RAF1", "ROLE": "Act", "CANCER_TYPE": "BLCA"},
        ])
        _run(seq, combined, census_roles, comp_roles, tmp_path / "out")
        comb = pd.read_csv(tmp_path / "out" / "cancer_driver.tsv", sep="\t", dtype=str)
        assert list(comb.columns) == ["Protein_ID", "Cancer Driver", "Role in Cancer", "Compendium Role"]
        row = comb.iloc[0]
        assert row["Role in Cancer"] == "oncogene, fusion"
        assert row["Compendium Role"] == "Act"


class TestEdgeCases:
    def test_no_roles_still_produces_membership(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        combined = _make_combined(tmp_path, [
            {"Protein_ID": "RAF1-201", "Cancer Driver": "Census"},
        ])
        no_file = tmp_path / "NO_FILE"
        no_file.write_text("")
        r = _run(seq, combined, no_file, no_file, tmp_path / "out")
        assert r.returncode == 0, r.stderr
        census = pd.read_csv(tmp_path / "out" / "census_driver.tsv", sep="\t", dtype=str)
        assert set(census["Protein_ID"]) == {"RAF1-201"}

    def test_empty_combined_writes_empty_outputs(self, tmp_path):
        seq = _make_seq(tmp_path, ["RAF1-201"])
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        r = _run(seq, empty, tmp_path / "NO_FILE", tmp_path / "NO_FILE", tmp_path / "out")
        assert r.returncode == 0
        for name in ("cancer_driver.tsv", "census_driver.tsv", "compendium_driver.tsv"):
            df = pd.read_csv(tmp_path / "out" / name, sep="\t", dtype=str)
            assert len(df) == 0
