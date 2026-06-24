"""
tests/test_subset_fasta.py

Tests for bin/subset_fasta.py

Run from the project root:
    pytest tests/test_subset_fasta.py -v
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
BIN          = PROJECT_ROOT / "bin" / "subset_fasta.py"

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINI_FASTA_CONTENT = """\
>sp|P04049|RAF1_HUMAN RAF proto-oncogene serine/threonine-protein kinase OS=Homo sapiens GN=RAF1 PE=1 SV=1
MEHIQGAWKTISNGFGFKDAVFDGSSCISPTIVQQFGYQRRASDDGKLTDPSKTSNTIRVFLPNKQRTVVNVR
NGMSLHDCLMKALKVRGLQPECCAVFRLLHEHKGKKARLDWNTDAASLIGEELQVDFLDHVPLTTHNFARKAFQ
>sp|P15056|BRAF_HUMAN Serine/threonine-protein kinase B-raf OS=Homo sapiens GN=BRAF PE=1 SV=1
MAALSGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPP
SIYLEAYEEFNRSYGKPSTELEEKFNMDKFTAIKVSQSQRTGQLPYPHESWMKPLVQIDPAEEEDSTFRDLASL
>sp|P01116|RASK_HUMAN GTPase KRas OS=Homo sapiens GN=KRAS PE=1 SV=1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY
>sp|Q13077|TRAF1_HUMAN TNF receptor-associated factor 1 OS=Homo sapiens GN=TRAF1 PE=1 SV=1
MAENSSFEEVFHFKPQNPIFPQPQEPPQEPAQEPKASGERASGMRSRTPGSGEAAGPVEGTQPGPQPGPQPQP
"""

MINI_FASTA_GZ_CONTENT = MINI_FASTA_CONTENT  # same content, will be written as .gz

@pytest.fixture
def mini_fasta(tmp_path):
    p = tmp_path / "mini.fasta"
    p.write_text(MINI_FASTA_CONTENT)
    return p

@pytest.fixture
def mini_fasta_gz(tmp_path):
    p = tmp_path / "mini.fasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(MINI_FASTA_GZ_CONTENT)
    return p


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def run_subset(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(BIN)] + args,
        capture_output=True, text=True
    )


def count_headers(fasta_path: str) -> int:
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def get_headers(fasta_path: str) -> list[str]:
    with open(fasta_path) as fh:
        return [line.rstrip() for line in fh if line.startswith(">")]


# ---------------------------------------------------------------------------
# Tests: pass-through mode
# ---------------------------------------------------------------------------

class TestPassThrough:

    def test_no_search_copies_all(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out)])
        assert r.returncode == 0
        assert count_headers(out) == 4

    def test_empty_search_copies_all(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", ""])
        assert r.returncode == 0
        assert count_headers(out) == 4

    def test_sequence_lines_preserved(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        run_subset(["--input", str(mini_fasta), "--output", str(out)])
        with open(out) as fh:
            lines = [l.rstrip() for l in fh if l.strip() and not l.startswith(">")]
        with open(mini_fasta) as fh:
            expected = [l.rstrip() for l in fh if l.strip() and not l.startswith(">")]
        assert lines == expected


# ---------------------------------------------------------------------------
# Tests: substring search
# ---------------------------------------------------------------------------

class TestSubstringSearch:

    def test_exact_gene_name_match(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", "GN=RAF1 "])
        assert r.returncode == 0
        headers = get_headers(out)
        assert len(headers) == 1
        assert "RAF1_HUMAN" in headers[0]

    def test_entry_name_match(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", "BRAF_HUMAN"])
        assert r.returncode == 0
        assert count_headers(out) == 1

    def test_accession_match(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", "P04049"])
        assert r.returncode == 0
        assert count_headers(out) == 1
        assert "RAF1_HUMAN" in get_headers(out)[0]

    def test_broad_search_matches_multiple(self, mini_fasta, tmp_path):
        """'RAF1' matches RAF1_HUMAN, BRAF (no), RASK (no), TRAF1 (yes=contains RAF1)."""
        out = tmp_path / "out.fasta"
        run_subset(["--input", str(mini_fasta), "--output", str(out),
                    "--search", "RAF1"])
        headers = get_headers(out)
        # RAF1_HUMAN and TRAF1_HUMAN both contain 'RAF1'
        assert any("RAF1_HUMAN" in h for h in headers)
        assert any("TRAF1_HUMAN" in h for h in headers)
        # BRAF and RASK do not contain 'RAF1'
        assert not any("BRAF_HUMAN" in h for h in headers)
        assert not any("RASK_HUMAN" in h for h in headers)

    def test_no_match_exits_nonzero(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", "ZZZZNOTEXIST"])
        assert r.returncode != 0

    def test_case_sensitive(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta), "--output", str(out),
                        "--search", "raf1"])
        # Headers use uppercase; lowercase search should match nothing
        assert r.returncode != 0


# ---------------------------------------------------------------------------
# Tests: gzip support
# ---------------------------------------------------------------------------

class TestGzipInput:

    def test_reads_gzipped_fasta(self, mini_fasta_gz, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta_gz), "--output", str(out)])
        assert r.returncode == 0
        assert count_headers(out) == 4

    def test_subset_from_gz(self, mini_fasta_gz, tmp_path):
        out = tmp_path / "out.fasta"
        r = run_subset(["--input", str(mini_fasta_gz), "--output", str(out),
                        "--search", "P04049"])
        assert r.returncode == 0
        assert count_headers(out) == 1


# ---------------------------------------------------------------------------
# Tests: invert mode
# ---------------------------------------------------------------------------

class TestInvert:

    def test_invert_excludes_match(self, mini_fasta, tmp_path):
        out = tmp_path / "out.fasta"
        run_subset(["--input", str(mini_fasta), "--output", str(out),
                    "--search", "P04049", "--invert"])
        headers = get_headers(out)
        # The P04049 (RAF1_HUMAN) entry must not appear; other 3 entries must
        assert all("P04049" not in h for h in headers)
        assert count_headers(out) == 3
