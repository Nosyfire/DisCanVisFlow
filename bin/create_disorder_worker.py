#!/usr/bin/env python3
"""
create_disorder_worker.py — Module 5b: Disorder & Binding Prediction

Per-residue scores and combined disorder annotation.

Data sources
------------
  IUPred3        — long-disorder scores  (local library: iupred3_lib.py)
  ANCHOR2        — binding-region scores  (local library: iupred3_lib.py)
  AIUPred        — attention-based disorder (local library: aiupred_lib.py)
  AIUPred-Binding— attention-based binding  (local AIUPred/binding library)
  AlphaFold      — per-residue pLDDT via EBI API download

Strategy for local libraries
----------------------------
  1. Try direct import (works if the current Python has scipy + torch)
  2. Fallback: subprocess via --aiupred_python (default: /opt/anaconda3/envs/aiupred/bin/python)
     The aiupred conda env has both scipy and torch.

Combined disorder rule
----------------------
  Per residue:
    • experimental = MobiDB curated (from --mobidb_tsv)
    • If pLDDT available:
        valid = (pseudo-RSA > 0.582) = (pLDDT < 41.8)
    • Else:
        valid = (IUPred3 > 0.4) OR experimental
    • Domain exclusion: Pfam Domain → not disordered
  Region ≥ 5 contiguous residues (gap ≤ 5 allowed).

Inputs
------
  --loc_chrom        loc_chrom_with_names_isoforms_with_seq.tsv
  --mobidb_tsv       mobidb_disorders.tsv (optional; for combined rule)
  --pfam_tsv         pfam_domains.tsv (optional; for domain exclusion)
  --ext_programs     path to DisCanVis_Data_Process/External_Programs/
  --aiupred_python   Python binary with scipy+torch (default: /opt/anaconda3/envs/aiupred/bin/python)
  --output_dir       output directory
  --request_delay    seconds between AlphaFold API requests (default: 0.5)
  --skip_alphafold   skip AlphaFold API download
  --skip_iupred      skip IUPred3/ANCHOR2
  --skip_aiupred     skip AIUPred/AIUPred-Binding

Outputs
-------
  iupred_scores.tsv          — Entry_Name | IUPredscores (comma-sep)
  anchor_scores.tsv          — Entry_Name | AnchorScore
  aiupred_scores.tsv         — Entry_Name | AIUPredscores
  aiupred_binding_scores.tsv — Entry_Name | AIUPredBinding
  alphafold_plddt.tsv        — Entry_Name | Plldtscores
  combined_disorder.tsv      — Entry_Name | Entry_Isoform | Gene | Start | End
  combined_disorder_pos.tsv  — Entry_Name | Position | CombinedDisorder (0/4)
"""

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
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

ALPHAFOLD_SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.pdb"
_ALPHAFOLD_FALLBACK_VERSIONS = (6, 5, 4)

# Default Python binary for subprocess-based predictions
DEFAULT_AIUPRED_PYTHON = "/opt/anaconda3/envs/aiupred/bin/python"


# ---------------------------------------------------------------------------
# Subprocess helpers — run IUPred3/AIUPred via a capable Python binary
# ---------------------------------------------------------------------------

def _subprocess_run(python_bin: str, script: str, input_text: str,
                    timeout: int = 600) -> str | None:
    """Execute a Python snippet via subprocess, passing input_text on stdin."""
    try:
        result = subprocess.run(
            [python_bin, "-c", script],
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode == 0:
            return result.stdout.strip()
        log.warning("Subprocess failed (rc=%d): %s", result.returncode, result.stderr[:600])
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        log.warning("Subprocess could not run (%s): %s", type(e).__name__, e)
        return None


def _filter_subprocess_output(raw: str | None) -> str | None:
    """Strip comment/warning lines (starting with '#') from subprocess stdout.

    Some prediction libraries print warnings like
    '# Warning: No GPU found...' to stdout instead of stderr,
    contaminating the parseable output.
    """
    if not raw:
        return None
    lines = [l for l in raw.split("\n") if l.strip() and not l.strip().startswith("#")]
    return "\n".join(lines) if lines else None


def _iupred3_subprocess(seq: str, python_bin: str, iupred3_dir: str
                        ) -> tuple[list[float], list[float]] | None:
    """Run IUPred3 + ANCHOR2 via subprocess. Returns (iupred_scores, anchor_scores)."""
    script = f"""
import sys
sys.path.insert(0, {repr(iupred3_dir)})
import iupred3_lib
seq = sys.stdin.read().strip()
iupred = iupred3_lib.iupred(seq)[0]
anchor = iupred3_lib.anchor2(seq)
print(",".join(f"{{v:.4f}}" for v in iupred))
print(",".join(f"{{v:.4f}}" for v in anchor))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, seq))
    if not out:
        return None
    lines = out.split("\n")
    if len(lines) < 2:
        return None
    try:
        iupred = [float(x) for x in lines[0].split(",") if x.strip()]
        anchor = [float(x) for x in lines[1].split(",") if x.strip()]
        return iupred, anchor
    except ValueError:
        log.warning("IUPred3: failed to parse output: %s", out[:200])
        return None


def _iupred3_batch_subprocess(seqs: dict, python_bin: str, iupred3_dir: str,
                               timeout: int = 3600) -> dict | None:
    """Run IUPred3 + ANCHOR2 for all sequences in one subprocess (import once).

    Returns {pid: {'iupred': [...], 'anchor': [...]}} or None on failure.
    """
    if not seqs:
        return {}
    import json as _json
    script = f"""
import sys, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, {repr(iupred3_dir)})
import iupred3_lib
data = json.loads(sys.stdin.read())
results = {{}}
for pid, seq in data.items():
    try:
        iupred_sc = iupred3_lib.iupred(seq)[0]
        anchor_sc = iupred3_lib.anchor2(seq)
        results[pid] = {{'iupred': [float(v) for v in iupred_sc],
                         'anchor': [float(v) for v in anchor_sc]}}
    except Exception as e:
        results[pid] = None
print(json.dumps(results))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, _json.dumps(seqs), timeout=timeout))
    if not out:
        return None
    try:
        return _json.loads(out)
    except Exception as e:
        log.warning("IUPred3 batch: failed to parse output: %s", str(e)[:200])
        return None


def _aiupred_disorder_subprocess(seq: str, python_bin: str, aiupred_dir: str
                                  ) -> list[float] | None:
    """Run AIUPred disorder via subprocess using aiupred-caid3 predict() API."""
    script = f"""
import sys
sys.path.insert(0, {repr(aiupred_dir)})
from aiupred_lib import init_models, predict, low_memory_predict
seq = sys.stdin.read().strip()
embedding, decoder, device = init_models('disorder')
if len(seq) < 1000:
    scores = predict(seq, embedding, decoder, device)
else:
    scores = low_memory_predict(seq, embedding, decoder, device)
print(",".join(f"{{v:.4f}}" for v in scores))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, seq))
    if not out:
        return None
    try:
        return [float(x) for x in out.split(",") if x.strip()]
    except ValueError:
        log.warning("AIUPred disorder: failed to parse output: %s", out[:200])
        return None


def _aiupred_disorder_batch_subprocess(seqs: dict, python_bin: str, aiupred_dir: str,
                                        timeout: int = 3600) -> dict | None:
    """Run AIUPred disorder for all sequences in one subprocess (load model once).

    Returns {pid: [scores]} or None on failure.
    """
    if not seqs:
        return {}
    import json as _json
    script = f"""
import sys, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, {repr(aiupred_dir)})
from aiupred_lib import init_models, predict, low_memory_predict
data = json.loads(sys.stdin.read())
embedding, decoder, device = init_models('disorder')
results = {{}}
for pid, seq in data.items():
    try:
        if len(seq) < 1000:
            sc = predict(seq, embedding, decoder, device)
        else:
            sc = low_memory_predict(seq, embedding, decoder, device)
        results[pid] = [float(v) for v in sc]
    except Exception as e:
        results[pid] = None
print(json.dumps(results))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, _json.dumps(seqs), timeout=timeout))
    if not out:
        return None
    try:
        return _json.loads(out)
    except Exception as e:
        log.warning("AIUPred disorder batch: failed to parse output: %s", str(e)[:200])
        return None


def _aiupred_binding_subprocess(seq: str, python_bin: str, aiupred_binding_dir: str
                                 ) -> list[float] | None:
    """Run AIUPred-Binding via subprocess."""
    script = f"""
import sys
sys.path.insert(0, {repr(aiupred_binding_dir)})
from aiupred_lib import init_models, predict_binding, low_memory_predict_binding
seq = sys.stdin.read().strip()
models = init_models("binding")
embedding, decoder, device = models
if len(seq) < 1000:
    sc = predict_binding(seq, embedding, decoder, device)
else:
    sc = low_memory_predict_binding(seq, embedding, decoder, device)
print(",".join(f"{{v:.4f}}" for v in sc))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, seq, timeout=900))
    if not out:
        return None
    try:
        return [float(x) for x in out.split(",") if x.strip()]
    except ValueError:
        log.warning("AIUPred binding: failed to parse output: %s", out[:200])
        return None


def _aiupred_binding_batch_subprocess(seqs: dict, python_bin: str, aiupred_binding_dir: str,
                                       timeout: int = 3600) -> dict | None:
    """Run AIUPred-Binding for all sequences in one subprocess (load model once).

    Returns {pid: [scores]} or None on failure.
    """
    if not seqs:
        return {}
    import json as _json
    script = f"""
import sys, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, {repr(aiupred_binding_dir)})
from aiupred_lib import init_models, predict_binding, low_memory_predict_binding
data = json.loads(sys.stdin.read())
embedding, decoder, device = init_models('binding')
results = {{}}
for pid, seq in data.items():
    try:
        if len(seq) < 1000:
            sc = predict_binding(seq, embedding, decoder, device)
        else:
            sc = low_memory_predict_binding(seq, embedding, decoder, device)
        results[pid] = [float(v) for v in sc]
    except Exception as e:
        results[pid] = None
print(json.dumps(results))
"""
    out = _filter_subprocess_output(_subprocess_run(python_bin, script, _json.dumps(seqs), timeout=timeout))
    if not out:
        return None
    try:
        return _json.loads(out)
    except Exception as e:
        log.warning("AIUPred binding batch: failed to parse output: %s", str(e)[:200])
        return None


# ---------------------------------------------------------------------------
# Direct-import helpers (used when scipy/torch available in current env)
# ---------------------------------------------------------------------------

_iupred3_lib = None

def _get_iupred3_direct(ext_programs: str):
    global _iupred3_lib
    if _iupred3_lib is not None:
        return _iupred3_lib
    iupred_path = str(Path(ext_programs) / "iupred3")
    if iupred_path not in sys.path:
        sys.path.insert(0, iupred_path)
    try:
        import iupred3_lib as lib
        _iupred3_lib = lib
        return lib
    except ImportError:
        return None


_aiupred_lib      = None
_aiupred_bind_lib = None

def _get_aiupred_bind_direct(ext_programs: str):
    """Try to directly import AIUPred-Binding (init_models + predict_binding)."""
    global _aiupred_bind_lib
    if _aiupred_bind_lib is not None:
        return _aiupred_bind_lib
    for subdir in ["AIUPred", "aiupred-caid3", "aiupred"]:
        apath = str(Path(ext_programs) / subdir)
        if apath not in sys.path:
            sys.path.insert(0, apath)
        try:
            import aiupred_lib as lib
            if hasattr(lib, 'init_models') and hasattr(lib, 'predict_binding'):
                _aiupred_bind_lib = lib
                return lib
        except ImportError:
            pass
    return None


def _get_aiupred_direct(ext_programs: str):
    """Try to directly import aiupred-caid3 predict() API (disorder)."""
    global _aiupred_lib
    if _aiupred_lib is not None:
        return _aiupred_lib
    # Prefer aiupred-caid3 which has init_models/predict API used for disorder
    for subdir in ["aiupred-caid3", "aiupred", "AIUPred"]:
        apath = str(Path(ext_programs) / subdir)
        if apath not in sys.path:
            sys.path.insert(0, apath)
        try:
            import aiupred_lib as lib
            # Verify it has the caid3-style predict() function (not just aiupred_disorder)
            if hasattr(lib, 'init_models') and hasattr(lib, 'predict'):
                _aiupred_lib = lib
                return lib
        except ImportError:
            pass
    return None


# ---------------------------------------------------------------------------
# AlphaFold pLDDT download
# ---------------------------------------------------------------------------

def _get_http(url: str, delay: float, timeout: int = 60):
    if not _HAS_REQUESTS:
        return None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(delay * (attempt + 1))
    return None


def _parse_pdb_plddt(text: str) -> list[float]:
    seen: set[int] = set()
    scores: list[float] = []
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            res_num  = int(line[22:26].strip())
            b_factor = float(line[60:66].strip())
        except (ValueError, IndexError):
            continue
        if res_num not in seen:
            seen.add(res_num)
            scores.append(b_factor)
    return scores


def fetch_alphafold_plddt(acc: str, delay: float) -> list[float] | None:
    # Strip isoform suffix — AlphaFold only has canonical entries (e.g. P04049, not P04049-2)
    acc_base = acc.split("-")[0] if "-" in acc else acc

    # Try summary API first to get the current pdbUrl (version-agnostic)
    if _HAS_REQUESTS:
        summary_r = _get_http(ALPHAFOLD_SUMMARY_URL.format(acc=acc_base), delay=0, timeout=15)
        if summary_r is not None:
            try:
                import json as _json
                data = _json.loads(summary_r.text)
                pdb_url = data[0].get("pdbUrl") if isinstance(data, list) and data else None
                if pdb_url:
                    r = _get_http(pdb_url, delay, timeout=60)
                    if r is not None:
                        scores = _parse_pdb_plddt(r.text)
                        return scores if scores else None
            except Exception:
                pass

    # Fallback: try known versioned URLs from newest to oldest
    for ver in _ALPHAFOLD_FALLBACK_VERSIONS:
        url = ALPHAFOLD_PDB_URL.format(acc=acc_base, ver=ver)
        r = _get_http(url, delay, timeout=60)
        if r is not None:
            scores = _parse_pdb_plddt(r.text)
            return scores if scores else None
    return None


# ---------------------------------------------------------------------------
# Score computation (tries direct import, falls back to subprocess)
# ---------------------------------------------------------------------------

def compute_scores(
    proteins:       pd.DataFrame,
    id_col:         str,
    ext_programs:   str,
    aiupred_python: str,
    delay:          float,
    skip_iupred:    bool,
    skip_aiupred:   bool,
    skip_alphafold: bool,
    plddt_local:    dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Returns (iupred_df, anchor_df, aiupred_df, aiupred_bind_df, plddt_df).

    plddt_local: optional dict {canonical_acc → [pLDDT scores]}.  When provided,
    AlphaFold API calls are replaced by dict lookups — no network needed.
    """

    ext = Path(ext_programs)
    iupred3_dir      = str(ext / "iupred3")
    aiupred_dir      = str(ext / "aiupred-caid3")   # disorder: caid3 API (predict())
    aiupred_bind_dir = str(ext / "AIUPred")          # binding: AIUPred API (predict_binding())

    # Probe direct import
    _direct_iupred  = None if skip_iupred  else _get_iupred3_direct(ext_programs)
    _direct_aiupred = None if skip_aiupred else _get_aiupred_direct(ext_programs)

    use_subprocess_iupred  = (not skip_iupred)  and (_direct_iupred  is None)
    use_subprocess_aiupred = (not skip_aiupred) and (_direct_aiupred is None)

    if use_subprocess_iupred or use_subprocess_aiupred:
        log.info("Direct library import failed — using subprocess via %s", aiupred_python)

    acc_col = "Entry_Isoform" if "Entry_Isoform" in proteins.columns else None

    # ── Batch subprocess pre-computation (load model once for all sequences) ──
    # When direct import is unavailable, run one subprocess per predictor that
    # loads the model once and scores every sequence, avoiding per-protein
    # subprocess startup + model-loading overhead (~10-13s → ~0.03s per protein).
    import json as _json

    def _valid_seqs(df):
        return {str(r.get(id_col, "")): str(r.get("Sequence", ""))
                for _, r in df.iterrows()
                if str(r.get("Sequence", "")) not in ("", "nan") and len(str(r.get("Sequence", ""))) >= 2}

    _batch_iupred   = {}   # pid → {'iupred': [...], 'anchor': [...]}
    _batch_aiu_dis  = {}   # pid → [scores]
    _batch_aiu_bind = {}   # pid → [scores]

    if use_subprocess_iupred and aiupred_python:
        seqs = _valid_seqs(proteins)
        log.info("Batch IUPred3+ANCHOR2 for %d sequences …", len(seqs))
        _res = _iupred3_batch_subprocess(seqs, aiupred_python, iupred3_dir)
        if _res is not None:
            _batch_iupred = _res
            n_ok = sum(1 for v in _res.values() if v)
            log.info("Batch IUPred3: %d/%d succeeded", n_ok, len(seqs))
        else:
            log.warning("Batch IUPred3 failed — falling back to per-protein subprocess")

    if use_subprocess_aiupred and aiupred_python:
        seqs = _valid_seqs(proteins)
        log.info("Batch AIUPred disorder for %d sequences …", len(seqs))
        _res = _aiupred_disorder_batch_subprocess(seqs, aiupred_python, aiupred_dir)
        if _res is not None:
            _batch_aiu_dis = _res
            n_ok = sum(1 for v in _res.values() if v)
            log.info("Batch AIUPred disorder: %d/%d succeeded", n_ok, len(seqs))
        else:
            log.warning("Batch AIUPred disorder failed — falling back to per-protein subprocess")

        log.info("Batch AIUPred-Binding for %d sequences …", len(seqs))
        _res = _aiupred_binding_batch_subprocess(seqs, aiupred_python, aiupred_bind_dir)
        if _res is not None:
            _batch_aiu_bind = _res
            n_ok = sum(1 for v in _res.values() if v)
            log.info("Batch AIUPred-Binding: %d/%d succeeded", n_ok, len(seqs))
        else:
            log.warning("Batch AIUPred-Binding failed — falling back to per-protein subprocess")

    iupred_rows, anchor_rows, aiupred_rows, aiu_bind_rows, plddt_rows = [], [], [], [], []

    try:
        from tqdm import tqdm as _tqdm
        _prot_iter = _tqdm(proteins.iterrows(), total=len(proteins), desc='Disorder scores', unit='isoform')
    except ImportError:
        _prot_iter = proteins.iterrows()

    for _, row in _prot_iter:
        pid = str(row.get(id_col, ""))
        seq = str(row.get("Sequence", ""))
        acc = str(row.get(acc_col, "")) if acc_col else ""

        if not seq or seq in ("nan", "") or len(seq) < 2:
            continue

        # ── IUPred3 + ANCHOR2 ─────────────────────────────────────────────────
        if not skip_iupred:
            iupred_sc, anchor_sc = None, None
            if _direct_iupred:
                try:
                    iupred_sc = _direct_iupred.iupred(seq)[0]
                    anchor_sc = _direct_iupred.anchor2(seq)
                except Exception as e:
                    log.debug("IUPred3 direct error (%s): %s", pid, e)

            if iupred_sc is None:
                # Prefer batch result; fall back to per-protein subprocess if batch missed this pid
                batch_result = _batch_iupred.get(pid)
                if batch_result:
                    iupred_sc = batch_result.get("iupred")
                    anchor_sc = batch_result.get("anchor")
                elif use_subprocess_iupred:
                    result = _iupred3_subprocess(seq, aiupred_python, iupred3_dir)
                    if result:
                        iupred_sc, anchor_sc = result

            if iupred_sc:
                iupred_rows.append({id_col: pid,
                    "IUPredscores": ", ".join(f"{v:.4f}" for v in iupred_sc)})
            if anchor_sc:
                anchor_rows.append({id_col: pid,
                    "AnchorScore": ", ".join(f"{v:.4f}" for v in anchor_sc)})

        # ── AIUPred disorder ──────────────────────────────────────────────────
        if not skip_aiupred:
            aiu_sc = None
            if _direct_aiupred:
                try:
                    if hasattr(_direct_aiupred, 'init_models') and hasattr(_direct_aiupred, 'predict'):
                        # aiupred-caid3 API: init_models once, then predict per sequence
                        if not hasattr(compute_scores, '_aiupred_models'):
                            compute_scores._aiupred_models = _direct_aiupred.init_models('disorder')
                        emb, dec, dev = compute_scores._aiupred_models
                        fn = (_direct_aiupred.predict if len(seq) < 1000
                              else _direct_aiupred.low_memory_predict)
                        aiu_sc = fn(seq, emb, dec, dev)
                    else:
                        aiu_sc = _direct_aiupred.aiupred_disorder(seq)[0]
                except Exception as e:
                    log.debug("AIUPred direct error (%s): %s", pid, e)

            if aiu_sc is None:
                # Prefer batch result; fall back to per-protein subprocess if batch missed this pid
                aiu_sc = _batch_aiu_dis.get(pid) or None
                if aiu_sc is None and use_subprocess_aiupred:
                    aiu_sc = _aiupred_disorder_subprocess(seq, aiupred_python, aiupred_dir)

            if aiu_sc:
                aiupred_rows.append({id_col: pid,
                    "AIUPredscores": ", ".join(f"{v:.4f}" for v in aiu_sc)})

        # ── AIUPred-Binding ───────────────────────────────────────────────────
        if not skip_aiupred:
            bind_sc = None
            _direct_aiupred_bind = _get_aiupred_bind_direct(aiupred_bind_dir
                                   if aiupred_bind_dir else ext_programs)
            if _direct_aiupred_bind:
                try:
                    if not hasattr(compute_scores, '_aiupred_bind_models'):
                        compute_scores._aiupred_bind_models = (
                            _direct_aiupred_bind.init_models('binding'))
                    emb_b, dec_b, dev_b = compute_scores._aiupred_bind_models
                    fn_b = (_direct_aiupred_bind.predict_binding if len(seq) < 1000
                            else _direct_aiupred_bind.low_memory_predict_binding)
                    bind_sc = fn_b(seq, emb_b, dec_b, dev_b)
                except Exception as e:
                    log.debug("AIUPred-Binding direct error (%s): %s", pid, e)
            if bind_sc is None:
                # Prefer batch result; fall back to per-protein subprocess if batch missed this pid
                bind_sc = _batch_aiu_bind.get(pid) or None
                if bind_sc is None and aiupred_python:
                    bind_sc = _aiupred_binding_subprocess(seq, aiupred_python, aiupred_bind_dir)
            if bind_sc:
                aiu_bind_rows.append({id_col: pid,
                    "AIUPredBinding": ", ".join(f"{v:.4f}" for v in bind_sc)})

        # ── AlphaFold pLDDT ───────────────────────────────────────────────────
        if not skip_alphafold and acc and acc not in ("nan", ""):
            acc_base = acc.split("-")[0] if "-" in acc else acc
            if plddt_local is not None:
                # Fast path: dict lookup from pre-extracted bulk tar
                plddt = plddt_local.get(acc_base)
            else:
                # Slow path: per-protein EBI API call (not recommended at full proteome scale)
                plddt = fetch_alphafold_plddt(acc, delay)
                time.sleep(delay)
            if plddt:
                if len(plddt) > len(seq):
                    plddt = plddt[:len(seq)]
                elif len(plddt) < len(seq):
                    plddt = plddt + [0.0] * (len(seq) - len(plddt))
                plddt_rows.append({id_col: pid,
                    "Plldtscores": ", ".join(f"{v:.1f}" for v in plddt)})

    def _df(rows, col):
        return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[id_col, col])

    return (
        _df(iupred_rows,    "IUPredscores"),
        _df(anchor_rows,    "AnchorScore"),
        _df(aiupred_rows,   "AIUPredscores"),
        _df(aiu_bind_rows,  "AIUPredBinding"),
        _df(plddt_rows,     "Plldtscores"),
    )


# ---------------------------------------------------------------------------
# Per-position combined disorder table
# ---------------------------------------------------------------------------

def _parse_scores(s: str) -> list[float]:
    if not s or str(s) in ("nan", ""):
        return []
    try:
        return [float(x.strip()) for x in str(s).split(",") if x.strip()]
    except ValueError:
        return []


def build_per_position(
    loc_df:    pd.DataFrame,
    id_col:    str,
    iupred_df: pd.DataFrame,
    plddt_df:  pd.DataFrame,
    mobidb_df: pd.DataFrame | None,
    pfam_df:   pd.DataFrame | None,
) -> pd.DataFrame:
    # Run on ALL isoforms (not just main_isoform) — matching legacy behaviour
    lc = loc_df.copy()

    iupred_map = (iupred_df.set_index(id_col)["IUPredscores"].to_dict()
                  if not iupred_df.empty and id_col in iupred_df.columns else {})
    plddt_map  = (plddt_df.set_index(id_col)["Plldtscores"].to_dict()
                  if not plddt_df.empty and id_col in plddt_df.columns else {})

    # MobiDB disorder positions {id → set of 1-based positions}
    mob_pos: dict[str, set] = {}
    if mobidb_df is not None and not mobidb_df.empty:
        pid_col = next((c for c in ["Protein_ID", "Entry_Name", id_col]
                        if c in mobidb_df.columns), None)
        if pid_col:
            for _, mr in mobidb_df.iterrows():
                pid = str(mr.get(pid_col, ""))
                try:
                    s = int(float(mr.get("Start", 0)))
                    e = int(float(mr.get("End",   0)))
                    mob_pos.setdefault(pid, set()).update(range(s, e + 1))
                except (ValueError, TypeError):
                    pass

    # Pfam Domain positions {id → set of 1-based positions}
    pfam_pos: dict[str, set] = {}
    if pfam_df is not None and not pfam_df.empty:
        pid_col = next((c for c in ["Protein_ID", "Accession", "Entry_Isoform"]
                        if c in pfam_df.columns), None)
        if pid_col:
            for _, pr in pfam_df.iterrows():
                if "Domain" not in str(pr.get("type", "")):
                    continue
                pid = str(pr.get(pid_col, ""))
                try:
                    s = int(float(pr.get("envelope_start", pr.get("Start", 0))))
                    e = int(float(pr.get("envelope_end",   pr.get("End",   0))))
                    pfam_pos.setdefault(pid, set()).update(range(s, e + 1))
                except (ValueError, TypeError):
                    pass

    rows = []
    for _, row in lc.iterrows():
        pid = str(row.get(id_col, ""))
        seq = str(row.get("Sequence", ""))
        if not seq or seq in ("nan", ""):
            continue

        iupred_lst = _parse_scores(iupred_map.get(pid, ""))
        plddt_lst  = _parse_scores(plddt_map.get(pid, ""))
        mob_set    = mob_pos.get(pid, set())
        pfam_set   = pfam_pos.get(pid, set())

        for i, _ in enumerate(seq):
            pos    = i + 1
            iupred = iupred_lst[i] if i < len(iupred_lst) else np.nan
            plddt  = plddt_lst[i]  if i < len(plddt_lst)  else np.nan
            rsa    = (100.0 - plddt) / 100.0 if not np.isnan(plddt) else np.nan
            mobidb = 1 if pos in mob_set else np.nan
            in_pfam = pos in pfam_set

            pfam_rule = not in_pfam
            if not np.isnan(rsa):
                valid = (rsa > 0.582) and pfam_rule
            else:
                valid = ((not np.isnan(iupred) and iupred > 0.4)
                         or (not np.isnan(mobidb))) and pfam_rule

            rows.append({
                id_col:              pid,
                "Position":          pos,
                "IUPredscores":      iupred,
                "RSA":               rsa,
                "MobiDB":            mobidb,
                "Pfam_Info":         "Domain" if in_pfam else np.nan,
                "ValidRegionStart":  valid,
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=[id_col, "Position", "IUPredscores", "RSA",
                 "MobiDB", "Pfam_Info", "ValidRegionStart"])


def find_regions(group: pd.DataFrame,
                 min_len: int = 5, allowed_gap: int = 5) -> list[tuple]:
    group = group.sort_values("Position").reset_index(drop=True)
    if group.empty:
        return []
    regions = []
    start = prev = group["Position"].iloc[0]
    for i in range(1, len(group)):
        cur = group["Position"].iloc[i]
        if cur - prev > allowed_gap:
            if prev - start + 1 >= min_len:
                regions.append((start, prev))
            start = cur
        prev = cur
    end = group["Position"].iloc[-1]
    if end - start + 1 >= min_len:
        regions.append((start, end))
    return regions


def compute_combined_disorder(
    per_pos_df: pd.DataFrame,
    loc_df:     pd.DataFrame,
    id_col:     str,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    if "ValidRegionStart" not in per_pos_df.columns or per_pos_df.empty:
        ec = pd.DataFrame(columns=[id_col, "Entry_Isoform", "Gene", "Start", "End"])
        ep = pd.DataFrame(columns=[id_col, "Position", "CombinedDisorder"])
        return ec, ep

    gene_col = next((c for c in ["Gene", "Gene_Gencode", "Gene_Uniprot"]
                     if c in loc_df.columns), None)
    meta_cols = [id_col, "Entry_Isoform"] + ([gene_col] if gene_col else [])
    meta = (loc_df[meta_cols].drop_duplicates()
            .rename(columns={gene_col: "Gene"} if gene_col else {}))

    valid_df = per_pos_df[per_pos_df["ValidRegionStart"]]
    region_rows, dis_rows = [], []

    # Use id_col that was embedded as column name in per_pos_df
    _pid_col = id_col if id_col in per_pos_df.columns else "Protein_ID"

    for pid, group in valid_df.groupby(_pid_col):
        gi = meta[meta[id_col] == pid]
        acc  = gi["Entry_Isoform"].iloc[0] if not gi.empty else ""
        gene = gi["Gene"].iloc[0]           if (not gi.empty and "Gene" in gi.columns) else ""
        for s, e in find_regions(group):
            region_rows.append({id_col: pid, "Entry_Isoform": acc,
                                 "Gene": gene, "Start": s, "End": e})
            dis_rows.extend({_pid_col: pid, "Position": p, "CombinedDisorder": 1}
                            for p in range(s, e + 1))

    region_df = pd.DataFrame(region_rows) if region_rows else pd.DataFrame(
        columns=[id_col, "Entry_Isoform", "Gene", "Start", "End"])

    all_pos = per_pos_df[[_pid_col, "Position"]].copy()
    if dis_rows:
        dis_df  = pd.DataFrame(dis_rows)
        pos_df  = all_pos.merge(dis_df, on=[_pid_col, "Position"], how="left")
    else:
        pos_df  = all_pos.copy()
        pos_df["CombinedDisorder"] = 0

    pos_df["CombinedDisorder"] = pos_df["CombinedDisorder"].fillna(0).astype(int)
    return region_df, pos_df[[_pid_col, "Position", "CombinedDisorder"]]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Module 5b: disorder + binding prediction")
    p.add_argument("--loc_chrom",      required=True)
    p.add_argument("--mobidb_tsv",     default=None)
    p.add_argument("--pfam_tsv",       default=None)
    p.add_argument("--ext_programs",   required=True,
                   help="Path to DisCanVis_Data_Process/External_Programs/")
    p.add_argument("--aiupred_python", default=DEFAULT_AIUPRED_PYTHON,
                   help="Python binary with scipy+torch for subprocess fallback")
    p.add_argument("--output_dir",     default=".")
    p.add_argument("--request_delay",  type=float, default=0.5)
    p.add_argument("--skip_alphafold",      action="store_true", default=False)
    p.add_argument("--skip_iupred",         action="store_true", default=False)
    p.add_argument("--skip_aiupred",        action="store_true", default=False)
    p.add_argument("--alphafold_plddt_tsv", default=None,
                   help="Pre-extracted pLDDT TSV (Accession, Plldtscores) — replaces EBI API calls")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading loc_chrom…")
    loc_df = pd.read_csv(args.loc_chrom, sep="\t", dtype=str)
    loc_df["Sequence"] = loc_df["Sequence"].fillna("")

    id_col = next((c for c in ["Protein_ID", "Entry_Name", "transcript_stable_id"]
                   if c in loc_df.columns), "Entry_Name")
    proteins = loc_df.copy()

    log.info("Computing disorder/binding scores for %d proteins (all isoforms, key=%s)…",
             len(proteins), id_col)

    # Load pre-extracted AlphaFold pLDDT TSV when provided (replaces per-protein API)
    plddt_local: dict | None = None
    _plddt_path = Path(args.alphafold_plddt_tsv) if args.alphafold_plddt_tsv else None
    _plddt_ok   = (_plddt_path is not None
                   and _plddt_path.exists()
                   and _plddt_path.name != "NO_FILE"
                   and _plddt_path.stat().st_size > 0)
    if _plddt_ok:
        af_df = pd.read_csv(args.alphafold_plddt_tsv, sep="\t", dtype=str)
        if not af_df.empty and "Accession" in af_df.columns and "Plldtscores" in af_df.columns:
            plddt_local = {}
            for _, row in af_df.iterrows():
                try:
                    scores = [float(x) for x in str(row["Plldtscores"]).split(",") if x.strip()]
                    plddt_local[str(row["Accession"])] = scores
                except (ValueError, AttributeError):
                    pass
            log.info("AlphaFold (local): loaded pLDDT for %d proteins (no API calls)",
                     len(plddt_local))

    iupred_df, anchor_df, aiupred_df, aiu_bind_df, plddt_df = compute_scores(
        proteins, id_col, args.ext_programs, args.aiupred_python,
        args.request_delay,
        args.skip_iupred, args.skip_aiupred, args.skip_alphafold,
        plddt_local=plddt_local)

    # Legacy output file names
    iupred_df.to_csv(outdir / "IUPredscores.tsv",      sep="\t", index=False)
    anchor_df.to_csv(outdir / "Anchorscores.tsv",       sep="\t", index=False)
    aiupred_df.to_csv(outdir / "AIUPredscores.tsv",     sep="\t", index=False)
    aiu_bind_df.to_csv(outdir / "AIUPredBinding.tsv",   sep="\t", index=False)
    plddt_df.to_csv(outdir / "AlphaFoldTable.tsv",      sep="\t", index=False)

    mobidb_df = None
    if args.mobidb_tsv and Path(args.mobidb_tsv).exists():
        mobidb_df = pd.read_csv(args.mobidb_tsv, sep="\t", dtype=str)

    pfam_df = None
    if args.pfam_tsv and Path(args.pfam_tsv).exists():
        pfam_df = pd.read_csv(args.pfam_tsv, sep="\t", dtype=str)

    log.info("Building per-position disorder table…")
    per_pos = build_per_position(loc_df, id_col, iupred_df, plddt_df, mobidb_df, pfam_df)
    region_df, pos_df = compute_combined_disorder(per_pos, loc_df, id_col)

    # Legacy output file names (CombinedDisorderNew = regions, CombinedDisorderNew_Pos = per-pos)
    region_df.to_csv(outdir / "CombinedDisorderNew.tsv",     sep="\t", index=False)
    pos_df.to_csv(outdir / "CombinedDisorderNew_Pos.tsv",    sep="\t", index=False)

    log.info(
        "Done — IUPred=%d  ANCHOR2=%d  AIUPred=%d  AIUPred-Binding=%d  pLDDT=%d  disorder_regions=%d",
        len(iupred_df), len(anchor_df), len(aiupred_df), len(aiu_bind_df),
        len(plddt_df), len(region_df),
    )


if __name__ == "__main__":
    main()
