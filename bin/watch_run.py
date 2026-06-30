#!/usr/bin/env python3
"""
bin/watch_run.py — Live progress monitor for a running DisCanVisFlow pipeline.

Shows all pipeline phases with task-by-task status, scatter-chunk progress,
BLAT chunk detail, and per-phase ETA estimates.

Data sources (in priority order):
  1. trace TSV  (results/<project>/reports/trace.tsv) — written by the REAL run,
     never overwritten by stub/test runs → used for completed/failed tasks
  2. work dir scan — detect currently-running tasks (dirs with .command.begin
     but no .exitcode)
  3. .nextflow.log — legacy fallback when no trace is found

Usage:
  python bin/watch_run.py                          # auto-detect trace, refresh every 30s
  python bin/watch_run.py --project discanvis      # explicit project
  python bin/watch_run.py --trace results/discanvis/reports/trace.tsv
  python bin/watch_run.py --interval 10            # refresh every 10s
  python bin/watch_run.py --once                   # print once and exit
  python bin/watch_run.py --work work/local
"""

import argparse
import csv
import re
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

# ── Terminal helpers ──────────────────────────────────────────────────────────

CLEAR  = "\033[2J\033[H"
BOLD   = "\033[1m"
DIM    = "\033[2m"
GREEN  = "\033[32m"
CYAN   = "\033[36m"
YELLOW = "\033[33m"
RED    = "\033[31m"
RESET  = "\033[0m"

def bar(frac: float, width: int = 28) -> str:
    frac = max(0.0, min(1.0, frac))
    filled = int(frac * width)
    return "█" * filled + "░" * (width - filled)

def human_time(seconds: float) -> str:
    if seconds < 0:
        return "?"
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        m, s = divmod(int(seconds), 60)
        return f"{m}m {s:02d}s"
    h, rem = divmod(int(seconds), 3600)
    m = rem // 60
    return f"{h}h {m:02d}m"

def elapsed_since(path: Path) -> str:
    try:
        s = time.time() - path.stat().st_mtime
        return human_time(s)
    except OSError:
        return "?"

# ── Phase definitions ─────────────────────────────────────────────────────────

PHASES = [
    ("Setup",          ["SUBSET_", "MAKEBLASTDB", "SPLIT_CDNA", "MERGE_UNIPROT",
                        "FETCH_", "DECOMPRESS_", "PPI_PREPROCESS", "ELM_CLASS_MAP",
                        "PARSE_UNIPROT", "PARSE_ALPHAFOLD", "NORMALISE_DEPMAP"]),
    ("BLAST + BLAT",   ["BLASTP", "MERGE_BLAST", "ID_MAP",
                        "BLAT_ALIGN", "MERGE_BLAT"]),
    ("Seq / Genome",   ["SEQUENCE_PROCESS", "SPLIT_SEQ_TABLE", "GENOME_MAP",
                        "GENOME_QUERY_MAP"]),
    ("Mutation",       ["MUTATION_MAP", "CLINVAR_DISEASE"]),
    ("Annotation",     ["ANNOTATION_MAP", "MOBIDB_MAP", "DISORDER_MAP",
                        "PDB_MAP", "PDB_BULK_MAP", "EXON_MAP", "GO_MAP",
                        "POLYMORPHISM_MAP", "PEM_MAP", "PEM_TRANSFER",
                        "COILEDCOILS_MAP", "PPI_MAP", "CONSERVATION_MAP",
                        "SCANSITE_MAP", "PROTEINGYM_MAP", "MAVEDB_MAP",
                        "DEPMAP_MAP", "DBNSFP_MAP", "PATHOGENICITY_MAP",
                        "OMIM_MAP", "ALPHAMISSENSE_MAP", "ELM_SWITCHES_MAP",
                        "FINCHES_MAP", "POSITION_BASED_MAP", "ISOFORM_ALIGN_MAP"]),
    ("Finalize",       ["TRANSCRIPT_MAP", "HOMOLOGY_MANIFEST", "MAPPING_REPORT"]),
]

def phase_of(name: str) -> str:
    for phase_name, keywords in PHASES:
        if any(kw in name.upper() for kw in keywords):
            return phase_name
    return "Other"

# ── Trace TSV parsing (primary source) ───────────────────────────────────────

def find_trace(project: str | None, outdir: str) -> Path | None:
    """Find the most recently modified trace.tsv under results/."""
    candidates = []
    if project:
        p = Path(outdir) / project / "reports" / "trace.tsv"
        if p.exists():
            return p
    # auto-detect: find all trace.tsv files and return newest
    for p in Path(outdir).glob("*/reports/trace.tsv"):
        candidates.append(p)
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)

def parse_trace(trace_path: Path) -> dict:
    """Parse trace.tsv → dict of task_id → task info for COMPLETED/FAILED tasks."""
    tasks = {}
    if not trace_path or not trace_path.exists():
        return tasks
    with open(trace_path, newline="", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name   = row.get("name", "").strip()
            status = row.get("status", "").strip()
            tid    = row.get("task_id", name)
            if not name:
                continue
            # Nextflow writes CACHED for -resume'd tasks; treat as COMPLETED
            if status == "CACHED":
                status = "COMPLETED"
            tasks[tid] = {
                "name":     name,
                "status":   status,
                "exit":     row.get("exit", "?").strip(),
                "phase":    phase_of(name),
                "duration": row.get("realtime", "").strip(),
            }
    return tasks

# ── Work-dir scan for RUNNING tasks ──────────────────────────────────────────

def find_running_in_workdir(work_dir: Path) -> list[dict]:
    """
    Scan work dir for tasks that have started (.command.begin) but not finished
    (.exitcode absent).  Returns list of {name, work_dir, elapsed_s}.
    """
    running = []
    if not work_dir.exists():
        return running
    for top in work_dir.iterdir():
        if not top.is_dir() or len(top.name) != 2:
            continue
        for sub in top.iterdir():
            if not sub.is_dir():
                continue
            begin = sub / ".command.begin"
            exitc = sub / ".exitcode"
            cmd_sh = sub / ".command.sh"
            if not begin.exists() or exitc.exists() or not cmd_sh.exists():
                continue
            # extract task name from .command.sh header comment: # NXF_TASK=NAME
            name = _task_name_from_sh(cmd_sh)
            if name is None:
                continue
            elapsed_s = time.time() - begin.stat().st_mtime
            running.append({
                "name":      name,
                "status":    "RUNNING",
                "exit":      "-",
                "phase":     phase_of(name),
                "elapsed_s": elapsed_s,
                "work_dir":  sub,
            })
    return running

_NAME_PATTERNS = [
    re.compile(r"#\s*NXF_TASK\s*=\s*(.+)"),           # newer Nextflow
    re.compile(r"nxf_main\(\)\s*\{.*?#\s*(.+?)\s*$",  # fallback
               re.DOTALL | re.MULTILINE),
]

def _task_name_from_sh(cmd_sh: Path) -> str | None:
    try:
        # Read only first 60 lines for speed
        lines = []
        with open(cmd_sh, errors="replace") as f:
            for i, line in enumerate(f):
                if i >= 60:
                    break
                lines.append(line)
        text = "".join(lines)
    except OSError:
        return None

    # Pattern: the process tag line looks like: # TASK: COILEDCOILS_MAP (coiledcoils chunk_003)
    m = re.search(r"#\s*TASK:\s*(.+)", text)
    if m:
        return m.group(1).strip()

    # Nextflow writes the process name as a comment near the top
    m = re.search(r"#\s*([\w_]+(?:\s+\([^)]+\))?)\s*$", text, re.MULTILINE)
    if m:
        candidate = m.group(1).strip()
        if re.match(r"^[A-Z][A-Z0-9_]+", candidate):
            return candidate

    # Last resort: grep for the python script call to infer module name
    script_map = {
        "create_coiledcoils_worker":  "COILEDCOILS_MAP",
        "create_disorder_worker":     "DISORDER_MAP",
        "create_annotation_worker":   "ANNOTATION_MAP",
        "create_transcript_map":      "TRANSCRIPT_MAP",
        "create_genome_map":          "GENOME_MAP",
        "create_mutation_map":        "MUTATION_MAP",
        "create_conservation_worker": "CONSERVATION_MAP",
        "create_pdb_worker":          "PDB_MAP",
        "create_go_worker":           "GO_MAP",
        "create_dbnsfp_map":          "DBNSFP_MAP",
        "create_ppi_worker":          "PPI_MAP",
        "create_blast_table":         "BLASTP",
        "deepcoil":                   "COILEDCOILS_MAP",
    }
    for key, process_name in script_map.items():
        if key in text:
            # Try to get the chunk tag
            m2 = re.search(r"(chunk_\d+)", text)
            tag = f" ({m2.group(1)})" if m2 else ""
            return f"{process_name}{tag}"
    return None

# ── Legacy .nextflow.log parsing (fallback) ───────────────────────────────────

def parse_nextflow_log(nf_log: Path) -> dict:
    """Parse .nextflow.log — used only when no trace.tsv is available."""
    if not nf_log.exists():
        return {}
    text = nf_log.read_text(errors="replace")
    tasks: dict[str, dict] = {}

    task_pattern = re.compile(
        r"TaskHandler\[id: (\d+); name: ([^;]+); status: (\w+); exit: ([^;]+);"
    )
    for m in task_pattern.finditer(text):
        tid, name, status, exit_code = m.groups()
        tasks[tid] = {
            "name": name.strip(), "status": status,
            "exit": exit_code.strip(), "phase": phase_of(name.strip()),
        }

    cached_pattern = re.compile(r"Cached process > ([^\n]+)")
    cached_id = 100_000
    for m in cached_pattern.finditer(text):
        name = m.group(1).strip()
        if not any(t["name"] == name for t in tasks.values()):
            tasks[str(cached_id)] = {
                "name": name, "status": "COMPLETED",
                "exit": "0", "phase": phase_of(name),
            }
            cached_id += 1

    submitted_pattern = re.compile(r"Submitted process > ([^\n]+)")
    submitted_id = 200_000
    for m in submitted_pattern.finditer(text):
        name = m.group(1).strip()
        existing = next((t for t in tasks.values() if t["name"] == name), None)
        if existing is None:
            tasks[str(submitted_id)] = {
                "name": name, "status": "RUNNING",
                "exit": "-", "phase": phase_of(name),
            }
            submitted_id += 1
        elif existing["status"] not in ("COMPLETED", "FAILED"):
            existing["status"] = "RUNNING"
    return tasks

def chunk_progress_from_err(work_dir: Path) -> list[dict]:
    """
    Parse .command.err of running COILEDCOILS_MAP / DISORDER_MAP tasks for
    DEEPCOIL_PROGRESS / AIUPRED_PROGRESS lines written by the workers.
    Returns list of {chunk, done, total, batch, n_batches, elapsed_s}.
    """
    results = []
    if not work_dir.exists():
        return results
    for top in work_dir.iterdir():
        if not top.is_dir() or len(top.name) != 2:
            continue
        for sub in top.iterdir():
            if not sub.is_dir():
                continue
            begin = sub / ".command.begin"
            exitc = sub / ".exitcode"
            cmd_sh = sub / ".command.sh"
            err_f  = sub / ".command.err"
            if not begin.exists() or exitc.exists() or not err_f.exists():
                continue
            # Only COILEDCOILS or DISORDER tasks
            if not cmd_sh.exists():
                continue
            try:
                sh_head = cmd_sh.read_text(errors="replace")[:300]
            except OSError:
                continue
            is_coil   = "deepcoil" in sh_head or "coiledcoils" in sh_head.lower()
            is_dis    = "create_disorder_worker" in sh_head
            if not (is_coil or is_dis):
                continue
            chunk = re.search(r"chunk_(\d+)", sh_head)
            chunk_id = chunk.group(0) if chunk else "?"
            elapsed_s = time.time() - begin.stat().st_mtime
            # Read last PROGRESS line from .command.err
            try:
                # Read last 4KB for speed
                err_text = err_f.read_text(errors="replace")[-4096:]
            except OSError:
                continue
            done = total = batch = n_batches = None
            for line in reversed(err_text.splitlines()):
                m = re.search(r"DEEPCOIL_PROGRESS (\d+)/(\d+) batch (\d+)/(\d+)", line)
                if m:
                    done, total, batch, n_batches = [int(x) for x in m.groups()]
                    break
                m = re.search(r"DEEPCOIL_INIT .* total=(\d+) batches=(\d+)", line)
                if m:
                    total, n_batches = int(m.group(1)), int(m.group(2))
                    done, batch = 0, 0
                    break
            if total:
                results.append({
                    "chunk": chunk_id, "done": done or 0, "total": total,
                    "batch": batch or 0, "n_batches": n_batches or "?",
                    "elapsed_s": elapsed_s, "is_coil": is_coil,
                })
    return sorted(results, key=lambda x: x["chunk"])


def start_time_from_log(run_log: Path):
    if not run_log.exists():
        return None
    text = run_log.read_text(errors="replace")
    matches = re.findall(
        r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*run START", text
    )
    for ts_str in reversed(matches):
        try:
            return datetime.strptime(ts_str, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            continue
    return None

# ── BLAT chunk details ────────────────────────────────────────────────────────

EXPECTED_PSL_PER_CHUNK = 41_000

def find_blat_dirs(work_dir: Path) -> list[dict]:
    by_chunk: dict[str, dict] = {}
    if not work_dir.exists():
        return []
    for top in sorted(work_dir.iterdir()):
        if not top.is_dir() or len(top.name) != 2:
            continue
        for sub in sorted(top.iterdir()):
            if not sub.is_dir():
                continue
            fasta_files = list(sub.glob("chunk_*.fasta"))
            if not fasta_files:
                continue
            fasta = fasta_files[0]
            chunk_name = fasta.stem
            psl_files = list(sub.glob("chunk_*.psl"))
            psl = psl_files[0] if psl_files else None
            exit_file = sub / ".exitcode"
            cmd_out   = sub / ".command.out"
            # Skip ghost dirs from killed/failed runs: old PSL (>2h) with
            # either no exitcode, or a non-zero exitcode (SIGTERM=143, etc.)
            if psl and psl.exists():
                psl_age = time.time() - psl.stat().st_mtime
                if psl_age > 7200:
                    if not exit_file.exists():
                        continue
                    try:
                        if exit_file.read_text().strip() != "0":
                            continue
                    except OSError:
                        continue
            entry = {
                "chunk": chunk_name, "fasta": fasta, "psl": psl,
                "exit_file": exit_file, "cmd_out": cmd_out,
            }
            existing = by_chunk.get(chunk_name)
            if existing is None:
                by_chunk[chunk_name] = entry
            else:
                new_done = exit_file.exists()
                old_done = existing["exit_file"].exists()
                if new_done and not old_done:
                    by_chunk[chunk_name] = entry
                elif not new_done and not old_done:
                    new_sz = psl.stat().st_size if psl and psl.exists() else 0
                    old_psl = existing["psl"]
                    old_sz = old_psl.stat().st_size if old_psl and old_psl.exists() else 0
                    if new_sz > old_sz:
                        by_chunk[chunk_name] = entry
    return sorted(by_chunk.values(), key=lambda x: x["chunk"])

def blat_status(info: dict) -> dict:
    exit_file = info["exit_file"]
    psl       = info["psl"]
    cmd_out   = info["cmd_out"]
    done = exit_file.exists() and exit_file.read_text().strip() == "0"
    final_hits = None
    if cmd_out.exists():
        t = cmd_out.read_text(errors="replace")
        m = re.search(r"Raw BLAT hits: (\d+)", t)
        if m:
            final_hits = int(m.group(1))
    try:
        with open(info["fasta"]) as f:
            n_seqs = sum(1 for line in f if line.startswith(">"))
    except OSError:
        n_seqs = 0
    psl_lines = 0
    psl_mtime = None
    if psl and psl.exists():
        try:
            psl_lines = int(
                subprocess.run(["wc", "-l", str(psl)], capture_output=True,
                               text=True).stdout.split()[0]
            )
            psl_mtime = datetime.fromtimestamp(psl.stat().st_mtime)
        except Exception:
            pass
    if done:
        pct = 100.0
    elif psl_lines > 0:
        pct = min(99.0, psl_lines / EXPECTED_PSL_PER_CHUNK * 100)
    else:
        pct = 0.0
    return {
        "chunk": info["chunk"], "done": done, "n_seqs": n_seqs,
        "psl_lines": psl_lines, "final_hits": final_hits,
        "pct": pct, "psl_mtime": psl_mtime,
    }

# ── Main display ──────────────────────────────────────────────────────────────

def print_status(trace_path: Path | None, nf_log: Path, work_dir: Path,
                 run_log: Path, _start: list, data_source: str):
    now = datetime.now()

    # ── Gather task data ─────────────────────────────────────────────────────
    if trace_path and trace_path.exists():
        # Primary: trace.tsv for completed/failed tasks in THIS run
        completed_tasks = parse_trace(trace_path)
        tasks = dict(completed_tasks)

        # Supplement: .nextflow.log has "Cached process >" lines for tasks
        # served from previous-session cache — these won't appear in trace yet.
        if nf_log.exists():
            cached_pattern = re.compile(r"Cached process > ([^\n]+)")
            cached_id = 100_000
            for m in cached_pattern.finditer(nf_log.read_text(errors="replace")):
                name = m.group(1).strip()
                if not any(t["name"] == name for t in tasks.values()):
                    tasks[str(cached_id)] = {
                        "name": name, "status": "COMPLETED",
                        "exit": "0", "phase": phase_of(name),
                    }
                    cached_id += 1

        # Work-dir scan for currently running tasks
        running_list = find_running_in_workdir(work_dir)
        running_names_in_tasks = {t["name"] for t in tasks.values()
                                   if t["status"] == "RUNNING"}
        rid = 900_000
        for rt in running_list:
            if rt["name"] not in running_names_in_tasks:
                tasks[str(rid)] = rt
                rid += 1
        source_label = f"trace: {trace_path}  +  cached from log"
    else:
        # Fallback: .nextflow.log
        tasks = parse_nextflow_log(nf_log)
        source_label = f"log: {nf_log}  {DIM}(no trace found — stub runs may pollute this){RESET}"

    # ── Start time ───────────────────────────────────────────────────────────
    if not _start:
        # Try trace mtime as approximate start, or run_log
        st = start_time_from_log(run_log)
        if st is None and trace_path and trace_path.exists():
            st = datetime.fromtimestamp(trace_path.stat().st_mtime - 3600)
        if st:
            _start.append(st)
    start_time = _start[0] if _start else None
    elapsed = (now - start_time).total_seconds() if start_time else 0

    done_tasks    = {t["name"]: t for t in tasks.values() if t["status"] == "COMPLETED"}
    running_tasks = {t["name"]: t for t in tasks.values() if t["status"] == "RUNNING"}
    failed_tasks  = {t["name"]: t for t in tasks.values() if t["status"] == "FAILED"}

    print(CLEAR, end="")
    print(f"{BOLD}DisCanVisFlow — Live Progress{RESET}  "
          f"{DIM}{now.strftime('%Y-%m-%d %H:%M:%S')}  elapsed: {human_time(elapsed)}{RESET}")
    print("═" * 72)

    # ── 1. Currently running ─────────────────────────────────────────────────
    print(f"\n{BOLD}▶ RUNNING NOW ({len(running_tasks)} tasks){RESET}")
    if running_tasks:
        for name, t in sorted(running_tasks.items()):
            elapsed_str = ""
            if "elapsed_s" in t:
                elapsed_str = f"  {DIM}({human_time(t['elapsed_s'])}){RESET}"
            print(f"    {CYAN}●{RESET} {name}{elapsed_str}")
    else:
        print(f"    {DIM}(none — pipeline idle, starting, or done){RESET}")

    if failed_tasks:
        print(f"\n{RED}{BOLD}✗ FAILED ({len(failed_tasks)}){RESET}")
        for name in failed_tasks:
            print(f"    {RED}✗{RESET} {name}")

    # ── 2. Phase progress table ───────────────────────────────────────────────
    print(f"\n{BOLD}Pipeline phases{RESET}")
    print(f"  {'Phase':<16} {'Done':>5}  {'Running':>7}  Status bar")
    print(f"  {'─'*60}")

    total_done = total_running = total_tasks = 0
    by_phase: dict[str, dict] = defaultdict(lambda: {"done": 0, "running": 0, "failed": 0, "total": 0})
    for t in tasks.values():
        ph = t["phase"]
        by_phase[ph]["total"] += 1
        if t["status"] == "COMPLETED":
            by_phase[ph]["done"] += 1
        elif t["status"] == "RUNNING":
            by_phase[ph]["running"] += 1
        elif t["status"] == "FAILED":
            by_phase[ph]["failed"] += 1

    for phase_name, _ in PHASES:
        d = by_phase.get(phase_name, {"done": 0, "running": 0, "failed": 0, "total": 0})
        if d["total"] == 0:
            continue
        n_done = d["done"]
        n_run  = d["running"]
        n_tot  = d["total"]
        total_done    += n_done
        total_running += n_run
        total_tasks   += n_tot
        frac = n_done / n_tot if n_tot else 0
        b = bar(frac, width=24)
        if n_run > 0:
            status = f"  {CYAN}{n_run} running{RESET}"
            icon   = f"{CYAN}▶{RESET}"
        elif n_done == n_tot:
            status = f"  {GREEN}complete{RESET}"
            icon   = f"{GREEN}✔{RESET}"
        else:
            status = ""
            icon   = " "
        fail_note = f"  {RED}{d['failed']} FAILED{RESET}" if d["failed"] else ""
        print(f"  {icon} {phase_name:<14}  {n_done:>4}/{n_tot:<4}  {n_run:>4} run  "
              f"[{b}] {frac*100:5.1f}%{status}{fail_note}")

    if by_phase.get("Other", {}).get("total", 0):
        d = by_phase["Other"]
        print(f"    Other: {d['done']}/{d['total']}")

    print(f"  {'─'*60}")
    print(f"  {'Total':<16}  {total_done:>4}/{total_tasks:<4}  {total_running:>4} run")

    # ── 3. Scatter module breakdown ───────────────────────────────────────────
    scatter_modules: dict[str, dict] = defaultdict(lambda: {"done": 0, "running": 0, "total": 0})
    scatter_keywords = [
        "GENOME_MAP", "DISORDER_MAP", "MUTATION_MAP", "TRANSCRIPT_MAP",
        "DBNSFP_MAP", "ALPHAMISSENSE_MAP", "DEPMAP_MAP", "PEM_TRANSFER",
        "MAVEDB_MAP", "PROTEINGYM_MAP", "SCANSITE_MAP", "CONSERVATION_MAP",
        "COILEDCOILS_MAP",
    ]
    for t in tasks.values():
        for kw in scatter_keywords:
            if kw in t["name"].upper():
                scatter_modules[kw]["total"] += 1
                if t["status"] == "COMPLETED":
                    scatter_modules[kw]["done"] += 1
                elif t["status"] == "RUNNING":
                    scatter_modules[kw]["running"] += 1
                break

    active_scatter = {k: v for k, v in scatter_modules.items() if v["total"] > 0}
    if active_scatter:
        print(f"\n{BOLD}Scatter module breakdown{RESET}  "
              f"{DIM}(modules that run N chunks in parallel){RESET}")
        print(f"  {'Module':<24} {'Done':>5}  {'Running':>7}  Progress")
        print(f"  {'─'*60}")
        for kw, d in sorted(active_scatter.items()):
            n, r, tot = d["done"], d["running"], d["total"]
            frac = n / tot if tot else 0
            b = bar(frac, width=20)
            icon = f"{GREEN}✔{RESET}" if n == tot else (f"{CYAN}▶{RESET}" if r else " ")
            print(f"  {icon} {kw:<22}  {n:>4}/{tot:<4}  {r:>4} run  [{b}] {frac*100:.0f}%")

    # ── 3b. Per-chunk DeepCoil / Disorder progress ───────────────────────────
    chunk_progs = chunk_progress_from_err(work_dir)
    if chunk_progs:
        print(f"\n{BOLD}DeepCoil / Disorder chunk detail{RESET}  "
              f"{DIM}(per-protein progress from .command.err){RESET}")
        print(f"  {'Chunk':<12} {'Proteins':>10}  {'Progress':<28}  {'Elapsed':>8}  Batch")
        print(f"  {'─'*72}")
        for cp in chunk_progs:
            done, total = cp["done"], cp["total"]
            frac = done / total if total else 0
            b    = bar(frac, width=24)
            elapsed_str = human_time(cp["elapsed_s"])
            batch_str = f"{cp['batch']}/{cp['n_batches']}" if cp["batch"] else "init"
            label = "coil" if cp["is_coil"] else "dis "
            print(f"  {CYAN}▶{RESET} {cp['chunk']:<10} [{label}]  {done:>5}/{total:<5}  "
                  f"[{b}] {frac*100:5.1f}%  {elapsed_str:>8}  batch {batch_str}")

    # ── 4. BLAT detail ────────────────────────────────────────────────────────
    blat_dirs = find_blat_dirs(work_dir)
    if blat_dirs:
        blat_statuses = [blat_status(d) for d in blat_dirs]
        n_blat_done   = sum(1 for s in blat_statuses if s["done"])
        blat_complete = n_blat_done == len(blat_statuses)

        print(f"\n{BOLD}BLAT_ALIGN detail{RESET}  "
              f"({len(blat_dirs)} chunks, {n_blat_done}/{len(blat_dirs)} done)")

        if not blat_complete:
            print(f"  {'Chunk':<14} {'Seqs':>5}  {'Progress':<30}  Status")
            print(f"  {'─'*58}")
            for st in blat_statuses:
                if st["done"]:
                    continue
                b = bar(st["pct"] / 100, width=26)
                psl_mtime = st["psl_mtime"]
                if psl_mtime:
                    idle = (now - psl_mtime).total_seconds()
                    status_str = f"stalled {idle:.0f}s" if idle > 90 else "writing PSL"
                else:
                    status_str = "starting"
                print(f"  {CYAN}▶{RESET} {st['chunk']:<12}  {st['n_seqs']:>5}  "
                      f"[{b}] {st['pct']:5.1f}%  {status_str}")
            print(f"\n  {GREEN}✔{RESET} {n_blat_done} done", end="")
            running_blat = [s for s in blat_statuses if not s["done"] and s["psl_lines"] > 0]
            if running_blat:
                avg_pct = sum(s["pct"] for s in running_blat) / len(running_blat)
                print(f"  |  {len(running_blat)} running (avg {avg_pct:.0f}%)", end="")
            print()
        else:
            total_hits = sum(s["final_hits"] or 0 for s in blat_statuses)
            print(f"  {GREEN}✔ All {len(blat_dirs)} chunks complete — "
                  f"{total_hits:,} total PSL hits{RESET}")

    # ── 5. ETA estimate ───────────────────────────────────────────────────────
    print(f"\n{BOLD}ETA estimate{RESET}  {DIM}(rough — based on observed timings){RESET}")

    running_names = set(running_tasks.keys())
    any_blast_blat = any(
        kw in n for n in running_names for kw in ("BLASTP", "BLAT_ALIGN", "ID_MAP", "MERGE_BLAT", "MERGE_BLAST")
    )
    any_genome = any("GENOME_MAP" in n or "SEQUENCE_PROCESS" in n for n in running_names)
    any_annotation = any(
        kw in n for n in running_names
        for kw in ("ANNOTATION_MAP", "DISORDER_MAP", "MUTATION_MAP",
                   "DBNSFP_MAP", "PDB_MAP", "GO_MAP", "ALPHAMISSENSE_MAP",
                   "COILEDCOILS_MAP")
    )
    any_final = any("TRANSCRIPT_MAP" in n or "MAPPING_REPORT" in n for n in running_names)

    lo = hi = timedelta(0)

    if any_blast_blat:
        print(f"  Current : BLAST/BLAT phase")
        print(f"  Seq/Genome after : ~5–20 min")
        print(f"  Annotation (20×) : ~60–120 min")
        print(f"  Final            : ~10–20 min")
        lo = timedelta(minutes=5+60+10)
        hi = timedelta(minutes=20+120+20)
    elif any_genome:
        print(f"  Current : Genome/Sequence processing")
        print(f"  Annotation (20×) : ~60–120 min")
        print(f"  Final            : ~10–20 min")
        lo = timedelta(minutes=60+10)
        hi = timedelta(minutes=120+20)
    elif any_annotation:
        coil_done = scatter_modules.get("COILEDCOILS_MAP", {}).get("done", 0)
        coil_tot  = scatter_modules.get("COILEDCOILS_MAP", {}).get("total", 0)
        coil_run  = scatter_modules.get("COILEDCOILS_MAP", {}).get("running", 0)
        dis_done  = scatter_modules.get("DISORDER_MAP", {}).get("done", 0)
        dis_tot   = scatter_modules.get("DISORDER_MAP", {}).get("total", 0)

        if coil_tot > 0 and coil_done < coil_tot:
            # COILEDCOILS is the bottleneck — estimate from deepcoil timing
            # Each chunk ~3h / 20 chunks; but running in parallel
            remaining_chunks = coil_tot - coil_done
            # Use elapsed time of longest-running chunk as proxy
            long_run = max((t.get("elapsed_s", 0) for t in running_tasks.values()), default=0)
            # DeepCoil on full proteome chunk ≈ 3-4h per chunk
            est_chunk_s = max(long_run, 7200)   # at least 2h remaining estimate
            print(f"  Current : COILEDCOILS_MAP scatter ({coil_done}/{coil_tot} done, {coil_run} running)")
            print(f"  Chunks running for: {human_time(long_run)}")
            print(f"  DeepCoil is CPU-only LSTM — ~3-4h per chunk on 64-CPU server")
            print(f"  Remaining chunks: {remaining_chunks}  (running {coil_run} in parallel)")
            # When all coil chunks finish, DISORDER + TRANSCRIPT take ~30 min
            print(f"  After coiledcoils: DISORDER_MAP + TRANSCRIPT ~30 min")
            lo = timedelta(seconds=est_chunk_s * 0.5)
            hi = timedelta(seconds=est_chunk_s * 1.5 + 1800)
        elif dis_tot > 0 and elapsed > 300 and dis_done > 0:
            rate = dis_done / elapsed
            remaining_tasks = dis_tot - dis_done
            est_s = remaining_tasks / rate if rate > 0 else 3600
            print(f"  Current : Annotation scatter ({dis_done}/{dis_tot} DISORDER chunks done)")
            print(f"  Est remaining annotation: ~{human_time(est_s)}")
            print(f"  Final TRANSCRIPT/REPORT: ~10–20 min")
            lo = timedelta(seconds=est_s + 600)
            hi = timedelta(seconds=est_s + 1200)
        else:
            print(f"  Current : Annotation scatter (still warming up)")
            lo = timedelta(minutes=60)
            hi = timedelta(minutes=120)
    elif any_final:
        print(f"  Current : TRANSCRIPT_MAP / MAPPING_REPORT (~10–20 min remaining)")
        lo = timedelta(minutes=5)
        hi = timedelta(minutes=20)
    else:
        print(f"  Pipeline appears idle — check log for errors or all done")

    if lo.total_seconds() > 0:
        eta_lo = now + lo
        eta_hi = now + hi
        print(f"  → ETA: {eta_lo.strftime('%H:%M')} – {eta_hi.strftime('%H:%M')}")

    print()
    print("═" * 72)
    print(f"{DIM}{source_label}  |  Ctrl+C to stop{RESET}")


def main():
    ap = argparse.ArgumentParser(description="Watch DisCanVisFlow pipeline progress")
    ap.add_argument("--project",  default=None,
                    help="Project name (e.g. discanvis). Auto-detected from results/ if omitted.")
    ap.add_argument("--outdir",   default="results",
                    help="Output directory root (default: results)")
    ap.add_argument("--trace",    default=None,
                    help="Explicit path to trace.tsv. Overrides --project/--outdir auto-detect.")
    ap.add_argument("--log",      default=".nextflow.log",
                    help="Nextflow log (fallback when no trace found)")
    ap.add_argument("--run_log",  default="logs/discanvis_full_run.log")
    ap.add_argument("--work",     default="work/local")
    ap.add_argument("--interval", type=int, default=30)
    ap.add_argument("--once",     action="store_true")
    args = ap.parse_args()

    # Resolve trace path
    if args.trace:
        trace_path = Path(args.trace)
    else:
        trace_path = find_trace(args.project, args.outdir)

    if trace_path and trace_path.exists():
        data_source = f"trace: {trace_path}"
    else:
        data_source = f"log (fallback): {args.log}"
        trace_path  = None

    nf_log   = Path(args.log)
    run_log  = Path(args.run_log)
    work_dir = Path(args.work)
    _start: list = []

    if args.once:
        print_status(trace_path, nf_log, work_dir, run_log, _start, data_source)
        return

    print(f"Monitoring {data_source} | work: {work_dir} | "
          f"refresh: {args.interval}s | Ctrl+C to stop")
    try:
        while True:
            print_status(trace_path, nf_log, work_dir, run_log, _start, data_source)
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nMonitoring stopped.")


if __name__ == "__main__":
    main()
