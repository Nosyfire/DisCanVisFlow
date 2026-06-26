#!/usr/bin/env python3
"""
bin/watch_run.py — Live progress monitor for a running DisCanVisFlow pipeline.

Shows all pipeline phases with task-by-task status, scatter-chunk progress,
BLAT chunk detail, and per-phase ETA estimates.

Usage:
  python bin/watch_run.py                   # auto-detect, refresh every 30s
  python bin/watch_run.py --interval 10     # refresh every 10s
  python bin/watch_run.py --once            # print once and exit
  python bin/watch_run.py --work work/local
"""

import argparse
import re
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

# ── Terminal helpers ──────────────────────────────────────────────────────────

CLEAR = "\033[2J\033[H"
BOLD  = "\033[1m"
DIM   = "\033[2m"
GREEN = "\033[32m"
CYAN  = "\033[36m"
YELLOW = "\033[33m"
RED   = "\033[31m"
RESET = "\033[0m"

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
    """Return human-readable time since file was last modified."""
    try:
        s = time.time() - path.stat().st_mtime
        return human_time(s)
    except OSError:
        return "?"

# ── Phase definitions ─────────────────────────────────────────────────────────
# Each phase = list of substrings that match task names belonging to that phase.

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
                        "FINCHES_MAP", "POSITION_BASED_MAP"]),
    ("Finalize",       ["TRANSCRIPT_MAP", "MAPPING_REPORT"]),
]

def phase_of(name: str) -> str:
    for phase_name, keywords in PHASES:
        if any(kw in name.upper() for kw in keywords):
            return phase_name
    return "Other"

# ── Nextflow log parsing ──────────────────────────────────────────────────────

def parse_nextflow_log(nf_log: Path):
    if not nf_log.exists():
        return {}
    text = nf_log.read_text(errors="replace")
    tasks: dict[str, dict] = {}

    # Format 1 (active tasks — all runs):
    # TaskHandler[id: N; name: X (tag); status: RUNNING/COMPLETED/FAILED; exit: 0; ...]
    task_pattern = re.compile(
        r"TaskHandler\[id: (\d+); name: ([^;]+); status: (\w+); exit: ([^;]+);"
    )
    for m in task_pattern.finditer(text):
        tid, name, status, exit_code = m.groups()
        name = name.strip()
        tasks[tid] = {
            "name":   name,
            "status": status,
            "exit":   exit_code.strip(),
            "phase":  phase_of(name),
        }

    # Format 2 (cached tasks on -resume):
    # INFO n.p.TaskProcessor - [ab/cdef12] Cached process > NAME (tag)
    cached_pattern = re.compile(r"Cached process > ([^\n]+)")
    cached_id = 100_000  # synthetic IDs starting high to avoid collision
    for m in cached_pattern.finditer(text):
        name = m.group(1).strip()
        already = any(t["name"] == name for t in tasks.values())
        if not already:
            tasks[str(cached_id)] = {
                "name": name, "status": "COMPLETED",
                "exit": "0", "phase": phase_of(name),
            }
            cached_id += 1

    # Format 3 (submitted but not yet done — shows as RUNNING):
    # INFO nextflow.Session - [8e/713ac6] Submitted process > NAME (tag)
    submitted_pattern = re.compile(r"Submitted process > ([^\n]+)")
    submitted_id = 200_000
    for m in submitted_pattern.finditer(text):
        name = m.group(1).strip()
        # Only add if not already seen as COMPLETED or in TaskHandler
        existing = next((t for t in tasks.values() if t["name"] == name), None)
        if existing is None:
            tasks[str(submitted_id)] = {
                "name": name, "status": "RUNNING",
                "exit": "-", "phase": phase_of(name),
            }
            submitted_id += 1
        elif existing["status"] not in ("COMPLETED", "FAILED"):
            # upgrade status to RUNNING if it's only been seen as submitted
            existing["status"] = "RUNNING"

    return tasks

def start_time_from_log(run_log: Path):
    """Return the MOST RECENT 'run START' timestamp — handles multiple restart entries."""
    if not run_log.exists():
        return None
    text = run_log.read_text(errors="replace")
    matches = re.findall(
        r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*run START", text
    )
    for ts_str in reversed(matches):  # latest first
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

def print_status(nf_log: Path, work_dir: Path, run_log: Path, _start: list):
    now = datetime.now()
    tasks = parse_nextflow_log(nf_log)

    if not _start:
        st = start_time_from_log(run_log)
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
            print(f"    {CYAN}●{RESET} {name}")
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

    # group tasks by phase
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
    # Group running/done tasks by base module name (strip parenthetical suffix)
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

    # ── 4. BLAT detail ────────────────────────────────────────────────────────
    blat_dirs = find_blat_dirs(work_dir)
    if blat_dirs:
        blat_statuses = [blat_status(d) for d in blat_dirs]
        n_blat_done   = sum(1 for s in blat_statuses if s["done"])
        blat_complete = n_blat_done == len(blat_statuses)

        print(f"\n{BOLD}BLAT_ALIGN detail{RESET}  "
              f"({len(blat_dirs)} chunks, {n_blat_done}/{len(blat_dirs)} done)")

        if not blat_complete:
            # show only running/incomplete chunks to keep display short
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
                   "DBNSFP_MAP", "PDB_MAP", "GO_MAP", "ALPHAMISSENSE_MAP")
    )
    any_final = any("TRANSCRIPT_MAP" in n or "MAPPING_REPORT" in n for n in running_names)

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
        # estimate from scatter progress
        ann_done = scatter_modules.get("DISORDER_MAP", {}).get("done", 0)
        ann_tot  = scatter_modules.get("DISORDER_MAP", {}).get("total", 20)
        if ann_tot > 0 and elapsed > 300 and ann_done > 0:
            rate = ann_done / elapsed  # tasks per second
            remaining_tasks = ann_tot - ann_done
            est_s = remaining_tasks / rate if rate > 0 else 3600
            print(f"  Current : Annotation scatter ({ann_done}/{ann_tot} DISORDER chunks done)")
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
        lo = hi = timedelta(0)
        print(f"  Pipeline appears idle — check log for errors")

    if lo.total_seconds() > 0:
        eta_lo = now + lo
        eta_hi = now + hi
        print(f"  → ETA: {eta_lo.strftime('%H:%M')} – {eta_hi.strftime('%H:%M')}")

    print()
    print("═" * 72)
    print(f"{DIM}nf log: {nf_log}  |  Ctrl+C to stop{RESET}")


def main():
    ap = argparse.ArgumentParser(description="Watch DisCanVisFlow pipeline progress")
    ap.add_argument("--log",      default=".nextflow.log")
    ap.add_argument("--run_log",  default="logs/discanvis_full_run.log")
    ap.add_argument("--work",     default="work/local")
    ap.add_argument("--interval", type=int, default=30)
    ap.add_argument("--once",     action="store_true")
    args = ap.parse_args()

    nf_log   = Path(args.log)
    run_log  = Path(args.run_log)
    work_dir = Path(args.work)
    _start: list = []

    if args.once:
        print_status(nf_log, work_dir, run_log, _start)
        return

    print(f"Monitoring {nf_log} | work: {work_dir} | "
          f"refresh: {args.interval}s | Ctrl+C to stop")
    try:
        while True:
            print_status(nf_log, work_dir, run_log, _start)
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nMonitoring stopped.")


if __name__ == "__main__":
    main()
