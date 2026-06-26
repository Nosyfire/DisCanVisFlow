#!/usr/bin/env python3
"""
bin/watch_run.py — Live progress monitor for a running DisCanVisFlow pipeline.

Reads .nextflow.log for task lifecycle events and tracks per-BLAT-chunk
progress from PSL file sizes so you can watch the run without staring at
the scrolling Nextflow output.

Usage:
  python bin/watch_run.py                         # auto-detect, refresh every 30s
  python bin/watch_run.py --interval 10           # refresh every 10s
  python bin/watch_run.py --once                  # print once and exit (good for cron)
  python bin/watch_run.py --log logs/discanvis_full_run.log --work work/local
"""

import argparse
import os
import re
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
import subprocess

# ── Helpers ──────────────────────────────────────────────────────────────────

def bar(frac: float, width: int = 35) -> str:
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

def human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"

# ── Nextflow log parsing ──────────────────────────────────────────────────────

def parse_nextflow_log(nf_log: Path) -> dict:
    """
    Parse .nextflow.log to get:
      - per-task status (RUNNING / COMPLETED / FAILED)
      - pipeline start time
    """
    if not nf_log.exists():
        return {"tasks": {}, "start_time": None, "raw_line_count": 0}

    text = nf_log.read_text(errors="replace")

    # Task lifecycle entries  e.g.
    # TaskHandler[id: 22; name: BLAT_ALIGN (blat_chunk_0001); status: COMPLETED; exit: 0; ...]
    task_pattern = re.compile(
        r"TaskHandler\[id: (\d+); name: ([^;]+); status: (\w+); exit: ([^;]+);"
    )
    tasks = {}
    for m in task_pattern.finditer(text):
        tid, name, status, exit_code = m.groups()
        tasks[tid] = {"name": name.strip(), "status": status, "exit": exit_code.strip()}

    # Pipeline start time from our wrapper log header
    start_match = re.search(
        r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*run START", text
    )
    start_time = None
    if start_match:
        try:
            start_time = datetime.strptime(start_match.group(1), "%Y-%m-%d %H:%M:%S")
        except ValueError:
            pass

    return {"tasks": tasks, "start_time": start_time, "raw_line_count": text.count("\n")}


# ── BLAT progress ─────────────────────────────────────────────────────────────

def find_blat_dirs(work_dir: Path) -> list[dict]:
    """
    Find BLAT_ALIGN task work directories under work_dir.
    A BLAT dir is any 2-level subdir containing chunk_NNNN.fasta.
    Deduplicates by chunk name: prefers completed > has-PSL > most-recent.
    """
    by_chunk: dict[str, dict] = {}
    for top in sorted(work_dir.iterdir()):
        if not top.is_dir() or len(top.name) != 2:
            continue
        for sub in sorted(top.iterdir()):
            if not sub.is_dir():
                continue
            fasta_files = list(sub.glob("chunk_*.fasta"))
            if fasta_files:
                fasta = fasta_files[0]
                chunk_name = fasta.stem          # e.g. chunk_0006
                psl_files = list(sub.glob("chunk_*.psl"))
                psl = psl_files[0] if psl_files else None
                cmd_out = sub / ".command.out"
                exit_file = sub / ".exitcode"
                entry = {
                    "chunk": chunk_name,
                    "fasta": fasta,
                    "psl": psl,
                    "cmd_out": cmd_out,
                    "exit_file": exit_file,
                    "work_dir": sub,
                }
                # Prefer: completed > has real PSL content > first-seen
                existing = by_chunk.get(chunk_name)
                if existing is None:
                    by_chunk[chunk_name] = entry
                else:
                    # Pick completed one; if both same, pick the one with more PSL content
                    new_done = exit_file.exists()
                    old_done = existing["exit_file"].exists()
                    if new_done and not old_done:
                        by_chunk[chunk_name] = entry
                    elif not new_done and not old_done:
                        # pick the one whose PSL file is bigger
                        new_psl_sz = psl.stat().st_size if psl and psl.exists() else 0
                        old_psl = existing["psl"]
                        old_psl_sz = old_psl.stat().st_size if old_psl and old_psl.exists() else 0
                        if new_psl_sz > old_psl_sz:
                            by_chunk[chunk_name] = entry
    return sorted(by_chunk.values(), key=lambda x: x["chunk"])


def blat_chunk_status(info: dict) -> dict:
    """
    Return status dict for one BLAT chunk:
      done, n_seqs, psl_lines, final_hits, pct_done, rate_seqs_per_min
    """
    chunk = info["chunk"]
    fasta = info["fasta"]
    psl = info["psl"]
    cmd_out = info["cmd_out"]
    exit_file = info["exit_file"]

    # Is it done?
    done = exit_file.exists() and exit_file.read_text().strip() == "0"

    # Final hit count from .command.out (written when BLAT finishes)
    final_hits = None
    if cmd_out.exists():
        out_text = cmd_out.read_text(errors="replace")
        m = re.search(r"Raw BLAT hits: (\d+)", out_text)
        if m:
            final_hits = int(m.group(1))
        seqs_m = re.search(r"Searched \d+ bases in (\d+) sequences", out_text)
        n_seqs_from_log = int(seqs_m.group(1)) if seqs_m else None
    else:
        n_seqs_from_log = None

    # Count sequences in the chunk FASTA (cached via fasta mtime)
    try:
        with open(fasta) as f:
            n_seqs = sum(1 for line in f if line.startswith(">"))
    except OSError:
        n_seqs = n_seqs_from_log or 0

    # PSL file stats (growing file = still running)
    psl_lines = 0
    psl_mtime = None
    if psl and psl.exists():
        try:
            psl_lines = int(
                subprocess.run(["wc", "-l", str(psl)],
                               capture_output=True, text=True).stdout.split()[0]
            )
            psl_mtime = datetime.fromtimestamp(psl.stat().st_mtime)
        except Exception:
            pass

    # Estimate progress: use calibration from completed chunks
    # Average raw PSL hits / chunk from empirical data ≈ 41 000 for 6 936 seqs
    # (varies by chunk content; used only when chunk is still running)
    EXPECTED_PSL_PER_CHUNK = 41_000
    if done and final_hits:
        pct = 100.0
    elif psl_lines > 0:
        pct = min(99.0, psl_lines / EXPECTED_PSL_PER_CHUNK * 100)
    else:
        pct = 0.0

    return {
        "chunk": chunk,
        "done": done,
        "n_seqs": n_seqs,
        "psl_lines": psl_lines,
        "final_hits": final_hits,
        "pct": pct,
        "psl_mtime": psl_mtime,
    }


# ── Display ───────────────────────────────────────────────────────────────────

CLEAR = "\033[2J\033[H"
BOLD  = "\033[1m"
DIM   = "\033[2m"
GREEN = "\033[32m"
CYAN  = "\033[36m"
RESET = "\033[0m"

def print_status(nf_log: Path, work_dir: Path, run_log: Path, start_time_cache: list):
    now = datetime.now()
    print(CLEAR, end="")
    print(f"{BOLD}DisCanVisFlow — Live Progress Monitor{RESET}  "
          f"{DIM}{now.strftime('%Y-%m-%d %H:%M:%S')}{RESET}")
    print("=" * 72)

    # ── Pipeline log parsing ──────────────────────────────────────────────────
    nf_data = parse_nextflow_log(nf_log)
    tasks = nf_data["tasks"]

    # Use cached start time if log didn't provide one this call
    if nf_data["start_time"]:
        start_time_cache[:] = [nf_data["start_time"]]
    start_time = start_time_cache[0] if start_time_cache else None

    elapsed = (now - start_time).total_seconds() if start_time else 0

    done_tasks    = [t for t in tasks.values() if t["status"] == "COMPLETED"]
    running_tasks = [t for t in tasks.values() if t["status"] == "RUNNING"]
    failed_tasks  = [t for t in tasks.values() if t["status"] == "FAILED"]

    print(f"\n{BOLD}Overall pipeline{RESET}   (log: {nf_log.name})")
    print(f"  Elapsed :  {human_time(elapsed)}")
    print(f"  Done    :  {GREEN}{len(done_tasks):4d}{RESET} tasks completed")
    print(f"  Running :  {CYAN}{len(running_tasks):4d}{RESET} tasks active   "
          f"{DIM}({', '.join(t['name'] for t in running_tasks[:4])}){RESET}")
    if failed_tasks:
        print(f"  FAILED  :  {len(failed_tasks)} tasks — check {nf_log}")

    # ── BLAT progress ─────────────────────────────────────────────────────────
    blat_dirs = find_blat_dirs(work_dir)
    if blat_dirs:
        print(f"\n{BOLD}BLAT_ALIGN chunks{RESET}   "
              f"({len(blat_dirs)} chunks × ~6 936 seqs against hg38.2bit)")
        print(f"  {'Chunk':<14} {'Seqs':>5}  {'PSL':>7}  {'Progress':<37}  Status")
        print(f"  {'-'*65}")

        all_done = True
        slow_pct = 100.0
        slow_chunk = None
        running_etas = []

        for info in blat_dirs:
            st = blat_chunk_status(info)
            icon = f"{GREEN}✔{RESET}" if st["done"] else f"{CYAN}▶{RESET}"
            pct_str = f"{st['pct']:5.1f}%"
            b = bar(st["pct"] / 100, width=30)
            psl_str = f"{st['psl_lines']:7,}"

            eta_str = ""
            if st["done"]:
                status_str = f"done ({st['final_hits']:,} hits)" if st["final_hits"] else "done"
            elif st["psl_mtime"]:
                idle_s = (now - st["psl_mtime"]).total_seconds()
                if idle_s > 60:
                    status_str = f"stalled? ({idle_s:.0f}s idle)"
                else:
                    status_str = "writing PSL"
                    all_done = False
                    if elapsed > 0 and st["pct"] > 5:
                        total_est = elapsed / (st["pct"] / 100)
                        eta_s = total_est - elapsed
                        running_etas.append(eta_s)
                        eta_str = f" ETA {human_time(eta_s)}"
                    if st["pct"] < slow_pct:
                        slow_pct = st["pct"]
                        slow_chunk = st["chunk"]
            else:
                status_str = "waiting"
                all_done = False

            print(f"  {icon} {st['chunk']:<12}  {st['n_seqs']:>5}  {psl_str}  "
                  f"[{b}] {pct_str}  {status_str}{eta_str}")

        n_done = sum(1 for info in blat_dirs if blat_chunk_status(info)["done"])
        print(f"\n  {n_done}/{len(blat_dirs)} BLAT chunks complete", end="")
        if running_etas:
            max_eta = max(running_etas)
            eta_time = now + timedelta(seconds=max_eta)
            print(f"  |  bottleneck ({slow_chunk}) ETA ≈ {human_time(max_eta)} "
                  f"(~{eta_time.strftime('%H:%M')})", end="")
        elif all_done:
            print(f"  {GREEN}— BLAT phase complete!{RESET}", end="")
        print()

    # ── What's next ───────────────────────────────────────────────────────────
    running_names = [t["name"] for t in running_tasks]
    blat_running = any("BLAT" in n for n in running_names)
    blast_running = any("BLASTP" in n for n in running_names)

    print(f"\n{BOLD}Pipeline phase{RESET}")
    if blat_running or blast_running:
        phase = "Phase 1: BLAT + BLASTP  (running in parallel)"
        next_ph = "→ MERGE_BLAT_PSL → GENOME_MAP → scatter annotation (20 chunks)"
    elif any("GENOME_MAP" in n or "MERGE_BLAT" in n or "SEQUENCE" in n for n in running_names):
        phase = "Phase 2: MERGE/GENOME_MAP/SEQUENCE_PROCESS"
        next_ph = "→ 20 scatter annotation chunks (ANNOTATION, DISORDER, MUTATION, ...)"
    elif running_names:
        phase = f"Phase 3: Annotation scatter  ({running_names[0]} ...)"
        next_ph = "→ TRANSCRIPT_MAP → MAPPING_REPORT → derivation scripts"
    else:
        phase = "Idle / checking..."
        next_ph = ""

    print(f"  {CYAN}{phase}{RESET}")
    if next_ph:
        print(f"  {DIM}{next_ph}{RESET}")

    # ── Rough ETA ─────────────────────────────────────────────────────────────
    print(f"\n{BOLD}Rough remaining estimate{RESET}  "
          f"{DIM}(heuristic — varies by disk I/O and annotation density){RESET}")
    if blat_running or blast_running:
        blat_rem = max(running_etas) if running_etas else 0
        print(f"  BLAT completion   : ~{human_time(blat_rem)}")
        print(f"  MERGE + GENOME_MAP: ~15–25 min after BLAT")
        print(f"  Annotation (20×)  : ~30–90 min  "
              f"{DIM}(DBNSFP/MUTATION/DISORDER are bottlenecks){RESET}")
        print(f"  TRANSCRIPT_MAP 20×: ~10–20 min")
        total_rem_lo = blat_rem + 15*60 + 40*60 + 10*60
        total_rem_hi = blat_rem + 25*60 + 90*60 + 20*60
        eta_lo = now + timedelta(seconds=total_rem_lo)
        eta_hi = now + timedelta(seconds=total_rem_hi)
        print(f"  → Total remaining : {human_time(total_rem_lo)} – {human_time(total_rem_hi)}")
        print(f"  → ETA             : {eta_lo.strftime('%H:%M')} – {eta_hi.strftime('%H:%M')}")
    else:
        print(f"  (update ETA once BLAT is done and annotation phase starts)")

    print()
    print("=" * 72)
    print(f"{DIM}tail -f {nf_log}  |  Ctrl+C to stop monitor{RESET}")


def main():
    ap = argparse.ArgumentParser(description="Watch DisCanVisFlow pipeline progress")
    ap.add_argument("--log",      default=".nextflow.log",
                    help="Nextflow internal log (default: .nextflow.log)")
    ap.add_argument("--run_log",  default="logs/discanvis_full_run.log",
                    help="run wrapper log file")
    ap.add_argument("--work",     default="work/local",
                    help="Nextflow work dir (default: work/local)")
    ap.add_argument("--interval", type=int, default=30,
                    help="refresh interval in seconds (default: 30)")
    ap.add_argument("--once",     action="store_true",
                    help="print once and exit")
    args = ap.parse_args()

    nf_log   = Path(args.log)
    run_log  = Path(args.run_log)
    work_dir = Path(args.work)
    start_time_cache: list = []

    # Try to get start time from run log first
    if run_log.exists():
        text = run_log.read_text(errors="replace")
        m = re.search(r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*run START", text)
        if m:
            try:
                start_time_cache.append(
                    datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
                )
            except ValueError:
                pass

    if args.once:
        print_status(nf_log, work_dir, run_log, start_time_cache)
        return

    print(f"Monitoring {nf_log} | work: {work_dir} | refresh: {args.interval}s | Ctrl+C to stop")
    try:
        while True:
            print_status(nf_log, work_dir, run_log, start_time_cache)
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nMonitoring stopped.")


if __name__ == "__main__":
    main()
