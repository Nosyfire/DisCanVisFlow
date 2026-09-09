#!/usr/bin/env python3
"""Locate a usable DisCanVisFlow checkout and report what still needs setting up.

Every discanvis-tools skill starts by answering the same three questions:
where is the pipeline, is its conda environment built, and is there already a
completed run whose results can just be sliced instead of recomputed. Answering
those by hand costs several tool calls and is easy to get subtly wrong (the
plugin cache looks like a checkout but is a bad place to write 27 GB of
reference data), so it lives here once.

Usage:
    python find_discanvis.py            # human-readable report
    python find_discanvis.py --json     # machine-readable, for scripting

Exit codes:
    0  a usable checkout was found
    1  no checkout found — the caller should offer to clone one
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# A directory is a DisCanVisFlow checkout only if all of these exist. Checking
# several markers rather than just main.nf avoids matching an unrelated
# Nextflow project that happens to sit in the search path.
MARKERS = ("main.nf", "nextflow.config", "bin", "modules")

CONDA_ENV_NAME = "discanvis"


def is_checkout(path: Path) -> bool:
    return all((path / m).exists() for m in MARKERS)


def walk_up(start: Path) -> Path | None:
    """Return the nearest ancestor of `start` that is a checkout."""
    for candidate in (start, *start.parents):
        if is_checkout(candidate):
            return candidate
    return None


def candidates() -> list[tuple[Path, str]]:
    """Places to look, best first, each with why it was considered.

    Order matters. An explicit $DISCANVIS_HOME is an instruction, so it wins.
    The working directory comes next because someone sitting inside a checkout
    means that one. The plugin's own copy is last among real candidates: it is
    a valid checkout, but writing results and a 27 GB reference cache into a
    versioned plugin directory means losing all of it on the next upgrade.
    """
    out: list[tuple[Path, str]] = []

    env_home = os.environ.get("DISCANVIS_HOME")
    if env_home:
        out.append((Path(env_home).expanduser(), "$DISCANVIS_HOME"))

    cwd_match = walk_up(Path.cwd())
    if cwd_match:
        out.append((cwd_match, "current directory"))

    for var in ("CLAUDE_PLUGIN_ROOT", "CURSOR_PLUGIN_ROOT"):
        plugin_root = os.environ.get(var)
        if plugin_root:
            out.append((Path(plugin_root), f"${var} (plugin cache)"))

    home = Path.home()
    for rel in (
        "DisCanVisFlow",
        "discanvisflow",
        "nextflow_discanvis",
        "PycharmProjects/nextflow_discanvis",
        "projects/DisCanVisFlow",
        "src/DisCanVisFlow",
    ):
        out.append((home / rel, f"~/{rel}"))

    return out


def conda_env_exists(name: str = CONDA_ENV_NAME) -> bool:
    conda = shutil.which("conda")
    if not conda:
        return False
    try:
        res = subprocess.run(
            [conda, "env", "list", "--json"],
            capture_output=True, text=True, timeout=60,
        )
        if res.returncode != 0:
            return False
        envs = json.loads(res.stdout).get("envs", [])
        return any(Path(e).name == name for e in envs)
    except (subprocess.SubprocessError, json.JSONDecodeError, OSError):
        return False


def completed_runs(root: Path) -> list[str]:
    """Project names under results/ that have a populated final/ directory.

    These are what makes extraction possible: slicing a finished run takes
    seconds, where recomputing the same genes takes minutes to hours.
    """
    results = root / "results"
    if not results.is_dir():
        return []
    found = []
    for project in sorted(results.iterdir()):
        final = project / "final"
        if final.is_dir() and any(final.iterdir()):
            found.append(project.name)
    return found


def local_config_usable(root: Path) -> bool:
    """Does config/data/local.config point at files that exist here?

    Existence alone is not enough. The file hard-codes absolute paths for one
    machine, so a copy that was moved between machines, or written against a
    layout that has since changed, exists but is useless. Reporting "present"
    for one of those sends the caller into --data local and a wall of
    missing-file errors, when --data discanvis_data would have just worked.
    Sample the paths it declares and require that most of them really resolve.
    """
    cfg = root / "config/data/local.config"
    if not cfg.is_file():
        return False
    try:
        text = cfg.read_text(errors="replace")
    except OSError:
        return False

    paths = re.findall(r"""['"](/[^'"]+)['"]""", text)
    if not paths:
        return False
    hits = sum(1 for p in paths if Path(p).exists())
    return hits >= max(1, len(paths) // 2)


def inspect(root: Path, why: str) -> dict:
    return {
        "root": str(root),
        "found_via": why,
        "is_plugin_cache": "plugin cache" in why,
        "conda_env": conda_env_exists(),
        "nextflow_on_path": shutil.which("nextflow") is not None,
        "completed_runs": completed_runs(root),
        "has_local_config": local_config_usable(root),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", action="store_true", help="emit JSON instead of prose")
    args = ap.parse_args()

    info = None
    for path, why in candidates():
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            continue
        if is_checkout(resolved):
            info = inspect(resolved, why)
            break

    if info is None:
        if args.json:
            print(json.dumps({"found": False}, indent=2))
        else:
            print("No DisCanVisFlow checkout found.")
            print()
            print("Looked in: $DISCANVIS_HOME, the current directory and its parents,")
            print("the plugin's own directory, and the usual spots under ~/.")
            print()
            print("To set one up:")
            print("  git clone https://github.com/Nosyfire/DisCanVisFlow ~/DisCanVisFlow")
            print("  cd ~/DisCanVisFlow && conda env create -f environment.yml")
            print("  export DISCANVIS_HOME=~/DisCanVisFlow")
        return 1

    if args.json:
        print(json.dumps({"found": True, **info}, indent=2))
        return 0

    print(f"Pipeline:      {info['root']}")
    print(f"Found via:     {info['found_via']}")
    print(f"conda env:     {'present' if info['conda_env'] else 'MISSING (conda env create -f environment.yml)'}")
    print(f"nextflow:      {'on PATH' if info['nextflow_on_path'] else 'NOT on PATH'}")
    print(f"local.config:  {'present' if info['has_local_config'] else 'absent — use --data discanvis_data'}")

    runs = info["completed_runs"]
    print(f"finished runs: {', '.join(runs) if runs else 'none — a request means a real pipeline run'}")

    if info["is_plugin_cache"]:
        print()
        print("NOTE: this is the plugin's own copy. It works, but a run writes tens of")
        print("gigabytes of reference data and results into a versioned plugin directory,")
        print("and a plugin upgrade discards all of it. Prefer a normal clone and point")
        print("DISCANVIS_HOME at it before running anything long.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
