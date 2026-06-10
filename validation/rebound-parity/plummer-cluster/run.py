"""Orchestrator for the Plummer cluster parity run.

Steps: softening smoke test -> 1PN registration gate (cargo test) -> apsis
side -> REBOUND side -> comparator. Per-step wall-time goes to stderr: the
phase-0 deliverable includes the cost measurement that sizes phase 1.

Run:
    python run.py --n 256            # phase 0 (informational)
    python run.py --n 1000           # phase 1 (blocked until gates freeze)

Protocol notebook:
    paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE_ROOT = HERE.parents[2]


def step(name: str, cmd: list[str], cwd: Path) -> None:
    print(f"-- {name} --", flush=True)
    t0 = time.perf_counter()
    result = subprocess.run(cmd, cwd=cwd)
    elapsed = time.perf_counter() - t0
    print(f"[time] {name}: {elapsed:.1f}s", file=sys.stderr)
    if result.returncode != 0:
        print(f"FAILED at step '{name}' (exit {result.returncode})", file=sys.stderr)
        sys.exit(result.returncode)


def main() -> int:
    args = parse_args()
    ics = f"ics_n{args.n}.csv"
    out = Path(args.output_dir)
    informational = args.n == 256

    py = sys.executable
    step("smoke: softening convention", [py, "smoke_pair.py", "--ics", ics], HERE)
    step(
        "gate: 1PN registration warning",
        ["cargo", "test", "-p", "apsis-1pn", "--test", "plummer_cluster_registration_gate"],
        WORKSPACE_ROOT,
    )
    step(
        "apsis side",
        [
            "cargo", "run", "--release", "--example", "rebound_parity_plummer_cluster",
            "-p", "apsis", "--",
            "--ics", str(HERE / ics),
            "--output", str(HERE / out / "apsis.csv"),
            "--stats-output", str(HERE / out / "apsis_stats.json"),
        ],
        WORKSPACE_ROOT,
    )
    step(
        "rebound side",
        [
            py, "rebound_side.py", "--ics", ics,
            "--apsis-csv", str(out / "apsis.csv"),
            "--output", str(out / "rebound.csv"),
            "--stats-output", str(out / "rebound_stats.json"),
        ],
        HERE,
    )
    compare_cmd = [
        py, "compare.py", "--ics", ics,
        "--apsis-csv", str(out / "apsis.csv"),
        "--rebound-csv", str(out / "rebound.csv"),
        "--apsis-stats", str(out / "apsis_stats.json"),
        "--rebound-stats", str(out / "rebound_stats.json"),
        "--output-dir", str(out),
    ]
    if informational:
        compare_cmd.append("--informational")
    step("comparator", compare_cmd, HERE)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plummer cluster parity orchestrator.")
    parser.add_argument("--n", type=int, choices=(256, 1000), default=256, help="body count")
    parser.add_argument("--output-dir", default="out", help="artefact directory")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
