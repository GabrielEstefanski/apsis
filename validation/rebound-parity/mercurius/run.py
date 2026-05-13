"""Orchestrator for REBOUND parity — Mercurius.

Runs the apsis side (Cargo example), then the REBOUND side (Python),
then the comparator. Exits with the comparator's exit code.

Usage from the scenario directory:

    python run.py

With explicit paths:

    python run.py --workspace-root /path/to/apsis --output-dir ./out

Exit codes:

- 0 — all metrics within tolerance.
- 1 — input file error in the comparator.
- 2 — at least one metric exceeded tolerance.
- non-zero (other) — apsis or REBOUND side failed.

Protocol notebook:
    docs/experiments/2026-05-13-rebound-parity-mercurius.md
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCENARIO_DIR: Path = Path(__file__).resolve().parent
# scenario dir is `<root>/validation/rebound-parity/mercurius`, so the
# workspace root is three levels up.
DEFAULT_WORKSPACE_ROOT: Path = SCENARIO_DIR.parent.parent.parent


def main() -> int:
    args = parse_args()
    workspace_root = Path(args.workspace_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    apsis_csv = output_dir / "apsis.csv"
    rebound_csv = output_dir / "rebound.csv"

    # ── 1. Apsis side via cargo ─────────────────────────────────────────── #
    print("=" * 72)
    print("[1/3] apsis Mercurius (Cargo example)")
    print("=" * 72)
    rc = subprocess.call(
        [
            "cargo", "run", "--release",
            "--example", "rebound_parity_mercurius",
            "-p", "apsis",
            "--",
            "--output", str(apsis_csv),
        ],
        cwd=workspace_root,
    )
    if rc != 0:
        print(f"ERROR: apsis side failed with exit code {rc}", file=sys.stderr)
        return rc

    # ── 2. REBOUND side ─────────────────────────────────────────────────── #
    print()
    print("=" * 72)
    print("[2/3] REBOUND MERCURIUS (Python)")
    print("=" * 72)
    rc = subprocess.call(
        [
            sys.executable, str(SCENARIO_DIR / "rebound_side.py"),
            "--apsis-csv", str(apsis_csv),
            "--output", str(rebound_csv),
        ],
    )
    if rc != 0:
        print(f"ERROR: REBOUND side failed with exit code {rc}", file=sys.stderr)
        return rc

    # ── 3. Comparator ───────────────────────────────────────────────────── #
    print()
    print("=" * 72)
    print("[3/3] Comparator")
    print("=" * 72)
    rc = subprocess.call(
        [
            sys.executable, str(SCENARIO_DIR / "compare.py"),
            "--apsis-csv", str(apsis_csv),
            "--rebound-csv", str(rebound_csv),
            "--output-dir", str(output_dir),
        ],
    )
    return rc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the full apsis/REBOUND parity comparison for Mercurius.",
    )
    parser.add_argument(
        "--workspace-root",
        default=str(DEFAULT_WORKSPACE_ROOT),
        help="Path to the apsis Cargo workspace root "
        "(default: inferred from script location).",
    )
    parser.add_argument(
        "--output-dir",
        default=str(SCENARIO_DIR / "out"),
        help="Directory for generated CSVs and reports (default: ./out).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
