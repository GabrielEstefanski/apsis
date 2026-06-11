"""Orchestrator for REBOUND parity — Pythagorean three-body (Burrau 1913).

Runs the apsis side (Cargo example), then the REBOUND side (Python),
then the comparator. Exits with the comparator's exit code.

Usage from the scenario directory:

    python run.py

With explicit paths or a longer-horizon stress run (NOT part of the
gated parity claim, see protocol notebook §Out of scope):

    python run.py --workspace-root /path/to/apsis --output-dir ./out
    python run.py --horizon 200 --output-dir ./out-200tu  # informational stress

Exit codes:

- 0 — all gated (Tier-1 + Tier-2) metrics within tolerance.
- 1 — input file error in the comparator.
- 2 — at least one gated metric exceeded tolerance.
- non-zero (other) — apsis or REBOUND side failed.

Protocol notebook:
    paper/notebooks/2026-04-30-rebound-parity-pythagorean.md
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCENARIO_DIR: Path = Path(__file__).resolve().parent
# scenario dir is `<root>/validation/rebound-parity/pythagorean`, so the
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
    print(f"[1/3] apsis IAS15 (Cargo example) — horizon {args.horizon}")
    print("=" * 72)
    apsis_cmd = [
        "cargo", "run", "--release",
        "--example", "rebound_parity_pythagorean",
        "-p", "apsis",
        "--",
        "--output", str(apsis_csv),
        "--horizon", str(args.horizon),
    ]
    rc = subprocess.call(apsis_cmd, cwd=workspace_root)
    if rc != 0:
        print(f"ERROR: apsis side failed with exit code {rc}", file=sys.stderr)
        return rc

    # ── 2. REBOUND side ─────────────────────────────────────────────────── #
    print()
    print("=" * 72)
    print("[2/3] REBOUND IAS15 (Python)")
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
        description="Run the full apsis/REBOUND parity comparison for the "
        "Pythagorean three-body problem.",
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
    parser.add_argument(
        "--horizon",
        type=float,
        default=70.0,
        help="Integration horizon in canonical time units. Default: 70 "
        "(the gated baseline declared in the protocol notebook). Larger "
        "values correspond to informational stress runs that are NOT "
        "part of the gated parity claim — see §Out of scope.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
