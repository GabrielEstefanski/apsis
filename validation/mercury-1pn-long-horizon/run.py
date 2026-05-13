"""Orchestrator for the long-horizon Mercury 1PN scenario.

Runs the IAS15+1PN side (Cargo example), optionally the Mercurius+1PN
side (also a Cargo example, available once PR #86 lands), then the
comparator. Exits with the comparator's exit code.

Usage from the scenario directory:

    python run.py                    # IAS15-only (Tier 1)
    python run.py --include-mercurius  # IAS15 + Mercurius (Tier 1 + 2 + 3)

Exit codes:

- 0 — all metrics within tolerance.
- 1 — input file error in the comparator.
- 2 — at least one metric exceeded tolerance.
- non-zero (other) — apsis side failed.

Protocol notebook:
    docs/experiments/2026-05-13-mercury-1pn-long-horizon.md
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCENARIO_DIR: Path = Path(__file__).resolve().parent
# scenario dir is `<root>/validation/mercury-1pn-long-horizon`, so the
# workspace root is two levels up.
DEFAULT_WORKSPACE_ROOT: Path = SCENARIO_DIR.parent.parent


def main() -> int:
    args = parse_args()
    workspace_root = Path(args.workspace_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ias15_csv = output_dir / "ias15.csv"
    mercurius_csv = output_dir / "mercurius.csv"

    # ── 1. IAS15+1PN side via cargo ──────────────────────────────────────
    print("=" * 72)
    print("[1/3] apsis IAS15 + apsis-1pn")
    print("=" * 72)
    rc = subprocess.call(
        [
            "cargo", "run", "--release",
            "--example", "mercury_1pn_long_horizon_ias15",
            "-p", "apsis-1pn",
            "--",
            "--output", str(ias15_csv),
        ],
        cwd=workspace_root,
    )
    if rc != 0:
        print(f"ERROR: IAS15 side failed with exit code {rc}", file=sys.stderr)
        return rc

    # ── 2. Mercurius+1PN side (optional) ────────────────────────────────
    if args.include_mercurius:
        print()
        print("=" * 72)
        print("[2/3] apsis Mercurius + apsis-1pn")
        print("=" * 72)
        rc = subprocess.call(
            [
                "cargo", "run", "--release",
                "--example", "mercury_1pn_long_horizon_mercurius",
                "-p", "apsis-1pn",
                "--",
                "--output", str(mercurius_csv),
            ],
            cwd=workspace_root,
        )
        if rc != 0:
            print(f"ERROR: Mercurius side failed with exit code {rc}", file=sys.stderr)
            return rc

    # ── 3. Comparator ────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("[3/3] Comparator")
    print("=" * 72)
    cmd = [
        sys.executable, str(SCENARIO_DIR / "compare.py"),
        "--ias15-csv", str(ias15_csv),
        "--output-dir", str(output_dir),
    ]
    if args.include_mercurius:
        cmd += ["--mercurius-csv", str(mercurius_csv)]
    return subprocess.call(cmd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the long-horizon Mercury 1PN comparison.",
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
        "--include-mercurius",
        action="store_true",
        help="Also run the Mercurius + apsis-1pn side and the cross-integrator "
        "Tier 2 + Tier 3 metrics. Available once PR #86 (Mercurius perturbation "
        "hole fix) merges.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
