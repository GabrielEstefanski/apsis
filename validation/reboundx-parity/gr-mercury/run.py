"""Orchestrator for REBOUNDx parity — Sun-Mercury 1PN.

Runs both sides and the comparator for the gr (1PN) case and the 1PN-off
control: apsis Cargo example -> reboundx_side.py -> compare.py.

Must run in a Linux environment (WSL / container / CI): REBOUNDx does not
build on Windows/MSVC (C99 VLAs in gr_full.c). The reboundx venv must have
reboundx compiled fresh against its co-installed rebound (`pip install
--no-cache-dir reboundx`) or the librebound RPATH breaks — see requirements.txt.

Run from the scenario directory with the reboundx venv's python:
    python run.py
    python run.py --workspace-root /path/to/apsis

Exit: 0 if both gr and control gate clean; non-zero otherwise.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCENARIO_DIR = Path(__file__).resolve().parent
DEFAULT_WORKSPACE_ROOT = SCENARIO_DIR.parents[2]


def run(cmd: list[str], cwd: Path) -> int:
    print(f"\n$ {' '.join(cmd)}  (cwd={cwd})", flush=True)
    return subprocess.run(cmd, cwd=cwd).returncode


def main() -> int:
    args = parse_args()
    root = Path(args.workspace_root).resolve()
    py = sys.executable
    out = SCENARIO_DIR / "out"

    steps: list[tuple[list[str], Path]] = [
        # apsis side (Cargo) — 1PN on, then control.
        (["cargo", "run", "--release", "--example", "reboundx_parity_gr", "-p", "apsis-1pn",
          "--", "--output", str(out / "apsis.csv")], root),
        (["cargo", "run", "--release", "--example", "reboundx_parity_gr", "-p", "apsis-1pn",
          "--", "--no-1pn", "--output", str(out / "apsis_kepler.csv")], root),
        # reboundx side.
        ([py, "reboundx_side.py", "--apsis-csv", "out/apsis.csv",
          "--output", "out/rebound.csv"], SCENARIO_DIR),
        ([py, "reboundx_side.py", "--no-1pn", "--apsis-csv", "out/apsis_kepler.csv",
          "--output", "out/rebound_kepler.csv"], SCENARIO_DIR),
    ]
    for cmd, cwd in steps:
        rc = run(cmd, cwd)
        if rc != 0:
            print(f"step failed (rc={rc})", file=sys.stderr)
            return rc

    rc_gr = run([py, "compare.py", "--apsis", "out/apsis.csv",
                 "--rebound", "out/rebound.csv", "--mode", "gr"], SCENARIO_DIR)
    rc_ctl = run([py, "compare.py", "--apsis", "out/apsis_kepler.csv",
                  "--rebound", "out/rebound_kepler.csv", "--mode", "control"], SCENARIO_DIR)
    return rc_gr or rc_ctl


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="REBOUNDx parity orchestrator (Sun-Mercury 1PN).")
    p.add_argument("--workspace-root", default=str(DEFAULT_WORKSPACE_ROOT))
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(main())
