"""Mercury 1PN error-budget ensemble orchestrator (Phase B).

Runs `crates/apsis-1pn/examples/error_budget_run.rs` via
`cargo run --release --example error_budget_run -p apsis-1pn -- ...`
from the workspace root.

Usage
-----
    python ensemble.py --phase smoke
    python ensemble.py --phase b1
    python ensemble.py --phase b3

Phases
------
smoke : 3 ULPs {0,1,2} x both constructors at N=500.
        Validates plumbing and prints a table.  Does not write CSV.

b1    : ULPs 0..24 (25 runs) x both constructors at N=500.
        Writes out/b1.csv.  Prints per-constructor mean rel_err,
        sigma_omega (std of rel_err), and the B5 central values
        (ulp=0 runs) for both constructors.

b3    : ULPs 0..24 x N in {100, 250, 500, 1000, 2000}, constructor raw_c only.
        Writes out/b3.csv.  Fits log(sigma) vs log(N) -> exponent alpha
        with standard error.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import time
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent.resolve()
# Walk up to the workspace root (contains Cargo.toml at the top level).
WORKSPACE_ROOT = SCRIPT_DIR.parent.parent.resolve()
OUT_DIR = SCRIPT_DIR / "out"

# ── Constants ─────────────────────────────────────────────────────────────────

CONSTRUCTORS = ["raw_c", "for_units"]
B1_ULPS = list(range(25))          # 0 .. 24
B1_N_ORBITS = 500
B3_ULPS = list(range(25))          # 0 .. 24
B3_N_VALUES = [100, 250, 500, 1000, 2000]


# ── Runner ────────────────────────────────────────────────────────────────────

def run_one(n_orbits: int, ulp: int, constructor: str) -> dict:
    """Execute one error_budget_run and return the parsed result dict."""
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "error_budget_run",
        "-p",
        "apsis-1pn",
        "--",
        "--orbits",
        str(n_orbits),
        "--ulp",
        str(ulp),
        "--constructor",
        constructor,
    ]
    t0 = time.perf_counter()
    result = subprocess.run(
        cmd,
        cwd=str(WORKSPACE_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"[ERROR] cargo run failed (exit {result.returncode}):", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)

    # Extract the one CSV line from stdout (cargo may print noise to stderr).
    csv_line = result.stdout.strip()
    if not csv_line:
        print("[ERROR] no output from error_budget_run", file=sys.stderr)
        print("stderr:", result.stderr, file=sys.stderr)
        sys.exit(1)

    parts = csv_line.split(",")
    if len(parts) != 6:
        print(
            f"[ERROR] expected 6 CSV fields, got {len(parts)}: {csv_line!r}",
            file=sys.stderr,
        )
        sys.exit(1)

    return {
        "orbits": int(parts[0]),
        "ulp": int(parts[1]),
        "constructor": parts[2],
        "measured_rad": float(parts[3]),
        "predicted_rad": float(parts[4]),
        "rel_err": float(parts[5]),
        "elapsed_s": elapsed,
    }


# ── Phase: smoke ──────────────────────────────────────────────────────────────

def phase_smoke() -> None:
    """3 ULPs {0,1,2} x both constructors at N=500."""
    print("Phase smoke: 3 ULPs x 2 constructors at N=500")
    print()
    # Header
    print(
        f"{'ulp':>4}  {'constructor':>12}  {'measured_rad':>22}  "
        f"{'predicted_rad':>22}  {'rel_err':>12}  {'time_s':>7}"
    )
    print("-" * 90)

    rows = []
    for constructor in CONSTRUCTORS:
        for ulp in [0, 1, 2]:
            row = run_one(B1_N_ORBITS, ulp, constructor)
            rows.append(row)
            print(
                f"{row['ulp']:>4}  {row['constructor']:>12}  "
                f"{row['measured_rad']:>22.17e}  "
                f"{row['predicted_rad']:>22.17e}  "
                f"{row['rel_err']:>12.6e}  "
                f"{row['elapsed_s']:>7.1f}s"
            )

    print()
    print("Gate check (ulp=0, raw_c rel_err < 1e-4):")
    raw_c_central = next(r for r in rows if r["constructor"] == "raw_c" and r["ulp"] == 0)
    status = "PASS" if raw_c_central["rel_err"] < 1e-4 else "BLOCKED"
    print(f"  raw_c ulp=0 rel_err = {raw_c_central['rel_err']:.6e}  [{status}]")
    if status == "BLOCKED":
        print("BLOCKED: raw_c central value exceeds gate bound 1e-4", file=sys.stderr)
        sys.exit(2)


# ── Phase: b1 ─────────────────────────────────────────────────────────────────

def phase_b1() -> None:
    """25 ULPs x both constructors at N=500. Writes out/b1.csv."""
    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / "b1.csv"

    print(f"Phase b1: {len(B1_ULPS)} ULPs x {len(CONSTRUCTORS)} constructors at N={B1_N_ORBITS}")
    print(f"Output: {out_path}")
    print()

    fieldnames = ["orbits", "ulp", "constructor", "measured_rad", "predicted_rad", "rel_err"]
    rows: list[dict] = []

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for constructor in CONSTRUCTORS:
            for ulp in B1_ULPS:
                row = run_one(B1_N_ORBITS, ulp, constructor)
                print(
                    f"  {constructor:>12}  ulp={ulp:>2}  rel_err={row['rel_err']:.6e}  "
                    f"({row['elapsed_s']:.1f}s)"
                )
                rows.append(row)
                writer.writerow({k: row[k] for k in fieldnames})
                fh.flush()

    print()
    for constructor in CONSTRUCTORS:
        subset = [r for r in rows if r["constructor"] == constructor]
        rel_errs = [r["rel_err"] for r in subset]
        mean_re = sum(rel_errs) / len(rel_errs)
        sigma = math.sqrt(
            sum((x - mean_re) ** 2 for x in rel_errs) / max(len(rel_errs) - 1, 1)
        )
        central = next(r for r in subset if r["ulp"] == 0)
        print(
            f"  {constructor:>12}:  mean_rel_err={mean_re:.6e}  "
            f"sigma_omega={sigma:.6e}  "
            f"B5_central(ulp=0)={central['rel_err']:.6e}"
        )


# ── Phase: b3 ─────────────────────────────────────────────────────────────────

def _fit_log_log(xs: list[float], ys: list[float]) -> tuple[float, float]:
    """Ordinary least squares on log(y) ~ alpha*log(x) + intercept.

    Returns (alpha, stderr_alpha).
    """
    n = len(xs)
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    mean_lx = sum(lx) / n
    mean_ly = sum(ly) / n
    sxx = sum((x - mean_lx) ** 2 for x in lx)
    sxy = sum((lx[i] - mean_lx) * (ly[i] - mean_ly) for i in range(n))
    alpha = sxy / sxx
    # Residual standard error
    ly_hat = [mean_ly + alpha * (lx[i] - mean_lx) for i in range(n)]
    sse = sum((ly[i] - ly_hat[i]) ** 2 for i in range(n))
    stderr = math.sqrt(sse / max(n - 2, 1) / sxx) if n > 2 else float("nan")
    return alpha, stderr


def phase_b3() -> None:
    """25 ULPs x N in B3_N_VALUES, raw_c only. Writes out/b3.csv."""
    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / "b3.csv"
    constructor = "raw_c"

    total = len(B3_ULPS) * len(B3_N_VALUES)
    print(
        f"Phase b3: {len(B3_ULPS)} ULPs x {len(B3_N_VALUES)} N-values, "
        f"constructor={constructor} ({total} runs total)"
    )
    print(f"Output: {out_path}")
    print()

    fieldnames = ["orbits", "ulp", "constructor", "measured_rad", "predicted_rad", "rel_err"]
    rows: list[dict] = []

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for n_orbits in B3_N_VALUES:
            for ulp in B3_ULPS:
                row = run_one(n_orbits, ulp, constructor)
                print(
                    f"  N={n_orbits:>5}  ulp={ulp:>2}  rel_err={row['rel_err']:.6e}  "
                    f"({row['elapsed_s']:.1f}s)"
                )
                rows.append(row)
                writer.writerow({k: row[k] for k in fieldnames})
                fh.flush()

    print()
    print("H3 fit: log(sigma_omega) ~ alpha * log(N)")
    n_values = []
    sigmas = []
    for n_orbits in B3_N_VALUES:
        subset = [r for r in rows if r["orbits"] == n_orbits]
        rel_errs = [r["rel_err"] for r in subset]
        mean_re = sum(rel_errs) / len(rel_errs)
        sigma = math.sqrt(
            sum((x - mean_re) ** 2 for x in rel_errs) / max(len(rel_errs) - 1, 1)
        )
        print(f"  N={n_orbits:>5}  sigma_omega={sigma:.6e}  (mean_rel_err={mean_re:.6e})")
        if sigma > 0:
            n_values.append(float(n_orbits))
            sigmas.append(sigma)

    if len(n_values) >= 3:
        alpha, stderr = _fit_log_log(n_values, sigmas)
        print()
        print(f"  Fitted exponent alpha = {alpha:.4f} ± {stderr:.4f}")
        print(
            "  (expected: alpha ~ 1/2 = random walk; "
            "~1 = truncation-dominated; ~0 = bounded)"
        )
    else:
        print("  [WARN] too few non-zero sigma values for a reliable fit")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mercury 1PN error-budget ensemble (Phase B)."
    )
    parser.add_argument(
        "--phase",
        required=True,
        choices=["smoke", "b1", "b3"],
        help="Which phase to run.",
    )
    args = parser.parse_args()

    if args.phase == "smoke":
        phase_smoke()
    elif args.phase == "b1":
        phase_b1()
    elif args.phase == "b3":
        phase_b3()


if __name__ == "__main__":
    main()
