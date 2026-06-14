"""Mercury 1PN error-budget ensemble orchestrator (Phase B / B').

Runs `crates/apsis-1pn/examples/error_budget_run.rs` via
`cargo run --release --example error_budget_run -p apsis-1pn -- ...`
from the workspace root.

Usage
-----
    python ensemble.py --phase smoke
    python ensemble.py --phase b1
    python ensemble.py --phase b3
    python ensemble.py --phase b4

Phases
------
smoke : 3 ULPs {0,1,2} x both constructors at N=500.
        Validates plumbing and prints a table.  Does not write CSV.

b1    : ULPs 0..24 (25 runs) x both constructors at N=500.
        Writes out/b1.csv.  Per constructor: raw mean/sigma of the signed
        rel_err, B5 central values, and the Phase-B' Q-corrected budget
        (residual - Q(nu_end), compared against the Phase-A floors).

b3    : ULPs 0..24 x N in {100, 250, 500, 1000, 2000}, constructor raw_c only.
        Writes out/b3.csv.  H3 fit on raw sigma AND on Q-corrected sigma
        (the raw exponent is contaminated by endpoint sampling; the
        corrected one measures the actual integration-noise regime).

b4    : ULPs 0..11 x eps_b in {1e-7..1e-11} at N=500, raw_c only.
        Writes out/b4.csv.  Per eps_b: mean overshoot, raw and corrected
        residuals.  Truncation-vs-round-off discrimination: the raw offset
        must track the endpoint step size; the corrected one must sit at
        the floors, eps_b-independent.

Phase-B' correction
-------------------
`error_budget_endpoint_symbolic.py` (gates GB-h, GB0..GB5) derives the
endpoint-offset function

    Q(nu) = eps * (3 nu - (3/e - e) sin nu - (5/2) sin 2 nu)

with eps recovered per-row as predicted_rad / (6 pi N) — no constants
duplicated here. `integrate_until` exits up to one adaptive IAS15 step
past t_end, so the endpoint samples osculating omega at true anomaly
nu_end != 0; Q(nu_end) is the deterministic part of that sampling error
(verified parameter-free in `error_budget_endpoint_numerical.py`,
gates GN1, GN2, GN4).
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
B4_ULPS = list(range(12))          # 0 .. 11
B4_EPS_B = [1e-7, 1e-8, 1e-9, 1e-10, 1e-11]
B4_N_ORBITS = 500

# Mirrors the gate constant (mercury_precession_gate.rs); consumed only by Q.
E_MERCURY = 0.20563

# Derivation floors at Mercury's e=0.20563, relative units. A1 = k_osc(e)*eps
# (closed form); A2 = C(e)*q, C measured at e by error_budget_a2_eccentricity.py
# (strongly e-dependent — not the e=0.2 value 8.291).
FLOOR_A1 = -3.46428 * 2.662484e-8
FLOOR_A2 = +8.0617 * 1.660114e-7

FIELDNAMES = [
    "orbits", "ulp", "constructor", "eps_b",
    "measured_rad", "predicted_rad", "rel_err", "t_overshoot", "nu_end",
]


# ── Phase-B' endpoint correction ──────────────────────────────────────────────

def q_of_nu(eps: float, nu: float) -> float:
    """Derived endpoint-offset (rad): error_budget_endpoint_symbolic.py."""
    e = E_MERCURY
    return eps * (3.0 * nu - (3.0 / e - e) * math.sin(nu) - 2.5 * math.sin(2.0 * nu))


def enrich(row: dict) -> dict:
    """Attach derived per-row quantities: eps, angle residual, Q, corrected."""
    n = row["orbits"]
    eps = row["predicted_rad"] / (6.0 * math.pi * n)
    r_angle = row["measured_rad"] - row["predicted_rad"]
    q = q_of_nu(eps, row["nu_end"])
    row["eps"] = eps
    row["r_angle"] = r_angle
    row["q_pred"] = q
    row["corr_angle"] = r_angle - q
    row["corr_rel"] = (r_angle - q) / row["predicted_rad"]
    return row


# ── Runner ────────────────────────────────────────────────────────────────────

def run_one(n_orbits: int, ulp: int, constructor: str, eps_b: float | None = None) -> dict:
    """Execute one error_budget_run and return the parsed result dict."""
    cmd = [
        "cargo", "run", "--release", "--example", "error_budget_run",
        "-p", "apsis-1pn", "--",
        "--orbits", str(n_orbits),
        "--ulp", str(ulp),
        "--constructor", constructor,
    ]
    if eps_b is not None:
        cmd += ["--eps-b", f"{eps_b:e}"]
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
    if len(parts) != 9:
        print(
            f"[ERROR] expected 9 CSV fields, got {len(parts)}: {csv_line!r}",
            file=sys.stderr,
        )
        sys.exit(1)

    row = {
        "orbits": int(parts[0]),
        "ulp": int(parts[1]),
        "constructor": parts[2],
        "eps_b": float(parts[3]),
        "measured_rad": float(parts[4]),
        "predicted_rad": float(parts[5]),
        "rel_err": float(parts[6]),       # signed
        "t_overshoot": float(parts[7]),
        "nu_end": float(parts[8]),
        "elapsed_s": elapsed,
    }
    return enrich(row)


# ── Statistics helpers ────────────────────────────────────────────────────────

def mean_sigma(xs: list[float]) -> tuple[float, float]:
    m = sum(xs) / len(xs)
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / max(len(xs) - 1, 1))
    return m, s


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
    ly_hat = [mean_ly + alpha * (lx[i] - mean_lx) for i in range(n)]
    sse = sum((ly[i] - ly_hat[i]) ** 2 for i in range(n))
    stderr = math.sqrt(sse / max(n - 2, 1) / sxx) if n > 2 else float("nan")
    return alpha, stderr


def budget_summary(rows: list[dict], label: str) -> None:
    """Raw + Q-corrected ensemble summary for one (constructor, N) subset."""
    n = rows[0]["orbits"]
    pred = rows[0]["predicted_rad"]
    raw_m, raw_s = mean_sigma([r["rel_err"] for r in rows])
    ang_m, ang_s = mean_sigma([r["r_angle"] for r in rows])
    cor_m, cor_s = mean_sigma([r["corr_angle"] for r in rows])
    over_m, _ = mean_sigma([r["t_overshoot"] for r in rows])
    nu_m, _ = mean_sigma([r["nu_end"] for r in rows])
    k = len(rows)
    floors_angle = (FLOOR_A1 + FLOOR_A2) * pred
    print(f"  {label}: K={k}  N={n}")
    print(f"    raw   : rel_err mean={raw_m:+.3e}  sigma={raw_s:.3e}")
    print(f"    angle : mean={ang_m:+.3e} rad  sigma={ang_s:.3e} rad")
    print(f"    endpoint: <overshoot>={over_m:.3e}  <nu_end>={nu_m:+.4e} rad")
    print(
        f"    Q-corrected: mean={cor_m:+.3e} rad  sigma={cor_s:.3e} rad  "
        f"sem={cor_s / math.sqrt(k):.3e}"
    )
    print(
        f"    floors (A1+A2): {floors_angle:+.3e} rad   "
        f"unexplained: {cor_m - floors_angle:+.3e} rad "
        f"({(cor_m - floors_angle) / (cor_s / math.sqrt(k)):+.2f} sem)"
    )


# ── Phase: smoke ──────────────────────────────────────────────────────────────

def phase_smoke() -> None:
    """3 ULPs {0,1,2} x both constructors at N=500."""
    print("Phase smoke: 3 ULPs x 2 constructors at N=500")
    print()
    print(
        f"{'ulp':>4}  {'constructor':>12}  {'rel_err':>12}  {'overshoot':>11}  "
        f"{'nu_end':>11}  {'Q(nu)':>11}  {'corr_angle':>11}  {'time_s':>7}"
    )
    print("-" * 92)

    rows = []
    for constructor in CONSTRUCTORS:
        for ulp in [0, 1, 2]:
            row = run_one(B1_N_ORBITS, ulp, constructor)
            rows.append(row)
            print(
                f"{row['ulp']:>4}  {row['constructor']:>12}  "
                f"{row['rel_err']:>12.3e}  "
                f"{row['t_overshoot']:>11.3e}  "
                f"{row['nu_end']:>11.3e}  "
                f"{row['q_pred']:>11.3e}  "
                f"{row['corr_angle']:>11.3e}  "
                f"{row['elapsed_s']:>7.1f}s"
            )

    print()
    print("Gate check (ulp=0, raw_c |rel_err| < 1e-4):")
    raw_c_central = next(r for r in rows if r["constructor"] == "raw_c" and r["ulp"] == 0)
    status = "PASS" if abs(raw_c_central["rel_err"]) < 1e-4 else "BLOCKED"
    print(f"  raw_c ulp=0 rel_err = {raw_c_central['rel_err']:+.6e}  [{status}]")
    if status == "BLOCKED":
        print("BLOCKED: raw_c central value exceeds gate bound 1e-4", file=sys.stderr)
        sys.exit(2)


# ── CSV writer shared by b1/b3/b4 ─────────────────────────────────────────────

def write_rows(out_path: Path, grid, label_fn) -> list[dict]:
    OUT_DIR.mkdir(exist_ok=True)
    rows: list[dict] = []
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for n_orbits, ulp, constructor, eps_b in grid:
            row = run_one(n_orbits, ulp, constructor, eps_b)
            print(f"  {label_fn(row)}  ({row['elapsed_s']:.1f}s)", flush=True)
            rows.append(row)
            writer.writerow({k: row[k] for k in FIELDNAMES})
            fh.flush()
    return rows


# ── Phase: b1 ─────────────────────────────────────────────────────────────────

def phase_b1() -> None:
    """25 ULPs x both constructors at N=500. Writes out/b1.csv."""
    out_path = OUT_DIR / "b1.csv"
    print(f"Phase b1: {len(B1_ULPS)} ULPs x {len(CONSTRUCTORS)} constructors at N={B1_N_ORBITS}")
    print(f"Output: {out_path}")
    print()

    grid = [
        (B1_N_ORBITS, ulp, constructor, None)
        for constructor in CONSTRUCTORS
        for ulp in B1_ULPS
    ]
    rows = write_rows(
        out_path, grid,
        lambda r: f"{r['constructor']:>12}  ulp={r['ulp']:>2}  rel_err={r['rel_err']:+.6e}",
    )

    print()
    for constructor in CONSTRUCTORS:
        subset = [r for r in rows if r["constructor"] == constructor]
        budget_summary(subset, constructor)
        central = next(r for r in subset if r["ulp"] == 0)
        print(f"    B5 central (ulp=0): rel_err={central['rel_err']:+.6e}")


# ── Phase: b3 ─────────────────────────────────────────────────────────────────

def phase_b3() -> None:
    """25 ULPs x N in B3_N_VALUES, raw_c only. Writes out/b3.csv."""
    out_path = OUT_DIR / "b3.csv"
    constructor = "raw_c"
    total = len(B3_ULPS) * len(B3_N_VALUES)
    print(
        f"Phase b3: {len(B3_ULPS)} ULPs x {len(B3_N_VALUES)} N-values, "
        f"constructor={constructor} ({total} runs total)"
    )
    print(f"Output: {out_path}")
    print()

    grid = [
        (n_orbits, ulp, constructor, None)
        for n_orbits in B3_N_VALUES
        for ulp in B3_ULPS
    ]
    rows = write_rows(
        out_path, grid,
        lambda r: f"N={r['orbits']:>5}  ulp={r['ulp']:>2}  rel_err={r['rel_err']:+.6e}",
    )

    print()
    print("Per-N summaries (raw vs Q-corrected):")
    sig_raw, sig_cor, n_vals = [], [], []
    for n_orbits in B3_N_VALUES:
        subset = [r for r in rows if r["orbits"] == n_orbits]
        budget_summary(subset, f"N={n_orbits}")
        _, s_raw = mean_sigma([r["r_angle"] for r in subset])
        _, s_cor = mean_sigma([r["corr_angle"] for r in subset])
        if s_raw > 0 and s_cor > 0:
            n_vals.append(float(n_orbits))
            sig_raw.append(s_raw)
            sig_cor.append(s_cor)

    print()
    print("H3 fit: log(sigma_angle) ~ alpha * log(N)   [angle units, rad]")
    a_raw, e_raw = _fit_log_log(n_vals, sig_raw)
    a_cor, e_cor = _fit_log_log(n_vals, sig_cor)
    print(f"  raw        : alpha = {a_raw:+.4f} ± {e_raw:.4f}   (endpoint-sampling contaminated)")
    print(f"  Q-corrected: alpha = {a_cor:+.4f} ± {e_cor:.4f}")
    print(
        "  (angle-units key: ~ +1 coherent/truncation, ~ +1/2 round-off random "
        "walk [Brouwer 1937], ~ 0 bounded)"
    )


# ── Phase: b4 ─────────────────────────────────────────────────────────────────

def phase_b4() -> None:
    """12 ULPs x eps_b in B4_EPS_B at N=500, raw_c only. Writes out/b4.csv."""
    out_path = OUT_DIR / "b4.csv"
    constructor = "raw_c"
    total = len(B4_ULPS) * len(B4_EPS_B)
    print(
        f"Phase b4: {len(B4_ULPS)} ULPs x {len(B4_EPS_B)} eps_b values at "
        f"N={B4_N_ORBITS}, constructor={constructor} ({total} runs total)"
    )
    print(f"Output: {out_path}")
    print()

    grid = [
        (B4_N_ORBITS, ulp, constructor, eps_b)
        for eps_b in B4_EPS_B
        for ulp in B4_ULPS
    ]
    rows = write_rows(
        out_path, grid,
        lambda r: (
            f"eps_b={r['eps_b']:.0e}  ulp={r['ulp']:>2}  "
            f"rel_err={r['rel_err']:+.6e}  overshoot={r['t_overshoot']:.3e}"
        ),
    )

    print()
    print("Per-eps_b summaries (raw must track the endpoint step; corrected must plateau):")
    for eps_b in B4_EPS_B:
        subset = [r for r in rows if r["eps_b"] == eps_b]
        budget_summary(subset, f"eps_b={eps_b:.0e}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mercury 1PN error-budget ensemble (Phase B / B')."
    )
    parser.add_argument(
        "--phase",
        required=True,
        choices=["smoke", "b1", "b3", "b4"],
        help="Which phase to run.",
    )
    args = parser.parse_args()

    if args.phase == "smoke":
        phase_smoke()
    elif args.phase == "b1":
        phase_b1()
    elif args.phase == "b3":
        phase_b3()
    elif args.phase == "b4":
        phase_b4()


if __name__ == "__main__":
    main()
