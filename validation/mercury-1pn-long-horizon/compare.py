"""Comparator for the long-horizon Mercury 1PN scenario.

Loads `out/ias15.csv` (and optionally `out/mercurius.csv` once the
Mercurius half lands), computes the cumulative perihelion advance
trajectory for each, and gates against the closed-form Schwarzschild
test-particle prediction `Δω = 6πGM/(c²a(1−e²))` per orbit.

## Gated metrics (declared *a priori* in the protocol notebook)

Per-side (Tier 1 for IAS15, Tier 2 for Mercurius):

- `|Δω_measured(end) − Δω_GR(end)| / |Δω_GR(end)|` ≤ 1 × 10⁻⁵
- per-orbit linearity `R² ≥ 0.99999` against `Δω_GR_per_orbit · k`

Cross-integrator (Tier 3, when both CSVs are present):

- `|Δω_IAS15(end) − Δω_Mercurius(end)| / |Δω_GR(end)|` ≤ 5 × 10⁻⁵

## Exit codes

- 0 — all gated metrics within tolerance.
- 1 — input file error (missing file, malformed schema).
- 2 — at least one gated metric exceeded tolerance.

## Run

    python compare.py
    python compare.py --ias15-csv path/to/ias15.csv --mercurius-csv path/to/mercurius.csv

## Protocol notebook

    docs/experiments/2026-05-13-mercury-1pn-long-horizon.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# ── Protocol constants (mirror the Rust example + apsis-1pn) ──────────── #

A: float = 0.387098
E: float = 0.20563

# Speed of light in canonical-Hénon units (mirrors `apsis_1pn::C_SOLAR_UNITS`).
# Computed at compile time from SI constants in the Rust crate; reproduced
# here using the same SI definitions so the Python derivation is the source
# of truth for this harness, not a hand-transcribed literal.
C_SI: float = 299_792_458.0
AU_SI: float = 149_597_870_700.0
YEAR_S: float = 365.25 * 86_400.0
TWO_PI: float = 2.0 * math.pi
C_SOLAR_UNITS: float = C_SI * (YEAR_S / TWO_PI) / AU_SI

# Per-orbit GR perihelion advance for the locked Mercury IC.
DELTA_OMEGA_GR_PER_ORBIT: float = 6.0 * math.pi / (
    C_SOLAR_UNITS * C_SOLAR_UNITS * A * (1.0 - E * E)
)

# ── Tolerances declared a priori ───────────────────────────────────────── #

TOL_REL_OMEGA: float = 1.0e-5
TOL_LINEARITY_R2: float = 0.99999
TOL_CROSS_INTEGRATOR: float = 5.0e-5


# ── Data records ───────────────────────────────────────────────────────── #


@dataclass
class Sample:
    orbit: int
    t: float
    x: float
    y: float
    vx: float
    vy: float
    a_osc: float
    e_osc: float
    omega_osc: float


@dataclass
class MetricResult:
    name: str
    observed: float
    tolerance: float
    passed: bool
    detail: dict = field(default_factory=dict)


# ── Trajectory analysis ────────────────────────────────────────────────── #


def unwrap(angles: list[float]) -> list[float]:
    """Remove `2π` jumps from an angle trajectory in-place fashion."""
    if not angles:
        return []
    out = [angles[0]]
    for a in angles[1:]:
        prev = out[-1]
        d = a - prev
        while d > math.pi:
            d -= 2.0 * math.pi
        while d < -math.pi:
            d += 2.0 * math.pi
        out.append(prev + d)
    return out


def linearity_r_squared(measured: list[float], slope: float) -> float:
    """`R²` of `measured[k]` against the linear model `slope · k` with
    intercept fixed at zero (since `measured[0] = 0` after subtracting
    the initial value). The fit is the GR prediction; `R²` measures how
    well the measured trajectory matches that prediction beyond the
    fit's intercept-free constraint."""
    n = len(measured)
    if n < 3:
        return 1.0
    mean_m = sum(measured) / n
    ss_tot = sum((m - mean_m) ** 2 for m in measured)
    if ss_tot < 1.0e-30:
        return 1.0
    ss_res = sum((m - slope * k) ** 2 for k, m in enumerate(measured))
    return 1.0 - ss_res / ss_tot


def evaluate_side(samples: list[Sample], label: str) -> tuple[MetricResult, MetricResult, list[float]]:
    """Returns (`Δω_end relative error`, `linearity R²`, unwrapped Δω trajectory)."""
    omegas_unwrapped = unwrap([s.omega_osc for s in samples])
    omega_initial = omegas_unwrapped[0]
    delta_omega = [w - omega_initial for w in omegas_unwrapped]

    n_orbits = samples[-1].orbit
    predicted_end = DELTA_OMEGA_GR_PER_ORBIT * n_orbits
    measured_end = delta_omega[-1]
    rel_err = abs(measured_end - predicted_end) / abs(predicted_end)

    r2 = linearity_r_squared(delta_omega, DELTA_OMEGA_GR_PER_ORBIT)

    m_rel = MetricResult(
        name=f"{label}: Δω relative error vs GR (end)",
        observed=rel_err,
        tolerance=TOL_REL_OMEGA,
        passed=rel_err <= TOL_REL_OMEGA,
        detail={
            "measured_rad": measured_end,
            "predicted_rad": predicted_end,
            "n_orbits": n_orbits,
            "Δω_per_orbit_rad": DELTA_OMEGA_GR_PER_ORBIT,
        },
    )
    m_lin = MetricResult(
        name=f"{label}: per-orbit linearity R²",
        observed=r2,
        tolerance=TOL_LINEARITY_R2,
        passed=r2 >= TOL_LINEARITY_R2,
    )
    return m_rel, m_lin, delta_omega


# ── Main ───────────────────────────────────────────────────────────────── #


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    args = parse_args()
    ias15_path = Path(args.ias15_csv).resolve()
    mercurius_path = Path(args.mercurius_csv).resolve() if args.mercurius_csv else None
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        ias15 = load_csv(ias15_path)
    except FileNotFoundError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1

    metrics: list[MetricResult] = []

    # ── Tier 1 — IAS15 + apsis-1pn ───────────────────────────────────── #
    m_rel_i, m_lin_i, delta_i = evaluate_side(ias15, "IAS15")
    metrics.extend([m_rel_i, m_lin_i])

    # ── Tier 2 + 3 — Mercurius (optional) ────────────────────────────── #
    delta_m: list[float] | None = None
    if mercurius_path is not None and mercurius_path.exists():
        mercurius = load_csv(mercurius_path)
        if len(mercurius) != len(ias15):
            print(
                f"ERROR: sample count mismatch — ias15 has {len(ias15)}, "
                f"mercurius has {len(mercurius)}",
                file=sys.stderr,
            )
            return 1
        m_rel_m, m_lin_m, delta_m = evaluate_side(mercurius, "Mercurius")
        metrics.extend([m_rel_m, m_lin_m])

        n_orbits = ias15[-1].orbit
        predicted_end = DELTA_OMEGA_GR_PER_ORBIT * n_orbits
        cross_err = abs(delta_i[-1] - delta_m[-1]) / abs(predicted_end)
        m_cross = MetricResult(
            name="cross-integrator: |Δω_IAS15 − Δω_Mercurius| / |Δω_GR|",
            observed=cross_err,
            tolerance=TOL_CROSS_INTEGRATOR,
            passed=cross_err <= TOL_CROSS_INTEGRATOR,
        )
        metrics.append(m_cross)

    all_passed = all(m.passed for m in metrics)
    print_report(metrics, len(ias15), ias15[-1].t, delta_m is not None)
    write_json_report(output_dir, metrics, len(ias15), ias15[-1].t, all_passed)
    return 0 if all_passed else 2


# ── I/O ─────────────────────────────────────────────────────────────────── #


def load_csv(path: Path) -> list[Sample]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found at {path}")
    samples: list[Sample] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            samples.append(
                Sample(
                    orbit=int(row["orbit"]),
                    t=float(row["t"]),
                    x=float(row["x"]),
                    y=float(row["y"]),
                    vx=float(row["vx"]),
                    vy=float(row["vy"]),
                    a_osc=float(row["a_osc"]),
                    e_osc=float(row["e_osc"]),
                    omega_osc=float(row["omega_osc"]),
                )
            )
    return samples


def print_report(
    metrics: list[MetricResult],
    n_samples: int,
    t_final: float,
    cross_done: bool,
) -> None:
    print()
    print("Long-horizon Mercury 1PN — comparison report")
    print(f"  samples per side : {n_samples}")
    print(f"  t_final          : {t_final:.6e}  (canonical Hénon)")
    print(f"  Δω_GR per orbit  : {DELTA_OMEGA_GR_PER_ORBIT:.6e} rad")
    if not cross_done:
        print("  Mercurius CSV    : not provided (Tier 2 / Tier 3 skipped)")
    print()
    print(f"  {'metric (gated)':<60} {'observed':>14} {'tolerance':>14}  verdict")
    print(f"  {'-' * 60} {'-' * 14} {'-' * 14}  {'-' * 7}")
    for m in metrics:
        verdict = "pass" if m.passed else "FAIL"
        print(f"  {m.name:<60} {m.observed:>14.3e} {m.tolerance:>14.3e}  {verdict}")
    print()


def write_json_report(
    output_dir: Path,
    metrics: list[MetricResult],
    n_samples: int,
    t_final: float,
    all_passed: bool,
) -> None:
    report = {
        "all_passed": all_passed,
        "n_samples": n_samples,
        "t_final": t_final,
        "delta_omega_gr_per_orbit": DELTA_OMEGA_GR_PER_ORBIT,
        "metrics": [asdict(m) for m in metrics],
    }
    (output_dir / "comparison.json").write_text(json.dumps(report, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Comparator for the long-horizon Mercury 1PN scenario.",
    )
    parser.add_argument(
        "--ias15-csv",
        default="out/ias15.csv",
        help="Path to the IAS15-side CSV. Default: out/ias15.csv.",
    )
    parser.add_argument(
        "--mercurius-csv",
        default="",
        help="Path to the Mercurius-side CSV (optional; "
        "when present, Tier 2 + Tier 3 metrics are computed). "
        "Default: empty (Tier 1 only).",
    )
    parser.add_argument(
        "--output-dir",
        default="out",
        help="Directory for the JSON report. Default: out/.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
