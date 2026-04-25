"""Comparator for REBOUND parity — Kepler e=0.5.

Loads `out/apsis.csv` and `out/rebound.csv`, computes the metrics defined
in the **revised** protocol (see §Revised Methodology in the notebook),
and exits with code 0 iff every gated metric is within its *a priori*
tolerance. On failure, exits non-zero with a structured report saved
to `out/comparison.json`.

## Why orbital invariants, not |Δr|

Adaptive high-order integrators (IAS15) are not bit-deterministic
across independent implementations: ULP-level differences in the
controller's `dt` selection accumulate as orbital *phase* drift over
many periods. This phase drift is *not* a numerical-correctness
signal — both implementations stay on the same Kepler orbit, just
advancing along it at slightly different rates. Sampling at fixed
times and computing `|r_apsis(t) − r_rebound(t)|` therefore conflates
phase drift (not invariant across implementations) with geometric
drift (the actual physical signal).

The revised protocol gates on **orbital invariants** instead — `a`,
`e`, `ω`, `h`, `E`. These are physical constants of pure Kepler
motion; both implementations should agree on them to within their
respective conservation-precision floors. `|Δr|` is reported but
flagged informational; see notebook §Pilot Results and §Revised
Methodology for the derivation.

## Gated metrics (declared *a priori* in §Revised Methodology)

- `|Δa| / a` — fractional cross-implementation drift in semi-major axis
- `|Δe|`    — eccentricity drift
- `|Δω|`    — periapsis-orientation drift (radians, NOT phase)
- `|Δh| / h` — fractional angular-momentum drift
- `|ΔE/E_0|` per side — energy conservation per implementation
- cross-impl `|ΔE|/|E_0|` — energy agreement between implementations

## Exit codes

- 0 — all gated metrics within tolerance.
- 1 — input file error (missing file, sample count mismatch).
- 2 — at least one gated metric exceeded tolerance.

## Run

    python compare.py
    python compare.py --apsis-csv path/to/apsis.csv --rebound-csv path/to/rebound.csv

## Protocol notebook

    docs/experiments/2026-04-25-rebound-parity-kepler.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# ── Protocol constants (mirror the Rust + Python sides) ────────────────── #

M_PRIMARY: float = 1.0
M_SECONDARY: float = 1.0e-6
# Reduced-mass parameter for relative motion: μ = G(m_primary + m_secondary)
MU: float = M_PRIMARY + M_SECONDARY

# ── Tolerances declared a priori (revised protocol) ────────────────────── #
#
# All four gated orbital invariants are bounded at ~50× f64 machine epsilon
# (≈ 1e-13), reflecting the ULP-level cross-implementation noise floor for
# two correct IAS15 implementations conserving Kepler invariants. The
# periapsis-orientation tolerance is one decade wider because the
# eccentricity-vector → atan2(ω) path has a 1/|e| condition factor; for
# e=0.5 this is benign but justifies the extra margin.
#
# Energy tolerances are unchanged from the original protocol — they
# already passed comfortably in the pilot run.

TOL_RELATIVE_SEMIAXIS: float = 1.0e-13
TOL_ECCENTRICITY: float = 1.0e-13
TOL_PERIAPSIS_OMEGA: float = 1.0e-12
TOL_RELATIVE_ANGULAR_MOMENTUM: float = 1.0e-13
TOL_RELATIVE_ENERGY_DRIFT: float = 1.0e-13
TOL_CROSS_IMPL_ENERGY: float = 1.0e-13


# ── Data records ───────────────────────────────────────────────────────── #


@dataclass
class Sample:
    orbit: int
    t: float
    x0: float
    y0: float
    vx0: float
    vy0: float
    x1: float
    y1: float
    vx1: float
    vy1: float
    e_total: float


@dataclass
class OrbitalState:
    """`(a, e, ω, h)` of the secondary's orbit relative to the primary."""

    a: float
    e: float
    omega: float
    h: float


@dataclass
class MetricResult:
    name: str
    observed: float
    tolerance: float
    passed: bool
    detail: dict = field(default_factory=dict)


# ── Orbital element extraction ─────────────────────────────────────────── #


def relative_elements(s: Sample, mu: float) -> OrbitalState:
    """Compute orbital elements of the secondary relative to the primary.

    2D Kepler reduction:
        r⃗ = r₁ − r₀,    v⃗ = v₁ − v₀,    μ = G(m₀ + m₁)

        ε = ½ v² − μ/r        (specific energy)
        a = −μ / (2ε)         (semi-major axis)
        h = x vy − y vx       (specific angular momentum, z-component)
        e² = 1 − h²/(μ a)     (eccentricity)
        e_vec = ((v² − μ/r) r⃗ − (r⃗·v⃗) v⃗) / μ
        ω = atan2(e_vec_y, e_vec_x)   (orientation of periapsis)

    All four are constants of pure Kepler motion, so both implementations
    must agree on them to within their respective conservation-precision
    floors regardless of orbital phase.
    """
    rx = s.x1 - s.x0
    ry = s.y1 - s.y0
    vrx = s.vx1 - s.vx0
    vry = s.vy1 - s.vy0

    r = math.sqrt(rx * rx + ry * ry)
    v_sq = vrx * vrx + vry * vry

    eps = 0.5 * v_sq - mu / r
    a = -mu / (2.0 * eps)
    h = rx * vry - ry * vrx
    # Clamp away from tiny negative values caused by f64 cancellation.
    e_sq = max(0.0, 1.0 - h * h / (mu * a))
    e = math.sqrt(e_sq)
    r_dot_v = rx * vrx + ry * vry
    ev_x = ((v_sq - mu / r) * rx - r_dot_v * vrx) / mu
    ev_y = ((v_sq - mu / r) * ry - r_dot_v * vry) / mu
    omega = math.atan2(ev_y, ev_x)

    return OrbitalState(a=a, e=e, omega=omega, h=h)


def angle_diff(a: float, b: float) -> float:
    """Smallest signed difference `a − b` wrapped into [−π, π]."""
    d = math.fmod(a - b + math.pi, 2.0 * math.pi)
    if d < 0.0:
        d += 2.0 * math.pi
    return d - math.pi


# ── Main ───────────────────────────────────────────────────────────────── #


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    args = parse_args()
    apsis_path = Path(args.apsis_csv).resolve()
    rebound_path = Path(args.rebound_csv).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        apsis = load_csv(apsis_path)
        rebound = load_csv(rebound_path)
    except FileNotFoundError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1

    if len(apsis) != len(rebound):
        print(
            f"ERROR: sample count mismatch — apsis has {len(apsis)}, "
            f"rebound has {len(rebound)}",
            file=sys.stderr,
        )
        return 1

    # ── Orbital elements per side per sample ────────────────────────────── #
    elem_apsis = [relative_elements(s, MU) for s in apsis]
    elem_rebound = [relative_elements(s, MU) for s in rebound]
    a0 = elem_apsis[0].a
    h0 = elem_apsis[0].h

    # ── Gated metric 1: Δa / a ──────────────────────────────────────────── #
    max_da_rel = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_da_rel = max(max_da_rel, abs(ea.a - er.a) / abs(a0))
    m_a = MetricResult(
        name="|Δa|/a (semi-major axis)",
        observed=max_da_rel,
        tolerance=TOL_RELATIVE_SEMIAXIS,
        passed=max_da_rel <= TOL_RELATIVE_SEMIAXIS,
    )

    # ── Gated metric 2: Δe ─────────────────────────────────────────────── #
    max_de = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_de = max(max_de, abs(ea.e - er.e))
    m_e = MetricResult(
        name="|Δe| (eccentricity)",
        observed=max_de,
        tolerance=TOL_ECCENTRICITY,
        passed=max_de <= TOL_ECCENTRICITY,
    )

    # ── Gated metric 3: Δω ─────────────────────────────────────────────── #
    max_domega = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_domega = max(max_domega, abs(angle_diff(ea.omega, er.omega)))
    m_om = MetricResult(
        name="|Δω| (periapsis orient.)",
        observed=max_domega,
        tolerance=TOL_PERIAPSIS_OMEGA,
        passed=max_domega <= TOL_PERIAPSIS_OMEGA,
    )

    # ── Gated metric 4: Δh / h ─────────────────────────────────────────── #
    max_dh_rel = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_dh_rel = max(max_dh_rel, abs(ea.h - er.h) / abs(h0))
    m_h = MetricResult(
        name="|Δh|/h (angular momentum)",
        observed=max_dh_rel,
        tolerance=TOL_RELATIVE_ANGULAR_MOMENTUM,
        passed=max_dh_rel <= TOL_RELATIVE_ANGULAR_MOMENTUM,
    )

    # ── Gated metrics 5–7: energy ──────────────────────────────────────── #
    e0_apsis = apsis[0].e_total
    e0_rebound = rebound[0].e_total
    max_drift_apsis = max(abs(a.e_total - e0_apsis) / abs(e0_apsis) for a in apsis)
    max_drift_rebound = max(
        abs(r.e_total - e0_rebound) / abs(e0_rebound) for r in rebound
    )
    max_cross_e = max(
        abs(a.e_total - r.e_total) / abs(e0_apsis)
        for a, r in zip(apsis, rebound)
    )
    m_ea = MetricResult(
        name="|ΔE/E_0| apsis",
        observed=max_drift_apsis,
        tolerance=TOL_RELATIVE_ENERGY_DRIFT,
        passed=max_drift_apsis <= TOL_RELATIVE_ENERGY_DRIFT,
    )
    m_er = MetricResult(
        name="|ΔE/E_0| rebound",
        observed=max_drift_rebound,
        tolerance=TOL_RELATIVE_ENERGY_DRIFT,
        passed=max_drift_rebound <= TOL_RELATIVE_ENERGY_DRIFT,
    )
    m_cross_e = MetricResult(
        name="cross-impl |ΔE|/|E_0|",
        observed=max_cross_e,
        tolerance=TOL_CROSS_IMPL_ENERGY,
        passed=max_cross_e <= TOL_CROSS_IMPL_ENERGY,
    )

    metrics = [m_a, m_e, m_om, m_h, m_ea, m_er, m_cross_e]
    all_passed = all(m.passed for m in metrics)

    # ── Informational (NOT gated): point-wise |Δr| on the secondary ──── #
    max_dr = 0.0
    max_dr_orbit = 0
    for a, r in zip(apsis, rebound):
        dx = a.x1 - r.x1
        dy = a.y1 - r.y1
        dr = (dx * dx + dy * dy) ** 0.5
        if dr > max_dr:
            max_dr = dr
            max_dr_orbit = a.orbit
    info = MetricResult(
        name="|Δr| (secondary, NOT gated)",
        observed=max_dr,
        tolerance=float("inf"),
        passed=True,
        detail={
            "argmax_orbit": max_dr_orbit,
            "note": "phase-drift contaminated; not invariant across adaptive integrators",
        },
    )

    print_report(metrics, info, len(apsis), apsis[-1].t, e0_apsis, e0_rebound)
    write_json_report(output_dir, metrics, info, len(apsis), apsis[-1].t, all_passed)

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
                    x0=float(row["x0"]),
                    y0=float(row["y0"]),
                    vx0=float(row["vx0"]),
                    vy0=float(row["vy0"]),
                    x1=float(row["x1"]),
                    y1=float(row["y1"]),
                    vx1=float(row["vx1"]),
                    vy1=float(row["vy1"]),
                    e_total=float(row["e_total"]),
                )
            )
    return samples


def print_report(
    metrics: list[MetricResult],
    info: MetricResult,
    n_samples: int,
    t_final: float,
    e0_apsis: float,
    e0_rebound: float,
) -> None:
    print()
    print("REBOUND parity — Kepler e=0.5 — comparison report (revised metrics)")
    print(f"  samples : {n_samples}")
    print(f"  t_final : {t_final:.6e}")
    print(f"  E_0 apsis   : {e0_apsis:+.18e}")
    print(f"  E_0 rebound : {e0_rebound:+.18e}")
    print()
    print(f"  {'metric (gated)':<32} {'observed':>14} {'tolerance':>14}  verdict")
    print(f"  {'-' * 32} {'-' * 14} {'-' * 14}  {'-' * 7}")
    for m in metrics:
        verdict = "pass" if m.passed else "FAIL"
        print(f"  {m.name:<32} {m.observed:>14.3e} {m.tolerance:>14.3e}  {verdict}")
    print()
    print(f"  {'informational (not gated)':<32} {'observed':>14}")
    print(f"  {'-' * 32} {'-' * 14}")
    print(f"  {info.name:<32} {info.observed:>14.3e}")
    argmax = info.detail.get("argmax_orbit", "n/a")
    print(f"    └── peak at orbit {argmax}; phase-drift contaminated, see protocol notebook")
    print()


def write_json_report(
    output_dir: Path,
    metrics: list[MetricResult],
    info: MetricResult,
    n_samples: int,
    t_final: float,
    all_passed: bool,
) -> None:
    report = {
        "all_passed": all_passed,
        "n_samples": n_samples,
        "t_final": t_final,
        "metrics": [asdict(m) for m in metrics],
        "informational": asdict(info),
    }
    (output_dir / "comparison.json").write_text(json.dumps(report, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Comparator for REBOUND vs apsis parity — Kepler e=0.5.",
    )
    parser.add_argument(
        "--apsis-csv",
        default="out/apsis.csv",
        help="Path to the apsis-side CSV. Default: out/apsis.csv.",
    )
    parser.add_argument(
        "--rebound-csv",
        default="out/rebound.csv",
        help="Path to the REBOUND-side CSV. Default: out/rebound.csv.",
    )
    parser.add_argument(
        "--output-dir",
        default="out",
        help="Directory for the JSON report. Default: out/.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
