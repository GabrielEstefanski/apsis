"""Comparator for REBOUND parity — Mercurius.

Loads `out/apsis.csv` and `out/rebound.csv`, computes the metrics defined
in §Hypothesis of the protocol notebook, and exits with code 0 iff every
gated metric is within its *a priori* tolerance. On failure, exits
non-zero with a structured report saved to `out/comparison.json`.

## Why orbital invariants for the test particle

Mercurius's encounter step engages IAS15 internally on the encountering
subset; IAS15 is not bit-deterministic across independent
implementations (the controller's `dt_next` selection branches on
ULP-scale truncation differences). At the end of a 10⁴-year integration
where the test particle has crossed Jupiter's orbit hundreds of times,
the cross-implementation phase drift on the encountering particles
dominates the point-wise position drift Δr. Orbital elements (`a`, `e`,
`i`) are constants of pure Kepler motion in the no-encounter regime
and change deterministically per encounter — both implementations
should agree on them to within their respective conservation-precision
floors. See the protocol notebook §Methodology for the full derivation.

## Gated metrics (declared *a priori* in the protocol notebook)

Tier 1 — conservation parity (3 gates):
- per-side ΔE/E₀ peak (apsis, REBOUND) — 2nd-order method floor
- cross-impl ΔE/E₀ peak — independent-implementations agreement
- cross-impl ΔLz/Lz₀ peak — angular momentum agreement

Tier 2 — test-particle orbital element parity (3 gates, end-of-integration):
- Δa relative
- Δe relative
- Δi relative

## Exit codes

- 0 — all gated metrics within tolerance.
- 1 — input file error (missing file, sample count mismatch).
- 2 — at least one gated metric exceeded tolerance.

## Run

    python compare.py
    python compare.py --apsis-csv path/to/apsis.csv --rebound-csv path/to/rebound.csv

## Protocol notebook

    paper/notebooks/2026-05-13-rebound-parity-mercurius.md
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

M_SUN: float = 1.0
M_TEST: float = 1.0e-9
N_BODIES: int = 6
TEST_PARTICLE_INDEX: int = 5  # last body in the IC list

# G in solar AU-year units (must match rebound_side.py exactly; same
# multiply-then-divide order as apsis's `UnitSystem::solar().g()` so
# the f64 representation agrees bit-for-bit across all three sides).
G_SI: float = 6.67430e-11
AU_M: float = 1.495978707e11
YR_S: float = 3.15576e7
GM_SUN_SI: float = 1.3271244e20  # IAU 2015 nominal; apsis units.rs GM_SUN_SI
MSUN_KG: float = GM_SUN_SI / G_SI
G_SOLAR: float = G_SI * MSUN_KG * YR_S * YR_S / (AU_M * AU_M * AU_M)

# Reduced-mass parameter for relative motion of TP around the Sun:
# μ = G(m_Sun + m_TP) ≈ G·m_Sun for m_TP ≪ M_Sun.
MU_TP: float = G_SOLAR * (M_SUN + M_TEST)

# ── Tolerances declared a priori (protocol notebook §Hypothesis) ───────── #

TOL_PER_SIDE_ENERGY: float = 1.0e-8
TOL_CROSS_IMPL_ENERGY: float = 5.0e-9
TOL_CROSS_IMPL_LZ: float = 1.0e-10
TOL_TP_SEMIAXIS_REL: float = 1.0e-5
TOL_TP_ECCENTRICITY_REL: float = 1.0e-5
TOL_TP_INCLINATION_REL: float = 1.0e-5


# ── Data records ───────────────────────────────────────────────────────── #


@dataclass
class Sample:
    """One row of the wide CSV. Body state is a flat list of length
    `6 · N_BODIES` in (x, y, z, vx, vy, vz) order per body."""

    year: int
    t: float
    body_state: list[float]
    e_total: float
    lz_total: float

    def body_pos(self, idx: int) -> tuple[float, float, float]:
        base = 6 * idx
        return (self.body_state[base], self.body_state[base + 1], self.body_state[base + 2])

    def body_vel(self, idx: int) -> tuple[float, float, float]:
        base = 6 * idx
        return (
            self.body_state[base + 3],
            self.body_state[base + 4],
            self.body_state[base + 5],
        )


@dataclass
class OrbitalState:
    """Osculating orbital elements `(a, e, i)` of the test particle
    relative to the Sun."""

    a: float
    e: float
    i: float


@dataclass
class MetricResult:
    name: str
    observed: float
    tolerance: float
    passed: bool
    detail: dict = field(default_factory=dict)


# ── Orbital element extraction (3D) ────────────────────────────────────── #


def test_particle_elements(s: Sample, mu: float) -> OrbitalState:
    """Compute osculating elements of body `TEST_PARTICLE_INDEX` relative
    to body 0 (the Sun) in 3D.

    Standard Kepler reduction (Murray & Dermott §2.8):
        r⃗ = r_TP − r_Sun,    v⃗ = v_TP − v_Sun,    μ = G(m_Sun + m_TP)
        ε = ½ v² − μ/r              (specific energy)
        a = −μ / (2ε)               (semi-major axis)
        h⃗ = r⃗ × v⃗                  (specific angular momentum)
        e⃗ = (v⃗ × h⃗)/μ − r⃗/r        (eccentricity vector)
        i = acos(h_z / |h|)         (inclination)
    """
    rsx, rsy, rsz = s.body_pos(0)
    vsx, vsy, vsz = s.body_vel(0)
    rpx, rpy, rpz = s.body_pos(TEST_PARTICLE_INDEX)
    vpx, vpy, vpz = s.body_vel(TEST_PARTICLE_INDEX)

    rx, ry, rz = rpx - rsx, rpy - rsy, rpz - rsz
    vx, vy, vz = vpx - vsx, vpy - vsy, vpz - vsz

    r = math.sqrt(rx * rx + ry * ry + rz * rz)
    v_sq = vx * vx + vy * vy + vz * vz

    eps = 0.5 * v_sq - mu / r
    a = -mu / (2.0 * eps)

    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    h_mag = math.sqrt(hx * hx + hy * hy + hz * hz)

    # Eccentricity vector e⃗ = (v⃗ × h⃗)/μ − r⃗/r
    evx = (vy * hz - vz * hy) / mu - rx / r
    evy = (vz * hx - vx * hz) / mu - ry / r
    evz = (vx * hy - vy * hx) / mu - rz / r
    e = math.sqrt(evx * evx + evy * evy + evz * evz)

    # Inclination from h_z / |h|, clamped against floating-point overshoot.
    cos_i = hz / h_mag if h_mag > 0.0 else 1.0
    cos_i = max(-1.0, min(1.0, cos_i))
    i = math.acos(cos_i)

    return OrbitalState(a=a, e=e, i=i)


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
        rebound_samples = load_csv(rebound_path)
    except FileNotFoundError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1

    if len(apsis) != len(rebound_samples):
        print(
            f"ERROR: sample count mismatch — apsis has {len(apsis)}, "
            f"rebound has {len(rebound_samples)}",
            file=sys.stderr,
        )
        return 1

    # ── Tier 1 — conservation parity ────────────────────────────────────── #
    e0_apsis = apsis[0].e_total
    e0_rebound = rebound_samples[0].e_total
    lz0_apsis = apsis[0].lz_total
    # Use the apsis-side baseline for all relative comparisons so the
    # report numbers are comparable across runs.
    e0_ref = abs(e0_apsis)
    lz0_ref = max(abs(lz0_apsis), 1.0e-30)

    max_drift_apsis = max(abs(s.e_total - e0_apsis) / e0_ref for s in apsis)
    max_drift_rebound = max(abs(s.e_total - e0_rebound) / abs(e0_rebound) for s in rebound_samples)
    max_cross_e = max(
        abs(a.e_total - r.e_total) / e0_ref for a, r in zip(apsis, rebound_samples)
    )
    max_cross_lz = max(
        abs(a.lz_total - r.lz_total) / lz0_ref for a, r in zip(apsis, rebound_samples)
    )

    m_e_apsis = MetricResult(
        name="ΔE/E₀ apsis (per side)",
        observed=max_drift_apsis,
        tolerance=TOL_PER_SIDE_ENERGY,
        passed=max_drift_apsis <= TOL_PER_SIDE_ENERGY,
    )
    m_e_rebound = MetricResult(
        name="ΔE/E₀ rebound (per side)",
        observed=max_drift_rebound,
        tolerance=TOL_PER_SIDE_ENERGY,
        passed=max_drift_rebound <= TOL_PER_SIDE_ENERGY,
    )
    m_cross_e = MetricResult(
        name="cross-impl ΔE/E₀",
        observed=max_cross_e,
        tolerance=TOL_CROSS_IMPL_ENERGY,
        passed=max_cross_e <= TOL_CROSS_IMPL_ENERGY,
    )
    m_cross_lz = MetricResult(
        name="cross-impl ΔLz/Lz₀",
        observed=max_cross_lz,
        tolerance=TOL_CROSS_IMPL_LZ,
        passed=max_cross_lz <= TOL_CROSS_IMPL_LZ,
    )

    # ── Tier 2 — test-particle orbital element parity (end-of-integration)
    elem_apsis_end = test_particle_elements(apsis[-1], MU_TP)
    elem_rebound_end = test_particle_elements(rebound_samples[-1], MU_TP)

    da_rel = abs(elem_apsis_end.a - elem_rebound_end.a) / max(abs(elem_apsis_end.a), 1.0e-30)
    de_rel = abs(elem_apsis_end.e - elem_rebound_end.e) / max(abs(elem_apsis_end.e), 1.0e-30)
    di_rel = abs(elem_apsis_end.i - elem_rebound_end.i) / max(abs(elem_apsis_end.i), 1.0e-30)

    m_a = MetricResult(
        name="Δa/a (TP, end of run)",
        observed=da_rel,
        tolerance=TOL_TP_SEMIAXIS_REL,
        passed=da_rel <= TOL_TP_SEMIAXIS_REL,
        detail={"apsis": elem_apsis_end.a, "rebound": elem_rebound_end.a},
    )
    m_e = MetricResult(
        name="Δe/e (TP, end of run)",
        observed=de_rel,
        tolerance=TOL_TP_ECCENTRICITY_REL,
        passed=de_rel <= TOL_TP_ECCENTRICITY_REL,
        detail={"apsis": elem_apsis_end.e, "rebound": elem_rebound_end.e},
    )
    m_i = MetricResult(
        name="Δi/i (TP, end of run)",
        observed=di_rel,
        tolerance=TOL_TP_INCLINATION_REL,
        passed=di_rel <= TOL_TP_INCLINATION_REL,
        detail={"apsis": elem_apsis_end.i, "rebound": elem_rebound_end.i},
    )

    metrics = [m_e_apsis, m_e_rebound, m_cross_e, m_cross_lz, m_a, m_e, m_i]
    all_passed = all(m.passed for m in metrics)

    # ── Informational (NOT gated): point-wise |Δr| on the TP ──────────── #
    max_dr = 0.0
    max_dr_year = 0
    for a, r in zip(apsis, rebound_samples):
        ax, ay, az = a.body_pos(TEST_PARTICLE_INDEX)
        rx, ry, rz = r.body_pos(TEST_PARTICLE_INDEX)
        dr = math.sqrt((ax - rx) ** 2 + (ay - ry) ** 2 + (az - rz) ** 2)
        if dr > max_dr:
            max_dr = dr
            max_dr_year = a.year
    info = MetricResult(
        name="|Δr| TP (NOT gated)",
        observed=max_dr,
        tolerance=float("inf"),
        passed=True,
        detail={
            "argmax_year": max_dr_year,
            "note": "phase-drift contaminated through encounter step; not invariant across IAS15 implementations",
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
            body_state = []
            for i in range(N_BODIES):
                body_state.extend([
                    float(row[f"x{i}"]),
                    float(row[f"y{i}"]),
                    float(row[f"z{i}"]),
                    float(row[f"vx{i}"]),
                    float(row[f"vy{i}"]),
                    float(row[f"vz{i}"]),
                ])
            samples.append(
                Sample(
                    year=int(row["year"]),
                    t=float(row["t"]),
                    body_state=body_state,
                    e_total=float(row["e_total"]),
                    lz_total=float(row["lz_total"]),
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
    print("REBOUND parity — Mercurius — comparison report")
    print(f"  samples : {n_samples}")
    print(f"  t_final : {t_final:.6e}  yr")
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
    argmax = info.detail.get("argmax_year", "n/a")
    print(f"    └── peak at year {argmax}; phase-drift contaminated, see protocol notebook")
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
        description="Comparator for REBOUND vs apsis parity — Mercurius.",
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
