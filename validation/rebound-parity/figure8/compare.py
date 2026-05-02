"""Comparator for REBOUND parity — Figure-8 choreography.

Loads `out/apsis.csv` and `out/rebound.csv`, computes the metrics defined
in the protocol notebook (see §Hypothesis), and exits with code 0 iff
every Tier-1 and Tier-2 gated metric is within its *a priori* tolerance.
On failure, exits non-zero with a structured report saved to
`out/comparison.json`.

## Metric hierarchy (declared *a priori* in the notebook)

The protocol organises metrics into three tiers reflecting evidentiary
weight. Verdict on the experiment is determined by Tier 1 and Tier 2
only. Tier 3 is reported for completeness, not as a gate.

### Tier 1 — Hard physical invariants (gated)
- `|ΔE/E_0|` per side          — energy conservation per implementation
- cross-impl `|ΔE|/|E_0|`      — energy agreement between implementations
- `|Δ𝐋|` per side              — angular-momentum conservation (vector norm; absolute)
- cross-impl `|Δ𝐋|`            — angular-momentum agreement (absolute)

### Tier 2 — Construction-level sanity (gated, weak)
- `|Δ𝐏|` per side              — linear-momentum conservation (absolute; 𝐏_0 ≈ 𝟎)
- cross-impl `|Δ𝐏|`            — linear-momentum agreement (absolute)
- `|Δ𝐫_COM|` per side          — COM drift from origin (absolute; 𝐫_COM(0) ≈ 𝟎)
- cross-impl `|Δ𝐫_COM|`        — COM agreement (absolute)

### Tier 3 — Geometric coherence (NOT gated)
- per-body `max |𝐫_apsis(t) − 𝐫_rebound(t)|`, max over bodies — phase-drift
  contaminated; reported as observational context only.

## Why the invariant set, not |Δr|

Adaptive high-order integrators (IAS15) are not bit-deterministic across
independent implementations: ULP-level differences in the controller's
`dt` selection accumulate as orbital *phase* drift over many periods.
Phase drift is not a numerical-correctness signal — both implementations
stay on the same dynamical trajectory, just advancing along it at
slightly different rates. Sampling at fixed times and computing
`|r_apsis(t) − r_rebound(t)|` therefore conflates phase drift (not
invariant across implementations) with geometric drift (the actual
physical signal). The protocol gates on the global integrals of motion
of the three-body system instead, which two correct implementations
must agree on regardless of phase. See the Kepler notebook
(`docs/experiments/2026-04-25-rebound-parity-kepler.md`) §Pilot
Interpretation for the full diagnostic narrative.

## Exit codes

- 0 — all Tier-1 and Tier-2 gated metrics within tolerance.
- 1 — input file error (missing file, sample count mismatch, t mismatch).
- 2 — at least one gated metric exceeded tolerance.

## Run

    python compare.py
    python compare.py --apsis-csv path/to/apsis.csv --rebound-csv path/to/rebound.csv

## Protocol notebook

    docs/experiments/2026-04-26-rebound-parity-figure8.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# ── Tolerances declared a priori (mirror the protocol notebook §Hypothesis) ─ #
#
# Tier 1 (hard physical invariants):
#   Energy and angular-momentum bounds at 50× f64 machine epsilon.
#   The angular-momentum reference scale is O(1) (per-body |m·𝐫×𝐯|);
#   absolute bound is appropriate because |𝐋_0| ≈ 0 by IC construction.
#
# Tier 2 (sanity):
#   Linear-momentum bound at 50× f64 ULP, absolute (|𝐏_0| ≈ 0).
#   COM-position bound one decade looser to account for position-scale
#   accumulation of ULP momentum noise over the 10T integration horizon.

# Tier 1
TOL_REL_ENERGY_PER_SIDE: float = 1.0e-13
TOL_REL_ENERGY_CROSS: float = 1.0e-13
TOL_ABS_ANGULAR_MOMENTUM_PER_SIDE: float = 1.0e-13
TOL_ABS_ANGULAR_MOMENTUM_CROSS: float = 1.0e-13

# Tier 2
TOL_ABS_LINEAR_MOMENTUM_PER_SIDE: float = 1.0e-13
TOL_ABS_LINEAR_MOMENTUM_CROSS: float = 1.0e-13
TOL_ABS_COM_PER_SIDE: float = 1.0e-12
TOL_ABS_COM_CROSS: float = 1.0e-12

# Number of bodies in the figure-8 system (fixed by protocol).
N_BODIES: int = 3
# Equal-mass IC; mirrored from the Rust example. Comparator uses this
# directly rather than parsing the CSV header to keep the metric formula
# transparent at the comparison site.
MASS: float = 1.0


# ── Data records ───────────────────────────────────────────────────────── #


@dataclass
class Sample:
    """One CSV row: per-body state at time `t`, plus total energy."""

    sample: int
    t: float
    bodies: list[tuple[float, float, float, float]]  # [(x, y, vx, vy)] × N_BODIES
    e_total: float


@dataclass
class Invariants:
    """Conserved quantities of the three-body system at one instant."""

    energy: float
    L: tuple[float, float, float]      # angular-momentum vector, full 3D
    P: tuple[float, float]             # linear momentum (planar; pz = 0)
    r_com: tuple[float, float]         # centre of mass (planar; z = 0)


@dataclass
class MetricResult:
    name: str
    tier: int
    observed: float
    tolerance: float
    passed: bool
    detail: dict = field(default_factory=dict)


# ── Invariant computation ──────────────────────────────────────────────── #


def physical_invariants(s: Sample, mass: float) -> Invariants:
    """Compute (E, 𝐋, 𝐏, 𝐫_COM) for one sample.

    Conventions match the apsis-side `total_energy` helper and REBOUND's
    `sim.energy()`: KE = ½ Σ m vᵢ², PE = −Σᵢ<ⱼ G mᵢ mⱼ / rᵢⱼ, with G = 1.

    Angular momentum is computed as the full 3-vector even though the
    system is planar by IC construction: this exposes any implementation-
    side breaking of planarity rather than baking the planar assumption
    into the metric. For the canonical figure-8 ICs, `L_x` and `L_y`
    must remain at the f64 round-off floor.
    """
    n = len(s.bodies)
    total_mass = mass * n

    # Energy
    ke = sum(0.5 * mass * (vx * vx + vy * vy) for (_, _, vx, vy) in s.bodies)
    pe = 0.0
    for i in range(n):
        xi, yi, _, _ = s.bodies[i]
        for j in range(i + 1, n):
            xj, yj, _, _ = s.bodies[j]
            dx = xi - xj
            dy = yi - yj
            r = math.sqrt(dx * dx + dy * dy)
            pe -= mass * mass / r
    energy = ke + pe

    # Angular momentum (full 3D vector; for planar ICs Lx = Ly = 0)
    Lx = 0.0
    Ly = 0.0
    Lz = 0.0
    for x, y, vx, vy in s.bodies:
        # 𝐫 × 𝐯 with 𝐫 = (x, y, 0), 𝐯 = (vx, vy, 0):
        #   (𝐫 × 𝐯)_x = y · 0 − 0 · vy = 0
        #   (𝐫 × 𝐯)_y = 0 · vx − x · 0 = 0
        #   (𝐫 × 𝐯)_z = x · vy − y · vx
        # Compute all three explicitly so any future migration to 3D
        # data flows through this formula unchanged.
        Lx += mass * (y * 0.0 - 0.0 * vy)
        Ly += mass * (0.0 * vx - x * 0.0)
        Lz += mass * (x * vy - y * vx)

    # Linear momentum
    Px = sum(mass * vx for (_, _, vx, _) in s.bodies)
    Py = sum(mass * vy for (_, _, _, vy) in s.bodies)

    # Centre-of-mass position
    com_x = sum(mass * x for (x, _, _, _) in s.bodies) / total_mass
    com_y = sum(mass * y for (_, y, _, _) in s.bodies) / total_mass

    return Invariants(
        energy=energy,
        L=(Lx, Ly, Lz),
        P=(Px, Py),
        r_com=(com_x, com_y),
    )


def vec_norm3(v: tuple[float, float, float]) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def vec_norm2(v: tuple[float, float]) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1])


def vec_diff_norm3(a, b) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    )


def vec_diff_norm2(a, b) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


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

    # Sanity check: t alignment. REBOUND should land on apsis's actual
    # sample times via exact_finish_time=1; significant divergence here
    # invalidates the comparison.
    max_dt = max(abs(a.t - r.t) for a, r in zip(apsis, rebound))
    if max_dt > 1.0e-12:
        print(
            f"ERROR: sample times disagree by up to {max_dt:.3e} — "
            "REBOUND's exact_finish_time may have failed to land on apsis's "
            "actual sample times.",
            file=sys.stderr,
        )
        return 1

    # ── Invariants per side per sample ─────────────────────────────────── #
    inv_apsis = [physical_invariants(s, MASS) for s in apsis]
    inv_rebound = [physical_invariants(s, MASS) for s in rebound]

    inv0_apsis = inv_apsis[0]
    inv0_rebound = inv_rebound[0]
    e0_apsis = inv0_apsis.energy
    e0_rebound = inv0_rebound.energy

    # ════════════════════════════════════════════════════════════════════ #
    # Tier 1 — Hard physical invariants (gated)
    # ════════════════════════════════════════════════════════════════════ #

    # 1a. |ΔE/E_0| per side
    max_de_apsis = max(abs(i.energy - e0_apsis) / abs(e0_apsis) for i in inv_apsis)
    max_de_rebound = max(
        abs(i.energy - e0_rebound) / abs(e0_rebound) for i in inv_rebound
    )
    m_e_apsis = MetricResult(
        name="|ΔE/E_0| apsis",
        tier=1,
        observed=max_de_apsis,
        tolerance=TOL_REL_ENERGY_PER_SIDE,
        passed=max_de_apsis <= TOL_REL_ENERGY_PER_SIDE,
    )
    m_e_rebound = MetricResult(
        name="|ΔE/E_0| rebound",
        tier=1,
        observed=max_de_rebound,
        tolerance=TOL_REL_ENERGY_PER_SIDE,
        passed=max_de_rebound <= TOL_REL_ENERGY_PER_SIDE,
    )

    # 1b. cross-impl |ΔE|/|E_0|
    max_de_cross = max(
        abs(a.energy - r.energy) / abs(e0_apsis)
        for a, r in zip(inv_apsis, inv_rebound)
    )
    m_e_cross = MetricResult(
        name="cross-impl |ΔE|/|E_0|",
        tier=1,
        observed=max_de_cross,
        tolerance=TOL_REL_ENERGY_CROSS,
        passed=max_de_cross <= TOL_REL_ENERGY_CROSS,
    )

    # 1c. |Δ𝐋| per side (vector norm of the drift from t=0)
    max_dL_apsis = max(vec_diff_norm3(i.L, inv0_apsis.L) for i in inv_apsis)
    max_dL_rebound = max(vec_diff_norm3(i.L, inv0_rebound.L) for i in inv_rebound)
    m_L_apsis = MetricResult(
        name="|Δ𝐋| apsis (abs)",
        tier=1,
        observed=max_dL_apsis,
        tolerance=TOL_ABS_ANGULAR_MOMENTUM_PER_SIDE,
        passed=max_dL_apsis <= TOL_ABS_ANGULAR_MOMENTUM_PER_SIDE,
    )
    m_L_rebound = MetricResult(
        name="|Δ𝐋| rebound (abs)",
        tier=1,
        observed=max_dL_rebound,
        tolerance=TOL_ABS_ANGULAR_MOMENTUM_PER_SIDE,
        passed=max_dL_rebound <= TOL_ABS_ANGULAR_MOMENTUM_PER_SIDE,
    )

    # 1d. cross-impl |Δ𝐋|
    max_dL_cross = max(
        vec_diff_norm3(a.L, r.L) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_L_cross = MetricResult(
        name="cross-impl |Δ𝐋| (abs)",
        tier=1,
        observed=max_dL_cross,
        tolerance=TOL_ABS_ANGULAR_MOMENTUM_CROSS,
        passed=max_dL_cross <= TOL_ABS_ANGULAR_MOMENTUM_CROSS,
    )

    # ════════════════════════════════════════════════════════════════════ #
    # Tier 2 — Construction-level sanity (gated, weak)
    # ════════════════════════════════════════════════════════════════════ #

    # 2a. |Δ𝐏| per side
    max_dP_apsis = max(vec_diff_norm2(i.P, inv0_apsis.P) for i in inv_apsis)
    max_dP_rebound = max(vec_diff_norm2(i.P, inv0_rebound.P) for i in inv_rebound)
    m_P_apsis = MetricResult(
        name="|Δ𝐏| apsis (abs)",
        tier=2,
        observed=max_dP_apsis,
        tolerance=TOL_ABS_LINEAR_MOMENTUM_PER_SIDE,
        passed=max_dP_apsis <= TOL_ABS_LINEAR_MOMENTUM_PER_SIDE,
    )
    m_P_rebound = MetricResult(
        name="|Δ𝐏| rebound (abs)",
        tier=2,
        observed=max_dP_rebound,
        tolerance=TOL_ABS_LINEAR_MOMENTUM_PER_SIDE,
        passed=max_dP_rebound <= TOL_ABS_LINEAR_MOMENTUM_PER_SIDE,
    )

    # 2b. cross-impl |Δ𝐏|
    max_dP_cross = max(
        vec_diff_norm2(a.P, r.P) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_P_cross = MetricResult(
        name="cross-impl |Δ𝐏| (abs)",
        tier=2,
        observed=max_dP_cross,
        tolerance=TOL_ABS_LINEAR_MOMENTUM_CROSS,
        passed=max_dP_cross <= TOL_ABS_LINEAR_MOMENTUM_CROSS,
    )

    # 2c. |Δ𝐫_COM| per side (drift from origin, since r_COM(0) ≈ 0)
    max_com_apsis = max(vec_diff_norm2(i.r_com, inv0_apsis.r_com) for i in inv_apsis)
    max_com_rebound = max(
        vec_diff_norm2(i.r_com, inv0_rebound.r_com) for i in inv_rebound
    )
    m_com_apsis = MetricResult(
        name="|Δ𝐫_COM| apsis (abs)",
        tier=2,
        observed=max_com_apsis,
        tolerance=TOL_ABS_COM_PER_SIDE,
        passed=max_com_apsis <= TOL_ABS_COM_PER_SIDE,
    )
    m_com_rebound = MetricResult(
        name="|Δ𝐫_COM| rebound (abs)",
        tier=2,
        observed=max_com_rebound,
        tolerance=TOL_ABS_COM_PER_SIDE,
        passed=max_com_rebound <= TOL_ABS_COM_PER_SIDE,
    )

    # 2d. cross-impl |Δ𝐫_COM|
    max_com_cross = max(
        vec_diff_norm2(a.r_com, r.r_com) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_com_cross = MetricResult(
        name="cross-impl |Δ𝐫_COM| (abs)",
        tier=2,
        observed=max_com_cross,
        tolerance=TOL_ABS_COM_CROSS,
        passed=max_com_cross <= TOL_ABS_COM_CROSS,
    )

    tier1 = [m_e_apsis, m_e_rebound, m_e_cross, m_L_apsis, m_L_rebound, m_L_cross]
    tier2 = [m_P_apsis, m_P_rebound, m_P_cross, m_com_apsis, m_com_rebound, m_com_cross]
    gated = tier1 + tier2
    all_passed = all(m.passed for m in gated)

    # ════════════════════════════════════════════════════════════════════ #
    # Tier 3 — Geometric coherence (informational, NOT gated)
    # ════════════════════════════════════════════════════════════════════ #

    max_dr_per_body: list[float] = [0.0] * N_BODIES
    argmax_sample_per_body: list[int] = [0] * N_BODIES
    for a, r in zip(apsis, rebound):
        for k in range(N_BODIES):
            ax, ay, _, _ = a.bodies[k]
            rx, ry, _, _ = r.bodies[k]
            dr = math.sqrt((ax - rx) ** 2 + (ay - ry) ** 2)
            if dr > max_dr_per_body[k]:
                max_dr_per_body[k] = dr
                argmax_sample_per_body[k] = a.sample
    max_dr_overall = max(max_dr_per_body)
    info_dr = MetricResult(
        name="|Δ𝐫| max over bodies (NOT gated)",
        tier=3,
        observed=max_dr_overall,
        tolerance=float("inf"),
        passed=True,
        detail={
            "per_body_max_dr": max_dr_per_body,
            "argmax_sample_per_body": argmax_sample_per_body,
            "note": "phase-drift contaminated; not invariant across adaptive integrators",
        },
    )

    # ── Report ──────────────────────────────────────────────────────────── #
    print_report(
        tier1, tier2, info_dr, len(apsis), apsis[-1].t, e0_apsis, e0_rebound,
    )
    write_json_report(
        output_dir, tier1, tier2, info_dr, len(apsis), apsis[-1].t, all_passed,
    )

    return 0 if all_passed else 2


# ── I/O ─────────────────────────────────────────────────────────────────── #


def load_csv(path: Path) -> list[Sample]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found at {path}")
    samples: list[Sample] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            bodies = [
                (
                    float(row[f"x{k}"]),
                    float(row[f"y{k}"]),
                    float(row[f"vx{k}"]),
                    float(row[f"vy{k}"]),
                )
                for k in range(N_BODIES)
            ]
            samples.append(
                Sample(
                    sample=int(row["sample"]),
                    t=float(row["t"]),
                    bodies=bodies,
                    e_total=float(row["e_total"]),
                )
            )
    return samples


def print_report(
    tier1: list[MetricResult],
    tier2: list[MetricResult],
    info: MetricResult,
    n_samples: int,
    t_final: float,
    e0_apsis: float,
    e0_rebound: float,
) -> None:
    print()
    print("REBOUND parity — Figure-8 — comparison report (tiered metrics)")
    print(f"  samples : {n_samples}")
    print(f"  t_final : {t_final:.6e}")
    print(f"  E_0 apsis   : {e0_apsis:+.18e}")
    print(f"  E_0 rebound : {e0_rebound:+.18e}")

    def _print_tier(label: str, metrics: list[MetricResult]) -> None:
        print()
        print(f"  ── {label} ──")
        print(f"  {'metric':<32} {'observed':>14} {'tolerance':>14}  verdict")
        print(f"  {'-' * 32} {'-' * 14} {'-' * 14}  {'-' * 7}")
        for m in metrics:
            verdict = "pass" if m.passed else "FAIL"
            print(
                f"  {m.name:<32} {m.observed:>14.3e} {m.tolerance:>14.3e}  {verdict}"
            )

    _print_tier("Tier 1 — hard physical invariants (gated)", tier1)
    _print_tier("Tier 2 — construction-level sanity (gated)", tier2)

    print()
    print("  ── Tier 3 — geometric coherence (NOT gated) ──")
    print(f"  {'metric':<32} {'observed':>14}")
    print(f"  {'-' * 32} {'-' * 14}")
    print(f"  {info.name:<32} {info.observed:>14.3e}")
    per_body = info.detail.get("per_body_max_dr", [])
    if per_body:
        per_body_fmt = ", ".join(f"{v:.3e}" for v in per_body)
        print(f"    └── per body: [{per_body_fmt}]")
    print(
        "    └── phase-drift contaminated; see protocol notebook §"
        "Why this metric set, not |Δr|"
    )
    print()


def write_json_report(
    output_dir: Path,
    tier1: list[MetricResult],
    tier2: list[MetricResult],
    info: MetricResult,
    n_samples: int,
    t_final: float,
    all_passed: bool,
) -> None:
    report = {
        "all_passed": all_passed,
        "n_samples": n_samples,
        "t_final": t_final,
        "tier1_hard": [asdict(m) for m in tier1],
        "tier2_sanity": [asdict(m) for m in tier2],
        "tier3_informational": asdict(info),
    }
    (output_dir / "comparison.json").write_text(json.dumps(report, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Comparator for REBOUND vs apsis parity — Figure-8.",
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
