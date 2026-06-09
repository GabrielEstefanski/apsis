"""Comparator for REBOUND parity — Pythagorean three-body (Burrau 1913).

Loads `out/apsis.csv` and `out/rebound.csv`, computes the metrics defined
in the protocol notebook (see §Hypothesis), and exits with code 0 iff
every Tier-1 and Tier-2 gated metric is within its *a priori* tolerance.
On failure, exits non-zero with a structured report saved to
`out/comparison.json`.

Mirror of `validation/rebound-parity/figure8/compare.py` with one
generalisation: `physical_invariants` accepts a `masses` tuple rather
than a single equal-mass scalar. The metric structure, tolerances, JSON
schema, and stdout report layout are byte-compatible with the figure-8
comparator so a downstream tool can consume any parity scenario's
output uniformly.

## Metric hierarchy (declared *a priori* in the notebook)

The protocol organises metrics into three tiers reflecting evidentiary
weight. Verdict on the experiment is determined by Tier 1 and Tier 2
only. Tier 3 is reported per-sample as observational context — never
aggregated into a pass/fail criterion.

### Tier 1 — Hard physical invariants (gated)
- `|ΔE/E_0|` per side          — energy conservation per implementation
- cross-impl `|ΔE|/|E_0|`      — energy agreement between implementations
- `|Δ𝐋|` per side              — angular-momentum conservation (vector norm; absolute)
- cross-impl `|Δ𝐋|`            — angular-momentum agreement (absolute)

### Tier 2 — Construction-level sanity (gated, weak)
- `|Δ𝐏|` per side              — linear-momentum conservation (absolute; 𝐏_0 = 𝟎 by IC)
- cross-impl `|Δ𝐏|`            — linear-momentum agreement (absolute)
- `|Δ𝐫_COM|` per side          — COM drift from origin (absolute; 𝐫_COM(0) = 𝟎 by IC)
- cross-impl `|Δ𝐫_COM|`        — COM agreement (absolute)

### Tier 3 — Geometric coherence (NOT gated)
- per-body `max |𝐫_apsis(t) − 𝐫_rebound(t)|`, max over bodies — phase-drift
  contaminated and Lyapunov-amplified for the Pythagorean dynamics;
  reported per-sample as observational context only. **Never aggregated
  into a pass/fail criterion.** Per-body magnitudes are expected to
  reach O(1) before the horizon.

## Why the invariant set, not |Δr|

For the Pythagorean three-body problem the trajectory is exponentially
sensitive to ULP-level controller decisions: two correct adaptive IAS15
implementations diverge in trajectory and agree in conserved quantities.
A protocol gated on `|Δr|` would be guaranteed to report failure
regardless of which implementation is correct. The protocol gates on
the global integrals of motion, which two correct implementations must
agree on regardless of the substep schedule that traces the flow. See
the Kepler notebook (`paper/notebooks/2026-04-25-rebound-parity-kepler.md`)
§Pilot Interpretation for the full diagnostic narrative on adaptive
high-order integrators.

## Exit codes

- 0 — all Tier-1 and Tier-2 gated metrics within tolerance.
- 1 — input file error (missing file, sample count mismatch, t mismatch).
- 2 — at least one gated metric exceeded tolerance.

## Run

    python compare.py
    python compare.py --apsis-csv path/to/apsis.csv --rebound-csv path/to/rebound.csv

## Protocol notebook

    paper/notebooks/2026-04-30-rebound-parity-pythagorean.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

EPS: float = 2.220446049250313e-16

# ── Tolerances ─────────────────────────────────────────────────────────── #
# Energy: the chaotic regime admits no tight a-priori floor — the deepest r_min
# is chaotic and shifts under ULP/platform differences. Sanity ceiling at the
# published f64 energy floor for the Pythagorean problem (Boekholt & Portegies
# Zwart 2015 §3.1, ~1e-8): ~1.4 decades above our IAS15 observed (~4e-10), ~2
# below the ~1e-6 gross-error level conventional in collisional N-body work.
# Deliberately the loose gate — parity is borne by L/P/COM (below) + the
# close-encounter alignment; energy conservation ≠ converged trajectory
# (Boekholt & PZ 2015 §3.3). See the protocol notebook.
TOL_REL_ENERGY_PER_SIDE: float = 1.0e-8
TOL_REL_ENERGY_CROSS: float = 1.0e-8
# L, P, COM = 0 by IC → Wilkinson cancellation floors, derived per-run in main().

# Number of bodies in the Pythagorean system (fixed by protocol).
N_BODIES: int = 3
# Burrau (1913) opposite-side mass convention. The comparator uses this
# directly rather than parsing the CSV header to keep the metric formula
# transparent at the comparison site.
MASSES: tuple[float, float, float] = (3.0, 4.0, 5.0)


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


def physical_invariants(s: Sample, masses: tuple[float, ...]) -> Invariants:
    """Compute (E, 𝐋, 𝐏, 𝐫_COM) for one sample.

    Conventions match the apsis-side `total_energy` helper and REBOUND's
    `sim.energy()`: KE = ½ Σ mᵢ vᵢ², PE = −Σᵢ<ⱼ G mᵢ mⱼ / rᵢⱼ, with G = 1.

    `masses` is per-body and may be heterogeneous (here: 3, 4, 5). The
    figure-8 comparator's equal-mass scalar form is the special case
    `masses = (M, M, ..., M)`; the formula here generalises it without
    altering the metric semantics.
    """
    assert len(masses) == len(s.bodies), "masses and bodies length must match"
    n = len(s.bodies)
    total_mass = sum(masses)

    # Energy
    ke = sum(
        0.5 * masses[i] * (vx * vx + vy * vy)
        for i, (_, _, vx, vy) in enumerate(s.bodies)
    )
    pe = 0.0
    for i in range(n):
        xi, yi, _, _ = s.bodies[i]
        for j in range(i + 1, n):
            xj, yj, _, _ = s.bodies[j]
            dx = xi - xj
            dy = yi - yj
            r = math.sqrt(dx * dx + dy * dy)
            pe -= masses[i] * masses[j] / r
    energy = ke + pe

    # Angular momentum: planar data, so only Lz is non-zero.
    Lx = Ly = 0.0
    Lz = sum(masses[i] * (x * vy - y * vx) for i, (x, y, vx, vy) in enumerate(s.bodies))

    # Linear momentum
    Px = sum(masses[i] * vx for i, (_, _, vx, _) in enumerate(s.bodies))
    Py = sum(masses[i] * vy for i, (_, _, _, vy) in enumerate(s.bodies))

    # Centre-of-mass position
    com_x = sum(masses[i] * x for i, (x, _, _, _) in enumerate(s.bodies)) / total_mass
    com_y = sum(masses[i] * y for i, (_, y, _, _) in enumerate(s.bodies)) / total_mass

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


def scale_P(bodies: list[tuple[float, float, float, float]], masses: tuple[float, ...]) -> float:
    """Wilkinson cancellation scale of 𝐏 = Σ mᵢ𝐯ᵢ (0 by IC): Σ|mᵢ𝐯ᵢ| per axis."""
    return math.hypot(
        sum(abs(masses[i] * vx) for i, (_, _, vx, _) in enumerate(bodies)),
        sum(abs(masses[i] * vy) for i, (_, _, _, vy) in enumerate(bodies)),
    )


def scale_L(bodies: list[tuple[float, float, float, float]], masses: tuple[float, ...]) -> float:
    """Wilkinson cancellation scale of L_z = Σ mᵢ(xᵢvyᵢ − yᵢvxᵢ) (0 by IC)."""
    return sum(
        abs(masses[i] * x * vy) + abs(masses[i] * y * vx)
        for i, (x, y, vx, vy) in enumerate(bodies)
    )


# ── Close-encounter detection (debug telemetry, never gated) ───────────── #


def pair_min_distance(s: Sample) -> float:
    """Smallest pairwise separation among the three bodies at sample `s`.

    Captures the proxy used by the IAS15 controller to drive `dt` shrinkage:
    when r_min collapses, the controller responds by tightening dt. Logging
    r_min(t) per side lets us check whether the two implementations are
    *responding to the same dynamical events* — the strongest available
    signal that the chaotic dynamics are being integrated equivalently.
    """
    n = len(s.bodies)
    r_min = float("inf")
    for i in range(n):
        xi, yi, _, _ = s.bodies[i]
        for j in range(i + 1, n):
            xj, yj, _, _ = s.bodies[j]
            dx = xi - xj
            dy = yi - yj
            r = math.sqrt(dx * dx + dy * dy)
            if r < r_min:
                r_min = r
    return r_min


@dataclass
class CloseEncounter:
    """One local minimum of r_min(t) — a close-encounter event."""

    t: float
    r_min: float
    sample: int


def detect_close_encounters(
    samples: list[Sample],
    half_window: int = 3,
    prominence_factor: float = 0.5,
) -> list[CloseEncounter]:
    """Detect prominent local minima of r_min(t).

    A sample `i` is an event if it is the minimum within a `[i - half_window,
    i + half_window]` window AND its `r_min` is below `prominence_factor`
    times the median `r_min` over the full series. The prominence filter
    rejects shallow background ripples; the window enforces locality.

    `half_window = 3` corresponds to ~0.1 t.u. at the 30-samples/t.u.
    cadence — small enough to resolve distinct close-approach events,
    large enough to suppress per-sample noise. Both parameters are
    debug-grade — they affect the diagnostic output, not any gated metric.
    """
    if len(samples) < 2 * half_window + 1:
        return []
    r_series = [pair_min_distance(s) for s in samples]
    sorted_r = sorted(r_series)
    median_r = sorted_r[len(sorted_r) // 2]
    threshold = prominence_factor * median_r

    events: list[CloseEncounter] = []
    for i in range(half_window, len(samples) - half_window):
        r_i = r_series[i]
        if r_i >= threshold:
            continue
        is_local_min = all(
            r_i <= r_series[j]
            for j in range(i - half_window, i + half_window + 1)
            if j != i
        )
        if is_local_min:
            events.append(CloseEncounter(t=samples[i].t, r_min=r_i, sample=samples[i].sample))
    return events


def match_encounters(
    apsis_events: list[CloseEncounter],
    rebound_events: list[CloseEncounter],
) -> tuple[list[tuple[CloseEncounter, CloseEncounter]], list[CloseEncounter], list[CloseEncounter]]:
    """Match apsis events to REBOUND events using nearest-neighbour-in-time
    pairing within a windowed tolerance.

    The window adapts to the typical inter-event interval on the apsis
    side: window = 0.5 × median(t_{i+1} - t_i). A pair (a, r) is matched
    iff |t_a - t_r| < window AND |log10(a.r_min / r.r_min)| < 0.5
    (within half a decade of the same magnitude).

    Returns `(matched, unmatched_apsis, unmatched_rebound)`. Used for
    diagnostic reporting; never enters a pass/fail criterion.
    """
    if len(apsis_events) < 2 or len(rebound_events) < 2:
        # Too few events for an inter-event timescale; fall back to
        # 1.0 t.u. as a conservative fixed window.
        window = 1.0
    else:
        intervals = sorted(
            apsis_events[i + 1].t - apsis_events[i].t
            for i in range(len(apsis_events) - 1)
        )
        median_interval = intervals[len(intervals) // 2]
        window = 0.5 * median_interval

    matched: list[tuple[CloseEncounter, CloseEncounter]] = []
    used_rebound: set[int] = set()

    for a_event in apsis_events:
        best_idx: int | None = None
        best_dt = float("inf")
        for j, r_event in enumerate(rebound_events):
            if j in used_rebound:
                continue
            dt = abs(a_event.t - r_event.t)
            if dt < best_dt:
                best_dt = dt
                best_idx = j
        if best_idx is not None and best_dt < window:
            r_event = rebound_events[best_idx]
            mag_ratio = abs(math.log10(a_event.r_min / r_event.r_min))
            if mag_ratio < 0.5:
                matched.append((a_event, r_event))
                used_rebound.add(best_idx)

    unmatched_apsis = [
        a for a in apsis_events
        if all(a is not pair[0] for pair in matched)
    ]
    unmatched_rebound = [
        r for j, r in enumerate(rebound_events)
        if j not in used_rebound
    ]
    return matched, unmatched_apsis, unmatched_rebound


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
    inv_apsis = [physical_invariants(s, MASSES) for s in apsis]
    inv_rebound = [physical_invariants(s, MASSES) for s in rebound]

    inv0_apsis = inv_apsis[0]
    inv0_rebound = inv_rebound[0]
    e0_apsis = inv0_apsis.energy
    e0_rebound = inv0_rebound.energy

    # 𝐏, 𝐋 = 0 by IC → Wilkinson floor 10·EPS·max_t Σ|terms|, each side on its
    # own scale; cross √2× the larger.
    p_scale_a = max(scale_P(s.bodies, MASSES) for s in apsis)
    p_scale_r = max(scale_P(s.bodies, MASSES) for s in rebound)
    l_scale_a = max(scale_L(s.bodies, MASSES) for s in apsis)
    l_scale_r = max(scale_L(s.bodies, MASSES) for s in rebound)
    tol_p_apsis = 10.0 * EPS * p_scale_a
    tol_p_rebound = 10.0 * EPS * p_scale_r
    tol_p_cross = 10.0 * math.sqrt(2.0) * EPS * max(p_scale_a, p_scale_r)
    tol_l_apsis = 10.0 * EPS * l_scale_a
    tol_l_rebound = 10.0 * EPS * l_scale_r
    tol_l_cross = 10.0 * math.sqrt(2.0) * EPS * max(l_scale_a, l_scale_r)

    # 𝐫_COM = 0 by IC → 1.5·EPS·P_scale·t_final/M (drift term; representation
    # subdominant), each side on its own scale; cross 2×.
    t_final = apsis[-1].t
    tol_com_apsis = 1.5 * EPS * p_scale_a * t_final / sum(MASSES)
    tol_com_rebound = 1.5 * EPS * p_scale_r * t_final / sum(MASSES)
    tol_com_cross = 1.5 * 2.0 * EPS * max(p_scale_a, p_scale_r) * t_final / sum(MASSES)

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
        tolerance=tol_l_apsis,
        passed=max_dL_apsis <= tol_l_apsis,
    )
    m_L_rebound = MetricResult(
        name="|Δ𝐋| rebound (abs)",
        tier=1,
        observed=max_dL_rebound,
        tolerance=tol_l_rebound,
        passed=max_dL_rebound <= tol_l_rebound,
    )

    # 1d. cross-impl |Δ𝐋|
    max_dL_cross = max(
        vec_diff_norm3(a.L, r.L) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_L_cross = MetricResult(
        name="cross-impl |Δ𝐋| (abs)",
        tier=1,
        observed=max_dL_cross,
        tolerance=tol_l_cross,
        passed=max_dL_cross <= tol_l_cross,
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
        tolerance=tol_p_apsis,
        passed=max_dP_apsis <= tol_p_apsis,
    )
    m_P_rebound = MetricResult(
        name="|Δ𝐏| rebound (abs)",
        tier=2,
        observed=max_dP_rebound,
        tolerance=tol_p_rebound,
        passed=max_dP_rebound <= tol_p_rebound,
    )

    # 2b. cross-impl |Δ𝐏|
    max_dP_cross = max(
        vec_diff_norm2(a.P, r.P) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_P_cross = MetricResult(
        name="cross-impl |Δ𝐏| (abs)",
        tier=2,
        observed=max_dP_cross,
        tolerance=tol_p_cross,
        passed=max_dP_cross <= tol_p_cross,
    )

    # 2c. |Δ𝐫_COM| per side (drift from origin, since r_COM(0) = 0 by IC)
    max_com_apsis = max(vec_diff_norm2(i.r_com, inv0_apsis.r_com) for i in inv_apsis)
    max_com_rebound = max(
        vec_diff_norm2(i.r_com, inv0_rebound.r_com) for i in inv_rebound
    )
    m_com_apsis = MetricResult(
        name="|Δ𝐫_COM| apsis (abs)",
        tier=2,
        observed=max_com_apsis,
        tolerance=tol_com_apsis,
        passed=max_com_apsis <= tol_com_apsis,
    )
    m_com_rebound = MetricResult(
        name="|Δ𝐫_COM| rebound (abs)",
        tier=2,
        observed=max_com_rebound,
        tolerance=tol_com_rebound,
        passed=max_com_rebound <= tol_com_rebound,
    )

    # 2d. cross-impl |Δ𝐫_COM|
    max_com_cross = max(
        vec_diff_norm2(a.r_com, r.r_com) for a, r in zip(inv_apsis, inv_rebound)
    )
    m_com_cross = MetricResult(
        name="cross-impl |Δ𝐫_COM| (abs)",
        tier=2,
        observed=max_com_cross,
        tolerance=tol_com_cross,
        passed=max_com_cross <= tol_com_cross,
    )

    tier1 = [m_e_apsis, m_e_rebound, m_e_cross, m_L_apsis, m_L_rebound, m_L_cross]
    tier2 = [m_P_apsis, m_P_rebound, m_P_cross, m_com_apsis, m_com_rebound, m_com_cross]
    gated = tier1 + tier2
    all_passed = all(m.passed for m in gated)

    # ════════════════════════════════════════════════════════════════════ #
    # Tier 3 — Geometric coherence (informational, NOT gated)
    # ════════════════════════════════════════════════════════════════════ #
    #
    # Per-sample, never aggregated into a pass/fail criterion. For the
    # Pythagorean dynamics, per-body |Δr| is expected to reach O(1) before
    # the horizon as Lyapunov instability amplifies ULP-level controller
    # decisions into trajectory-scale separation. The number is reported
    # for shape; it has no tolerance and never affects the verdict.

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
            "note": "phase-drift contaminated and Lyapunov-amplified; not a gate, never aggregated into pass/fail",
        },
    )

    # ── Close-encounter detection (debug telemetry) ─────────────────────── #
    #
    # Per-side prominent local minima of r_min(t), then nearest-in-time
    # pairing within a windowed tolerance. Reports match fraction and any
    # unmatched events. Never enters the verdict — see the protocol notebook
    # §Verdict criterion (Tier 1 + Tier 2 only gate).
    apsis_events = detect_close_encounters(apsis)
    rebound_events = detect_close_encounters(rebound)
    matched, unmatched_apsis, unmatched_rebound = match_encounters(
        apsis_events, rebound_events,
    )
    encounters_summary = {
        "n_apsis": len(apsis_events),
        "n_rebound": len(rebound_events),
        "n_matched": len(matched),
        "match_fraction_apsis": (
            len(matched) / len(apsis_events) if apsis_events else 0.0
        ),
        "match_fraction_rebound": (
            len(matched) / len(rebound_events) if rebound_events else 0.0
        ),
        "matched": [
            {
                "t_apsis": a.t, "r_apsis": a.r_min, "sample_apsis": a.sample,
                "t_rebound": r.t, "r_rebound": r.r_min, "sample_rebound": r.sample,
                "dt": a.t - r.t,
                "log10_r_ratio": math.log10(a.r_min / r.r_min),
            }
            for a, r in matched
        ],
        "unmatched_apsis": [
            {"t": e.t, "r_min": e.r_min, "sample": e.sample}
            for e in unmatched_apsis
        ],
        "unmatched_rebound": [
            {"t": e.t, "r_min": e.r_min, "sample": e.sample}
            for e in unmatched_rebound
        ],
    }

    # ── Report ──────────────────────────────────────────────────────────── #
    print_report(
        tier1, tier2, info_dr, encounters_summary,
        len(apsis), apsis[-1].t, e0_apsis, e0_rebound,
    )
    write_json_report(
        output_dir, tier1, tier2, info_dr, encounters_summary,
        len(apsis), apsis[-1].t, all_passed,
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
    encounters: dict,
    n_samples: int,
    t_final: float,
    e0_apsis: float,
    e0_rebound: float,
) -> None:
    print()
    print("REBOUND parity — Pythagorean — comparison report (tiered metrics)")
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
    print("  ── Tier 3 — geometric coherence (NOT gated, per-sample) ──")
    print(f"  {'metric':<32} {'observed':>14}")
    print(f"  {'-' * 32} {'-' * 14}")
    print(f"  {info.name:<32} {info.observed:>14.3e}")
    per_body = info.detail.get("per_body_max_dr", [])
    if per_body:
        per_body_fmt = ", ".join(f"{v:.3e}" for v in per_body)
        print(f"    └── per body: [{per_body_fmt}]")
    print(
        "    └── phase-drift contaminated and Lyapunov-amplified; "
        "see protocol notebook §Why this metric set"
    )

    print()
    print("  ── Close-encounter alignment (debug telemetry, NEVER gated) ──")
    print(
        f"  events apsis  : {encounters['n_apsis']} prominent local minima of r_min(t)"
    )
    print(
        f"  events rebound: {encounters['n_rebound']} prominent local minima of r_min(t)"
    )
    print(
        f"  matched pairs : {encounters['n_matched']} "
        f"(apsis fraction {encounters['match_fraction_apsis']:.0%}, "
        f"rebound fraction {encounters['match_fraction_rebound']:.0%})"
    )
    if encounters["matched"]:
        max_dt = max(abs(p["dt"]) for p in encounters["matched"])
        max_log_ratio = max(
            abs(p["log10_r_ratio"]) for p in encounters["matched"]
        )
        print(
            f"  worst |Δt| in matched set: {max_dt:.3e} t.u.   "
            f"worst |log10(r_apsis / r_rebound)|: {max_log_ratio:.3f}"
        )
    n_unmatched = (
        len(encounters["unmatched_apsis"]) + len(encounters["unmatched_rebound"])
    )
    if n_unmatched:
        print(
            f"  unmatched: {len(encounters['unmatched_apsis'])} apsis + "
            f"{len(encounters['unmatched_rebound'])} rebound (chaos-driven phase drift)"
        )
    print()


def write_json_report(
    output_dir: Path,
    tier1: list[MetricResult],
    tier2: list[MetricResult],
    info: MetricResult,
    encounters: dict,
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
        "close_encounters_debug": encounters,
    }
    (output_dir / "comparison.json").write_text(json.dumps(report, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Comparator for REBOUND vs apsis parity — Pythagorean three-body.",
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
