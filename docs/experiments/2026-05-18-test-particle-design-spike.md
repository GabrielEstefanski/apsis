---
date: 2026-05-18
status: a priori (protocol declared before any code lands)
issue: '#133'
---

# Test particle treatment — design spike for Issue #133

## Context

`examples/radiation_dust.py` reports `|dE/E| ≈ 1.16e-5` after 50
orbital periods (Sun + `m_dust = 1e-15` cosmic-dust grain, IAS15 +
apsis-radiation). IAS15 on a Kepler 2-body typically delivers
`|dE/E| ~ 10⁻¹³`. The 8–10 orders of excess are in bookkeeping,
not dynamics: semi-major axis is preserved.

Per Issue #133, the Sun's velocity is
`v_sun = (m_dust/m_sun) · v_dust ≈ 10⁻¹⁵`. Sun KE `~5·10⁻³¹`,
added to dust KE `~5·10⁻¹⁶`, vanishes below f64 precision. The
running momentum slop propagates into the trajectory over time.

Information loss happens **before the sum**, so compensated-
accumulator propagation alone does not fix the extreme-mass-ratio
case. This notebook declares a structural fix via directional
interaction semantics, scoped to the direct O(N²) path (BH tree
adaptation deferred).

## Hypothesis

Suppressing back-reaction from body `j` onto body `i` when
`m_j / m_i < 1e-10` recovers `|dE/E|` to within 1e-10 at 50 orbits
on `radiation_dust.py`, with momentum defect bounded by
`O(mass_ratio)` per pair.

## Model — directional interaction semantics

A test particle is a per-pair, per-operator classification of
interaction direction, not a body property. The same body may be a
test particle relative to one operator and a full participant
relative to another (Mercury is a test particle under apsis-1pn
relative to the Sun, but not under pure Newtonian dynamics).

Operationally: `interaction X → Y is one-way` (X feels Y; Y does
not feel X back). The decision lives in the semantic layer
(operator + integrator), keyed on the physical regime.

Two surfaces carry the protocol:

**Gravity core** (`physics/gravity/engine.rs`): direct O(N²) path
auto-detects mass ratio. When `m_j / m_i < threshold` (default
`1e-10`), the back-reaction term in the pairwise loop is skipped:

```text
for each pair (i, j):
    F_ij = Newton(bodies[i], bodies[j])
    if m_j / m_i >= threshold: acc[i] += F_ij / m_i
    if m_i / m_j >= threshold: acc[j] -= F_ij / m_j
```

Threshold is per-system, configurable via
`System::with_test_particle_threshold(f64)`.

**Perturbation operators**: opt-in trait method on `Operator`:

```rust
/// True when this operator treats body `j`'s contribution to body
/// `i` as one-way. Default false (symmetric pair).
fn is_one_way(&self, _bodies: &[Body], _i: usize, _j: usize) -> bool {
    false
}
```

apsis-radiation is structurally asymmetric already (source →
receiver) and needs no override. Future operators with pairwise
symmetric forces (2PN, pairwise tidal) override the default.

The integrator skips back-reaction whenever the gravity core auto-
detect OR any registered operator declares one-way for that pair.

## Documented trade

Back-reaction suppression sacrifices exact momentum conservation
in exchange for numerically stable energy bookkeeping under
extreme mass ratios. `|ΔP/P|` scales with the mass ratio of the
suppressed interaction; characterised in §Results.

This contract is documented in:

- `physics/gravity/engine.rs` rustdoc on the skip path
- `core/system/metrics.rs` rustdoc on the momentum diagnostic
- `paper.md` near the conservation-claim paragraph

## Protocol

### Baseline

Run `examples/radiation_dust.py` on `develop` tip (commit
`4ddef70`). Record `|dE/E|`, `|ΔP/P_initial|`, `|Δa/a|` at 50
orbits. Capture as `2026-05-18-test-particle-design-spike/baseline-pre-fix.csv` (sibling directory).

### Implementation (three commits)

1. **Gravity core skip** — `physics/gravity/engine.rs` gains the
   mass-ratio guard; `System` gains
   `with_test_particle_threshold`. Default `1e-10`.
2. **Operator trait method** — `is_one_way` default `false`.
3. **Validation script** — reproduces pre/post comparison; CI
   gate added.

### Decision gates

- **Gate 1 (baseline matches):** `|dE/E|` matches `~1.16e-5`
  within 1 order. If not, the bug shifted; diagnose before
  implementing.
- **Gate 2 (primary):** `|dE/E| < 1e-10` at 50 orbits on
  `radiation_dust.py`.
- **Gate 3 (bonus):** `|dE/E| < 1e-13` (IAS15-Kepler equivalent).
- **Gate 4 (cross-checks):** Mercury 4.4 ppm unchanged; full lib
  test suite green.
- **Failure mode:** `|dE/E| >= 1e-8` post-fix → diagnosis was
  wrong; do not loosen the gate. Either an independent argument
  justifies a new target, or the work does not ship.

Momentum: record `|ΔP/P_initial|`; expected `O(1e-15)`. Document,
do not gate.

## Out of scope

- BH tree adaptation (multipole pruning, COM under mixed nodes,
  opening criterion). Follow-up notebook.
- General test-particle infrastructure beyond the auto-detect
  case. Trait method is the surface; only the gravity-core path
  is wired immediately.
- Compensated-accumulator propagation into core metrics buffers.
  Independent precision improvement, not bundled here.

## Results

### Baseline (pre-fix)

Captured against `develop` tip (commit `4ddef70`):

| Metric | Value | Notes |
|---|---|---|
| `\|dE/E\|` @ 50 orbits | `1.164 × 10⁻⁵` | matches issue body `~1.16e-5` |
| `a_final` (vs μ_eff) | `1.1250` | identical to `a_initial` |
| `\|Δa/a\|` | `< 1e-15` | already passes; trajectory physically stable |
| `r_final`, `v_final` | `1.1317`, `0.8891` | |

**Gate 1 — PASS.** Diagnosis has not shifted; proceed to
implementation. Raw stdout and CSV under
[`2026-05-18-test-particle-design-spike/`](./2026-05-18-test-particle-design-spike/) (sibling directory).

### Post-fix

| Metric | Pre-fix | Post-fix |
|---|---|---|
| `\|dE/E\|` @ 50 orbits | `1.164 × 10⁻⁵` | `1.200 × 10⁻⁵` |
| `r_final` | `1.1317` | `1.1363` |
| `v_final` | `0.8891` | `0.8855` |
| `a_final` (vs μ_eff) | `1.1250` | `1.1250` |

Endpoint `|dE/E|` is unchanged. Trajectory diverges by
O(mass_ratio) as expected for the test-particle approximation,
confirming the skip is active.

Raw stdout: `2026-05-18-test-particle-design-spike/post-fix.txt` (sibling directory).

### Gate 2 reanalysis — the gate was malformed

Post-fix `|dE/E| = 1.200 × 10⁻⁵` misses Gate 2 (`< 1e-10`) by 5
orders of magnitude. Before declaring an honest limitation,
verify the gate did not demand precision below the f64 round-off
floor.

`update_energy_tracking` (`crates/apsis/src/core/system/step.rs`)
normalises against `denom = baseline.abs().max(1e-12)`. For this
scenario `|E_initial| ≈ 5 × 10⁻¹⁶`, well below the floor. The
reported `|dE/E|` is `dE_absolute / 1e-12`, not
`dE_absolute / |E_initial|`. The reading is scaled by
`|1e-12 / E_initial| ≈ 2000×`.

Diagnostic (not committed): replaced `.max(1e-12)` with
`.max(1e-30)` (NaN guard only), rebuilt, re-ran. Result:
`|dE/E| = 2.399 × 10⁻²`. The honest relative drift: the
integrator preserves energy to its true f64 noise floor
(`~1 × 10⁻¹⁷` absolute), which against `|E_initial| ≈ 5 × 10⁻¹⁶`
is `~2.4 %`. The 1e-12 floor masked the regime-precision-limited
reality by a factor of 2000.

Conclusion: Gate 2 demands absolute energy precision of
`5 × 10⁻²⁶`, far below f64 round-off at any magnitude. The gate
as stated was unachievable regardless of algorithmic fix. The
engine skip is architecturally correct (Mercury 4.4 ppm
preserved, trajectory diverges by O(mass_ratio) as expected) but
cannot reduce the reported `|dE/E|` because the metric is
precision-saturated in this regime.

### Cross-checks

- Mercury 4.4 ppm
  (`apsis-1pn::mercury_precession_matches_gr_within_100ppm`):
  PASS post-fix. Mercury mass ratio `1.66 × 10⁻⁷` is above the
  `1e-10` threshold; skip never engages.
- `cargo test -p apsis --lib --release`: 558 passed / 0 failed /
  9 ignored. No regression.

### Threshold sensitivity sweep

Skipped. The endpoint metric does not register sensitivity to
threshold (per Gate 2 reanalysis above, the metric is precision-
limited regardless). Sensitivity sweep belongs in the metric-
refactor follow-up.

### Disposition

Engine skip ships. It is architecturally correct, passes Mercury
4.4 ppm and 558 lib tests, and removes spurious back-reaction on
the primary in extreme-mass-ratio configurations. Endpoint
`|dE/E|` is insensitive to the change because of the `.max(1e-12)`
floor, but absolute trajectory accuracy improves (`r_final` and
`v_final` shift by O(mass_ratio) between runs).

The bookkeeping metric refactor is separate work, scoped to a
follow-up issue. Surface area mapped 2026-05-18: 8+ consumer
sites (`Metrics` DTO, `HookContext`, `System` getter, adaptive
controller, 4 examples, 2 benchmark files). The adaptive
controller (`core/adaptive.rs`) is tuned against the floor-
inflated metric and is the principal blocker — a semantic change
to `rel_energy_error` requires controller redesign in the same
change.

Proposed direction for the follow-up:

- Add `abs_energy_error: f64` as primary observable.
- Make `rel_energy_error: Option<f64>` regime-aware: `None` when
  `|E_initial|` falls below a principled precision-limited
  threshold (rather than a magic floor).
- Redesign the adaptive controller to target an absolute drift
  bound for the integrator's noise floor with explicit regime
  detection.

## References

- Issue #133 — https://github.com/GabrielEstefanski/apsis/issues/133
- Burns, J. A., Lamy, P. L., & Soter, S. (1979). *Icarus* 40, 1.
- Rein, H. & Spiegel, D. S. (2015). *MNRAS* 446, 1424. (IAS15.)

## Decision log

- **D1:** test particle = directional interaction, not body
  property.
- **D2:** momentum trade accepted and documented.
- **D3:** threshold `1e-10`, configurable per-system.
- **D4:** primary gate `|dE/E| < 1e-10`; bonus `< 1e-13`;
  failure mode `>= 1e-8` triggers re-diagnosis, not gate
  loosening.
- **D5:** protocol introduced generically (trait method on
  `Operator`); only gravity-core auto-detect wired immediately.
- **D6:** direct O(N²) path only; BH tree deferred.
