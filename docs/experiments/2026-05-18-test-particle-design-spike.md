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
orbits. Capture as `validation/test-particle/baseline-pre-fix.csv`.

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
`validation/test-particle/`.

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
