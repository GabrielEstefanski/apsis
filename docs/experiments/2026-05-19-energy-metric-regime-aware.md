---
date: 2026-05-19
status: a priori (protocol declared before any code lands)
issue: '#146'
closes: '#133'
---

# Energy bookkeeping metric — regime-aware refactor

## Context

`update_energy_tracking` (`crates/apsis/src/core/system/step.rs`)
normalises against `denom = baseline.abs().max(1e-12)`. When
`|E_initial|` is smaller than the `1e-12` floor (extreme mass
ratio, near-cancellation, low total energy), the reported
`rel_energy_error` is `dE_absolute / 1e-12`, not
`dE_absolute / |E_initial|`. The reading is artificially scaled
by `|1e-12 / E_initial|` and stops measuring physics.

The same pattern exists in `update_angular_momentum_tracking`
(`step.rs`).

Issue #133 surfaced this on `radiation_dust.py`:
`|E_initial| ≈ 5 × 10⁻¹⁶`, reported `|dE/E| = 1.16 × 10⁻⁵`,
true honest relative drift `≈ 2.4 %`. The floor masked the
regime-precision-limited reality by a factor of 2000. The Issue
#133 engine fix (PR #147, directional back-reaction suppression)
removed spurious back-reaction on the primary but did not change
the reported `|dE/E|` because the metric itself was the bug.

This refactor replaces the floor with a regime-aware metric.
Closes #133.

## Hypothesis

Relative drift `(X - X₀) / |X₀|` is well-defined only when
`|X₀| ≥ √ε_machine`. Below that threshold the f64 round-off
noise floor dominates any signal, and the metric reports
amplified rounding artefact rather than physical drift.

Exposing the metric as `Option<f64>` (`None` in the
precision-limited regime) and adding `abs_energy_error: f64`
as a primary observable produces a metric semantically valid
in every regime. The adaptive controller, switched to operate
on `rel` when available and to disable feedback (fixed dt)
otherwise, stops generating spurious dt jitter driven by
rounding noise.

## Model

### Public types

```rust
/// Conditioning floor below which a relative metric (X - X₀)/|X₀| is
/// undefined: f64 round-off dominates signal when |X₀| < √ε.
pub const MIN_RELATIVE_DENOMINATOR: f64 = f64::EPSILON.sqrt();
// ≈ 1.49e-8

/// Energy / angular-momentum bookkeeping on Metrics.
pub struct Metrics {
    pub energy: f64,
    pub initial_energy: f64,
    pub abs_energy_error: f64,
    pub rel_energy_error: Option<f64>,

    pub angular_momentum_z: f64,
    pub initial_angular_momentum_z: f64,
    pub abs_angular_momentum_error: f64,
    pub rel_angular_momentum_error: Option<f64>,
    // ... existing fields ...
}
```

### Shared helper

```rust
// crates/apsis/src/core/system/regime.rs
pub(crate) fn regime_aware_rel(delta: f64, baseline: f64) -> Option<f64> {
    if baseline.abs() < MIN_RELATIVE_DENOMINATOR {
        None
    } else {
        Some(delta / baseline.abs())
    }
}
```

Signature takes `delta` (the already-computed `current - initial`)
rather than `value`, so the helper cannot be miscalled as
`(value - initial - initial) / initial`. Used by both
`update_energy_tracking` and `update_angular_momentum_tracking`.

### Controller mode

```rust
// crates/apsis/src/core/adaptive.rs
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeedbackMode {
    /// Relative-error feedback is well-conditioned;
    /// dt is adjusted against `target_rel_energy_error`.
    Active,
    /// |E_initial| is in the precision-limited regime;
    /// controller disables feedback and runs fixed dt.
    DisabledPrecisionLimited,
}
```

Controller decision each step:

```text
match rel_energy_error {
    Some(rel) => adjust dt against target_rel_energy_error;
                 feedback_mode = Active;
    None      => keep dt = user_dt; feedback_mode =
                 DisabledPrecisionLimited;
}
```

`FeedbackMode` is exposed via `AdaptiveStats.feedback_mode` for
telemetry. Default behaviour requires no user configuration:
regime is auto-detected from `|E_initial|`.

### Raw-observable accessors

`System` gains public accessors for the raw energy and angular-
momentum observables. Power users who want a custom normalisation
compose `energy()` and `initial_energy()` themselves; the
convenience metric remains for the default 80 % case:

```rust
impl System {
    pub fn energy(&self) -> f64;             // exists
    pub fn initial_energy(&self) -> f64;     // NEW
    pub fn abs_energy_drift(&self) -> f64;   // NEW
    pub fn energy_delta(&self) -> Option<f64>;  // CHANGED sig

    pub fn lz(&self) -> f64;                       // exists
    pub fn initial_lz(&self) -> f64;               // NEW
    pub fn abs_lz_drift(&self) -> f64;             // NEW
    pub fn lz_delta(&self) -> Option<f64>;         // CHANGED sig
}
```

## Documented contract

Relative drift `(X - X₀)/|X₀|` is defined only when
`|X₀| ≥ √ε_machine`. Below that threshold the metric is
undefined and the library reports `None`. Absolute drift
`X - X₀` is always reported and remains the primary
observable across regimes.

This goes into:

- `core/system/regime.rs` rustdoc on the helper.
- `core/system/metrics.rs` rustdoc on the Metrics fields.
- `paper.md` near the conservation-claim paragraph (paper voice
  audit; not blocking this PR).

## Protocol

### Baseline

Run the following on `develop` tip and record current outputs:

1. `examples/radiation_dust.py` — `|dE/E|` reading (expect
   `1.200e-5`, the post-#147 value).
2. `examples/kepler_2body.rs` — `|dE/E|` reading (non-degenerate;
   expect existing IAS15-Kepler typical value).
3. `examples/figure_eight.rs` — `|dE/E|` reading (non-degenerate).
4. Mercury 1PN precession gate
   (`apsis-1pn::mercury_precession_matches_gr_within_100ppm`) —
   PASS reading.

Capture as `validation/energy-metric/baseline-pre-fix.{csv,txt}`.

### Implementation gate

Six commits:

1. **`regime` module + helper + threshold constant.**
   `core/system/regime.rs` with `regime_aware_rel` +
   `MIN_RELATIVE_DENOMINATOR`. Unit tests at module level.
2. **Metrics + accumulators.**
   `update_energy_tracking` / `update_angular_momentum_tracking`
   call the helper. `Metrics` DTO gains `abs_*_error` and
   changes `rel_*_error` to `Option<f64>`. `System` field
   `initial_energy: Option<f64>` becomes accessible via new
   `initial_energy() -> f64` (panics if `None`; callers go
   through `try_initial_energy()` for the pre-first-step case).
3. **Adaptive controller + `FeedbackMode`.**
   `core/adaptive.rs` consumes `Option<f64>` from
   `rel_energy_error`. `AdaptiveStats.feedback_mode` exposed.
   `target_rel_energy_error: 1e-6` setpoint retained; no
   `target_abs_energy_error` introduced (controller disables in
   degenerate regime).
4. **HookContext + Stats propagation.**
   `HookContext.rel_*_error: Option<f64>`. `Stats` DTO mirrors
   Metrics changes. Python binding `Record.stats` re-exports
   updated shape.
5. **Examples + benchmarks.**
   `figure_eight.rs`, `kepler_2body.rs`,
   `pythagorean_close_encounter.rs`, `solar_system_long.rs` —
   `Option<f64>` handling in `println!`. Benchmark
   `runner.rs:79` switches to `abs_energy_error`.
6. **Python binding + tests + release note.**
   PyO3 accessors mirror shape. `apsis/_native/__init__.pyi`
   stubs updated. Smoke tests cover `None` case. Release note
   under `docs/releases/`.

### Decision gates

- **Gate 1 (baseline matches).** All four baseline values
  reproducible within 1 ULP.
- **Gate 2 (regime-aware semantics).** New regime tests pass:
  - dust scenario returns `None` for `rel_energy_error`,
    finite value for `abs_energy_error`.
  - kepler scenario returns `Some(_)` for `rel_energy_error`,
    matches existing value within 1 ULP.
  - dust scenario reports `FeedbackMode::DisabledPrecisionLimited`;
    kepler reports `Active`.
- **Gate 3 (no regression in conditioned regimes).** Mercury
  4.4 ppm gate PASS, REBOUND parity portfolio (Kepler /
  figure8 / pythagorean / retrograde / mercurius) unchanged at
  1 ULP for all gated invariant metrics, `cargo test -p apsis
  --lib --release` green.
- **Gate 4 (controller behaviour).** Non-degenerate runs show
  unchanged dt selection (within 1 ULP of pre-fix values);
  degenerate runs (dust) show `dt = user_dt` constant (no
  jitter). Telemetry exposes `feedback_mode` correctly per
  regime.
- **Failure mode.** Any regression in Gate 3 → diagnose, do
  not loosen the regression check. Any unexpected `None` in
  non-degenerate scenarios → threshold or helper has a bug.

## Out of scope

- Records v0.2 Diagnostic frames (`d_energy_rel`, `d_lz_rel`).
  Fold into #143 (v0.3 FORMAT_VER 3 schema work). NaN markers
  rejected as an encoding choice — use explicit tagged
  optional encoding when the format bump happens.
- Linear momentum tracking. Degeneracy is the default case for
  COM-centered systems (`|P_initial| ≈ 0` by construction);
  adding a third regime-aware metric with always-on weirdness
  is low value.
- `paper.md` conservation-claim wording update. Defer to paper
  voice audit. Mention in release note.

## Expected outcomes

### Default (non-degenerate) systems

User-visible behaviour unchanged at the convenience level.
`sys.energy_delta` returns `Some(x)` with the same `x` as
pre-fix. `Stats` reports same numbers. Benchmarks unchanged.

### Degenerate systems (dust, near-cancellation)

- `sys.energy_delta` returns `None`. Python users handle via
  `if (d := sys.energy_delta) is not None: ...` or treat as
  `nan` if they want REBOUND-style behaviour.
- `sys.abs_energy_drift` returns the honest absolute drift
  (e.g. `~1e-17` for dust scenario at noise floor).
- `sys.feedback_mode` returns `DisabledPrecisionLimited`.
- Adaptive controller does not jitter dt; fixed at `user_dt`.

### Paper claim posture

"Energy conservation to N × 10⁻ᵏ" claims become
regime-qualified. Default IAS15-Kepler claims unchanged. Dust
or near-cancellation regimes documented as
"precision-limited; absolute drift only".

## Results

### Baseline (pre-fix)

Captured against `develop` tip (commit `4ddef70` post-#147 merge):

| Scenario | Regime | Metric | Value | Expected post-fix |
|---|---|---|---|---|
| `radiation_dust.py` | degenerate (`\|E_initial\| ≈ 5e-16`) | `\|dE/E\|` | `1.200 × 10⁻⁵` | `None` (regime degenerate); abs drift `~1e-17` |
| `kepler_2body` | well-conditioned (`\|E_initial\| ≈ 0.5`) | `\|dE/E\|` | `3.775 × 10⁻¹⁵` | `Some(3.775e-15)` ±1 ULP |
| `mercury_precession_gate` | well-conditioned | 4.4 ppm vs GR | PASS | PASS |

**Gate 1 — PASS.** Dust value matches post-#147 develop tip;
Kepler value matches IAS15 noise floor expectation; Mercury gate
green. Diagnosis is stable; ready for implementation.

Raw stdout under `validation/energy-metric/`.

## References

- Issue #146 — `https://github.com/GabrielEstefanski/apsis/issues/146`
- Issue #133 — original report; closed by this work.
- PR #147 — directional back-reaction suppression (companion
  architectural fix landed separately).
- Higham, N. J. (2002). *Accuracy and Stability of Numerical
  Algorithms*, 2nd ed. §1.4, conditioning of division.

## Decision log

- **D1:** `abs_energy_error: f64` (always) +
  `rel_energy_error: Option<f64>` (regime-aware). Mirror for
  angular momentum.
- **D2:** `MIN_RELATIVE_DENOMINATOR = √f64::EPSILON ≈ 1.49 × 10⁻⁸`.
  Named for conditioning role, not physics threshold.
- **D3:** Controller uses `rel` when `Some`; disables feedback
  (fixed `user_dt`) when `None`. `FeedbackMode` enum exposed
  via `AdaptiveStats`. No `target_abs_energy_error` introduced.
- **D4:** Shared `regime_aware_rel(delta, baseline)` helper.
  Same threshold for energy and angular momentum.
- **D5:** Pre-beta breaking changes accepted. No backward
  compat shims. Migration touches `Metrics`, `Stats`,
  `HookContext`, `System` accessors, adaptive controller,
  4 Rust examples, 1 benchmark file, Python binding,
  Python stubs, smoke tests.
- **D6:** Market-standard escape hatch: `System.initial_energy()`
  and `System.initial_lz()` exposed publicly so power users
  can compose their own normalisation.
- **D7:** Records v0.2 Diagnostic frames out of scope; tracked
  in #143 for FORMAT_VER 3.
- **D8:** Linear momentum tracking out of scope (degeneracy is
  default for COM-centered systems).
