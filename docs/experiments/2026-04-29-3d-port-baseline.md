# 3D port — physics regression baseline

**Date:** 2026-04-29
**Branch:** `feat/3d-port` at commit `df7a0df`
**Platform:** Windows 11, MSVC 14.43, libm via std, x86_64-pc-windows-msvc
**Toolchain:** stable Rust (workspace `Cargo.lock` resolved)

## Purpose

The 3D port migrates the physics stack from `(f64, f64)` to `Vec3`
across four phased commits (Body fields → kernel/acc API → IAS15 buffers
→ orbital elements). Every commit adds the third spatial component
without altering algorithmic structure, and every commit is staged
under the same regression contract:

> For any input with `z = vz = 0`, every numeric measurement reported by
> the test suite must be **less than or equal to** its pre-port value.
> Errors may decrease (FMA fusion or autovectorisation may give
> incidental gains in ULPs), but **may not increase**, even within the
> assertion threshold.

The reduction-order regression caught during commit 3 (Mercury 1PN
went from ~1.08 ppm to ~122 ppm when `m * dx * fac` was re-associated
to `(m * fac) * dx`) is exactly the failure mode this contract is
designed to detect. Threshold-passing alone would have masked it; a
locked baseline catches it.

## Captured values

Numbers are the actual measurements emitted by each test in release
or debug mode as run by the workspace gate suite. They are NOT
thresholds — thresholds are visible in the cited assertion sites.

| Gate | Mode | Measured | Assertion threshold | Source |
|------|------|----------|---------------------|--------|
| Mercury 1PN perihelion precession (`rel_err`) | release | **1.075879e−6** | 1e−4 | [crates/apsis-1pn/tests/mercury_precession_gate.rs](../../crates/apsis-1pn/tests/mercury_precession_gate.rs) |
| Newtonian Kepler closure (`drift` of ω over 300 orbits) | release | **1.049457e−14** | 1e−9 | same file, `baseline_newtonian_kepler_is_closed` |
| IAS15 Kepler `e=0.5` energy peak (`peak |δE/E₀|`) | debug | **2.664535e−15** | 1e−12 | [crates/apsis/src/physics/integrator/ias15.rs](../../crates/apsis/src/physics/integrator/ias15.rs) `ias15_kepler_energy_peak_and_monotonic_drift` |
| IAS15 Kepler `e=0.5` monotonic drift (`slope · t_final`) | debug | **6.145135e−16** | 1e−13 | same |
| IAS15 Kepler `e=0.5` linear-regression slope | debug | **−4.890143e−19** | (informative) | same |
| IAS15 high-eccentricity `e=0.9` peak energy | debug | **1.065814e−14** | 1e−11 | `ias15_kepler_high_eccentricity` |
| IAS15 Pythagorean (Burrau 1913) peak energy through close encounters | debug | **1.122362e−12** | 1e−11 | `ias15_pythagorean_energy_through_close_encounters` |
| Barnes–Hut vs direct (4-body symmetric layout) | debug | **0.0** for every component of every body | 1e−2 | [crates/apsis/src/physics/gravity/engine.rs](../../crates/apsis/src/physics/gravity/engine.rs) `barnes_hut_matches_exact_with_small_error` |

The BH-vs-exact `0.0` is genuine: the symmetric four-body layout
folds the BH traversal into pairwise interactions that match the
direct kernel bit-for-bit, with no ULP-level disagreement to
report. Any non-zero value after the IAS15 migration would indicate
a reduction-order shift in the BH path, even though commit 4 does
not touch BH.

## Acceptance contract for commit 4

The IAS15 buffer migration (`b`, `e`, `g`, `csb`, `csx`, `csv`,
`pic_*`, `snap_*` from `(f64, f64)` to `Vec3`; `BodyCoeffs` from
`[(f64, f64); 7]` to `[Vec3; 7]`; `predict_ias15` /
`predict_v_ias15` / `predict_order2` rewritten with `Vec3` inputs
and outputs) is accepted iff:

1. Every value above is reproduced **bit-equivalent or smaller** in
   absolute magnitude, on the same hardware, after the migration.
2. No new test fails. Every existing test continues to pass against
   its current threshold.
3. `cargo fmt --all -- --check`, `cargo clippy --workspace
   --all-targets -- -D warnings`, and the full `cargo test
   --workspace` are clean.

A measurement that **exceeds** the captured baseline — even by one
ULP, even within the assertion threshold — is treated as a
regression. The expected diagnostic path is to identify the
re-associated operation and restore the original reduction order
(per the playbook in
[2026-04-28-ias15-velocity-prediction-bug.md](2026-04-28-ias15-velocity-prediction-bug.md)
and the post-mortem comments inside `engine.rs::exact_eval`).

## How to reproduce

```bash
# Workspace unit tests (debug)
CARGO_TARGET_DIR=target-3d cargo test --workspace --lib --nocapture --test-threads=1 \
  2>&1 | grep BASELINE

# Release-mode physics gates
CARGO_TARGET_DIR=target-3d cargo test --release -p apsis-1pn --tests -- \
  --ignored --nocapture --test-threads=1 2>&1 | grep BASELINE
```

The `BASELINE` lines printed by the temporary `eprintln!` calls in
each assertion site are not committed to the repository — this
document is the evidence of record. After commit 4 lands, the
instrumentation is re-applied transiently to verify, then stripped.

## Post-migration verification (commit 4)

Re-running the same instrumentation on the same hardware after the
IAS15 buffer migration produced identical values to the baseline at
every measurement site:

| Gate | Baseline | Post commit 4 | Δ |
|------|----------|---------------|---|
| Mercury 1PN `rel_err` | 1.075879e−6 | 1.075879e−6 | bit-exact |
| Newtonian Kepler closure `drift` | 1.049457e−14 | 1.049457e−14 | bit-exact |
| IAS15 Kepler `e=0.5` peak | 2.664535e−15 | 2.664535e−15 | bit-exact |
| IAS15 Kepler `e=0.5` drift | 6.145135e−16 | 6.145135e−16 | bit-exact |
| IAS15 Kepler `e=0.5` slope | −4.890143e−19 | −4.890143e−19 | bit-exact |
| IAS15 high-e `e=0.9` peak | 1.065814e−14 | 1.065814e−14 | bit-exact |
| IAS15 Pythagorean peak | 1.122362e−12 | 1.122362e−12 | bit-exact |
| BH-vs-direct (8 components) | 0.0 | 0.0 | bit-exact |

Every quantitative gate produced its baseline value to the last
printed digit. The acceptance contract holds: the IAS15 buffer
migration introduced zero numerical drift for planar input.

The mechanism that delivered this is the discipline applied at every
re-association site:

* `b₆.length_squared()` is **not** used for the Picard residual or
  truncation-error norms. The hand-written `b.x*b.x + b.y*b.y +
  b.z*b.z` reproduces the pre-port `b.x*b.x + b.y*b.y` reduction
  followed by an addition of an exactly-zero `b.z*b.z` (IEEE-754
  exact additive identity) for `z = 0` inputs.
* `update_g_and_b`, `advance_state`, `warmstart_b`, and
  `recompute_g_from_b` write their three-axis arithmetic in the same
  scalar form as the original two-axis code, with the third axis
  appended; they do **not** consolidate into Vec3 ops. The cost is
  ~70 lines of explicit per-axis algebra; the benefit is the
  bit-exact preservation table above.
* The persistent integrator buffers (`b/e/g/csb/csx/csv/pic_*/snap_*`)
  start as `Vec3::ZERO` from `ensure_capacity`. With planar input
  (`z = vz = 0` on every body), every `.z` slot they accumulate is
  also exactly zero, so memory layout changes (24-byte rows where
  there used to be 16-byte rows) do not propagate into observable
  state.

## Not captured (rationale)

Tests that assert booleans or use `assert_relative_eq!(_, _, epsilon
= …)` rather than reporting a numeric magnitude do not contribute to
this baseline. They pass or fail; there is no value to compare
against. The 269 unit tests in `apsis::*` that gate symmetry,
sign, additivity, layout, and contract semantics are governed by
their own assertions and are checked for binary pass/fail.

The truncated-kernel continuity counter
(`kernel_continuity_counter_test.rs`) is a bijection-of-events test
without a free-running magnitude; it is binary-checked.
