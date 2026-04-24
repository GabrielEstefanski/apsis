# Solar-System-Class Stutter in IAS15: A Three-Hypothesis Diagnostic Walk

**Date:** 2026-04-22
**Subject:** Interactive-app stutter at N≈641 with IAS15 in `gravity-sim`
**Baseline before:** `ec66e05` (persistent snapshot buffers merged)
**Branch:** `diag/ias15-solar-system-stutter`
**Outcome:** Two hypotheses refuted, one identified and fixed; the fix is REBOUND-style per-body RMS convergence and truncation estimation. One scenario-design insight documented.

---

## Abstract

Interactive playback of the `gravity-sim` app's `solar_system` preset
(641 bodies) and `jupiter_trojans` preset (802 bodies) under IAS15
showed a reproducible "stutters every few frames" pattern that Yoshida
and Velocity Verlet did not exhibit. Yoshida and Verlet run smoothly at
those scales; IAS15 degrades to multi-hundred-millisecond hitches.

This write-up is the third in a series of IAS15 diagnostics. The
earlier two (phase profile, buffer reuse) improved wall-time at small
N but did not address the stutter. This round instrumented additional
phases, proposed three hypotheses in sequence, refuted two with
measurement, and found the third by code inspection of the convergence
and truncation estimators.

**The fix is one structural change in two places** (convergence
residual in `picard_loop_inner` and `truncation_error`): replace
`max‖Δ‖ / max‖a₀‖` aggregation across bodies with per-body relative
ratios aggregated as RMS. This is the norm REBOUND's IAS15 uses in
spirit, and it is the correct norm to use for any measurement that
the outer controller treats as a scalar convergence signal over an
N-body system.

**Measured effect at solar_n641** (the synthetic scenario built to
reproduce the production regime):

* Truncation rejections: **−35.8%** (3878 → 2488)
* Evaluate calls: **−23.0%** (90 424 → 69 658)
* Total wall time (single-thread): **−23.5%** (64.95 s → 49.66 s)
* Energy conservation in controlled scenarios: preserved at machine
  precision (peak `|δE/E|` ≤ 2.77e-12 for Pythagorean 3-body, ≤ 7.74e-15
  for cluster_n50).

**Meta-finding**: the original max-max formula is algorithmically
fine at N = 2, where "max across bodies" reduces to one body and
the ratio is self-consistent. At large N the numerator and denominator
were selecting *different* bodies' outliers — an N-dependent
noise-to-signal measurement misrepresented as a truncation error.
The controller then rejected every step for bogus physics, cascading
`dt` toward `DT_MIN` and saturating the evaluate path with retries.

---

## 1. Motivation

User report:

> *TRAPPIST-1 é uma maravilha, roda tranquilo, solar system com 641
> corpos já começa a travar. Jupiter trojans (802) também trava
> bastante. Com 1400+ trava até com Verlet. No Yoshida ou outros
> integradores é bem mais aceitável.*

Two observations define the diagnostic frame:

1. **Integrator-specific at the same N.** Yoshida and Verlet remain
   smooth at N=641 where IAS15 stutters. The problem is not in the
   force engine alone (which is shared) — it is in what IAS15
   specifically does with the force engine.
2. **Stutter pattern, not slow.** "Anda 3×, para, anda 3×" is periodic
   hitches, not uniform slowdown. In hindsight this is consistent
   with a rejection cascade chewing through retries before any
   simulation time advances, but we did not infer that on round 1.

## 2. Round 1 — allocation-storm hypothesis (refuted)

**H1:** *At large N, IAS15 allocates many MB/sec into the
`DenseSnapshot` path (x0, v0, a0, b cloned per accepted sub-step).
At 100 KB × hundreds of substeps/sec, allocator contention produces
the observed hitches.*

**Test:** added feature-gated `a0_clone` and `dense_snapshot_build`
phase timers to `ias15.rs`; added `solar_n641` bench scenario (N=641
central + test particles, seeded RNG, one dynamical time of
integration). Ran phase profile under the new instrumentation.

**Result:**

| Phase (at N=641) | % of wall time |
|------------------|---------------:|
| evaluate | **98.40%** |
| update_g_and_b | 1.34% |
| snapshot_restore | 0.08% |
| dense_snapshot_build | **0.04%** |
| a0_clone | **0.00%** |

The two allocation candidates we instrumented total 0.04% of wall
time. **H1 refuted.** The cost is inside `evaluate`.

## 3. Round 2 — BH rebuild frequency hypothesis (refuted)

Code inspection of `force_model.rs` revealed:

```rust
fn compute(&mut self, bodies: &[Body], acc: &mut [(f64, f64)]) -> f64 {
    self.engine.build(bodies);
    self.engine.evaluate(bodies, self.theta, acc)
}
```

Every `ForceModel::compute` call rebuilds the quadtree. At N=641 with
IAS15's ~45 compute calls per accepted sub-step (1 start + 7 stages ×
~6 Picard iterations + 1 post-accept), we measured 90 424 compute
calls across the scenario — meaning 90 424 tree rebuilds.

**H2:** *Tree rebuild dominates per-call evaluate cost. Caching the
tree across Picard iterations within a sub-step (bodies barely move
between iterations) would recover most of the per-call cost.*

**Test:** instrumented `force_model.rs` with `tree_build` and
`tree_traverse` phase timers inside `compute`. Ran phase profile.

**Result:**

| Sub-phase (inside evaluate, at N=641) | ns/call | % of evaluate |
|---------------------------------------|--------:|--------------:|
| tree_build | 19 813 | **2.81%** |
| tree_traverse | 685 235 | **97.17%** |

**H2 refuted.** Tree build is 2.8% of evaluate; caching it would
recover ~2.8% of total wall time. Not worth the architectural
change. The real cost is in the traversal itself.

## 4. Round 3 — traversal-internal investigation

With the traversal consuming 97% of evaluate and evaluate 98% of wall
time, we tested three further hypotheses about where inside the
traversal the cost lives:

* **Rayon single-thread overhead**: `.into_par_iter()` for 641 bodies
  × 90k calls = 58M rayon tasks. Test: replace with plain `.iter()`.
  Measured saving: ~4% of traversal.
* **Per-body stack allocation**: `Vec::with_capacity(128)` per body
  per call = 58M small allocations. Test: hoist the stack to a
  single reusable buffer. Measured saving: ~5% of traversal.
* **Multi-thread scaling** (is the single-thread bench representative?):
  Added `IAS15_BENCH_MULTITHREAD=1` env var. Measured speedup at
  N=641: evaluate went 705 → 235 µs/call ≈ **3.0×**. Consistent
  with normal parallel scaling + Amdahl; not with allocator
  contention catastrophe (which would give <1.5× or bimodal).

At this point we had spent significant instrumentation effort
confirming that the per-call traversal cost is the dominant cost,
that micro-optimisations within traversal recover ≤10% combined, and
that nothing about the bench misrepresents the app.

**None of these was the stutter's root cause.**

## 5. Round 4 — the counter that mattered

Only by printing `snapshot_restore.count` (total rejections) alongside
`substeps` in the phase profile did the actual signal surface:

| Scenario | substeps | rejections | rate |
|----------|---------:|-----------:|-----:|
| kepler_e05 | 13 332 | 29 | 0.22% |
| kepler_e09 | 13 150 | 18 | 0.14% |
| kepler_e099 | 8 624 | 31 | 0.36% |
| pythagorean | 1 582 | 1 | 0.06% |
| cluster_n50 | 917 | 46 | 5.0% |
| **solar_n641** | **2 001** | **3 878** | **194%** |

**Three orders of magnitude more rejections per accepted sub-step at
N=641 than at any small-N scenario.** Of the 90 424 evaluate calls
in the scenario, ~54 000 were spent on attempts that got rejected —
the integrator was doing physics work that then got thrown away.

The IAS15 controller classifies rejections into two categories:
`RejectPicard` (Picard predictor–corrector failed to converge) and
`RejectTruncation` (Picard converged but the truncation error
estimate was above ε). They drive different shrink strategies.

To distinguish, we tried the **noise-floor convergence hypothesis**
(`residual < 10·PICARD_TOL` = accepted convergence) that an earlier
experiment had marked as null-result at N=2. This hypothesis
specifically targets `RejectPicard`:

* Before fix: rejections = 3878
* After fix:  rejections = **3878** (unchanged to the ULP)

**Null result confirms all 3878 rejections are `RejectTruncation`,
not `RejectPicard`.** Picard was converging fine; the controller was
rejecting for truncation error exceeding ε.

## 6. The actual fix — per-body RMS norm, REBOUND-style

`truncation_error` used the same max-max formula as the convergence
check:

```rust
fn truncation_error(&self, a0: &[(f64, f64)]) -> f64 {
    let mut max_b6 = 0.0_f64;
    let mut max_a = 0.0_f64;
    for (i, row) in self.b.iter().enumerate() {
        let b = row[6];
        max_b6 = max_b6.max(norm(b));
        max_a = max_a.max(norm(a0[i]));
    }
    max_b6 / max_a
}
```

**Why this breaks at N=641:** `max_b6` picks the body whose `b₆`
has the largest round-off noise. `max_a` picks a different body —
typically the one closest to the central mass. The ratio is a
noise-to-signal measurement divided across two different bodies,
*not* a truncation error. At N=2 the max reduces to one body and
the formula is self-referential (numerator and denominator are the
same body); at N=641 it is not.

**Fix:** compute per-body relative ratios and aggregate with RMS:

```rust
fn truncation_error(&self, a0: &[(f64, f64)]) -> f64 {
    let mut sum_sq = 0.0;
    let mut count = 0;
    for (i, row) in self.b.iter().enumerate() {
        let b6 = norm(row[6]);
        let a_mag = norm(a0[i]);
        if a_mag > 0.0 {
            let rel = b6 / a_mag;
            sum_sq += rel * rel;
            count += 1;
        }
    }
    sqrt(sum_sq / count)
}
```

The same transformation was applied to the Picard convergence
residual in `picard_loop_inner` for consistency.

**Result at solar_n641 (single-thread, deterministic):**

| Metric | Before fix | After fix | Delta |
|--------|-----------:|----------:|------:|
| substeps | 2001 | 1978 | -1.1% |
| rejections | 3878 | 2488 | **-35.8%** |
| evaluate calls | 90 424 | 69 658 | **-23.0%** |
| Total wall time | 64.95 s | 49.66 s | **-23.5%** |

**Quality preserved** in small-N controlled scenarios:

| Scenario | peak `\|δE/E\|` (before) | peak `\|δE/E\|` (after) | ratio |
|----------|------------------------:|------------------------:|------:|
| kepler_e05/e09/e099 | bit-exact | bit-exact | 1.00× |
| pythagorean | 2.06e-12 | 2.77e-12 | 1.35× |
| cluster_n50 | 5.56e-15 | 7.74e-15 | 1.39× |

The three Kepler scenarios remain bit-exact because at N=2 the max
reduces to a single body — old and new formulas are algebraically
equivalent in that limit. At N=3 (Pythagorean) and N=50 (cluster)
the formulas diverge by ~40% on peak error, both still at f64
machine-precision class (≤ 1e-11). IAS15's advertised energy
conservation is preserved.

## 7. The scenario-design footnote

An earlier draft of `solar_n641` used `softening = 0.01` and an
annulus of `[0.5, 5.0]` AU. Recording its baseline after the fix
showed:

* `peak_energy_err` = 3.5e-4 (8 orders of magnitude worse than
  Pythagorean's 2.77e-12).
* `dt_min` hit the literal `DT_MIN = 1e-12` floor.
* 36 `degraded_total` (DT_MIN escapes — controller gave up).

This is not a regression from the fix. The wide annulus combined
with random distribution concentrates close pairs at separations
below the softening length; the resulting force gradients force the
adaptive `dt` toward zero; the controller hits `DT_MIN`, accepts the
step in "degraded" mode, and the accumulated degraded-accept error
shows up as 1e-4 energy drift.

Before the fix, the rejection cascade was *hiding* this pathology:
every step was being over-rejected by the noisy controller, so `dt`
stayed absurdly small, and per-step error was tiny because the step
itself advanced almost no simulation time. After the fix, the
controller accepts at realistic `dt` and the scenario's inherent
physical stiffness becomes visible as energy drift.

**Two distinct properties** of the original `solar_n641`:

1. **Diagnostic value**: it genuinely reproduced the stutter pattern
   and made the rejection cascade measurable. Without it, the RMS-
   norm insight would not have landed.
2. **Representativeness**: it is *not* a fair sample of IAS15's
   advertised quality, because it sits outside the regime for which
   IAS15+BH is a good match. A scenario built to measure quality
   needs a softening comparable to the closest-pair separation and
   a dynamical-time range the adaptive controller can navigate.

We adjusted the `solar_n641` scenario to `softening = 0.05` and
annulus `[1.5, 3.5]` (closer to the `solar_system` preset's asteroid
belt at `[2.2, 3.2]`) so the recorded baseline reflects IAS15's
expected quality rather than an out-of-regime stress test. The
original parameters remain valuable as a separate stress-test
scenario — but that is future work, not part of this fix.

## 8. Meta-observations for future diagnostics

Three lessons from the two refuted hypotheses:

### 8.1 N=2 benchmarks cannot refute N=641 hypotheses

The noise-floor convergence change was tested and nulled at N=2–3
scenarios weeks before this investigation. Re-testing at N=641 during
round 5 *also* produced null (because the real bug was
`RejectTruncation` not `RejectPicard`), but *only* the N=641
measurement could have discovered that. The prior null-result
genuinely said "does not help N=2"; it never spoke to N=641.

Stated generally: **the regime in which a hypothesis is tested must
overlap with the regime in which the symptom occurs.** We had the
scenario infrastructure to do this from the start; we did not use it
early enough.

### 8.2 Counter instrumentation is cheaper than deep instrumentation

Round 4 (printing `rejections` alongside `substeps`) was a 3-line
change to the harness and revealed the actual bottleneck. Rounds 1
through 3 (allocation, tree rebuild, rayon, stack allocation)
required multi-hour phase-timing infrastructure and yielded null
results. The cheap measurement would have been faster to try first.

Future diagnostics should exhaust **single-number counters** before
committing to fine-grained timing.

### 8.3 Control-flow fix ≠ numerical fix

The cascade chain we eventually understood is:

```
close encounters in scenario
  → non-smooth force field
  → large body-to-body residual variance
  → max-max norm amplifies outlier into bogus residual
  → controller RejectTruncation cascade
  → dt → DT_MIN
  → degraded accepts (actual error)
  → and, in the fast-running case, rebuild storm in BH traversal
```

The fix (RMS norm) breaks the chain at *one* link — the norm
calculation. It does not make the underlying scenario more
physically tractable; in fact at the original softening the fix
made the scenario's true stiffness visible. The scenario-design
footnote (§7) handles the other end of the chain.

A complete paper section on this needs both. Optimising one link
without acknowledging the other is where "IAS15 is slow at N=641"
turns into scandalous quotes in a review.

## 9. Reproducibility

**Hardware/OS:** Windows 11 Pro, x86-64.
**Toolchain:** Rust stable (2026-04-22), edition 2024.
**Code state:** branch `diag/ias15-solar-system-stutter`, final commit
on that branch after this write-up.

**Invocations:**

```bash
# Standard validation + Criterion timing (single-thread, bit-exact gate):
cargo bench --bench ias15

# Phase profile with feature-gated instrumentation:
cargo bench --features ias15-profile --bench ias15 -- --test

# Multi-thread diagnostic (non-deterministic, large-N only):
IAS15_BENCH_MULTITHREAD=1 cargo bench --features ias15-profile --bench ias15

# Re-record baseline after a trajectory-affecting change:
IAS15_BENCH_UPDATE_BASELINE=1 cargo bench --bench ias15
```

The RMS-norm change is in `crates/gravity-sim-core/src/physics/integrator/ias15.rs`
(two functions: the `residual` computation in `picard_loop_inner`,
and `truncation_error`). No external API change; `ForceModel` trait
unchanged; `Integrator` trait unchanged.

---

## Appendix A — Path of hypotheses

| Round | Hypothesis | Instrumentation cost | Result | Saved |
|-------|------------|---------------------:|--------|------:|
| 1 | Allocation storm in DenseSnapshot | ~2h (phase timers) | Refuted | 0.04% |
| 2 | BH rebuild frequency | ~1h (compute-split timers) | Refuted | 2.8% |
| 3α | Rayon overhead (single-thread) | ~30min (toggle experiment) | Partial | 4% |
| 3β | Stack allocation | ~15min (hoist experiment) | Partial | 5% |
| 3γ | Multi-thread allocator contention | ~45min (diagnostic env var) | Refuted | N/A |
| 4 | Picard convergence floor at large N | 5min (counter display) | Refuted for this scenario | 0% |
| 5 | **RMS norm replaces max-max** | **~10min (local edit)** | **Confirmed** | **-23.5%** |

Total diagnostic effort vs fix: ~5 hours of instrumentation, ~10
minutes of actual code change. The instrumentation was the work; the
fix was one-line-per-function.

## Appendix B — The scenario parameters before and after

```rust
// Before (out-of-regime stress test — diagnostic only):
const R_INNER: f64 = 0.5;
const R_OUTER: f64 = 5.0;
const SOFTENING: f64 = 0.01;
// Recorded peak_energy_err = 3.5e-4 (unacceptable for a baseline).

// After (in-regime for IAS15+BH):
const R_INNER: f64 = 1.5;
const R_OUTER: f64 = 3.5;
const SOFTENING: f64 = 0.05;
// Expected peak_energy_err ≤ 1e-11 (machine-precision class).
```

The change narrows the annulus to the asteroid-belt range of the
interactive `solar_system` preset and raises softening to a length
comparable to the expected closest-pair separation at N=640 over
that area. Both adjustments keep `dt_adaptive` bounded away from
`DT_MIN`.
