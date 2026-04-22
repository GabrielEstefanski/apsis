# Picard Noise-Floor Convergence: Experiment Log

**Date:** 2026-04-22
**Subject:** IAS15 integrator in `gravity-sim`
**Baseline commit:** `ffc3489` (drift-metrics baseline merged into `develop`)
**Experimental branch:** `perf/ias15-picard-noise-floor` (not merged; null result)
**Outcome:** Hypothesis falsified. Change abandoned; documentation retained as design record.

---

## Abstract

We tested a single localised modification to the Picard predictor–corrector inside IAS15:
accept convergence when the iteration residual drops below `10 · PICARD_TOL` (the round-off
noise floor of the divided-difference recurrence), rather than requiring the strict
`< PICARD_TOL` threshold or the two-consecutive stagnation heuristic. The hypothesis, based
on REBOUND's implicit "stop when the machine can't do better" convergence criterion, was that
this would reduce `picard_iters_total` by 10–25% on typical scenarios without degrading
numerical quality.

The hypothesis was falsified. Across four scenarios (Kepler at *e* ∈ {0.5, 0.9, 0.99} and
the Pythagorean three-body problem), the change reduced Picard iterations by 0 to 15
(0.00% to 0.04% of total) while leaving every quality metric bit-identical to the baseline.
Criterion-measured wall time was unchanged in three of four scenarios and showed a small
(−2.5%, *p* < 0.05) reduction in kepler\_e05 that is ambiguously attributable to the change
versus run-to-run noise.

The falsification is diagnostic: the steady-state IAS15 Picard residual *reaches exactly
`0.0` in f64* through the combination of warmstart and the compensated-summation `b`
update, rather than hovering in the interval `(PICARD_TOL, 10·PICARD_TOL)`. Strict
convergence therefore already handles the vast majority of substeps and the
noise-floor branch has almost nothing to catch. This redirects future optimisation
work away from iteration-count reduction and toward **cost per iteration** — the
7-stage divided-difference update in `update_g_and_b` is the real hot path.

---

## 1. Motivation

IAS15 (Rein & Spiegel 2015) is a 15th-order adaptive Gauss-Radau integrator used in
`gravity-sim` as the high-fidelity integrator of choice for gravitational dynamics.
Its inner loop is a predictor–corrector Picard iteration that refines the power-series
coefficients `b[0..6]` at each adaptive sub-step until a convergence criterion is met.

Profiling inspection of the integrator showed that the Picard loop was running
approximately **4 iterations per accepted sub-step** across all scenarios, while the
original Rein & Spiegel paper reports typical values of 2–3 and REBOUND's reference
implementation is similarly economical. The gap raised the question: is our loop doing
unnecessary work at the round-off floor?

The immediate motivator was to match REBOUND-level efficiency without compromising the
energy-conservation quality that IAS15 is chosen for.

## 2. Experimental Infrastructure

All measurements were produced by the `gravity-sim` benchmark harness introduced in
commits `f5e24e1` through `ffc3489`. Its design is a methodological contribution in its
own right and is documented here as part of the experimental setup.

### 2.1 Determinism Enforcement

The bench entry point (`benches/ias15.rs::main`) forces the rayon global thread pool to
one worker before any scenario runs:

```rust
rayon::ThreadPoolBuilder::new()
    .num_threads(1)
    .build_global()
    .expect("rayon single-thread init");
```

This eliminates reduction-order non-determinism from parallel force evaluation, a
prerequisite for bit-exact metric comparison across runs.

### 2.2 Tiered Baseline Gate

The harness stores a versioned baseline in `benches/baselines/ias15.toml` and compares
every recorded metric against it on each bench run. Tolerances are typed:

| Tier | Metrics | Default tolerance | Semantics |
|------|---------|-------------------|-----------|
| Counter | substeps, rejections\_picard, rejections\_truncation, picard\_iters\_total, degraded\_total | `tol_abs = 0` | Integer quantities; expected bit-deterministic. Any change is a real behavioural shift. |
| Float | dt\_min, dt\_max, dt\_mean, dt\_p05, dt\_p50, dt\_p95, peak\_energy\_err, rel\_energy\_err\_rms, energy\_drift\_slope | `tol_factor` calibrated from observed jitter | Real-valued quantities; allow ULP-level drift from FMA / compiler reassociation. |

A hard upper cap `MAX_ALLOWED_TOL_FACTOR = 1.5` rejects any baseline with a tolerance
above 1.5× at load time, preventing silent drift through successive relaxations.

### 2.3 Recording Procedure

`IAS15_BENCH_UPDATE_BASELINE=1 cargo bench` triggers recording mode: each scenario is
run 10 times. The resulting per-metric distributions drive tolerance calibration:

* **Counter metrics** must have `min == max` across runs; otherwise the harness aborts
  with a determinism-violation error rather than widening the tolerance band.
* **Float metrics** get `tol_factor = max(1.0, 1 + 2·(max − min)/mean)`, capped at
  `MAX_ALLOWED_TOL_FACTOR`.

Baseline updates are intended to be explicit, reviewable actions: the resulting TOML
diff is committed alongside the code change that motivated it, so a reviewer sees
exactly which metrics moved and how much.

### 2.4 Scenarios

Four scenarios cover the IAS15 regimes relevant to `gravity-sim`:

| Scenario | N | Description | Intended stress |
|----------|---|-------------|-----------------|
| kepler\_e05 | 2 | Kepler orbit at *e* = 0.5, 100 orbits | Steady-state; warmstart maximally effective |
| kepler\_e09 | 2 | Kepler orbit at *e* = 0.9, 50 orbits | Moderate eccentricity shrink/grow cycle |
| kepler\_e099 | 2 | Kepler orbit at *e* = 0.99, 20 orbits | Controlled close encounter (pericentre at 0.01·a) |
| pythagorean | 3 | Burrau 1913 three-body, *t* ∈ [0, 10] | Chaotic multi-body with close encounters |

All use `ε = 1e-9` (IAS15 target relative error on the dominant truncation term).

### 2.5 Metrics Captured

Per-scenario metrics recorded on every run:

**Controller counters (Tier 1):**
- `substeps` — accepted adaptive sub-steps
- `rejections_picard`, `rejections_truncation` — rejection counts by cause
- `picard_iters_total` — cumulative Picard iterations across all attempts
- `degraded_total` — DT\_MIN-floor or deadline escapes

**`dt` distribution (Tier 2):**
- `dt_min`, `dt_max`, `dt_mean` — distribution extremes and centroid
- `dt_p05`, `dt_p50`, `dt_p95` — lower-tail shape, median, upper-tail shape

**Numerical quality (Tier 2):**
- `peak_energy_err` — `max |δE / E₀|` over the run
- `rel_energy_err_rms` — `sqrt(mean(|δE/E₀|²))` over all sampled sub-steps; penalises
  sustained error rather than isolated spikes
- `energy_drift_slope` — least-squares slope of `|δE/E₀|(t)` vs `t`; positive slope
  indicates secular envelope growth (drift), near-zero indicates oscillatory error

Criterion timing is recorded separately through its own baseline system
(`target/criterion/`), independent of the metric gate.

## 3. Hypothesis

**H0:** A convergence branch inserted in the Picard loop at `residual < 10 · PICARD_TOL`,
between the strict check (`residual < PICARD_TOL`) and the stagnation guard, will
reduce `picard_iters_total` by 10–25% across all four scenarios without degrading
any Tier 2 quality metric.

**Rationale:** The strict threshold `PICARD_TOL = 1e-16` sits at f64 machine epsilon.
If the typical residual at convergence is dominated by round-off in the Newton
divided-difference recurrence, it will oscillate in the band `(PICARD_TOL, 10·PICARD_TOL)`
— never reaching strict convergence, never triggering stagnation (two-consecutive-worse
rule), and costing an extra 1–2 iterations per sub-step waiting for one of the two
termination conditions to fire. REBOUND's implementation treats this band as convergent
implicitly.

**Falsification criterion:** if the change produces a measurable reduction in
`picard_iters_total` (> 5%) and holds all Tier 2 metrics within calibrated tolerance,
H0 is supported and the change is merged. If iteration reduction is < 1% or any quality
metric degrades, H0 is falsified and the change is abandoned.

## 4. Methodology

### 4.1 The Change

A single constant and a single branch in `Ias15::picard_loop_inner`:

```rust
/// Round-off noise floor of the divided-difference recurrence.
const PICARD_NOISE_FLOOR: f64 = 10.0 * PICARD_TOL;

// (inside the Picard iteration loop, between strict convergence and stagnation:)

if residual < PICARD_TOL {
    // Strict convergence (pre-existing)
    restore_xv(bodies, x0, v0);
    return (true, residual, iters);
}

// NEW: effective convergence at the round-off noise floor.
if residual < PICARD_NOISE_FLOOR {
    restore_xv(bodies, x0, v0);
    return (true, residual, iters);
}

// Stagnation guard (pre-existing, unchanged)
if iter >= 2 && residual > last_residual {
    no_improve += 1;
    if no_improve >= 2 {
        restore_xv(bodies, x0, v0);
        return (false, residual, iters);
    }
} else {
    no_improve = 0;
}
```

No other modifications. `decide_dt`, the outer retry loop, `optimal_dt`, and
`warmstart_b` were untouched.

### 4.2 Evaluation Procedure

1. On clean `develop` at commit `ffc3489`, record a fresh baseline with
   `IAS15_BENCH_UPDATE_BASELINE=1 cargo bench`. Confirm the harness reports
   bit-exact determinism across all 10 runs of all 4 scenarios.
2. Apply the change described in §4.1.
3. Run `cargo test --release ias15` — IAS15 unit and regression tests must pass
   (11 tests, covering Kepler energy drift over 100 orbits, high-eccentricity
   Kepler, Pythagorean close encounters, and the decide\_dt decision table).
4. Run `cargo bench --bench ias15 -- --test` (validation mode). The gate should
   report any per-metric change relative to the baseline.
5. Run `IAS15_BENCH_UPDATE_BASELINE=1 cargo bench` to regenerate the baseline
   under the new code. Inspect the TOML diff to read the before/after for every
   metric in every scenario.
6. Run `cargo bench` (timing mode) to measure Criterion-level wall-time impact.
7. Evaluate against the falsification criterion in §3.

## 5. Results

### 5.1 Metric Changes

Tier 1 counters (full table, absolute values):

| Metric | Scenario | Baseline | Post-change | Δ | % change |
|--------|----------|---------:|------------:|---:|--------:|
| `picard_iters_total` | kepler\_e05 | 59505 | 59505 | **0** | 0.00% |
| `picard_iters_total` | kepler\_e09 | 53885 | 53882 | **−3** | −0.006% |
| `picard_iters_total` | kepler\_e099 | 35042 | 35027 | **−15** | −0.04% |
| `picard_iters_total` | pythagorean | 6362 | 6361 | **−1** | −0.016% |
| `substeps` | (all) | — | — | **0** | 0.00% |
| `rejections_picard` | (all) | — | — | **0** | 0.00% |
| `rejections_truncation` | (all) | — | — | **0** | 0.00% |
| `degraded_total` | (all) | — | — | **0** | 0.00% |

Tier 2 floats: **every single metric bit-identical** to baseline across all four
scenarios. Same u64 bit pattern. No degradation, no improvement — zero trajectory
change in 99.96% of sub-steps, and the remaining 0.04% produced acceleration
profiles equivalent to f64 precision.

### 5.2 Criterion Timing

| Scenario | Median time (post-change) | Change vs pre-change | *p* value | Verdict |
|----------|--------------------------:|--------------------:|---------:|---------|
| kepler\_e05 | 338.94 µs | −2.51% [−4.02, −0.93] | 0.00 | **Improvement** |
| kepler\_e09 | 342.16 µs | −0.13% [−1.85, +1.68] | 0.89 | No change |
| kepler\_e099 | 335.08 µs | +0.10% [−1.61, +1.81] | 0.91 | No change |
| pythagorean | 427.55 µs | +1.57% [−0.16, +3.36] | 0.07 | No change |

The kepler\_e05 improvement is statistically significant but inconsistent with the
accompanying Tier 1 data — kepler\_e05 saw *zero* change in `picard_iters_total`,
`substeps`, or any Tier 2 metric. The most defensible interpretation is that the
2.5% reading reflects Criterion's sample variance rather than the code change. A
re-run on cold system state would be needed to attribute it conclusively, but the
other three scenarios show no effect and collectively constitute stronger evidence.

### 5.3 Falsification

By the criteria of §3: iteration reduction is 0 to 0.04% (vs 10–25% predicted);
quality is preserved (bit-exact); timing impact is within noise in three of four
scenarios. **H0 is falsified** on the iteration-reduction dimension. The change is
numerically safe but has no measurable impact on the quantity it was designed to
affect.

## 6. Interpretation — Why the Noise-Floor Band Is Empty

The result is initially surprising but mechanically explainable once the
compensated-summation path is considered.

Inside `update_g_and_b`, the propagation of Δg into `b` uses Neumaier summation
(`add_cs`). When the warmstart-produced `b` is near steady state, the divided
difference of the accelerations at the 7 Gauss-Radau nodes yields a Δg that is
numerically exactly zero (not "small") because:

1. The node accelerations computed at consecutive Picard iterations differ only
   in bits already below the compensated-summation residual.
2. `add_cs(&mut b, &mut csb, 0.0)` updates neither `b` nor `csb` — the new `b[i][6]`
   bit pattern equals the old `b[i][6]` bit pattern.
3. Therefore `residual = max|Δb[6]| / max|a₀| = 0.0 / max|a₀| = 0.0`.
4. `0.0 < PICARD_TOL` — strict convergence fires.

The inhabitants of the band `(PICARD_TOL, 10·PICARD_TOL)` are only those sub-steps
where warmstart is imperfect *and* the first iteration has not yet cancelled the
imbalance to exact zero. This is a narrow regime — empirically 0.04% of sub-steps
in the most demanding scenario (kepler\_e099) and 0% in the smooth one
(kepler\_e05).

The 4-iterations-per-sub-step average is therefore *not* waste at the noise floor
— it is the natural quadratic/cubic convergence of Picard on a 15th-order
implicit scheme, from an initial state that the warmstart places "close but not
identical" to the fixed point. An integrator that converged in 2–3 iterations
either started from a better initial guess, had a less stringent convergence
criterion, or operated under different arithmetic semantics.

## 7. Conclusion

Picard iterations in our IAS15 implementation terminate at the arithmetically
strongest possible condition (`residual = 0.0` exactly) for the majority of
sub-steps under realistic warmstart. The hypothesised round-off band between
strict convergence and stagnation detection is populated at the 0.04% level — too
sparse to yield a useful optimisation.

The change tested here is numerically correct and matches the REBOUND
"stop-when-the-machine-can't-do-better" philosophy, but in our concrete
implementation it has no measurable effect because the machine *already* can't do
better — the strict test catches that case first.

## 8. Implications for Future Work

Iteration count is not the optimisation target. The hot path in IAS15 at
`gravity-sim`'s current scenarios is **cost per iteration**, not iteration count:

* **`update_g_and_b`** (7 stages × N bodies × compensated-summation updates to `b`)
  runs approximately 4 × 7 × N times per accepted sub-step. For N = 2 at
  4 iters/sub-step this is 56 `add_cs` calls plus 7 Newton divided-difference
  updates per sub-step. This is where any real per-sub-step speedup has to come
  from.
* **Layout**: `b[i][k]` is currently `Vec<[(f64, f64); 7]>` (struct-of-arrays
  along the outer axis, tuple-of-arrays on the inner). SIMD across the body
  axis would require either a transposed layout or a `(bx, by)` split; worth
  measuring.
* **Force evaluation scaling**: at N = 2–3 the force model is a trivial share of
  per-sub-step time. To identify the crossover at which force dominates, add a
  100-body or 1000-body cluster scenario to the harness. Until then the
  integrator-vs-force split is speculative.
* **REBOUND cross-check**: the 2–3 iterations-per-sub-step figure quoted in
  Rein & Spiegel 2015 may reflect a different convergence measurement
  (e.g. relative to `b` magnitude rather than `a₀` magnitude) or a different
  compensated-summation discipline. A direct cross-run of identical initial
  conditions in both implementations would resolve the apparent discrepancy
  without further speculation.

## 9. Threats to Validity

* **Platform-specific.** All measurements were taken on a single machine
  (Windows 11, single-threaded rayon, default Rust toolchain). Fused-multiply-add
  availability, compiler reassociation choices, and allocator behaviour differ
  across platforms; the bit-exact baseline is machine-specific by design.
  Cross-platform reproducibility would require either disabling FMA explicitly
  or relaxing Tier 2 tolerances — both deferred to when platform portability
  becomes a concrete requirement.
* **Scenario coverage.** The four scenarios cover small-N (N ∈ {2, 3}) and short
  integration windows (up to 100 orbits for the smoothest case). Regimes not
  covered include: large-N Barnes-Hut-dominated force evaluation, extreme
  long-duration integrations (> 10⁶ sub-steps where f64 noise accumulation
  matters), and scenarios with body spawn/despawn (which would invalidate
  warmstart and potentially populate the noise-floor band more densely). The
  conclusions of this experiment apply confidently within the tested regime and
  hypothetically outside it.
* **Criterion variance.** The kepler\_e05 −2.5% reading that could not be
  attributed to any Tier 1 or Tier 2 change suggests Criterion's inherent
  sample variance at the few-hundred-microsecond scale is roughly 2–3%. This
  is adequate for detecting optimisations > 5% and insufficient for finer
  effects — treat sub-3% differences here as inconclusive unless corroborated
  by a Tier 1 change.
* **Null hypothesis scope.** The experiment falsifies "adding the noise-floor
  branch reduces iteration count meaningfully", not "the branch is incorrect".
  The branch would be correct to adopt under different numerical conditions
  (looser `PICARD_TOL`, different recurrence, different warmstart) where the
  strict check is less frequently reachable.

## 10. Reproducibility

**Hardware/OS:** Windows 11 Pro for Workstations, x86-64.
**Toolchain:** Rust stable (as of 2026-04-22), edition 2024.
**Code state:**
* Baseline `develop` HEAD: `ffc3489` (drift metrics merged).
* Experimental diff: two insertions in `src/physics/integrator/ias15.rs`
  (one constant, one `if` block) — reproduced verbatim in §4.1.

**Invocation:**
```bash
# Baseline recording (10 runs per scenario, writes benches/baselines/ias15.toml)
IAS15_BENCH_UPDATE_BASELINE=1 cargo bench --bench ias15

# Validation (compares a single run against the recorded baseline)
cargo bench --bench ias15 -- --test

# Timing (Criterion's statistical harness)
cargo bench --bench ias15
```

**Determinism checks:** All counter metrics recorded with `tol_abs = 0` across
both the baseline and the post-change run. All float metrics recorded with
`tol_factor = 1.0` (bit-exact) from 10-run calibration. This confirms rayon
single-thread enforcement and the integrator's arithmetic determinism were both
in effect for the measurements reported here.

---

## Appendix A — Full Baseline Diff

The `benches/baselines/ias15.toml` diff between the two states, showing exactly
which metrics moved:

```diff
@@ kepler_e09 @@
-picard_iters_total = { value = 53885.0, tol_abs = 0.0 }
+picard_iters_total = { value = 53882.0, tol_abs = 0.0 }

@@ kepler_e099 @@
-picard_iters_total = { value = 35042.0, tol_abs = 0.0 }
+picard_iters_total = { value = 35027.0, tol_abs = 0.0 }

@@ pythagorean @@
-picard_iters_total = { value = 6362.0, tol_abs = 0.0 }
+picard_iters_total = { value = 6361.0, tol_abs = 0.0 }
```

kepler\_e05 produced zero diff. All other metrics (dt distribution,
`peak_energy_err`, `rel_energy_err_rms`, `energy_drift_slope`, rejection counts,
`substeps`) produced zero diff across all four scenarios.

## Appendix B — Related Observations

**RMS vs peak ratio as a regime indicator.** The baseline data captured in this
experiment illustrates a property useful for future diagnosis:

| Scenario | peak_energy_err | rel_energy_err_rms | peak / RMS |
|----------|----------------:|-------------------:|-----------:|
| kepler\_e05 | 2.66e-15 | 1.08e-15 | 2.5× |
| kepler\_e09 | 1.07e-14 | 1.86e-15 | 5.7× |
| kepler\_e099 | 8.53e-14 | 2.25e-14 | 3.8× |
| pythagorean | 2.06e-12 | 4.45e-13 | 4.6× |

Ratios near 1 indicate error that is near-constant across the integration;
larger ratios indicate error concentrated in transient events (close encounters,
chaotic bursts). kepler\_e099 and pythagorean both show moderate concentration
despite very different dynamical regimes. This ratio is a candidate summary
metric for future comparative work across integrators or tunings.

**`dt` percentile spread as a close-encounter intensity indicator.** The ratio
`dt_p50 / dt_p05` scales with how much the controller compresses `dt` in rare
regimes:

| Scenario | dt\_p50 / dt\_p05 |
|----------|------------------:|
| kepler\_e05 | 1.95× |
| kepler\_e09 | 6.22× |
| kepler\_e099 | 31.1× |
| pythagorean | 97.5× |

This captures the progression from smooth (kepler\_e05) through controlled
pericentre (kepler\_e099) to chaotic multi-encounter (pythagorean) with an
interpretable single number. Not formally proposed as a published metric — but
useful as an informal health-check when evaluating a controller change.

## Appendix C — Why This Document Exists

Null results in numerical software engineering are rare in published literature
but valuable: they document which optimisation hypotheses are *not* productive,
save future workers from re-running the same dead ends, and surface diagnostic
insights (like §6 above) that would otherwise live only in tribal knowledge.

The gravity-sim project treats this document as the permanent artefact of the
experiment. The code change itself was not merged; the methodology and
numerical findings are what we keep.
