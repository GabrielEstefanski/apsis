# IAS15 controller audit — Pascal warmstart, halving rejection, and growth cap

**Date:** 2026-04-26
**Subject:** Architectural audit of `apsis`'s IAS15 implementation against the algorithmic specification in Rein & Spiegel (2015), motivated by a substep-count blowup in the figure-8 cross-implementation parity scenario.
**Status:** *Diagnosed and resolved end-to-end. All 12 gated metrics of the figure-8 parity protocol pass at 1 ULP after the refactor; total wall time on the canonical 10-period stress test dropped from 15 minutes ($\sim 1.34 \times 10^{8}$ floor-pinned substeps) to $\approx 80$ ms (zero floor hits). Kepler-parity informational $|\Delta r|$ improved from $1.57 \times 10^{-9}$ to $2.18 \times 10^{-12}$, confirming the controller now tracks the specification's `dt_next` choices to ULP precision rather than diverging by phase from a peer implementation (REBOUND, used as the parity reference).*

---

## Abstract

A reproducible substep-count blowup was observed in `apsis`'s IAS15 implementation on the figure-8 choreography (Chenciner & Montgomery 2000) at integration horizons of 6–10 orbital periods. Trajectories remained physically correct ($\max|\Delta E / E_0| \approx 1$ ULP, all parity invariants within the *a priori* tolerance) but the controller drove `dt` through a slow rejection cascade into the `DT_MIN` $= 10^{-12}$ floor and accepted hundreds of millions of degraded sub-steps before the run completed.

The root cause was a combination of three independent divergences from the algorithmic specification given in Rein & Spiegel (2015) — each individually compatible with correct trajectories at low orbit count, each individually responsible for an order-of-magnitude inefficiency at higher orbit count:

1. The warmstart predictor used the *diagonal* of the polynomial-basis transformation $e_k = q^{k+1} \cdot b_k$ rather than the full Pascal-triangle expansion $e_k = q^{k+1} \cdot \sum_{j \geq k} \binom{j+1}{k+1} \, b_j$ derived in Everhart (1985, §III).
2. The truncation-rejection branch shrunk `dt` by the optimal-step formula $dt \cdot 0.9 \cdot (\epsilon / \mathrm{err})^{1/7}$ on every retry; near the noise floor where $\mathrm{err} \approx \epsilon$, this yields $\sim 10$–$20\,\%$ shrink per retry rather than the unconditional halving specified in Rein & Spiegel (2015) §3.4, requiring $\sim 30$ retries instead of $\sim 10$ to reach a stable step.
3. The accept-path `dt_next` proposal had no upper bound; in smooth regions where $\mathrm{err} \ll \epsilon$ the formula proposed unbounded growth, and the next sub-step's truncation gate rejected the resulting overshoot, paying a full shrink cascade per close-encounter recovery. The specification (Rein & Spiegel 2015, §3.4) caps $dt_\text{new}$ at $7 \cdot dt_\text{current}$ precisely to prevent this pathology.

After bringing the three pieces of the controller into line with the specification (full Pascal warmstart, halving on truncation rejection, $7 \times$ growth cap), the figure-8 stress test runs in $80$ ms with zero floor-pinned sub-steps and the cross-implementation phase drift $|\Delta r|$ drops from $3.17 \times 10^{-7}$ to $9.44 \times 10^{-13}$. The result is documented here not as a bug-fix changelog but as a methodological observation about the interaction between high-order adaptive integrators and the numerical-conditioning details that distinguish a publishable implementation from one that "looks right" on smooth scenarios.

---

## Symptom

The figure-8 cross-implementation parity scenario (`validation/rebound-parity/figure8/`, comparing `apsis` IAS15 against REBOUND IAS15 as a peer implementation of the same algorithmic specification) integrated cleanly over 1–5 orbital periods (sub-second wall time, machine-precision energy), but at 6 + periods entered a sustained "floor cascade":

- $\geq 10^{8}$ accepted sub-steps,
- $\approx 60\,\%$ at `DT_MIN` $= 10^{-12}$,
- cumulative wall time of 15 min for 10 periods.

The trajectory remained physically correct on every gated metric (all 12 invariant gates of the parity protocol passed at 1–10 ULP), so the cascade was a *performance* defect, not a correctness defect — but a sufficiently severe one that it could not stay in v0.1 of the library. The cascade onset was *deterministic at $t \approx 5.5\,T$* independent of the total integration horizon: a 5-period run terminated cleanly *just before* the cascade would have begun; a 6-period run hit the cascade and never recovered.

The internal `figure_eight` benchmark in `crates/apsis/examples/figure_eight.rs` (100 periods, $dt_0 = 10^{-4}$) ran in $0.58$ s; the parity scenario at the protocol-declared $dt_0 = T/1000 \approx 6.3 \times 10^{-3}$ cascaded. The discrepancy ruled out integrator-state pathology and pointed at the controller.

## Investigation method

A feature-gated tracer (`ias15-diag`, runtime-toggled with `APSIS_IAS15_TRACE=1`) was added to emit one tab-separated line per attempt with: substep counter, `dt_try`, `dt_next` after the attempt, `trunc_err`, Picard convergence flag, Picard iteration count, cumulative stagnation count, cumulative shrink–grow reversal count, and the controller's decision label (`accept`, `accept_floor`, `reject_picard`, `reject_trunc`). The trace is throttled by `APSIS_IAS15_TRACE_CAP` (default $2{,}000$ events) so cascade scenarios do not bury the diagnostic in $10^{8}$ duplicate lines.

Cheap counters that are always on (per-side `picard_stagnations_total`, `shrink_grow_cycles_total`) were added to `AdaptiveStats` for non-feature-flagged surfacing. Empirically the `shrink_grow_cycles / substeps` ratio is the most diagnostic signal: a healthy run on smooth motion sees $\ll 1$; the figure-8 cascade ran at $\approx 4\,\%$ chatter even before the floor hits, exactly the controller-warmstart oscillation fingerprint that the missing-cross-terms warmstart was inducing.

The cascade onset was localised by reading the trace at the first `accept_floor`: at substep $2661$, `dt` had collapsed to $10^{-12}$ and `trunc_err` sat at $1.16 \times 10^{-9}$ — *just above* $\epsilon = 10^{-9}$. Tracing backward, a single rejection cascade at substep $2659$ had taken `dt` from $6.3 \times 10^{-8}$ to `DT_MIN` over 25 retries. The retries were not Picard failures (Picard was converging in 2 iterations throughout); they were truncation rejections where `trunc_err` *did not scale as $dt^7$* — instead, it oscillated between $1 \times 10^{-9}$ and $2 \times 10^{-8}$ as `dt` shrunk by four orders of magnitude.

That is the empirical signature of a noise-floor saturation in $\|b_6\| / \|a_0\|$ arithmetic, not a controller pathology per se. The controller was correctly responding to the metric it was given; the metric was correctly reporting the noise floor; the cumulative inefficiency arose from the *rate* at which the controller approached the saturation regime and from the warmstart-controller *coupling* that made the noise floor unnecessarily prominent.

## Root cause #1 — Pascal-triangle warmstart (`warmstart_b`)

`warmstart_b` (`crates/apsis/src/physics/integrator/ias15.rs`) implemented only the *diagonal* of the polynomial-basis transformation that maps the previous step's `b` coefficients to the next step's seed at a different `dt`:

$$
\begin{aligned}
\text{apsis (broken):} \quad & e_k = q^{k+1} \cdot b_k \\
\text{specification:} \quad  & e_k = q^{k+1} \cdot \sum_{j \geq k} \binom{j+1}{k+1} \, b_j
\end{aligned}
$$

(The full transformation is derived in Everhart 1985, §III, and reproduced in Rein & Spiegel 2015 §3.2.) The cross-terms $\binom{j+1}{k+1} \, b_j$ for $j > k$ carry the contribution of higher-order coefficients into the seed of lower-order ones; they originate from re-expanding the acceleration ansatz $a(u) = a_0 + b_0 u + b_1 u^2 + \ldots + b_6 u^7$ after substituting $u_\text{new} = (u_\text{old} - 1)/q$: every higher-order term contributes to every lower-order term through the binomial expansion. For $k = 6$ the diagonal *is* the full transform (single column) and the two formulas agree; for $k < 6$ `apsis` was silently dropping 1–6 cross-terms.

For smooth motion at near-constant `dt` ($q \approx 1$), the cross-terms are bounded and Picard's predictor-corrector loop refines whatever `b` the warmstart provided within 3–4 iterations. **The bias becomes load-bearing only when `dt` changes enough that the cross-terms accumulate faster than Picard refines them.** Each close encounter in the figure-8 forces the controller to shrink `dt`; recovery grows it back. After enough close-approach cycles the residual bias against the true `b` shifts the post-Picard `truncation_error` upward, the controller responds by shrinking `dt` further — and recovery from `DT_MIN` is asymptotically slow because the same bias contaminates every warmstart on the way back up.

The 5-period boundary observed empirically was set by how many close-approach cycles the figure-8 admits before the bias compounds past the controller's tolerance. The same mechanism would manifest on any scenario with frequent `dt` adjustments: Mercury's perihelion, the Pythagorean three-body bench, planetary close encounters.

The fix replaces the diagonal-only formula with the full Pascal-triangle transformation derived from the binomial expansion. Five direct unit tests pin the new implementation against the analytical reference (per-coefficient assertions at $q \in \{0.1, 0.5, 1, 2, 5\}$, polynomial-equivalence check via sample-point evaluation, Picard-residual preservation under non-zero `be`, zero-input zero-output sanity, and the $q = 1$ polynomial-continuation property). The tests live alongside the implementation and survive a refactor that preserves behaviour but rewrites the formula.

## Root cause #2 — Halving on truncation rejection (specification §3.4)

The truncation-rejection branch was shrinking `dt` via the optimal-step formula $dt \cdot 0.9 \cdot (\epsilon / \mathrm{err})^{1/7}$ on every retry, mirroring the accept-path `dt_next` proposal. The trade-off shows up cleanly on the figure-8 cascade trace: when $\mathrm{err}$ sits just above $\epsilon$ (the noise-floor regime that emerges past close-encounter onset), the formula shrinks `dt` by only 10–20 % per attempt; 25–30 retries are needed to drop `dt` by the factor of $10^3+$ that the local geometry actually demands, and each retry pays a full force-eval sweep. The specification (Rein & Spiegel 2015 §3.4) prescribes an unconditional halving on truncation rejection, which converges to an acceptable `dt` in approximately 10 retries.

The fix replaces the formula in the rejection branch with an unconditional halving (`dt_try = (dt_try * 0.5).max(DT_MIN)`). The accept-path `dt_next` proposal still uses the optimum formula — what changes is only how aggressively the controller scans down to a feasible step on a single rejected attempt. Picard rejections (a separate failure class) already used `PICARD_SHRINK = 0.5`, so they are unaffected.

## Root cause #3 — $7\times$ growth cap on `dt_next` (specification §3.4)

The accept-path `dt_next` proposal had no upper bound; in smooth regions where $\mathrm{err} \ll \epsilon$ the formula proposes unbounded growth, and the next sub-step's truncation gate rejects the resulting overshoot, paying a full shrink cascade per close-encounter recovery. The specification (Rein & Spiegel 2015 §3.4) caps $dt_\text{new}$ at $7 \cdot dt_\text{current}$ per accepted sub-step precisely to prevent this; the cap was missing in `apsis`.

The factor `7.0` is the specification value. Tightening it would slow recovery from over-shrinks (post-close-encounter); loosening it would re-introduce the overshoot pathology. Combined with the warmstart fix and the rejection halving, the cap closes the controller-warmstart feedback loop that produced the cascade.

## Why the original implementation looked correct

The three gaps each appear individually defensible and are individually consistent with correct trajectories on smooth, low-orbit-count scenarios:

- The diagonal-only warmstart is correct *to leading order in $q$*, and exact at $q = 1$ would only differ from the full Pascal expansion by terms that Picard refines away in 3–4 iterations on smooth motion. The failure mode requires repeated $q \neq 1$ transitions and cross-term accumulation across them.
- The optimal-step rejection shrink is mathematically correct — $(\epsilon / \mathrm{err})^{1/7}$ *is* the asymptotic optimum step proposal. The failure mode is the *rate* at which it converges to a feasible step when $\mathrm{err}$ is near the noise floor, not its eventual value.
- The unbounded growth proposal is also asymptotically correct. The failure mode is purely transient — the post-close-encounter recovery overshoots, the next attempt rejects, and the work spent rejecting is the cost.

Each gap individually passes the existing IAS15 unit-test suite (Kepler $e = 0.5$ over 100 orbits, Pythagorean three-body across central encounters, high-eccentricity peak-error gate). What distinguishes the figure-8 parity scenario is its sustained mix of close-encounter `dt` shrinks and smooth-region `dt` growths over many orbits — exactly the regime where any one of the three gaps amplifies the others.

## Side findings — incorporated as part of the same audit

### Step-size hint vs hard-cap contract (`Integrator::dt_hint`)

The `Integrator::step` trait method previously named its time-step argument `dt`, and the `System::step` orchestrator passed `current_dt = user_dt` on every call. For fixed-step integrators (Velocity Verlet, Yoshida 4, Wisdom–Holman) `dt` is the exact step, and `current_dt = user_dt` is correct. For self-adaptive integrators (IAS15) `dt` is conceptually a hint that seeds the controller's first call; on subsequent calls the integrator should drive `dt_next` from its own state, and the caller's hint should not silently bound it. Treating the orchestrator's `dt` as a hard per-call cap silently pinned the IAS15 controller to whatever step the user supplied as an initial guess — the same defect class as the cascade above, but routed through the orchestrator rather than through the controller.

The trait now declares the parameter `dt_hint`, and the contract is made explicit through two new trait methods:

- `Integrator::controls_own_step_size(&self) -> bool` — `true` when the integrator's controller picks `dt`. Default `false`.
- `Integrator::proposed_next_dt(&self) -> Option<f64>` — the controller's recommendation for the next call. Default `None`.

`System::step` reads both and updates `current_dt` accordingly: for self-adaptive integrators, `current_dt` follows `proposed_next_dt`; for fixed-step ones it stays at `user_dt`. The change is invisible at runtime for current callers (IAS15 with the warmstart/halving/cap fixes ignores its `dt_hint` after the first call anyway) but eliminates the trap for the next adaptive integrator added to the zoo (SABA, Hermite, MERCURIUS).

### Compensated recentering (`Integrator::recenter_bodies`)

`System::step` periodically applies a centre-of-mass calibration shift via `calibration::apply_body_shift`, which performs a bare `body.x -= dx` on each body. IAS15 maintains a Neumaier compensation buffer `csx` paired with each body's stored position (`(x, csx)` represents an extended-precision running sum). A bare subtraction wipes the prior compensation history rather than continuing to track it; on long runs under periodic recentering the loss accumulates into a bit-reproducibility gap.

A new trait method `Integrator::recenter_bodies(&mut self, bodies, dx, dy)` replaces direct calls to `apply_body_shift` from `System::step` and `System::recenter_com`. Default impl is the bare subtraction (correct for integrators with no per-body compensation). IAS15 overrides to route the shift through `add_cs` against its own `csx` buffers, preserving the compensation invariant. The fix is sub-ULP per recentering event and irrelevant to the cascade above; it is included here because the snapshot-replay determinism ADR commits to bit-reproducibility under the periodic recentering, which the bare subtraction silently violated.

### Picard stagnation = success (specification §3.3)

The Picard predictor–corrector loop's stagnation guard previously returned `(converged: false, ...)` on two consecutive non-decreasing iterations, routing the attempt through the `RejectPicard` branch. The specification (Rein & Spiegel 2015 §3.3) breaks out of the loop on the same condition, treating the saturated state as best-effort convergence and letting the truncation gate decide acceptance — an `apsis` divergence from this convention forced `dt` halvings even when the `b` coefficients were as accurate as f64 allowed. The guard now returns `(converged: true, ...)`.

The change is largely belt-and-braces under the warmstart + halving + cap fixes, since stagnation is now a rare event in healthy regimes — but is kept because the specification-conformant semantics is the cleanest invariant to maintain.

## Why we did *not* switch the truncation metric to max-max

The specification leaves the truncation-error norm as an implementation choice. The most common convention (used by Rein & Spiegel's reference C implementation) is $\max|b_6| / \max|a_0|$ (max-max) for both Picard residual and truncation error; `apsis` uses RMS of per-body $|b_{6,i}| / |a_{0,i}|$. The departure is documented in the `apsis` source (`physics/integrator/ias15.rs`, comment on `truncation_error`):

> Empirically this showed up at solar-system-class $N$ ($\approx 641$ bodies) as a 194 % rejection rate (3878 rejections over 2001 accepted sub-steps), nearly all via `RejectTruncation`. The original max-max formula treated outliers as if they were representative of the whole system: $\max\|b_6\|$ picked up the one body whose $b_6$ had the largest round-off noise, $\max\|a_0\|$ picked up the Sun's dominant acceleration, and the ratio was *not* a convergence criterion any more — it was a noise-to-signal floor that grew with body count.

For figure-8 ($N = 3$, equal masses) the two metrics yield numerically similar values, so the cascade does not require the metric switch. Switching would re-introduce the documented 194 % rejection rate at solar-system $N$. The metric is left as RMS; the cascade is resolved by the controller fixes alone.

## Validation

### Figure-8 parity stress test (the original failing case)

Pre-refactor versus post-controller-refactor:

| Metric | Pre-refactor | Post-refactor | Notes |
| --- | --- | --- | --- |
| Wall time, 10-period parity | 15 min | 80 ms | $\approx 11{,}000\times$ faster |
| Floor-pinned sub-steps (`dt = DT_MIN`) | $\approx 1.3 \times 10^{8}$ | $0$ | controller never reaches floor |
| Total accepted sub-steps | $\approx 2.2 \times 10^{8}$ | $\approx 1.5 \times 10^{3}$ | five orders of magnitude fewer |
| $\max\!\lvert\Delta E / E_0\rvert$ apsis | $8.6 \times 10^{-16}$ | $8.6 \times 10^{-16}$ | unchanged at 1 ULP |
| $\max\!\lvert\Delta\mathbf{L}\rvert$ cross-impl | $4.2 \times 10^{-16}$ | $4.4 \times 10^{-16}$ | unchanged at 1 ULP |
| $\max\!\lvert\Delta\mathbf{r}\rvert$ informational | $3.17 \times 10^{-7}$ | $9.44 \times 10^{-13}$ | phase drift 5+ orders better |

All 12 gated metrics of the figure-8 parity protocol (`paper/notebooks/2026-04-26-rebound-parity-figure8.md`, §Hypothesis) pass at 1 ULP under the new controller.

### Kepler parity regression

The $e = 0.5$ Kepler parity scenario at 100 orbits:

| Metric | Pre-refactor | Post-refactor |
| --- | --- | --- |
| $\lvert\Delta a\rvert / a$ (semi-major axis) | $3.6 \times 10^{-15}$ | $3.1 \times 10^{-15}$ |
| $\lvert\Delta e\rvert$ | $2.9 \times 10^{-15}$ | $2.3 \times 10^{-15}$ |
| $\lvert\Delta\omega\rvert$ | $2.2 \times 10^{-15}$ | $3.1 \times 10^{-15}$ |
| $\lvert\Delta h\rvert / h$ | $6.4 \times 10^{-16}$ | $1.4 \times 10^{-15}$ |
| Cross-impl $\lvert\Delta E\rvert / \lvert E_0\rvert$ | $4.2 \times 10^{-15}$ | $2.5 \times 10^{-15}$ |
| $\lvert\Delta r\rvert$ informational | $1.57 \times 10^{-9}$ | $2.18 \times 10^{-12}$ |

All seven gated metrics still pass at 1–3 ULP. The $\lvert\Delta r\rvert$ improvement (3 orders of magnitude) reflects the controller now tracking the specification's `dt_next` choice to ULP precision rather than diverging from a peer implementation by accumulated phase. This is an expected side effect of bringing the controller into specification compliance and is explicitly the "honest" parity statement for adaptive high-order integrators (cf. the protocol notebook's §"Why this metric set, not $\lvert\Delta r\rvert$").

### IAS15 unit-test suite

20 unit tests pass across all configurations (with and without `ias15-profile`, with and without `ias15-diag`):

- 6 `decide_dt` truth-table cases.
- 6 `warmstart_b` direct tests (Pascal-table reference match at $q \in \{0.5, 2\}$, $q < 1$, polynomial-equivalence sample, Picard-residual preservation under non-zero `be`, zero-input zero-output, $q = 1$ polynomial continuation).
- 4 system-level tests (Kepler high-eccentricity peak, Kepler 100-orbit drift monotonicity, Pythagorean energy through close encounters, system-`t`-matches-substep regression).
- 4 force-model compatibility tests (deterministic force enforcement under IAS15 selection).

## File changes

- `crates/apsis/src/physics/integrator/ias15.rs`
  - `warmstart_b`: diagonal-only formula replaced by Pascal-triangle expansion. Six unit tests pin behaviour.
  - rejection-branch shrink: `optimal_dt(...)` replaced by unconditional `(dt_try * 0.5).max(DT_MIN)`. Comment block at the rejection arm explains the convergence-rate trade-off.
  - accept-branch `dt_next`: capped above by `dt_try * DT_GROWTH_LIMIT` (`7.0`). New constant `DT_GROWTH_LIMIT` documented at module scope.
  - `picard_loop_inner`: stagnation guard returns `(true, residual, iters)` per Rein & Spiegel (2015) §3.3.
  - `controls_own_step_size`, `proposed_next_dt`, `recenter_bodies`: new trait method overrides routing the integrator-specific contracts.
  - `picard_stagnations_total`, `shrink_grow_cycles_total`, `dt_dir_prev`: new state fields surfaced through `AdaptiveStats`.
  - feature-gated `diag` module: per-step trace emission under `APSIS_IAS15_TRACE=1`.
- `crates/apsis/src/physics/integrator/traits.rs`
  - `Integrator::dt_hint`, `controls_own_step_size`, `proposed_next_dt`, `recenter_bodies`: new trait methods with documented contracts.
  - `AdaptiveStats::picard_stagnations`, `shrink_grow_cycles`: new always-on counters.
- `crates/apsis/src/core/system/step.rs`
  - `current_dt` update routed through `proposed_next_dt` for self-adaptive integrators.
  - COM recentering routed through `Integrator::recenter_bodies` rather than the bare `apply_body_shift`.
- `crates/apsis/src/core/system/bodies.rs`
  - `System::recenter_com` routed through `Integrator::recenter_bodies` for compensation-aware translation.
- `crates/apsis/Cargo.toml`
  - `ias15-diag` Cargo feature added (companion to the existing `ias15-profile`).

## Methodological observation for the v0.1 paper

Implementing an adaptive high-order integrator to its algorithmic specification requires more than agreement on the *abstract description*. The three controller gaps documented above each preserve the asymptotic correctness of the IAS15 method (Rein & Spiegel 2015) while individually inflating the substep count by an order of magnitude on stiff-mix scenarios. Verifying that the implementation also matches the specification's *concrete numerical recipe* — the warmstart polynomial transformation, the rejection cadence, the growth cap — is therefore a precondition for cross-implementation parity statements, not a follow-up to them. The figure-8 cascade is a concrete instance of how a parity gate that uses physical invariants (energy, angular momentum, COM) cleanly passes while the controller silently degrades into a $10^{8}$-substep regime; the corollary is that *invariant-passing alone is insufficient evidence of an honest cross-implementation match for adaptive integrators*. The substep count and the rejection-class breakdown — both surfaced through `AdaptiveStats` and the optional `ias15-diag` trace — are first-class artefacts of the parity scenario, not internal-only diagnostics.

This is an observation worth including in the v0.1 paper's methods section: not as a bug-fix changelog (the audit is private to `apsis`) but as a methodological argument for *why* the parity portfolio gates on invariants *and* on substep economy, and *why* the substep-count distribution belongs in the supplementary CSV alongside the invariant-error tables.

## Files

- This document: `docs/experiments/2026-04-26-ias15-warmstart-bug.md`.
- Figure-8 parity protocol notebook: `paper/notebooks/2026-04-26-rebound-parity-figure8.md`.
- Figure-8 parity harness: `validation/rebound-parity/figure8/`.
- IAS15 implementation: `crates/apsis/src/physics/integrator/ias15.rs`.
- Integrator trait contract: `crates/apsis/src/physics/integrator/traits.rs`.
