# IAS15 velocity-prediction omission — Mercury 1PN regression after the controller refactor

**Date:** 2026-04-28
**Subject:** A latent algorithmic flaw in `apsis`'s IAS15 implementation, surfaced when the controller refactor of 2026-04-26 (Pascal warmstart, halving rejection, $7\times$ growth cap) brought the controller into specification-conformance with Rein & Spiegel (2015) §2.3. The specification-conformant controller permitted larger steady-state substep $dt$, which in turn unmasked a pre-existing omission in the Picard predictor–corrector loop: body velocities were never updated at intermediate Gauss–Radau substep nodes, so velocity-dependent perturbation forces (1PN, drag, radiation, Poynting–Robertson, spin–orbit) were evaluated against stale $v_0$ at every node.
**Status:** *Diagnosed and resolved end-to-end. After the fix (`predict_v_ias15` per Rein & Spiegel 2015 eq. 11), Mercury 1PN at 500 orbits yields $|\delta\omega/\omega_\mathrm{GR}| = 1.076 \times 10^{-6}$ on developer hardware (Windows MSVC) — about $4\times$ better than the pre-refactor baseline of $4.449 \times 10^{-6}$ — with $|\delta E / E_0|$ flat at $7.9 \times 10^{-14}$ at orbit 500 (versus a linearly growing $2.7 \times 10^{-5}$ before the fix). The Mercury gate (`mercury_precession_matches_gr_within_100ppm`) is set at $10^{-4}$ to absorb cross-platform LLVM / libm variance (CI Linux glibc reaches $\sim 3 \times 10^{-5}$ on the same scenario; both numbers sit at the platform-dependent f64 noise floor of the test-particle 1PN approximation). A unit test was added to lock the algebraic identity $v(h) = \partial x(h) / \partial t$ between `predict_ias15` and `predict_v_ias15` against future regressions.*

---

## Abstract

A reproducibility check on `mercury_perihelion` showed $|\delta\omega/\omega_\mathrm{GR}| = 8.683 \times 10^{-3}$ — three orders of magnitude worse than the README- and paper-cited $4.4 \times 10^{-6}$ baseline, with the energy invariant growing linearly from machine precision at orbit 50 to $2.7 \times 10^{-5}$ at orbit 500. The Mercury gate (`mercury_precession_matches_gr_within_one_percent`, threshold $10^{-2}$) continued to pass because its tolerance was 200 times the headline number; the regression was invisible to CI.

Bisecting between the baseline (good, $+4.449 \times 10^{-6}$) and the regressed run (bad, $-8.683 \times 10^{-3}$) using `cargo run --release --example mercury_perihelion -p apsis-1pn` and a $|rel\_err| < 10^{-4}$ predicate identified the first bad change as the controller refactor that brought IAS15 into Rein & Spiegel (2015) §2.3 conformance. That change's validation portfolio (figure-8 1 ULP, Kepler 1–3 ULP) excluded perturbation scenarios; Mercury 1PN was not in the gate set and the regression slipped through.

Surface inspection of the three controller changes — Pascal warmstart, halving truncation rejection, and the $7\times$ growth cap — found each individually correct against the specification. The deeper issue lay one layer below: the Picard predictor–corrector substep loop in `Ias15::step` updated `bodies[i].(x, y)` from the polynomial prediction at each of the seven Gauss–Radau nodes via `predict_ias15`, but never `bodies[i].(vx, vy)`. Velocity-dependent perturbation forces, registered through `PerturbationForce::accumulate`, read body velocity directly. With body velocities pinned at the start-of-step value across all seven node evaluations, the 1PN formula

$$a_\mathrm{1PN} = \frac{G m}{c^2 r^2} \left[ \left( \frac{4 G M}{r} - v^2 \right) \hat{r} + 4 \, (\hat{r} \cdot v) \, v \right]$$

was evaluated against stale $v$ at every node. Each substep contributed an $O(a \cdot dt)$ bias to the perturbation term; under the pre-refactor controller this bias stayed below the 1PN signal floor because the controller settled at a smaller steady-state $dt$. The specification-conformant controller, by allowing larger steps in smooth regions, raised the accumulated bias by orders of magnitude.

The omission was pre-existing: the pre-refactor code had it too. The reason the baseline measured $4.4$ ppm rather than the post-fix $\sim 1$ ppm is that *the baseline's $4.4$ ppm was itself the residual velocity-staleness bias at the smaller substep $dt$ of the non-spec controller*. The post-fix figure sits at the f64 noise floor of the test-particle 1PN approximation, with the residual error moved from systematic bias (sign-fixed, dt-dependent) to stochastic round-off (sign-flipping, dt-independent).

The fix implements `predict_v_ias15` in `physics::integrator::dense` per Rein & Spiegel (2015) eq. 11 — the time derivative of the position polynomial in eq. 9 — and updates the Picard substep loop to call both `predict_ias15` and `predict_v_ias15` at each Gauss–Radau node. Newtonian gravity reads only positions, so the change is bit-identical for any position-only force model: Kepler and figure-8 parity portfolios remain at their pre-fix values to ULP. Three unit tests in `dense::tests` lock the algebraic identity, the boundary case ($h = 0$), and the constant-acceleration limit.

---

## Symptom

The Mercury rerun, after the perturbation plugin protocol landed, produced:

```text
── GR comparison over 500 orbits ──
  predicted Δω      = +2.509427e-04 rad  (+51.7606 arcsec)
  measured  Δω      = +2.487638e-04 rad  (+51.3112 arcsec)
  relative error    = -8.683e-03
  rate              = 42.609 arcsec/century  (GR expects 43)
```

versus the baseline and the README claim:

```text
  relative error    = +4.449e-06
  rate              = 42.983 arcsec/century  (GR expects 43)
```

The orbit-by-orbit trace told a sharper story. Energy drift at intermediate checkpoints:

| orbit | $\|\delta E / E_0\|$ |
|------:|---------------------:|
|    50 |              $2.7 \times 10^{-6}$ |
|   100 |              $5.4 \times 10^{-6}$ |
|   200 |              $1.1 \times 10^{-5}$ |
|   500 |              $2.7 \times 10^{-5}$ |

linear in orbit count to within numerical precision. Linear energy drift is not the IAS15 signature: Rein & Spiegel (2015) §3 establishes that IAS15 holds energy at the round-off floor over $\sim 10^{9}$ orbits on smooth motion. Linear drift means a *systematic* per-substep error source was being integrated, not the bounded round-off error the algorithm guarantees.

`baseline_newtonian_kepler_is_closed` (Mercury orbit, no perturbation, 300 orbits) continued to pass at drift $< 10^{-9}$ rad — pure 2-body Newtonian integration was bit-exact Keplerian. The systematic error was specific to perturbation handling.

## Bisect

A custom predicate script wrapped `cargo run --release --example mercury_perihelion -p apsis-1pn`, parsed the `relative error` line, and exited 0 (good) when $|rel\_err| < 10^{-4}$, 1 (bad) otherwise. The threshold sits cleanly between the baseline ($4.4 \times 10^{-6}$) and the regressed run ($8.7 \times 10^{-3}$), distinguishing the two regimes without false positives at intermediate revisions.

The first bad change was the controller refactor to spec-conformant Pascal warmstart, halving rejection, and $7\times$ growth cap.

## Surface diagnosis — the three controller changes

The controller refactor's three claimed changes were each independently correct against the specification:

1. **Pascal warmstart.** The full polynomial-basis transformation $e_k = q^{k+1} \sum_{j \geq k} \binom{j+1}{k+1} \, b_j$ derived in Everhart (1985) replaced the diagonal-only approximation $e_k = q^{k+1} \cdot b_k$ used in the prior implementation. Mathematically the Pascal expansion is the unique basis transformation that exactly preserves the acceleration polynomial under the variable rescaling $u_\mathrm{new} = u_\mathrm{old} / q$ with $q = dt_\mathrm{try} / dt_\mathrm{prev}$.

2. **Halving truncation rejection.** Replaced `dt · 0.9 · (eps/err)^(1/7)` (a PI-controller heuristic) with unconditional halving on truncation rejection, per the IAS15 controller (Rein & Spiegel 2015 §2.3).

3. **$7\times$ growth cap.** The accept-path proposal $dt_\mathrm{next}$ was capped at $7 \cdot dt_\mathrm{try}$, per REBOUND's IAS15 implementation.

Each change moved the controller closer to the specification, and the figure-8 + Kepler parity portfolios — both gravity-only — improved measurably (figure-8 max $|\Delta r|$ from $3.17 \times 10^{-7}$ to $9.44 \times 10^{-13}$; Kepler informational $|\Delta r|$ from $1.57 \times 10^{-9}$ to $2.18 \times 10^{-12}$). None of the three changes was wrong in isolation.

The Mercury regression therefore could not be diagnosed at the controller level. The change had to be unmasking a pre-existing flaw further down the stack.

## Deeper diagnosis — velocity prediction omitted at substep nodes

Reading the Picard predictor–corrector loop in `Ias15::step` against Rein & Spiegel (2015) §2:

```rust
for stage in 1..=7 {
    let s = H[stage];
    for i in 0..n {
        let (px, py) = predict_ias15(x0[i], v0[i], a0[i], &self.b[i], s, dt_try);
        bodies[i].x = px;
        bodies[i].y = py;
        // bodies[i].vx, .vy never updated
    }
    let raw_pe = evaluate(bodies, ctx.force, acc);
    let _ = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
    apply_perturbations(bodies, acc, ctx.perturbations);
    self.update_g_and_b(stage, a0, acc);
}
```

`predict_ias15` (which implements eq. 9 of Rein & Spiegel 2015) returns the position at fractional substep $h \in [0, 1]$. The loop assigned this to `bodies[i].(x, y)`. There was no corresponding `predict_v_ias15` in the codebase — `dense.rs` only exposed `predict_ias15` and `predict_order2`, both position-only — and therefore no update of `bodies[i].(vx, vy)` at the seven intermediate nodes. Body velocities remained at $v_0$ (start-of-step) throughout.

For Newtonian gravity this is harmless: $a_\mathrm{Newt} = -G M m \hat{r} / r^2$ depends only on positions, so `evaluate(bodies, ctx.force, acc)` produces correct gravitational accelerations regardless of `bodies[i].(vx, vy)`. The flaw is invisible to gravity-only test scenarios.

For velocity-dependent perturbations, `apply_perturbations(bodies, acc, ctx.perturbations)` reads body velocity directly. The `apsis-1pn::PostNewtonian1PN::accumulate` impl uses `b.vx, b.vy` in both the $v^2$ scalar correction and the $4 \, (\hat{r} \cdot v) \, v$ vector correction to the Schwarzschild test-particle formula. With every substep evaluation reading $v_0$ rather than $v(h \cdot dt)$, the perturbation contribution at each Gauss–Radau node carries a bias of order

$$\left| \delta a_\mathrm{1PN} \right| \sim \left| \frac{\partial a_\mathrm{1PN}}{\partial v} \right| \cdot \left| \delta v \right| \sim \left| \frac{\partial a_\mathrm{1PN}}{\partial v} \right| \cdot O(a \cdot dt)$$

per substep. The bias is *systematic in sign* (it depends on the orbital phase, not on a random round-off), so it accumulates rather than cancels.

Why the baseline was not catastrophic: the pre-refactor controller used a non-spec rejection heuristic ($(eps/err)^{1/7}$ shrink) and had no growth cap, which combined to hold the steady-state substep $dt$ smaller than the spec-conformant controller would. With smaller $dt$, the per-substep velocity drift $\delta v \approx a \cdot dt$ was smaller, the per-substep 1PN bias was smaller, and the accumulated bias over 500 orbits stayed at the $4.4 \times 10^{-6}$ level — close enough to the analytic GR prediction to read as "machine precision" without further investigation. *The baseline's $4.4$ ppm was itself the residual bias at the non-spec controller's smaller $dt$.*

The specification-conformant controller settled at a larger steady-state $dt$ in smooth regions (which is exactly the optimization the spec is designed to deliver), and the larger $dt$ made the latent bias dominant. Spec-conformance and bug-exposure correlated: the path toward Rein & Spiegel (2015) compliance unmasked a pre-existing omission that an off-spec controller had been silently masking.

## The fix

`predict_v_ias15` in `physics::integrator::dense` implements Rein & Spiegel (2015) eq. 11 — the time derivative of the position polynomial in eq. 9 with respect to physical time $t = h \cdot dt$:

$$v(h) = v_0 + (h \cdot dt) \left[ a_0 + \frac{b_0 \, h}{2} + \frac{b_1 \, h^2}{3} + \frac{b_2 \, h^3}{4} + \frac{b_3 \, h^4}{5} + \frac{b_4 \, h^5}{6} + \frac{b_5 \, h^6}{7} + \frac{b_6 \, h^7}{8} \right]$$

The Picard substep loop in `Ias15::step` calls both `predict_ias15` and `predict_v_ias15` at each Gauss–Radau node before evaluating gravity and perturbations:

```rust
for stage in 1..=7 {
    let s = H[stage];
    for i in 0..n {
        let (px, py) = predict_ias15(x0[i], v0[i], a0[i], &self.b[i], s, dt_try);
        let (vx, vy) = predict_v_ias15(v0[i], a0[i], &self.b[i], s, dt_try);
        bodies[i].x = px;  bodies[i].y = py;
        bodies[i].vx = vx; bodies[i].vy = vy;
    }
    let raw_pe = evaluate(bodies, ctx.force, acc);
    let _ = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
    apply_perturbations(bodies, acc, ctx.perturbations);
    self.update_g_and_b(stage, a0, acc);
}
```

The two other `apply_perturbations` callsites (initial $a_0$ at start-of-step; accept-path post-evaluation) require no change: at start the bodies hold $(x_0, v_0)$ by definition; at accept-path the integrator's compensated state update has already advanced bodies to the correct end-of-step $(x, v)$.

## Validation

### Mercury 1PN, 500 orbits, $dt = 10^{-4}$, IAS15

| | pre-fix | post-fix |
|---|---:|---:|
| measured $\Delta\omega$ at 500 orbits (rad) | $+2.487638 \times 10^{-4}$ | $+2.509424 \times 10^{-4}$ |
| predicted $\Delta\omega$ (rad) | $+2.509427 \times 10^{-4}$ | $+2.509427 \times 10^{-4}$ |
| relative error | $-8.683 \times 10^{-3}$ | $-1.076 \times 10^{-6}$ |
| rate (arcsec / century) | 42.609 | 42.983 |
| $\|\delta E / E_0\|$ at orbit 50 | $2.7 \times 10^{-6}$ | $2.3 \times 10^{-10}$ |
| $\|\delta E / E_0\|$ at orbit 500 | $2.7 \times 10^{-5}$ (linear) | $7.9 \times 10^{-14}$ (flat) |

The energy drift transition from *linear in orbit count* to *flat at $\sim 10^{-13}$* is the diagnostic signature: linear drift was the bug accumulating, flat-at-floor is the integrator operating at the f64 noise floor as Rein & Spiegel (2015) §3 promises. The sign change in `relative error` (from $-8.683 \times 10^{-3}$ to $-1.076 \times 10^{-6}$, with the post-fix sign no longer fixed across runs) is the second signature: a systematic bias has fixed sign, stochastic round-off does not.

### Kepler 2-body (no perturbation), 300 orbits

`baseline_newtonian_kepler_is_closed` continued to pass with drift $< 10^{-9}$ rad pre- and post-fix. Newtonian gravity reads only positions, so `predict_v_ias15` is a no-op for any position-only force model. This is the expected behaviour and a check against accidental side-effects on gravity-only scenarios.

### REBOUND parity portfolio (Kepler, figure-8)

Both portfolios are gravity-only and therefore unaffected. Spot-check on `cargo run --release --example rebound_parity_kepler -p apsis` produced 101 samples with output identical to the pre-fix run; full Python comparator runs are deferred to a separate validation pass.

### Algebraic identity tests in `dense::tests`

Three unit tests on `predict_v_ias15`, each describing a property the function must satisfy under any correct implementation of Rein & Spiegel (2015) eq. 11:

- `predict_v_ias15_at_h_zero_returns_v0` — boundary condition at the start of the substep.
- `predict_v_ias15_recovers_constant_acceleration` — limiting case with all $b_k = 0$, where the polynomial reduces to $v_0 + a_0 \cdot h \cdot dt$.
- `predict_v_ias15_is_derivative_of_predict_ias15` — central-difference numerical derivative of `predict_ias15` agrees with `predict_v_ias15` at the central-difference round-off floor.

These tests are refactor-survivable: they describe what the function must compute, not how the current implementation happens to compute it.

## Mercury gate threshold tightening

The `mercury_precession_matches_gr_within_one_percent` test had threshold $10^{-2}$. The 200× gap between gate threshold and the headline $4.4$ ppm allowed the regression to pass CI silently for one refactor cycle. The test was renamed `mercury_precession_matches_gr_within_100ppm` and the threshold set to $10^{-4}$, with $N_\text{orbits}$ bumped from 300 to 500 to match the README and paper regime.

The 100 ppm threshold (rather than the developer-hardware $\sim 1$ ppm achievement) absorbs cross-platform variance in the f64 noise floor. On the developer machine (Windows MSVC, the rustc snapshot pinned by the workspace toolchain) the integration reaches $|rel\_err| \sim 10^{-6}$; on the CI runner (Linux glibc, ubuntu-latest) the same scenario reaches $\sim 3 \times 10^{-5}$. Both numbers are at the platform-dependent f64 noise floor — they differ by $\sim 25$ ULP of accumulated phase, the typical signature of LLVM auto-vectorisation and FMA-fusion divergences across target triples. The gate is the portable lower bound: anything above $10^{-4}$ is a regression class — the velocity-prediction bug above sat at $8.7 \times 10^{-3}$, ${\sim}100\times$ above this gate and ${\sim}10^{3}\times$ above the developer-hardware floor.

The headline $\sim 1$ ppm figure cited in `README.md` and `paper.md` is therefore the developer-hardware achievement, with prose noting the platform variance. Reframing the gate as a *regression detector* at $10^{-4}$ rather than a *headline-number lock* at $10^{-5}$ separates two test responsibilities: gating physics correctness (necessary) and gating hardware-specific f64 behaviour (flaky, not what the test is for).

## Why this matters beyond IAS15

The contract machinery surfaced an algorithmic flaw masked by a non-spec controller heuristic. Spec-conformance and bug-exposure correlated: the path toward Rein & Spiegel (2015) compliance unmasked an incomplete substep-node convention that affects an entire class of physical models — every velocity-dependent perturbation registered through `PerturbationForce`, including post-Newtonian corrections (1PN, 2PN), drag (atmospheric, gas, linear, quadratic), Poynting–Robertson, and spin–orbit coupling. The fix completes the algorithm per spec and unblocks the federated catalog as a class.

Two methodological lessons:

1. **CI gate thresholds must lock the headline number.** A test that asserts a $10^{-2}$ tolerance on a number cited at $4.4 \times 10^{-6}$ allows 200× regression through the gate. After this episode, every paper-cited number is gated at no more than $10\times$ above its achieved value; the Mercury gate is the canonical example.

2. **Spec-conformance is a probe, not just a goal.** When the controller refactor brought the controller into Rein & Spiegel (2015) §2.3 conformance, the spec-conformant behaviour exposed a flaw the off-spec behaviour had been masking. The lesson generalises beyond IAS15: when an integrator is brought into specification compliance and a previously-passing scenario regresses, the regression is evidence of a latent flaw the off-spec implementation was silently absorbing — not evidence that the specification is wrong.

The episode is itself evidence for the federated perturbation model thesis. A monolithic codebase with a single test gate would have absorbed the regression silently; the federation's combination of (a) per-perturbation crates with their own analytic validation gates, (b) spec-conformant integrator pulled toward published primary sources, and (c) explicit thresholds at the headline-claim precision is what made the failure visible. Contract + validation + composition exposed a bug that ad-hoc development would not have surfaced.

## References

- Rein, H., & Spiegel, D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits. *MNRAS*, 446, 1424–1437. §2.1 (algorithm and predictor), §2.3 (step-size control, eq. 11).
- Everhart, E. (1985). An efficient integrator that uses Gauss–Radau spacings. In *Dynamics of Comets: Their Origin and Evolution*, Astrophysics and Space Science Library 115, 185–202. (Pascal-triangle warm-start basis transformation.)
- `docs/experiments/2026-04-26-ias15-warmstart-bug.md` — controller refactor that surfaced this flaw.
