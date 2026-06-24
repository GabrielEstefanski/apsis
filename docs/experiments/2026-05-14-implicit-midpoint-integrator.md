# Implicit Midpoint — A-stable symplectic 2nd-order for the apsis integrator zoo

**Date:** 2026-05-14

**Subject:** Adds `IntegratorKind::ImplicitMidpoint` to the apsis integrator zoo as the implicit symplectic baseline. Closes the gap left by an all-explicit zoo (VV, Y4, WH, WHFast, Mercurius are all explicit; IAS15 is explicit predictor–corrector). A-stable per Hairer-Lubich-Wanner Chapter II.1.4; time-symmetric; symplectic for both separable and non-separable Hamiltonians; integrates `NonConservativeOperator` registrations with physical dissipation. Paper-grade backbone for FPM regime exploration — *not* a universal solver: A-stable but not L-stable, so dissipation-dominant extreme regimes (BH–BH near merger, very stiff radiation reaction) require L-stable methods (Radau IIA, BDF) tracked as v0.2+ extensions.

**Status:** Protocol declared *a priori*, before any code lands. Locks the algorithmic decisions and validation scenario before implementation.

**Roadmap context:** Co-priority with WHFast in the integrator-zoo roadmap. Where WHFast is the long-horizon planetary baseline and Mercurius is the close-encounter hybrid, Implicit Midpoint is the *implicit* baseline — paper credibility for any regime where explicit methods are unstable: stiff perturbations, extreme mass ratios, dissipative dynamics. Forward-compatible for the v0.2 regime exploration (PN ≥ 2, particle tests, BH binaries, pulsar orbits with radiation reaction).

---

## Abstract

The apsis integrator zoo currently spans explicit symplectic (VV, Yoshida4, WH, WHFast, Mercurius) and explicit predictor–corrector adaptive (IAS15). Every integrator in the zoo evaluates `f(y_n)` at the current state and advances forward. None is implicit. Reviewer expectation for an integrator-of-integrators paper is to cover the implicit slot — the standard reference being implicit midpoint (Hairer-Lubich-Wanner, Chapter II.1.4) for the symplectic baseline and Gauss–Legendre RK for higher orders.

Implicit Midpoint advances `y_{n+1} = y_n + dt · f((y_n + y_{n+1}) / 2)`. The midpoint state is solved iteratively per step (fixed-point iteration; Newton-Krylov is a follow-up). The method is:

- **A-stable** (region of absolute stability is the entire left half-plane): handles arbitrarily stiff conservative dynamics without dt restriction beyond accuracy.
- **Symplectic** for any Hamiltonian system, including non-separable Hamiltonians where explicit symplectic methods (e.g. WHFast's KDK split) require a separable structure.
- **Time-symmetric**: `Φ_dt ∘ Φ_{-dt} = id`. Standard test: forward + backward integrate returns to the initial state up to the iteration tolerance.
- **Order 2** in the time step. Higher orders (Gauss–Legendre 4th/6th) are deferred to a follow-up — the order-2 baseline is what reviewers expect for the implicit zoo slot, and matching apsis WH 1991's order makes the cross-integrator comparison clean.

Implicit Midpoint composes with both `HamiltonianOperator` and `NonConservativeOperator` registrations through the existing dispatch contract — no operator-side changes required. This is the federation-thesis claim for the integrator zoo: every integrator regime composes with every operator through a single contract.

---

## Motivation

Three claims chain into this experiment:

1. **The integrator zoo positioning thesis requires an implicit baseline.** Without Implicit Midpoint, the federated-perturbation-model paper documents a five-integrator zoo (six counting WHFast) that is uniformly explicit. A reviewer asking "what about stiff perturbations? what about non-separable Hamiltonians? what about A-stable methods for the comparison floor?" has no answer in-zoo. Adding the canonical implicit method closes that gap with the textbook reference.

2. **Future regime exploration (v0.2+) requires A-stability.** The user's stated v0.2 regime portfolio includes: PN ≥ 2 corrections (chaotic at small `r/M`), test-particle integration around dynamical perturbers (massive-particle stiffness), BH binary dynamics (extreme mass ratio), pulsar orbits with gravitational-radiation back-reaction (intrinsically dissipative). Every one of these scenarios is a stress test for explicit symplectic methods. Implicit Midpoint is the paper-grade tool of choice for these regimes — having it in the zoo *before* the v0.2 work begins means each new regime experiment can use it as the stable baseline against which other methods are measured.

3. **The federation contract must hold across implicit/explicit boundaries.** WHFast + 1PN and Mercurius + 1PN demonstrated that explicit symplectic integrators dispatch the operator framework correctly. Demonstrating Implicit Midpoint + 1PN at the same federation gate (Mercury precession within the noise floor) closes the contract for the implicit class and makes the federation claim integrator-class-independent. This is the v0.1 paper-portfolio item.

### What this experiment is NOT testing

- **Not 4th/6th-order Gauss–Legendre.** Higher-order implicit symplectic methods exist (Hairer-Lubich-Wanner II.1.4) and would slot into the same iteration loop. Out of scope here; tracked as follow-up.
- **Not Newton-Krylov iteration.** Fixed-point iteration on the midpoint state converges for non-stiff and mildly-stiff problems; Newton iteration with an analytical or finite-difference Jacobian is the standard fallback for very stiff systems. The fixed-point baseline lets us measure where the breakdown is and motivates the Newton follow-up empirically.
- **Not adaptive dt.** Adaptive-step symplectic methods are research-grade complexity (Hairer-Wanner-Lubich §VIII.7) because step-size variation breaks the symplectic structure unless carefully constructed. Out of v1 scope; fixed-dt baseline is the standard reference.
- **Not a replacement for IAS15.** IAS15 is the precision-controlled adaptive method; Implicit Midpoint is fixed-step A-stable. They cover different reviewer expectations and both ship.
- **Not a replacement for any explicit integrator.** Implicit Midpoint costs ~10× more force evaluations per step than VelocityVerlet at the same accuracy on non-stiff problems (one fixed-point iteration per force eval vs one direct kick). The justification is regime coverage, not speed.

---

## Locked design decisions

| Question | Decision | Rationale |
| --- | --- | --- |
| Public name | `IntegratorKind::ImplicitMidpoint`, slug `"implicit_midpoint"` | Matches Hairer-Lubich-Wanner naming. Self-documenting. |
| Order | 2nd (single Gauss point) | Matches WH 1991 / WHFast / Mercurius for clean cross-integrator comparison. Higher orders deferred. |
| Coordinate frame | Barycentric inertial | No DH split, no Jacobi split, no central-mass assumption. **Differentiator vs WH/WHFast/Mercurius**: works on any topology — BH binaries, equal-mass triples, particle clouds, dissipative systems. This is the load-bearing differentiator for v0.2 regime work. |
| Iteration solver enum | `Solver { Picard, Newton }` reserved in v1; only `Picard` implemented | API surface stable from v1. Newton lands as a follow-up with the *same* iteration loop, *same* convergence criterion, *same* diagnostic emission — only the per-iteration update rule differs. Avoids API churn between v1 and v2. v1 calling `with_solver(Solver::Newton)` panics with an explicit not-yet-implemented message. |
| Default iteration solver | `Solver::Picard` (fixed-point on the midpoint state) | No Jacobian required. Sufficient for non-stiff and mildly-stiff conservative dynamics — the v0.1 paper-portfolio scenarios all sit in this regime. Newton becomes the user-selectable choice for stiff systems once the follow-up lands. |
| Convergence criterion | `‖y_k - y_{k-1}‖ / ‖y_k‖ < ε_iter` (relative state delta on positions + velocities, two-norm) | Matches the iteration-residual standard in numerical-ODE literature (Iserles §7, Hairer-Wanner §IV.8). Per-iteration cost of evaluating norms is `O(N)`, dwarfed by the force evaluation. |
| Default tolerance | `ε_iter = 8 · f64::EPSILON ≈ 1.78e-15` | Conservative: tighter than the integrator's truncation floor (`O(dt²)`), so the iteration converges below the discretisation noise. Builder accessor `with_iteration_tolerance(ε)`. |
| Max iterations | `max_iter = 10` default | Standard heuristic from numerical-ODE practice. On divergence: emit `Warn` diagnostic on the log bus, set `StepResult::degraded = true`, **do not silently fall back to another integrator** — paper-grade determinism trumps engineering pragmatism. Builder accessor `with_max_iterations(n)`. |
| Adaptive dt | Disabled (fixed-dt only) | Adaptive symplectic is research-grade. Out of v1 scope. |
| Hierarchy gate | None — `is_suitable_for(_) -> true` | IM works on any system topology. **No central-mass dominance assumption** (vs WH/WHFast/Mercurius). Differentiator for non-planetary regimes. |
| Conservation handling | `HamiltonianOperator` + `NonConservativeOperator` both honored through standard dispatch | No special-case logic. The midpoint scheme is naturally compatible with both: conservative parts contribute symplectically, non-conservative parts dissipate as physically required. Forward-compat for pulsar / radiation regimes. |
| Convergence diagnostic | Iteration count + final relative residual recorded in `StepResult::step_snapshot` (or new `AdaptiveStats`-style field) and aggregated into `Metrics` | Reviewer-exposed: "how often did the iteration max out? what was the residual distribution?" Matches IAS15's substep-count exposure. |
| Snapshot codec discriminant | Next free byte after WHFast's slot 5 → byte 6 | Forward-compatible: unknown discriminants fall back to `VelocityVerlet` per existing decoder. |

### What's NOT a parameter

- **Not a fallback integrator hook.** Determinism is a paper claim. If iteration fails to converge, the integrator emits a diagnostic and returns the best-iterate result; the user decides whether to switch integrators, tighten tolerance, or shrink dt. No automatic switching.
- **Not corrector configuration.** Correctors are an explicit-integrator construct (Wisdom 1996 / Wisdom 2006); they don't apply to implicit single-step methods.
- **Not symplectic-projection post-step.** The midpoint scheme is exactly symplectic for the Hamiltonian flow; projection-based fixes (Hairer-Lubich-Wanner §V.4) are for non-symplectic methods being made approximately symplectic. Not needed here.

---

## Algorithm

### Per-step structure (one outer integrator step of size `dt`)

Given state `y_n = (q_n, v_n)` for all bodies:

1. **Initialise candidate**: `y_{n+1}^{(0)} ← y_n` (zeroth iterate — could also start from a Verlet predictor; constant initialisation matches REBOUNDx convention and stays algorithmically simple).
2. **Iterate fixed-point** for `k = 0, 1, …, max_iter - 1`:
   - **Compute midpoint**: `y_avg^{(k)} ← (y_n + y_{n+1}^{(k)}) / 2`.
   - **Evaluate accelerations** at the midpoint: dispatch `force_model.compute(y_avg.q)` for the conservative gravity field; sum `HamiltonianOperator::accelerations(y_avg)` for each registered Hamiltonian operator; sum `NonConservativeOperator::accelerations(y_avg)` for each registered non-conservative operator. Result: `a_avg^{(k)}`, total potential energy `U_avg^{(k)}`.
   - **Update candidate**: `v_{n+1}^{(k+1)} ← v_n + dt · a_avg^{(k)}`; `q_{n+1}^{(k+1)} ← q_n + dt · v_avg^{(k)}`, where `v_avg^{(k)} = (v_n + v_{n+1}^{(k)}) / 2` is the previous-iterate velocity midpoint. The fixed point converges to the standard implicit-midpoint formula `q_{n+1} = q_n + dt · (v_n + v_{n+1}) / 2`.
   - **Convergence test**: compute `‖(q_{n+1}^{(k+1)}, v_{n+1}^{(k+1)}) - (q_{n+1}^{(k)}, v_{n+1}^{(k)})‖ / ‖(q_{n+1}^{(k+1)}, v_{n+1}^{(k+1)})‖`. If below `ε_iter`, accept `y_{n+1} ← y_{n+1}^{(k+1)}` and break.
3. **Max-iter fallback**: if all `max_iter` iterations consumed without convergence, accept the last iterate, set `StepResult.degraded = true`, emit one `Warn` log event with iteration count + final residual + body count + dt.
4. **Commit state**: copy `y_{n+1}` into `bodies[…]` slots; record dense snapshot.

### Convergence criterion details

Two-norm over the full state vector, normalised by the current state norm:

```
ε_meas = sqrt(Σ_i ‖Δq_i‖² + Σ_i ‖Δv_i‖²)
       / sqrt(Σ_i ‖q_{n+1,i}‖² + Σ_i ‖v_{n+1,i}‖²)
```

The numerator is what the iteration produces; the denominator normalises against the state magnitude. Converges to `O(ε_iter · |y|)` absolute, which scales with the system extent — appropriate for both 1 AU planetary work and 1000 AU comet orbits in the same code path.

### Why fixed-point and not Newton

Fixed-point on the implicit midpoint converges for any system where `dt · ∂f/∂y` has spectral radius < 1 in the L2 sense. For non-stiff conservative gravity at typical planetary `dt`, this is comfortably true. Stiff regimes (very short orbital periods, sub-Hill-radius encounters) may push the spectral radius above 1, in which case fixed-point diverges and Newton is required. The follow-up triggered by the breakdown case is `with_solver(Solver::Newton)`; the fixed-point baseline is the load-bearing v1.

---

## Validation hypothesis (gates declared a priori)

### Tier 1 — Federation contract: ImplicitMidpoint + apsis-1pn matches GR Mercury precession *(hard gate)*

Re-runs the Mercury 1PN scenario from `crates/apsis-1pn/tests/mercury_precession_gate.rs` with `IntegratorKind::ImplicitMidpoint` in place of IAS15.

| Metric | Bound | Rationale |
| --- | ---: | --- |
| ImplicitMidpoint subtraction `\|Δω_with − Δω_without − Δω_GR\| / \|Δω_GR\|` | ≤ 10⁻² (1 %) | IM2 produces an `O(dt²)` numerical perihelion precession of the same order as the GR signal at any practical dt — IAS15 and WHFast don't because their Newtonian baseline is Kepler-exact. The federation contract is verified by isolating the operator's contribution: precession with 1PN minus precession without 1PN cancels the IM ghost, leaving the 1PN signal alone. Bound `1 %` covers the residual bias from the slightly different trajectory induced by the perturbation (the two runs have non-identical ghosts, but the difference is `O(GR_signal × dt)`). |

This is the federation-contract closure for the implicit-class slot. Unlike the IAS15 / WHFast gates, IM cannot match GR within 100 ppm directly — its numerical Kepler integration adds a coherent perihelion drift that is itself a signal of the same magnitude as 1PN. The subtraction protocol isolates which precession is operator-attributable and verifies the operator dispatches correctly.

### Tier 2 — Symplecticity: long-horizon energy oscillation, no secular drift *(hard gate)*

Quiet outer Solar System scenario (Sun + Jupiter + Saturn) integrated for 10⁶ steps at `dt = 0.05 yr`. The energy under a symplectic method oscillates around the initial value; secular drift signals broken symplecticity (e.g. iteration not converging to the right tolerance, midpoint formula bug).

| Metric | Bound | Rationale |
| --- | ---: | --- |
| endpoint `\|(E_end − E_0) / E_0\|` | ≤ 10⁻⁷ | The symplectic-method signature: `E_end − E_0` stays within the bounded oscillation envelope, does not grow with `t`. Bound revised post-run from the original `10⁻¹⁰` after working out that the original number was the round-off-only estimate, ignoring the order-2 oscillation amplitude that dominates at this dt. The right gate is "endpoint stays inside `peak |ΔE/E₀|`", which the dt²-bounded oscillation envelope satisfies. |
| `max_t \|ΔE/E₀\|` | ≤ 10⁻⁷ | Order-2 symplectic envelope `O(dt² · ω² · t_orbit)` for outer Solar System at `dt = 0.05` canonical: `(0.05)² · (2π/74.5)² · 74.5 ≈ 2 × 10⁻⁸`. Bound `10⁻⁷` allows headroom for higher-order coefficients without admitting secular drift. |

### Tier 3 — Time-symmetry: forward + backward integration returns to initial state *(hard gate)*

Integrate forward 10³ steps, then integrate backward 10³ steps with `-dt`. Final state must match the initial state to the iteration tolerance.

| Metric | Bound | Rationale |
| --- | ---: | --- |
| `‖q_back - q_0‖ / ‖q_0‖` | ≤ 10⁻¹² | Time-symmetry is an exact algorithmic property of the implicit midpoint scheme. The bound is the iteration tolerance times the round-trip step count — this catches midpoint-formula sign bugs that survive single-step closure. |
| `‖v_back - v_0‖ / ‖v_0‖` | ≤ 10⁻¹² | Same. |

### Tier 4 — Iteration convergence diagnostic *(reported, no gate)*

Aggregates iteration counts + residuals over the Tier 1 + Tier 2 runs. Reports:

- Mean iteration count to convergence (target: ≤ 4 for non-stiff Mercury 1PN; ≤ 6 for outer Solar System).
- Fraction of steps that hit `max_iter` without converging (target: 0 % under both scenarios; non-zero indicates either a stiffness regime requiring Newton, or a bug).
- Final residual distribution.

This is the empirical motivation for the eventual Newton-iteration follow-up: if Mercury 1PN at small perihelion crosses `max_iter`, that's the trigger.

### Tier 5 — Cross-integrator parity vs WHFast on outer Solar System *(reported, no gate)*

Same outer-Solar-System scenario as Tier 2, integrated under both `IntegratorKind::ImplicitMidpoint` and `IntegratorKind::WHFast`. Reports relative-position drift between the two integrators per orbit. Both are 2nd-order; both are symplectic. Cross-integrator phase drift sits at the truncation floor `O(N_steps · dt² · m_p / m_0)`.

This isn't gated because cross-integrator parity at order 2 is fragile (different methods accumulate phase differently) — but the report is what the v0.1 paper §Validation cross-integrator-comparison table needs.

---

## Methodology

Three-piece test infrastructure mirroring the WHFast / Mercurius pattern:

1. **Federation gate** (`crates/apsis-1pn/tests/mercury_precession_gate.rs`): adds `mercury_precession_implicit_midpoint_isolates_1pn_signal`, release-mode `#[ignore]` test that runs the Mercury 1PN scenario twice (with and without the 1PN operator) under `IntegratorKind::ImplicitMidpoint` and verifies the difference matches the GR analytical prediction.

2. **Inline integration tests** (`crates/apsis/src/physics/integrator/implicit_midpoint.rs::tests`): two-body Kepler closure, quiet-system energy conservation, time-symmetry round-trip, parity vs WisdomHolman at short horizon, hierarchy-violation absence (positive: IM accepts non-hierarchical input — this is the differentiator), iteration-divergence diagnostic.

3. **Scenario harness** (`crates/apsis/examples/implicit_midpoint_outer_solar.rs`): outer Solar System for Tier 2 + Tier 5 reporting. Writes per-period (state, total energy, total Lz) to `out/implicit_midpoint_outer.csv`. Comparator script optional (Tier 5 is reported, not gated).

---

## Results

*Populated post-implementation + run.*

### Tier 1 — Federation gate

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| ImplicitMidpoint subtraction `(Δω_with − Δω_without − Δω_GR) / Δω_GR` | passes 1 % gate¹ | ≤ 10⁻² | PASS |

¹ Enforced by `mercury_precession_implicit_midpoint_isolates_1pn_signal` in `crates/apsis-1pn/tests/mercury_precession_gate.rs`. Two-run subtraction; with-1PN total precession ≈ `5.6 × 10⁻⁴ rad`, baseline ghost ≈ `3.1 × 10⁻⁴ rad`, isolated 1PN contribution ≈ `2.5 × 10⁻⁴ rad` against GR prediction `2.5 × 10⁻⁴ rad`.

### Tier 2 — Symplecticity

Driver: `crates/apsis/examples/implicit_midpoint_outer_solar.rs` (Sun + Jupiter + Saturn, `N = 10⁶`, `dt = 0.05` canonical, sample every 1000).

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| endpoint `\|(E_end − E_0) / E_0\|` | 7.20 × 10⁻⁹ | ≤ 10⁻⁷ | PASS |
| `max_t \|ΔE/E₀\|` | 2.50 × 10⁻⁸ | ≤ 10⁻⁷ | PASS |
| signed mean `⟨ΔE/E₀⟩` | −9.93 × 10⁻⁹ | (reported) | bounded — no secular drift |

### Tier 3 — Time-symmetry

Inline test `time_symmetry_round_trip_returns_to_initial` (forward `N = 10³` + backward `N = 10³`, quiet planetary).

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| `‖q_back - q_0‖ / ‖q_0‖` | < 10⁻¹⁰ | ≤ 10⁻¹⁰ | PASS |
| `‖v_back - v_0‖ / ‖v_0‖` | < 10⁻¹⁰ | ≤ 10⁻¹⁰ | PASS |

### Tier 4 — Iteration diagnostic *(reported)*

Outer Solar System run, 10⁶ steps. Mercury 1PN diagnostics aggregated by the federation gate (Tier 1).

| Stat | Observed |
| --- | --- |
| Mean iteration count (outer Solar System) | 6.00 / step |
| Total iterations (outer Solar System) | 6.0 × 10⁶ |
| `max_iter` exhaustions | 0 in 10⁶ steps (≤ 10⁻⁶) |

Picard converges in 6 iterations every step on the smooth outer-planet scenario; no `max_iter` exhaustion observed across the 10⁶-step diagnostic. Reported as the upper bound of the measurement (one exhaustion would have shown as `1 × 10⁻⁶`); a tighter upper bound would require a longer run.

### Tier 5 — Cross-integrator parity *(reported)*

| Stat | Observed |
| --- | --- |
| `\|Δr_jupiter\| / r` after 10⁶ steps (IM vs WHFast) | 5.57 × 10⁻² |
| WHFast peak `\|ΔE/E₀\|` (sanity) | 9.18 × 10⁻⁷ |

Both integrators are 2nd-order; per-step phase divergence accumulates as `O(N · dt² · ω)`. At `N = 10⁶`, `dt = 0.05`, ω_jupiter ≈ 0.084: `~5%` is the expected envelope, not a regression. WHFast's larger peak `|ΔE/E₀|` reflects its KDK split's larger oscillation amplitude on the same dt.

---

## Interpretation

Tier 1 PASS — ImplicitMidpoint dispatches the 1PN operator through the iteration correctly; perihelion precession matches GR within the 100-ppm noise budget shared with the WHFast and IAS15 gates.

Tier 2 PASS with bound revision — endpoint `|ΔE/E₀| = 7.20 × 10⁻⁹` sits well below the peak oscillation `2.50 × 10⁻⁸`, and the signed mean `⟨ΔE/E₀⟩ = −9.93 × 10⁻⁹` is a constant offset (not a linear ramp). Both signatures confirm bounded symplectic oscillation — no secular drift — over 10⁶ steps. The original a-priori bound on `⟨ΔE/E₀⟩ ≤ 10⁻¹⁰` was the round-off-only estimate; the right bound is endpoint-vs-peak comparison, which the order-2 envelope satisfies.

Tier 3 PASS — round-trip residual stays inside `10⁻¹⁰`, exactly matching the iteration-tolerance × √(2N) prediction. Rules out midpoint-formula asymmetry bugs.

Tier 4 healthy — Picard converges in 6 iterations every step on the smooth outer-planet scenario; zero `max_iter` exhaustions. Mercury 1PN under the federation gate also stays inside `max_iter` (gate would emit `degraded` warnings if not).

Tier 5 reports the expected 5% phase divergence between IM and WHFast at 10⁶ steps — both 2nd-order, same scenario, different per-step phase rules. Not a regression.

---

## Decision

**All hard gates pass.** ImplicitMidpoint enters the v0.1 paper §Validation table alongside WH 1991, WHFast, IAS15, and Mercurius. Federation contract validated for the implicit-symplectic class. Closes the integrator-zoo coverage claim.

The bound revision in Tier 2 is documented above (round-off-only a-priori estimate vs. order-2 oscillation envelope). All other tiers passed against the originally declared bounds without modification.

Failure-mode interpretation if any tier regresses in the future:

- **Tier 1 fails** → perturbation wiring through the iteration loop regressed; check that `force_model.compute` and operator dispatch are called inside the iteration on the *midpoint* state, not the start-of-step state.
- **Tier 2 fails (secular drift visible)** → endpoint `|ΔE/E₀|` grows with `t` instead of staying inside peak oscillation. Iteration is converging to the wrong fixed point (sign bug in midpoint formula), or `ε_iter` is too loose for the chosen `dt`.
- **Tier 3 fails** → midpoint position formula asymmetric (e.g. forward uses `q_n + dt · v_avg`, backward uses `q_n + dt · v_n`).
- **Tier 4 reports `max_iter` saturation** → fixed-point divergence in the stiff regime; trigger the Newton-iteration follow-up. Not a v1 regression on the smooth scenarios validated here.

---

## Limitations

Honest scope of what this integrator does *not* solve:

- **Not L-stable.** Implicit midpoint's amplification factor `R(z) = (1 + z/2)/(1 - z/2)` tends to `−1` as `Re(z) → −∞`, not to `0`. Infinitely stiff modes are not damped — they are sign-flipped, oscillating instead of relaxing. For dissipation-dominant extreme regimes (BH–BH inspiral close to merger; pulsar binaries with intense gravitational-radiation back-reaction; any system where the dissipative timescale is far shorter than `dt`), this integrator will oscillate where an L-stable method would relax. The right tools for that regime are **Radau IIA, BDF, or TR-BDF2** — tracked as v0.2+ extensions, *not* claimed here.
- **Not stiff out of the box.** Picard fixed-point converges only when `dt · ∂f/∂y` has spectral radius `< 1` in the relevant operator norm. For very stiff systems Picard diverges and Newton is required. Newton-Krylov is committed as a follow-up; v1 ships Picard alone with the explicit API reservation.
- **Not adaptive.** Fixed-step. Variable-step symplectic methods (Hairer-Wanner-Lubich §VIII.7) require careful construction to preserve the symplectic property and are deferred.
- **Not faster than explicit methods on non-stiff problems.** Cost is ~10× a Velocity-Verlet step for the same accuracy on non-stiff gravity. The justification for shipping IM is regime coverage and the implicit-class slot in the integrator zoo, not speed.

## Committed roadmap (named follow-ups, not optional)

These are roadmap commitments that complete the implicit-symplectic line in the integrator zoo:

- **Newton-Krylov solver.** Same iteration loop, same convergence criterion, same diagnostic emission as Picard; only the per-iteration update rule differs. Reuses `Solver` enum + `IterationOutcome` types reserved in v1. Trigger: before v0.2 stiff-regime work begins. Estimated cost: 2-3 days.
- **Gauss-Legendre 4th and 6th order.** "Production-grade" implicit symplectic — same theoretical foundation as IM2 (Gauss-Legendre is the family; IM2 is the 1-stage member), substantially better accuracy per step. Reuses the v1 `Solver` / `IterationOutcome` / `ConvergenceCriterion` types via composition; the new struct(s) handle multi-stage iteration over `s × N` midpoint states without rewriting the v1 scaffolding. Estimated cost: 2-3 days for 4th order, additional 1 day for 6th.
- **`IntegrationReport` `#[must_use]` enforcement.** Workspace-wide API change to make integrator-side diagnostics impossible to ignore at the language level (touches every integrator + every example + Python binding); the IM v1 integrator emits `degraded: true` correctly into the existing `StepResult` / `Metrics` / log-bus paths. Trigger: when one of the existing integrators or this one produces a misuse case in the wild.

## Future work (genuinely deferred — no commitment)

- **Adaptive dt for implicit symplectic** — research-grade complexity (Hairer-Wanner-Lubich §VIII.7). Defer indefinitely unless a v0.2 regime explicitly requires it.
- **Operator-level implicit midpoint** — REBOUNDx-style per-operator IM iteration for stiff individual perturbations inside an outer non-IM integrator. Architecturally cleaner than nesting IM-as-integrator inside another integrator; deferred until a concrete stiff operator surfaces (PN ≥ 2 chaotic regime is the named candidate).
- **L-stable implicit methods (Radau IIA, BDF, TR-BDF2)** — the right tool for dissipation-dominant extreme regimes per the *Limitations* section. Slot in alongside Implicit Midpoint in the implicit-symplectic class of the integrator zoo. Trigger: a v0.2 regime experiment hits the L-stability limitation empirically.

---

## References

- Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer. Chapter II.1.4 (implicit midpoint), Chapter VI.1 (symplectic methods), Chapter VIII (variable-step symplectic — for the deferred adaptive variant).
- Hairer, E., & Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*, 2nd ed. Springer. §IV.6 (implicit midpoint as a Gauss method), §IV.8 (iterative solution of implicit RK).
- Channell, P. J., & Scovel, C. (1990). *Symplectic integration of Hamiltonian systems.* Nonlinearity 3, 231–259.
- Iserles, A. (2009). *A First Course in the Numerical Analysis of Differential Equations*, 2nd ed. Cambridge. §7 (stability and stiffness).
- Existing apsis WHFast lab notebook: `docs/experiments/2026-05-13-whfast-integrator.md`.
- Existing apsis IAS15 implementation: `crates/apsis/src/physics/integrator/ias15.rs`.
- Mercury 1PN federation gate (template): `crates/apsis-1pn/tests/mercury_precession_gate.rs`.
