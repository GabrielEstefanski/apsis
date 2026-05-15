//! Scenario specifications for the IAS15 benchmark harness.
//!
//! Each builder returns a [`ScenarioSpec`] with fully-deterministic
//! initial conditions and a fixed integration window. There is no
//! dependency on the runtime side (no `System`, no force model); that
//! wiring is done in [`super::runner`], keeping the scenario catalog
//! pure data — easy to reason about, trivially testable, and safe to
//! read from both the validation and recording code paths.
//!
//! # Adding a scenario
//!
//! 1. Add a builder fn returning `ScenarioSpec`.
//! 2. Register it in [`all`]. The order there determines the order
//!    entries appear in the baseline file.
//! 3. Run `IAS15_BENCH_UPDATE_BASELINE=1 cargo bench` to record an
//!    initial entry in `benches/baselines/ias15.toml`.
//!
//! # On close-encounter coverage
//!
//! We intentionally include Kepler e=0.99 (pericenter at 0.01·a) as a
//! *controlled* close-encounter scenario in addition to the chaotic
//! Pythagorean three-body. The two cover different failure modes:
//! Kepler e=0.99 stresses the controller's `dt` shrink/grow cycle
//! around a reproducible pericenter passage each orbit; Pythagorean
//! exercises multiple close encounters in a row without reference
//! trajectory — more robust to subtle trajectory-altering bugs, but
//! harder to interpret when a regression fires.

use apsis::domain::body::Body;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

/// Fully-specified benchmark scenario: initial conditions + the size of
/// the integration window used for baseline validation.
///
/// The `name` is the stable identifier used as the section key in
/// `benches/baselines/ias15.toml`. Renaming a scenario therefore
/// requires a baseline update — by design, so the regression gate
/// cannot silently ignore a renamed scenario.
pub struct ScenarioSpec {
    pub name: &'static str,
    pub bodies: Vec<Body>,
    pub dt_budget: f64,
    pub duration: f64,
    /// Whether this scenario participates in the Criterion wall-time
    /// benchmark group. `false` for scenarios whose single sub-step is
    /// expensive enough that Criterion's sample batching (typically 100
    /// iterations × `STEPS_PER_ITER` sub-steps) would push total bench
    /// wall time beyond a few minutes. Such scenarios still run through
    /// the validation + phase-profile path, which is the diagnostic
    /// signal we care about for them — the Criterion-specific
    /// statistical comparisons are not worth the wait at large N.
    pub criterion_bench: bool,

    /// Whether mismatches against the stored baseline are a regression
    /// (exit 2) or an advisory (scenario runs, phase profile prints,
    /// exit stays 0). `false` for scenarios still under active
    /// diagnosis — the metrics are expected to shift as the
    /// investigation progresses, and a hard gate would block progress
    /// without adding signal. Flip back to `true` once the scenario
    /// has a known-good baseline that a PR is supposed to preserve.
    pub gate_on_baseline: bool,
}

/// Two-body Keplerian orbit at moderate eccentricity (e=0.5).
///
/// The steady-state baseline: warmstart is maximally effective and
/// Picard should converge in 2–3 iterations per sub-step. A regression
/// here points at a change in the controller's quiescent behaviour.
pub fn kepler_e05() -> ScenarioSpec {
    kepler_two_body("kepler_e05", 2.0, 0.5, 2.0, 100)
}

/// Two-body Keplerian orbit at e=0.9.
///
/// Pericenter at 0.1·a — the controller must shrink `dt` aggressively
/// at each pericenter passage and re-grow it over the apocenter arc.
/// Tests the full shrink/grow cycle every orbit.
pub fn kepler_e09() -> ScenarioSpec {
    kepler_two_body("kepler_e09", 1.0, 0.9, 2.0, 50)
}

/// Two-body Keplerian orbit at e=0.99 (controlled close encounter).
///
/// Pericenter at 0.01·a drives the step size down by ~2 orders of
/// magnitude during passage. This is where Picard tends to stagger
/// (stagnation guard + retry path are exercised heavily) and where
/// round-off accumulation can leak into energy conservation.
pub fn kepler_e099() -> ScenarioSpec {
    kepler_two_body("kepler_e099", 1.0, 0.99, 2.0, 20)
}

/// Pythagorean three-body problem (Burrau 1913): masses 3-4-5 at rest
/// on the vertices of a 3-4-5 triangle.
///
/// Chaotic with violent close encounters between t≈2 and t≈5 followed
/// by further encounters; stress-tests rejection rollback across
/// multiple concurrent close approaches. The tight integration window
/// (t∈[0,10]) covers the strongest encounter phase without letting
/// chaos amplify bit-level noise into macroscopic trajectory divergence.
pub fn pythagorean() -> ScenarioSpec {
    let bodies = vec![
        Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0),
        Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0),
        Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0),
    ];
    ScenarioSpec {
        name: "pythagorean",
        bodies,
        dt_budget: 0.1,
        duration: 10.0,
        criterion_bench: true,
        gate_on_baseline: true,
    }
}

/// Uniform disk cluster of 50 equal-mass bodies with circular velocities
/// around the centre of mass. Seeded random layout under the default
/// exact `NewtonKernel`.
///
/// This is the first scenario in the harness with non-trivial N. Its
/// purpose is to expose how the IAS15 phases scale with body count —
/// specifically, the transition where `evaluate` (O(N²) pairwise sum
/// below `EXACT_THRESHOLD = 64`) overtakes the other phases, and to
/// give `update_g_and_b` (linear in N) a large enough body axis to
/// make SIMD / SoA layout arguments testable rather than speculative.
///
/// The seed is pinned (`0xc1a55e1` — "cluster seed") so the scenario
/// is bit-deterministic across runs on the same machine. N=50 sits
/// just under the Barnes-Hut crossover so the force path is pure
/// O(N²) — useful both as a baseline for future BH comparisons and
/// because tree-build overhead would otherwise dominate a single
/// short scenario.
pub fn cluster_n50() -> ScenarioSpec {
    const N: usize = 50;
    const R_DISK: f64 = 1.0;
    const SEED: u64 = 0xc1a55e1;

    // Total mass 1 → each body mass 1/N.
    let mass_per_body = 1.0 / N as f64;

    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut bodies = Vec::with_capacity(N);

    for _ in 0..N {
        // Uniform disk sampling: θ ~ U(0, 2π), r ~ sqrt(U(0,1)) · R.
        // The square root on r is the inverse CDF for uniform area
        // density (a linear-in-U sampling would concentrate bodies
        // near the centre).
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        let r = R_DISK * rng.random::<f64>().sqrt();
        let x = r * theta.cos();
        let y = r * theta.sin();

        // Circular velocity magnitude for an enclosed-mass-proportional
        // potential (which a uniform disk approximates at leading
        // order): v_circ(r) = sqrt(M_enc · G / r) where M_enc ∝ r²
        // gives v_circ ∝ r. Direction is tangential (CCW).
        //
        // This is not the self-consistent solution — the real
        // uniform disk has a harmonic potential, not Keplerian — but
        // it produces a bound, non-pathological initial state whose
        // dynamics exercise the integrator across a mix of tight and
        // loose pairs. The benchmark target is timing, not orbital
        // perfection.
        let v_mag = r; // sqrt(G · r · M_tot / R²) with G = M = R = 1
        let vx = -v_mag * theta.sin();
        let vy = v_mag * theta.cos();

        let b = Body::rocky(mass_per_body).at(x, y).with_velocity(vx, vy);
        bodies.push(b);
    }

    ScenarioSpec {
        name: "cluster_n50",
        bodies,
        // dt_budget comfortably above the controller's natural step
        // (confirmed empirically after recording the baseline —
        // dt_p95 << dt_budget, so the budget acts only as an upper
        // cap during warm-up, never constrains steady state).
        dt_budget: 0.05,
        // Short duration by design: chaotic N-body dynamics amplify
        // bit-level noise into trajectory divergence on long
        // integrations. Bench relevance — cost per sub-step, phase
        // distribution — is fully exercised in a fraction of a
        // dynamical time.
        duration: 0.5,
        criterion_bench: true,
        gate_on_baseline: true,
    }
}

/// Central-body + 640 test-particle disk approximating the interactive
/// app's `solar_system` preset (Sun + planets + asteroid belt +
/// comets ≈ 641 bodies). Purpose: diagnose IAS15 stutter reported at
/// this scale in normal playback mode.
///
/// Simplifications relative to the full preset:
///
///   * All test particles are equal-mass (1e-10 solar masses) rather
///     than a mix of planets + asteroids + comets. The goal is
///     allocator pressure and phase distribution at N≈641, not
///     trajectory fidelity.
///   * Single annulus [0.5, 5] AU instead of the preset's
///     planets-and-belt structure. Keplerian circular velocities
///     around the central mass.
///   * `G = 1` rather than `G = 4π²` (solar-AU-year). The preset's
///     physical correctness is preserved elsewhere; the bench only
///     cares about wall-time / alloc behaviour, which is invariant
///     under uniform velocity / time scaling.
///
/// Marked `criterion_bench: false` — at N=641 a single Criterion
/// iteration (`STEPS_PER_ITER = 100` accepted sub-steps) takes on
/// the order of seconds; a 100-sample Criterion run would consume
/// tens of minutes. The validation + phase-profile path provides
/// the diagnostic signal we need (per-phase breakdown with the new
/// `a0_clone` and `dense_snapshot_build` timers).
pub fn solar_n641() -> ScenarioSpec {
    const N_TEST: usize = 640; // + 1 central body = 641 total
    const M_CENTRAL: f64 = 1.0;
    const M_TEST: f64 = 1e-10;
    // Annulus chosen to approximate the interactive app's
    // `solar_system` preset asteroid belt ([2.2, 3.2] AU), with a bit
    // of margin on each side so the controller exercises a range of
    // dynamical periods rather than a single narrow shell.
    //
    // ## On the regime and its limits
    //
    // Under the default exact `NewtonKernel` (ε = 0),
    // N=640 random-uniform test particles in a thin annulus concentrate
    // enough close pairs that the adaptive dt
    // shrinks toward DT_MIN around close-encounter events, yielding a
    // baseline `peak_energy_err` of order 10â»â´ — NOT representative of
    // IAS15's machine-precision regime.
    //
    // This is a deliberate *regime stress test*, not a quality
    // benchmark. It was introduced to reproduce the stutter reported
    // in the interactive app's 641-body preset and to make the
    // rejection cascade (pre-RMS-norm-fix: 194% rejection rate)
    // measurable. The stress test served its purpose — the RMS-norm
    // fix documented in `docs/experiments/
    // 2026-04-22-solar-system-stutter-diagnosis.md` reduced rejections
    // by 36% and wall time by 23% on this scenario. A separate
    // controlled-quality scenario at N=641 (e.g. structured rings or
    // a well-spaced Keplerian disk) is future work.
    //
    // Because the scenario sits outside IAS15's efficient regime,
    // `gate_on_baseline` is set to `false`: metric shifts are
    // expected across algorithmic changes and should be reviewed as
    // advisories, not treated as regressions.
    const R_INNER: f64 = 1.5;
    const R_OUTER: f64 = 3.5;
    const SEED: u64 = 0x501a5; // "solaš"

    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut bodies = Vec::with_capacity(N_TEST + 1);

    // Central body: massive star at origin, at rest.
    bodies.push(Body::star(M_CENTRAL).at(0.0, 0.0).with_velocity(0.0, 0.0));

    for _ in 0..N_TEST {
        // Uniform in [R_INNER, R_OUTER] via rejection-free sampling.
        let r = R_INNER + (R_OUTER - R_INNER) * rng.random::<f64>();
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        let x = r * theta.cos();
        let y = r * theta.sin();

        // Circular Keplerian velocity around the central body:
        // v_circ = sqrt(G · M_central / r). With G = 1, v = 1/sqrt(r).
        let v_circ = (M_CENTRAL / r).sqrt();
        let vx = -v_circ * theta.sin();
        let vy = v_circ * theta.cos();

        let b = Body::asteroid(M_TEST).at(x, y).with_velocity(vx, vy);
        bodies.push(b);
    }

    ScenarioSpec {
        name: "solar_n641",
        bodies,
        // dt_budget generous: the controller will shrink it at the
        // inner-radius bodies (shortest orbital periods). Chosen to
        // leave headroom for the adaptive mechanism rather than
        // clipping it.
        dt_budget: 0.05,
        // Very short duration: at N=641 with IAS15 and BH each
        // accepted sub-step takes on the order of tens of ms, so
        // even 0.1 time units produces enough sub-steps (~hundreds)
        // for reliable phase statistics. Longer windows only add
        // wall time without improving the diagnostic signal.
        duration: 0.1,
        criterion_bench: false,
        gate_on_baseline: false,
    }
}

/// Central body + 640 test particles arranged as 20 concentric rings of
/// 32 bodies each — structured, well-separated geometry intended to
/// exercise IAS15 at solar-system-class N without the close-encounter
/// cascade that `solar_n641` exhibits.
///
/// ## Status: kept out of the default catalog
///
/// This scenario is **not** registered in [`all`]. Building it landed
/// an unrelated architectural finding instead: IAS15 paired with
/// Barnes-Hut cascades at large N regardless of the scenario's
/// scenario-level quality, because BH's tree approximation is not a
/// deterministic function of state across Picard iterations. After
/// that finding, the integrator/force pairing became a first-class
/// concern (see `System::set_integrator` and the integrator
/// execution-profile ADR), and the app's default integrator moved
/// from IAS15 to Yoshida 4.
///
/// The function is kept because the scenario itself remains valuable:
///
///   * It is the cleanest N=641 stress the harness has — regular
///     geometry, bounded pair separations, no birthday-problem close
///     encounters. Useful for future work that wants to isolate
///     integrator behaviour at high N from scenario stiffness.
///   * Running it under IAS15 with the auto-switch to direct O(N²)
///     gives a clean reading of IAS15's in-regime quality at N=641.
///   * Comparing its metrics between force models (direct vs a
///     hypothetical smooth multipole method in the future) is
///     exactly the experiment the scenario was designed for.
///
/// To use: call [`structured_rings_n641`] directly from a diagnostic
/// script or add it temporarily to [`all`] for a one-off recording.
///
/// ## Geometry
///
/// * Radial range: R_INNER = 1.5, R_OUTER = 3.5 (matches
///   `solar_n641` for an apples-to-apples comparison at N=641).
/// * 20 rings, linearly spaced in `r` → Î”r ≈ 0.105.
/// * 32 equally-spaced bodies per ring → minimum intra-ring chord at
///   R_INNER = 2·R·sin(π/32) ≈ 0.294.
/// * Per-ring pseudo-random phase offset (seeded) to break radial
///   alignment between rings without destroying their co-rotating
///   structure.
/// * Softening = 0.03 → intra-ring margin 9.8×, inter-ring margin
///   3.5×. Co-rotating rings preserve intra-ring separation to
///   leading order; radial ordering preserves inter-ring separation.
///
/// ## Expected in-regime targets (when IAS15+direct is used)
///
/// * `peak_energy_err` ≤ 1e-11 (machine-precision class at f64).
/// * `degraded_total` = 0 (controller never saturates DT_MIN).
/// * `dt_min` ≫ DT_MIN = 1e-12.
/// * Rejection rate < 10% of accepted substeps.
#[allow(dead_code)] // kept out of the default `all()` catalog on purpose — see doc.
pub fn structured_rings_n641() -> ScenarioSpec {
    const N_RINGS: usize = 20;
    const BODIES_PER_RING: usize = 32;
    const N_TEST: usize = N_RINGS * BODIES_PER_RING; // 640
    const M_CENTRAL: f64 = 1.0;
    const M_TEST: f64 = 1e-10;
    const R_INNER: f64 = 1.5;
    const R_OUTER: f64 = 3.5;
    const SEED: u64 = 0x1166_5EED; // "rings seed"

    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut bodies = Vec::with_capacity(N_TEST + 1);

    bodies.push(Body::star(M_CENTRAL).at(0.0, 0.0).with_velocity(0.0, 0.0));

    let dr = (R_OUTER - R_INNER) / (N_RINGS as f64 - 1.0);
    let dtheta = std::f64::consts::TAU / BODIES_PER_RING as f64;

    for k in 0..N_RINGS {
        let r = R_INNER + dr * k as f64;
        // Per-ring phase offset breaks radial alignment between rings
        // while preserving uniform angular spacing within each ring.
        let phase = rng.random::<f64>() * std::f64::consts::TAU;
        let v_circ = (M_CENTRAL / r).sqrt();

        for j in 0..BODIES_PER_RING {
            let theta = phase + j as f64 * dtheta;
            let x = r * theta.cos();
            let y = r * theta.sin();
            let vx = -v_circ * theta.sin();
            let vy = v_circ * theta.cos();

            let b = Body::asteroid(M_TEST).at(x, y).with_velocity(vx, vy);
            bodies.push(b);
        }
    }

    ScenarioSpec {
        name: "structured_rings_n641",
        bodies,
        // Budget generous — the controller's natural dt at this scale
        // is set by the innermost orbit (T_inner ≈ 11.5); budget acts
        // only as an upper cap during warm-up.
        dt_budget: 0.1,
        // ~half the innermost orbital period (T_inner/2 ≈ 5.77) —
        // long enough for the controller to settle and for a credible
        // `energy_drift_slope` fit.
        duration: 6.0,
        criterion_bench: false,
        gate_on_baseline: false,
    }
}

/// Ordered list of scenarios registered in the default bench catalog.
/// The order determines the order of sections in the baseline file.
///
/// Scenarios defined in this module but *not* included here
/// (currently: [`structured_rings_n641`]) are kept for ad-hoc
/// diagnostic runs and future investigations — they remain callable
/// directly by name from scripts or temporary additions to this list.
pub fn all() -> Vec<ScenarioSpec> {
    vec![kepler_e05(), kepler_e09(), kepler_e099(), pythagorean(), cluster_n50(), solar_n641()]
}

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Build a two-body Keplerian orbit with the specified semi-major axis,
/// eccentricity, and total mass parameter μ. The bodies are placed on
/// either side of the origin at pericenter with equal-magnitude opposite
/// velocities (symmetric about the centre of momentum).
fn kepler_two_body(name: &'static str, a: f64, e: f64, mu: f64, n_orbits: u64) -> ScenarioSpec {
    let r_peri = a * (1.0 - e);
    let v_peri = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();
    let period = 2.0 * std::f64::consts::PI * (a.powi(3) / mu).sqrt();

    let b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
    let b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);

    ScenarioSpec {
        name,
        bodies: vec![b1, b2],
        // Budget of period/20 lets the controller settle near its
        // natural step size (~period/30 at ε=1e-9 for e=0.5) without
        // the budget acting as a cap.
        dt_budget: period / 20.0,
        duration: n_orbits as f64 * period,
        criterion_bench: true,
        gate_on_baseline: true,
    }
}
