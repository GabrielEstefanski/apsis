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

use gravity_sim::domain::body::Body;
use gravity_sim::domain::materials::Material;
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
    let mut bodies = vec![
        Body::new(1.0, 3.0, 0.0, 0.0, 3.0, Material::Rocky),
        Body::new(-2.0, -1.0, 0.0, 0.0, 4.0, Material::Rocky),
        Body::new(1.0, -1.0, 0.0, 0.0, 5.0, Material::Rocky),
    ];
    for b in &mut bodies {
        b.softening = 0.0;
    }
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
/// around the centre of mass. Seeded random layout, softening applied
/// to keep pairwise forces bounded when bodies pass close to each other.
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
    const SOFTENING: f64 = 0.02;
    const SEED: u64 = 0xc1a55e1;

    // Total mass 1 → each body mass 1/N. Keeps pairwise force at
    // softened close approach bounded by m_i·m_j / ε² = (1/N²) / ε²
    // which is well within f64 dynamic range for N=50, ε=0.02.
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
        // order): v_circ(r) = sqrt(M_enc · G / r) where M_enc ∝ r²
        // gives v_circ ∝ r. Direction is tangential (CCW).
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

        let mut b = Body::new(x, y, vx, vy, mass_per_body, Material::Rocky);
        b.softening = SOFTENING;
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
    // Even at this reduced spread and with softening = 0.05,
    // N=640 random-uniform test particles in a thin annulus concentrate
    // some pairs below the softening length. The adaptive dt therefore
    // shrinks toward DT_MIN around close-encounter events, yielding a
    // baseline `peak_energy_err` of order 10⁻⁴ — NOT representative of
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
    const SOFTENING: f64 = 0.05;
    const SEED: u64 = 0x501a5; // "solaš"

    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut bodies = Vec::with_capacity(N_TEST + 1);

    // Central body: massive star at origin, at rest.
    bodies.push(Body::new(0.0, 0.0, 0.0, 0.0, M_CENTRAL, Material::Star));

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

        let mut b = Body::new(x, y, vx, vy, M_TEST, Material::Asteroid);
        b.softening = SOFTENING;
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

/// Ordered list of all scenarios. The order determines the order of

/// Ordered list of all scenarios. The order determines the order of
/// sections in the baseline file.
pub fn all() -> Vec<ScenarioSpec> {
    vec![
        kepler_e05(),
        kepler_e09(),
        kepler_e099(),
        pythagorean(),
        cluster_n50(),
        solar_n641(),
    ]
}

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Build a two-body Keplerian orbit with the specified semi-major axis,
/// eccentricity, and total mass parameter μ. The bodies are placed on
/// either side of the origin at pericenter with equal-magnitude opposite
/// velocities (symmetric about the centre of momentum).
fn kepler_two_body(
    name: &'static str,
    a: f64,
    e: f64,
    mu: f64,
    n_orbits: u64,
) -> ScenarioSpec {
    let r_peri = a * (1.0 - e);
    let v_peri = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();
    let period = 2.0 * std::f64::consts::PI * (a.powi(3) / mu).sqrt();

    let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
    b1.softening = 0.0;
    let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
    b2.softening = 0.0;

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
