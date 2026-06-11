//! Scenario specifications for the IAS15 benchmark harness: pure data
//! (deterministic ICs + integration window); the runtime wiring lives
//! in [`super::runner`].
//!
//! Adding a scenario: builder fn → register in [`all`] (order = order
//! in the baseline file) → `IAS15_BENCH_UPDATE_BASELINE=1 cargo bench`
//! to record the initial entry.
//!
//! Close-encounter coverage is deliberately double: Kepler e=0.99
//! stresses a *reproducible* per-orbit dt shrink/grow cycle; the
//! Pythagorean stresses repeated encounters with no reference
//! trajectory — broader coverage, harder to interpret on regression.

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
    /// Whether this scenario joins the Criterion wall-time group.
    /// `false` where sample batching would cost minutes at large N;
    /// the validation + phase-profile path still runs.
    pub criterion_bench: bool,

    /// Whether baseline mismatches are a regression (exit 2) or an
    /// advisory. `false` while a scenario is under active diagnosis;
    /// flip to `true` once it has a known-good baseline to preserve.
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

/// Seeded uniform-disk cluster, 50 equal masses on circular orbits.
/// Exposes how the IAS15 phases scale with N (where the O(N²)
/// `evaluate` overtakes the rest). N=50 sits under the Barnes–Hut
/// crossover, so the force path stays pure O(N²).
pub fn cluster_n50() -> ScenarioSpec {
    const N: usize = 50;
    const R_DISK: f64 = 1.0;
    const SEED: u64 = 0xc1a55e1;

    // Total mass 1 → each body mass 1/N.
    let mass_per_body = 1.0 / N as f64;

    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut bodies = Vec::with_capacity(N);

    for _ in 0..N {
        // sqrt on r = inverse CDF for uniform area density.
        let theta = rng.random::<f64>() * std::f64::consts::TAU;
        let r = R_DISK * rng.random::<f64>().sqrt();
        let x = r * theta.cos();
        let y = r * theta.sin();

        // v_circ grows linearly with r for an enclosed-mass-proportional
        // potential. Not the self-consistent disk solution; the target is
        // a bound, non-pathological state for timing, not orbital fidelity.
        let v_mag = r;
        let vx = -v_mag * theta.sin();
        let vy = v_mag * theta.cos();

        let b = Body::rocky(mass_per_body).at(x, y).with_velocity(vx, vy);
        bodies.push(b);
    }

    ScenarioSpec {
        name: "cluster_n50",
        bodies,
        dt_budget: 0.05, // above the controller's natural step; caps warm-up only
        // Short on purpose: per-sub-step cost and phase distribution are
        // fully exercised in a fraction of a dynamical time.
        duration: 0.5,
        criterion_bench: true,
        gate_on_baseline: true,
    }
}

/// Central body + 640 test particles in a single annulus,
/// approximating the interactive app's 641-body `solar_system`
/// preset. Built to diagnose the IAS15 stutter reported at this
/// scale; equal masses, plain annulus, and G = 1 on purpose — the
/// bench measures wall-time/alloc behaviour, not trajectory
/// fidelity. `criterion_bench: false`: one Criterion iteration
/// takes seconds at this N; the validation + phase-profile path
/// carries the diagnostic signal.
pub fn solar_n641() -> ScenarioSpec {
    const N_TEST: usize = 640; // + 1 central body = 641 total
    const M_CENTRAL: f64 = 1.0;
    const M_TEST: f64 = 1e-10;
    // Annulus around the preset's asteroid belt with margin so the
    // controller sees a range of dynamical periods. Deliberate regime
    // stress test: random close pairs under the exact kernel drive dt
    // toward DT_MIN, so energy error is far above IAS15's smooth-flow
    // floor and `gate_on_baseline` stays false. See
    // docs/experiments/2026-04-22-solar-system-stutter-diagnosis.md.
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
        dt_budget: 0.05,
        // Hundreds of sub-steps — enough for phase statistics.
        duration: 0.1,
        criterion_bench: false,
        gate_on_baseline: false,
    }
}

/// Central body + 640 test particles as 20 concentric rings of 32 —
/// structured, well-separated N=641 without `solar_n641`'s
/// close-encounter cascade. Kept out of [`all`]: building it exposed
/// the IAS15 + Barnes-Hut non-determinism finding that led to the
/// integrator/force pairing rule (ADR-003) and the Yoshida-4 default,
/// but the scenario itself remains the cleanest high-N stress for
/// isolating integrator behaviour from scenario stiffness. Geometry
/// matches `solar_n641`'s radial range for apples-to-apples N=641;
/// per-ring phase offsets break radial alignment. In-regime targets
/// when run under IAS15 + direct: peak energy error at the f64 class
/// (~1e-11), zero degraded accepts, rejection rate under 10%.
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
        dt_budget: 0.1,
        // ~half the innermost orbital period: controller settles and
        // the energy_drift_slope fit gets a credible window.
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
