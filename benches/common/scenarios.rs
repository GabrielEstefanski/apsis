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
    }
}

/// Ordered list of all scenarios. The order determines the order of
/// sections in the baseline file.
pub fn all() -> Vec<ScenarioSpec> {
    vec![kepler_e05(), kepler_e09(), kepler_e099(), pythagorean()]
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
    }
}
