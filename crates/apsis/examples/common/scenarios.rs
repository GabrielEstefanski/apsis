//! Scenario catalogue shared across benchmark examples.
//!
//! Each scenario is a triple: a deterministic body-list builder, a
//! natural integration timestep, and a characteristic dynamical
//! timescale. Consumers (speed benchmarks, conservation benchmarks,
//! integrator-Pareto sweeps, …) compose these into cell-by-cell
//! measurements without needing to duplicate scenario definitions.
//!
//! ## Why the catalogue lives here, not in `apsis`
//!
//! The scenario definitions are *benchmark scaffolding* — they pick
//! specific masses, distributions, and timestep conventions that are
//! neither public API nor physical primitives the library should
//! expose. Extracting them to a shared example-side module keeps them
//! out of the public API of `apsis` while still avoiding per-example
//! duplication.

use std::f64::consts::TAU;

use apsis::domain::body::Body;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

/// A self-contained scenario specification.
///
/// Fields are plain data plus three function pointers:
/// [`build`](Scenario::build) materialises the initial body list,
/// [`dt_hint`](Scenario::dt_hint) returns the natural timestep at a
/// given N (chosen per scenario to avoid the hidden bias of a
/// one-size-fits-all `dt`), and [`t_characteristic`](Scenario::t_characteristic)
/// returns the dynamical timescale used by long-integration
/// benchmarks to pick a window over which to measure drift.
pub struct Scenario {
    pub name: &'static str,
    pub description: &'static str,
    pub build: fn(usize, u64) -> Vec<Body>,
    pub dt_hint: fn(usize) -> f64,
    pub t_characteristic: fn(usize) -> f64,
    pub min_n: usize,
    pub max_n: usize,
}

/// The full scenario catalogue.
///
/// Benchmarks that only care about a subset should filter by
/// [`Scenario::name`] or intersect with [`Scenario::min_n`] /
/// [`Scenario::max_n`]. New scenarios are added by pushing to this
/// list in one place; every consumer picks them up automatically.
pub fn all() -> Vec<Scenario> {
    vec![
        Scenario {
            name: "friendly_cluster",
            description: "Uniform 2D disk, virial-scale velocities",
            build: build_friendly_cluster,
            dt_hint: dt_hint_cluster,
            t_characteristic: t_char_cluster,
            min_n: 128,
            max_n: 65_536,
        },
        Scenario {
            name: "hierarchical_kepler",
            description: "1 M☉ + (N−1) test bodies on log-uniform Kepler orbits",
            build: build_hierarchical_kepler,
            dt_hint: dt_hint_hierarchical,
            t_characteristic: t_char_hierarchical,
            min_n: 128,
            max_n: 65_536,
        },
        Scenario {
            name: "clustered_substructure",
            description: "Main Plummer-like cluster + sub-clumps",
            build: build_clustered_substructure,
            dt_hint: dt_hint_cluster,
            t_characteristic: t_char_cluster,
            min_n: 128,
            max_n: 32_768,
        },
        Scenario {
            name: "multiple_binaries",
            description: "N/2 isolated two-body systems, far-separated",
            build: build_multiple_binaries,
            dt_hint: dt_hint_binaries,
            t_characteristic: t_char_binaries,
            min_n: 128,
            max_n: 32_768,
        },
    ]
}

// ── dt_hint functions ────────────────────────────────────────────────────── //

/// Natural timestep for a uniform-disk cluster of radius 10, total mass
/// scaling as M ∝ N (per [`build_friendly_cluster`]). Picked as
/// `t_dyn / 1000`: ~1000 steps per crossing time. A dense sub-clump scenario
/// with the same mass profile uses the same hint (tree depth and leaf
/// capacity dominate, not force magnitude).
pub fn dt_hint_cluster(n: usize) -> f64 {
    t_char_cluster(n) / 1000.0
}

/// Natural timestep for the hierarchical Kepler scenario: the innermost
/// orbit sits at `a = 0.3 AU` with period `T = 2π · 0.3^{3/2}`; we use
/// `T / 1000 ≈ 1e-3` as `dt`. Independent of N because the innermost
/// orbit does not scale with body count.
pub fn dt_hint_hierarchical(_n: usize) -> f64 {
    t_char_hierarchical(0) / 1000.0
}

/// Natural timestep for the multiple-binaries scenario: each binary has
/// period `T = 2π`; `dt = T / 1000`.
pub fn dt_hint_binaries(_n: usize) -> f64 {
    t_char_binaries(0) / 1000.0
}

// ── characteristic timescales ───────────────────────────────────────────── //

/// Dynamical time of the uniform-disk cluster: `t_dyn = √(R³ / M)`.
pub fn t_char_cluster(n: usize) -> f64 {
    let r_disk = 10.0_f64;
    let m_total = (n as f64) * 1e-4;
    (r_disk.powi(3) / m_total).sqrt()
}

/// Inner-orbit period of the hierarchical scenario: `T = 2π · a_min^{3/2}`
/// with `a_min = 0.3 AU` and `GM = 1`.
pub fn t_char_hierarchical(_n: usize) -> f64 {
    let a_min = 0.3_f64;
    TAU * a_min.powf(1.5)
}

/// Orbital period of each binary in the multiple-binaries scenario:
/// `T = 2π · √(a³ / M)` with `a = 1`, `M = 1` → `T = 2π`.
pub fn t_char_binaries(_n: usize) -> f64 {
    TAU
}

// ── Scenario builders ────────────────────────────────────────────────────── //

pub fn build_friendly_cluster(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let r_disk: f64 = 10.0;
    let m_total: f64 = (n as f64) * 1e-4;
    let m_each: f64 = m_total / n as f64;
    let v_scale = (m_total / r_disk).sqrt() * 0.2;

    (0..n)
        .map(|_| {
            let theta = rng.random::<f64>() * TAU;
            let r = r_disk * rng.random::<f64>().sqrt();
            let x = r * theta.cos();
            let y = r * theta.sin();
            let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_scale;
            let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_scale;
            Body::rocky(m_each).at(x, y).with_velocity(vx, vy)
        })
        .collect()
}

/// Sun at origin plus N−1 bodies on circular Kepler orbits with
/// log-uniform semi-major axes in [0.3, 30] AU.
pub fn build_hierarchical_kepler(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut bodies = Vec::with_capacity(n);
    bodies.push(Body::star(1.0).at(0.0, 0.0).unsoftened());

    let log_a_min = 0.3_f64.log10();
    let log_a_max = 30.0_f64.log10();
    let log_m_min = -6.0_f64;
    let log_m_max = -4.0_f64;

    for _ in 1..n {
        let a = 10f64.powf(log_a_min + (log_a_max - log_a_min) * rng.random::<f64>());
        let theta = rng.random::<f64>() * TAU;
        let x = a * theta.cos();
        let y = a * theta.sin();
        let v = (1.0 / a).sqrt();
        let vx = -v * theta.sin();
        let vy = v * theta.cos();
        let m = 10f64.powf(log_m_min + (log_m_max - log_m_min) * rng.random::<f64>());
        bodies.push(Body::rocky(m).at(x, y).with_velocity(vx, vy));
    }
    bodies
}

/// 80% of bodies form a diffuse main disk; 20% form `K_sub` dense
/// sub-clumps placed at random positions around the main centre.
pub fn build_clustered_substructure(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let r_main: f64 = 10.0;
    let r_sub: f64 = 1.0;
    let m_total: f64 = (n as f64) * 1e-4;
    let m_each: f64 = m_total / n as f64;
    let v_main = (m_total / r_main).sqrt() * 0.2;
    let v_sub = (m_total / r_sub).sqrt() * 0.2;

    let n_main = (n as f64 * 0.8) as usize;
    let n_sub_bodies = n - n_main;
    let k_sub = ((n_sub_bodies as f64).sqrt().ceil() as usize).max(2);
    let per_clump = n_sub_bodies / k_sub;

    let mut bodies = Vec::with_capacity(n);

    for _ in 0..n_main {
        let theta = rng.random::<f64>() * TAU;
        let r = r_main * rng.random::<f64>().sqrt();
        let (x, y) = (r * theta.cos(), r * theta.sin());
        let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_main;
        let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_main;
        bodies.push(Body::rocky(m_each).at(x, y).with_velocity(vx, vy));
    }

    for clump in 0..k_sub {
        let theta_c = rng.random::<f64>() * TAU;
        let r_c = r_main * 0.5 * (1.0 + rng.random::<f64>());
        let (cx, cy) = (r_c * theta_c.cos(), r_c * theta_c.sin());
        let members = if clump == k_sub - 1 { n_sub_bodies - per_clump * clump } else { per_clump };
        for _ in 0..members {
            let theta = rng.random::<f64>() * TAU;
            let r = r_sub * rng.random::<f64>().sqrt();
            let (x, y) = (cx + r * theta.cos(), cy + r * theta.sin());
            let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_sub;
            let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_sub;
            bodies.push(Body::rocky(m_each).at(x, y).with_velocity(vx, vy));
        }
    }
    bodies
}

/// `N/2` equal-mass two-body systems, each with separation 1 AU and
/// circular velocity, placed at random positions in a large box so that
/// inter-binary distances dominate intra-binary distances.
pub fn build_multiple_binaries(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let n_binaries = n / 2;
    let field_size = 100.0_f64;
    let sep = 1.0_f64;
    let m_each = 0.5_f64;
    let v_each = 0.5_f64;

    let mut bodies = Vec::with_capacity(n_binaries * 2);
    for _ in 0..n_binaries {
        let cx = (rng.random::<f64>() - 0.5) * field_size;
        let cy = (rng.random::<f64>() - 0.5) * field_size;
        let phi = rng.random::<f64>() * TAU;
        bodies.push(
            Body::rocky(m_each)
                .at(cx - 0.5 * sep * phi.cos(), cy - 0.5 * sep * phi.sin())
                .with_velocity(v_each * phi.sin(), -v_each * phi.cos()),
        );
        bodies.push(
            Body::rocky(m_each)
                .at(cx + 0.5 * sep * phi.cos(), cy + 0.5 * sep * phi.sin())
                .with_velocity(-v_each * phi.sin(), v_each * phi.cos()),
        );
    }
    bodies
}
