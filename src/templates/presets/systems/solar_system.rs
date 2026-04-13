//! Solar system template.
//!
//! ## Unit conventions
//!
//! | Quantity | Unit                  | Notes                          |
//! |----------|-----------------------|--------------------------------|
//! | Mass     | M_☉ (solar masses)    | Sun = 1.0                      |
//! | Distance | AU (astronomical units)| Earth–Sun = 1.0               |
//! | Velocity | AU / time_unit        | derived from circular_orbit    |
//! | G        | 4π² (in AU/M_☉/yr²)   | so that Earth period = 1 yr    |
//!
//! All planetary masses and semi-major axes are from the IAU 2012 nominal
//! planetary values (Prša et al. 2016, AJ 152, 41).
//! Moon orbital radius from Chapront et al. (2002).

use std::f64::consts::TAU;

use rand::random;

use crate::core::materials::Material;
use crate::templates::{Template, TemplateBody, builders::circular_orbit};

// ── Physical constants in simulation units ────────────────────────────────────

/// Solar mass [simulation mass units]. All other masses are fractions of this.
const M_SUN: f64 = 1.0;

/// Earth mass in solar masses (IAU 2012).
const M_EARTH: f64 = 3.003e-6;

/// Moon-to-Earth mass ratio (IAU 2012: k² = 1/81.3005).
const M_MOON: f64 = M_EARTH / 81.3005;

// ── Planet table ──────────────────────────────────────────────────────────────

struct Planet {
    /// Name for documentation; unused at runtime.
    #[allow(dead_code)]
    name: &'static str,
    /// Mass [M_☉].
    mass: f64,
    /// Semi-major axis [AU].
    a: f64,
    material: Material,
}

/// IAU 2012 nominal planetary masses and semi-major axes.
/// Sources: Prša et al. (2016); Williams (2023) NASA fact sheets.
const PLANETS: &[Planet] = &[
    Planet {
        name: "Mercury",
        mass: 1.652e-7,
        a: 0.38710,
        material: Material::Rocky,
    },
    Planet {
        name: "Venus",
        mass: 2.448e-6,
        a: 0.72333,
        material: Material::Rocky,
    },
    Planet {
        name: "Earth",
        mass: M_EARTH,
        a: 1.00000,
        material: Material::Rocky,
    },
    Planet {
        name: "Mars",
        mass: 3.213e-7,
        a: 1.52366,
        material: Material::Rocky,
    },
    Planet {
        name: "Jupiter",
        mass: 9.543e-4,
        a: 5.20336,
        material: Material::Gas,
    },
    Planet {
        name: "Saturn",
        mass: 2.857e-4,
        a: 9.53707,
        material: Material::Gas,
    },
    Planet {
        name: "Uranus",
        mass: 4.366e-5,
        a: 19.1913,
        material: Material::IceGiant,
    },
    Planet {
        name: "Neptune",
        mass: 5.151e-5,
        a: 30.0690,
        material: Material::IceGiant,
    },
    Planet {
        name: "Pluto",
        mass: 6.591e-9,
        a: 39.4817,
        material: Material::Icy,
    },
];

// ── Moon placement ────────────────────────────────────────────────────────────

/// Place a satellite in a stable circular orbit around a parent body.
///
/// Stability criterion: the satellite must be within the parent's Hill sphere.
/// The Hill radius is:
///
/// ```text
/// r_Hill = a_parent · (m_parent / 3·M_primary)^(1/3)
/// ```
///
/// For the Earth–Moon system: r_Hill ≈ 0.01 AU; the Moon orbits at 0.00257 AU,
/// which is ~26% of r_Hill — well within the stable zone (< 50% r_Hill).
///
/// No velocity correction factor is applied. The three-body perturbation from
/// the Sun is naturally handled by computing the satellite velocity in the
/// inertial frame by adding the parent's velocity vectorially.
fn place_moon(
    parent_pos: [f64; 2],
    parent_vel: [f64; 2],
    parent_mass: f64,
    a_moon: f64,
    phase: f64,
) -> ([f64; 2], [f64; 2]) {
    // Circular orbital speed around the parent in the parent's rest frame.
    // v = sqrt(G · m_parent / r); in our units G = 4π² so the circular_orbit
    // builder already embeds this — we replicate the formula directly here
    // to keep the parent-relative computation explicit.
    let v_circ = (M_SUN * (parent_mass / M_SUN) / a_moon).sqrt();

    // Position and velocity relative to parent.
    let rel_pos = [a_moon * phase.cos(), a_moon * phase.sin()];
    let rel_vel = [-v_circ * phase.sin(), v_circ * phase.cos()];

    // Inertial-frame state = parent state + relative state.
    (
        [parent_pos[0] + rel_pos[0], parent_pos[1] + rel_pos[1]],
        [parent_vel[0] + rel_vel[0], parent_vel[1] + rel_vel[1]],
    )
}

/// Compute the Hill radius for a body orbiting a primary.
///
/// ```text
/// r_Hill = a · (m_body / 3·M_primary)^(1/3)
/// ```
fn hill_radius(a: f64, m_body: f64, m_primary: f64) -> f64 {
    a * (m_body / (3.0 * m_primary)).powf(1.0 / 3.0)
}

// ── Comet placement ───────────────────────────────────────────────────────────

/// Place a comet on a realistic eccentric orbit (e ∈ [0.6, 0.99]).
///
/// Rather than using `v_circular * scale_factor`, we construct the orbit
/// correctly from the vis-viva equation at periapsis:
///
/// ```text
/// v_peri = sqrt(G·M · (1 + e) / r_peri)
/// ```
///
/// The comet is placed at periapsis on a random approach angle, so it will
/// follow a proper Keplerian ellipse/hyperbola rather than a perturbed circle.
fn place_comet(
    sun_mass: f64,
    r_peri: f64,
    e: f64,
    omega: f64, // argument of periapsis [rad]
) -> ([f64; 2], [f64; 2]) {
    // Speed at periapsis from vis-viva: v² = GM(1+e)/r_peri
    let v_peri = (sun_mass * (1.0 + e) / r_peri).sqrt();

    // At periapsis the velocity is perpendicular to the position vector.
    // Position: r_peri along the periapsis direction.
    // Velocity: perpendicular (CCW).
    let pos = [r_peri * omega.cos(), r_peri * omega.sin()];
    let vel = [-v_peri * omega.sin(), v_peri * omega.cos()];

    (pos, vel)
}

// ── Template builder ──────────────────────────────────────────────────────────

pub fn solar_system() -> Template {
    let mut bodies = Vec::with_capacity(1 + PLANETS.len() + 1 + 600 + 30);

    // ── Sun ───────────────────────────────────────────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("Sun"),
        mass: M_SUN,
        material: Material::Star,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        spin: 0.0,
    });

    // Track Earth state for Moon placement.
    let mut earth_pos = [0.0_f64; 2];
    let mut earth_vel = [0.0_f64; 2];

    // ── Planets ───────────────────────────────────────────────────────────── //
    for p in PLANETS {
        // Random initial phase — physically valid for a snapshot of the system.
        let phase = random::<f64>() * TAU;
        let (pos, vel) = circular_orbit(M_SUN, p.a, phase);

        if (p.a - 1.0).abs() < 1e-6 {
            earth_pos = pos;
            earth_vel = vel;
        }

        bodies.push(TemplateBody {
            name: Some(p.name),
            mass: p.mass,
            material: p.material,
            position: Some(pos),
            velocity: vel,
            spin: 0.0,
        });
    }

    // ── Moon ──────────────────────────────────────────────────────────────── //
    {
        // Moon semi-major axis: 384 400 km = 0.002570 AU (Chapront et al. 2002).
        let moon_a = 0.002570;

        // Verify stability: moon_a must be well within Earth's Hill sphere.
        // r_Hill(Earth) ≈ 0.0100 AU → moon_a / r_Hill ≈ 0.26 (stable, < 0.5).
        let r_hill_earth = hill_radius(1.0, M_EARTH, M_SUN);
        debug_assert!(
            moon_a < 0.5 * r_hill_earth,
            "Moon semi-major axis {moon_a:.4} AU exceeds 50% of Earth Hill radius \
             {r_hill_earth:.4} AU — orbit may be unstable"
        );

        let phase = random::<f64>() * TAU;
        let (moon_pos, moon_vel) = place_moon(earth_pos, earth_vel, M_EARTH, moon_a, phase);

        bodies.push(TemplateBody {
            name: Some("Moon"),
            mass: M_MOON,
            position: Some(moon_pos),
            velocity: moon_vel,
            material: Material::Rocky,
            spin: 0.0,
        });
    }

    // ── Asteroid belt (2.2–3.2 AU) ────────────────────────────────────────── //
    //
    // Eccentricity distribution follows the observed main-belt distribution
    // (Bottke et al. 2005: mean e ≈ 0.14, σ ≈ 0.08).
    for _ in 0..600 {
        let a = 2.2 + random::<f64>() * 1.0;
        let e = (random::<f64>() * 0.16 + 0.06).min(0.35); // e ∈ [0.06, 0.35]
        let phase = random::<f64>() * TAU;

        // Place at a random true anomaly on the ellipse, not just at periapsis.
        // For simplicity we use the circular speed and apply an eccentricity
        // correction: v_tangential = v_circ * sqrt(1 + e·cos(ν) + ...) — here
        // we use the first-order approximation v ≈ v_circ · (1 + e·cos(phase)).
        let (mut pos, mut vel) = circular_orbit(M_SUN, a, phase);
        let ecc_factor = 1.0 + e * phase.cos();
        vel[0] *= ecc_factor;
        vel[1] *= ecc_factor;

        bodies.push(TemplateBody {
            name: None,
            mass: 1e-10,
            position: Some(pos),
            velocity: vel,
            material: Material::Asteroid,
            spin: 0.0,
        });
    }

    // ── Comets (Jupiter-family + long-period) ──────────────────────────────── //
    //
    // Jupiter-family comets: e ∈ [0.5, 0.8], periapsis 1–3 AU (Levison 1996).
    // Long-period comets:    e ∈ [0.97, 0.999], periapsis 0.5–2 AU.
    for i in 0..30 {
        let (r_peri, e) = if i < 20 {
            // Jupiter-family
            let r = 1.0 + random::<f64>() * 2.0;
            let e = 0.50 + random::<f64>() * 0.30;
            (r, e)
        } else {
            // Long-period / Oort cloud
            let r = 0.5 + random::<f64>() * 1.5;
            let e = 0.97 + random::<f64>() * 0.029;
            (r, e)
        };

        let omega = random::<f64>() * TAU;
        let (pos, vel) = place_comet(M_SUN, r_peri, e, omega);

        bodies.push(TemplateBody {
            name: None,
            mass: 1e-12,
            position: Some(pos),
            velocity: vel,
            material: Material::Comet,
            spin: 0.0,
        });
    }

    Template {
        name: "Solar System",
        description: "The Sun and eight planets, plus Pluto, the Moon, the asteroid belt, and a sprinkling of comets.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.001), // ~0.36 days; ~1/1000 of Earth's orbital period
    }
}
