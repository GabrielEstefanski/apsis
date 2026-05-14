//! Lagrange's equilateral three-body solution (equal-mass case).
//!
//! Three equal masses at the vertices of an equilateral triangle,
//! rotating rigidly about the common centre of mass. With `r` denoting
//! the distance from each vertex to the COM (side length `L = r·√3`),
//! the gravitational force on each body sums to `G·m²·√3 / L²` toward
//! the COM. Balancing centripetal `m·ω²·r` against this gives
//! `ω² = G·m / (r³·√3)` and tangential speed `v = ω·r`.
//!
//! Lagrange's equilateral solution is *only* linearly stable when one
//! mass strongly dominates the other two — Routh's criterion requires
//!
//! ```text
//! 27·(m₁·m₂ + m₂·m₃ + m₃·m₁) < (m₁ + m₂ + m₃)²
//! ```
//!
//! which fails for equal masses. This preset deliberately exhibits the
//! unstable case: rounding-noise perturbations grow exponentially and
//! the configuration disintegrates within a few orbits. It is the
//! pedagogical counterpart to the stable Sun–Earth L4/L5 trojan preset
//! (where `m_test ≪ m_sun` satisfies Routh).

use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_lagrange_equilateral(_seed: u64) -> Template {
    let m = 1.0_f64;
    let r = 1.0_f64;
    let omega = (m / (r * r * r * 3.0_f64.sqrt())).sqrt();
    let v = omega * r;

    let p1 = [r, 0.0, 0.0];
    let p2 = [-0.5 * r, (3.0_f64).sqrt() / 2.0 * r, 0.0];
    let p3 = [-0.5 * r, -(3.0_f64).sqrt() / 2.0 * r, 0.0];

    Template {
        name: "Lagrange Equilateral (unstable)",
        description: "Three equal-mass bodies at the vertices of an equilateral triangle. \
                      Lagrange's analytic equilibrium; linearly unstable for equal masses \
                      (Routh's criterion violated), so the configuration drifts apart visibly \
                      within a few orbits.",
        bodies: vec![
            TemplateBody {
                name: Some("Body 1"),
                mass: m,
                position: Some(p1),
                velocity: tangential(p1, v),
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("Body 2"),
                mass: m,
                position: Some(p2),
                velocity: tangential(p2, v),
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("Body 3"),
                mass: m,
                position: Some(p3),
                velocity: tangential(p3, v),
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
        ],
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        suggested_dt: Some(0.001),
        suggested_integrator: None,
        units: UnitSystem::dimensionless(),
    }
}

fn tangential(position: [f64; 3], speed: f64) -> [f64; 3] {
    let r = (position[0] * position[0] + position[1] * position[1]).sqrt();
    [-position[1] * speed / r, position[0] * speed / r, 0.0]
}
