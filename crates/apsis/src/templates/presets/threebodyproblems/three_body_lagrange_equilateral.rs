//! Lagrange's equilateral three-body solution (equal-mass case).
//!
//! Three equal masses at the vertices of an equilateral triangle,
//! rotating rigidly about the common centre of mass. For three equal
//! masses at distance `r` from the COM, the angular velocity that
//! balances gravity is `ω = √(3·G·m / r³)` and each body's velocity
//! is tangential with magnitude `ω·r`.
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
    let omega = (3.0 * m / (r * r * r)).sqrt();
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
                preset: &body_preset::ROCKY,
            },
            TemplateBody {
                name: Some("Body 2"),
                mass: m,
                position: Some(p2),
                velocity: tangential(p2, v),
                preset: &body_preset::ROCKY,
            },
            TemplateBody {
                name: Some("Body 3"),
                mass: m,
                position: Some(p3),
                velocity: tangential(p3, v),
                preset: &body_preset::ROCKY,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}

fn tangential(position: [f64; 3], speed: f64) -> [f64; 3] {
    let r = (position[0] * position[0] + position[1] * position[1]).sqrt();
    [-position[1] * speed / r, position[0] * speed / r, 0.0]
}
