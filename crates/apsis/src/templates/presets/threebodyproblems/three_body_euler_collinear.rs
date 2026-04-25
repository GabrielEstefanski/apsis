//! Euler's collinear three-body solution (1767).
//!
//! Three bodies on a straight line through the centre of mass, rotating
//! rigidly about that axis. For three equal masses at `-d`, `0`, `+d`
//! the angular velocity that balances gravity on each outer body is
//!
//!     ω² = G·m / d³ + G·m / (4d³) = 5·G·m / (4·d³)
//!
//! and the outer-body tangential speed is `v = ω·d = √(5·G·m / (4·d))`.
//! The middle body sits on the rotation axis with zero velocity.
//!
//! Linearly unstable for any mass ratio (Euler's collinear configurations
//! are saddle points of the rotating-frame potential). The configuration
//! disintegrates within a few orbits — pedagogical contrast to the
//! Lagrange equilateral solution, which is stable for sufficiently
//! asymmetric masses.

use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_euler_collinear(_seed: u64) -> Template {
    let m = 1.0_f64;
    let d = 1.0_f64;
    let v = (5.0 * m / (4.0 * d)).sqrt();

    Template {
        name: "Euler Collinear (unstable)",
        description: "Three equal-mass bodies on a rotating straight line — Euler's 1767 \
                      collinear solution. Linearly unstable for any mass ratio; the \
                      configuration disintegrates within a few orbits.",
        bodies: vec![
            TemplateBody {
                name: Some("Outer body 1"),
                mass: m,
                position: Some([-d, 0.0]),
                velocity: [0.0, -v],
                material: Material::Rocky,
            },
            TemplateBody {
                name: Some("Centre body"),
                mass: m,
                position: Some([0.0, 0.0]),
                velocity: [0.0, 0.0],
                material: Material::Rocky,
            },
            TemplateBody {
                name: Some("Outer body 2"),
                mass: m,
                position: Some([d, 0.0]),
                velocity: [0.0, v],
                material: Material::Rocky,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}
