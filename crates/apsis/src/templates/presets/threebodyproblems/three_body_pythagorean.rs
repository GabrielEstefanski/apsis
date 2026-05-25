//! Burrau (1913) Pythagorean three-body problem.
//!
//! Three masses (3, 4, 5) at the vertices of a 3-4-5 right triangle,
//! released from rest. The opposite-side convention is preserved:
//! the side opposite mass `m_i` has length `m_i`. The centre of mass
//! coincides with the origin by construction.
//!
//! The system is chaotic; one body is typically ejected after a long
//! sequence of close encounters near `t ≈ 46` ($G = 1$, dimensionless).
//! Historically used as a benchmark for adaptive integrators.

use crate::{
    domain::body_preset,
    physics::integrator::IntegratorKind,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_pythagorean(_seed: u64) -> Template {
    Template {
        name: "Pythagorean (Burrau 1913)",
        description: "Masses 3, 4, 5 at the vertices of a 3-4-5 right triangle, released from \
                      rest. Canonical chaotic three-body problem; one body is eventually \
                      ejected after a complex sequence of close encounters.",
        bodies: vec![
            TemplateBody {
                name: Some("m = 3"),
                mass: 3.0,
                position: Some([1.0, 3.0, 0.0]),
                velocity: [0.0, 0.0, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("m = 4"),
                mass: 4.0,
                position: Some([-2.0, -1.0, 0.0]),
                velocity: [0.0, 0.0, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("m = 5"),
                mass: 5.0,
                position: Some([1.0, -1.0, 0.0]),
                velocity: [0.0, 0.0, 0.0],
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
        // Pythagorean has tight close encounters; Mercurius's hybrid
        // changeover handles them; Newton + fixed-step explicit blows up.
        suggested_integrator: Some(IntegratorKind::Mercurius),
        units: UnitSystem::dimensionless(),
    }
}
