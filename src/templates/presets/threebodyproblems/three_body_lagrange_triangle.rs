use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody},
};

pub fn three_body_lagrange_triangle() -> Template {
    let m: f64 = 1.0;
    let r: f64 = 1.0;
    // Correct three-body equilateral formula: net centripetal force from two
    // bodies at side-length d = r·√3 gives v² = m / (r·√3).
    // (Using v = sqrt(m/r) — the two-body Kepler formula — is 32% too high
    // and causes the triangle to immediately fly apart.)
    let v = (m / (r * 3.0_f64.sqrt())).sqrt();

    let p1 = [r, 0.0];
    let p2 = [-0.5 * r, (3.0_f64).sqrt() / 2.0 * r];
    let p3 = [-0.5 * r, -(3.0_f64).sqrt() / 2.0 * r];

    fn perp([x, y]: [f64; 2]) -> [f64; 2] {
        [-y, x]
    }

    fn norm(v: [f64; 2]) -> [f64; 2] {
        let len = (v[0] * v[0] + v[1] * v[1]).sqrt();
        [v[0] / len, v[1] / len]
    }

    fn scale(v: [f64; 2], s: f64) -> [f64; 2] {
        [v[0] * s, v[1] * s]
    }

    let bodies = vec![
        TemplateBody {
            mass: m,
            radius: 0.006,
            position: Some(p1),
            velocity: scale(norm(perp(p1)), v),
            material: Material::Rocky,
        },
        TemplateBody {
            mass: m,
            radius: 0.006,
            position: Some(p2),
            velocity: scale(norm(perp(p2)), v),
            material: Material::Rocky,
        },
        TemplateBody {
            mass: m,
            radius: 0.006,
            position: Some(p3),
            velocity: scale(norm(perp(p3)), v),
            material: Material::Rocky,
        },
    ];

    Template {
        name: "Three Body - Lagrange (Newtonian)",
        bodies,
        scale: 1.0,
    }
}
