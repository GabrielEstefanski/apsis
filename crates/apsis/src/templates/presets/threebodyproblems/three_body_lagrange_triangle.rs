use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_lagrange_triangle(_seed: u64) -> Template {
    let m: f64 = 1.0;
    let r: f64 = 1.0;

    // Para triângulo equilátero estável:
    // ω² = (G * m_total) / R³  com R = distância ao CoM = r
    // m_total = 3m  →  ω = sqrt(3m / r³)
    let omega = (3.0 * m / (r * r * r)).sqrt();
    let v = omega * r;

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
            name: Some("Body 1"),
            mass: m,
            position: Some(p1),
            velocity: scale(norm(perp(p1)), v),
            material: Material::Rocky,
        },
        TemplateBody {
            name: Some("Body 2"),
            mass: m,
            position: Some(p2),
            velocity: scale(norm(perp(p2)), v),
            material: Material::Rocky,
        },
        TemplateBody {
            name: Some("Body 3"),
            mass: m,
            position: Some(p3),
            velocity: scale(norm(perp(p3)), v),
            material: Material::Rocky,
        },
    ];

    Template {
        name: "Three Body - Lagrange (Newtonian)",
        description: "Three equal-mass bodies in a stable rotating equilateral triangle configuration.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}
