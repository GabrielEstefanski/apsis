use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody, UnitSystem};

pub fn binary_star() -> Template {
    let a = 1.0_f64;
    let m = 1.0_f64;

    let v = (m / a).sqrt();

    let r = 0.5;

    let v_body = v * 0.5;

    let omega = v_body / r;

    Template {
        name: "Binary Stars",
        description: "Two equal-mass stars in a circular orbit.",
        bodies: vec![
            TemplateBody {
                name: Some("Star A"),
                mass: m,
                position: Some([-r, 0.0]),
                velocity: [0.0, -v_body],
                material: Material::Star,
                spin: omega,
            },
            TemplateBody {
                name: Some("Star B"),
                mass: m,
                position: Some([r, 0.0]),
                velocity: [0.0, v_body],
                material: Material::Star,
                spin: omega,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::solar_au(),
    }
}

/// Unequal binary: a solar-mass primary with a sub-stellar companion (0.08 M_sun).
///
/// Both orbit the shared CoM. The companion's orbit is computed via reduced-mass
/// Kepler:  v_2 = sqrt(G·(m₁ + m₂) · m₁ / (m₁ + m₂) / a)
///               = sqrt(G · m₁ / a)  in mass-ratio approximation for m₁ ≫ m₂.
pub fn star_companion() -> Template {
    let m1 = 1.0_f64;
    let m2 = 0.08_f64;
    let a = 1.0_f64;

    let m_total = m1 + m2;

    // distâncias ao CoM
    let r1 = m2 * a / m_total;
    let r2 = m1 * a / m_total;

    // velocidade orbital base
    let v_orb = (m_total / a).sqrt();

    let v1 = v_orb * r1 / a;
    let v2 = v_orb * r2 / a;

    // rotações sincronizadas (tidal lock aproximado)
    let omega = v_orb / a;

    Template {
        name: "Star + Companion",
        description: "A solar-mass star with a brown dwarf companion.",
        bodies: vec![
            TemplateBody {
                name: Some("Primary Star"),
                mass: m1,
                position: Some([-r1, 0.0]),
                velocity: [0.0, -v1],
                material: Material::Star,
                spin: omega,
            },
            TemplateBody {
                name: Some("Companion"),
                mass: m2,
                position: Some([r2, 0.0]),
                velocity: [0.0, v2],
                material: Material::BrownDwarf,
                spin: omega,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::solar_au(),
    }
}
