use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody};

/// Two equal-mass stars in a circular orbit.
///
/// Setup: m₁ = m₂ = 1.0, separation a = 1.0.
/// Each orbits the CoM at r = a/2 = 0.5 with
///   v = sqrt(G·m / a) = 1.0  (G = 1 units).
pub fn binary_star() -> Template {
    let v = (1.0_f64 / 1.0_f64).sqrt(); // sqrt(m/a) with G=1

    Template {
        name: "Binary Stars",
        bodies: vec![
            TemplateBody {
                mass: 1.0,
                radius: 0.045,
                position: Some([-0.5, 0.0]),
                velocity: [0.0, -v * 0.5],
                material: Material::Star,
            },
            TemplateBody {
                mass: 1.0,
                radius: 0.045,
                position: Some([0.5, 0.0]),
                velocity: [0.0, v * 0.5],
                material: Material::Star,
            },
        ],
        scale: 1.0,
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

    // CoM position = m2·a / m_total  from m1 (along +x)
    let r1 = m2 * a / m_total; // m1 distance from CoM
    let r2 = m1 * a / m_total; // m2 distance from CoM

    // Orbital velocity: v = 2π·r / T, using Kepler T² = 4π²a³/(G·m_total) → v = sqrt(G·m_total/a)·(r/a)
    let v_orb = (m_total / a).sqrt();
    let v1 = v_orb * r1 / a;
    let v2 = v_orb * r2 / a;

    Template {
        name: "Star + Companion",
        bodies: vec![
            TemplateBody {
                mass: m1,
                radius: 0.05,
                position: Some([-r1, 0.0]),
                velocity: [0.0, -v1],
                material: Material::Star,
            },
            TemplateBody {
                mass: m2,
                radius: 0.025,
                position: Some([r2, 0.0]),
                velocity: [0.0, v2],
                material: Material::BrownDwarf,
            },
        ],
        scale: 1.0,
    }
}
