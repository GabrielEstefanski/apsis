use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn sun_earth_lagrange() -> Template {
    let m_sun = 1.0;
    let m_earth = 3.0e-6;

    let a = 1.0; // AU

    // Earth orbit
    let (earth_pos, earth_vel) = circular_orbit(m_sun, a, 0.0);

    // L4 (+60°)
    let (l4_pos, l4_vel) = circular_orbit(m_sun, a, std::f64::consts::PI / 3.0);

    // L5 (-60°)
    let (l5_pos, l5_vel) = circular_orbit(m_sun, a, -std::f64::consts::PI / 3.0);

    Template {
        name: "Sun–Earth L4/L5",
        bodies: vec![
            // Sun
            TemplateBody {
                mass: m_sun,
                radius: 0.02,
                position: Some([0.0, 0.0]),
                velocity: [0.0, 0.0],
                material: Material::Star,
            },
            // Earth
            TemplateBody {
                mass: m_earth,
                radius: 0.002,
                position: Some(earth_pos),
                velocity: earth_vel,
                material: Material::Rocky,
            },
            // Trojan L4
            TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(l4_pos),
                velocity: l4_vel,
                material: Material::Asteroid,
            },
            // Trojan L5
            TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(l5_pos),
                velocity: l5_vel,
                material: Material::Asteroid,
            },
        ],
        scale: 1.0,
    }
}
