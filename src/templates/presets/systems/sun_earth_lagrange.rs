use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn sun_earth_lagrange() -> Template {
    let m_sun = 1.0;
    let m_earth = 3.0e-6;

    let a = 1.0; // AU

    // Earth orbit around the Sun
    let (earth_pos, earth_vel) = circular_orbit(m_sun, a, 0.0);

    // L4 (+60° ahead of Earth)
    let (l4_pos, l4_vel) = circular_orbit(m_sun, a, std::f64::consts::PI / 3.0);

    // L5 (-60° behind Earth)
    let (l5_pos, l5_vel) = circular_orbit(m_sun, a, -std::f64::consts::PI / 3.0);

    Template {
        name: "Sun–Earth L4/L5",
        description: "Sun–Earth system with Trojan test particles at the L4 and L5 Lagrange points.",
        bodies: vec![
            // Sun
            TemplateBody {
                name: Some("Sun"),
                mass: m_sun,
                position: Some([0.0, 0.0]),
                velocity: [0.0, 0.0],
                material: Material::Star,
                spin: 0.0,
            },
            // Earth
            TemplateBody {
                name: Some("Earth"),
                mass: m_earth,
                position: Some(earth_pos),
                velocity: earth_vel,
                material: Material::Rocky,
                spin: 0.0,
            },
            // Trojan at L4
            TemplateBody {
                name: Some("L4 Probe"),
                mass: 1e-12,
                position: Some(l4_pos),
                velocity: l4_vel,
                material: Material::Asteroid,
                spin: 0.0,
            },
            // Trojan at L5
            TemplateBody {
                name: Some("L5 Probe"),
                mass: 1e-12,
                position: Some(l5_pos),
                velocity: l5_vel,
                material: Material::Asteroid,
                spin: 0.0,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.002),
    }
}
