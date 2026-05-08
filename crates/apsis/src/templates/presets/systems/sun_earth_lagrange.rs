use crate::{
    domain::body_preset,
    templates::{
        Template, TemplateBody, UnitSystem,
        builders::{KG_M3_TO_SOLAR_AU3, circular_orbit},
    },
};

pub fn sun_earth_lagrange(_seed: u64) -> Template {
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
                position: Some([0.0, 0.0, 0.0]),
                velocity: [0.0, 0.0, 0.0],
                class_override: None,
                preset: &body_preset::STAR,
                density: Some(1408.0 * KG_M3_TO_SOLAR_AU3),
                albedo: None,
            },
            // Earth
            TemplateBody {
                name: Some("Earth"),
                mass: m_earth,
                position: Some(earth_pos),
                velocity: earth_vel,
                class_override: None,
                preset: &body_preset::ROCKY,
                density: Some(5514.0 * KG_M3_TO_SOLAR_AU3),
                albedo: None,
            },
            // Trojan at L4
            TemplateBody {
                name: Some("L4 Probe"),
                mass: 1e-12,
                position: Some(l4_pos),
                velocity: l4_vel,
                class_override: None,
                preset: &body_preset::ASTEROID,
                density: None,
                albedo: None,
            },
            // Trojan at L5
            TemplateBody {
                name: Some("L5 Probe"),
                mass: 1e-12,
                position: Some(l5_pos),
                velocity: l5_vel,
                class_override: None,
                preset: &body_preset::ASTEROID,
                density: None,
                albedo: None,
            },
        ],
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        suggested_dt: Some(0.002),
        units: UnitSystem::solar_au(),
    }
}
