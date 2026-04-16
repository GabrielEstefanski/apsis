use crate::{
    core::materials::Material,
    templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit},
};

pub fn sun_earth_unstable_lagrange() -> Template {
    let m_sun = 1.0;
    let m_earth = 3.0e-6;

    let a = 1.0;

    // Earth orbit around the Sun
    let (earth_pos, earth_vel) = circular_orbit(m_sun, a, 0.0);

    // ── L1 / L2 distance approximation ── //
    let mu = m_earth / (m_sun + m_earth);
    let r_l = a * (mu / 3.0).cbrt(); // ≈ 0.01 AU

    // L1 (between Sun and Earth)
    let l1_pos = [earth_pos[0] - r_l, earth_pos[1]];
    let l1_vel = earth_vel;

    // L2 (beyond Earth)
    let l2_pos = [earth_pos[0] + r_l, earth_pos[1]];
    let l2_vel = earth_vel;

    // L3 (opposite side of the Sun)
    let (l3_pos, l3_vel) = circular_orbit(m_sun, a, std::f64::consts::PI);

    Template {
        name: "Sun–Earth L1/L2/L3 (Unstable)",
        description: "Sun–Earth system with test particles at unstable Lagrange points (L1, L2, L3).",
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
            // L1 (unstable)
            TemplateBody {
                name: Some("L1 Probe"),
                mass: 1e-12,
                position: Some(l1_pos),
                velocity: l1_vel,
                material: Material::Asteroid,
                spin: 0.0,
            },
            // L2 (unstable)
            TemplateBody {
                name: Some("L2 Probe"),
                mass: 1e-12,
                position: Some(l2_pos),
                velocity: l2_vel,
                material: Material::Asteroid,
                spin: 0.0,
            },
            // L3 (unstable)
            TemplateBody {
                name: Some("L3 Probe"),
                mass: 1e-12,
                position: Some(l3_pos),
                velocity: l3_vel,
                material: Material::Asteroid,
                spin: 0.0,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.002),
        units: UnitSystem::solar_au(),
    }
}
