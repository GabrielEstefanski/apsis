use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn sun_earth_unstable_lagrange() -> Template {
    let m_sun = 1.0;
    let m_earth = 3.0e-6;

    let a = 1.0;

    let (earth_pos, earth_vel) = circular_orbit(m_sun, a, 0.0);

    // ── L1 / L2 distance approximation ── //
    let mu = m_earth / (m_sun + m_earth);
    let r_l = a * (mu / 3.0).cbrt(); // ~0.01

    // L1 (between Sun and Earth)
    let l1_pos = [earth_pos[0] - r_l, earth_pos[1]];
    let l1_vel = earth_vel;

    // L2 (beyond Earth)
    let l2_pos = [earth_pos[0] + r_l, earth_pos[1]];
    let l2_vel = earth_vel;

    // L3 (opposite side of Sun)
    let (l3_pos, l3_vel) = circular_orbit(m_sun, a, std::f64::consts::PI);

    Template {
        name: "Sun–Earth L1/L2/L3 (Unstable)",
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
            // L1
            TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(l1_pos),
                velocity: l1_vel,
                material: Material::Asteroid,
            },
            // L2
            TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(l2_pos),
                velocity: l2_vel,
                material: Material::Asteroid,
            },
            // L3
            TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(l3_pos),
                velocity: l3_vel,
                material: Material::Asteroid,
            },
        ],
        scale: 1.0,
    }
}
