use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn jupiter_trojans() -> Template {
    use rand::random;
    use std::f64::consts::TAU;

    let m_sun = 1.0;
    let m_jupiter = 9.5e-4;

    let a = 5.204;

    // Jupiter orbit
    let (j_pos, j_vel) = circular_orbit(m_sun, a, 0.0);

    let l4_angle = std::f64::consts::PI / 3.0;
    let l5_angle = -std::f64::consts::PI / 3.0;

    let mut bodies = vec![
        // Sun
        TemplateBody {
            mass: m_sun,
            radius: 0.04,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Star,
        },
        // Jupiter
        TemplateBody {
            mass: m_jupiter,
            radius: 0.01,
            position: Some(j_pos),
            velocity: j_vel,
            material: Material::Gas,
        },
    ];

    // ── Trojan cloud ── //
    let n = 400;
    let spread = 0.15; // maior que Terra → sistema mais “vivo”
    let vel_disp = 0.08;

    for &center in &[l4_angle, l5_angle] {
        for _ in 0..n {
            let angle = center + (random::<f64>() - 0.5) * spread;
            let r = a * (1.0 + (random::<f64>() - 0.5) * 0.05);

            let (mut pos, mut vel) = circular_orbit(m_sun, r, angle);

            // velocity dispersion
            vel[0] *= 1.0 + vel_disp * (random::<f64>() - 0.5);
            vel[1] *= 1.0 + vel_disp * (random::<f64>() - 0.5);

            bodies.push(TemplateBody {
                mass: 1e-12,
                radius: 0.0003,
                position: Some(pos),
                velocity: vel,
                material: Material::Asteroid,
            });
        }
    }

    Template {
        name: "Jupiter Trojans",
        bodies,
        scale: 1.0,
    }
}
