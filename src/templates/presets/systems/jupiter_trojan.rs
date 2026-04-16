use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit},
};

pub fn jupiter_trojans() -> Template {
    use rand::random;
    use std::f64::consts::{PI, TAU};

    let m_sun = 1.0;
    let m_jupiter = 9.5e-4;

    let a = 5.204;

    // Jupiter orbit around the Sun
    let (j_pos, j_vel) = circular_orbit(m_sun, a, 0.0);

    let l4_angle = PI / 3.0;
    let l5_angle = -PI / 3.0;

    let mut bodies = vec![
        // Sun
        TemplateBody {
            name: Some("Sun"),
            mass: m_sun,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Star,
            spin: 0.0,
        },
        // Jupiter
        TemplateBody {
            name: Some("Jupiter"),
            mass: m_jupiter,
            position: Some(j_pos),
            velocity: j_vel,
            material: Material::Gas,
            spin: 0.0,
        },
    ];

    // ── Trojan clouds (L4 and L5) ── //
    let n = 400;
    let spread = 0.15; // angular spread around Lagrange points
    let vel_disp = 0.08; // velocity dispersion

    for &center in &[l4_angle, l5_angle] {
        for _ in 0..n {
            let angle = center + (random::<f64>() - 0.5) * spread;
            let r = a * (1.0 + (random::<f64>() - 0.5) * 0.05);

            let (pos, mut vel) = circular_orbit(m_sun, r, angle);

            // Small velocity dispersion to keep the cloud dynamically alive
            vel[0] *= 1.0 + vel_disp * (random::<f64>() - 0.5);
            vel[1] *= 1.0 + vel_disp * (random::<f64>() - 0.5);

            bodies.push(TemplateBody {
                name: None,
                mass: 1e-12,
                position: Some(pos),
                velocity: vel,
                material: Material::Asteroid,
                spin: 0.0,
            });
        }
    }

    Template {
        name: "Jupiter Trojans",
        description: "Sun–Jupiter system with Trojan asteroid clouds at L4 and L5.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.002),
        units: UnitSystem::solar_au(),
    }
}
