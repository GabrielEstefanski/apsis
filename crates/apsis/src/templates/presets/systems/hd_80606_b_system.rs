use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem},
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub fn hd_80606(seed: u64) -> Template {
    let mut bodies = Vec::with_capacity(2);

    let m_star = 1.0;
    let m_planet = 4.0e-3;

    let a = 0.455;
    let e = 0.93;

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("HD 80606"),
        mass: m_star,
        preset: &body_preset::STAR,
        density: None,
        position: Some([0.0, 0.0, 0.0]),
        velocity: [0.0, 0.0, 0.0],
    });

    // ── Planet (placed at periapsis) ───── //
    let r_peri = a * (1.0 - e);

    let v_peri = (m_star * (1.0 + e) / r_peri).sqrt();

    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };
    let phase = rng.random::<f64>() * std::f64::consts::TAU;

    let pos = [r_peri * phase.cos(), r_peri * phase.sin(), 0.0];
    let vel = [-v_peri * phase.sin(), v_peri * phase.cos(), 0.0];

    bodies.push(TemplateBody {
        name: Some("HD 80606 b"),
        mass: m_planet,
        preset: &body_preset::GAS,
        density: None,
        position: Some(pos),
        velocity: vel,
    });

    Template {
        name: "HD 80606",
        description: "Extreme eccentric exoplanet orbit.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0001),
        units: UnitSystem::solar_au(),
    }
}
