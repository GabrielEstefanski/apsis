use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn hd_80606() -> Template {
    let mut bodies = Vec::with_capacity(2);

    let m_star = 1.0;
    let m_planet = 4.0e-3;

    let a = 0.455;
    let e = 0.93;

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("HD 80606"),
        mass: m_star,
        material: Material::Star,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        spin: 0.0,
    });

    // ── Planet (placed at periapsis) ───── //
    let r_peri = a * (1.0 - e);

    let v_peri = (m_star * (1.0 + e) / r_peri).sqrt();

    let phase = rand::random::<f64>() * std::f64::consts::TAU;

    let pos = [r_peri * phase.cos(), r_peri * phase.sin()];
    let vel = [-v_peri * phase.sin(), v_peri * phase.cos()];

    bodies.push(TemplateBody {
        name: Some("HD 80606 b"),
        mass: m_planet,
        material: Material::Gas,
        position: Some(pos),
        velocity: vel,
        spin: 0.0,
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
