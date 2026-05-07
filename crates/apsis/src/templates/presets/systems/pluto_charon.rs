//! Pluto–Charon binary configuration.
//!
//! Two bodies on a circular orbit about their common barycentre with
//! the mass ratio of the real Pluto–Charon system (≈ 8.7 : 1). With
//! that ratio the barycentre lies at ~10% of the orbital separation
//! from the primary's centre, outside the primary's bulk for any
//! plausible body radius — the canonical visual demonstration of
//! "barycentre outside the larger body".
//!
//! Units are dimensionless: only the mass ratio and the resulting
//! barycentre offset are visually meaningful here, so an absolute
//! scale would be misleading.

use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn pluto_charon(_seed: u64) -> Template {
    let m_pluto = 1.0_f64;
    let m_charon = 0.115_f64;
    let a = 1.0_f64;

    let m_total = m_pluto + m_charon;
    let r_pluto = a * m_charon / m_total;
    let r_charon = a * m_pluto / m_total;

    let v_relative = (m_total / a).sqrt();
    let v_pluto = v_relative * m_charon / m_total;
    let v_charon = v_relative * m_pluto / m_total;

    Template {
        name: "Pluto–Charon",
        description: "Dwarf-planet binary with mass ratio ~8.7:1; the barycentre lies outside \
                      Pluto, so both bodies trace circles around an empty point in space.",
        bodies: vec![
            TemplateBody {
                name: Some("Pluto"),
                mass: m_pluto,
                position: Some([-r_pluto, 0.0, 0.0]),
                velocity: [0.0, -v_pluto, 0.0],
                class_override: None,
                preset: &body_preset::ICY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("Charon"),
                mass: m_charon,
                position: Some([r_charon, 0.0, 0.0]),
                velocity: [0.0, v_charon, 0.0],
                class_override: None,
                preset: &body_preset::ICY,
                density: None,
                albedo: None,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}
