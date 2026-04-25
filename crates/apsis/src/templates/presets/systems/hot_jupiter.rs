//! Hot Jupiter — Sun-like star with a Jupiter-mass companion on a
//! short-period circular orbit.
//!
//! Reference parameters approximate observed hot Jupiters: a Jovian
//! companion at ~0.05 AU around a Sun-like primary, with an orbital
//! period of a few days. The post-migration state is integrated
//! directly; this preset is a contrast to the standard solar-system
//! cadence rather than a study of migration itself.

use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit},
};

pub fn hot_jupiter(_seed: u64) -> Template {
    let m_star = 1.0;
    let m_jupiter = 9.55e-4;
    let a_jupiter = 0.05;

    let (jupiter_pos, jupiter_vel) = circular_orbit(m_star, a_jupiter, 0.0);

    Template {
        name: "Hot Jupiter",
        description: "A Sun-like star with a Jupiter-mass companion on a close (~0.05 AU) \
                      circular orbit, period ~4 days; canonical post-migration configuration.",
        bodies: vec![
            TemplateBody {
                name: Some("Star"),
                mass: m_star,
                position: Some([0.0, 0.0]),
                velocity: [0.0, 0.0],
                material: Material::Star,
            },
            TemplateBody {
                name: Some("Hot Jupiter"),
                mass: m_jupiter,
                position: Some(jupiter_pos),
                velocity: jupiter_vel,
                material: Material::Gas,
            },
        ],
        display_scale: 50.0,
        suggested_dt: Some(0.0001),
        units: UnitSystem::solar_au(),
    }
}
