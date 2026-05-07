//! Alpha Centauri AB binary system.
//!
//! The two main components of the closest stellar system to the Sun,
//! both Sun-like (G2V and K1V), on a 79.91-year mutual orbit at high
//! eccentricity (e ≈ 0.52, semi-major axis ≈ 23.4 AU).
//!
//! Proxima Centauri orbits the AB pair at ~8 700 AU with a period
//! of ~547 000 years; including it in an interactive simulator
//! produces motion indistinguishable from escape on any reasonable
//! timescale, so this preset omits it.

use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem, builders::KG_M3_TO_SOLAR_AU3},
};

pub fn alpha_centauri_ab(_seed: u64) -> Template {
    let m_a = 1.10_f64;
    let m_b = 0.91_f64;
    let a = 23.4_f64;
    let e = 0.52_f64;

    let m_total = m_a + m_b;
    let r_peri = a * (1.0 - e);
    let v_peri = (m_total * (1.0 + e) / r_peri).sqrt();

    let r1 = r_peri * m_b / m_total;
    let r2 = r_peri * m_a / m_total;
    let v1 = v_peri * m_b / m_total;
    let v2 = v_peri * m_a / m_total;

    Template {
        name: "Alpha Centauri AB",
        description: "The closest binary star system to the Sun: two Sun-like stars on a \
                      79.91-year mutual orbit at eccentricity 0.52. Proxima omitted (period \
                      too long for interactive playback).",
        bodies: vec![
            // Densities from Bigot et al. (2006) interferometric radii
            // combined with the dynamical masses: ρ_A ≈ 1450 kg/m³,
            // ρ_B ≈ 1900 kg/m³. Both denser than the Sun (1408 kg/m³)
            // because they're slightly smaller relative to mass.
            TemplateBody {
                name: Some("Alpha Centauri A"),
                mass: m_a,
                position: Some([-r1, 0.0, 0.0]),
                velocity: [0.0, -v1, 0.0],
                preset: &body_preset::STAR,
                density: Some(1450.0 * KG_M3_TO_SOLAR_AU3),
            },
            TemplateBody {
                name: Some("Alpha Centauri B"),
                mass: m_b,
                position: Some([r2, 0.0, 0.0]),
                velocity: [0.0, v2, 0.0],
                preset: &body_preset::STAR,
                density: Some(1900.0 * KG_M3_TO_SOLAR_AU3),
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.002),
        units: UnitSystem::solar_au(),
    }
}
