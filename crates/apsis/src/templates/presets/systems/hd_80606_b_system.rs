//! HD 80606 — single-star system with the most eccentric known
//! transiting hot Jupiter.
//!
//! HD 80606 b is a 4 M_Jup gas giant on an `e = 0.93` orbit. Periapsis
//! takes the planet to 0.03 AU from the star (closer than Mercury);
//! apoapsis sits past 0.88 AU. Surface temperature swings ~800 K in
//! six hours during the periastron passage — the canonical example
//! of an unrelaxed, dynamically heated giant.
//!
//! Reference: Naef et al. (2001) discovery; Laughlin et al. (2009)
//! Spitzer thermal observations of the periastron passage.

use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem, builders::KG_M3_TO_SOLAR_AU3},
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub fn hd_80606(seed: u64) -> Template {
    let mut bodies = Vec::with_capacity(2);

    // HD 80606: G6V, M = 0.97 M_☉ (Pepe et al. 2002).
    let m_star = 0.97;
    // HD 80606 b: 4.115 M_Jup ≈ 3.93e-3 M_☉ (Hébrard et al. 2010).
    let m_planet = 3.93e-3;

    let a = 0.455;
    let e = 0.93;

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("HD 80606"),
        mass: m_star,
        preset: &body_preset::STAR,
        density: Some(1500.0 * KG_M3_TO_SOLAR_AU3),
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
        // Bulk density 980 kg/m³ from transit + RV (Hébrard et al.
        // 2010). Slightly less dense than Saturn — typical for an
        // inflated tidally heated giant.
        density: Some(980.0 * KG_M3_TO_SOLAR_AU3),
        position: Some(pos),
        velocity: vel,
    });

    Template {
        name: "HD 80606",
        description: "G6V star with the most eccentric known transiting hot Jupiter \
                      (HD 80606 b, e = 0.93). The planet swings from 0.03 AU at periapsis \
                      to 0.88 AU at apoapsis on a 111-day orbit, undergoing extreme tidal \
                      heating during the periastron passage.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0001),
        units: UnitSystem::solar_au(),
    }
}
