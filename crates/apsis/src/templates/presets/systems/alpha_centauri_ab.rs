//! Alpha Centauri triple-star system.
//!
//! The closest stellar system to the Sun. The inner pair (Alpha
//! Centauri A and B) form a 79.91-year eccentric binary at
//! a = 23.4 AU. Proxima Centauri orbits the AB pair at ~8 700 AU
//! with a period of ~547 000 yr (Kervella et al. 2017).
//!
//! ## Orbit-overlay caveat for Proxima
//!
//! The orbit-overlay path picks a single primary per body via Hill
//! sphere / dominant-mass heuristics. Proxima is gravitationally
//! bound to the *combined* mass of A + B (M ≈ 2.01 M_☉), but the
//! detector typically pairs it with A alone — and the vis-viva
//! relation against `M_A = 1.10 M_☉` puts Proxima past the escape
//! velocity for that single body, so the overlay draws a hyperbola.
//! The N-body simulation itself is unaffected: gravity sees the full
//! mass distribution and Proxima stays bound.
//!
//! Fixing the overlay requires extending the hierarchical-primary
//! detector to recognise tight binaries as a single barycentric
//! primary. Tracked separately.

use crate::{
    domain::body_preset,
    templates::{
        Template, TemplateBody, UnitSystem, builders::KG_M3_TO_SOLAR_AU3,
        keplerian::state_from_elements,
    },
};

pub fn alpha_centauri_ab(_seed: u64) -> Template {
    // ── Inner binary: A + B ─────────────────────────────────────────
    let m_a = 1.10_f64;
    let m_b = 0.91_f64;
    let a_ab = 23.4_f64;
    let e_ab = 0.52_f64;

    let m_ab_total = m_a + m_b;
    let r_peri_ab = a_ab * (1.0 - e_ab);
    let v_peri_ab = (m_ab_total * (1.0 + e_ab) / r_peri_ab).sqrt();

    let r1 = r_peri_ab * m_b / m_ab_total;
    let r2 = r_peri_ab * m_a / m_ab_total;
    let v1 = v_peri_ab * m_b / m_ab_total;
    let v2 = v_peri_ab * m_a / m_ab_total;

    // ── Proxima Centauri ────────────────────────────────────────────
    // Mass and orbital elements from Kervella, Thévenin & Lovis
    // (2017) A&A 598 — long-baseline radial-velocity + astrometric
    // solution that established Proxima as gravitationally bound to
    // AB (probability > 99.99%).
    let m_proxima = 0.1221_f64;
    // Orbit relative to the AB barycenter (origin in our frame
    // after A + B were placed about it). Inclination 107.6° is
    // measured from the AB orbital plane and produces a near-
    // perpendicular, slightly retrograde orbit; we route it through
    // `state_from_elements` so the inclination shows up in 3D.
    let (proxima_pos, proxima_vel) = state_from_elements(
        m_ab_total,
        8700.0,                 // a [AU]
        0.50,                   // e
        107.6_f64.to_radians(), // i (relative to AB plane)
        126.0_f64.to_radians(), // RAAN (ω̃ = 72°, Ω chosen
        72.0_f64.to_radians(),  // ω so the orbit is
        0.0,                    // visually distinct)
    );

    Template {
        name: "Alpha Centauri",
        description: "The closest stellar system to the Sun: an inner binary (Alpha Centauri \
                      A + B, 79.91-year orbit at e = 0.52) with Proxima Centauri orbiting at \
                      ~8 700 AU. Proxima is a low-mass M5.5V red dwarf with the only known \
                      Earth-mass habitable-zone planet (Proxima b, omitted here for \
                      simulation cadence reasons). Note: the orbit overlay may render \
                      Proxima as hyperbolic because the hierarchy detector pairs it with A \
                      alone rather than the AB barycenter — the underlying N-body \
                      simulation is correct.",
        bodies: vec![
            // Densities for A and B from Bigot et al. (2006)
            // interferometric radii + dynamical masses.
            TemplateBody {
                name: Some("Alpha Centauri A"),
                mass: m_a,
                position: Some([-r1, 0.0, 0.0]),
                velocity: [0.0, -v1, 0.0],
                class_override: None,
                preset: &body_preset::STAR,
                density: Some(1450.0 * KG_M3_TO_SOLAR_AU3),
                albedo: None,
            },
            TemplateBody {
                name: Some("Alpha Centauri B"),
                mass: m_b,
                position: Some([r2, 0.0, 0.0]),
                velocity: [0.0, v2, 0.0],
                class_override: None,
                preset: &body_preset::STAR,
                density: Some(1900.0 * KG_M3_TO_SOLAR_AU3),
                albedo: None,
            },
            // Proxima — M5.5V red dwarf, ρ ≈ 56 800 kg/m³ from
            // Demory et al. (2009) interferometric radius.
            TemplateBody {
                name: Some("Proxima Centauri"),
                mass: m_proxima,
                position: Some(proxima_pos),
                velocity: proxima_vel,
                preset: &body_preset::RED_DWARF,
                class_override: None,
                density: Some(56_800.0 * KG_M3_TO_SOLAR_AU3),
                albedo: None,
            },
        ],
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        suggested_dt: Some(0.002),
        units: UnitSystem::solar_au(),
    }
}
