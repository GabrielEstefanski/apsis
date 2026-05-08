//! TRAPPIST-1 — seven Earth-sized planets in a chain of mean-motion
//! resonances around an ultra-cool M8V red dwarf.
//!
//! Star and planet parameters follow Agol et al. (2021) AJ 161:
//! 7-planet TTV solution with refined masses and bulk densities.
//! All seven planets are rocky / partly volatile-rich; planets e
//! through h orbit in the habitable zone.

use crate::{
    domain::body_preset,
    templates::{
        Template, TemplateBody, UnitSystem,
        builders::{KG_M3_TO_SOLAR_AU3, circular_orbit},
    },
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub fn trappist_1(seed: u64) -> Template {
    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };
    let mut bodies = Vec::with_capacity(8);

    // ── Star ─────────────────────────────────────────────────────── //
    // TRAPPIST-1: M8V ultra-cool red dwarf. Mass 0.0898 M_☉ (Mann et al.
    // 2019), bulk density 51 100 kg/m³ from interferometric radius
    // 0.1192 R_☉ (Van Grootel et al. 2018).
    let m_star = 0.0898;
    bodies.push(TemplateBody {
        name: Some("TRAPPIST-1"),
        mass: m_star,
        preset: &body_preset::RED_DWARF,
        density: Some(51_100.0 * KG_M3_TO_SOLAR_AU3),
        albedo: None,
        position: Some([0.0, 0.0, 0.0]),
        velocity: [0.0, 0.0, 0.0],
        class_override: None,
    });

    // ── Planets ──────────────────────────────────────────────────── //
    // Mass in M_☉ (Earth-mass × 3.003e-6); density in kg/m³ from
    // Agol et al. (2021) Table 4 + earlier Grimm et al. (2018) for h.
    // Semi-major axes from the same ephemerides; orbits assumed
    // circular for visualisation (real e are 10⁻³–10⁻², invisible
    // at template scale).
    let planets: [(&str, f64, f64, f64); 7] = [
        // (name, a [AU], mass [M_☉],         ρ [kg/m³])
        ("b", 0.01154, 1.374 * 3.003e-6, 5360.0),
        ("c", 0.01580, 1.308 * 3.003e-6, 5630.0),
        ("d", 0.02227, 0.388 * 3.003e-6, 4720.0),
        ("e", 0.02925, 0.692 * 3.003e-6, 5550.0),
        ("f", 0.03849, 1.039 * 3.003e-6, 5200.0),
        ("g", 0.04683, 1.321 * 3.003e-6, 4900.0),
        ("h", 0.06189, 0.326 * 3.003e-6, 4200.0),
    ];

    for (name, a, mass, density_si) in planets {
        let phase = rng.random::<f64>() * std::f64::consts::TAU;

        let (pos, vel) = circular_orbit(m_star, a, phase);

        bodies.push(TemplateBody {
            name: Some(name),
            mass,
            preset: &body_preset::ROCKY,
            density: Some(density_si * KG_M3_TO_SOLAR_AU3),
            albedo: None,
            position: Some(pos),
            velocity: vel,
            class_override: None,
        });
    }

    Template {
        name: "TRAPPIST-1",
        description: "Ultra-cool M8V red dwarf with seven Earth-sized planets in a chain \
                      of mean-motion resonances. Compact (innermost orbital period 1.51 d, \
                      outermost 18.77 d) and densely interacting; planets e–g sit in the \
                      habitable zone. Masses and densities from the Agol et al. (2021) \
                      TTV refinement.",
        bodies,
        // Inner planet b has period 1.51 d ≈ 0.00414 yr; integrator
        // needs ~period/100 ≈ 4×10⁻⁵ in T_AU units (1 yr = 2π T_AU)
        // to resolve close encounters.
        suggested_dt: Some(0.0001),
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        units: UnitSystem::solar_au(),
    }
}
