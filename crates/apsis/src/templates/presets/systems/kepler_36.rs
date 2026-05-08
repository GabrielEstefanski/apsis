//! Kepler-36 — two planets at extreme close approach.
//!
//! Kepler-36 b (rocky super-Earth) and Kepler-36 c (puffy
//! mini-Neptune) orbit a slightly-evolved G-type star with semi-major
//! axes 0.115 AU and 0.128 AU — separations at conjunction reach
//! 0.013 AU (~five Earth-Moon distances). The mutual perturbations
//! drive observable transit-timing variations and an exotic 6:7
//! mean-motion resonance signature.
//!
//! Densities differ by 8× between b and c despite their orbital
//! proximity — among the cleanest known examples of compositional
//! diversity in a tightly-packed system.
//!
//! Reference: Carter et al. (2012) Science 337, 556; Vissapragada
//! et al. (2020) refines b's density.

use crate::{
    domain::body_preset::{self, BodyPreset},
    templates::{
        Template, TemplateBody, UnitSystem,
        builders::{KG_M3_TO_SOLAR_AU3, circular_orbit},
    },
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub fn kepler_36(seed: u64) -> Template {
    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };
    let mut bodies = Vec::with_capacity(3);

    // Kepler-36: G1IV subgiant, 1.071 M_☉, ρ ≈ 1290 kg/m³ (slightly
    // less dense than the Sun — it has expanded onto the subgiant
    // branch).
    let m_star = 1.071;
    bodies.push(TemplateBody {
        name: Some("Kepler-36"),
        mass: m_star,
        preset: &body_preset::STAR,
        density: Some(1290.0 * KG_M3_TO_SOLAR_AU3),
        albedo: None,
        position: Some([0.0, 0.0, 0.0]),
        velocity: [0.0, 0.0, 0.0],
        class_override: None,
    });

    // Per-body data. Mass in M_☉ (M_⊕ × 3.003e-6); density in kg/m³.
    // Source: Carter et al. (2012) + Vissapragada et al. (2020).
    //
    // The large density gap (b ~7.46 g/cm³, c ~0.89 g/cm³) is the
    // headline feature — b is iron-rich rocky, c is a puffy
    // hydrogen-helium envelope on a small core.
    let planets: [(&str, f64, f64, &'static BodyPreset, f64); 2] = [
        // (name, a [AU], mass [M_☉],            preset,                ρ [kg/m³])
        ("b", 0.1153, 4.45 * 3.003e-6, &body_preset::ROCKY, 7460.0),
        ("c", 0.1283, 8.08 * 3.003e-6, &body_preset::GAS, 890.0),
    ];

    for (name, a, mass, preset, density_si) in planets {
        let phase = rng.random::<f64>() * std::f64::consts::TAU;

        let (pos, vel) = circular_orbit(m_star, a, phase);

        bodies.push(TemplateBody {
            name: Some(name),
            mass,
            preset,
            position: Some(pos),
            velocity: vel,
            class_override: None,
            density: Some(density_si * KG_M3_TO_SOLAR_AU3),
            albedo: None,
        });
    }

    Template {
        name: "Kepler-36",
        description: "Two planets in close, dynamically interacting orbits around a G1IV \
                      subgiant. Kepler-36 b (rocky super-Earth, 7.46 g/cm³) and Kepler-36 c \
                      (puffy mini-Neptune, 0.89 g/cm³) sit only 0.013 AU apart at \
                      conjunction. The 8× density gap and tight packing make this a \
                      benchmark for transit-timing-variation studies of compositional \
                      diversity.",
        bodies,
        // Inner planet b period 13.84 d → use ~period/100 ≈ 1×10⁻⁴
        // in T_AU units to resolve the close encounters cleanly.
        suggested_dt: Some(0.0001),
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        units: UnitSystem::solar_au(),
    }
}
