use crate::{
    domain::body_preset::{self, BodyPreset},
    templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit},
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

pub fn kepler_36(seed: u64) -> Template {
    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };
    let mut bodies = Vec::with_capacity(3);

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("Kepler-36"),
        mass: 1.07,
        preset: &body_preset::STAR,
        position: Some([0.0, 0.0, 0.0]),
        velocity: [0.0, 0.0, 0.0],
        class_override: None,
    });

    // Planets — `&'static BodyPreset` carries the construction defaults
    // (density, colour, q_pr) per body without polluting the runtime.
    let planets: [(&str, f64, f64, &'static BodyPreset); 2] =
        [("b", 0.115, 4.0e-6, &body_preset::ROCKY), ("c", 0.128, 2.0e-5, &body_preset::GAS)];

    for (name, a, mass, preset) in planets {
        let phase = rng.random::<f64>() * std::f64::consts::TAU;

        let (pos, vel) = circular_orbit(1.07, a, phase);

        bodies.push(TemplateBody {
            name: Some(name),
            mass,
            preset,
            position: Some(pos),
            velocity: vel,
            class_override: None,
        });
    }

    Template {
        name: "Kepler-36",
        description: "Two planets in close, interacting orbits.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0002),
        units: UnitSystem::solar_au(),
    }
}
