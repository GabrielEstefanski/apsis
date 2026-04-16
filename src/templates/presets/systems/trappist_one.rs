use crate::{
    core::materials::Material,
    templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit},
};

pub fn trappist_1() -> Template {
    let mut bodies = Vec::with_capacity(8);

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("TRAPPIST-1"),
        mass: 0.089,
        material: Material::Star,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        spin: 0.0,
    });

    let planets = [
        ("b", 0.0115, 3.0e-6),
        ("c", 0.0158, 3.5e-6),
        ("d", 0.0223, 1.0e-6),
        ("e", 0.0292, 2.5e-6),
        ("f", 0.0385, 3.0e-6),
        ("g", 0.0469, 3.5e-6),
        ("h", 0.0619, 1.0e-6),
    ];

    for (name, a, mass) in planets {
        let phase = rand::random::<f64>() * std::f64::consts::TAU;

        let (pos, vel) = circular_orbit(0.089, a, phase);

        bodies.push(TemplateBody {
            name: Some(name),
            mass,
            material: Material::Rocky,
            position: Some(pos),
            velocity: vel,
            spin: 0.0,
        });
    }

    Template {
        name: "TRAPPIST-1",
        description: "Compact resonant exoplanetary system.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0001), // CRÍTICO
        units: UnitSystem::solar_au(),
    }
}
