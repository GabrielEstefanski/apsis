use crate::{
    core::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn kepler_36() -> Template {
    let mut bodies = Vec::with_capacity(3);

    // ── Star ───────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("Kepler-36"),
        mass: 1.07,
        material: Material::Star,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        spin: 0.0,
    });

    // Planets
    let planets = [
        ("b", 0.115, 4.0e-6, Material::Rocky),
        ("c", 0.128, 2.0e-5, Material::Gas),
    ];

    for (name, a, mass, material) in planets {
        let phase = rand::random::<f64>() * std::f64::consts::TAU;

        let (pos, vel) = circular_orbit(1.07, a, phase);

        bodies.push(TemplateBody {
            name: Some(name),
            mass,
            material,
            position: Some(pos),
            velocity: vel,
            spin: 0.0,
        });
    }

    Template {
        name: "Kepler-36",
        description: "Two planets in close, interacting orbits.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0002),
    }
}
