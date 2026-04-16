use rand::random;
use std::f64::consts::TAU;

use crate::core::materials::Material;
use crate::templates::{Template, TemplateBody, UnitSystem, builders::circular_orbit};

pub fn simple_three_body() -> Template {
    let mut bodies = Vec::with_capacity(3);

    // ── Central star (Sun-like) ── //
    bodies.push(TemplateBody {
        name: Some("Primary Star"),
        mass: 1.0,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        material: Material::Star,
        spin: 0.0,
    });

    // ── Planet (Earth-like) ── //
    // Random orbital phase for variability
    let earth_phase = random::<f64>() * TAU;
    let (earth_pos, earth_vel) = circular_orbit(1.0, 1.0, earth_phase);

    bodies.push(TemplateBody {
        name: Some("Earth-like Planet"),
        mass: 3.0e-6,
        position: Some(earth_pos),
        velocity: earth_vel,
        material: Material::Rocky,
        spin: 0.0,
    });

    // ── Moon (hierarchical orbit around the planet) ── //
    let moon_phase = random::<f64>() * TAU;
    let moon_a = 0.00257;

    // Relative position in the planet frame
    let rel_pos = [moon_a * moon_phase.cos(), moon_a * moon_phase.sin()];

    // Orbital velocity relative to the planet (uses planet mass)
    let v = (3.0e-6 / moon_a).sqrt();
    let rel_vel = [-v * moon_phase.sin(), v * moon_phase.cos()];

    // Convert to inertial frame by adding planet state
    let moon_pos = [earth_pos[0] + rel_pos[0], earth_pos[1] + rel_pos[1]];
    let moon_vel = [earth_vel[0] + rel_vel[0], earth_vel[1] + rel_vel[1]];

    bodies.push(TemplateBody {
        name: Some("Moon"),
        mass: 3.7e-8,
        position: Some(moon_pos),
        velocity: moon_vel,
        material: Material::Icy,
        spin: 0.0,
    });

    Template {
        name: "3-body hierarchical",
        description: "Star–planet–moon hierarchical system with correct inertial composition.",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.0005),
        units: UnitSystem::solar_au(),
    }
}
