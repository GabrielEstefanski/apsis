use rand::random;
use std::f64::consts::TAU;

use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody, builders::circular_orbit};

pub fn simple_three_body() -> Template {
    let mut bodies = Vec::with_capacity(3);

    // ── Sun ── //
    bodies.push(TemplateBody {
        mass: 1.0,
        radius: 0.02,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        material: Material::Star,
    });

    // ── Earth ── //
    let earth_phase = random::<f64>() * TAU;
    let (earth_pos, earth_vel) = circular_orbit(1.0, 1.0, earth_phase);

    bodies.push(TemplateBody {
        mass: 3.0e-6,
        radius: 0.002,
        position: Some(earth_pos),
        velocity: earth_vel,
        material: Material::Rocky,
    });

    // ── Moon (hierárquico correto) ── //
    let moon_phase = random::<f64>() * TAU;
    let moon_a = 0.00257;

    // posição relativa à Terra
    let rel_pos = [moon_a * moon_phase.cos(), moon_a * moon_phase.sin()];

    // velocidade relativa (importante: usa massa da Terra)
    let v = (3.0e-6 / moon_a).sqrt();

    let rel_vel = [-v * moon_phase.sin(), v * moon_phase.cos()];

    // composição absoluta
    let moon_pos = [earth_pos[0] + rel_pos[0], earth_pos[1] + rel_pos[1]];

    let moon_vel = [earth_vel[0] + rel_vel[0], earth_vel[1] + rel_vel[1]];

    bodies.push(TemplateBody {
        mass: 3.7e-8,
        radius: 0.0007,
        position: Some(moon_pos),
        velocity: moon_vel,
        material: Material::Icy,
    });

    Template {
        name: "3-body hierarchical",
        bodies,
        scale: 1.0,
    }
}
