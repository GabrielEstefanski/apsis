use std::f64::consts::TAU;

use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody, builders::circular_orbit};

/// Inner solar system in dimensionless G = 1 units.
///
/// Unit convention (Gaussian-like):
///   1 mass  unit = M_sun
///   1 length unit = 1 AU
///   G = 1  (so orbital period of Earth ≈ 2π time units)
///
/// Planet masses (M_sun):
///   Mercury  1.66e-7   Venus 2.45e-6   Earth 3.00e-6   Mars 3.23e-7
pub fn solar_system() -> Template {
    use rand::random;

    struct Planet {
        mass: f64,
        a: f64,
        material: Material,
    }

    let planets = [
        Planet {
            mass: 1.66e-7,
            a: 0.387,
            material: Material::Rocky,
        }, // Mercury
        Planet {
            mass: 2.45e-6,
            a: 0.723,
            material: Material::Rocky,
        }, // Venus
        Planet {
            mass: 3.00e-6,
            a: 1.000,
            material: Material::Rocky,
        }, // Earth
        Planet {
            mass: 3.23e-7,
            a: 1.524,
            material: Material::Rocky,
        }, // Mars
        Planet {
            mass: 9.5e-4,
            a: 5.204,
            material: Material::Gas,
        }, // Jupiter
    ];

    let mut bodies = vec![TemplateBody {
        mass: 1.0,
        radius: 0.02,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        material: Material::Star,
    }];

    // ── Planets ── //
    for p in &planets {
        let phase = random::<f64>() * TAU;

        let (mut pos, mut vel) = circular_orbit(1.0, p.a, phase);

        let eps = 0.02;
        pos[0] *= 1.0 + eps * (random::<f64>() - 0.5);
        pos[1] *= 1.0 + eps * (random::<f64>() - 0.5);

        vel[0] *= 1.0 + eps * (random::<f64>() - 0.5);
        vel[1] *= 1.0 + eps * (random::<f64>() - 0.5);

        bodies.push(TemplateBody {
            mass: p.mass,
            radius: 0.002,
            position: Some(pos),
            velocity: vel,
            material: p.material,
        });
    }

    // ── Asteroid belt ── //
    let belt_inner = 2.2;
    let belt_outer = 3.2;
    let belt_count = 600;

    for _ in 0..belt_count {
        let a = belt_inner + random::<f64>() * (belt_outer - belt_inner);
        let phase = random::<f64>() * TAU;

        let (mut pos, mut vel) = circular_orbit(1.0, a, phase);

        let jitter = 0.15;
        vel[0] *= 1.0 + jitter * (random::<f64>() - 0.5);
        vel[1] *= 1.0 + jitter * (random::<f64>() - 0.5);

        bodies.push(TemplateBody {
            mass: 1e-10,
            radius: 0.0003,
            position: Some(pos),
            velocity: vel,
            material: Material::Asteroid,
        });
    }

    // ── Comets ── //
    for _ in 0..20 {
        let a = 5.0 + random::<f64>() * 10.0;
        let phase = random::<f64>() * TAU;

        let (mut pos, mut vel) = circular_orbit(1.0, a, phase);

        vel[0] *= 0.6;
        vel[1] *= 0.6;

        bodies.push(TemplateBody {
            mass: 1e-12,
            radius: 0.0002,
            position: Some(pos),
            velocity: vel,
            material: Material::Comet,
        });
    }

    Template {
        name: "Solar System",
        bodies,
        scale: 1.0,
    }
}

/// Star with a single gas giant companion at 5 AU.
pub fn star_gas_giant() -> Template {
    let a = 5.0_f64;
    let (pos, vel) = circular_orbit(1.0, a, 0.0);

    Template {
        name: "Star + Gas Giant",
        bodies: vec![
            TemplateBody {
                mass: 1.0,
                radius: 0.04,
                position: Some([0.0, 0.0]),
                velocity: [0.0, 0.0],
                material: Material::Star,
            },
            TemplateBody {
                mass: 9.5e-4,
                radius: 0.006,
                position: Some(pos),
                velocity: vel,
                material: Material::Gas,
            },
        ],
        scale: 1.0,
    }
}
