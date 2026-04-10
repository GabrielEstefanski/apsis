use std::f64::consts::TAU;

use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody, builders::circular_orbit};

pub fn solar_system() -> Template {
    use rand::random;

    struct Planet {
        mass: f64,
        a: f64,
        material: Material,
    }

    // helper: órbita relativa estável
    fn orbit_relative(
        parent_pos: [f64; 2],
        parent_vel: [f64; 2],
        parent_mass: f64,
        radius: f64,
        phase: f64,
    ) -> ([f64; 2], [f64; 2]) {
        let rel_pos = [radius * phase.cos(), radius * phase.sin()];

        let v = (parent_mass / radius).sqrt();

        let rel_vel = [-v * phase.sin(), v * phase.cos()];

        (
            [parent_pos[0] + rel_pos[0], parent_pos[1] + rel_pos[1]],
            [parent_vel[0] + rel_vel[0], parent_vel[1] + rel_vel[1]],
        )
    }

    let planets = [
        Planet {
            mass: 1.66e-7,
            a: 0.387,
            material: Material::Rocky,
        },
        Planet {
            mass: 2.45e-6,
            a: 0.723,
            material: Material::Rocky,
        },
        Planet {
            mass: 3.00e-6,
            a: 1.000,
            material: Material::Rocky,
        }, // Earth
        Planet {
            mass: 3.23e-7,
            a: 1.524,
            material: Material::Rocky,
        },
        Planet {
            mass: 9.5e-4,
            a: 5.204,
            material: Material::Gas,
        },
        Planet {
            mass: 2.86e-4,
            a: 9.58,
            material: Material::Gas,
        },
        Planet {
            mass: 4.37e-5,
            a: 19.2,
            material: Material::IceGiant,
        },
        Planet {
            mass: 5.15e-5,
            a: 30.05,
            material: Material::IceGiant,
        },
        Planet {
            mass: 6.5e-9,
            a: 39.48,
            material: Material::Icy,
        },
    ];

    let mut bodies = Vec::with_capacity(1 + planets.len() + 700);

    // ── Sun ── //
    bodies.push(TemplateBody {
        mass: 1.0,
        radius: 0.02,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        material: Material::Star,
    });

    let mut earth_state: Option<([f64; 2], [f64; 2])> = None;

    // ── Planets ── //
    for p in &planets {
        let phase = random::<f64>() * TAU;
        let (mut pos, mut vel) = circular_orbit(1.0, p.a, phase);

        // jitter leve (ok, não destrutivo)
        let eps = 1e-4;
        pos[0] *= 1.0 + eps * (random::<f64>() - 0.5);
        pos[1] *= 1.0 + eps * (random::<f64>() - 0.5);
        vel[0] *= 1.0 + eps * (random::<f64>() - 0.5);
        vel[1] *= 1.0 + eps * (random::<f64>() - 0.5);

        if (p.a - 1.0).abs() < 1e-6 {
            earth_state = Some((pos, vel));
        }

        bodies.push(TemplateBody {
            mass: p.mass,
            radius: 0.002,
            position: Some(pos),
            velocity: vel,
            material: p.material,
        });
    }

    // ── Moon (AGORA ESTÁVEL) ── //
    if let Some((earth_pos, earth_vel)) = earth_state {
        let moon_a = 0.001; // 🔥 dentro da zona estável (< ~0.002)

        let phase = random::<f64>() * TAU;

        let (mut moon_pos, mut moon_vel) =
            orbit_relative(earth_pos, earth_vel, 3.0e-6, moon_a, phase);

        // 🔧 pequeno ajuste pró-estabilidade (compensa o Sol)
        let correction = 0.98;
        moon_vel[0] *= correction;
        moon_vel[1] *= correction;

        bodies.push(TemplateBody {
            mass: 3.7e-8,
            radius: 0.0007,
            position: Some(moon_pos),
            velocity: moon_vel,
            material: Material::Icy,
        });
    }

    // ── Asteroid belt ── //
    for _ in 0..600 {
        let a = 2.2 + random::<f64>() * (3.2 - 2.2);
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
        name: "Solar System (stable hierarchical moon)",
        bodies,
        scale: 1.0,
    }
}
