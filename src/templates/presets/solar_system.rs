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
pub fn inner_solar_system() -> Template {
    struct Planet {
        mass: f64,
        radius: f64,
        a: f64,
        phase: f64,
        material: Material,
    }

    let planets = [
        Planet { mass: 1.66e-7, radius: 0.0007, a: 0.387, phase: 0.0,           material: Material::Rocky },
        Planet { mass: 2.45e-6, radius: 0.0018, a: 0.723, phase: TAU / 4.0,     material: Material::Rocky },
        Planet { mass: 3.00e-6, radius: 0.0018, a: 1.000, phase: TAU * 2.0/5.0, material: Material::Rocky },
        Planet { mass: 3.23e-7, radius: 0.0010, a: 1.524, phase: TAU * 3.0/5.0, material: Material::Rocky },
    ];

    let mut bodies = vec![TemplateBody {
        mass: 1.0,
        radius: 0.02,
        position: Some([0.0, 0.0]),
        velocity: [0.0, 0.0],
        material: Material::Star,
    }];

    for p in &planets {
        let (pos, vel) = circular_orbit(1.0, p.a, p.phase);
        bodies.push(TemplateBody {
            mass: p.mass,
            radius: p.radius,
            position: Some(pos),
            velocity: vel,
            material: p.material,
        });
    }

    Template {
        name: "Inner Solar System",
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
