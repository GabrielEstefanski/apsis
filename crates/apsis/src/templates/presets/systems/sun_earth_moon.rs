//! Sun–Earth–Moon hierarchical three-body system.
//!
//! The Earth orbits the Sun on a unit-radius circular orbit; the Moon
//! orbits the Earth at the real-world Earth–Moon distance ratio
//! (`a_moon ≈ 2.57e-3 AU`). Both phases are seeded so successive runs
//! with different seeds explore different starting configurations.
//!
//! Velocities are composed in the inertial frame: the Moon carries the
//! Earth's instantaneous velocity plus the local Earth-frame circular
//! velocity. This preserves the hierarchy under any integrator that
//! conserves total momentum.

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use std::f64::consts::TAU;

use crate::{
    domain::body_preset,
    templates::{
        Template, TemplateBody, UnitSystem,
        builders::{KG_M3_TO_SOLAR_AU3, circular_orbit},
    },
};

pub fn sun_earth_moon(seed: u64) -> Template {
    const M_SUN: f64 = 1.0;
    const M_EARTH: f64 = 3.0e-6;
    const M_MOON: f64 = 3.7e-8;
    const A_EARTH: f64 = 1.0;
    const A_MOON: f64 = 2.57e-3;

    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };

    let earth_phase = rng.random::<f64>() * TAU;
    let (earth_pos, earth_vel) = circular_orbit(M_SUN, A_EARTH, earth_phase);

    let moon_phase = rng.random::<f64>() * TAU;
    let (moon_rel_pos, moon_rel_vel) = circular_orbit(M_EARTH, A_MOON, moon_phase);
    let moon_pos = [
        earth_pos[0] + moon_rel_pos[0],
        earth_pos[1] + moon_rel_pos[1],
        earth_pos[2] + moon_rel_pos[2],
    ];
    let moon_vel = [
        earth_vel[0] + moon_rel_vel[0],
        earth_vel[1] + moon_rel_vel[1],
        earth_vel[2] + moon_rel_vel[2],
    ];

    Template {
        name: "Sun–Earth–Moon",
        description: "Hierarchical three-body system: Earth in a circular orbit around the Sun, \
                      Moon in a circular orbit around the Earth. Velocities are composed in the \
                      inertial frame so total linear momentum is exactly zero.",
        bodies: vec![
            TemplateBody {
                name: Some("Sun"),
                mass: M_SUN,
                position: Some([0.0, 0.0, 0.0]),
                velocity: [0.0, 0.0, 0.0],
                class_override: None,
                preset: &body_preset::STAR,
                density: Some(1408.0 * KG_M3_TO_SOLAR_AU3),
            },
            TemplateBody {
                name: Some("Earth"),
                mass: M_EARTH,
                position: Some(earth_pos),
                velocity: earth_vel,
                class_override: None,
                preset: &body_preset::ROCKY,
                density: Some(5514.0 * KG_M3_TO_SOLAR_AU3),
            },
            TemplateBody {
                name: Some("Moon"),
                mass: M_MOON,
                position: Some(moon_pos),
                velocity: moon_vel,
                class_override: None,
                preset: &body_preset::ICY,
                density: Some(3344.0 * KG_M3_TO_SOLAR_AU3),
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.0005),
        units: UnitSystem::solar_au(),
    }
}
