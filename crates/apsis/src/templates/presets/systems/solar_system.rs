//! Solar System template — heliocentric ecliptic J2000 frame.
//!
//! ## Sources
//!
//! All orbital elements are heliocentric and quoted in the J2000 mean
//! ecliptic frame, except moons of giant planets, whose elements are
//! given relative to the parent body's equator (IAU 2009 WGCCRE pole
//! orientation) and rotated into the ecliptic frame at template
//! instantiation by [`crate::templates::keplerian::parent_equator_basis`].
//!
//! | Quantity                              | Source                               |
//! |---------------------------------------|--------------------------------------|
//! | Planet `(a, e, i, Ω, ω̃)`              | NASA JPL `approx_pos.html` (J2000)   |
//! | Planet masses                         | IAU 2012 nominal (Prša et al. 2016)  |
//! | Dwarf planet elements                 | JPL Horizons body pages              |
//! | Moon elements                         | NASA SSD body fact sheets            |
//! | Parent pole orientations              | IAU 2009 WGCCRE Report               |
//! | Earth Moon ecliptic inclination       | Williams (2023) Moon fact sheet      |
//! | Asteroid main-belt i distribution     | Bottke et al. (2005)                 |
//! | Jupiter-family / long-period comet i  | Levison & Duncan (1997)              |
//!
//! ## Phase convention
//!
//! `(a, e, i, Ω, ω)` are real for every named body so orbital planes
//! and apsidal lines line up with the literature. The mean anomaly
//! at epoch `M₀` is randomised per seed: an interactive simulator
//! benefits from variety, and any specific J2000 phase would be
//! immediately obsolete by the next integration step anyway.
//!
//! ## Mass units
//!
//! All masses are in solar masses (`M_☉`); the planet masses translate
//! the IAU 2012 nominal values via `1 M_☉ = 1.989 × 10³⁰ kg`. Dwarf-
//! planet masses come from JPL's published kg values.
//!
//! ## Unit system
//!
//! `M_☉ / AU / T_AU` (with `G = 1`), so `Earth period = 2π T_AU = 1 yr`.

use std::f64::consts::TAU;

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use crate::domain::body_preset::{self, BodyPreset};
use crate::templates::keplerian::{
    parent_equator_basis, state_from_elements, state_from_elements_in_basis,
};
use crate::templates::{Template, TemplateBody, UnitSystem};

// ── Fundamental constants ─────────────────────────────────────────────────────

/// Sun mass in simulation units (canonical anchor).
const M_SUN: f64 = 1.0;

/// Earth mass in solar masses (IAU 2012 nominal).
const M_EARTH: f64 = 3.003e-6;

/// Earth–Moon mass ratio: k² = 1 / 81.3005 (IAU 2012).
const M_MOON: f64 = M_EARTH / 81.3005;

// ── Major-planet table ────────────────────────────────────────────────────────

struct PlanetData {
    name: &'static str,
    mass: f64,
    a: f64,
    e: f64,
    i_deg: f64,
    raan_deg: f64,
    argp_deg: f64,
    preset: &'static BodyPreset,
}

/// J2000 ecliptic mean orbital elements (NASA JPL `approx_pos.html`,
/// 1800-2050 fit). Argument of periapsis is computed from the published
/// longitude of periapsis ω̃ as `ω = ω̃ − Ω`.
const PLANETS: &[PlanetData] = &[
    PlanetData {
        name: "Mercury",
        mass: 1.660e-7,
        a: 0.38710,
        e: 0.20564,
        i_deg: 7.005,
        raan_deg: 48.331,
        argp_deg: 29.130,
        preset: &body_preset::ROCKY,
    },
    PlanetData {
        name: "Venus",
        mass: 2.448e-6,
        a: 0.72333,
        e: 0.00678,
        i_deg: 3.395,
        raan_deg: 76.680,
        argp_deg: 54.852,
        preset: &body_preset::ROCKY,
    },
    PlanetData {
        name: "Earth",
        mass: M_EARTH,
        a: 1.00000,
        e: 0.01671,
        i_deg: 0.000,
        raan_deg: 0.000,
        argp_deg: 102.938,
        preset: &body_preset::ROCKY,
    },
    PlanetData {
        name: "Mars",
        mass: 3.213e-7,
        a: 1.52366,
        e: 0.09339,
        i_deg: 1.850,
        raan_deg: 49.558,
        argp_deg: 286.502,
        preset: &body_preset::ROCKY,
    },
    PlanetData {
        name: "Jupiter",
        mass: 9.543e-4,
        a: 5.20336,
        e: 0.04839,
        i_deg: 1.304,
        raan_deg: 100.464,
        argp_deg: 274.255,
        preset: &body_preset::GAS,
    },
    PlanetData {
        name: "Saturn",
        mass: 2.857e-4,
        a: 9.53707,
        e: 0.05415,
        i_deg: 2.486,
        raan_deg: 113.665,
        argp_deg: 338.766,
        preset: &body_preset::GAS,
    },
    PlanetData {
        name: "Uranus",
        mass: 4.366e-5,
        a: 19.1913,
        e: 0.04717,
        i_deg: 0.773,
        raan_deg: 74.006,
        argp_deg: 96.999,
        preset: &body_preset::ICE_GIANT,
    },
    PlanetData {
        name: "Neptune",
        mass: 5.151e-5,
        a: 30.0690,
        e: 0.00859,
        i_deg: 1.770,
        raan_deg: 131.784,
        argp_deg: 273.187,
        preset: &body_preset::ICE_GIANT,
    },
];

// ── Dwarf planets and large TNOs ──────────────────────────────────────────────

/// Reuses [`PlanetData`]: same shape, different physical category.
/// Pluto is included here so the IAU dwarf set sits together.
const DWARFS: &[PlanetData] = &[
    PlanetData {
        name: "Ceres",
        mass: 4.72e-10,
        a: 2.7691,
        e: 0.0760,
        i_deg: 10.594,
        raan_deg: 80.305,
        argp_deg: 73.598,
        preset: &body_preset::ASTEROID,
    },
    PlanetData {
        name: "Pluto",
        mass: 6.55e-9,
        a: 39.4817,
        e: 0.2488,
        i_deg: 17.16,
        raan_deg: 110.299,
        argp_deg: 113.834,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Haumea",
        mass: 2.014e-9,
        a: 43.218,
        e: 0.19501,
        i_deg: 28.214,
        raan_deg: 121.79,
        argp_deg: 240.20,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Makemake",
        mass: 1.56e-9,
        a: 45.430,
        e: 0.16126,
        i_deg: 29.007,
        raan_deg: 79.382,
        argp_deg: 294.834,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Eris",
        mass: 8.28e-9,
        a: 67.781,
        e: 0.43607,
        i_deg: 44.044,
        raan_deg: 35.951,
        argp_deg: 151.639,
        preset: &body_preset::ICY,
    },
    PlanetData {
        // Mass uncertain (~1×10²¹ kg); JPL lists no determined value.
        // The estimate keeps Sedna in the dwarf-planet mass band and
        // makes its trajectory observable; it is not a Voyager-grade
        // physical claim.
        name: "Sedna",
        mass: 5.0e-10,
        a: 506.8,
        e: 0.8496,
        i_deg: 11.931,
        raan_deg: 144.248,
        argp_deg: 311.352,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Quaoar",
        mass: 7.04e-10,
        a: 43.694,
        e: 0.0392,
        i_deg: 7.989,
        raan_deg: 188.83,
        argp_deg: 147.479,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Orcus",
        mass: 3.18e-10,
        a: 39.387,
        e: 0.22701,
        i_deg: 20.592,
        raan_deg: 268.799,
        argp_deg: 72.310,
        preset: &body_preset::ICY,
    },
    PlanetData {
        name: "Gonggong",
        mass: 8.80e-10,
        a: 67.485,
        e: 0.5063,
        i_deg: 30.627,
        raan_deg: 336.866,
        argp_deg: 207.628,
        preset: &body_preset::ICY,
    },
];

// ── Moon table ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Parent {
    Earth,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
}

struct MoonData {
    name: &'static str,
    mass: f64,
    parent: Parent,
    a: f64,
    e: f64,
    /// Inclination relative to the parent's equator [deg], except for
    /// Earth's Moon, whose orbital plane is referenced to the ecliptic.
    i_deg: f64,
    raan_deg: f64,
    argp_deg: f64,
    preset: &'static BodyPreset,
}

const MOONS: &[MoonData] = &[
    // Earth
    MoonData {
        name: "Moon",
        mass: M_MOON,
        parent: Parent::Earth,
        a: 0.002569,
        e: 0.0549,
        i_deg: 5.145,
        raan_deg: 125.08,
        argp_deg: 318.15,
        preset: &body_preset::ROCKY,
    },
    // Jupiter — Galilean moons (Jupiter equatorial frame)
    MoonData {
        name: "Io",
        mass: 4.491e-8,
        parent: Parent::Jupiter,
        a: 0.002819,
        e: 0.0041,
        i_deg: 0.050,
        raan_deg: 43.977,
        argp_deg: 84.129,
        preset: &body_preset::ROCKY,
    },
    MoonData {
        name: "Europa",
        mass: 2.413e-8,
        parent: Parent::Jupiter,
        a: 0.004485,
        e: 0.0094,
        i_deg: 0.471,
        raan_deg: 219.106,
        argp_deg: 88.970,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Ganymede",
        mass: 7.452e-8,
        parent: Parent::Jupiter,
        a: 0.007155,
        e: 0.0013,
        i_deg: 0.204,
        raan_deg: 63.552,
        argp_deg: 192.417,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Callisto",
        mass: 5.410e-8,
        parent: Parent::Jupiter,
        a: 0.012585,
        e: 0.0074,
        i_deg: 0.205,
        raan_deg: 298.848,
        argp_deg: 52.643,
        preset: &body_preset::ICY,
    },
    // Saturn (Saturn equatorial frame)
    MoonData {
        name: "Mimas",
        mass: 1.886e-11,
        parent: Parent::Saturn,
        a: 0.001240,
        e: 0.0202,
        i_deg: 1.566,
        raan_deg: 173.027,
        argp_deg: 332.499,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Enceladus",
        mass: 5.430e-11,
        parent: Parent::Saturn,
        a: 0.001591,
        e: 0.0047,
        i_deg: 0.009,
        raan_deg: 169.506,
        argp_deg: 0.000,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Dione",
        mass: 5.508e-10,
        parent: Parent::Saturn,
        a: 0.002523,
        e: 0.0022,
        i_deg: 0.019,
        raan_deg: 290.415,
        argp_deg: 168.820,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Rhea",
        mass: 1.160e-9,
        parent: Parent::Saturn,
        a: 0.003524,
        e: 0.001,
        i_deg: 0.345,
        raan_deg: 351.042,
        argp_deg: 256.609,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Titan",
        mass: 6.764e-8,
        parent: Parent::Saturn,
        a: 0.008168,
        e: 0.0288,
        i_deg: 0.349,
        raan_deg: 28.058,
        argp_deg: 78.371,
        preset: &body_preset::ICY,
    },
    MoonData {
        name: "Iapetus",
        mass: 9.082e-10,
        parent: Parent::Saturn,
        a: 0.023803,
        e: 0.0286,
        i_deg: 15.470,
        raan_deg: 75.831,
        argp_deg: 271.606,
        preset: &body_preset::ICY,
    },
    // Uranus (Uranus equatorial frame — pole tilted 97.77° from
    // ecliptic, so satellites orbit ~perpendicular to it).
    MoonData {
        name: "Titania",
        mass: 1.774e-9,
        parent: Parent::Uranus,
        a: 0.002914,
        e: 0.0011,
        i_deg: 0.340,
        raan_deg: 99.771,
        argp_deg: 165.522,
        preset: &body_preset::ICY,
    },
    // Neptune — Triton orbits retrograde (i > 90° relative to parent
    // equator), the only large moon known to do so. The result in
    // ecliptic coordinates is a strongly inclined retrograde orbit.
    MoonData {
        name: "Triton",
        mass: 1.075e-8,
        parent: Parent::Neptune,
        a: 0.002372,
        e: 0.000016,
        i_deg: 156.865,
        raan_deg: 177.608,
        argp_deg: 234.412,
        preset: &body_preset::ICY,
    },
];

// ── Parent body lookup ────────────────────────────────────────────────────────

struct ParentInfo {
    mass: f64,
    /// `(α, δ)` of the spin axis in J2000 equatorial coords [deg].
    /// IAU 2009 WGCCRE Report, Table 1. Earth's pole is reported as
    /// `(0, 90)` purely so the lookup is total — the Moon never uses
    /// this rotation; its elements are already in the ecliptic frame.
    pole_ra_deg: f64,
    pole_dec_deg: f64,
}

fn parent_info(parent: Parent) -> ParentInfo {
    match parent {
        Parent::Earth => ParentInfo { mass: M_EARTH, pole_ra_deg: 0.0, pole_dec_deg: 90.0 },
        Parent::Jupiter => {
            ParentInfo { mass: 9.543e-4, pole_ra_deg: 268.057, pole_dec_deg: 64.495 }
        },
        Parent::Saturn => ParentInfo { mass: 2.857e-4, pole_ra_deg: 40.589, pole_dec_deg: 83.537 },
        Parent::Uranus => {
            ParentInfo { mass: 4.366e-5, pole_ra_deg: 257.311, pole_dec_deg: -15.175 }
        },
        Parent::Neptune => ParentInfo { mass: 5.151e-5, pole_ra_deg: 299.36, pole_dec_deg: 43.46 },
    }
}

// ── Distribution helpers (asteroid belt, comets) ──────────────────────────────

/// Box–Muller transform: convert two uniforms in `(0, 1]` into a
/// standard normal sample. Used for asteroid-belt inclination drawn
/// from `N(8°, 6°)` — Bottke et al. (2005) main-belt fit.
fn normal_sample(rng: &mut SmallRng) -> f64 {
    let u1: f64 = rng.random::<f64>().max(1e-12);
    let u2: f64 = rng.random::<f64>();
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

// ── Build ─────────────────────────────────────────────────────────────────────

/// Construct the heliocentric Solar System template.
pub fn solar_system(seed: u64) -> Template {
    let mut rng: SmallRng =
        if seed == 0 { rand::make_rng() } else { SmallRng::seed_from_u64(seed) };

    // Capacity: Sun + 8 planets + 9 dwarfs + 13 moons + 600 asteroids + 30 comets.
    let mut bodies = Vec::with_capacity(1 + PLANETS.len() + DWARFS.len() + MOONS.len() + 630);

    // ── Sun ───────────────────────────────────────────────────────────────── //
    bodies.push(TemplateBody {
        name: Some("Sun"),
        mass: M_SUN,
        preset: &body_preset::STAR,
        position: Some([0.0, 0.0, 0.0]),
        velocity: [0.0, 0.0, 0.0],
        class_override: None,
    });

    // ── Planets and dwarfs share the same heliocentric construction ──────── //
    // We capture each parent body's state as it is built, then use it to
    // place its satellites (planet + relative offset = inertial state).
    let mut parent_state: [Option<([f64; 3], [f64; 3])>; 5] = [None; 5];
    let parent_slot = |p: Parent| -> usize {
        match p {
            Parent::Earth => 0,
            Parent::Jupiter => 1,
            Parent::Saturn => 2,
            Parent::Uranus => 3,
            Parent::Neptune => 4,
        }
    };

    let push_helio = |bodies: &mut Vec<TemplateBody>,
                      parent_state: &mut [Option<([f64; 3], [f64; 3])>; 5],
                      rng: &mut SmallRng,
                      data: &PlanetData| {
        let mean_anom = rng.random::<f64>() * TAU;
        let (pos, vel) = state_from_elements(
            M_SUN,
            data.a,
            data.e,
            data.i_deg.to_radians(),
            data.raan_deg.to_radians(),
            data.argp_deg.to_radians(),
            mean_anom,
        );
        // Cache state for any parent that has moons in the table.
        let parent = match data.name {
            "Earth" => Some(Parent::Earth),
            "Jupiter" => Some(Parent::Jupiter),
            "Saturn" => Some(Parent::Saturn),
            "Uranus" => Some(Parent::Uranus),
            "Neptune" => Some(Parent::Neptune),
            _ => None,
        };
        if let Some(p) = parent {
            parent_state[parent_slot(p)] = Some((pos, vel));
        }
        bodies.push(TemplateBody {
            name: Some(data.name),
            mass: data.mass,
            preset: data.preset,
            position: Some(pos),
            velocity: vel,
            class_override: None,
        });
    };

    for p in PLANETS {
        push_helio(&mut bodies, &mut parent_state, &mut rng, p);
    }
    for d in DWARFS {
        push_helio(&mut bodies, &mut parent_state, &mut rng, d);
    }

    // ── Moons ─────────────────────────────────────────────────────────────── //
    for m in MOONS {
        let Some((parent_pos, parent_vel)) = parent_state[parent_slot(m.parent)] else {
            continue;
        };
        let info = parent_info(m.parent);

        // Earth's Moon uses ecliptic-frame elements directly; everything
        // else is given in the parent's equator and rotated into the
        // ecliptic via the IAU pole orientation.
        let mean_anom = rng.random::<f64>() * TAU;
        let (rel_pos, rel_vel) = if matches!(m.parent, Parent::Earth) {
            state_from_elements(
                info.mass,
                m.a,
                m.e,
                m.i_deg.to_radians(),
                m.raan_deg.to_radians(),
                m.argp_deg.to_radians(),
                mean_anom,
            )
        } else {
            let basis = parent_equator_basis(info.pole_ra_deg, info.pole_dec_deg);
            state_from_elements_in_basis(
                info.mass,
                m.a,
                m.e,
                m.i_deg.to_radians(),
                m.raan_deg.to_radians(),
                m.argp_deg.to_radians(),
                mean_anom,
                basis,
            )
        };

        bodies.push(TemplateBody {
            name: Some(m.name),
            mass: m.mass,
            preset: m.preset,
            position: Some([
                parent_pos[0] + rel_pos[0],
                parent_pos[1] + rel_pos[1],
                parent_pos[2] + rel_pos[2],
            ]),
            velocity: [
                parent_vel[0] + rel_vel[0],
                parent_vel[1] + rel_vel[1],
                parent_vel[2] + rel_vel[2],
            ],
            // Moons orbit a planet regardless of which preset (ROCKY for the
            // Moon, ICY for the Galileans/Saturnians, etc.) supplies their
            // density. The class filter groups them all under Moon so the
            // user can hide moons en masse without touching planets.
            class_override: Some(crate::domain::body_preset::BodyClass::Moon),
        });
    }

    // ── Main-belt asteroid swarm ──────────────────────────────────────────── //
    //
    // a ∈ [2.2, 3.2] AU, e ∈ [0.06, 0.35], i ∼ N(8°, 6°) clipped to [0°, 30°].
    // Inclination distribution from Bottke et al. (2005); the random Ω/ω
    // give a thick disc rather than a coplanar sheet.
    for _ in 0..600 {
        let a = 2.2 + rng.random::<f64>();
        let e = (0.06 + rng.random::<f64>() * 0.16).min(0.35);
        let inc_deg = (8.0_f64 + 6.0 * normal_sample(&mut rng)).clamp(0.0, 30.0);
        let raan_deg = rng.random::<f64>() * 360.0;
        let argp_deg = rng.random::<f64>() * 360.0;
        let mean_anom = rng.random::<f64>() * TAU;

        let (pos, vel) = state_from_elements(
            M_SUN,
            a,
            e,
            inc_deg.to_radians(),
            raan_deg.to_radians(),
            argp_deg.to_radians(),
            mean_anom,
        );

        bodies.push(TemplateBody {
            name: None,
            mass: 1e-10,
            preset: &body_preset::ASTEROID,
            position: Some(pos),
            velocity: vel,
            class_override: None,
        });
    }

    // ── Comets ────────────────────────────────────────────────────────────── //
    //
    // 20 Jupiter-family (e ∈ [0.5, 0.8], i ∈ [0°, 30°], periapsis 1–3 AU)
    // and 10 long-period (e ∈ [0.97, 0.999], periapsis 0.5–2 AU,
    // isotropic inclination including retrograde — Levison & Duncan 1997).
    for i in 0..30 {
        let (a, e, inc_deg, raan_deg, argp_deg) = if i < 20 {
            // Jupiter-family: prograde, low to moderate inclination.
            let r_peri = 1.0 + rng.random::<f64>() * 2.0;
            let e = 0.50 + rng.random::<f64>() * 0.30;
            let a = r_peri / (1.0 - e);
            let inc = rng.random::<f64>() * 30.0;
            let raan = rng.random::<f64>() * 360.0;
            let argp = rng.random::<f64>() * 360.0;
            (a, e, inc, raan, argp)
        } else {
            // Long-period: extreme eccentricity, inclination isotropic
            // on the sphere (uniform in cos i over [-1, 1]).
            let r_peri = 0.5 + rng.random::<f64>() * 1.5;
            let e = 0.97 + rng.random::<f64>() * 0.029;
            let a = r_peri / (1.0 - e);
            let cos_i = 1.0 - 2.0 * rng.random::<f64>();
            let inc = cos_i.acos().to_degrees();
            let raan = rng.random::<f64>() * 360.0;
            let argp = rng.random::<f64>() * 360.0;
            (a, e, inc, raan, argp)
        };

        let mean_anom = rng.random::<f64>() * TAU;
        let (pos, vel) = state_from_elements(
            M_SUN,
            a,
            e,
            inc_deg.to_radians(),
            raan_deg.to_radians(),
            argp_deg.to_radians(),
            mean_anom,
        );

        bodies.push(TemplateBody {
            name: None,
            mass: 1e-12,
            preset: &body_preset::COMET,
            position: Some(pos),
            velocity: vel,
            class_override: None,
        });
    }

    Template {
        name: "Solar System",
        description: "The Sun, eight planets, nine IAU and candidate dwarf planets (Ceres, \
                      Pluto, Haumea, Makemake, Eris, Sedna, Quaoar, Orcus, Gonggong), thirteen \
                      large moons (Earth's Moon plus the Galilean four, six major Saturnian \
                      moons, Titania, and Triton), the asteroid main belt, and a sample of \
                      Jupiter-family and long-period comets. Heliocentric ecliptic J2000 \
                      frame; orbital elements from NASA JPL; mean anomaly randomised per seed.",
        bodies,
        display_scale: 1.0,
        // Inner planets need ~Mercury-period / 1000 ≈ 0.001 yr per step
        // for stable Velocity Verlet. The 2π factor sets the implicit
        // period to one year for Earth.
        suggested_dt: Some(0.001),
        units: UnitSystem::solar_au(),
    }
}
