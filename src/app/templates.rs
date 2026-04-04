use std::f64::consts::PI;

use crate::domain::body::Body;
use crate::physics::gravity::G;

// ── Category ──────────────────────────────────────────────────────────────── //

#[derive(PartialEq, Clone, Copy)]
pub enum TemplateCategory {
    Bodies,
    Formations,
    Collisions,
}

impl TemplateCategory {
    pub fn label(self) -> &'static str {
        match self {
            Self::Bodies => "CELESTIAL BODIES",
            Self::Formations => "FORMATIONS",
            Self::Collisions => "COLLISIONS",
        }
    }
    pub fn grid_id(self) -> &'static str {
        match self {
            Self::Bodies => "tpl_bodies",
            Self::Formations => "tpl_formations",
            Self::Collisions => "tpl_collisions",
        }
    }
}

// ── Catalog ───────────────────────────────────────────────────────────────── //

pub struct TemplateEntry {
    pub key: &'static str,
    pub label: &'static str,
    pub category: TemplateCategory,
}

/// Add one entry here + one `match` arm in `template_bodies` to register a new scenario.
pub const TEMPLATE_CATALOG: &[TemplateEntry] = &[
    // Celestial bodies
    TemplateEntry {
        key: "inner_solar",
        label: "Inner Solar",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "outer_solar",
        label: "Outer Solar",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "earth_moon",
        label: "Earth-Moon",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "jupiter_system",
        label: "Jupiter Moons",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "hot_jupiter",
        label: "Hot Jupiter",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "trappist1",
        label: "TRAPPIST-1",
        category: TemplateCategory::Bodies,
    },
    // Formations
    TemplateEntry {
        key: "binary",
        label: "Binary Star",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "figure8",
        label: "Figure-8",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "pythagorean",
        label: "Pythagorean",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "belt",
        label: "Asteroid Belt",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "galaxies",
        label: "Galaxies",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l1",
        label: "Lagrange L1",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l2",
        label: "Lagrange L2",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l3",
        label: "Lagrange L3",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l4",
        label: "Lagrange L4",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l5",
        label: "Lagrange L5",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "hierarchical",
        label: "Hierarchical",
        category: TemplateCategory::Formations,
    },
    // Collisions
    TemplateEntry {
        key: "merge_head_on",
        label: "Head-on Merge",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "giant_impact",
        label: "Giant Impact",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "scatter_flyby",
        label: "Flyby Scatter",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "shattering",
        label: "Shattering",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "chain_reaction",
        label: "Chain React.",
        category: TemplateCategory::Collisions,
    },
];

// ── Initial conditions ────────────────────────────────────────────────────── //

/// Returns the initial body list for a given catalog key.
/// Circular orbit speed (G=1): v = sqrt(M_central / r).
pub fn template_bodies(key: &str) -> Vec<Body> {
    match key {
        // ── Celestial bodies ─────────────────────────────────────────────── //

        // Inner solar system: Sun + Mercury, Venus, Earth, Mars
        "inner_solar" => {
            let m_star = 100.0f64;
            let planets: &[(f64, f64)] = &[
                (2.5, 0.0002),  // Mercury-like
                (4.5, 0.0025),  // Venus-like
                (6.5, 0.0030),  // Earth-like
                (10.0, 0.0003), // Mars-like
            ];
            let mut v = vec![Body::new(0.0, 0.0, 0.0, 0.0, m_star)];
            for &(r, m) in planets {
                let vc = (G * m_star / r).sqrt();
                v.push(Body::new(r, 0.0, 0.0, vc, m));
            }
            v
        }

        // Outer solar system: Sun + Jupiter, Saturn, Uranus, Neptune
        "outer_solar" => {
            let m_star = 100.0f64;
            let planets: &[(f64, f64)] = &[
                (10.0, 0.30), // Jupiter-like
                (16.0, 0.10), // Saturn-like
                (22.0, 0.03), // Uranus-like
                (30.0, 0.03), // Neptune-like
            ];
            let mut v = vec![Body::new(0.0, 0.0, 0.0, 0.0, m_star)];
            for &(r, m) in planets {
                let vc = (G * m_star / r).sqrt();
                v.push(Body::new(r, 0.0, 0.0, vc, m));
            }
            v
        }

        // Earth-Moon system
        "earth_moon" => {
            let m_earth = 1.0f64;
            let m_moon = 0.012f64;
            let r_moon = 0.5f64;
            let vc = (G * m_earth / r_moon).sqrt();
            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_earth),
                Body::new(r_moon, 0.0, 0.0, vc, m_moon),
            ]
        }

        // Jupiter + four Galilean moons
        "jupiter_system" => {
            let m_jup = 5.0f64;
            // (r, mass) — Io, Europa, Ganymede, Callisto
            let moons: &[(f64, f64)] = &[(1.0, 0.008), (1.6, 0.015), (2.6, 0.020), (4.0, 0.010)];
            let mut v = vec![Body::new(0.0, 0.0, 0.0, 0.0, m_jup)];
            for &(r, m) in moons {
                let vc = (G * m_jup / r).sqrt();
                v.push(Body::new(r, 0.0, 0.0, vc, m));
            }
            v
        }

        // Hot Jupiter: close-in gas giant on tight orbit
        "hot_jupiter" => {
            let m_star = 50.0f64;
            let r = 1.5f64;
            let vc = (G * m_star / r).sqrt();
            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, vc, 2.0),
            ]
        }

        // TRAPPIST-1 simplified: compact M-dwarf + 7 small planets
        "trappist1" => {
            let m_star = 8.0f64;
            // (r, mass) — planets b through h
            let planets: &[(f64, f64)] = &[
                (0.40, 0.00055),
                (0.60, 0.00070),
                (0.80, 0.00045),
                (1.00, 0.00062),
                (1.30, 0.00068),
                (1.65, 0.00128),
                (2.20, 0.00033),
            ];
            let mut v = vec![Body::new(0.0, 0.0, 0.0, 0.0, m_star)];
            for &(r, m) in planets {
                let vc = (G * m_star / r).sqrt();
                v.push(Body::new(r, 0.0, 0.0, vc, m));
            }
            v
        }

        // ── Formations ───────────────────────────────────────────────────── //
        "binary" => vec![
            Body::new(-1.0, 0.0, 0.0, -0.5, 1.0),
            Body::new(1.0, 0.0, 0.0, 0.5, 1.0),
        ],

        "figure8" => vec![
            Body::new(0.97000436, -0.24308753, 0.46620369, 0.43236573, 1.0),
            Body::new(0.0, 0.0, -0.93240737, -0.86473146, 1.0),
            Body::new(-0.97000436, 0.24308753, 0.46620369, 0.43236573, 1.0),
        ],

        "pythagorean" => vec![
            Body::new(-1.5, 0.0, 0.0, 0.0, 3.0),
            Body::new(1.5, 0.0, 0.0, 0.0, 4.0),
            Body::new(0.0, 2.0, 0.0, 0.0, 5.0),
        ],

        "belt" => {
            let m_star = 100.0f64;
            let mut v = vec![Body::new(0.0, 0.0, 0.0, 0.0, m_star)];
            let n = 120usize;
            let r = 12.0f64;
            let v_orb = (G * m_star / r).sqrt();
            for i in 0..n {
                let a = 2.0 * PI * i as f64 / n as f64;
                v.push(Body::new(
                    r * a.cos(),
                    r * a.sin(),
                    -v_orb * a.sin(),
                    v_orb * a.cos(),
                    0.001,
                ));
            }
            v
        }

        "galaxies" => {
            let mut v = Vec::new();
            for g in 0..2i32 {
                let (cx, vxg) = if g == 0 {
                    (-22.0f64, 0.7f64)
                } else {
                    (22.0f64, -0.7f64)
                };
                let m_core = 60.0f64;
                v.push(Body::new(cx, 0.0, vxg, 0.0, m_core));
                let nr = 24usize;
                let r = 5.0f64;
                let vo = (G * m_core / r).sqrt();
                for i in 0..nr {
                    let a = 2.0 * PI * i as f64 / nr as f64;
                    v.push(Body::new(
                        cx + r * a.cos(),
                        r * a.sin(),
                        vxg - vo * a.sin(),
                        vo * a.cos(),
                        0.1,
                    ));
                }
            }
            v
        }

        "lagrange_l5" => {
            let m_star = 100.0f64;
            let r = 8.0f64;
            let v_p = (G * m_star / r).sqrt();
            let a = -PI / 3.0; // -60°

            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, v_p, 0.5),
                Body::new(
                    r * a.cos(),
                    r * a.sin(),
                    -v_p * a.sin(),
                    v_p * a.cos(),
                    0.01,
                ),
            ]
        }

        "lagrange_l1" => {
            let m_star = 100.0f64;
            let m_planet = 0.5f64;
            let r = 8.0f64;

            let mu = m_planet / (m_star + m_planet);
            let omega = (G * (m_star + m_planet) / r.powi(3)).sqrt();

            let x = r * (1.0 - (mu / 3.0).cbrt());
            let v = omega * x;

            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, (G * m_star / r).sqrt(), m_planet),
                Body::new(x, 0.0, 0.0, v, 0.01),
            ]
        }

        "lagrange_l2" => {
            let m_star = 100.0f64;
            let m_planet = 0.5f64;
            let r = 8.0f64;

            let mu = m_planet / (m_star + m_planet);
            let omega = (G * (m_star + m_planet) / r.powi(3)).sqrt();

            let x = r * (1.0 + (mu / 3.0).cbrt());
            let v = omega * x;

            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, (G * m_star / r).sqrt(), m_planet),
                Body::new(x, 0.0, 0.0, v, 0.01),
            ]
        }

        "lagrange_l3" => {
            let m_star = 100.0f64;
            let m_planet = 0.5f64;
            let r = 8.0f64;

            let mu = m_planet / (m_star + m_planet);
            let omega = (G * (m_star + m_planet) / r.powi(3)).sqrt();

            let x = -r * (1.0 + 5.0 * mu / 12.0);
            let v = omega * x.abs();

            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, (G * m_star / r).sqrt(), m_planet),
                Body::new(x, 0.0, 0.0, -v, 0.01),
            ]
        }

        // Lagrange L4 Trojan: star + planet + trojan asteroid at 60° ahead
        "lagrange_l4" => {
            let m_star = 100.0f64;
            let r = 8.0f64;
            let v_p = (G * m_star / r).sqrt();
            let a = PI / 3.0; // 60° ahead
            vec![
                Body::new(0.0, 0.0, 0.0, 0.0, m_star),
                Body::new(r, 0.0, 0.0, v_p, 0.5),
                Body::new(
                    r * a.cos(),
                    r * a.sin(),
                    -v_p * a.sin(),
                    v_p * a.cos(),
                    0.01,
                ),
            ]
        }

        // Hierarchical triple: tight inner binary + distant outer companion
        "hierarchical" => {
            let m = 1.0f64;
            let r_in = 1.5f64;
            let r_out = 10.0f64;
            // Binary: two masses m at ±r_in; v = sqrt(G*m / (4*r_in))
            let v_b = (G * m / (4.0 * r_in)).sqrt();
            // Outer companion orbiting total mass 2m at r_out
            let v_o = (G * 2.0 * m / r_out).sqrt();
            vec![
                Body::new(-r_in, 0.0, 0.0, -v_b, m),
                Body::new(r_in, 0.0, 0.0, v_b, m),
                Body::new(0.0, r_out, v_o, 0.0, 0.5),
            ]
        }

        // ── Collisions ────────────────────────────────────────────────────── //

        // Head-on merge: two equal planets on direct collision path
        "merge_head_on" => vec![
            Body::new(-5.0, 0.0, 0.3, 0.0, 2.0),
            Body::new(5.0, 0.0, -0.3, 0.0, 2.0),
        ],

        // Giant impact: large planet + smaller oblique impactor (Moon-forming analog)
        "giant_impact" => vec![
            Body::new(-3.0, 0.4, 0.7, 0.0, 3.0),
            Body::new(3.0, 0.0, -0.7, 0.06, 0.5),
        ],

        // Flyby scatter: fast hyperbolic pass — gravitational slingshot / scattering
        "scatter_flyby" => vec![
            Body::new(-12.0, 2.5, 2.8, 0.0, 2.0),
            Body::new(0.0, 0.0, 0.0, 0.0, 2.0),
        ],

        // Shattering: very high-velocity collision → strong fragmentation
        "shattering" => vec![
            Body::new(-6.0, 0.0, 3.5, 0.0, 2.0),
            Body::new(6.0, 0.0, -3.5, 0.0, 2.0),
        ],

        // Chain reaction: four bodies — first triggers a cascade
        "chain_reaction" => vec![
            Body::new(-12.0, 0.0, 2.5, 0.0, 1.5),
            Body::new(-3.0, 0.0, 0.0, 0.0, 1.5),
            Body::new(3.5, 0.0, 0.0, 0.0, 1.5),
            Body::new(12.0, 0.0, -1.8, 0.0, 1.5),
        ],

        _ => vec![],
    }
}

// ── Procedural spawners ───────────────────────────────────────────────────── //

pub fn spawn_ring(
    cx: f64,
    cy: f64,
    radius: f64,
    count: usize,
    mass: f64,
    orbit_vel: f64,
) -> Vec<Body> {
    (0..count)
        .map(|i| {
            let a = 2.0 * PI * i as f64 / count as f64;
            Body::new(
                cx + radius * a.cos(),
                cy + radius * a.sin(),
                -orbit_vel * a.sin(),
                orbit_vel * a.cos(),
                mass,
            )
        })
        .collect()
}

pub fn spawn_cluster(
    cx: f64,
    cy: f64,
    radius: f64,
    count: usize,
    mass: f64,
    vel_disp: f64,
) -> Vec<Body> {
    (0..count)
        .map(|_| {
            let r: f64 = rand::random::<f64>().sqrt() * radius;
            let a: f64 = rand::random::<f64>() * 2.0 * PI;
            let va: f64 = rand::random::<f64>() * 2.0 * PI;
            Body::new(
                cx + r * a.cos(),
                cy + r * a.sin(),
                vel_disp * va.cos(),
                vel_disp * va.sin(),
                mass,
            )
        })
        .collect()
}
