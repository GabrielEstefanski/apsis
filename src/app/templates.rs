use std::f64::consts::PI;

use crate::domain::body::{Body, density_from_mass_radius};
use crate::domain::materials::Material;
use crate::physics::gravity::G;

const SOLAR_KG: f64 = 1.988_47e30;
const EARTH_KG: f64 = 5.972_2e24;
const JUPITER_KG: f64 = 1.898_13e27;
const AU_KM: f64 = 149_597_870.7;
const SOLAR_RADIUS_KM: f64 = 695_700.0;
const EARTH_RADIUS_KM: f64 = 6_371.0;
const JUPITER_RADIUS_KM: f64 = 69_911.0;

// With G = 1, choosing mass units as 4*pi^2 solar masses makes circular
// orbital speeds in AU come out naturally in AU / year.
const KEPLER_MASS_SCALE: f64 = 4.0 * PI * PI;

#[derive(PartialEq, Clone, Copy)]
pub enum TemplateCategory {
    Bodies,
    Formations,
    Collisions,
}

impl TemplateCategory {
    pub fn label(self) -> &'static str {
        match self {
            Self::Bodies => "REAL SYSTEMS",
            Self::Formations => "REAL CONFIGURATIONS",
            Self::Collisions => "IMPACT ANALOGS",
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

pub struct TemplateEntry {
    pub key: &'static str,
    pub label: &'static str,
    pub category: TemplateCategory,
}

pub const TEMPLATE_CATALOG: &[TemplateEntry] = &[
    TemplateEntry {
        key: "inner_solar",
        label: "Inner Solar",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "outer_solar",
        label: "Giant Planets",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "earth_moon",
        label: "Earth-Moon",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "jupiter_system",
        label: "Galilean Moons",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "hot_jupiter",
        label: "51 Pegasi b",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "trappist1",
        label: "TRAPPIST-1",
        category: TemplateCategory::Bodies,
    },
    TemplateEntry {
        key: "binary",
        label: "Pluto-Charon",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "belt",
        label: "Main Belt",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l1",
        label: "Sun-Earth L1",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l2",
        label: "Sun-Earth L2",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l4",
        label: "Jupiter Trojan L4",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "lagrange_l5",
        label: "Jupiter Trojan L5",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "hierarchical",
        label: "Kepler-16 AB-b",
        category: TemplateCategory::Formations,
    },
    TemplateEntry {
        key: "merge_head_on",
        label: "Arrokoth Analog",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "giant_impact",
        label: "Earth-Theia",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "scatter_flyby",
        label: "Earth-Apophis",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "shattering",
        label: "Vesta Family",
        category: TemplateCategory::Collisions,
    },
    TemplateEntry {
        key: "chain_reaction",
        label: "SL9 @ Jupiter",
        category: TemplateCategory::Collisions,
    },
];

#[inline]
fn deg(value: f64) -> f64 {
    value.to_radians()
}

#[inline]
fn solar_mass(mass_solar: f64) -> f64 {
    KEPLER_MASS_SCALE * mass_solar
}

#[inline]
fn kg(mass_kg: f64) -> f64 {
    solar_mass(mass_kg / SOLAR_KG)
}

#[inline]
fn au_from_km(distance_km: f64) -> f64 {
    distance_km / AU_KM
}

#[inline]
fn au_per_year_from_km_s(speed_km_s: f64) -> f64 {
    speed_km_s * 86_400.0 * 365.25 / AU_KM
}

fn body_real(
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    mass_kg: f64,
    radius_km: f64,
    material: Material,
) -> Body {
    let mut b = Body::new(x, y, vx, vy, kg(mass_kg), material);
    let radius_au = au_from_km(radius_km);
    b.density = density_from_mass_radius(b.mass, radius_au);
    b.sync_physical_properties();
    b.radius = b.physical_radius;
    b.softening = b.softening.max(b.physical_radius * 2.0);
    b
}

fn circular_orbit_real(
    primary: Body,
    orbital_radius: f64,
    phase: f64,
    mass_kg: f64,
    radius_km: f64,
    material: Material,
) -> Body {
    let mass = kg(mass_kg);
    let omega = (G * (primary.mass + mass) / orbital_radius.powi(3)).sqrt();
    let (c, s) = (phase.cos(), phase.sin());
    body_real(
        primary.x + orbital_radius * c,
        primary.y + orbital_radius * s,
        primary.vx - orbital_radius * omega * s,
        primary.vy + orbital_radius * omega * c,
        mass_kg,
        radius_km,
        material,
    )
}

fn barycentric_binary_real(
    separation: f64,
    phase: f64,
    mass_a_kg: f64,
    radius_a_km: f64,
    material_a: Material,
    mass_b_kg: f64,
    radius_b_km: f64,
    material_b: Material,
) -> [Body; 2] {
    let mass_a = kg(mass_a_kg);
    let mass_b = kg(mass_b_kg);
    let total = mass_a + mass_b;
    let omega = (G * total / separation.powi(3)).sqrt();
    let ra = separation * mass_b / total;
    let rb = separation * mass_a / total;
    let (c, s) = (phase.cos(), phase.sin());

    let a = body_real(
        -ra * c,
        -ra * s,
        ra * omega * s,
        -ra * omega * c,
        mass_a_kg,
        radius_a_km,
        material_a,
    );
    let b = body_real(
        rb * c,
        rb * s,
        -rb * omega * s,
        rb * omega * c,
        mass_b_kg,
        radius_b_km,
        material_b,
    );

    [a, b]
}

pub fn template_bodies(key: &str) -> Vec<Body> {
    match key {
        "inner_solar" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            vec![
                sun,
                circular_orbit_real(
                    sun,
                    0.387_098,
                    deg(252.251),
                    3.301_1e23,
                    2_439.7,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    sun,
                    0.723_332,
                    deg(181.979),
                    4.867_5e24,
                    6_051.8,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    sun,
                    1.0,
                    deg(100.464),
                    EARTH_KG,
                    EARTH_RADIUS_KM,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    sun,
                    1.523_679,
                    deg(355.453),
                    6.417_1e23,
                    3_389.5,
                    Material::Rocky,
                ),
            ]
        }

        "outer_solar" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            vec![
                sun,
                circular_orbit_real(
                    sun,
                    5.2044,
                    deg(34.404),
                    JUPITER_KG,
                    JUPITER_RADIUS_KM,
                    Material::Gas,
                ),
                circular_orbit_real(
                    sun,
                    9.5826,
                    deg(49.944),
                    5.683_4e26,
                    58_232.0,
                    Material::Gas,
                ),
                circular_orbit_real(
                    sun,
                    19.2184,
                    deg(313.232),
                    8.681_3e25,
                    25_362.0,
                    Material::IceGiant,
                ),
                circular_orbit_real(
                    sun,
                    30.11,
                    deg(304.880),
                    1.024_13e26,
                    24_622.0,
                    Material::IceGiant,
                ),
            ]
        }

        "earth_moon" => barycentric_binary_real(
            au_from_km(384_400.0),
            deg(135.0),
            EARTH_KG,
            EARTH_RADIUS_KM,
            Material::Rocky,
            7.346e22,
            1_737.4,
            Material::Rocky,
        )
        .into(),

        "jupiter_system" => {
            let jupiter = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                JUPITER_KG,
                JUPITER_RADIUS_KM,
                Material::Gas,
            );
            vec![
                jupiter,
                circular_orbit_real(
                    jupiter,
                    au_from_km(421_700.0),
                    deg(20.0),
                    8.931_9e22,
                    1_821.6,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    jupiter,
                    au_from_km(671_100.0),
                    deg(135.0),
                    4.799_8e22,
                    1_560.8,
                    Material::Icy,
                ),
                circular_orbit_real(
                    jupiter,
                    au_from_km(1_070_400.0),
                    deg(235.0),
                    1.481_9e23,
                    2_634.1,
                    Material::Icy,
                ),
                circular_orbit_real(
                    jupiter,
                    au_from_km(1_882_700.0),
                    deg(320.0),
                    1.075_9e23,
                    2_410.3,
                    Material::Icy,
                ),
            ]
        }

        "hot_jupiter" => {
            let star = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                1.11 * SOLAR_KG,
                1.237 * SOLAR_RADIUS_KM,
                Material::Star,
            );
            vec![
                star,
                circular_orbit_real(
                    star,
                    0.0520,
                    deg(40.0),
                    0.46 * JUPITER_KG,
                    1.9 * JUPITER_RADIUS_KM,
                    Material::Gas,
                ),
            ]
        }

        "trappist1" => {
            let star = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                0.0898 * SOLAR_KG,
                0.121 * SOLAR_RADIUS_KM,
                Material::BrownDwarf,
            );
            vec![
                star,
                circular_orbit_real(
                    star,
                    0.01154,
                    deg(15.0),
                    1.374 * EARTH_KG,
                    1.116 * EARTH_RADIUS_KM,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    star,
                    0.01580,
                    deg(120.0),
                    1.308 * EARTH_KG,
                    1.097 * EARTH_RADIUS_KM,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    star,
                    0.02227,
                    deg(215.0),
                    0.388 * EARTH_KG,
                    0.788 * EARTH_RADIUS_KM,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    star,
                    0.02928,
                    deg(300.0),
                    0.692 * EARTH_KG,
                    0.920 * EARTH_RADIUS_KM,
                    Material::Rocky,
                ),
                circular_orbit_real(
                    star,
                    0.03853,
                    deg(65.0),
                    1.039 * EARTH_KG,
                    1.045 * EARTH_RADIUS_KM,
                    Material::Icy,
                ),
                circular_orbit_real(
                    star,
                    0.04688,
                    deg(170.0),
                    1.321 * EARTH_KG,
                    1.129 * EARTH_RADIUS_KM,
                    Material::Icy,
                ),
                circular_orbit_real(
                    star,
                    0.06193,
                    deg(260.0),
                    0.326 * EARTH_KG,
                    0.773 * EARTH_RADIUS_KM,
                    Material::Icy,
                ),
            ]
        }

        "binary" => barycentric_binary_real(
            au_from_km(19_573.0),
            deg(20.0),
            1.303e22,
            1_188.3,
            Material::Icy,
            1.586e21,
            606.0,
            Material::Icy,
        )
        .into(),

        "belt" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            vec![
                sun,
                circular_orbit_real(
                    sun,
                    2.3615,
                    deg(45.0),
                    2.590_76e20,
                    262.7,
                    Material::Asteroid,
                ),
                circular_orbit_real(
                    sun,
                    2.7675,
                    deg(150.0),
                    9.383_5e20,
                    473.0,
                    Material::Asteroid,
                ),
                circular_orbit_real(sun, 2.773, deg(255.0), 2.14e20, 256.0, Material::Asteroid),
                circular_orbit_real(sun, 3.1415, deg(330.0), 8.67e19, 217.0, Material::Asteroid),
            ]
        }

        "lagrange_l1" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            let earth = circular_orbit_real(
                sun,
                1.0,
                deg(0.0),
                EARTH_KG,
                EARTH_RADIUS_KM,
                Material::Rocky,
            );
            let l1 = circular_orbit_real(
                sun,
                1.0 - au_from_km(1_500_000.0),
                deg(0.0),
                1_850.0,
                0.005,
                Material::Asteroid,
            );
            vec![sun, earth, l1]
        }

        "lagrange_l2" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            let earth = circular_orbit_real(
                sun,
                1.0,
                deg(0.0),
                EARTH_KG,
                EARTH_RADIUS_KM,
                Material::Rocky,
            );
            let l2 = circular_orbit_real(
                sun,
                1.0 + au_from_km(1_500_000.0),
                deg(0.0),
                6_200.0,
                0.010,
                Material::Asteroid,
            );
            vec![sun, earth, l2]
        }

        "lagrange_l4" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            let jupiter = circular_orbit_real(
                sun,
                5.2044,
                deg(0.0),
                JUPITER_KG,
                JUPITER_RADIUS_KM,
                Material::Gas,
            );
            let trojan =
                circular_orbit_real(sun, 5.2044, deg(60.0), 7.9e18, 57.0, Material::Asteroid);
            vec![sun, jupiter, trojan]
        }

        "lagrange_l5" => {
            let sun = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                SOLAR_KG,
                SOLAR_RADIUS_KM,
                Material::Star,
            );
            let jupiter = circular_orbit_real(
                sun,
                5.2044,
                deg(0.0),
                JUPITER_KG,
                JUPITER_RADIUS_KM,
                Material::Gas,
            );
            let trojan =
                circular_orbit_real(sun, 5.2044, deg(-60.0), 1.3e18, 32.0, Material::Asteroid);
            vec![sun, jupiter, trojan]
        }

        "hierarchical" => {
            let [star_a, star_b] = barycentric_binary_real(
                0.224_31,
                deg(10.0),
                0.6897 * SOLAR_KG,
                0.6489 * SOLAR_RADIUS_KM,
                Material::Star,
                0.202_55 * SOLAR_KG,
                0.22623 * SOLAR_RADIUS_KM,
                Material::Star,
            );
            let total_mass = star_a.mass + star_b.mass;
            let omega = (G * total_mass / 0.7048_f64.powi(3)).sqrt();
            let phase = deg(210.0);
            let (c, s) = (phase.cos(), phase.sin());
            let planet = body_real(
                0.7048 * c,
                0.7048 * s,
                -0.7048 * omega * s,
                0.7048 * omega * c,
                0.333 * JUPITER_KG,
                0.7538 * JUPITER_RADIUS_KM,
                Material::Gas,
            );
            vec![star_a, star_b, planet]
        }

        "merge_head_on" => {
            let left = body_real(
                -au_from_km(28_000.0),
                0.0,
                au_per_year_from_km_s(0.010),
                0.0,
                8.0e17,
                10.0,
                Material::Icy,
            );
            let right = body_real(
                au_from_km(28_000.0),
                0.0,
                -au_per_year_from_km_s(0.010),
                0.0,
                6.0e17,
                7.0,
                Material::Icy,
            );
            vec![left, right]
        }

        "giant_impact" => {
            let earth = body_real(
                -au_from_km(22_000.0),
                au_from_km(3_000.0),
                au_per_year_from_km_s(4.5),
                0.0,
                EARTH_KG,
                EARTH_RADIUS_KM,
                Material::Rocky,
            );
            let theia = body_real(
                au_from_km(28_000.0),
                0.0,
                -au_per_year_from_km_s(4.9),
                au_per_year_from_km_s(0.7),
                0.107 * EARTH_KG,
                0.53 * EARTH_RADIUS_KM,
                Material::Rocky,
            );
            vec![earth, theia]
        }

        "scatter_flyby" => {
            let earth = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                EARTH_KG,
                EARTH_RADIUS_KM,
                Material::Rocky,
            );
            let apophis = body_real(
                -au_from_km(90_000.0),
                au_from_km(38_000.0),
                au_per_year_from_km_s(7.4),
                -au_per_year_from_km_s(1.2),
                6.1e10,
                0.17,
                Material::Asteroid,
            );
            vec![earth, apophis]
        }

        "shattering" => {
            let vesta = body_real(
                -au_from_km(4_500.0),
                0.0,
                au_per_year_from_km_s(2.6),
                0.0,
                2.590_76e20,
                262.7,
                Material::Asteroid,
            );
            let impactor = body_real(
                au_from_km(6_500.0),
                au_from_km(500.0),
                -au_per_year_from_km_s(5.1),
                -au_per_year_from_km_s(0.3),
                3.0e18,
                35.0,
                Material::Asteroid,
            );
            vec![vesta, impactor]
        }

        "chain_reaction" => {
            let jupiter = body_real(
                0.0,
                0.0,
                0.0,
                0.0,
                JUPITER_KG,
                JUPITER_RADIUS_KM,
                Material::Gas,
            );
            let fragments = [
                (-180_000.0, 74_000.0, 20.8, -5.2, 2.5e13),
                (-150_000.0, 61_000.0, 20.7, -5.0, 2.0e13),
                (-120_000.0, 48_000.0, 20.6, -4.8, 1.7e13),
                (-90_000.0, 35_000.0, 20.5, -4.6, 1.4e13),
            ];
            let mut bodies = vec![jupiter];
            for (x_km, y_km, vx_km_s, vy_km_s, mass_kg) in fragments {
                bodies.push(body_real(
                    au_from_km(x_km),
                    au_from_km(y_km),
                    au_per_year_from_km_s(vx_km_s),
                    au_per_year_from_km_s(vy_km_s),
                    mass_kg,
                    1.0,
                    Material::Comet,
                ));
            }
            bodies
        }

        _ => vec![],
    }
}

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
                Material::Asteroid,
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
                Material::Asteroid,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{TEMPLATE_CATALOG, template_bodies};

    #[test]
    fn every_catalog_entry_resolves_to_bodies() {
        for entry in TEMPLATE_CATALOG {
            assert!(
                !template_bodies(entry.key).is_empty(),
                "template {} returned no bodies",
                entry.key
            );
        }
    }

    #[test]
    fn non_collision_templates_do_not_start_overlapping() {
        for key in [
            "inner_solar",
            "outer_solar",
            "earth_moon",
            "jupiter_system",
            "hot_jupiter",
            "trappist1",
            "binary",
            "belt",
            "lagrange_l1",
            "lagrange_l2",
            "lagrange_l4",
            "lagrange_l5",
            "hierarchical",
        ] {
            let bodies = template_bodies(key);
            for i in 0..bodies.len() {
                for j in (i + 1)..bodies.len() {
                    let dx = bodies[i].x - bodies[j].x;
                    let dy = bodies[i].y - bodies[j].y;
                    let distance = (dx * dx + dy * dy).sqrt();
                    let min_distance = bodies[i].physical_radius + bodies[j].physical_radius;
                    assert!(
                        distance > min_distance,
                        "template {key} starts with overlap between {i} and {j}"
                    );
                }
            }
        }
    }
}
