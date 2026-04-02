use std::f64::consts::PI;

use crate::domain::body::Body;
use crate::physics::gravity::G;

pub fn template_bodies(key: &str) -> Vec<Body> {
    match key {
        "binary" => vec![
            Body { x: -1.0, y: 0.0, vx: 0.0, vy: -0.5, mass: 1.0 },
            Body { x: 1.0, y: 0.0, vx: 0.0, vy: 0.5, mass: 1.0 },
        ],
        "figure8" => vec![
            Body { x: 0.97000436, y: -0.24308753, vx: 0.46620369, vy: 0.43236573, mass: 1.0 },
            Body { x: 0.0, y: 0.0, vx: -0.93240737, vy: -0.86473146, mass: 1.0 },
            Body { x: -0.97000436, y: 0.24308753, vx: 0.46620369, vy: 0.43236573, mass: 1.0 },
        ],
        "solar" => vec![
            Body { x: 0.0, y: 0.0, vx: 0.0, vy: 0.0, mass: 100.0 },
            Body { x: 5.0, y: 0.0, vx: 0.0, vy: 4.47, mass: 0.1 },
            Body { x: 10.0, y: 0.0, vx: 0.0, vy: 3.16, mass: 0.2 },
            Body { x: 18.0, y: 0.0, vx: 0.0, vy: 2.36, mass: 0.5 },
            Body { x: 28.0, y: 0.0, vx: 0.0, vy: 1.89, mass: 1.5 },
        ],
        "pythagorean" => vec![
            Body { x: -1.5, y: 0.0, vx: 0.0, vy: 0.0, mass: 3.0 },
            Body { x: 1.5, y: 0.0, vx: 0.0, vy: 0.0, mass: 4.0 },
            Body { x: 0.0, y: 2.0, vx: 0.0, vy: 0.0, mass: 5.0 },
        ],
        "belt" => {
            let mut v = vec![Body { x: 0.0, y: 0.0, vx: 0.0, vy: 0.0, mass: 100.0 }];
            let n = 120usize;
            let r = 12.0f64;
            let v_orb = (G * 100.0 / r).sqrt();
            for i in 0..n {
                let a = 2.0 * PI * i as f64 / n as f64;
                v.push(Body {
                    x: r * a.cos(),
                    y: r * a.sin(),
                    vx: -v_orb * a.sin(),
                    vy: v_orb * a.cos(),
                    mass: 0.001,
                });
            }
            v
        }
        "galaxies" => {
            let mut v = Vec::new();
            for g in 0..2i32 {
                let (cx, vxg) = if g == 0 { (-22.0f64, 0.7f64) } else { (22.0f64, -0.7f64) };
                v.push(Body { x: cx, y: 0.0, vx: vxg, vy: 0.0, mass: 60.0 });
                let nr = 24usize;
                let r = 5.0f64;
                let vo = (G * 60.0 / r).sqrt();
                for i in 0..nr {
                    let a = 2.0 * PI * i as f64 / nr as f64;
                    v.push(Body {
                        x: cx + r * a.cos(),
                        y: r * a.sin(),
                        vx: vxg - vo * a.sin(),
                        vy: vo * a.cos(),
                        mass: 0.1,
                    });
                }
            }
            v
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
            Body {
                x: cx + radius * a.cos(),
                y: cy + radius * a.sin(),
                vx: -orbit_vel * a.sin(),
                vy: orbit_vel * a.cos(),
                mass,
            }
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
            Body {
                x: cx + r * a.cos(),
                y: cy + r * a.sin(),
                vx: vel_disp * va.cos(),
                vy: vel_disp * va.sin(),
                mass,
            }
        })
        .collect()
}
