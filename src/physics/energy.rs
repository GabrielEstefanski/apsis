use crate::domain::body::Body;

pub fn kinetic_energy(bodies: &[Body]) -> f64 {
    let mut k = 0.0;
    for b in bodies {
        k += 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy);
    }
    k
}

pub fn angular_momentum_z(bodies: &[Body]) -> f64 {
    let mut lz = 0.0;
    for b in bodies {
        lz += b.mass * (b.x * b.vy - b.y * b.vx);
    }
    lz
}

pub fn total_energy(kinetic: f64, potential: f64) -> f64 {
    kinetic + potential
}

pub fn center_of_mass_state(bodies: &[Body]) -> (f64, f64, f64, f64) {
    let mut m = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;

    for b in bodies {
        m += b.mass;
        x += b.mass * b.x;
        y += b.mass * b.y;
        vx += b.mass * b.vx;
        vy += b.mass * b.vy;
    }

    if m == 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    (x / m, y / m, vx / m, vy / m)
}

