//! Velocity Verlet (Stormer-Verlet) integrator for N-body gravity.
//!
//! ## Algorithm - two-stage symplectic leapfrog
//!
//! Given state (**x**(t), **v**(t)) and accelerations **a**(**x**(t)):
//!
//! 1. **Half-kick**   `v(t + 1/2 dt) = v(t) + 1/2 a(x(t)) * dt`
//! 2. **Drift**       `x(t + dt)     = x(t) + v(t + 1/2 dt) * dt`
//! 3. **Force eval**  `a(t + dt)     = grad Phi(x(t + dt)) / m`
//! 4. **Half-kick**   `v(t + dt)     = v(t + 1/2 dt) + 1/2 a(t + dt) * dt`
//!
//! ## Symplectic property
//!
//! Velocity Verlet is a 2nd-order symplectic integrator: it conserves a
//! modified (shadow) Hamiltonian, so total energy oscillates around its initial
//! value without secular drift. Explicit Euler, by contrast, always loses or
//! gains energy monotonically.
//!
//! ## Time-reversal symmetry
//!
//! For any conservative (time-symmetric) force law, running N steps with `+dt`
//! then N steps with `-dt` exactly recovers the initial state (up to
//! floating-point rounding). This is a direct consequence of the algorithm's
//! structure: negating dt inverts each half-kick and drift in order.
//!
//! ## References
//! - Verlet (1967). *Phys. Rev.* 159, 98.
//! - Swope, Andersen, Berens, Wilson (1982). *J. Chem. Phys.* 76, 637.

use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;

/// Output produced by a single integration step.
pub struct StepResult<'a> {
    /// Accelerations at the end of the step (= a(t + dt)).
    pub acc1: &'a [(f64, f64)],

    /// Gravitational potential energy evaluated at the new positions x(t+dt).
    pub potential: f64,
}

/// Advance all bodies by one velocity-Verlet step of size `dt`.
pub fn step<'a>(
    bodies: &mut [Body],
    dt: f64,
    theta: f64,
    engine: &mut BarnesHutEngine,
    scratch_acc: &'a mut Vec<(f64, f64)>,
) -> StepResult<'a> {
    let n = bodies.len();

    if scratch_acc.len() != n {
        scratch_acc.resize(n, (0.0, 0.0));
    }

    evaluate_accelerations(bodies, theta, engine, scratch_acc);
    half_kick(bodies, scratch_acc, 0.5 * dt);
    drift(bodies, dt);

    let potential = evaluate_accelerations(bodies, theta, engine, scratch_acc);
    half_kick(bodies, scratch_acc, 0.5 * dt);

    StepResult {
        acc1: scratch_acc.as_slice(),
        potential,
    }
}

/// Rebuild the gravity structure and fill `scratch_acc` with current accelerations.
pub fn evaluate_accelerations(
    bodies: &[Body],
    theta: f64,
    engine: &mut BarnesHutEngine,
    scratch_acc: &mut Vec<(f64, f64)>,
) -> f64 {
    if scratch_acc.len() != bodies.len() {
        scratch_acc.resize(bodies.len(), (0.0, 0.0));
    }

    engine.build(bodies);
    engine.evaluate(bodies, theta, scratch_acc)
}

/// Apply a velocity kick `v += a * dt`.
pub fn half_kick(bodies: &mut [Body], acc: &[(f64, f64)], dt: f64) {
    for (body, &(ax, ay)) in bodies.iter_mut().zip(acc.iter()) {
        body.vx += ax * dt;
        body.vy += ay * dt;
    }
}

/// Advance all positions using the current velocities.
pub fn drift(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::energy::{angular_momentum_z, kinetic_energy, total_energy};

    fn make_engine() -> BarnesHutEngine {
        BarnesHutEngine::new(16)
    }

    /// Two-body circular orbit initial conditions for G=1, m=1, orbital radius=1.
    fn circular_orbit_bodies() -> Vec<Body> {
        vec![
            Body::new(-1.0, 0.0, 0.0, -0.5, 1.0),
            Body::new(1.0, 0.0, 0.0, 0.5, 1.0),
        ]
    }

    #[test]
    fn single_body_moves_uniformly_under_zero_force() {
        let mut bodies = vec![Body::new(1.0, 2.0, 3.0, -1.0, 1.0)];
        let mut engine = make_engine();
        let mut acc = Vec::new();
        let dt = 0.1;

        step(&mut bodies, dt, 0.5, &mut engine, &mut acc);

        assert!((bodies[0].x - 1.3).abs() < 1e-12, "x = {}", bodies[0].x);
        assert!((bodies[0].y - 1.9).abs() < 1e-12, "y = {}", bodies[0].y);
        assert!((bodies[0].vx - 3.0).abs() < 1e-12, "vx = {}", bodies[0].vx);
        assert!(
            (bodies[0].vy - (-1.0)).abs() < 1e-12,
            "vy = {}",
            bodies[0].vy
        );
    }

    #[test]
    fn time_reversal_recovers_initial_state() {
        let mut bodies = circular_orbit_bodies();
        let x0 = bodies[0].x;
        let y0 = bodies[0].y;
        let vx0 = bodies[0].vx;
        let vy0 = bodies[0].vy;

        let mut engine = make_engine();
        let mut acc = Vec::new();
        let dt = 0.01;
        let n_steps = 50;

        for _ in 0..n_steps {
            step(&mut bodies, dt, 0.5, &mut engine, &mut acc);
        }
        for _ in 0..n_steps {
            step(&mut bodies, -dt, 0.5, &mut engine, &mut acc);
        }

        assert!(
            (bodies[0].x - x0).abs() < 1e-10,
            "x drift = {}",
            bodies[0].x - x0
        );
        assert!(
            (bodies[0].y - y0).abs() < 1e-10,
            "y drift = {}",
            bodies[0].y - y0
        );
        assert!(
            (bodies[0].vx - vx0).abs() < 1e-10,
            "vx drift = {}",
            bodies[0].vx - vx0
        );
        assert!(
            (bodies[0].vy - vy0).abs() < 1e-10,
            "vy drift = {}",
            bodies[0].vy - vy0
        );
    }

    #[test]
    fn angular_momentum_z_conserved_in_two_body_orbit() {
        let mut bodies = circular_orbit_bodies();
        let lz0 = angular_momentum_z(&bodies);

        let mut engine = make_engine();
        let mut acc = Vec::new();

        for _ in 0..200 {
            step(&mut bodies, 0.005, 0.5, &mut engine, &mut acc);
        }

        let lz1 = angular_momentum_z(&bodies);
        let rel_err = (lz1 - lz0).abs() / lz0.abs();
        assert!(rel_err < 1e-8, "Lz relative drift = {:.2e}", rel_err);
    }

    #[test]
    fn total_energy_has_no_secular_drift() {
        let mut bodies = circular_orbit_bodies();
        let mut engine = make_engine();
        let mut acc = Vec::new();
        let dt = 0.005;

        let e0 = {
            let r = step(&mut bodies, dt, 0.5, &mut engine, &mut acc);
            total_energy(kinetic_energy(&bodies), r.potential)
        };

        for _ in 0..499 {
            step(&mut bodies, dt, 0.5, &mut engine, &mut acc);
        }

        let e1 = {
            let r = step(&mut bodies, dt, 0.5, &mut engine, &mut acc);
            total_energy(kinetic_energy(&bodies), r.potential)
        };

        let rel_err = (e1 - e0).abs() / e0.abs();
        assert!(rel_err < 1e-3, "energy relative drift = {:.2e}", rel_err);
    }
}
