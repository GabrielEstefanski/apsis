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
