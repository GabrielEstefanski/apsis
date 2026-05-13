//! Helpers integrators call to apply operators at canonical positions
//! in the integration step. See
//! [`crate::physics::integrator::operator`] §Dispatch-contract.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Operator,
};

/// Sum the force contributions of every Hamiltonian and
/// non-conservative operator into `acc`. The caller owns `acc`'s
/// initialisation contract — typically the gravity force model writes
/// first and perturbations accumulate on top.
pub fn accumulate_perturbation_forces(
    bodies: &[Body],
    acc: &mut [Vec3],
    hamiltonian: &[Box<dyn HamiltonianOperator>],
    non_conservative: &[Box<dyn NonConservativeOperator>],
) {
    for op in hamiltonian {
        op.accumulate_force(bodies, acc);
    }
    for op in non_conservative {
        op.accumulate_force(bodies, acc);
    }
}

/// Sum `Σᵢ V_i` across Hamiltonian operators for inclusion in
/// [`crate::core::system::System::total_energy`].
pub fn total_hamiltonian_contribution(
    bodies: &[Body],
    hamiltonian: &[Box<dyn HamiltonianOperator>],
) -> f64 {
    hamiltonian.iter().map(|op| op.energy_contribution(bodies)).sum()
}

/// Dispatch the boundary `observe` hook on every registered operator
/// after the integrator has synchronised body state.
pub fn dispatch_observers(
    bodies: &[Body],
    t: f64,
    dt: f64,
    hamiltonian: &mut [Box<dyn HamiltonianOperator>],
    non_conservative: &mut [Box<dyn NonConservativeOperator>],
    observers: &mut [Box<dyn Operator>],
) {
    for op in hamiltonian.iter_mut() {
        op.observe(bodies, t, dt);
    }
    for op in non_conservative.iter_mut() {
        op.observe(bodies, t, dt);
    }
    for op in observers.iter_mut() {
        op.observe(bodies, t, dt);
    }
}
