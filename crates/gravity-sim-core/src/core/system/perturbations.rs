//! Non-gravitational perturbation force registration.

use crate::core::system::System;
use crate::physics::integrator::PerturbationForce;

impl System {
    /// Register a non-gravitational perturbation force.
    ///
    /// Applied at every subsequent integration step, additively with other
    /// registered perturbations.  Remove with [`clear_perturbations`].
    pub fn add_perturbation(&mut self, p: Box<dyn PerturbationForce>) {
        self.perturbations.push(p);
    }

    /// Remove all registered perturbation forces.
    pub fn clear_perturbations(&mut self) {
        self.perturbations.clear();
    }

    pub fn perturbation_count(&self) -> usize {
        self.perturbations.len()
    }
}
