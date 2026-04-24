//! Non-gravitational perturbation force registration.

use crate::core::log::Source;
use crate::core::system::System;
use crate::physics::integrator::PerturbationForce;

impl System {
    /// Register a non-gravitational perturbation force.
    ///
    /// Applied at every subsequent integration step, additively with other
    /// registered perturbations. Remove with [`clear_perturbations`].
    ///
    /// # Softening-compatibility check
    ///
    /// If `p.requires_exact_gravity()` returns `true` (declared by the
    /// perturbation author — e.g. `PostNewtonian1PN` does) and any body in
    /// the system currently carries nonzero Plummer softening, this method
    /// emits a [`warn_diag!`](crate::warn_diag) diagnostic. The simulator's
    /// default softening (`EPS_BASE · mass^(1/3)`, ε ≈ 0.02 AU on the Sun)
    /// introduces a numerical apsidal precession that can silently dominate
    /// the perturbation's physical signal by orders of magnitude.
    ///
    /// Dismiss the warning by unsoftening the relevant bodies before or
    /// after this call:
    ///
    /// ```ignore
    /// sys = sys.with_exact_gravity();                  // whole system
    /// // or per-body at construction:
    /// let sun = Body::star(1.0).unsoftened();
    /// ```
    pub fn add_perturbation(&mut self, p: Box<dyn PerturbationForce>) {
        if p.requires_exact_gravity() {
            let softened = self.bodies.iter().filter(|b| b.softening != 0.0).count();
            if softened > 0 {
                let max_softening =
                    self.bodies.iter().map(|b| b.softening.abs()).fold(0.0_f64, f64::max);
                crate::warn_diag!(
                    Source::System,
                    "perturbation requires exact 1/r gravity but bodies are softened; \
                     call System::with_exact_gravity() or Body::unsoftened() — numerical \
                     apsidal precession from Plummer softening will otherwise swamp the signal",
                    softened_bodies = softened,
                    total_bodies = self.bodies.len(),
                    max_softening = max_softening,
                );
            }
        }
        self.perturbations.push(p);
    }

    /// Remove all registered perturbation forces.
    pub fn clear_perturbations(&mut self) {
        self.perturbations.clear();
    }

    pub fn perturbation_count(&self) -> usize {
        self.perturbations.len()
    }

    #[cfg(test)]
    pub(crate) fn softening_scale_value(&self) -> f64 {
        self.softening_scale
    }

    /// Zero the Plummer softening on every currently-registered body, and
    /// set the system's `softening_scale` to `0` so any body added later
    /// inherits exact `1/r` gravity too.
    ///
    /// Use at construction time for fine-physics experiments — post-Newtonian
    /// corrections, J2 oblateness, tidal dissipation — where the default
    /// material-scaled softening would contribute a numerical apsidal
    /// precession that dominates the measurement.
    ///
    /// ```ignore
    /// let mut sys = System::from_template(TemplateKind::SolarSystem)
    ///     .with_exact_gravity()
    ///     .with_integrator(IntegratorKind::Ias15);
    /// sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));
    /// ```
    #[must_use]
    pub fn with_exact_gravity(mut self) -> Self {
        for b in &mut self.bodies {
            b.softening = 0.0;
        }
        self.softening_scale = 0.0;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    #[test]
    fn with_exact_gravity_zeroes_existing_bodies() {
        let bodies = vec![Body::star(1.0), Body::rocky(3e-6)];
        // Pre-condition: material-scaled softening is nonzero.
        assert!(bodies.iter().all(|b| b.softening > 0.0));

        let sys = System::new(bodies).with_exact_gravity();
        assert!(sys.bodies().iter().all(|b| b.softening == 0.0));
        assert_eq!(sys.softening_scale_value(), 0.0);
    }

    #[test]
    fn with_exact_gravity_persists_for_later_added_bodies() {
        // Bodies added *after* `with_exact_gravity` must also end up
        // unsoftened — otherwise the guarantee is leaky.
        let mut sys = System::new(vec![Body::star(1.0)]).with_exact_gravity();
        sys.add_body(Body::rocky(3e-6));
        assert!(sys.bodies().iter().all(|b| b.softening == 0.0));
    }
}
