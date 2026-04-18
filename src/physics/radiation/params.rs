//! Per-body radiation interaction parameters.

use crate::physics::radiation::source::RadiationSource;
use std::f64::consts::PI;

/// Radiation interaction parameters intrinsic to a small body.
///
/// These are physical properties of the particle, independent of any
/// particular source or simulation frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiationParams {
    /// Geometric cross-section presented to the radiation field.
    /// For a sphere of radius `s`: `π s²`.
    pub area: f64,
    /// Mass of the body in internal mass units.
    pub mass: f64,
    /// Radiation pressure efficiency factor `Q_pr ∈ (0, 2]`.
    ///
    /// `Q_pr = 1` for a perfect absorber; `Q_pr = 2` for a perfect
    /// back-reflector; intermediate values for realistic grains.
    /// Wavelength-averaged values from Mie theory are appropriate here.
    pub q_pr: f64,
}

impl RadiationParams {
    /// Returns params that produce zero radiation force.
    ///
    /// Useful as a default for bodies that should not be affected
    /// (stars, massive planets).
    pub fn inert() -> Self {
        Self { area: 0.0, mass: 1.0, q_pr: 0.0 }
    }

    /// Computes the dimensionless β parameter for a given source.
    ///
    /// ```text
    /// β = F_rad / F_grav = (L · Q_pr · A) / (4π · μ · c · m)
    /// ```
    ///
    /// Requires `mu_source = G · M_source` from the gravitational subsystem.
    /// β ≈ 1 for micron-sized dust; β ≪ 1 for macroscopic bodies.
    pub fn beta(&self, source: &RadiationSource, mu_source: f64) -> f64 {
        if mu_source <= 0.0 || self.mass <= 0.0 {
            return 0.0;
        }
        (source.luminosity * self.q_pr * self.area) / (4.0 * PI * mu_source * source.c * self.mass)
    }
}
