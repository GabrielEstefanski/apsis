//! Osculating orbital element cache.

use crate::core::system::System;
use crate::physics::orbital::{self, OrbitalElements};

impl System {
    /// Recompute osculating orbital elements for all bodies and cache the result.
    ///
    /// O(N²) — call once per rendered frame, not every physics step.
    pub fn update_orbital_elements(&mut self) {
        self.orbital_cache = orbital::compute_all(&self.bodies, self.g_factor);
    }

    /// Cached osculating orbital elements (one slot per body).
    ///
    /// Call [`update_orbital_elements`] first to get fresh values.
    pub fn orbital_elements(&self) -> &[Option<OrbitalElements>] {
        &self.orbital_cache
    }
}
