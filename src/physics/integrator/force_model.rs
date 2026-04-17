//! Trait abstraction for N-body force evaluation.
//!
//! [`ForceModel`] decouples integrators from the concrete force engine
//! (Barnes-Hut, direct O(N²), GPU, etc.).  Integrators call
//! [`ForceModel::compute`] without knowing which algorithm or data structure
//! produces the accelerations.
//!
//! [`GravityForceModel`] is the default implementation, wrapping a
//! [`BarnesHutEngine`] with a fixed opening angle θ.

use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;

// ── Trait ─────────────────────────────────────────────────────────────────────

/// A force model that can compute accelerations for a set of bodies.
///
/// Implementations own whatever internal state they need (quadtrees, GPU
/// buffers, neighbour lists, …).  The only contract is:
///
/// 1. After `compute()` returns, `acc[i]` holds the acceleration of body `i`.
/// 2. The return value is the **raw** (unscaled) potential energy.
///
/// `compute()` takes `&mut self` because most implementations need to rebuild
/// internal structures (e.g. the Barnes-Hut tree) every evaluation.
pub trait ForceModel: Send {
    /// Compute accelerations for `bodies` and write them into `acc`.
    ///
    /// `acc` is guaranteed to have length ≥ `bodies.len()`.
    /// Returns the raw gravitational potential energy (before any g_factor
    /// scaling).
    fn compute(&mut self, bodies: &[Body], acc: &mut [(f64, f64)]) -> f64;

    /// Barnes-Hut opening angle θ, if the model uses one.
    ///
    /// Defaults to `0.5`. Force models that do not use a tree-based
    /// approximation may ignore this and return any constant.
    fn theta(&self) -> f64 {
        0.5
    }

    /// Update the opening angle θ.
    ///
    /// No-op for models that do not use a hierarchical tree approximation.
    fn set_theta(&mut self, _theta: f64) {}

    /// Access the underlying Barnes-Hut engine for spatial queries.
    ///
    /// Returns `None` for force models that do not use a Barnes-Hut tree
    /// (e.g. direct O(N²), GPU kernels, relativistic corrections).
    /// Used by the adaptive-θ controller; the controller is automatically
    /// disabled when this returns `None`.
    fn bh_engine(&self) -> Option<&BarnesHutEngine> {
        None
    }
}

// ── GravityForceModel ─────────────────────────────────────────────────────────

/// Default force model: Barnes-Hut / exact O(N²) gravity with Plummer
/// softening, parameterised by the opening angle θ.
pub struct GravityForceModel {
    engine: BarnesHutEngine,
    theta: f64,
}

impl GravityForceModel {
    /// Create a new gravity force model.
    ///
    /// - `theta`:     Barnes-Hut opening angle (controls accuracy vs speed).
    /// - `max_depth`: Maximum quadtree depth (16 is sufficient for all
    ///                practical particle counts).
    pub fn new(theta: f64, max_depth: usize) -> Self {
        Self {
            engine: BarnesHutEngine::new(max_depth),
            theta,
        }
    }

}

impl ForceModel for GravityForceModel {
    fn compute(&mut self, bodies: &[Body], acc: &mut [(f64, f64)]) -> f64 {
        self.engine.build(bodies);
        self.engine.evaluate(bodies, self.theta, acc)
    }

    fn theta(&self) -> f64 {
        self.theta
    }

    fn set_theta(&mut self, theta: f64) {
        self.theta = theta.clamp(0.05, 1.5);
    }

    fn bh_engine(&self) -> Option<&BarnesHutEngine> {
        Some(&self.engine)
    }
}
