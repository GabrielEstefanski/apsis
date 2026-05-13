//! Trait abstraction for N-body force evaluation.
//!
//! [`ForceModel`] decouples integrators from the concrete force engine
//! (Barnes-Hut, direct O(N²), GPU, etc.).  Integrators call
//! [`ForceModel::compute`] without knowing which algorithm or data structure
//! produces the accelerations.
//!
//! [`GravityForceModel`] is the default implementation, wrapping a
//! [`BarnesHutEngine`] with a fixed opening angle θ.

use std::sync::Arc;

use crate::domain::body::Body;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::gravity::{BarnesHutEngine, Kernel, PlummerKernel};

// ── Trait ─────────────────────────────────────────────────────────────────────

/// A force model that can compute accelerations for a set of bodies.
///
/// Implementations own whatever internal state they need (octrees, GPU
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
    fn compute(&mut self, bodies: &[Body], acc: &mut [Vec3]) -> f64;

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

    /// N threshold below which exact O(N²) evaluation is used instead of BH.
    ///
    /// Defaults to 64. Force models without a BH tree may return any constant.
    fn exact_threshold(&self) -> usize {
        64
    }

    /// Set the exact-evaluation threshold.
    ///
    /// No-op for force models that do not use a BH tree.
    fn set_exact_threshold(&mut self, _n: usize) {}

    /// Access the underlying Barnes-Hut engine for spatial queries.
    ///
    /// Returns `None` for force models that do not use a Barnes-Hut tree
    /// (e.g. direct O(N²), GPU kernels, relativistic corrections).
    /// Used by the adaptive-θ controller; the controller is automatically
    /// disabled when this returns `None`.
    fn bh_engine(&self) -> Option<&BarnesHutEngine> {
        None
    }

    /// Handle to the gravitational kernel this force model dispatches through.
    ///
    /// The default returns [`PlummerKernel`] for force models that do not
    /// have an explicit kernel concept — preserving the simulator's
    /// canonical Plummer-softened semantics for consumers that query
    /// kernel properties via [`Kernel::properties`].
    fn kernel(&self) -> Arc<dyn Kernel> {
        Arc::new(PlummerKernel::new())
    }

    /// Swap the active kernel.
    ///
    /// Default: no-op. Force models that carry a configurable kernel
    /// (e.g., [`GravityForceModel`] with its Barnes-Hut engine) override
    /// this to propagate the new kernel down.
    fn set_kernel(&mut self, _kernel: Arc<dyn Kernel>) {}

    /// Whether this force model is a deterministic function of state
    /// — i.e. `compute(bodies)` returns the same accelerations (to
    /// within f64 ULP) on two calls with identical `bodies`.
    ///
    /// Read by `System::set_integrator` to enforce the pairing rule
    /// with the integrator's
    /// [`requires_deterministic_force`](crate::physics::integrator::traits::Integrator::requires_deterministic_force).
    ///
    /// Default `true` covers the direct O(N²) case and force models
    /// with no hierarchical approximation. Implementations whose
    /// internal structure is position-dependent (BH tree rebuild,
    /// neighbour lists refreshed per-step) should override this to
    /// return `false` *when the approximation is active* — the
    /// determinism is a property of the current configuration, not
    /// only of the type.
    ///
    /// # Future evolution
    ///
    /// Will be upgraded to a `DeterminismLevel` enum once a second
    /// non-trivial force model (FMM, GPU) makes the `Strict` /
    /// `Approximate` distinction load-bearing. See
    /// [`Integrator::requires_deterministic_force`](crate::physics::integrator::traits::Integrator::requires_deterministic_force)
    /// for the corresponding evolution on the integrator side.
    fn is_deterministic(&self) -> bool {
        true
    }
}

// ── GravityForceModel ─────────────────────────────────────────────────────────

/// Default force model: Barnes-Hut / exact O(N²) gravity with Plummer
/// softening, parameterised by the opening angle θ.
///
/// Owns the per-`compute()` SoA snapshot buffer
/// ([`BodyArrays`](crate::domain::body_arrays::BodyArrays)). Each call to
/// [`compute`](Self::compute) re-packs the buffer from the integrator's
/// `Vec<Body>` and passes it to the engine. See
/// `docs/experiments/2026-05-10-soa-layout.md` §Design constraint for the
/// lifecycle contract.
pub struct GravityForceModel {
    engine: BarnesHutEngine,
    theta: f64,
    /// Per-`compute()` SoA snapshot. Reused across calls (allocator
    /// hits zero after warmup); written by `pack_from`, read by the
    /// engine, conceptually discarded when `compute` returns.
    body_arrays: BodyArrays,
}

impl GravityForceModel {
    /// Create a new gravity force model.
    ///
    /// - `theta`:     Barnes-Hut opening angle (controls accuracy vs speed).
    /// - `max_depth`: Maximum octree depth (16 is sufficient for all
    ///   practical particle counts).
    pub fn new(theta: f64, max_depth: usize) -> Self {
        Self { engine: BarnesHutEngine::new(max_depth), theta, body_arrays: BodyArrays::new() }
    }
}

impl ForceModel for GravityForceModel {
    fn kernel(&self) -> Arc<dyn Kernel> {
        self.engine.kernel()
    }

    fn set_kernel(&mut self, kernel: Arc<dyn Kernel>) {
        self.engine.set_kernel(kernel);
    }

    fn compute(&mut self, bodies: &[Body], acc: &mut [Vec3]) -> f64 {
        // Pack the SoA snapshot for this compute() call. Body positions
        // mutate between integrator substeps (IAS15's Picard iterations,
        // VV's drift), so the snapshot is rebuilt each call from the
        // authoritative AoS state. `pack_from` reuses buffer capacity;
        // allocations occur only when N grows past a previous high-water
        // mark.
        self.body_arrays.pack_from(bodies);

        // Phase-split instrumentation for the IAS15 diagnostic harness:
        // separate the tree-build half from the traversal half of the
        // evaluate work so an optimisation that caches the tree across
        // Picard iterations can be sized against real data rather than
        // speculation. Entirely compiled out when `ias15-profile` is
        // off; accesses thread-local storage owned by `ias15::profile`.
        #[cfg(feature = "ias15-profile")]
        let build_start = std::time::Instant::now();
        self.engine.build(&self.body_arrays);
        #[cfg(feature = "ias15-profile")]
        crate::physics::integrator::ias15::profile::record_tree_build(build_start.elapsed());

        #[cfg(feature = "ias15-profile")]
        let traverse_start = std::time::Instant::now();
        let pe = self.engine.evaluate(&self.body_arrays, self.theta, acc);
        #[cfg(feature = "ias15-profile")]
        crate::physics::integrator::ias15::profile::record_tree_traverse(traverse_start.elapsed());

        pe
    }

    fn theta(&self) -> f64 {
        self.theta
    }

    fn set_theta(&mut self, theta: f64) {
        self.theta = theta.clamp(0.05, 1.5);
    }

    fn exact_threshold(&self) -> usize {
        self.engine.exact_threshold()
    }

    fn set_exact_threshold(&mut self, n: usize) {
        self.engine.set_exact_threshold(n);
    }

    fn bh_engine(&self) -> Option<&BarnesHutEngine> {
        Some(&self.engine)
    }

    /// Delegates to [`BarnesHutEngine::is_direct_mode`]. The engine is
    /// a deterministic function of state iff it is configured so the
    /// BH branch is unreachable — i.e. `exact_threshold ≥ DIRECT_MODE_THRESHOLD`.
    ///
    /// Note the state-sensitive nature: `set_exact_threshold(usize::MAX)`
    /// flips this model to deterministic; any threshold below the
    /// clamp ceiling flips it back. `System::set_integrator` uses this
    /// to enforce the integrator/force-model compatibility rule.
    fn is_deterministic(&self) -> bool {
        self.engine.is_direct_mode()
    }
}
