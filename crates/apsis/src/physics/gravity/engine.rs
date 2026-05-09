//! Barnes-Hut force engine — orchestrates the octree and the Plummer kernel.
//!
//! ## Two evaluation strategies
//!
//! | N | Strategy | Complexity |
//! |---|---|---|
//! | N ≤ [`EXACT_THRESHOLD`] | Direct O(N²) pairwise sum | exact |
//! | N > [`EXACT_THRESHOLD`] | Barnes-Hut tree traversal | O(N log N), approximate |
//!
//! For small N the tree overhead exceeds the savings, so `exact_eval` is always
//! used.  For large N the BH approximation is controlled by the opening angle θ:
//! a cell of width s at distance d is accepted as a point mass when `s/d < θ`.
//!
//! ## Barnes-Hut criterion
//!
//! Given a node with aggregated mass M and COM at distance d from the target
//! body, the node is treated as a single pseudo-body when:
//!
//! ```text
//! s / d < θ   (s = cell side length)
//! ```
//!
//! - θ = 0: forces exact evaluation (recurse to all leaves)
//! - θ → ∞: monopole only (fast, inaccurate)
//! - Typical production values: 0.5 – 0.9
//!
//! ## References
//! - Barnes & Hut (1986). *Nature* 324, 446–449.
//! - Dehnen (2014). *Comput. Astrophys. Cosmol.* 1, 1.

use std::sync::Arc;

use crate::domain::body::Body;
use crate::math::Vec3;
use rayon::prelude::*;

use super::kernel::{G, Kernel, PlummerKernel, pair_eps2};
use super::tree::{DIRECT_MODE_THRESHOLD, EXACT_THRESHOLD, MultipoleOrder, NO_CHILD, Node, Octree};

// ── BarnesHutEngine ───────────────────────────────────────────────────────── //

/// N-body force engine using a Barnes-Hut octree.
///
/// Each call to [`build`](Self::build) reconstructs the octree from the
/// current body positions.  [`evaluate`](Self::evaluate) then computes
/// gravitational accelerations and total potential energy using that tree.
///
/// The engine contains no body state — it is safe to rebuild and re-evaluate
/// every step without any carry-over from previous steps.
///
/// The pair potential is supplied by a [`Kernel`] held behind [`Arc`]; the
/// default (from [`BarnesHutEngine::new`]) is [`PlummerKernel`], which
/// reproduces the Plummer-softened force law used throughout the library.
pub struct BarnesHutEngine {
    tree: Octree,
    /// N ≤ this → exact O(N²); N > this → Barnes-Hut traversal.
    exact_threshold: usize,
    kernel: Arc<dyn Kernel>,
    /// Multipole expansion order. See [`MultipoleOrder`] for toggle scope.
    multipole_order: MultipoleOrder,
}

impl BarnesHutEngine {
    /// Create a new engine with the default [`PlummerKernel`].
    ///
    /// `max_depth` bounds the octree depth; 16 is sufficient for all
    /// practical particle counts.
    pub fn new(max_depth: usize) -> Self {
        Self::with_kernel(max_depth, Arc::new(PlummerKernel::new()))
    }

    /// Create a new engine with a caller-supplied [`Kernel`] implementation.
    ///
    /// Use this to run the BH traversal against a non-Plummer kernel — for
    /// example, a kernel that demonstrates or tests a different Exactness
    /// or Continuity class.
    pub fn with_kernel(max_depth: usize, kernel: Arc<dyn Kernel>) -> Self {
        Self {
            tree: Octree::new(max_depth),
            exact_threshold: EXACT_THRESHOLD,
            kernel,
            multipole_order: MultipoleOrder::Monopole,
        }
    }

    /// Handle to the kernel this engine dispatches through.
    ///
    /// Used by [`System::add_perturbation`](crate::core::system::System::add_perturbation)
    /// to query the active kernel's
    /// [`KernelProperties`](crate::physics::gravity::kernel::KernelProperties)
    /// against each perturbation's
    /// [`KernelRequirements`](crate::physics::gravity::kernel::KernelRequirements).
    pub fn kernel(&self) -> Arc<dyn Kernel> {
        Arc::clone(&self.kernel)
    }

    /// Swap the active kernel.
    ///
    /// Used by
    /// [`System::with_kernel`](crate::core::system::System::with_kernel) to
    /// let researchers configure a non-default kernel for experiments such
    /// as the continuity counter-test (see
    /// [`TruncatedPlummerKernel`](crate::physics::gravity::kernel::TruncatedPlummerKernel)).
    pub fn set_kernel(&mut self, kernel: Arc<dyn Kernel>) {
        self.kernel = kernel;
    }

    /// Set the N threshold below which exact O(N²) evaluation is used.
    ///
    /// Range is clamped to `[1, DIRECT_MODE_THRESHOLD]`. Passing
    /// `usize::MAX` (or any value at or above `DIRECT_MODE_THRESHOLD`)
    /// forces the engine into "direct mode" — BH is never used
    /// regardless of body count. See [`is_direct_mode`](Self::is_direct_mode).
    pub fn set_exact_threshold(&mut self, n: usize) {
        self.exact_threshold = n.clamp(1, DIRECT_MODE_THRESHOLD);
    }

    /// Current exact-evaluation threshold.
    pub fn exact_threshold(&self) -> usize {
        self.exact_threshold
    }

    /// `true` iff the engine is configured so direct O(N²) summation
    /// is used for any practical body count — i.e.
    /// `exact_threshold() >= DIRECT_MODE_THRESHOLD`.
    ///
    /// This is the canonical way to ask "is the BH branch
    /// unreachable here?". Callers that need to reason about
    /// determinism (notably `ForceModel::is_deterministic`) should
    /// go through this rather than hard-coding the clamp ceiling.
    pub fn is_direct_mode(&self) -> bool {
        self.exact_threshold >= DIRECT_MODE_THRESHOLD
    }

    // ── Multipole order — experiment toggle ────────────────────────────────
    //
    // Removed in the final commit of the perf 2×2 experiment
    // (`docs/experiments/2026-05-08-octree-perf-2x2.md`) once §Decision
    // is written and the chosen multipole order is baked-in.

    /// Switch between [`MultipoleOrder::Monopole`] and
    /// [`MultipoleOrder::Quadrupole`] for subsequent [`build`] calls. The
    /// next [`build`] re-aggregates the tree under the new order; an
    /// already-built tree retains the order it was last built with until
    /// the engine rebuilds.
    ///
    /// [`build`]: Self::build
    #[allow(dead_code)] // perf 2x2 harness only; allow removed when bench lands
    pub(crate) fn set_multipole_order(&mut self, order: MultipoleOrder) {
        self.multipole_order = order;
    }

    /// Currently active multipole expansion order.
    #[allow(dead_code)] // read by perf 2x2 harness only; lib path passes the field directly
    pub(crate) fn multipole_order(&self) -> MultipoleOrder {
        self.multipole_order
    }

    /// Rebuild the octree from the current body positions.
    ///
    /// Must be called before [`evaluate`](Self::evaluate) whenever bodies have moved.
    pub fn build(&mut self, bodies: &[Body]) {
        self.tree.build(bodies, self.multipole_order);
    }

    /// Compute gravitational accelerations and return total potential energy.
    ///
    /// Fills `acc[i] = (aₓ, aᵧ, a_z)` for each body.
    /// Returns `PE = Σᵢ<ⱼ −G mᵢ mⱼ / r_ij` (softened).
    ///
    /// - N ≤ `exact_threshold`: uses exact O(N²) pairwise sum.
    /// - N > `exact_threshold`: uses parallel BH traversal.
    ///
    /// Spatial partition is the 3D octree (`Octree`) and the kernel
    /// arithmetic is fully 3D — `r² = Δx² + Δy² + Δz²` at every site.
    pub fn evaluate(&self, bodies: &[Body], theta: f64, acc: &mut [Vec3]) -> f64 {
        let n = bodies.len();
        acc.fill(Vec3::ZERO);

        if n == 0 {
            return 0.0;
        }

        let kernel: &dyn Kernel = &*self.kernel;

        if n <= self.exact_threshold {
            return exact_eval(bodies, kernel, acc);
        }

        let nodes = self.tree.nodes();

        let results: Vec<(Vec3, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                bh_eval_body(nodes, i, &bodies[i], bodies, theta, kernel, &mut stack)
            })
            .collect();

        let mut potential = 0.0_f64;
        for (i, (a, phi)) in results.into_iter().enumerate() {
            acc[i] = a;
            // phi is the specific potential at body i; multiply by mass for energy
            potential += bodies[i].mass * phi;
        }

        // Each pair counted once from each side → divide by 2
        0.5 * potential
    }

    /// Approximate θ-error proxy for a single body.
    ///
    /// Computes a mass-weighted RMS of `(s/d)²` over all nodes that would be
    /// accepted by the BH criterion at the given `theta`.  Used by the adaptive
    /// θ controller to estimate the current force truncation error.
    pub fn theta_error_proxy(&self, body_idx: usize, bodies: &[Body], theta: f64) -> f64 {
        if self.tree.nodes().is_empty() {
            return 0.0;
        }

        let body = &bodies[body_idx];
        let eps2 = body.softening * body.softening;
        let mut violation_sum = 0.0_f64;
        let mut weight_sum = 0.0_f64;

        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);

        while let Some(raw) = stack.pop() {
            let node = &self.tree.nodes()[raw as usize];

            if node.mass <= 0.0 || node.is_leaf() {
                continue;
            }

            let dx = node.com_x - body.x;
            let dy = node.com_y - body.y;
            let dz = node.com_z - body.z;
            let d = (dx * dx + dy * dy + dz * dz + eps2).sqrt();
            let ratio = node.size() / d;

            if ratio < theta {
                violation_sum += node.mass * ratio * ratio;
                weight_sum += node.mass;
            } else {
                for &c in &node.children {
                    if c != NO_CHILD {
                        stack.push(c);
                    }
                }
            }
        }

        if weight_sum > 0.0 { (violation_sum / weight_sum).sqrt() } else { 0.0 }
    }

    fn node_density(&self, node: &Node, x: f64, y: f64, z: f64, theta: f64) -> f64 {
        let dx = node.com_x - x;
        let dy = node.com_y - y;
        let dz = node.com_z - z;
        let dist2 = dx * dx + dy * dy + dz * dz + 1e-6;

        let size = node.size();

        if size * size / dist2 < theta * theta || node.is_leaf() {
            let dist = dist2.sqrt();
            return node.mass / dist;
        }

        let mut sum = 0.0;

        for &c in &node.children {
            if c != NO_CHILD {
                let child = &self.tree.nodes()[c as usize];
                sum += self.node_density(child, x, y, z, theta);
            }
        }

        sum
    }

    pub fn estimate_local_density(&self, x: f64, y: f64, z: f64, theta: f64) -> f64 {
        if self.tree.nodes().is_empty() {
            return 0.0;
        }

        let root = &self.tree.nodes()[0];
        self.node_density(root, x, y, z, theta)
    }

    pub fn query_neighbors(&self, x: f64, y: f64, z: f64, radius: f64, out: &mut Vec<usize>) {
        out.clear();

        let nodes = self.tree.nodes();
        if nodes.is_empty() {
            return;
        }

        self.query_node(nodes, 0, x, y, z, radius * radius, out);

        out.sort_unstable();
        out.dedup();
    }

    fn query_node(
        &self,
        nodes: &[Node],
        node_idx: u32,
        x: f64,
        y: f64,
        z: f64,
        radius2: f64,
        out: &mut Vec<usize>,
    ) {
        let node = &nodes[node_idx as usize];

        if node.mass <= 0.0 {
            return;
        }

        if !self.node_intersects(node, x, y, z, radius2) {
            return;
        }

        if node.is_leaf() {
            for k in 0..node.body_len as usize {
                out.push(node.bodies[k] as usize);
            }
            return;
        }

        for &c in &node.children {
            if c != NO_CHILD {
                self.query_node(nodes, c, x, y, z, radius2, out);
            }
        }
    }

    fn node_intersects(&self, node: &Node, x: f64, y: f64, z: f64, radius2: f64) -> bool {
        let half = node.half;

        let dx = ((x - node.cx).abs() - half).max(0.0);
        let dy = ((y - node.cy).abs() - half).max(0.0);
        let dz = ((z - node.cz).abs() - half).max(0.0);

        dx * dx + dy * dy + dz * dz <= radius2
    }
}

// ── Private evaluation strategies ─────────────────────────────────────────── //

/// Direct O(N²) pairwise force evaluation — exact for any N.
///
/// Iterates over all unique pairs (i, j).  For each pair, applies Newton's
/// 3rd law by updating both `acc[i]` and `acc[j]` from the same kernel
/// evaluation, using pairwise softening ε²_ij = (ε²_i + ε²_j)/2.
///
/// Returns the total gravitational potential energy PE = Σᵢ<ⱼ mᵢ Φᵢⱼ.
fn exact_eval(bodies: &[Body], kernel: &dyn Kernel, acc: &mut [Vec3]) -> f64 {
    let n = bodies.len();
    let mut potential = 0.0_f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = bodies[j].x - bodies[i].x;
            let dy = bodies[j].y - bodies[i].y;
            let dz = bodies[j].z - bodies[i].z;
            let eps2 = pair_eps2(bodies[i].softening, bodies[j].softening);
            let r_sq = dx * dx + dy * dy + dz * dz;

            // Shared geometric factor: G · f(r², ε²) = G / (r² + ε²)^(3/2)
            let fac = G * kernel.acceleration_factor(r_sq, eps2);

            // Newton's 3rd law: F_ij = −F_ji
            // acc_i += m_j · (dx, dy, dz) · fac
            // acc_j += m_i · (−dx, −dy, −dz) · fac
            //
            // The component-by-component `m · d · fac` chain is
            // load-bearing: re-associating into a shared
            // `m_fac = mass * fac` factor shifts ULPs and is observable
            // on the Mercury 1PN gate, which sits at the f64 noise
            // floor.
            acc[i].x += bodies[j].mass * dx * fac;
            acc[i].y += bodies[j].mass * dy * fac;
            acc[i].z += bodies[j].mass * dz * fac;
            acc[j].x -= bodies[i].mass * dx * fac;
            acc[j].y -= bodies[i].mass * dy * fac;
            acc[j].z -= bodies[i].mass * dz * fac;

            // Pair potential energy: E_ij = m_i · Φ_ij,  Φ_ij = −G · m_j · K
            let phi_ij = -G * bodies[j].mass * kernel.potential(r_sq, eps2);
            potential += bodies[i].mass * phi_ij;
        }
    }

    potential
}

/// Barnes-Hut force evaluation for a single body — O(log N) per body.
///
/// Returns `(a, φ)` where `a` is the acceleration vector and `φ` is the
/// specific gravitational potential (potential per unit mass) at the body's
/// position.  Multiply by `body.mass` to get the contribution to total PE.
///
/// Node interactions use the target body's own ε² — the tree stores only
/// aggregated mass and 3D COM, not per-body softening in internal nodes.
fn bh_eval_body(
    nodes: &[Node],
    body_idx: usize,
    body: &Body,
    bodies: &[Body],
    theta: f64,
    kernel: &dyn Kernel,
    stack: &mut Vec<u32>,
) -> (Vec3, f64) {
    let mut a = Vec3::ZERO;
    let mut phi = 0.0_f64;

    stack.clear();
    if !nodes.is_empty() {
        stack.push(0);
    }

    while let Some(raw) = stack.pop() {
        let node = &nodes[raw as usize];
        if node.mass <= 0.0 {
            continue;
        }

        if node.is_leaf() {
            // Exact pairwise kernel for all bodies in this leaf
            for k in 0..node.body_len as usize {
                let bi = node.bodies[k] as usize;
                if bi == body_idx {
                    continue;
                }
                let other = bodies[bi];
                let dx = other.x - body.x;
                let dy = other.y - body.y;
                let dz = other.z - body.z;
                let eps2 = pair_eps2(body.softening, other.softening);
                let r_sq = dx * dx + dy * dy + dz * dz;

                let fac = G * other.mass * kernel.acceleration_factor(r_sq, eps2);
                a.x += dx * fac;
                a.y += dy * fac;
                a.z += dz * fac;
                phi += -G * other.mass * kernel.potential(r_sq, eps2);
            }
            continue;
        }

        // BH criterion: accept this node as a pseudo-body when s/d < θ.
        let dx = node.com_x - body.x;
        let dy = node.com_y - body.y;
        let dz = node.com_z - body.z;
        let eps2 = body.softening * body.softening;
        let d = (dx * dx + dy * dy + dz * dz + eps2).sqrt();

        if node.size() / d < theta {
            let r_sq = dx * dx + dy * dy + dz * dz;
            let fac = G * node.mass * kernel.acceleration_factor(r_sq, eps2);
            a.x += dx * fac;
            a.y += dy * fac;
            a.z += dz * fac;
            phi += -G * node.mass * kernel.potential(r_sq, eps2);
        } else {
            for &c in &node.children {
                if c != NO_CHILD {
                    stack.push(c);
                }
            }
        }
    }

    (a, phi)
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    use approx::assert_relative_eq;

    fn eval(bodies: &[Body]) -> (Vec<Vec3>, f64) {
        let mut engine = BarnesHutEngine::new(16);
        engine.build(bodies);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        let potential = engine.evaluate(bodies, 0.5, &mut acc);
        (acc, potential)
    }

    fn body(x: f64, y: f64, m: f64) -> Body {
        Body::rocky(m).at(x, y).with_velocity(0.0, 0.0)
    }

    // ── Newton's 3rd law ───────────────────────────────────────────── //

    #[test]
    fn total_force_on_system_is_zero() {
        let bodies = vec![body(0.0, 0.0, 1.0), body(3.0, 0.0, 2.0)];

        let (acc, _) = eval(&bodies);

        let fx: f64 = acc.iter().zip(&bodies).map(|(a, b)| b.mass * a.x).sum();
        let fy: f64 = acc.iter().zip(&bodies).map(|(a, b)| b.mass * a.y).sum();
        let fz: f64 = acc.iter().zip(&bodies).map(|(a, b)| b.mass * a.z).sum();

        assert_relative_eq!(fx, 0.0, epsilon = 1e-12);
        assert_relative_eq!(fy, 0.0, epsilon = 1e-12);
        assert_relative_eq!(fz, 0.0, epsilon = 1e-12);
    }

    // ── Force direction ───────────────────────────────────────────── //

    #[test]
    fn force_direction_is_attractive() {
        let bodies = vec![body(0.0, 0.0, 1.0), body(4.0, 0.0, 1.0)];

        let (acc, _) = eval(&bodies);

        assert!(acc[0].x > 0.0);
        assert!(acc[1].x < 0.0);
    }

    // ── Superposition ───────────────────────────────────────────── //

    #[test]
    fn symmetric_configuration_has_zero_net_x_force_on_center() {
        let bodies = vec![body(-5.0, 0.0, 1.0), body(0.0, 0.0, 1.0), body(5.0, 0.0, 1.0)];

        let (acc, _) = eval(&bodies);

        assert_relative_eq!(acc[1].x, 0.0, epsilon = 1e-12);
    }

    // ── Potential sign ───────────────────────────────────────────── //

    #[test]
    fn gravitational_potential_is_negative() {
        let bodies = vec![body(0.0, 0.0, 1.0), body(2.0, 0.0, 1.0)];

        let (_, potential) = eval(&bodies);

        assert!(potential < 0.0);
    }

    // ── Barnes-Hut vs Exact ─────────────────────────────────────── //

    #[test]
    fn barnes_hut_matches_exact_with_small_error() {
        fn rel_err(a: f64, b: f64) -> f64 {
            (a - b).abs() / b.abs().max(1e-12)
        }

        let bodies = vec![
            body(-2.0, 0.0, 1.0),
            body(2.0, 0.0, 1.0),
            body(0.0, 3.0, 2.0),
            body(0.0, -3.0, 2.0),
        ];

        // Exato
        let mut engine_exact = BarnesHutEngine::new(16);
        engine_exact.set_exact_threshold(usize::MAX);
        engine_exact.build(&bodies);

        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        engine_exact.evaluate(&bodies, 0.5, &mut acc_exact);

        // BH
        let mut engine_bh = BarnesHutEngine::new(16);
        engine_bh.set_exact_threshold(1);
        engine_bh.build(&bodies);

        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        engine_bh.evaluate(&bodies, 0.5, &mut acc_bh);

        for i in 0..bodies.len() {
            let ex = acc_exact[i];
            let bh = acc_bh[i];

            assert!(rel_err(bh.x, ex.x) < 1e-2);
            assert!(rel_err(bh.y, ex.y) < 1e-2);
        }
    }

    // ── Octree validation (lab notebook 2026-05-08-octree-port) ────────── //
    //
    // Bounds declared a priori in
    // `docs/experiments/2026-05-08-octree-port.md`. Failure here means the
    // implementation is wrong, not the bound.

    /// Tier 1 — Barnes-Hut force accuracy on a general 3D distribution.
    ///
    /// 100 bodies sampled in the unit sphere with log-normal masses, fixed
    /// seed for reproducibility. Per-body max relative force error against
    /// exact O(N²) must stay under the Salmon-Warren 5% bound at θ = 0.5.
    #[test]
    fn tier1_octree_bh_force_error_under_5pct_at_theta_0_5() {
        let bodies = sphere_distribution_lognormal(100, 0x6F637472);

        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&bodies);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&bodies, 0.5, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&bodies);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&bodies, 0.5, &mut acc_bh);

        let max_rel = body_max_rel_error(&acc_bh, &acc_exact);
        eprintln!("[octree-tier1] θ=0.5 max rel-err = {max_rel:.4e}");
        assert!(
            max_rel <= 5e-2,
            "max per-body rel-err = {max_rel:.4e} exceeds 5e-2 (Salmon-Warren) at θ = 0.5"
        );
    }

    /// Tier 1 (exact-mode sanity) — Newton's third law at the round-off
    /// floor. Exact pairwise evaluation IS symmetric by construction, so
    /// `Σ m_i a_i` accumulates only floating-point summation noise. BH
    /// mode is not gated on this — the monopole approximation breaks
    /// pairwise symmetry by design (body A sees a far node at its COM,
    /// the bodies inside that node see A individually; the action and
    /// reaction sums are not algebraically equal). Failure of this test
    /// would indicate a defect in the exact pairwise kernel, not in the
    /// BH traversal.
    #[test]
    fn tier1_exact_mode_preserves_newton_third_law_at_roundoff() {
        let bodies = sphere_distribution_lognormal(100, 0x6F637472);

        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&bodies);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&bodies, 0.5, &mut acc);

        let net: Vec3 = acc.iter().zip(&bodies).fold(Vec3::ZERO, |s, (a, b)| s + b.mass * *a);
        eprintln!("[octree-tier1] exact mode |Σ m a| = {:.4e}", net.length());
        assert!(
            net.length() < 1e-12,
            "exact-mode Σ m_i a_i = {} exceeds round-off floor 1e-12",
            net.length(),
        );
    }

    /// Tier 1 — Loose-θ regime. Same bodies, θ = 0.9 widens the bound to 10 %.
    #[test]
    fn tier1_octree_bh_force_error_under_10pct_at_theta_0_9() {
        let bodies = sphere_distribution_lognormal(100, 0x6F637472);

        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&bodies);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&bodies, 0.9, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&bodies);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&bodies, 0.9, &mut acc_bh);

        let max_rel = body_max_rel_error(&acc_bh, &acc_exact);
        eprintln!("[octree-tier1] θ=0.9 max rel-err = {max_rel:.4e}");
        assert!(max_rel <= 1e-1, "max per-body rel-err = {max_rel:.4e} exceeds 1e-1 at θ = 0.9");
    }

    /// Tier 2 — Inclined Kepler. Two-body at i = 30° padded to N > EXACT_THRESHOLD
    /// so the BH branch is exercised. Integrate 100 orbital periods with
    /// Velocity Verlet at `dt = T/200`. Specific angular momentum `|L|/|L₀|`
    /// must drift no more than 1 × 10⁻³ — the Bug #4 bound from the WH
    /// refactor (`docs/experiments/2026-05-03-wh-refactor.md`), reused here
    /// because Lz was the diagnostic that originally caught 2D-only defects.
    #[test]
    fn tier2_octree_inclined_kepler_lz_below_1e_minus_3() {
        // 2-body Kepler, mass ratio 1:1e-3, e = 0.3, i = 30°.
        let m_central = 1.0_f64;
        let m_planet = 1.0e-3_f64;
        let a = 1.0_f64;
        let e = 0.3_f64;
        let inc = 30.0_f64.to_radians();
        let mu = m_central + m_planet;
        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / mu).sqrt();

        // Periapsis state in orbital plane, then rotated by i around x-axis.
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) * mu / (a * (1.0 - e))).sqrt();
        let (sin_i, cos_i) = inc.sin_cos();

        let mut bodies = vec![
            // Central body at origin
            Body::rocky(m_central).at(0.0, 0.0).with_velocity(0.0, 0.0),
            // Planet at periapsis with inclined velocity
            {
                let mut b = Body::rocky(m_planet).at(r_peri, 0.0).with_velocity(0.0, 0.0);
                b.vy = v_peri * cos_i;
                b.vz = v_peri * sin_i;
                b
            },
        ];

        // Pad to N > EXACT_THRESHOLD with massless test particles far from
        // the binary so they don't perturb the Kepler orbit. Using mass 1e-30
        // keeps Newton's law evaluating but the contribution to Lz drift is
        // negligible.
        let n_pad = 100;
        for i in 0..n_pad {
            let phi = 2.0 * std::f64::consts::PI * (i as f64) / (n_pad as f64);
            let r_far = 1.0e6;
            bodies.push(
                Body::rocky(1.0e-30)
                    .at(r_far * phi.cos(), r_far * phi.sin())
                    .with_velocity(0.0, 0.0),
            );
        }

        // Initial angular momentum |L_0| over the binary only.
        let l0 = orbital_angular_momentum(&bodies[0], &bodies[1]);
        let l0_mag = l0.length();

        let dt = period / 200.0;
        let n_steps = (100.0 * period / dt).ceil() as usize;

        let mut engine = BarnesHutEngine::new(16);
        let mut acc = vec![Vec3::ZERO; bodies.len()];

        // Velocity Verlet driver — minimal in-test loop, avoids pulling
        // in System orchestration just for this measurement.
        engine.build(&bodies);
        engine.evaluate(&bodies, 0.5, &mut acc);

        let mut peak_rel_drift = 0.0_f64;

        for _ in 0..n_steps {
            // kick (½dt)
            for (b, a) in bodies.iter_mut().zip(&acc) {
                b.vx += 0.5 * dt * a.x;
                b.vy += 0.5 * dt * a.y;
                b.vz += 0.5 * dt * a.z;
            }
            // drift (dt)
            for b in bodies.iter_mut() {
                b.x += dt * b.vx;
                b.y += dt * b.vy;
                b.z += dt * b.vz;
            }
            // recompute forces
            engine.build(&bodies);
            engine.evaluate(&bodies, 0.5, &mut acc);
            // kick (½dt)
            for (b, a) in bodies.iter_mut().zip(&acc) {
                b.vx += 0.5 * dt * a.x;
                b.vy += 0.5 * dt * a.y;
                b.vz += 0.5 * dt * a.z;
            }

            let l = orbital_angular_momentum(&bodies[0], &bodies[1]);
            let drift = (l - l0).length() / l0_mag;
            if drift > peak_rel_drift {
                peak_rel_drift = drift;
            }
        }

        eprintln!("[octree-tier2] inclined-kepler peak |ΔL|/|L₀| = {peak_rel_drift:.4e}");
        assert!(
            peak_rel_drift <= 1.0e-3,
            "peak |ΔL|/|L₀| = {peak_rel_drift:.4e} exceeds 1e-3 over 100 inclined Kepler periods"
        );
    }

    /// Tier 1 — Larger-N variant that genuinely exercises the BH approximation.
    ///
    /// At N = 100 the tree is shallow enough that θ = 0.5 opens most internal
    /// nodes down to leaves — the traversal effectively does exact pairwise
    /// work and the per-body error sits at the round-off floor (which meets
    /// the bound but doesn't probe the algorithm). At N = 1000 the tree
    /// reaches depth ≈ log₈(1000) ≈ 3-4 and the BH criterion accepts a
    /// meaningful number of distant nodes as monopoles. The 5 % bound is
    /// the same Salmon-Warren value; if it holds here, the algorithm
    /// approximates correctly under the load it was designed for.
    #[test]
    fn tier1_octree_bh_force_error_under_5pct_at_theta_0_5_n_1000() {
        let bodies = sphere_distribution_lognormal(1000, 0x6F637472);

        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&bodies);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&bodies, 0.5, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&bodies);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&bodies, 0.5, &mut acc_bh);

        let max_rel = body_max_rel_error(&acc_bh, &acc_exact);
        eprintln!("[octree-tier1] N=1000 θ=0.5 max rel-err = {max_rel:.4e}");
        assert!(
            max_rel <= 5e-2,
            "max per-body rel-err = {max_rel:.4e} exceeds 5e-2 (Salmon-Warren) at N=1000, θ=0.5"
        );
    }

    /// Tier 3 — Empirical wall-time scaling. Builds + evaluates the octree at
    /// a range of body counts in BH mode (θ = 0.5), reports the mean wall
    /// time per `evaluate` call after a warm-up iteration. Output goes to
    /// stderr (visible with `cargo test ... -- --nocapture`).
    ///
    /// The gate here is weak by design — absolute numbers vary by hardware
    /// and Rayon thread count — but the **growth ratio** between consecutive
    /// N values is intrinsic to the algorithm. O(N²) gives a 4× ratio when
    /// N doubles; O(N log N) gives ~2.1-2.3×. The assert at the end checks
    /// the worst observed ratio stays under 4× (i.e. better than O(N²)),
    /// which is the bare minimum for "BH is doing its job".
    ///
    /// `#[ignore]`d from the default unit-test loop because per-evaluate
    /// timings at N ∈ [100, 2500] sit in the sub-millisecond range, where
    /// run-to-run variance from OS scheduling, allocator warm-up, and CPU
    /// frequency scaling routinely pushes the worst observed ratio across
    /// the 4× gate even when the algorithm is healthy. Opt-in with
    /// `cargo test --release -p apsis tier3_octree_evaluate -- --ignored
    /// --nocapture`.
    #[test]
    #[ignore = "wall-time gate: opt-in via --ignored, run in --release on a quiet machine"]
    fn tier3_octree_evaluate_scaling_better_than_n_squared() {
        let ns = [100, 250, 500, 1000, 2500];
        let theta = 0.5;
        let warmup = 1;
        let measured = 5;

        let mut times_ms = Vec::with_capacity(ns.len());
        for &n in &ns {
            let bodies = sphere_distribution_lognormal(n, 0x6F637472);
            let mut bh = BarnesHutEngine::new(16);
            bh.set_exact_threshold(1);
            bh.build(&bodies);
            let mut acc = vec![Vec3::ZERO; bodies.len()];

            for _ in 0..warmup {
                bh.evaluate(&bodies, theta, &mut acc);
            }
            let start = std::time::Instant::now();
            for _ in 0..measured {
                bh.evaluate(&bodies, theta, &mut acc);
            }
            let mean_ms = start.elapsed().as_secs_f64() * 1000.0 / (measured as f64);
            times_ms.push(mean_ms);
            eprintln!("[octree-tier3] N={n:5} θ={theta} mean evaluate = {mean_ms:.3} ms");
        }

        // Worst growth ratio across consecutive N pairs whose N ratio is
        // approximately 2× (used to compare against the O(N²) reference of
        // 4×). Pairs in `ns` with N ratios: 250/100=2.5, 500/250=2,
        // 1000/500=2, 2500/1000=2.5 — all ~2-2.5×, consistent with each
        // other for the ratio test.
        let worst_ratio = times_ms
            .windows(2)
            .zip(ns.windows(2))
            .map(|(t, n)| {
                let n_ratio = n[1] as f64 / n[0] as f64;
                // Normalise the time ratio to a 2× N-doubling so all pairs
                // are compared on the same scale.
                let t_ratio = t[1] / t[0];
                t_ratio.powf((2.0_f64.ln()) / n_ratio.ln())
            })
            .fold(0.0_f64, f64::max);
        eprintln!("[octree-tier3] worst N-doubling time ratio = {worst_ratio:.2}× (O(N²) = 4×)");

        assert!(
            worst_ratio < 4.0,
            "worst N-doubling time ratio {worst_ratio:.2}× ≥ 4× — BH not pruning effectively, \
             traversal degraded to O(N²)-class behaviour"
        );
    }

    // ── Helpers for the validation tests ─────────────────────────────────── //

    /// Sample N bodies uniformly inside the unit sphere with log-normal
    /// masses (μ = 0, σ = 1). Deterministic for a given `seed`.
    fn sphere_distribution_lognormal(n: usize, seed: u64) -> Vec<Body> {
        // Linear congruential generator — sufficient for reproducible test
        // initial conditions; not cryptographic.
        let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut next_u64 = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            state
        };
        let mut next_unit = || (next_u64() >> 11) as f64 / (1u64 << 53) as f64;

        let mut bodies = Vec::with_capacity(n);
        while bodies.len() < n {
            // Rejection sampling for uniform-in-ball.
            let x = 2.0 * next_unit() - 1.0;
            let y = 2.0 * next_unit() - 1.0;
            let z = 2.0 * next_unit() - 1.0;
            if x * x + y * y + z * z > 1.0 {
                continue;
            }
            // Log-normal mass via Box-Muller standard normal then exp.
            let u1 = next_unit().max(1e-12);
            let u2 = next_unit();
            let normal = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            let mass = normal.exp();

            let mut b = Body::rocky(mass).at(x, y).with_velocity(0.0, 0.0);
            b.z = z;
            bodies.push(b);
        }
        bodies
    }

    /// Per-body maximum relative force error: `max_i |a_i − a_ref_i| / |a_ref_i|`,
    /// with a small absolute-magnitude floor so bodies with near-zero
    /// reference acceleration don't blow up the relative metric.
    fn body_max_rel_error(acc: &[Vec3], reference: &[Vec3]) -> f64 {
        acc.iter().zip(reference).fold(0.0_f64, |peak, (a, r)| {
            let r_mag = r.length().max(1e-30);
            let err = (*a - *r).length() / r_mag;
            peak.max(err)
        })
    }

    /// Specific angular momentum `r × v` of the relative orbit between
    /// two bodies. Returns `m_planet · (r_planet − r_central) ×
    /// (v_planet − v_central)` so the magnitude is dimensionally
    /// `mass · length² / time`.
    fn orbital_angular_momentum(central: &Body, planet: &Body) -> Vec3 {
        let r = Vec3::new(planet.x - central.x, planet.y - central.y, planet.z - central.z);
        let v = Vec3::new(planet.vx - central.vx, planet.vy - central.vy, planet.vz - central.vz);
        let cross = Vec3::new(r.y * v.z - r.z * v.y, r.z * v.x - r.x * v.z, r.x * v.y - r.y * v.x);
        planet.mass * cross
    }

    // ── MultipoleOrder toggle scaffold ────────────────────────────────── //
    //
    // Wired up at the BarnesHutEngine surface in the perf 2×2 experiment
    // (`docs/experiments/2026-05-08-octree-perf-2x2.md`). The Quadrupole
    // branch becomes physically active in the subsequent commit that adds
    // the tensor aggregation; until then both variants must produce
    // identical forces because the evaluation path is multipole-agnostic.

    #[test]
    fn multipole_order_default_is_monopole() {
        let engine = BarnesHutEngine::new(16);
        assert_eq!(engine.multipole_order(), MultipoleOrder::Monopole);
    }

    #[test]
    fn multipole_order_setter_round_trips() {
        let mut engine = BarnesHutEngine::new(16);

        engine.set_multipole_order(MultipoleOrder::Quadrupole);
        assert_eq!(engine.multipole_order(), MultipoleOrder::Quadrupole);

        engine.set_multipole_order(MultipoleOrder::Monopole);
        assert_eq!(engine.multipole_order(), MultipoleOrder::Monopole);
    }
}
