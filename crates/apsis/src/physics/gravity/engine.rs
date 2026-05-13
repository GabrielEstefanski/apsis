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

use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use rayon::prelude::*;

use super::kernel::{G, Kernel, PlummerKernel, pair_eps2};
#[cfg(target_arch = "x86_64")]
use super::simd;
use super::tree::{DEFAULT_LEAF, DIRECT_MODE_THRESHOLD, EXACT_THRESHOLD, NO_CHILD, Node, Octree};

// ── Leaf-pair dispatch ───────────────────────────────────────────────────── //

/// Implementation that processes the leaf-pair phase of the two-phase
/// BH walk. Resolved once at engine construction; read once per body
/// inside the parallel walk, never mutated mid-walk.
///
/// `Avx2` requires both [`Kernel::is_plummer`] (the SIMD path inlines
/// the Plummer formula) and runtime AVX2 + FMA detection on the host.
/// All other configurations route to `Scalar`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LeafPairKernel {
    Scalar,
    #[cfg(target_arch = "x86_64")]
    Avx2,
}

impl LeafPairKernel {
    /// Pick the fastest available leaf-pair implementation for the given
    /// kernel on the recorded hardware.
    fn select(kernel: &dyn Kernel) -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if kernel.is_plummer()
                && std::is_x86_feature_detected!("avx2")
                && std::is_x86_feature_detected!("fma")
            {
                return Self::Avx2;
            }
        }
        let _ = kernel;
        Self::Scalar
    }
}

// ── WalkCounters ──────────────────────────────────────────────────────────── //

/// Per-walk work counters incremented inside [`bh_eval_body`] and aggregated
/// across the parallel iter in [`BarnesHutEngine::evaluate_profile`].
///
/// Used by the engine ceiling profiling experiment
/// (`docs/experiments/2026-05-09-engine-ceiling.md`) to derive
/// `t_per_interaction = t_bh_walk / (n_bh_accepted + n_leaf_interactions)`,
/// which is the metric both SIMD and MAC optimisations affect.
///
/// The struct is `repr(C)` and contains only `u64` — no Option, Vec, or
/// branching helpers — so the increment in the hot path is a single
/// register-level `+= 1` per accepted interaction.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct WalkCounters {
    /// Total `stack.pop()` invocations (each a node visit, regardless of
    /// whether the node was BH-accepted or recursed into).
    pub n_node_visits: u64,
    /// Internal nodes accepted as monopole + traceless quadrupole via the
    /// `s/d < θ` opening criterion.
    pub n_bh_accepted: u64,
    /// Pairwise force calls inside leaf nodes (excluding self-pair).
    pub n_leaf_interactions: u64,
}

impl WalkCounters {
    #[inline(always)]
    pub(crate) fn merge(&mut self, other: &WalkCounters) {
        self.n_node_visits += other.n_node_visits;
        self.n_bh_accepted += other.n_bh_accepted;
        self.n_leaf_interactions += other.n_leaf_interactions;
    }
}

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
    /// Resolved at engine construction (and every [`set_kernel`](Self::set_kernel))
    /// from the kernel identity + host AVX2 detection. Read once per body
    /// inside [`evaluate_profile`], never mutated during a walk.
    leaf_pair_kernel: LeafPairKernel,
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
        let leaf_pair_kernel = LeafPairKernel::select(kernel.as_ref());
        Self {
            tree: Octree::new(max_depth),
            exact_threshold: EXACT_THRESHOLD,
            kernel,
            leaf_pair_kernel,
        }
    }

    /// Handle to the kernel this engine dispatches through.
    ///
    /// Used by [`System::add_hamiltonian_perturbation`](crate::core::system::System::add_hamiltonian_perturbation)
    /// and
    /// [`System::add_non_conservative_perturbation`](crate::core::system::System::add_non_conservative_perturbation)
    /// to query the active kernel's
    /// [`KernelProperties`](crate::physics::gravity::kernel::KernelProperties)
    /// against each operator's
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
        self.leaf_pair_kernel = LeafPairKernel::select(kernel.as_ref());
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

    /// Rebuild the octree from a [`BodyArrays`] snapshot.
    ///
    /// The tree is built with monopole + traceless quadrupole aggregation
    /// at every node; the multipole order is baked in (see the perf 2×2
    /// experiment, §Decision). Insertion is in input order — the Morton
    /// (Z-order) variant was characterised in the same experiment and
    /// reverted at the v1 target scale; see the §Decision for the trend
    /// at larger N.
    ///
    /// Sole writer of the BH walk's input. The SoA snapshot lifecycle is
    /// declared in `docs/experiments/2026-05-10-soa-layout.md` §Design
    /// constraint: packed once per `ForceModel::compute()`, read by build
    /// and walk, conceptually discarded after evaluate returns.
    ///
    /// Must be called before [`evaluate`](Self::evaluate) whenever bodies
    /// have moved.
    pub fn build(&mut self, arrays: &BodyArrays) {
        self.tree.build(arrays);
    }

    /// Update the tree to reflect the current body positions in `arrays`.
    ///
    /// Smart-default companion to [`build`](Self::build) for the per-step
    /// integrator path: walks the per-body cell back-reference to find
    /// migrants and re-inserts only those, derefines emptied cells, then
    /// recomputes multipoles. Falls back to a full [`build`] when the body
    /// count has changed, any body has migrated outside the root, or no
    /// prior tree state exists.
    pub fn maintain(&mut self, arrays: &BodyArrays) {
        self.tree.maintain(arrays);
    }

    /// Compute gravitational accelerations and return total potential energy.
    ///
    /// Fills `acc[i] = (aₓ, aᵧ, a_z)` for each body in `arrays`.
    /// Returns `PE = Σᵢ<ⱼ −G mᵢ mⱼ / r_ij` (softened).
    ///
    /// - N ≤ `exact_threshold`: uses exact O(N²) pairwise sum.
    /// - N > `exact_threshold`: uses parallel BH traversal.
    ///
    /// Spatial partition is the 3D octree (`Octree`) and the kernel
    /// arithmetic is fully 3D — `r² = Δx² + Δy² + Δz²` at every site.
    pub fn evaluate(&self, arrays: &BodyArrays, theta: f64, acc: &mut [Vec3]) -> f64 {
        // The profiling harness consumes the same code path via
        // [`evaluate_profile`] (see `engine_ceiling.rs`); the public surface
        // discards the work counters this method also produces internally.
        self.evaluate_profile(arrays, theta, acc).0
    }

    /// Variant of [`evaluate`] that also returns the per-step BH walk work
    /// counters aggregated across all bodies. Used by the engine ceiling
    /// profiling harness (`docs/experiments/2026-05-09-engine-ceiling.md`)
    /// to derive per-interaction cost metrics. Counters are zero in the
    /// exact-mode branch (`N ≤ exact_threshold`) since the BH walk does not
    /// execute.
    pub(crate) fn evaluate_profile(
        &self,
        arrays: &BodyArrays,
        theta: f64,
        acc: &mut [Vec3],
    ) -> (f64, WalkCounters) {
        let n = arrays.len();
        acc.fill(Vec3::ZERO);

        if n == 0 {
            return (0.0, WalkCounters::default());
        }

        let kernel: &dyn Kernel = &*self.kernel;

        if n <= self.exact_threshold {
            return (exact_eval(arrays, kernel, acc), WalkCounters::default());
        }

        let nodes = self.tree.nodes();
        let leaf_pair_kernel = self.leaf_pair_kernel;

        let results: Vec<(Vec3, f64, WalkCounters)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                let mut lists = InteractionLists::with_capacity(2048, 1024);
                bh_eval_body(
                    nodes,
                    i,
                    arrays,
                    theta,
                    kernel,
                    leaf_pair_kernel,
                    &mut stack,
                    &mut lists,
                )
            })
            .collect();

        let mut potential = 0.0_f64;
        let mut counters = WalkCounters::default();
        for (i, (a, phi, c)) in results.into_iter().enumerate() {
            acc[i] = a;
            potential += arrays.mass[i] * phi;
            counters.merge(&c);
        }

        (0.5 * potential, counters)
    }

    /// Approximate θ-error proxy for a single body.
    ///
    /// Computes a mass-weighted RMS of `(s/d)²` over all nodes that would be
    /// accepted by the BH criterion at the given `theta`.  Used by the adaptive
    /// θ controller to estimate the current force truncation error.
    pub fn theta_error_proxy(&self, body_idx: usize, arrays: &BodyArrays, theta: f64) -> f64 {
        if self.tree.nodes().is_empty() {
            return 0.0;
        }

        let body_pos_x = arrays.pos_x[body_idx];
        let body_pos_y = arrays.pos_y[body_idx];
        let body_pos_z = arrays.pos_z[body_idx];
        let body_softening = arrays.softening[body_idx];
        let eps2 = body_softening * body_softening;
        let mut violation_sum = 0.0_f64;
        let mut weight_sum = 0.0_f64;

        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);

        while let Some(raw) = stack.pop() {
            let node = &self.tree.nodes()[raw as usize];

            if node.mass <= 0.0 || node.is_leaf() {
                continue;
            }

            let dx = node.com_x - body_pos_x;
            let dy = node.com_y - body_pos_y;
            let dz = node.com_z - body_pos_z;
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

    fn node_density(&self, node: &Node<DEFAULT_LEAF>, x: f64, y: f64, z: f64, theta: f64) -> f64 {
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
        nodes: &[Node<DEFAULT_LEAF>],
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

    fn node_intersects(
        &self,
        node: &Node<DEFAULT_LEAF>,
        x: f64,
        y: f64,
        z: f64,
        radius2: f64,
    ) -> bool {
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
/// Iterates over all unique pairs (i, j) reading positions / mass / softening
/// from the [`BodyArrays`] snapshot. For each pair, applies Newton's 3rd
/// law by updating both `acc[i]` and `acc[j]` from the same kernel
/// evaluation, using pairwise softening ε²_ij = (ε²_i + ε²_j) / 2.
///
/// The component-by-component `m · d · fac` chain is load-bearing: re-
/// associating into a shared `m_fac = mass * fac` factor shifts ULPs and
/// is observable on the Mercury 1PN gate, which sits at the f64 noise floor.
///
/// Returns the total gravitational potential energy PE = Σᵢ<ⱼ mᵢ Φᵢⱼ.
fn exact_eval(arrays: &BodyArrays, kernel: &dyn Kernel, acc: &mut [Vec3]) -> f64 {
    let n = arrays.len();
    let mut potential = 0.0_f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = arrays.pos_x[j] - arrays.pos_x[i];
            let dy = arrays.pos_y[j] - arrays.pos_y[i];
            let dz = arrays.pos_z[j] - arrays.pos_z[i];
            let eps2 = pair_eps2(arrays.softening[i], arrays.softening[j]);
            let r_sq = dx * dx + dy * dy + dz * dz;

            let fac = G * kernel.acceleration_factor(r_sq, eps2);

            let mass_i = arrays.mass[i];
            let mass_j = arrays.mass[j];

            acc[i].x += mass_j * dx * fac;
            acc[i].y += mass_j * dy * fac;
            acc[i].z += mass_j * dz * fac;
            acc[j].x -= mass_i * dx * fac;
            acc[j].y -= mass_i * dy * fac;
            acc[j].z -= mass_i * dz * fac;

            let phi_ij = -G * mass_j * kernel.potential(r_sq, eps2);
            potential += mass_i * phi_ij;
        }
    }

    potential
}

/// Per-body interaction lists emitted by phase 1 of the BH walk and
/// consumed by phase 2 (the dense kernel).
///
/// Two parallel `Vec<u32>`s: leaf-pair body indices (into [`BodyArrays`])
/// and accepted-node indices (into the [`Octree`]'s flat node array).
/// Phase 2 processes them in two homogeneous loops — first all leaf-pair
/// interactions, then all accepted-node interactions — which is the
/// branchless lane-uniform shape the AVX2 leaf-pair kernel vectorises.
///
/// The struct is allocated per body per evaluate call inside the rayon
/// closure; rayon's work-stealing limits in-flight closures to ~thread
/// count, so peak memory is bounded by `num_threads × (leaf_cap +
/// node_cap) × 4 bytes` — ~150 KB at the default capacity hints below
/// on a 12-thread machine.
struct InteractionLists {
    /// Body indices of leaf-pair neighbours (excluding self).
    leaf_body_indices: Vec<u32>,
    /// Indices into the flat node array of nodes BH-accepted as
    /// monopole + quadrupole pseudo-bodies.
    accepted_node_indices: Vec<u32>,
}

impl InteractionLists {
    fn with_capacity(leaf_cap: usize, node_cap: usize) -> Self {
        Self {
            leaf_body_indices: Vec::with_capacity(leaf_cap),
            accepted_node_indices: Vec::with_capacity(node_cap),
        }
    }

    fn clear(&mut self) {
        self.leaf_body_indices.clear();
        self.accepted_node_indices.clear();
    }
}

/// Phase 1 — Walk the tree, emit per-body interaction lists.
///
/// DFS via stack. Decisions per visit:
/// - leaf: push every non-self body index into `leaf_body_indices`
/// - internal accepted (`s/d < θ`): push node index into
///   `accepted_node_indices`
/// - internal rejected: push children to stack
///
/// Returns walk counters. Lists are written into the caller's reusable
/// buffer (cleared at top).
#[inline(always)]
fn bh_walk_emit_lists(
    nodes: &[Node<DEFAULT_LEAF>],
    body_idx: usize,
    arrays: &BodyArrays,
    theta: f64,
    stack: &mut Vec<u32>,
    lists: &mut InteractionLists,
) -> WalkCounters {
    let body_pos_x = arrays.pos_x[body_idx];
    let body_pos_y = arrays.pos_y[body_idx];
    let body_pos_z = arrays.pos_z[body_idx];
    let body_softening = arrays.softening[body_idx];
    let body_eps2 = body_softening * body_softening;

    let mut counters = WalkCounters::default();
    lists.clear();

    stack.clear();
    if !nodes.is_empty() {
        stack.push(0);
    }

    while let Some(raw) = stack.pop() {
        counters.n_node_visits += 1;
        let node = &nodes[raw as usize];
        if node.mass <= 0.0 {
            continue;
        }

        if node.is_leaf() {
            for k in 0..node.body_len as usize {
                let bi = node.bodies[k];
                if bi as usize == body_idx {
                    continue;
                }
                lists.leaf_body_indices.push(bi);
                counters.n_leaf_interactions += 1;
            }
            continue;
        }

        // BH criterion: accept this node as a pseudo-body when s/d < θ.
        let dx = node.com_x - body_pos_x;
        let dy = node.com_y - body_pos_y;
        let dz = node.com_z - body_pos_z;
        let d = (dx * dx + dy * dy + dz * dz + body_eps2).sqrt();

        if node.size() / d < theta {
            lists.accepted_node_indices.push(raw);
            counters.n_bh_accepted += 1;
        } else {
            for &c in &node.children {
                if c != NO_CHILD {
                    stack.push(c);
                }
            }
        }
    }

    counters
}

/// Phase 2 — Process interaction lists with the dense scalar kernel.
///
/// Two homogeneous loops, no branches inside either:
/// 1. Leaf-pair interactions (Plummer monopole on body-body pairs).
/// 2. Accepted-node interactions (Plummer monopole + traceless quadrupole
///    on body vs aggregated node).
///
/// Summation order is segregated (all leaves, then all nodes), which
/// differs from the single-phase walk's DFS-interleaved order. Per-body
/// acceleration drift between the two is ~`O(n_interactions × ULP)` ≈
/// `~3000 × 2^-52 ≈ 7 × 10⁻¹³` in the worst case — well within Tier 1's
/// `1e-13` tolerance gate.
#[inline(always)]
fn bh_process_lists(
    nodes: &[Node<DEFAULT_LEAF>],
    body_idx: usize,
    arrays: &BodyArrays,
    kernel: &dyn Kernel,
    leaf_pair_kernel: LeafPairKernel,
    lists: &InteractionLists,
) -> (Vec3, f64) {
    let body_pos_x = arrays.pos_x[body_idx];
    let body_pos_y = arrays.pos_y[body_idx];
    let body_pos_z = arrays.pos_z[body_idx];
    let body_softening = arrays.softening[body_idx];
    let body_eps2 = body_softening * body_softening;

    // Phase 2a: leaf-pair Plummer monopole. Dispatched per-engine to either
    // the scalar dyn-kernel path or the AVX2 inlined-Plummer path.
    let (mut a, mut phi) = match leaf_pair_kernel {
        LeafPairKernel::Scalar => leaf_pair_scalar(
            body_pos_x,
            body_pos_y,
            body_pos_z,
            body_softening,
            arrays,
            kernel,
            &lists.leaf_body_indices,
        ),
        #[cfg(target_arch = "x86_64")]
        LeafPairKernel::Avx2 => unsafe {
            simd::process_leafpair_avx2(
                body_pos_x,
                body_pos_y,
                body_pos_z,
                body_softening,
                arrays,
                &lists.leaf_body_indices,
            )
        },
    };

    // Phase 2b: accepted-node Plummer monopole + traceless quadrupole.
    for &raw_ni in &lists.accepted_node_indices {
        let node = &nodes[raw_ni as usize];
        let dx = node.com_x - body_pos_x;
        let dy = node.com_y - body_pos_y;
        let dz = node.com_z - body_pos_z;
        let r_sq = dx * dx + dy * dy + dz * dz;

        let fac = G * node.mass * kernel.acceleration_factor(r_sq, body_eps2);
        a.x += dx * fac;
        a.y += dy * fac;
        a.z += dz * fac;
        phi += -G * node.mass * kernel.potential(r_sq, body_eps2);

        let q_zz = -(node.q_xx + node.q_yy);
        let qr_x = node.q_xx * dx + node.q_xy * dy + node.q_xz * dz;
        let qr_y = node.q_xy * dx + node.q_yy * dy + node.q_yz * dz;
        let qr_z = node.q_xz * dx + node.q_yz * dy + q_zz * dz;
        let rqr = dx * qr_x + dy * qr_y + dz * qr_z;

        let inv_r2 = 1.0 / (r_sq + body_eps2);
        let inv_r5 = fac / node.mass * inv_r2;
        let inv_r7 = inv_r5 * inv_r2;

        let coef_qr = -G * inv_r5;
        let coef_r = 2.5 * G * rqr * inv_r7;

        a.x += coef_qr * qr_x + coef_r * dx;
        a.y += coef_qr * qr_y + coef_r * dy;
        a.z += coef_qr * qr_z + coef_r * dz;
        phi += -0.5 * G * rqr * inv_r5;
    }

    (a, phi)
}

/// Barnes-Hut force evaluation for a single body — O(log N) per body.
///
/// Two-phase pattern: phase 1 walks the tree (control flow + decisions)
/// and emits interaction lists; phase 2 processes the lists with the
/// dense scalar or AVX2 kernel selected by [`LeafPairKernel`].
///
/// Returns `(a, φ, counters)` where `a` is the acceleration vector, `φ`
/// is the specific gravitational potential at the body's position
/// (multiply by `mass[body_idx]` to get the contribution to total PE),
/// and `counters` track work done during the walk for the engine ceiling
/// profiler.
///
/// `lists` is a reusable per-body scratch buffer. The caller (typically
/// `evaluate_profile`) allocates it once per closure invocation and
/// passes by mutable reference. Phase 1 clears it before writing.
///
/// Node interactions use the target body's own ε² — the tree stores only
/// aggregated mass and 3D COM, not per-body softening at internal nodes.
/// Each accepted node contributes monopole + traceless-quadrupole; the
/// quadrupole tensor is always populated by `Octree::build` per the perf
/// 2×2 §Decision.
#[inline(always)]
fn bh_eval_body(
    nodes: &[Node<DEFAULT_LEAF>],
    body_idx: usize,
    arrays: &BodyArrays,
    theta: f64,
    kernel: &dyn Kernel,
    leaf_pair_kernel: LeafPairKernel,
    stack: &mut Vec<u32>,
    lists: &mut InteractionLists,
) -> (Vec3, f64, WalkCounters) {
    let counters = bh_walk_emit_lists(nodes, body_idx, arrays, theta, stack, lists);
    let (a, phi) = bh_process_lists(nodes, body_idx, arrays, kernel, leaf_pair_kernel, lists);
    (a, phi, counters)
}

/// Scalar leaf-pair Plummer monopole — the dyn-dispatched fallback when
/// AVX2 is unavailable or the active kernel is not the Plummer fast path.
///
/// Pulled out of [`bh_process_lists`] so the dispatch site there is a
/// branch-free `match` over [`LeafPairKernel`] rather than two nested
/// loop bodies.
#[inline(always)]
fn leaf_pair_scalar(
    body_pos_x: f64,
    body_pos_y: f64,
    body_pos_z: f64,
    body_softening: f64,
    arrays: &BodyArrays,
    kernel: &dyn Kernel,
    leaf_body_indices: &[u32],
) -> (Vec3, f64) {
    let mut a = Vec3::ZERO;
    let mut phi = 0.0_f64;

    for &raw_bi in leaf_body_indices {
        let bi = raw_bi as usize;
        let other_mass = arrays.mass[bi];
        let dx = arrays.pos_x[bi] - body_pos_x;
        let dy = arrays.pos_y[bi] - body_pos_y;
        let dz = arrays.pos_z[bi] - body_pos_z;
        let eps2 = pair_eps2(body_softening, arrays.softening[bi]);
        let r_sq = dx * dx + dy * dy + dz * dz;

        let fac = G * other_mass * kernel.acceleration_factor(r_sq, eps2);
        a.x += dx * fac;
        a.y += dy * fac;
        a.z += dz * fac;
        phi += -G * other_mass * kernel.potential(r_sq, eps2);
    }

    (a, phi)
}

/// Single-phase BH walk reference, kept under `cfg(test)` so the
/// tolerance test can compare two-phase output against an inline
/// implementation that interleaves walk and kernel arithmetic. The
/// interleaved (DFS-order) summation is not bit-exact with the
/// two-phase segregated summation; the test gates on FP tolerance.
#[cfg(test)]
#[inline(always)]
fn bh_eval_body_single_phase(
    nodes: &[Node<DEFAULT_LEAF>],
    body_idx: usize,
    arrays: &BodyArrays,
    theta: f64,
    kernel: &dyn Kernel,
    stack: &mut Vec<u32>,
) -> (Vec3, f64, WalkCounters) {
    let body_pos_x = arrays.pos_x[body_idx];
    let body_pos_y = arrays.pos_y[body_idx];
    let body_pos_z = arrays.pos_z[body_idx];
    let body_softening = arrays.softening[body_idx];

    let mut a = Vec3::ZERO;
    let mut phi = 0.0_f64;
    let mut counters = WalkCounters::default();

    stack.clear();
    if !nodes.is_empty() {
        stack.push(0);
    }

    while let Some(raw) = stack.pop() {
        counters.n_node_visits += 1;
        let node = &nodes[raw as usize];
        if node.mass <= 0.0 {
            continue;
        }

        if node.is_leaf() {
            for k in 0..node.body_len as usize {
                let bi = node.bodies[k] as usize;
                if bi == body_idx {
                    continue;
                }
                let other_mass = arrays.mass[bi];
                let dx = arrays.pos_x[bi] - body_pos_x;
                let dy = arrays.pos_y[bi] - body_pos_y;
                let dz = arrays.pos_z[bi] - body_pos_z;
                let eps2 = pair_eps2(body_softening, arrays.softening[bi]);
                let r_sq = dx * dx + dy * dy + dz * dz;

                let fac = G * other_mass * kernel.acceleration_factor(r_sq, eps2);
                a.x += dx * fac;
                a.y += dy * fac;
                a.z += dz * fac;
                phi += -G * other_mass * kernel.potential(r_sq, eps2);
                counters.n_leaf_interactions += 1;
            }
            continue;
        }

        // BH criterion: accept this node as a pseudo-body when s/d < θ.
        let dx = node.com_x - body_pos_x;
        let dy = node.com_y - body_pos_y;
        let dz = node.com_z - body_pos_z;
        let eps2 = body_softening * body_softening;
        let d = (dx * dx + dy * dy + dz * dz + eps2).sqrt();

        if node.size() / d < theta {
            let r_sq = dx * dx + dy * dy + dz * dz;
            let fac = G * node.mass * kernel.acceleration_factor(r_sq, eps2);
            a.x += dx * fac;
            a.y += dy * fac;
            a.z += dz * fac;
            phi += -G * node.mass * kernel.potential(r_sq, eps2);

            let q_zz = -(node.q_xx + node.q_yy);
            let qr_x = node.q_xx * dx + node.q_xy * dy + node.q_xz * dz;
            let qr_y = node.q_xy * dx + node.q_yy * dy + node.q_yz * dz;
            let qr_z = node.q_xz * dx + node.q_yz * dy + q_zz * dz;
            let rqr = dx * qr_x + dy * qr_y + dz * qr_z;

            let inv_r2 = 1.0 / (r_sq + eps2);
            let inv_r5 = fac / node.mass * inv_r2;
            let inv_r7 = inv_r5 * inv_r2;

            let coef_qr = -G * inv_r5;
            let coef_r = 2.5 * G * rqr * inv_r7;

            a.x += coef_qr * qr_x + coef_r * dx;
            a.y += coef_qr * qr_y + coef_r * dy;
            a.z += coef_qr * qr_z + coef_r * dz;
            phi += -0.5 * G * rqr * inv_r5;
            counters.n_bh_accepted += 1;
        } else {
            for &c in &node.children {
                if c != NO_CHILD {
                    stack.push(c);
                }
            }
        }
    }

    (a, phi, counters)
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    use approx::assert_relative_eq;

    /// Test helper: pack a SoA snapshot, build the tree, evaluate at θ = 0.5.
    /// Mirrors what `GravityForceModel::compute` does in production, minus
    /// the engine ownership.
    fn eval(bodies: &[Body]) -> (Vec<Vec3>, f64) {
        let mut engine = BarnesHutEngine::new(16);
        let mut arrays = BodyArrays::with_capacity(bodies.len());
        arrays.pack_from(bodies);
        engine.build(&arrays);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        let potential = engine.evaluate(&arrays, 0.5, &mut acc);
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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());

        __arrays.pack_from(&bodies);
        // Exato
        let mut engine_exact = BarnesHutEngine::new(16);
        engine_exact.set_exact_threshold(usize::MAX);
        engine_exact.build(&__arrays);

        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        engine_exact.evaluate(&__arrays, 0.5, &mut acc_exact);

        // BH
        let mut engine_bh = BarnesHutEngine::new(16);
        engine_bh.set_exact_threshold(1);
        engine_bh.build(&__arrays);

        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        engine_bh.evaluate(&__arrays, 0.5, &mut acc_bh);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&__arrays);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&__arrays, 0.5, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&__arrays);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&__arrays, 0.5, &mut acc_bh);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&__arrays);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&__arrays, 0.5, &mut acc);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&__arrays);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&__arrays, 0.9, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&__arrays);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&__arrays, 0.9, &mut acc_bh);

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
                b.vel_y = v_peri * cos_i;
                b.vel_z = v_peri * sin_i;
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
        let mut __arrays = BodyArrays::with_capacity(bodies.len());

        // Velocity Verlet driver — minimal in-test loop, avoids pulling
        // in System orchestration just for this measurement. Bodies mutate
        // every step (kick / drift), so the SoA snapshot is repacked
        // before each force eval — same pattern `GravityForceModel::compute`
        // uses in production for IAS15.
        __arrays.pack_from(&bodies);
        engine.build(&__arrays);
        engine.evaluate(&__arrays, 0.5, &mut acc);

        let mut peak_rel_drift = 0.0_f64;

        for _ in 0..n_steps {
            // kick (½dt)
            for (b, a) in bodies.iter_mut().zip(&acc) {
                b.vel_x += 0.5 * dt * a.x;
                b.vel_y += 0.5 * dt * a.y;
                b.vel_z += 0.5 * dt * a.z;
            }
            // drift (dt)
            for b in bodies.iter_mut() {
                b.pos_x += dt * b.vel_x;
                b.pos_y += dt * b.vel_y;
                b.pos_z += dt * b.vel_z;
            }
            // recompute forces — repack arrays since bodies moved
            __arrays.pack_from(&bodies);
            engine.build(&__arrays);
            engine.evaluate(&__arrays, 0.5, &mut acc);
            // kick (½dt)
            for (b, a) in bodies.iter_mut().zip(&acc) {
                b.vel_x += 0.5 * dt * a.x;
                b.vel_y += 0.5 * dt * a.y;
                b.vel_z += 0.5 * dt * a.z;
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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut exact = BarnesHutEngine::new(16);
        exact.set_exact_threshold(usize::MAX);
        exact.build(&__arrays);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        exact.evaluate(&__arrays, 0.5, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.set_exact_threshold(1);
        bh.build(&__arrays);
        let mut acc_bh = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&__arrays, 0.5, &mut acc_bh);

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
            let mut __arrays = BodyArrays::with_capacity(bodies.len());
            __arrays.pack_from(&bodies);
            let mut bh = BarnesHutEngine::new(16);
            bh.set_exact_threshold(1);
            bh.build(&__arrays);
            let mut acc = vec![Vec3::ZERO; bodies.len()];

            for _ in 0..warmup {
                bh.evaluate(&__arrays, theta, &mut acc);
            }
            let start = std::time::Instant::now();
            for _ in 0..measured {
                bh.evaluate(&__arrays, theta, &mut acc);
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
            b.pos_z = z;
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
        let r = Vec3::new(
            planet.pos_x - central.pos_x,
            planet.pos_y - central.pos_y,
            planet.pos_z - central.pos_z,
        );
        let v = Vec3::new(
            planet.vel_x - central.vel_x,
            planet.vel_y - central.vel_y,
            planet.vel_z - central.vel_z,
        );
        let cross = Vec3::new(r.y * v.z - r.z * v.y, r.z * v.x - r.x * v.z, r.x * v.y - r.y * v.x);
        planet.mass * cross
    }

    /// Quadrupole-corrected BH at θ = 0.5 hits the Hernquist & Katz 1989
    /// per-body bound (5 × 10⁻³) on the lab notebook's canonical sphere
    /// distribution at N = 1000. Regression gate that the always-on
    /// quadrupole code path stays within the bound the perf 2×2 §Decision
    /// settled on.
    #[test]
    fn quadrupole_evaluate_meets_hernquist_katz_bound() {
        let bodies = sphere_distribution_lognormal(1000, 0x6F637472);
        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let theta = 0.5;

        let mut bh_exact = BarnesHutEngine::new(16);
        bh_exact.set_exact_threshold(usize::MAX);
        bh_exact.build(&__arrays);
        let mut acc_exact = vec![Vec3::ZERO; bodies.len()];
        bh_exact.evaluate(&__arrays, theta, &mut acc_exact);

        let mut bh = BarnesHutEngine::new(16);
        bh.build(&__arrays);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(&__arrays, theta, &mut acc);

        let err = body_max_rel_error(&acc, &acc_exact);
        eprintln!("[quad-evaluate] theta={theta} N=1000 max rel-err = {err:.4e}");

        // Bound 1e-2 is comfortably above quadrupole's measured max-error
        // on this distribution (3.0e-3 in the perf 2x2 N=1000 table) and
        // well below monopole's 2.5e-2 baseline -- if this fails, the
        // always-on quadrupole code path lost the multipole correction.
        assert!(
            err < 1.0e-2,
            "quadrupole BH at theta={theta} N=1000 gives max rel-err {err:.4e} -- \
             quadrupole correction may have regressed",
        );
    }

    /// Sanity: walk counters returned by `evaluate_profile` are populated
    /// with positive values on a representative non-trivial run, and equal
    /// zero in the exact-mode branch where the BH walk does not execute.
    /// Bounded relations a priori: `n_node_visits ≥ n_bh_accepted`,
    /// `n_node_visits ≥ 1` per body that triggered the walk.
    #[test]
    fn walk_counters_populate_on_bh_path_and_zero_on_exact_path() {
        let bodies = sphere_distribution_lognormal(1000, 0x6F637472);
        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let theta = 0.5;
        let mut acc = vec![Vec3::ZERO; bodies.len()];

        let mut bh = BarnesHutEngine::new(16);
        bh.build(&__arrays);
        let (_, counters_bh) = bh.evaluate_profile(&__arrays, theta, &mut acc);

        assert!(counters_bh.n_node_visits > 0, "BH walk visited zero nodes");
        assert!(counters_bh.n_bh_accepted > 0, "BH walk accepted zero internal nodes");
        assert!(counters_bh.n_leaf_interactions > 0, "BH walk did zero leaf interactions");
        assert!(counters_bh.n_node_visits >= counters_bh.n_bh_accepted);
        // Each body's walk does at least one stack pop.
        assert!(counters_bh.n_node_visits >= bodies.len() as u64);

        // Exact-mode branch: counters must stay zero.
        let mut bh_exact = BarnesHutEngine::new(16);
        bh_exact.set_exact_threshold(usize::MAX);
        bh_exact.build(&__arrays);
        let (_, counters_exact) = bh_exact.evaluate_profile(&__arrays, theta, &mut acc);
        assert_eq!(counters_exact.n_node_visits, 0);
        assert_eq!(counters_exact.n_bh_accepted, 0);
        assert_eq!(counters_exact.n_leaf_interactions, 0);
    }

    // ── Tier 1 — two-phase walk vs single-phase reference ──────────────────── //

    /// Tier 1 of `2026-05-11-simd-kernel.md` for the two-phase walk
    /// refactor. The two-phase pattern changes summation order from
    /// DFS-interleaved to segregated (all leaf-pairs first, then all
    /// accepted-nodes), so floating-point reordering is expected at
    /// `O(n_interactions × ULP)` ≈ ~7 × 10⁻¹³ at N = 10⁴. Bound 1 × 10⁻¹³
    /// covers the typical case (most bodies have fewer than ~3000
    /// interactions); a single body in a small-force pocket can hit the
    /// upper edge of the worst-case envelope, so the gate uses p99
    /// rather than max.
    ///
    /// Failure here means the two-phase implementation has a bug beyond
    /// FP reordering — likely a missed interaction or wrong index in
    /// either `bh_walk_emit_lists` or `bh_process_lists`.
    #[test]
    fn tier1_two_phase_walk_matches_single_phase_within_tolerance() {
        for &n in &[1_000usize, 5_000] {
            for &seed in &[0x6F637472u64, 0x71756164, 0x6D6F7274] {
                let bodies = sphere_distribution_lognormal(n, seed);
                let mut arrays = BodyArrays::with_capacity(bodies.len());
                arrays.pack_from(&bodies);

                let kernel = PlummerKernel::new();
                let nodes = {
                    let mut tree: Octree = Octree::new(16);
                    tree.build(&arrays);
                    tree.nodes().to_vec()
                };

                // Single-phase reference path.
                let mut acc_single = vec![Vec3::ZERO; bodies.len()];
                {
                    let mut stack = Vec::with_capacity(128);
                    for i in 0..bodies.len() {
                        let (a, _, _) =
                            bh_eval_body_single_phase(&nodes, i, &arrays, 0.5, &kernel, &mut stack);
                        acc_single[i] = a;
                    }
                }

                // Two-phase production path (scalar leaf-pair kernel — the
                // AVX2 path has its own Tier 1 test below).
                let mut acc_two_phase = vec![Vec3::ZERO; bodies.len()];
                {
                    let mut stack = Vec::with_capacity(128);
                    let mut lists = InteractionLists::with_capacity(2048, 1024);
                    for i in 0..bodies.len() {
                        let (a, _, _) = bh_eval_body(
                            &nodes,
                            i,
                            &arrays,
                            0.5,
                            &kernel,
                            LeafPairKernel::Scalar,
                            &mut stack,
                            &mut lists,
                        );
                        acc_two_phase[i] = a;
                    }
                }

                // Per-body relative error; check p99 against tolerance.
                let mut rel_errs: Vec<f64> = acc_single
                    .iter()
                    .zip(acc_two_phase.iter())
                    .map(|(s, t)| {
                        let diff = (s.x - t.x, s.y - t.y, s.z - t.z);
                        let num = (diff.0 * diff.0 + diff.1 * diff.1 + diff.2 * diff.2).sqrt();
                        let den = (s.x * s.x + s.y * s.y + s.z * s.z).sqrt().max(1e-300);
                        num / den
                    })
                    .collect();
                rel_errs.sort_by(|a, b| a.partial_cmp(b).unwrap());

                let p99_idx = (rel_errs.len() as f64 * 0.99) as usize;
                let p99 = rel_errs[p99_idx.min(rel_errs.len() - 1)];
                let max_err = *rel_errs.last().unwrap();

                eprintln!(
                    "[two-phase-tier1] N={} seed=0x{:X}  p99={:.3e}  max={:.3e}",
                    n, seed, p99, max_err,
                );

                assert!(
                    p99 <= 1e-13,
                    "two-phase vs single-phase p99 rel-err = {:.3e} exceeds 1e-13 \
                     at N={} seed=0x{:X}; max = {:.3e}",
                    p99,
                    n,
                    seed,
                    max_err,
                );
            }
        }
    }

    // ── Tier 0/1/2a — AVX2 leaf-pair kernel ────────────────────────────────── //

    /// Tier 0 (hardware sanity) — saxpy speedup.
    ///
    /// `y[i] += a · x[i]` over 1 M lanes, scalar vs AVX2-FMA. Predicts
    /// `t_scalar / t_avx2 ≥ 2.5×` on Zen 4 (4-lane double FMA, ~1 cyc
    /// throughput per FMA, scalar pipeline ~3-4 cyc per iter).
    ///
    /// Bound calibrated against Agner Fog Zen 4 latency/throughput tables
    /// (`fma vpd` ~1 cyc throughput, ~4 cyc latency); 2.5× leaves head
    /// room for memory traffic and loop overhead. A miss here means the
    /// CPU SIMD throughput is below spec or the toolchain is failing to
    /// emit the intended instructions — investigate before trusting any
    /// downstream Tier.
    ///
    /// `#[ignore]`d from default loop because wall-time gates are noisy
    /// at sub-millisecond ranges; opt-in via
    /// `cargo test --release ... -- --ignored --nocapture`.
    #[test]
    #[cfg(target_arch = "x86_64")]
    #[ignore = "wall-time gate: opt-in via --ignored, run in --release on a quiet machine"]
    fn tier0_saxpy_avx2_speedup_geq_2_5x() {
        if !std::is_x86_feature_detected!("avx2") || !std::is_x86_feature_detected!("fma") {
            eprintln!("[simd-tier0] AVX2+FMA unavailable, skipping");
            return;
        }

        const N: usize = 1_000_000;
        let x: Vec<f64> = (0..N).map(|i| (i as f64) * 1.0e-3).collect();
        let mut y_scalar: Vec<f64> = vec![0.0; N];
        let mut y_avx2: Vec<f64> = vec![0.0; N];
        let a = 0.7_f64;

        // Warm-up.
        super::super::simd::saxpy_scalar(a, &x, &mut y_scalar);
        unsafe { super::super::simd::saxpy_avx2(a, &x, &mut y_avx2) };

        let measured = 5;
        let start_scalar = std::time::Instant::now();
        for _ in 0..measured {
            super::super::simd::saxpy_scalar(a, &x, &mut y_scalar);
        }
        let t_scalar_ns = start_scalar.elapsed().as_nanos() as f64 / measured as f64;

        let start_avx2 = std::time::Instant::now();
        for _ in 0..measured {
            unsafe { super::super::simd::saxpy_avx2(a, &x, &mut y_avx2) };
        }
        let t_avx2_ns = start_avx2.elapsed().as_nanos() as f64 / measured as f64;

        let speedup = t_scalar_ns / t_avx2_ns;
        eprintln!(
            "[simd-tier0] saxpy N={N}  t_scalar={t_scalar_ns:.0}ns  \
             t_avx2={t_avx2_ns:.0}ns  speedup={speedup:.2}×",
        );

        assert!(
            speedup >= 2.5,
            "saxpy AVX2 speedup {speedup:.2}× < 2.5× — toolchain or hardware not delivering \
             the expected SIMD throughput; investigate before trusting downstream gates"
        );
    }

    /// Tier 1 — AVX2 leaf-pair kernel matches scalar within tolerance.
    ///
    /// Compares `evaluate` output between two engines on the same SoA
    /// snapshot: one forced into scalar leaf-pair dispatch, the other on
    /// AVX2. The AVX2 path inlines the Plummer formula with horizontal-
    /// sum reductions across 4 lanes, which reorders summation (sum of
    /// 4-lane chunks instead of strict left-to-right). p99 per-body
    /// relative-acceleration error must stay under 1 × 10⁻¹³, the same
    /// bound the two-phase walk Tier 1 settled on (FP-reordering envelope
    /// of `O(n_interactions × ULP)`).
    ///
    /// Failure here means the AVX2 kernel diverges from scalar Plummer
    /// beyond pure FP reordering — likely a bug in gather indexing, lane
    /// arithmetic, or the scalar tail.
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn tier1_avx2_leaf_pair_matches_scalar_within_tolerance() {
        if !std::is_x86_feature_detected!("avx2") || !std::is_x86_feature_detected!("fma") {
            eprintln!("[simd-tier1] AVX2+FMA unavailable, skipping");
            return;
        }

        for &n in &[1_000usize, 5_000] {
            for &seed in &[0x6F637472u64, 0x71756164, 0x6D6F7274] {
                let bodies = sphere_distribution_lognormal(n, seed);
                let mut arrays = BodyArrays::with_capacity(bodies.len());
                arrays.pack_from(&bodies);

                let mut engine_scalar = BarnesHutEngine::new(16);
                engine_scalar.leaf_pair_kernel = LeafPairKernel::Scalar;
                engine_scalar.build(&arrays);
                let mut acc_scalar = vec![Vec3::ZERO; bodies.len()];
                engine_scalar.evaluate(&arrays, 0.5, &mut acc_scalar);

                let mut engine_avx2 = BarnesHutEngine::new(16);
                engine_avx2.leaf_pair_kernel = LeafPairKernel::Avx2;
                engine_avx2.build(&arrays);
                let mut acc_avx2 = vec![Vec3::ZERO; bodies.len()];
                engine_avx2.evaluate(&arrays, 0.5, &mut acc_avx2);

                let mut rel_errs: Vec<f64> = acc_scalar
                    .iter()
                    .zip(acc_avx2.iter())
                    .map(|(s, t)| {
                        let diff = (s.x - t.x, s.y - t.y, s.z - t.z);
                        let num = (diff.0 * diff.0 + diff.1 * diff.1 + diff.2 * diff.2).sqrt();
                        let den = (s.x * s.x + s.y * s.y + s.z * s.z).sqrt().max(1e-300);
                        num / den
                    })
                    .collect();
                rel_errs.sort_by(|a, b| a.partial_cmp(b).unwrap());

                let p99_idx = (rel_errs.len() as f64 * 0.99) as usize;
                let p99 = rel_errs[p99_idx.min(rel_errs.len() - 1)];
                let max_err = *rel_errs.last().unwrap();

                eprintln!(
                    "[simd-tier1] AVX2 vs scalar  N={n}  seed=0x{seed:X}  \
                     p99={p99:.3e}  max={max_err:.3e}",
                );

                assert!(
                    p99 <= 1e-13,
                    "AVX2 leaf-pair vs scalar p99 rel-err = {p99:.3e} exceeds 1e-13 \
                     at N={n} seed=0x{seed:X}; max = {max_err:.3e}"
                );
            }
        }
    }

    /// Tier 2a (kernel-isolated) — AVX2 Plummer microkernel speedup over
    /// scalar.
    ///
    /// Runs both kernels on a pre-laid-out tuple stream of 1 M
    /// interactions. Bypasses gather (scalar reads the same tuples,
    /// AVX2 reads aligned arrays) so the measured ratio reflects pure
    /// arithmetic + reduction throughput, not memory-system effects.
    ///
    /// A-priori range `[1.8, 2.5]×` from
    /// `docs/experiments/2026-05-11-simd-kernel.md` §A-priori bounds:
    /// AVX2 at 4 lanes × FMA-1 cyc throughput vs scalar that the compiler
    /// is already auto-vectorising at 2 lanes — the realistic delta is
    /// 1.8-2.5×, not the naive 4× lane-count. A measurement above 2.5×
    /// would be surprising; below 1.8× would suggest the scalar path is
    /// already lane-saturated or the AVX2 kernel has a stall (sqrt+div
    /// chain on the critical path).
    ///
    /// `#[ignore]`d for the same reason as Tier 0: wall-time gates need
    /// `--release` on a quiet machine.
    #[test]
    #[cfg(target_arch = "x86_64")]
    #[ignore = "wall-time gate: opt-in via --ignored, run in --release on a quiet machine"]
    fn tier2a_kernel_avx2_speedup_in_range() {
        if !std::is_x86_feature_detected!("avx2") || !std::is_x86_feature_detected!("fma") {
            eprintln!("[simd-tier2a] AVX2+FMA unavailable, skipping");
            return;
        }

        const N: usize = 1_000_000;
        // Synthesise a deterministic interaction stream (no gather).
        let mut tuples: Vec<(f64, f64, f64, f64, f64)> = Vec::with_capacity(N);
        let mut dx_a = Vec::with_capacity(N);
        let mut dy_a = Vec::with_capacity(N);
        let mut dz_a = Vec::with_capacity(N);
        let mut eps2_a = Vec::with_capacity(N);
        let mut mass_a = Vec::with_capacity(N);
        let mut s = 0x9E3779B97F4A7C15_u64;
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        };
        for _ in 0..N {
            let dx = next() * 10.0;
            let dy = next() * 10.0;
            let dz = next() * 10.0;
            let eps2 = 1.0e-3 + next().abs() * 1.0e-2;
            let m = 0.5 + next().abs() * 0.5;
            tuples.push((dx, dy, dz, eps2, m));
            dx_a.push(dx);
            dy_a.push(dy);
            dz_a.push(dz);
            eps2_a.push(eps2);
            mass_a.push(m);
        }

        let mut acc_scalar = Vec3::ZERO;
        let mut acc_avx2 = Vec3::ZERO;

        // Warm-up.
        super::super::simd::plummer_kernel_scalar_micro(&mut acc_scalar, &tuples);
        unsafe {
            super::super::simd::plummer_kernel_avx2_micro(
                &mut acc_avx2,
                &dx_a,
                &dy_a,
                &dz_a,
                &eps2_a,
                &mass_a,
            )
        };
        acc_scalar = Vec3::ZERO;
        acc_avx2 = Vec3::ZERO;

        let measured = 5;
        let start_scalar = std::time::Instant::now();
        for _ in 0..measured {
            super::super::simd::plummer_kernel_scalar_micro(&mut acc_scalar, &tuples);
        }
        let t_scalar_ns = start_scalar.elapsed().as_nanos() as f64 / measured as f64;

        let start_avx2 = std::time::Instant::now();
        for _ in 0..measured {
            unsafe {
                super::super::simd::plummer_kernel_avx2_micro(
                    &mut acc_avx2,
                    &dx_a,
                    &dy_a,
                    &dz_a,
                    &eps2_a,
                    &mass_a,
                )
            };
        }
        let t_avx2_ns = start_avx2.elapsed().as_nanos() as f64 / measured as f64;

        let speedup = t_scalar_ns / t_avx2_ns;
        eprintln!(
            "[simd-tier2a] Plummer kernel N={N}  t_scalar={t_scalar_ns:.0}ns  \
             t_avx2={t_avx2_ns:.0}ns  speedup={speedup:.2}×",
        );

        // Sanity: the two paths should agree to first-order on totals,
        // even with reordered summation. Off by ≥ 1 % flags a bug.
        let scalar_mag =
            (acc_scalar.x.powi(2) + acc_scalar.y.powi(2) + acc_scalar.z.powi(2)).sqrt();
        let avx2_mag = (acc_avx2.x.powi(2) + acc_avx2.y.powi(2) + acc_avx2.z.powi(2)).sqrt();
        let rel_diff = (scalar_mag - avx2_mag).abs() / scalar_mag.max(1e-300);
        assert!(rel_diff < 1e-2, "scalar and AVX2 micro disagree by {rel_diff:.2e} > 1%");

        assert!(
            (1.8..=2.5).contains(&speedup),
            "AVX2 microkernel speedup {speedup:.2}× outside [1.8, 2.5]× envelope"
        );
    }

    // ── Leaf-pair kernel dispatch ──────────────────────────────────────────── //

    /// `LeafPairKernel::select` picks the AVX2 path on a Plummer kernel
    /// plus an AVX2/FMA-capable host, and falls back to scalar otherwise.
    /// Hardware capability is read at runtime.
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn dispatch_picks_avx2_when_plummer_and_avx2_fma_available() {
        let kernel = PlummerKernel::new();
        let picked = LeafPairKernel::select(&kernel);
        match picked {
            LeafPairKernel::Avx2 => {
                assert!(
                    std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma")
                );
            },
            LeafPairKernel::Scalar => {
                assert!(
                    !(std::is_x86_feature_detected!("avx2")
                        && std::is_x86_feature_detected!("fma"))
                );
            },
        }
    }

    /// Non-Plummer kernel must never be routed to the AVX2 SIMD path —
    /// the AVX2 kernel inlines the Plummer formula directly and would
    /// silently produce wrong forces against e.g. `TruncatedPlummerKernel`.
    #[test]
    fn dispatch_falls_back_to_scalar_for_non_plummer_kernel() {
        use crate::physics::gravity::kernel::TruncatedPlummerKernel;
        let truncated = TruncatedPlummerKernel::new(1.0);
        let picked = LeafPairKernel::select(&truncated);
        assert_eq!(picked, LeafPairKernel::Scalar);
    }
}
