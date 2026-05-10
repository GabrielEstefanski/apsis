//! Barnes-Hut octree — pure spatial data structure, no force physics.
//!
//! ## Structure
//!
//! Flat `Vec<Node>`. Internal nodes subdivide their cubic cell into 8
//! octants; leaf nodes store up to [`LEAF`] body indices directly.
//!
//! ## Octant numbering
//!
//! Bit-packed relative to the cell centre:
//!
//! ```text
//! octant = (z >= cz) << 2 | (y >= cy) << 1 | (x >= cx)
//! ```
//!
//! Canonical Morton-like ordering (Salmon 1991; Warren & Salmon 1993).
//! Deterministic across runs, which the bit-stable replay tests rely on.
//!
//! ## Mass aggregation
//!
//! After [`Octree::build`] every node holds the aggregated total mass and
//! 3D centre of mass (`com_x`, `com_y`, `com_z`) of all bodies in its
//! subtree. These are the quantities inspected by the Barnes-Hut criterion
//! during force evaluation.
//!
//! ## Invariants
//! - `nodes[0]` is always the root after a successful `build`.
//! - A leaf has `children == [NO_CHILD; 8]`.
//! - `node.body_count` equals the sum of `body_count` of all children
//!   (or `body_len` for leaves).

use crate::domain::body::Body;

// ── Constants ─────────────────────────────────────────────────────────────── //

/// Production default for the [`Octree`] leaf-capacity generic parameter.
/// Matches GADGET-2 / PKDGRAV3 defaults; chosen at the speed end of the
/// `{4, 8, 16, 32}` Pareto trade-off characterised by the perf 2×2
/// experiment (`docs/experiments/2026-05-08-octree-perf-2x2.md`, §Results
/// leaf-sensitivity sweep).
pub(crate) const DEFAULT_LEAF: usize = 8;

/// Sentinel value for an absent child pointer.
pub(crate) const NO_CHILD: u32 = u32::MAX;

/// For N ≤ this threshold the engine falls back to exact O(N²) evaluation,
/// avoiding tree overhead that dominates at small particle counts.
pub(crate) const EXACT_THRESHOLD: usize = 64;

/// Upper clamp on the configurable exact-evaluation threshold, and the
/// canonical "direct mode" threshold. When
/// [`BarnesHutEngine::set_exact_threshold`] is called with a value
/// ≥ this constant, the engine's BH branch becomes unreachable for any
/// practical N and the force computation is guaranteed deterministic
/// (see [`BarnesHutEngine::is_direct_mode`]).
pub(crate) const DIRECT_MODE_THRESHOLD: usize = 10_000;

/// Small padding added to the root bounding cube so that no body ever sits
/// exactly on a cell boundary (which would cause ambiguous octant assignment).
const TREE_PAD: f64 = 1e-2;

// ── Node ──────────────────────────────────────────────────────────────────── //

/// One node in the Barnes-Hut octree, generic over leaf capacity.
///
/// Leaf nodes hold up to `LEAF` body indices directly. Internal nodes store
/// only aggregated mass / COM and eight child pointers. The `LEAF` generic
/// is propagated through [`Octree`] and ultimately pinned by
/// [`BarnesHutEngine`] to [`DEFAULT_LEAF`] in production; the perf 2×2
/// harness instantiates other values for the sensitivity sweep.
#[derive(Clone, Copy)]
pub(crate) struct Node<const LEAF: usize> {
    /// Cell centre, world coordinates.
    pub cx: f64,
    pub cy: f64,
    pub cz: f64,
    /// Half the side length of the cube cell (cell side = `2 * half`).
    pub half: f64,

    /// Aggregated mass of all bodies in this subtree.
    pub mass: f64,
    /// Aggregated centre-of-mass.
    pub com_x: f64,
    pub com_y: f64,
    pub com_z: f64,

    /// Symmetric traceless quadrupole tensor about this node's COM. Five
    /// independent components stored; `q_zz = -(q_xx + q_yy)` is reconstructed
    /// at the evaluation site. Populated by [`Octree::build`]'s second-pass
    /// `aggregate_quadrupole` traversal.
    pub q_xx: f64,
    pub q_xy: f64,
    pub q_xz: f64,
    pub q_yy: f64,
    pub q_yz: f64,

    /// Total body count in this subtree (leaf + all descendants).
    pub body_count: u32,

    /// Number of body indices stored in `bodies` (leaf nodes only).
    pub body_len: u8,
    /// Body index buffer for leaf nodes.  Valid range: `0..body_len`.
    pub bodies: [u32; LEAF],

    /// Child node indices in the flat node array. [`NO_CHILD`] means absent.
    /// Index follows the bit-pack `(z >= cz) << 2 | (y >= cy) << 1 | (x >= cx)`.
    pub children: [u32; 8],
}

impl<const LEAF: usize> Node<LEAF> {
    fn new(cx: f64, cy: f64, cz: f64, half: f64) -> Self {
        Self {
            cx,
            cy,
            cz,
            half,
            mass: 0.0,
            com_x: 0.0,
            com_y: 0.0,
            com_z: 0.0,
            q_xx: 0.0,
            q_xy: 0.0,
            q_xz: 0.0,
            q_yy: 0.0,
            q_yz: 0.0,
            body_count: 0,
            body_len: 0,
            bodies: [0u32; LEAF],
            children: [NO_CHILD; 8],
        }
    }

    /// True iff this node is a leaf (no child pointers set).
    #[inline]
    pub(crate) fn is_leaf(&self) -> bool {
        self.children[0] == NO_CHILD
    }

    /// Side length of this cell: `2 * half`.
    #[inline]
    pub(crate) fn size(&self) -> f64 {
        self.half * 2.0
    }
}

// ── Octree ────────────────────────────────────────────────────────────────── //

/// Flat-array Barnes-Hut octree generic over leaf capacity.
///
/// Call [`build`](Self::build) to (re)construct the tree from a body slice.
/// The resulting [`Node`] array is accessed by the force engine for both the
/// BH traversal and the `theta_error_proxy` heuristic.
///
/// `LEAF` defaults to [`DEFAULT_LEAF`] = 8, which preserves the production
/// callsite ergonomics (`Octree::new(max_depth)` continues to work unchanged
/// from PR-perf-1). The perf 2×2 leaf-sensitivity sweep instantiates other
/// values directly (`Octree::<4>::new(max_depth)` etc.) without going through
/// `BarnesHutEngine`.
pub(crate) struct Octree<const LEAF: usize = DEFAULT_LEAF> {
    pub(crate) nodes: Vec<Node<LEAF>>,
    max_depth: usize,
}

impl<const LEAF: usize> Octree<LEAF> {
    pub(crate) fn new(max_depth: usize) -> Self {
        Self { nodes: Vec::new(), max_depth }
    }

    /// Rebuild the tree from scratch for the given body slice.
    ///
    /// After this call `nodes[0]` is the root covering an axis-aligned cubic
    /// cell that contains all bodies with a small pad. Every node's `mass`,
    /// `com_{x,y,z}`, and traceless quadrupole tensor (`q_xx, q_xy, q_xz,
    /// q_yy, q_yz`) fields reflect the aggregated state of its subtree.
    /// Quadrupole aggregation is always performed — the perf 2×2 §Decision
    /// settled it as the production multipole order.
    pub(crate) fn build(&mut self, bodies: &[Body]) {
        self.nodes.clear();

        if bodies.is_empty() {
            return;
        }

        // ── Compute padded cubic AABB ────────────────────────────────── //
        let mut min_x = bodies[0].x;
        let mut max_x = bodies[0].x;
        let mut min_y = bodies[0].y;
        let mut max_y = bodies[0].y;
        let mut min_z = bodies[0].z;
        let mut max_z = bodies[0].z;
        for b in &bodies[1..] {
            min_x = min_x.min(b.x);
            max_x = max_x.max(b.x);
            min_y = min_y.min(b.y);
            max_y = max_y.max(b.y);
            min_z = min_z.min(b.z);
            max_z = max_z.max(b.z);
        }
        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let cz = 0.5 * (min_z + max_z);
        let extent = (max_x - min_x).max(max_y - min_y).max(max_z - min_z);
        let mut half = 0.5 * extent;
        half = if half <= 0.0 { TREE_PAD } else { half * 1.0001 + TREE_PAD };

        self.nodes.push(Node::new(cx, cy, cz, half));

        // ── Insert all bodies in input order ─────────────────────────── //
        for i in 0..bodies.len() {
            self.insert(0, i, bodies, 0);
        }

        // ── Aggregate mass / COM bottom-up, then quadrupole tensor ──── //
        self.aggregate_mass(0, bodies);
        self.aggregate_quadrupole(0, bodies);
    }

    /// Read-only access to the flat node array.
    #[inline]
    pub(crate) fn nodes(&self) -> &[Node<LEAF>] {
        &self.nodes
    }

    // ── Private tree-building helpers ─────────────────────────────────── //

    fn insert(&mut self, mut node_idx: usize, body_idx: usize, bodies: &[Body], mut depth: usize) {
        loop {
            // Hard depth cap: just store in current node and skip.
            if depth > self.max_depth {
                let node = &mut self.nodes[node_idx];
                if (node.body_len as usize) < LEAF {
                    node.bodies[node.body_len as usize] = body_idx as u32;
                    node.body_len += 1;
                }
                return;
            }

            if self.nodes[node_idx].is_leaf() {
                let len = self.nodes[node_idx].body_len as usize;

                // Leaf has room, or we've hit the depth cap — store here.
                if len < LEAF || depth == self.max_depth {
                    if (self.nodes[node_idx].body_len as usize) < LEAF {
                        self.nodes[node_idx].bodies[len] = body_idx as u32;
                        self.nodes[node_idx].body_len += 1;
                    }
                    return;
                }

                // Leaf is full — split into eight children, reinsert existing bodies.
                let existing_len = self.nodes[node_idx].body_len as usize;
                let existing = self.nodes[node_idx].bodies;
                self.nodes[node_idx].body_len = 0;

                self.subdivide(node_idx);

                for &bi in &existing[..existing_len] {
                    let child = self.child_octant(node_idx, bi as usize, bodies);
                    self.insert(child, bi as usize, bodies, depth + 1);
                }
            }

            node_idx = self.child_octant(node_idx, body_idx, bodies);
            depth += 1;
        }
    }

    fn subdivide(&mut self, idx: usize) {
        let (cx, cy, cz, half) = {
            let n = &self.nodes[idx];
            (n.cx, n.cy, n.cz, n.half)
        };
        let h = half * 0.5;
        // Octant order matches the bit-pack convention:
        //   bit 0 = x sign, bit 1 = y sign, bit 2 = z sign.
        // Iteration order [0..8] therefore produces:
        //   0:(−,−,−)  1:(+,−,−)  2:(−,+,−)  3:(+,+,−)
        //   4:(−,−,+)  5:(+,−,+)  6:(−,+,+)  7:(+,+,+)
        let mut children = [NO_CHILD; 8];
        for octant in 0..8 {
            let sx = if octant & 0b001 != 0 { h } else { -h };
            let sy = if octant & 0b010 != 0 { h } else { -h };
            let sz = if octant & 0b100 != 0 { h } else { -h };
            children[octant] = self.push_node(cx + sx, cy + sy, cz + sz, h) as u32;
        }
        self.nodes[idx].children = children;
    }

    fn push_node(&mut self, cx: f64, cy: f64, cz: f64, half: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node::new(cx, cy, cz, half));
        idx
    }

    /// Returns the index of the child node covering the octant that
    /// contains `bodies[body_idx]`.
    fn child_octant(&self, node_idx: usize, body_idx: usize, bodies: &[Body]) -> usize {
        let n = &self.nodes[node_idx];
        let b = bodies[body_idx];
        let octant =
            ((b.z >= n.cz) as usize) << 2 | ((b.y >= n.cy) as usize) << 1 | (b.x >= n.cx) as usize;
        self.nodes[node_idx].children[octant] as usize
    }

    /// Recursively aggregate mass and 3D centre-of-mass bottom-up.
    /// Returns `(mass, com_x, com_y, com_z)` for the subtree rooted at `idx`.
    fn aggregate_mass(&mut self, idx: usize, bodies: &[Body]) -> (f64, f64, f64, f64) {
        if self.nodes[idx].is_leaf() {
            let len = self.nodes[idx].body_len as usize;
            let mut m = 0.0_f64;
            let mut wx = 0.0_f64;
            let mut wy = 0.0_f64;
            let mut wz = 0.0_f64;

            for k in 0..len {
                let b = bodies[self.nodes[idx].bodies[k] as usize];
                m += b.mass;
                wx += b.mass * b.x;
                wy += b.mass * b.y;
                wz += b.mass * b.z;
            }

            self.nodes[idx].body_count = len as u32;
            self.nodes[idx].mass = m;

            if m > 0.0 {
                self.nodes[idx].com_x = wx / m;
                self.nodes[idx].com_y = wy / m;
                self.nodes[idx].com_z = wz / m;
                return (m, self.nodes[idx].com_x, self.nodes[idx].com_y, self.nodes[idx].com_z);
            }
            return (0.0, 0.0, 0.0, 0.0);
        }

        let children = self.nodes[idx].children;
        let mut m = 0.0_f64;
        let mut wx = 0.0_f64;
        let mut wy = 0.0_f64;
        let mut wz = 0.0_f64;
        let mut cnt = 0u32;

        for &c in &children {
            if c == NO_CHILD {
                continue;
            }
            let (cm, cx, cy, cz) = self.aggregate_mass(c as usize, bodies);
            m += cm;
            wx += cm * cx;
            wy += cm * cy;
            wz += cm * cz;
            cnt += self.nodes[c as usize].body_count;
        }

        self.nodes[idx].body_count = cnt;
        self.nodes[idx].mass = m;
        if m > 0.0 {
            self.nodes[idx].com_x = wx / m;
            self.nodes[idx].com_y = wy / m;
            self.nodes[idx].com_z = wz / m;
        }

        (self.nodes[idx].mass, self.nodes[idx].com_x, self.nodes[idx].com_y, self.nodes[idx].com_z)
    }

    /// Bottom-up second pass: populate the symmetric traceless quadrupole
    /// tensor `Q` at every node, expressed about that node's COM.
    ///
    /// Leaves: `Q = Σ_b m_b · (3 d_b ⊗ d_b − I |d_b|²)` with `d_b = pos_b − COM`.
    /// Internals: parallel-axis theorem on children — `Q = Σ_c [Q_c + M_c (3 D_c ⊗ D_c − I |D_c|²)]`
    /// with `D_c = COM_c − COM_node`. Reference: Goldstein, Poole & Safko §11.3;
    /// Hernquist & Katz 1989 eq. (2.7–2.10).
    fn aggregate_quadrupole(&mut self, idx: usize, bodies: &[Body]) {
        if self.nodes[idx].is_leaf() {
            let cmx = self.nodes[idx].com_x;
            let cmy = self.nodes[idx].com_y;
            let cmz = self.nodes[idx].com_z;
            let len = self.nodes[idx].body_len as usize;

            let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);

            for k in 0..len {
                let b = bodies[self.nodes[idx].bodies[k] as usize];
                let dx = b.x - cmx;
                let dy = b.y - cmy;
                let dz = b.z - cmz;
                let d2 = dx * dx + dy * dy + dz * dz;
                q_xx += b.mass * (3.0 * dx * dx - d2);
                q_xy += b.mass * 3.0 * dx * dy;
                q_xz += b.mass * 3.0 * dx * dz;
                q_yy += b.mass * (3.0 * dy * dy - d2);
                q_yz += b.mass * 3.0 * dy * dz;
            }

            let n = &mut self.nodes[idx];
            n.q_xx = q_xx;
            n.q_xy = q_xy;
            n.q_xz = q_xz;
            n.q_yy = q_yy;
            n.q_yz = q_yz;
            return;
        }

        let children = self.nodes[idx].children;
        for &c in &children {
            if c != NO_CHILD {
                self.aggregate_quadrupole(c as usize, bodies);
            }
        }

        let pcom_x = self.nodes[idx].com_x;
        let pcom_y = self.nodes[idx].com_y;
        let pcom_z = self.nodes[idx].com_z;

        let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);

        for &c in &children {
            if c == NO_CHILD {
                continue;
            }
            let child = &self.nodes[c as usize];
            let dx = child.com_x - pcom_x;
            let dy = child.com_y - pcom_y;
            let dz = child.com_z - pcom_z;
            let d2 = dx * dx + dy * dy + dz * dz;
            let m = child.mass;

            q_xx += child.q_xx + m * (3.0 * dx * dx - d2);
            q_xy += child.q_xy + m * 3.0 * dx * dy;
            q_xz += child.q_xz + m * 3.0 * dx * dz;
            q_yy += child.q_yy + m * (3.0 * dy * dy - d2);
            q_yz += child.q_yz + m * 3.0 * dy * dz;
        }

        let n = &mut self.nodes[idx];
        n.q_xx = q_xx;
        n.q_xy = q_xy;
        n.q_xz = q_xz;
        n.q_yy = q_yy;
        n.q_yz = q_yz;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    use approx::assert_relative_eq;
    use proptest::prelude::*;

    fn make_tree() -> Octree {
        Octree::new(16)
    }

    fn body_xyz(x: f64, y: f64, z: f64, m: f64) -> Body {
        let mut b = Body::rocky(m).at(x, y).with_velocity(0.0, 0.0);
        b.z = z;
        b
    }

    fn body_xy(x: f64, y: f64, m: f64) -> Body {
        body_xyz(x, y, 0.0, m)
    }

    /// After build, root COM must equal the mass-weighted average:
    /// r_com = Σ mᵢ rᵢ / M
    #[test]
    fn root_com_equals_mass_weighted_average() {
        let bodies = vec![body_xy(0.0, 0.0, 1.0), body_xy(4.0, 0.0, 3.0)];

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];

        assert_relative_eq!(root.mass, 4.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_x, 3.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_z, 0.0, epsilon = 1e-12);
    }

    /// Same invariant in 3D with non-zero z.
    #[test]
    fn root_com_3d_equals_mass_weighted_average() {
        let bodies = vec![body_xyz(0.0, 0.0, 0.0, 1.0), body_xyz(4.0, 2.0, -2.0, 3.0)];

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];

        assert_relative_eq!(root.mass, 4.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_x, 3.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_y, 1.5, epsilon = 1e-12);
        assert_relative_eq!(root.com_z, -1.5, epsilon = 1e-12);
    }

    /// Root must contain all bodies
    #[test]
    fn root_body_count_equals_n() {
        let bodies: Vec<Body> = (0..10).map(|i| body_xy(i as f64, 0.0, 1.0)).collect();

        let mut tree = make_tree();
        tree.build(&bodies);

        assert_eq!(tree.nodes[0].body_count, 10);
    }

    /// Single body => no subdivision
    #[test]
    fn single_body_root_is_leaf_with_no_children() {
        let bodies = vec![body_xy(1.0, 2.0, 5.0)];

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];

        assert!(root.is_leaf());
        assert_eq!(root.body_len, 1);
        assert_eq!(root.body_count, 1);
    }

    /// Edge case: multiple bodies at same position
    #[test]
    fn bodies_same_position() {
        let bodies = vec![body_xy(0.0, 0.0, 1.0), body_xy(0.0, 0.0, 2.0)];

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];

        assert_eq!(root.body_count, 2);
        assert_relative_eq!(root.com_x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(root.com_z, 0.0, epsilon = 1e-12);
    }

    /// Octant numbering must satisfy the bit-pack contract that the rest
    /// of the engine reads from (`(z >= cz) << 2 | (y >= cy) << 1 | (x >= cx)`).
    /// A regression here silently re-assigns bodies to wrong cells and
    /// poisons every BH force computation.
    #[test]
    fn octant_numbering_matches_bit_pack_contract() {
        // 16 bodies: 2 distinct positions per octant. 16 > LEAF = 8
        // forces the root to subdivide. After build, each child[octant]
        // should contain exactly the two bodies whose corner sign pattern
        // matches that octant's bit-pack index.
        let mut bodies = Vec::new();
        for octant in 0..8 {
            let sx = if octant & 0b001 != 0 { 1.0 } else { -1.0 };
            let sy = if octant & 0b010 != 0 { 1.0 } else { -1.0 };
            let sz = if octant & 0b100 != 0 { 1.0 } else { -1.0 };
            // Two bodies per octant at slightly different positions to
            // avoid same-position degeneracy.
            bodies.push(body_xyz(sx * 1.0, sy * 1.0, sz * 1.0, 1.0));
            bodies.push(body_xyz(sx * 0.5, sy * 0.5, sz * 0.5, 1.0));
        }

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert!(!root.is_leaf(), "root should subdivide after 16 inserts");
        for octant in 0..8 {
            let child_idx = root.children[octant] as usize;
            let child = &tree.nodes[child_idx];
            assert_eq!(
                child.body_count, 2,
                "octant {octant} should contain exactly its two bodies"
            );
            // The two bodies inserted at this octant were indices
            // (2*octant, 2*octant+1). Verify the child holds those two.
            let in_subtree = collect_subtree_bodies(&tree, child_idx);
            assert!(
                in_subtree.contains(&(2 * octant)) && in_subtree.contains(&(2 * octant + 1)),
                "octant {octant} subtree should hold body indices {} and {}; got {:?}",
                2 * octant,
                2 * octant + 1,
                in_subtree,
            );
        }
    }

    /// Walk a subtree and gather every body index its leaves carry. Used by
    /// octant-assignment tests to verify child placement after subdivision.
    fn collect_subtree_bodies(tree: &Octree, idx: usize) -> Vec<usize> {
        let mut out = Vec::new();
        let mut stack = vec![idx];
        while let Some(i) = stack.pop() {
            let node = &tree.nodes[i];
            if node.is_leaf() {
                for k in 0..node.body_len as usize {
                    out.push(node.bodies[k] as usize);
                }
            } else {
                for &c in &node.children {
                    if c != NO_CHILD {
                        stack.push(c as usize);
                    }
                }
            }
        }
        out
    }

    /// Tree built from a body cloud with non-zero z must aggregate the z
    /// coordinate through every level. A regression where com_z stays 0
    /// would silently make the BH branch read the wrong COM in inclined
    /// systems.
    #[test]
    fn aggregate_propagates_com_z_through_subtree() {
        // 16 bodies at random-ish 3D positions with non-zero z. Forces
        // multi-level subdivision (16 > LEAF) and exercises the
        // recursive aggregate path.
        let bodies: Vec<Body> = (0..16)
            .map(|i| {
                let t = i as f64;
                body_xyz(t * 0.7 - 5.0, (t * 1.3).sin() * 3.0, (t * 0.5).cos() * 2.0, 1.0 + 0.1 * t)
            })
            .collect();

        let m_total: f64 = bodies.iter().map(|b| b.mass).sum();
        let com_x_expected: f64 = bodies.iter().map(|b| b.mass * b.x).sum::<f64>() / m_total;
        let com_y_expected: f64 = bodies.iter().map(|b| b.mass * b.y).sum::<f64>() / m_total;
        let com_z_expected: f64 = bodies.iter().map(|b| b.mass * b.z).sum::<f64>() / m_total;

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert_relative_eq!(root.mass, m_total, epsilon = 1e-12);
        assert_relative_eq!(root.com_x, com_x_expected, epsilon = 1e-12);
        assert_relative_eq!(root.com_y, com_y_expected, epsilon = 1e-12);
        assert_relative_eq!(root.com_z, com_z_expected, epsilon = 1e-12);
    }

    proptest! {
        #[test]
        fn total_mass_is_preserved(
            xs in prop::collection::vec(-100.0..100.0, 1..50),
            ys in prop::collection::vec(-100.0..100.0, 1..50),
            zs in prop::collection::vec(-100.0..100.0, 1..50),
            masses in prop::collection::vec(0.1..10.0, 1..50)
        ) {
            let len = xs.len().min(ys.len()).min(zs.len()).min(masses.len());

            let bodies: Vec<Body> = (0..len)
                .map(|i| body_xyz(xs[i], ys[i], zs[i], masses[i]))
                .collect();

            let expected_mass: f64 = bodies.iter().map(|b| b.mass).sum();

            let mut tree = make_tree();
            tree.build(&bodies);

            let root = &tree.nodes[0];

            prop_assert!(
                (root.mass - expected_mass).abs() < 1e-6
            );
        }
    }

    // ── Quadrupole aggregation ─────────────────────────────────────────── //

    /// Direct closed-form check on a 2-body leaf (no subdivision).
    /// Two equal-mass bodies on the x-axis at ±1 give COM at the origin
    /// and `Q_xx = 4`, `Q_yy = Q_zz = −2`, off-diagonal = 0.
    #[test]
    fn quadrupole_leaf_two_bodies_matches_closed_form() {
        let bodies = vec![body_xyz(1.0, 0.0, 0.0, 1.0), body_xyz(-1.0, 0.0, 0.0, 1.0)];

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert!(root.is_leaf(), "2 bodies must fit in the root leaf");
        assert_relative_eq!(root.com_x, 0.0, epsilon = 1e-12);

        assert_relative_eq!(root.q_xx, 4.0, epsilon = 1e-12);
        assert_relative_eq!(root.q_yy, -2.0, epsilon = 1e-12);
        // Q_zz reconstructed from traceless invariant.
        let q_zz = -(root.q_xx + root.q_yy);
        assert_relative_eq!(q_zz, -2.0, epsilon = 1e-12);
        assert_relative_eq!(root.q_xy, 0.0, epsilon = 1e-12);
        assert_relative_eq!(root.q_xz, 0.0, epsilon = 1e-12);
        assert_relative_eq!(root.q_yz, 0.0, epsilon = 1e-12);
    }

    /// End-to-end invariant: regardless of how the tree subdivides, the
    /// root's aggregated `Q` (computed bottom-up via the parallel-axis
    /// theorem) must equal `Q` computed directly from all bodies relative
    /// to the root COM. Validates the leaf path, the recursive
    /// parallel-axis combination, and that the two paths agree byte-wise.
    #[test]
    fn quadrupole_root_matches_direct_sum_under_subdivision() {
        // 16 bodies arranged so the root must subdivide (LEAF = 8),
        // log-normal masses, asymmetric positions to exercise every cross
        // term (q_xy, q_xz, q_yz all nonzero).
        let positions: Vec<(f64, f64, f64, f64)> = (0..16)
            .map(|i| {
                let t = i as f64;
                let x = (t * 0.31).sin();
                let y = (t * 0.47).cos();
                let z = (t * 0.19).sin() * (t * 0.71).cos();
                let m = (t * 0.13).sin().abs() + 0.5;
                (x, y, z, m)
            })
            .collect();
        let bodies: Vec<Body> =
            positions.iter().map(|&(x, y, z, m)| body_xyz(x, y, z, m)).collect();

        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert!(!root.is_leaf(), "16 bodies must force the root to subdivide");

        let cmx = root.com_x;
        let cmy = root.com_y;
        let cmz = root.com_z;

        let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for b in &bodies {
            let dx = b.x - cmx;
            let dy = b.y - cmy;
            let dz = b.z - cmz;
            let d2 = dx * dx + dy * dy + dz * dz;
            q_xx += b.mass * (3.0 * dx * dx - d2);
            q_xy += b.mass * 3.0 * dx * dy;
            q_xz += b.mass * 3.0 * dx * dz;
            q_yy += b.mass * (3.0 * dy * dy - d2);
            q_yz += b.mass * 3.0 * dy * dz;
        }

        // Bound covers the FP-reorder drift between the bottom-up
        // recursion's accumulation order and the direct sum's order:
        // ≈ 16 leaf accumulations + log₂(16) = 4 levels of internal
        // combinations gives ~20 floating-point adds along each diagonal,
        // bounded by 20 · ε ≈ 5 × 10⁻¹⁵; 1e-12 has 200× headroom.
        assert_relative_eq!(root.q_xx, q_xx, epsilon = 1e-12);
        assert_relative_eq!(root.q_xy, q_xy, epsilon = 1e-12);
        assert_relative_eq!(root.q_xz, q_xz, epsilon = 1e-12);
        assert_relative_eq!(root.q_yy, q_yy, epsilon = 1e-12);
        assert_relative_eq!(root.q_yz, q_yz, epsilon = 1e-12);
    }

    // ── Generic LEAF parameter ─────────────────────────────────────────── //

    /// Sanity check that the generic `Octree<const LEAF: usize>` parameter
    /// actually changes leaf capacity. With LEAF = 4 a 5-body distribution
    /// must subdivide the root; with LEAF = 16 the same distribution stays
    /// in a single root leaf. Aggregated mass is independent of LEAF.
    #[test]
    fn generic_leaf_parameter_changes_split_threshold() {
        let bodies: Vec<Body> = (0..5)
            .map(|i| {
                let t = i as f64 * 0.31;
                body_xyz(t.sin(), t.cos(), (t * 0.7).sin(), 1.0)
            })
            .collect();

        let mut tight: Octree<4> = Octree::new(16);
        tight.build(&bodies);
        assert!(!tight.nodes[0].is_leaf(), "5 > LEAF=4 must subdivide the root");

        let mut loose: Octree<16> = Octree::new(16);
        loose.build(&bodies);
        assert!(loose.nodes[0].is_leaf(), "5 ≤ LEAF=16 keeps the root as a leaf");

        let total_mass: f64 = bodies.iter().map(|b| b.mass).sum();
        assert_relative_eq!(tight.nodes[0].mass, total_mass, epsilon = 1e-12);
        assert_relative_eq!(loose.nodes[0].mass, total_mass, epsilon = 1e-12);
    }

    /// Default `Octree::new` instantiates `Octree<DEFAULT_LEAF>` so existing
    /// callers (the engine, every test in this module) continue to compile
    /// unchanged. Spot-checked via `make_tree()` which still resolves to
    /// the default.
    #[test]
    fn default_octree_uses_default_leaf() {
        // If this changes silently, the production tree-build path's
        // measured costs in the lab notebook stop matching the deployed
        // binary's behaviour.
        assert_eq!(DEFAULT_LEAF, 8);
        let _tree: Octree = Octree::new(16);
    }
}
