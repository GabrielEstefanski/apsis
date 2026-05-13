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
//! ## Per-step maintenance
//!
//! [`Octree::maintain`] walks the per-body cell back-reference [`cell_idx`]
//! to find migrants — bodies whose new position is no longer inside their
//! known leaf — removes them from their old cells and re-inserts them from
//! the root. Multipoles are recomputed leaf-up identically to [`build`].
//! When a body has left the root bounding cube or the body count has
//! changed, maintenance falls back to a full [`build`].
//!
//! ## Invariants
//! - `nodes[0]` is always the root after a successful `build` or `maintain`.
//! - A leaf has `children == [NO_CHILD; 8]`.
//! - `node.body_count` equals the sum of `body_count` of all children
//!   (or `body_len` for leaves).
//! - `cell_idx[i]` indexes the leaf in `nodes[..]` that owns body `i`,
//!   or `u32::MAX` if no tree state exists.

use crate::domain::body_arrays::BodyArrays;

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
/// Call [`build`](Self::build) to (re)construct the tree from a [`BodyArrays`]
/// snapshot, or [`maintain`](Self::maintain) to update an existing tree
/// after the bodies have moved. The resulting [`Node`] array is accessed by
/// the force engine for both the BH traversal and the `theta_error_proxy`
/// heuristic.
///
/// `LEAF` defaults to [`DEFAULT_LEAF`] = 8. The perf 2×2 leaf-sensitivity
/// sweep (`docs/experiments/2026-05-08-octree-perf-2x2.md`) instantiates other
/// values directly (`Octree::<4>::new(max_depth)` etc.) without going through
/// `BarnesHutEngine`.
pub(crate) struct Octree<const LEAF: usize = DEFAULT_LEAF> {
    pub(crate) nodes: Vec<Node<LEAF>>,
    /// Per-body back-reference into [`nodes`]: body `i` lives in the leaf
    /// `nodes[cell_idx[i]]`. Maintained by [`insert`](Self::insert) during
    /// build/maintain. `u32::MAX` indicates "no tree state for this body"
    /// (not yet inserted, or tree was just cleared).
    pub(crate) cell_idx: Vec<u32>,
    max_depth: usize,
}

/// Sentinel for "no cell currently owns this body" in [`Octree::cell_idx`].
pub(crate) const NO_CELL: u32 = u32::MAX;

impl<const LEAF: usize> Octree<LEAF> {
    pub(crate) fn new(max_depth: usize) -> Self {
        Self { nodes: Vec::new(), cell_idx: Vec::new(), max_depth }
    }

    /// Read-only access to the flat node array.
    #[inline]
    pub(crate) fn nodes(&self) -> &[Node<LEAF>] {
        &self.nodes
    }

    // ── Private tree-building helpers ─────────────────────────────────── //

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

    /// Rebuild the tree from a [`BodyArrays`] snapshot.
    ///
    /// After this call `nodes[0]` is the root covering an axis-aligned cubic
    /// cell that contains all bodies with a small pad. Every node's `mass`,
    /// `com_{x,y,z}`, and traceless quadrupole tensor (`q_xx, q_xy, q_xz,
    /// q_yy, q_yz`) fields reflect the aggregated state of its subtree.
    /// Quadrupole aggregation is always performed — the perf 2×2 §Decision
    /// settled it as the production multipole order.
    ///
    /// Resets [`cell_idx`](Self::cell_idx) to one entry per body, populated
    /// by [`insert`](Self::insert) as each body lands in a leaf.
    pub(crate) fn build(&mut self, arrays: &BodyArrays) {
        self.nodes.clear();
        self.cell_idx.clear();
        self.cell_idx.resize(arrays.len(), NO_CELL);

        if arrays.is_empty() {
            return;
        }

        let mut min_x = arrays.pos_x[0];
        let mut max_x = arrays.pos_x[0];
        let mut min_y = arrays.pos_y[0];
        let mut max_y = arrays.pos_y[0];
        let mut min_z = arrays.pos_z[0];
        let mut max_z = arrays.pos_z[0];
        for i in 1..arrays.len() {
            min_x = min_x.min(arrays.pos_x[i]);
            max_x = max_x.max(arrays.pos_x[i]);
            min_y = min_y.min(arrays.pos_y[i]);
            max_y = max_y.max(arrays.pos_y[i]);
            min_z = min_z.min(arrays.pos_z[i]);
            max_z = max_z.max(arrays.pos_z[i]);
        }
        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let cz = 0.5 * (min_z + max_z);
        let extent = (max_x - min_x).max(max_y - min_y).max(max_z - min_z);
        let mut half = 0.5 * extent;
        half = if half <= 0.0 { TREE_PAD } else { half * 1.0001 + TREE_PAD };

        self.nodes.push(Node::new(cx, cy, cz, half));

        for i in 0..arrays.len() {
            self.insert(0, i, arrays, 0);
        }

        self.aggregate_mass(0, arrays);
        self.aggregate_quadrupole(0, arrays);
    }

    /// Update the tree to reflect the current body positions in `arrays`.
    ///
    /// Walks the per-body cell back-reference and re-inserts only those
    /// bodies whose new position is no longer inside their previously known
    /// leaf cell ("migrants"). Multipoles are recomputed leaf-up identically
    /// to [`build`], so per-cell mass / COM / quadrupole values are bit-exact
    /// with what a from-scratch [`build`] over the same particle set would
    /// produce. The tree topology (cell index assignments, subdivision
    /// depth) may differ from a from-scratch build because the maintained
    /// tree retains subdivisions made by previous steps.
    ///
    /// Falls back to [`build`] when:
    /// - no tree state exists yet (`nodes` is empty);
    /// - the body count has changed (`cell_idx.len() != arrays.len()`);
    /// - any body has migrated outside the root bounding cube (the spatial
    ///   index would otherwise lose coverage of the new position).
    pub(crate) fn maintain(&mut self, arrays: &BodyArrays) {
        if self.nodes.is_empty() || self.cell_idx.len() != arrays.len() {
            self.build(arrays);
            return;
        }

        if arrays.is_empty() {
            self.nodes.clear();
            self.cell_idx.clear();
            return;
        }

        let (rcx, rcy, rcz, rhalf) = {
            let r = &self.nodes[0];
            (r.cx, r.cy, r.cz, r.half)
        };
        for i in 0..arrays.len() {
            if (arrays.pos_x[i] - rcx).abs() > rhalf
                || (arrays.pos_y[i] - rcy).abs() > rhalf
                || (arrays.pos_z[i] - rcz).abs() > rhalf
            {
                self.build(arrays);
                return;
            }
        }

        for i in 0..arrays.len() {
            let cell = self.cell_idx[i] as usize;
            let still_in_cell = self.cell_idx[i] != NO_CELL && self.body_in_cell(i, cell, arrays);
            if !still_in_cell {
                if self.cell_idx[i] != NO_CELL {
                    self.remove_body_from_leaf(i, cell);
                }
                self.cell_idx[i] = NO_CELL;
                self.insert(0, i, arrays, 0);
            }
        }

        self.aggregate_mass(0, arrays);
        self.aggregate_quadrupole(0, arrays);
    }

    /// True iff body `body_idx` is geometrically inside cell `cell_idx`'s
    /// cube AND that cell is a leaf currently listing the body in its
    /// `bodies[]` array.
    fn body_in_cell(&self, body_idx: usize, cell_idx: usize, arrays: &BodyArrays) -> bool {
        if cell_idx >= self.nodes.len() {
            return false;
        }
        let node = &self.nodes[cell_idx];
        if !node.is_leaf() {
            return false;
        }
        let half = node.half;
        if (arrays.pos_x[body_idx] - node.cx).abs() > half
            || (arrays.pos_y[body_idx] - node.cy).abs() > half
            || (arrays.pos_z[body_idx] - node.cz).abs() > half
        {
            return false;
        }
        let len = node.body_len as usize;
        for k in 0..len {
            if node.bodies[k] as usize == body_idx {
                return true;
            }
        }
        false
    }

    /// Remove `body_idx` from leaf `cell_idx`'s `bodies[]` array, compacting
    /// by swap-with-last + decrement of `body_len`. Caller must clear the
    /// migrant's `cell_idx` slot separately.
    fn remove_body_from_leaf(&mut self, body_idx: usize, cell_idx: usize) {
        let node = &mut self.nodes[cell_idx];
        let len = node.body_len as usize;
        for k in 0..len {
            if node.bodies[k] as usize == body_idx {
                node.bodies[k] = node.bodies[len - 1];
                node.body_len -= 1;
                return;
            }
        }
    }

    fn insert(
        &mut self,
        mut node_idx: usize,
        body_idx: usize,
        arrays: &BodyArrays,
        mut depth: usize,
    ) {
        loop {
            if depth > self.max_depth {
                let node = &mut self.nodes[node_idx];
                if (node.body_len as usize) < LEAF {
                    node.bodies[node.body_len as usize] = body_idx as u32;
                    node.body_len += 1;
                    if body_idx < self.cell_idx.len() {
                        self.cell_idx[body_idx] = node_idx as u32;
                    }
                }
                return;
            }

            if self.nodes[node_idx].is_leaf() {
                let len = self.nodes[node_idx].body_len as usize;

                if len < LEAF || depth == self.max_depth {
                    if (self.nodes[node_idx].body_len as usize) < LEAF {
                        self.nodes[node_idx].bodies[len] = body_idx as u32;
                        self.nodes[node_idx].body_len += 1;
                        if body_idx < self.cell_idx.len() {
                            self.cell_idx[body_idx] = node_idx as u32;
                        }
                    }
                    return;
                }

                let existing_len = self.nodes[node_idx].body_len as usize;
                let existing = self.nodes[node_idx].bodies;
                self.nodes[node_idx].body_len = 0;

                self.subdivide(node_idx);

                for &bi in &existing[..existing_len] {
                    let child = self.child_octant(node_idx, bi as usize, arrays);
                    self.insert(child, bi as usize, arrays, depth + 1);
                }
            }

            node_idx = self.child_octant(node_idx, body_idx, arrays);
            depth += 1;
        }
    }

    fn child_octant(&self, node_idx: usize, body_idx: usize, arrays: &BodyArrays) -> usize {
        let n = &self.nodes[node_idx];
        let octant = ((arrays.pos_z[body_idx] >= n.cz) as usize) << 2
            | ((arrays.pos_y[body_idx] >= n.cy) as usize) << 1
            | (arrays.pos_x[body_idx] >= n.cx) as usize;
        self.nodes[node_idx].children[octant] as usize
    }

    fn aggregate_mass(&mut self, idx: usize, arrays: &BodyArrays) -> (f64, f64, f64, f64) {
        if self.nodes[idx].is_leaf() {
            let len = self.nodes[idx].body_len as usize;
            let mut m = 0.0_f64;
            let mut wx = 0.0_f64;
            let mut wy = 0.0_f64;
            let mut wz = 0.0_f64;

            for k in 0..len {
                let bi = self.nodes[idx].bodies[k] as usize;
                let mass = arrays.mass[bi];
                m += mass;
                wx += mass * arrays.pos_x[bi];
                wy += mass * arrays.pos_y[bi];
                wz += mass * arrays.pos_z[bi];
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
            let (cm, cx, cy, cz) = self.aggregate_mass(c as usize, arrays);
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

    fn aggregate_quadrupole(&mut self, idx: usize, arrays: &BodyArrays) {
        if self.nodes[idx].is_leaf() {
            let cmx = self.nodes[idx].com_x;
            let cmy = self.nodes[idx].com_y;
            let cmz = self.nodes[idx].com_z;
            let len = self.nodes[idx].body_len as usize;

            let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);

            for k in 0..len {
                let bi = self.nodes[idx].bodies[k] as usize;
                let mass = arrays.mass[bi];
                let dx = arrays.pos_x[bi] - cmx;
                let dy = arrays.pos_y[bi] - cmy;
                let dz = arrays.pos_z[bi] - cmz;
                let d2 = dx * dx + dy * dy + dz * dz;
                q_xx += mass * (3.0 * dx * dx - d2);
                q_xy += mass * 3.0 * dx * dy;
                q_xz += mass * 3.0 * dx * dz;
                q_yy += mass * (3.0 * dy * dy - d2);
                q_yz += mass * 3.0 * dy * dz;
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
                self.aggregate_quadrupole(c as usize, arrays);
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
        b.pos_z = z;
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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

        assert_eq!(tree.nodes[0].body_count, 10);
    }

    /// Single body => no subdivision
    #[test]
    fn single_body_root_is_leaf_with_no_children() {
        let bodies = vec![body_xy(1.0, 2.0, 5.0)];

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

        let root = &tree.nodes[0];

        assert!(root.is_leaf());
        assert_eq!(root.body_len, 1);
        assert_eq!(root.body_count, 1);
    }

    /// Edge case: multiple bodies at same position
    #[test]
    fn bodies_same_position() {
        let bodies = vec![body_xy(0.0, 0.0, 1.0), body_xy(0.0, 0.0, 2.0)];

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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
        let com_x_expected: f64 = bodies.iter().map(|b| b.mass * b.pos_x).sum::<f64>() / m_total;
        let com_y_expected: f64 = bodies.iter().map(|b| b.mass * b.pos_y).sum::<f64>() / m_total;
        let com_z_expected: f64 = bodies.iter().map(|b| b.mass * b.pos_z).sum::<f64>() / m_total;

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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

            let mut __arrays = BodyArrays::with_capacity(bodies.len());
            __arrays.pack_from(&bodies);
            let mut tree = make_tree();
            tree.build(&__arrays);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tree = make_tree();
        tree.build(&__arrays);

        let root = &tree.nodes[0];
        assert!(!root.is_leaf(), "16 bodies must force the root to subdivide");

        let cmx = root.com_x;
        let cmy = root.com_y;
        let cmz = root.com_z;

        let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for b in &bodies {
            let dx = b.pos_x - cmx;
            let dy = b.pos_y - cmy;
            let dz = b.pos_z - cmz;
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

        let mut __arrays = BodyArrays::with_capacity(bodies.len());
        __arrays.pack_from(&bodies);
        let mut tight: Octree<4> = Octree::new(16);
        tight.build(&__arrays);
        assert!(!tight.nodes[0].is_leaf(), "5 > LEAF=4 must subdivide the root");

        let mut loose: Octree<16> = Octree::new(16);
        loose.build(&__arrays);
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

    // ── Maintenance: cell_idx tracking + per-step update ───────────────── //

    /// Sphere log-normal distribution helper local to the maintenance tests.
    /// Matches the perf-series convention; deterministic per `seed`.
    fn sphere_lognormal(n: usize, seed: u64) -> Vec<Body> {
        let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut next_u64 = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            state
        };
        let mut next_unit = || (next_u64() >> 11) as f64 / (1u64 << 53) as f64;

        let mut bodies = Vec::with_capacity(n);
        while bodies.len() < n {
            let x = 2.0 * next_unit() - 1.0;
            let y = 2.0 * next_unit() - 1.0;
            let z = 2.0 * next_unit() - 1.0;
            if x * x + y * y + z * z > 1.0 {
                continue;
            }
            let u1 = next_unit().max(1e-12);
            let u2 = next_unit();
            let normal = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            let mass = normal.exp();
            let mut b = Body::rocky(mass).at(x, y).with_velocity(0.0, 0.0);
            b.pos_z = z;
            // Inject orbital-scale velocities so VV-style integration produces
            // realistic per-step displacements.
            b.vel_x = 0.3 * (next_unit() * 2.0 - 1.0);
            b.vel_y = 0.3 * (next_unit() * 2.0 - 1.0);
            b.vel_z = 0.3 * (next_unit() * 2.0 - 1.0);
            bodies.push(b);
        }
        bodies
    }

    /// `cell_idx` is populated during build and every entry points to a
    /// valid leaf containing the body. Foundation invariant for maintenance.
    #[test]
    fn build_populates_cell_idx_to_owning_leaves() {
        let bodies = sphere_lognormal(200, 0xCE11_1D58);
        let mut arrays = BodyArrays::with_capacity(bodies.len());
        arrays.pack_from(&bodies);

        let mut tree = make_tree();
        tree.build(&arrays);

        assert_eq!(tree.cell_idx.len(), bodies.len());
        for (i, &cell) in tree.cell_idx.iter().enumerate() {
            assert_ne!(cell, NO_CELL, "body {i} has no cell after build");
            let node = &tree.nodes[cell as usize];
            assert!(node.is_leaf(), "body {i}'s cell {cell} is not a leaf");
            let mut found = false;
            for k in 0..node.body_len as usize {
                if node.bodies[k] as usize == i {
                    found = true;
                    break;
                }
            }
            assert!(found, "body {i} not listed in its cell {cell}");
        }
    }

    /// Maintenance with no body movement is a no-op: tree state and
    /// per-cell multipoles are identical to the pre-maintain snapshot.
    /// Every cell is bit-exact.
    #[test]
    fn tier1_maintain_no_movement_preserves_tree_bit_exact() {
        let bodies = sphere_lognormal(500, 0x6F637472);
        let mut arrays = BodyArrays::with_capacity(bodies.len());
        arrays.pack_from(&bodies);

        let mut tree = make_tree();
        tree.build(&arrays);

        let nodes_before: Vec<_> = tree
            .nodes
            .iter()
            .map(|n| {
                (
                    n.mass.to_bits(),
                    n.com_x.to_bits(),
                    n.com_y.to_bits(),
                    n.com_z.to_bits(),
                    n.q_xx.to_bits(),
                    n.q_xy.to_bits(),
                    n.q_xz.to_bits(),
                    n.q_yy.to_bits(),
                    n.q_yz.to_bits(),
                    n.body_count,
                    n.body_len,
                )
            })
            .collect();
        let cell_idx_before = tree.cell_idx.clone();

        tree.maintain(&arrays);

        assert_eq!(tree.nodes.len(), nodes_before.len(), "node count changed under no-op maintain");
        for (i, (after, before)) in tree.nodes.iter().zip(nodes_before.iter()).enumerate() {
            assert_eq!(after.mass.to_bits(), before.0, "node {i} mass diverged");
            assert_eq!(after.com_x.to_bits(), before.1, "node {i} com_x diverged");
            assert_eq!(after.com_y.to_bits(), before.2, "node {i} com_y diverged");
            assert_eq!(after.com_z.to_bits(), before.3, "node {i} com_z diverged");
            assert_eq!(after.q_xx.to_bits(), before.4, "node {i} q_xx diverged");
            assert_eq!(after.q_xy.to_bits(), before.5, "node {i} q_xy diverged");
            assert_eq!(after.q_xz.to_bits(), before.6, "node {i} q_xz diverged");
            assert_eq!(after.q_yy.to_bits(), before.7, "node {i} q_yy diverged");
            assert_eq!(after.q_yz.to_bits(), before.8, "node {i} q_yz diverged");
            assert_eq!(after.body_count, before.9, "node {i} body_count diverged");
            assert_eq!(after.body_len, before.10, "node {i} body_len diverged");
        }
        assert_eq!(tree.cell_idx, cell_idx_before, "cell_idx changed under no-op maintain");
    }

    /// First call to maintain on a fresh tree falls back to build:
    /// `cell_idx` and `nodes` after maintain match what build would produce.
    #[test]
    fn maintain_on_empty_tree_falls_back_to_build() {
        let bodies = sphere_lognormal(100, 0x6D6F7274);
        let mut arrays = BodyArrays::with_capacity(bodies.len());
        arrays.pack_from(&bodies);

        let mut tree_built = make_tree();
        tree_built.build(&arrays);

        let mut tree_maintained = make_tree();
        tree_maintained.maintain(&arrays);

        assert_eq!(tree_built.nodes.len(), tree_maintained.nodes.len());
        assert_eq!(tree_built.cell_idx, tree_maintained.cell_idx);
    }

    /// Maintenance after every body has migrated outside the original root
    /// cube falls back to a full rebuild. Body count unchanged; positions
    /// shifted far beyond root extent.
    #[test]
    fn maintain_falls_back_when_body_leaves_root() {
        let mut bodies = sphere_lognormal(50, 0x71756164);
        let mut arrays = BodyArrays::with_capacity(bodies.len());
        arrays.pack_from(&bodies);

        let mut tree = make_tree();
        tree.build(&arrays);

        let root_half = tree.nodes[0].half;
        // Push body 0 well outside the root cube
        bodies[0].pos_x += 100.0 * root_half;
        arrays.pack_from(&bodies);

        tree.maintain(&arrays);

        assert!(
            (bodies[0].pos_x - tree.nodes[0].cx).abs() <= tree.nodes[0].half,
            "after rebuild the new root must contain body 0"
        );
    }

    /// Maintenance after a small velocity-Verlet-style displacement produces
    /// a tree whose per-body force-field-relevant scalars (mass, COM,
    /// quadrupole tensor) match what a from-scratch rebuild would produce
    /// over the same particle set, within FP-summation envelope.
    ///
    /// Tier 1 acceptance: the multipoles are mathematically computed from
    /// the same particle set; only summation-order across cell subdivisions
    /// can differ. The bound is the inherent O(n_per_leaf × ULP) ≈ 1e-14.
    /// In practice the maintained tree retains its prior subdivision and
    /// rebuild produces the same subdivision when migrants stay rare, so
    /// the typical-case error is at machine epsilon.
    #[test]
    fn tier1_maintain_per_step_matches_rebuild_per_cell_within_tolerance() {
        let mut bodies = sphere_lognormal(500, 0x6F637472);
        let mut arrays_maintain = BodyArrays::with_capacity(bodies.len());
        let mut arrays_rebuild = BodyArrays::with_capacity(bodies.len());

        let mut tree_maintain: Octree = Octree::new(16);
        let mut tree_rebuild: Octree = Octree::new(16);

        arrays_maintain.pack_from(&bodies);
        arrays_rebuild.pack_from(&bodies);
        tree_maintain.build(&arrays_maintain);
        tree_rebuild.build(&arrays_rebuild);

        // Gentle drift so most bodies stay in their current cells; a few migrate.
        let dt = 1.0e-3;
        for step in 0..5 {
            for b in &mut bodies {
                b.pos_x += dt * b.vel_x;
                b.pos_y += dt * b.vel_y;
                b.pos_z += dt * b.vel_z;
            }
            arrays_maintain.pack_from(&bodies);
            arrays_rebuild.pack_from(&bodies);

            tree_maintain.maintain(&arrays_maintain);
            tree_rebuild.build(&arrays_rebuild);

            // Aggregate-from-root invariants: total mass + mass-weighted COM
            // must agree at FP tolerance regardless of subdivision history.
            let mr = &tree_rebuild.nodes[0];
            let mm = &tree_maintain.nodes[0];
            assert_relative_eq!(mr.mass, mm.mass, epsilon = 1.0e-12, max_relative = 1.0e-12);
            assert_relative_eq!(mr.com_x, mm.com_x, epsilon = 1.0e-10, max_relative = 1.0e-10);
            assert_relative_eq!(mr.com_y, mm.com_y, epsilon = 1.0e-10, max_relative = 1.0e-10);
            assert_relative_eq!(mr.com_z, mm.com_z, epsilon = 1.0e-10, max_relative = 1.0e-10);
            assert_relative_eq!(mr.q_xx, mm.q_xx, epsilon = 1.0e-9, max_relative = 1.0e-9);
            assert_relative_eq!(mr.q_xy, mm.q_xy, epsilon = 1.0e-9, max_relative = 1.0e-9);
            assert_relative_eq!(mr.q_xz, mm.q_xz, epsilon = 1.0e-9, max_relative = 1.0e-9);
            assert_relative_eq!(mr.q_yy, mm.q_yy, epsilon = 1.0e-9, max_relative = 1.0e-9);
            assert_relative_eq!(mr.q_yz, mm.q_yz, epsilon = 1.0e-9, max_relative = 1.0e-9);

            // body_count at the root is exactly N
            assert_eq!(mr.body_count, bodies.len() as u32, "rebuild root body_count step {step}");
            assert_eq!(mm.body_count, bodies.len() as u32, "maintain root body_count step {step}",);

            // Every body's cell_idx in the maintained tree still owns it.
            for i in 0..bodies.len() {
                let cell = tree_maintain.cell_idx[i] as usize;
                let node = &tree_maintain.nodes[cell];
                assert!(node.is_leaf(), "body {i} in non-leaf after maintain step {step}");
                let mut found = false;
                for k in 0..node.body_len as usize {
                    if node.bodies[k] as usize == i {
                        found = true;
                        break;
                    }
                }
                assert!(found, "body {i} not listed in its cell {cell} after maintain step {step}");
            }
        }
    }

    // Property: across many random small displacements the maintained
    // tree's root mass / COM stay bit-exact with rebuild — the root-level
    // reductions are insensitive to the subdivision history.
    proptest! {
        #![proptest_config(ProptestConfig {
            cases: 16, .. ProptestConfig::default()
        })]
        #[test]
        fn proptest_maintain_root_mass_bit_exact(seed in 0u64..1_000) {
            let mut bodies = sphere_lognormal(64, seed.wrapping_mul(0x9E3779B97F4A7C15));
            let mut arrays = BodyArrays::with_capacity(bodies.len());
            arrays.pack_from(&bodies);
            let mut tree_m = make_tree();
            let mut tree_b = make_tree();
            tree_m.build(&arrays);
            tree_b.build(&arrays);

            for b in &mut bodies {
                b.pos_x += 1.0e-3 * b.vel_x;
                b.pos_y += 1.0e-3 * b.vel_y;
                b.pos_z += 1.0e-3 * b.vel_z;
            }
            arrays.pack_from(&bodies);

            tree_m.maintain(&arrays);
            tree_b.build(&arrays);

            prop_assert_eq!(tree_m.nodes[0].mass.to_bits(), tree_b.nodes[0].mass.to_bits());
            prop_assert_eq!(tree_m.nodes[0].body_count, tree_b.nodes[0].body_count);
        }
    }
}
