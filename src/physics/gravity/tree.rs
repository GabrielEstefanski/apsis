//! Barnes-Hut quadtree — pure spatial data structure, no force physics.
//!
//! ## Structure
//!
//! The tree is stored as a flat `Vec<Node>`.  Each node covers a square
//! sub-region of the simulation domain.  Internal nodes subdivide into four
//! quadrants (NW, NE, SW, SE); leaf nodes store up to [`LEAF_CAPACITY`] body
//! indices directly.
//!
//! ## Mass aggregation
//!
//! After [`QuadTree::build`] every node holds the aggregated total mass and
//! center-of-mass (`com_x`, `com_y`) of all bodies in its subtree.  These are
//! the quantities inspected by the Barnes-Hut criterion during force evaluation.
//!
//! ## Invariants
//! - `nodes[0]` is always the root after a successful `build`.
//! - A leaf has `children == [NO_CHILD; 4]`.
//! - `node.body_count` equals the sum of `body_count` of all children (or
//!   `body_len` for leaves).

use crate::domain::body::Body;

// ── Constants ─────────────────────────────────────────────────────────────── //

/// Maximum number of body indices stored directly in a leaf node before it
/// is split into four children.
pub(crate) const LEAF_CAPACITY: usize = 8;

/// Sentinel value for an absent child pointer.
pub(crate) const NO_CHILD: u32 = u32::MAX;

/// For N ≤ this threshold the engine falls back to exact O(N²) evaluation,
/// avoiding tree overhead that dominates at small particle counts.
pub(crate) const EXACT_THRESHOLD: usize = 64;

/// Small padding added to the root bounding box so that no body ever sits
/// exactly on a cell boundary (which would cause ambiguous quadrant assignment).
const TREE_PAD: f64 = 1e-2;

// ── Node ──────────────────────────────────────────────────────────────────── //

/// One node in the Barnes-Hut quadtree.
///
/// Leaf nodes hold up to [`LEAF_CAPACITY`] body indices directly.
/// Internal nodes store only aggregated mass / COM and four child pointers.
#[derive(Clone, Copy)]
pub(crate) struct Node {
    /// x-coordinate of the cell centre.
    pub cx: f64,
    /// y-coordinate of the cell centre.
    pub cy: f64,
    /// Half the side length of the cell (cell side = `2 * half`).
    pub half: f64,

    /// Aggregated mass of all bodies in this subtree.
    pub mass: f64,
    /// x-coordinate of the aggregated center-of-mass.
    pub com_x: f64,
    /// y-coordinate of the aggregated center-of-mass.
    pub com_y: f64,
    /// Total body count in this subtree (leaf + all descendants).
    pub body_count: u32,

    /// Number of body indices stored in `bodies` (leaf nodes only).
    pub body_len: u8,
    /// Body index buffer for leaf nodes.  Valid range: `0..body_len`.
    pub bodies: [u32; LEAF_CAPACITY],

    /// Child node indices in the flat node array.  [`NO_CHILD`] means absent.
    /// Order: [SW, SE, NW, NE] — see [`QuadTree::child_quadrant`].
    pub children: [u32; 4],
}

impl Node {
    fn new(cx: f64, cy: f64, half: f64) -> Self {
        Self {
            cx,
            cy,
            half,
            mass: 0.0,
            com_x: 0.0,
            com_y: 0.0,
            body_count: 0,
            body_len: 0,
            bodies: [0u32; LEAF_CAPACITY],
            children: [NO_CHILD; 4],
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

// ── QuadTree ──────────────────────────────────────────────────────────────── //

/// Flat-array Barnes-Hut quadtree.
///
/// Call [`build`](Self::build) to (re)construct the tree from a body slice.
/// The resulting [`Node`] array is accessed by the force engine for both the
/// BH traversal and the `theta_error_proxy` heuristic.
pub(crate) struct QuadTree {
    pub(crate) nodes: Vec<Node>,
    max_depth: usize,
}

impl QuadTree {
    pub(crate) fn new(max_depth: usize) -> Self {
        Self {
            nodes: Vec::new(),
            max_depth,
        }
    }

    /// Rebuild the tree from scratch for the given body slice.
    ///
    /// After this call `nodes[0]` is the root covering a bounding box that
    /// contains all bodies with a small pad.  Every node's `mass` and
    /// `com_{x,y}` fields reflect the aggregated state of its subtree.
    pub(crate) fn build(&mut self, bodies: &[Body]) {
        self.nodes.clear();

        if bodies.is_empty() {
            return;
        }

        // ── Compute bounding box ───────────────────────────────────────── //
        let mut min_x = bodies[0].x;
        let mut max_x = bodies[0].x;
        let mut min_y = bodies[0].y;
        let mut max_y = bodies[0].y;

        for b in &bodies[1..] {
            min_x = min_x.min(b.x);
            max_x = max_x.max(b.x);
            min_y = min_y.min(b.y);
            max_y = max_y.max(b.y);
        }

        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let mut half = 0.5 * (max_x - min_x).max(max_y - min_y);
        half = if half <= 0.0 {
            TREE_PAD
        } else {
            half * 1.0001 + TREE_PAD
        };

        self.nodes.push(Node::new(cx, cy, half));

        // ── Insert all bodies ──────────────────────────────────────────── //
        for i in 0..bodies.len() {
            self.insert(0, i, bodies, 0);
        }

        // ── Aggregate mass / COM bottom-up ─────────────────────────────── //
        self.aggregate_mass(0, bodies);
    }

    /// Read-only access to the flat node array.
    #[inline]
    pub(crate) fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    // ── Private tree-building helpers ─────────────────────────────────── //

    fn insert(&mut self, mut node_idx: usize, body_idx: usize, bodies: &[Body], mut depth: usize) {
        loop {
            // Hard depth cap: just store in current node and skip.
            if depth > self.max_depth {
                let node = &mut self.nodes[node_idx];
                if (node.body_len as usize) < LEAF_CAPACITY {
                    node.bodies[node.body_len as usize] = body_idx as u32;
                    node.body_len += 1;
                }
                return;
            }

            if self.nodes[node_idx].is_leaf() {
                let len = self.nodes[node_idx].body_len as usize;

                // Leaf has room, or we've hit the depth cap — store here.
                if len < LEAF_CAPACITY || depth == self.max_depth {
                    if (self.nodes[node_idx].body_len as usize) < LEAF_CAPACITY {
                        self.nodes[node_idx].bodies[len] = body_idx as u32;
                        self.nodes[node_idx].body_len += 1;
                    }
                    return;
                }

                // Leaf is full — split into four children, reinsert existing bodies.
                let existing_len = self.nodes[node_idx].body_len as usize;
                let existing = self.nodes[node_idx].bodies;
                self.nodes[node_idx].body_len = 0;

                self.subdivide(node_idx);

                for &bi in &existing[..existing_len] {
                    let child = self.child_quadrant(node_idx, bi as usize, bodies);
                    self.insert(child, bi as usize, bodies, depth + 1);
                }
            }

            node_idx = self.child_quadrant(node_idx, body_idx, bodies);
            depth += 1;
        }
    }

    fn subdivide(&mut self, idx: usize) {
        let (cx, cy, half) = {
            let n = &self.nodes[idx];
            (n.cx, n.cy, n.half)
        };
        let h = half * 0.5;
        // Quadrant order: [SW, SE, NW, NE]
        let c = [
            self.push_node(cx - h, cy - h, h),
            self.push_node(cx + h, cy - h, h),
            self.push_node(cx - h, cy + h, h),
            self.push_node(cx + h, cy + h, h),
        ];
        self.nodes[idx].children = [c[0] as u32, c[1] as u32, c[2] as u32, c[3] as u32];
    }

    fn push_node(&mut self, cx: f64, cy: f64, half: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node::new(cx, cy, half));
        idx
    }

    /// Returns the index of the child node covering the quadrant that
    /// contains `bodies[body_idx]`.
    fn child_quadrant(&self, node_idx: usize, body_idx: usize, bodies: &[Body]) -> usize {
        let n = &self.nodes[node_idx];
        let b = bodies[body_idx];
        let q = match (b.x >= n.cx, b.y >= n.cy) {
            (false, false) => 0, // SW
            (true, false) => 1,  // SE
            (false, true) => 2,  // NW
            (true, true) => 3,   // NE
        };
        self.nodes[node_idx].children[q] as usize
    }

    /// Recursively aggregate mass and center-of-mass bottom-up.
    /// Returns `(mass, com_x, com_y)` for the subtree rooted at `idx`.
    fn aggregate_mass(&mut self, idx: usize, bodies: &[Body]) -> (f64, f64, f64) {
        if self.nodes[idx].is_leaf() {
            let len = self.nodes[idx].body_len as usize;
            let mut m = 0.0_f64;
            let mut wx = 0.0_f64;
            let mut wy = 0.0_f64;

            for k in 0..len {
                let b = bodies[self.nodes[idx].bodies[k] as usize];
                m += b.mass;
                wx += b.mass * b.x;
                wy += b.mass * b.y;
            }

            self.nodes[idx].body_count = len as u32;
            self.nodes[idx].mass = m;

            if m > 0.0 {
                self.nodes[idx].com_x = wx / m;
                self.nodes[idx].com_y = wy / m;
                return (m, self.nodes[idx].com_x, self.nodes[idx].com_y);
            }
            return (0.0, 0.0, 0.0);
        }

        let children = self.nodes[idx].children;
        let mut m = 0.0_f64;
        let mut wx = 0.0_f64;
        let mut wy = 0.0_f64;
        let mut cnt = 0u32;

        for &c in &children {
            if c == NO_CHILD {
                continue;
            }
            let (cm, cx, cy) = self.aggregate_mass(c as usize, bodies);
            m += cm;
            wx += cm * cx;
            wy += cm * cy;
            cnt += self.nodes[c as usize].body_count;
        }

        self.nodes[idx].body_count = cnt;
        self.nodes[idx].mass = m;
        if m > 0.0 {
            self.nodes[idx].com_x = wx / m;
            self.nodes[idx].com_y = wy / m;
        }

        (
            self.nodes[idx].mass,
            self.nodes[idx].com_x,
            self.nodes[idx].com_y,
        )
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    fn make_tree() -> QuadTree {
        QuadTree::new(16)
    }

    /// After build, root COM must equal the mass-weighted average of all bodies:
    /// r_com = Σ mᵢ rᵢ / M.
    #[test]
    fn root_com_equals_mass_weighted_average() {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 1.0),
            Body::new(4.0, 0.0, 0.0, 0.0, 3.0),
        ];
        // COM_x = (1·0 + 3·4) / 4 = 3.0
        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert!((root.mass - 4.0).abs() < 1e-12);
        assert!((root.com_x - 3.0).abs() < 1e-12);
        assert!(root.com_y.abs() < 1e-12);
    }

    /// `body_count` in the root must equal N after build — no body lost.
    #[test]
    fn root_body_count_equals_n() {
        let bodies: Vec<Body> = (0..10)
            .map(|i| Body::new(i as f64, 0.0, 0.0, 0.0, 1.0))
            .collect();
        let mut tree = make_tree();
        tree.build(&bodies);
        assert_eq!(tree.nodes[0].body_count, 10);
    }

    /// A single body produces a root that is already a leaf (no subdivision).
    #[test]
    fn single_body_root_is_leaf_with_no_children() {
        let bodies = vec![Body::new(1.0, 2.0, 0.0, 0.0, 5.0)];
        let mut tree = make_tree();
        tree.build(&bodies);

        let root = &tree.nodes[0];
        assert!(root.is_leaf());
        assert_eq!(root.body_len, 1);
        assert_eq!(root.body_count, 1);
    }
}
