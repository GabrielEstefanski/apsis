//! Barnes-Hut force engine — orchestrates the quadtree and the Plummer kernel.
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

use crate::domain::body::Body;
use rayon::prelude::*;

use super::kernel::{G, pair_eps2, plummer_acc, plummer_phi};
use super::tree::{EXACT_THRESHOLD, NO_CHILD, Node, QuadTree};

// ── BarnesHutEngine ───────────────────────────────────────────────────────── //

/// N-body force engine using a Barnes-Hut quadtree.
///
/// Each call to [`build`](Self::build) reconstructs the quadtree from the
/// current body positions.  [`evaluate`](Self::evaluate) then computes
/// gravitational accelerations and total potential energy using that tree.
///
/// The engine contains no body state — it is safe to rebuild and re-evaluate
/// every step without any carry-over from previous steps.
pub struct BarnesHutEngine {
    tree: QuadTree,
}

impl BarnesHutEngine {
    /// Create a new engine.
    ///
    /// `max_depth` bounds the quadtree depth; 16 is sufficient for all
    /// practical particle counts.
    pub fn new(max_depth: usize) -> Self {
        Self {
            tree: QuadTree::new(max_depth),
        }
    }

    /// Rebuild the quadtree from the current body positions.
    ///
    /// Must be called before [`evaluate`](Self::evaluate) whenever bodies have moved.
    pub fn build(&mut self, bodies: &[Body]) {
        self.tree.build(bodies);
    }

    /// Compute gravitational accelerations and return total potential energy.
    ///
    /// Fills `acc[i] = (aₓ, aᵧ)` for each body.
    /// Returns `PE = Σᵢ<ⱼ −G mᵢ mⱼ / r_ij` (softened).
    ///
    /// - N ≤ [`EXACT_THRESHOLD`]: uses exact O(N²) pairwise sum.
    /// - N > [`EXACT_THRESHOLD`]: uses parallel BH traversal.
    pub fn evaluate(&self, bodies: &[Body], theta: f64, acc: &mut [(f64, f64)]) -> f64 {
        let n = bodies.len();
        acc.fill((0.0, 0.0));

        if n == 0 {
            return 0.0;
        }

        if n <= EXACT_THRESHOLD {
            return exact_eval(bodies, acc);
        }

        let nodes = self.tree.nodes();

        let results: Vec<(f64, f64, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                bh_eval_body(nodes, i, &bodies[i], bodies, theta, &mut stack)
            })
            .collect();

        let mut potential = 0.0_f64;
        for (i, (ax, ay, phi)) in results.into_iter().enumerate() {
            acc[i] = (ax, ay);
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
            let d = (dx * dx + dy * dy + eps2).sqrt();
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

        if weight_sum > 0.0 {
            (violation_sum / weight_sum).sqrt()
        } else {
            0.0
        }
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
fn exact_eval(bodies: &[Body], acc: &mut [(f64, f64)]) -> f64 {
    let n = bodies.len();
    let mut potential = 0.0_f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = bodies[j].x - bodies[i].x;
            let dy = bodies[j].y - bodies[i].y;
            let eps2 = pair_eps2(bodies[i].softening, bodies[j].softening);

            // Shared geometric factor: G / (r² + ε²)^(3/2)
            let d2 = dx * dx + dy * dy + eps2;
            let fac = G * d2.sqrt().recip().powi(3);

            // Newton's 3rd law: F_ij = −F_ji
            // acc_i += m_j · (dx, dy) · fac
            // acc_j += m_i · (−dx, −dy) · fac
            acc[i].0 += bodies[j].mass * dx * fac;
            acc[i].1 += bodies[j].mass * dy * fac;
            acc[j].0 -= bodies[i].mass * dx * fac;
            acc[j].1 -= bodies[i].mass * dy * fac;

            // Pair potential energy: E_ij = m_i · Φ_ij
            potential += bodies[i].mass * plummer_phi(dx, dy, bodies[j].mass, eps2);
        }
    }

    potential
}

/// Barnes-Hut force evaluation for a single body — O(log N) per body.
///
/// Returns `(aₓ, aᵧ, φ)` where `φ` is the specific gravitational potential
/// (potential per unit mass) at the body's position.  Multiply by `body.mass`
/// to get the contribution to total PE.
///
/// Node interactions use the target body's own ε² (the tree stores only
/// aggregated mass/COM, not per-body softening in internal nodes).
fn bh_eval_body(
    nodes: &[Node],
    body_idx: usize,
    body: &Body,
    bodies: &[Body],
    theta: f64,
    stack: &mut Vec<u32>,
) -> (f64, f64, f64) {
    let mut ax = 0.0_f64;
    let mut ay = 0.0_f64;
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
                let eps2 = pair_eps2(body.softening, other.softening);

                let (dax, day) = plummer_acc(dx, dy, other.mass, eps2);
                ax += dax;
                ay += day;
                phi += plummer_phi(dx, dy, other.mass, eps2);
            }
            continue;
        }

        // BH criterion: accept this node as a pseudo-body when s/d < θ
        let dx = node.com_x - body.x;
        let dy = node.com_y - body.y;
        let eps2 = body.softening * body.softening;
        let d = (dx * dx + dy * dy + eps2).sqrt();

        if node.size() / d < theta {
            let (dax, day) = plummer_acc(dx, dy, node.mass, eps2);
            ax += dax;
            ay += day;
            phi += plummer_phi(dx, dy, node.mass, eps2);
        } else {
            for &c in &node.children {
                if c != NO_CHILD {
                    stack.push(c);
                }
            }
        }
    }

    (ax, ay, phi)
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    fn eval(bodies: &[Body]) -> (Vec<(f64, f64)>, f64) {
        let mut engine = BarnesHutEngine::new(16);
        engine.build(bodies);
        let mut acc = vec![(0.0, 0.0); bodies.len()];
        let potential = engine.evaluate(bodies, 0.5, &mut acc);
        (acc, potential)
    }

    // ── Newton's 3rd law ────────────────────────────────────────────────── //

    /// Σᵢ mᵢ aᵢ = 0 — the total force on an isolated system is zero.
    /// This is a consequence of Newton's 3rd law: every action has an equal
    /// and opposite reaction, so internal forces sum to zero.
    #[test]
    fn total_force_on_system_is_zero() {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 1.0),
            Body::new(3.0, 0.0, 0.0, 0.0, 2.0),
        ];
        let (acc, _) = eval(&bodies);

        let fx: f64 = acc.iter().zip(&bodies).map(|(a, b)| b.mass * a.0).sum();
        let fy: f64 = acc.iter().zip(&bodies).map(|(a, b)| b.mass * a.1).sum();

        assert!(fx.abs() < 1e-12, "net Fx = {fx}");
        assert!(fy.abs() < 1e-12, "net Fy = {fy}");
    }

    // ── Force direction ─────────────────────────────────────────────────── //

    /// Gravity is attractive: each body accelerates toward the other.
    #[test]
    fn force_direction_is_attractive() {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 1.0),
            Body::new(4.0, 0.0, 0.0, 0.0, 1.0),
        ];
        let (acc, _) = eval(&bodies);
        assert!(acc[0].0 > 0.0, "b0 should accelerate toward b1 (+x)");
        assert!(acc[1].0 < 0.0, "b1 should accelerate toward b0 (−x)");
    }

    // ── Superposition ───────────────────────────────────────────────────── //

    /// A body at the midpoint between two equal masses has zero net x-force:
    /// the forces from left and right cancel exactly (superposition principle).
    #[test]
    fn symmetric_configuration_has_zero_net_x_force_on_center() {
        let bodies = vec![
            Body::new(-5.0, 0.0, 0.0, 0.0, 1.0),
            Body::new(0.0, 0.0, 0.0, 0.0, 1.0), // center
            Body::new(5.0, 0.0, 0.0, 0.0, 1.0),
        ];
        let (acc, _) = eval(&bodies);
        assert!(acc[1].0.abs() < 1e-12, "net Fx on center = {}", acc[1].0);
    }

    // ── Potential sign ──────────────────────────────────────────────────── //

    /// PE = Σᵢ<ⱼ −G mᵢ mⱼ / r_ij < 0 for any configuration of positive masses.
    #[test]
    fn gravitational_potential_is_negative() {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 1.0),
            Body::new(2.0, 0.0, 0.0, 0.0, 1.0),
        ];
        let (_, potential) = eval(&bodies);
        assert!(potential < 0.0, "PE = {potential}");
    }
}
