//! Morton (Z-order) spatial encoding for body insertion + walk ordering.
//!
//! ## Why Morton
//!
//! In a Barnes-Hut walk, two bodies that are spatially adjacent visit
//! overlapping subsets of the tree. If consecutive bodies in the parallel
//! iterator are also spatially adjacent, the second body's walk finds the
//! shared tree nodes warm in cache from the first body's walk. Morton
//! ordering — sorting bodies by interleaved 3D coordinate bits — is the
//! standard way to achieve that property without a heuristic.
//!
//! ## Bit layout
//!
//! 21 bits per axis (`MORTON_BITS`), interleaved into a 63-bit word with
//! the same bit-pack convention the octree uses for child indexing
//! (`(z >= cz) << 2 | (y >= cy) << 1 | (x >= cx)`). The k-th triple of
//! Morton bits therefore selects the same octant the tree's k-th
//! subdivision visits, so high-order bits cluster bodies by root octant,
//! next bits cluster within child octants, and so on.
//!
//! 21 bits gives 2²¹ ≈ 2.1 M cells per axis. The maximum tree depth of 16
//! consumes only 16 bits, so the encoding has 5-level headroom: bodies at
//! the deepest tree level still encode to distinct Morton codes unless
//! their coordinates collide at sub-cell precision.
//!
//! ## References
//! - Morton (1966). *A computer oriented geodetic data base*. IBM Tech. Rep.
//! - Salmon (1991). *Parallel hierarchical N-body methods*. Caltech PhD.
//! - Warren & Salmon (1993). *A parallel hashed Oct-Tree N-body algorithm*.

// `dead_code` allowed for the scaffolding window. The module's items are
// exercised by tests (below) and by the perf 2×2 harness, but the next
// commit (Morton-aware build in `Octree::build`) is what brings them into
// the lib's non-test code path; the allow is removed there.
#![allow(dead_code)]

use crate::domain::body::Body;

use super::tree::TREE_PAD;

// ── Aabb ─────────────────────────────────────────────────────────────────── //

/// Padded cubic axis-aligned bounding box used by both the octree's root
/// cell and Morton normalisation. Sharing the same struct guarantees the
/// two views of "the simulation domain" never drift.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Aabb {
    pub center: [f64; 3],
    /// Half the side length of the cubic cell (full side = `2 * half`).
    pub half: f64,
}

/// Compute the padded cubic AABB enclosing every body, using the same
/// padding rule the octree applies to its root cell.
///
/// Returns an arbitrary unit-half AABB at the origin if `bodies` is empty
/// — caller is responsible for skipping the empty case if it matters.
pub(crate) fn compute_aabb(bodies: &[Body]) -> Aabb {
    if bodies.is_empty() {
        return Aabb { center: [0.0; 3], half: 1.0 };
    }
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
    let center = [0.5 * (min_x + max_x), 0.5 * (min_y + max_y), 0.5 * (min_z + max_z)];
    let extent = (max_x - min_x).max(max_y - min_y).max(max_z - min_z);
    let mut half = 0.5 * extent;
    half = if half <= 0.0 { TREE_PAD } else { half * 1.0001 + TREE_PAD };
    Aabb { center, half }
}

// ── Morton encoding ──────────────────────────────────────────────────────── //

/// Bits per axis. 21 × 3 = 63, fits in `u64` with one bit to spare.
const MORTON_BITS: u32 = 21;

/// Spread a 21-bit integer across every third bit position of a 63-bit
/// word. Output bit positions are 0, 3, 6, …, 60 (every third bit, low to
/// high). Standard "magic-number" sequence; widely cited inline
/// implementation (Holcomb 2010; rust-numerics-style variants).
#[inline]
fn expand_bits_21(v: u64) -> u64 {
    let mut x = v & 0x1F_FFFF;
    x = (x | (x << 32)) & 0x001F_0000_0000_FFFF;
    x = (x | (x << 16)) & 0x001F_0000_FF00_00FF;
    x = (x | (x << 8)) & 0x100F_00F0_0F00_F00F;
    x = (x | (x << 4)) & 0x10C3_0C30_C30C_30C3;
    x = (x | (x << 2)) & 0x1249_2492_4924_9249;
    x
}

/// 63-bit Morton code interleaving three 21-bit quantised coordinates.
/// Bit ordering matches the octree's child-index convention so the k-th
/// triple of code bits gives the octant the k-th subdivision would pick.
#[inline]
pub(crate) fn morton_encode_3d(qx: u64, qy: u64, qz: u64) -> u64 {
    expand_bits_21(qx) | (expand_bits_21(qy) << 1) | (expand_bits_21(qz) << 2)
}

// ── Permutation ──────────────────────────────────────────────────────────── //

/// Compute the permutation of body indices that orders bodies by ascending
/// Morton code in the supplied AABB.
///
/// Bodies sharing a Morton code (only possible when both fall in the same
/// 2²¹ sub-cell along every axis) are stably ordered by their original
/// index, which keeps determinism intact for the bit-stable replay tests.
pub(crate) fn compute_morton_permutation(bodies: &[Body], aabb: &Aabb) -> Vec<u32> {
    if bodies.is_empty() {
        return Vec::new();
    }
    let inv_extent = 1.0 / (2.0 * aabb.half);
    let q_max = (1u64 << MORTON_BITS) - 1;
    let q_max_f = q_max as f64;

    let mut codes: Vec<(u64, u32)> = bodies
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let nx = ((b.x - aabb.center[0]) * inv_extent + 0.5).clamp(0.0, 1.0);
            let ny = ((b.y - aabb.center[1]) * inv_extent + 0.5).clamp(0.0, 1.0);
            let nz = ((b.z - aabb.center[2]) * inv_extent + 0.5).clamp(0.0, 1.0);
            let qx = (nx * q_max_f) as u64;
            let qy = (ny * q_max_f) as u64;
            let qz = (nz * q_max_f) as u64;
            (morton_encode_3d(qx, qy, qz), i as u32)
        })
        .collect();

    codes.sort_by_key(|&(code, _)| code);
    codes.into_iter().map(|(_, idx)| idx).collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;

    fn body_at(x: f64, y: f64, z: f64) -> Body {
        let mut b = Body::rocky(1.0).at(x, y).with_velocity(0.0, 0.0);
        b.z = z;
        b
    }

    #[test]
    fn morton_encode_zero_is_zero() {
        assert_eq!(morton_encode_3d(0, 0, 0), 0);
    }

    /// The lowest bit of each axis lands at the position the octant index
    /// would occupy: x → bit 0, y → bit 1, z → bit 2.
    #[test]
    fn morton_encode_low_bits_match_octant_pack() {
        assert_eq!(morton_encode_3d(1, 0, 0), 0b001);
        assert_eq!(morton_encode_3d(0, 1, 0), 0b010);
        assert_eq!(morton_encode_3d(0, 0, 1), 0b100);
        assert_eq!(morton_encode_3d(1, 1, 1), 0b111);
    }

    /// The k-th bit of each axis lands at code position 3k + axis_offset.
    #[test]
    fn morton_encode_high_bits_separated_by_strides() {
        assert_eq!(morton_encode_3d(0b10, 0, 0), 1u64 << 3);
        assert_eq!(morton_encode_3d(0, 0b10, 0), 1u64 << 4);
        assert_eq!(morton_encode_3d(0, 0, 0b10), 1u64 << 5);
    }

    #[test]
    fn morton_encode_max_21_bit_round_trip() {
        let q_max = (1u64 << MORTON_BITS) - 1;
        // Every bit set in every axis → every bit of the 63-bit code.
        let code = morton_encode_3d(q_max, q_max, q_max);
        assert_eq!(code, (1u64 << 63) - 1);
    }

    #[test]
    fn aabb_empty_slice_returns_unit_at_origin() {
        let aabb = compute_aabb(&[]);
        assert_eq!(aabb.center, [0.0; 3]);
        assert_eq!(aabb.half, 1.0);
    }

    #[test]
    fn aabb_is_cubic_padded_around_extreme_axis() {
        // y axis has the largest extent (4 units); cubic AABB picks that.
        let bodies = vec![body_at(-1.0, -2.0, -1.5), body_at(1.0, 2.0, 1.5)];
        let aabb = compute_aabb(&bodies);
        assert!((aabb.center[0] - 0.0).abs() < 1e-12);
        assert!((aabb.center[1] - 0.0).abs() < 1e-12);
        assert!((aabb.center[2] - 0.0).abs() < 1e-12);
        // half is at least extent/2 = 2.0, plus pad
        assert!(aabb.half > 2.0);
        assert!(aabb.half < 2.05);
    }

    #[test]
    fn permutation_is_a_bijection_on_indices() {
        let bodies: Vec<Body> = (0..32)
            .map(|i| {
                let t = i as f64 * 0.31;
                body_at(t.sin(), t.cos(), (t * 0.7).sin())
            })
            .collect();
        let perm = compute_morton_permutation(&bodies, &compute_aabb(&bodies));
        assert_eq!(perm.len(), bodies.len());
        let mut seen = vec![false; bodies.len()];
        for &i in &perm {
            assert!(!seen[i as usize], "body {i} appears twice");
            seen[i as usize] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn permutation_orders_octants_by_pack_index() {
        // One body in each of the 8 octants of a unit cube. Morton sort
        // must visit them in the same order the tree's child-index pack
        // would: octant 0 (---) first, octant 7 (+++) last.
        let bodies = vec![
            body_at(-0.5, -0.5, -0.5), // 0 (---)
            body_at(0.5, -0.5, -0.5),  // 1 (+--)
            body_at(-0.5, 0.5, -0.5),  // 2 (-+-)
            body_at(0.5, 0.5, -0.5),   // 3 (++-)
            body_at(-0.5, -0.5, 0.5),  // 4 (--+)
            body_at(0.5, -0.5, 0.5),   // 5 (+-+)
            body_at(-0.5, 0.5, 0.5),   // 6 (-++)
            body_at(0.5, 0.5, 0.5),    // 7 (+++)
        ];
        let perm = compute_morton_permutation(&bodies, &compute_aabb(&bodies));
        // Each input body's octant index equals its input index (set up
        // that way above), so the sorted permutation must be [0..8].
        assert_eq!(perm, (0u32..8).collect::<Vec<_>>());
    }

    #[test]
    fn permutation_clusters_bodies_within_same_octant() {
        // 4 bodies in the (---) octant + 4 in the (+++) octant should
        // come out as [4 from ---, then 4 from +++], with the within-
        // octant sub-order itself Morton-sorted.
        let mut bodies = Vec::new();
        for k in 0..4 {
            let d = -0.4 + 0.05 * k as f64;
            bodies.push(body_at(d, d, d)); // all in octant 0 (---)
        }
        for k in 0..4 {
            let d = 0.1 + 0.05 * k as f64;
            bodies.push(body_at(d, d, d)); // all in octant 7 (+++)
        }
        let perm = compute_morton_permutation(&bodies, &compute_aabb(&bodies));
        // First half of perm must be the indices [0, 1, 2, 3] (octant ---),
        // second half must be [4, 5, 6, 7] (octant +++).
        let first_half: std::collections::HashSet<u32> = perm[..4].iter().copied().collect();
        let second_half: std::collections::HashSet<u32> = perm[4..].iter().copied().collect();
        assert_eq!(first_half, [0u32, 1, 2, 3].into_iter().collect());
        assert_eq!(second_half, [4u32, 5, 6, 7].into_iter().collect());
    }

    #[test]
    fn permutation_handles_empty_slice() {
        let aabb = Aabb { center: [0.0; 3], half: 1.0 };
        let perm = compute_morton_permutation(&[], &aabb);
        assert!(perm.is_empty());
    }
}
