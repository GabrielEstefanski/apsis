//! Per-step Struct-of-Arrays snapshot of the body state read by the
//! gravity hot path.
//!
//! [`BodyArrays`] holds the four fields the BH walk and the tree build
//! actually read — `pos_x`, `pos_y`, `pos_z`, `mass` — laid out as four
//! contiguous `Vec<f64>`. The integrator and the public API continue to
//! read [`Body`]; SoA is execution state for one phase, not a domain
//! type. Force-law softening (when the active kernel uses it) is read
//! from the kernel itself, not the bodies.
//!
//! ## Lifecycle (per step)
//!
//! ```text
//! pack_from(&[Body])    once at step start
//! tree build reads SoA
//! BH walk reads SoA
//! force eval returns
//! SoA stale until the next pack_from overwrites
//! ```
//!
//! Never mutated during the step. No invalidation flags, no dirty bits —
//! [`pack_from`] is the synchronisation point. Reusable buffer: callers
//! allocate one [`BodyArrays`] at startup and call [`pack_from`] each step.
//!
//! ## Field selection
//!
//! Four fields, exactly one BH-walk leaf-pair payload per body (32 bytes).
//! Velocity is deliberately excluded — including `vel_x/y/z` would inflate
//! the row size and halve the per-cache-line density that motivates the
//! refactor. The integrator reads/writes velocity through [`Body`] and is
//! compute-bound in its own internal coefficient arrays, not Body-bound.
//!
//! Reference: `docs/experiments/2026-05-10-soa-layout.md` §Design constraint.

use crate::domain::body::Body;

/// Four-field Struct-of-Arrays snapshot of the body state for the gravity
/// hot path. See module docs for the lifecycle contract.
///
/// All fields are public so the kernel and tree-build inner loops can read
/// `arrays.pos_x[i]` directly without an accessor — keeps the loop in a
/// shape SIMD can unroll later without structural refactor.
#[derive(Debug, Clone, Default)]
pub struct BodyArrays {
    pub pos_x: Vec<f64>,
    pub pos_y: Vec<f64>,
    pub pos_z: Vec<f64>,
    pub mass: Vec<f64>,
}

impl BodyArrays {
    /// Empty arrays. Call [`pack_from`] before reading.
    pub fn new() -> Self {
        Self::default()
    }

    /// Empty arrays with `capacity` reserved per field. Useful at startup
    /// when the maximum body count is known.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            pos_x: Vec::with_capacity(capacity),
            pos_y: Vec::with_capacity(capacity),
            pos_z: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
        }
    }

    /// Body count after the last [`pack_from`]. All field arrays are
    /// invariantly the same length.
    #[inline]
    pub fn len(&self) -> usize {
        self.pos_x.len()
    }

    /// `true` iff the snapshot carries zero bodies.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.pos_x.is_empty()
    }

    /// Discard the current snapshot. Useful for callers that retain the
    /// buffer across steps and want to assert clean state in tests; routine
    /// use prefers [`pack_from`], which clears as part of its contract.
    pub fn clear(&mut self) {
        self.pos_x.clear();
        self.pos_y.clear();
        self.pos_z.clear();
        self.mass.clear();
    }

    /// Overwrite the snapshot with the four hot fields of `bodies`.
    ///
    /// Existing buffer capacity is reused; only allocations occur if the
    /// new body count exceeds previous capacity. Field arrays are always
    /// pushed in lockstep, so [`len`](Self::len) is well-defined post-call.
    ///
    /// **Contract.** This is the sole writer of the snapshot. The arrays
    /// must not be mutated between calls — the SoA snapshot is read-only
    /// during a step (the experiment's design constraint).
    pub fn pack_from(&mut self, bodies: &[Body]) {
        self.pos_x.clear();
        self.pos_y.clear();
        self.pos_z.clear();
        self.mass.clear();

        let n = bodies.len();
        self.pos_x.reserve(n);
        self.pos_y.reserve(n);
        self.pos_z.reserve(n);
        self.mass.reserve(n);

        for b in bodies {
            self.pos_x.push(b.pos_x);
            self.pos_y.push(b.pos_y);
            self.pos_z.push(b.pos_z);
            self.mass.push(b.mass);
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_bodies() -> Vec<Body> {
        vec![
            Body::rocky(1.0).at_3d(1.0, 2.0, 3.0).with_velocity_3d(0.1, 0.2, 0.3),
            Body::star(2.0).at_3d(-4.0, 5.0, -6.0).with_velocity_3d(-0.4, -0.5, -0.6),
            Body::asteroid(0.001).at_3d(7.0, -8.0, 9.0).with_velocity_3d(0.7, -0.8, 0.9),
        ]
    }

    /// `pack_from` preserves every hot field bit-for-bit. If a field is
    /// reordered, coerced, or arithmetically transformed during packing,
    /// the accelerations the BH walk computes downstream will diverge
    /// from those computed against the equivalent `&[Body]`.
    #[test]
    fn pack_from_preserves_fields_bit_exact() {
        let bodies = sample_bodies();
        let mut arrays = BodyArrays::new();
        arrays.pack_from(&bodies);

        assert_eq!(arrays.len(), bodies.len());
        for (i, b) in bodies.iter().enumerate() {
            assert_eq!(arrays.pos_x[i].to_bits(), b.pos_x.to_bits());
            assert_eq!(arrays.pos_y[i].to_bits(), b.pos_y.to_bits());
            assert_eq!(arrays.pos_z[i].to_bits(), b.pos_z.to_bits());
            assert_eq!(arrays.mass[i].to_bits(), b.mass.to_bits());
        }
    }

    /// All four field arrays must always be the same length. A divergence
    /// would mean a partial pack happened (panic mid-loop, etc.) and any
    /// indexed read of `arrays.pos_x[i]` paired with `arrays.mass[i]`
    /// silently reads stale data.
    #[test]
    fn pack_from_keeps_field_arrays_in_lockstep() {
        let bodies = sample_bodies();
        let mut arrays = BodyArrays::new();
        arrays.pack_from(&bodies);

        let n = arrays.len();
        assert_eq!(arrays.pos_x.len(), n);
        assert_eq!(arrays.pos_y.len(), n);
        assert_eq!(arrays.pos_z.len(), n);
        assert_eq!(arrays.mass.len(), n);
    }

    /// Empty input produces empty arrays (not panic, not stale carryover).
    #[test]
    fn pack_from_empty_slice_yields_empty_arrays() {
        let mut arrays = BodyArrays::with_capacity(8);
        arrays.pack_from(&[]);
        assert!(arrays.is_empty());
        assert_eq!(arrays.len(), 0);
    }

    /// Reusing the buffer across calls overwrites cleanly — the
    /// canonical use pattern (one buffer allocated at startup, packed
    /// every step).
    #[test]
    fn pack_from_reuses_buffer_without_carryover() {
        let mut arrays = BodyArrays::with_capacity(8);

        let first = sample_bodies();
        arrays.pack_from(&first);
        assert_eq!(arrays.len(), 3);

        let second = vec![Body::rocky(5.0).at_3d(10.0, 20.0, 30.0)];
        arrays.pack_from(&second);
        assert_eq!(arrays.len(), 1);
        assert_eq!(arrays.pos_x[0].to_bits(), second[0].pos_x.to_bits());
        assert_eq!(arrays.pos_y[0].to_bits(), second[0].pos_y.to_bits());
        assert_eq!(arrays.pos_z[0].to_bits(), second[0].pos_z.to_bits());
    }

    /// `clear` brings the buffer back to empty without dropping capacity —
    /// matches `Vec::clear` semantics. Used by tests to assert clean
    /// state; production callers prefer `pack_from` which clears as part
    /// of its contract.
    #[test]
    fn clear_empties_without_dropping_capacity() {
        let mut arrays = BodyArrays::with_capacity(16);
        arrays.pack_from(&sample_bodies());
        let cap_before = arrays.pos_x.capacity();
        arrays.clear();
        assert!(arrays.is_empty());
        assert_eq!(arrays.pos_x.capacity(), cap_before);
    }

    /// Constructed via `with_capacity`, the arrays start empty but with
    /// the requested capacity reserved.
    #[test]
    fn with_capacity_reserves_without_pushing() {
        let arrays = BodyArrays::with_capacity(64);
        assert!(arrays.is_empty());
        assert!(arrays.pos_x.capacity() >= 64);
        assert!(arrays.pos_y.capacity() >= 64);
        assert!(arrays.pos_z.capacity() >= 64);
        assert!(arrays.mass.capacity() >= 64);
    }
}
