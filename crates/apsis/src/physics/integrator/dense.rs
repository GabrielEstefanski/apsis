//! Dense output — sub-step position interpolation for smooth rendering.
//!
//! Each completed integration step (or sub-step for IAS15) records a
//! [`DenseSnapshot`] that captures the state needed to evaluate body
//! positions at any time `t ∈ [t₀, t₀ + dt]` without re-running physics.
//!
//! # Interpolation formulas
//!
//! | Integrator | Formula | Order |
//! |------------|---------|-------|
//! | IAS15      | Rein & Spiegel (2015) polynomial via b-coefficients | ≤ 15 |
//! | VV / Y4 / WH | 2nd-order Taylor: `x₀ + v₀·h·dt + ½·a₀·(h·dt)²` | 2 |
//!
//! The IAS15 polynomial is exact to the precision of the accepted b-coefficients.
//! The Order-2 fallback is sufficient for smooth visual rendering between steps —
//! the physical error is O(dt²), below the pixel resolution for typical step sizes.
//!
//! # Usage
//!
//! ```ignore
//! let h = (t_render - snap.t0) / snap.dt;   // ∈ [0, 1]
//! let (x, y) = snap.interpolate(body_idx, h.clamp(0.0, 1.0));
//! ```

use crate::physics::integrator::IntegratorKind;

// ── DenseSnapshot ─────────────────────────────────────────────────────────────

/// Per-body IAS15 b-coefficients captured at the end of an accepted sub-step.
/// Laid out as `[(bx, by); 7]` matching the seven Gauss-Radau nodes.
pub type DenseCoeffs = [(f64, f64); 7];

/// Snapshot of the state needed to interpolate body positions within one step.
///
/// The caller (physics thread) sets [`t0`] to `system.t() - dt` after the step
/// completes; the integrator only needs to fill the shape-of-trajectory fields
/// (`x0`, `v0`, `a0`, `b`).
#[derive(Clone)]
pub struct DenseSnapshot {
    /// Absolute sim time at the start of this step.
    pub t0: f64,

    /// Duration of this step (sub-step dt for IAS15, full system dt for others).
    pub dt: f64,

    /// World positions at `t0`, one per body.
    pub x0: Vec<(f64, f64)>,

    /// Velocities at `t0`, one per body.
    pub v0: Vec<(f64, f64)>,

    /// Accelerations at `t0`, one per body.
    pub a0: Vec<(f64, f64)>,

    /// IAS15 b-coefficients, one [`DenseCoeffs`] per body.
    /// **Empty for non-IAS15 integrators** — the [`interpolate`](Self::interpolate)
    /// method falls back to the 2nd-order Taylor formula automatically.
    pub b: Vec<DenseCoeffs>,

    /// Identifies the integrator that produced this snapshot.
    pub kind: IntegratorKind,
}

impl DenseSnapshot {
    /// Interpolated world position for body `i` at normalised time `h ∈ [0, 1]`.
    ///
    /// Panics in debug mode if `i >= self.x0.len()`.
    #[inline]
    pub fn interpolate(&self, i: usize, h: f64) -> (f64, f64) {
        debug_assert!(i < self.x0.len(), "body index out of range");

        let x0 = self.x0[i];
        let v0 = self.v0[i];
        let a0 = self.a0[i];
        let dt = self.dt;

        if !self.b.is_empty() {
            predict_ias15(x0, v0, a0, &self.b[i], h, dt)
        } else {
            predict_order2(x0, v0, a0, h, dt)
        }
    }

    /// Number of bodies in this snapshot.
    #[inline]
    pub fn n_bodies(&self) -> usize {
        self.x0.len()
    }
}

// ── Interpolation kernels ─────────────────────────────────────────────────────

/// IAS15 degree-15 polynomial interpolation (Rein & Spiegel 2015, eq. 9).
///
/// Evaluates position at substep fraction `h ∈ [0, 1]` given the start-of-step
/// kinematics and the seven Gauss-Radau b-coefficients.
///
/// `x(h) = x₀ + v₀·h·dt + (h·dt)² · [a₀/2 + b₀·h/6 + b₁·h²/12 + ··· + b₆·h⁷/72]`
#[inline]
pub fn predict_ias15(
    x0: (f64, f64),
    v0: (f64, f64),
    a0: (f64, f64),
    b: &DenseCoeffs,
    h: f64,
    dt: f64,
) -> (f64, f64) {
    let h2 = h * h;
    let h3 = h2 * h;
    let h4 = h3 * h;
    let h5 = h4 * h;
    let h6 = h5 * h;
    let h7 = h6 * h;

    let dt2 = dt * dt;

    let ax = a0.0 * 0.5
        + b[0].0 * h / 6.0
        + b[1].0 * h2 / 12.0
        + b[2].0 * h3 / 20.0
        + b[3].0 * h4 / 30.0
        + b[4].0 * h5 / 42.0
        + b[5].0 * h6 / 56.0
        + b[6].0 * h7 / 72.0;

    let ay = a0.1 * 0.5
        + b[0].1 * h / 6.0
        + b[1].1 * h2 / 12.0
        + b[2].1 * h3 / 20.0
        + b[3].1 * h4 / 30.0
        + b[4].1 * h5 / 42.0
        + b[5].1 * h6 / 56.0
        + b[6].1 * h7 / 72.0;

    (x0.0 + v0.0 * h * dt + h2 * dt2 * ax, x0.1 + v0.1 * h * dt + h2 * dt2 * ay)
}

/// 2nd-order Taylor interpolation: `x₀ + v₀·h·dt + ½·a₀·(h·dt)²`.
///
/// Used for VV, Yoshida-4, and Wisdom–Holman.  Accurate to O(dt²) which is
/// sufficient for visual smoothness at typical interactive step sizes.
#[inline]
pub fn predict_order2(
    x0: (f64, f64),
    v0: (f64, f64),
    a0: (f64, f64),
    h: f64,
    dt: f64,
) -> (f64, f64) {
    let s = h * dt;
    (x0.0 + v0.0 * s + 0.5 * a0.0 * s * s, x0.1 + v0.1 * s + 0.5 * a0.1 * s * s)
}
