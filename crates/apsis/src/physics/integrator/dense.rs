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

/// IAS15 degree-15 velocity at substep fraction `h ∈ [0, 1]` (Rein & Spiegel
/// 2015, eq. 11).
///
/// Differentiating the position polynomial in [`predict_ias15`] (eq. 9) once
/// with respect to physical time `t = h · dt` gives:
///
/// `v(h) = v₀ + (h·dt) · [a₀ + b₀·h/2 + b₁·h²/3 + b₂·h³/4 + b₃·h⁴/5 + b₄·h⁵/6 + b₅·h⁶/7 + b₆·h⁷/8]`
///
/// Required at every Gauss–Radau substep node when forces are evaluated
/// inside Picard predictor–corrector iteration: any velocity-dependent
/// perturbation registered through
/// [`PerturbationForce::accumulate`](crate::physics::integrator::PerturbationForce::accumulate)
/// reads `body.(vx, vy)` directly, so leaving the body velocities at
/// their start-of-step values biases every node evaluation by `O(a · dt)`.
/// On a Mercury 1PN integration the bias accumulates linearly to
/// ~10⁻³ relative precession error over 500 orbits — see
/// `docs/experiments/2026-04-28-ias15-velocity-prediction-bug.md`.
#[inline]
pub fn predict_v_ias15(
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

    let inner_x = a0.0
        + b[0].0 * h / 2.0
        + b[1].0 * h2 / 3.0
        + b[2].0 * h3 / 4.0
        + b[3].0 * h4 / 5.0
        + b[4].0 * h5 / 6.0
        + b[5].0 * h6 / 7.0
        + b[6].0 * h7 / 8.0;

    let inner_y = a0.1
        + b[0].1 * h / 2.0
        + b[1].1 * h2 / 3.0
        + b[2].1 * h3 / 4.0
        + b[3].1 * h4 / 5.0
        + b[4].1 * h5 / 6.0
        + b[5].1 * h6 / 7.0
        + b[6].1 * h7 / 8.0;

    (v0.0 + h * dt * inner_x, v0.1 + h * dt * inner_y)
}

#[cfg(test)]
mod tests {
    use super::{DenseCoeffs, predict_ias15, predict_v_ias15};

    fn sample_b() -> DenseCoeffs {
        [
            (0.11, 0.21),
            (0.12, 0.22),
            (0.13, 0.23),
            (0.14, 0.24),
            (0.15, 0.25),
            (0.16, 0.26),
            (0.17, 0.27),
        ]
    }

    #[test]
    fn predict_v_ias15_at_h_zero_returns_v0() {
        let v0 = (1.5, -0.7);
        let a0 = (0.3, 0.2);
        let b = sample_b();
        assert_eq!(predict_v_ias15(v0, a0, &b, 0.0, 1e-3), v0);
    }

    #[test]
    fn predict_v_ias15_recovers_constant_acceleration() {
        let v0 = (1.5, -0.7);
        let a0 = (0.3, 0.2);
        let b: DenseCoeffs = [(0.0, 0.0); 7];
        let dt = 1e-3;
        for h in [0.1, 0.3, 0.5, 0.7, 1.0] {
            let (vx, vy) = predict_v_ias15(v0, a0, &b, h, dt);
            let expected_vx = v0.0 + a0.0 * h * dt;
            let expected_vy = v0.1 + a0.1 * h * dt;
            assert!(
                (vx - expected_vx).abs() < 1e-15,
                "vx at h={h}: got {vx}, expected {expected_vx}"
            );
            assert!(
                (vy - expected_vy).abs() < 1e-15,
                "vy at h={h}: got {vy}, expected {expected_vy}"
            );
        }
    }

    #[test]
    fn predict_v_ias15_is_derivative_of_predict_ias15() {
        // Tolerance reflects the central-difference round-off floor
        // (O(eps²) + O(ε_mach / eps)).
        let x0 = (0.5, 0.3);
        let v0 = (1.5, -0.7);
        let a0 = (0.3, 0.2);
        let b = sample_b();
        let dt = 1e-3;
        let eps = 1e-5;
        for h in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let (xp, yp) = predict_ias15(x0, v0, a0, &b, h + eps, dt);
            let (xm, ym) = predict_ias15(x0, v0, a0, &b, h - eps, dt);
            // Central difference in `h` then convert to derivative in
            // physical time: `dx/dt = (1/dt) · dx/dh`.
            let v_num_x = (xp - xm) / (2.0 * eps * dt);
            let v_num_y = (yp - ym) / (2.0 * eps * dt);
            let (vx, vy) = predict_v_ias15(v0, a0, &b, h, dt);
            assert!(
                (vx - v_num_x).abs() < 1e-7,
                "vx at h={h}: analytical {vx}, numerical {v_num_x}"
            );
            assert!(
                (vy - v_num_y).abs() < 1e-7,
                "vy at h={h}: analytical {vy}, numerical {v_num_y}"
            );
        }
    }
}
