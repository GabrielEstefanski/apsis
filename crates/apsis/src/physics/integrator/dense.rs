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
//! let p = snap.interpolate(body_idx, h.clamp(0.0, 1.0));
//! ```

use crate::math::Vec3;
use crate::physics::integrator::IntegratorKind;

// ── DenseSnapshot ─────────────────────────────────────────────────────────────

/// Per-body IAS15 b-coefficients captured at the end of an accepted sub-step.
/// Laid out as `[Vec3; 7]` matching the seven Gauss-Radau nodes.
pub type DenseCoeffs = [Vec3; 7];

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
    pub x0: Vec<Vec3>,

    /// Velocities at `t0`, one per body.
    pub v0: Vec<Vec3>,

    /// Accelerations at `t0`, one per body.
    pub a0: Vec<Vec3>,

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
    pub fn interpolate(&self, i: usize, h: f64) -> Vec3 {
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

    /// Whether `x0`, `v0`, `a0`, and (when populated) `b` all carry the
    /// same body count.
    ///
    /// `interpolate(i, h)` indexes every internal vector at the same `i`,
    /// so a snapshot whose internal arrays disagree on length will panic
    /// when the consumer's loop runs past the shortest. Producers that
    /// build a snapshot from heterogeneous sources (for example, the
    /// Order-2 fallback in `System::step`, which captures `x0` / `v0`
    /// from `bodies` but `a0` from `scratch_acc`) must verify shape
    /// consistency at construction time, and consumers that hold a
    /// `DenseSnapshot` across mutation of `bodies` should re-check
    /// before each render.
    #[inline]
    pub fn is_shape_consistent(&self) -> bool {
        let n = self.x0.len();
        self.v0.len() == n && self.a0.len() == n && (self.b.is_empty() || self.b.len() == n)
    }
}

// ── Interpolation kernels ─────────────────────────────────────────────────────

/// IAS15 degree-15 polynomial interpolation (Rein & Spiegel 2015, eq. 9).
///
/// Evaluates position at substep fraction `h ∈ [0, 1]` given the start-of-step
/// kinematics and the seven Gauss-Radau b-coefficients.
///
/// `x(h) = x₀ + v₀·h·dt + (h·dt)² · [a₀/2 + b₀·h/6 + b₁·h²/12 + ··· + b₆·h⁷/72]`
///
/// Component-by-component scalar form: `(b·h^k)/c + a·0.5` is computed
/// per axis. Re-associating into `Vec3` ops would shift ULPs and is
/// therefore avoided — the IAS15 module sits at the f64 noise floor
/// where reduction order is observable downstream
/// (cf. `docs/experiments/2026-04-29-3d-port-baseline.md`).
#[inline]
pub fn predict_ias15(x0: Vec3, v0: Vec3, a0: Vec3, b: &DenseCoeffs, h: f64, dt: f64) -> Vec3 {
    let h2 = h * h;
    let h3 = h2 * h;
    let h4 = h3 * h;
    let h5 = h4 * h;
    let h6 = h5 * h;
    let h7 = h6 * h;

    let dt2 = dt * dt;

    let ax = a0.x * 0.5
        + b[0].x * h / 6.0
        + b[1].x * h2 / 12.0
        + b[2].x * h3 / 20.0
        + b[3].x * h4 / 30.0
        + b[4].x * h5 / 42.0
        + b[5].x * h6 / 56.0
        + b[6].x * h7 / 72.0;

    let ay = a0.y * 0.5
        + b[0].y * h / 6.0
        + b[1].y * h2 / 12.0
        + b[2].y * h3 / 20.0
        + b[3].y * h4 / 30.0
        + b[4].y * h5 / 42.0
        + b[5].y * h6 / 56.0
        + b[6].y * h7 / 72.0;

    let az = a0.z * 0.5
        + b[0].z * h / 6.0
        + b[1].z * h2 / 12.0
        + b[2].z * h3 / 20.0
        + b[3].z * h4 / 30.0
        + b[4].z * h5 / 42.0
        + b[5].z * h6 / 56.0
        + b[6].z * h7 / 72.0;

    Vec3::new(
        x0.x + v0.x * h * dt + h2 * dt2 * ax,
        x0.y + v0.y * h * dt + h2 * dt2 * ay,
        x0.z + v0.z * h * dt + h2 * dt2 * az,
    )
}

/// 2nd-order Taylor interpolation: `x₀ + v₀·h·dt + ½·a₀·(h·dt)²`.
///
/// Used for VV, Yoshida-4, and Wisdom–Holman.  Accurate to O(dt²) which is
/// sufficient for visual smoothness at typical interactive step sizes.
#[inline]
pub fn predict_order2(x0: Vec3, v0: Vec3, a0: Vec3, h: f64, dt: f64) -> Vec3 {
    let s = h * dt;
    Vec3::new(
        x0.x + v0.x * s + 0.5 * a0.x * s * s,
        x0.y + v0.y * s + 0.5 * a0.y * s * s,
        x0.z + v0.z * s + 0.5 * a0.z * s * s,
    )
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
/// reads `body.(vx, vy, vz)` directly, so leaving the body velocities at
/// their start-of-step values biases every node evaluation by `O(a · dt)`.
/// On a Mercury 1PN integration the bias accumulates linearly to
/// ~10⁻³ relative precession error over 500 orbits — see
/// `docs/experiments/2026-04-28-ias15-velocity-prediction-bug.md`.
#[inline]
pub fn predict_v_ias15(v0: Vec3, a0: Vec3, b: &DenseCoeffs, h: f64, dt: f64) -> Vec3 {
    let h2 = h * h;
    let h3 = h2 * h;
    let h4 = h3 * h;
    let h5 = h4 * h;
    let h6 = h5 * h;
    let h7 = h6 * h;

    let inner_x = a0.x
        + b[0].x * h / 2.0
        + b[1].x * h2 / 3.0
        + b[2].x * h3 / 4.0
        + b[3].x * h4 / 5.0
        + b[4].x * h5 / 6.0
        + b[5].x * h6 / 7.0
        + b[6].x * h7 / 8.0;

    let inner_y = a0.y
        + b[0].y * h / 2.0
        + b[1].y * h2 / 3.0
        + b[2].y * h3 / 4.0
        + b[3].y * h4 / 5.0
        + b[4].y * h5 / 6.0
        + b[5].y * h6 / 7.0
        + b[6].y * h7 / 8.0;

    let inner_z = a0.z
        + b[0].z * h / 2.0
        + b[1].z * h2 / 3.0
        + b[2].z * h3 / 4.0
        + b[3].z * h4 / 5.0
        + b[4].z * h5 / 6.0
        + b[5].z * h6 / 7.0
        + b[6].z * h7 / 8.0;

    Vec3::new(v0.x + h * dt * inner_x, v0.y + h * dt * inner_y, v0.z + h * dt * inner_z)
}

#[cfg(test)]
mod tests {
    use super::{DenseCoeffs, predict_ias15, predict_v_ias15};
    use crate::math::Vec3;

    fn sample_b() -> DenseCoeffs {
        [
            Vec3::new(0.11, 0.21, 0.31),
            Vec3::new(0.12, 0.22, 0.32),
            Vec3::new(0.13, 0.23, 0.33),
            Vec3::new(0.14, 0.24, 0.34),
            Vec3::new(0.15, 0.25, 0.35),
            Vec3::new(0.16, 0.26, 0.36),
            Vec3::new(0.17, 0.27, 0.37),
        ]
    }

    #[test]
    fn predict_v_ias15_at_h_zero_returns_v0() {
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let b = sample_b();
        assert_eq!(predict_v_ias15(v0, a0, &b, 0.0, 1e-3), v0);
    }

    #[test]
    fn predict_v_ias15_recovers_constant_acceleration() {
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let b: DenseCoeffs = [Vec3::ZERO; 7];
        let dt = 1e-3;
        for h in [0.1, 0.3, 0.5, 0.7, 1.0] {
            let v = predict_v_ias15(v0, a0, &b, h, dt);
            let expected =
                Vec3::new(v0.x + a0.x * h * dt, v0.y + a0.y * h * dt, v0.z + a0.z * h * dt);
            assert!(
                (v.x - expected.x).abs() < 1e-15,
                "vx at h={h}: got {} expected {}",
                v.x,
                expected.x
            );
            assert!(
                (v.y - expected.y).abs() < 1e-15,
                "vy at h={h}: got {} expected {}",
                v.y,
                expected.y
            );
            assert!(
                (v.z - expected.z).abs() < 1e-15,
                "vz at h={h}: got {} expected {}",
                v.z,
                expected.z
            );
        }
    }

    #[test]
    fn predict_v_ias15_is_derivative_of_predict_ias15() {
        // Tolerance reflects the central-difference round-off floor
        // (O(eps²) + O(ε_mach / eps)).
        let x0 = Vec3::new(0.5, 0.3, -0.2);
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let b = sample_b();
        let dt = 1e-3;
        let eps = 1e-5;
        for h in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let xp = predict_ias15(x0, v0, a0, &b, h + eps, dt);
            let xm = predict_ias15(x0, v0, a0, &b, h - eps, dt);
            // Central difference in `h` then convert to derivative in
            // physical time: `dx/dt = (1/dt) · dx/dh`.
            let v_num = Vec3::new(
                (xp.x - xm.x) / (2.0 * eps * dt),
                (xp.y - xm.y) / (2.0 * eps * dt),
                (xp.z - xm.z) / (2.0 * eps * dt),
            );
            let v = predict_v_ias15(v0, a0, &b, h, dt);
            assert!(
                (v.x - v_num.x).abs() < 1e-7,
                "vx at h={h}: analytical {} numerical {}",
                v.x,
                v_num.x
            );
            assert!(
                (v.y - v_num.y).abs() < 1e-7,
                "vy at h={h}: analytical {} numerical {}",
                v.y,
                v_num.y
            );
            assert!(
                (v.z - v_num.z).abs() < 1e-7,
                "vz at h={h}: analytical {} numerical {}",
                v.z,
                v_num.z
            );
        }
    }
}
