//! Keplerian two-body propagator used by the Wisdom–Holman integrator.
//!
//! Operates in the **relative frame** (body relative to a fixed central mass).
//! Units are internal simulation units where `G = 1`; `mu` is the gravitational
//! parameter in those units.
//!
//! # Formulation
//!
//! The propagator uses the **universal variable** `χ` and **Stumpff functions**
//! `c(ψ)`, `s(ψ)` throughout, following Battin (1987) §4.5. The formulation is
//! continuously valid for elliptic (`e < 1`), parabolic (`e = 1`), and
//! hyperbolic (`e > 1`) trajectories without branching on the orbit type. The
//! universal-variable solver depends only on scalar invariants of the state
//! (`|r|`, `|v|²`, `r · v`), so the same algorithm propagates state in any
//! number of spatial dimensions; this implementation uses [`Vec3`].
//!
//! # References
//!
//! - Battin, R. H. (1987). *An Introduction to the Mathematics and Methods
//!   of Astrodynamics*. AIAA Education Series. §4.5.
//! - Bate, R. R., Mueller, D. D., & White, J. E. (1971).
//!   *Fundamentals of Astrodynamics*. Dover. §2.8.
//! - Wisdom, J. & Holman, M. (1991). *Astron. J.* 102, 1528–1538.

use crate::math::Vec3;

// ── Stumpff functions ─────────────────────────────────────────────────────────

/// Stumpff function `c(ψ) = (1 − cos √ψ) / ψ`, continuously defined for all
/// `ψ ∈ ℝ` via the closed-form identities for `ψ ≷ 0` and a Taylor expansion
/// near zero.
fn stumpff_c(psi: f64) -> f64 {
    if psi.abs() < 1e-6 {
        return 0.5 - psi / 24.0 + psi * psi / 720.0;
    }
    if psi > 0.0 { (1.0 - psi.sqrt().cos()) / psi } else { ((-psi).sqrt().cosh() - 1.0) / (-psi) }
}

/// Stumpff function `s(ψ) = (√ψ − sin √ψ) / ψ^(3/2)`, continuously defined
/// for all `ψ ∈ ℝ` via the closed-form identities for `ψ ≷ 0` and a Taylor
/// expansion near zero.
fn stumpff_s(psi: f64) -> f64 {
    if psi.abs() < 1e-6 {
        return 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
    }
    if psi > 0.0 {
        let sq = psi.sqrt();
        (sq - sq.sin()) / (psi * sq)
    } else {
        let sq = (-psi).sqrt();
        (sq.sinh() - sq) / ((-psi) * sq)
    }
}

// ── Universal-variable time equation ─────────────────────────────────────────

/// Evaluates the universal Kepler time equation and its derivative.
///
/// Returns `(f(χ) − √μ · dt,  df/dχ)` where
///
/// ```text
/// f(χ) = (r₀ · v_r₀ / √μ) · χ² c(ψ)  +  (1 − r₀ α) · χ³ s(ψ)  +  r₀ χ
/// ```
///
/// and `ψ = α χ²`, `α = 1/a` (negative for hyperbolic orbits).
/// The derivative `df/dχ = r(χ) / √μ` is the instantaneous radius, which is
/// always positive and guarantees monotone convergence of Newton's method.
#[inline]
fn kepler_universal_fd(chi: f64, alpha: f64, r0: f64, vr0: f64, mu: f64, dt: f64) -> (f64, f64) {
    let sqrt_mu = mu.sqrt();
    let psi = alpha * chi * chi;
    let c = stumpff_c(psi);
    let s = stumpff_s(psi);
    let chi2 = chi * chi;
    let chi3 = chi2 * chi;

    let r = (r0 * vr0 / sqrt_mu) * chi * (1.0 - psi * s) + (1.0 - r0 * alpha) * chi2 * c + r0;

    let f_val =
        (r0 * vr0 / sqrt_mu) * chi2 * c + (1.0 - r0 * alpha) * chi3 * s + r0 * chi / sqrt_mu
            - sqrt_mu * dt;

    (f_val, r / sqrt_mu)
}

// ── Public propagator ─────────────────────────────────────────────────────────

/// Advances a Keplerian two-body orbit by time `dt` using the universal
/// variable formulation.
///
/// Given the relative state vector `(r0, v0)` of a body with respect to a
/// fixed central mass and the gravitational parameter `mu`, returns the
/// propagated state `(r1, v1)` via the Lagrange coefficient equations
///
/// ```text
/// r⃗' = f · r⃗₀  +  g · v⃗₀
/// v⃗' = ḟ · r⃗₀  +  ġ · v⃗₀
/// ```
///
/// where `f`, `g`, `ḟ`, `ġ` are scalar coefficients expressed in terms of `χ`
/// and Stumpff functions, valid for all orbit types without branching. Because
/// the coefficients are scalar and the Lagrange relations are linear in the
/// input vectors, the propagator is dimensionally agnostic — the same code
/// path handles planar and inclined orbits identically.
///
/// # Convergence
///
/// Newton–Raphson on the universal time equation converges quadratically. The
/// derivative `df/dχ = r/√μ > 0` is always positive, guaranteeing monotone
/// convergence. The tolerance is `|Δχ| < 10⁻¹²`, with a safety cap at 50
/// iterations.
///
/// # Degenerate inputs
///
/// A zero or near-zero radius `|r₀| → 0` indicates a collision singularity;
/// the function returns the input state unchanged in that case.
///
/// # References
///
/// - Battin (1987) §4.5; Bate–Mueller–White (1971) §2.8.
/// - Wisdom & Holman (1991) *Astron. J.* 102, 1528–1538.
pub fn kepler_step(r0: Vec3, v0: Vec3, dt: f64, mu: f64) -> (Vec3, Vec3) {
    let r0_norm = r0.length();
    if r0_norm < 1e-30 {
        return (r0, v0);
    }

    let v_sq = v0.length_squared();
    let vr0 = r0.dot(v0) / r0_norm;
    let alpha = 2.0 / r0_norm - v_sq / mu; // = 1/a; negative iff hyperbolic
    let sqrt_mu = mu.sqrt();

    // Initial guess for χ (Bate–Mueller–White §2.8, adapted for sign of dt)
    let chi0 = if alpha > 0.0 {
        // Elliptic: χ ≈ n · dt · a,  n = √(μ / a³) = √μ · α^(3/2)
        sqrt_mu * dt * alpha
    } else {
        // Hyperbolic: use the asymptotic chord length as seed
        let a = 1.0 / alpha; // a < 0
        let s = dt.signum() * (2.0 * mu * alpha.abs()).sqrt() * dt.abs()
            / (r0_norm * vr0 + dt.signum() * (mu * (-a)).sqrt() * (1.0 - r0_norm * alpha));
        (2.0 * mu * alpha.abs()).sqrt() * s.tanh() / alpha.abs().sqrt()
    };

    // Newton–Raphson on the universal Kepler equation
    let mut chi = chi0;
    for _ in 0..50 {
        let (f_val, df) = kepler_universal_fd(chi, alpha, r0_norm, vr0, mu, dt);
        let dchi = f_val / df;
        chi -= dchi;

        if !chi.is_finite() {
            break;
        }

        if dchi.abs() < 1e-12 {
            break;
        }
    }

    // Final Stumpff evaluation for Lagrange coefficients
    let psi = alpha * chi * chi;
    let c = stumpff_c(psi);
    let s = stumpff_s(psi);
    let chi2 = chi * chi;
    let chi3 = chi2 * chi;

    let r1_norm = (r0_norm * vr0 / sqrt_mu) * chi * (1.0 - psi * s)
        + (1.0 - r0_norm * alpha) * chi2 * c
        + r0_norm;

    let f_lag = 1.0 - chi2 * c / r0_norm;
    let g_lag = dt - chi3 * s / sqrt_mu;
    let df_lag = -sqrt_mu / (r1_norm * r0_norm) * chi * (1.0 - psi * s);
    let dg_lag = 1.0 - chi2 * c / r1_norm;

    let r1 = f_lag * r0 + g_lag * v0;
    let v1 = df_lag * r0 + dg_lag * v0;

    (r1, v1)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Circular orbit at `r = 1`, `v = 1`, `mu = 1` should return to its
    /// starting point after one full period `T = 2π`.
    #[test]
    fn circular_orbit_returns_to_start_after_one_period() {
        let r0 = Vec3::new(1.0, 0.0, 0.0);
        let v0 = Vec3::new(0.0, 1.0, 0.0);
        let mu = 1.0;
        let period = 2.0 * std::f64::consts::PI;

        let (r1, v1) = kepler_step(r0, v0, period, mu);

        assert!((r1 - r0).length() < 1e-10, "position drift {}", (r1 - r0).length());
        assert!((v1 - v0).length() < 1e-10, "velocity drift {}", (v1 - v0).length());
    }

    /// Energy is exactly conserved by the analytical propagator (within f64
    /// round-off accumulated through the Newton iteration).
    #[test]
    fn elliptic_orbit_conserves_energy() {
        let r0 = Vec3::new(1.0, 0.0, 0.0);
        let v0 = Vec3::new(0.0, 1.5, 0.0); // bound: v² < 2μ/r
        let mu = 1.0;
        let energy0 = 0.5 * v0.length_squared() - mu / r0.length();

        let mut r = r0;
        let mut v = v0;
        for _ in 0..1000 {
            let (rn, vn) = kepler_step(r, v, 0.05, mu);
            r = rn;
            v = vn;
        }

        let energy1 = 0.5 * v.length_squared() - mu / r.length();
        let drift = ((energy1 - energy0) / energy0).abs();
        assert!(drift < 1e-12, "energy drifted by {drift}");
    }

    /// Specific angular momentum is exactly conserved by the analytical
    /// propagator.
    #[test]
    fn elliptic_orbit_conserves_angular_momentum_3d() {
        // Inclined orbit (out-of-plane component); the propagator should treat
        // it identically to the planar case.
        let r0 = Vec3::new(1.0, 0.0, 0.5);
        let v0 = Vec3::new(0.0, 1.2, 0.3);
        let mu = 1.0;

        let h0 = Vec3::new(
            r0.y * v0.z - r0.z * v0.y,
            r0.z * v0.x - r0.x * v0.z,
            r0.x * v0.y - r0.y * v0.x,
        );

        let mut r = r0;
        let mut v = v0;
        for _ in 0..1000 {
            let (rn, vn) = kepler_step(r, v, 0.05, mu);
            r = rn;
            v = vn;
        }

        let h1 = Vec3::new(r.y * v.z - r.z * v.y, r.z * v.x - r.x * v.z, r.x * v.y - r.y * v.x);
        let drift = ((h1 - h0).length() / h0.length()).abs();
        assert!(drift < 1e-12, "angular momentum drifted by {drift}");
    }

    /// State with `z = 0`, `vz = 0` confined to the `xy`-plane should
    /// produce a planar trajectory through the [`Vec3`] API. This guards
    /// against a regression where the propagator would inadvertently leak
    /// `z` motion through f64 round-off.
    #[test]
    fn planar_state_stays_planar_through_vec3_api() {
        let r0 = Vec3::new(1.0, 0.0, 0.0);
        let v0 = Vec3::new(0.0, 1.5, 0.0);
        let mu = 1.0;

        let mut r = r0;
        let mut v = v0;
        for _ in 0..500 {
            let (rn, vn) = kepler_step(r, v, 0.07, mu);
            r = rn;
            v = vn;
        }

        assert!(r.z.abs() < 1e-14, "z leaked into planar trajectory: {}", r.z);
        assert!(v.z.abs() < 1e-14, "vz leaked into planar trajectory: {}", v.z);
    }

    /// A reverse step `kepler_step(r, v, dt) → kepler_step(·, ·, −dt)`
    /// reproduces the original state to f64 precision (time-reversal
    /// symmetry of the analytical Kepler flow).
    #[test]
    fn time_reversal_returns_initial_state() {
        let r0 = Vec3::new(1.5, 0.4, -0.2);
        let v0 = Vec3::new(-0.3, 1.1, 0.2);
        let mu = 1.0;
        let dt = 1.5;

        let (r1, v1) = kepler_step(r0, v0, dt, mu);
        let (r0_back, v0_back) = kepler_step(r1, v1, -dt, mu);

        assert!(
            (r0_back - r0).length() < 1e-12,
            "position reversal drift {}",
            (r0_back - r0).length()
        );
        assert!(
            (v0_back - v0).length() < 1e-12,
            "velocity reversal drift {}",
            (v0_back - v0).length()
        );
    }
}
