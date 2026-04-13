//! Keplerian two-body propagator used by the Wisdom–Holman integrator.
//!
//! All functions operate in the **relative frame** (planet minus central body).
//! Units are internal simulation units where `G = 1`; `mu = G * M_central`.
//!
//! # Formulation
//!
//! The propagator uses the **universal variable** `χ` and **Stumpff functions**
//! `c(ψ)`, `s(ψ)` throughout, following Battin (1987) §4.5.  This single
//! formulation is continuously valid for elliptic (`e < 1`), parabolic
//! (`e = 1`), and hyperbolic (`e > 1`) trajectories without branching on the
//! orbit type.  The classical elliptic Kepler equation is recovered as the
//! special case `ψ = χ²/a > 0`.
//!
//! # References
//!
//! - Battin, R. H. (1987). *An Introduction to the Mathematics and Methods
//!   of Astrodynamics*. AIAA Education Series. §4.5.
//! - Bate, R. R., Mueller, D. D., & White, J. E. (1971).
//!   *Fundamentals of Astrodynamics*. Dover. §2.8.
//! - Wisdom, J. & Holman, M. (1991). *Astron. J.* 102, 1528–1538.

// ── Stumpff functions ─────────────────────────────────────────────────────────

/// Stumpff function `c(ψ) = (1 − cos √ψ) / ψ`.
///
/// Continuously defined for all `ψ ∈ ℝ`:
///
/// | Domain   | Identity used |
/// |----------|---------------|
/// | `ψ > 0`  | `(1 − cos √ψ) / ψ` |
/// | `ψ < 0`  | `(cosh √(−ψ) − 1) / (−ψ)` |
/// | `ψ ≈ 0`  | Taylor series `½ − ψ/24 + ψ²/720 − …` |
fn stumpff_c(psi: f64) -> f64 {
    if psi.abs() < 1e-6 {
        return 0.5 - psi / 24.0 + psi * psi / 720.0;
    }
    if psi > 0.0 {
        (1.0 - psi.sqrt().cos()) / psi
    } else {
        ((-psi).sqrt().cosh() - 1.0) / (-psi)
    }
}

/// Stumpff function `s(ψ) = (√ψ − sin √ψ) / ψ^(3/2)`.
///
/// Continuously defined for all `ψ ∈ ℝ`:
///
/// | Domain   | Identity used |
/// |----------|---------------|
/// | `ψ > 0`  | `(√ψ − sin √ψ) / ψ^(3/2)` |
/// | `ψ < 0`  | `(sinh √(−ψ) − √(−ψ)) / (−ψ)^(3/2)` |
/// | `ψ ≈ 0`  | Taylor series `1/6 − ψ/120 + ψ²/5040 − …` |
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
/// The derivative `df/dχ = r(χ) / √μ` is the instantaneous radius, which
/// is always positive and guarantees monotone convergence of Newton's method.
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
/// Given the relative state vector `(x, y, vx, vy)` of a body with respect
/// to the central mass and the gravitational parameter `mu = G * M_central`,
/// returns the propagated state `(x', y', vx', vy')` via the Lagrange
/// coefficient equations
///
/// ```text
/// r⃗'  = f · r⃗₀  +  g · v⃗₀
/// v⃗'  = ḟ · r⃗₀  +  ġ · v⃗₀
/// ```
///
/// where `f`, `g`, `ḟ`, `ġ` are expressed in terms of `χ` and Stumpff
/// functions, valid for all orbit types without branching.
///
/// # Convergence
///
/// Newton–Raphson on the universal time equation converges quadratically.
/// The derivative `df/dχ = r/√μ > 0` is always positive, guaranteeing
/// monotone convergence. The tolerance is `|Δχ| < 10⁻¹²`.
///
/// # Degenerate inputs
///
/// A zero or near-zero radius `r → 0` indicates a collision singularity;
/// the function returns the input state unchanged in that case.
///
/// # References
///
/// - Battin (1987) §4.5; Bate–Mueller–White (1971) §2.8.
/// - Wisdom & Holman (1991) *Astron. J.* 102, 1528–1538.
pub fn kepler_step(x: f64, y: f64, vx: f64, vy: f64, dt: f64, mu: f64) -> (f64, f64, f64, f64) {
    let r0 = (x * x + y * y).sqrt();
    if r0 < 1e-30 {
        return (x, y, vx, vy);
    }

    let v2 = vx * vx + vy * vy;
    let vr0 = (x * vx + y * vy) / r0;
    let alpha = 2.0 / r0 - v2 / mu; // = 1/a; negative iff hyperbolic
    let sqrt_mu = mu.sqrt();

    // Initial guess for χ (Bate–Mueller–White §2.8, adapted for sign of dt)
    let chi0 = if alpha > 0.0 {
        // Elliptic: χ ≈ n · dt · a,  n = √(μ / a³) = √μ · α^(3/2)
        sqrt_mu * dt * alpha
    } else {
        // Hyperbolic: use the asymptotic chord length as seed
        let a = 1.0 / alpha; // a < 0
        let s = dt.signum() * (2.0 * mu * alpha.abs()).sqrt() * dt.abs()
            / (r0 * vr0 + dt.signum() * (mu * (-a)).sqrt() * (1.0 - r0 * alpha));
        (2.0 * mu * alpha.abs()).sqrt() * s.tanh() / alpha.abs().sqrt()
    };

    // Newton–Raphson on the universal Kepler equation
    let mut chi = chi0;
    for _ in 0..50 {
        let (f_val, df) = kepler_universal_fd(chi, alpha, r0, vr0, mu, dt);
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

    let r1 = (r0 * vr0 / sqrt_mu) * chi * (1.0 - psi * s) + (1.0 - r0 * alpha) * chi2 * c + r0;

    let f_lag = 1.0 - chi2 * c / r0;
    let g_lag = dt - chi3 * s / sqrt_mu;
    let df_lag = -sqrt_mu / (r1 * r0) * chi * (1.0 - psi * s);
    let dg_lag = 1.0 - chi2 * c / r1;

    let x1 = f_lag * x + g_lag * vx;
    let y1 = f_lag * y + g_lag * vy;
    let vx1 = df_lag * x + dg_lag * vx;
    let vy1 = df_lag * y + dg_lag * vy;

    (x1, y1, vx1, vy1)
}
