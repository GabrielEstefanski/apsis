//! Euler's collinear three-body solution (1767).
//!
//! Three bodies on a rotating straight line through the COM, all three
//! orbiting at the same angular velocity. For arbitrary masses
//! `(m₁, m₂, m₃)` the position ratio `λ = |x₃ − x₂| / |x₂ − x₁|` is the
//! unique positive real root of Euler's quintic
//!
//! ```text
//! (m₁ + m₂)·λ⁵ + (3m₁ + 2m₂)·λ⁴ + (3m₁ + m₂)·λ³
//!   − (m₂ + 3m₃)·λ² − (2m₂ + 3m₃)·λ − (m₂ + m₃) = 0
//! ```
//!
//! Once `λ` is known, all three bodies orbit the COM at radii determined
//! by the COM constraint `Σ mᵢ·rᵢ = 0`. Each body's tangential speed is
//! `vᵢ = ω·|rᵢ|`. The lightest body sits farthest from the COM (largest
//! orbit), the heaviest closest.
//!
//! The preset uses `(m₁, m₂, m₃) = (0.1, 1.0, 0.5)`, which gives `λ ≈ 1.5`
//! and produces a visible asymmetric configuration: the light body traces
//! the largest orbit, the heavy middle body the smallest, the medium-mass
//! body in between.
//!
//! Linearly unstable for any mass ratio (Routh's criterion). The
//! configuration disintegrates after a few orbits.

use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem},
};

/// Solve Euler's quintic for `λ` via bisection. The polynomial is
/// monotone increasing in `λ ∈ (0, ∞)` (positive leading coefficient,
/// negative at `λ = 0`), so it has exactly one positive real root.
fn euler_lambda(m1: f64, m2: f64, m3: f64) -> f64 {
    let a5 = m1 + m2;
    let a4 = 3.0 * m1 + 2.0 * m2;
    let a3 = 3.0 * m1 + m2;
    let a2 = -(m2 + 3.0 * m3);
    let a1 = -(2.0 * m2 + 3.0 * m3);
    let a0 = -(m2 + m3);
    let f = |l: f64| a5 * l.powi(5) + a4 * l.powi(4) + a3 * l.powi(3) + a2 * l * l + a1 * l + a0;

    let mut lo = 1.0e-6_f64;
    let mut hi = 100.0_f64;
    let f_lo_sign = f(lo).signum();
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid).signum() == f_lo_sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn three_body_euler_collinear(_seed: u64) -> Template {
    // Asymmetric Euler configuration: light + heavy + medium.
    let m1 = 0.1_f64;
    let m2 = 1.0_f64;
    let m3 = 0.5_f64;

    let lambda = euler_lambda(m1, m2, m3);

    // Place bodies along x with body 2 at origin, body 1 at -1, body 3 at +λ.
    let a = 1.0_f64;
    let x1_raw = -a;
    let x2_raw = 0.0_f64;
    let x3_raw = lambda * a;

    // Shift so the COM is at the origin.
    let m_total = m1 + m2 + m3;
    let x_com = (m1 * x1_raw + m2 * x2_raw + m3 * x3_raw) / m_total;
    let r1 = x1_raw - x_com;
    let r2 = x2_raw - x_com;
    let r3 = x3_raw - x_com;

    // Angular velocity from force balance on body 1: net gravitational
    // pull (toward COM) equals m₁·ω²·|r₁|.
    let pull_on_1 = m2 / (a * a) + m3 / ((1.0 + lambda).powi(2) * a * a);
    let omega = (pull_on_1 / r1.abs()).sqrt();

    // Tangential velocities (counter-clockwise rotation around COM).
    let v1 = omega * r1;
    let v2 = omega * r2;
    let v3 = omega * r3;

    Template {
        name: "Euler Collinear (unstable)",
        description: "Three bodies on a rotating straight line — Euler's 1767 \
                      collinear solution with masses (0.1, 1.0, 0.5). All three \
                      orbit the COM at the same ω; the light body (A) traces the \
                      largest orbit, the heavy middle body (B) the smallest, the \
                      medium body (C) in between. Linearly unstable for any mass \
                      ratio; the configuration disintegrates after a few orbits.",
        bodies: vec![
            TemplateBody {
                name: Some("A (light)"),
                mass: m1,
                position: Some([r1, 0.0, 0.0]),
                velocity: [0.0, v1, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("B (heavy)"),
                mass: m2,
                position: Some([r2, 0.0, 0.0]),
                velocity: [0.0, v2, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("C (medium)"),
                mass: m3,
                position: Some([r3, 0.0, 0.0]),
                velocity: [0.0, v3, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
        ],
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        suggested_dt: Some(0.001),
        suggested_integrator: None,
        units: UnitSystem::dimensionless(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Quintic root for `(0.1, 1.0, 0.5)` is `λ ≈ 1.5`. Verifies the
    /// solver finds a positive root that satisfies the polynomial.
    #[test]
    fn quintic_root_satisfies_polynomial() {
        let lambda = euler_lambda(0.1, 1.0, 0.5);
        assert!(lambda > 0.0 && lambda < 5.0);
        let m1 = 0.1;
        let m2 = 1.0;
        let m3 = 0.5;
        let f = (m1 + m2) * lambda.powi(5)
            + (3.0 * m1 + 2.0 * m2) * lambda.powi(4)
            + (3.0 * m1 + m2) * lambda.powi(3)
            - (m2 + 3.0 * m3) * lambda * lambda
            - (2.0 * m2 + 3.0 * m3) * lambda
            - (m2 + m3);
        assert!(f.abs() < 1.0e-10, "quintic residual {f}");
    }

    /// COM lies at the origin after the shift; force balance on each
    /// body reproduces the same `ω²` (collinearity-preserving condition).
    #[test]
    fn force_balance_yields_consistent_omega() {
        let template = three_body_euler_collinear(0);
        let bodies = &template.bodies;
        let m: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
        let r: Vec<f64> = bodies.iter().map(|b| b.position.unwrap()[0]).collect();
        let m_total: f64 = m.iter().sum();
        let r_com: f64 = m.iter().zip(&r).map(|(mi, ri)| mi * ri).sum::<f64>() / m_total;
        assert!(r_com.abs() < 1.0e-12, "COM not at origin: {r_com}");

        // ω² for each body from |F_net| / (mᵢ · |rᵢ|).
        let mut omega_sq = [0.0_f64; 3];
        for i in 0..3 {
            let mut force = 0.0_f64;
            for j in 0..3 {
                if i == j {
                    continue;
                }
                let dx = r[j] - r[i];
                let f_ij = m[i] * m[j] / (dx * dx) * dx.signum();
                force += f_ij;
            }
            // Net force on body i must point toward COM (origin); sign
            // matches -r[i].
            omega_sq[i] = -force / (m[i] * r[i]);
        }
        // All three ω² must agree to high precision.
        let max = omega_sq.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min = omega_sq.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!((max - min) / max < 1.0e-9, "ω² varies across bodies: {omega_sq:?}");
    }
}
