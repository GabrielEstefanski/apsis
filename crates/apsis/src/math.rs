//! Mathematical primitives shared across the physics stack.

pub mod compensated;
pub mod vec3;

pub use compensated::CompensatedF64;
pub use vec3::Vec3;

/// Wraps an angle into the canonical `[-π, π]` range.
///
/// Robust under repeated application: `wrap_pi(wrap_pi(x)) == wrap_pi(x)`
/// for any finite `x`. Apply only to results of *composite* angular
/// arithmetic (sums, differences) — `atan2` already returns in
/// `[-π, π]` and re-wrapping it adds a `(angle + π).rem_euclid(2π) − π`
/// round-trip that perturbs the result by 1 ULP.
#[inline]
pub fn wrap_pi(angle: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let wrapped = (angle + PI).rem_euclid(TAU);
    wrapped - PI
}

#[cfg(test)]
mod tests {
    use super::wrap_pi;
    use std::f64::consts::PI;

    #[test]
    fn wrap_pi_is_identity_inside_range() {
        for a in [-PI + 1e-9, -1.0, 0.0, 1.0, PI - 1e-9] {
            assert!((wrap_pi(a) - a).abs() < 1e-15, "wrap_pi({a}) drifted from itself");
        }
    }

    #[test]
    fn wrap_pi_is_idempotent() {
        for a in [-10.0, -PI, -2.5, 0.0, 2.5, PI, 10.0] {
            let once = wrap_pi(a);
            let twice = wrap_pi(once);
            assert_eq!(twice, once, "wrap_pi not idempotent at a = {a}");
        }
    }

    #[test]
    fn wrap_pi_brings_large_angles_into_range() {
        for a in [-100.0, -3.0 * PI, 3.0 * PI, 100.0] {
            let w = wrap_pi(a);
            assert!(w > -PI - 1e-9 && w <= PI + 1e-9, "wrap_pi({a}) = {w} out of range");
        }
    }

    #[test]
    fn wrap_pi_preserves_angle_modulo_tau() {
        // Wrapping cannot shift the angle by anything other than a
        // multiple of 2π — its purpose is range, not value.
        for a in [-7.5, -2.0, 0.7, 4.2, 9.0] {
            let w = wrap_pi(a);
            let diff = (a - w) / std::f64::consts::TAU;
            assert!(
                (diff - diff.round()).abs() < 1e-12,
                "wrap_pi({a}) shifted by non-multiple of 2π"
            );
        }
    }
}
