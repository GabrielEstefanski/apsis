//! Three-dimensional Cartesian vector for positions, velocities, and
//! accelerations in physical space.
//!
//! `Vec3` is the storage and transport type for kinematic quantities
//! throughout the physics stack. Layout is three plain `f64` fields with
//! `#[repr(C)]` so a `&[Vec3]` is bitwise compatible with a `&[[f64; 3]]`
//! and with C/FFI consumers. No SIMD intrinsics — the integrator hot
//! loops are written against the scalars and rely on LLVM
//! autovectorisation.

use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[repr(C)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };

    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// `(x, y, 0)` — bridge for 2D initialisation paths preserved during
    /// the 3D port. Equivalent to `Vec3::new(x, y, 0.0)`.
    #[inline]
    pub const fn from_xy(x: f64, y: f64) -> Self {
        Self { x, y, z: 0.0 }
    }

    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    #[inline]
    pub fn cross(self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    #[inline]
    pub fn length_squared(self) -> f64 {
        self.dot(self)
    }

    #[inline]
    pub fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    /// Squared Euclidean distance to `rhs`.
    ///
    /// Equivalent to `(self - rhs).length_squared()` but spelled directly
    /// at the call site. Prefer this over `.length()` whenever the
    /// comparison only needs `r²` (e.g. neighbour cutoffs, BH opening
    /// criterion) — saves one `sqrt`.
    #[inline]
    pub fn distance_squared(self, rhs: Self) -> f64 {
        (self - rhs).length_squared()
    }

    /// Euclidean distance to `rhs`.
    ///
    /// Equivalent to `(self - rhs).length()`.
    #[inline]
    pub fn distance(self, rhs: Self) -> f64 {
        self.distance_squared(rhs).sqrt()
    }

    /// Unit vector along `self`. **Caller must ensure the vector is non-zero.**
    ///
    /// Returns NaN components when `self == Vec3::ZERO`. Use
    /// [`try_normalize`](Self::try_normalize) when the input is not
    /// already validated by the caller (e.g. force kernels that gate on
    /// `r² < 1e-30` already protect this path; collision normals from
    /// arbitrary input do not).
    #[inline]
    pub fn normalize(self) -> Self {
        self / self.length()
    }

    /// Unit vector along `self`, or `None` when `self` is at or below the
    /// f64 normalisation floor.
    ///
    /// The threshold is chosen so the division `1 / length` cannot
    /// overflow into infinity for any finite, non-subnormal input.
    #[inline]
    pub fn try_normalize(self) -> Option<Self> {
        let len_sq = self.length_squared();
        if len_sq > f64::MIN_POSITIVE { Some(self / len_sq.sqrt()) } else { None }
    }

    #[inline]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, z: -self.z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self { x: self.x * rhs, y: self.y * rhs, z: self.z * rhs }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, rhs: Vec3) -> Vec3 {
        rhs * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f64) -> Self {
        Self { x: self.x / rhs, y: self.y / rhs, z: self.z / rhs }
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl MulAssign<f64> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl From<(f64, f64, f64)> for Vec3 {
    #[inline]
    fn from((x, y, z): (f64, f64, f64)) -> Self {
        Self { x, y, z }
    }
}

impl From<[f64; 3]> for Vec3 {
    #[inline]
    fn from([x, y, z]: [f64; 3]) -> Self {
        Self { x, y, z }
    }
}

impl From<Vec3> for [f64; 3] {
    #[inline]
    fn from(v: Vec3) -> Self {
        [v.x, v.y, v.z]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_is_additive_identity() {
        let v = Vec3::new(1.5, -2.0, 3.25);
        assert_eq!(v + Vec3::ZERO, v);
        assert_eq!(Vec3::ZERO + v, v);
    }

    #[test]
    fn addition_is_componentwise() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(0.5, -1.0, 4.0);
        assert_eq!(a + b, Vec3::new(1.5, 1.0, 7.0));
    }

    #[test]
    fn subtraction_is_componentwise() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(0.5, -1.0, 4.0);
        assert_eq!(a - b, Vec3::new(0.5, 3.0, -1.0));
    }

    #[test]
    fn scalar_mul_commutes() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        assert_eq!(v * 2.5, 2.5 * v);
        assert_eq!(v * 2.5, Vec3::new(2.5, -5.0, 7.5));
    }

    #[test]
    fn dot_matches_definition() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, -5.0, 6.0);
        assert_eq!(a.dot(b), 1.0 * 4.0 + 2.0 * (-5.0) + 3.0 * 6.0);
    }

    #[test]
    fn cross_satisfies_right_hand_rule() {
        let ex = Vec3::new(1.0, 0.0, 0.0);
        let ey = Vec3::new(0.0, 1.0, 0.0);
        let ez = Vec3::new(0.0, 0.0, 1.0);
        assert_eq!(ex.cross(ey), ez);
        assert_eq!(ey.cross(ez), ex);
        assert_eq!(ez.cross(ex), ey);
    }

    #[test]
    fn cross_anticommutes() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!(a.cross(b), -b.cross(a));
    }

    #[test]
    fn cross_with_self_is_zero() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.cross(v), Vec3::ZERO);
    }

    #[test]
    fn length_squared_matches_dot_self() {
        let v = Vec3::new(3.0, 4.0, 12.0);
        assert_eq!(v.length_squared(), v.dot(v));
        assert_eq!(v.length_squared(), 169.0);
        assert_eq!(v.length(), 13.0);
    }

    #[test]
    fn from_xy_zeroes_z() {
        assert_eq!(Vec3::from_xy(1.0, 2.0), Vec3::new(1.0, 2.0, 0.0));
    }

    #[test]
    fn from_tuple_and_array_roundtrip() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(Vec3::from((1.0, 2.0, 3.0)), v);
        assert_eq!(Vec3::from([1.0, 2.0, 3.0]), v);
        assert_eq!(<[f64; 3]>::from(v), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn negation_inverts_addition() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        assert_eq!(v + (-v), Vec3::ZERO);
    }

    #[test]
    fn add_assign_matches_add() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(0.5, -1.0, 4.0);
        let mut c = a;
        c += b;
        assert_eq!(c, a + b);
    }

    #[test]
    fn distance_matches_length_of_difference() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 6.0, 15.0);
        // Δ = (3, 4, 12) → length 13, length_squared 169.
        assert_eq!(a.distance_squared(b), 169.0);
        assert_eq!(a.distance(b), 13.0);
        assert_eq!(a.distance_squared(b), (a - b).length_squared());
        assert_eq!(a.distance(b), (a - b).length());
    }

    #[test]
    fn distance_is_symmetric() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-2.0, 5.0, 1.5);
        assert_eq!(a.distance(b), b.distance(a));
        assert_eq!(a.distance_squared(b), b.distance_squared(a));
    }

    #[test]
    fn normalize_returns_unit_vector() {
        let v = Vec3::new(3.0, 4.0, 12.0); // length 13
        let u = v.normalize();
        assert!((u.length() - 1.0).abs() < 1e-15);
        assert!((u.x - 3.0 / 13.0).abs() < 1e-15);
        assert!((u.y - 4.0 / 13.0).abs() < 1e-15);
        assert!((u.z - 12.0 / 13.0).abs() < 1e-15);
    }

    #[test]
    fn normalize_axis_aligned_recovers_basis() {
        let ex = Vec3::new(2.5, 0.0, 0.0).normalize();
        assert_eq!(ex, Vec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn try_normalize_returns_some_for_nonzero() {
        let v = Vec3::new(3.0, 4.0, 12.0);
        let u = v.try_normalize().unwrap();
        assert!((u.length() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn try_normalize_returns_none_for_zero() {
        assert_eq!(Vec3::ZERO.try_normalize(), None);
    }

    #[test]
    fn try_normalize_returns_none_for_subnormal_length() {
        // length_squared at or below f64::MIN_POSITIVE rejects.
        let tiny = Vec3::new(f64::MIN_POSITIVE.sqrt() * 0.5, 0.0, 0.0);
        assert!(tiny.try_normalize().is_none());
    }

    #[test]
    fn is_finite_rejects_nan_and_inf() {
        assert!(Vec3::new(1.0, 2.0, 3.0).is_finite());
        assert!(!Vec3::new(f64::NAN, 0.0, 0.0).is_finite());
        assert!(!Vec3::new(0.0, f64::INFINITY, 0.0).is_finite());
        assert!(!Vec3::new(0.0, 0.0, f64::NEG_INFINITY).is_finite());
    }

    #[test]
    fn repr_c_layout_matches_array() {
        assert_eq!(std::mem::size_of::<Vec3>(), 3 * std::mem::size_of::<f64>());
        assert_eq!(std::mem::align_of::<Vec3>(), std::mem::align_of::<f64>());
        let v = Vec3::new(1.0, 2.0, 3.0);
        let arr: [f64; 3] = v.into();
        assert_eq!(arr, [1.0, 2.0, 3.0]);
    }
}
