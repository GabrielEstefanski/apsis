//! 3D orbit camera state and math.
//!
//! [`OrbitCamera`] follows the conventional simulator-grade camera
//! used by Universe Sandbox, KSP map view, REBOUND-viz: a pivot point
//! in world space plus three scalar controls (azimuth, elevation,
//! distance) that locate the eye on a sphere centred on the pivot.
//!
//! Two state copies — `current` and `target` — drive a critically
//! damped spring approach that gives smooth motion under any input
//! cadence. Inertia falls out of this naturally: a sharp gesture
//! injects a target offset that the spring decays to zero over
//! roughly `1 / OMEGA_N` seconds without a separate "fling" path.

pub mod input;

use glam::{DMat4, DVec3};

/// Vertical field of view used by the body pass, in radians.
/// Shared so canvas projection and any framing helper that has to
/// reason about on-screen pixel sizes cannot drift apart.
pub const FOV_Y_RAD: f32 = 0.698_131_7; // 40°.to_radians()
/// Near plane of the reverse-Z perspective. 0.001 AU ≈ 150 000 km —
/// past the surface of any planet, well inside the camera's typical
/// orbit; framing helpers floor distances at a small multiple of this.
pub const NEAR_PLANE: f32 = 0.001;

/// Singularity guard for elevation: at exactly ±π/2 the up-vector
/// degenerates and azimuth becomes ill-defined. Clamping at this
/// margin (≈ 0.057°) is invisible in practice.
const ELEVATION_LIMIT: f64 = std::f64::consts::FRAC_PI_2 - 1e-3;
/// Linear distance lower bound. Keeps the eye outside numerical
/// epsilons of the pivot and prevents division-by-zero in `view_matrix`.
const MIN_DISTANCE: f64 = 1e-6;

/// One camera pose: pivot, spherical-coordinate orientation, and
/// distance. Lives twice in [`OrbitCamera`] (current vs. target).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CameraPose {
    pub pivot: DVec3,
    /// Rotation around the world Y axis, in radians. Positive =
    /// counter-clockwise looking down (right-hand rule about +Y).
    pub azimuth: f64,
    /// Elevation above the XZ plane, in radians. Positive = looking
    /// down on the pivot. Clamped to ±([`ELEVATION_LIMIT`]).
    pub elevation: f64,
    /// Eye-to-pivot distance, in world units. Always ≥ [`MIN_DISTANCE`].
    pub distance: f64,
}

impl CameraPose {
    pub fn new(pivot: DVec3, azimuth: f64, elevation: f64, distance: f64) -> Self {
        Self {
            pivot,
            azimuth,
            elevation: elevation.clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT),
            distance: distance.max(MIN_DISTANCE),
        }
    }

    /// Eye position in world space.
    pub fn eye(&self) -> DVec3 {
        let (sa, ca) = self.azimuth.sin_cos();
        let (se, ce) = self.elevation.sin_cos();
        // Eye sits on the pivot-centred sphere of radius `distance`.
        // Identity orientation (azimuth=elevation=0) places the eye on +Z
        // looking toward −Z, matching wgpu's NDC convention.
        let dir = DVec3::new(ce * sa, se, ce * ca);
        self.pivot + dir * self.distance
    }

    /// Unit vector from eye toward pivot.
    pub fn forward(&self) -> DVec3 {
        (self.pivot - self.eye()).normalize()
    }

    /// World-space right vector for the camera frame. Computed from
    /// world up (+Y) and forward; orthonormal as long as the
    /// elevation clamp holds.
    pub fn right(&self) -> DVec3 {
        self.forward().cross(DVec3::Y).normalize()
    }

    /// Camera-frame up. Recomputed from forward × right so the basis
    /// stays orthonormal regardless of accumulated float error in
    /// azimuth/elevation.
    pub fn up(&self) -> DVec3 {
        self.right().cross(self.forward()).normalize()
    }

    /// Right-handed look-at matrix from eye toward pivot.
    pub fn view_matrix(&self) -> DMat4 {
        DMat4::look_at_rh(self.eye(), self.pivot, DVec3::Y)
    }
}

impl Default for CameraPose {
    fn default() -> Self {
        Self::new(DVec3::ZERO, 0.0, 0.3, 50.0)
    }
}

/// Spring-damper around a target pose. `current` chases `target`; a
/// gesture writes `target` and the spring carries `current` along.
#[derive(Debug, Clone)]
pub struct OrbitCamera {
    pub current: CameraPose,
    pub target: CameraPose,
    /// Natural angular frequency of the spring (rad/s). Higher =
    /// snappier follow. Default tuned to settle visually in roughly
    /// 200 ms (ω_n = 12, ζ = 1, settles in ~4/ω_n).
    pub omega_n: f64,
}

impl OrbitCamera {
    pub fn new(initial: CameraPose) -> Self {
        Self { current: initial, target: initial, omega_n: 12.0 }
    }

    /// Replace the pose immediately, snapping both `current` and
    /// `target`. Used on first frame and on hard recentre.
    pub fn snap(&mut self, pose: CameraPose) {
        self.current = pose;
        self.target = pose;
    }

    /// Advance `current` toward `target` by `dt` seconds using a
    /// critically damped exponential approach. Frame-rate independent.
    pub fn integrate(&mut self, dt: f64) {
        let alpha = 1.0 - (-self.omega_n * dt).exp();
        let c = &mut self.current;
        let t = &self.target;

        c.pivot = c.pivot.lerp(t.pivot, alpha);
        c.azimuth = lerp_angle(c.azimuth, t.azimuth, alpha);
        c.elevation = c.elevation + alpha * (t.elevation - c.elevation);
        c.distance = c.distance + alpha * (t.distance - c.distance);

        c.elevation = c.elevation.clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT);
        c.distance = c.distance.max(MIN_DISTANCE);
    }

    /// Apply gesture deltas to `target`. Caller is responsible for
    /// scaling these by their own sensitivity factors before calling.
    pub fn rotate(&mut self, d_azimuth: f64, d_elevation: f64) {
        self.target.azimuth += d_azimuth;
        self.target.elevation =
            (self.target.elevation + d_elevation).clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT);
    }

    /// Multiplicative zoom. `factor > 1` zooms out. Distance is
    /// updated geometrically so wheel ticks feel uniform across
    /// scales.
    pub fn zoom(&mut self, factor: f64) {
        self.target.distance = (self.target.distance * factor).max(MIN_DISTANCE);
    }

    /// Translate the pivot in the camera's screen plane. `dx` is
    /// along `right`, `dy` is along `up`. Magnitudes are world units.
    pub fn pan_screen(&mut self, dx: f64, dy: f64) {
        let r = self.target.right();
        let u = self.target.up();
        self.target.pivot += r * dx + u * dy;
    }

    /// `true` when `current` is close enough to `target` that further
    /// integration cannot move the eye by a visible amount. Used by
    /// the canvas adapter to stop requesting repaints when the
    /// camera has settled.
    pub fn is_at_rest(&self) -> bool {
        const ANGLE_EPS: f64 = 1e-5;
        const DIST_EPS_REL: f64 = 1e-5;
        const PIVOT_EPS_REL: f64 = 1e-5;

        let c = &self.current;
        let t = &self.target;
        let dist_eps = t.distance.max(1.0) * DIST_EPS_REL;
        let pivot_eps = t.distance.max(1.0) * PIVOT_EPS_REL;

        angle_diff(c.azimuth, t.azimuth).abs() < ANGLE_EPS
            && (c.elevation - t.elevation).abs() < ANGLE_EPS
            && (c.distance - t.distance).abs() < dist_eps
            && (c.pivot - t.pivot).length() < pivot_eps
    }
}

/// Shortest signed angular distance from `a` to `b`, in `(-π, π]`.
/// Sibling of [`lerp_angle`] — both treat radian inputs as living on
/// the unit circle so wrap-around at ±τ is invisible.
fn angle_diff(a: f64, b: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let mut d = (b - a) % two_pi;
    if d > std::f64::consts::PI {
        d -= two_pi;
    } else if d < -std::f64::consts::PI {
        d += two_pi;
    }
    d
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self::new(CameraPose::default())
    }
}

/// Linear interpolation on the unit circle, choosing the shortest
/// arc. Plain lerp on radians wraps the wrong way across ±π.
fn lerp_angle(a: f64, b: f64, t: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let mut d = (b - a) % two_pi;
    if d > std::f64::consts::PI {
        d -= two_pi;
    } else if d < -std::f64::consts::PI {
        d += two_pi;
    }
    a + d * t
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }
    fn vec_approx_eq(a: DVec3, b: DVec3, eps: f64) -> bool {
        (a - b).length() < eps
    }

    #[test]
    fn identity_pose_eye_on_positive_z() {
        let p = CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0);
        assert!(vec_approx_eq(p.eye(), DVec3::new(0.0, 0.0, 10.0), 1e-12));
        assert!(vec_approx_eq(p.forward(), DVec3::NEG_Z, 1e-12));
    }

    #[test]
    fn azimuth_quarter_turn_eye_on_positive_x() {
        let p = CameraPose::new(DVec3::ZERO, FRAC_PI_2, 0.0, 10.0);
        assert!(vec_approx_eq(p.eye(), DVec3::new(10.0, 0.0, 0.0), 1e-12));
    }

    #[test]
    fn elevation_top_down_eye_on_positive_y() {
        let p = CameraPose::new(DVec3::ZERO, 0.0, ELEVATION_LIMIT, 10.0);
        let eye = p.eye();
        assert!(approx_eq(eye.y, 10.0, 1e-2));
    }

    #[test]
    fn frame_is_orthonormal() {
        for &az in &[0.0, 0.7, 1.5, -2.3, std::f64::consts::PI] {
            for &el in &[-1.2, -0.4, 0.0, 0.5, 1.3] {
                let p = CameraPose::new(DVec3::new(1.0, 2.0, 3.0), az, el, 7.5);
                let f = p.forward();
                let r = p.right();
                let u = p.up();
                assert!(approx_eq(f.length(), 1.0, 1e-12));
                assert!(approx_eq(r.length(), 1.0, 1e-12));
                assert!(approx_eq(u.length(), 1.0, 1e-12));
                assert!(approx_eq(f.dot(r), 0.0, 1e-12));
                assert!(approx_eq(f.dot(u), 0.0, 1e-12));
                assert!(approx_eq(r.dot(u), 0.0, 1e-12));
            }
        }
    }

    #[test]
    fn elevation_is_clamped() {
        let p = CameraPose::new(DVec3::ZERO, 0.0, 100.0, 10.0);
        assert!(p.elevation <= ELEVATION_LIMIT);
        let p = CameraPose::new(DVec3::ZERO, 0.0, -100.0, 10.0);
        assert!(p.elevation >= -ELEVATION_LIMIT);
    }

    #[test]
    fn distance_is_floored() {
        let p = CameraPose::new(DVec3::ZERO, 0.0, 0.0, -5.0);
        assert!(p.distance >= MIN_DISTANCE);
    }

    #[test]
    fn integrate_with_zero_dt_is_noop() {
        let initial = CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0);
        let mut cam = OrbitCamera::new(initial);
        cam.target = CameraPose::new(DVec3::ZERO, 1.0, 0.5, 20.0);
        cam.integrate(0.0);
        assert_eq!(cam.current, initial);
    }

    #[test]
    fn integrate_converges_to_target() {
        let initial = CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0);
        let target = CameraPose::new(DVec3::new(1.0, 2.0, 3.0), 0.6, -0.4, 25.0);
        let mut cam = OrbitCamera::new(initial);
        cam.target = target;
        // At ω_n = 12 rad/s, 1 second of integration is e^-12 ≈ 6e-6 of
        // the way out — close enough for sub-mm tolerance on 25-unit
        // distances.
        for _ in 0..600 {
            cam.integrate(1.0 / 60.0);
        }
        assert!(approx_eq(cam.current.azimuth, target.azimuth, 1e-4));
        assert!(approx_eq(cam.current.elevation, target.elevation, 1e-4));
        assert!(approx_eq(cam.current.distance, target.distance, 1e-4));
        assert!(vec_approx_eq(cam.current.pivot, target.pivot, 1e-4));
    }

    #[test]
    fn lerp_angle_takes_shortest_arc() {
        // From −0.1 toward +6.18 (≈ −0.1 in wrapped coords): the
        // shortest path is the negative direction, not +6.28.
        let a = -0.1_f64;
        let b = std::f64::consts::TAU - 0.1; // same point, wrapped
        let mid = lerp_angle(a, b, 0.5);
        // Shortest path stays near the original point.
        assert!(approx_eq(mid, a, 1e-12));
    }

    #[test]
    fn zoom_factor_scales_distance() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.zoom(2.0);
        assert!(approx_eq(cam.target.distance, 20.0, 1e-12));
        cam.zoom(0.25);
        assert!(approx_eq(cam.target.distance, 5.0, 1e-12));
    }

    #[test]
    fn is_at_rest_after_convergence_through_full_turn() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        // Target one full turn past current. Mathematically the same
        // point on the circle; the rest predicate must collapse the
        // 2π gap.
        cam.target.azimuth = std::f64::consts::TAU;
        for _ in 0..600 {
            cam.integrate(1.0 / 60.0);
        }
        assert!(cam.is_at_rest());
    }

    #[test]
    fn pan_screen_translates_pivot_in_camera_frame() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        // At identity pose: right = +X, up = +Y, forward = −Z.
        cam.pan_screen(3.0, 4.0);
        assert!(vec_approx_eq(cam.target.pivot, DVec3::new(3.0, 4.0, 0.0), 1e-12));
    }
}
