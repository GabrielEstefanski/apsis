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

/// Default near-plane ceiling. 0.001 AU ≈ 150 000 km — past the surface
/// of any planet, well inside the typical solar-system orbit. The
/// runtime [`adaptive_near_plane`] caps at this value when the camera
/// is far from its pivot, and steps down for close-up work; framing
/// helpers floor distances at a small multiple of it.
pub const NEAR_PLANE: f32 = 0.001;

/// Smallest near-plane the depth buffer can carry without losing
/// per-body resolution at solar-system far distances. Reverse-Z with
/// `perspective_infinite_reverse_rh` resolves depths to roughly
/// `near × 2^23` (f32 mantissa); 1e-7 AU near with bodies out to
/// ~10⁴ AU keeps the ratio inside what float depth can distinguish.
const NEAR_PLANE_FLOOR: f32 = 1e-7;

/// Compute the perspective near plane to use for a given camera-to-
/// pivot distance.
///
/// Fixed at [`NEAR_PLANE`] (0.001 AU) every projection clipped bodies
/// at ~150 000 km from the camera. For solar-system overview that is
/// fine, but for close-up work on small bodies the body itself is
/// well inside the clip — Earth's radius is 4.3·10⁻⁵ AU, so any
/// camera distance < 0.005 AU (Earth radius × ~100) shows the body
/// disappearing before the camera even arrives.
///
/// The adaptive plane scales as `distance / 1000`, so the body always
/// sits at least 1000× the near plane away — well inside the clip
/// on every frame regardless of zoom level. Capped above by the
/// original [`NEAR_PLANE`] so distant views keep the depth precision
/// they had before; floored at [`NEAR_PLANE_FLOOR`] so extreme zoom
/// can't collapse the depth buffer.
pub fn adaptive_near_plane(camera_distance: f32) -> f32 {
    (camera_distance * 1e-3).clamp(NEAR_PLANE_FLOOR, NEAR_PLANE)
}

/// Singularity guard for elevation: at exactly ±π/2 the up-vector
/// degenerates and azimuth becomes ill-defined. Clamping at this
/// margin (≈ 0.057°) is invisible in practice.
const ELEVATION_LIMIT: f64 = std::f64::consts::FRAC_PI_2 - 1e-3;
/// Linear distance lower bound. Keeps the eye outside numerical
/// epsilons of the pivot and prevents division-by-zero in `view_matrix`.
const MIN_DISTANCE: f64 = 1e-6;

/// Phase-locked transition between two camera viewpoints.
///
/// Captured on click-to-focus. Each frame the follow loop derives a
/// single fraction `t = 1 - alpha_remaining` and feeds it to
/// [`CameraPose::vanwijk_to`], which evaluates the
/// (pivot, distance, azimuth, elevation) pose along the
/// van Wijk & Nuij (2003) smooth zoom-and-pan path — so the body
/// stays monotonically on its way to centre on screen rather than
/// drifting off-frame mid-zoom (the artefact a separable
/// `linear pivot, log distance` lerp produces).
///
/// The pivot endpoint tracks the body's current position each frame
/// (instead of being pinned to where the body was at click time), so
/// fast-moving targets stay in frame throughout the transition.
#[derive(Debug, Clone, Copy)]
pub struct FollowTransition {
    pub body_idx: usize,
    pub initial: CameraPose,
    pub target_azimuth: f64,
    pub target_elevation: f64,
    pub target_distance: f64,
    /// Decays from 1 (untransitioned) to 0 (settled).
    pub alpha_remaining: f64,
    pub tau: f64,
}

impl FollowTransition {
    /// `alpha_remaining` falls by `1/e` in 150 ms; visually settled
    /// in roughly 3·τ.
    pub const DEFAULT_TAU: f64 = 0.15;
    /// Below this remaining alpha the transition is considered done
    /// and the caller hands off to steady-state feedforward.
    pub const SETTLED_EPS: f64 = 1e-3;

    pub fn capture(
        body_idx: usize,
        initial: CameraPose,
        target_azimuth: f64,
        target_elevation: f64,
        target_distance: f64,
    ) -> Self {
        Self {
            body_idx,
            initial,
            target_azimuth,
            target_elevation,
            target_distance,
            alpha_remaining: 1.0,
            tau: Self::DEFAULT_TAU,
        }
    }

    /// Advance the phase by `dt`. Returns `true` once `alpha_remaining`
    /// drops below [`SETTLED_EPS`].
    pub fn step(&mut self, dt: f64) -> bool {
        self.alpha_remaining *= (-dt / self.tau).exp();
        if self.alpha_remaining < Self::SETTLED_EPS {
            self.alpha_remaining = 0.0;
            true
        } else {
            false
        }
    }

    /// Lerp fraction `t = 1 - alpha_remaining` ∈ [0, 1].
    pub fn t(&self) -> f64 {
        1.0 - self.alpha_remaining
    }

    /// Live re-target: shift the orientation endpoint by gesture deltas
    /// without cancelling the in-flight transition. The remaining
    /// portion of the van Wijk path lerps toward the new endpoint, so
    /// a user starting to orbit mid-transition gets the camera
    /// arriving smoothly at the new orientation rather than the
    /// transition snapping to None and the steady-state pivot-snap
    /// jumping into place.
    pub fn rotate_target(&mut self, d_azimuth: f64, d_elevation: f64) {
        self.target_azimuth += d_azimuth;
        self.target_elevation =
            (self.target_elevation + d_elevation).clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT);
    }

    /// Live re-target: scale the distance endpoint by `factor` while
    /// preserving the configured min-distance floor. Same intent as
    /// [`rotate_target`](Self::rotate_target) — scroll-zoom mid-
    /// transition smoothly modifies where the camera lands instead
    /// of cancelling the transition.
    pub fn zoom_target(&mut self, factor: f64, min_distance: f64) {
        self.target_distance = (self.target_distance * factor).max(min_distance);
    }
}

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

    /// Rotation-only view matrix — `look_at_rh` rebuilt as if the eye
    /// sat at the world origin. The render path runs under Floating
    /// Origin (geometry shifted by `render_origin = eye()` on the CPU
    /// side), so the camera is effectively at the origin and only the
    /// orientation transforms matter.
    pub fn view_rotation_only(&self) -> DMat4 {
        let look_dir = self.pivot - self.eye();
        DMat4::look_at_rh(DVec3::ZERO, look_dir, DVec3::Y)
    }

    /// Phase-locked interpolation between two poses by fraction `t`.
    /// Distance lerps in log space; angles take the shortest arc;
    /// pivot lerps linearly. `t = 0` returns `self`, `t = 1` returns
    /// `other`.
    pub fn lerp_to(&self, other: &CameraPose, t: f64) -> CameraPose {
        let log_a = self.distance.max(MIN_DISTANCE).ln();
        let log_b = other.distance.max(MIN_DISTANCE).ln();
        CameraPose::new(
            self.pivot.lerp(other.pivot, t),
            lerp_angle(self.azimuth, other.azimuth, t),
            self.elevation + t * (other.elevation - self.elevation),
            (log_a + t * (log_b - log_a)).exp(),
        )
    }

    /// van Wijk & Nuij (2003) smooth zoom-and-pan interpolation.
    ///
    /// A separable `(linear pivot, log distance)` lerp produces the
    /// "screen-space bowing" artefact during click-to-focus reframes:
    /// the body's apparent on-screen offset is proportional to
    /// `|pivot - body| / distance`, and the log-distance term collapses
    /// faster than the linear pivot term at intermediate `t`. The body
    /// visibly drifts off-centre by an order of magnitude before
    /// snapping back at `t = 1`.
    ///
    /// The van Wijk path parameterises `(pivot, distance)` jointly so
    /// that the perceived velocity (combined screen-space pan + zoom)
    /// stays constant. With `ρ = √2` (the paper's perception-tuned
    /// default) the artefact disappears entirely. Orientation
    /// (azimuth, elevation) still lerps linearly — angles don't
    /// participate in the bowing.
    ///
    /// Reference: van Wijk & Nuij, "Smooth and Efficient Zooming and
    /// Panning", IEEE InfoVis 2003.
    pub fn vanwijk_to(&self, other: &CameraPose, t: f64) -> CameraPose {
        let azimuth = lerp_angle(self.azimuth, other.azimuth, t);
        let elevation = self.elevation + t * (other.elevation - self.elevation);
        let (pivot, distance) =
            vanwijk_pivot_distance(self.pivot, self.distance, other.pivot, other.distance, t);
        CameraPose::new(pivot, azimuth, elevation, distance)
    }
}

/// `ρ²` for van Wijk's smooth zoom-and-pan path. The paper derives
/// `ρ = √2` from a perception study minimising "perceived
/// instantaneous velocity" along the path.
const VANWIJK_RHO_SQ: f64 = 2.0;

/// van Wijk & Nuij 2003 closed-form path between two viewpoints.
///
/// Returns `(pivot, distance)` at progress `t ∈ [0, 1]` such that
/// `t = 0 → (p0, w0)` and `t = 1 → (p1, w1)`. The path is a
/// hyperbolic curve in the `(u, ln w)` plane that keeps the perceived
/// screen-space velocity constant, eliminating the bowing artefact of
/// the separable lerp.
///
/// Degenerate cases:
/// - `|p1 - p0| → 0` (pure zoom, no pan): falls back to log-lerp on
///   distance, holds pivot. Necessary because `b_i` has `|u_d|` in the
///   denominator.
/// - Path length `S` non-finite or near zero (initial and target
///   numerically coincident): falls back to separable lerp.
fn vanwijk_pivot_distance(p0: DVec3, w0: f64, p1: DVec3, w1: f64, t: f64) -> (DVec3, f64) {
    let u_vec = p1 - p0;
    let u_d = u_vec.length();

    let log_lerp_dist = |t: f64| -> f64 {
        let log_w0 = w0.max(MIN_DISTANCE).ln();
        let log_w1 = w1.max(MIN_DISTANCE).ln();
        (log_w0 + t * (log_w1 - log_w0)).exp()
    };

    if u_d < 1e-9 {
        return (p0, log_lerp_dist(t));
    }

    let rho_sq = VANWIJK_RHO_SQ;
    let rho_4 = rho_sq * rho_sq;

    let b0 = (w1 * w1 - w0 * w0 + rho_4 * u_d * u_d) / (2.0 * w0 * rho_sq * u_d);
    let b1 = (w1 * w1 - w0 * w0 - rho_4 * u_d * u_d) / (2.0 * w1 * rho_sq * u_d);

    let r0 = (-b0 + (b0 * b0 + 1.0).sqrt()).ln();
    let r1 = (-b1 + (b1 * b1 + 1.0).sqrt()).ln();

    let s_total = (r1 - r0) / rho_sq;
    if !s_total.is_finite() || s_total.abs() < 1e-12 {
        return (p0.lerp(p1, t), log_lerp_dist(t));
    }

    let s = t * s_total;
    let arg = rho_sq * s + r0;

    let cosh_r0 = r0.cosh();
    let sinh_r0 = r0.sinh();
    let cosh_arg = arg.cosh();
    let tanh_arg = arg.tanh();

    let w = w0 * cosh_r0 / cosh_arg;
    let u_progress = (w0 / rho_sq) * (cosh_r0 * tanh_arg - sinh_r0);

    let pivot = p0 + u_vec * (u_progress / u_d);
    (pivot, w)
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
    /// snappier follow. Default settles visually in ~170 ms
    /// (ω_n = 24, ζ = 1, ~4/ω_n) — KSP / Universe Sandbox feel.
    /// Per-frame progress at 60 fps: 33 %.
    pub omega_n: f64,
    /// Lower bound on `target.distance`, in world units. The
    /// canvas writes the selected body's physical radius (with a
    /// small safety margin) so the camera can't scroll into / past
    /// the body it's following. `None` falls back to the global
    /// [`MIN_DISTANCE`] floor — appropriate when no body is
    /// selected, since the global floor exists only to keep
    /// `view_matrix` numerically well-defined.
    pub min_distance_floor: Option<f64>,
}

impl OrbitCamera {
    pub fn new(initial: CameraPose) -> Self {
        Self { current: initial, target: initial, omega_n: 24.0, min_distance_floor: None }
    }

    /// Effective minimum distance for the current frame. Combines
    /// the global numeric floor with the optional context-aware
    /// floor (selected body radius) supplied by the canvas.
    #[inline]
    pub fn effective_min_distance(&self) -> f64 {
        match self.min_distance_floor {
            Some(floor) => floor.max(MIN_DISTANCE),
            None => MIN_DISTANCE,
        }
    }

    /// Replace the pose immediately, snapping both `current` and
    /// `target`. Used on first frame and on hard recentre.
    pub fn snap(&mut self, pose: CameraPose) {
        self.current = pose;
        self.target = pose;
    }

    /// Advance `current` toward `target` by `dt` seconds.
    /// Distance lerps in log space so a click reframe feels uniform
    /// across the AU-to-light-year dynamic range.
    pub fn integrate(&mut self, dt: f64) {
        let alpha = 1.0 - (-self.omega_n * dt).exp();
        if alpha == 0.0 {
            return;
        }
        let min_dist = self.effective_min_distance();
        let c = &mut self.current;
        let t = &self.target;

        c.pivot = c.pivot.lerp(t.pivot, alpha);
        c.azimuth = lerp_angle(c.azimuth, t.azimuth, alpha);
        c.elevation += alpha * (t.elevation - c.elevation);

        let log_c = c.distance.max(min_dist).ln();
        let log_t = t.distance.max(min_dist).ln();
        c.distance = (log_c + alpha * (log_t - log_c)).exp();

        c.elevation = c.elevation.clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT);
        c.distance = c.distance.max(min_dist);
    }

    /// Apply rotation gesture deltas. Direct manipulation: `current`
    /// snaps to `target` for the rotated axes so the cursor doesn't
    /// lead the camera. The spring still owns `pivot` (shared with
    /// the follow loop, which needs damped pivot motion to stay
    /// centred under feedforward).
    ///
    /// Industry convention: KSP map-view, Universe Sandbox, Cinemachine
    /// "Damping = 0" on direct gesture axes.
    pub fn rotate(&mut self, d_azimuth: f64, d_elevation: f64) {
        self.target.azimuth += d_azimuth;
        self.target.elevation =
            (self.target.elevation + d_elevation).clamp(-ELEVATION_LIMIT, ELEVATION_LIMIT);
        self.current.azimuth = self.target.azimuth;
        self.current.elevation = self.target.elevation;
    }

    /// Multiplicative zoom. `factor > 1` zooms out. Distance is
    /// updated geometrically so wheel ticks feel uniform across
    /// scales. Direct manipulation: `current.distance` snaps to
    /// `target.distance` (see [`rotate`](Self::rotate) for rationale).
    ///
    /// Floors at [`effective_min_distance`](Self::effective_min_distance)
    /// so the camera can't scroll into / past the body it's
    /// following — the canvas writes the selected body's physical
    /// radius into [`min_distance_floor`](Self::min_distance_floor)
    /// so a wheel tick that would land inside the body lands at its
    /// surface instead.
    pub fn zoom(&mut self, factor: f64) {
        let min_dist = self.effective_min_distance();
        self.target.distance = (self.target.distance * factor).max(min_dist);
        self.current.distance = self.target.distance;
    }

    /// Translate the pivot along the camera's right/up axes by world-
    /// unit deltas.
    pub fn pan_pivot(&mut self, dx: f64, dy: f64) {
        let r = self.target.right();
        let u = self.target.up();
        self.target.pivot += r * dx + u * dy;
    }

    /// Pivot target that cancels the spring's discrete-step lag for
    /// motion at most quadratic in time. The factor `(1/α − 1)`
    /// (with `α = 1 - exp(-ω·dt)`) reduces to `1/ω` in the
    /// `dt → 0` limit while staying exact at the actual frame rate.
    /// Inputs are wall-time derivatives sampled *after* the
    /// integrator step that produced `body_pos` (matches what the
    /// physics thread publishes), hence the negative sign on the
    /// `½·a·dt²` term: it removes the `a·dt²` contribution that
    /// `wall_vel` already carries vs. the velocity that drove the
    /// last position step.
    pub fn feedforward_pivot(
        &self,
        dt: f64,
        body_pos: DVec3,
        wall_vel: DVec3,
        wall_acc: DVec3,
    ) -> DVec3 {
        let omega = self.omega_n.max(1e-9);
        let alpha = 1.0 - (-omega * dt).exp();
        if alpha < 1e-12 {
            return body_pos;
        }
        let k = 1.0 / alpha - 1.0;
        body_pos + (wall_vel * dt - wall_acc * (0.5 * dt * dt)) * k
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
        // 10 s of integration at ω_n = 24 leaves `exp(-240)` of the
        // initial gap — well below any sensible tolerance.
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
    fn pan_pivot_translates_along_camera_right_and_up() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.pan_pivot(3.0, 4.0);
        assert!(vec_approx_eq(cam.target.pivot, DVec3::new(3.0, 4.0, 0.0), 1e-12));
    }

    // ── Direct gesture snap ──────────────────────────────────────────────────

    #[test]
    fn rotate_snaps_current_to_target_for_azimuth_and_elevation() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.rotate(0.4, -0.2);
        assert!(approx_eq(cam.current.azimuth, cam.target.azimuth, 1e-12));
        assert!(approx_eq(cam.current.elevation, cam.target.elevation, 1e-12));
        assert!(approx_eq(cam.current.azimuth, 0.4, 1e-12));
        assert!(approx_eq(cam.current.elevation, -0.2, 1e-12));
    }

    #[test]
    fn zoom_snaps_current_to_target_distance() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.zoom(0.5);
        assert!(approx_eq(cam.current.distance, cam.target.distance, 1e-12));
        assert!(approx_eq(cam.current.distance, 5.0, 1e-12));
    }

    #[test]
    fn pan_pivot_does_not_snap_current() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.pan_pivot(3.0, 4.0);
        assert_eq!(cam.current.pivot, DVec3::ZERO);
        assert_eq!(cam.target.pivot, DVec3::new(3.0, 4.0, 0.0));
    }

    // ── Log-space distance lerp ──────────────────────────────────────────────

    #[test]
    fn distance_log_lerp_uniform_per_step_ratio() {
        // 10 → 0.001 in log-space: each integration step covers the
        // same RATIO of the remaining log-distance, so 60 frames
        // cover most of the 4-decade gap.
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.target.distance = 0.001;
        for _ in 0..60 {
            cam.integrate(1.0 / 60.0);
        }
        // After 1 s at ω_n = 24 the spring covers 1 - exp(-24) ≈ 1.
        // In log space that means 4 decades of progress.
        let progress_decades = (10.0_f64.log10() - cam.current.distance.log10())
            / (10.0_f64.log10() - 0.001_f64.log10());
        assert!(progress_decades > 0.999, "log progress = {progress_decades}");
    }

    #[test]
    fn distance_log_lerp_no_overshoot() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 100.0));
        cam.target.distance = 1.0;
        for _ in 0..600 {
            cam.integrate(1.0 / 60.0);
            assert!(
                cam.current.distance >= 1.0 - 1e-6,
                "distance {} undershot target",
                cam.current.distance
            );
        }
    }

    // ── PD-with-feedforward ──────────────────────────────────────────────────

    #[test]
    fn feedforward_zero_for_static_body() {
        let cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        let pos = DVec3::new(1.0, 2.0, 3.0);
        let ff = cam.feedforward_pivot(1.0 / 60.0, pos, DVec3::ZERO, DVec3::ZERO);
        assert!(vec_approx_eq(ff, pos, 1e-12));
    }

    #[test]
    fn feedforward_cancels_constant_velocity_lag() {
        // Discrete-aware feedforward should pin the spring to the
        // body for any constant-velocity motion at the actual dt.
        let dt = 1.0 / 60.0;
        let v = DVec3::new(2.0, 0.0, 0.0);
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));

        let mut body_pos = DVec3::ZERO;
        for _ in 0..120 {
            body_pos += v * dt;
            let ff = cam.feedforward_pivot(dt, body_pos, v, DVec3::ZERO);
            cam.target.pivot = ff;
            cam.integrate(dt);
        }
        assert!(
            vec_approx_eq(cam.current.pivot, body_pos, 1e-9),
            "current {:?} vs body {:?}",
            cam.current.pivot,
            body_pos,
        );
    }

    #[test]
    fn feedforward_cancels_constant_acceleration_lag() {
        // Body advances Verlet-style — `body += v·dt + ½·a·dt²` —
        // matching the discrete model the feedforward inverts.
        let dt = 1.0 / 60.0;
        let a = DVec3::new(0.0, 1.0, 0.0);
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));

        let mut body_pos = DVec3::ZERO;
        let mut body_vel = DVec3::ZERO;
        for _ in 0..120 {
            body_pos += body_vel * dt + a * (0.5 * dt * dt);
            body_vel += a * dt;
            let ff = cam.feedforward_pivot(dt, body_pos, body_vel, a);
            cam.target.pivot = ff;
            cam.integrate(dt);
        }
        assert!(
            vec_approx_eq(cam.current.pivot, body_pos, 1e-9),
            "current {:?} vs body {:?}",
            cam.current.pivot,
            body_pos,
        );
    }

    // ── FollowTransition (phase-locked) ──────────────────────────────────────

    fn pose_at(distance: f64) -> CameraPose {
        CameraPose::new(DVec3::ZERO, 0.0, 0.0, distance)
    }

    #[test]
    fn follow_transition_alpha_decays_at_expected_rate() {
        // After one τ ≈ 9 frames at 60 fps, alpha_remaining ≈ e⁻¹.
        let mut t = FollowTransition::capture(0, pose_at(10.0), 0.0, 0.0, 1.0);
        let dt = 1.0 / 60.0;
        let frames = (FollowTransition::DEFAULT_TAU / dt).round() as u32;
        for _ in 0..frames {
            t.step(dt);
        }
        assert!(
            (t.alpha_remaining - (-1.0_f64).exp()).abs() < 0.02,
            "alpha_remaining = {}",
            t.alpha_remaining,
        );
    }

    #[test]
    fn follow_transition_settles_within_two_seconds() {
        let mut t = FollowTransition::capture(0, pose_at(10.0), 0.0, 0.0, 1.0);
        let dt = 1.0 / 60.0;
        let mut steps = 0;
        while !t.step(dt) {
            steps += 1;
            assert!(steps < 120, "transition failed to settle in 2 s");
        }
        assert_eq!(t.alpha_remaining, 0.0);
    }

    #[test]
    fn follow_transition_t_is_zero_at_capture_one_at_settle() {
        let mut t = FollowTransition::capture(0, pose_at(10.0), 0.0, 0.0, 1.0);
        assert!((t.t() - 0.0).abs() < 1e-12);
        while !t.step(1.0 / 60.0) {}
        assert!((t.t() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pose_lerp_distance_midpoint_is_geometric_mean() {
        // 100 → 0.01 in log space gives a midpoint of 1.0
        // (geometric mean), not 50.005 (arithmetic).
        let initial = CameraPose::new(DVec3::new(10.0, 0.0, 0.0), 0.0, 0.0, 100.0);
        let target = CameraPose::new(DVec3::new(0.0, 5.0, 0.0), 1.0, 0.5, 0.01);
        let mid = initial.lerp_to(&target, 0.5);
        assert!((mid.distance - 1.0).abs() < 1e-12, "midpoint distance = {}", mid.distance);
    }

    #[test]
    fn pose_lerp_endpoints_match_inputs() {
        let a = CameraPose::new(DVec3::new(1.0, 2.0, 3.0), 0.4, -0.2, 5.0);
        let b = CameraPose::new(DVec3::new(7.0, 8.0, 9.0), -0.6, 1.0, 50.0);
        let at0 = a.lerp_to(&b, 0.0);
        let at1 = a.lerp_to(&b, 1.0);
        assert!(vec_approx_eq(at0.pivot, a.pivot, 1e-9));
        assert!((at0.distance - a.distance).abs() < 1e-9);
        assert!(vec_approx_eq(at1.pivot, b.pivot, 1e-9));
        assert!((at1.distance - b.distance).abs() < 1e-9);
    }

    // ── van Wijk smooth zoom-and-pan ─────────────────────────────────────────

    #[test]
    fn vanwijk_endpoints_match_inputs() {
        let a = CameraPose::new(DVec3::new(10.0, 0.0, 0.0), 0.3, -0.1, 100.0);
        let b = CameraPose::new(DVec3::new(0.0, 0.0, 0.0), 0.3, -0.1, 1.0);
        let at0 = a.vanwijk_to(&b, 0.0);
        let at1 = a.vanwijk_to(&b, 1.0);
        assert!(vec_approx_eq(at0.pivot, a.pivot, 1e-9), "t=0 pivot");
        assert!((at0.distance - a.distance).abs() < 1e-9, "t=0 distance");
        assert!(vec_approx_eq(at1.pivot, b.pivot, 1e-9), "t=1 pivot");
        assert!((at1.distance - b.distance).abs() < 1e-9, "t=1 distance");
    }

    #[test]
    fn vanwijk_pure_zoom_falls_back_to_log_lerp() {
        // No pan: pivot identical at endpoints. Path collapses to log
        // distance lerp (geometric mean at midpoint).
        let p = DVec3::new(1.0, 2.0, 3.0);
        let a = CameraPose::new(p, 0.0, 0.0, 100.0);
        let b = CameraPose::new(p, 0.0, 0.0, 0.01);
        let mid = a.vanwijk_to(&b, 0.5);
        assert!(vec_approx_eq(mid.pivot, p, 1e-9), "pivot held");
        assert!((mid.distance - 1.0).abs() < 1e-9, "geometric mean = {}", mid.distance);
    }

    #[test]
    fn vanwijk_screen_offset_is_monotonic_under_aggressive_zoom() {
        // The exact scenario the separable lerp mishandles: 100 AU →
        // 1 AU zoom while panning by 1 unit. Screen offset is
        // proportional to `|pivot - body| / distance`. With the
        // separable lerp this peaks ~5–6× the initial value at
        // intermediate `t`. The van Wijk path keeps it monotonically
        // decreasing.
        let body = DVec3::new(0.0, 0.0, 0.0);
        let initial = CameraPose::new(DVec3::new(1.0, 0.0, 0.0), 0.0, 0.0, 100.0);
        let target = CameraPose::new(body, 0.0, 0.0, 1.0);

        let screen_offset = |pose: &CameraPose| (pose.pivot - body).length() / pose.distance;

        let initial_offset = screen_offset(&initial);
        let mut prev = initial_offset;
        for k in 1..=100 {
            let t = k as f64 / 100.0;
            let pose = initial.vanwijk_to(&target, t);
            let off = screen_offset(&pose);
            assert!(off <= prev + 1e-9, "screen offset grew at t={t}: {prev} → {off}");
            assert!(
                off <= initial_offset + 1e-9,
                "screen offset exceeded initial at t={t}: {off} > {initial_offset}",
            );
            prev = off;
        }
    }

    #[test]
    fn vanwijk_separable_lerp_does_bow_for_comparison() {
        // Documents the artefact `vanwijk_to` exists to fix.
        // Same scenario as the monotonicity test above, evaluated
        // under the separable `lerp_to`. Screen offset peaks well
        // above the initial value — that peak is what the user sees
        // as the body "fugindo do centro" mid-zoom.
        let body = DVec3::new(0.0, 0.0, 0.0);
        let initial = CameraPose::new(DVec3::new(1.0, 0.0, 0.0), 0.0, 0.0, 100.0);
        let target = CameraPose::new(body, 0.0, 0.0, 1.0);

        let screen_offset = |pose: &CameraPose| (pose.pivot - body).length() / pose.distance;
        let initial_offset = screen_offset(&initial);

        let mut peak = 0.0_f64;
        for k in 1..=100 {
            let t = k as f64 / 100.0;
            let pose = initial.lerp_to(&target, t);
            peak = peak.max(screen_offset(&pose));
        }
        // Sanity: the separable lerp peak is at least 3× the initial
        // offset for this scenario (analytic ≈ 5–6×).
        assert!(
            peak > initial_offset * 3.0,
            "expected separable-lerp bow > 3× initial; peak = {peak}, initial = {initial_offset}",
        );
    }

    // ── Adaptive near plane ──────────────────────────────────────────────────

    #[test]
    fn adaptive_near_plane_caps_at_legacy_value_for_far_views() {
        // Distant overview (camera 30 AU from pivot): retain the
        // original 0.001 AU near plane so the depth buffer keeps
        // the precision distant-orbit body rendering relies on.
        for distance in [10.0, 30.0, 100.0, 1000.0] {
            assert_eq!(adaptive_near_plane(distance), NEAR_PLANE);
        }
    }

    #[test]
    fn adaptive_near_plane_scales_down_for_close_views() {
        // Close-up on a small body (camera ~1 Earth radius from
        // pivot): near plane shrinks proportionally so the body
        // doesn't get clipped before the camera arrives.
        let earth_radius_au = 4.3e-5_f32;
        let near = adaptive_near_plane(earth_radius_au * 100.0);
        assert!(near < NEAR_PLANE, "should be below the legacy near at {earth_radius_au}*100 AU");
        assert!(near >= NEAR_PLANE_FLOOR, "should clamp at the depth-precision floor");
    }

    #[test]
    fn adaptive_near_plane_clamps_at_floor_for_pathologically_close_views() {
        // Below the floor the depth buffer collapses — clamp to
        // preserve usable f32 depth precision.
        for distance in [1e-12, 1e-10, NEAR_PLANE_FLOOR / 2.0] {
            assert_eq!(adaptive_near_plane(distance), NEAR_PLANE_FLOOR);
        }
    }

    // ── Min-distance floor (selected-body radius) ────────────────────────────

    #[test]
    fn effective_min_distance_falls_back_to_global_when_unset() {
        let cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        assert_eq!(cam.effective_min_distance(), MIN_DISTANCE);
    }

    #[test]
    fn effective_min_distance_uses_floor_when_set() {
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.min_distance_floor = Some(0.5);
        assert_eq!(cam.effective_min_distance(), 0.5);
    }

    #[test]
    fn effective_min_distance_never_below_global_floor() {
        // A floor smaller than `MIN_DISTANCE` would let the eye
        // collapse past the numerical guard. Guard wins.
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.min_distance_floor = Some(1e-12);
        assert_eq!(cam.effective_min_distance(), MIN_DISTANCE);
    }

    #[test]
    fn zoom_clamps_at_min_distance_floor() {
        // Wheel-zoom that would otherwise pull the eye into the body
        // lands at the body's surface (5 % beyond the radius).
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 10.0));
        cam.min_distance_floor = Some(2.0);
        cam.zoom(0.001); // would naively go to 0.01
        assert_eq!(cam.target.distance, 2.0, "zoom should clamp at floor");
        assert_eq!(cam.current.distance, 2.0, "current should match target post-snap");
    }

    // ── FollowTransition live re-target ──────────────────────────────────────

    fn sample_transition() -> FollowTransition {
        FollowTransition::capture(
            0,
            CameraPose::new(DVec3::new(5.0, 0.0, 0.0), 0.1, 0.0, 100.0),
            0.5,
            0.3,
            10.0,
        )
    }

    #[test]
    fn rotate_target_shifts_endpoint_and_clamps_elevation() {
        let mut t = sample_transition();
        t.rotate_target(0.2, 0.1);
        assert!(approx_eq(t.target_azimuth, 0.7, 1e-12));
        assert!(approx_eq(t.target_elevation, 0.4, 1e-12));

        // Push elevation past the singularity guard; clamp must hold.
        let mut t = sample_transition();
        t.rotate_target(0.0, 100.0);
        assert!(t.target_elevation <= ELEVATION_LIMIT);
        let mut t = sample_transition();
        t.rotate_target(0.0, -100.0);
        assert!(t.target_elevation >= -ELEVATION_LIMIT);
    }

    #[test]
    fn rotate_target_does_not_modify_alpha_remaining() {
        // Live re-target must keep the transition's progress fraction
        // intact — only the endpoint moves. Otherwise the camera
        // would either snap forward (alpha drops) or restart (alpha
        // climbs back to 1) on every gesture frame.
        let mut t = sample_transition();
        let alpha_before = t.alpha_remaining;
        t.rotate_target(0.5, 0.2);
        assert_eq!(t.alpha_remaining, alpha_before);
    }

    #[test]
    fn zoom_target_scales_distance_with_floor() {
        let mut t = sample_transition();
        t.zoom_target(0.5, 1e-6);
        assert!(approx_eq(t.target_distance, 5.0, 1e-12));

        // Floor binds: factor that would land below the floor is
        // clamped at floor.
        let mut t = sample_transition();
        t.zoom_target(0.001, 0.5);
        assert_eq!(t.target_distance, 0.5);
    }

    #[test]
    fn zoom_target_does_not_modify_alpha_remaining() {
        let mut t = sample_transition();
        let alpha_before = t.alpha_remaining;
        t.zoom_target(0.5, 0.0);
        assert_eq!(t.alpha_remaining, alpha_before);
    }

    #[test]
    fn integrate_keeps_distance_above_floor() {
        // Spring chase from a low-distance starting point: floor
        // must hold even when the log-lerp would land below it
        // (e.g., target was already clamped by an earlier zoom).
        let mut cam = OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, 5.0));
        cam.target.distance = 1.0;
        cam.min_distance_floor = Some(2.0);
        // Big dt so the spring covers most of the gap.
        cam.integrate(1.0);
        assert!(cam.current.distance >= 2.0, "spring landed at {}", cam.current.distance);
    }
}
