//! Dense output — sub-step position, velocity, and acceleration interpolation
//! for smooth rendering.
//!
//! Each completed integration step (or sub-step for IAS15) records a
//! [`DenseSnapshot`] that captures the state needed to evaluate body
//! positions, velocities, and accelerations at any time
//! `t ∈ [t₀, t₀ + dt]` without re-running physics.
//!
//! # Interpolation formulas
//!
//! | Integrator | Position | Velocity | Acceleration |
//! |------------|----------|----------|--------------|
//! | IAS15      | Rein & Spiegel (2015) polynomial via b-coefficients (eq. 9) | derivative of position polynomial (eq. 11) | second derivative — the b-coefficients ARE the higher-order acceleration terms |
//! | WH         | Kepler-analytical: each non-central body propagated via `kepler_step(q₀, v₀, h·dt, μ)` in the rest frame; central body reconstructed from barycenter conservation; Galilean shift back | from the same Kepler propagation (returns the velocity at `h`) | inertial gravitational acceleration at the propagated configuration |
//! | VV / Y4    | 2nd-order Taylor: `x₀ + v₀·h·dt + ½·a₀·(h·dt)²` | analytical derivative: `v₀ + a₀·h·dt` | constant within step: `a₀` |
//!
//! The IAS15 polynomial is exact to the precision of the accepted b-coefficients.
//! The WH Kepler-analytical kernel is exact for the central-force motion;
//! deviation from the integrator's post-step state at `h = 1` is bounded by
//! the perturbation-kick magnitude (O(planet-planet acc · dt)), which is
//! the same truncation order WH itself uses to apportion work between H_K
//! and H_I, and visually invisible at solar-system-scale dt. The Order-2
//! fallback for VV / Y4 is sufficient for smooth visual rendering — those
//! integrators have no analytical drift kernel of their own and the
//! quadratic prediction matches the leapfrog post-step exactly at `h = 1`.
//!
//! # Usage
//!
//! ```ignore
//! let h = (t_render - snap.t0) / snap.dt;   // ∈ [0, 1]
//! let p = snap.interpolate(body_idx, h.clamp(0.0, 1.0));
//! let v = snap.velocity_at(body_idx, h.clamp(0.0, 1.0));
//! let a = snap.acceleration_at(body_idx, h.clamp(0.0, 1.0));
//! ```

use crate::math::Vec3;
use crate::physics::integrator::IntegratorKind;
use crate::physics::integrator::kepler::kepler_step;

// ── DenseSnapshot ─────────────────────────────────────────────────────────────

/// Per-body IAS15 b-coefficients captured at the end of an accepted sub-step.
/// Laid out as `[Vec3; 7]` matching the seven Gauss-Radau nodes.
pub type DenseCoeffs = [Vec3; 7];

/// Wisdom–Holman dense-output state carried alongside the standard
/// snapshot fields. When present, sub-step interpolation evaluates each
/// planet's heliocentric trajectory analytically via [`kepler_step`]
/// rather than via the order-2 Taylor fallback that the [`v0`, `a0`]
/// triple supports.
///
/// Order-2 Taylor on a Keplerian orbit is faithful when the per-step
/// orbital fraction is `≪ 1` (Earth at solar dt). It bows visibly for
/// bodies that cover an appreciable arc per step (Galilean moons,
/// Phobos at default dt) — the curvature of the orbit is not in the
/// model. Replacing the planet-position kernel with the same Kepler
/// propagator the integrator's drift step uses removes that bow at
/// the cost of one Newton-Raphson universal-variable iteration per
/// planet per sample, which renders flat in practice.
///
/// The interpolation reproduces the Wisdom–Holman split structure
/// minus the symmetric kicks and indirect drift — those are
/// perturbation-scale and dominated by the Kepler motion at the dt
/// range WH typically operates in. The error at `h = 1` against the
/// post-step state is bounded by the same perturbation magnitude
/// times `dt`, and visually invisible for hierarchical Solar-System-
/// scale configurations.
#[derive(Clone)]
pub struct WhDenseData {
    /// `G · m_central`. Read by `kepler_step`.
    pub mu: f64,
    /// Central body mass (rest frame); needed for total-momentum
    /// reconstruction of the central-body velocity at sub-step time.
    pub m_sun: f64,
    /// `Σ m_i` over all bodies. Drives the barycenter constraint that
    /// derives the central body's position from the planet
    /// configuration at sub-step time.
    pub m_total: f64,
    /// Centre-of-mass velocity in the original (pre-WH) frame. Each
    /// body's inertial state is the rest-frame Kepler propagation
    /// shifted by `v_com · h · dt`.
    pub v_com: Vec3,
    /// Central body's position in the rest frame at step entry.
    /// `r_0(h)` derives from this via the barycenter constraint
    /// `r_0(h) = r_0(0) + (m_q_in − m_q_h) / m_total`.
    pub r0_sun_rest: Vec3,
    /// `Σ m_i q_i_helio` evaluated at step entry. Constant for the
    /// life of the snapshot — frozen value enters the barycenter
    /// reconstruction at every `h`.
    pub m_q_in: Vec3,
    /// Heliocentric positions of the non-central bodies at step entry,
    /// in the rest frame. One entry per planet, indexed `0..N-1` to
    /// match the `bodies[1..N]` slice WH operates on.
    pub q0_helio_rest: Vec<Vec3>,
    /// Inertial velocities of the non-central bodies at step entry,
    /// in the rest frame. Same indexing convention as
    /// [`q0_helio_rest`](Self::q0_helio_rest).
    pub v0_inertial_rest: Vec<Vec3>,
    /// Masses of the non-central bodies in declaration order. Carried
    /// alongside the kinematic state so the barycenter reconstruction
    /// can sum mass-weighted positions without an external lookup.
    pub planet_masses: Vec<f64>,
}

/// Bulk interpolation result: every body's `(position, velocity,
/// acceleration)` triple at one sub-step `h`, indexed `0..N` in the
/// system's body order.
///
/// Returned by [`WhDenseData::interpolate_kinematics`] so the renderer
/// can pay the O(N) Kepler-propagation cost once per render frame
/// instead of repeating it inside per-body queries.
#[derive(Debug, Clone)]
pub struct WhKinematics {
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    pub accelerations: Vec<Vec3>,
}

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

    /// Wisdom–Holman dense-output state. `Some` only when WH produced
    /// the snapshot; the per-body interpolation methods route through
    /// the analytical Kepler kernel when present, falling back to the
    /// `(x0, v0, a0)` order-2 Taylor otherwise.
    pub wh_data: Option<WhDenseData>,
}

impl WhDenseData {
    /// Compute every body's `(position, velocity, acceleration)` triple
    /// at normalised time `h ∈ [0, 1]` along the Wisdom–Holman step.
    ///
    /// Each non-central body is propagated by the same universal-
    /// variable Kepler solver the integrator's drift step uses
    /// (`kepler_step`), evaluated at `h · dt` from the rest-frame
    /// step-entry state. The central body's position is reconstructed
    /// from the barycenter constraint at the propagated planet
    /// configuration; its velocity from total-momentum conservation in
    /// the rest frame. A final Galilean shift carries the rest-frame
    /// state back to the original (pre-WH) frame.
    ///
    /// Accelerations are evaluated as the inertial gravitational acc
    /// at the propagated configuration:
    ///   * non-central: `−μ q_i(h) / |q_i(h)|³` (Sun's pull only;
    ///     planet-planet perturbation is the kick term and is omitted
    ///     for sub-step interpolation, which mirrors the same
    ///     truncation the position kernel makes by skipping the
    ///     KDK kick + indirect drift).
    ///   * central: `Σ G m_i q_i(h) / |q_i(h)|³`, satisfying
    ///     `Σ m_i a_i = 0` by construction.
    ///
    /// O(N) total work per call. Renderers that need every body's
    /// state at the same `h` should call this once per frame and
    /// index the returned vectors, rather than dispatching through
    /// the per-body methods on [`DenseSnapshot`] (which would each
    /// repeat the full O(N) sweep).
    pub fn interpolate_kinematics(&self, h: f64, dt: f64) -> WhKinematics {
        let n_planets = self.q0_helio_rest.len();
        let n = n_planets + 1;
        let dt_sub = h * dt;

        // Step 1: Kepler propagation per planet. q_h is heliocentric
        // (rest frame); v_h is inertial (rest frame).
        let mut q_helio = Vec::with_capacity(n_planets);
        let mut v_inertial_rest = Vec::with_capacity(n_planets);
        for (q0, v0) in self.q0_helio_rest.iter().zip(&self.v0_inertial_rest) {
            let (q_h, v_h) = kepler_step(*q0, *v0, dt_sub, self.mu);
            q_helio.push(q_h);
            v_inertial_rest.push(v_h);
        }

        // Step 2: barycenter reconstruction for the central body's
        // rest-frame position. In the rest frame Q_0 (mass-weighted
        // total) is invariant, so r_0_rest(h) = r_0(0) + (m_q_in
        // − m_q_h) / m_total.
        let m_q_h: Vec3 =
            q_helio.iter().zip(&self.planet_masses).fold(Vec3::ZERO, |s, (q, m)| s + *m * *q);
        let r0_rest = self.r0_sun_rest + (self.m_q_in - m_q_h) / self.m_total;

        // Step 3: total-momentum conservation in the rest frame
        // (Σ m_i v_i = 0) ⇒ v_0_rest = −(1/m_sun) Σ m_i v_i.
        let p_planets: Vec3 = v_inertial_rest
            .iter()
            .zip(&self.planet_masses)
            .fold(Vec3::ZERO, |s, (v, m)| s + *m * *v);
        let v0_rest = -p_planets / self.m_sun;

        // Step 4: Galilean shift back to the original frame.
        let dr_com = self.v_com * dt_sub;

        let mut positions = Vec::with_capacity(n);
        let mut velocities = Vec::with_capacity(n);
        let mut accelerations = Vec::with_capacity(n);

        // Central body first (index 0 by system convention).
        positions.push(r0_rest + dr_com);
        velocities.push(v0_rest + self.v_com);

        // Sun's gravitational acceleration at the propagated config.
        // a_0 = G Σ m_i q_i(h) / |q_i(h)|³ — placeholder, filled after
        // the planet loop populates the inertial-frame helio vectors.
        let mut sun_acc = Vec3::ZERO;
        for (i, q) in q_helio.iter().enumerate() {
            let r2 = q.length_squared().max(1e-60);
            let inv_r3 = 1.0 / (r2 * r2.sqrt());
            let g = self.mu / self.m_sun; // G = μ / m_central
            sun_acc += g * self.planet_masses[i] * *q * inv_r3;

            positions.push(*q + r0_rest + dr_com);
            velocities.push(v_inertial_rest[i] + self.v_com);
            // Planet's inertial gravitational acceleration: Sun's pull
            // alone. Perturbation kicks are omitted (see method docs).
            let kepler_pull = -self.mu * *q * inv_r3;
            accelerations.push(kepler_pull);
        }
        // Sun's slot was reserved with a default; fill in.
        if n > 0 {
            accelerations.insert(0, sun_acc);
        }

        WhKinematics { positions, velocities, accelerations }
    }
}

impl DenseSnapshot {
    /// Interpolated world position for body `i` at normalised time `h ∈ [0, 1]`.
    ///
    /// Panics in debug mode if `i >= self.x0.len()`.
    #[inline]
    pub fn interpolate(&self, i: usize, h: f64) -> Vec3 {
        debug_assert!(i < self.x0.len(), "body index out of range");

        if let Some(wh) = &self.wh_data {
            return wh.interpolate_kinematics(h, self.dt).positions[i];
        }

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

    /// Interpolated velocity for body `i` at normalised time `h ∈ [0, 1]`.
    ///
    /// Returned alongside [`interpolate`](Self::interpolate) so a render
    /// frame's `(position, velocity)` pair lives at the same point in
    /// the integrator step. Sub-step consumers (camera follow's
    /// feedforward predictor, perturbation forces that read body.vel
    /// inside Picard iteration) need this consistency to avoid biasing
    /// `O(a · h · dt)` per evaluation.
    ///
    /// Panics in debug mode if `i >= self.x0.len()`.
    #[inline]
    pub fn velocity_at(&self, i: usize, h: f64) -> Vec3 {
        debug_assert!(i < self.x0.len(), "body index out of range");

        if let Some(wh) = &self.wh_data {
            return wh.interpolate_kinematics(h, self.dt).velocities[i];
        }

        let v0 = self.v0[i];
        let a0 = self.a0[i];
        let dt = self.dt;

        if !self.b.is_empty() {
            predict_v_ias15(v0, a0, &self.b[i], h, dt)
        } else {
            predict_v_order2(v0, a0, h, dt)
        }
    }

    /// Interpolated acceleration for body `i` at normalised time `h ∈ [0, 1]`.
    ///
    /// Companion to [`interpolate`](Self::interpolate) and
    /// [`velocity_at`](Self::velocity_at) — render consumers that
    /// combine the kinematic triple (camera follow's feedforward
    /// predictor; field queries that paint by `|a|`) need it sampled
    /// at the same point inside the step rather than pinned to the
    /// step boundary.
    ///
    /// For VV / Y4 / WH the order-2 Taylor model treats acceleration
    /// as constant within a step, so this returns `a0` directly. For
    /// IAS15 it evaluates the polynomial second derivative — the
    /// b-coefficients ARE the higher-order acceleration terms in
    /// Gauss–Radau form, so this is the lightest of the three IAS15
    /// kernels.
    ///
    /// Panics in debug mode if `i >= self.x0.len()`.
    #[inline]
    pub fn acceleration_at(&self, i: usize, h: f64) -> Vec3 {
        debug_assert!(i < self.x0.len(), "body index out of range");

        if let Some(wh) = &self.wh_data {
            return wh.interpolate_kinematics(h, self.dt).accelerations[i];
        }

        if self.b.is_empty() { self.a0[i] } else { predict_a_ias15(self.a0[i], &self.b[i], h) }
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

/// Analytical derivative of [`predict_order2`]: `v₀ + a₀·h·dt`.
///
/// Companion to [`predict_order2`] for VV, Yoshida-4, and Wisdom–Holman.
/// Render consumers that read interpolated position must read
/// interpolated velocity from the same `h` to keep their
/// `(position, velocity)` pair consistent at the same point inside the step.
#[inline]
pub fn predict_v_order2(v0: Vec3, a0: Vec3, h: f64, dt: f64) -> Vec3 {
    let s = h * dt;
    Vec3::new(v0.x + a0.x * s, v0.y + a0.y * s, v0.z + a0.z * s)
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
/// operator registered through
/// [`HamiltonianOperator::accumulate_force`](crate::physics::integrator::HamiltonianOperator::accumulate_force)
/// or
/// [`NonConservativeOperator::accumulate_force`](crate::physics::integrator::NonConservativeOperator::accumulate_force)
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

/// IAS15 degree-15 acceleration at substep fraction `h ∈ [0, 1]`.
///
/// Differentiating the velocity polynomial in [`predict_v_ias15`]
/// once more with respect to physical time `t = h · dt` collapses
/// the divisors and gives:
///
/// `a(h) = a₀ + b₀·h + b₁·h² + b₂·h³ + b₃·h⁴ + b₄·h⁵ + b₅·h⁶ + b₆·h⁷`
///
/// The Gauss–Radau b-coefficients are precisely the higher-order
/// acceleration terms, so `predict_a_ias15` is the lightest of the
/// three IAS15 kernels (no `dt` factor, no division).
#[inline]
pub fn predict_a_ias15(a0: Vec3, b: &DenseCoeffs, h: f64) -> Vec3 {
    let h2 = h * h;
    let h3 = h2 * h;
    let h4 = h3 * h;
    let h5 = h4 * h;
    let h6 = h5 * h;
    let h7 = h6 * h;

    let ax = a0.x
        + b[0].x * h
        + b[1].x * h2
        + b[2].x * h3
        + b[3].x * h4
        + b[4].x * h5
        + b[5].x * h6
        + b[6].x * h7;

    let ay = a0.y
        + b[0].y * h
        + b[1].y * h2
        + b[2].y * h3
        + b[3].y * h4
        + b[4].y * h5
        + b[5].y * h6
        + b[6].y * h7;

    let az = a0.z
        + b[0].z * h
        + b[1].z * h2
        + b[2].z * h3
        + b[3].z * h4
        + b[4].z * h5
        + b[5].z * h6
        + b[6].z * h7;

    Vec3::new(ax, ay, az)
}

#[cfg(test)]
mod tests {
    use super::{
        DenseCoeffs, DenseSnapshot, predict_a_ias15, predict_ias15, predict_order2,
        predict_v_ias15, predict_v_order2,
    };
    use crate::math::Vec3;
    use crate::physics::integrator::IntegratorKind;

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

    #[test]
    fn predict_v_order2_at_h_zero_returns_v0() {
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        assert_eq!(predict_v_order2(v0, a0, 0.0, 1e-3), v0);
    }

    #[test]
    fn predict_v_order2_is_derivative_of_predict_order2() {
        // Central difference of position polynomial against analytical
        // velocity. Same tolerance class as the IAS15 derivative test.
        let x0 = Vec3::new(0.5, 0.3, -0.2);
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let dt = 1e-3;
        let eps = 1e-5;
        for h in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let xp = predict_order2(x0, v0, a0, h + eps, dt);
            let xm = predict_order2(x0, v0, a0, h - eps, dt);
            let v_num = Vec3::new(
                (xp.x - xm.x) / (2.0 * eps * dt),
                (xp.y - xm.y) / (2.0 * eps * dt),
                (xp.z - xm.z) / (2.0 * eps * dt),
            );
            let v = predict_v_order2(v0, a0, h, dt);
            assert!((v.x - v_num.x).abs() < 1e-7, "vx at h={h}");
            assert!((v.y - v_num.y).abs() < 1e-7, "vy at h={h}");
            assert!((v.z - v_num.z).abs() < 1e-7, "vz at h={h}");
        }
    }

    #[test]
    fn snapshot_velocity_at_dispatches_to_order2_when_b_empty() {
        // VV / Y4 / WH leave `b` empty; `velocity_at` must fall back to
        // the Order-2 predictor (matching `interpolate`'s position fallback).
        let snap = DenseSnapshot {
            t0: 0.0,
            dt: 1e-3,
            x0: vec![Vec3::new(0.0, 0.0, 0.0)],
            v0: vec![Vec3::new(2.0, -1.0, 0.5)],
            a0: vec![Vec3::new(0.4, 0.0, -0.1)],
            b: Vec::new(),
            kind: IntegratorKind::Yoshida4,
            wh_data: None,
        };
        let v = snap.velocity_at(0, 0.5);
        let expected = predict_v_order2(snap.v0[0], snap.a0[0], 0.5, snap.dt);
        assert_eq!(v, expected);
    }

    #[test]
    fn snapshot_velocity_at_uses_ias15_when_b_present() {
        let snap = DenseSnapshot {
            t0: 0.0,
            dt: 1e-3,
            x0: vec![Vec3::new(0.0, 0.0, 0.0)],
            v0: vec![Vec3::new(2.0, -1.0, 0.5)],
            a0: vec![Vec3::new(0.4, 0.0, -0.1)],
            b: vec![sample_b()],
            kind: IntegratorKind::Ias15,
            wh_data: None,
        };
        let v = snap.velocity_at(0, 0.5);
        let expected = predict_v_ias15(snap.v0[0], snap.a0[0], &snap.b[0], 0.5, snap.dt);
        assert_eq!(v, expected);
    }

    // ── Acceleration interpolation ───────────────────────────────────────────

    #[test]
    fn predict_a_ias15_at_h_zero_returns_a0() {
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let b = sample_b();
        assert_eq!(predict_a_ias15(a0, &b, 0.0), a0);
    }

    #[test]
    fn predict_a_ias15_is_derivative_of_predict_v_ias15() {
        // Central difference of velocity polynomial against analytical
        // acceleration. Same tolerance class as the velocity-derivative
        // test in this module.
        let v0 = Vec3::new(1.5, -0.7, 0.4);
        let a0 = Vec3::new(0.3, 0.2, -0.1);
        let b = sample_b();
        let dt = 1e-3;
        let eps = 1e-5;
        for h in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let vp = predict_v_ias15(v0, a0, &b, h + eps, dt);
            let vm = predict_v_ias15(v0, a0, &b, h - eps, dt);
            let a_num = Vec3::new(
                (vp.x - vm.x) / (2.0 * eps * dt),
                (vp.y - vm.y) / (2.0 * eps * dt),
                (vp.z - vm.z) / (2.0 * eps * dt),
            );
            let a = predict_a_ias15(a0, &b, h);
            assert!((a.x - a_num.x).abs() < 1e-6, "ax at h={h}: {} vs {}", a.x, a_num.x);
            assert!((a.y - a_num.y).abs() < 1e-6, "ay at h={h}: {} vs {}", a.y, a_num.y);
            assert!((a.z - a_num.z).abs() < 1e-6, "az at h={h}: {} vs {}", a.z, a_num.z);
        }
    }

    #[test]
    fn snapshot_acceleration_at_returns_a0_when_b_empty() {
        // VV / Y4 / WH: order-2 model has constant acceleration within
        // a step.
        let snap = DenseSnapshot {
            t0: 0.0,
            dt: 1e-3,
            x0: vec![Vec3::new(0.0, 0.0, 0.0)],
            v0: vec![Vec3::new(2.0, -1.0, 0.5)],
            a0: vec![Vec3::new(0.4, 0.0, -0.1)],
            b: Vec::new(),
            kind: IntegratorKind::Yoshida4,
            wh_data: None,
        };
        for h in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_eq!(snap.acceleration_at(0, h), snap.a0[0]);
        }
    }

    #[test]
    fn snapshot_acceleration_at_uses_ias15_when_b_present() {
        let snap = DenseSnapshot {
            t0: 0.0,
            dt: 1e-3,
            x0: vec![Vec3::new(0.0, 0.0, 0.0)],
            v0: vec![Vec3::new(2.0, -1.0, 0.5)],
            a0: vec![Vec3::new(0.4, 0.0, -0.1)],
            b: vec![sample_b()],
            kind: IntegratorKind::Ias15,
            wh_data: None,
        };
        let a = snap.acceleration_at(0, 0.5);
        let expected = predict_a_ias15(snap.a0[0], &snap.b[0], 0.5);
        assert_eq!(a, expected);
    }
}
