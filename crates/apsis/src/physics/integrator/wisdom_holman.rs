//! Wisdom–Holman mixed-variable symplectic integrator (Wisdom & Holman 1991).
//!
//! Second-order kick-drift-kick symplectic split in heliocentric Cartesian
//! coordinates with inertial momenta. The Hamiltonian decomposes as
//!
//! ```text
//! H = H_K + H_I
//! H_K = Σ ( p_i² / 2 m_i  −  G m_0 m_i / |q_i| )      (Keplerian, analytical)
//! H_I = Σ_{i<j} ( −G m_i m_j / |q_i − q_j| )          (planet-planet)
//!     + ( Σ p_i )² / (2 m_0)                          (indirect)
//! ```
//!
//! where `q_i = r_i − r_0` is the heliocentric position of planet `i` and
//! `p_i = m_i v_i` is the inertial momentum (Wisdom & Holman 1991 §III).
//!
//! Per-step structure:
//!
//! 1. Snapshot the central body's inertial state and the system's total
//!    inertial momentum (conserved exactly through the step).
//! 2. Translate planet positions to the heliocentric frame.
//! 3. Half-kick on planet velocities from `H_I` (planet-planet potential
//!    gradient + indirect-acceleration term).
//! 4. Drift on planet states under `H_K`: each planet propagates
//!    analytically along its Keplerian orbit around a fixed-center
//!    potential at the origin via the universal-variable Kepler
//!    propagator in [`super::kepler::kepler_step`]. Planet inertial
//!    velocities are also updated by this step.
//! 5. Drift on planet positions under the indirect-momentum term of
//!    `H_I`: each `q_i` shifts by `(Σ m_j v_j) / m_0 · dt`, where the
//!    sum is over post-Kepler planet momenta. The two drift operators
//!    commute (one depends only on `q`, the other only on `p`); applying
//!    them in either order yields the same KDK truncation.
//! 6. Half-kick again.
//! 7. Translate planets back to inertial coordinates. The central body's
//!    new inertial position is recovered from the barycenter conservation
//!    constraint (`Q_0` = total-mass-weighted system position is invariant
//!    in the rest frame, and otherwise advances at `P_0 / M`); its new
//!    inertial velocity is recovered from total-momentum conservation.
//!
//! The integrator is permissive on the force-model pairing: any `ForceModel`
//! satisfies the requirements, since the planet-planet kick is computed
//! through the standard force-evaluation interface.
//!
//! # References
//! - Wisdom, J. & Holman, M. (1991). *Astron. J.* 102, 1528–1538.
//! - Duncan, M., Levison, H., & Lee, M. H. (1998). *Astron. J.* 116, 2067.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::dense::{DenseSnapshot, WhDenseData};
use crate::physics::integrator::helpers::{
    apply_perturbations_planets, evaluate, scale_acc_and_pe,
};
use crate::physics::integrator::kepler::kepler_step;
use crate::physics::integrator::traits::{
    HierarchySignal, Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Minimum ratio `M_central / Σ m_i (i > 0)` for which the WH split derivation
/// holds without a perturbation expansion that would dominate the integrator
/// truncation error.
const WH_DOMINANCE_RATIO: f64 = 10.0;

/// Wisdom–Holman mixed-variable symplectic map.
pub struct WisdomHolman;

impl Default for WisdomHolman {
    fn default() -> Self {
        Self::new()
    }
}

impl WisdomHolman {
    pub fn new() -> Self {
        Self
    }

    /// Returns `true` if `bodies[0]` dominates the system mass distribution
    /// to the threshold required for the Wisdom–Holman perturbation expansion.
    ///
    /// Two conditions must hold: the central body must be at least as massive
    /// as any other single body, and the central-to-rest mass ratio must be
    /// at least `WH_DOMINANCE_RATIO`.
    pub fn is_suitable_for(bodies: &[Body]) -> bool {
        if bodies.len() < 2 {
            return false;
        }
        let m0 = bodies[0].mass;
        let m_rest: f64 = bodies[1..].iter().map(|b| b.mass).sum();
        let max_other = bodies[1..].iter().map(|b| b.mass).fold(0.0_f64, f64::max);
        m0 >= max_other && m0 >= WH_DOMINANCE_RATIO * m_rest
    }
}

impl Integrator for WisdomHolman {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        if bodies.len() < 2 {
            return StepResult {
                consumed_dt: dt,
                potential_energy: 0.0,
                used_fallback: false,
                step_snapshot: None,
                degraded: false,
                hierarchy_signal: Some(HierarchySignal::Violated),
            };
        }

        let m0 = bodies[0].mass;
        let m_total: f64 = bodies.iter().map(|b| b.mass).sum();
        let mu = ctx.g_factor * m0;

        // Capture pre-step inertial state for the dense-output snapshot.
        // Done before any frame transformation so the snapshot's `x0`,
        // `v0`, `a0` fields carry the body-aligned, original-frame
        // kinematics every consumer of `DenseSnapshot` expects (the
        // wh_data path uses its own rest-frame state separately).
        let pre_x0_inertial: Vec<Vec3> = bodies.iter().map(|b| Vec3::new(b.x, b.y, b.z)).collect();
        let pre_v0_inertial: Vec<Vec3> =
            bodies.iter().map(|b| Vec3::new(b.vx, b.vy, b.vz)).collect();
        let pre_a0_inertial: Vec<Vec3> = if acc.len() == bodies.len() {
            acc.clone()
        } else {
            // First step (no scratch_acc populated yet) or post-resize
            // mismatch: ZERO is a safe placeholder — `wh_data` is the
            // primary interpolation path and does not read `a0`. The
            // fallback only triggers if `wh_data` itself is None.
            vec![Vec3::ZERO; bodies.len()]
        };

        // The Wisdom-Holman canonical formulation uses heliocentric positions
        // and inertial momenta in the rest frame. To extend to arbitrary
        // initial frames without altering the algorithm, this implementation
        // performs a Galilean transformation to the centre-of-mass rest frame
        // at step entry, runs the symplectic split there, and applies the
        // inverse Galilean transformation at step exit. Total momentum is
        // exactly conserved by the symplectic split, so `v_com` is unchanged
        // through the step.
        let p_total =
            bodies.iter().fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vx, b.vy, b.vz));
        let v_com = p_total / m_total;
        for b in bodies.iter_mut() {
            b.vx -= v_com.x;
            b.vy -= v_com.y;
            b.vz -= v_com.z;
        }

        // Step-entry snapshots in the rest frame.
        let r0_in = Vec3::new(bodies[0].x, bodies[0].y, bodies[0].z);

        // The barycenter-constraint reconstruction at step exit needs the
        // step-entry value of `Σ m_i q_i_in` (i ≥ 1) to derive `r_0_out` from
        // the post-step planet positions. In the rest frame `Q_0` is invariant.
        let m_q_in: Vec3 = bodies[1..]
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * (Vec3::new(b.x, b.y, b.z) - r0_in));

        // ── Translate planets to heliocentric (positions only) ────────────
        for b in bodies[1..].iter_mut() {
            b.x -= r0_in.x;
            b.y -= r0_in.y;
            b.z -= r0_in.z;
        }

        // Capture the rest-frame (q_helio, v_inertial) pair for every
        // non-central body. Done after the heliocentric translation
        // (so `q` is what `kepler_step` expects) and before the first
        // half-kick mutates velocities. This is the input the dense-
        // output Kepler interpolator replays at sub-step times.
        let q0_helio_rest: Vec<Vec3> =
            bodies[1..].iter().map(|b| Vec3::new(b.x, b.y, b.z)).collect();
        let v0_inertial_rest: Vec<Vec3> =
            bodies[1..].iter().map(|b| Vec3::new(b.vx, b.vy, b.vz)).collect();
        let planet_masses: Vec<f64> = bodies[1..].iter().map(|b| b.mass).collect();

        // ── First half-kick ───────────────────────────────────────────────
        let pe = wh_kick(bodies, ctx, 0.5 * dt, acc, mu);

        // ── Drift: H_K (analytical Kepler around fixed origin, per planet)
        // The Hamiltonian `H_K = p_i² / 2 m_i − G m_0 m_i / |q_i|` propagates
        // each planet around a fixed central potential of strength `mu = G m_0`
        // (Wisdom & Holman 1991 §III). The leading O(m_i / m_0) correction
        // relative to the true two-body problem is absorbed by the H_indirect
        // drift below.
        for b in bodies[1..].iter_mut() {
            let q = Vec3::new(b.x, b.y, b.z);
            let v = Vec3::new(b.vx, b.vy, b.vz);
            let (q_new, v_new) = kepler_step(q, v, dt, mu);
            b.x = q_new.x;
            b.y = q_new.y;
            b.z = q_new.z;
            b.vx = v_new.x;
            b.vy = v_new.y;
            b.vz = v_new.z;
        }

        // ── Drift: H_indirect (uniform shift on all heliocentric positions)
        // The indirect kinetic cross-term `(Σ p_i)² / (2 m_0)` depends only on
        // momenta and so generates a position drift; under H_K + H_indirect
        // applied sequentially the indirect drift uses the post-Kepler planet
        // momenta. The shift is identical for all planets.
        let p_planets_post_kepler: Vec3 =
            bodies[1..].iter().fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vx, b.vy, b.vz));
        let indirect_shift = (p_planets_post_kepler / m0) * dt;
        for b in bodies[1..].iter_mut() {
            b.x += indirect_shift.x;
            b.y += indirect_shift.y;
            b.z += indirect_shift.z;
        }

        // ── Second half-kick ──────────────────────────────────────────────
        let _ = wh_kick(bodies, ctx, 0.5 * dt, acc, mu);

        // ── Reconstruct central body inertial state in the rest frame ─────
        // Barycenter constraint in the rest frame: `Q_0` is invariant, so
        //   r_0_out = r_0_in + (m_q_in − m_q_out) / M
        // where m_q_out is `Σ m_i q_i_post` evaluated on the post-step
        // heliocentric positions.
        let m_q_out: Vec3 =
            bodies[1..].iter().fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.x, b.y, b.z));
        let r0_out = r0_in + (m_q_in - m_q_out) / m_total;

        // ── Translate planets back to rest-frame inertial coordinates ────
        for b in bodies[1..].iter_mut() {
            b.x += r0_out.x;
            b.y += r0_out.y;
            b.z += r0_out.z;
        }

        // Rest-frame total-momentum conservation: `Σ_all m_i v_i = 0`, so
        //   v_0_out = −(1/m_0) Σ_{i≥1} m_i v_i_out.
        let p_planets_out: Vec3 =
            bodies[1..].iter().fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vx, b.vy, b.vz));
        let v0_out_rest = -p_planets_out / m0;

        bodies[0].x = r0_out.x;
        bodies[0].y = r0_out.y;
        bodies[0].z = r0_out.z;
        bodies[0].vx = v0_out_rest.x;
        bodies[0].vy = v0_out_rest.y;
        bodies[0].vz = v0_out_rest.z;

        // ── Inverse Galilean transformation back to the original frame ────
        // Advance all bodies by `v_com · dt` (centre-of-mass position drift
        // over the step) and restore each body's `v_com` velocity component.
        let dr_com = v_com * dt;
        for b in bodies.iter_mut() {
            b.x += dr_com.x;
            b.y += dr_com.y;
            b.z += dr_com.z;
            b.vx += v_com.x;
            b.vy += v_com.y;
            b.vz += v_com.z;
        }

        // Populate `acc` with the total inertial acceleration of every
        // body (size N) so consumers downstream of the integrator see a
        // body-aligned acceleration vector. Three of them depend on this
        // contract:
        //
        //   * The dense-output snapshot path in `System::step` requires
        //     `scratch_acc.len() == bodies.len()`; without it, no Order-2
        //     snapshot is built for WH and the renderer cannot
        //     interpolate body positions within a step (visible as
        //     bodies "freezing" between sparse publish ticks at slow
        //     achieved sim rates).
        //   * `RenderState::accelerations` is indexed by body in the UI;
        //     a size-(N-1) buffer published from WH would silently
        //     off-by-one every consumer (camera feedforward, |a| field).
        //   * `Metrics::last_accelerations` exposes the buffer to
        //     external readers (CSV export, hooks); a misaligned shape
        //     is a contract leak.
        //
        // At this point `acc` holds planet-planet + perturbation
        // accelerations on planets (size N-1) from the second half-kick.
        // The Sun's gravitational pull on each planet was applied
        // analytically by `kepler_step` and is not in the buffer; the
        // Sun's own inertial acceleration was never computed (its state
        // was reconstructed from barycenter conservation, not
        // integrated). The total inertial acceleration of body `i` is:
        //
        //   planet i (i ≥ 1):  a_i = acc[i−1]  −  μ q_i / |q_i|³
        //   central body 0:    a_0 = G Σ_{i≥1} m_i q_i / |q_i|³
        //
        // where `q_i = bodies[i].pos − bodies[0].pos`. Both expressions
        // are evaluated in the original inertial frame the bodies are
        // now in (post Galilean inverse shift above).
        let n = bodies.len();
        let r0_inertial = Vec3::new(bodies[0].x, bodies[0].y, bodies[0].z);
        let mut planet_total_acc = Vec::with_capacity(n - 1);
        let mut sun_acc = Vec3::ZERO;
        for (i, b) in bodies[1..].iter().enumerate() {
            let q = Vec3::new(b.x, b.y, b.z) - r0_inertial;
            let r2 = q.length_squared().max(1e-60);
            let inv_r3 = 1.0 / (r2 * r2.sqrt());
            let kepler_pull_on_planet = -mu * q * inv_r3;
            planet_total_acc.push(acc[i] + kepler_pull_on_planet);
            sun_acc += ctx.g_factor * b.mass * q * inv_r3;
        }
        acc.resize(n, Vec3::ZERO);
        acc[0] = sun_acc;
        for (i, &a) in planet_total_acc.iter().enumerate() {
            acc[i + 1] = a;
        }

        // Classify the system's current mass distribution against the WH
        // dominance criterion and surface it through `hierarchy_signal`.
        // Observability only — the integrator has already run; downstream
        // consumers (Metrics, UI, logging) read the signal to detect when
        // the system has left the validated regime.
        let masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
        let signal = HierarchySignal::classify(&masses);

        // Build the dense-output snapshot. The Kepler-analytical kernel
        // in `WhDenseData::interpolate_kinematics` replays each
        // non-central body's drift at sub-step times via the same
        // universal-variable solver the integrator's drift step used,
        // so the renderer sees curved orbital trajectories within a
        // step instead of order-2 Taylor straight-line approximations.
        // The standard `x0 / v0 / a0` fields are populated as a
        // defensive fallback path if `wh_data` ever fails to dispatch.
        let wh_data = WhDenseData {
            mu,
            m_sun: m0,
            m_total,
            v_com,
            r0_sun_rest: r0_in,
            m_q_in,
            q0_helio_rest,
            v0_inertial_rest,
            planet_masses,
        };
        let step_snapshot = DenseSnapshot {
            t0: 0.0,
            dt,
            x0: pre_x0_inertial,
            v0: pre_v0_inertial,
            a0: pre_a0_inertial,
            b: Vec::new(),
            kind: IntegratorKind::WisdomHolman,
            wh_data: Some(wh_data),
        };

        StepResult {
            consumed_dt: dt,
            potential_energy: pe,
            used_fallback: false,
            step_snapshot: Some(step_snapshot),
            degraded: false,
            hierarchy_signal: Some(signal),
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::WisdomHolman
    }
}

/// Apply a perturbation kick of duration `dt` to all planet velocities.
///
/// The kick comprises three contributions, all in inertial-velocity units:
/// the planet-planet pair forces, the indirect-acceleration response
/// `−Σ m_j a_j / m_0` accounting for the non-inertial heliocentric frame, and
/// any registered non-gravitational perturbation forces. Returns the
/// gravitational potential energy (planet-planet plus planet-central) for
/// the post-kick configuration, scaled by `g_factor`.
fn wh_kick(
    bodies: &mut [Body],
    ctx: &mut IntegratorContext<'_>,
    dt: f64,
    acc: &mut Vec<Vec3>,
    mu: f64,
) -> f64 {
    let m0 = bodies[0].mass;

    let raw_pe = evaluate(&bodies[1..], ctx.force, acc);

    let bary_acc_raw =
        acc.iter().zip(bodies[1..].iter()).fold(Vec3::ZERO, |a, (a_i, b)| a + b.mass * *a_i);
    let indirect_raw = -bary_acc_raw / m0;

    let pe_inter = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
    let pe_central = central_potential(bodies, mu);

    apply_perturbations_planets(&bodies[1..], acc, ctx.perturbations);

    let indirect = indirect_raw * ctx.g_factor;
    for (i, ai) in acc.iter().enumerate() {
        let kick = (*ai + indirect) * dt;
        bodies[i + 1].vx += kick.x;
        bodies[i + 1].vy += kick.y;
        bodies[i + 1].vz += kick.z;
    }

    pe_inter + pe_central
}

/// Central Keplerian potential `−μ Σ m_i / |q_i|` evaluated in the
/// heliocentric frame (planet positions relative to central body).
fn central_potential(bodies: &[Body], mu: f64) -> f64 {
    bodies[1..]
        .iter()
        .map(|b| {
            let r = (b.x * b.x + b.y * b.y + b.z * b.z).sqrt().max(1e-30);
            -mu * b.mass / r
        })
        .sum()
}
