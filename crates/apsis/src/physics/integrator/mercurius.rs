//! Mercurius — close-encounter hybrid symplectic integrator.
//!
//! Implementation of the algorithm specified in:
//!
//!   * Rein H., Hernandez D. M., Tamayo D. & Brown G. (2019). *Hybrid
//!     symplectic integrators for planetary dynamics*, MNRAS **489**,
//!     4632–4640. [arXiv:1908.03468](https://arxiv.org/abs/1908.03468)
//!
//! Faithful port of REBOUND's `integrator_mercurius.c` to apsis. The
//! algorithm is a **rewind hybrid**: a Wisdom-Holman outer step
//! advances every particle, an encounter detector identifies pairs
//! that came within their critical radius during the step, those
//! pairs rewind to their pre-Kepler state, and IAS15 re-integrates
//! them over the same outer window with the full Sun pull plus the
//! (1−K)-weighted planet-planet residual. Non-encountering particles
//! keep their analytical-Kepler positions.
//!
//! The "Hamiltonian split with IAS15 integrating only the (1−K)
//! residual" is mathematically incomplete: a pure (1−K)·V Hamiltonian
//! generates a kick (no position evolution), so feeding it to IAS15
//! (which integrates ẍ = a) introduces a v·τ free drift that
//! double-counts Kepler. The rewind structure avoids this by giving
//! IAS15 the *full* dynamical equations on the encountering subset.
//!
//! # Per-step structure (democratic-heliocentric coordinates)
//!
//! 1. `interaction(τ/2)` — K-weighted half-kick on planet velocities.
//! 2. `jump(τ/2)` — uniform position drift on planets, accounting for
//!    Sun's recoil from planet momenta.
//! 3. `com(τ)` — advance the COM by `τ · v_com`.
//! 4. `backup` — snapshot every particle's pre-Kepler state.
//! 5. `kepler(τ)` — analytical Kepler drift around the Sun, every planet.
//! 6. `encounter_predict` — for every pair, fit `r²(t)` over the
//!    interval via cubic-Hermite interpolation between pre- and
//!    post-Kepler `(r, dr/dt)` and find `r_min`. Flag pairs with
//!    `r_min < dcrit_ij`.
//! 7. `encounter_step(τ)` — rewind every flagged particle to its
//!    pre-Kepler state. IAS15 integrates the rewound system over `τ`
//!    with the close-field force model (full Sun pull on planets +
//!    (1−K)-weighted planet-planet for encountering pairs). After
//!    IAS15, restore non-encountering planets to their post-Kepler
//!    state (encountering ones keep the IAS15 result).
//! 8. `jump(τ/2)` — repeat stage 2.
//! 9. `interaction(τ/2)` — repeat stage 1.
//!
//! # Changeover function `L_mercury`
//!
//! REBOUND's default `L_mercury(d, dcrit)`:
//!
//! ```text
//!     y = (d - 0.1·dcrit) / (0.9·dcrit)
//!     L = 0                                  for y ≤ 0   (deep encounter)
//!     L = 10·y³ - 15·y⁴ + 6·y⁵               for 0 < y < 1
//!     L = 1                                  for y ≥ 1   (no encounter)
//! ```
//!
//! C² quintic Hermite polynomial. The 0.1·dcrit deadband ensures
//! `L ≡ 0` deep inside the encounter so IAS15 carries the full
//! responsibility without leakage from the K-kick.
//!
//! # Critical radius `dcrit` (REBOUND `dcrit_for_particle`)
//!
//! Per-particle, max of four criteria:
//!
//! ```text
//!     dcrit_i = max( v_c · 0.4·τ ,                  # average velocity
//!                    |v_i| · 0.4·τ ,                # current velocity
//!                    α · |a_i| · ((m_i)/(3 m_0))^(1/3) ,  # mutual Hill
//!                    2 · r_i^physical )                   # physical radius
//! ```
//!
//! with `a_i` the osculating semi-major axis,
//! `v_c = √(G m_0 / |a_i|)` the circular velocity at that distance.
//! Pair-wise reduction: `dcrit_ij = max(dcrit_i, dcrit_j)`.
//!
//! # Hierarchical-system requirement
//!
//! Mercurius assumes a dominant central body. On non-hierarchical
//! configurations the analytical Kepler drift around `bodies[0]` is
//! ill-posed; `step` returns `used_fallback = true` and emits a
//! `warn_diag!` event without advancing time.
//!
//! # References
//!
//! - Rein et al. (2019). MNRAS 489, 4632–4640.
//! - Wisdom & Holman (1991). AJ 102, 1528.
//! - Rein & Spiegel (2015). MNRAS 446, 1424.
//! - Lab notebook: `docs/experiments/2026-05-13-mercurius-hybrid.md`.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::ias15::Ias15;
use crate::physics::integrator::kepler::kepler_step;
use crate::physics::integrator::traits::{
    HierarchySignal, Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Default Hill-radius multiplier (REBOUND's `r_crit_hill`).
///
/// Matches REBOUND's MERCURIUS default and the value validated by
/// Rein et al. (2019) §3 across several planetary scattering scenarios.
pub const DEFAULT_HILL_FACTOR: f64 = 3.0;

/// Tiny floor on `dcrit` to prevent `0/0` in the changeover.
const DCRIT_FLOOR: f64 = 1.0e-30;

/// Mercurius hybrid integrator — see module documentation for the
/// algorithm.
pub struct Mercurius {
    /// Hill-radius multiplier α (REBOUND's `r_crit_hill`).
    alpha: f64,

    /// Internal IAS15 sub-integrator used inside the encounter step.
    ias15: Ias15,

    /// Pre-Kepler state of every particle, captured each outer step
    /// before the analytical Kepler drift. Read by `encounter_predict`
    /// and `encounter_step`.
    particles_backup: Vec<Body>,

    /// Per-particle critical radius for this outer step. Recomputed
    /// inside `step` after the inertial→DH conversion.
    dcrit: Vec<f64>,

    /// Encounter map: `[i] = true` if particle `i` participated in any
    /// close-encounter pair during the outer step.
    encounter_map: Vec<bool>,

    /// COM position in inertial coords, updated each step.
    com_pos: Vec3,
    /// COM velocity in inertial coords; constant across the outer step
    /// (the algorithm is momentum-conserving by construction).
    com_vel: Vec3,

    /// Working buffer for the K-weighted interaction acceleration.
    acc_int: Vec<Vec3>,

    /// Working buffer for the close-field IAS15 acceleration.
    acc_close: Vec<Vec3>,

    /// Bodies count from the most recent outer step. Used to resize
    /// `particles_backup`, `dcrit`, `encounter_map` lazily.
    last_n: usize,
}

impl Default for Mercurius {
    fn default() -> Self {
        Self::new()
    }
}

impl Mercurius {
    /// Create a Mercurius integrator with the canonical Rein et al.
    /// (2019) defaults: `α = 3`, embedded IAS15 with default tolerance.
    pub fn new() -> Self {
        Self::with_alpha(DEFAULT_HILL_FACTOR)
    }

    /// Create a Mercurius integrator with a custom Hill-radius
    /// multiplier `α`.
    ///
    /// Larger `α` enlarges the changeover band (more pairs spend time
    /// in the smooth-K transition; encounter step engages more often).
    /// Smaller `α` pushes the K → 0 boundary closer to the pair,
    /// localising the encounter-step engagement at the cost of risking
    /// missed encounters when a fast pair zips through the band in
    /// less than one outer step.
    pub fn with_alpha(alpha: f64) -> Self {
        Self {
            alpha: alpha.max(0.0),
            ias15: Ias15::new(),
            particles_backup: Vec::new(),
            dcrit: Vec::new(),
            encounter_map: Vec::new(),
            com_pos: Vec3::ZERO,
            com_vel: Vec3::ZERO,
            acc_int: Vec::new(),
            acc_close: Vec::new(),
            last_n: 0,
        }
    }

    /// Active Hill-radius multiplier α.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Override the Hill-radius multiplier at runtime.
    pub fn set_alpha(&mut self, alpha: f64) {
        self.alpha = alpha.max(0.0);
    }

    // ── Changeover function (REBOUND `L_mercury`) ──────────────────────

    /// REBOUND's `L_mercury` C² quintic Hermite changeover with a
    /// `0.1·dcrit` deadband.
    #[inline]
    fn changeover_l(d: f64, dcrit: f64) -> f64 {
        let y = (d - 0.1 * dcrit) / (0.9 * dcrit);
        if y <= 0.0 {
            0.0
        } else if y >= 1.0 {
            1.0
        } else {
            10.0 * y.powi(3) - 15.0 * y.powi(4) + 6.0 * y.powi(5)
        }
    }

    // ── Coordinate transformations (REBOUND `inertial_to_dh` / `dh_to_inertial`) ──

    /// Convert inertial coordinates to democratic-heliocentric: planet
    /// positions become Sun-relative; every body's velocity becomes
    /// COM-relative. The Sun's position becomes zero in DH and its
    /// velocity becomes `-v_com`.
    fn inertial_to_dh(&mut self, bodies: &mut [Body]) {
        let n = bodies.len();
        let mut com_pos = Vec3::ZERO;
        let mut com_vel = Vec3::ZERO;
        let mut mtot = 0.0_f64;
        for b in bodies.iter() {
            let m = b.mass;
            com_pos += m * Vec3::new(b.pos_x, b.pos_y, b.pos_z);
            com_vel += m * Vec3::new(b.vel_x, b.vel_y, b.vel_z);
            mtot += m;
        }
        com_pos = com_pos / mtot;
        com_vel = com_vel / mtot;

        let r0 = Vec3::new(bodies[0].pos_x, bodies[0].pos_y, bodies[0].pos_z);
        // Reverse iteration so bodies[0] is shifted last (its position is
        // the reference for everyone else's translation).
        for i in (0..n).rev() {
            bodies[i].pos_x -= r0.x;
            bodies[i].pos_y -= r0.y;
            bodies[i].pos_z -= r0.z;
            bodies[i].vel_x -= com_vel.x;
            bodies[i].vel_y -= com_vel.y;
            bodies[i].vel_z -= com_vel.z;
        }

        self.com_pos = com_pos;
        self.com_vel = com_vel;
    }

    /// Inverse of `inertial_to_dh`. Reconstructs the central body's
    /// inertial position from the COM constraint and restores every
    /// other body's COM-relative velocity to inertial velocity.
    fn dh_to_inertial(&self, bodies: &mut [Body]) {
        let n = bodies.len();
        let mut tmp = Vec3::ZERO;
        let mut tmp_v = Vec3::ZERO;
        let mut tmp_m = 0.0_f64;
        for b in bodies.iter().skip(1) {
            let m = b.mass;
            tmp += m * Vec3::new(b.pos_x, b.pos_y, b.pos_z);
            tmp_v += m * Vec3::new(b.vel_x, b.vel_y, b.vel_z);
            tmp_m += m;
        }
        let m_total = tmp_m + bodies[0].mass;
        tmp = tmp / m_total;
        tmp_v = tmp_v / bodies[0].mass;

        // Sun's inertial position: COM_pos − Σ_{i≥1} m_i q_i / M_total
        let r0 = self.com_pos - tmp;
        bodies[0].pos_x = r0.x;
        bodies[0].pos_y = r0.y;
        bodies[0].pos_z = r0.z;

        for b in bodies.iter_mut().skip(1) {
            b.pos_x += r0.x;
            b.pos_y += r0.y;
            b.pos_z += r0.z;
            b.vel_x += self.com_vel.x;
            b.vel_y += self.com_vel.y;
            b.vel_z += self.com_vel.z;
        }

        // Sun's inertial velocity from total-momentum conservation:
        //   v_0 = v_com − Σ_{i≥1} m_i v_i / m_0
        bodies[0].vel_x = self.com_vel.x - tmp_v.x;
        bodies[0].vel_y = self.com_vel.y - tmp_v.y;
        bodies[0].vel_z = self.com_vel.z - tmp_v.z;
        let _ = n;
    }

    // ── dcrit computation (REBOUND `dcrit_for_particle`) ──────────────

    /// Per-particle critical radius from REBOUND's four criteria.
    /// Operates on a planet in DH coords (`q_i` heliocentric, `v_i`
    /// COM-relative; we approximate `v_i − v_0 ≈ v_i` because Sun's
    /// COM-relative velocity is small).
    fn compute_dcrit_for(&self, planet: &Body, m_central: f64, dt: f64) -> f64 {
        let q = Vec3::new(planet.pos_x, planet.pos_y, planet.pos_z);
        let v = Vec3::new(planet.vel_x, planet.vel_y, planet.vel_z);
        let r = q.length();
        let v2 = v.length_squared();

        let g_m = m_central + planet.mass;
        // Osculating semi-major axis from the vis-viva relation:
        //   v² = GM (2/r − 1/a)  ⇒  a = GM·r / (2GM − r·v²)
        let denom = 2.0 * g_m - r * v2;
        let a = if denom.abs() > 1.0e-30 { g_m * r / denom } else { r };
        let vc = if a.abs() > 1.0e-30 { (g_m / a.abs()).sqrt() } else { 0.0 };

        let mut dcrit = 0.0_f64;
        // 1. Average velocity criterion.
        dcrit = dcrit.max(vc * 0.4 * dt.abs());
        // 2. Current velocity criterion.
        dcrit = dcrit.max(v2.sqrt() * 0.4 * dt.abs());
        // 3. Hill-radius criterion (α-multiplied).
        if m_central > 0.0 {
            let mass_ratio = libm::cbrt((planet.mass / (3.0 * m_central)).max(0.0));
            dcrit = dcrit.max(self.alpha * a.abs() * mass_ratio);
        }
        // 4. Physical radius criterion.
        dcrit = dcrit.max(2.0 * planet.physical_radius);

        dcrit.max(DCRIT_FLOOR)
    }

    /// Recompute `self.dcrit` for the current particles. Sun's `dcrit`
    /// uses the physical-radius criterion only.
    fn rebuild_dcrit(&mut self, bodies: &[Body], dt: f64) {
        let n = bodies.len();
        self.dcrit.resize(n, 0.0);
        if n == 0 {
            return;
        }
        self.dcrit[0] = (2.0 * bodies[0].physical_radius).max(DCRIT_FLOOR);
        let m0 = bodies[0].mass;
        for i in 1..n {
            self.dcrit[i] = self.compute_dcrit_for(&bodies[i], m0, dt);
        }
    }

    // ── Operators (REBOUND `interaction_step`, `jump_step`, `com_step`, `kepler_step`) ──

    /// K-weighted planet-planet acceleration on every planet.
    /// In DH coords the Sun's pull is excluded — handled analytically
    /// by the Kepler drift.
    fn evaluate_interaction(&mut self, bodies: &[Body], g_factor: f64) {
        let n = bodies.len();
        self.acc_int.clear();
        self.acc_int.resize(n, Vec3::ZERO);
        for i in 1..n {
            let qi = Vec3::new(bodies[i].pos_x, bodies[i].pos_y, bodies[i].pos_z);
            for j in 1..n {
                if i == j {
                    continue;
                }
                let qj = Vec3::new(bodies[j].pos_x, bodies[j].pos_y, bodies[j].pos_z);
                let dq = qj - qi;
                let r2 = dq.length_squared().max(DCRIT_FLOOR);
                let r = r2.sqrt();
                let dcrit_pair = self.dcrit[i].max(self.dcrit[j]);
                let l = Self::changeover_l(r, dcrit_pair);
                if l <= 0.0 {
                    continue;
                }
                let inv_r3 = 1.0 / (r * r2);
                self.acc_int[i] += dq * (l * g_factor * bodies[j].mass * inv_r3);
            }
        }
    }

    /// K-weighted half-kick + perturbation accumulation. The encounter
    /// step's `CloseFieldForceModel` deliberately skips perturbations
    /// to avoid double-counting the contribution already folded in
    /// here.
    fn interaction_step(
        &mut self,
        bodies: &mut [Body],
        g_factor: f64,
        dt: f64,
        hamiltonian: &[Box<dyn crate::physics::integrator::HamiltonianOperator>],
        non_conservative: &[Box<dyn crate::physics::integrator::NonConservativeOperator>],
    ) {
        self.evaluate_interaction(bodies, g_factor);
        for op in hamiltonian {
            op.accumulate_force(bodies, &mut self.acc_int);
        }
        for op in non_conservative {
            op.accumulate_force(bodies, &mut self.acc_int);
        }
        for (i, b) in bodies.iter_mut().enumerate().skip(1) {
            let kick = self.acc_int[i] * dt;
            b.vel_x += kick.x;
            b.vel_y += kick.y;
            b.vel_z += kick.z;
        }
    }

    /// `jump_step(dt)` — uniform position drift on every planet by
    /// `dt · (Σ_j m_j v_j) / m_0`. Captures the Sun's recoil from
    /// planet momenta.
    fn jump_step(&mut self, bodies: &mut [Body], dt: f64) {
        let m0 = bodies[0].mass;
        if m0 <= 0.0 {
            return;
        }
        let mut p = Vec3::ZERO;
        for b in bodies.iter().skip(1) {
            p += b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z);
        }
        let shift = p * (dt / m0);
        for b in bodies.iter_mut().skip(1) {
            b.pos_x += shift.x;
            b.pos_y += shift.y;
            b.pos_z += shift.z;
        }
    }

    /// `com_step(dt)` — advance the inertial COM by `dt · v_com`.
    fn com_step(&mut self, dt: f64) {
        self.com_pos += self.com_vel * dt;
    }

    /// `kepler_step(dt)` — analytical Kepler drift around the Sun for
    /// every planet.
    fn kepler_step_all(&self, bodies: &mut [Body], mu: f64, dt: f64) {
        for b in bodies.iter_mut().skip(1) {
            let q = Vec3::new(b.pos_x, b.pos_y, b.pos_z);
            let v = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
            let (q_new, v_new) = kepler_step(q, v, dt, mu);
            b.pos_x = q_new.x;
            b.pos_y = q_new.y;
            b.pos_z = q_new.z;
            b.vel_x = v_new.x;
            b.vel_y = v_new.y;
            b.vel_z = v_new.z;
        }
    }

    // ── Encounter prediction (REBOUND `encounter_predict`) ─────────────

    /// Find the minimum pairwise separation over the step using cubic-
    /// Hermite interpolation between pre- and post-Kepler `(r², dr²/dt)`.
    /// Returns `r²_min` (squared, to avoid an extra sqrt in the caller).
    #[inline]
    fn rmin_squared_over_step(ro: f64, rn: f64, drodt: f64, drndt: f64, dt: f64) -> f64 {
        // REBOUND-faithful Hermite minimisation. Find the t ∈ (0, 1)
        // minimum of the cubic Hermite of r²(t) using its quadratic
        // derivative; sample the cubic at any minimum that lies inside
        // the unit interval and reduce r²_min accordingly.
        let a = 6.0 * (ro - rn) + 3.0 * dt * (drodt + drndt);
        let b = 6.0 * (rn - ro) - 2.0 * dt * (2.0 * drodt + drndt);
        let c = dt * drodt;
        let s = b * b - 4.0 * a * c;
        let sr = s.max(0.0).sqrt();
        let mut rmin2 = ro.min(rn);
        if a.abs() > 1.0e-30 {
            let two_a = 2.0 * a;
            for tmin in [(-b + sr) / two_a, (-b - sr) / two_a] {
                if tmin > 0.0 && tmin < 1.0 {
                    let omt = 1.0 - tmin;
                    let rmin_t = omt * omt * (1.0 + 2.0 * tmin) * ro
                        + tmin * tmin * (3.0 - 2.0 * tmin) * rn
                        + tmin * omt * omt * dt * drodt
                        - tmin * tmin * omt * dt * drndt;
                    rmin2 = rmin2.min(rmin_t.max(0.0));
                }
            }
        }
        rmin2
    }

    /// Populate `self.encounter_map` from current vs backup particle
    /// positions / velocities.
    fn encounter_predict(&mut self, bodies: &[Body], dt: f64) {
        let n = bodies.len();
        self.encounter_map.clear();
        self.encounter_map.resize(n, false);

        for i in 0..n {
            for j in (i + 1)..n {
                let (xn, yn, zn, vxn, vyn, vzn) = (
                    bodies[i].pos_x - bodies[j].pos_x,
                    bodies[i].pos_y - bodies[j].pos_y,
                    bodies[i].pos_z - bodies[j].pos_z,
                    bodies[i].vel_x - bodies[j].vel_x,
                    bodies[i].vel_y - bodies[j].vel_y,
                    bodies[i].vel_z - bodies[j].vel_z,
                );
                let bi = &self.particles_backup[i];
                let bj = &self.particles_backup[j];
                let (xo, yo, zo, vxo, vyo, vzo) = (
                    bi.pos_x - bj.pos_x,
                    bi.pos_y - bj.pos_y,
                    bi.pos_z - bj.pos_z,
                    bi.vel_x - bj.vel_x,
                    bi.vel_y - bj.vel_y,
                    bi.vel_z - bj.vel_z,
                );
                let rn = xn * xn + yn * yn + zn * zn;
                let ro = xo * xo + yo * yo + zo * zo;
                let drndt = 2.0 * (xn * vxn + yn * vyn + zn * vzn);
                let drodt = 2.0 * (xo * vxo + yo * vyo + zo * vzo);

                let rmin2 = Self::rmin_squared_over_step(ro, rn, drodt, drndt, dt);
                let dcrit_pair = self.dcrit[i].max(self.dcrit[j]);
                // REBOUND's 1.21× safety factor (≈ (1.1)²) — be slightly
                // conservative on the trigger to absorb interpolation error.
                let trigger2 = 1.21 * dcrit_pair * dcrit_pair;
                if rmin2 < trigger2 {
                    self.encounter_map[i] = true;
                    self.encounter_map[j] = true;
                }
            }
        }
    }

    // ── Encounter step (REBOUND `encounter_step`) ─────────────────────

    /// Returns `true` iff at least one planet (index ≥ 1) is in the
    /// encounter map. The Sun (index 0) is set to `true` unconditionally
    /// in REBOUND's encounter map, so the count of "real" encounters is
    /// the number of `true` entries beyond the Sun.
    fn any_planet_in_encounter(&self) -> bool {
        self.encounter_map.iter().skip(1).any(|&b| b)
    }

    /// Re-integrate the encountering subset over `[0, dt]` with IAS15.
    /// Encountering particles are rewound to their pre-Kepler state
    /// before IAS15 starts; non-encountering particles keep their
    /// post-Kepler positions but are restored at the end of the
    /// session (IAS15's free-drift on them is discarded).
    fn encounter_step(&mut self, bodies: &mut [Body], g_factor: f64, dt: f64) {
        if !self.any_planet_in_encounter() {
            return;
        }

        let n = bodies.len();
        // Snapshot the full post-Kepler state so we can restore non-
        // encountering particles after IAS15.
        let post_kepler: Vec<Body> = bodies.to_vec();

        // Rewind encountering particles to their pre-Kepler state.
        for i in 0..n {
            if self.encounter_map[i] {
                bodies[i] = self.particles_backup[i];
            }
        }
        // The Sun is always rewound (and pinned) — it does not move
        // in DH coords during the encounter sub-integration.
        bodies[0] = self.particles_backup[0];
        bodies[0].vel_x = 0.0;
        bodies[0].vel_y = 0.0;
        bodies[0].vel_z = 0.0;

        // Sub-context with empty operator slices — the outer
        // interaction half-kicks have already folded the perturbation
        // contribution; passing it through here would double-count.
        let hamiltonian_empty: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> =
            Vec::new();
        let non_conservative_empty: Vec<
            Box<dyn crate::physics::integrator::NonConservativeOperator>,
        > = Vec::new();
        let mut observers_empty: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        let mut close_force =
            CloseFieldForceModel { dcrit: &self.dcrit, encounter_map: &self.encounter_map };
        let mut ctx_close = IntegratorContext {
            force: &mut close_force,
            g_factor,
            hamiltonian_perturbations: &hamiltonian_empty,
            non_conservative_perturbations: &non_conservative_empty,
            observers: &mut observers_empty,
        };

        // Drive IAS15 over the outer window [0, dt]. The controller's
        // dt_next is clamped per-call so it cannot overshoot the
        // boundary; FSAL is invalidated at the start because the
        // surrounding outer-step operators have moved bodies since the
        // last IAS15 invocation.
        self.ias15.invalidate_force_cache();
        let mut consumed = 0.0_f64;
        let exit_tol = dt.abs() * 1.0e-12;
        while consumed + exit_tol < dt {
            let remaining = dt - consumed;
            // Pin the Sun's velocity to zero before each substep so
            // IAS15 cannot free-drift it.
            bodies[0].vel_x = 0.0;
            bodies[0].vel_y = 0.0;
            bodies[0].vel_z = 0.0;
            self.ias15.cap_proposed_dt(remaining);
            let result = self.ias15.step(bodies, &mut ctx_close, remaining, &mut self.acc_close);
            consumed += result.consumed_dt;
            if result.consumed_dt <= 0.0 {
                break;
            }
        }

        // Restore non-encountering planets to their post-Kepler state.
        // The Sun also returns to its pre-encounter DH state (origin,
        // zero velocity). Encountering particles keep IAS15's result.
        bodies[0] = self.particles_backup[0];
        bodies[0].vel_x = 0.0;
        bodies[0].vel_y = 0.0;
        bodies[0].vel_z = 0.0;
        for i in 1..n {
            if !self.encounter_map[i] {
                bodies[i] = post_kepler[i];
            }
        }

        // FSAL stays invalid for the next outer step's encounter
        // session: the surrounding interaction/jump kicks will move
        // bodies before IAS15 sees them again.
        self.ias15.invalidate_force_cache();
    }
}

impl Integrator for Mercurius {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        // ── Validate hierarchy ─────────────────────────────────────────
        let masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
        let signal = HierarchySignal::classify(&masses);
        if !matches!(signal, HierarchySignal::Hierarchical | HierarchySignal::Borderline) {
            crate::warn_diag!(
                crate::core::log::Source::Integrator,
                "Mercurius selected on non-hierarchical configuration; step refused",
                regime = signal.label(),
                hint = "Mercurius requires a dominant central body; switch to ias15 directly",
            );
            return StepResult {
                consumed_dt: 0.0,
                potential_energy: 0.0,
                used_fallback: true,
                step_snapshot: None,
                degraded: true,
                hierarchy_signal: Some(signal),
            };
        }

        let n = bodies.len();
        if n < 2 {
            return StepResult {
                consumed_dt: dt,
                potential_energy: 0.0,
                used_fallback: false,
                step_snapshot: None,
                degraded: false,
                hierarchy_signal: Some(signal),
            };
        }

        // Resize per-step buffers if N changed.
        if self.last_n != n {
            self.particles_backup.resize(n, bodies[0]);
            self.dcrit.resize(n, 0.0);
            self.encounter_map.resize(n, false);
            self.acc_int.resize(n, Vec3::ZERO);
            self.acc_close.resize(n, Vec3::ZERO);
            self.last_n = n;
        }

        let g_factor = ctx.g_factor;
        let mu = g_factor * bodies[0].mass;

        // ── Convert inertial → DH ──────────────────────────────────────
        self.inertial_to_dh(bodies);

        // ── dcrit table ────────────────────────────────────────────────
        self.rebuild_dcrit(bodies, dt);

        // ── 1. interaction(τ/2) ────────────────────────────────────────
        self.interaction_step(
            bodies,
            g_factor,
            0.5 * dt,
            ctx.hamiltonian_perturbations,
            ctx.non_conservative_perturbations,
        );

        // ── 2. jump(τ/2) ───────────────────────────────────────────────
        self.jump_step(bodies, 0.5 * dt);

        // ── 3. com(τ) ──────────────────────────────────────────────────
        self.com_step(dt);

        // ── 4. backup ──────────────────────────────────────────────────
        self.particles_backup.copy_from_slice(bodies);

        // ── 5. kepler(τ) ───────────────────────────────────────────────
        self.kepler_step_all(bodies, mu, dt);

        // ── 6. encounter_predict ──────────────────────────────────────
        self.encounter_predict(bodies, dt);

        // ── 7. encounter_step(τ) ──────────────────────────────────────
        self.encounter_step(bodies, g_factor, dt);

        // ── 8. jump(τ/2) ───────────────────────────────────────────────
        self.jump_step(bodies, 0.5 * dt);

        // ── 9. interaction(τ/2) ────────────────────────────────────────
        // dcrit may have grown / shrunk because particle positions and
        // velocities have changed; recompute before the closing kick.
        self.rebuild_dcrit(bodies, dt);
        self.interaction_step(
            bodies,
            g_factor,
            0.5 * dt,
            ctx.hamiltonian_perturbations,
            ctx.non_conservative_perturbations,
        );

        // ── Convert DH → inertial ──────────────────────────────────────
        self.dh_to_inertial(bodies);

        // ── Populate `acc` with the total inertial acceleration ────────
        // The K-weighted decomposition is algorithmic, not physical;
        // diagnostics consumers see the real Newtonian acceleration on
        // the post-step positions.
        acc.clear();
        acc.resize(n, Vec3::ZERO);
        for i in 0..n {
            let qi = Vec3::new(bodies[i].pos_x, bodies[i].pos_y, bodies[i].pos_z);
            for j in 0..n {
                if i == j {
                    continue;
                }
                let qj = Vec3::new(bodies[j].pos_x, bodies[j].pos_y, bodies[j].pos_z);
                let dq = qj - qi;
                let r2 = dq.length_squared().max(DCRIT_FLOOR);
                let inv_r3 = 1.0 / (r2 * r2.sqrt());
                acc[i] += dq * (g_factor * bodies[j].mass * inv_r3);
            }
        }

        StepResult {
            consumed_dt: dt,
            potential_energy: 0.0,
            used_fallback: false,
            step_snapshot: None,
            degraded: false,
            hierarchy_signal: Some(signal),
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Mercurius
    }

    fn requires_deterministic_force(&self) -> bool {
        false
    }

    fn set_hill_factor(&mut self, alpha: f64) {
        self.set_alpha(alpha);
    }

    fn hill_factor(&self) -> Option<f64> {
        Some(self.alpha())
    }

    fn resume_state(&self) -> Vec<u8> {
        mercurius_resume::encode(self)
    }

    fn restore_resume_state(
        &mut self,
        bytes: &[u8],
    ) -> Result<(), crate::physics::integrator::traits::ResumeError> {
        mercurius_resume::decode_into(self, bytes)
    }
}

mod mercurius_resume {
    use super::Mercurius;
    use crate::physics::integrator::traits::{Integrator, ResumeError};

    /// Layout: `magic(b"MER")` ‖ `version(u8 = 1)` ‖ `alpha(f64 LE)` ‖
    /// `inner_len(u32 LE)` ‖ embedded `Ias15::resume_state()` bytes.
    /// Particles_backup / dcrit / encounter_map / acc_int / acc_close
    /// are intra-step scratch (rebuilt each outer step) and excluded;
    /// COM state is reconstructed from body positions at restore time.
    const MAGIC: &[u8; 3] = b"MER";
    const VERSION: u8 = 1;

    pub fn encode(s: &Mercurius) -> Vec<u8> {
        let inner = s.ias15.resume_state();
        let mut out = Vec::with_capacity(3 + 1 + 8 + 4 + inner.len());
        out.extend_from_slice(MAGIC);
        out.push(VERSION);
        out.extend_from_slice(&s.alpha.to_le_bytes());
        out.extend_from_slice(&(inner.len() as u32).to_le_bytes());
        out.extend_from_slice(&inner);
        out
    }

    pub fn decode_into(s: &mut Mercurius, bytes: &[u8]) -> Result<(), ResumeError> {
        if bytes.len() < 16 || &bytes[..3] != MAGIC || bytes[3] != VERSION {
            return Err(ResumeError::UnsupportedFormat);
        }
        s.alpha = f64::from_le_bytes(bytes[4..12].try_into().unwrap());
        let inner_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        if bytes.len() < 16 + inner_len {
            return Err(ResumeError::Truncated);
        }
        s.ias15.restore_resume_state(&bytes[16..16 + inner_len])
    }
}

// ── CloseFieldForceModel ──────────────────────────────────────────────────────

/// Force model passed to the embedded IAS15 during the encounter step.
/// Computes the close-field acceleration: full Sun pull on every
/// non-Sun particle, plus (1−L)-weighted planet-planet from every
/// other particle (encountering or not).
struct CloseFieldForceModel<'a> {
    /// Per-particle critical radius; pair-wise reduction
    /// `dcrit_ij = max(dcrit_i, dcrit_j)`.
    dcrit: &'a [f64],
    /// True for particles in the encounter map. The Sun (index 0) is
    /// set to `true` unconditionally so the Sun pull is computed for
    /// every body; other entries gate which planets receive force.
    encounter_map: &'a [bool],
}

impl<'a> crate::physics::integrator::ForceModel for CloseFieldForceModel<'a> {
    fn compute(&mut self, bodies: &[Body], acc: &mut [Vec3]) -> f64 {
        let n = bodies.len();
        for a in acc.iter_mut().take(n) {
            *a = Vec3::ZERO;
        }
        if n < 2 {
            return 0.0;
        }

        let mut pe = 0.0_f64;
        let q0 = Vec3::new(bodies[0].pos_x, bodies[0].pos_y, bodies[0].pos_z);
        let m0 = bodies[0].mass;

        // Sun pull on each planet (full strength). Sun receives no
        // acceleration in DH coords — its velocity is pinned to zero
        // by the encounter step driver.
        for i in 1..n {
            // Compute force on every planet, but only honour it for
            // particles in the encounter map. Non-encountering
            // particles have their post-Kepler state restored at the
            // end of the encounter session, so their post-IAS15 state
            // is discarded.
            if !self.encounter_map[i] {
                continue;
            }
            let qi = Vec3::new(bodies[i].pos_x, bodies[i].pos_y, bodies[i].pos_z);
            let dq = q0 - qi;
            let r2 = dq.length_squared().max(DCRIT_FLOOR);
            let inv_r3 = 1.0 / (r2 * r2.sqrt());
            acc[i] += dq * (m0 * inv_r3);
            pe -= m0 * bodies[i].mass / r2.sqrt();
        }

        // (1−L)-weighted planet-planet pulls.
        for i in 1..n {
            if !self.encounter_map[i] {
                continue;
            }
            let qi = Vec3::new(bodies[i].pos_x, bodies[i].pos_y, bodies[i].pos_z);
            let dcrit_i = self.dcrit[i];
            for j in 1..n {
                if i == j {
                    continue;
                }
                let qj = Vec3::new(bodies[j].pos_x, bodies[j].pos_y, bodies[j].pos_z);
                let dq = qj - qi;
                let r2 = dq.length_squared().max(DCRIT_FLOOR);
                let r = r2.sqrt();
                let dcrit_pair = dcrit_i.max(self.dcrit[j]);
                let one_minus_l = 1.0 - Mercurius::changeover_l(r, dcrit_pair);
                if one_minus_l <= 0.0 {
                    continue;
                }
                let inv_r3 = 1.0 / (r * r2);
                acc[i] += dq * (one_minus_l * bodies[j].mass * inv_r3);
                if i < j {
                    pe -= one_minus_l * bodies[i].mass * bodies[j].mass / r;
                }
            }
        }

        pe
    }

    fn is_deterministic(&self) -> bool {
        true
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::integrator::WisdomHolman;
    use crate::physics::integrator::force_model::GravityForceModel;

    fn quiet_planetary() -> Vec<Body> {
        // Sun + 2 widely-separated planets on circular Keplerian
        // orbits. r ∈ {1, 2}; v chosen for circular Kepler at G·M = 1.
        vec![
            Body::star(1.0),
            Body::rocky(1.0e-6).at(1.0, 0.0).with_velocity(0.0, 1.0),
            Body::rocky(1.0e-6).at(2.0, 0.0).with_velocity(0.0, std::f64::consts::FRAC_1_SQRT_2),
        ]
    }

    fn step_via(integrator: &mut dyn Integrator, bodies: &mut [Body], dt: f64, n_steps: usize) {
        let mut force = GravityForceModel::new(0.5, 16);
        let mut acc: Vec<Vec3> = vec![Vec3::ZERO; bodies.len()];
        let hamiltonian: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> = Vec::new();
        let non_conservative: Vec<Box<dyn crate::physics::integrator::NonConservativeOperator>> =
            Vec::new();
        let mut observers: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        for _ in 0..n_steps {
            let mut ctx = IntegratorContext {
                force: &mut force,
                g_factor: 1.0,
                hamiltonian_perturbations: &hamiltonian,
                non_conservative_perturbations: &non_conservative,
                observers: &mut observers,
            };
            integrator.step(bodies, &mut ctx, dt, &mut acc);
        }
    }

    #[test]
    fn changeover_l_endpoints() {
        // L_mercury(d=0, dcrit=1) → y = -0.111... → 0
        assert_eq!(Mercurius::changeover_l(0.0, 1.0), 0.0);
        // L_mercury(d=0.1·dcrit, dcrit=1) → y = 0 → 0
        assert!(Mercurius::changeover_l(0.1, 1.0).abs() < 1.0e-15);
        // L_mercury(d=dcrit, dcrit=1) → y = 1 → 1
        assert_eq!(Mercurius::changeover_l(1.0, 1.0), 1.0);
        // L_mercury(d > dcrit) → y > 1 → clamped 1
        assert_eq!(Mercurius::changeover_l(2.0, 1.0), 1.0);
    }

    #[test]
    fn changeover_l_is_monotone_on_active_band() {
        // Band: 0.1 ≤ d ≤ 1.0 with dcrit = 1.
        let mut prev = -1.0;
        for k in 0..=20 {
            let d = 0.1 + 0.9 * (k as f64 / 20.0);
            let l = Mercurius::changeover_l(d, 1.0);
            assert!(l >= prev, "L should be non-decreasing; L({d}) = {l}, prev = {prev}");
            prev = l;
        }
    }

    #[test]
    fn changeover_l_is_c2_at_endpoints() {
        // L'(y) = 30 y² - 60 y³ + 30 y⁴ = 30 y²(1 - y)²
        // L'(0) = L'(1) = 0; L''(0) = L''(1) = 0.
        let h = 1.0e-5;
        let dcrit = 1.0;
        // First derivative near y = 0 (d = 0.1).
        let dl_at_in =
            (Mercurius::changeover_l(0.1 + h, dcrit) - Mercurius::changeover_l(0.1, dcrit)) / h;
        let dl_at_out =
            (Mercurius::changeover_l(1.0, dcrit) - Mercurius::changeover_l(1.0 - h, dcrit)) / h;
        assert!(dl_at_in.abs() < 1e-3, "L'(at 0.1 dcrit) ≈ 0; got {dl_at_in}");
        assert!(dl_at_out.abs() < 1e-3, "L'(at dcrit) ≈ 0; got {dl_at_out}");
    }

    #[test]
    fn refuses_non_hierarchical_step() {
        let mut bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut force = GravityForceModel::new(0.5, 16);
        let mut acc: Vec<Vec3> = vec![Vec3::ZERO; bodies.len()];
        let hamiltonian: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> = Vec::new();
        let non_conservative: Vec<Box<dyn crate::physics::integrator::NonConservativeOperator>> =
            Vec::new();
        let mut observers: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        let mut ctx = IntegratorContext {
            force: &mut force,
            g_factor: 1.0,
            hamiltonian_perturbations: &hamiltonian,
            non_conservative_perturbations: &non_conservative,
            observers: &mut observers,
        };
        let mut merc = Mercurius::new();
        let result = merc.step(&mut bodies, &mut ctx, 0.01, &mut acc);
        assert_eq!(result.consumed_dt, 0.0);
        assert!(result.used_fallback);
    }

    #[test]
    fn quiet_system_takes_no_encounters_at_alpha_zero() {
        // α = 0 collapses the Hill criterion; only physical-radius and
        // velocity criteria contribute. On the quiet system at dt =
        // 1e-3, dcrit ≈ max(velocity·0.4·dt, 2·r_phys) ≈ O(1e-3) ≪
        // separation O(1). The encounter detector should never fire.
        let dt = 1.0e-3;
        let mut bodies = quiet_planetary();
        let mut merc = Mercurius::with_alpha(0.0);
        step_via(&mut merc, &mut bodies, dt, 100);
        // After 100 steps, the encounter map should remain entirely
        // false on planet entries (it's only populated when at least
        // one pair triggers in the most recent step).
        assert!(
            !merc.any_planet_in_encounter(),
            "quiet system at α=0 should not have triggered any encounters"
        );
    }

    #[test]
    fn tier1_no_encounters_matches_wh_to_split_error_floor() {
        // No-encounter limit (α = 0, quiet system): Mercurius reduces
        // to a 5-stage symplectic split (int(τ/2) jmp(τ/2) com(τ) +
        // kep(τ) + jmp(τ/2) int(τ/2)). Apsis WH uses (kick(τ/2)
        // kep(τ) jmp(τ) kick(τ/2)). Both are 2nd-order symplectic but
        // the jump placement differs, producing an O(τ²·m_p/m_0) per-
        // step truncation difference.
        //
        // Bound estimate at dt = 1e-3, m_p/m_0 = 1e-6, n_steps = 200:
        //   200 · (1e-3)² · 1e-6 ≈ 2e-13 cumulative trajectory drift,
        //   plus round-off accumulation ~ N·ε ≈ 1e-14.
        // Total expected: well below 1e-5 on |Δr|/r.
        let dt = 1.0e-3;
        let n_steps = 200;
        let mut bodies_merc = quiet_planetary();
        let mut bodies_wh = quiet_planetary();
        let mut merc = Mercurius::with_alpha(0.0);
        let mut wh = WisdomHolman::new();

        step_via(&mut merc, &mut bodies_merc, dt, n_steps);
        step_via(&mut wh, &mut bodies_wh, dt, n_steps);

        for (bm, bw) in bodies_merc.iter().zip(bodies_wh.iter()) {
            let dx = bm.pos_x - bw.pos_x;
            let dy = bm.pos_y - bw.pos_y;
            let dz = bm.pos_z - bw.pos_z;
            let r = (bm.pos_x.powi(2) + bm.pos_y.powi(2) + bm.pos_z.powi(2)).sqrt().max(1.0e-30);
            let rel = (dx * dx + dy * dy + dz * dz).sqrt() / r;
            assert!(
                rel < 1.0e-5,
                "no-encounter limit vs WisdomHolman: |Δr|/r = {rel:.3e}, expected < 1e-5"
            );
        }
    }

    #[test]
    fn close_pair_engages_encounter_step() {
        // Force a close encounter via large α so the Hill-radius
        // criterion dominates: planets at separation 1.0 with
        // dcrit = α · a · (m/3M)^(1/3) ≈ 100 · 1.0 · (1e-3/3)^(1/3) ≈ 7.
        // Both planets fall inside each other's dcrit on the first
        // step. The encounter step must engage and consume the outer
        // dt without crashing.
        let bodies = vec![
            Body::star(1.0),
            Body::rocky(1.0e-3).at(1.0, 0.0).with_velocity(0.0, 1.0),
            Body::rocky(1.0e-3).at(2.0, 0.0).with_velocity(0.0, std::f64::consts::FRAC_1_SQRT_2),
        ];
        let mut bs = bodies;
        let mut merc = Mercurius::with_alpha(100.0);
        step_via(&mut merc, &mut bs, 1.0e-3, 5);
        assert!(
            merc.any_planet_in_encounter(),
            "large α should engage the encounter step on the first call"
        );
        // Sanity: bodies should still be finite after IAS15 sub-integration.
        for b in &bs {
            assert!(
                b.pos_x.is_finite()
                    && b.pos_y.is_finite()
                    && b.pos_z.is_finite()
                    && b.vel_x.is_finite()
                    && b.vel_y.is_finite()
                    && b.vel_z.is_finite(),
                "encounter step left a body with non-finite kinematics"
            );
        }
    }

    #[test]
    fn no_encounter_limit_conserves_energy() {
        // Independent sanity: even without comparing to WH, a quiet
        // system under Mercurius should not visibly drift in energy.
        let dt = 1.0e-3;
        let bodies0 = quiet_planetary();
        let mut bodies = bodies0.clone();
        let mut merc = Mercurius::with_alpha(0.0);

        let energy = |bs: &[Body]| -> f64 {
            let mut ke = 0.0;
            let mut pe = 0.0;
            for (i, b) in bs.iter().enumerate() {
                ke += 0.5 * b.mass * (b.vel_x.powi(2) + b.vel_y.powi(2) + b.vel_z.powi(2));
                for j in (i + 1)..bs.len() {
                    let dx = b.pos_x - bs[j].pos_x;
                    let dy = b.pos_y - bs[j].pos_y;
                    let dz = b.pos_z - bs[j].pos_z;
                    let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-30);
                    pe -= b.mass * bs[j].mass / r;
                }
            }
            ke + pe
        };

        let e0 = energy(&bodies);
        step_via(&mut merc, &mut bodies, dt, 1000);
        let e1 = energy(&bodies);
        let drift = ((e1 - e0) / e0.abs().max(1.0e-30)).abs();
        assert!(
            drift < 1.0e-8,
            "1000 quiet steps: |ΔE/E0| = {drift:.3e}, expected < 1e-8 (symplectic)"
        );
    }

    struct ConstantYKickOnPlanets {
        a_y: f64,
    }

    impl crate::physics::integrator::Operator for ConstantYKickOnPlanets {}

    impl crate::physics::integrator::HamiltonianOperator for ConstantYKickOnPlanets {
        fn accumulate_force(&self, bodies: &[crate::domain::body::Body], acc: &mut [Vec3]) {
            for i in 1..bodies.len() {
                acc[i].y += self.a_y;
            }
        }

        fn potential(
            &self,
            bodies: &[crate::domain::body::Body],
        ) -> crate::physics::integrator::Potential {
            // V = -a_y * Σ y_i so that −∂V/∂y_i = a_y.
            crate::physics::integrator::Potential::Value(
                -self.a_y * bodies.iter().skip(1).map(|b| b.pos_y).sum::<f64>(),
            )
        }
    }

    #[test]
    fn registered_perturbations_are_honored_by_interaction_step() {
        let dt = 1.0e-3;
        let n_steps = 100;
        let a_y = 1.0e-6;

        let baseline_bodies = quiet_planetary();
        let mut bodies_no_pert = baseline_bodies.clone();
        let mut bodies_with_pert = baseline_bodies;

        let mut merc_a = Mercurius::with_alpha(0.0);
        let mut merc_b = Mercurius::with_alpha(0.0);
        let mut force_a = GravityForceModel::new(0.5, 16);
        let mut force_b = GravityForceModel::new(0.5, 16);

        let no_h: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> = Vec::new();
        let with_h: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> =
            vec![Box::new(ConstantYKickOnPlanets { a_y })];
        let nc: Vec<Box<dyn crate::physics::integrator::NonConservativeOperator>> = Vec::new();
        let mut obs_a: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        let mut obs_b: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();

        let mut acc_a: Vec<Vec3> = vec![Vec3::ZERO; 3];
        let mut acc_b: Vec<Vec3> = vec![Vec3::ZERO; 3];

        for _ in 0..n_steps {
            let mut ctx_a = IntegratorContext {
                force: &mut force_a,
                g_factor: 1.0,
                hamiltonian_perturbations: &no_h,
                non_conservative_perturbations: &nc,
                observers: &mut obs_a,
            };
            merc_a.step(&mut bodies_no_pert, &mut ctx_a, dt, &mut acc_a);

            let mut ctx_b = IntegratorContext {
                force: &mut force_b,
                g_factor: 1.0,
                hamiltonian_perturbations: &with_h,
                non_conservative_perturbations: &nc,
                observers: &mut obs_b,
            };
            merc_b.step(&mut bodies_with_pert, &mut ctx_b, dt, &mut acc_b);
        }

        let total_t = (n_steps as f64) * dt;
        let expected_dvy = a_y * total_t;
        for i in 1..bodies_no_pert.len() {
            let dvy = bodies_with_pert[i].vel_y - bodies_no_pert[i].vel_y;
            assert!(
                (dvy - expected_dvy).abs() / expected_dvy < 5.0e-2,
                "planet {i}: Δvy = {dvy:.6e}, expected ~{expected_dvy:.6e} (within 5%)",
            );
        }
    }
}
