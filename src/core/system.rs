//! Simulation orchestrator for an N-body gravitational system.
//!
//! This module defines the [`System`] type, responsible for advancing
//! the state of a set of massive bodies interacting via gravity.
//!
//! ## Design goals
//!
//! - Preserve physical consistency with a symplectic integrator (Velocity Verlet)
//! - Maintain deterministic and reproducible evolution
//! - Provide diagnostics for energy and angular momentum conservation
//! - Support Barnes–Hut acceleration for scalable simulations
//!
//! ## Important assumptions
//!
//! - Time step (`dt`) is constant (required for symplectic behavior)
//! - The force field is evaluated consistently at well-defined points
//! - No discrete events (e.g., collisions) are applied within integration steps
//!
//! ## Notes
//!
//! This system is intended for scientific and numerical experiments in
//! gravitational dynamics, not for general-purpose physics engines.

use crate::core::adaptive::{
    AccelerationStats, DtAdaptationConfig, DtController, DtMode, ThetaController,
};
use crate::core::body::{Body, NamedBody};
use crate::core::calibration;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::metrics::Metrics;
use crate::core::trail_buffer::{TrailBuffer, adaptive_capacity};

const MASS_TO_SOLAR: f64 = 1.0;
const RADIUS_TO_SOLAR: f64 = 1.0 / 0.00465;
const L_SUN: f64 = 1.0;

/// Minimum ratio `M_central / Σ m_i (i > 0)` required for the Wisdom–Holman
/// integrator to be considered valid.
///
/// WH is formally correct for any mass ratio, but its accuracy as an
/// *integrator* is proportional to the perturbation parameter
/// `ε = m_perturber / M_central`.  Below a ratio of 10 the perturbations
/// are large enough that the Keplerian splitting breaks down and the method
/// produces silently wrong trajectories.
const WH_DOMINANCE_RATIO: f64 = 10.0;

/// Number of bodies that actually need individual trail rendering.
///
/// Belt members and sub-threshold bodies are excluded because their trails are
/// suppressed by the renderer anyway. Using this count for ring-buffer capacity
/// allocation keeps GPU memory proportional to what's actually rendered, not
/// to the total body count (which can be dominated by asteroid belts).
/// Generate an auto-name for a new body given existing names.
/// Counts existing names that start with the material prefix and appends N+1.
fn auto_name(material: crate::core::materials::Material, existing: &[String]) -> String {
    let prefix = material.display_name();
    let count = existing.iter().filter(|n| n.starts_with(prefix)).count() + 1;
    format!("{prefix} {count}")
}

fn resolved_name(
    explicit: Option<String>,
    material: crate::core::materials::Material,
    existing: &[String],
) -> String {
    explicit
        .map(|name| name.trim().to_owned())
        .filter(|name| !name.is_empty())
        .unwrap_or_else(|| auto_name(material, existing))
}

fn trail_body_count(bodies: &[Body]) -> usize {
    if bodies.is_empty() {
        return 0;
    }
    let max_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);
    if max_mass <= 0.0 {
        return bodies.len();
    }
    bodies
        .iter()
        .filter(|b| b.mass / max_mass > 1e-6)
        .count()
        .max(1)
}
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{
    Integrator, PerturbationForce, Y4_C, Y4_D, drift, evaluate_accelerations, kick,
};
use crate::physics::orbital::{self, OrbitalElements};

/// Central simulation state for an N-body gravitational system.
pub struct System {
    /// Bodies participating in the simulation.
    bodies: Vec<Body>,

    /// GPU-ready ring buffer of trail positions and colours.
    trail_buf: TrailBuffer,
    trail_every: usize,

    /// Total mass of the system (used for COM recentering).
    total_mass: f64,

    /// Last computed energies.
    last_kinetic: f64,
    last_potential: f64,

    /// Initial total energy (used as reference).
    initial_energy: Option<f64>,

    /// Relative energy error (diagnostic only).
    rel_energy_error: f64,

    /// Barnes–Hut engine for approximate force computation.
    engine: BarnesHutEngine,

    /// Scratch buffer for accelerations.
    scratch_acc: Vec<(f64, f64)>,

    /// Barnes–Hut opening angle parameter (θ).
    theta: f64,

    /// Active integration algorithm.
    integrator: Integrator,

    /// Cached osculating orbital elements — one slot per body.
    /// Updated on demand via [`System::update_orbital_elements`], not every step.
    orbital_cache: Vec<Option<OrbitalElements>>,

    /// Global Plummer softening scale applied on top of the per-body
    /// mass-proportional default: `ε = EPS_BASE · m^(1/3) · softening_scale`.
    softening_scale: f64,

    /// Diagnostics subsystem.
    diagnostics: DiagnosticsComputer,
    last_diag: SimulationDiagnostics,

    /// Step counter.
    steps: u64,

    /// Total simulated time elapsed (t = steps × dt, but tracked as f64
    /// so it remains correct even if dt changes mid-run).
    t: f64,

    /// Timestep currently used by the integrator.
    ///
    /// In [`DtMode::Fixed`] this always equals `user_dt`.
    /// In [`DtMode::Adaptive`] it is the output of [`DtController`] and may
    /// differ from `user_dt`.
    current_dt: f64,

    /// User-requested timestep — the value set via [`set_dt`] and used as the
    /// proposed baseline by the adaptive controller.
    user_dt: f64,

    /// Timestep management policy.  Default: [`DtMode::Fixed`].
    ///
    /// See [`DtMode`] for the full scientific rationale and the consequences
    /// of each choice.
    dt_mode: DtMode,

    /// Adaptive timestep controller.  Only consulted when
    /// `dt_mode == DtMode::Adaptive`.
    dt_ctrl: DtController,

    /// Adaptive Barnes–Hut opening-angle controller.  Only active when
    /// `adaptive_theta == true`.
    theta_ctrl: ThetaController,

    /// Whether the adaptive θ controller is active.  Default: `false`.
    ///
    /// Note: varying θ between steps changes force accuracy per step but does
    /// **not** break symplecticity.  For reproducible force-accuracy budgets in
    /// published runs, keep this `false` and set θ manually.
    adaptive_theta: bool,

    /// Gravitational scaling factor (G multiplier).
    g_factor: f64,

    /// Initial angular momentum (z-component) used as reference.
    initial_angular_momentum: Option<f64>,

    /// Relative angular momentum error (diagnostic only).
    rel_angular_momentum_error: f64,

    /// Absolute angular momentum error (always meaningful).
    abs_angular_momentum_error: f64,

    /// Human-readable label for each body, parallel to `bodies`.
    /// Kept separate because `Body` is `Copy` and cannot own a `String`.
    names: Vec<String>,

    /// Minimum pairwise separation cached from the most recent step.
    r_min: f64,

    /// Maximum effective pairwise softening length cached from the most recent step.
    softening_max: f64,

    perturbations: Vec<Box<dyn PerturbationForce>>,
}

impl System {
    /// Creates a new simulation system.
    ///
    /// # Parameters
    ///
    /// - `bodies`: Initial set of bodies
    /// - `theta`: Barnes–Hut opening angle (controls accuracy vs performance)
    /// - `dt`: Fixed time step
    /// - `max_depth`: Maximum tree depth for Barnes–Hut
    /// - `trail_every`: Sampling interval for trails (ring-buffer depth is
    ///   chosen automatically via [`adaptive_capacity`])
    ///
    /// # Notes
    ///
    /// - Smaller `theta` increases accuracy (approaches O(N²))
    /// - Smaller `dt` improves stability and energy conservation
    pub fn new(
        bodies: Vec<Body>,
        theta: f64,
        dt: f64,
        max_depth: usize,
        trail_every: usize,
    ) -> Self {
        let n = bodies.len();
        let trail_n = trail_body_count(&bodies);
        let cap = adaptive_capacity(trail_n.max(1));
        let mut trail_buf = TrailBuffer::new_with_capacity(n, cap);
        trail_buf.update_colors(&bodies);

        let total_mass = bodies.iter().map(|b| b.mass).sum();
        let names = bodies
            .iter()
            .map(|b| auto_name(b.material, &[]))
            .collect::<Vec<_>>();
        // Re-generate with correct counters (so Star 1, Star 2 … instead of all "Star 1")
        let names = {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for b in &bodies {
                acc.push(auto_name(b.material, &acc));
            }
            acc
        };

        let (r_min, softening_max) = Self::compute_closeness(&bodies);

        Self {
            bodies,
            trail_buf,
            trail_every: trail_every.max(1),
            total_mass,
            last_kinetic: 0.0,
            last_potential: 0.0,
            initial_energy: None,
            rel_energy_error: 0.0,
            engine: BarnesHutEngine::new(max_depth),
            scratch_acc: Vec::new(),
            theta,
            integrator: Integrator::VelocityVerlet,
            orbital_cache: Vec::new(),
            softening_scale: 1.0,
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            steps: 0,
            t: 0.0,
            current_dt: dt,
            user_dt: dt,
            dt_mode: DtMode::Fixed,
            dt_ctrl: DtController::new(DtAdaptationConfig {
                // Always ready; `dt_mode` is the master switch.
                enabled: true,
                min_dt: 1e-9,
                max_dt: 1e6,
                target_rel_energy_error: 1e-6,
                accel_epsilon: 0.1,
                grow_limit: 1.2,
                shrink_limit: 0.5,
                dt_slew_fraction: 0.1,
            }),
            theta_ctrl: ThetaController::new(1e-3, 0.05, 1.5).with_initial_theta(theta),
            adaptive_theta: false,
            g_factor: 1.0,
            initial_angular_momentum: None,
            rel_angular_momentum_error: 0.0,
            abs_angular_momentum_error: 0.0,
            names,
            r_min,
            softening_max,
            perturbations: Vec::new(),
        }
    }
}

impl System {
    /// Advance the simulation by one time step using the configured integrator.
    pub fn step(&mut self) {
        match self.integrator {
            Integrator::VelocityVerlet => self.step_vv(),
            Integrator::Yoshida4 => self.step_yoshida4(),
            Integrator::WisdomHolman => self.step_wisdom_holman(),
        }

        let dt = self.current_dt;
        self.steps += 1;
        self.t += dt;

        self.last_diag = self
            .diagnostics
            .compute(&self.scratch_acc, &self.bodies, dt);

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        // ── Adaptive controllers ──────────────────────────────────────────────
        // Run after diagnostics and energy tracking so all inputs are fresh.
        //
        // DtMode::Fixed: current_dt is always user_dt; the DtController is
        //   never consulted.  The symplectic guarantee is fully preserved.
        //
        // DtMode::Adaptive: DtController modulates current_dt each step.
        //   The symplectic guarantee is broken — see DtMode documentation.
        self.current_dt = match self.dt_mode {
            DtMode::Fixed => self.user_dt,
            DtMode::Adaptive => {
                let stats = AccelerationStats::new(self.last_diag.max_acc, self.last_diag.jerk);
                self.dt_ctrl
                    .update(self.user_dt, self.rel_energy_error, stats)
            }
        };

        // Adaptive θ: does not break symplecticity but does vary force accuracy
        // per step.  Only active when explicitly enabled.
        // `theta_error_proxy` requires the quadtree to be built — guaranteed
        // here because the final force eval of the step leaves the tree populated.
        // For N ≤ EXACT_THRESHOLD the tree is empty and the proxy returns 0.0.
        if self.adaptive_theta && !self.bodies.is_empty() {
            let e_theta = self.engine.theta_error_proxy(0, &self.bodies, self.theta);
            self.theta = self.theta_ctrl.update(e_theta, self.current_dt);
        }

        // Periodically remove COM drift.  The trail buffer is translated by
        // the same vector so stored positions remain consistent.
        if self.steps % 97 == 0 {
            if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
                calibration::apply_body_shift(&mut self.bodies, dx, dy);
                self.trail_buf.translate(-dx as f32, -dy as f32);
            }
        }

        // Update softening diagnostics every step (O(N²) but bounded by threshold).
        let (r_min, soft_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = soft_max;
    }

    // ── Velocity Verlet (KDK leapfrog, 2nd-order symplectic) ─────────────────────
    //
    // Scheme: F(t) → kick(½dt) → drift(dt) → F(t+dt) → kick(½dt)
    //
    // The two half-kicks bracketing the drift are equivalent to a single
    // full kick at the midpoint, giving 2nd-order accuracy with one force
    // evaluation per amortised step.

    fn step_vv(&mut self) {
        let dt = self.current_dt;
        let theta = self.theta;

        let raw_pe =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);
        self.last_potential = self.scale_acc_and_pe(raw_pe);
        self.apply_perturbations();

        kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);
        drift(&mut self.bodies, dt);

        let raw_pe =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);
        self.last_potential = self.scale_acc_and_pe(raw_pe);
        self.apply_perturbations();

        kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);
    }

    // ── Yoshida 4th-order (Forest–Ruth DKD composition) ──────────────────────────
    //
    // Scheme: drift(c₀) → F → kick(d₀) → drift(c₁) → F → kick(d₁) → drift(c₂) → F → kick(d₂) → drift(c₃)
    //
    // The middle kick coefficient d₁ = w₀ ≈ −1.70 is negative, meaning the
    // second sub-step is a backward kick in time. This is not a bug — it is the
    // mechanism by which leading error terms cancel to achieve 4th-order accuracy.
    //
    // References:
    //   Forest & Ruth (1990). Nucl. Instrum. Methods Phys. Res. A 290, 395–400.
    //   Yoshida (1990). Phys. Lett. A 150, 262–268.

    fn step_yoshida4(&mut self) {
        let dt = self.current_dt;
        let theta = self.theta;

        for i in 0..3 {
            drift(&mut self.bodies, Y4_C[i] * dt);

            let raw_pe = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe);
            self.apply_perturbations();

            kick(&mut self.bodies, &self.scratch_acc, Y4_D[i] * dt);
        }

        drift(&mut self.bodies, Y4_C[3] * dt);

        // ── Consistent energy snapshot ────────────────────────────────────────
        // After the final drift the phase-space point is (q(t+dt), v(t+dt)).
        // `last_potential` still holds PE(q‴) — the potential at the positions
        // BEFORE the drift — which is inconsistent with the current body state.
        // Without this correction, `update_energy_tracking` computes
        //   E_shadow = KE(v(t+dt)) + PE(q‴)
        // which oscillates at O(dt) rather than O(dt⁴), making Y4 appear
        // dramatically worse than VV in the metrics panel.
        //
        // Re-evaluating the potential at q(t+dt) costs one additional force call
        // per step (3 → 4 total).  The accelerations are also updated so that
        // `scratch_acc` is consistent with the final positions, which improves
        // jerk diagnostics on the next step.
        {
            let raw_pe = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe);
        }
    }

    /// Evaluates inter-planetary perturbation accelerations (excluding the central
    /// body), computes the heliocentric indirect-term correction, and applies a
    /// velocity kick of magnitude `dt` to all planets.
    ///
    /// # Heliocentric indirect term
    ///
    /// In heliocentric coordinates the perturbation Hamiltonian contains a
    /// momentum-dependent cross term that contributes an additional acceleration
    ///
    /// ```text
    /// a_indirect,i = −(Σ_j m_j a_j) / M₀
    /// ```
    ///
    /// where `a_j` are the **raw** (pre-scaled) inter-planetary accelerations and
    /// the sum runs over all planets. The indirect term shares the same `g_factor`
    /// scaling as the direct perturbation and must be computed before
    /// [`scale_acc_and_pe`] is called to avoid a spurious double-application.
    ///
    /// # Returns
    ///
    /// The total gravitational potential at the evaluated positions:
    /// inter-planetary interaction energy plus the central `−μ Σ mᵢ/rᵢ` term.
    fn wh_kick(&mut self, dt: f64, mu: f64) -> f64 {
        let theta = self.theta;
        let total_m0 = self.bodies[0].mass;

        let raw_pe = evaluate_accelerations(
            &self.bodies[1..],
            theta,
            &mut self.engine,
            &mut self.scratch_acc,
        );

        let (ax_bary_raw, ay_bary_raw) = self
            .scratch_acc
            .iter()
            .zip(self.bodies[1..].iter())
            .fold((0.0_f64, 0.0_f64), |(ax, ay), (&(axi, ayi), b)| {
                (ax + b.mass * axi, ay + b.mass * ayi)
            });

        let indirect_x_raw = -ax_bary_raw / total_m0;
        let indirect_y_raw = -ay_bary_raw / total_m0;

        let potential = self.scale_acc_and_pe(raw_pe) + self.central_potential(mu);

        // Non-gravitational perturbations act on planets only (bodies[1..]),
        // aligned with scratch_acc which has length N-1.
        self.apply_perturbations_planets();

        let indirect_x = indirect_x_raw * self.g_factor;
        let indirect_y = indirect_y_raw * self.g_factor;
        for (i, &(ax, ay)) in self.scratch_acc.iter().enumerate() {
            self.bodies[i + 1].vx += (ax + indirect_x) * dt;
            self.bodies[i + 1].vy += (ay + indirect_y) * dt;
        }

        potential
    }

    // ── Wisdom–Holman mixed-variable symplectic (2nd-order) ──────────────────────
    //
    // Scheme (heliocentric frame):
    //
    //   kick_pert(½dt)  →  drift_Kepler(dt)  →  kick_pert(½dt)
    //
    // The Hamiltonian is split as H = H_Kepler + H_pert. H_Kepler is integrated
    // exactly via the analytic universal-variable propagator; H_pert contributes
    // velocity kicks that include the heliocentric indirect term (momentum
    // cross-term) required to preserve symplecticity.
    //
    // The integration is performed entirely in heliocentric coordinates and
    // converted back to the inertial barycentric frame at the end of each step
    // via total-momentum conservation.
    //
    // Assumptions:
    //   - `bodies[0]` is the dominant central mass.
    //   - The system is hierarchical: M_central ≫ all other masses.
    //   - Close encounters between planets degrade accuracy; switch to Yoshida4
    //     if any separation approaches the mutual Hill radius.
    //
    // References:
    //   Wisdom, J. & Holman, M. (1991). Astron. J. 102, 1528–1538.

    fn step_wisdom_holman(&mut self) {
        // Safety guard: WH is only valid when bodies[0] dominates the system.
        // If the criterion is not met, fall back to Yoshida4 silently so the
        // simulation produces physically correct (if slower) results rather than
        // a silently wrong trajectory.
        if !self.is_wh_suitable() {
            self.step_yoshida4();
            return;
        }

        let dt = self.current_dt;
        let mu = self.g_factor * self.bodies[0].mass;
        let total_m0 = self.bodies[0].mass;

        // ── To heliocentric frame ─────────────────────────────────────────────
        let (cx0, cy0, cvx0, cvy0) = (
            self.bodies[0].x,
            self.bodies[0].y,
            self.bodies[0].vx,
            self.bodies[0].vy,
        );
        for b in &mut self.bodies[1..] {
            b.x -= cx0;
            b.y -= cy0;
            b.vx -= cvx0;
            b.vy -= cvy0;
        }

        // ── First half-kick (perturbations + indirect term) ───────────────────
        // Potential at t is recorded but immediately overwritten; diagnostics
        // always report the end-of-step value at x(t + dt).
        let _ = self.wh_kick(0.5 * dt, mu);

        // ── Exact Keplerian drift ─────────────────────────────────────────────
        // bodies[0] is the origin and remains at rest in this frame.
        for i in 1..self.bodies.len() {
            let b = &self.bodies[i];
            let (nx, ny, nvx, nvy) =
                crate::physics::kepler::kepler_step(b.x, b.y, b.vx, b.vy, dt, mu);
            self.bodies[i].x = nx;
            self.bodies[i].y = ny;
            self.bodies[i].vx = nvx;
            self.bodies[i].vy = nvy;
        }

        // ── Second half-kick (perturbations + indirect term) ──────────────────
        // Potential at x(t + dt) — this is the value exposed by metrics().
        self.last_potential = self.wh_kick(0.5 * dt, mu);

        // ── Back to inertial (barycentric) frame ──────────────────────────────
        // Recover the central-body velocity from total-momentum conservation:
        //   M₀ v₀ = −Σᵢ mᵢ vᵢ
        // then shift all positions by the central body's inertial displacement.
        let (px, py) = self.bodies[1..]
            .iter()
            .fold((0.0_f64, 0.0_f64), |(px, py), b| {
                (px + b.mass * b.vx, py + b.mass * b.vy)
            });

        self.bodies[0].vx = -px / total_m0;
        self.bodies[0].vy = -py / total_m0;
        self.bodies[0].x += self.bodies[0].vx * dt;
        self.bodies[0].y += self.bodies[0].vy * dt;

        let (cx1, cy1, cvx1, cvy1) = (
            self.bodies[0].x,
            self.bodies[0].y,
            self.bodies[0].vx,
            self.bodies[0].vy,
        );
        for b in &mut self.bodies[1..] {
            b.x += cx1;
            b.y += cy1;
            b.vx += cvx1;
            b.vy += cvy1;
        }
    }
}

impl System {
    /// Compute the minimum pairwise separation and maximum effective softening
    /// length over all body pairs.
    ///
    /// Skipped (returns sentinels) when N < 2 or N > [`N_CLOSENESS_THRESHOLD`],
    /// to keep overhead bounded for large asteroid-belt simulations.
    fn compute_closeness(bodies: &[Body]) -> (f64, f64) {
        const N_CLOSENESS_THRESHOLD: usize = 512;

        if bodies.len() < 2 || bodies.len() > N_CLOSENESS_THRESHOLD {
            return (f64::MAX, 0.0);
        }

        let mut r_min = f64::MAX;
        let mut soft_max = 0.0_f64;

        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                let dx = bodies[i].x - bodies[j].x;
                let dy = bodies[i].y - bodies[j].y;
                let r = (dx * dx + dy * dy).sqrt();
                if r < r_min {
                    r_min = r;
                }
                let eps2_ij = (bodies[i].softening * bodies[i].softening
                    + bodies[j].softening * bodies[j].softening)
                    * 0.5;
                let eps_ij = eps2_ij.sqrt();
                if eps_ij > soft_max {
                    soft_max = eps_ij;
                }
            }
        }

        (r_min, soft_max)
    }

    /// Multiply every acceleration in `scratch_acc` and the raw potential
    /// by `g_factor`, then return the scaled potential.
    ///
    /// The engine always uses the hard-coded `G₀ = 1.0`; multiplying the
    /// output is equivalent to running with `G_eff = G₀ · g_factor`.
    fn scale_acc_and_pe(&mut self, raw_pe: f64) -> f64 {
        if (self.g_factor - 1.0).abs() > 1e-15 {
            for a in &mut self.scratch_acc {
                a.0 *= self.g_factor;
                a.1 *= self.g_factor;
            }
        }
        raw_pe * self.g_factor
    }
}

impl System {
    /// Updates energy diagnostics for the current simulation state.
    fn update_energy_tracking(&mut self) {
        let kinetic = kinetic_energy(&self.bodies);
        self.last_kinetic = kinetic;

        let total = total_energy(kinetic, self.last_potential);

        let baseline = match self.initial_energy {
            Some(v) => v,
            None => {
                self.initial_energy = Some(total);
                total
            }
        };

        let denom = baseline.abs().max(1e-12);
        self.rel_energy_error = (total - baseline) / denom;
    }

    /// Updates angular momentum diagnostics.
    fn update_angular_momentum_tracking(&mut self) {
        let lz = angular_momentum_z(&self.bodies);

        let baseline = match self.initial_angular_momentum {
            Some(v) => v,
            None => {
                self.initial_angular_momentum = Some(lz);
                lz
            }
        };

        self.abs_angular_momentum_error = (lz - baseline).abs();

        let denom = baseline.abs().max(1e-12);
        self.rel_angular_momentum_error = (lz - baseline) / denom;
    }

    fn central_potential(&self, mu: f64) -> f64 {
        self.bodies[1..]
            .iter()
            .map(|b| {
                let r = (b.x * b.x + b.y * b.y).sqrt().max(1e-30);
                -mu * b.mass / r
            })
            .sum()
    }
}

impl System {
    /// Registers a non-gravitational perturbation force.
    ///
    /// The force is applied at every subsequent integration step until
    /// removed via [`clear_perturbations`]. Multiple perturbations are
    /// applied in registration order and accumulate additively.
    ///
    /// # Example
    ///
    /// ```rust
    /// use physics::radiation::perturbation::RadiationField;
    ///
    /// system.add_perturbation(Box::new(RadiationField::new(source, n, true)));
    /// ```
    pub fn add_perturbation(&mut self, p: Box<dyn PerturbationForce>) {
        self.perturbations.push(p);
    }

    /// Removes all registered perturbation forces.
    pub fn clear_perturbations(&mut self) {
        self.perturbations.clear();
    }

    /// Returns the number of currently registered perturbations.
    pub fn perturbation_count(&self) -> usize {
        self.perturbations.len()
    }

    /// Accumulates all registered perturbation forces into `scratch_acc`.
    ///
    /// Must be called **after** [`scale_acc_and_pe`] so gravitational and
    /// non-gravitational contributions are separable in diagnostics.
    /// Perturbation forces are independent of `g_factor`.
    fn apply_perturbations(&mut self) {
        if self.perturbations.is_empty() {
            return;
        }
        for p in &self.perturbations {
            p.accumulate(&self.bodies, &mut self.scratch_acc);
        }
    }

    /// Variant of [`apply_perturbations`] for use inside [`wh_kick`].
    ///
    /// During the Wisdom–Holman sub-step the Barnes–Hut tree is built from
    /// `bodies[1..]` only, so `scratch_acc` has length `N − 1`.
    /// This helper passes the matching slice of `bodies` to each perturbation
    /// so indices remain aligned.
    fn apply_perturbations_planets(&mut self) {
        if self.perturbations.is_empty() {
            return;
        }
        let bodies_planets = &self.bodies[1..];
        for p in &self.perturbations {
            p.accumulate_offset(bodies_planets, &mut self.scratch_acc, 1);
        }
    }
}

impl System {
    /// Adds a new body to the simulation.
    ///
    /// The trail buffer is reset to accommodate the new body count; trail
    /// history is lost.  Energy baseline is reset because the system
    /// topology has changed.
    pub fn add_body(&mut self, mut body: Body) {
        use crate::core::body::default_softening;
        body.sync_physical_properties();
        if (self.softening_scale - 1.0).abs() > 1e-15 {
            body.softening = default_softening(body.mass) * self.softening_scale;
        }
        self.total_mass += body.mass;
        self.names.push(auto_name(body.material, &self.names));
        body.update_luminosity(MASS_TO_SOLAR, RADIUS_TO_SOLAR, L_SUN);
        self.bodies.push(body);

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
    }

    /// Adds a single body while preserving an explicit display name when given.
    pub fn add_named_body(&mut self, named_body: NamedBody) {
        self.add_named_bodies(vec![named_body]);
    }

    /// Add multiple bodies in a single batch.
    ///
    /// More efficient than calling [`add_body`] in a loop: the trail buffer is
    /// reset only once and the energy baseline is invalidated once.
    pub fn add_bodies(&mut self, new_bodies: Vec<Body>) {
        use crate::core::body::default_softening;
        for mut body in new_bodies {
            body.sync_physical_properties();
            if (self.softening_scale - 1.0).abs() > 1e-15 {
                body.softening = default_softening(body.mass) * self.softening_scale;
            }
            self.total_mass += body.mass;
            self.names.push(auto_name(body.material, &self.names));
            body.update_luminosity(MASS_TO_SOLAR, RADIUS_TO_SOLAR, L_SUN);
            self.bodies.push(body);
        }

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);
        self.initial_energy = None;
    }

    /// Add multiple bodies in a single batch while preserving explicit names.
    ///
    /// Each `NamedBody` may provide a pre-authored display name. Bodies without
    /// an explicit name fall back to the standard material-based naming scheme.
    pub fn add_named_bodies(&mut self, new_bodies: Vec<NamedBody>) {
        use crate::core::body::default_softening;
        for mut named_body in new_bodies {
            let mut body = named_body.body;
            body.sync_physical_properties();
            if (self.softening_scale - 1.0).abs() > 1e-15 {
                body.softening = default_softening(body.mass) * self.softening_scale;
            }
            self.total_mass += body.mass;
            let name = resolved_name(named_body.name.take(), body.material, &self.names);
            body.update_luminosity(MASS_TO_SOLAR, RADIUS_TO_SOLAR, L_SUN);
            self.names.push(name);
            self.bodies.push(body);
        }

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);
        self.initial_energy = None;
    }

    /// Removes the centre-of-mass velocity so the system is in its rest frame.
    pub fn zero_com_velocity(&mut self) {
        calibration::zero_com_velocity(&mut self.bodies, self.total_mass);
    }

    /// Recenters the system so that the centre of mass is at the origin.
    ///
    /// The trail buffer is translated by the same vector so stored positions
    /// remain visually consistent.
    pub fn recenter_com(&mut self) {
        if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
            calibration::apply_body_shift(&mut self.bodies, dx, dy);
            self.trail_buf.translate(-dx as f32, -dy as f32);
        }
    }

    /// All body names (parallel to `bodies()`).
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Rename body `idx`. Silently ignores out-of-range indices.
    pub fn set_name(&mut self, idx: usize, name: String) {
        if let Some(slot) = self.names.get_mut(idx) {
            *slot = name;
        }
    }

    /// Replaces the entire set of bodies in the simulation.
    ///
    /// All previous state is cleared, the trail buffer is reset, and the
    /// system is normalised to its COM rest frame.
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.scratch_acc.clear();
        self.names.clear();
        self.total_mass = 0.0;

        for mut b in bodies {
            b.sync_physical_properties();
            b.update_luminosity(MASS_TO_SOLAR, RADIUS_TO_SOLAR, L_SUN);
            self.total_mass += b.mass;
            self.names.push(auto_name(b.material, &self.names));
            self.bodies.push(b);
        }

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
        self.rel_energy_error = 0.0;
        self.steps = 0;
        self.t = 0.0;
        self.last_potential = 0.0;
        self.last_kinetic = 0.0;
        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();
        // Reset controllers: the system topology changed; slew history from
        // a previous run is stale.  dt/theta settings are preserved since
        // load_bodies doesn't alter those parameters.
        self.dt_ctrl.reset();
        self.theta_ctrl.set(self.theta);

        self.zero_com_velocity();
        self.recenter_com();

        let (r_min, softening_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }

    /// Removes a body from the simulation.
    ///
    /// Uses `swap_remove` for O(1) removal.  The trail buffer is reset
    /// because body indices change.
    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            let removed = self.bodies.swap_remove(index);
            self.total_mass -= removed.mass;
            if index < self.names.len() {
                self.names.swap_remove(index);
            }

            let n = self.bodies.len();
            let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
            self.trail_buf.reset(n, cap);
            self.trail_buf.update_colors(&self.bodies);

            self.initial_energy = None;
            self.rel_energy_error = 0.0;
        }
    }

    /// Updates a body in-place, recomputing derived physical properties.
    ///
    /// If the body colour changes, the trail colour buffer is re-uploaded on
    /// the next render frame.
    pub fn update_body(&mut self, index: usize, mut body: Body) {
        if let Some(slot) = self.bodies.get_mut(index) {
            let mass_changed = (slot.mass - body.mass).abs() > 1e-15;

            body.sync_physical_properties();

            if mass_changed {
                self.total_mass += body.mass - slot.mass;
            }

            body.update_luminosity(MASS_TO_SOLAR, RADIUS_TO_SOLAR, L_SUN);
            *slot = body;

            if mass_changed {
                self.initial_energy = None;
                self.rel_energy_error = 0.0;
            }

            self.trail_buf.update_colors(&self.bodies);
        }
    }
}

impl System {
    // ── Timestep guidance ────────────────────────────────────────────────────

    /// Computes a physics-justified recommended timestep from the current system
    /// state, using two complementary N-body criteria:
    ///
    /// 1. **Power et al. (2003) acceleration criterion:**
    ///    `dt_acc = η · √(ε_min / a_max)`
    ///    Ensures no body moves more than ~η softening lengths per step.
    ///
    /// 2. **Aarseth jerk criterion:**
    ///    `dt_jerk = η · √(a_max / j_max)`
    ///    Limits the fractional change in acceleration per step.
    ///    Only used after the first integration step (requires computed jerk).
    ///
    /// Returns the **minimum** of both estimates, clamped to `[1e-9, 1e6]`.
    /// Returns `None` before the first force evaluation or when no bodies exist.
    ///
    /// `η = 0.05` is a conservative default suitable for publication-quality
    /// runs.  For exploratory work η = 0.1 is commonly used.
    ///
    /// # References
    /// - Power et al. (2003). MNRAS 338, 14–34. §3.
    /// - Aarseth, S. J. (2003). *Gravitational N-Body Simulations*. Cambridge. §2.
    fn compute_recommended_dt(&self) -> Option<f64> {
        if self.bodies.is_empty() || self.last_diag.max_acc <= 1e-30 {
            return None;
        }

        // `body.softening` already incorporates `softening_scale`
        // (applied during add_body / load_bodies / set_softening_scale).
        let eps_min = self
            .bodies
            .iter()
            .map(|b| b.softening)
            .fold(f64::MAX, f64::min);

        if eps_min >= f64::MAX || eps_min <= 0.0 {
            return None;
        }

        const ETA: f64 = 0.05;

        let dt_acc = ETA * (eps_min / self.last_diag.max_acc).sqrt();

        let dt_jerk = if self.last_diag.jerk > 1e-30 {
            ETA * (self.last_diag.max_acc / self.last_diag.jerk).sqrt()
        } else {
            f64::MAX
        };

        Some(dt_acc.min(dt_jerk).clamp(1e-9, 1e6))
    }

    /// Returns an immutable slice of all bodies in the simulation.
    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    pub fn dt(&self) -> f64 {
        self.current_dt
    }

    /// Total simulated time elapsed.
    pub fn t(&self) -> f64 {
        self.t
    }

    /// Number of integration steps completed.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// Returns a shared reference to the GPU-ready trail ring buffer.
    pub fn trail_buf(&self) -> &TrailBuffer {
        &self.trail_buf
    }

    /// Returns a mutable reference to the GPU-ready trail ring buffer.
    ///
    /// Required by the trail renderer to drain dirty flags each frame.
    pub fn trail_buf_mut(&mut self) -> &mut TrailBuffer {
        &mut self.trail_buf
    }

    /// Returns the total mass of the system.
    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    /// Sets the gravitational scaling factor.
    pub fn set_g_factor(&mut self, g: f64) {
        self.g_factor = g.max(0.0);
    }

    /// Returns the current gravitational scaling factor.
    pub fn g_factor(&self) -> f64 {
        self.g_factor
    }

    pub fn set_dt(&mut self, dt: f64) {
        self.user_dt = dt;
        self.current_dt = dt;
        // Discard the controller's slew history so the next step starts
        // fresh from the new user-requested value rather than sleweing from
        // the previous adapted value.
        self.dt_ctrl.reset();
    }

    /// Returns the active integrator.
    pub fn integrator(&self) -> Integrator {
        self.integrator
    }

    /// Switches the integration algorithm.  Takes effect on the next [`step`].
    ///
    /// When switching to [`Integrator::WisdomHolman`] on a system that does not
    /// satisfy the dominance criterion, [`step`] will fall back to Yoshida4
    /// automatically.  Check [`is_wh_suitable`] beforehand if you need to know.
    pub fn set_integrator(&mut self, i: Integrator) {
        self.integrator = i;
    }

    /// Returns `true` if the system satisfies the Wisdom–Holman dominance
    /// criterion and the integrator is safe to use.
    ///
    /// The criterion is:
    ///
    /// 1. At least two bodies are present.
    /// 2. `bodies[0]` is the most massive body.
    /// 3. `bodies[0].mass ≥ WH_DOMINANCE_RATIO × Σ mᵢ (i > 0)`.
    ///
    /// The threshold [`WH_DOMINANCE_RATIO`] is 10.  This is intentionally
    /// conservative: WH is formally valid for any ε = m_perturber/M_central < 1,
    /// but accuracy degrades rapidly once the ratio drops below ~10.  For
    /// the Solar System (Jupiter/Sun ≈ 1/1047) the criterion is met with large
    /// margin; equal-mass or figure-8 systems fail immediately.
    pub fn is_wh_suitable(&self) -> bool {
        if self.bodies.len() < 2 {
            return false;
        }
        let m0 = self.bodies[0].mass;
        let m_rest: f64 = self.bodies[1..].iter().map(|b| b.mass).sum();
        // bodies[0] must also be the heaviest individual body.
        let max_other = self.bodies[1..]
            .iter()
            .map(|b| b.mass)
            .fold(0.0_f64, f64::max);
        m0 >= max_other && m0 >= WH_DOMINANCE_RATIO * m_rest
    }

    /// Returns the current Barnes–Hut opening angle θ.
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Sets the Barnes–Hut opening angle θ (clamped to [0.05, 1.5]).
    ///
    /// Smaller θ → more accurate (approaches O(N²) as θ → 0).
    /// Larger θ → faster but less accurate.
    ///
    /// Also syncs the adaptive controller's internal state so that if
    /// adaptation is later enabled, it starts from the user-set value.
    pub fn set_theta(&mut self, theta: f64) {
        let t = theta.clamp(0.05, 1.5);
        self.theta = t;
        self.theta_ctrl.set(t);
    }

    /// Returns the user-requested timestep.
    ///
    /// When adaptive dt is disabled this equals `dt()`.  When enabled, `dt()`
    /// may differ as the controller modulates it step-by-step.
    pub fn user_dt(&self) -> f64 {
        self.user_dt
    }

    /// Set the timestep management policy.
    ///
    /// # Scientific implications
    ///
    /// Setting [`DtMode::Adaptive`] breaks the symplectic structure of the
    /// integrator and may produce secular energy drift.  See [`DtMode`] for
    /// the full rationale.  **Use [`DtMode::Fixed`] for any run whose results
    /// will be analysed or cited.**
    ///
    /// Switching to [`DtMode::Fixed`] immediately restores `current_dt` to
    /// `user_dt` and resets the controller's slew history so no adapted state
    /// bleeds into the fixed-dt run.
    pub fn set_dt_mode(&mut self, mode: DtMode) {
        self.dt_mode = mode;
        if mode == DtMode::Fixed {
            self.current_dt = self.user_dt;
            self.dt_ctrl.reset();
        }
    }

    /// Returns the active timestep management policy.
    pub fn dt_mode(&self) -> DtMode {
        self.dt_mode
    }

    /// Enable or disable the adaptive Barnes–Hut θ controller.
    ///
    /// **Disabled by default.**  When enabled, θ is adjusted each step based
    /// on the BH force-truncation error proxy, targeting the controller's
    /// configured error tolerance.  When disabled, θ is fixed at the
    /// user-set value.
    ///
    /// Note: varying θ does **not** break symplecticity, but it does change
    /// the force accuracy budget per step, making quantitative error analysis
    /// harder.  For reproducible accuracy in published runs, keep θ fixed.
    ///
    /// Disabling re-syncs the controller's internal state to the current θ
    /// so that re-enabling starts from a consistent baseline.
    pub fn set_adaptive_theta(&mut self, enabled: bool) {
        self.adaptive_theta = enabled;
        if !enabled {
            self.theta_ctrl.set(self.theta);
        }
    }

    /// Returns `true` if the adaptive θ controller is active.
    pub fn adaptive_theta_enabled(&self) -> bool {
        self.adaptive_theta
    }

    /// Returns the current global softening scale factor.
    pub fn softening_scale(&self) -> f64 {
        self.softening_scale
    }

    /// Sets a global Plummer softening scale applied on top of the
    /// per-body mass-proportional default (`ε = ε_default · scale`).
    ///
    /// Also rescales all existing body softenings immediately.
    pub fn set_softening_scale(&mut self, scale: f64) {
        use crate::core::body::default_softening;
        self.softening_scale = scale.max(0.0);
        for b in &mut self.bodies {
            b.softening = default_softening(b.mass) * self.softening_scale;
        }
    }

    pub fn trail_every(&self) -> usize {
        self.trail_every
    }

    pub fn set_trail_every(&mut self, n: usize) {
        self.trail_every = n.max(1);
    }

    /// Records the current body positions into the trail ring buffer.
    ///
    /// Call this **once per rendered frame** (not per physics step) so the
    /// trail density is proportional to the amount of simulated time per
    /// frame rather than to a fixed physics-step count.
    pub fn push_trail(&mut self) {
        self.trail_buf.push(&self.bodies);
    }

    /// Returns a reference to the Barnes–Hut engine.
    pub fn engine(&self) -> &BarnesHutEngine {
        &self.engine
    }

    /// Returns diagnostic metrics for the current simulation state.
    pub fn metrics(&self) -> Metrics {
        let kinetic = self.last_kinetic;
        let potential = self.last_potential;
        let total = total_energy(kinetic, potential);

        let lz = angular_momentum_z(&self.bodies);
        let (com_x, com_y, com_vx, com_vy) = center_of_mass_state(&self.bodies);

        Metrics {
            kinetic,
            potential,
            total_energy: total,
            rel_energy_error: self.rel_energy_error,

            angular_momentum_z: lz,
            rel_angular_momentum_error: self.rel_angular_momentum_error,
            abs_angular_momentum_error: self.abs_angular_momentum_error,

            com_x,
            com_y,
            com_vx,
            com_vy,

            t: self.t,
            steps: self.steps,

            integrator: self.integrator,
            g_factor: self.g_factor,
            theta: self.theta,
            dt: self.current_dt,
            user_dt: self.user_dt,
            dt_mode: self.dt_mode,
            adaptive_theta: self.adaptive_theta,

            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,

            r_min: self.r_min,
            softening_max: self.softening_max,

            recommended_dt: self.compute_recommended_dt(),
        }
    }

    /// Returns accelerations computed during the last integration step.
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }

    // ── Orbital elements ─────────────────────────────────────────────────────

    /// Recomputes osculating orbital elements for all bodies and caches the result.
    ///
    /// This is O(N²) and should be called **once per rendered frame**, not every
    /// physics step. The result is available via [`orbital_elements`].
    pub fn update_orbital_elements(&mut self) {
        self.orbital_cache = orbital::compute_all(&self.bodies, self.g_factor);
    }

    /// Returns the cached osculating orbital elements (one slot per body).
    ///
    /// Call [`update_orbital_elements`] first to get fresh values.
    pub fn orbital_elements(&self) -> &[Option<OrbitalElements>] {
        &self.orbital_cache
    }

    // ── Snapshot (save / load) ───────────────────────────────────────────────

    /// Capture the minimal state required for deterministic reproduction.
    pub fn to_snapshot(&self) -> crate::core::snapshot::SimSnapshot {
        use crate::core::snapshot::{BodyRecord, SimSnapshot};
        SimSnapshot {
            save_id: 0,
            t: self.t,
            steps: self.steps,
            dt: self.current_dt,
            theta: self.theta,
            softening_scale: self.softening_scale,
            g_factor: self.g_factor,
            integrator: self.integrator,
            trail_every: self.trail_every,
            sim_name: String::new(),
            seed: 0,
            trail: None,
            bodies: self.bodies.iter().map(BodyRecord::from_body).collect(),
            names: self.names.clone(),
        }
    }

    /// Replace the current simulation state with a saved snapshot.
    ///
    /// The trail buffer is cleared (it is cosmetic and cannot be restored).
    /// Energy / angular-momentum references are reset so the first post-load
    /// step establishes new baselines.
    pub fn restore_from_snapshot(&mut self, snap: &crate::core::snapshot::SimSnapshot) {
        let bodies: Vec<Body> = snap.bodies.iter().map(|r| r.into_body()).collect();
        self.names = if snap.names.len() == bodies.len() {
            snap.names.clone()
        } else {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for b in &bodies {
                acc.push(auto_name(b.material, &acc));
            }
            acc
        };
        let n = bodies.len();

        self.bodies = bodies;
        self.total_mass = self.bodies.iter().map(|b| b.mass).sum();
        self.scratch_acc.clear();

        // Rebuild trail buffer — restore saved trail if dimensions match,
        // otherwise start empty.
        let cap =
            crate::core::trail_buffer::adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);
        if let Some(trail_snap) = &snap.trail {
            if trail_snap.n_bodies == n as u32
                && trail_snap.positions.len()
                    == (trail_snap.n_bodies * trail_snap.capacity) as usize
            {
                self.trail_buf.restore_from_snapshot(trail_snap);
            }
        }

        self.t = snap.t;
        self.steps = snap.steps;
        self.current_dt = snap.dt;
        self.user_dt = snap.dt;
        self.theta = snap.theta;
        self.softening_scale = snap.softening_scale;
        self.g_factor = snap.g_factor;
        self.integrator = snap.integrator;
        self.trail_every = snap.trail_every.max(1);

        self.initial_energy = None;
        self.initial_angular_momentum = None;
        self.rel_energy_error = 0.0;
        self.rel_angular_momentum_error = 0.0;
        self.abs_angular_momentum_error = 0.0;
        self.last_kinetic = 0.0;
        self.last_potential = 0.0;
        self.diagnostics = crate::core::diagnostics::DiagnosticsComputer::new();
        self.last_diag = crate::core::diagnostics::SimulationDiagnostics::default();
        self.orbital_cache.clear();
        // Reset controllers: discard stale slew history from the previous run.
        self.dt_ctrl.reset();
        self.theta_ctrl.set(snap.theta);

        let (r_min, softening_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }
}

// ── End-to-end integration tests ─────────────────────────────────────────────
//
// These tests verify that the full simulation pipeline — force evaluation,
// integrator, energy tracking — correctly conserves the Hamiltonian over many
// orbital periods.  They test the *integrated system*, not individual
// primitives.
//
// Physical scenario
// ─────────────────
// Two equal-mass bodies in a circular orbit about their common centre of mass.
//
//   G = 1 (simulation units), M₁ = M₂ = 1
//   Positions: (−1, 0) and (+1, 0) — separation d = 2, orbital radius r = 1
//   Velocities: (0, −0.5) and (0, +0.5) — counter-clockwise orbit
//
// Derivation of initial conditions:
//   Circular orbit requires centripetal = gravitational acceleration:
//     v²/r = G·M_partner/d²  →  v²/1 = 1·1/4  →  v = 0.5   ✓
//   Orbital period:
//     T = 2πr/v = 2π·1/0.5 = 4π ≈ 12.566
//   Centre-of-mass velocity: (m₁·v₁ + m₂·v₂)/(m₁+m₂) = (−0.5+0.5)/2 = 0  ✓
//   Angular momentum: Lz = m(x₁vy₁ − y₁vx₁) + m(x₂vy₂ − y₂vx₂)
//                       = 1·(−1·(−0.5)) + 1·(1·0.5) = 1.0 > 0 (CCW)  ✓
//
// Energy-error measurement
// ────────────────────────
// `System::update_energy_tracking` sets E₀ on the first step, then tracks
//   δE/E₀ = (Eₙ − E₀) / |E₀|
// continuously.  `metrics().rel_energy_error` returns the instantaneous value.
// We record the *maximum* over all steps: for symplectic integrators the error
// oscillates rather than drifts, so the final sample can underestimate the
// true peak.
//
// Tolerance derivation
// ────────────────────
// With dt = 0.01 and T = 4π, the ratio dt/T ≈ 7.96 × 10⁻⁴.
//
// Velocity Verlet (2nd order):
//   Amplitude of energy oscillation ~ (dt/T)² ≈ 6.3 × 10⁻⁷
//   Tolerance 1 × 10⁻⁴ gives a factor-of-160 safety margin.
//
// Yoshida 4th-order:
//   Amplitude ~ (dt/T)⁴ ≈ 4 × 10⁻¹³
//   Tolerance 1 × 10⁻⁷ gives ≥ 10⁶× safety margin (also validates that
//   Yoshida4 is *observably* more accurate than VV for the same dt).
//
// Both tolerances accommodate:
//   - Plummer softening (ε = 0.02 per unit-mass body): modifies potential
//     smoothly but the modified system still has a symplectic structure
//   - Floating-point rounding: accumulated over ~10⁵ steps, remains ≪ 1 × 10⁻⁷
//   - Bounded oscillatory transients from the first few orbits
#[cfg(test)]
mod integration_tests {
    use super::*;
    use crate::core::materials::Material;
    use crate::physics::integrator::Integrator;

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Constructs the two-body circular-orbit scenario described above.
    ///
    /// θ = 0.5: exact O(N²) forces are used regardless (N = 2 < EXACT_THRESHOLD
    /// = 64), so θ has no effect on accuracy here.
    /// `DtMode::Fixed` is the default — integration is fully symplectic.
    fn two_body_circular_system(integrator: Integrator, dt: f64) -> System {
        let bodies = vec![
            Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
            Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
        ];
        let mut sys = System::new(bodies, 0.5, dt, 10, 1);
        sys.set_integrator(integrator);
        sys
    }

    /// Runs `sys` for `n_periods` orbital periods and returns the maximum
    /// relative energy error |δE/E₀| observed over all steps.
    ///
    /// Tracking the maximum is essential for symplectic integrators: their
    /// energy error oscillates, so sampling only the final step would miss the
    /// true peak and could give a falsely optimistic result.
    fn max_rel_energy_error(sys: &mut System, n_periods: u64, dt: f64) -> f64 {
        // T = 4π (derived above)
        const PERIOD: f64 = 4.0 * std::f64::consts::PI;
        let total_steps = (n_periods as f64 * PERIOD / dt).ceil() as u64;

        let mut peak: f64 = 0.0;
        for _ in 0..total_steps {
            sys.step();
            peak = peak.max(sys.metrics().rel_energy_error.abs());
        }
        peak
    }

    // ── Tests ─────────────────────────────────────────────────────────────────

    /// Velocity Verlet energy conservation over 100 circular orbits.
    ///
    /// VV is a 2nd-order symplectic integrator: energy error is bounded and
    /// oscillatory with amplitude O((dt/T)²).
    ///
    /// With dt = 0.01, T = 4π:
    ///   (dt/T)² ≈ 6.3 × 10⁻⁷ → tolerance 10⁻⁴ (factor-of-160 safety margin)
    #[test]
    fn energy_conservation_velocity_verlet() {
        const DT: f64 = 0.01;
        const N_PERIODS: u64 = 100;
        const TOLERANCE: f64 = 1e-4;

        let mut sys = two_body_circular_system(Integrator::VelocityVerlet, DT);
        let peak_err = max_rel_energy_error(&mut sys, N_PERIODS, DT);

        assert!(
            peak_err < TOLERANCE,
            "VelocityVerlet: peak |δE/E₀| = {:.3e} exceeds {:.0e} \
             after {} periods (dt = {}, T = 4π ≈ 12.566)",
            peak_err,
            TOLERANCE,
            N_PERIODS,
            DT,
        );
    }

    /// Yoshida 4th-order energy conservation over 100 circular orbits.
    ///
    /// Yoshida4 is a 4th-order symplectic integrator: energy error amplitude
    /// is O((dt/T)⁴), far smaller than VV for the same dt.
    ///
    /// With dt = 0.01, T = 4π:
    ///   (dt/T)⁴ ≈ 4 × 10⁻¹³ → tolerance 10⁻⁷
    ///
    /// The tighter tolerance vs. VV validates that the higher-order method
    /// delivers its theoretical accuracy advantage in the full pipeline.
    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect raw peak errors"]
    fn print_peak_errors_diagnostic() {
        for &(label, integrator, dt) in &[
            ("VV    dt=0.01 ", Integrator::VelocityVerlet, 0.01_f64),
            ("VV    dt=0.001", Integrator::VelocityVerlet, 0.001_f64),
            ("Y4    dt=0.01 ", Integrator::Yoshida4, 0.01_f64),
            ("Y4    dt=0.001", Integrator::Yoshida4, 0.001_f64),
        ] {
            let mut sys = two_body_circular_system(integrator, dt);
            let peak = max_rel_energy_error(&mut sys, 10, dt);
            println!("{label}  peak |δE/E₀| = {peak:.3e}");
        }
    }

    #[test]
    fn energy_conservation_yoshida4() {
        const DT: f64 = 0.01;
        const N_PERIODS: u64 = 100;
        const TOLERANCE: f64 = 1e-7;

        let mut sys = two_body_circular_system(Integrator::Yoshida4, DT);
        let peak_err = max_rel_energy_error(&mut sys, N_PERIODS, DT);

        assert!(
            peak_err < TOLERANCE,
            "Yoshida4: peak |δE/E₀| = {:.3e} exceeds {:.0e} \
             after {} periods (dt = {}, T = 4π ≈ 12.566)",
            peak_err,
            TOLERANCE,
            N_PERIODS,
            DT,
        );
    }
}

// ── Guard: Wisdom–Holman dominance criterion ─────────────────────────────────
//
// Tests that `is_wh_suitable` correctly identifies hierarchical and
// non-hierarchical systems, and that `step_wisdom_holman` falls back to
// Yoshida4 (producing physically correct output) when the criterion is not met.
#[cfg(test)]
mod wh_guard_tests {
    use super::*;
    use crate::core::materials::Material;
    use crate::physics::integrator::Integrator;

    /// A Sun-like body at the origin plus one lightweight planet.
    /// Sun mass = 1000, planet mass = 1 → ratio = 1000 ≥ 10.
    fn hierarchical_system() -> System {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 1000.0, Material::Star),
            Body::new(10.0, 0.0, 0.0, 10.0, 1.0, Material::Rocky),
        ];
        let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
        sys.set_integrator(Integrator::WisdomHolman);
        sys
    }

    /// Two equal-mass bodies — no dominant central body.
    fn equal_mass_system() -> System {
        let bodies = vec![
            Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
            Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
        ];
        let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
        sys.set_integrator(Integrator::WisdomHolman);
        sys
    }

    /// Three equal-mass bodies (figure-8 topology) — clearly non-hierarchical.
    fn three_equal_mass_system() -> System {
        let bodies = vec![
            Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
            Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
            Body::new(0.0, 1.0, 0.5, 0.0, 1.0, Material::Rocky),
        ];
        let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
        sys.set_integrator(Integrator::WisdomHolman);
        sys
    }

    /// Central body at the boundary: mass = exactly 10 × rest → ratio = 10 → suitable.
    fn boundary_system_just_above() -> System {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 10.0, Material::Star),
            Body::new(10.0, 0.0, 0.0, 1.0, 1.0, Material::Rocky),
        ];
        System::new(bodies, 0.5, 0.01, 10, 1)
    }

    /// Central body just below: mass = 9.9 × rest → ratio < 10 → not suitable.
    fn boundary_system_just_below() -> System {
        let bodies = vec![
            Body::new(0.0, 0.0, 0.0, 0.0, 9.9, Material::Star),
            Body::new(10.0, 0.0, 0.0, 1.0, 1.0, Material::Rocky),
        ];
        System::new(bodies, 0.5, 0.01, 10, 1)
    }

    #[test]
    fn hierarchical_system_is_suitable() {
        assert!(hierarchical_system().is_wh_suitable());
    }

    #[test]
    fn equal_mass_system_is_not_suitable() {
        assert!(!equal_mass_system().is_wh_suitable());
    }

    #[test]
    fn three_equal_mass_is_not_suitable() {
        assert!(!three_equal_mass_system().is_wh_suitable());
    }

    #[test]
    fn boundary_at_exactly_10x_is_suitable() {
        assert!(boundary_system_just_above().is_wh_suitable());
    }

    #[test]
    fn boundary_below_10x_is_not_suitable() {
        assert!(!boundary_system_just_below().is_wh_suitable());
    }

    #[test]
    fn single_body_is_not_suitable() {
        let bodies = vec![Body::new(0.0, 0.0, 0.0, 0.0, 1.0, Material::Rocky)];
        let sys = System::new(bodies, 0.5, 0.01, 10, 1);
        assert!(!sys.is_wh_suitable());
    }

    /// When WH is selected on a non-hierarchical system, `step` must not panic
    /// and must produce finite positions (the fallback to Yoshida4 is active).
    #[test]
    fn wh_on_non_hierarchical_does_not_panic_and_stays_finite() {
        let mut sys = equal_mass_system();
        for _ in 0..100 {
            sys.step();
        }
        for b in sys.bodies() {
            assert!(b.x.is_finite() && b.y.is_finite(), "body left finite domain");
            assert!(b.vx.is_finite() && b.vy.is_finite(), "velocity left finite domain");
        }
    }

    /// The fallback must conserve energy as well as plain Yoshida4 would.
    /// We run both WH-on-equal-mass (falls back to Y4) and explicit Y4 for
    /// 100 steps and confirm the energy errors are identical.
    #[test]
    fn wh_fallback_energy_matches_yoshida4_directly() {
        let mut sys_wh = equal_mass_system(); // WH selected → Y4 fallback
        let mut sys_y4 = equal_mass_system();
        sys_y4.set_integrator(Integrator::Yoshida4);

        for _ in 0..100 {
            sys_wh.step();
            sys_y4.step();
        }

        let err_wh = sys_wh.metrics().rel_energy_error.abs();
        let err_y4 = sys_y4.metrics().rel_energy_error.abs();
        // Both paths call the same step_yoshida4 kernel; errors should be identical.
        assert!(
            (err_wh - err_y4).abs() < 1e-15,
            "WH fallback energy error {err_wh:.3e} ≠ direct Y4 {err_y4:.3e}"
        );
    }
}

// ── Benchmark: Kepler vs. analytical solution ────────────────────────────────
//
// Validates that the simulated position of a two-body elliptical orbit matches
// the exact Keplerian solution at a given time t.
//
// Setup
// -----
// Two equal-mass bodies (m₁ = m₂ = 1), zero Plummer softening so the kernel
// reduces to the exact Newtonian 1/r² force.
//
//   - Relative semi-major axis  a = 2.0
//   - Eccentricity              e = 0.5
//   - μ = G·(m₁+m₂) = 2.0  →  T = 2π√(a³/μ) = 4π ≈ 12.566
//   - Periapsis separation      r_peri = a·(1−e) = 1.0
//   - Relative speed at peri    v_peri = √(μ·(1+e)/(a·(1−e))) = √3
//
// By symmetry the centre of mass is at the origin throughout, so the
// absolute position of each body is ±r_rel/2 where r_rel = r₂ − r₁.
//
// Comparison
// ----------
// After N steps (total time t = N·dt), the simulated relative position is
// compared against the Kepler analytical prediction at exactly that time.
// Using the simulated time t = N·dt (instead of the nominal period T) removes
// any period-discretisation error: the only remaining error is the integrator's
// local truncation error accumulated over N steps.
//
// Expected accuracy (two-body elliptical orbit, a = 2, e = 0.5, T = 4π)
// -----------------------------------------------------------------------
//   VV   (2nd order, dt = 0.01): O(dt²/T) × 2π ≈ 5 × 10⁻⁵  →  tol 10⁻³
//   Y4   (4th order, dt = 0.01): O(dt⁴/T) × 2π ≈ 5 × 10⁻⁹  →  tol 10⁻⁶
//
// References
// ----------
//   - Yoshida (1990). Phys. Lett. A 150, 262–268.
//   - Forest & Ruth (1990). Nucl. Instrum. Methods Phys. Res. A 290, 395–400.
#[cfg(test)]
mod benchmark_kepler {
    use super::*;
    use crate::core::materials::Material;
    use crate::physics::integrator::Integrator;

    // ── Kepler helpers ────────────────────────────────────────────────────────

    /// Solve Kepler's equation  M = E − e·sin(E)  for the eccentric anomaly E
    /// using Newton–Raphson iteration (converges in ≤ 10 iterations for e < 1).
    fn solve_kepler(mean_anomaly: f64, e: f64) -> f64 {
        let mut ea = mean_anomaly;
        for _ in 0..60 {
            let d = (mean_anomaly - ea + e * ea.sin()) / (1.0 - e * ea.cos());
            ea += d;
            if d.abs() < 1e-14 {
                break;
            }
        }
        ea
    }

    /// Analytical relative position **r = r₂ − r₁** at time `t` for a
    /// two-body Keplerian orbit starting at periapsis.
    ///
    /// - `mu`:  gravitational parameter G·(m₁ + m₂)
    /// - `a`:   semi-major axis of the relative orbit
    /// - `e`:   eccentricity (0 ≤ e < 1)
    ///
    /// Returns `(x, y)` in the orbital plane with the periapsis along +x.
    fn kepler_relative_pos(t: f64, mu: f64, a: f64, e: f64) -> (f64, f64) {
        let n = (mu / a.powi(3)).sqrt(); // mean motion ω = √(μ/a³)
        let ea = solve_kepler(n * t, e); // eccentric anomaly
        let x = a * (ea.cos() - e);
        let y = a * (1.0 - e * e).sqrt() * ea.sin();
        (x, y)
    }

    /// Two equal-mass bodies (m = 1, ε = 0) placed at periapsis of an
    /// elliptical orbit with a_rel = 2, e = 0.5.
    fn kepler_two_body(integrator: Integrator, dt: f64) -> System {
        const A: f64 = 2.0;
        const E: f64 = 0.5;
        const MU: f64 = 2.0; // G·(1 + 1)
        let r_peri = A * (1.0 - E); // = 1.0
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt(); // = √3

        let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
        b1.softening = 0.0; // exact Newtonian kernel
        let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], 0.5, dt, 10, 1);
        sys.set_integrator(integrator);
        sys
    }

    /// Run `n_steps` and return the Euclidean error between the simulated
    /// relative position and the Kepler prediction at time `t = n_steps · dt`.
    fn kepler_position_error(integrator: Integrator, dt: f64, n_steps: u64) -> f64 {
        const A: f64 = 2.0;
        const E: f64 = 0.5;
        const MU: f64 = 2.0;

        let mut sys = kepler_two_body(integrator, dt);
        for _ in 0..n_steps {
            sys.step();
        }
        let t = n_steps as f64 * dt;

        let bodies = sys.bodies();
        let rx = bodies[1].x - bodies[0].x;
        let ry = bodies[1].y - bodies[0].y;

        let (ex, ey) = kepler_relative_pos(t, MU, A, E);
        ((rx - ex).powi(2) + (ry - ey).powi(2)).sqrt()
    }

    // ── Tests ─────────────────────────────────────────────────────────────────

    /// VV relative position matches the Kepler analytical orbit to within 10⁻²
    /// after ≈ 1 orbital period (T = 4π ≈ 12.566, dt = 0.01, 1257 steps).
    ///
    /// VV is 2nd-order symplectic: phase error accumulates as O(dt²) per orbit.
    /// Measured error at dt = 0.01, e = 0.5: ≈ 2.2 × 10⁻³  →  tol 10⁻².
    /// The tighter dt = 0.001 result (diagnostic) confirms O(dt²) scaling.
    #[test]
    fn kepler_position_accuracy_velocity_verlet() {
        const DT: f64 = 0.01;
        const N: u64 = 1257; // round(4π / 0.01) ≈ 1256.6 → 1257 steps
        const TOL: f64 = 1e-2;

        let err = kepler_position_error(Integrator::VelocityVerlet, DT, N);
        assert!(
            err < TOL,
            "VV Kepler: |Δr| = {:.3e} exceeds {:.0e} \
             after {N} steps (dt={DT}, t≈{:.4}, T=4π≈12.566)",
            err,
            TOL,
            N as f64 * DT,
        );
    }

    /// Yoshida4 relative position matches the Kepler analytical orbit to within
    /// 10⁻⁶ after ≈ 1 orbital period (same setup as VV benchmark).
    ///
    /// Y4 is 4th-order symplectic: phase error O(dt⁴/T) per orbit.
    /// With dt = 0.01, T = 4π:  expected error ≈ 5 × 10⁻⁹  →  tol 10⁻⁶.
    ///
    /// The gap between VV (10⁻³) and Y4 (10⁻⁶) validates the order advantage
    /// of the higher-order integrator on a non-trivial, non-circular orbit.
    #[test]
    fn kepler_position_accuracy_yoshida4() {
        const DT: f64 = 0.01;
        const N: u64 = 1257;
        const TOL: f64 = 1e-6;

        let err = kepler_position_error(Integrator::Yoshida4, DT, N);
        assert!(
            err < TOL,
            "Y4 Kepler: |Δr| = {:.3e} exceeds {:.0e} \
             after {N} steps (dt={DT}, t≈{:.4}, T=4π≈12.566)",
            err,
            TOL,
            N as f64 * DT,
        );
    }

    /// Diagnostic: print actual position errors for both integrators and
    /// several dt values.  Run with `cargo test -- --ignored` to inspect.
    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect raw Kepler errors"]
    fn print_kepler_errors_diagnostic() {
        for &(label, integrator, dt, n) in &[
            ("VV  dt=0.01  ", Integrator::VelocityVerlet, 0.01_f64, 1257u64),
            ("VV  dt=0.001 ", Integrator::VelocityVerlet, 0.001_f64, 12567u64),
            ("Y4  dt=0.01  ", Integrator::Yoshida4, 0.01_f64, 1257u64),
            ("Y4  dt=0.001 ", Integrator::Yoshida4, 0.001_f64, 12567u64),
        ] {
            let err = kepler_position_error(integrator, dt, n);
            println!("{label}  |Δr| = {err:.3e}");
        }
    }
}

// ── Benchmark: figure-8 three-body orbit closure ─────────────────────────────
//
// Validates that the Chenciner–Montgomery figure-8 choreography closes after
// one period T ≈ 6.3259 (G = 1, m = 1 units).
//
// Setup
// -----
// Three equal masses (m = 1) with zero Plummer softening, initialised at the
// published Chenciner & Montgomery (2000) initial conditions.  The centre of
// mass is exactly at the origin with zero total momentum.
//
// After STEPS = round(T/dt) integration steps the simulated positions of all
// three bodies are compared against their initial positions.  Orbit closure is
// the defining property of a choreography: if closure fails, either the orbit
// is non-periodic (wrong ICs), the integrator has too large an error, or the
// force kernel has been altered.
//
// Error budget (Yoshida4, dt = 0.001, T ≈ 6.3259, 6326 steps)
// -------------------------------------------------------------
//   Timing discretisation:  |t_actual − T| = |6.326 − 6.3259| ≈ 8.6 × 10⁻⁵
//     → positional floor   ≈ 8.6 × 10⁻⁵ × v_max ≈ 8.6 × 10⁻⁵
//   Y4 integration error:   O(dt³ · T) ≈ 6 × 10⁻⁹   (negligible)
//   Tolerance:  10⁻³  (factor-of-12 safety over the timing floor)
//
// References
// ----------
//   - Chenciner & Montgomery (2000). Ann. Math. 152, 881–901.
//   - Simó (2002). Celest. Mech. Dyn. Astron. 83, 85–100.
#[cfg(test)]
mod benchmark_figure8 {
    use super::*;
    use crate::core::materials::Material;
    use crate::physics::integrator::Integrator;

    // Initial conditions: Chenciner & Montgomery (2000), G = 1, m = 1.
    // Layout: (x, y, vx, vy)
    const IC: [(f64, f64, f64, f64); 3] = [
        (-0.97000436,  0.24308753,  0.46620369,  0.43236573),
        ( 0.97000436, -0.24308753,  0.46620369,  0.43236573),
        ( 0.0,         0.0,        -0.93240737, -0.86473146),
    ];

    /// Period from Simó (2002) in G = 1, m = 1 units.
    const T: f64 = 6.32591398;

    fn figure8_system(integrator: Integrator, dt: f64) -> System {
        let bodies = IC
            .iter()
            .map(|&(x, y, vx, vy)| {
                let mut b = Body::new(x, y, vx, vy, 1.0, Material::Rocky);
                b.softening = 0.0; // exact Newtonian kernel
                b
            })
            .collect();
        let mut sys = System::new(bodies, 0.5, dt, 10, 1);
        sys.set_integrator(integrator);
        sys
    }

    // ── Tests ─────────────────────────────────────────────────────────────────

    /// Figure-8 orbit closes to within 10⁻³ after one period (Yoshida4).
    ///
    /// STEPS = 6326 ≈ T/0.001 = 6325.9.  Actual simulation time = 6.326,
    /// giving a timing floor of ≈ 8.6 × 10⁻⁵.  Tolerance 10⁻³ leaves a
    /// factor-of-12 margin above that floor.
    #[test]
    fn figure8_orbit_closure_yoshida4() {
        const DT: f64 = 0.001;
        const STEPS: u64 = 6326; // round(T / DT)
        const TOL: f64 = 1e-3;

        let mut sys = figure8_system(Integrator::Yoshida4, DT);

        for _ in 0..STEPS {
            sys.step();
        }

        let bodies = sys.bodies();
        let max_err = IC
            .iter()
            .zip(bodies.iter())
            .map(|(&(x0, y0, _, _), b)| ((b.x - x0).powi(2) + (b.y - y0).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);

        assert!(
            max_err < TOL,
            "Figure-8 (Y4, dt={DT}): max |Δr| = {:.3e} > {:.0e} \
             after {STEPS} steps (t={:.6}, T≈{T:.6})",
            max_err,
            TOL,
            STEPS as f64 * DT,
        );
    }

    /// Diagnostic: print per-body closure errors and peak energy error.
    /// Run with `cargo test -- --ignored` to inspect raw values.
    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect figure-8 closure errors"]
    fn print_figure8_closure_diagnostic() {
        for &(label, integrator, dt, steps) in &[
            ("Y4  dt=0.001 ", Integrator::Yoshida4, 0.001_f64, 6326u64),
            ("Y4  dt=0.0001", Integrator::Yoshida4, 0.0001_f64, 63259u64),
            ("VV  dt=0.001 ", Integrator::VelocityVerlet, 0.001_f64, 6326u64),
        ] {
            let mut sys = figure8_system(integrator, dt);
            for _ in 0..steps {
                sys.step();
            }
            let bodies = sys.bodies();
            let t_actual = steps as f64 * dt;
            println!("{label}  t={t_actual:.6}  T={T:.6}");
            for (i, (&(x0, y0, _, _), b)) in IC.iter().zip(bodies.iter()).enumerate() {
                let err = ((b.x - x0).powi(2) + (b.y - y0).powi(2)).sqrt();
                println!("  body {i}: |Δr| = {err:.3e}");
            }
        }
    }
}
