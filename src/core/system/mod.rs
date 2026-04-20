//! Simulation orchestrator for an N-body gravitational system.
//!
//! [`System`] advances a set of massive bodies under gravity. It owns the
//! integration loop, trail buffer, and all diagnostic state, but delegates
//! force evaluation and integration to pluggable trait objects so that new
//! physics (relativistic corrections, drag, radiation, alternative solvers)
//! can be wired in without touching the core.
//!
//! ## Pluggable contracts
//!
//! | Seam | Trait | Default |
//! |---|---|---|
//! | Force engine | [`ForceModel`] | [`GravityForceModel`] (Barnes-Hut) |
//! | Integrator | [`Integrator`] | [`VelocityVerlet`] |
//! | Extra forces | [`PerturbationForce`] | none |
//!
//! ## Module layout
//!
//! | Module | Responsibility |
//! |---|---|
//! | `mod.rs` | [`System`] struct, constructors |
//! | `step` | `step()` and conservation-law tracking |
//! | `bodies` | body CRUD, names, COM calibration |
//! | `config` | getters/setters (θ, dt, integrator, softening, …) |
//! | `metrics` | [`Metrics`] assembly and recommended-dt |
//! | `orbital` | osculating-element cache |
//! | `snapshot` | save/load via [`SimSnapshot`] |
//! | `perturbations` | non-gravitational force registration |
//! | `helpers` | free functions (naming, closeness, trail count) |

pub(crate) mod bodies;
pub(crate) mod config;
pub(crate) mod helpers;
pub(crate) mod metrics;
pub(crate) mod orbital;
pub(crate) mod perturbations;
pub(crate) mod snapshot;
pub(crate) mod step;
#[cfg(test)]
mod tests;

use crate::core::adaptive::{
    AccelerationStats, DtAdaptationConfig, DtController, DtMode, ThetaController,
};
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::hooks::HookRegistry;
use crate::core::metrics::Metrics;
use crate::domain::body::Body;
use crate::physics::integrator::{
    ForceModel, GravityForceModel, Integrator, IntegratorKind, PerturbationForce, make_integrator,
};
use crate::physics::orbital::OrbitalElements;

/// Central simulation state for an N-body gravitational system.
pub struct System {
    /// Bodies participating in the simulation.
    pub(crate) bodies: Vec<Body>,

    /// Total mass of the system (used for COM recentering).
    pub(crate) total_mass: f64,

    /// Last computed energies.
    pub(crate) last_kinetic: f64,
    pub(crate) last_potential: f64,

    /// Initial total energy (used as reference for relative error).
    pub(crate) initial_energy: Option<f64>,

    /// Relative energy error (diagnostic only).
    pub(crate) rel_energy_error: f64,

    /// Pluggable force model — default: Barnes-Hut gravity.
    ///
    /// Swap via [`System::set_force_model`] to use any [`ForceModel`]
    /// implementation (direct O(N²), GPU kernels, PN corrections, …).
    pub(crate) force_model: Box<dyn ForceModel>,

    /// Scratch buffer for accelerations — reused every step.
    pub(crate) scratch_acc: Vec<(f64, f64)>,

    /// Active integration algorithm (trait object).
    pub(crate) integrator: Box<dyn Integrator>,

    /// Cached osculating orbital elements — one slot per body.
    pub(crate) orbital_cache: Vec<Option<OrbitalElements>>,

    /// Global Plummer softening scale: `ε = ε_default · softening_scale`.
    pub(crate) softening_scale: f64,

    /// Diagnostics subsystem.
    pub(crate) diagnostics: DiagnosticsComputer,
    pub(crate) last_diag: SimulationDiagnostics,

    /// Step counter.
    pub(crate) steps: u64,

    /// Total simulated time elapsed.
    pub(crate) t: f64,

    /// Timestep currently used by the integrator.
    pub(crate) current_dt: f64,

    /// User-requested timestep (baseline for the adaptive controller).
    pub(crate) user_dt: f64,

    /// Timestep management policy.
    pub(crate) dt_mode: DtMode,

    /// Adaptive timestep controller.
    pub(crate) dt_ctrl: DtController,

    /// Adaptive Barnes-Hut opening-angle controller.
    pub(crate) theta_ctrl: ThetaController,

    /// Whether the adaptive θ controller is active.
    pub(crate) adaptive_theta: bool,

    /// Gravitational scaling factor (G multiplier).
    pub(crate) g_factor: f64,

    /// Initial angular momentum (z-component) — conservation baseline.
    pub(crate) initial_angular_momentum: Option<f64>,

    /// Relative angular momentum error.
    pub(crate) rel_angular_momentum_error: f64,

    /// Absolute angular momentum error.
    pub(crate) abs_angular_momentum_error: f64,

    /// Human-readable label for each body, parallel to `bodies`.
    /// Separate because `Body` is `Copy` and cannot own a `String`.
    pub(crate) names: Vec<String>,

    /// Minimum pairwise separation from the most recent step.
    pub(crate) r_min: f64,

    /// Maximum effective pairwise softening length from the most recent step.
    pub(crate) softening_max: f64,

    /// Registered non-gravitational perturbation forces.
    pub(crate) perturbations: Vec<Box<dyn PerturbationForce>>,

    /// Reproducibility seed. Consumed by preset builders and cluster spawners.
    /// Persisted in snapshots so a run can be replayed exactly.
    pub(crate) seed: u64,

    /// Registered observer/command hooks. Dispatched from [`System::step`].
    pub(crate) hooks: HookRegistry,

    /// Set by a [`Command::Stop`](crate::core::hooks::Command::Stop) request.
    /// Headless runners honour this; the GUI may ignore it.
    pub(crate) stop_requested: bool,

    /// Accumulated world-space COM translation since the last call to
    /// [`take_com_shift`](System::take_com_shift). The render-side
    /// [`TrailRecorder`](crate::render::TrailRecorder) reads and clears this
    /// each frame to keep trail positions aligned with the shifted bodies.
    pub(crate) pending_com_shift: (f32, f32),
}

impl System {
    /// Create a simulation with the default Barnes-Hut force model.
    ///
    /// - `theta`:      Barnes-Hut opening angle (accuracy vs speed).
    /// - `dt`:         Fixed timestep.
    /// - `max_depth`:  Maximum quadtree depth.
    /// - `trail_every`: Sampling interval for trail ring-buffer.
    pub fn new(
        bodies: Vec<Body>,
        theta: f64,
        dt: f64,
        max_depth: usize,
        trail_every: usize,
    ) -> Self {
        Self::with_force_model(
            bodies,
            Box::new(GravityForceModel::new(theta, max_depth)),
            dt,
            trail_every,
        )
    }

    /// Create a simulation with an arbitrary pluggable force model.
    ///
    /// Use this constructor to inject alternative gravity engines (direct
    /// O(N²), GPU kernel, post-Newtonian, …) at construction time.
    /// The force model can also be swapped at runtime via [`set_force_model`].
    pub fn with_force_model(
        bodies: Vec<Body>,
        force_model: Box<dyn ForceModel>,
        dt: f64,
        _trail_every: usize,
    ) -> Self {
        let total_mass = bodies.iter().map(|b| b.mass).sum();
        let names = {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for b in &bodies {
                acc.push(helpers::auto_name(b.material, &acc));
            }
            acc
        };

        let theta = force_model.theta();
        let (r_min, softening_max) = helpers::compute_closeness(&bodies);

        Self {
            bodies,
            total_mass,
            last_kinetic: 0.0,
            last_potential: 0.0,
            initial_energy: None,
            rel_energy_error: 0.0,
            force_model,
            scratch_acc: Vec::new(),
            integrator: make_integrator(IntegratorKind::VelocityVerlet),
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
            seed: 0,
            hooks: HookRegistry::new(),
            stop_requested: false,
            pending_com_shift: (0.0, 0.0),
        }
    }
}
