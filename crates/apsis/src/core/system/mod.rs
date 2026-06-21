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
//! | Extra forces | [`HamiltonianOperator`] / [`NonConservativeOperator`] | none |
//!
//! ## Module layout
//!
//! | Module | Responsibility |
//! |---|---|
//! | `mod.rs` | [`System`] struct, constructors |
//! | `step` | `step()` and conservation-law tracking |
//! | `bodies` | body CRUD, names, COM calibration |
//! | `config` | getters/setters (θ, dt, integrator, kernel, …) |
//! | `metrics` | [`Metrics`] assembly and recommended-dt |
//! | `orbital` | osculating-element cache |
//! | `perturbations` | non-gravitational force registration |
//! | `helpers` | free functions (naming, closeness) |

pub(crate) mod bodies;
pub(crate) mod config;
pub(crate) mod helpers;
pub(crate) mod metrics;
pub(crate) mod orbital;
pub(crate) mod perturbations;
pub(crate) mod regime;
pub(crate) mod step;
#[cfg(test)]
mod tests;

pub use regime::MIN_RELATIVE_DENOMINATOR;

use crate::core::adaptive::{DtAdaptationConfig, DtController, DtMode, ThetaController};
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::hooks::HookRegistry;
use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::{
    ForceModel, GravityForceModel, HamiltonianOperator, Integrator, IntegratorKind,
    NonConservativeOperator, Operator, make_integrator,
};
use crate::physics::orbital::OrbitalElements;
use crate::templates::instantiate::instantiate;
use crate::templates::kind::TemplateKind;
use crate::units::UnitSystem;

// ── Default parameters (used by System::new) ──────────────────────────────────

/// Default Barnes-Hut opening angle. Standard in the literature for accuracy
/// vs speed on mixed scenes.
const DEFAULT_THETA: f64 = 0.6;

/// Default fixed timestep. Safe for unit-scale scenarios and Yoshida-4; users
/// with stiff scenes should call [`System::with_dt`] explicitly or enable
/// adaptive dt via [`System::set_dt_mode`].
const DEFAULT_DT: f64 = 1e-4;

/// Default maximum octree depth. Covers scenes up to ~10⁹ spatial extent
/// before degrading to leaf splits; rarely touched.
const DEFAULT_MAX_DEPTH: usize = 32;

/// Central simulation state for an N-body gravitational system.
pub struct System {
    /// Bodies participating in the simulation.
    pub(crate) bodies: Vec<Body>,

    /// Total mass of the system (used for COM recentering).
    pub(crate) total_mass: f64,

    /// Last computed energies.
    pub(crate) last_kinetic: f64,
    pub(crate) last_potential: f64,

    /// Initial total energy (used as reference for drift).
    pub(crate) initial_energy: Option<f64>,

    /// Absolute energy drift `E - E_initial` (signed).
    pub(crate) abs_energy_error: f64,

    /// Relative energy drift; `None` in precision-limited regime.
    pub(crate) rel_energy_error: Option<f64>,

    /// Pluggable force model — default: Barnes-Hut gravity.
    ///
    /// Swap via [`System::set_force_model`] to use any [`ForceModel`]
    /// implementation (direct O(N²), GPU kernels, PN corrections, …).
    pub(crate) force_model: Box<dyn ForceModel>,

    /// Scratch buffer for accelerations — reused every step.
    pub(crate) scratch_acc: Vec<Vec3>,

    /// Active integration algorithm (trait object).
    pub(crate) integrator: Box<dyn Integrator>,

    /// Cached osculating orbital elements — one slot per body.
    pub(crate) orbital_cache: Vec<Option<OrbitalElements>>,

    /// Diagnostics subsystem.
    pub(crate) diagnostics: DiagnosticsComputer,
    pub(crate) last_diag: SimulationDiagnostics,

    /// `true` if the most recent step was accepted under duress (an
    /// IAS15 sub-step that hit the `DT_MIN` floor without satisfying
    /// the tolerance). Mirrors [`StepResult::degraded`]; surfaced via
    /// [`Metrics::last_step_degraded`].
    pub(crate) last_step_degraded: bool,

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

    /// Active unit system. Frozen post-construction; see [`crate::units`].
    pub(crate) units: UnitSystem,

    /// Effective `G` in canonical units. Seeded from `units.g()` at
    /// construction; callers may rescale it independently via
    /// [`set_g_factor`](Self::set_g_factor).
    pub(crate) g_factor: f64,

    /// Initial Lz — conservation baseline.
    pub(crate) initial_angular_momentum: Option<f64>,

    /// Relative Lz drift; `None` in precision-limited regime.
    pub(crate) rel_angular_momentum_error: Option<f64>,

    /// Absolute Lz drift `|Lz - Lz_initial|`.
    pub(crate) abs_angular_momentum_error: f64,

    /// Minimum pairwise separation from the most recent step.
    pub(crate) r_min: f64,

    /// Optional close-encounter advisory threshold. When `Some(t)`, the
    /// step loop classifies `r_min` against `t` via
    /// [`crate::physics::encounter::EncounterFlag`] and emits a
    /// `warn_diag!` event on the `Far`/`Approaching` → `Close`
    /// transition. `None` (the default) disables the diagnostic.
    pub(crate) close_encounter_threshold: Option<f64>,

    /// Encounter flag from the most recent step. Tracked across steps so
    /// the warn-on-transition rule fires exactly once per descent into
    /// the `Close` band; stays observable to external readers between
    /// steps.
    pub(crate) last_encounter_flag: crate::physics::encounter::EncounterFlag,

    /// Hamiltonian-class perturbations (force = −∇V derivable, with
    /// energy contribution summed into [`total_energy`](Self::total_energy)).
    /// Symplectic integrators preserve invariants when only operators
    /// of this class are registered.
    pub(crate) hamiltonian_perturbations: Vec<Box<dyn HamiltonianOperator>>,

    /// Non-conservative perturbations (force without a Hamiltonian:
    /// drag, radiation reaction, dissipative coupling). Symplectic
    /// integrators degrade silently when one of these is registered;
    /// the registration site emits a `warn_diag` so the broken
    /// invariant is documented.
    pub(crate) non_conservative_perturbations: Vec<Box<dyn NonConservativeOperator>>,

    /// Pure observers — read-only operators called at synchronized
    /// step boundaries (Shadow Hamiltonian tracker, audit trail
    /// emitters). Contribute no force, no energy.
    pub(crate) observers: Vec<Box<dyn Operator>>,

    /// Set of `(operator, bound)` pairs already reported via
    /// `warn_diag` for regime-of-validity violations. The dedup state
    /// is per-`System`-instance, lifetime-of-process: a violation that
    /// is true at registration AND persists through dynamic checks
    /// emits exactly one warning. Cleared by
    /// [`reset_regime_warnings`](Self::reset_regime_warnings) when the
    /// caller wants to re-arm the bus (e.g. after a deliberate
    /// scenario change).
    pub(crate) regime_warnings_emitted: std::collections::HashSet<(&'static str, &'static str)>,

    /// Reproducibility seed. Consumed by preset builders and cluster spawners.
    /// Persisted in snapshots so a run can be replayed exactly.
    pub(crate) seed: u64,

    /// When the system was built via [`System::from_template`], remembers the
    /// preset so [`with_seed`](System::with_seed) can rebuild the body list
    /// with a new seed without forcing a separate `from_template_with_seed`
    /// entry point. `None` after manual construction, snapshot restore, or
    /// any mutation that invalidates the "bodies equal `kind.build(seed)`"
    /// invariant.
    pub(crate) template_source: Option<TemplateKind>,

    /// Registered observer/command hooks. Dispatched from [`System::step`].
    pub(crate) hooks: HookRegistry,

    /// Set by a [`Command::Stop`](crate::core::hooks::Command::Stop) request.
    /// Honoured by the run loop; downstream callers may inspect it.
    pub(crate) stop_requested: bool,

    /// Land [`integrate_until`](System::integrate_until) exactly on
    /// `t_end` by clipping the final step (REBOUND's
    /// `exact_finish_time=1`). On by default; disable via
    /// [`set_exact_finish_time`](System::set_exact_finish_time) to run
    /// whole steps past `t_end`, e.g. to preserve a fixed-step
    /// symplectic rhythm across sampling boundaries.
    pub(crate) exact_finish_time: bool,

    /// Dense-output snapshot from the most recent integration step.
    /// Produced each step; consumed by downstream interpolators (e.g.
    /// trail samplers, sub-step position renderers) that need a smooth
    /// curve between integrator step boundaries.
    pub(crate) last_dense_snapshot: Option<crate::physics::integrator::DenseSnapshot>,

    /// Set on first [`System::finish`] call so subsequent ones are no-ops.
    pub(crate) finished: bool,
}

impl System {
    /// Create a simulation from a body list and an explicit unit system.
    ///
    /// The `units` argument is **mandatory and immutable**. Every body
    /// coordinate, velocity, mass, and `dt` passed in (now or later)
    /// is interpreted in the canonical units of this [`UnitSystem`];
    /// passing a value in the wrong unit is a silent physical error,
    /// not a runtime error. The unit system cannot be changed after
    /// construction — the only way to "change units" is to rebuild
    /// the `System`.
    ///
    /// Defaults for everything else (integrator, `dt`, θ, kernel,
    /// max octree depth) match the conventions of small-N research
    /// scripts:
    ///
    /// | Parameter              | Default                     |
    /// |------------------------|-----------------------------|
    /// | Integrator             | Yoshida 4th order (symplectic) |
    /// | dt                     | `1e-4` simulation time units |
    /// | Barnes-Hut θ           | `0.6`                        |
    /// | Max octree depth     | `32`                         |
    /// | Softening scale        | `1.0`                        |
    ///
    /// Override any of these with the fluent [`with_*`](Self::with_dt)
    /// builder methods.
    ///
    /// ```ignore
    /// use apsis::core::system::System;
    /// use apsis::domain::body::Body;
    /// use apsis::physics::integrator::IntegratorKind;
    /// use apsis::units::UnitSystem;
    ///
    /// let sun = Body::star(1.0);
    /// let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
    ///
    /// let mut sys = System::new(vec![sun, earth], UnitSystem::canonical())
    ///     .with_integrator(IntegratorKind::Ias15)
    ///     .with_dt(1e-4);
    ///
    /// sys.integrate_for(100.0);
    /// println!("dE/E = {:.3e}", sys.energy_delta());
    /// ```
    pub fn new(bodies: Vec<Body>, units: UnitSystem) -> Self {
        Self::with_force_model_inner(
            bodies,
            Box::new(GravityForceModel::new(DEFAULT_THETA, DEFAULT_MAX_DEPTH)),
            DEFAULT_DT,
            units,
        )
    }

    /// Construct with a pluggable force model (direct O(N²), GPU kernel,
    /// post-Newtonian, …). Escape hatch for advanced users; most callers
    /// prefer [`new`](Self::new) followed by builder methods.
    pub fn with_force_model(
        bodies: Vec<Body>,
        force_model: Box<dyn ForceModel>,
        units: UnitSystem,
    ) -> Self {
        Self::with_force_model_inner(bodies, force_model, DEFAULT_DT, units)
    }

    /// Construct a system from a built-in preset.
    ///
    /// Defaults match [`System::new`]; override any with `.with_*` builder
    /// methods. For randomised presets (e.g. [`TemplateKind::JupiterTrojans`])
    /// the initial seed is `0`; change it via `.with_seed(42)` — which
    /// rebuilds the body list with the new seed automatically, keeping a
    /// single builder entry point for the whole construction chain.
    ///
    /// ```ignore
    /// use apsis::core::system::System;
    /// use apsis::physics::integrator::IntegratorKind;
    /// use apsis::templates::TemplateKind;
    ///
    /// let mut sys = System::from_template(TemplateKind::JupiterTrojans)
    ///     .with_seed(42)        // rebuilds the trojan layout with seed=42
    ///     .with_integrator(IntegratorKind::Ias15)
    ///     .with_dt(1e-4);
    /// ```
    ///
    /// For string-keyed lookup (config files, CLI input), resolve the name
    /// first via [`TemplateKind::from_name`]:
    ///
    /// ```ignore
    /// let kind = TemplateKind::from_name(&config.preset)?;
    /// let sys  = System::from_template(kind);
    /// ```
    pub fn from_template(kind: TemplateKind) -> Self {
        let template = kind.build(0);
        let named = instantiate(&template);
        // Presets are calibrated for G = 1 (Hénon).
        let mut sys = Self::new(Vec::new(), UnitSystem::canonical());
        sys.add_bodies(named);
        // Restore the template handle that add_bodies cleared, so a
        // follow-up .with_seed(...) can rebuild from the same preset.
        sys.template_source = Some(kind);
        sys
    }

    fn with_force_model_inner(
        bodies: Vec<Body>,
        force_model: Box<dyn ForceModel>,
        dt: f64,
        units: UnitSystem,
    ) -> Self {
        let total_mass = bodies.iter().map(|b| b.mass).sum();
        // Resolve each body's name against the running set so explicit
        // names (e.g. "Sun", "Mercury") ride through and `None` slots
        // receive monotonic `"Body N"` placeholders.
        let mut acc_names: Vec<String> = Vec::with_capacity(bodies.len());
        let bodies: Vec<Body> = bodies
            .into_iter()
            .map(|mut b| {
                let resolved =
                    helpers::resolved_name(b.name.take(), helpers::DEFAULT_NAME_PREFIX, &acc_names);
                acc_names.push(resolved.clone());
                b.name = Some(resolved);
                b
            })
            .collect();

        let theta = force_model.theta();
        let r_min = helpers::compute_closeness(&bodies);

        Self {
            bodies,
            total_mass,
            last_kinetic: 0.0,
            last_potential: 0.0,
            initial_energy: None,
            abs_energy_error: 0.0,
            rel_energy_error: None,
            force_model,
            scratch_acc: Vec::new(),
            integrator: make_integrator(IntegratorKind::Ias15),
            orbital_cache: Vec::new(),
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            last_step_degraded: false,
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
            }),
            theta_ctrl: ThetaController::new(1e-3, 0.05, 1.5).with_initial_theta(theta),
            adaptive_theta: false,
            g_factor: units.g(),
            units,
            initial_angular_momentum: None,
            rel_angular_momentum_error: None,
            abs_angular_momentum_error: 0.0,
            r_min,
            close_encounter_threshold: None,
            last_encounter_flag: crate::physics::encounter::EncounterFlag::Far,
            hamiltonian_perturbations: Vec::new(),
            non_conservative_perturbations: Vec::new(),
            observers: Vec::new(),
            regime_warnings_emitted: std::collections::HashSet::new(),
            seed: 0,
            hooks: HookRegistry::new(),
            stop_requested: false,
            exact_finish_time: true,
            last_dense_snapshot: None,
            template_source: None,
            finished: false,
        }
    }

    // ── Fluent construction builder ───────────────────────────────────────────
    // Consume-and-return-Self chain. Runtime mutation lives in the `set_*`
    // counterparts in [`crate::core::system::config`].

    /// Fixed timestep for the integrator.
    #[inline]
    #[must_use]
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.set_dt(dt);
        self
    }

    /// Barnes-Hut opening angle `θ` (accuracy ↔ speed trade-off).
    #[inline]
    #[must_use]
    pub fn with_theta(mut self, theta: f64) -> Self {
        self.set_theta(theta);
        self
    }

    /// Integrator choice (see [`IntegratorKind`]).
    #[inline]
    #[must_use]
    pub fn with_integrator(mut self, kind: IntegratorKind) -> Self {
        self.set_integrator(kind);
        self
    }

    /// Maximum Barnes-Hut octree depth.
    ///
    /// Most scenes do not need to touch this; the default (32) covers a
    /// spatial extent of ~10⁹ before degrading to forced leaf splits.
    #[inline]
    #[must_use]
    pub fn with_max_depth(mut self, max_depth: usize) -> Self {
        let theta = self.force_model.theta();
        self.set_force_model(Box::new(GravityForceModel::new(theta, max_depth)));
        self
    }

    /// RNG seed.
    ///
    /// For systems built from scratch with [`System::new`], this only sets
    /// the runtime seed forwarded to seeded integrators and samplers.
    ///
    /// For systems built via [`System::from_template`], this **also rebuilds
    /// the body list** using the preset's builder with the new seed —
    /// keeping a single fluent chain for the whole construction, no
    /// separate `from_template_with_seed` entry point.
    ///
    /// The rebuild is a no-op for deterministic presets (most of them) and
    /// regenerates only the randomised ones (Jupiter Trojans, cluster
    /// layouts, …). `template_source` is preserved, so subsequent
    /// `.with_seed(...)` calls still rebuild.
    #[must_use]
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.set_seed(seed);
        if let Some(kind) = self.template_source {
            self.rebuild_from_template(kind, seed);
        }
        self
    }

    fn rebuild_from_template(&mut self, kind: TemplateKind, seed: u64) {
        self.bodies.clear();
        self.total_mass = 0.0;
        self.initial_energy = None;
        self.initial_angular_momentum = None;
        self.rel_energy_error = None;
        self.abs_energy_error = 0.0;
        self.rel_angular_momentum_error = None;
        self.abs_angular_momentum_error = 0.0;
        let template = kind.build(seed);
        let named = instantiate(&template);
        self.add_bodies(named);
        // add_bodies cleared template_source; this is an internal
        // rebuild path, so restore the invariant.
        self.template_source = Some(kind);
    }
}
