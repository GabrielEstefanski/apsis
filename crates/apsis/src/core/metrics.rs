use crate::core::adaptive::DtMode;
use crate::physics::integrator::{AdaptiveStats, IntegratorKind};

/// Snapshot of the physical state after the most recent integration
/// step. Pure data, no logic.
///
/// Absolute conservation drifts are always defined; the relative forms
/// are `None` when the baseline is below
/// [`crate::core::system::MIN_RELATIVE_DENOMINATOR`] (round-off
/// dominates the metric there).
#[derive(Debug, Clone, Copy)]
pub struct Metrics {
    // ── Energetics ────────────────────────────────────────────────────────── //
    pub kinetic: f64,
    pub potential: f64,
    pub total_energy: f64,

    /// Initial total energy (set after the first force evaluation).
    pub initial_energy: f64,

    /// Absolute energy drift `E_now − E_initial`.
    pub abs_energy_error: f64,

    /// Relative energy drift `(E_now − E_initial) / |E_initial|`,
    /// or `None` when `|E_initial|` is below
    /// [`crate::core::system::MIN_RELATIVE_DENOMINATOR`].
    pub rel_energy_error: Option<f64>,

    // ── Angular momentum & COM ─────────────────────────────────────────────── //
    pub angular_momentum_z: f64,

    /// Initial Lz (set after the first state evaluation).
    pub initial_angular_momentum_z: f64,

    /// Absolute angular momentum drift `|Lz_now − Lz_initial|`.
    pub abs_angular_momentum_error: f64,

    /// Relative angular momentum drift, or `None` when `|Lz_initial|`
    /// is below [`crate::core::system::MIN_RELATIVE_DENOMINATOR`].
    pub rel_angular_momentum_error: Option<f64>,

    pub com_x: f64,
    pub com_y: f64,
    pub com_z: f64,
    pub com_vx: f64,
    pub com_vy: f64,
    pub com_vz: f64,

    // ── Time ──────────────────────────────────────────────────────────────── //
    pub t: f64,
    pub steps: u64,

    // ── Simulation parameters ─────────────────────────────────────────────── //
    pub integrator_kind: IntegratorKind,

    /// Effective gravitational multiplier (G_eff = G₀ · g_factor).
    pub g_factor: f64,

    /// Barnes–Hut opening angle θ. Meaningful only when
    /// [`force_is_direct`](Self::force_is_direct) is `false`.
    pub theta: f64,

    /// `true` when the force model skips Barnes–Hut entirely
    /// ([`ForceModel::is_deterministic`] at snapshot time); θ has no
    /// effect then.
    pub force_is_direct: bool,

    /// Current integration timestep (may differ from `user_dt` when
    /// `dt_mode == DtMode::Adaptive`).
    pub dt: f64,

    /// User-requested timestep (the adaptation target; equals `dt` when
    /// `dt_mode == DtMode::Fixed`).
    pub user_dt: f64,

    /// Active timestep policy. `Adaptive` breaks symplecticity.
    pub dt_mode: DtMode,

    /// Whether the adaptive Barnes–Hut θ controller is active. Varying θ
    /// changes per-step force accuracy (unlike adaptive dt, it does not
    /// break symplecticity).
    pub adaptive_theta: bool,

    // ── Diagnostics ───────────────────────────────────────────────────────── //
    pub max_acc: f64,
    pub jerk: f64,
    pub max_vel: f64,

    /// `true` if the last step was accepted under duress (adaptive
    /// integrator at its minimum step without meeting tolerance);
    /// always `false` for fixed-step integrators.
    pub last_step_degraded: bool,

    // ── Geometry diagnostics ─────────────────────────────────────────────── //
    /// Minimum pairwise separation observed at the last step (simulation units).
    ///
    /// Set to `f64::MAX` when fewer than 2 bodies are present or when N is too
    /// large for O(N²) computation (> [`N_CLOSENESS_THRESHOLD`]).
    pub r_min: f64,

    /// Squared softening of the active gravity kernel. Zero for the
    /// default Newton kernel; positive when a Plummer (or other softened)
    /// kernel is selected.
    pub kernel_epsilon_squared: f64,

    // ── Timestep guidance ─────────────────────────────────────────────────── //
    /// Recommended timestep derived from the current system state.
    ///
    /// Computed as the minimum of two standard N-body criteria:
    ///
    /// 1. **Power et al. (2003) acceleration criterion:**
    ///    `dt_acc = η · √(ε_min / a_max)`, η = 0.05
    ///    — ensures each body moves less than one softening length per step.
    ///
    /// 2. **Aarseth jerk criterion:**
    ///    `dt_jerk = η · √(a_max / j_max)`
    ///    — limits the fractional change in acceleration per step.
    ///    Available only after the first step (jerk is zero before).
    ///
    /// `None` when no bodies are present or before the first force evaluation.
    ///
    /// # References
    /// - Power et al. (2003). MNRAS 338, 14–34. §3.
    /// - Aarseth, S. J. (2003). *Gravitational N-Body Simulations*. Cambridge. §2.
    pub recommended_dt: Option<f64>,

    /// Cumulative adaptive-integrator counters (sub-steps, rejections,
    /// Picard iterations, degraded accepts). `None` for fixed-step
    /// integrators; `Some(..)` only when the active integrator is adaptive
    /// (currently: IAS15). See [`AdaptiveStats`] for field semantics.
    pub adaptive_stats: Option<AdaptiveStats>,
}
