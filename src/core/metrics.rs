/// A snapshot of every observable simulation quantity at a single instant.
///
/// `Metrics` is a pure data-transfer object: it carries no logic, performs
/// no computation, and holds no mutable state.  It is produced each frame by
/// [`crate::core::system::System::metrics`] and consumed by the UI and
/// diagnostics layers.
///
/// All values reflect the state **after** the most recent integration step,
/// including any collision merges and COM re-centring that occurred in that
/// step.
#[derive(Debug, Clone, Copy)]
pub struct Metrics {
    // ── Energetics ────────────────────────────────────────────────────────── //
    pub kinetic: f64,
    pub potential: f64,
    pub total_energy: f64,

    /// Relative energy drift since the first step (or since the last baseline
    /// reset caused by a COM-velocity correction or a collision merge):
    ///   δE = (E_now − E_baseline) / |E_baseline_scale|
    pub rel_energy_error: f64,

    /// Largest |δE| ever recorded in this run.
    pub max_rel_energy_error: f64,

    // ── Angular momentum & COM ─────────────────────────────────────────────── //
    /// Total orbital angular momentum about the (re-centred) origin.
    pub angular_momentum_z: f64,

    /// Relative angular momentum drift since the first step:
    ///   δLz = (Lz_now − Lz_baseline) / |Lz_scale|
    pub rel_angular_momentum_error: f64,

    /// Largest |δLz| ever recorded in this run.
    pub max_rel_angular_momentum_error: f64,

    pub com_x: f64,
    pub com_y: f64,
    pub com_vx: f64,
    pub com_vy: f64,

    // ── Physics config ────────────────────────────────────────────────────── //
    /// Effective gravitational multiplier (G_eff = G₀ · g_factor).
    pub g_factor: f64,

    // ── Adaptive integrator state ─────────────────────────────────────────── //
    /// Current Barnes-Hut opening angle θ (dimensionless).
    pub theta: f64,
    /// Actual time-step used in the last integration step.
    pub dt: f64,

    /// Smoothed relative error attributed to θ (theta-controller input).
    pub theta_fixed_rel_error: f64,
    /// Relative truncation error attributed to dt (dt-controller input).
    pub dt_fixed_rel_error: f64,

    pub last_theta_error_norm: f64,
    pub theta_error_smoothed_norm: f64,
    pub dt_controller_state: f64,

    // ── Diagnostics ───────────────────────────────────────────────────────── //
    pub max_acc: f64,
    pub jerk: f64,

    // ── Collision bookkeeping ─────────────────────────────────────────────── //
    /// Number of pairwise merges that occurred during the last integration step.
    pub merges_this_step: usize,
    /// Number of elastic / partial-restitution bounces during the last step.
    pub bounces_this_step: usize,
    /// Number of pairs within 2× contact distance but NOT resolved.
    /// A non-zero value means the dt_controller should reduce dt next step to
    /// prevent future tunneling.
    pub near_miss_count: usize,
    /// Number of fragment bodies spawned during the last step (debris events).
    pub fragments_spawned_this_step: usize,
    /// Number of hit-and-run events during the last step.
    pub hit_and_runs_this_step: usize,
    /// Cumulative dust mass (ejecta below tracking threshold) since simulation start.
    pub total_dust_mass: f64,
}
