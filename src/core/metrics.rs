use crate::physics::integrator::Integrator;

/// Snapshot of the physical state of the simulation at a single instant.
///
/// `Metrics` is a pure data-transfer object with no logic or side effects.
/// It represents the system state **after the most recent integration step**.
///
/// # Design principles
///
/// - No historical accumulation (no "max" tracking)
/// - No heuristic scaling
/// - All values are directly interpretable physically
///
/// # Notes
///
/// - Energy and angular momentum errors are measured relative to
///   their initial values
/// - These errors should remain bounded in a correct symplectic simulation
#[derive(Debug, Clone, Copy)]
pub struct Metrics {
    // ── Energetics ────────────────────────────────────────────────────────── //
    pub kinetic: f64,
    pub potential: f64,
    pub total_energy: f64,

    /// Relative energy drift:
    ///   δE = (E_now − E_initial) / |E_initial|
    pub rel_energy_error: f64,

    // ── Angular momentum & COM ─────────────────────────────────────────────── //
    /// Total angular momentum (z-component).
    pub angular_momentum_z: f64,

    /// Relative angular momentum drift:
    ///   δLz = (Lz_now − Lz_initial) / |Lz_initial|
    pub rel_angular_momentum_error: f64,

    /// Absolute angular momentum drift:
    ///   |Lz_now − Lz_initial|
    ///
    /// This is always meaningful, even when Lz ≈ 0.
    pub abs_angular_momentum_error: f64,

    pub com_x: f64,
    pub com_y: f64,
    pub com_vx: f64,
    pub com_vy: f64,

    // ── Time ──────────────────────────────────────────────────────────────── //
    /// Total simulated time elapsed.
    pub t: f64,

    /// Number of integration steps completed.
    pub steps: u64,

    // ── Simulation parameters ─────────────────────────────────────────────── //
    /// Active integration algorithm.
    pub integrator: Integrator,

    /// Effective gravitational multiplier (G_eff = G₀ · g_factor).
    pub g_factor: f64,

    /// Barnes–Hut opening angle θ.
    pub theta: f64,

    /// Fixed integration time step.
    pub dt: f64,

    // ── Diagnostics ───────────────────────────────────────────────────────── //
    pub max_acc: f64,
    pub jerk: f64,
    pub max_vel: f64,

    // ── Softening diagnostics ─────────────────────────────────────────────── //
    /// Minimum pairwise separation observed at the last step (simulation units).
    ///
    /// Set to `f64::MAX` when fewer than 2 bodies are present or when N is too
    /// large for O(N²) computation (> [`N_CLOSENESS_THRESHOLD`]).
    pub r_min: f64,

    /// Maximum effective pairwise softening length at the last step:
    ///   `ε_ij = √((ε²_i + ε²_j) / 2)`,  maximised over all pairs.
    ///
    /// Zero when no pairs exist.
    pub softening_max: f64,
}
