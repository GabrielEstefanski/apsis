use crate::physics::integrator::Integrator;

/// Runtime physics configuration — tunable constants and display unit labels.
#[derive(Clone)]
pub struct PhysicsConfig {
    /// Active integration algorithm.
    pub integrator: Integrator,

    // ── Force accuracy ────────────────────────────────────────────────────────

    /// Barnes–Hut opening angle θ ∈ [0.05, 1.5].
    ///
    /// Controls the accuracy/performance trade-off of the tree force solver:
    /// - `θ → 0` : every pair evaluated exactly → O(N²), maximum accuracy
    /// - `θ = 0.5`: standard compromise (default)
    /// - `θ → 1.5`: aggressively approximated → O(N log N), less accurate
    pub theta: f64,

    /// Global Plummer softening scale applied on top of the per-body default.
    ///
    /// Per-body default: `ε = EPS_BASE · m^(1/3)` (Plummer mass-proportional).
    /// This scale multiplies that value: `ε_eff = ε_default · softening_scale`.
    ///
    /// - `1.0` — default; physically motivated size
    /// - `< 1.0` — sharper forces, closer to point-mass singularity
    /// - `> 1.0` — smoother forces, suppress close encounters
    pub softening_scale: f64,

    // ── Gravity ───────────────────────────────────────────────────────────────

    /// Effective gravitational strength multiplier.
    ///
    /// Scales all accelerations: `G_eff = G₀ · g_factor` (G₀ = 1.0).
    pub g_factor: f64,

    // ── Trails ────────────────────────────────────────────────────────────────

    /// Trail sampling interval: record a trail point every N rendered frames.
    pub trail_every: usize,

    // ── Unit labels (cosmetic) ────────────────────────────────────────────────

    /// Display label appended to mass values, e.g. `"M"`, `"M☉"`, `"kg"`.
    pub mass_label: String,

    /// Display label appended to distance values, e.g. `"u"`, `"AU"`, `"pc"`.
    pub dist_label: String,

    /// Display label appended to time values, e.g. `"t"`, `"yr"`, `"Myr"`.
    pub time_label: String,
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            integrator: Integrator::VelocityVerlet,
            theta: 0.5,
            softening_scale: 1.0,
            g_factor: 1.0,
            trail_every: 1,
            mass_label: "M".into(),
            dist_label: "u".into(),
            time_label: "t".into(),
        }
    }
}
