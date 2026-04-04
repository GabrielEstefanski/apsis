/// Runtime physics configuration — tunable constants and display unit labels.
///
/// `g_factor` is the only value that changes simulation behaviour; the label
/// fields are cosmetic and only affect how values are annotated in the UI.
#[derive(Clone)]
pub struct PhysicsConfig {
    /// Effective gravitational strength multiplier.
    ///
    /// All gravitational accelerations are scaled by this factor, so the
    /// simulation behaves as if `G_eff = G₀ · g_factor` (where `G₀ = 1.0`
    /// is the natural simulation unit).
    ///
    /// - `1.0` — default; matches template initial conditions exactly
    /// - `< 1.0` — weaker gravity; orbits need less velocity to stay bound
    /// - `> 1.0` — stronger gravity; faster collapse, tighter orbits
    pub g_factor: f64,

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
            g_factor: 1.0,
            mass_label: "M".into(),
            dist_label: "u".into(),
            time_label: "t".into(),
        }
    }
}
