use crate::domain::body_preset::{self, BodyClass, BodyPreset};
use std::f64::consts::PI;

/// Base softening length for a body of mass 1.0.
/// Per-body softening scales as `EPS_BASE * mass^(1/3)`, so each body's
/// softening volume is proportional to its mass — physically motivated by
/// the Plummer-equivalent equal-mass softening criterion.
///
/// Exposed as `#[doc(hidden)] pub` so calibration and softening-scaling
/// tools can read the baseline without committing it as a stable API
/// constant — future calibration work may change the numeric value.
#[doc(hidden)]
pub const EPS_BASE: f64 = 0.02;

/// Point-mass body: kinematics, mass, and the small set of physical
/// properties read by gravitational and radiation force evaluators.
///
/// Body owns its physical state directly. It carries no taxonomy
/// reference — material classifications are construction-time presets
/// (see [`crate::domain::body_preset`]), not runtime fields. The
/// resulting struct matches the REBOUND `reb_particle` shape: every
/// field is a quantity a force evaluator might read, and every default
/// is sensible for a non-emitting non-receiving point mass.
#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub pos_x: f64,
    pub pos_y: f64,
    pub pos_z: f64,
    pub vel_x: f64,
    pub vel_y: f64,
    pub vel_z: f64,
    pub mass: f64,

    /// Gravitational softening length ε for this body.
    ///
    /// Pairwise: ε²_ij = (ε²_i + ε²_j) / 2.
    ///
    /// Calibrated by `System::calibrate_softening`.
    pub softening: f64,

    /// True physical radius derived from mass and density.
    ///
    /// Used by:
    /// - energy calculations (e.g. disruption threshold Q*)
    /// - radiation pressure (cross section ∝ r²)
    /// - rendering and collision geometry
    pub physical_radius: f64,

    /// Bulk density of the body: ρ = m / V, V = 4/3 π r³.
    ///
    /// Primary size property — physical radius is recomputed from
    /// this via [`Body::sync_physical_properties`] whenever mass or
    /// density change.
    pub density: f64,

    /// Display colour `[R, G, B]`. Set by the construction preset
    /// (e.g. brown for [`body_preset::ROCKY`], blue-grey for
    /// [`body_preset::ICY`]) and overridable per body via
    /// [`Body::with_color`]. Never read by the physics path.
    pub color: [u8; 3],

    /// Bolometric luminosity in internal energy · time⁻¹ units.
    ///
    /// Set by the construction preset for luminous classes
    /// ([`body_preset::STAR`], [`body_preset::BROWN_DWARF`],
    /// [`body_preset::WHITE_DWARF`]); zero for everything else.
    /// Static after construction — mutating mass or density does not
    /// recompute luminosity automatically. Override manually with
    /// [`Body::with_luminosity`] if a sim layer needs to refresh it.
    pub luminosity: f64,

    /// Radiation-pressure receiver coefficient `Q_pr` (Burns,
    /// Lamy & Soter 1979). Zero on stars and large planets;
    /// positive on dust grains, asteroid surfaces, comets.
    ///
    /// Read by [`crate::physics::radiation`] when computing per-body
    /// radiation forces. Bodies with `q_pr == 0.0` are silently
    /// skipped by the radiation pipeline.
    pub q_pr: f64,

    /// UX taxonomy bucket — Star, Planet, Moon, Asteroid, Comet, or
    /// Unknown. Set from the construction preset's
    /// [`BodyPreset::default_class`] (overridable per body via
    /// [`Body::with_class`]); never read by physics. The render and
    /// inspector layers use it for category-level filters and
    /// grouping.
    pub class: BodyClass,

    /// Bond albedo — fraction of incident bolometric radiation that
    /// the body reflects, integrated over all wavelengths and phase
    /// angles. Dimensionless, range `[0, 1]`. Read by the
    /// [`crate::physics::photometry`] pipeline to compute the
    /// reflected flux at the observer; ignored by force evaluators.
    ///
    /// Defaults are taken from [`BodyPreset::default_albedo`] (which
    /// are class-typical placeholders, e.g. `0.10` for asteroid,
    /// `0.30` for rocky). Templates that quote real bodies override
    /// to the published Bond value (Earth `0.306`, Moon `0.11`,
    /// Vesta `0.42`).
    ///
    /// Stars carry `0.0` — their flux is emitted, not reflected.
    pub albedo: f64,
}

/// Body payload with an optional explicit display name.
///
/// Returned by [`Body::named`] and by the template catalog. Consumed by
/// [`System::add_named_body`](crate::core::system::System::add_named_body)
/// to preserve authored names; otherwise the system derives a stable
/// preset-based fallback.
#[derive(Clone, Debug)]
pub struct NamedBody {
    pub body: Body,
    pub name: Option<String>,
}

impl Body {
    // ── Low-level constructor ────────────────────────────────────────────────
    //
    // Explicit mass and density; everything else takes a sensible default
    // (zero kinematics, mid-grey colour, no luminosity, no radiation
    // receiver). This is the form documented in the paper.

    /// Construct a body with explicit mass and density. Physical
    /// radius is derived as `r = (3m / 4πρ)^(1/3)`; softening uses
    /// the project default `EPS_BASE * m^(1/3)`. Position and
    /// velocity start at the origin with no velocity; tune via the
    /// fluent builder methods or the preset-based factories.
    pub fn new(mass: f64, density: f64) -> Self {
        let physical_radius = radius_from_density_mass(density, mass);
        let softening = default_softening(mass);
        Self {
            pos_x: 0.0,
            pos_y: 0.0,
            pos_z: 0.0,
            vel_x: 0.0,
            vel_y: 0.0,
            vel_z: 0.0,
            mass,
            softening,
            physical_radius,
            density,
            color: [180, 180, 180],
            luminosity: 0.0,
            q_pr: 0.0,
            class: BodyClass::Unknown,
            albedo: 0.30,
        }
    }

    /// Construct a body from a [`BodyPreset`] at the requested mass.
    ///
    /// Density, colour, `q_pr`, and luminosity are all taken from the
    /// preset; physical radius and softening are derived. Mass is
    /// assumed to be in solar units — the preset's luminosity model
    /// runs with that assumption. For non-solar unit systems, override
    /// the resulting body's `luminosity` field manually after
    /// construction.
    pub fn from_preset(preset: &BodyPreset, mass: f64) -> Self {
        // Preset density models live in human-readable kg/m³ so the
        // source reads like the literature; cross to coherent
        // M_☉/AU³ here, once, before any geometry runs.
        let density_kg_m3 = preset.density.density_at(mass);
        let density = density_kg_m3 * crate::templates::builders::KG_M3_TO_SOLAR_AU3;
        let physical_radius = radius_from_density_mass(density, mass);
        let softening = default_softening(mass);
        let luminosity = preset
            .luminosity
            .map(|src| src.compute(mass, physical_radius / SOLAR_RADIUS_AU))
            .unwrap_or(0.0);

        Self {
            pos_x: 0.0,
            pos_y: 0.0,
            pos_z: 0.0,
            vel_x: 0.0,
            vel_y: 0.0,
            vel_z: 0.0,
            mass,
            softening,
            physical_radius,
            density,
            color: preset.default_color,
            luminosity,
            q_pr: preset.default_q_pr,
            class: preset.default_class,
            albedo: preset.default_albedo,
        }
    }

    // ── Preset-named constructors (high-level ergonomic factories) ──────────
    //
    // Each is a one-liner over `from_preset` with the corresponding built-in
    // preset. Preserves the `Body::rocky(mass)` etc. spelling without
    // re-exposing the preset enum on the body itself.

    /// Star — main-sequence luminous body. See [`body_preset::STAR`].
    pub fn star(mass: f64) -> Self {
        Self::from_preset(&body_preset::STAR, mass)
    }

    /// Brown dwarf — sub-stellar deuterium burner.
    pub fn brown_dwarf(mass: f64) -> Self {
        Self::from_preset(&body_preset::BROWN_DWARF, mass)
    }

    /// White dwarf — degenerate stellar remnant.
    pub fn white_dwarf(mass: f64) -> Self {
        Self::from_preset(&body_preset::WHITE_DWARF, mass)
    }

    /// Gas giant — H/He envelope (Jupiter, Saturn, hot Jupiters).
    pub fn gas_giant(mass: f64) -> Self {
        Self::from_preset(&body_preset::GAS, mass)
    }

    /// Ice giant — water/methane envelope (Uranus, Neptune).
    pub fn ice_giant(mass: f64) -> Self {
        Self::from_preset(&body_preset::ICE_GIANT, mass)
    }

    /// Rocky body — silicate terrestrial planet or large rocky moon.
    pub fn rocky(mass: f64) -> Self {
        Self::from_preset(&body_preset::ROCKY, mass)
    }

    /// Icy body — water-ice dominated (icy moons, KBOs).
    pub fn icy(mass: f64) -> Self {
        Self::from_preset(&body_preset::ICY, mass)
    }

    /// Asteroid — rocky minor body.
    pub fn asteroid(mass: f64) -> Self {
        Self::from_preset(&body_preset::ASTEROID, mass)
    }

    /// Comet — volatile-rich minor body.
    pub fn comet(mass: f64) -> Self {
        Self::from_preset(&body_preset::COMET, mass)
    }

    // ── Fluent builder ────────────────────────────────────────────────────────
    //
    // Each method consumes and returns `Self`, so they chain naturally:
    //
    //     Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0)

    /// Position in simulation coordinates. `z` is left at its current value
    /// (default `0.0`); use [`at_3d`](Self::at_3d) for explicit 3D placement.
    #[inline]
    #[must_use]
    pub fn at(mut self, x: f64, y: f64) -> Self {
        self.pos_x = x;
        self.pos_y = y;
        self
    }

    /// Position in 3D simulation coordinates.
    #[inline]
    #[must_use]
    pub fn at_3d(mut self, x: f64, y: f64, z: f64) -> Self {
        self.pos_x = x;
        self.pos_y = y;
        self.pos_z = z;
        self
    }

    /// Velocity in simulation coordinates. `vz` is left at its current value
    /// (default `0.0`); use [`with_velocity_3d`](Self::with_velocity_3d) for
    /// explicit 3D motion.
    #[inline]
    #[must_use]
    pub fn with_velocity(mut self, vx: f64, vy: f64) -> Self {
        self.vel_x = vx;
        self.vel_y = vy;
        self
    }

    /// Velocity in 3D simulation coordinates.
    #[inline]
    #[must_use]
    pub fn with_velocity_3d(mut self, vx: f64, vy: f64, vz: f64) -> Self {
        self.vel_x = vx;
        self.vel_y = vy;
        self.vel_z = vz;
        self
    }

    /// Override the preset-default density. Physical radius is
    /// recomputed to match.
    #[inline]
    #[must_use]
    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
        self
    }

    /// Override the preset-default colour. Visualisation only — never
    /// affects the physics path.
    #[inline]
    #[must_use]
    pub fn with_color(mut self, color: [u8; 3]) -> Self {
        self.color = color;
        self
    }

    /// Set the bolometric luminosity directly. Use after construction
    /// when running outside the canonical solar unit system, or to
    /// install a non-default luminosity (a tabulated pulsar curve,
    /// for example).
    #[inline]
    #[must_use]
    pub fn with_luminosity(mut self, l: f64) -> Self {
        self.luminosity = l;
        self
    }

    /// Override the preset-default radiation-pressure receiver
    /// coefficient `Q_pr`. Set to `0.0` to opt the body out of the
    /// radiation pipeline.
    #[inline]
    #[must_use]
    pub fn with_q_pr(mut self, q_pr: f64) -> Self {
        self.q_pr = q_pr;
        self
    }

    /// Override the preset-default UX class. Use to tag a body whose
    /// preset does not match its role: Earth's Moon is constructed
    /// from [`body_preset::ROCKY`] (default class
    /// [`BodyClass::Planet`]) but should render under
    /// [`BodyClass::Moon`] for filter purposes.
    #[inline]
    #[must_use]
    pub fn with_class(mut self, class: BodyClass) -> Self {
        self.class = class;
        self
    }

    /// Override the preset-default Bond albedo. Templates of named
    /// bodies (Earth, Sun, Moon, Mercury, …) call this with their
    /// published value; bulk-anonymous bodies inherit the preset's
    /// class-typical placeholder.
    #[inline]
    #[must_use]
    pub fn with_albedo(mut self, albedo: f64) -> Self {
        self.albedo = albedo.clamp(0.0, 1.0);
        self
    }

    /// Attach an explicit display name, producing a [`NamedBody`] consumable
    /// by [`System::add_named_body`](crate::core::system::System::add_named_body).
    #[inline]
    #[must_use]
    pub fn named(self, name: impl Into<String>) -> NamedBody {
        NamedBody { body: self, name: Some(name.into()) }
    }

    /// Zero this body's Plummer softening length, producing the exact `1/r`
    /// potential for every interaction involving it.
    ///
    /// The simulator defaults every body to a preset-scaled Plummer
    /// softening (`EPS_BASE · mass^(1/3)`). For a solar-mass body this
    /// gives ε ≈ 0.02 AU — about 5 % of Mercury's perihelion distance —
    /// which introduces a numerical apsidal precession that can easily
    /// dominate a fine-physics signal (post-Newtonian precession, J2
    /// oblateness, tidal dissipation). Call this when the body participates
    /// in a measurement of such a deviation-from-Kepler effect.
    ///
    /// Equivalent to `body.softening = 0.0;` but expresses intent:
    ///
    /// ```ignore
    /// let sun     = Body::star(1.0).unsoftened();
    /// let mercury = Body::rocky(3e-6)
    ///     .at(0.307, 0.0)
    ///     .with_velocity(0.0, 1.98)
    ///     .unsoftened();
    /// ```
    ///
    /// See also [`System::with_exact_gravity`](crate::core::system::System::with_exact_gravity)
    /// to unsoften the whole system in one call.
    #[inline]
    #[must_use]
    pub fn unsoftened(mut self) -> Self {
        self.softening = 0.0;
        self
    }

    // ── Mutators ──────────────────────────────────────────────────────────────

    /// Recompute the physical radius from the current mass and density.
    ///
    /// Must be called whenever `mass` or `density` is mutated in place
    /// (e.g. via direct field assignment on a `&mut Body`). It intentionally
    /// does **not** touch the calibrated contact radius, which belongs to the
    /// numerical collision model rather than the body's physical geometry.
    pub fn sync_physical_properties(&mut self) {
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
    }

    /// Returns `true` if this body emits radiation.
    ///
    /// Equivalent to `self.luminosity > 0.0` but more expressive at call sites.
    #[inline]
    pub fn is_luminous(&self) -> bool {
        self.luminosity > 0.0
    }
}

/// Default softening before system-scale calibration.
pub(crate) fn default_softening(mass: f64) -> f64 {
    EPS_BASE * mass.abs().cbrt()
}

/// Radius of a uniform sphere from its bulk density and mass.
///
/// **Pure geometric.** `density` and `mass` must be in coherent
/// units of the active simulation system (e.g. for `solar_au`:
/// mass in M_☉, density in M_☉/AU³, returned radius in AU).
///
/// Callers receiving SI-format physical data (kg/m³ from NASA fact
/// sheets, GADGET parameter files, …) are responsible for
/// converting to coherent units at the boundary. The two boundary
/// crossings in this codebase are:
///
/// * `Body::from_preset`, which converts the preset's kg/m³ output
///   via [`crate::templates::builders::KG_M3_TO_SOLAR_AU3`]; and
/// * `templates::*`, which write per-body density as
///   `kg_m3_value * KG_M3_TO_SOLAR_AU3` directly in the source.
///
/// Keeping the conversion at the construction boundary rather than
/// inside this function lets `Body::density` mean exactly one thing
/// (coherent M_☉/AU³ in solar_au) regardless of how the body was
/// built.
pub fn radius_from_density_mass(density: f64, mass: f64) -> f64 {
    let vol = mass / density.max(1e-30);
    sphere_radius_from_volume(vol)
}

#[inline]
pub(crate) fn sphere_radius_from_volume(volume: f64) -> f64 {
    ((3.0 * volume) / (4.0 * PI)).cbrt()
}

/// Solar radius in astronomical units. Used internally by
/// [`Body::from_preset`] to convert `physical_radius` (in simulation
/// length units, i.e. AU under `solar_au`) to solar radii for the
/// luminosity model.
const SOLAR_RADIUS_AU: f64 = 0.00465047;

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unsoftened_zeroes_softening() {
        let b = Body::rocky(1.0);
        assert!(b.softening > 0.0, "default softening should be nonzero");
        let b = b.unsoftened();
        assert_eq!(b.softening, 0.0);
    }

    #[test]
    fn fluent_builder_produces_expected_state() {
        let b = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
        assert_eq!(b.pos_x, 1.0);
        assert_eq!(b.pos_y, 0.0);
        assert_eq!(b.pos_z, 0.0);
        assert_eq!(b.vel_x, 0.0);
        assert_eq!(b.vel_y, 1.0);
        assert_eq!(b.vel_z, 0.0);
        assert_eq!(b.mass, 3e-6);
    }

    #[test]
    fn material_constructors_default_z_and_vz_to_zero() {
        let b = Body::star(1.0);
        assert_eq!(b.pos_z, 0.0);
        assert_eq!(b.vel_z, 0.0);
    }

    #[test]
    fn at_leaves_z_untouched() {
        let b = Body::rocky(1.0).at_3d(0.0, 0.0, 5.0).at(1.0, 2.0);
        assert_eq!(b.pos_x, 1.0);
        assert_eq!(b.pos_y, 2.0);
        assert_eq!(b.pos_z, 5.0);
    }

    #[test]
    fn at_3d_sets_all_three_components() {
        let b = Body::rocky(1.0).at_3d(1.0, 2.0, 3.0);
        assert_eq!(b.pos_x, 1.0);
        assert_eq!(b.pos_y, 2.0);
        assert_eq!(b.pos_z, 3.0);
    }

    #[test]
    fn with_velocity_leaves_vz_untouched() {
        let b = Body::rocky(1.0).with_velocity_3d(0.0, 0.0, 7.0).with_velocity(1.0, 2.0);
        assert_eq!(b.vel_x, 1.0);
        assert_eq!(b.vel_y, 2.0);
        assert_eq!(b.vel_z, 7.0);
    }

    #[test]
    fn with_velocity_3d_sets_all_three_components() {
        let b = Body::rocky(1.0).with_velocity_3d(1.0, 2.0, 3.0);
        assert_eq!(b.vel_x, 1.0);
        assert_eq!(b.vel_y, 2.0);
        assert_eq!(b.vel_z, 3.0);
    }

    // ── Density / luminosity propagation from presets ────────────────────────

    #[test]
    fn rocky_preset_uses_material_default_density() {
        // Anchor mass is 1 M_⊕ ≈ 3.0034e-6 M_☉. The preset's EOS pivot
        // lands on ρ₀ = 5514 kg/m³ (Earth's bulk density); `from_preset`
        // crosses to coherent M_☉/AU³ at construction, so the stored
        // value is `5514 × KG_M3_TO_SOLAR_AU3 ≈ 9.28 × 10⁶`.
        let earth = Body::rocky(3.0034e-6);
        let expected = 5514.0 * crate::templates::builders::KG_M3_TO_SOLAR_AU3;
        assert!((earth.density - expected).abs() / expected < 1e-3);
    }

    #[test]
    fn icy_density_differs_from_rocky() {
        let rocky = Body::rocky(1.0);
        let icy = Body::icy(1.0);
        assert_ne!(rocky.density, icy.density);
    }

    #[test]
    fn star_preset_sets_luminosity() {
        let sun = Body::star(1.0);
        // Mass=1 in solar units → luminosity ≈ 1 solar luminosity.
        assert!((sun.luminosity - 1.0).abs() < 0.05, "L = {}", sun.luminosity);
    }

    #[test]
    fn rocky_preset_has_zero_luminosity() {
        assert_eq!(Body::rocky(1.0).luminosity, 0.0);
    }

    #[test]
    fn is_luminous_false_for_planets() {
        assert!(!Body::rocky(1.0).is_luminous());
        assert!(!Body::gas_giant(317.0).is_luminous());
    }

    #[test]
    fn is_luminous_true_for_main_sequence_star() {
        assert!(Body::star(1.0).is_luminous());
    }

    // ── Q_pr propagation ─────────────────────────────────────────────────────

    #[test]
    fn dust_class_presets_have_positive_q_pr() {
        assert!(Body::asteroid(1e-4).q_pr > 0.0);
        assert!(Body::comet(1e-6).q_pr > 0.0);
        assert!(Body::icy(1.0).q_pr > 0.0);
    }

    #[test]
    fn massive_classes_have_zero_q_pr() {
        assert_eq!(Body::rocky(1.0).q_pr, 0.0);
        assert_eq!(Body::star(1.0).q_pr, 0.0);
    }

    #[test]
    fn with_q_pr_overrides_preset_default() {
        let custom = Body::rocky(1.0).with_q_pr(0.5);
        assert_eq!(custom.q_pr, 0.5);
    }

    // ── Low-level Body::new ──────────────────────────────────────────────────

    #[test]
    fn body_new_uses_explicit_density_and_zero_optionals() {
        let b = Body::new(1.0, 2000.0);
        assert_eq!(b.mass, 1.0);
        assert_eq!(b.density, 2000.0);
        assert_eq!(b.luminosity, 0.0);
        assert_eq!(b.q_pr, 0.0);
        assert_eq!(b.color, [180, 180, 180]);
    }

    #[test]
    fn body_new_recovers_radius_from_density() {
        let b = Body::new(1.0, 2000.0);
        let expected = radius_from_density_mass(2000.0, 1.0);
        assert_eq!(b.physical_radius, expected);
    }

    #[test]
    fn with_density_recomputes_radius() {
        let b = Body::rocky(1.0).with_density(10_000.0);
        let expected = radius_from_density_mass(10_000.0, b.mass);
        assert_eq!(b.physical_radius, expected);
    }
}
