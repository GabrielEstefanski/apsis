use crate::domain::materials::{Material, density};
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

#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub mass: f64,

    /// Gravitational softening length ε for this body.
    ///
    /// Pairwise: ε²_ij = (ε²_i + ε²_j) / 2.
    ///
    /// Calibrated by `System::calibrate_softening`.
    pub softening: f64,

    /// True physical radius derived from mass and density.
    ///
    /// This represents the actual size of the body and is used for:
    /// - energy calculations (e.g. disruption threshold Q*)
    /// - physically meaningful scaling
    ///
    /// Unlike `radius`, this value is **never modified by calibration**.
    pub physical_radius: f64,

    /// Bulk density of the body: ρ = m / V, V = 4/3 π r³.
    ///
    /// This is the **primary size property** of a body — the physical radius
    /// is derived from it via `r = (3m / 4πρ)^(1/3)`.
    ///
    /// This value is invariant during simulation except for merge/fragmentation
    /// events where material composition changes.
    pub density: f64,

    /// Astrophysical material class.
    pub material: Material,

    /// Display colour [R, G, B].
    pub color: [u8; 3],

    /// Bolometric luminosity in internal energy · time⁻¹ units.
    ///
    /// This field is **not** computed automatically on construction — it starts
    /// at `0.0` and must be populated explicitly via [`update_luminosity`]
    /// before any radiation calculation.
    ///
    /// [`RadiationField::from_bodies`] reads this field directly; calling
    /// `from_bodies` without first calling `update_luminosity` on all luminous
    /// bodies will silently produce no radiation sources.
    pub luminosity: f64,
}

/// Body payload with an optional explicit display name.
///
/// Returned by [`Body::named`] and by the template catalog. Consumed by
/// [`System::add_named_body`](crate::core::system::System::add_named_body)
/// to preserve authored names; otherwise the system derives a stable
/// material-based fallback.
#[derive(Clone, Debug)]
pub struct NamedBody {
    pub body: Body,
    pub name: Option<String>,
}

impl Body {
    // ── Material constructors ─────────────────────────────────────────────────
    //
    // Each constructor creates a body at the origin, at rest, with physical
    // properties derived from `mass` and the canonical material profile.
    // Position and velocity are set via the fluent builder methods below.

    /// Star — main-sequence luminous body. Default density, luminous material.
    pub fn star(mass: f64) -> Self {
        Self::from_material(mass, Material::Star)
    }

    /// Brown dwarf — sub-stellar, deuterium-burning regime.
    pub fn brown_dwarf(mass: f64) -> Self {
        Self::from_material(mass, Material::BrownDwarf)
    }

    /// White dwarf — compact stellar remnant.
    pub fn white_dwarf(mass: f64) -> Self {
        Self::from_material(mass, Material::WhiteDwarf)
    }

    /// Gas giant — Jupiter-class hydrogen/helium envelope.
    pub fn gas_giant(mass: f64) -> Self {
        Self::from_material(mass, Material::Gas)
    }

    /// Ice giant — Neptune-class water/methane envelope.
    pub fn ice_giant(mass: f64) -> Self {
        Self::from_material(mass, Material::IceGiant)
    }

    /// Rocky body — terrestrial planet or large rocky satellite.
    pub fn rocky(mass: f64) -> Self {
        Self::from_material(mass, Material::Rocky)
    }

    /// Icy body — water-dominated composition (outer satellites, KBOs).
    pub fn icy(mass: f64) -> Self {
        Self::from_material(mass, Material::Icy)
    }

    /// Asteroid — rocky minor body.
    pub fn asteroid(mass: f64) -> Self {
        Self::from_material(mass, Material::Asteroid)
    }

    /// Comet — volatile-rich minor body.
    pub fn comet(mass: f64) -> Self {
        Self::from_material(mass, Material::Comet)
    }

    /// Body with an explicit material.
    ///
    /// Prefer the material-named constructors ([`star`](Self::star),
    /// [`rocky`](Self::rocky), …) for readability; this is the escape hatch
    /// when the material is computed programmatically.
    pub fn of(mass: f64, material: Material) -> Self {
        Self::from_material(mass, material)
    }

    fn from_material(mass: f64, material: Material) -> Self {
        let density = density(material, mass);
        let physical_radius = radius_from_density_mass(density, mass);
        let softening = default_softening(mass);

        Self {
            x: 0.0,
            y: 0.0,
            vx: 0.0,
            vy: 0.0,
            mass,
            softening,
            physical_radius,
            density,
            material,
            color: material.props().base_color,
            luminosity: 0.0,
        }
    }

    // ── Fluent builder ────────────────────────────────────────────────────────
    //
    // Each method consumes and returns `Self`, so they chain naturally:
    //
    //     Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0)

    /// Position in simulation coordinates.
    #[inline]
    #[must_use]
    pub fn at(mut self, x: f64, y: f64) -> Self {
        self.x = x;
        self.y = y;
        self
    }

    /// Velocity in simulation coordinates.
    #[inline]
    #[must_use]
    pub fn with_velocity(mut self, vx: f64, vy: f64) -> Self {
        self.vx = vx;
        self.vy = vy;
        self
    }

    /// Override the material-default density. Radius is recomputed to match.
    #[inline]
    #[must_use]
    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
        self
    }

    /// Attach an explicit display name, producing a [`NamedBody`] consumable
    /// by [`System::add_named_body`](crate::core::system::System::add_named_body).
    #[inline]
    #[must_use]
    pub fn named(self, name: impl Into<String>) -> NamedBody {
        NamedBody {
            body: self,
            name: Some(name.into()),
        }
    }

    // ── Mutators ──────────────────────────────────────────────────────────────

    /// Recompute physical-only quantities from the current mass and density.
    ///
    /// Must be used whenever `mass` or `density` is mutated in place
    /// (e.g. via direct field assignment on a `&mut Body`). It intentionally
    /// does **not** touch the calibrated contact radius, which belongs to the
    /// numerical collision model rather than the body's physical geometry.
    pub fn sync_physical_properties(&mut self) {
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
    }

    /// Updates the cached [`luminosity`](Self::luminosity) field from the
    /// current mass, radius, and supplied unit conversion factors.
    ///
    /// Must be called explicitly after:
    /// - construction if radiation is enabled
    /// - any change to `mass` or `density` (after [`sync_physical_properties`])
    /// - any change to `material`
    ///
    /// `l_sun` is L☉ expressed in internal energy · time⁻¹ units.
    pub(crate) fn update_luminosity(
        &mut self,
        mass_to_solar: f64,
        radius_to_solar: f64,
        l_sun: f64,
    ) {
        self.luminosity = self.luminosity_solar(mass_to_solar, radius_to_solar) * l_sun;
    }

    pub(crate) fn luminosity_solar(&self, mass_to_solar: f64, radius_to_solar: f64) -> f64 {
        let m = self.mass * mass_to_solar;
        match self.material {
            Material::Star => main_sequence_luminosity_smooth(m),
            Material::BrownDwarf => brown_dwarf_luminosity(m),
            Material::WhiteDwarf => white_dwarf_luminosity(self.physical_radius * radius_to_solar),
            _ => 0.0,
        }
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

/// Radius from density and mass: r = (3m / 4πρ)^(1/3).
pub fn radius_from_density_mass(density: f64, mass: f64) -> f64 {
    let vol = mass / density.max(1e-30);
    sphere_radius_from_volume(vol)
}

#[inline]
pub(crate) fn sphere_radius_from_volume(volume: f64) -> f64 {
    ((3.0 * volume) / (4.0 * PI)).cbrt()
}

// ── Luminosity models ─────────────────────────────────────────────────────────

/// Logistic sigmoid: smooth step from 0 → 1 centred at `m0` with width `w`.
#[inline]
fn logistic(m: f64, m0: f64, w: f64) -> f64 {
    1.0 / (1.0 + ((m0 - m) / w).exp())
}

/// Main-sequence mass–luminosity relation with a continuously differentiable
/// exponent.
///
/// The exponent `α(M)` blends smoothly across two physical regimes:
///
/// | Regime       | Mass range  | α   | Dominant physics          |
/// |--------------|-------------|-----|---------------------------|
/// | Low-mass     | M ≲ 0.43 M☉ | 2.3 | Fully convective interior |
/// | Solar-type   | M ≳ 0.43 M☉ | 3.5 | Radiative core            |
///
/// References: Salaris & Cassisi (2005) §5.3; Tout et al. (1996) *MNRAS* 281.
fn main_sequence_luminosity_smooth(m: f64) -> f64 {
    if m <= 0.0 {
        return 0.0;
    }

    let alpha = 2.3 + (3.5 - 2.3) * logistic(m, 0.43, 0.15);

    m.powf(alpha)
}

/// Deuterium-burning luminosity for sub-stellar objects.
///
/// Reference: Burrows et al. (1997) *ApJ* 491, 856.
fn brown_dwarf_luminosity(m: f64) -> f64 {
    if m <= 0.013 {
        return 0.0;
    }
    let onset = logistic(m, 0.013, 0.002);
    1e-3 * (m / 0.05).powi(2) * onset
}

/// Stefan–Boltzmann cooling luminosity for a white dwarf.
///
/// Reference: Koester & Chanmugam (1990) *Rep. Prog. Phys.* 53, 837.
fn white_dwarf_luminosity(r_solar: f64) -> f64 {
    const T_EFF: f64 = 10_000.0;
    const T_SUN: f64 = 5_778.0;
    let t_ratio = T_EFF / T_SUN;
    r_solar * r_solar * t_ratio.powi(4)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_mass_star_has_unit_luminosity() {
        let l = main_sequence_luminosity_smooth(1.0);
        assert!((l - 1.0).abs() < 0.01, "L(1 M☉) = {l}, expected ~1");
    }

    #[test]
    fn luminosity_increases_with_mass() {
        let l1 = main_sequence_luminosity_smooth(1.0);
        let l2 = main_sequence_luminosity_smooth(2.0);
        assert!(l2 > l1);
    }

    #[test]
    fn luminosity_continuous_at_low_mass_boundary() {
        let eps = 1e-4;
        let l_lo = main_sequence_luminosity_smooth(0.43 - eps);
        let l_hi = main_sequence_luminosity_smooth(0.43 + eps);
        let slope = (l_hi - l_lo) / (2.0 * eps);
        assert!(slope.is_finite(), "discontinuity at 0.43 M☉");
        assert!(slope > 0.0, "luminosity must increase with mass");
    }

    #[test]
    fn luminosity_continuous_at_eddington_boundary() {
        let eps = 1e-3;
        let l_lo = main_sequence_luminosity_smooth(50.0 - eps);
        let l_hi = main_sequence_luminosity_smooth(50.0 + eps);
        let slope = (l_hi - l_lo) / (2.0 * eps);
        assert!(slope.is_finite());
        assert!(slope > 0.0);
    }

    #[test]
    fn brown_dwarf_below_threshold_is_zero() {
        assert_eq!(brown_dwarf_luminosity(0.01), 0.0);
    }

    #[test]
    fn brown_dwarf_increases_with_mass() {
        let l1 = brown_dwarf_luminosity(0.03);
        let l2 = brown_dwarf_luminosity(0.07);
        assert!(l2 > l1);
    }

    #[test]
    fn white_dwarf_typical_luminosity_in_range() {
        let l = white_dwarf_luminosity(0.01);
        assert!(l > 1e-4 && l < 0.1, "WD luminosity out of expected range: {l}");
    }

    #[test]
    fn fluent_builder_produces_expected_state() {
        let b = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
        assert_eq!(b.x, 1.0);
        assert_eq!(b.y, 0.0);
        assert_eq!(b.vx, 0.0);
        assert_eq!(b.vy, 1.0);
        assert_eq!(b.mass, 3e-6);
        assert_eq!(b.material, Material::Rocky);
    }

    #[test]
    fn material_constructors_use_material_default_density() {
        let rocky = Body::rocky(1.0);
        let icy = Body::icy(1.0);
        assert!(rocky.density > 0.0);
        assert!(icy.density > 0.0);
        assert!(rocky.density != icy.density);
    }

    #[test]
    fn non_luminous_materials_return_zero() {
        let body = Body::rocky(1.0);
        assert_eq!(body.luminosity_solar(1.0, 1.0), 0.0);
    }

    #[test]
    fn is_luminous_false_before_update() {
        let body = Body::star(1.0);
        assert!(!body.is_luminous());
    }

    #[test]
    fn is_luminous_true_after_update() {
        let mut body = Body::star(1.0);
        body.update_luminosity(1.0, 1.0 / 0.00465, 1.0);
        assert!(body.is_luminous());
    }
}
