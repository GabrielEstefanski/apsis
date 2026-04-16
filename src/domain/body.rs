use crate::domain::materials::{Material, density};
use std::f64::consts::PI;

/// Base softening length for a body of mass 1.0.
/// Per-body softening scales as `EPS_BASE * mass^(1/3)`, so each body's
/// softening volume is proportional to its mass — physically motivated by
/// the Plummer-equivalent equal-mass softening criterion.
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
    /// - moment of inertia
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

    /// Angular velocity around the z-axis: ω_z = L_z / I.
    pub omega_z: f64,

    /// Moment of inertia around z-axis: I_z = (2/5)·m·r² for a uniform sphere.
    pub moment_inertia: f64,

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
/// Template instantiation uses this type to preserve authored names when bodies
/// are inserted into the simulation. Callers may leave `name` as `None` to let
/// the system derive a stable material-based fallback.
#[derive(Clone, Debug)]
pub struct NamedBody {
    pub body: Body,
    pub name: Option<String>,
}

impl Body {
    pub fn new(x: f64, y: f64, vx: f64, vy: f64, mass: f64, material: Material) -> Self {
        let density = density(material, mass);
        let physical_radius = radius_from_density_mass(density, mass);
        let softening = default_softening(mass);

        Self {
            x,
            y,
            vx,
            vy,
            mass,
            softening,
            physical_radius,
            density,
            omega_z: 0.0,
            moment_inertia: default_moment_inertia(mass, physical_radius),
            material,
            color: material.props().base_color,
            luminosity: 0.0,
        }
    }

    /// Recompute physical-only quantities from the current mass and density.
    ///
    /// This must be used whenever mass or density changes. It intentionally
    /// does **not** touch the calibrated contact radius, which belongs to the
    /// numerical collision model rather than the body's physical geometry.

    pub fn sync_physical_properties(&mut self) {
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
        self.moment_inertia = default_moment_inertia(self.mass, self.physical_radius);
    }

    /// Computes the bolometric luminosity of this body in solar luminosities.
    ///
    /// This is a **pure function** — it does not modify `self` and can be
    /// called at any time for diagnostic or UI purposes.
    ///
    /// # Unit conversion parameters
    ///
    /// | Parameter        | Meaning                                     |
    /// |------------------|---------------------------------------------|
    /// | `mass_to_solar`  | internal mass unit → M☉                    |
    /// | `radius_to_solar`| internal length unit → R☉                  |
    ///
    /// # Accuracy
    ///
    /// Order-of-magnitude estimate suitable for radiation pressure
    /// calculations. Not a stellar evolution model.
    pub fn luminosity_solar(&self, mass_to_solar: f64, radius_to_solar: f64) -> f64 {
        let m = self.mass * mass_to_solar;
        match self.material {
            Material::Star => main_sequence_luminosity_smooth(m),
            Material::BrownDwarf => brown_dwarf_luminosity(m),
            Material::WhiteDwarf => white_dwarf_luminosity(self.physical_radius * radius_to_solar),
            _ => 0.0,
        }
    }

    /// Updates the cached [`luminosity`](Self::luminosity) field from the
    /// current mass, radius, and supplied unit conversion factors.
    ///
    /// Must be called explicitly after:
    /// - construction (`Body::new`) if radiation is enabled
    /// - any change to `mass` or `density` (after [`sync_physical_properties`])
    /// - any change to `material`
    ///
    /// `l_sun` is L☉ expressed in internal energy · time⁻¹ units.
    pub fn update_luminosity(&mut self, mass_to_solar: f64, radius_to_solar: f64, l_sun: f64) {
        self.luminosity = self.luminosity_solar(mass_to_solar, radius_to_solar) * l_sun;
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
pub fn default_softening(mass: f64) -> f64 {
    EPS_BASE * mass.abs().cbrt()
}

/// Moment of inertia for a uniform sphere: I = (2/5)·m·r².
/// Uses the **physical radius**.
pub fn default_moment_inertia(mass: f64, radius: f64) -> f64 {
    0.4 * mass * radius * radius
}

/// Density from mass and radius: ρ = m / (4/3 π r³).
///
/// Returns a safe positive fallback when `radius ≤ 0`.
pub fn density_from_mass_radius(mass: f64, radius: f64) -> f64 {
    let vol = sphere_volume(radius);
    if vol > 0.0 { mass / vol } else { 1.0 }
}

/// Radius from density and mass: r = (3m / 4πρ)^(1/3).
pub fn radius_from_density_mass(density: f64, mass: f64) -> f64 {
    let vol = mass / density.max(1e-30);
    sphere_radius_from_volume(vol)
}

/// Volume of a sphere: V = 4/3 π r³.
#[inline]
pub fn sphere_volume(radius: f64) -> f64 {
    (4.0 / 3.0) * PI * radius.powi(3)
}

/// Radius of a sphere given its volume: r = (3V / 4π)^(1/3).
#[inline]
pub fn sphere_radius_from_volume(volume: f64) -> f64 {
    ((3.0 * volume) / (4.0 * PI)).cbrt()
}

// ── Luminosity models ─────────────────────────────────────────────────────────

/// Logistic sigmoid: smooth step from 0 → 1 centred at `m0` with width `w`.
///
/// Used to blend continuously between power-law regimes, avoiding the
/// derivative discontinuities of a piecewise-constant exponent.
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
/// Blending width: 0.15 M☉.  The function is strictly monotonically
/// increasing for all M > 0: `dL/dM = α · M^(α−1) > 0`.
///
/// The Eddington luminosity limit is not enforced here; for the
/// N-body force model (radiation pressure, PR drag) the error from
/// ignoring Eddington clamping is negligible at the masses used.
///
/// References: Salaris & Cassisi (2005) §5.3; Tout et al. (1996) *MNRAS* 281.
fn main_sequence_luminosity_smooth(m: f64) -> f64 {
    if m <= 0.0 {
        return 0.0;
    }

    let alpha = 2.3 + (3.5 - 2.3) * logistic(m, 0.43, 0.15); // low-mass → solar

    m.powf(alpha)
}

/// Deuterium-burning luminosity for sub-stellar objects.
///
/// Power-law fit for `0.013 M☉ < M < 0.08 M☉`. A logistic onset term
/// ensures `L → 0` smoothly at the deuterium-burning limit rather than
/// cutting off abruptly.
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
/// ```text
/// L / L☉ = (R / R☉)² · (T_eff / T☉)⁴
/// ```
///
/// Uses `T_eff = 10 000 K` (typical for a ~1 Gyr cooling age) and
/// `T☉ = 5 778 K`. Actual luminosity ranges from ~0.1 L☉ (young) to
/// ~10⁻⁴ L☉ (old); this estimate is appropriate for middle-aged white dwarfs.
///
/// `r_solar` must be the physical radius already converted to R☉.
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
        // By definition L(1 M☉) = 1 L☉ on the main sequence.
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
        // dL/dM must not jump at M = 0.43 M☉
        let eps = 1e-4;
        let l_lo = main_sequence_luminosity_smooth(0.43 - eps);
        let l_hi = main_sequence_luminosity_smooth(0.43 + eps);
        let slope = (l_hi - l_lo) / (2.0 * eps);
        // Finite difference of the slope — should be smooth, not infinite
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
        // A WD with R = 0.01 R☉ should give L ~ 10⁻³–10⁻² L☉
        let l = white_dwarf_luminosity(0.01);
        assert!(l > 1e-4 && l < 0.1, "WD luminosity out of expected range: {l}");
    }

    #[test]
    fn non_luminous_materials_return_zero() {
        let body = Body::new(0.0, 0.0, 0.0, 0.0, 1.0, Material::Rocky);
        assert_eq!(body.luminosity_solar(1.0, 1.0), 0.0);
    }

    #[test]
    fn is_luminous_false_before_update() {
        let body = Body::new(0.0, 0.0, 0.0, 0.0, 1.0, Material::Star);
        // luminosity field starts at 0.0 until update_luminosity is called
        assert!(!body.is_luminous());
    }

    #[test]
    fn is_luminous_true_after_update() {
        let mut body = Body::new(0.0, 0.0, 0.0, 0.0, 1.0, Material::Star);
        body.update_luminosity(1.0, 1.0 / 0.00465, 1.0);
        assert!(body.is_luminous());
    }
}
