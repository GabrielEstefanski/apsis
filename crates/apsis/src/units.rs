//! Physical unit systems for the simulation.
//!
//! [`UnitSystem`] is the closed contract that turns dimensionless body
//! state (`x`, `vy`, `mass`, `dt`) into a physically meaningful run.
//! Every [`crate::core::system::System`] is constructed against an
//! explicit `UnitSystem`; there is no implicit default and no setter
//! after construction.
//!
//! `G` is always derived from the three SI scales by
//! `G_code = G_SI · M[kg] · T[s]² / L[m]³`; nothing here hardcodes
//! `4π²` or `1`. There is no dimensional checking — passing a value
//! in the wrong unit is a silent physical error, matching REBOUND's
//! trade-off (no per-multiplication overhead in the integrator loop).

use std::fmt;

// ── SI constants ──────────────────────────────────────────────────────────────

/// Newtonian gravitational constant in SI (m³ kg⁻¹ s⁻²). CODATA 2018.
pub const G_SI: f64 = 6.674_30e-11;

/// Astronomical unit in metres. Exact since the IAU 2012 redefinition.
pub const AU_M: f64 = 1.495_978_707e11;

/// Julian year in seconds (365.25 × 86 400).
pub const YR_S: f64 = 3.155_76e7;

/// Solar mass in kilograms (IAU 2015 nominal).
pub const MSUN_KG: f64 = 1.988_92e30;

/// One centimetre in metres.
pub const CM_M: f64 = 1.0e-2;

/// One gram in kilograms.
pub const G_KG: f64 = 1.0e-3;

// ── UnitSystem ────────────────────────────────────────────────────────────────

/// A closed system of units for length, time, and mass.
///
/// Constructable only via the named factories
/// ([`UnitSystem::canonical`], [`UnitSystem::si`], [`UnitSystem::solar`],
/// [`UnitSystem::cgs`]) or [`UnitSystem::custom`]. There is no `Default`,
/// no public field access, and no public setter — the contract is that
/// once a `UnitSystem` is chosen at [`System::new`] time, it cannot be
/// changed without rebuilding the `System` from scratch.
///
/// [`System::new`]: crate::core::system::System::new
#[derive(Clone, Copy, Debug)]
pub struct UnitSystem {
    length_m: f64,
    time_s: f64,
    mass_kg: f64,

    // Display-only labels — ignored by `PartialEq`.
    length_label: &'static str,
    time_label: &'static str,
    mass_label: &'static str,
}

impl UnitSystem {
    // ── Named factories ──────────────────────────────────────────────────

    /// Hénon-style canonical N-body units: `G = 1` by construction.
    ///
    /// Length and time scales are nominally `1` in SI; the mass scale
    /// absorbs the `1/G_SI` factor required to satisfy
    /// `G = G_SI · M · T² / L³ = 1`. This matches the implicit
    /// default of REBOUND when no units are specified, and the
    /// convention of stellar-dynamics literature (Aarseth, Hénon).
    pub const fn canonical() -> Self {
        Self {
            length_m: 1.0,
            time_s: 1.0,
            mass_kg: 1.0 / G_SI,
            length_label: "L",
            time_label: "T",
            mass_label: "M",
        }
    }

    /// Alias for [`canonical`](Self::canonical) using the literature name.
    /// Use this in code that targets stellar-dynamics readers; the
    /// values and behaviour are identical.
    pub const fn henon() -> Self {
        Self::canonical()
    }

    /// SI units: metre, second, kilogram. `G ≈ 6.674e-11`.
    pub const fn si() -> Self {
        Self {
            length_m: 1.0,
            time_s: 1.0,
            mass_kg: 1.0,
            length_label: "m",
            time_label: "s",
            mass_label: "kg",
        }
    }

    /// Solar-system IAU units: astronomical unit, Julian year, solar
    /// mass. The derived `G` is `≈ 39.478`, the IAU approximation to
    /// `4π²` that satisfies Kepler's third law for Earth's orbit by
    /// construction.
    ///
    /// Distinct from [`solar_canonical`](Self::solar_canonical), which
    /// uses `T = year/(2π)` to make `G = 1` exactly. Pick `solar` for
    /// IAU compatibility (REBOUND default with `G = 4π²`); pick
    /// `solar_canonical` for Hénon-style normalisation (REBOUND with
    /// `G = 1`, the apsis-1pn validation portfolio convention).
    pub const fn solar() -> Self {
        Self {
            length_m: AU_M,
            time_s: YR_S,
            mass_kg: MSUN_KG,
            length_label: "AU",
            time_label: "yr",
            mass_label: "Msun",
        }
    }

    /// Solar-system canonical (Hénon-normalised) units: astronomical
    /// unit, **Gaussian time** = `sqrt(AU³ / (G_SI · M_sun))`, solar
    /// mass. The time scale is chosen so the derived `G` equals `1`
    /// **exactly** by construction — `G_code = G_SI · M · T² / L³ = 1`.
    /// Standard convention for stellar-dynamics literature (Aarseth)
    /// when fixing both physical units AND `G = 1`, and the unit
    /// system the apsis-1pn validation portfolio (Mercury 1PN gate,
    /// long-horizon experiments) runs in.
    ///
    /// The Gaussian time unit numerically differs from `YR_S/(2π)`
    /// (IAU julian year over 2π) by ~0.009 % — the historical
    /// astrodynamics gap between the IAU-defined year and the year
    /// implied by the Gaussian gravitational constant. The Gaussian
    /// definition is what yields G = 1 exactly; the IAU julian
    /// definition only approximately. Distinct from [`solar`](Self::solar),
    /// which uses the IAU julian year directly and yields `G ≈ 4π²`.
    ///
    /// Not `const fn` because `f64::sqrt` is not stable in const
    /// context; the value is otherwise immutable.
    pub fn solar_canonical() -> Self {
        // T = sqrt(L³ / (G · M)) so that G_code = G · M · T² / L³ = 1.
        let l3 = AU_M * AU_M * AU_M;
        let gm = G_SI * MSUN_KG;
        let t_gaussian = (l3 / gm).sqrt();
        Self {
            length_m: AU_M,
            time_s: t_gaussian,
            mass_kg: MSUN_KG,
            length_label: "AU",
            time_label: "T_G",
            mass_label: "Msun",
        }
    }

    /// CGS units: centimetre, second, gram. `G ≈ 6.674e-8` (cm³ g⁻¹ s⁻²).
    pub const fn cgs() -> Self {
        Self {
            length_m: CM_M,
            time_s: 1.0,
            mass_kg: G_KG,
            length_label: "cm",
            time_label: "s",
            mass_label: "g",
        }
    }

    /// Build a `UnitSystem` from explicit SI scales.
    ///
    /// All three scales must be strictly positive and finite. Zero,
    /// negative, infinite, and NaN values are rejected at the boundary
    /// because they would cause `g()` to return a non-finite value
    /// that would only manifest as a numerical explosion deep inside
    /// the integrator — far from the cause.
    ///
    /// Labels are generic (`"length"`, `"time"`, `"mass"`); the named
    /// factories are the right tool when a literature shorthand
    /// (`AU`, `yr`, …) is meaningful.
    pub fn custom(length_m: f64, time_s: f64, mass_kg: f64) -> Result<Self, UnitError> {
        if !length_m.is_finite() || length_m <= 0.0 {
            return Err(UnitError::InvalidLength(length_m));
        }
        if !time_s.is_finite() || time_s <= 0.0 {
            return Err(UnitError::InvalidTime(time_s));
        }
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(UnitError::InvalidMass(mass_kg));
        }
        Ok(Self {
            length_m,
            time_s,
            mass_kg,
            length_label: "length",
            time_label: "time",
            mass_label: "mass",
        })
    }

    // ── Scale accessors ──────────────────────────────────────────────────

    /// SI metres per code-unit length.
    #[inline]
    pub fn length_scale_si(&self) -> f64 {
        self.length_m
    }

    /// SI seconds per code-unit time.
    #[inline]
    pub fn time_scale_si(&self) -> f64 {
        self.time_s
    }

    /// SI kilograms per code-unit mass.
    #[inline]
    pub fn mass_scale_si(&self) -> f64 {
        self.mass_kg
    }

    // ── Derived gravitational constant ───────────────────────────────────

    /// Newtonian gravitational constant in this system's canonical units.
    ///
    /// Computed from the SI scales by
    /// `G_code = G_SI · mass_scale · time_scale² / length_scale³`.
    /// Never read from a hardcoded literature value — `solar().g()`
    /// returns `≈ 39.478` (the IAU approximation to `4π²`) by
    /// derivation, not by definition. `solar_canonical().g()`
    /// returns `1.0 ± 1 ULP` rather than literal `1.0` because
    /// `time_s = sqrt(AU³/(G·M))` is itself a derivation.
    #[inline]
    pub fn g(&self) -> f64 {
        G_SI * self.mass_kg * self.time_s * self.time_s
            / (self.length_m * self.length_m * self.length_m)
    }

    // ── Explicit conversions ─────────────────────────────────────────────

    /// Convert a length expressed in this system's canonical units to SI metres.
    #[inline]
    pub fn length_to_si(&self, x: f64) -> f64 {
        x * self.length_m
    }

    /// Convert a length expressed in SI metres to this system's canonical units.
    #[inline]
    pub fn length_from_si(&self, x: f64) -> f64 {
        x / self.length_m
    }

    /// Convert a duration expressed in this system's canonical units to SI seconds.
    #[inline]
    pub fn time_to_si(&self, x: f64) -> f64 {
        x * self.time_s
    }

    /// Convert a duration expressed in SI seconds to this system's canonical units.
    #[inline]
    pub fn time_from_si(&self, x: f64) -> f64 {
        x / self.time_s
    }

    /// Convert a mass expressed in this system's canonical units to SI kilograms.
    #[inline]
    pub fn mass_to_si(&self, x: f64) -> f64 {
        x * self.mass_kg
    }

    /// Convert a mass expressed in SI kilograms to this system's canonical units.
    #[inline]
    pub fn mass_from_si(&self, x: f64) -> f64 {
        x / self.mass_kg
    }

    // ── Display labels ───────────────────────────────────────────────────

    /// Display symbol for the length axis (`"AU"`, `"m"`, `"cm"`, ...).
    #[inline]
    pub fn length_label(&self) -> &'static str {
        self.length_label
    }

    /// Display symbol for the time axis (`"yr"`, `"s"`, ...).
    #[inline]
    pub fn time_label(&self) -> &'static str {
        self.time_label
    }

    /// Display symbol for the mass axis (`"Msun"`, `"kg"`, `"g"`, ...).
    #[inline]
    pub fn mass_label(&self) -> &'static str {
        self.mass_label
    }
}

// Manual impl: only SI scales count; labels are display metadata so a
// `custom(AU_M, YR_S, MSUN_KG)` compares equal to `solar()`.
impl PartialEq for UnitSystem {
    fn eq(&self, other: &Self) -> bool {
        self.length_m == other.length_m
            && self.time_s == other.time_s
            && self.mass_kg == other.mass_kg
    }
}

impl fmt::Display for UnitSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "UnitSystem(L=1 {} = {:.4e} m, T=1 {} = {:.4e} s, M=1 {} = {:.4e} kg, G={:.4e})",
            self.length_label,
            self.length_m,
            self.time_label,
            self.time_s,
            self.mass_label,
            self.mass_kg,
            self.g(),
        )
    }
}

// ── Errors ────────────────────────────────────────────────────────────────────

/// Failure to construct a [`UnitSystem`] from explicit SI scales.
///
/// Each variant carries the offending value so a downstream binding
/// (Python, CLI parser, config loader) can include it verbatim in the
/// error message it surfaces to the user.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitError {
    /// `length_m` was zero, negative, infinite, or NaN.
    InvalidLength(f64),
    /// `time_s` was zero, negative, infinite, or NaN.
    InvalidTime(f64),
    /// `mass_kg` was zero, negative, infinite, or NaN.
    InvalidMass(f64),
}

impl fmt::Display for UnitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(v) => {
                write!(f, "length scale must be a strictly positive finite f64, got {v}")
            },
            Self::InvalidTime(v) => {
                write!(f, "time scale must be a strictly positive finite f64, got {v}")
            },
            Self::InvalidMass(v) => {
                write!(f, "mass scale must be a strictly positive finite f64, got {v}")
            },
        }
    }
}

impl std::error::Error for UnitError {}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// `canonical()` must satisfy `G = 1` exactly. Any drift here
    /// silently rescales every gravitational interaction, so this is
    /// the cheapest invariant to assert at the unit-system level.
    #[test]
    fn canonical_g_is_exactly_one() {
        let u = UnitSystem::canonical();
        assert_eq!(u.g(), 1.0, "canonical units must have G = 1 exactly");
    }

    /// `henon()` is documented as a literature alias for `canonical()`.
    /// They must compare equal under `PartialEq` so callers can use
    /// either name interchangeably without surprises.
    #[test]
    fn henon_equals_canonical() {
        assert_eq!(UnitSystem::henon(), UnitSystem::canonical());
    }

    /// `si()` must reproduce `G_SI` exactly — the round-trip identity
    /// for the trivial system where the SI scales are all `1`.
    #[test]
    fn si_g_equals_g_si() {
        let u = UnitSystem::si();
        assert_eq!(u.g(), G_SI);
    }

    /// `solar()`'s derived `G` must land near `4π²` (the dimensionless
    /// form of Kepler's third law for the Sun-Earth system). Tolerance
    /// is 1% — `4π²` is an idealisation that assumes the Earth's orbit
    /// is exactly circular at exactly `1 AU` with exactly `1 yr`
    /// period. Real CODATA `G_SI` × IAU `MSUN_KG` × Julian `YR_S` /
    /// `AU_M³` gives `≈ 39.487`, about 0.02% off from `4π² ≈ 39.478`.
    /// The test pins the order-of-magnitude / sign / shape, not the
    /// idealised value.
    #[test]
    fn solar_g_approximates_four_pi_squared() {
        let u = UnitSystem::solar();
        let four_pi_sq = 4.0 * std::f64::consts::PI * std::f64::consts::PI;
        let rel_err = (u.g() - four_pi_sq).abs() / four_pi_sq;
        assert!(
            rel_err < 1e-2,
            "solar G ({}) deviates too far from 4π² ({four_pi_sq}); rel_err = {rel_err:.3e}",
            u.g()
        );
    }

    /// Conversions must round-trip to f64 round-off. Fundamental
    /// invariant: applying `to_si` then `from_si` (and vice versa)
    /// returns the original value within ULP.
    #[test]
    fn conversions_round_trip() {
        let u = UnitSystem::solar();
        let original = 1.234_567_8;
        for converted in [
            u.length_from_si(u.length_to_si(original)),
            u.time_from_si(u.time_to_si(original)),
            u.mass_from_si(u.mass_to_si(original)),
        ] {
            assert!((converted - original).abs() < 1e-15);
        }
    }

    /// `length_to_si(1.0)` is the length scale itself — semantic check.
    #[test]
    fn length_to_si_unit_input_is_scale() {
        let u = UnitSystem::solar();
        assert_eq!(u.length_to_si(1.0), AU_M);
        assert_eq!(u.time_to_si(1.0), YR_S);
        assert_eq!(u.mass_to_si(1.0), MSUN_KG);
    }

    /// `custom()` must reject zero, negative, infinite, and NaN scales
    /// at the boundary — the failure mode otherwise is a `NaN` or `Inf`
    /// `g()` that explodes inside the integrator far from the cause.
    #[test]
    fn custom_rejects_invalid_scales() {
        for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(UnitSystem::custom(bad, 1.0, 1.0), Err(UnitError::InvalidLength(_))));
            assert!(matches!(UnitSystem::custom(1.0, bad, 1.0), Err(UnitError::InvalidTime(_))));
            assert!(matches!(UnitSystem::custom(1.0, 1.0, bad), Err(UnitError::InvalidMass(_))));
        }
    }

    /// `custom()` with the same SI scales as a named factory must
    /// compare equal to it. Labels diverge (named vs. generic) but
    /// `PartialEq` ignores labels by design.
    #[test]
    fn custom_with_solar_scales_equals_solar() {
        let custom = UnitSystem::custom(AU_M, YR_S, MSUN_KG).unwrap();
        let named = UnitSystem::solar();
        assert_eq!(custom, named);
        // But labels diverge — this is the metadata path.
        assert_ne!(custom.length_label(), named.length_label());
    }

    /// `Display` includes the chosen labels, the SI scales, and `G` —
    /// the minimum information required to reproduce a run from a log line.
    #[test]
    fn display_carries_labels_scales_and_g() {
        let s = format!("{}", UnitSystem::solar());
        assert!(s.contains("AU"), "display should mention AU label: {s}");
        assert!(s.contains("yr"), "display should mention yr label: {s}");
        assert!(s.contains("Msun"), "display should mention Msun label: {s}");
        assert!(s.contains("G="), "display should mention G value: {s}");
    }
}
