//! Numeric formatters for the Inspector.
//!
//! Single dispatcher [`format_value`] takes a value and a [`QuantityType`]
//! and returns the displayable `(value, unit)` pair. Each formatter is
//! deterministic and produces a fixed-precision output: the
//! [`crate::app::design::primitives::BloombergFlash`] mechanism compares
//! consecutive **formatted** outputs rather than raw values, so a
//! sub-display change does not trigger a flash.
//!
//! Conventions follow Vallado §3.5 and the REBOUND/SPICE output style:
//! mantissa width fixed, exponent always signed and zero-padded to two
//! digits, units abbreviated in lowercase except where standard demands
//! otherwise (`AU`).

use super::unit_strategy::{distance, time};

const NAN_DASH: &str = "—";

/// Quantity tag — drives unit selection and precision in [`format_value`].
///
/// The tag captures *what* a value represents (vector component vs scalar
/// distance, dynamic vs static angle, etc.). Reference frame information
/// (intrinsic vs camera-relative vs primary-relative) lives on the section
/// header / label, not here — see `feedback_scientific_app_idiom.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantityType {
    /// Position vector component (x, y, z). Always SI scientific in metres
    /// so cross-component comparison is honest.
    DistanceVector,
    /// Scalar distance (semi-major axis, pericenter, apocenter, body
    /// radius, distance to camera). Auto-unit km/AU/ly per magnitude.
    DistanceScalar,
    /// Velocity vector component (vx, vy, vz). SI scientific in m/s.
    VelocityVector,
    /// Mass — kg in scientific notation.
    Mass,
    /// Time / period — auto-unit days/years/scientific years.
    Time,
    /// Static angle (i, Ω, ω) — 2 decimal degrees.
    AngleStatic,
    /// Dynamic angle (ν, M, E) — 3 decimal degrees.
    AngleDynamic,
    /// Energy — J in scientific notation.
    Energy,
    /// Eccentricity (or any other dimensionless quantity needing 5
    /// decimals of precision).
    Eccentricity,
    /// Arcsecond — 1 decimal above 10″, 2 decimals below.
    Arcsecond,
}

/// Render a value through the [`QuantityType`] dispatcher. Returns the
/// formatted display string plus the unit string (empty when the unit
/// glyph is part of the value or the quantity is dimensionless).
pub fn format_value(value: f64, quantity: QuantityType) -> (String, &'static str) {
    match quantity {
        QuantityType::DistanceVector => (fmt_sci(value), "m"),
        QuantityType::DistanceScalar => fmt_distance_auto(value),
        QuantityType::VelocityVector => (fmt_sci(value), "m/s"),
        QuantityType::Mass => (fmt_sci(value), "kg"),
        QuantityType::Time => fmt_period_auto(value),
        QuantityType::AngleStatic => (fmt_deg_static(value), "°"),
        QuantityType::AngleDynamic => (fmt_deg_dynamic(value), "°"),
        QuantityType::Energy => (fmt_sci(value), "J"),
        QuantityType::Eccentricity => (fmt_decimal(value, 5), ""),
        QuantityType::Arcsecond => (fmt_arcsec(value), "″"),
    }
}

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Scientific notation with a 4-significant-figure mantissa and an
/// always-signed two-digit exponent. Output is constant 10 characters for
/// finite inputs (`+1.234e+10`, `-9.100e-03`); `NaN`/`±∞` collapse to `—`.
fn fmt_sci(x: f64) -> String {
    if !x.is_finite() {
        return NAN_DASH.to_owned();
    }
    if x == 0.0 {
        return "+0.000e+00".to_owned();
    }
    let sign = if x.is_sign_negative() { '-' } else { '+' };
    let abs = x.abs();
    let exp = abs.log10().floor() as i32;
    let mantissa = abs / 10f64.powi(exp);
    let (mantissa, exp) =
        if mantissa >= 9.9995 { (mantissa / 10.0, exp + 1) } else { (mantissa, exp) };
    let exp_sign = if exp < 0 { '-' } else { '+' };
    let exp_abs = exp.unsigned_abs();
    format!("{sign}{mantissa:.3}e{exp_sign}{exp_abs:02}")
}

/// Static angle — 2 decimals (no degree glyph; dispatcher attaches °).
fn fmt_deg_static(rad: f64) -> String {
    if !rad.is_finite() {
        return NAN_DASH.to_owned();
    }
    format!("{:.2}", rad.to_degrees())
}

/// Dynamic angle — 3 decimals (no degree glyph; dispatcher attaches °).
fn fmt_deg_dynamic(rad: f64) -> String {
    if !rad.is_finite() {
        return NAN_DASH.to_owned();
    }
    format!("{:.3}", rad.to_degrees())
}

/// Fixed-decimal formatter for dimensionless quantities (e, mass ratio, etc.).
fn fmt_decimal(x: f64, places: usize) -> String {
    if !x.is_finite() {
        return NAN_DASH.to_owned();
    }
    format!("{x:.*}", places)
}

/// Period — selects unit via [`time::select`] and applies the matching
/// display rule (3 decimals for days/years; scientific for large years).
fn fmt_period_auto(seconds: f64) -> (String, &'static str) {
    if !seconds.is_finite() {
        return (NAN_DASH.to_owned(), "");
    }
    let unit = time::select(seconds);
    let value = time::convert(seconds, unit);
    let s = match unit {
        time::TimeUnit::Day | time::TimeUnit::Year => format!("{value:.3}"),
        time::TimeUnit::ScientificYear => fmt_sci(value),
    };
    (s, unit.label())
}

/// Distance — selects unit via [`distance::select`] (km below 0.01 AU,
/// AU up to 1 ly, ly above) with per-unit precision (km: 2 decimals,
/// AU: 4 sig figs, ly: scientific).
fn fmt_distance_auto(meters: f64) -> (String, &'static str) {
    if !meters.is_finite() {
        return (NAN_DASH.to_owned(), "");
    }
    let unit = distance::select(meters);
    let value = distance::convert(meters, unit);
    let s = match unit {
        distance::DistanceUnit::Km => format!("{value:.2}"),
        distance::DistanceUnit::Au => fmt_au_fixed(value),
        distance::DistanceUnit::Ly => fmt_sci(value),
    };
    (s, unit.label())
}

/// Arcsecond — ≥ 10″ uses 1 decimal, < 10″ uses 2 decimals (no glyph).
fn fmt_arcsec(arcsec: f64) -> String {
    if !arcsec.is_finite() {
        return NAN_DASH.to_owned();
    }
    if arcsec.abs() >= 10.0 { format!("{arcsec:.1}") } else { format!("{arcsec:.2}") }
}

/// AU-specific fixed-precision: 4 significant figures, never scientific.
fn fmt_au_fixed(value: f64) -> String {
    if value == 0.0 {
        return "0.000".to_owned();
    }
    let abs = value.abs();
    let exp = abs.log10().floor() as i32;
    let places = (3 - exp).max(0) as usize;
    format!("{value:.*}", places)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dispatcher integration ──

    #[test]
    fn dispatcher_distance_vector_emits_si_scientific() {
        // Mercury position x ~ 5.79e10 m
        let (s, u) = format_value(5.79e10, QuantityType::DistanceVector);
        assert_eq!(s, "+5.790e+10");
        assert_eq!(u, "m");
    }

    #[test]
    fn dispatcher_distance_scalar_picks_unit_per_magnitude() {
        // Mercury semi-major axis 0.387 AU
        let (s, u) = format_value(0.387 * distance::M_PER_AU, QuantityType::DistanceScalar);
        assert_eq!(u, "AU");
        assert_eq!(s, "0.3870");

        // Lunar distance 384400 km — under 0.01 AU threshold, stays in km.
        let (s, u) = format_value(3.844e8, QuantityType::DistanceScalar);
        assert_eq!(u, "km");
        assert_eq!(s, "384400.00");

        // Proxima Centauri ~ 4.24 ly
        let (_s, u) = format_value(4.24 * distance::M_PER_LY, QuantityType::DistanceScalar);
        assert_eq!(u, "ly");
    }

    #[test]
    fn dispatcher_velocity_vector_emits_si_scientific() {
        let (s, u) = format_value(-12_300.0, QuantityType::VelocityVector);
        assert_eq!(s, "-1.230e+04");
        assert_eq!(u, "m/s");
    }

    #[test]
    fn dispatcher_mass_emits_si_scientific() {
        let (s, u) = format_value(3.302e23, QuantityType::Mass);
        assert_eq!(s, "+3.302e+23");
        assert_eq!(u, "kg");
    }

    #[test]
    fn dispatcher_time_picks_unit_per_magnitude() {
        // Mercury period ~ 87.97 d
        let (s, u) = format_value(87.97 * time::S_PER_DAY, QuantityType::Time);
        assert_eq!(u, "d");
        assert_eq!(s, "87.970");

        // Pluto period ~ 248 yr
        let (_s, u) = format_value(248.0 * time::S_PER_YEAR, QuantityType::Time);
        assert_eq!(u, "yr");

        // Sedna period ~ 11400 yr — scientific.
        let (s, u) = format_value(11_400.0 * time::S_PER_YEAR, QuantityType::Time);
        assert_eq!(u, "yr");
        assert!(s.contains('e'));
    }

    #[test]
    fn dispatcher_static_angle_two_decimals() {
        let (s, u) = format_value(7.005_f64.to_radians(), QuantityType::AngleStatic);
        assert_eq!(s, "7.00");
        assert_eq!(u, "°");
    }

    #[test]
    fn dispatcher_dynamic_angle_three_decimals() {
        let (s, u) = format_value(174.790_f64.to_radians(), QuantityType::AngleDynamic);
        assert_eq!(s, "174.790");
        assert_eq!(u, "°");
    }

    #[test]
    fn dispatcher_energy_si_scientific() {
        let (s, u) = format_value(2.736e32, QuantityType::Energy);
        assert_eq!(s, "+2.736e+32");
        assert_eq!(u, "J");
    }

    #[test]
    fn dispatcher_eccentricity_five_decimals_no_unit() {
        let (s, u) = format_value(0.20563, QuantityType::Eccentricity);
        assert_eq!(s, "0.20563");
        assert_eq!(u, "");
    }

    #[test]
    fn dispatcher_arcsec_picks_precision_per_magnitude() {
        let (s, u) = format_value(13.21, QuantityType::Arcsecond);
        assert_eq!(s, "13.2");
        assert_eq!(u, "″");

        let (s, _u) = format_value(1.30, QuantityType::Arcsecond);
        assert_eq!(s, "1.30");
    }

    #[test]
    fn dispatcher_renders_dash_for_nan_across_quantities() {
        for q in [
            QuantityType::DistanceVector,
            QuantityType::DistanceScalar,
            QuantityType::VelocityVector,
            QuantityType::Mass,
            QuantityType::Time,
            QuantityType::AngleStatic,
            QuantityType::AngleDynamic,
            QuantityType::Energy,
            QuantityType::Eccentricity,
            QuantityType::Arcsecond,
        ] {
            let (s, _) = format_value(f64::NAN, q);
            assert_eq!(s, NAN_DASH, "NaN must render as dash for {q:?}");
        }
    }

    // ── fmt_sci internals (width + stability) ──

    #[test]
    fn fmt_sci_width_is_constant_10_chars() {
        for x in [1e-30, -1e-5, 1.23, -100.0, 1.234e15, -9.999e99] {
            let s = fmt_sci(x);
            assert_eq!(s.chars().count(), 10, "fmt_sci({x}) = {s}");
        }
    }

    #[test]
    fn fmt_sci_zero_renders_as_explicit_positive() {
        assert_eq!(fmt_sci(0.0), "+0.000e+00");
    }

    #[test]
    fn fmt_sci_handles_rounding_carry_into_exponent() {
        assert_eq!(fmt_sci(9.9999), "+1.000e+01");
    }

    // ── Visual stability (the Bloomberg-flash contract) ──

    #[test]
    fn fmt_sci_sub_display_change_yields_identical_string() {
        // Two values within 1e-5 relative — must format identically so
        // BloombergFlash does not trigger on sub-precision noise.
        let a = 1.23435e10;
        let b = 1.23440e10;
        assert_eq!(fmt_sci(a), fmt_sci(b));
    }

    #[test]
    fn fmt_sci_above_display_change_yields_different_string() {
        let a = 1.23435e10;
        let b = 1.23500e10;
        assert_ne!(fmt_sci(a), fmt_sci(b));
    }

    // ── AU precision preserves orbital identity ──

    #[test]
    fn au_format_preserves_four_sig_figs() {
        // Mercury 0.3871 AU should not collapse to "0.39".
        assert_eq!(fmt_au_fixed(0.3871), "0.3871");
        // Sedna apohelion 937 AU.
        assert_eq!(fmt_au_fixed(937.0), "937.0");
        // 39.482 AU (Pluto) — 4 sig figs.
        assert_eq!(fmt_au_fixed(39.482), "39.48");
    }
}
