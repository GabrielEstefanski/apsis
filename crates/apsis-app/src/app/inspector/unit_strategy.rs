//! Unit selection thresholds for the Inspector.
//!
//! Centralises the decision of which unit to display a quantity in. Any
//! Inspector field rendering a distance or duration consults the strategy
//! here rather than implementing its own threshold logic, so the
//! km/AU/ly and d/yr boundaries stay consistent across rows.
//!
//! The thresholds match the convention used in REBOUND, SPICE, and
//! astrodynamics texts (Vallado §3.5; Murray–Dermott §2.1): SI base
//! everywhere except where readability collapses on solar-system or
//! interstellar scales.

/// Distance unit selection.
pub mod distance {
    /// Conversion factors to SI metres.
    pub const M_PER_KM: f64 = 1.0e3;
    pub const M_PER_AU: f64 = 1.495_978_707e11;
    pub const M_PER_LY: f64 = 9.460_730_472_580_8e15;

    /// Boundary thresholds (in metres). A value at exactly the boundary
    /// rounds up to the larger unit.
    ///
    /// `AU_FROM_KM = 0.01 AU` (≈ 1.5×10⁶ km) separates the planetary
    /// regime (km — orbital around a planet, lunar distances, GEO/LEO)
    /// from the heliocentric regime (AU — solar-system bodies including
    /// Mercury at 0.387 AU). This keeps the orbital-element identity of
    /// inner planets visible (`0.387 AU` over `5.79×10⁷ km`) without
    /// putting near-Earth artificial-satellite scales into AU.
    pub const KM_FROM_M: f64 = 1.0; // floor is km — the m unit is not used
    pub const AU_FROM_KM: f64 = 0.01 * M_PER_AU;
    pub const LY_FROM_AU: f64 = M_PER_LY; // ≥ 1 ly expressed in ly

    /// Selected unit for a metres value.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum DistanceUnit {
        Km,
        Au,
        Ly,
    }

    impl DistanceUnit {
        pub fn label(self) -> &'static str {
            match self {
                Self::Km => "km",
                Self::Au => "AU",
                Self::Ly => "ly",
            }
        }
    }

    /// Choose the unit for displaying a distance in metres. Negative or
    /// non-finite inputs select [`DistanceUnit::Km`] as a neutral default —
    /// formatters render the special value (`NaN`/`±∞`) directly.
    pub fn select(meters: f64) -> DistanceUnit {
        let m = meters.abs();
        if !m.is_finite() {
            return DistanceUnit::Km;
        }
        if m >= LY_FROM_AU {
            DistanceUnit::Ly
        } else if m >= AU_FROM_KM {
            DistanceUnit::Au
        } else {
            DistanceUnit::Km
        }
    }

    /// Convert metres into the chosen unit's natural value.
    pub fn convert(meters: f64, unit: DistanceUnit) -> f64 {
        match unit {
            DistanceUnit::Km => meters / M_PER_KM,
            DistanceUnit::Au => meters / M_PER_AU,
            DistanceUnit::Ly => meters / M_PER_LY,
        }
    }
}

/// Time unit selection.
pub mod time {
    /// Conversion factors to SI seconds.
    pub const S_PER_DAY: f64 = 86_400.0;
    pub const S_PER_YEAR: f64 = 365.25 * S_PER_DAY;

    /// Boundary thresholds (in days). `< 1000 d` displays in days,
    /// `< 1e6 d` in years with fixed decimals, `≥ 1e6 d` in years with
    /// scientific notation.
    pub const DAYS_TO_YEARS_THRESHOLD: f64 = 1_000.0;
    pub const YEARS_FIXED_TO_SCI_THRESHOLD: f64 = 1_000_000.0;

    /// Selected unit for a seconds value.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum TimeUnit {
        /// Days, fixed decimal display.
        Day,
        /// Years, fixed decimal display.
        Year,
        /// Years, scientific notation display.
        ScientificYear,
    }

    impl TimeUnit {
        pub fn label(self) -> &'static str {
            match self {
                Self::Day => "d",
                Self::Year | Self::ScientificYear => "yr",
            }
        }
    }

    /// Choose the unit for displaying a duration in seconds. Negative or
    /// non-finite inputs select [`TimeUnit::Day`] as a neutral default.
    pub fn select(seconds: f64) -> TimeUnit {
        let s = seconds.abs();
        if !s.is_finite() {
            return TimeUnit::Day;
        }
        let days = s / S_PER_DAY;
        if days < DAYS_TO_YEARS_THRESHOLD {
            TimeUnit::Day
        } else if days < YEARS_FIXED_TO_SCI_THRESHOLD {
            TimeUnit::Year
        } else {
            TimeUnit::ScientificYear
        }
    }

    /// Convert seconds into the chosen unit's natural value.
    pub fn convert(seconds: f64, unit: TimeUnit) -> f64 {
        match unit {
            TimeUnit::Day => seconds / S_PER_DAY,
            TimeUnit::Year | TimeUnit::ScientificYear => seconds / S_PER_YEAR,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Distance ──

    #[test]
    fn distance_below_planetary_threshold_selects_km() {
        // Geosynchronous orbit ~ 4.2e7 m = 4.2e4 km — under 0.01 AU.
        assert_eq!(distance::select(4.2e7), distance::DistanceUnit::Km);
        // Lunar distance ~ 384400 km ≈ 0.0026 AU — still km.
        assert_eq!(distance::select(3.844e8), distance::DistanceUnit::Km);
    }

    #[test]
    fn distance_above_planetary_threshold_selects_au() {
        // Mercury distance ~ 0.387 AU — over the 0.01 AU threshold.
        assert_eq!(distance::select(0.387 * distance::M_PER_AU), distance::DistanceUnit::Au);
        // 1 AU exact — AU.
        assert_eq!(distance::select(distance::M_PER_AU), distance::DistanceUnit::Au);
        // Sedna's apohelion ~ 937 AU.
        assert_eq!(distance::select(937.0 * distance::M_PER_AU), distance::DistanceUnit::Au);
    }

    #[test]
    fn distance_above_one_ly_selects_ly() {
        // 1 ly exact -> Ly
        assert_eq!(distance::select(distance::M_PER_LY), distance::DistanceUnit::Ly);
        // Proxima Centauri ~ 4.24 ly
        assert_eq!(distance::select(4.24 * distance::M_PER_LY), distance::DistanceUnit::Ly);
    }

    #[test]
    fn distance_handles_negative_and_nonfinite() {
        // `select` reads magnitude; negative AU still selects AU.
        assert_eq!(distance::select(-distance::M_PER_AU * 5.0), distance::DistanceUnit::Au);
        // Non-finite inputs fall back to Km neutrally.
        assert_eq!(distance::select(f64::NAN), distance::DistanceUnit::Km);
        assert_eq!(distance::select(f64::INFINITY), distance::DistanceUnit::Km);
    }

    #[test]
    fn distance_conversion_inverts_constants() {
        let m = 0.387 * distance::M_PER_AU;
        assert!(
            (distance::convert(m, distance::DistanceUnit::Au) - 0.387).abs() < 1e-12,
            "AU conversion off",
        );
    }

    // ── Time ──

    #[test]
    fn time_short_period_selects_day() {
        // Mercury period ~ 87.97 d
        assert_eq!(time::select(87.97 * time::S_PER_DAY), time::TimeUnit::Day);
        // ISS period ~ 0.064 d
        assert_eq!(time::select(0.064 * time::S_PER_DAY), time::TimeUnit::Day);
    }

    #[test]
    fn time_medium_period_selects_year() {
        // Pluto period ~ 248 yr ≈ 90581 d
        assert_eq!(time::select(248.0 * time::S_PER_YEAR), time::TimeUnit::Year);
        // Sedna period ~ 11400 yr ≈ 4.16e6 d → over the year-fixed threshold
        // (≥ 1e6 d) so it should flip into scientific.
        assert_eq!(time::select(11_400.0 * time::S_PER_YEAR), time::TimeUnit::ScientificYear,);
    }

    #[test]
    fn time_threshold_boundaries() {
        // Just below 1000 d → Day.
        assert_eq!(time::select(999.999 * time::S_PER_DAY), time::TimeUnit::Day,);
        // 1000 d exactly → Year.
        assert_eq!(time::select(1_000.0 * time::S_PER_DAY), time::TimeUnit::Year);
        // 1e6 d exactly → ScientificYear.
        assert_eq!(time::select(1_000_000.0 * time::S_PER_DAY), time::TimeUnit::ScientificYear,);
    }

    #[test]
    fn time_conversion_inverts_constants() {
        let s = 87.969 * time::S_PER_DAY;
        let days = time::convert(s, time::TimeUnit::Day);
        assert!((days - 87.969).abs() < 1e-9, "day conversion off");
    }
}
