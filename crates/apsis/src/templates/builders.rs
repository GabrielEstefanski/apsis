//! Helpers for constructing initial states from orbital parameters.
//!
//! Templates pull from these instead of computing trig and rotations inline,
//! so a Keplerian setup reads as one line per body and the heavy math is
//! audited in one place.

/// Convert a published density from SI (kg/m³) to the canonical
/// solar-AU simulation units (M_☉/AU³).
///
/// Templates that quote real bodies use this to write per-body density
/// in source as the literal NASA fact-sheet number multiplied by this
/// constant — e.g. `5514.0 * KG_M3_TO_SOLAR_AU3` for Earth. The
/// resulting value goes through [`Body::with_density`](crate::domain::body::Body::with_density),
/// which recomputes `physical_radius` so the rendered geometry tracks
/// the published radius to within a few percent.
///
/// # Derivation
///
/// `1 M_☉ / AU³ = (1.495 978 707 × 10¹¹ m)⁻³ · 1.988 92 × 10³⁰ kg`
/// (IAU 2012 nominal AU; IAU 2015 nominal solar-mass parameter
/// divided by CODATA 2018 G):
///
/// ```text
/// AU³        = 3.347 928 9 × 10³³ m³
/// M_☉ / AU³  = 1.988 92e30 / 3.347 928 9e33 kg/m³
///            = 5.940 84 × 10⁻⁴ kg/m³
/// kg/m³ → M_☉/AU³ multiplier
///            = 1 / 5.940 84e-4
///            = 1683.26
/// ```
///
/// This value is unit-system specific (canonical solar-AU). Templates
/// that run in different unit systems must derive their own conversion
/// from `length_to_m` / `mass_to_kg` on the [`UnitSystem`](crate::templates::UnitSystem).
pub const KG_M3_TO_SOLAR_AU3: f64 = 1683.26;

/// Circular orbit in the XY plane around a body of mass `center_mass`
/// fixed at the origin. Returns inertial-frame `(position, velocity)`
/// in 3D, with `z = vz = 0`.
///
/// `radius` and `center_mass` are in simulation units; the embedded
/// implicit `G = 1` makes velocity drop out of `v = sqrt(GM/r)`. Phase
/// `phase` is the true anomaly at t=0 in radians, measured CCW from
/// the +X axis.
pub fn circular_orbit(center_mass: f64, radius: f64, phase: f64) -> ([f64; 3], [f64; 3]) {
    let x = radius * phase.cos();
    let y = radius * phase.sin();

    let v = (center_mass / radius).sqrt();

    let vx = -v * phase.sin();
    let vy = v * phase.cos();

    ([x, y, 0.0], [vx, vy, 0.0])
}
