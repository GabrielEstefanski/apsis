//! Yoshida-4 (Forest–Ruth) composition coefficients.
//!
//! ```text
//! cbrt2  = 2^(1/3)
//! w₁     = 1 / (2 − cbrt2)  ≈  1.3512071919596578
//! w₀     = 1 − 2 w₁          ≈ −1.7024143839193156
//!
//! Drift coefficients:  c = [w₁/2,  (w₁+w₀)/2,  (w₀+w₁)/2,  w₁/2]
//! Kick  coefficients:  d = [w₁,     w₀,          w₁        ]
//! ```
//!
//! The middle kick coefficient `w₀` is negative — the second sub-step runs
//! "backwards" in time. This is correct and essential to the 4th-order
//! cancellation.
//!
//! # References
//!
//! - Forest & Ruth (1990). *Nucl. Instrum. Methods Phys. Res.* A 290, 395–400.
//! - Yoshida (1990). *Phys. Lett. A* 150, 262–268.

/// 2^(1/3)
const CBRT2: f64 = 1.2599210498948732_f64;

/// w₁ = 1 / (2 − 2^(1/3))
pub const Y4_W1: f64 = 1.0_f64 / (2.0_f64 - CBRT2);

/// w₀ = 1 − 2 w₁  (negative — middle sub-step goes backwards)
pub const Y4_W0: f64 = 1.0_f64 - 2.0_f64 * Y4_W1;

/// Drift (position) coefficients: c[i] applied before force eval i, plus final drift.
pub const Y4_C: [f64; 4] = [
    Y4_W1 * 0.5,
    (Y4_W1 + Y4_W0) * 0.5,
    (Y4_W0 + Y4_W1) * 0.5, // == Y4_C[1] by symmetry
    Y4_W1 * 0.5,           // == Y4_C[0] by symmetry
];

/// Kick (velocity) coefficients: d[i] applied after force eval i.
pub const Y4_D: [f64; 3] = [Y4_W1, Y4_W0, Y4_W1];
