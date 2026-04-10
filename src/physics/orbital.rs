//! Osculating orbital element computation for 2D N-body systems.
//!
//! ## Osculating elements
//!
//! In an N-body system, orbital elements are only strictly defined for isolated
//! two-body problems. The **osculating** elements are the Keplerian elements of
//! the equivalent two-body orbit that shares the body's current position and
//! velocity. They are a snapshot, not a conserved quantity, but they are the
//! standard way to characterise instantaneous orbits in N-body codes.
//!
//! ## Primary selection
//!
//! For each body `i`, the **primary** is the body `j ≠ i` that produces the
//! largest gravitational acceleration at body `i`'s location:
//!
//! ```text
//! dominant j = argmax_j  G·m_j / r_ij²
//! ```
//!
//! This correctly handles hierarchical systems: the Moon's primary is Earth
//! even though the Sun is more massive, because Earth is closer.
//!
//! ## Computed quantities
//!
//! | Symbol | Name | Valid when |
//! |--------|------|-----------|
//! | `a`    | semi-major axis | bound orbit (e < 1) |
//! | `e`    | eccentricity | always |
//! | `T`    | period | bound orbit |
//! | `h`    | specific angular momentum (z) | always |
//! | `ε`    | specific orbital energy | always |
//! | `ω`    | argument of periapsis | e > 1e-6 |
//!
//! ## References
//! - Murray & Dermott (1999). *Solar System Dynamics*. Cambridge.
//! - Bate, Mueller & White (1971). *Fundamentals of Astrodynamics*. Dover.

use std::f64::consts::TAU;

use crate::domain::body::Body;

// ── Orbit classification ──────────────────────────────────────────────────────

/// Keplerian orbit type derived from specific orbital energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitType {
    /// Bound elliptical orbit (ε < 0, e < 1).
    Elliptical,
    /// Near-parabolic transition (|ε| < threshold, e ≈ 1).
    Parabolic,
    /// Unbound hyperbolic flyby (ε > 0, e > 1).
    Hyperbolic,
}

impl OrbitType {
    pub fn label(self) -> &'static str {
        match self {
            Self::Elliptical => "elliptical",
            Self::Parabolic => "parabolic",
            Self::Hyperbolic => "hyperbolic",
        }
    }

    pub fn is_bound(self) -> bool {
        matches!(self, Self::Elliptical | Self::Parabolic)
    }
}

// ── OrbitalElements ───────────────────────────────────────────────────────────

/// Osculating Keplerian orbital elements for one body at one instant.
///
/// All quantities are in simulation units. `a`, `T`, and `ω` are only
/// physically meaningful for bound (elliptical) orbits.
#[derive(Debug, Clone, Copy)]
pub struct OrbitalElements {
    /// Index of the dominant primary body this orbit is computed relative to.
    pub primary_idx: usize,

    /// Semi-major axis.  `f64::INFINITY` for parabolic/hyperbolic orbits.
    pub a: f64,

    /// Eccentricity (dimensionless). `0` = circular, `1` = parabolic, `>1` = hyperbolic.
    pub e: f64,

    /// Orbital period.  `f64::INFINITY` for unbound orbits.
    pub period: f64,

    /// Specific angular momentum (z-component, scalar in 2D).
    /// Positive = counter-clockwise.
    pub h: f64,

    /// Specific orbital energy (`ε = v²/2 − GM/r`).
    /// Negative = bound, positive = unbound.
    pub energy: f64,

    /// Argument of periapsis ω ∈ [−π, π] (radians).
    /// Undefined when `e < 1e-6` (circular); returns 0.
    pub omega: f64,

    /// Orbit classification derived from `energy`.
    pub orbit_type: OrbitType,
}

impl OrbitalElements {
    /// Returns `true` for elliptical and parabolic orbits.
    pub fn is_bound(self) -> bool {
        self.orbit_type.is_bound()
    }

    /// CSV header row matching [`Self::to_csv_row`].
    pub fn csv_header() -> &'static str {
        "t,body_idx,primary_idx,a,e,period,h,energy,omega_deg,orbit_type"
    }

    /// Serialise to a CSV data row.  `t` and `body_idx` are injected by the caller.
    pub fn to_csv_row(self, t: f64, body_idx: usize) -> String {
        format!(
            "{t:.6e},{body_idx},{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.4},{:?}",
            self.primary_idx,
            self.a,
            self.e,
            self.period,
            self.h,
            self.energy,
            self.omega.to_degrees(),
            self.orbit_type.label(),
        )
    }
}

// ── Primary selection ─────────────────────────────────────────────────────────

/// Returns the index of the gravitationally dominant body for body `idx`.
///
/// Dominant = maximises `G·m_j / r_ij²` (i.e. largest acceleration contribution).
/// Returns `None` when the system has fewer than 2 bodies.
pub fn dominant_primary(bodies: &[Body], idx: usize) -> Option<usize> {
    if bodies.len() < 2 {
        return None;
    }

    let bi = &bodies[idx];

    bodies
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != idx)
        .max_by(|(_, bj), (_, bk)| {
            let rj2 = (bj.x - bi.x).powi(2) + (bj.y - bi.y).powi(2);
            let rk2 = (bk.x - bi.x).powi(2) + (bk.y - bi.y).powi(2);
            // Compare m_j/r_j² vs m_k/r_k²  (G cancels)
            let score_j = if rj2 > 0.0 { bj.mass / rj2 } else { f64::INFINITY };
            let score_k = if rk2 > 0.0 { bk.mass / rk2 } else { f64::INFINITY };
            score_j.partial_cmp(&score_k).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(j, _)| j)
}

// ── Element computation ───────────────────────────────────────────────────────

/// Compute osculating orbital elements for body `idx` relative to `primary_idx`.
///
/// Returns `None` if the bodies are co-located (r < 1e-15) or have zero combined mass.
pub fn compute_elements(
    bodies: &[Body],
    idx: usize,
    primary_idx: usize,
    g_factor: f64,
) -> Option<OrbitalElements> {
    let b = &bodies[idx];
    let p = &bodies[primary_idx];

    // Relative state vector
    let rx = b.x - p.x;
    let ry = b.y - p.y;
    let vrx = b.vx - p.vx;
    let vry = b.vy - p.vy;

    let r = (rx * rx + ry * ry).sqrt();
    let v2 = vrx * vrx + vry * vry;
    let gm = g_factor * (b.mass + p.mass);

    if r < 1e-15 || gm < 1e-30 {
        return None;
    }

    // Specific orbital energy  ε = ½v² − GM/r
    let energy = 0.5 * v2 - gm / r;

    // Specific angular momentum (z-component)   h = r × v
    let h = rx * vry - ry * vrx;

    // Eccentricity via Laplace–Runge–Lenz vector
    //   e_vec = (v × h) / GM  −  r̂
    // In 2D:  (v × h̄) = (vry·h,  −vrx·h)
    let ex = vry * h / gm - rx / r;
    let ey = -vrx * h / gm - ry / r;
    let e = (ex * ex + ey * ey).sqrt();

    // Argument of periapsis ω = angle of eccentricity vector
    let omega = if e > 1e-6 { ey.atan2(ex) } else { 0.0 };

    // Semi-major axis and period
    const ENERGY_THRESH: f64 = 1e-12;

    let (a, period, orbit_type) = if energy < -ENERGY_THRESH {
        // Bound elliptical:  a = −GM / (2ε)
        let a = -gm / (2.0 * energy);
        // Kepler III:  T = 2π √(a³/GM)
        let period = TAU * (a * a * a / gm).sqrt();
        (a, period, OrbitType::Elliptical)
    } else if energy < ENERGY_THRESH {
        // Near-parabolic transition
        (f64::INFINITY, f64::INFINITY, OrbitType::Parabolic)
    } else {
        // Hyperbolic:  a is negative by convention (distance to focus)
        let a = -gm / (2.0 * energy);
        (a, f64::INFINITY, OrbitType::Hyperbolic)
    };

    Some(OrbitalElements {
        primary_idx,
        a,
        e,
        period,
        h,
        energy,
        omega,
        orbit_type,
    })
}

// ── Batch computation ─────────────────────────────────────────────────────────

/// Compute osculating elements for every body.
///
/// Returns one `Option<OrbitalElements>` per body. The dominant primary is
/// determined automatically via [`dominant_primary`].
/// Returns `None` for a body when its primary cannot be found or the geometry
/// is degenerate.
pub fn compute_all(bodies: &[Body], g_factor: f64) -> Vec<Option<OrbitalElements>> {
    bodies
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let p = dominant_primary(bodies, i)?;
            compute_elements(bodies, i, p, g_factor)
        })
        .collect()
}
