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

use std::f64::consts::{PI, TAU};

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
///
/// # 3D readiness
///
/// The struct carries `inclination` and `lon_ascending_node` even though
/// the simulation is currently 2D. In 2D both are exactly zero and the
/// perifocal → world rotation collapses to a single rotation by `ω`. When
/// the engine moves to 3D, [`compute_elements`] will populate them from
/// the angular-momentum vector and the ascending-node geometry; every
/// consumer of [`sample_orbit`] keeps working unchanged.
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

    /// Inclination i ∈ [0, π] (radians). Always `0` in 2D.
    pub inclination: f64,

    /// Longitude of the ascending node Ω ∈ [−π, π] (radians). Always `0` in 2D.
    pub lon_ascending_node: f64,

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

    /// Samples the predicted Keplerian orbit in **world coordinates**.
    ///
    /// Returns `steps + 1` points.
    ///
    /// # Parametrisation
    ///
    /// **Elliptical** — uniform eccentric anomaly `E ∈ [0, 2π]`:
    ///
    /// ```text
    /// x_pf = a (cos E − e)
    /// y_pf = a √(1 − e²) sin E
    /// ```
    ///
    /// The polyline **closes exactly** (first = last) and avoids
    /// periastro-clustering for any eccentricity up to ~0.95.
    ///
    /// **Hyperbolic** — uniform hyperbolic anomaly `H ∈ [−H_max, +H_max]`:
    ///
    /// ```text
    /// x_pf = |a| (e − cosh H)
    /// y_pf = |a| √(e² − 1) sinh H
    /// ```
    ///
    /// At `H = 0` the sample lies on the periapsis `(|a|(e−1), 0)`. The
    /// polyline is **open** — it traces one branch of the hyperbola
    /// clipped at `|H| ≤ H_max` (see [`Self::H_MAX_HYPER`]), since the
    /// true branch extends to infinity along the asymptotes.
    ///
    /// **Parabolic** — no finite semi-major axis; returns empty.
    ///
    /// # Frame transform
    ///
    /// Perifocal → world via `R₃(Ω) · R₁(i) · R₃(ω)`, then translated by
    /// the primary's position (the conic is centred on the **focus**,
    /// not the geometric centre). In 2D (i = 0, Ω = 0) this collapses to
    /// a rotation by `ω` alone.
    ///
    /// # Return value
    ///
    /// * `Vec<[f64; 3]>` in world coordinates, z ≡ 0 in 2D.
    /// * Empty vector for parabolic orbits, degenerate geometry, or
    ///   `steps < 2`.
    pub fn sample_orbit(&self, primary_pos: [f64; 3], steps: usize) -> Vec<[f64; 3]> {
        if steps < 2 {
            return Vec::new();
        }

        // Composite rotation R₃(Ω) · R₁(i) · R₃(ω) applied to a perifocal
        // point (x_pf, y_pf, 0). Only the first two columns matter because
        // z_pf = 0; we pre-compute those as six scalars.
        let (sw, cw) = self.omega.sin_cos();
        let (si, ci) = self.inclination.sin_cos();
        let (so, co) = self.lon_ascending_node.sin_cos();

        let r11 = co * cw - so * sw * ci;
        let r12 = -co * sw - so * cw * ci;
        let r21 = so * cw + co * sw * ci;
        let r22 = -so * sw + co * cw * ci;
        let r31 = sw * si;
        let r32 = cw * si;

        let project = |x_pf: f64, y_pf: f64| -> [f64; 3] {
            [
                r11 * x_pf + r12 * y_pf + primary_pos[0],
                r21 * x_pf + r22 * y_pf + primary_pos[1],
                r31 * x_pf + r32 * y_pf + primary_pos[2],
            ]
        };

        match self.orbit_type {
            OrbitType::Elliptical => {
                if !self.a.is_finite() || self.a <= 0.0 {
                    return Vec::new();
                }
                let a = self.a;
                let e = self.e.clamp(0.0, 0.999);
                let b = a * (1.0 - e * e).sqrt();

                let mut out = Vec::with_capacity(steps + 1);
                for k in 0..=steps {
                    let ek = TAU * (k as f64) / (steps as f64);
                    let (s_ek, c_ek) = ek.sin_cos();
                    let x_pf = a * (c_ek - e);
                    let y_pf = b * s_ek;
                    out.push(project(x_pf, y_pf));
                }
                out
            },
            OrbitType::Hyperbolic => {
                if !self.a.is_finite() || self.e <= 1.0 {
                    return Vec::new();
                }
                // Hyperbolic: a < 0 by the codebase convention
                // (compute_elements uses a = -GM/(2ε) with ε > 0).
                let a_h = self.a.abs();
                let e = self.e;
                let b_h = a_h * (e * e - 1.0).sqrt();
                let h_max = Self::H_MAX_HYPER;

                let mut out = Vec::with_capacity(steps + 1);
                for k in 0..=steps {
                    // H uniformly spans [-h_max, +h_max].
                    let t = (k as f64) / (steps as f64); // 0..=1
                    let h = -h_max + 2.0 * h_max * t;
                    let ch = h.cosh();
                    let sh = h.sinh();
                    let x_pf = a_h * (e - ch); // H=0 → a_h(e−1) = r_peri
                    let y_pf = b_h * sh;
                    out.push(project(x_pf, y_pf));
                }
                out
            },
            OrbitType::Parabolic => Vec::new(),
        }
    }

    /// Hyperbolic anomaly clip used by [`Self::sample_orbit`]. At `H = 3`
    /// the sampled arm reaches ≈ 10×e periapsis distances (`cosh 3 ≈ 10`),
    /// which covers the interesting part of a flyby for any `e > 1` while
    /// staying well short of the asymptotes (ν_∞ = arccos(−1/e)).
    pub const H_MAX_HYPER: f64 = 3.0;

    /// Periapsis position in **world coordinates**.
    ///
    /// Returns `None` for degenerate / unresolvable geometry:
    /// * parabolic (`a = ∞`) — the parabola has a periapsis but the
    ///   element set here does not encode it separately;
    /// * non-finite semi-major axis;
    /// * eccentricity ≥ 1 with `a ≥ 0` (malformed element).
    pub fn periapsis_world(&self, primary_pos: [f64; 3]) -> Option<[f64; 3]> {
        if !self.a.is_finite() {
            return None;
        }
        // r_peri handles both sign conventions: ellipse (a > 0) → a(1−e),
        // hyperbola (a < 0) → |a|(e−1). Both reduce to |a(1−e)|.
        let r_peri = match self.orbit_type {
            OrbitType::Elliptical => self.a * (1.0 - self.e),
            OrbitType::Hyperbolic => self.a.abs() * (self.e - 1.0),
            OrbitType::Parabolic => return None,
        };
        if !r_peri.is_finite() || r_peri <= 0.0 {
            return None;
        }
        // Perifocal periapsis: (r_peri, 0, 0), then rotate by (Ω, i, ω).
        Some(self.rotate_perifocal([r_peri, 0.0, 0.0], primary_pos))
    }

    /// Apoapsis position in **world coordinates**.
    ///
    /// Returns `None` for unbound orbits (hyperbolic has no apoapsis; the
    /// body escapes to infinity) and for parabolic / degenerate cases.
    pub fn apoapsis_world(&self, primary_pos: [f64; 3]) -> Option<[f64; 3]> {
        if !matches!(self.orbit_type, OrbitType::Elliptical) {
            return None;
        }
        if !self.a.is_finite() || self.a <= 0.0 {
            return None;
        }
        // Apoapsis in perifocal: x = −a(1+e), y = 0 (opposite side of focus
        // from periapsis).
        let x_pf = -self.a * (1.0 + self.e);
        Some(self.rotate_perifocal([x_pf, 0.0, 0.0], primary_pos))
    }

    /// Apply the perifocal → world rotation `R₃(Ω)·R₁(i)·R₃(ω)` and
    /// translate by `primary_pos`. Shared by apsis accessors so they stay
    /// in sync with [`Self::sample_orbit`].
    fn rotate_perifocal(&self, pf: [f64; 3], primary_pos: [f64; 3]) -> [f64; 3] {
        let (sw, cw) = self.omega.sin_cos();
        let (si, ci) = self.inclination.sin_cos();
        let (so, co) = self.lon_ascending_node.sin_cos();

        let r11 = co * cw - so * sw * ci;
        let r12 = -co * sw - so * cw * ci;
        let r21 = so * cw + co * sw * ci;
        let r22 = -so * sw + co * cw * ci;
        let r31 = sw * si;
        let r32 = cw * si;

        // z_pf ≡ 0 for all perifocal apsis points.
        [
            r11 * pf[0] + r12 * pf[1] + primary_pos[0],
            r21 * pf[0] + r22 * pf[1] + primary_pos[1],
            r31 * pf[0] + r32 * pf[1] + primary_pos[2],
        ]
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
/// Linear-in-state primitives from which all Keplerian elements derive.
///
/// Exposing these lets external code (e.g. the render-side osculating-element
/// smoother in `apsis-app`) apply transformations *before* the non-linear
/// reconstruction into (a, e, ω). EMA on (a, e, sin ω, cos ω) is dimensionally
/// noisier than EMA on these four scalars: ε is linear in (v², 1/r), the
/// Laplace–Runge–Lenz components (ex, ey) are linear in (r, v), and h is
/// bilinear in (r, v). They have no angular wraparound and no singularity at
/// the parabolic limit.
#[derive(Debug, Clone, Copy)]
pub struct OrbitInvariants {
    /// Specific orbital energy ε = ½v² − GM/r.
    pub energy: f64,
    /// Specific angular momentum z-component h = r × v.
    pub h: f64,
    /// Laplace–Runge–Lenz vector x-component (eccentricity vector).
    pub ex: f64,
    /// Laplace–Runge–Lenz vector y-component (eccentricity vector).
    pub ey: f64,
}

/// Computes the linear primitives (ε, h, ex, ey) for body `idx` relative to
/// `primary_idx` under gravitational parameter `g_factor`.
///
/// Returns `None` for degenerate inputs (coincident bodies, zero GM).
pub fn compute_invariants(
    bodies: &[Body],
    idx: usize,
    primary_idx: usize,
    g_factor: f64,
) -> Option<OrbitInvariants> {
    let b = &bodies[idx];
    let p = &bodies[primary_idx];

    let rx = b.x - p.x;
    let ry = b.y - p.y;
    let vrx = b.vx - p.vx;
    let vry = b.vy - p.vy;

    let r = (rx * rx + ry * ry).sqrt();
    let gm = g_factor * (b.mass + p.mass);

    if r < 1e-15 || gm < 1e-30 {
        return None;
    }

    let v2 = vrx * vrx + vry * vry;
    let energy = 0.5 * v2 - gm / r;
    let h = rx * vry - ry * vrx;
    let ex = vry * h / gm - rx / r;
    let ey = -vrx * h / gm - ry / r;

    Some(OrbitInvariants { energy, h, ex, ey })
}

/// Reconstructs `OrbitalElements` from the linear primitives plus the
/// gravitational parameter and the chosen primary index.
///
/// Pure function — no body list lookup. This is the inverse of
/// [`compute_invariants`]: `elements_from_invariants(compute_invariants(...).
/// unwrap(), primary_idx, gm)` is equal (up to floating-point) to
/// [`compute_elements`].
pub fn elements_from_invariants(
    inv: &OrbitInvariants,
    primary_idx: usize,
    gm: f64,
) -> OrbitalElements {
    let e = (inv.ex * inv.ex + inv.ey * inv.ey).sqrt();
    let omega = if e > 1e-6 { inv.ey.atan2(inv.ex) } else { 0.0 };

    const ENERGY_THRESH: f64 = 1e-12;

    let (a, period, orbit_type) = if inv.energy < -ENERGY_THRESH {
        let a = -gm / (2.0 * inv.energy);
        let period = TAU * (a * a * a / gm).sqrt();
        (a, period, OrbitType::Elliptical)
    } else if inv.energy < ENERGY_THRESH {
        (f64::INFINITY, f64::INFINITY, OrbitType::Parabolic)
    } else {
        let a = -gm / (2.0 * inv.energy);
        (a, f64::INFINITY, OrbitType::Hyperbolic)
    };

    OrbitalElements {
        primary_idx,
        a,
        e,
        period,
        h: inv.h,
        energy: inv.energy,
        omega,
        inclination: 0.0,
        lon_ascending_node: 0.0,
        orbit_type,
    }
}

pub fn compute_elements(
    bodies: &[Body],
    idx: usize,
    primary_idx: usize,
    g_factor: f64,
) -> Option<OrbitalElements> {
    let inv = compute_invariants(bodies, idx, primary_idx, g_factor)?;
    let gm = g_factor * (bodies[idx].mass + bodies[primary_idx].mass);
    Some(elements_from_invariants(&inv, primary_idx, gm))
}

/// Eccentricity below which [`elements_anchored_to_body`] skips the
/// geometric anchor and returns the smoothed-direction ellipse instead.
///
/// At this threshold the body sits at most ~5% of `a` off the displayed
/// orbit — sub-pixel at typical viewing zoom for solar-system-scale
/// scenes — which is invisible compared to the angular wobble the
/// anchor produces in the same regime.
pub const ANCHOR_MIN_E: f64 = 0.05;

/// Reconstructs `OrbitalElements` whose ellipse passes geometrically
/// through the body's current position.
///
/// Uses the **smoothed** invariants for shape (a, e via energy, e_vec
/// magnitude) but recomputes the argument of periapsis ω from the
/// instantaneous (r⃗, v⃗) so the rendered ellipse contains the body
/// exactly. This is the rendering contract the EMA smoother needs to
/// preserve consistency between an averaged orbit shape and the body's
/// instantaneous position.
///
/// # Algorithm
///
/// Given smoothed (a, e) and body's current state vector (r⃗, v⃗) in
/// the primary's frame:
///
/// 1. From the focus-conic equation r = p / (1 + e·cos ν) with p = a(1−e²):
///    `cos ν = (p/r − 1) / e`. Clamp to [−1, 1] to absorb cases where
///    the body's distance falls slightly outside the smoothed ellipse's
///    [a(1−e), a(1+e)] range — the body is then placed at the closer
///    apsis of the displayed orbit, which is the safe fallback.
///
/// 2. Sign of sin ν from the analytic identity
///    `r_x·ey − r_y·ex = −h·(r⃗·v⃗)/GM`, giving
///    `sign(sin ν) = sign(h_inst · (r⃗·v⃗))`. This handles both prograde
///    and retrograde orbits correctly. The factor `h_inst` is computed
///    from the *current* state (orbit's observed chirality), not the
///    smoothed h — the displayed orbit must match what the user is
///    seeing right now.
///
/// 3. ω = atan2(r_y, r_x) − ν.
///
/// # Edge cases
///
/// * Hyperbolic / parabolic / unbound (`energy ≥ 0`): falls back to
///   [`elements_from_invariants`]. The smoother already short-circuits
///   these, but defensive.
/// * Near-circular (`e < ANCHOR_MIN_E`): falls back to
///   [`elements_from_invariants`]. Anchoring divides by e to recover ν,
///   so a small smoothed residual `|δe_vec|` produces an angular wobble
///   `δω ≈ |δe_vec|/e` that becomes visible as the orbit polyline
///   "rotating" frame-to-frame. The smoothed eccentricity vector
///   `atan2(ey, ex)` already gives a stable orientation; the body sits
///   at most ~e·a off the displayed circle, which is sub-pixel for the
///   regime this branch catches.
/// * `r ≈ 0` (coincident with primary): falls back to non-anchored.
pub fn elements_anchored_to_body(
    inv: &OrbitInvariants,
    primary_idx: usize,
    gm: f64,
    body: &Body,
    primary: &Body,
) -> OrbitalElements {
    const ENERGY_THRESH: f64 = 1e-12;

    // Unbound (or near-parabolic): no closed conic to anchor — defer.
    if inv.energy >= -ENERGY_THRESH {
        return elements_from_invariants(inv, primary_idx, gm);
    }

    let e = (inv.ex * inv.ex + inv.ey * inv.ey).sqrt();

    // Near-circular: anchor's (p/r − 1)/e branch amplifies smoothed
    // residual into visible ω jitter. Use the smoothed eccentricity
    // vector's direction directly — stable, with at most ~e·a body
    // offset from the displayed orbit (sub-pixel for e < this threshold).
    if e < ANCHOR_MIN_E {
        return elements_from_invariants(inv, primary_idx, gm);
    }

    let a = -gm / (2.0 * inv.energy);
    let period = TAU * (a * a * a / gm).sqrt();

    let rx = body.x - primary.x;
    let ry = body.y - primary.y;
    let r = (rx * rx + ry * ry).sqrt();
    if r < 1e-15 {
        return elements_from_invariants(inv, primary_idx, gm);
    }
    let vrx = body.vx - primary.vx;
    let vry = body.vy - primary.vy;

    // Solve focus-conic for cos ν, clamping to absorb body distances
    // slightly outside the smoothed ellipse's apsidal range.
    let p = a * (1.0 - e * e);
    let cos_nu = ((p / r - 1.0) / e).clamp(-1.0, 1.0);

    // sign(sin ν) = sign(h_inst · (r⃗·v⃗)). Both factors are needed:
    // h_inst alone gives orbit chirality (prograde/retrograde);
    // (r⃗·v⃗) alone is wrong for retrograde. The product encodes both
    // the leg (outbound/inbound) and the chirality.
    //
    // Near peri/apo the (r⃗·v⃗) term passes through zero and its sign
    // is FP-fragile, but sin ν ≈ 0 there too — both ν = +ε and ν = −ε
    // map to nearly identical positions, so a flicker is sub-pixel.
    let h_inst = rx * vry - ry * vrx;
    let r_dot_v = rx * vrx + ry * vry;
    let s = (h_inst * r_dot_v).signum();
    // signum() returns 0 only when h·(r·v) = 0 exactly; pick +1 then.
    let sin_nu_sign = if s == 0.0 { 1.0 } else { s };

    let nu = sin_nu_sign * cos_nu.acos();
    let theta_body = ry.atan2(rx);
    let omega_raw = theta_body - nu;
    // Wrap to [−π, π] to match the convention of elements_from_invariants.
    let omega = (omega_raw + PI).rem_euclid(TAU) - PI;

    OrbitalElements {
        primary_idx,
        a,
        e,
        period,
        h: inv.h,
        energy: inv.energy,
        omega,
        inclination: 0.0,
        lon_ascending_node: 0.0,
        orbit_type: OrbitType::Elliptical,
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;
    use std::f64::consts::{PI, TAU};

    // ── Helpers ───────────────────────────────────────────────────────────────

    const G: f64 = 1.0;
    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::rocky(mass).at(x, y).with_velocity(vx, vy)
    }

    /// Órbita circular perfeita ao redor da origem.
    ///
    /// Para r e M dados, a velocidade circular é v = sqrt(GM/r).
    /// Coloca o corpo em (r, 0) movendo-se em (0, v_c) → CCW.
    fn circular_orbit(r: f64, primary_mass: f64) -> (Body, Body) {
        let v_c = (G * primary_mass / r).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let satellite = body(r, 0.0, 0.0, v_c, 1e-10); // massa negligível
        (primary, satellite)
    }

    fn elements(primary: Body, satellite: Body) -> OrbitalElements {
        let bodies = vec![primary, satellite];
        compute_elements(&bodies, 1, 0, G).expect("elementos devem ser computáveis")
    }

    // ── 0. Refactor invariants (linear primitives) ───────────────────────────
    //
    // compute_elements is now defined as the composition
    // elements_from_invariants ∘ compute_invariants. This must hold exactly
    // (modulo floating-point fma differences) for every regime.

    fn assert_elements_equiv(a: &OrbitalElements, b: &OrbitalElements, label: &str) {
        let rel = |x: f64, y: f64| {
            if x.is_finite() && y.is_finite() {
                (x - y).abs() / x.abs().max(y.abs()).max(1e-30)
            } else if x.is_infinite() && y.is_infinite() && x.signum() == y.signum() {
                0.0
            } else {
                f64::INFINITY
            }
        };
        let tol = 1e-12;
        assert!(rel(a.a, b.a) < tol, "{label}: a {} vs {}", a.a, b.a);
        assert!(rel(a.e, b.e) < tol, "{label}: e {} vs {}", a.e, b.e);
        assert!(rel(a.energy, b.energy) < tol, "{label}: ε {} vs {}", a.energy, b.energy);
        assert!(rel(a.h, b.h) < tol, "{label}: h {} vs {}", a.h, b.h);
        assert!(rel(a.omega, b.omega) < tol || (a.e < 1e-6 && b.e < 1e-6),
            "{label}: ω {} vs {}", a.omega, b.omega);
        assert_eq!(a.orbit_type, b.orbit_type, "{label}: orbit_type");
    }

    #[test]
    fn invariants_roundtrip_circular() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let bodies = vec![p, s];
        let direct = compute_elements(&bodies, 1, 0, G).unwrap();
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let via_inv = elements_from_invariants(&inv, 0, gm);
        assert_elements_equiv(&direct, &via_inv, "circular");
    }

    #[test]
    fn invariants_roundtrip_eccentric() {
        // e = 0.5 ellipse
        let r_peri = 10.0;
        let m = 1e6;
        let gm = G * m;
        let v_peri = (gm * 1.5 / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let bodies = vec![primary, satellite];
        let direct = compute_elements(&bodies, 1, 0, G).unwrap();
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let via_inv = elements_from_invariants(&inv, 0, gm);
        assert_elements_equiv(&direct, &via_inv, "eccentric");
    }

    /// Anchoring through the body must produce an ellipse that contains
    /// the body geometrically: at the body's azimuth from focus, the
    /// orbit's r equals the body's |r| exactly.
    #[test]
    fn anchored_ellipse_contains_body_prograde() {
        let r_peri = 10.0;
        let m = 1e6;
        let gm = G * m;
        let v_peri = (gm * 1.5 / r_peri).sqrt();
        let bodies =
            vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, v_peri, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].x - bodies[0].x;
        let ry = bodies[1].y - bodies[0].y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        assert!((r - r_orbit).abs() / r < 1e-12,
            "body off orbit: r={r}, r_orbit={r_orbit}");
    }

    /// Retrograde body: sign(h·(r·v)) governs which side of the orbit
    /// the body sits. The check is identical to prograde.
    #[test]
    fn anchored_ellipse_contains_body_retrograde() {
        let r_peri = 10.0;
        let m = 1e6;
        let gm = G * m;
        let v_peri = (gm * 1.5 / r_peri).sqrt();
        // Sign of v_y inverted → retrograde
        let bodies =
            vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, -v_peri, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].x - bodies[0].x;
        let ry = bodies[1].y - bodies[0].y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        assert!((r - r_orbit).abs() / r < 1e-12,
            "retrograde body off orbit: r={r}, r_orbit={r_orbit}");
    }

    /// When fed *smoothed* invariants whose (a, e) drift slightly from
    /// the body's instantaneous orbit, anchoring still places the body
    /// exactly on the displayed ellipse — that is the whole point of
    /// the anchor: shape is averaged, orientation is exact.
    #[test]
    fn anchored_ellipse_contains_body_under_invariant_drift() {
        let r_peri = 10.0;
        let m = 1e6;
        let gm_kep = G * m;
        let v_peri = (gm_kep * 1.5 / r_peri).sqrt();
        let bodies =
            vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, v_peri, 1e-10)];
        // Real invariants of the body's orbit.
        let inv_true = compute_invariants(&bodies, 1, 0, G).unwrap();
        // Simulate smoothed invariants drifted by 1% in energy and h.
        let inv_smooth = OrbitInvariants {
            energy: inv_true.energy * 1.01,
            h: inv_true.h * 0.99,
            ex: inv_true.ex * 1.01,
            ey: inv_true.ey * 0.99,
        };
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv_smooth, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].x - bodies[0].x;
        let ry = bodies[1].y - bodies[0].y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        // Body should still lie on the (drifted-shape) ellipse to FP
        // precision, because anchoring depends only on r and the
        // smoothed (a, e) — it does not require the orbit to be the
        // body's *true* orbit.
        assert!((r - r_orbit).abs() / r < 1e-12,
            "drifted-invariant anchored orbit: r={r}, r_orbit={r_orbit}");
    }

    /// Near-circular bypass: when e < ANCHOR_MIN_E, ω must depend only
    /// on the smoothed eccentricity vector — *not* on the body's
    /// instantaneous position. Moving the body around its orbit while
    /// holding the smoothed invariants fixed must yield byte-identical
    /// elements. This is the contract that kills the low-e wobble seen
    /// on default-circular outer planets in the solar preset.
    #[test]
    fn anchor_low_e_omega_independent_of_body_position() {
        let m = 1e6;
        let r = 10.0;
        let gm = G * m;
        let v_c = (gm / r).sqrt();
        // Slightly perturb a circular orbit to e ≈ 0.01 (well below
        // ANCHOR_MIN_E = 0.05). Anchor must NOT engage here.
        let v = 1.005 * v_c;
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let bodies = vec![primary, body(r, 0.0, 0.0, v, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        assert!(
            (inv.ex * inv.ex + inv.ey * inv.ey).sqrt() < ANCHOR_MIN_E,
            "test setup: e must be below threshold to exercise the bypass"
        );

        let gm_pair = G * (bodies[1].mass + primary.mass);
        let el_at_peri = elements_anchored_to_body(&inv, 0, gm_pair, &bodies[1], &primary);

        // Move the body to a different point along its orbit (any
        // arbitrary location); invariants stay the same since they are
        // the *smoothed* state passed in, not recomputed from the body.
        let displaced = body(0.0, r * 1.1, -v, 0.0, 1e-10);
        let el_displaced = elements_anchored_to_body(&inv, 0, gm_pair, &displaced, &primary);

        assert_eq!(el_at_peri.omega, el_displaced.omega,
            "low-e bypass must produce constant ω across body positions");
        assert_eq!(el_at_peri.a, el_displaced.a);
        assert_eq!(el_at_peri.e, el_displaced.e);
    }

    /// Symmetric check: above the threshold, anchoring DOES engage and
    /// ω moves with the body. Catches accidental over-broadening of the
    /// bypass.
    #[test]
    fn anchor_high_e_omega_tracks_body_position() {
        let m = 1e6;
        let r = 10.0;
        let gm = G * m;
        let v_c = (gm / r).sqrt();
        let v = 1.2 * v_c; // e ≈ 0.44, well above threshold
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let bodies = vec![primary, body(r, 0.0, 0.0, v, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm_pair = G * (bodies[1].mass + primary.mass);
        let el_a = elements_anchored_to_body(&inv, 0, gm_pair, &bodies[1], &primary);

        let displaced = body(0.0, r, -v, 0.0, 1e-10);
        let el_b = elements_anchored_to_body(&inv, 0, gm_pair, &displaced, &primary);

        assert_ne!(el_a.omega, el_b.omega,
            "above threshold, anchor must rotate ω to keep body on orbit");
    }

    #[test]
    fn invariants_roundtrip_hyperbolic() {
        let r_peri = 10.0;
        let m = 1e6;
        let gm = G * m;
        let v_peri = 1.5 * (2.0 * gm / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let bodies = vec![primary, satellite];
        let direct = compute_elements(&bodies, 1, 0, G).unwrap();
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let via_inv = elements_from_invariants(&inv, 0, gm);
        assert_elements_equiv(&direct, &via_inv, "hyperbolic");
    }

    // ── 1. Órbita circular ────────────────────────────────────────────────────
    //
    // Resultado analítico: e = 0, a = r, T = 2π√(r³/GM)

    #[test]
    fn circular_orbit_has_zero_eccentricity() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        assert!(el.e < 1e-6, "e = {}, esperado ≈ 0", el.e);
    }

    #[test]
    fn circular_orbit_semimajor_axis_equals_radius() {
        let r = 10.0;
        let (p, s) = circular_orbit(r, 1e6);
        let el = elements(p, s);
        let err = (el.a - r).abs() / r;
        assert!(err < 1e-6, "a = {}, esperado {r}", el.a);
    }

    #[test]
    fn circular_orbit_period_matches_kepler_iii() {
        // T = 2π √(a³ / GM)
        let r = 10.0;
        let m = 1e6;
        let (p, s) = circular_orbit(r, m);
        let el = elements(p, s);
        let t_expected = TAU * (r.powi(3) / (G * m)).sqrt();
        let err = (el.period - t_expected).abs() / t_expected;
        assert!(err < 1e-6, "T = {}, esperado {t_expected}", el.period);
    }

    #[test]
    fn circular_orbit_is_classified_elliptical() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        assert_eq!(el.orbit_type, OrbitType::Elliptical);
    }

    #[test]
    fn circular_orbit_energy_is_negative() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        assert!(el.energy < 0.0, "energia = {}, deve ser negativa", el.energy);
    }

    // ── 2. Momento angular ────────────────────────────────────────────────────
    //
    // Para CCW: h = r × v > 0
    // Para CW:  h = r × v < 0

    #[test]
    fn ccw_orbit_has_positive_angular_momentum() {
        let (p, s) = circular_orbit(10.0, 1e6); // CCW por construção
        let el = elements(p, s);
        assert!(el.h > 0.0, "h = {}, deve ser positivo (CCW)", el.h);
    }

    #[test]
    fn cw_orbit_has_negative_angular_momentum() {
        let (p, mut s) = circular_orbit(10.0, 1e6);
        s.vy = -s.vy; // inverte para CW
        let el = elements(p, s);
        assert!(el.h < 0.0, "h = {}, deve ser negativo (CW)", el.h);
    }

    #[test]
    fn angular_momentum_magnitude_equals_r_cross_v() {
        let r = 10.0;
        let m = 1e6;
        let v_c = (G * m / r).sqrt();
        let (p, s) = circular_orbit(r, m);
        let el = elements(p, s);
        // h = r × v = r * v_c para órbita circular em (r,0) com v=(0,v_c)
        let h_expected = r * v_c;
        let err = (el.h - h_expected).abs() / h_expected;
        assert!(err < 1e-6, "h = {}, esperado {h_expected}", el.h);
    }

    // ── 3. Kepler III — scaling ───────────────────────────────────────────────
    //
    // T² ∝ a³: dobrar o semi-eixo → período escala por 2^(3/2)

    #[test]
    fn period_scales_with_semimajor_axis_cubed() {
        let m = 1e6;
        let (p1, s1) = circular_orbit(10.0, m);
        let (p2, s2) = circular_orbit(20.0, m);
        let t1 = elements(p1, s1).period;
        let t2 = elements(p2, s2).period;
        let ratio = t2 / t1;
        let expected = 2_f64.powf(1.5); // (20/10)^(3/2) = 2√2
        let err = (ratio - expected).abs() / expected;
        assert!(err < 1e-6, "T2/T1 = {ratio}, esperado {expected}");
    }

    // ── 4. Energia orbital ────────────────────────────────────────────────────
    //
    // ε = -GM / (2a) para qualquer órbita kepleriana

    #[test]
    fn energy_equals_minus_gm_over_2a() {
        let r = 15.0;
        let m = 1e6;
        let (p, s) = circular_orbit(r, m);
        let el = elements(p, s);
        // Para órbita circular a = r
        let energy_expected = -(G * m) / (2.0 * r);
        let err = (el.energy - energy_expected).abs() / energy_expected.abs();
        assert!(err < 1e-6, "ε = {}, esperado {energy_expected}", el.energy);
    }

    // ── 5. Órbitas elípticas ──────────────────────────────────────────────────
    //
    // Órbita elíptica com e conhecido: v no periapsis = √(GM(1+e)/r_peri)

    #[test]
    fn elliptical_orbit_eccentricity_matches_construction() {
        // Construção: no periapsis (r_peri), velocidade v_peri = sqrt(GM(1+e)/r_peri)
        // Para e = 0.5, r_peri = 10 → a = r_peri/(1-e) = 20
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let err = (el.e - e_target).abs();
        assert!(err < 1e-6, "e = {}, esperado {e_target}", el.e);
    }

    #[test]
    fn elliptical_orbit_semimajor_axis_from_periapsis() {
        // a = r_peri / (1 - e)
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let a_expected = r_peri / (1.0 - e_target); // = 20
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let err = (el.a - a_expected).abs() / a_expected;
        assert!(err < 1e-6, "a = {}, esperado {a_expected}", el.a);
    }

    // ── 6. Órbita hiperbólica ─────────────────────────────────────────────────
    //
    // v > v_escape = sqrt(2GM/r) → energia positiva, e > 1

    #[test]
    fn hyperbolic_orbit_has_positive_energy() {
        let r = 10.0;
        let m = 1e6;
        let v_escape = (2.0 * G * m / r).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r, 0.0, 0.0, v_escape * 1.5, 1e-10);
        let el = elements(primary, satellite);
        assert!(el.energy > 0.0, "ε = {}, deve ser positivo", el.energy);
        assert_eq!(el.orbit_type, OrbitType::Hyperbolic);
    }

    #[test]
    fn hyperbolic_orbit_eccentricity_greater_than_one() {
        let r = 10.0;
        let m = 1e6;
        let v_escape = (2.0 * G * m / r).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r, 0.0, 0.0, v_escape * 1.5, 1e-10);
        let el = elements(primary, satellite);
        assert!(el.e > 1.0, "e = {}, deve ser > 1", el.e);
    }

    #[test]
    fn escape_velocity_gives_near_zero_energy() {
        // At exactly v_escape the specific energy is mathematically zero.
        // In f64, (sqrt(2GM/r))^2 ≠ 2GM/r by ~1 ULP, so the energy may land
        // just below or just above −ENERGY_THRESH.  The orbit_type is therefore
        // Elliptical or Parabolic depending on the rounding; what we can test
        // reliably is that |energy| < ENERGY_THRESH × safety_factor.
        let r = 10.0;
        let m = 1e6;
        let v_escape = (2.0 * G * m / r).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r, 0.0, 0.0, v_escape, 1e-10);
        let el = elements(primary, satellite);
        // energy = ½v² − GM/r; at v = v_escape this is ~0 up to FP error.
        assert!(
            el.energy.abs() < 1e-8,
            "energy = {:.3e}, expected |energy| < 1e-8 at escape velocity",
            el.energy
        );
        // The orbit must be bound (Elliptical or Parabolic), never Hyperbolic.
        assert!(
            el.orbit_type != OrbitType::Hyperbolic,
            "orbit should not be Hyperbolic at exactly escape velocity"
        );
    }

    // ── 7. Argumento do periapsis ─────────────────────────────────────────────
    //
    // Órbita no eixo x → ω = 0; rotacionada 90° → ω = π/2

    #[test]
    fn periapsis_on_x_axis_gives_omega_zero() {
        // Periapsis em (r, 0): ω deve ser 0
        let e = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let v_peri = (G * m * (1.0 + e) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        assert!(el.omega.abs() < 1e-6, "ω = {}, esperado 0", el.omega);
    }

    #[test]
    fn periapsis_on_y_axis_gives_omega_pi_over_2() {
        // Periapsis em (0, r): ω deve ser π/2
        let e = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let v_peri = (G * m * (1.0 + e) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        // Rotaciona 90°: posição (0, r_peri), velocidade (-v_peri, 0) para CCW
        let satellite = body(0.0, r_peri, -v_peri, 0.0, 1e-10);
        let el = elements(primary, satellite);
        let err = (el.omega - PI / 2.0).abs();
        assert!(err < 1e-6, "ω = {}, esperado π/2", el.omega);
    }

    // ── 8. Seleção de primário dominante ──────────────────────────────────────
    //
    // Corpo mais próximo deve vencer mesmo sendo menos massivo (se m/r² for maior)

    #[test]
    fn dominant_primary_prefers_closer_body_by_acceleration() {
        // Sol em (-1000, 0), Terra em (1, 0), satélite em (0, 0)
        // m_sol/r² = 1e6/1e6 = 1.0   vs   m_terra/r² = 1e3/1 = 1000 → Terra vence
        let sun = body(-1000.0, 0.0, 0.0, 0.0, 1e6);
        let earth = body(1.0, 0.0, 0.0, 0.0, 1e3);
        let satellite = body(0.0, 0.0, 0.0, 0.0, 1.0);
        let bodies = vec![sun, earth, satellite];
        let primary = dominant_primary(&bodies, 2).unwrap();
        assert_eq!(primary, 1, "primário deve ser a Terra (índice 1), não o Sol");
    }

    #[test]
    fn dominant_primary_prefers_massive_body_when_equidistant() {
        // Dois corpos equidistantes: o mais massivo deve vencer
        let heavy = body(-10.0, 0.0, 0.0, 0.0, 1e6);
        let light = body(10.0, 0.0, 0.0, 0.0, 1e3);
        let probe = body(0.0, 0.0, 0.0, 0.0, 1.0);
        let bodies = vec![heavy, light, probe];
        let primary = dominant_primary(&bodies, 2).unwrap();
        assert_eq!(primary, 0, "primário deve ser o corpo mais massivo");
    }

    #[test]
    fn dominant_primary_returns_none_for_single_body() {
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, 1e6)];
        assert!(dominant_primary(&bodies, 0).is_none());
    }

    // ── 9. Conservação — vis-viva ─────────────────────────────────────────────
    //
    // v² = GM(2/r − 1/a) em qualquer ponto da órbita

    // ── 10. sample_orbit — geometry ───────────────────────────────────────────

    #[test]
    fn sample_orbit_closes_the_loop() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        assert_eq!(pts.len(), 65, "should return steps+1 points");
        let first = pts.first().unwrap();
        let last = pts.last().unwrap();
        for k in 0..3 {
            assert!(
                (first[k] - last[k]).abs() < 1e-10,
                "loop must close exactly on axis {k}: first={:?}, last={:?}",
                first,
                last,
            );
        }
    }

    #[test]
    fn sample_orbit_circular_all_points_equidistant_from_focus() {
        let r = 10.0;
        let (p, s) = circular_orbit(r, 1e6);
        let el = elements(p, s);
        // Primary at origin in this fixture.
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 128);
        for pt in &pts {
            let d = (pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt();
            let err = (d - r).abs() / r;
            assert!(err < 1e-6, "distance {d} should equal {r}");
        }
    }

    #[test]
    fn sample_orbit_ellipse_extremes_match_r_peri_r_apo() {
        // e = 0.5, r_peri = 10 → a = 20, r_apo = 30
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 256);
        let dists: Vec<f64> = pts.iter().map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt()).collect();
        let r_min = dists.iter().cloned().fold(f64::INFINITY, f64::min);
        let r_max = dists.iter().cloned().fold(0.0_f64, f64::max);
        let r_apo_expected = r_peri * (1.0 + e_target) / (1.0 - e_target); // 30
        assert!((r_min - r_peri).abs() / r_peri < 1e-4, "r_min = {r_min}, expected {r_peri}",);
        assert!(
            (r_max - r_apo_expected).abs() / r_apo_expected < 1e-4,
            "r_max = {r_max}, expected {r_apo_expected}",
        );
    }

    #[test]
    fn sample_orbit_is_focus_centred_not_geometry_centred() {
        // For e > 0, the primary (focus) is offset from the geometric centre
        // by a·e. If we (incorrectly) centred on the geometric centre, the
        // focus-distance extremes would become (a − a·e) and (a + a·e) only
        // when measured from the geometric centre — not from the focus.
        //
        // Direct test: with primary at origin, the closest sampled point
        // must sit at distance a(1-e) from the origin, NOT at distance a.
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 128);
        let r_min =
            pts.iter().map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt()).fold(f64::INFINITY, f64::min);
        let a = el.a; // 20
        // Focus-centred: r_min = a(1-e) = 10.
        // Geometry-centred: r_min would be a(1-e) relative to geometric
        // centre, but geometric centre is at focus + a·e on major axis,
        // giving r_min from origin = a(1-e) + a·e = a = 20.
        // So a passing test at r_peri = 10 rules out the geometric-centre bug.
        assert!(
            (r_min - r_peri).abs() / r_peri < 1e-4,
            "r_min = {r_min}, should equal r_peri = {r_peri} (focus-centred), \
             would be {a} if geometry-centred",
        );
    }

    #[test]
    fn sample_orbit_translates_with_primary() {
        let r = 10.0;
        let (p, s) = circular_orbit(r, 1e6);
        let el = elements(p, s);
        let shift = [42.5_f64, -17.25, 0.0];
        let a = el.sample_orbit([0.0, 0.0, 0.0], 32);
        let b = el.sample_orbit(shift, 32);
        assert_eq!(a.len(), b.len());
        for (pa, pb) in a.iter().zip(b.iter()) {
            for k in 0..3 {
                let expected = pa[k] + shift[k];
                assert!(
                    (pb[k] - expected).abs() < 1e-9,
                    "axis {k}: got {}, expected {expected}",
                    pb[k],
                );
            }
        }
    }

    // ── Hyperbolic sampling ───────────────────────────────────────────────────

    /// Builds a canonical hyperbolic flyby: body passing periapsis on +x,
    /// moving +y at speed v > v_escape. Returns (primary, satellite).
    fn hyperbolic_flyby(r_peri: f64, v_multiplier: f64, primary_mass: f64) -> (Body, Body) {
        let gm = G * primary_mass;
        let v_peri = v_multiplier * (2.0 * gm / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        (primary, satellite)
    }

    #[test]
    fn sample_orbit_hyperbolic_returns_open_polyline() {
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let el = elements(p, s);
        assert_eq!(el.orbit_type, OrbitType::Hyperbolic);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        assert_eq!(pts.len(), 65, "should return steps+1 points");
        // Open polyline — first and last are NOT equal.
        let first = pts.first().unwrap();
        let last = pts.last().unwrap();
        let sep = ((first[0] - last[0]).powi(2) + (first[1] - last[1]).powi(2)).sqrt();
        assert!(sep > 1.0, "hyperbolic polyline must not close (got sep = {sep})");
    }

    #[test]
    fn sample_orbit_hyperbolic_middle_point_is_at_periapsis() {
        let r_peri = 10.0;
        let m = 1e6;
        let (p, s) = hyperbolic_flyby(r_peri, 1.5, m);
        let el = elements(p, s);
        // H = 0 at k = steps/2. Use even steps so the midpoint exists.
        let steps = 64usize;
        let pts = el.sample_orbit([0.0, 0.0, 0.0], steps);
        let mid = pts[steps / 2];
        let d = (mid[0] * mid[0] + mid[1] * mid[1]).sqrt();
        let err = (d - r_peri).abs() / r_peri;
        assert!(err < 1e-6, "midpoint distance = {d}, expected r_peri = {r_peri}");
        // And the point lies on +x (ω = 0 for this fixture).
        assert!(mid[0] > 0.0, "periapsis should have +x sign, got {}", mid[0]);
        assert!(mid[1].abs() < 1e-6, "periapsis should have y ≈ 0, got {}", mid[1]);
    }

    #[test]
    fn sample_orbit_hyperbolic_distance_grows_monotonically_from_midpoint() {
        // Because r = |a|(e cosh H − 1), distance to focus strictly
        // increases as |H| increases from 0.
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let el = elements(p, s);
        let steps = 64usize;
        let pts = el.sample_orbit([0.0, 0.0, 0.0], steps);
        let mid = steps / 2;
        let dist = |k: usize| (pts[k][0].powi(2) + pts[k][1].powi(2)).sqrt();
        // Backward half: dist(0) > dist(1) > … > dist(mid)
        for k in 0..mid {
            assert!(
                dist(k) > dist(k + 1),
                "backward monotonicity broken at k={k}: {} vs {}",
                dist(k),
                dist(k + 1),
            );
        }
        // Forward half: dist(mid) < dist(mid+1) < …
        for k in mid..steps {
            assert!(
                dist(k) < dist(k + 1),
                "forward monotonicity broken at k={k}: {} vs {}",
                dist(k),
                dist(k + 1),
            );
        }
    }

    #[test]
    fn sample_orbit_hyperbolic_respects_omega_rotation() {
        // Rotate the fixture 90°: periapsis lands on +y, not +x.
        let r_peri = 10.0;
        let m = 1e6;
        let gm = G * m;
        let v_peri = 1.5 * (2.0 * gm / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        // Position (0, r_peri), velocity (-v_peri, 0) → CCW, periapsis on +y.
        let satellite = body(0.0, r_peri, -v_peri, 0.0, 1e-10);
        let el = elements(primary, satellite);
        let steps = 64usize;
        let pts = el.sample_orbit([0.0, 0.0, 0.0], steps);
        let mid = pts[steps / 2];
        assert!(mid[0].abs() < 1e-6, "rotated peri.x should be 0, got {}", mid[0]);
        assert!((mid[1] - r_peri).abs() < 1e-6, "rotated peri.y = {}", mid[1]);
    }

    #[test]
    fn sample_orbit_hyperbolic_translates_with_primary() {
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let el = elements(p, s);
        let shift = [42.5_f64, -17.25, 0.0];
        let a = el.sample_orbit([0.0, 0.0, 0.0], 32);
        let b = el.sample_orbit(shift, 32);
        assert_eq!(a.len(), b.len());
        for (pa, pb) in a.iter().zip(b.iter()) {
            for k in 0..3 {
                let expected = pa[k] + shift[k];
                assert!(
                    (pb[k] - expected).abs() < 1e-9,
                    "axis {k}: got {}, expected {expected}",
                    pb[k],
                );
            }
        }
    }

    // ── Apsis accessors ────────────────────────────────────────────────────────

    #[test]
    fn periapsis_world_ellipse_lies_on_plus_x() {
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let p = el.periapsis_world([0.0, 0.0, 0.0]).unwrap();
        assert!((p[0] - r_peri).abs() < 1e-6, "peri.x = {}, expected {r_peri}", p[0]);
        assert!(p[1].abs() < 1e-6, "peri.y = {}", p[1]);
    }

    #[test]
    fn apoapsis_world_ellipse_matches_r_apo_on_minus_x() {
        // e = 0.5, r_peri = 10 → a = 20, r_apo = 30 on −x side of focus.
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let r_apo = r_peri * (1.0 + e_target) / (1.0 - e_target); // = 30
        let p = el.apoapsis_world([0.0, 0.0, 0.0]).unwrap();
        assert!((p[0] + r_apo).abs() < 1e-6, "apo.x = {}, expected {}", p[0], -r_apo);
        assert!(p[1].abs() < 1e-6, "apo.y = {}", p[1]);
    }

    #[test]
    fn periapsis_world_hyperbola_matches_r_peri() {
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let el = elements(p, s);
        let pt = el.periapsis_world([0.0, 0.0, 0.0]).unwrap();
        let d = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
        assert!((d - 10.0).abs() < 1e-6, "hyper peri distance = {d}");
        assert!(pt[0] > 0.0, "hyper peri should sit on +x");
    }

    #[test]
    fn apoapsis_world_is_none_for_hyperbola() {
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let el = elements(p, s);
        assert!(el.apoapsis_world([0.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn apsides_translate_with_primary_position() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        let shift = [3.0_f64, -4.0, 0.0];
        let peri_a = el.periapsis_world([0.0, 0.0, 0.0]).unwrap();
        let peri_b = el.periapsis_world(shift).unwrap();
        for k in 0..3 {
            assert!((peri_b[k] - peri_a[k] - shift[k]).abs() < 1e-9);
        }
    }

    #[test]
    fn apsides_respect_omega_rotation() {
        // Periastro rotated 90° onto +y; apoapsis onto −y.
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(0.0, r_peri, -v_peri, 0.0, 1e-10);
        let el = elements(primary, satellite);
        let peri = el.periapsis_world([0.0, 0.0, 0.0]).unwrap();
        assert!(peri[0].abs() < 1e-6);
        assert!((peri[1] - r_peri).abs() < 1e-6);
        let apo = el.apoapsis_world([0.0, 0.0, 0.0]).unwrap();
        let r_apo = r_peri * (1.0 + e_target) / (1.0 - e_target);
        assert!(apo[0].abs() < 1e-6);
        assert!((apo[1] + r_apo).abs() < 1e-6);
    }

    #[test]
    fn sample_orbit_parabolic_returns_empty() {
        // Near-parabolic: v ≈ v_escape puts energy near zero; type may
        // round to Elliptical or Parabolic depending on ULP. Force the
        // type explicitly on a computed element to test the branch.
        let (p, s) = hyperbolic_flyby(10.0, 1.5, 1e6);
        let mut el = elements(p, s);
        el.orbit_type = OrbitType::Parabolic;
        el.a = f64::INFINITY;
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        assert!(pts.is_empty(), "parabolic must not sample");
    }

    #[test]
    fn sample_orbit_rejects_too_few_steps() {
        let (p, s) = circular_orbit(10.0, 1e6);
        let el = elements(p, s);
        assert!(el.sample_orbit([0.0, 0.0, 0.0], 0).is_empty());
        assert!(el.sample_orbit([0.0, 0.0, 0.0], 1).is_empty());
    }

    #[test]
    fn sample_orbit_periastro_at_omega_zero_lies_on_positive_x() {
        // ω = 0: periastro in the +x direction from the focus.
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        let el = elements(primary, satellite);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        // E = 0 is the first sample → periastro.
        let peri = pts[0];
        assert!((peri[0] - r_peri).abs() < 1e-6, "peri.x = {}", peri[0]);
        assert!(peri[1].abs() < 1e-6, "peri.y = {}", peri[1]);
        assert!(peri[2].abs() < 1e-12, "peri.z = {}", peri[2]);
    }

    #[test]
    fn sample_orbit_rotates_with_omega() {
        // Rotate the fixture by 90°: periastro must land on +y.
        let e_target = 0.5_f64;
        let r_peri = 10.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let v_peri = (gm * (1.0 + e_target) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(0.0, r_peri, -v_peri, 0.0, 1e-10);
        let el = elements(primary, satellite);
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        let peri = pts[0];
        assert!(peri[0].abs() < 1e-6, "peri.x = {}", peri[0]);
        assert!((peri[1] - r_peri).abs() < 1e-6, "peri.y = {}", peri[1]);
    }

    #[test]
    fn sample_orbit_inclination_pi_over_2_puts_orbit_in_xz_plane() {
        // Build a circular 2D orbit, then force i = π/2 and re-sample.
        // With Ω = 0, ω = 0, rotating by i about the line of nodes (x-axis)
        // sends y → z. Result: every sampled point must have y ≈ 0.
        let (p, s) = circular_orbit(10.0, 1e6);
        let mut el = elements(p, s);
        el.inclination = PI / 2.0;
        el.lon_ascending_node = 0.0;
        el.omega = 0.0;
        let pts = el.sample_orbit([0.0, 0.0, 0.0], 64);
        for pt in &pts {
            assert!(pt[1].abs() < 1e-6, "y = {} should be 0 for orbit in xz-plane", pt[1],);
        }
    }

    // ── 11. vis-viva at apoapsis (pre-existing) ───────────────────────────────

    #[test]
    fn vis_viva_holds_at_apoapsis() {
        // Construção via periapsis, depois coloca o satélite no apoapsis
        // r_apo = a(1+e), v_apo = sqrt(GM(1-e)/r_apo)
        let e = 0.6_f64;
        let r_peri = 5.0_f64;
        let m = 1e6_f64;
        let gm = G * m;
        let a = r_peri / (1.0 - e);
        let r_apo = a * (1.0 + e);
        let v_apo = (gm * (1.0 - e) / r_apo).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = body(-r_apo, 0.0, 0.0, -v_apo, 1e-10); // apoapsis em (-r_apo, 0)
        let el = elements(primary, satellite);
        // vis-viva: v² = GM(2/r - 1/a)
        let v2 = v_apo * v_apo;
        let v2_visviva = gm * (2.0 / r_apo - 1.0 / el.a);
        let err = (v2 - v2_visviva).abs() / v2;
        assert!(err < 1e-6, "vis-viva: v² = {v2}, esperado {v2_visviva}");
    }
}
