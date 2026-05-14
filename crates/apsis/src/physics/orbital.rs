//! Osculating orbital element computation for 3D N-body systems.
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
//! Two distinct notions of "primary" coexist for an N-body snapshot:
//!
//! * [`dominant_primary`] — the body `j ≠ i` whose instantaneous
//!   gravitational pull on body `i` is largest, `argmax_j G·m_j / r_ij²`.
//!   This is the natural reference for osculating two-body elements and
//!   is what [`compute_elements`] uses.
//!
//! * [`hierarchical_primary`] — the body in whose Hill sphere body `i`
//!   currently sits. For a Moon–Earth–Sun configuration the dominant
//!   primary of the Moon is the **Sun** (`G·M☉ / r_⊙^2 ≈ 5.9 × 10⁻³ m/s²`
//!   exceeds Earth's `≈ 2.7 × 10⁻³ m/s²`), while the hierarchical primary
//!   is **Earth**. The two coincide for non-hierarchical configurations
//!   such as Mercury orbiting the Sun directly.
//!
//! Both are pure snapshot helpers: nothing in the integrator consumes
//! either, and neither is cached across steps. Renderers and analysis
//! consumers may compute them on demand.
//!
//! ## Computed quantities
//!
//! | Symbol | Name | Valid when |
//! |--------|------|-----------|
//! | `a`    | semi-major axis | bound orbit (e < 1) |
//! | `e`    | eccentricity | always |
//! | `T`    | period | bound orbit |
//! | `h_vec`| specific angular momentum vector | always |
//! | `ε`    | specific orbital energy | always |
//! | `i`    | inclination relative to the simulation's `ẑ` axis | always |
//! | `Ω`    | longitude of ascending node | `\|n\| > N_DEGENERATE_EPS` |
//! | `ω`    | argument of periapsis | `e > E_CIRCULAR_EPS` (and inclined branch when `\|n\| > N_DEGENERATE_EPS`) |
//! | `ν`    | true anomaly | always (degenerate to argument of latitude when `e < E_CIRCULAR_EPS`) |
//! | `E`/`H`| eccentric anomaly (`E` elliptical, hyperbolic `H` unbound) | bound or hyperbolic; NaN parabolic |
//! | `M`    | mean anomaly (linearly time-evolving) | bound or hyperbolic; NaN parabolic |
//! | `q`,`Q`| pericenter / apocenter distance | helper methods on [`OrbitalElements`] |
//! | `n`    | mean motion `2π / T` | helper method on [`OrbitalElements`] |
//!
//! ## Frame and convention alignment
//!
//! The reference frame is the **inertial right-handed Cartesian** axes
//! `{x̂, ŷ, ẑ}` of the simulation. `ẑ` is the inclination axis: `i = 0`
//! means the orbital plane coincides with the `xy`-plane and the orbit
//! is prograde with respect to `+ẑ`; `i = π` is retrograde in the same
//! plane. The line of nodes is `n = ẑ × h_vec`, lying in the `xy`-plane.
//!
//! Element definitions, frame, and singularity handling follow the
//! standard astrodynamics convention as implemented by REBOUND
//! (`tools.c::reb_tools_orbit_to_particle` and inverse). The single
//! deliberate divergence is the angular range — `[-π, π]` here vs
//! `[0, 2π]` in REBOUND — required to keep the locked planar baseline
//! bit-equivalent across the 3D port.
//!
//! Apsis-specific surface (kernel `Exactness`/`Continuity` invariants,
//! per-operator `KernelRequirements`, federated extension crates)
//! lives in `physics::gravity::kernel` and
//! `physics::integrator::operator`. This module is intentionally vanilla.
//!
//! ## Singularity contracts
//!
//! Two thresholds gate the angular calculations and are documented at
//! their constants:
//!
//! * `E_CIRCULAR_EPS` — eccentricity below which `ω` is undefined
//!   (circular orbit). Convention: `ω = 0`.
//! * `N_DEGENERATE_EPS` — node-vector magnitude below which `Ω` and the
//!   inclined `ω` formulation are undefined. Triggers for both `i ≈ 0`
//!   (planar prograde) and `i ≈ π` (planar retrograde) since `\|ẑ × h\| =
//!   \|h\|·\|sin i\|` vanishes in both. Convention: `Ω = 0`,
//!   `ω = atan2(e_vec.y, e_vec.x)`.
//!
//! Angular continuity across step boundaries — the `π → −π` wrap that
//! visualisations need to unwrap for smooth animation — is the
//! responsibility of the presentation layer (e.g. `apsis-app`'s
//! `orbit_smoother`), not this module. The element values returned here
//! are always principal-value angles in `[-π, π]` enforced by
//! [`crate::math::wrap_pi`].
//!
//! ## References
//! - Murray & Dermott (1999). *Solar System Dynamics*. Cambridge.
//! - Bate, Mueller & White (1971). *Fundamentals of Astrodynamics*. Dover.
//! - Vallado (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed.

use std::f64::consts::TAU;

use crate::domain::body::Body;
use crate::math::{Vec3, wrap_pi};

/// Eccentricity below which `ω` is undefined (circular orbit). Returns 0.
const E_CIRCULAR_EPS: f64 = 1e-6;

/// Node-vector magnitude below which `Ω` and the inclined `ω` form are
/// undefined. Triggers for both prograde-planar (`i ≈ 0`) and
/// retrograde-planar (`i ≈ π`) orbits, since `n = ẑ × h_vec` vanishes
/// in both. `1e-12` sits well above the f64 round-off floor for the
/// cross-product evaluation while still capturing genuine planar input
/// where `h_vec.x` and `h_vec.y` are exactly zero.
const N_DEGENERATE_EPS: f64 = 1e-12;

/// CSV schema version for [`OrbitalElements::csv_header`] and
/// [`OrbitalElements::to_csv_row`].
///
/// **Version 1** (legacy): `t, body_idx, primary_idx, a, e, period, h,
/// energy, omega_deg, orbit_type` — 10 columns.
///
/// **Version 2** (current): adds `inclination_deg, lon_asc_node_deg,
/// true_anomaly_deg, eccentric_anomaly_deg, mean_anomaly_deg, pericenter,
/// apocenter` — 17 columns total. External tooling parsing the apsis
/// CSV format must check this constant before consuming the file; v1
/// readers will silently misalign on v2 output.
pub const CSV_SCHEMA_VERSION: u32 = 2;

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
/// # Frame
///
/// Elements live in the simulation's inertial Cartesian frame `{x̂, ŷ, ẑ}`.
/// `inclination` is measured against `+ẑ`; `lon_ascending_node` against
/// `+x̂`; both are populated from `h_vec` and the ascending-node geometry.
/// For planar configurations (`|n| < N_DEGENERATE_EPS`) `Ω` collapses to
/// 0 by convention and `ω` falls back to `atan2(e_vec.y, e_vec.x)`.
///
/// # Anomalies
///
/// `true_anomaly`, `eccentric_anomaly`, and `mean_anomaly` characterise
/// the body's instantaneous phase on the orbit. They satisfy Kepler's
/// equation `M = E − e sin E` for elliptical orbits, with the analogous
/// hyperbolic identity `M_h = e sinh H − H` for unbound flybys. The mean
/// anomaly is the linearly time-evolving angle `M = n (t − t_peri)` where
/// `n` is the mean motion (see [`Self::mean_motion`]).
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

    /// Specific angular momentum vector `r × v`.
    ///
    /// Magnitude `|h_vec|` is the conserved quantity in two-body
    /// motion; direction defines the orbital plane. The signed
    /// z-component `h_vec.z` is positive for orbits prograde with
    /// respect to `+ẑ`.
    pub h_vec: Vec3,

    /// Specific orbital energy (`ε = v²/2 − GM/r`).
    /// Negative = bound, positive = unbound.
    pub energy: f64,

    /// Argument of periapsis ω ∈ [−π, π] (radians).
    /// Undefined when `e < 1e-6` (circular); returns 0.
    pub omega: f64,

    /// Inclination i ∈ [0, π] (radians).
    pub inclination: f64,

    /// Longitude of the ascending node Ω ∈ [−π, π] (radians).
    /// Returns 0 when the line of nodes is degenerate (`i ≈ 0` or `i ≈ π`).
    pub lon_ascending_node: f64,

    /// True anomaly ν ∈ [−π, π] (radians) — body's angular position on the
    /// orbit measured from periapsis, in the direction of motion. For
    /// circular orbits (`e < E_CIRCULAR_EPS`), reduces to the argument of
    /// latitude (angle from ascending node). NaN for parabolic.
    pub true_anomaly: f64,

    /// Eccentric anomaly. For elliptical orbits, `E ∈ [−π, π]` (radians)
    /// related to `ν` by `tan(E/2) = √((1−e)/(1+e)) tan(ν/2)`. For
    /// hyperbolic orbits, hyperbolic anomaly `H ∈ ℝ` related by
    /// `tanh(H/2) = √((e−1)/(e+1)) tan(ν/2)`. NaN for parabolic.
    ///
    /// The half-angle recovery is numerically robust for `e < 0.999`.
    /// Beyond that the `(1−e).sqrt()` factor enters the round-off floor
    /// of the source ellipse and an asymptotic branch (parabolic limit)
    /// would be needed; no current caller exercises that regime.
    pub eccentric_anomaly: f64,

    /// Mean anomaly.
    ///
    /// **Elliptical:** `M = E − e sin E`, naturally landing in
    /// `[−π − e, π + e]`. The value is **not** wrapped to `[−π, π]` —
    /// `M` is semantically a temporal variable (`M(t) = M₀ + n t`) and an
    /// artificial principal-value collapse would corrupt consumers
    /// tracking `dM/dt` across frames. Snapshot recovery from
    /// instantaneous `(r, v)` cannot supply absolute temporal context;
    /// the `2π` discontinuity at orbit closure must be unwrapped by the
    /// consumer with a continuity-tracking layer if monotonicity in time
    /// is required.
    ///
    /// **Hyperbolic:** `M_h = e sinh H − H ∈ ℝ`. Already unbounded.
    ///
    /// **Parabolic:** NaN. Detection threshold is `|energy| < 1e-12`
    /// (see `ENERGY_THRESH` inside [`elements_from_invariants`]); orbits
    /// inside that band are flagged but not characterised. Barker's
    /// variable `D` (the parabolic anomaly) is reserved for the day a
    /// caller renders parabolic flybys.
    pub mean_anomaly: f64,

    /// Orbit classification derived from `energy`.
    pub orbit_type: OrbitType,
}

impl OrbitalElements {
    /// Returns `true` for elliptical and parabolic orbits.
    pub fn is_bound(self) -> bool {
        self.orbit_type.is_bound()
    }

    /// Pericenter distance `q`.
    ///
    /// Elliptical: `q = a(1 − e)`. Hyperbolic: `q = |a|(e − 1)` (positive,
    /// the closest-approach distance on the incoming branch). Parabolic
    /// returns NaN — the canonical parabolic parameter is `q = h² / (2GM)`,
    /// which this struct does not carry.
    pub fn pericenter(self) -> f64 {
        match self.orbit_type {
            OrbitType::Elliptical => self.a * (1.0 - self.e),
            OrbitType::Hyperbolic => self.a.abs() * (self.e - 1.0),
            OrbitType::Parabolic => f64::NAN,
        }
    }

    /// Apocenter distance `Q = a(1 + e)`. Defined only for elliptical
    /// orbits; hyperbolic and parabolic flybys have no bounded apoapsis
    /// and return NaN.
    pub fn apocenter(self) -> f64 {
        match self.orbit_type {
            OrbitType::Elliptical => self.a * (1.0 + self.e),
            _ => f64::NAN,
        }
    }

    /// Mean motion `n = 2π / T`. Returns NaN for unbound orbits.
    pub fn mean_motion(self) -> f64 {
        if self.period.is_finite() && self.period > 0.0 { TAU / self.period } else { f64::NAN }
    }

    /// CSV header row matching [`Self::to_csv_row`].
    ///
    /// Schema is versioned — see [`CSV_SCHEMA_VERSION`]. The current
    /// header reflects v2 (17 columns); external parsers built against
    /// the v1 format (10 columns) will silently misalign and must check
    /// the version constant before consuming the file.
    pub fn csv_header() -> &'static str {
        "t,body_idx,primary_idx,a,e,period,h,energy,\
         inclination_deg,lon_asc_node_deg,omega_deg,\
         true_anomaly_deg,eccentric_anomaly_deg,mean_anomaly_deg,\
         pericenter,apocenter,orbit_type"
    }

    /// Serialise to a CSV data row.  `t` and `body_idx` are injected by the caller.
    pub fn to_csv_row(self, t: f64, body_idx: usize) -> String {
        format!(
            "{t:.6e},{body_idx},{},\
             {:.6e},{:.6e},{:.6e},{:.6e},{:.6e},\
             {:.4},{:.4},{:.4},\
             {:.4},{:.4},{:.4},\
             {:.6e},{:.6e},{:?}",
            self.primary_idx,
            self.a,
            self.e,
            self.period,
            self.h_vec.z,
            self.energy,
            self.inclination.to_degrees(),
            self.lon_ascending_node.to_degrees(),
            self.omega.to_degrees(),
            self.true_anomaly.to_degrees(),
            self.eccentric_anomaly.to_degrees(),
            self.mean_anomaly.to_degrees(),
            self.pericenter(),
            self.apocenter(),
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
            let rj2 = (bj.pos_x - bi.pos_x).powi(2)
                + (bj.pos_y - bi.pos_y).powi(2)
                + (bj.pos_z - bi.pos_z).powi(2);
            let rk2 = (bk.pos_x - bi.pos_x).powi(2)
                + (bk.pos_y - bi.pos_y).powi(2)
                + (bk.pos_z - bi.pos_z).powi(2);
            // Compare m_j/r_j² vs m_k/r_k²  (G cancels)
            let score_j = if rj2 > 0.0 { bj.mass / rj2 } else { f64::INFINITY };
            let score_k = if rk2 > 0.0 { bk.mass / rk2 } else { f64::INFINITY };
            score_j.partial_cmp(&score_k).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(j, _)| j)
}

/// Returns `true` when no body in the snapshot is strictly more massive
/// than `bodies[idx]`.
///
/// System root has no Keplerian orbit; rendering one would misrepresent
/// N-body dynamics — the most massive body's motion is the integrated
/// reflex of every other body, not a closed conic around any single
/// reference. Inspector and canvas overlays both consult this predicate
/// to skip drawing an orbit they have no honest way to construct.
pub fn is_system_root(bodies: &[Body], idx: usize) -> bool {
    if idx >= bodies.len() {
        return false;
    }
    let m = bodies[idx].mass;
    !bodies.iter().enumerate().any(|(j, b)| j != idx && b.mass > m)
}

/// How the hierarchical relationship was established — surfaced to the
/// caller so consumers can show *why* the parent was chosen rather than
/// asserting it as fact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HierarchicalRelation {
    /// Selected body sits inside the candidate's Hill sphere — the
    /// canonical capture criterion.
    HillSphere,
    /// Hill-sphere check failed, but the selected body is gravitationally
    /// bound (`ε < 0`) to the candidate; smallest such candidate wins.
    Energy,
}

/// Hierarchical primary — the body in whose Hill sphere `idx` currently
/// resides, paired with the criterion that established the relationship.
/// Distinct from [`dominant_primary`] (the strongest attractor): for the
/// Earth-Moon-Sun configuration the dominant primary of the Moon is the
/// Sun while the hierarchical primary is Earth.
///
/// # Algorithm
///
/// Hill-sphere check, smallest-containing wins. For each candidate
/// `cand` with `mass > bodies[idx].mass`, sorted ascending by mass:
///
/// * Treat `cand` as a "system root" when no body in the snapshot is
///   strictly more massive than it. The root's Hill sphere is infinite —
///   anything not contained by a smaller candidate falls back to it.
/// * Otherwise, compute the Hill radius
///   `r_H = a · ∛(m_cand / 3·m_parent)` using `cand`'s **current**
///   distance to its [`dominant_primary`] as `a` (a snapshot proxy for
///   the semi-major axis). Return `cand` if `idx` lies inside that Hill
///   radius.
///
/// When no candidate's Hill sphere contains `idx`, fall back to an
/// energy test: pick the smallest-mass candidate around which `idx` is
/// gravitationally bound (`ε < 0`). This catches partial captures and
/// chaotic configurations where Hill-sphere reasoning fails.
///
/// # Pure snapshot
///
/// This function consumes only the bodies' current `(r, v, m)`. It is
/// not memoised, not consulted by the integrator, and carries no
/// implication for dynamics — apsis remains a direct N-body code, with
/// the hierarchy used solely as a render/analysis lens.
///
/// Returns `None` when there is no body more massive than `idx` and the
/// energy fallback finds no bound candidate.
pub fn hierarchical_primary(bodies: &[Body], idx: usize) -> Option<(usize, HierarchicalRelation)> {
    if bodies.len() < 2 {
        return None;
    }
    let bi = &bodies[idx];

    let mut candidates: Vec<usize> =
        (0..bodies.len()).filter(|&j| j != idx && bodies[j].mass > bi.mass).collect();
    candidates.sort_by(|&a, &b| {
        bodies[a].mass.partial_cmp(&bodies[b].mass).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Hill-sphere primary — smallest-containing wins.
    for &cand in &candidates {
        let bc = &bodies[cand];
        let is_root = bodies.iter().enumerate().all(|(j, b)| j == cand || b.mass <= bc.mass);

        let r_hill = if is_root {
            f64::INFINITY
        } else {
            let parent_idx = match dominant_primary(bodies, cand) {
                Some(p) => p,
                None => continue,
            };
            let bp = &bodies[parent_idx];
            let a_cp = ((bc.pos_x - bp.pos_x).powi(2)
                + (bc.pos_y - bp.pos_y).powi(2)
                + (bc.pos_z - bp.pos_z).powi(2))
            .sqrt();
            if a_cp < 1e-15 || bp.mass < 1e-30 {
                continue;
            }
            a_cp * (bc.mass / (3.0 * bp.mass)).cbrt()
        };

        let r_self = ((bi.pos_x - bc.pos_x).powi(2)
            + (bi.pos_y - bc.pos_y).powi(2)
            + (bi.pos_z - bc.pos_z).powi(2))
        .sqrt();
        if r_self < r_hill {
            return Some((cand, HierarchicalRelation::HillSphere));
        }
    }

    // Energy fallback — pick the smallest bound candidate.
    for &cand in &candidates {
        let bc = &bodies[cand];
        let dx = bi.pos_x - bc.pos_x;
        let dy = bi.pos_y - bc.pos_y;
        let dz = bi.pos_z - bc.pos_z;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-15 {
            continue;
        }
        let dvx = bi.vel_x - bc.vel_x;
        let dvy = bi.vel_y - bc.vel_y;
        let dvz = bi.vel_z - bc.vel_z;
        let v2 = dvx * dvx + dvy * dvy + dvz * dvz;
        // Sign of the orbital energy proxy `½v² − m_cand/r`. The constant
        // `G` would scale both terms identically and is irrelevant to the
        // sign comparison; bound iff this is negative.
        if 0.5 * v2 - bc.mass / r < 0.0 {
            return Some((cand, HierarchicalRelation::Energy));
        }
    }

    None
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
    /// Specific angular momentum vector `h_vec = r × v`.
    pub h_vec: Vec3,
    /// Eccentricity (Laplace–Runge–Lenz) vector
    /// `e_vec = (v × h)/μ − r̂`. Magnitude is the orbit's eccentricity;
    /// direction points from the focus toward the periapsis.
    pub e_vec: Vec3,
    /// Body's position relative to the primary, `r_body − r_primary`.
    /// Required by anomaly recovery in [`elements_from_invariants`].
    pub r_rel: Vec3,
    /// Body's velocity relative to the primary, `v_body − v_primary`.
    /// Required by anomaly recovery in [`elements_from_invariants`].
    pub v_rel: Vec3,
}

/// Computes the linear primitives (ε, h_vec, e_vec) for body `idx`
/// relative to `primary_idx` under gravitational parameter `g_factor`.
///
/// All arithmetic uses fixed `(x, y, z)` reduction order; re-associating
/// any inner sum into a Vec3-level expression would shift ULPs and is
/// observable on the energy / momentum / Mercury 1PN gates. Specifically:
///
/// * `r²` and `v²` use [`Vec3::length_squared`] (`((x² + y²) + z²)`
///   left-to-right via the `dot(self, self)` form).
/// * `h_vec = r × v` uses the canonical [`Vec3::cross`] component formula.
/// * `e_vec` is assembled element-wise from `(v×h)/μ − r̂` with three
///   separate `(v_cross_h.{x,y,z} / μ) − r̂.{x,y,z}` lines, never
///   consolidated into something like `(v.cross(h) − r̂·μ) / μ`.
///
/// `r̂ = r_vec / r` is computed by bare division after the `r < 1e-15`
/// gate guarantees a safe denominator. [`Vec3::try_normalize`] is
/// deliberately **not** used here: it applies a different threshold
/// (`f64::MIN_POSITIVE` ≈ `2 × 10⁻³⁰⁸`) and changes the rounding of
/// the result.
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

    let r_vec = Vec3::new(b.pos_x - p.pos_x, b.pos_y - p.pos_y, b.pos_z - p.pos_z);
    let v_vec = Vec3::new(b.vel_x - p.vel_x, b.vel_y - p.vel_y, b.vel_z - p.vel_z);

    let r = r_vec.length();
    let gm = g_factor * (b.mass + p.mass);

    if r < 1e-15 || gm < 1e-30 {
        return None;
    }

    let v2 = v_vec.length_squared();
    let energy = 0.5 * v2 - gm / r;

    let h_vec = r_vec.cross(v_vec);
    let v_cross_h = v_vec.cross(h_vec);

    // r̂ = r_vec / r — safe after the `r < 1e-15` gate above. Element-wise
    // assembly of `e_vec` preserves the per-axis reduction order.
    let inv_r = 1.0 / r;
    let r_hat = Vec3::new(r_vec.x * inv_r, r_vec.y * inv_r, r_vec.z * inv_r);

    let inv_gm = 1.0 / gm;
    let e_vec = Vec3::new(
        v_cross_h.x * inv_gm - r_hat.x,
        v_cross_h.y * inv_gm - r_hat.y,
        v_cross_h.z * inv_gm - r_hat.z,
    );

    Some(OrbitInvariants { energy, h_vec, e_vec, r_rel: r_vec, v_rel: v_vec })
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
    let e = inv.e_vec.length();
    let h_mag = inv.h_vec.length();

    // Inclination from `h_vec.z / |h_vec|`. Clamped because
    // |h_vec.z| ≤ |h_vec| holds analytically but a +1 ULP overshoot
    // out of `acos`'s domain would NaN the result.
    let inclination = if h_mag > 0.0 { (inv.h_vec.z / h_mag).clamp(-1.0, 1.0).acos() } else { 0.0 };

    // Node line `n = ẑ × h_vec = (-h_vec.y, h_vec.x, 0)`. `|n|` measures
    // the planar degeneracy: vanishes for both prograde-planar
    // (`h_vec.z > 0, |h_vec.x| = |h_vec.y| = 0`) and retrograde-planar
    // (`h_vec.z < 0, |h_vec.x| = |h_vec.y| = 0`) orbits.
    let n_x = -inv.h_vec.y;
    let n_y = inv.h_vec.x;
    let n_mag = (n_x * n_x + n_y * n_y).sqrt();

    // `atan2` results are already in `(-π, π]` — no `wrap_pi` is
    // applied here. Reserving `wrap_pi` for sites where two angles
    // are summed or differenced (e.g. `θ_body − ν` in the anchored
    // path) keeps the planar baseline bit-exact: `wrap_pi` of an
    // in-range angle is the identity in exact arithmetic but
    // introduces ULP-level rounding through the
    // `rem_euclid(2π) − π` round-trip, observable on the
    // `baseline_newtonian_kepler_is_closed` drift gate.
    let (lon_ascending_node, omega) = if n_mag < N_DEGENERATE_EPS {
        // Planar fallback: `Ω` is undefined → 0; `ω = atan2(e_y, e_x)`.
        let omega = if e > E_CIRCULAR_EPS { inv.e_vec.y.atan2(inv.e_vec.x) } else { 0.0 };
        (0.0, omega)
    } else {
        // Inclined branch. `Ω = atan2(n_y, n_x)` from the node line.
        let lon_asc = n_y.atan2(n_x);

        let omega = if e > E_CIRCULAR_EPS {
            // ω = atan2(e · (ĥ × n̂), e · n̂)
            //
            // Both atan2 args carry an identical scaling factor of
            // `|h||n|` (see header doc): the perpendicular axis
            // `ĥ × n̂` has magnitude `|h||n|` because ĥ ⊥ n̂; dotting
            // `e_vec` against the un-normalised `h_vec × n_vec` and
            // against `n_vec` and dividing the first by `|h|`
            // recovers two quantities scaled identically by `|n|`.
            // atan2 is invariant under common positive scaling, so
            // the explicit normalisation by `|n|` is skipped.
            let n_vec = Vec3::new(n_x, n_y, 0.0);
            let h_cross_n = inv.h_vec.cross(n_vec);
            let sin_arg = inv.e_vec.dot(h_cross_n) / h_mag;
            let cos_arg = inv.e_vec.dot(n_vec);
            sin_arg.atan2(cos_arg)
        } else {
            0.0
        };

        (lon_asc, omega)
    };

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

    let (true_anomaly, eccentric_anomaly, mean_anomaly) = anomalies(orbit_type, e, omega, inv);

    OrbitalElements {
        primary_idx,
        a,
        e,
        period,
        h_vec: inv.h_vec,
        energy: inv.energy,
        omega,
        inclination,
        lon_ascending_node,
        true_anomaly,
        eccentric_anomaly,
        mean_anomaly,
        orbit_type,
    }
}

/// Recover the true / eccentric / mean anomalies from the orbit invariants.
///
/// For the elliptical and hyperbolic cases the routine returns the canonical
/// triple `(ν, E, M)` (or `(ν, H, M_h)` for hyperbolic). Parabolic orbits
/// return `NaN` for all three — the parabolic anomaly system is governed
/// by Barker's equation, which is intentionally out of scope while no caller
/// renders parabolic flybys.
fn anomalies(orbit_type: OrbitType, e: f64, omega: f64, inv: &OrbitInvariants) -> (f64, f64, f64) {
    let r = inv.r_rel.length();
    if r < 1e-15 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }

    let nu = if e > E_CIRCULAR_EPS {
        // ν from focus-conic: cos ν = (e_vec · r) / (e · r); sign from r·v
        // (positive when receding from primary, i.e. moving from peri- to apo-centre).
        let cos_nu = (inv.e_vec.dot(inv.r_rel) / (e * r)).clamp(-1.0, 1.0);
        let raw = cos_nu.acos();
        if inv.r_rel.dot(inv.v_rel) >= 0.0 { raw } else { -raw }
    } else {
        // Circular: ν is undefined. Convention: return the argument of
        // latitude `u = atan2(r·m̂, r·n̂)` (with `n̂` = ascending node
        // direction, `m̂ = ĥ × n̂`), shifted by ω. Since ω = 0 for circular
        // orbits by convention, this returns the body's angle from the
        // ascending node — preserving position information for visualisation.
        let u = argument_of_latitude(inv.r_rel, inv.h_vec);
        wrap_pi(u - omega)
    };

    match orbit_type {
        OrbitType::Elliptical => {
            // E from ν via half-angle identity, evaluated as `atan2` to
            // avoid the tan(π/2) branch when ν approaches ±π.
            let half_nu = 0.5 * nu;
            let (s_half, c_half) = half_nu.sin_cos();
            let big_e = 2.0 * ((1.0 - e).sqrt() * s_half).atan2((1.0 + e).sqrt() * c_half);
            // M = E − e sin E. Not wrapped to `[−π, π]` because M is
            // semantically a temporal variable (`M(t) = M₀ + n t`); a
            // principal-value collapse would corrupt consumers tracking
            // dM/dt across frames. Snapshot recovery from instantaneous
            // (r, v) cannot supply absolute temporal context regardless,
            // so the value lands naturally in `[−π − e, π + e]` and the
            // `2π` jump at orbit closure must be unwrapped by the
            // consumer with a continuity-tracking layer.
            let m = big_e - e * big_e.sin();
            (nu, big_e, m)
        },
        OrbitType::Hyperbolic => {
            // H from ν via half-angle: tanh(H/2) = √((e−1)/(e+1)) tan(ν/2).
            // The argument can grow beyond ±1 if ν is past the asymptote
            // ν_∞ = ±acos(−1/e); clamp to keep atanh defined and let the
            // saturated value mark the body sitting on the asymptote.
            let half_nu = 0.5 * nu;
            let arg = ((e - 1.0) / (e + 1.0)).sqrt() * half_nu.tan();
            let half_h = arg.clamp(-1.0 + 1e-15, 1.0 - 1e-15).atanh();
            let big_h = 2.0 * half_h;
            let m_h = e * big_h.sinh() - big_h;
            (nu, big_h, m_h)
        },
        OrbitType::Parabolic => (f64::NAN, f64::NAN, f64::NAN),
    }
}

/// Argument of latitude — angle from the ascending node to the body's
/// current position, measured in the orbital plane in the direction of
/// motion. Reduces to `atan2(r.y, r.x)` for planar prograde orbits.
fn argument_of_latitude(r_rel: Vec3, h_vec: Vec3) -> f64 {
    let h_mag = h_vec.length();
    if h_mag < 1e-15 {
        return 0.0;
    }
    let h_hat = h_vec / h_mag;

    // Node line `n = ẑ × h`. Vanishes for planar prograde or retrograde
    // orbits — fall back to `+x̂` so a planar circular orbit reports the
    // body's argument as `atan2(r.y, r.x)`.
    let n_x = -h_vec.y;
    let n_y = h_vec.x;
    let n_mag = (n_x * n_x + n_y * n_y).sqrt();
    let n_hat = if n_mag > N_DEGENERATE_EPS {
        Vec3::new(n_x / n_mag, n_y / n_mag, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };

    let m_hat = h_hat.cross(n_hat);
    let along_node = r_rel.dot(n_hat);
    let along_perp = r_rel.dot(m_hat);
    along_perp.atan2(along_node)
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

    let e = inv.e_vec.length();

    // Near-circular: anchor's (p/r − 1)/e branch amplifies smoothed
    // residual into visible ω jitter. Use the smoothed eccentricity
    // vector's direction directly — stable, with at most ~e·a body
    // offset from the displayed orbit (sub-pixel for e < this threshold).
    if e < ANCHOR_MIN_E {
        return elements_from_invariants(inv, primary_idx, gm);
    }

    let a = -gm / (2.0 * inv.energy);
    let period = TAU * (a * a * a / gm).sqrt();

    // 3D state vectors in the primary's frame.
    let r_vec = Vec3::new(
        body.pos_x - primary.pos_x,
        body.pos_y - primary.pos_y,
        body.pos_z - primary.pos_z,
    );
    let v_vec = Vec3::new(
        body.vel_x - primary.vel_x,
        body.vel_y - primary.vel_y,
        body.vel_z - primary.vel_z,
    );

    let h_mag = inv.h_vec.length();

    // Out-of-plane projection breaks down without a defined orbital
    // normal; fall back to the unanchored elements.
    if h_mag < 1e-15 {
        return elements_from_invariants(inv, primary_idx, gm);
    }
    let h_hat = inv.h_vec / h_mag;

    // Project the body's position onto the smoothed orbital plane to
    // get an in-plane radius `r_eff`. Any out-of-plane component is
    // discarded — the displayed ellipse always lives in that plane,
    // so the body's projection is the closest point that can possibly
    // sit on it. For a body whose true plane matches the smoothed
    // plane, `r_eff ≈ |r_vec|` and the anchor reduces to the planar
    // case below.
    let r_in_plane = r_vec - h_hat * r_vec.dot(h_hat);
    let r_eff = r_in_plane.length();
    if r_eff < 1e-15 {
        return elements_from_invariants(inv, primary_idx, gm);
    }

    // Reference frame inside the orbital plane. `node_hat` is the
    // standard ascending-node direction (ẑ × ĥ). When the orbit is
    // planar (h_hat ≈ ±ẑ) this collapses; pick the world x-axis as a
    // stable fallback so atan2 keeps a meaningful zero.
    let node_raw = Vec3::new(-h_hat.y, h_hat.x, 0.0);
    let node_mag = node_raw.length();
    let node_hat = if node_mag > 1e-9 {
        node_raw / node_mag
    } else {
        let proj = Vec3::new(1.0, 0.0, 0.0);
        let proj = proj - h_hat * proj.dot(h_hat);
        let pl = proj.length();
        if pl > 1e-15 {
            proj / pl
        } else {
            return elements_from_invariants(inv, primary_idx, gm);
        }
    };
    let perp_node_hat = h_hat.cross(node_hat);

    // Solve focus-conic for cos ν, clamping to absorb body distances
    // slightly outside the smoothed ellipse's apsidal range.
    let p = a * (1.0 - e * e);
    let cos_nu = ((p / r_eff - 1.0) / e).clamp(-1.0, 1.0);

    // sign(sin ν) = sign((r⃗ × v⃗)·ĥ · (r⃗·v⃗)). The triple product gives
    // the orbit chirality relative to the smoothed plane; (r⃗·v⃗) gives
    // outbound/inbound. Their product disambiguates the ±ν branch.
    let h_inst_signed = r_vec.cross(v_vec).dot(h_hat);
    let r_dot_v = r_vec.dot(v_vec);
    let s = (h_inst_signed * r_dot_v).signum();
    let sin_nu_sign = if s == 0.0 { 1.0 } else { s };
    let nu = sin_nu_sign * cos_nu.acos();

    // Body's azimuthal angle inside the smoothed orbital plane,
    // measured from the line of nodes. ω anchors the perifocal
    // x-axis so the body lands at parametric angle ν on the ellipse.
    let theta_body = r_in_plane.dot(perp_node_hat).atan2(r_in_plane.dot(node_hat));
    let omega = wrap_pi(theta_body - nu);

    let inclination = h_hat.z.clamp(-1.0, 1.0).acos();
    let lon_ascending_node = if node_mag > 1e-9 { node_hat.y.atan2(node_hat.x) } else { 0.0 };

    let half_nu = 0.5 * nu;
    let (s_half, c_half) = half_nu.sin_cos();
    let big_e = 2.0 * ((1.0 - e).sqrt() * s_half).atan2((1.0 + e).sqrt() * c_half);
    let mean_anomaly = big_e - e * big_e.sin();

    OrbitalElements {
        primary_idx,
        a,
        e,
        period,
        h_vec: inv.h_vec,
        energy: inv.energy,
        omega,
        inclination,
        lon_ascending_node,
        true_anomaly: nu,
        eccentric_anomaly: big_e,
        mean_anomaly,
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
        assert!(rel(a.h_vec.z, b.h_vec.z) < tol, "{label}: h_z {} vs {}", a.h_vec.z, b.h_vec.z);
        assert!(
            rel(a.omega, b.omega) < tol || (a.e < 1e-6 && b.e < 1e-6),
            "{label}: ω {} vs {}",
            a.omega,
            b.omega
        );
        // Anomalies — only compare when both are defined and the orbit is
        // non-circular. Circular orbits have ill-conditioned anomalies and
        // ω jitter dominates the residual.
        if a.e >= 1e-6 && b.e >= 1e-6 {
            for (name, x, y) in [
                ("ν", a.true_anomaly, b.true_anomaly),
                ("E", a.eccentric_anomaly, b.eccentric_anomaly),
                ("M", a.mean_anomaly, b.mean_anomaly),
            ] {
                if x.is_finite() && y.is_finite() {
                    assert!(rel(x, y) < tol, "{label}: {name} {} vs {}", x, y);
                }
            }
        }
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
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, v_peri, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].pos_x - bodies[0].pos_x;
        let ry = bodies[1].pos_y - bodies[0].pos_y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        assert!((r - r_orbit).abs() / r < 1e-12, "body off orbit: r={r}, r_orbit={r_orbit}");
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
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, -v_peri, 1e-10)];
        let inv = compute_invariants(&bodies, 1, 0, G).unwrap();
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].pos_x - bodies[0].pos_x;
        let ry = bodies[1].pos_y - bodies[0].pos_y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        assert!(
            (r - r_orbit).abs() / r < 1e-12,
            "retrograde body off orbit: r={r}, r_orbit={r_orbit}"
        );
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
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, v_peri, 1e-10)];
        // Real invariants of the body's orbit.
        let inv_true = compute_invariants(&bodies, 1, 0, G).unwrap();
        // Simulate smoothed invariants drifted by 1% in energy and h.
        let inv_smooth = OrbitInvariants {
            energy: inv_true.energy * 1.01,
            h_vec: inv_true.h_vec * 0.99,
            e_vec: Vec3::new(inv_true.e_vec.x * 1.01, inv_true.e_vec.y * 0.99, inv_true.e_vec.z),
            r_rel: inv_true.r_rel,
            v_rel: inv_true.v_rel,
        };
        let gm = G * (bodies[1].mass + bodies[0].mass);
        let el = elements_anchored_to_body(&inv_smooth, 0, gm, &bodies[1], &bodies[0]);

        let rx = bodies[1].pos_x - bodies[0].pos_x;
        let ry = bodies[1].pos_y - bodies[0].pos_y;
        let r = (rx * rx + ry * ry).sqrt();
        let nu = ry.atan2(rx) - el.omega;
        let r_orbit = el.a * (1.0 - el.e * el.e) / (1.0 + el.e * nu.cos());
        // Body should still lie on the (drifted-shape) ellipse to FP
        // precision, because anchoring depends only on r and the
        // smoothed (a, e) — it does not require the orbit to be the
        // body's *true* orbit.
        assert!(
            (r - r_orbit).abs() / r < 1e-12,
            "drifted-invariant anchored orbit: r={r}, r_orbit={r_orbit}"
        );
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
            inv.e_vec.length() < ANCHOR_MIN_E,
            "test setup: e must be below threshold to exercise the bypass"
        );

        let gm_pair = G * (bodies[1].mass + primary.mass);
        let el_at_peri = elements_anchored_to_body(&inv, 0, gm_pair, &bodies[1], &primary);

        // Move the body to a different point along its orbit (any
        // arbitrary location); invariants stay the same since they are
        // the *smoothed* state passed in, not recomputed from the body.
        let displaced = body(0.0, r * 1.1, -v, 0.0, 1e-10);
        let el_displaced = elements_anchored_to_body(&inv, 0, gm_pair, &displaced, &primary);

        assert_eq!(
            el_at_peri.omega, el_displaced.omega,
            "low-e bypass must produce constant ω across body positions"
        );
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

        assert_ne!(
            el_a.omega, el_b.omega,
            "above threshold, anchor must rotate ω to keep body on orbit"
        );
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
        assert!(el.h_vec.z > 0.0, "h_z = {}, deve ser positivo (CCW)", el.h_vec.z);
    }

    #[test]
    fn cw_orbit_has_negative_angular_momentum() {
        let (p, mut s) = circular_orbit(10.0, 1e6);
        s.vel_y = -s.vel_y; // inverte para CW
        let el = elements(p, s);
        assert!(el.h_vec.z < 0.0, "h_z = {}, deve ser negativo (CW)", el.h_vec.z);
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
        let err = (el.h_vec.z - h_expected).abs() / h_expected;
        assert!(err < 1e-6, "h_z = {}, esperado {h_expected}", el.h_vec.z);
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

    // ── 5. Property tests for the 3D element pipeline ────────────────────────
    //
    // Two empirical contracts that protect against a future stylistic
    // refactor that quietly reorders a Vec3 expression or returns the
    // wrong quadrant from atan2:
    //
    //   (a) Planar input (`z = vz = 0`) produces the explicit planar
    //       formulas — `Ω = 0`, `ω = atan2(e_y, e_x)` — to f64 precision.
    //   (b) Inclined input recovers `inclination`, `Ω`, `ω` correctly,
    //       with angles in `[-π, π]` and stable across small
    //       perturbations.

    /// Singularity contract: `e < E_CIRCULAR_EPS` ⇒ `ω = 0` exactly.
    ///
    /// A circular orbit has no defined argument of periapsis. The
    /// convention pinned by `elements_from_invariants` is to return 0
    /// rather than the noise-amplified `atan2` of a near-zero
    /// eccentricity vector. This is the contract the smoother and the
    /// orbit-overlay UI rely on to avoid frame-to-frame ω jitter on
    /// near-circular outer planets.
    #[test]
    fn omega_is_zero_for_circular_orbit() {
        let (primary, satellite) = circular_orbit(10.0, 1e6);
        let el = elements(primary, satellite);
        assert!(el.e < E_CIRCULAR_EPS, "test setup: e must be below the circular threshold");
        assert_eq!(el.omega, 0.0, "ω must be exactly zero for a circular orbit");
    }

    /// Singularity contract: planar orbit (`|n| < N_DEGENERATE_EPS`)
    /// ⇒ `Ω = 0`, `inclination = 0` (prograde) or `π` (retrograde).
    ///
    /// For a body confined to `z = vz = 0`, the angular momentum is
    /// purely along ±ẑ and the node line `n = ẑ × h_vec` is exactly
    /// zero. Returning `Ω = 0` by convention prevents the inclined
    /// branch's `atan2` from operating on a numerically degenerate
    /// vector.
    #[test]
    fn lon_ascending_node_is_zero_for_planar_orbit() {
        // Prograde planar orbit (inclination ≈ 0).
        let (primary, satellite) = circular_orbit(5.0, 1e6);
        let el = elements(primary, satellite);
        assert_eq!(el.lon_ascending_node, 0.0, "Ω must be 0 for prograde planar input");
        assert!(
            el.inclination.abs() < 1e-12,
            "inclination must be ≈ 0 for prograde planar input, got {}",
            el.inclination
        );

        // Retrograde planar orbit (inclination ≈ π). Same input but
        // velocity flipped on y-axis → angular momentum along −ẑ.
        let v_c = (G * 1e6 / 5.0).sqrt();
        let primary_retro = body(0.0, 0.0, 0.0, 0.0, 1e6);
        let satellite_retro = body(5.0, 0.0, 0.0, -v_c, 1e-10);
        let el_retro = elements(primary_retro, satellite_retro);
        assert_eq!(el_retro.lon_ascending_node, 0.0, "Ω must be 0 for retrograde planar input");
        assert!(
            (el_retro.inclination - PI).abs() < 1e-12,
            "inclination must be ≈ π for retrograde planar input, got {}",
            el_retro.inclination
        );
    }

    /// Bit-exact contract: the planar branch in
    /// `elements_from_invariants` produces `ω = atan2(e_vec.y, e_vec.x)`.
    ///
    /// For `z = vz = 0` input, `h_vec.x = h_vec.y = 0` exactly and the
    /// node vector `n` is the additive identity, so the branch is
    /// taken structurally (not by threshold luck). Equality, not
    /// tolerance — protects against any reduction-order shift.
    #[test]
    fn planar_input_omega_matches_atan2_e_y_e_x() {
        // Eccentric planar orbit: e ≈ 0.3, periapsis displaced from
        // the +x axis so `ω` is non-trivial.
        let primary = body(0.0, 0.0, 0.0, 0.0, 1e6);
        let satellite = body(7.0, 4.0, -0.6, 0.4, 1e-10);

        let inv = compute_invariants(&[primary, satellite], 1, 0, G).unwrap();
        let gm = G * (primary.mass + satellite.mass);
        let el = elements_from_invariants(&inv, 0, gm);

        let omega_ref = inv.e_vec.y.atan2(inv.e_vec.x);
        assert_eq!(el.omega, omega_ref, "planar fallback ω must equal atan2(e_y, e_x) bit-for-bit");
    }

    /// Inclined elliptical orbit recovers the textbook inclination.
    ///
    /// Configures a 30°-inclined orbit by rotating a planar circular
    /// orbit's velocity into the `xz`-plane (rotation around `ŷ`)
    /// and checks that the recovered `inclination` matches `π/6`
    /// to f64 precision. This is the simplest test that fails under
    /// any sign or component error in the cross-product chain.
    #[test]
    fn inclined_orbit_recovers_inclination() {
        let i_target = PI / 6.0; // 30°
        let r = 10.0;
        let m = 1e6;
        let v_c = (G * m / r).sqrt();
        // Rotate the planar velocity (0, v_c, 0) by `i` around ŷ:
        // v_rot = (v_c · sin i, v_c · cos i · 0 + …) — but we placed
        // the body on +x axis, so we rotate the velocity vector
        // (0, v_c, 0) → (0, v_c · cos i, v_c · sin i).
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = Body::rocky(1e-10).at_3d(r, 0.0, 0.0).with_velocity_3d(
            0.0,
            v_c * i_target.cos(),
            v_c * i_target.sin(),
        );

        let bodies = vec![primary, satellite];
        let el = compute_elements(&bodies, 1, 0, G).expect("inclined orbit must yield elements");

        assert!(
            (el.inclination - i_target).abs() < 1e-12,
            "inclination: got {}, expected {} ({} arcsec error)",
            el.inclination,
            i_target,
            (el.inclination - i_target).abs() * 206_265.0,
        );
        // Sanity: angular range pinned.
        assert!(
            el.lon_ascending_node >= -PI && el.lon_ascending_node <= PI,
            "Ω must lie in [-π, π], got {}",
            el.lon_ascending_node,
        );
        assert!(el.omega >= -PI && el.omega <= PI, "ω must lie in [-π, π], got {}", el.omega,);
    }

    /// Hyperbolic flyby in 3D: continuity of angular elements under a
    /// small kinematic perturbation.
    ///
    /// Builds a hyperbolic orbit (e > 1) at 30° inclination and then
    /// perturbs the body's velocity by 1% in each component. The
    /// angular elements (`Ω`, `ω`) must respond by amounts of the
    /// same order — never flip by `±π`, which would indicate a
    /// quadrant misclassification in the `atan2` projections. This
    /// is the regression test for the trap that "ω = atan2(e_vec.y,
    /// e_vec.x)" would walk into for inclined orbits.
    #[test]
    fn hyperbolic_inclined_continuity_under_perturbation() {
        let i_target = PI / 6.0;
        let r = 5.0;
        let m = 1e6;
        // Hyperbolic excess: v > v_escape.
        let v_esc = (2.0 * G * m / r).sqrt();
        let v_hyp = 1.5 * v_esc;
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite = Body::rocky(1e-10).at_3d(r, 0.0, 0.0).with_velocity_3d(
            0.0,
            v_hyp * i_target.cos(),
            v_hyp * i_target.sin(),
        );

        let el_a = compute_elements(&[primary, satellite], 1, 0, G).unwrap();

        // 1% perturbation on each velocity component.
        let perturbed = Body::rocky(1e-10).at_3d(r, 0.0, 0.0).with_velocity_3d(
            0.01 * v_hyp,
            v_hyp * i_target.cos() * 1.01,
            v_hyp * i_target.sin() * 0.99,
        );
        let el_b = compute_elements(&[primary, perturbed], 1, 0, G).unwrap();

        assert!(matches!(el_a.orbit_type, OrbitType::Hyperbolic));
        assert!(matches!(el_b.orbit_type, OrbitType::Hyperbolic));

        // Continuity: changes are of the same order as the
        // perturbation (~1%), not π.
        assert!(
            (el_a.inclination - el_b.inclination).abs() < 0.05,
            "inclination jumped by {} rad under 1% perturbation — likely a quadrant flip",
            (el_a.inclination - el_b.inclination).abs(),
        );
        assert!(
            (el_a.lon_ascending_node - el_b.lon_ascending_node).abs() < 0.5,
            "Ω jumped by {} rad under 1% perturbation",
            (el_a.lon_ascending_node - el_b.lon_ascending_node).abs(),
        );
        assert!(
            (el_a.omega - el_b.omega).abs() < 0.5,
            "ω jumped by {} rad under 1% perturbation — sign of a quadrant misclassification",
            (el_a.omega - el_b.omega).abs(),
        );
    }

    /// Invariance contract: for planar input, `|h_vec|` exactly equals
    /// `|h_z|`.
    ///
    /// Asserts the precondition `h_vec.x == 0 && h_vec.y == 0` first
    /// so the contract cannot be re-used out of context — anyone
    /// reading the test sees that the equality only holds when the
    /// orbit is in the `xy`-plane.
    #[test]
    fn h_vec_magnitude_equals_h_z_when_planar() {
        let (primary, satellite) = circular_orbit(8.0, 1e6);
        let inv = compute_invariants(&[primary, satellite], 1, 0, G).unwrap();

        // Precondition: orbit is in the `xy`-plane, so `h_vec` is along ±ẑ.
        assert_eq!(inv.h_vec.x, 0.0, "test setup: planar input must have h_vec.x == 0");
        assert_eq!(inv.h_vec.y, 0.0, "test setup: planar input must have h_vec.y == 0");

        // Under that precondition, the magnitude reduces to |h_z|.
        assert_eq!(
            inv.h_vec.length(),
            inv.h_vec.z.abs(),
            "|h_vec| must equal |h_z| when h_vec.x = h_vec.y = 0",
        );
    }

    // ── 6. 3D validation portfolio (algebraic) ────────────────────────────────
    //
    // The planar tests above cover regression: they confirm that `z = vz = 0`
    // input still produces the same numbers it did before the 3D port. The
    // tests below cover the complementary direction: that `z ≠ 0` /
    // `vz ≠ 0` input is *handled correctly*, not just propagated as zeros.
    //
    // These three tests are deliberately structural rather than statistical:
    //   - cross-axis isolation: motion in z must NOT contaminate x / y
    //   - h_vec direction: r × v orientation must match the right-hand rule
    //   - (i, Ω, ω) round-trip: orbital elements must invert correctly,
    //     including at quadrant boundaries where atan2 is fragile
    //
    // The (i, Ω, ω) test deliberately constructs orbits via a test-local
    // Euler rotation written from scratch — NOT through `sample_orbit` —
    // so the reference comes from an independent implementation rather
    // than from the same code path the test exercises.

    /// Apply the standard 3-1-3 Euler rotation `R_z(Ω) R_x(i) R_z(ω)` to a
    /// perifocal-frame vector, mapping it to the inertial frame. Written
    /// from scratch so the reference is independent of any rotation in
    /// `apsis::physics::orbital`. Convention: right-handed, i ∈ [0, π],
    /// Ω, ω ∈ [-π, π].
    fn perifocal_to_inertial(perifocal: Vec3, omega: f64, i: f64, lon_asc: f64) -> Vec3 {
        let (sw, cw) = omega.sin_cos();
        let (si, ci) = i.sin_cos();
        let (so, co) = lon_asc.sin_cos();

        // First R_z(ω): rotates within the orbital plane around ẑ_pf.
        let x1 = perifocal.x * cw - perifocal.y * sw;
        let y1 = perifocal.x * sw + perifocal.y * cw;
        let z1 = perifocal.z;

        // Then R_x(i): tilts the plane around x̂.
        let x2 = x1;
        let y2 = y1 * ci - z1 * si;
        let z2 = y1 * si + z1 * ci;

        // Finally R_z(Ω): rotates the line of nodes into place around ẑ.
        let x3 = x2 * co - y2 * so;
        let y3 = x2 * so + y2 * co;
        let z3 = z2;

        Vec3::new(x3, y3, z3)
    }

    /// Cross-axis isolation: a body moving purely along ẑ must keep `x` and
    /// `y` exactly at their initial values across the entire integration.
    ///
    /// This is the cheapest test that catches accidental coupling between
    /// the z-axis acceleration and the in-plane components — a `dx` term
    /// receiving a contribution from `dz`, an `acc.x` accumulator
    /// receiving the wrong index, or a Vec3 swizzle bug. Two equal-mass
    /// bodies on the z-axis fall toward each other under gravity that
    /// must be purely along ẑ; any non-zero `x` or `y` after integration
    /// is a coupling leak.
    #[test]
    fn pure_z_motion_preserves_xy_components() {
        use crate::core::system::System;
        use crate::physics::integrator::IntegratorKind;
        use crate::units::UnitSystem;

        let m = 1.0;
        let z0 = 1.0;
        let a = Body::rocky(m).at_3d(0.0, 0.0, z0).with_velocity_3d(0.0, 0.0, -0.05);
        let b = Body::rocky(m).at_3d(0.0, 0.0, -z0).with_velocity_3d(0.0, 0.0, 0.05);

        let mut sys = System::new(vec![a, b], UnitSystem::canonical())
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(1e-3);

        for _ in 0..200 {
            sys.step();
        }

        for (k, body) in sys.bodies().iter().enumerate() {
            // The integrator must not introduce any in-plane motion. The
            // bound here is the f64 round-off floor — anything above 1e-14
            // indicates a structural coupling, not a tolerance failure.
            assert!(
                body.pos_x.abs() < 1e-14,
                "body {k}: x drifted to {} from a pure-z initial condition",
                body.pos_x,
            );
            assert!(
                body.pos_y.abs() < 1e-14,
                "body {k}: y drifted to {} from a pure-z initial condition",
                body.pos_y,
            );
            assert!(
                body.vel_x.abs() < 1e-14,
                "body {k}: vx drifted to {} from a pure-z initial condition",
                body.vel_x,
            );
            assert!(
                body.vel_y.abs() < 1e-14,
                "body {k}: vy drifted to {} from a pure-z initial condition",
                body.vel_y,
            );
        }
    }

    /// `h_vec = r × v` must point in the direction prescribed by the
    /// right-hand rule for a circular orbit tilted out of the `xy`-plane.
    ///
    /// Setup: body at `(R, 0, 0)` with velocity rotated 60° around `x̂`
    /// from the planar `(0, v_c, 0)`. The orbital plane is then tilted
    /// 60° from the `xy`-plane and `h_vec` should sit perpendicular to
    /// the velocity direction within that plane:
    ///
    ///   `h_vec / |h_vec| = (0, -sin(i), cos(i))`
    ///
    /// A bug that swaps a cross-product component (e.g.,
    /// `(r × v).x = ry · vz` written as `rx · vz`) flips this direction
    /// detectably, even though the magnitude `|h_vec|` would be unchanged.
    #[test]
    fn h_vec_direction_matches_right_hand_rule_for_inclined_circular() {
        let m = 1e6;
        let r = 4.0;
        let v_c = (G * m / r).sqrt();
        let i = std::f64::consts::FRAC_PI_3; // 60°
        let (sin_i, cos_i) = i.sin_cos();

        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        let satellite =
            Body::rocky(1e-10).at_3d(r, 0.0, 0.0).with_velocity_3d(0.0, v_c * cos_i, v_c * sin_i);

        let inv = compute_invariants(&[primary, satellite], 1, 0, G).unwrap();
        let h_mag = inv.h_vec.length();

        // Analytic prediction: r × v with r = (R, 0, 0) and
        // v = (0, v_c·cos i, v_c·sin i) is (0, -R·v_c·sin i, R·v_c·cos i).
        let expected = Vec3::new(0.0, -sin_i, cos_i);
        let unit = inv.h_vec / h_mag;

        // Component-wise comparison rather than dot-product — the test
        // localises a swizzle bug to the offending axis.
        assert!(
            (unit.x - expected.x).abs() < 1e-12,
            "h.pos_x direction off: {} vs {}",
            unit.x,
            expected.x
        );
        assert!(
            (unit.y - expected.y).abs() < 1e-12,
            "h.pos_y direction off: {} vs {}",
            unit.y,
            expected.y
        );
        assert!(
            (unit.z - expected.z).abs() < 1e-12,
            "h.pos_z direction off: {} vs {}",
            unit.z,
            expected.z
        );
    }

    /// `compute_invariants → elements_from_invariants` round-trip: build
    /// an orbit with prescribed `(i, Ω, ω)` via the independent Euler
    /// rotation, recover the elements, and assert the recovered angles
    /// match the input.
    ///
    /// Two sub-cases:
    ///
    ///   - **mid-quadrant** `(i, Ω, ω) = (30°, 45°, 60°)` — sanity that
    ///     the algebra works in a generic regime.
    ///   - **quadrant-boundary** `(i, Ω, ω) = (60°, π − ε, −π + ε)` —
    ///     `Ω` and `ω` sit on opposite sides of the wrap; recovery must
    ///     stay close to the input value (no spurious flip by `±π`,
    ///     no atan2 quadrant misclassification).
    ///
    /// The quadrant-boundary case is what catches sign errors in the
    /// `n` vector or in `e_vec · (h × n)` that the mid-quadrant case
    /// would silently absorb.
    #[test]
    fn orbital_elements_round_trip_under_3d_rotation() {
        // Reference Keplerian orbit in the perifocal frame: periapsis on
        // +x̂_pf, body at the periapsis (true anomaly 0). For a circular
        // orbit `e ≈ 0` makes `ω` undefined, so we use a moderately
        // eccentric configuration where `ω` is meaningful.
        let m = 1.0;
        let a = 1.0;
        let e = 0.3;
        let r_peri = a * (1.0 - e);
        let v_peri = (G * m * (1.0 + e) / (a * (1.0 - e))).sqrt();

        let cases = [
            ("mid-quadrant", PI / 6.0, PI / 4.0, PI / 3.0),
            ("quadrant-boundary", PI / 3.0, PI - 0.01, -PI + 0.01),
        ];

        for (label, i_target, lon_asc_target, omega_target) in cases {
            // Perifocal-frame state at periapsis: position on +x̂, velocity on +ŷ.
            let r_pf = Vec3::new(r_peri, 0.0, 0.0);
            let v_pf = Vec3::new(0.0, v_peri, 0.0);

            // Map to inertial via R_z(Ω) R_x(i) R_z(ω) — independent of
            // anything in `physics::orbital`.
            let r_inertial = perifocal_to_inertial(r_pf, omega_target, i_target, lon_asc_target);
            let v_inertial = perifocal_to_inertial(v_pf, omega_target, i_target, lon_asc_target);

            let primary = body(0.0, 0.0, 0.0, 0.0, m);
            let satellite = Body::rocky(1e-10)
                .at_3d(r_inertial.x, r_inertial.y, r_inertial.z)
                .with_velocity_3d(v_inertial.x, v_inertial.y, v_inertial.z);

            let el = compute_elements(&[primary, satellite], 1, 0, G)
                .expect("elements must be computable");

            // Tolerance: 1e-9 rad ≈ 0.2 milli-arcsec at 1 AU — captures
            // rotation-algebra correctness without being defeated by the
            // trig chain ULPs (sin/cos × cross × dot accumulates ~1 ULP).
            // A real bug surfaces as residuals ≥ 1e-3 (atan2 quadrant flip
            // ⇒ ~π) or ≥ 1 (sign error).
            assert!(
                (el.inclination - i_target).abs() < 1e-9,
                "{label}: i recovered as {} rad, expected {i_target} rad",
                el.inclination,
            );

            // Ω and ω are in `[-π, π]`; near the boundary, a recovered
            // value `+π − ε` and an input `-π + ε` are almost-identical
            // points but their plain difference is ~2π. Compare via the
            // wrapped difference so a true π-flip surfaces as a > 1
            // residual while the boundary equivalence reads as ~0.
            let omega_diff = ((el.omega - omega_target + PI).rem_euclid(TAU) - PI).abs();
            assert!(
                omega_diff < 1e-9,
                "{label}: ω recovered as {} rad, expected {omega_target} rad (wrapped diff {})",
                el.omega,
                omega_diff,
            );
            let lon_asc_diff =
                ((el.lon_ascending_node - lon_asc_target + PI).rem_euclid(TAU) - PI).abs();
            assert!(
                lon_asc_diff < 1e-9,
                "{label}: Ω recovered as {} rad, expected {lon_asc_target} rad (wrapped diff {})",
                el.lon_ascending_node,
                lon_asc_diff,
            );

            // Eccentricity is rotation-invariant — sanity check that the
            // 3D round-trip preserves the magnitude as well as the angles.
            assert!(
                (el.e - e).abs() < 1e-9,
                "{label}: e recovered as {} (expected {e}); rotation should not change e",
                el.e,
            );
        }
    }

    // ── Anomalies — Kepler equation, half-angle identities, time linearity ──

    /// Build a body at periapsis on an elliptical orbit of given `(a, e)`,
    /// with primary at the origin and orbit in the xy-plane (prograde).
    /// At periapsis: r = a(1−e), v ⊥ r, v = √(GM(1+e)/(a(1−e))).
    fn ellipse_at_periapsis(a: f64, e: f64, primary_mass: f64) -> (Body, Body) {
        let gm = G * primary_mass;
        let r_peri = a * (1.0 - e);
        let v_peri = (gm * (1.0 + e) / r_peri).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let satellite = body(r_peri, 0.0, 0.0, v_peri, 1e-10);
        (primary, satellite)
    }

    /// Body on an elliptical orbit at given true anomaly ν, with `(a, e)`
    /// in the xy-plane. Uses focus-conic and the in-plane frame, then
    /// rotates into the world frame via ω = 0 (no rotation needed since
    /// periapsis is on +x̂).
    fn ellipse_at_true_anomaly(a: f64, e: f64, nu: f64, primary_mass: f64) -> (Body, Body) {
        let gm = G * primary_mass;
        let p = a * (1.0 - e * e);
        let r = p / (1.0 + e * nu.cos());
        let (s_nu, c_nu) = nu.sin_cos();
        let x = r * c_nu;
        let y = r * s_nu;
        // Velocity in perifocal frame (Murray & Dermott eq. 2.36):
        // v_r = √(GM/p) · e sin ν,  v_θ = √(GM/p) · (1 + e cos ν)
        let h = (gm * p).sqrt();
        let inv_r = 1.0 / r;
        // Cartesian velocity = v_r · r̂ + v_θ · θ̂. With r̂ = (cos ν, sin ν)
        // and θ̂ = (-sin ν, cos ν).
        let v_r = (gm / p).sqrt() * e * s_nu;
        let v_theta = h * inv_r;
        let vx = v_r * c_nu - v_theta * s_nu;
        let vy = v_r * s_nu + v_theta * c_nu;
        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let satellite = body(x, y, vx, vy, 1e-10);
        (primary, satellite)
    }

    /// Kepler's equation `M = E − e sin E` must hold to round-off for any
    /// elliptical orbit at any phase.
    #[test]
    fn kepler_equation_holds_across_eccentricity_range() {
        let nu_samples = [-2.5, -1.3, -0.5, 0.0, 0.4, 1.1, 2.2, 3.0];
        for &e in &[0.05, 0.21, 0.45, 0.7, 0.85, 0.95] {
            for &nu in &nu_samples {
                let (p, s) = ellipse_at_true_anomaly(10.0, e, nu, 1e6);
                let el = elements(p, s);
                let big_e = el.eccentric_anomaly;
                let m = el.mean_anomaly;
                let lhs = m;
                let rhs = wrap_pi(big_e - e * big_e.sin());
                let diff = ((lhs - rhs + PI).rem_euclid(TAU) - PI).abs();
                assert!(
                    diff < 1e-12,
                    "Kepler eq fails at e={e}, ν={nu}: M={m}, E−e·sinE={rhs}, diff={diff}",
                );
            }
        }
    }

    /// Half-angle identity: `tan(E/2) = √((1−e)/(1+e)) · tan(ν/2)`.
    /// Exercised across the eccentricity range covering Mercury (0.21),
    /// Eris (0.44), and Sedna (0.85).
    #[test]
    fn true_to_eccentric_anomaly_half_angle_identity() {
        let nu_samples = [-2.0, -0.7, 0.3, 1.5, 2.8];
        for &e in &[0.05, 0.21, 0.44, 0.7, 0.85] {
            for &nu in &nu_samples {
                let (p, s) = ellipse_at_true_anomaly(10.0, e, nu, 1e6);
                let el = elements(p, s);
                let lhs = (0.5 * el.eccentric_anomaly).tan();
                let rhs = ((1.0 - e) / (1.0 + e)).sqrt() * (0.5 * el.true_anomaly).tan();
                let diff = (lhs - rhs).abs();
                assert!(
                    diff < 1e-10,
                    "ν↔E identity fails at e={e}, ν={nu}: tan(E/2)={lhs}, √(...)·tan(ν/2)={rhs}",
                );
            }
        }
    }

    /// At periapsis (ν = 0) all three anomalies must vanish.
    #[test]
    fn anomalies_zero_at_periapsis() {
        for &e in &[0.05, 0.21, 0.5, 0.85] {
            let (p, s) = ellipse_at_periapsis(10.0, e, 1e6);
            let el = elements(p, s);
            assert!(el.true_anomaly.abs() < 1e-12, "ν({e}) = {}", el.true_anomaly);
            assert!(el.eccentric_anomaly.abs() < 1e-12, "E({e}) = {}", el.eccentric_anomaly);
            assert!(el.mean_anomaly.abs() < 1e-12, "M({e}) = {}", el.mean_anomaly);
        }
    }

    /// Mean motion `n = 2π / T = √(GM / a³)`. Pure algebra on stored fields.
    #[test]
    fn mean_motion_matches_kepler_third_law() {
        for &e in &[0.0, 0.21, 0.7] {
            let a = 10.0;
            let m = 1e6;
            let (p, s) =
                if e < 1e-10 { circular_orbit(a, m) } else { ellipse_at_periapsis(a, e, m) };
            let el = elements(p, s);
            let expected = (G * m / (a * a * a)).sqrt();
            let got = el.mean_motion();
            assert!(
                (got - expected).abs() / expected < 1e-12,
                "n mismatch at e={e}: got {got}, expected {expected}",
            );
        }
    }

    /// `pericenter() = a(1−e)` and `apocenter() = a(1+e)` for elliptical orbits.
    #[test]
    fn pericenter_apocenter_match_ae_formulas() {
        for &(a, e) in &[(10.0, 0.0), (10.0, 0.5), (50.0, 0.85), (1.0, 0.21)] {
            let (p, s) = ellipse_at_periapsis(a, e, 1e6);
            let el = elements(p, s);
            assert!(
                (el.pericenter() - a * (1.0 - e)).abs() < 1e-9,
                "q mismatch at (a,e)=({a},{e}): got {} expected {}",
                el.pericenter(),
                a * (1.0 - e),
            );
            assert!(
                (el.apocenter() - a * (1.0 + e)).abs() < 1e-9,
                "Q mismatch at (a,e)=({a},{e}): got {} expected {}",
                el.apocenter(),
                a * (1.0 + e),
            );
        }
    }

    /// Hyperbolic identity: `M_h = e sinh H − H`, with `H` related to ν by
    /// `tanh(H/2) = √((e−1)/(e+1)) · tan(ν/2)` for ν inside the asymptote.
    #[test]
    fn hyperbolic_kepler_identity_holds() {
        // Hyperbolic flyby: e = 1.5, q = 5 (pericenter distance).
        let e = 1.5;
        let q = 5.0;
        let m = 1e6;
        let gm = G * m;
        // At periapsis: r = q, v = √(GM(1+e)/q)
        let v_peri = (gm * (1.0 + e) / q).sqrt();
        let primary = body(0.0, 0.0, 0.0, 0.0, m);
        // Place body at periapsis on +x̂, velocity in +ŷ direction.
        let satellite = body(q, 0.0, 0.0, v_peri, 1e-10);
        let bodies = vec![primary, satellite];
        let el = compute_elements(&bodies, 1, 0, G).unwrap();
        assert_eq!(el.orbit_type, OrbitType::Hyperbolic);
        // At periapsis: ν = H = M_h = 0.
        assert!(el.true_anomaly.abs() < 1e-12);
        assert!(el.eccentric_anomaly.abs() < 1e-12);
        assert!(el.mean_anomaly.abs() < 1e-12);

        // Off-periapsis sample: place body at ν = 0.5 rad, verify identity.
        let nu: f64 = 0.5;
        let p = q * (1.0 + e); // semi-latus rectum = a(1−e²); for hyperbolic q = a(e−1) so a = q/(e−1) and p = q(e+1)/1 ... actually p = q(1+e). check.
        // For hyperbolic: q = a(e−1) with a < 0 by convention or |a| convention.
        // Using p = h²/GM = q(1+e), valid for any conic.
        let r = p / (1.0 + e * nu.cos());
        let (s_nu, c_nu) = nu.sin_cos();
        let x = r * c_nu;
        let y = r * s_nu;
        let h = (gm * p).sqrt();
        let v_r = (gm / p).sqrt() * e * s_nu;
        let v_theta = h / r;
        let vx = v_r * c_nu - v_theta * s_nu;
        let vy = v_r * s_nu + v_theta * c_nu;
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(x, y, vx, vy, 1e-10)];
        let el = compute_elements(&bodies, 1, 0, G).unwrap();
        assert_eq!(el.orbit_type, OrbitType::Hyperbolic);
        let big_h = el.eccentric_anomaly;
        let m_h = el.mean_anomaly;
        let kepler_residual = (m_h - (e * big_h.sinh() - big_h)).abs();
        assert!(
            kepler_residual < 1e-10,
            "Hyperbolic Kepler fails: M_h={m_h}, e·sinh(H)−H={}, residual={kepler_residual}",
            e * big_h.sinh() - big_h,
        );
    }

    /// Inclined elliptical orbit (Eris-like, i = 44°): inclination must be
    /// recovered correctly and anomalies must satisfy Kepler's equation just
    /// the same as in the planar case.
    #[test]
    fn anomalies_correct_on_inclined_orbit_eris_like() {
        // Build a body on an ellipse in the xy-plane, then rotate the
        // entire (r, v) by inclination i around the +x̂ axis. The
        // ascending node sits on +x̂.
        let a = 50.0;
        let e = 0.44;
        let nu = 1.2;
        let i = 44.0_f64.to_radians();
        let m = 1e6;

        let (p_planar, s_planar) = ellipse_at_true_anomaly(a, e, nu, m);
        // Rotate (r, v) around x̂: y' = y cos i, z' = y sin i; same for v.
        let (s_i, c_i) = i.sin_cos();
        let mut s =
            body(s_planar.pos_x, s_planar.pos_y * c_i, s_planar.vel_x, s_planar.vel_y * c_i, 1e-10);
        s.pos_z = s_planar.pos_y * s_i;
        s.vel_z = s_planar.vel_y * s_i;
        s.sync_physical_properties();

        let bodies = vec![p_planar, s];
        let el = compute_elements(&bodies, 1, 0, G).unwrap();

        // Inclination matches.
        assert!(
            (el.inclination - i).abs() < 1e-10,
            "i: got {} expected {}",
            el.inclination.to_degrees(),
            i.to_degrees(),
        );
        // True anomaly preserved through rotation (rotation does not move
        // the body along the orbit).
        assert!((el.true_anomaly - nu).abs() < 1e-9, "ν: got {} expected {}", el.true_anomaly, nu,);
        // Kepler's equation still holds.
        let kep = wrap_pi(el.eccentric_anomaly - e * el.eccentric_anomaly.sin());
        let diff = ((el.mean_anomaly - kep + PI).rem_euclid(TAU) - PI).abs();
        assert!(diff < 1e-12, "Kepler eq on inclined orbit: M={}, E−e·sinE={kep}", el.mean_anomaly);
    }

    // ── Reversibility: state → elements → state ──────────────────────────────
    //
    // Reconstruct (r, v) from a full Keplerian element set via the
    // perifocal-to-world rotation `R₃(Ω) · R₁(i) · R₃(ω)`. Test helper —
    // promote to public API in a follow-up if Playback mode or the Python
    // binding ever needs to construct bodies from elements.

    /// Build a (r_world, v_world) pair from `(a, e, i, Ω, ω, ν, μ)` via
    /// the standard astrodynamics perifocal-to-inertial rotation. For
    /// hyperbolic input, `a` should be negative (codebase convention) and
    /// the same conic formulas apply with `p = a(1 − e²)` retaining its
    /// sign so `r = p / (1 + e cos ν)` is positive on the active branch.
    fn reconstruct_state(
        a: f64,
        e: f64,
        i: f64,
        lon_asc: f64,
        omega: f64,
        nu: f64,
        mu: f64,
    ) -> (Vec3, Vec3) {
        let p = a * (1.0 - e * e);
        let (s_nu, c_nu) = nu.sin_cos();
        let r = p / (1.0 + e * c_nu);

        let r_pf = Vec3::new(r * c_nu, r * s_nu, 0.0);
        let factor = (mu / p.abs()).sqrt();
        let v_pf = Vec3::new(-factor * s_nu, factor * (e + c_nu), 0.0);

        // Composite rotation R₃(Ω) · R₁(i) · R₃(ω) — same matrix as in
        // `OrbitalElements::sample_orbit`, full 3×3 form here.
        let (s_o, c_o) = omega.sin_cos();
        let (s_i, c_i) = i.sin_cos();
        let (s_w, c_w) = lon_asc.sin_cos();

        let r11 = c_w * c_o - s_w * s_o * c_i;
        let r12 = -c_w * s_o - s_w * c_o * c_i;
        let r21 = s_w * c_o + c_w * s_o * c_i;
        let r22 = -s_w * s_o + c_w * c_o * c_i;
        let r31 = s_o * s_i;
        let r32 = c_o * s_i;

        let rotate = |p: Vec3| -> Vec3 {
            Vec3::new(r11 * p.x + r12 * p.y, r21 * p.x + r22 * p.y, r31 * p.x + r32 * p.y)
        };
        (rotate(r_pf), rotate(v_pf))
    }

    /// Run a state → elements → state round-trip on a chosen regime.
    /// Asserts the recovered position and velocity match the original to
    /// `1e-9` relative.
    fn assert_state_roundtrip(
        label: &str,
        a: f64,
        e: f64,
        i_deg: f64,
        lon_asc_deg: f64,
        omega_deg: f64,
        nu_deg: f64,
    ) {
        let primary_mass = 1e6;
        let mu = G * primary_mass;
        let i = i_deg.to_radians();
        let lon_asc = lon_asc_deg.to_radians();
        let omega = omega_deg.to_radians();
        let nu = nu_deg.to_radians();

        let (r0, v0) = reconstruct_state(a, e, i, lon_asc, omega, nu, mu);

        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let mut sat = body(r0.x, r0.y, v0.x, v0.y, 1e-10);
        sat.pos_z = r0.z;
        sat.vel_z = v0.z;
        sat.sync_physical_properties();

        let bodies = vec![primary, sat];
        let el = compute_elements(&bodies, 1, 0, G).expect("regime is computable");

        let (r1, v1) = reconstruct_state(
            el.a,
            el.e,
            el.inclination,
            el.lon_ascending_node,
            el.omega,
            el.true_anomaly,
            mu,
        );

        let dr = (r1 - r0).length() / r0.length().max(1e-30);
        let dv = (v1 - v0).length() / v0.length().max(1e-30);
        assert!(dr < 1e-9, "{label}: |Δr|/|r| = {dr:.3e} (r0={r0:?}, r1={r1:?})");
        assert!(dv < 1e-9, "{label}: |Δv|/|v| = {dv:.3e} (v0={v0:?}, v1={v1:?})");
    }

    /// Mercury-like: low eccentricity, moderate inclination, full Ω/ω.
    #[test]
    fn state_roundtrip_mercury_like() {
        assert_state_roundtrip("Mercury-like", 0.387, 0.21, 7.0, 48.331, 29.124, 174.79);
    }

    /// Eris-like: extreme inclination (44°), moderate eccentricity. The
    /// regime where 2D-projection of the orbit collapses scientific
    /// information.
    #[test]
    fn state_roundtrip_eris_like() {
        assert_state_roundtrip("Eris-like", 68.0, 0.44, 44.04, 35.95, 151.0, 90.0);
    }

    /// Sedna-like: extreme eccentricity (0.85). Stresses the
    /// `(1−e).sqrt()` factor inside the half-angle E recovery.
    #[test]
    fn state_roundtrip_sedna_like() {
        assert_state_roundtrip("Sedna-like", 507.0, 0.85, 11.93, 144.0, 311.0, 45.0);
    }

    /// Hyperbolic flyby: `e > 1` with non-trivial inclination.
    #[test]
    fn state_roundtrip_hyperbolic_flyby() {
        // a < 0 by convention for hyperbolic; ν = 0.5 rad keeps the body
        // well inside the asymptote `ν_∞ = acos(−1/e) ≈ 2.30` for e = 1.5.
        let nu_deg = 0.5_f64.to_degrees();
        assert_state_roundtrip("Hyperbolic", -10.0, 1.5, 20.0, 60.0, 45.0, nu_deg);
    }

    /// Circular planar: degenerate fallback path (`e ≈ 0`, `Ω = 0`,
    /// `ω = 0`, ν taken as argument of latitude).
    #[test]
    fn state_roundtrip_circular_planar() {
        let primary_mass = 1e6;
        let mu = G * primary_mass;
        let a = 10.0;
        let e = 0.0;
        let nu = std::f64::consts::FRAC_PI_4;

        let (r0, v0) = reconstruct_state(a, e, 0.0, 0.0, 0.0, nu, mu);
        let primary = body(0.0, 0.0, 0.0, 0.0, primary_mass);
        let mut sat = body(r0.x, r0.y, v0.x, v0.y, 1e-10);
        sat.pos_z = r0.z;
        sat.vel_z = v0.z;
        sat.sync_physical_properties();
        let bodies = vec![primary, sat];
        let el = compute_elements(&bodies, 1, 0, G).expect("circular regime is computable");

        let (r1, v1) = reconstruct_state(
            el.a,
            el.e,
            el.inclination,
            el.lon_ascending_node,
            el.omega,
            el.true_anomaly,
            mu,
        );

        let dr = (r1 - r0).length() / r0.length();
        let dv = (v1 - v0).length() / v0.length();
        assert!(dr < 1e-9, "circular: |Δr|/|r| = {dr:.3e}");
        assert!(dv < 1e-9, "circular: |Δv|/|v| = {dv:.3e}");
    }

    // ── hierarchical_primary — Hill-sphere parent, distinct from
    //                          the strongest-attractor `dominant_primary` ──

    /// Make a body at world position (x, y, z) with given mass. Velocity is
    /// zero — `hierarchical_primary` does not consult velocity except in
    /// the energy fallback, which is exercised separately.
    fn body3(x: f64, y: f64, z: f64, mass: f64) -> Body {
        let mut b = Body::rocky(mass).at(x, y);
        b.pos_z = z;
        b.sync_physical_properties();
        b
    }

    /// Sun + Earth + Moon at SI distances. Sun pulls the Moon ~2× harder
    /// than Earth does (`G·M☉/r_⊙² > G·M⊕/r⊕²`), so `dominant_primary`
    /// returns the Sun for the Moon. Hierarchical primary must instead
    /// return Earth — the body whose Hill sphere contains the Moon.
    #[test]
    fn hierarchical_primary_recovers_earth_for_the_moon() {
        // Indices: 0 = Sun, 1 = Earth, 2 = Moon.
        let m_sun = 1.989e30;
        let m_earth = 5.972e24;
        let m_moon = 7.342e22;
        let au = 1.495_978_707e11;
        let r_em = 3.844e8;
        let bodies = vec![
            body3(0.0, 0.0, 0.0, m_sun),
            body3(au, 0.0, 0.0, m_earth),
            body3(au + r_em, 0.0, 0.0, m_moon),
        ];

        // Sanity: confirm dominant_primary picks the Sun (the documented
        // 2× ratio holds — anyone reading this test should find the
        // motivation in the function-level rustdoc).
        assert_eq!(dominant_primary(&bodies, 2), Some(0), "dominant_primary(Moon) = Sun");

        // Hierarchical primary recovers Earth — established via Hill sphere.
        assert_eq!(
            hierarchical_primary(&bodies, 2),
            Some((1, HierarchicalRelation::HillSphere)),
            "hierarchical_primary(Moon) = Earth via Hill sphere",
        );
    }

    /// For a body whose dominant primary already coincides with its
    /// hierarchical primary (no Hill-sphere divergence), both functions
    /// return the same answer.
    #[test]
    fn hierarchical_primary_matches_dominant_for_inner_planets() {
        let m_sun = 1.989e30;
        let m_mercury = 3.302e23;
        let m_earth = 5.972e24;
        let au = 1.495_978_707e11;
        let bodies = vec![
            body3(0.0, 0.0, 0.0, m_sun),
            body3(0.387 * au, 0.0, 0.0, m_mercury),
            body3(au, 0.0, 0.0, m_earth),
        ];

        for idx in [1, 2] {
            let h = hierarchical_primary(&bodies, idx).map(|(i, _)| i);
            assert_eq!(
                h,
                dominant_primary(&bodies, idx),
                "non-hierarchical body {idx}: both primaries should coincide on the Sun",
            );
        }
    }

    /// Inclined Earth–Moon configuration — the Hill-sphere check must
    /// reduce a 3D distance, not a planar projection. Place the Moon
    /// straight above Earth (along `+ẑ`) and confirm Earth still wins.
    #[test]
    fn hierarchical_primary_handles_inclined_separation() {
        let m_sun = 1.989e30;
        let m_earth = 5.972e24;
        let m_moon = 7.342e22;
        let au = 1.495_978_707e11;
        let r_em = 3.844e8;
        let bodies = vec![
            body3(0.0, 0.0, 0.0, m_sun),
            body3(au, 0.0, 0.0, m_earth),
            // Moon directly above Earth — the planar (x, y) projection
            // collapses; only the 3D distance recovers the geometry.
            body3(au, 0.0, r_em, m_moon),
        ];

        assert_eq!(
            hierarchical_primary(&bodies, 2),
            Some((1, HierarchicalRelation::HillSphere)),
            "Moon stacked along ẑ above Earth must still resolve to Earth",
        );
    }

    /// A truly isolated body — the only body more massive than itself is
    /// the system root, captured via the infinite-Hill-sphere branch.
    #[test]
    fn hierarchical_primary_returns_root_for_unfiltered_planet() {
        let m_sun = 1.989e30;
        let m_earth = 5.972e24;
        let au = 1.495_978_707e11;
        let bodies = vec![body3(0.0, 0.0, 0.0, m_sun), body3(au, 0.0, 0.0, m_earth)];
        assert_eq!(
            hierarchical_primary(&bodies, 1),
            Some((0, HierarchicalRelation::HillSphere)),
            "Earth's hierarchical primary is the Sun (system root)",
        );
    }

    /// The most massive body has no parent; result is `None`.
    #[test]
    fn hierarchical_primary_returns_none_for_heaviest_body() {
        let bodies =
            vec![body3(0.0, 0.0, 0.0, 1.989e30), body3(1.495_978_707e11, 0.0, 0.0, 5.972e24)];
        assert!(hierarchical_primary(&bodies, 0).is_none());
    }

    /// Single-body and empty inputs short-circuit to `None` without
    /// inspecting any state.
    #[test]
    fn hierarchical_primary_returns_none_for_trivial_input() {
        let bodies = vec![body3(0.0, 0.0, 0.0, 1.0e30)];
        assert!(hierarchical_primary(&bodies, 0).is_none());
    }

    // ── is_system_root — shared predicate consumed by inspector and
    //                    canvas to skip orbit rendering for the heaviest
    //                    body ──

    #[test]
    fn is_system_root_marks_the_heaviest_body() {
        let bodies = vec![
            body3(0.0, 0.0, 0.0, 1.989e30),     // Sun
            body3(1.5e11, 0.0, 0.0, 5.972e24),  // Earth
            body3(7.78e11, 0.0, 0.0, 1.898e27), // Jupiter
        ];
        assert!(is_system_root(&bodies, 0), "Sun is system root");
        assert!(!is_system_root(&bodies, 1), "Earth has the Sun above it");
        assert!(!is_system_root(&bodies, 2), "Jupiter has the Sun above it");
    }

    #[test]
    fn is_system_root_handles_a_single_body() {
        let bodies = vec![body3(0.0, 0.0, 0.0, 1.0e30)];
        assert!(is_system_root(&bodies, 0));
    }

    #[test]
    fn is_system_root_returns_false_for_out_of_range_index() {
        let bodies = vec![body3(0.0, 0.0, 0.0, 1.0e30)];
        assert!(!is_system_root(&bodies, 5));
    }

    /// Equal-mass siblings: neither has a body strictly heavier, so both
    /// register as roots. The application layer can disambiguate (e.g.
    /// pick a deterministic representative) when this matters.
    #[test]
    fn is_system_root_treats_equal_mass_pair_as_co_roots() {
        let bodies = vec![body3(0.0, 0.0, 0.0, 1.0e30), body3(1.0e11, 0.0, 0.0, 1.0e30)];
        assert!(is_system_root(&bodies, 0));
        assert!(is_system_root(&bodies, 1));
    }
}
