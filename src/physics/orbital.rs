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

use crate::core::body::Body;

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
            let score_j = if rj2 > 0.0 {
                bj.mass / rj2
            } else {
                f64::INFINITY
            };
            let score_k = if rk2 > 0.0 {
                bk.mass / rk2
            } else {
                f64::INFINITY
            };
            score_j
                .partial_cmp(&score_k)
                .unwrap_or(std::cmp::Ordering::Equal)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{body::Body, materials::Material};
    use std::f64::consts::{PI, TAU};

    // ── Helpers ───────────────────────────────────────────────────────────────

    const G: f64 = 1.0;
    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::new(x, y, vx, vy, mass, Material::Rocky)
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
        assert!(
            el.energy < 0.0,
            "energia = {}, deve ser negativa",
            el.energy
        );
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
        assert_eq!(
            primary, 1,
            "primário deve ser a Terra (índice 1), não o Sol"
        );
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
