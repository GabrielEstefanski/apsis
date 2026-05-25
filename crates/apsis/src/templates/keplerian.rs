//! Keplerian elements → Cartesian state vector.
//!
//! Standard astrodynamics conversion. Solves Kepler's equation by
//! Newton–Raphson, computes the state in the perifocal frame, and
//! rotates by the 3-1-3 Euler sequence (Ω, i, ω) into the inertial
//! frame.
//!
//! For heliocentric bodies the inertial frame is the ecliptic J2000
//! frame and `state_from_elements` is enough. For satellites whose
//! published elements are given relative to a parent's equator,
//! [`parent_equator_basis`] returns the rotation that takes a vector
//! from the parent's equatorial frame to the ecliptic frame; the
//! caller composes this with the moon's orbital rotation through
//! [`state_from_elements_in_basis`].
//!
//! All angles are in **radians**. `mu = G · M_central` in simulation
//! units (with the simulator's implicit `G = 1`, this is just the
//! central mass).
//!
//! Reference: Curtis, *Orbital Mechanics for Engineering Students*,
//! 3e, §4.4 (state from orbital elements) and §11.5 (Newton–Raphson
//! solver for Kepler's equation).

/// Obliquity of the ecliptic at J2000 [rad]. From IAU 2006 precession
/// model: 23°26′21.40″ ≈ 23.439279°.
pub const ECLIPTIC_OBLIQUITY: f64 = 0.40909280422232897; // 23.4393° in radians

/// Maximum Newton–Raphson iterations for Kepler's equation. Ten is
/// already more than enough for any e ≤ 0.99 to converge to f64 epsilon.
const KEPLER_MAX_ITER: usize = 32;
/// Convergence tolerance for the eccentric anomaly residual.
const KEPLER_TOL: f64 = 1e-14;

/// Solve Kepler's equation `M = E − e·sin(E)` for the eccentric
/// anomaly `E`, given mean anomaly `M` and eccentricity `e ∈ [0, 1)`.
///
/// Uses Newton–Raphson with the standard `E₀ = M + e·sin(M)` seed,
/// which converges quadratically across the elliptic regime.
pub fn solve_kepler(mean_anom: f64, e: f64) -> f64 {
    // Wrap M into [-π, π] first so the Newton iteration starts in a
    // monotone region of the residual.
    let m = mean_anom.rem_euclid(std::f64::consts::TAU);
    let m = if m > std::f64::consts::PI { m - std::f64::consts::TAU } else { m };

    let mut e_anom = m + e * m.sin();
    for _ in 0..KEPLER_MAX_ITER {
        let f = e_anom - e * e_anom.sin() - m;
        let fp = 1.0 - e * e_anom.cos();
        let delta = f / fp;
        e_anom -= delta;
        if delta.abs() < KEPLER_TOL {
            break;
        }
    }
    e_anom
}

/// Convert classical orbital elements to a Cartesian state vector
/// `(position, velocity)` in the central body's inertial frame.
///
/// * `mu` — `G · M_central` (with `G = 1` in simulation units, this is
///   just `M_central`).
/// * `a` — semi-major axis (length).
/// * `e` — eccentricity, `0 ≤ e < 1` (elliptic orbits only).
/// * `inc` — inclination relative to the reference plane [rad].
/// * `raan` — longitude of ascending node Ω [rad].
/// * `argp` — argument of periapsis ω [rad].
/// * `mean_anom` — mean anomaly at epoch M [rad].
///
/// The reference plane is whatever the elements are quoted in; for
/// heliocentric planet elements that's the ecliptic, so the returned
/// state is heliocentric ecliptic.
pub fn state_from_elements(
    mu: f64,
    a: f64,
    e: f64,
    inc: f64,
    raan: f64,
    argp: f64,
    mean_anom: f64,
) -> ([f64; 3], [f64; 3]) {
    // ── Solve Kepler for the eccentric anomaly ────────────────────────────
    let e_anom = solve_kepler(mean_anom, e);

    // True anomaly via the half-angle identity: avoids the sign ambiguity
    // of atan2 and stays well-conditioned for high e.
    let half = 0.5 * e_anom;
    let true_anom = 2.0 * (((1.0 + e).sqrt() * half.sin()).atan2((1.0 - e).sqrt() * half.cos()));

    // ── Perifocal frame ────────────────────────────────────────────────────
    let r = a * (1.0 - e * e_anom.cos());
    let p_xy = [r * true_anom.cos(), r * true_anom.sin()];

    // Velocity from the vis-viva relation in perifocal coords.
    let p_param = a * (1.0 - e * e); // semi-latus rectum
    let coef = (mu / p_param).sqrt();
    let v_xy = [-coef * true_anom.sin(), coef * (e + true_anom.cos())];

    // ── Rotate perifocal → reference frame via 3-1-3 Euler (Ω, i, ω) ──────
    let (cos_o, sin_o) = (raan.cos(), raan.sin());
    let (cos_i, sin_i) = (inc.cos(), inc.sin());
    let (cos_w, sin_w) = (argp.cos(), argp.sin());

    // R = Rz(Ω) · Rx(i) · Rz(ω). Expanded form (column-major, applied
    // to a column vector):
    let r11 = cos_o * cos_w - sin_o * sin_w * cos_i;
    let r12 = -cos_o * sin_w - sin_o * cos_w * cos_i;
    let r21 = sin_o * cos_w + cos_o * sin_w * cos_i;
    let r22 = -sin_o * sin_w + cos_o * cos_w * cos_i;
    let r31 = sin_w * sin_i;
    let r32 = cos_w * sin_i;

    let pos = [
        r11 * p_xy[0] + r12 * p_xy[1],
        r21 * p_xy[0] + r22 * p_xy[1],
        r31 * p_xy[0] + r32 * p_xy[1],
    ];
    let vel = [
        r11 * v_xy[0] + r12 * v_xy[1],
        r21 * v_xy[0] + r22 * v_xy[1],
        r31 * v_xy[0] + r32 * v_xy[1],
    ];

    (pos, vel)
}

/// Same as [`state_from_elements`] but the resulting state is then
/// rotated through `basis` — a 3×3 column-major matrix that takes a
/// vector from the elements' reference frame to the desired output
/// frame (typically heliocentric ecliptic).
///
/// Use [`parent_equator_basis`] to build `basis` for satellites whose
/// elements are quoted relative to a parent's equator.
pub fn state_from_elements_in_basis(
    mu: f64,
    a: f64,
    e: f64,
    inc: f64,
    raan: f64,
    argp: f64,
    mean_anom: f64,
    basis: [[f64; 3]; 3],
) -> ([f64; 3], [f64; 3]) {
    let (pos, vel) = state_from_elements(mu, a, e, inc, raan, argp, mean_anom);
    (apply_basis(&basis, pos), apply_basis(&basis, vel))
}

#[inline]
fn apply_basis(basis: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        basis[0][0] * v[0] + basis[0][1] * v[1] + basis[0][2] * v[2],
        basis[1][0] * v[0] + basis[1][1] * v[1] + basis[1][2] * v[2],
        basis[2][0] * v[0] + basis[2][1] * v[1] + basis[2][2] * v[2],
    ]
}

/// Unit vector pointing along a body's spin axis, expressed in the
/// heliocentric ecliptic J2000 frame.
///
/// Inputs are the body's pole right ascension and declination in the
/// equatorial J2000 frame [degrees] — the standard form in which IAU
/// publishes pole orientations (WGCCRE reports).
///
/// The conversion is one rotation by `−ε` around the X axis: the J2000
/// equatorial and ecliptic frames share their X axis (vernal equinox),
/// and the ecliptic is tilted `+ε` from the equator about that axis.
pub fn ecliptic_pole_from_ra_dec(pole_ra_deg: f64, pole_dec_deg: f64) -> [f64; 3] {
    let ra = pole_ra_deg.to_radians();
    let dec = pole_dec_deg.to_radians();

    // Equatorial unit vector.
    let xe = dec.cos() * ra.cos();
    let ye = dec.cos() * ra.sin();
    let ze = dec.sin();

    // Rotate by −ε around X to get the ecliptic-frame components.
    let cos_eps = ECLIPTIC_OBLIQUITY.cos();
    let sin_eps = ECLIPTIC_OBLIQUITY.sin();
    [xe, cos_eps * ye + sin_eps * ze, -sin_eps * ye + cos_eps * ze]
}

/// Build the 3×3 rotation that takes a vector in a parent body's
/// equatorial frame to the heliocentric ecliptic J2000 frame.
///
/// The parent's equator-z axis is its spin axis (`pole`); the
/// equator-x axis is the line of nodes between the parent equator and
/// the ecliptic (so that `equator-x` lies in the ecliptic plane and
/// the prime meridian is taken as the ascending node by convention);
/// equator-y completes the right-handed basis.
///
/// Returned as a row-major 3×3 with columns `[equator_x, equator_y, equator_z]`,
/// suitable for [`state_from_elements_in_basis`].
pub fn parent_equator_basis(pole_ra_deg: f64, pole_dec_deg: f64) -> [[f64; 3]; 3] {
    let pole = ecliptic_pole_from_ra_dec(pole_ra_deg, pole_dec_deg);

    // ascending node direction = ecliptic_z × pole, normalised
    let ecl_z = [0.0_f64, 0.0, 1.0];
    let mut x = cross(ecl_z, pole);
    let len = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    if len > 1e-12 {
        x = [x[0] / len, x[1] / len, x[2] / len];
    } else {
        // Degenerate: pole parallel to ecliptic z (equatorial frame ==
        // ecliptic frame). Pick the +X axis as a fallback.
        x = [1.0, 0.0, 0.0];
    }
    let y = cross(pole, x); // already unit since pole ⟂ x

    // Row-major matrix: row i, column j. Column j = j-th basis vector.
    [[x[0], y[0], pole[0]], [x[1], y[1], pole[1]], [x[2], y[2], pole[2]]]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn solve_kepler_zero_mean_anomaly_returns_zero() {
        assert!(approx_eq(solve_kepler(0.0, 0.5), 0.0, 1e-12));
    }

    #[test]
    fn solve_kepler_circular_returns_mean_anomaly() {
        // e = 0 ⇒ E = M exactly.
        for &m in &[0.1_f64, 1.0, 2.5, 4.7] {
            let e_anom = solve_kepler(m, 0.0);
            let wrapped = m.rem_euclid(std::f64::consts::TAU);
            let wrapped = if wrapped > PI { wrapped - std::f64::consts::TAU } else { wrapped };
            assert!(approx_eq(e_anom, wrapped, 1e-12));
        }
    }

    #[test]
    fn solve_kepler_high_eccentricity_satisfies_equation() {
        let e = 0.9;
        for &m in &[0.0_f64, 0.5, 1.5, 3.1] {
            let e_anom = solve_kepler(m, e);
            let m_check = e_anom - e * e_anom.sin();
            let wrapped = m.rem_euclid(std::f64::consts::TAU);
            let wrapped = if wrapped > PI { wrapped - std::f64::consts::TAU } else { wrapped };
            assert!(
                approx_eq(m_check, wrapped, 1e-10),
                "Kepler residual: M={m} e={e} → E={e_anom}, back-substituted M={m_check}"
            );
        }
    }

    #[test]
    fn circular_equatorial_orbit_lies_in_xy_plane() {
        // a=1, e=0, i=0, all angles 0 → planar circular orbit, body at +X
        // moving in +Y at speed sqrt(mu/a) = 1.
        let (pos, vel) = state_from_elements(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(approx_eq(pos[0], 1.0, 1e-12));
        assert!(approx_eq(pos[1], 0.0, 1e-12));
        assert!(approx_eq(pos[2], 0.0, 1e-12));
        assert!(approx_eq(vel[0], 0.0, 1e-12));
        assert!(approx_eq(vel[1], 1.0, 1e-12));
        assert!(approx_eq(vel[2], 0.0, 1e-12));
    }

    #[test]
    fn polar_orbit_lies_in_xz_plane() {
        // i = π/2, all other angles 0 → orbit goes through +X and +Z.
        let (pos, vel) = state_from_elements(1.0, 1.0, 0.0, PI / 2.0, 0.0, 0.0, 0.0);
        assert!(approx_eq(pos[0], 1.0, 1e-12));
        assert!(approx_eq(pos[1], 0.0, 1e-12));
        assert!(approx_eq(pos[2], 0.0, 1e-12));
        // Velocity rotated by 90° around X: was (0, 1, 0), becomes (0, 0, 1).
        assert!(approx_eq(vel[0], 0.0, 1e-12));
        assert!(approx_eq(vel[1], 0.0, 1e-12));
        assert!(approx_eq(vel[2], 1.0, 1e-12));
    }

    #[test]
    fn ecliptic_pole_from_north_ecliptic_pole_is_z() {
        // The north ecliptic pole, in equatorial J2000 coords, sits at
        // RA = 270°, DEC = 90° − ε ≈ 66.56°.
        let pole = ecliptic_pole_from_ra_dec(270.0, 90.0 - ECLIPTIC_OBLIQUITY.to_degrees());
        assert!(approx_eq(pole[0], 0.0, 1e-3));
        assert!(approx_eq(pole[1], 0.0, 1e-3));
        assert!(approx_eq(pole[2], 1.0, 1e-3));
    }

    #[test]
    fn parent_equator_basis_is_orthonormal() {
        let basis = parent_equator_basis(268.057, 64.495); // Jupiter
        let col = |j: usize| [basis[0][j], basis[1][j], basis[2][j]];
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

        let x = col(0);
        let y = col(1);
        let z = col(2);
        assert!(approx_eq(dot(x, x), 1.0, 1e-12));
        assert!(approx_eq(dot(y, y), 1.0, 1e-12));
        assert!(approx_eq(dot(z, z), 1.0, 1e-12));
        assert!(approx_eq(dot(x, y), 0.0, 1e-12));
        assert!(approx_eq(dot(x, z), 0.0, 1e-12));
        assert!(approx_eq(dot(y, z), 0.0, 1e-12));
    }
}
