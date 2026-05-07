//! Helpers for constructing initial states from orbital parameters.
//!
//! Templates pull from these instead of computing trig and rotations inline,
//! so a Keplerian setup reads as one line per body and the heavy math is
//! audited in one place.

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
