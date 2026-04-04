use std::f64::consts::PI;

/// Base softening length for a body of mass 1.0.
/// Per-body softening scales as `EPS_BASE * mass^(1/3)`, so each body's
/// softening volume is proportional to its mass — physically motivated by
/// the Plummer-equivalent equal-mass softening criterion.
pub const EPS_BASE: f64 = 0.02;

/// Placeholder radius for a body of mass 1.0 before `System::calibrate_radii`
/// is called.  Uses the same 3-D mass scaling as softening.
const R_PLACEHOLDER: f64 = EPS_BASE * 0.5;

#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub mass: f64,

    /// Gravitational softening length ε for this body.
    /// Pairwise: ε²_ij = (ε²_i + ε²_j) / 2.
    /// Calibrated by `System::calibrate_softening`.
    pub softening: f64,

    /// Physical collision radius.  A pair collides when their separation falls
    /// below r_i + r_j.  Calibrated by `System::calibrate_radii`.
    ///
    /// Invariant: radius ≤ softening (kept by calibration so that the force
    /// is already in the softened regime when two bodies touch).
    pub radius: f64,

    /// Bulk density of the body: ρ = m / V, V = 4/3 π r³.
    ///
    /// This is the **primary size property** of a body — radius is derived from
    /// it via `r = (3m / 4πρ)^(1/3)`.  Storing density instead of radius alone
    /// lets merges correctly compute the combined volume and thus the merged body's
    /// radius, instead of assuming constant density for all bodies.
    pub density: f64,

    /// Angular velocity around the z-axis: ω_z = L_z / I.
    /// Updated after inelastic collisions to conserve the internal angular
    /// momentum lost in the merge.  Used for visualization and collision
    /// response in elastic/partial-restitution collisions (Phase 2).
    pub omega_z: f64,

    /// Moment of inertia around z-axis: I_z = (2/5)·m·r² for a uniform sphere.
    /// Used to compute ω_z from the internal angular momentum after a merge.
    /// For a point mass in 2-D, this is a notional measure; the important
    /// physics is the conversion L_internal → ω_z for visual effects.
    pub moment_inertia: f64,

    /// Optional custom display color [R, G, B].  When `None` the renderer falls
    /// back to the index-based palette in `theme::body_color`.
    pub color: Option<[u8; 3]>,
}

impl Body {
    pub fn new(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Self {
        let softening = default_softening(mass);
        let radius = default_radius(mass);
        let density = density_from_mass_radius(mass, radius);
        Self {
            x,
            y,
            vx,
            vy,
            mass,
            softening,
            radius,
            density,
            omega_z: 0.0,
            moment_inertia: default_moment_inertia(mass, radius),
            color: None,
        }
    }
}

/// Default softening before system-scale calibration.
pub fn default_softening(mass: f64) -> f64 {
    EPS_BASE * mass.abs().cbrt()
}

/// Default collision radius before system-scale calibration.
/// Kept smaller than softening so the force is already softened at contact.
pub fn default_radius(mass: f64) -> f64 {
    R_PLACEHOLDER * mass.abs().cbrt()
}

/// Moment of inertia for a uniform sphere: I = (2/5)·m·r².
/// Used to convert internal angular momentum to rotational velocity.
pub fn default_moment_inertia(mass: f64, radius: f64) -> f64 {
    0.4 * mass * radius * radius
}

/// Density from mass and radius: ρ = m / (4/3 π r³).
///
/// Returns a safe positive fallback when `radius ≤ 0`.
pub fn density_from_mass_radius(mass: f64, radius: f64) -> f64 {
    let vol = sphere_volume(radius);
    if vol > 0.0 {
        mass / vol
    } else {
        // Degenerate body — assign unit density so radius can be recovered later
        1.0
    }
}

/// Radius from density and mass: r = (3m / 4πρ)^(1/3).
pub fn radius_from_density_mass(density: f64, mass: f64) -> f64 {
    let vol = mass / density.max(1e-30);
    sphere_radius_from_volume(vol)
}

/// Volume of a sphere: V = 4/3 π r³.
#[inline]
pub fn sphere_volume(radius: f64) -> f64 {
    (4.0 / 3.0) * PI * radius.powi(3)
}

/// Radius of a sphere given its volume: r = (3V / 4π)^(1/3).
#[inline]
pub fn sphere_radius_from_volume(volume: f64) -> f64 {
    ((3.0 * volume) / (4.0 * PI)).cbrt()
}
