use crate::domain::materials::{Material, density};
use std::f64::consts::PI;

/// Base softening length for a body of mass 1.0.
/// Per-body softening scales as `EPS_BASE * mass^(1/3)`, so each body's
/// softening volume is proportional to its mass — physically motivated by
/// the Plummer-equivalent equal-mass softening criterion.
pub const EPS_BASE: f64 = 0.02;

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

    /// True physical radius derived from mass and density.
    ///
    /// This represents the actual size of the body and is used for:
    /// - energy calculations (e.g. disruption threshold Q*)
    /// - moment of inertia
    /// - physically meaningful scaling
    ///
    /// Unlike `radius`, this value is **never modified by calibration**.
    pub physical_radius: f64,

    /// Bulk density of the body: ρ = m / V, V = 4/3 π r³.
    ///
    /// This is the **primary size property** of a body — the physical radius
    /// is derived from it via `r = (3m / 4πρ)^(1/3)`.
    ///
    /// This value is invariant during simulation except for merge/fragmentation
    /// events where material composition changes.
    pub density: f64,

    /// Angular velocity around the z-axis: ω_z = L_z / I.
    pub omega_z: f64,

    /// Moment of inertia around z-axis: I_z = (2/5)·m·r² for a uniform sphere.
    ///
    /// NOTE: This is computed using the **physical radius**, not the collision radius.
    pub moment_inertia: f64,

    /// Astrophysical material class.
    pub material: Material,

    /// Display colour [R, G, B].
    pub color: [u8; 3],
}

impl Body {
    pub fn new(x: f64, y: f64, vx: f64, vy: f64, mass: f64, material: Material) -> Self {
        let density = density(material, mass);

        // True physical radius
        let physical_radius = radius_from_density_mass(density, mass);

        let softening = default_softening(mass);

        Self {
            x,
            y,
            vx,
            vy,
            mass,
            softening,
            physical_radius,
            density,
            omega_z: 0.0,
            moment_inertia: default_moment_inertia(mass, physical_radius),
            material,
            color: material.props().base_color,
        }
    }

    /// Recompute physical-only quantities from the current mass and density.
    ///
    /// This must be used whenever mass or density changes. It intentionally
    /// does **not** touch the calibrated contact radius, which belongs to the
    /// numerical collision model rather than the body's physical geometry.
    pub fn sync_physical_properties(&mut self) {
        self.physical_radius = radius_from_density_mass(self.density, self.mass);
        self.moment_inertia = default_moment_inertia(self.mass, self.physical_radius);
    }

    #[inline]
    pub fn is_diffuse_cloud(&self) -> bool {
        self.material.is_diffuse()
    }
}

/// Default softening before system-scale calibration.
pub fn default_softening(mass: f64) -> f64 {
    EPS_BASE * mass.abs().cbrt()
}

/// Moment of inertia for a uniform sphere: I = (2/5)·m·r².
/// Uses the **physical radius**.
pub fn default_moment_inertia(mass: f64, radius: f64) -> f64 {
    0.4 * mass * radius * radius
}

/// Density from mass and radius: ρ = m / (4/3 π r³).
///
/// Returns a safe positive fallback when `radius ≤ 0`.
pub fn density_from_mass_radius(mass: f64, radius: f64) -> f64 {
    let vol = sphere_volume(radius);
    if vol > 0.0 { mass / vol } else { 1.0 }
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
