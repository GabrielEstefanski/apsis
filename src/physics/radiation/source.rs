/// A radiating point source (star or compact object).
///
/// Carries only intrinsic source properties. The gravitational parameter
/// `μ = G M` is deliberately excluded: it belongs to the gravitational
/// subsystem, not to the radiation model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiationSource {
    /// Position in the inertial frame.
    pub x: f64,
    pub y: f64,
    /// Velocity in the inertial frame.
    ///
    /// Required for Poynting–Robertson drag: the aberrated flux direction
    /// depends on the relative velocity between source and particle.
    pub vx: f64,
    pub vy: f64,
    /// Bolometric luminosity in internal units (energy · time⁻¹).
    pub luminosity: f64,
    /// Speed of light in internal units.
    ///
    /// Carrying `c` explicitly decouples the module from any particular
    /// choice of unit system.
    pub c: f64,
}
