//! Event payloads emitted after integration.
//!
//! Events are detected **inside** [`System::step`] on the integrated state, but
//! dispatched to hooks **after** the physics advance completes. Pre-computed
//! quantities (reduced mass, pairwise kinetic energy, momenta) save hooks from
//! recomputing invariants and keep the dispatch order deterministic.

/// Pairwise contact/merge candidate detected after a step.
#[derive(Debug, Clone, Copy)]
pub struct CollisionEvent {
    pub i: usize,
    pub j: usize,
    pub t: f64,
    /// Centre-to-centre separation at detection.
    pub separation: f64,
    /// Sum of physical radii at detection.
    pub contact_radius: f64,
    /// Reduced mass: μ = m_i m_j / (m_i + m_j).
    pub reduced_mass: f64,
    /// Kinetic energy in the pair's centre-of-momentum frame.
    pub relative_ke: f64,
    /// Pairwise linear momentum (p_i + p_j).
    pub total_momentum: (f64, f64),
    /// Pairwise angular momentum about the origin (z-component).
    pub total_angular_momentum: f64,
}

/// Body that escaped the system's bound region.
#[derive(Debug, Clone, Copy)]
pub struct EscapeEvent {
    pub body: usize,
    pub t: f64,
    /// Radial distance from COM at detection.
    pub radius: f64,
    /// Speed (‖v‖) at detection.
    pub speed: f64,
    /// Specific orbital energy (v²/2 − GM/r).
    pub specific_energy: f64,
}
