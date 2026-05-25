//! System-level radiation perturbation.
//!
//! Bridges [`force`](super::force) and the integrator's operator surface.
//! This is the only file in the radiation module that depends on
//! [`Body`] or the integrator interface.
//!
//! # Responsibility boundary
//!
//! - [`Body`] carries `q_pr` directly as a physical property (set by the
//!   construction preset for receiver classes — asteroids, comets, icy
//!   bodies — and zero on emitters and large planets).
//! - This module reads `(Body::q_pr, Body::physical_radius, Body::mass)`
//!   and packs them into [`RadiationParams`] for the force kernels.
//!
//! # Operator classification
//!
//! `RadiationField` registers as
//! [`NonConservativeOperator`](crate::physics::integrator::NonConservativeOperator)
//! to accommodate the Poynting–Robertson drag branch
//! (`include_pr_drag = true`). The pure radiation-pressure branch is
//! Hamiltonian (`V_rad = −F_rad · x`); a future split into
//! `RadiationPressure` + `PoyntingRobertsonDrag` would let symplectic
//! integrators preserve invariants in the pressure-only configuration.
//! Tracking issue: split deferred until a downstream consumer asks for
//! the conservation guarantee.

use std::f64::consts::PI;

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::{NonConservativeOperator, Operator};
use crate::physics::radiation::force::{pr_drag_acceleration, radiation_acceleration};
use crate::physics::radiation::params::RadiationParams;
use crate::physics::radiation::source::RadiationSource;

// ── RadiationField ────────────────────────────────────────────────────────────

/// Radiation pressure and optional Poynting–Robertson drag from a single
/// radiating source, applied to an arbitrary number of bodies.
///
/// Each body carries its own [`RadiationParams`] (area, mass, Q_pr).
/// The source is shared across all bodies and evaluated once per call.
///
/// Bodies with `None` in `body_params` are silently skipped — zero force,
/// zero cost. This is the mechanism by which massive bodies (stars, giant
/// planets) are exempted without special-casing the integrator.
pub struct RadiationField {
    /// The radiating source (star or compact object).
    pub source: RadiationSource,
    /// Per-body radiation parameters, indexed by **global** body index,
    /// parallel to `System::bodies`. `None` entries are treated as inert.
    pub body_params: Vec<Option<RadiationParams>>,
    /// If `true`, includes the Poynting–Robertson drag term alongside
    /// direct radiation pressure.
    pub include_pr_drag: bool,
}

impl RadiationField {
    /// Constructs an empty field for `n_bodies` bodies.
    ///
    /// All slots are initialised to `None` (inert). Use [`set_params`] to
    /// assign radiation parameters to individual bodies manually, or
    /// [`from_bodies`] to populate automatically from materials.
    pub fn new(source: RadiationSource, n_bodies: usize, include_pr_drag: bool) -> Self {
        Self { source, body_params: vec![None; n_bodies], include_pr_drag }
    }

    /// Assigns radiation parameters to the body at `index`.
    /// Pass `None` to make the body inert.
    pub fn set_params(&mut self, index: usize, params: Option<RadiationParams>) {
        if let Some(slot) = self.body_params.get_mut(index) {
            *slot = params;
        }
    }

    /// Constructs one [`RadiationField`] per luminous body found in `bodies`.
    ///
    /// For each body whose `luminosity` field is positive, a field is created
    /// with that body as the source. Receiver params are derived automatically
    /// from each other body's material via [`body_radiation_params`]. A source
    /// is never its own receiver.
    ///
    /// Returns an empty `Vec` if no body emits radiation — registering an
    /// empty vec has zero cost.
    ///
    /// # Precondition
    ///
    /// [`Body::update_luminosity`] must have been called on all luminous bodies
    /// before invoking this function. Bodies with `luminosity == 0.0` are
    /// silently skipped regardless of their material.
    ///
    /// # Unit conversion parameters
    ///
    /// | Parameter         | Meaning                                     |
    /// |-------------------|---------------------------------------------|
    /// | `c`               | speed of light in internal length · time⁻¹  |
    pub fn from_bodies(bodies: &[Body], c: f64, include_pr_drag: bool) -> Vec<Self> {
        let mut fields = Vec::new();

        for (src_idx, src_body) in bodies.iter().enumerate() {
            let luminosity = src_body.luminosity;
            if luminosity <= 0.0 {
                continue;
            }

            let source = RadiationSource {
                x: src_body.pos_x,
                y: src_body.pos_y,
                z: src_body.pos_z,
                vx: src_body.vel_x,
                vy: src_body.vel_y,
                vz: src_body.vel_z,
                luminosity,
                c,
            };

            let body_params = bodies
                .iter()
                .enumerate()
                .map(|(i, b)| if i == src_idx { None } else { body_radiation_params(b) })
                .collect();

            fields.push(Self { source, body_params, include_pr_drag });
        }

        fields
    }
}

// ── Operator impls ────────────────────────────────────────────────────────────

impl Operator for RadiationField {}

impl NonConservativeOperator for RadiationField {
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]) {
        for (i, (body, a_slot)) in bodies.iter().zip(acc.iter_mut()).enumerate() {
            let Some(params) = self.body_params.get(i).and_then(|p| p.as_ref()) else {
                continue;
            };

            let pos = Vec3::new(body.pos_x, body.pos_y, body.pos_z);
            let vel = Vec3::new(body.vel_x, body.vel_y, body.vel_z);

            let a = if self.include_pr_drag {
                pr_drag_acceleration(pos, vel, params, &self.source)
            } else {
                radiation_acceleration(pos, params, &self.source)
            };

            *a_slot += a;
        }
    }
}

// ── Private helper ────────────────────────────────────────────────────────────

/// Derives [`RadiationParams`] for a body from its geometry and the
/// `q_pr` it carries.
///
/// Returns `None` for non-receiver bodies (`q_pr == 0`), which covers
/// stars, planets, and any user body that opted out via
/// [`Body::with_q_pr(0.0)`].
fn body_radiation_params(body: &Body) -> Option<RadiationParams> {
    if body.q_pr <= 0.0 {
        return None;
    }
    Some(RadiationParams {
        area: PI * body.physical_radius * body.physical_radius,
        mass: body.mass,
        q_pr: body.q_pr,
    })
}
