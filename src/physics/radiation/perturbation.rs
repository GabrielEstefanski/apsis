//! System-level radiation perturbation.
//!
//! Bridges [`force`](super::force) and the integrator's
//! [`PerturbationForce`] trait. This is the only file in the radiation
//! module that depends on [`Body`] or the integrator interface.
//!
//! # Responsibility boundary
//!
//! - [`Material`] owns `q_pr()` — a scalar property of surface composition.
//! - This module owns the translation from `(Body, Material)` → [`RadiationParams`],
//!   because that requires knowing both the geometry (`physical_radius`) and
//!   the composition (`q_pr`). Neither `Body` nor `Material` alone has enough
//!   information.

use std::f64::consts::PI;

use crate::core::body::Body;
use crate::physics::integrator::PerturbationForce;
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
                x: src_body.x,
                y: src_body.y,
                vx: src_body.vx,
                vy: src_body.vy,
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

// ── PerturbationForce impl ────────────────────────────────────────────────────

impl PerturbationForce for RadiationField {
    /// Accumulates radiation accelerations for the full body slice (`offset = 0`).
    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)]) {
        self.accumulate_offset(bodies, scratch_acc, 0);
    }

    /// Accumulates radiation accelerations for a sub-slice of bodies.
    ///
    /// `offset` is the global index of `bodies[0]` within `System::bodies`.
    /// Used by [`System::apply_perturbations_planets`] during the
    /// Wisdom–Holman sub-step, where `scratch_acc` covers only `bodies[1..]`
    /// and the global index of each planet is `local_index + 1`.
    fn accumulate_offset(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)], offset: usize) {
        for (local_i, (body, acc)) in bodies.iter().zip(scratch_acc.iter_mut()).enumerate() {
            let global_i = local_i + offset;
            let Some(params) = self.body_params.get(global_i).and_then(|p| p.as_ref()) else {
                continue;
            };

            let pos = [body.x, body.y];
            let vel = [body.vx, body.vy];

            let [ax, ay] = if self.include_pr_drag {
                pr_drag_acceleration(pos, vel, params, &self.source)
            } else {
                radiation_acceleration(pos, params, &self.source)
            };

            acc.0 += ax;
            acc.1 += ay;
        }
    }
}

// ── Private helper ────────────────────────────────────────────────────────────

/// Derives [`RadiationParams`] for a body from its geometry and material.
///
/// Returns `None` if the material is not a radiation receiver (`q_pr == 0`),
/// which covers stars, planets, and all massive bodies. This is the single
/// site where [`Material::q_pr`] and [`Body::physical_radius`] are combined —
/// neither type needs to know about the other.
fn body_radiation_params(body: &Body) -> Option<RadiationParams> {
    let q_pr = body.material.q_pr();
    if q_pr <= 0.0 {
        return None;
    }
    Some(RadiationParams {
        area: PI * body.physical_radius * body.physical_radius,
        mass: body.mass,
        q_pr,
    })
}
