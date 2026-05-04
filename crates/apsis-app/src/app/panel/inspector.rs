//! Inspector right-panel content — adapter from live `SimulationApp`
//! state into [`InspectorData`] consumed by [`crate::app::inspector::view`].
//!
//! The adapter is the only place that knows about both the simulation
//! domain (`Body`, `OrbitalElements`, `System`) and the inspector data
//! shape. The view itself is domain-agnostic by contract.

use crate::app::inspector::{
    self, ActionData, ActionKind, EnergyData, Header, Identity, InspectorData, KinematicState,
    OrbitData,
};
use crate::app::ui::{SimulationApp, UndoRecord};
use eframe::egui::{self, Color32};

const ACTION_FOCUS: usize = 0;
const ACTION_HIDE_TRAIL: usize = 1;
const ACTION_DELETE: usize = 2;

impl SimulationApp {
    pub(super) fn inspector_content(&mut self, ui: &mut egui::Ui, idx: usize) {
        if self.is_editing_locked() {
            ui.disable();
        }

        let data = self.build_inspector_data(idx);
        let clicked = inspector::show(ui, &data, &mut self.inspector_state);

        match clicked {
            Some(ACTION_FOCUS) => self.toggle_follow_selected(),
            Some(ACTION_HIDE_TRAIL) => self.toggle_hide_trail(idx),
            Some(ACTION_DELETE) => self.delete_selected_body(idx),
            _ => {},
        }
    }

    /// Build the inspector payload for body `idx` from the current
    /// system snapshot. Pure projection — nothing here mutates state or
    /// triggers domain computation beyond what the cached snapshot already
    /// publishes.
    fn build_inspector_data(&self, idx: usize) -> InspectorData {
        let bodies = self.system.bodies();
        let body = bodies[idx];
        let elements = self.system.orbital_elements();
        let elem = elements.get(idx).and_then(|e| *e);

        let raw_name = self.system.name(idx);
        let display_name =
            if raw_name.is_empty() { format!("body #{idx}") } else { raw_name.to_owned() };
        let total = bodies.len();
        let breadcrumb = format!("body {} of {total}", idx + 1);

        let [cr, cg, cb] = body.color;

        let orbit = elem.map(|el| {
            let primary_name = self.system.name(el.primary_idx).to_owned();
            let primary = if primary_name.is_empty() {
                format!("body #{}", el.primary_idx)
            } else {
                primary_name
            };
            OrbitData {
                primary_name: primary,
                semi_major_axis_m: el.a,
                eccentricity: el.e,
                inclination_rad: el.inclination,
                period_s: el.period,
                lon_ascending_node_rad: el.lon_ascending_node,
                argument_of_pericenter_rad: el.omega,
                true_anomaly_rad: el.true_anomaly,
                mean_anomaly_rad: el.mean_anomaly,
                eccentric_anomaly_rad: el.eccentric_anomaly,
                pericenter_m: el.pericenter(),
                apocenter_m: el.apocenter(),
            }
        });

        let energy = build_energy(bodies, idx, elem, self.system.g_factor());

        InspectorData {
            header: Header {
                name: display_name,
                breadcrumb,
                swatch: Color32::from_rgb(cr, cg, cb),
            },
            identity: Identity { mass_kg: body.mass, radius_m: Some(body.physical_radius) },
            state: KinematicState {
                position_m: [body.x, body.y, body.z],
                velocity_m_s: [body.vx, body.vy, body.vz],
            },
            orbit,
            energy: Some(energy),
            // Perturbation enumeration by name is not yet exposed by
            // `System`; the `PerturbationForce` trait carries no human
            // identifier. Section auto-hides while `perturbations` is
            // empty.
            perturbations: Vec::new(),
            // Camera-relative readouts assume a 3D camera with a known
            // world-space position. The current 2D-pan camera does not
            // expose one cleanly; the section auto-hides until the 3D
            // canvas + floating-origin work lands.
            camera_relative: None,
            actions: vec![
                ActionData {
                    label: if self.follow_selected_body {
                        "Unfollow camera".to_owned()
                    } else {
                        "Focus camera".to_owned()
                    },
                    icon: None,
                    shortcut: Some("F".to_owned()),
                    kind: ActionKind::Neutral,
                },
                ActionData {
                    label: "Hide trail".to_owned(),
                    icon: None,
                    shortcut: Some("H".to_owned()),
                    kind: ActionKind::Neutral,
                },
                ActionData {
                    label: "Delete".to_owned(),
                    icon: None,
                    shortcut: Some("Del".to_owned()),
                    kind: ActionKind::Destructive,
                },
            ],
        }
    }

    fn toggle_follow_selected(&mut self) {
        self.follow_selected_body = !self.follow_selected_body;
    }

    fn toggle_hide_trail(&mut self, idx: usize) {
        if let Some(hint) = self.render_hints.get_mut(idx) {
            hint.show_trail = !hint.show_trail;
        }
    }

    fn delete_selected_body(&mut self, idx: usize) {
        let body = self.system.bodies()[idx];
        let name = self.system.name(idx).to_owned();
        let old_last = self.system.bodies().len().saturating_sub(1);
        self.push_undo(UndoRecord::RemovedBody { body, name });
        self.system.remove_body(idx);
        self.pins_after_swap_remove(idx, old_last);
        self.selected_body = None;
        self.follow_selected_body = false;
        self.selection_form = None;
    }
}

/// Per-body kinetic and potential plus the orbital-element specific
/// energy (when defined). Potential is a pairwise sum over all other
/// bodies — O(N) per inspector frame, fine for the regimes the inspector
/// is used in (selection of one body in a system of < 10³).
fn build_energy(
    bodies: &[apsis::domain::body::Body],
    idx: usize,
    orbit: Option<apsis::physics::orbital::OrbitalElements>,
    g: f64,
) -> EnergyData {
    let body = bodies[idx];
    let v2 = body.vx * body.vx + body.vy * body.vy + body.vz * body.vz;
    let kinetic = 0.5 * body.mass * v2;

    let mut potential = 0.0;
    for (j, other) in bodies.iter().enumerate() {
        if j == idx {
            continue;
        }
        let dx = body.x - other.x;
        let dy = body.y - other.y;
        let dz = body.z - other.z;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r > 1e-15 {
            potential -= g * body.mass * other.mass / r;
        }
    }

    let specific = orbit.map(|el| el.energy).unwrap_or(f64::NAN);
    EnergyData { kinetic_j: kinetic, potential_j: potential, specific_j: specific }
}
