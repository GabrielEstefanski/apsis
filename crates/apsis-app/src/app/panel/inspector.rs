//! Inspector right-panel content — adapter from live `SimulationApp`
//! state into [`InspectorData`] consumed by [`crate::app::inspector::view`].
//!
//! The adapter is the only place that knows about both the simulation
//! domain (`Body`, `OrbitalElements`, `System`) and the inspector data
//! shape. The view itself is domain-agnostic by contract.

use crate::app::inspector::{
    self, ActionData, ActionKind, AggregateData, EnergyData, Header, Identity, InspectorData,
    KinematicState, OrbitData, PerturbationData, RelationKind, RelationsData,
};
use crate::app::ui::{BodySelection, SimulationApp, UndoRecord};
use apsis::physics::orbital::{self as physics_orbital, HierarchicalRelation, is_system_root};
use eframe::egui::{self, Color32};
use std::collections::BTreeSet;

const ACTION_FOCUS: usize = 0;
const ACTION_HIDE_TRAIL: usize = 1;
const ACTION_DELETE: usize = 2;

const ACTION_AGG_DELETE: usize = 0;
const ACTION_AGG_DESELECT: usize = 1;

impl SimulationApp {
    pub(super) fn aggregate_content(&mut self, ui: &mut egui::Ui, indices: &BTreeSet<usize>) {
        if self.is_editing_locked() {
            ui.disable();
        }

        let data = self.build_aggregate_data(indices);
        let clicked = inspector::show_aggregate(ui, &data);

        match clicked {
            Some(ACTION_AGG_DELETE) => self.delete_selected_bodies(indices.clone()),
            Some(ACTION_AGG_DESELECT) => {
                self.selection = BodySelection::default();
                self.selection_form = None;
            },
            _ => {},
        }
    }

    fn build_aggregate_data(&self, indices: &BTreeSet<usize>) -> AggregateData {
        let bodies = self.system.bodies();

        let mut total_mass = 0.0_f64;
        let mut com = [0.0_f64; 3];
        let mut v_com = [0.0_f64; 3];

        for &i in indices {
            let b = bodies[i];
            total_mass += b.mass;
            com[0] += b.mass * b.x;
            com[1] += b.mass * b.y;
            com[2] += b.mass * b.z;
            v_com[0] += b.mass * b.vx;
            v_com[1] += b.mass * b.vy;
            v_com[2] += b.mass * b.vz;
        }

        if total_mass > 0.0 {
            com[0] /= total_mass;
            com[1] /= total_mass;
            com[2] /= total_mass;
            v_com[0] /= total_mass;
            v_com[1] /= total_mass;
            v_com[2] /= total_mass;
        }

        // bounding_radius = max(|r_i − COM|)
        let bounding_radius = indices
            .iter()
            .map(|&i| {
                let b = bodies[i];
                let dx = b.x - com[0];
                let dy = b.y - com[1];
                let dz = b.z - com[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .fold(0.0_f64, f64::max);

        let body_names: Vec<String> = indices
            .iter()
            .map(|&i| {
                let n = self.system.name(i);
                if n.is_empty() { format!("body #{i}") } else { n.to_owned() }
            })
            .collect();

        AggregateData {
            count: indices.len(),
            body_names,
            total_mass_kg: total_mass,
            com_m: com,
            v_com_m_s: v_com,
            bounding_radius_m: bounding_radius,
            actions: vec![
                ActionData {
                    label: "Delete selected".to_owned(),
                    icon: None,
                    shortcut: Some("Del".to_owned()),
                    kind: ActionKind::Destructive,
                },
                ActionData {
                    label: "Deselect all".to_owned(),
                    icon: None,
                    shortcut: Some("Esc".to_owned()),
                    kind: ActionKind::Neutral,
                },
            ],
        }
    }

    /// Delete all bodies in `indices` in descending index order so that
    /// swap-remove semantics don't invalidate the remaining indices.
    fn delete_selected_bodies(&mut self, indices: BTreeSet<usize>) {
        // BTreeSet iterates ascending; reverse gives descending.
        let sorted: Vec<usize> = indices.into_iter().rev().collect();
        for idx in sorted {
            if idx < self.system.bodies().len() {
                let body = self.system.bodies()[idx];
                let name = self.system.name(idx).to_owned();
                let old_last = self.system.bodies().len().saturating_sub(1);
                self.push_undo(UndoRecord::RemovedBody { body, name });
                self.system.remove_body(idx);
                self.pins_after_swap_remove(idx, old_last);
            }
        }
        self.selection = BodySelection::default();
        self.follow_selected_body = false;
        self.selection_form = None;
    }

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
        // System root has no Keplerian orbit; rendering one would
        // misrepresent N-body dynamics. Suppress both ORBIT and RELATIONS
        // sections — there is no meaningful primary to report.
        let elem =
            if is_system_root(bodies, idx) { None } else { elements.get(idx).and_then(|e| *e) };

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
        let relations = build_relations(bodies, idx, elem, &display_name, |i| {
            let n = self.system.name(i);
            if n.is_empty() { format!("body #{i}") } else { n.to_owned() }
        });

        InspectorData {
            header: Header {
                name: display_name.clone(),
                breadcrumb,
                swatch: Color32::from_rgb(cr, cg, cb),
            },
            identity: Identity { mass_kg: body.mass, radius_m: Some(body.physical_radius) },
            state: KinematicState {
                position_m: [body.x, body.y, body.z],
                velocity_m_s: [body.vx, body.vy, body.vz],
            },
            orbit,
            relations,
            energy: Some(energy),
            perturbations: self
                .perturbation_catalog
                .iter()
                .filter(|e| e.enabled)
                .map(|e| PerturbationData {
                    name: e.descriptor.name().to_owned(),
                    active: true,
                    readouts: Vec::new(),
                })
                .collect(),
            // Needs 3D camera world-pose; section auto-hides until the canvas lands.
            camera_relative: None,
            actions: vec![
                ActionData {
                    label: if self.follow_selected_body {
                        "Unfollow camera"
                    } else {
                        "Focus camera"
                    }
                    .to_owned(),
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
        self.selection = BodySelection::default();
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

/// Compute the hierarchical relationship for the selected body. Returns
/// `Some(RelationsData)` only when the hierarchical primary differs from
/// the strongest attractor used by `orbit` — the section auto-hides
/// otherwise so the inspector stays compact for non-hierarchical cases.
///
/// `name_of` resolves a body index to a human-readable name, falling back
/// to `body #N` when the simulation didn't supply one.
fn build_relations(
    bodies: &[apsis::domain::body::Body],
    idx: usize,
    orbit: Option<apsis::physics::orbital::OrbitalElements>,
    secondary_name: &str,
    name_of: impl Fn(usize) -> String,
) -> Option<RelationsData> {
    let (primary_idx, kind) = physics_orbital::hierarchical_primary(bodies, idx)?;
    let strongest_idx = orbit.map(|el| el.primary_idx);
    if strongest_idx == Some(primary_idx) {
        return None;
    }
    let primary_name = name_of(primary_idx);
    let frame_label = format!("Barycentric ({primary_name}–{secondary_name})");
    Some(RelationsData {
        primary_name,
        secondary_name: secondary_name.to_owned(),
        kind: match kind {
            HierarchicalRelation::HillSphere => RelationKind::BoundHillSphere,
            HierarchicalRelation::Energy => RelationKind::BoundEnergy,
        },
        frame_label,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use apsis::domain::body::Body;

    fn make_body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        let mut b = Body::rocky(mass).at(x, y).with_velocity(vx, vy);
        b.z = 0.0;
        b.vz = 0.0;
        b
    }

    fn aggregate(bodies: &[Body], indices: &[usize]) -> AggregateData {
        let mut total_mass = 0.0_f64;
        let mut com = [0.0_f64; 3];
        let mut v_com = [0.0_f64; 3];

        for &i in indices {
            let b = bodies[i];
            total_mass += b.mass;
            com[0] += b.mass * b.x;
            com[1] += b.mass * b.y;
            com[2] += b.mass * b.z;
            v_com[0] += b.mass * b.vx;
            v_com[1] += b.mass * b.vy;
            v_com[2] += b.mass * b.vz;
        }

        if total_mass > 0.0 {
            com[0] /= total_mass;
            com[1] /= total_mass;
            com[2] /= total_mass;
            v_com[0] /= total_mass;
            v_com[1] /= total_mass;
            v_com[2] /= total_mass;
        }

        let bounding_radius = indices
            .iter()
            .map(|&i| {
                let b = bodies[i];
                let dx = b.x - com[0];
                let dy = b.y - com[1];
                let dz = b.z - com[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .fold(0.0_f64, f64::max);

        AggregateData {
            count: indices.len(),
            body_names: indices.iter().map(|i| format!("body #{i}")).collect(),
            total_mass_kg: total_mass,
            com_m: com,
            v_com_m_s: v_com,
            bounding_radius_m: bounding_radius,
            actions: vec![],
        }
    }

    #[test]
    fn aggregate_com_equal_masses() {
        let bodies = vec![make_body(-1.0, 0.0, 0.0, 0.0, 1.0), make_body(1.0, 0.0, 0.0, 0.0, 1.0)];
        let d = aggregate(&bodies, &[0, 1]);
        assert!((d.com_m[0]).abs() < 1e-12, "COM x should be 0 for symmetric pair");
        assert!((d.bounding_radius_m - 1.0).abs() < 1e-12);
    }

    #[test]
    fn aggregate_com_unequal_masses() {
        // Body A: mass 1 at x=0; Body B: mass 3 at x=4 → COM = 3.0
        let bodies = vec![make_body(0.0, 0.0, 0.0, 0.0, 1.0), make_body(4.0, 0.0, 0.0, 0.0, 3.0)];
        let d = aggregate(&bodies, &[0, 1]);
        assert!((d.com_m[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn aggregate_translational_invariance() {
        let shift = 1e11_f64;
        let bodies_a = vec![make_body(1.0, 2.0, 0.0, 0.0, 5.0), make_body(3.0, 4.0, 0.0, 0.0, 5.0)];
        let bodies_b = vec![
            make_body(1.0 + shift, 2.0 + shift, 0.0, 0.0, 5.0),
            make_body(3.0 + shift, 4.0 + shift, 0.0, 0.0, 5.0),
        ];
        let da = aggregate(&bodies_a, &[0, 1]);
        let db = aggregate(&bodies_b, &[0, 1]);
        // Bounding radius is translation-invariant.
        assert!((da.bounding_radius_m - db.bounding_radius_m).abs() < 1e-6);
        // COM shifts by exactly `shift`.
        assert!((db.com_m[0] - da.com_m[0] - shift).abs() < 1e-6);
    }

    #[test]
    fn aggregate_zero_mass_no_panic() {
        // All-zero mass bodies — total_mass = 0; guard prevents division; COM stays 0.
        let bodies = vec![make_body(1.0, 0.0, 0.0, 0.0, 0.0), make_body(-1.0, 0.0, 0.0, 0.0, 0.0)];
        let d = aggregate(&bodies, &[0, 1]);
        assert_eq!(d.total_mass_kg, 0.0);
        assert!(d.com_m[0].is_finite());
    }

    #[test]
    fn body_selection_toggle_invariant() {
        // None → Single → Multi → Single → None
        let sel = BodySelection::default();
        let sel = sel.toggle(0);
        assert!(matches!(sel, BodySelection::Single(0)));

        let sel = sel.toggle(1);
        assert!(matches!(sel, BodySelection::Multi(_)));
        if let BodySelection::Multi(ref s) = sel {
            assert_eq!(s.len(), 2);
        }

        let sel = sel.toggle(1);
        assert!(matches!(sel, BodySelection::Single(0)));

        let sel = sel.toggle(0);
        assert!(matches!(sel, BodySelection::None));
    }

    #[test]
    fn body_selection_multi_never_has_one_element() {
        // Toggling the second body back out must collapse to Single, not Multi({0}).
        let sel = BodySelection::None.toggle(5).toggle(7).toggle(7);
        assert!(matches!(sel, BodySelection::Single(5)));
    }
}
