//! Zone orchestration for the main UI.
//!
//! Layout (see `project_ui_redesign.md` in memory for the full spec):
//!
//! ```text
//! ┌ top bar (session) ──────────────────────────────────────────────┐
//! ├──┬──────────────┬──────────────────────────────┬────────────────┤
//! │🔍│              │                              │                │
//! │➕│  contextual  │         canvas (wgpu)        │   inspector    │
//! │★ │  panel       │                              │   (auto-show)  │
//! │👁│  (280 px)    │                              │   (320 px)     │
//! │⎆ │              │                              │                │
//! │⚙ │              │                              │                │
//! ├──┴──────────────┴──────────────────────────────┴────────────────┤
//! │ playbar (play/pause · dt · SPF · integrator · reset)            │
//! └─────────────────────────────────────────────────────────────────┘
//!   48 px   280 px                                   320 px
//! ```
//!
//! Panel registration order matters in egui: each panel carves space out of
//! the remaining rect. Top → bottom → left-most → next-left → right →
//! central. Keep this order stable so the canvas always occupies the
//! residual rectangle.

mod inspector;
mod metrics;
mod notifications_panel;
mod playbar;
mod precision_confirmation_modal;
mod precision_panel;
mod save_modal;
mod settings_modal;
mod shortcuts_modal;
mod tabs;
mod tool_rail;
mod toolbar;

use crate::app::theme::{BORDER, PANEL_BG};
use crate::app::ui::{BodySelection, SimulationApp};
use eframe::egui::{self, Stroke};
use std::collections::BTreeSet;

const CONTEXTUAL_MIN: f32 = 240.0;
const CONTEXTUAL_DEFAULT: f32 = 300.0;
const INSPECTOR_WIDTH: f32 = 320.0;

impl SimulationApp {
    // ── Top bar (session) ──────────────────────────────────────────────────
    pub(super) fn draw_toolbar(&mut self, ctx: &egui::Context) {
        egui::Panel::top("toolbar")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14, 10))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .show(ctx, |ui| self.toolbar_content(ui));
    }

    // ── Contextual panel (driven by tool rail selection) ───────────────────
    pub(super) fn draw_panel(&mut self, ctx: &egui::Context) {
        egui::Panel::left("contextual")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14, 12))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .default_size(CONTEXTUAL_DEFAULT)
            .min_size(CONTEXTUAL_MIN)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    self.panel_tab_dispatch(ui);
                });
            });
    }

    // ── Inspector (right, auto-show when a body is selected) ───────────────
    pub(super) fn draw_inspector(&mut self, ctx: &egui::Context) {
        if matches!(self.selection, BodySelection::None) {
            self.system.set_orbital_elements_needed(false);
            return;
        }

        // Stale-index guard: prune selected indices that have gone out of range.
        let n_bodies = self.system.bodies().len();
        let stale = match &self.selection {
            BodySelection::Single(i) if *i >= n_bodies => Some(BodySelection::default()),
            BodySelection::Multi(set) => {
                let valid: BTreeSet<usize> =
                    set.iter().copied().filter(|&i| i < n_bodies).collect();
                (valid.len() < set.len()).then(|| match valid.len() {
                    0 => BodySelection::default(),
                    1 => BodySelection::Single(*valid.iter().next().unwrap()),
                    _ => BodySelection::Multi(valid),
                })
            },
            _ => None,
        };
        if let Some(sel) = stale {
            if matches!(sel, BodySelection::None) {
                self.follow_selected_body = false;
                self.selection_form = None;
                self.system.set_orbital_elements_needed(false);
            }
            self.selection = sel;
            if matches!(self.selection, BodySelection::None) {
                return;
            }
        }

        // Extract dispatch info before the closure to avoid borrow conflicts.
        let single_idx = self.selection.single();
        let multi_set: Option<BTreeSet<usize>> = match &self.selection {
            BodySelection::Multi(s) => Some(s.clone()),
            _ => None,
        };

        self.system.set_orbital_elements_needed(single_idx.is_some());

        egui::Panel::right("inspector")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(16, 14))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .default_size(INSPECTOR_WIDTH)
            .min_size(INSPECTOR_WIDTH)
            .max_size(INSPECTOR_WIDTH)
            .resizable(false)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    if let Some(idx) = single_idx {
                        self.inspector_content(ui, idx);
                    } else if let Some(ref indices) = multi_set {
                        self.aggregate_content(ui, indices);
                    }
                });
            });
    }

    /// Apply the freshest pending view: template-supplied hints when
    /// present, bounding-sphere fit otherwise. Both consume
    /// [`SimulationApp::pending_fit`].
    pub(in crate::app) fn apply_pending_view(&mut self) {
        if let Some(hints) = self.pending_camera_hints.take() {
            self.apply_template_view_hints(hints);
        } else {
            self.fit_to_view();
        }
    }

    fn apply_template_view_hints(&mut self, hints: crate::app::ui::TemplateCameraHints) {
        let bodies = self.system.bodies();
        if bodies.is_empty() {
            return;
        }
        let centroid = body_centroid(bodies);
        let distance = hints.distance.unwrap_or_else(|| bounding_sphere_distance(bodies, centroid));
        let up = hints
            .up
            .map(|[x, y, z]| glam::DVec3::new(x, y, z))
            .filter(|v| v.length_squared() > 1e-12)
            .map(|v| v.normalize())
            .unwrap_or(glam::DVec3::Z);

        let pose = camera_pose_from_orbital_up(centroid, distance, up);
        self.camera.target = pose;
        self.camera.current = pose;
        self.orbital_plane_up = up;
        self.follow_selected_body = false;
    }

    // ── Fit view (shared by Camera tool and F shortcut) ────────────────────
    pub(in crate::app) fn fit_to_view(&mut self) {
        let bodies = self.system.bodies();
        if bodies.is_empty() {
            return;
        }
        let centroid = body_centroid(bodies);
        let dist = bounding_sphere_distance(bodies, centroid);
        self.camera.target.pivot = centroid;
        self.camera.target.distance = dist;
        self.follow_selected_body = false;
    }
}

/// Centroid of the AABB enclosing every body. Single-pass, O(N).
fn body_centroid(bodies: &[apsis::domain::body::Body]) -> glam::DVec3 {
    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_z, mut max_z) = (f64::INFINITY, f64::NEG_INFINITY);
    for b in bodies {
        min_x = min_x.min(b.x);
        max_x = max_x.max(b.x);
        min_y = min_y.min(b.y);
        max_y = max_y.max(b.y);
        min_z = min_z.min(b.z);
        max_z = max_z.max(b.z);
    }
    glam::DVec3::new((min_x + max_x) * 0.5, (min_y + max_y) * 0.5, (min_z + max_z) * 0.5)
}

/// Distance such that the half-AABB-extent projects to half the vertical
/// FOV with a 1.2× margin. Floored at 5× the near plane so single-body
/// systems don't clip into the camera.
fn bounding_sphere_distance(bodies: &[apsis::domain::body::Body], centroid: glam::DVec3) -> f64 {
    let mut extent: f64 = 0.0;
    for b in bodies {
        extent = extent
            .max((b.x - centroid.x).abs())
            .max((b.y - centroid.y).abs())
            .max((b.z - centroid.z).abs());
    }
    let extent = extent.max(1e-9);
    let half_fov = (crate::app::camera::FOV_Y_RAD as f64) * 0.5;
    (extent / half_fov.tan() * 1.2).max(crate::app::camera::NEAR_PLANE as f64 * 5.0)
}

/// Build a [`CameraPose`] looking at `pivot` from above the orbital
/// plane (`up` direction in world coords) with a ~28° tilt toward the
/// world-Y axis — matches the NASA-Eyes / Universe-Sandbox default
/// scene-load feel.
fn camera_pose_from_orbital_up(
    pivot: glam::DVec3,
    distance: f64,
    up: glam::DVec3,
) -> crate::app::camera::CameraPose {
    // Tilt direction lives in the plane perpendicular to `up`. Project
    // world-Y into that plane; if `up` is parallel to Y, fall back to
    // world-Z so the camera always lands at a defined pose.
    let proj_y = glam::DVec3::Y - up * up.dot(glam::DVec3::Y);
    let tilt_axis = if proj_y.length() > 1e-6 { proj_y.normalize() } else { glam::DVec3::Z };
    let tilt = 0.5_f64;
    let dir = (up * tilt.cos() + tilt_axis * tilt.sin()).normalize();
    // Camera convention: dir = (cos(el)·sin(az), sin(el), cos(el)·cos(az)).
    let el = dir.y.asin();
    let az = dir.x.atan2(dir.z);
    crate::app::camera::CameraPose::new(pivot, az, el, distance)
}
