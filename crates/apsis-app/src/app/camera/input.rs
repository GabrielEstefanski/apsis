//! Gesture → camera-target translation.
//!
//! This module is deliberately egui-free. The canvas adapter (CP3)
//! reads `egui::Response` and `egui::Input`, packs the relevant
//! values into [`DragInput`] / scroll scalars, and calls the
//! `apply_*` entry points here. Keeping the translator pure makes
//! sensitivity, sign conventions, and modifier mappings unit-testable
//! without spinning up a windowing context.
//!
//! # Sign conventions
//!
//! Following the modeller idiom (Blender, Maya, Universe Sandbox):
//!
//! - **Pan**: scene follows the cursor. Dragging right slides the
//!   world right under the cursor — the pivot moves left in world
//!   space.
//! - **Orbit**: eye follows the cursor. Dragging right orbits the
//!   eye toward the scene's right; the scene visually rotates left.
//! - **Zoom**: scrolling up zooms in (distance shrinks).
//!
//! Pan magnitude scales with `distance` so a one-pixel drag covers
//! the same fraction of the visible frame at every zoom level —
//! standard professional behaviour.

use super::OrbitCamera;
use glam::DVec2;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PointerButton {
    #[default]
    Primary,
    Secondary,
    Middle,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Modifiers {
    pub shift: bool,
    pub alt: bool,
    pub ctrl: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct DragInput {
    /// Pixel delta since the last frame, in egui screen coordinates
    /// (y grows downward).
    pub delta_px: DVec2,
    pub button: PointerButton,
    pub modifiers: Modifiers,
}

#[derive(Debug, Clone, Copy)]
pub struct CameraInputConfig {
    /// Radians of azimuth/elevation per pixel of orbit drag.
    pub rotate_rad_per_px: f64,
    /// Pan sensitivity expressed as a fraction of `distance` per
    /// pixel. Zoom-aware: a 1-px drag pans the same on-screen
    /// fraction at any distance.
    pub pan_per_px_per_distance: f64,
    /// Exponential zoom rate. `factor = exp(-amount · rate)`, so
    /// scrolling by `amount = 1` multiplies distance by `e^-rate`.
    pub zoom_rate: f64,
}

impl Default for CameraInputConfig {
    fn default() -> Self {
        // Tuned for a 16-px-tall scrollwheel notch on a 1080p display:
        // one notch ≈ 12% distance change, a full screen-width orbit
        // drag ≈ a quarter turn.
        Self { rotate_rad_per_px: 0.005, pan_per_px_per_distance: 0.0015, zoom_rate: 0.12 }
    }
}

/// Resolved gesture intent. The orbit camera does not care which
/// button or modifier produced the gesture, only what to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GestureKind {
    Orbit,
    Pan,
}

fn classify(input: DragInput) -> GestureKind {
    match input.button {
        PointerButton::Middle => GestureKind::Pan,
        PointerButton::Primary if input.modifiers.shift => GestureKind::Pan,
        PointerButton::Primary => GestureKind::Orbit,
        PointerButton::Secondary => GestureKind::Orbit,
    }
}

pub fn apply_drag(camera: &mut OrbitCamera, input: DragInput, config: &CameraInputConfig) {
    match classify(input) {
        GestureKind::Orbit => {
            // Eye-follows-cursor: drag right (+Δx) increases azimuth so
            // the eye orbits toward the scene's right; drag down (+Δy in
            // egui's y-down screen) increases elevation so the scene
            // tilts forward.
            camera.rotate(
                input.delta_px.x * config.rotate_rad_per_px,
                input.delta_px.y * config.rotate_rad_per_px,
            );
        },
        GestureKind::Pan => {
            // Scene-follows-cursor: drag right slides scene right, so the
            // pivot moves left along the camera's right vector. The
            // y-axis sign also inverts because egui screen-y is down
            // while our `up` is world up.
            let world_per_px = camera.target.distance * config.pan_per_px_per_distance;
            camera.pan_pivot(-input.delta_px.x * world_per_px, -input.delta_px.y * world_per_px);
        },
    }
}

/// `amount` is a unitless scroll measure: positive = wheel rolled up
/// (zoom in), negative = wheel rolled down (zoom out). The canvas
/// adapter normalises egui's pixel-valued `smooth_scroll_delta.y`
/// into ticks before calling.
pub fn apply_scroll(camera: &mut OrbitCamera, amount: f64, config: &CameraInputConfig) {
    if amount == 0.0 {
        return;
    }
    let factor = (-amount * config.zoom_rate).exp();
    camera.zoom(factor);
}

#[cfg(test)]
mod tests {
    use super::super::CameraPose;
    use super::*;
    use glam::DVec3;

    fn camera_at(distance: f64) -> OrbitCamera {
        OrbitCamera::new(CameraPose::new(DVec3::ZERO, 0.0, 0.0, distance))
    }

    #[test]
    fn scroll_zero_is_noop() {
        let mut cam = camera_at(10.0);
        apply_scroll(&mut cam, 0.0, &CameraInputConfig::default());
        assert_eq!(cam.target.distance, 10.0);
    }

    #[test]
    fn scroll_up_zooms_in() {
        let mut cam = camera_at(10.0);
        apply_scroll(&mut cam, 1.0, &CameraInputConfig::default());
        assert!(cam.target.distance < 10.0);
    }

    #[test]
    fn scroll_down_zooms_out() {
        let mut cam = camera_at(10.0);
        apply_scroll(&mut cam, -1.0, &CameraInputConfig::default());
        assert!(cam.target.distance > 10.0);
    }

    #[test]
    fn scroll_up_then_down_returns_to_start() {
        let mut cam = camera_at(10.0);
        let cfg = CameraInputConfig::default();
        apply_scroll(&mut cam, 1.5, &cfg);
        apply_scroll(&mut cam, -1.5, &cfg);
        assert!((cam.target.distance - 10.0).abs() < 1e-12);
    }

    #[test]
    fn primary_drag_orbits() {
        let mut cam = camera_at(10.0);
        let initial = cam.target;
        let drag = DragInput {
            delta_px: DVec2::new(50.0, 30.0),
            button: PointerButton::Primary,
            modifiers: Modifiers::default(),
        };
        apply_drag(&mut cam, drag, &CameraInputConfig::default());
        assert!(cam.target.azimuth > initial.azimuth);
        assert!(cam.target.elevation > initial.elevation);
        assert_eq!(cam.target.pivot, initial.pivot);
    }

    #[test]
    fn middle_drag_pans() {
        let mut cam = camera_at(10.0);
        let initial = cam.target;
        let drag = DragInput {
            delta_px: DVec2::new(50.0, 0.0),
            button: PointerButton::Middle,
            modifiers: Modifiers::default(),
        };
        apply_drag(&mut cam, drag, &CameraInputConfig::default());
        assert_eq!(cam.target.azimuth, initial.azimuth);
        assert_eq!(cam.target.elevation, initial.elevation);
        // Drag right → pivot moves left along world +X.
        assert!(cam.target.pivot.x < 0.0);
    }

    #[test]
    fn shift_primary_drag_pans() {
        let mut cam = camera_at(10.0);
        let drag = DragInput {
            delta_px: DVec2::new(50.0, 0.0),
            button: PointerButton::Primary,
            modifiers: Modifiers { shift: true, ..Modifiers::default() },
        };
        apply_drag(&mut cam, drag, &CameraInputConfig::default());
        assert!(cam.target.pivot.x < 0.0);
    }

    #[test]
    fn pan_scales_with_distance() {
        let cfg = CameraInputConfig::default();
        let drag = DragInput {
            delta_px: DVec2::new(10.0, 0.0),
            button: PointerButton::Middle,
            modifiers: Modifiers::default(),
        };

        let mut near = camera_at(1.0);
        apply_drag(&mut near, drag, &cfg);

        let mut far = camera_at(100.0);
        apply_drag(&mut far, drag, &cfg);

        let near_pan = near.target.pivot.length();
        let far_pan = far.target.pivot.length();
        assert!((far_pan / near_pan - 100.0).abs() < 1e-9);
    }

    #[test]
    fn pan_y_sign_matches_screen_y_down() {
        let mut cam = camera_at(10.0);
        let drag = DragInput {
            delta_px: DVec2::new(0.0, 50.0),
            button: PointerButton::Middle,
            modifiers: Modifiers::default(),
        };
        apply_drag(&mut cam, drag, &CameraInputConfig::default());
        // Drag down on a y-down screen → scene slides down → pivot
        // moves down along world up = pivot.y < 0.
        assert!(cam.target.pivot.y < 0.0);
    }

    #[test]
    fn zero_drag_is_noop() {
        let mut cam = camera_at(10.0);
        let initial = cam.target;
        let drag = DragInput {
            delta_px: DVec2::ZERO,
            button: PointerButton::Primary,
            modifiers: Modifiers::default(),
        };
        apply_drag(&mut cam, drag, &CameraInputConfig::default());
        assert_eq!(cam.target, initial);
    }
}
