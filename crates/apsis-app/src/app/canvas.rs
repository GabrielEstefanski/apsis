use crate::app::camera::{FOV_Y_RAD, NEAR_PLANE};
use crate::app::ui::{BodySelection, SelectionForm, SemanticScaleMode, SimulationApp, UndoRecord};
use crate::render::CallbackFn;
use crate::render::lighting::{LightSpec, SceneLighting};
use crate::render::orbit_overlay::{
    OrbitOverlayStyle, closest_sample_index, draw_orbit_apsides, draw_orbit_polyline,
    draw_orbit_polyline_with_halo,
};
use crate::render::render_relative::RenderRelativeVec3;
use apsis::physics::orbital::{compute_elements, dominant_primary, is_system_root};
use apsis::templates::instantiate_at;
use eframe::egui::{self, Color32, FontId, Pos2, Stroke};
use eframe::egui_wgpu;

// ── Tunables ──────────────────────────────────────────────────────────────────

/// Minimum hit-test radius in pixels. Makes small/distant bodies easier to click.
const MIN_HIT_PX: f32 = 10.0;

/// Maximum distance in pixels from the body centre to place its name label.
/// Prevents labels drifting far from large bodies at high zoom.
const MAX_LABEL_OFFSET_PX: f32 = 48.0;

/// Label font size.
const LABEL_FONT_SIZE: f32 = 10.5;

/// Selection ring: base gap outside the body disc.
const RING_GAP: f32 = 5.0;

/// Camera pan animation: fraction of remaining distance applied each frame.
const CAM_LERP: f32 = 0.16;

/// Body framing: target on-screen pixel size for a body when click-to-focus
/// reframes the camera. The chosen distance puts `physical_radius` at
/// roughly this many pixels; mode-specific min-px floors take over when the
/// body would otherwise be sub-pixel.
const FRAME_TARGET_PX: f32 = 80.0;

/// Eccentricity window around 1.0 within which an orbit is considered
/// numerically degenerate (near-parabolic). Chosen deliberately tight so
/// high-eccentricity cometary orbits (e ≈ 0.95 – 0.99) remain visible —
/// those are exactly the orbits worth looking at.
const DEGEN_ECC_WINDOW: f64 = 0.005;

// ── Camera-triad helpers ──────────────────────────────────────────────────────

/// Project a camera-space direction onto unit screen coordinates.
/// Returns `(dx, dy, depth, colour, label)` where `dy` flips because
/// egui screen-y points down, and positive `depth` means behind the
/// camera (right-handed view space).
fn project(
    cam: glam::DVec3,
    color: eframe::egui::Color32,
    label: &'static str,
) -> (f32, f32, f32, eframe::egui::Color32, &'static str) {
    (cam.x as f32, -cam.y as f32, cam.z as f32, color, label)
}

// ── World → screen projection ────────────────────────────────────────────────

/// Render-frame projection: `projection × rotation_only_view`. Consumed
/// by shaders and CPU helpers whose geometry has already been shifted
/// by the current `render_origin`, so the camera sits at `(0,0,0)` and
/// only the view rotation matters.
fn camera_view_proj_relative(
    camera: &crate::app::camera::OrbitCamera,
    rect: egui::Rect,
) -> glam::Mat4 {
    let aspect = (rect.width() / rect.height().max(1.0)).max(0.001);
    let proj = glam::Mat4::perspective_infinite_reverse_rh(FOV_Y_RAD, aspect, NEAR_PLANE);
    let view = camera.current.view_rotation_only().as_mat4();
    proj * view
}

/// Ray-cast through the camera onto the absolute-world ecliptic plane
/// (`z = 0`). Operates entirely in the render frame (camera at the
/// origin, plane at `z = -render_origin.z`); the absolute hit point is
/// recovered by adding `render_origin` back at `f64` precision.
///
/// Returns `None` when the camera ray doesn't hit that plane in front
/// of the eye — looking up, grazing, or with the pivot above the
/// horizon.
fn screen_to_world_on_z_plane(
    screen_pos: egui::Pos2,
    view_proj_relative: glam::Mat4,
    rect: egui::Rect,
    render_origin: glam::DVec3,
) -> Option<glam::DVec3> {
    let inv = view_proj_relative.inverse();
    let ndc_x = ((screen_pos.x - rect.min.x) / rect.width()) * 2.0 - 1.0;
    let ndc_y = -(((screen_pos.y - rect.min.y) / rect.height()) * 2.0 - 1.0);
    // Reverse-Z infinite-far: clip z = 0 lands on the far plane, so
    // unprojecting that gives a stable far-ray endpoint in the render
    // frame. Camera origin is `(0,0,0)` here, so the ray direction is
    // just the normalised endpoint.
    let far_clip = glam::Vec4::new(ndc_x, ndc_y, 0.0, 1.0);
    let far_relative = inv * far_clip;
    if far_relative.w.abs() < 1e-12 {
        return None;
    }
    let far_pos = far_relative.truncate() / far_relative.w;
    let ray_dir = far_pos.normalize();
    if ray_dir.z.abs() < 1e-6 {
        return None;
    }
    // Absolute `z = 0` plane sits at `z_relative = -render_origin.z`
    // in the render frame.
    let target_z_relative = -render_origin.z as f32;
    let t = target_z_relative / ray_dir.z;
    if t <= 0.0 {
        return None;
    }
    let hit_relative = ray_dir * t;
    // Recover absolute world coordinates at `f64` precision.
    let hit_x = render_origin.x + hit_relative.x as f64;
    let hit_y = render_origin.y + hit_relative.y as f64;
    Some(glam::DVec3::new(hit_x, hit_y, 0.0))
}

/// Project an absolute-world point onto canvas screen coordinates.
/// The world position is shifted into the render frame at `f64`
/// precision (via [`RenderRelativeVec3::from_world`]) before the `f32`
/// cast, so the projection stays stable for bodies at AU-scale absolute
/// positions even when the camera is close to them.
///
/// Returns `None` when the point is at or behind the near plane
/// (clip-space `w ≤ 0`).
fn world_to_screen(
    world_pos: glam::DVec3,
    render_origin: glam::DVec3,
    view_proj_relative: glam::Mat4,
    rect: egui::Rect,
) -> Option<egui::Pos2> {
    let rel = RenderRelativeVec3::from_world(world_pos, render_origin);
    let clip = view_proj_relative * rel.as_vec3().extend(1.0);
    if clip.w <= 0.0 {
        return None;
    }
    let ndc = clip.truncate() / clip.w;
    let sx = rect.min.x + (ndc.x * 0.5 + 0.5) * rect.width();
    // Egui screen-y grows downward; flip the y-NDC sign.
    let sy = rect.min.y + (-ndc.y * 0.5 + 0.5) * rect.height();
    Some(egui::pos2(sx, sy))
}

/// Pixel radius of a sphere with physical size `radius_world` whose
/// centre sits at `center_relative` in the render frame. Matches the
/// formula used by the body shader so hit-test and ring overlays land
/// on the rendered silhouette.
///
/// The camera is at the render-frame origin, so `center_relative.length()`
/// is the camera-to-body distance — small magnitude when the camera is
/// close, with full `f32` precision because both sides of the original
/// subtraction were taken in `f64`.
fn projected_radius_px(
    center_relative: RenderRelativeVec3,
    radius_world: f32,
    canvas_height_px: f32,
) -> f32 {
    let dist = center_relative.as_vec3().length().max(1e-6);
    let focal_y = canvas_height_px / (2.0 * (FOV_Y_RAD * 0.5).tan());
    radius_world * focal_y / dist
}

/// World-space render radius for a body under the active semantic scale mode.
///
/// `Comparative` floors at `min_px`-equivalent world units; `Illustrative`
/// compresses physical size to a 2.5–20 px window so a Solar-System view
/// keeps both Mercury and the Sun on screen without one becoming a dot or
/// the other a wall. `Physical` returns the raw radius — visibility is the
/// caller's problem at that point.
fn radius_world_3d(
    physical_radius: f64,
    mode: SemanticScaleMode,
    min_px: f32,
    camera_dist: f32,
    focal_y: f32,
) -> f32 {
    let r = physical_radius as f32;
    let px_to_world = camera_dist / focal_y;

    match mode {
        SemanticScaleMode::Physical => r,

        SemanticScaleMode::Comparative => r.max(min_px * px_to_world),

        SemanticScaleMode::Illustrative => {
            let physical_px = r * focal_y / camera_dist;
            let k = 0.15;
            let scaled_px = (1.0 - (-k * physical_px).exp()) * 20.0;
            scaled_px.max(min_px).max(2.5) * px_to_world
        },
    }
}

// ── Orbit overlay filter helpers ──────────────────────────────────────────────

/// Maps a hierarchy level to an index into the user's per-level toggle
/// array. Levels 0..=2 map directly; level 3 and deeper fold into slot 3
/// ("L3+") so the UI has a finite, fixed set of toggles.
fn level_is_visible(level: u8, toggles: &[bool; 4]) -> bool {
    let slot = (level as usize).min(3);
    toggles[slot]
}

/// Returns `true` when the orbit is *geometrically* degenerate enough to
/// hide even though `sample_orbit` would still draw it. Two cases:
///
/// * Near-parabolic (|1 − e| < DEGEN_ECC_WINDOW). The curve is numerically
///   unstable and visually indistinguishable from a straight line at
///   typical zoom levels.
/// * Periapsis inside the primary's physical radius — the body would
///   intersect its own primary. Physically impossible as a sustained
///   orbit; usually indicates a glitched state that should not be
///   advertised as a prediction.
fn is_degenerate_orbit(
    el: &apsis::physics::orbital::OrbitalElements,
    primary: &apsis::domain::body::Body,
) -> bool {
    if (el.e - 1.0).abs() < DEGEN_ECC_WINDOW {
        return true;
    }
    if !el.a.is_finite() {
        return true; // parabolic (a = ∞) — sample_orbit returns empty
    }
    // Elliptical: a > 0, r_peri = a(1−e). Hyperbolic: a < 0 by convention,
    // r_peri = |a|(e−1). Both collapse to a*(1−e) with the right sign.
    let r_peri = if el.a > 0.0 { el.a * (1.0 - el.e) } else { el.a.abs() * (el.e - 1.0) };
    r_peri < primary.physical_radius
}

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ui: &mut egui::Ui) {
        let rect = ui.max_rect();
        let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
        let ctx = ui.ctx();

        let (scroll_y, hover_pos, dt, time) = ctx.input(|i| {
            (
                i.smooth_scroll_delta.y,
                i.pointer.hover_pos(),
                i.stable_dt.min(0.05_f32),
                i.time as f32,
            )
        });

        // ── Camera pan animation (spring toward target) ───────────────────────
        if let Some(target) = self.camera_anim_target {
            let delta = target - self.offset;
            if delta.length_sq() < 0.1 {
                self.offset = target;
                self.camera_anim_target = None;
            } else {
                self.offset += delta * CAM_LERP;
                ctx.request_repaint();
            }
        }

        // centre is derived AFTER any animation update
        let center = rect.center() + self.offset;

        // ── 3D camera gestures ────────────────────────────────────────────────
        {
            use crate::app::camera::input::{
                DragInput, Modifiers as CamMods, PointerButton as CamBtn, apply_drag, apply_scroll,
            };
            use glam::DVec2;

            let pointer_in_canvas = hover_pos.is_some_and(|p| rect.contains(p));
            let (rmb, mmb, ptr_delta, mods) = ctx.input(|i| {
                (
                    i.pointer.button_down(egui::PointerButton::Secondary),
                    i.pointer.button_down(egui::PointerButton::Middle),
                    i.pointer.delta(),
                    CamMods {
                        shift: i.modifiers.shift,
                        alt: i.modifiers.alt,
                        ctrl: i.modifiers.ctrl,
                    },
                )
            });

            if pointer_in_canvas && (rmb || mmb) && ptr_delta != egui::Vec2::ZERO {
                // Pan moves the pivot in world space — incompatible with the
                // body-follow loop, which clobbers the pivot each frame. Drop
                // follow on pan so the user's manual drag is what they see.
                // Orbit (RMB) and zoom keep follow alive, matching the
                // Universe Sandbox / KSP map-view idiom.
                if mmb {
                    self.follow_selected_body = false;
                }
                // Any user drag overrides an in-flight cinematic
                // transition — gesture intent wins immediately.
                self.follow_transition = None;
                apply_drag(
                    &mut self.camera,
                    DragInput {
                        delta_px: DVec2::new(ptr_delta.x as f64, ptr_delta.y as f64),
                        button: if mmb { CamBtn::Middle } else { CamBtn::Secondary },
                        modifiers: mods,
                    },
                    &self.camera_input_config,
                );
            }

            if pointer_in_canvas && scroll_y.abs() > 0.0 {
                self.follow_transition = None;
                // smooth_scroll_delta is in pixels; normalise to wheel ticks.
                apply_scroll(&mut self.camera, scroll_y as f64 / 60.0, &self.camera_input_config);
            }
        }

        // ── Pan with inertia ──────────────────────────────────────────────────
        if response.dragged() {
            let delta = response.drag_delta();
            self.offset += delta;
            // Cancel smooth-pan if the user takes manual control
            self.camera_anim_target = None;
            self.follow_selected_body = false;
            self.follow_transition = None;
            let frame_vel = delta / dt.max(1.0 / 120.0);
            self.pan_vel = self.pan_vel * 0.4 + frame_vel * 0.6;
            ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        } else {
            if self.pan_vel.length_sq() > 1.0 {
                self.offset += self.pan_vel * dt;
                self.pan_vel *= (-dt * 2.8_f32).exp();
                if self.pan_vel.length_sq() < 1.0 {
                    self.pan_vel = egui::Vec2::ZERO;
                }
                ctx.request_repaint();
            }
        }

        // ── Body-follow tracking (PD + feedforward) ──────────────────────────
        // Pivot target = body.pos + v/ω + a/ω², so the spring settles
        // on the body for any motion at most quadratic in time. While
        // a transition is active, a body-relative offset decays on
        // top so the new target stays in frame from frame one. Runs
        // after pan/zoom so manual gestures override follow on the
        // same frame.
        if self.follow_selected_body {
            if let Some(idx) = self.selection.single() {
                if let Some(body) = self.system.bodies().get(idx).copied() {
                    let body_pos = glam::DVec3::new(body.x, body.y, body.z);
                    let body_vel = glam::DVec3::new(body.vx, body.vy, body.vz);
                    let body_acc = self
                        .system
                        .accelerations()
                        .get(idx)
                        .copied()
                        .map(|a| glam::DVec3::new(a.x, a.y, a.z))
                        .unwrap_or(glam::DVec3::ZERO);

                    let sim_rate = self.system.sim_rate().max(0.0);
                    let wall_vel = body_vel * sim_rate;
                    let wall_acc = body_acc * sim_rate * sim_rate;

                    let ff = self.camera.feedforward_pivot(dt as f64, body_pos, wall_vel, wall_acc);

                    // Drop any transition pointing at a different body
                    // — happens when the user clicks a new target before
                    // the previous transition settled.
                    if let Some(t) = &self.follow_transition {
                        if t.body_idx != idx {
                            self.follow_transition = None;
                        }
                    }

                    if let Some(state) = self.follow_transition.as_mut() {
                        let settled = state.step(dt as f64);
                        if settled {
                            self.follow_transition = None;
                        } else {
                            // Body-anchored target pose; pivot endpoint
                            // tracks the live body so a fast target stays
                            // in frame throughout the transition.
                            let body_target = crate::app::camera::CameraPose::new(
                                body_pos,
                                state.target_azimuth,
                                state.target_elevation,
                                state.target_distance,
                            );
                            let lerped = state.initial.lerp_to(&body_target, state.t());
                            self.camera.current = lerped;
                            self.camera.target = body_target;
                            ctx.request_repaint();
                        }
                    }

                    if self.follow_transition.is_none() {
                        // Steady state: spring + feedforward pins the
                        // body centred without lag.
                        self.camera.target.pivot = ff;
                        ctx.request_repaint();
                    }
                } else {
                    self.follow_selected_body = false;
                    self.follow_transition = None;
                    self.selection = BodySelection::default();
                    self.selection_form = None;
                }
            } else {
                self.follow_selected_body = false;
                self.follow_transition = None;
            }
        } else if self.follow_transition.is_some() {
            self.follow_transition = None;
        }

        // Spring chase against this frame's `target`. Runs after the
        // gesture and follow blocks write `target`, before matrices
        // and `render_origin` are read below.
        self.camera.integrate(dt as f64);
        if !self.camera.is_at_rest() {
            ctx.request_repaint();
        }

        // World → screen transform for the 3D body pass. Built once per
        // frame and threaded through hit-test, hover, label drawing and
        // the render callback so a single source of truth governs where
        // each body lives on screen.
        let view_proj_relative = camera_view_proj_relative(&self.camera, rect);
        let render_origin = self.camera.current.eye();

        // ── Hover detection (before click, drives cursor + ring) ──────────────
        self.hovered_body = hover_pos
            .filter(|p| rect.contains(*p))
            .and_then(|p| self.find_body_at(p, view_proj_relative, rect));

        // Cursor: pointer over bodies, grab over empty space, grabbing while dragging
        if response.dragged() {
            ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        } else if self.hovered_body.is_some() {
            ctx.set_cursor_icon(egui::CursorIcon::PointingHand);
        } else if hover_pos.map_or(false, |p| rect.contains(p)) {
            ctx.set_cursor_icon(egui::CursorIcon::Grab);
        }

        // ── Data-driven colour pipeline (SPLASH / yt-style) ───────────────────
        // Resolve the active ColorView once per frame, producing one RGB
        // triple per body. `None` means "use material colours" — bodies and
        // trails both fall back to Body::color downstream.
        let body_colors_override: Option<Vec<[u8; 3]>> = {
            let sel_clone = self.color_view.clone();
            sel_clone.and_then(|sel| self.evaluate_color_view(&sel))
        };

        // ── GPU body rendering ────────────────────────────────────────────────
        {
            let mut backend = self.backend.lock().unwrap();
            backend.begin();
            backend.show_grid = self.show_grid;
            backend.user_ev_stops = self.exposure_ev;
            // Floating Origin anchor for the frame. Renderers that
            // upload geometry consume this to produce camera-relative
            // `f32` positions.
            backend.render_origin = self.camera.current.eye();

            let bodies = self.system.bodies();

            // Build the per-frame scene lighting from every luminous body.
            // Intensities are normalised by the brightest source so the
            // primary star sits at 1.0 and companions scale relatively —
            // this keeps the shader's attenuation knob (`r_ref`) tuned to
            // scene scale rather than absolute luminosity units.
            let max_lum = bodies
                .iter()
                .filter(|b| b.is_luminous())
                .map(|b| b.luminosity)
                .fold(0.0_f64, f64::max);
            let lights: Vec<LightSpec> = if max_lum > 0.0 {
                bodies
                    .iter()
                    .filter(|b| b.is_luminous())
                    .map(|b| {
                        let rel = RenderRelativeVec3::from_world(
                            glam::DVec3::new(b.x, b.y, b.z),
                            render_origin,
                        );
                        LightSpec {
                            pos_relative: rel.as_array(),
                            intensity: (b.luminosity / max_lum) as f32,
                        }
                    })
                    .collect()
            } else {
                Vec::new()
            };

            // Characteristic distance for the 1/r² falloff. Without this,
            // scenes spanning >10× in radius (Solar System: Mercury at 0.4 AU
            // vs Neptune at 30 AU → 75× range → 5600× flux ratio) drop outer
            // planets to ~0.1% brightness — physically correct but visually
            // invisible against HDR black.
            //
            // Using the RMS distance from the primary light to all non-lit
            // bodies lands r_ref mid-scene: inner planets saturate (good, they
            // are already fully lit in reality), outer planets keep a usable
            // fraction. The RMS (vs. plain mean) weights toward outer bodies
            // deliberately — it's the outer planets that need the help.
            //
            // Computed in absolute `f64` space against the primary
            // luminous body — the LightSpec we built above is already
            // in render-relative coordinates and has lost a tiny amount
            // of precision, irrelevant here but worth taking from the
            // raw body for clarity.
            let primary_luminous = bodies.iter().find(|b| b.is_luminous());
            let r_ref = if let Some(primary) = primary_luminous {
                let (lx, ly, lz) = (primary.x, primary.y, primary.z);
                let (sum_sq, n) = bodies.iter().filter(|b| !b.is_luminous()).fold(
                    (0.0_f64, 0usize),
                    |(acc, k), b| {
                        let dx = b.x - lx;
                        let dy = b.y - ly;
                        let dz = b.z - lz;
                        (acc + dx * dx + dy * dy + dz * dz, k + 1)
                    },
                );
                if n > 0 { (sum_sq / n as f64).sqrt().max(1e-3) as f32 } else { 1.0 }
            } else {
                1.0
            };

            backend.set_scene_lighting(SceneLighting {
                lights,
                r_ref,
                // Backstop for bodies far beyond r_ref whose attenuation still
                // rounds to zero (e.g. Pluto in a Solar-System view). Without
                // this floor they collapse into pure black in HDR space.
                ambient_floor: 0.10,
                ..Default::default()
            });

            let body_canvas_h = rect.height().max(1.0);
            let body_focal_y = body_canvas_h / (2.0 * (FOV_Y_RAD * 0.5).tan());
            let body_mode = self.semantic_scale_mode;
            let body_min_px: f32 = match body_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            };

            // Photometry lights: every luminous body's position +
            // intrinsic luminosity, consumed by the bolometric flux
            // pipeline that drives sub-pixel reflective point sprites.
            let photometry_lights: Vec<(apsis::math::Vec3, f64)> = bodies
                .iter()
                .filter(|b| b.is_luminous())
                .map(|b| (apsis::math::Vec3::new(b.x, b.y, b.z), b.luminosity))
                .collect();

            // Reference magnitude that maps to linear pixel = 1.0
            // before auto-exposure metering compresses around the peak.
            // Shifted by the user EV offset (positive = brighter scene).
            let m_ref_base = -4.0_f64;
            let m_ref = m_ref_base - self.exposure_ev as f64;

            // Cross-fade band between disc and point sprite, in pixels
            // of projected radius.
            const DISK_OFF_PX: f32 = 0.5;
            const DISK_FULL_PX: f32 = 2.5;

            for (i, b) in bodies.iter().enumerate() {
                let rgb = match body_colors_override.as_ref() {
                    Some(colors) => colors[i],
                    None => b.color,
                };
                let body_world = glam::DVec3::new(b.x, b.y, b.z);
                let center_rel = RenderRelativeVec3::from_world(body_world, render_origin);
                let body_dist = center_rel.as_vec3().length().max(1e-6);
                let r_world = radius_world_3d(
                    b.physical_radius,
                    body_mode,
                    body_min_px,
                    body_dist,
                    body_focal_y,
                );
                let r_px = projected_radius_px(center_rel, r_world, body_canvas_h);

                let base = [rgb[0] as f32 / 255.0, rgb[1] as f32 / 255.0, rgb[2] as f32 / 255.0];
                let luminous = b.is_luminous();

                let (albedo_full, emissive_full) = if luminous {
                    // Emissive scales with luminosity so HDR has dynamic
                    // range for the bloom to grade against. Saturating
                    // sigmoid: brown dwarf ≈ 1×, Sun ≈ 7×, O-star ≤ 10×.
                    let lum_solar = b.luminosity as f32;
                    let scale = 1.0 + 9.0 * (1.0 - (-lum_solar).exp());
                    let e = [base[0] * scale, base[1] * scale, base[2] * scale, 1.0];
                    ([0.0, 0.0, 0.0, 1.0], e)
                } else {
                    ([base[0], base[1], base[2], 1.0], [0.0, 0.0, 0.0, 1.0])
                };

                // Cross-fade disc ↔ point applies to both branches —
                // sub-pixel luminous bodies need the sprite path or
                // they vanish from the rasteriser entirely. Point
                // sprites for reflective accumulate into HDR_R; for
                // luminous, into HDR_L (so bloom picks them up).
                let t = ((r_px - DISK_OFF_PX) / (DISK_FULL_PX - DISK_OFF_PX)).clamp(0.0, 1.0);
                let weight_disk = t * t * (3.0 - 2.0 * t);
                let weight_point = 1.0 - weight_disk;

                if weight_disk > 0.001 {
                    let albedo = [
                        albedo_full[0] * weight_disk,
                        albedo_full[1] * weight_disk,
                        albedo_full[2] * weight_disk,
                        albedo_full[3],
                    ];
                    let emissive = [
                        emissive_full[0] * weight_disk,
                        emissive_full[1] * weight_disk,
                        emissive_full[2] * weight_disk,
                        emissive_full[3],
                    ];
                    backend.draw_body(
                        center_rel,
                        r_world,
                        albedo,
                        emissive,
                        luminous,
                        b.albedo as f32,
                    );
                }

                if weight_point > 0.001 {
                    if let Some(pix) =
                        world_to_screen(body_world, render_origin, view_proj_relative, rect)
                    {
                        let body_pos = apsis::math::Vec3::new(b.x, b.y, b.z);
                        let observer = apsis::math::Vec3::new(
                            render_origin.x,
                            render_origin.y,
                            render_origin.z,
                        );
                        let mag = apsis::physics::photometry::apparent_magnitude(
                            b,
                            body_pos,
                            observer,
                            &photometry_lights,
                        );
                        let intensity =
                            apsis::physics::photometry::magnitude_to_linear_intensity(mag, m_ref)
                                as f32
                                * weight_point;
                        if intensity.is_finite() && intensity > 0.0 {
                            if luminous {
                                backend.draw_point_luminous([pix.x, pix.y], intensity, base);
                            } else {
                                backend.draw_point([pix.x, pix.y], intensity, base);
                            }
                        }
                    }
                }
            }

            // Predicted Keplerian orbits — pure visual annotation, never
            // feeds back into physics. Filter pipeline (Lote 3.2):
            //
            //   1. Hierarchy gate: primary known + level visible
            //   2. Degeneracy gate: reject near-parabolic and periapsis
            //      inside the primary's body radius
            //   3. Influence ranking: log(mass) × viewport-proximity
            //   4. Top-N truncation
            //
            // Background pass draws the survivors with faint alpha; a
            // final foreground pass paints the selected body on top.
            {
                let g_factor = self.system.g_factor();
                let t_sim = self.system.t();
                // Orbit-overlay sample projection: take each `f64` world
                // sample, shift into the render frame at full precision,
                // then apply the relative view-proj. Same pattern as the
                // body draw loop so an orbit polyline lands on the
                // pixels its body draws on.
                let project = |p: [f64; 3]| -> ([f32; 2], f32) {
                    let world = glam::DVec3::new(p[0], p[1], p[2]);
                    let rel = RenderRelativeVec3::from_world(world, render_origin);
                    let clip = view_proj_relative * rel.as_vec3().extend(1.0);
                    let w = clip.w;
                    if w > 0.0 {
                        let ndc = clip.truncate() / w;
                        let sx = rect.min.x + (ndc.x * 0.5 + 0.5) * rect.width();
                        let sy = rect.min.y + (-ndc.y * 0.5 + 0.5) * rect.height();
                        ([sx, sy], w)
                    } else {
                        ([f32::NAN, f32::NAN], w)
                    }
                };

                self.orbit_hierarchy.tick(bodies, g_factor);

                // If sim time jumped backwards (snapshot restore, system
                // reset), the cached invariants are stale: drop everything
                // and let the smoother cold-start. NaN sentinel = first
                // frame; skip the comparison.
                if !self.orbit_smoother_last_t.is_nan() && t_sim < self.orbit_smoother_last_t - 1e-9
                {
                    self.orbit_smoother.clear();
                }

                // Prune stale pins (e.g. after a collision-merge swap_remove
                // on the sim thread invalidated the index).
                let n_bodies = bodies.len();
                self.pinned_orbits.retain(|&i| i < n_bodies);

                // Clone names once for the smoother (cache key is the body's
                // stable name — survives swap_remove because index identity
                // is meaningless across collisions).
                let names_owned: Vec<String> = self.system.names().to_vec();

                // Pre-compute siblings for every primary used as an overlay
                // anchor, once per frame. The smoother re-uses these slices
                // to compute the tidal-vector weighted τ; doing this once
                // avoids re-scanning the body list per overlay candidate.
                let mut siblings_by_primary: std::collections::HashMap<usize, Vec<usize>> =
                    std::collections::HashMap::new();
                for j in 0..bodies.len() {
                    if let Some(p) = self.orbit_hierarchy.primary(j) {
                        siblings_by_primary.entry(p).or_default().push(j);
                    }
                }
                let empty_siblings: Vec<usize> = Vec::new();
                let siblings_for = |p: usize| -> &[usize] {
                    siblings_by_primary.get(&p).map(|v| v.as_slice()).unwrap_or(&empty_siblings)
                };

                if self.show_orbit_ellipses {
                    let bg_style = OrbitOverlayStyle::background_default();
                    let vp_center = rect.center();
                    let vp_half_diag = (rect.width().powi(2) + rect.height().powi(2)).sqrt() * 0.5;

                    // Filter pipeline runs on raw `compute_elements` because
                    // the degeneracy gate is a binary geometry test (e ≈ 1
                    // window, periapsis-inside-primary) and smoothing has
                    // negligible effect on those classifications. Only the
                    // top-N survivors get drawn through the smoother — a
                    // body that wasn't shown last frame cold-starts on
                    // first display, which is correct.
                    let mut candidates: Vec<(usize, usize, f32)> =
                        Vec::with_capacity(bodies.len().min(self.orbit_top_n * 2));

                    for i in 0..bodies.len() {
                        if self.selection.contains(i) {
                            continue;
                        }
                        if self.pinned_orbits.contains(&i) {
                            continue;
                        }
                        // System root has no Keplerian orbit; rendering
                        // one would misrepresent N-body dynamics.
                        if is_system_root(bodies, i) {
                            continue;
                        }
                        let Some(level) = self.orbit_hierarchy.level(i) else {
                            continue;
                        };
                        if !level_is_visible(level, &self.orbit_visible_levels) {
                            continue;
                        }
                        let Some(primary_idx) = self.orbit_hierarchy.primary(i) else {
                            continue;
                        };
                        let Some(el) = compute_elements(bodies, i, primary_idx, g_factor) else {
                            continue;
                        };
                        if self.orbit_hide_degenerate
                            && is_degenerate_orbit(&el, &bodies[primary_idx])
                        {
                            continue;
                        }

                        let b = &bodies[i];
                        // Viewport-proximity weight from the actual 3D
                        // projection of the body. Off-screen and
                        // behind-camera bodies skip the candidate
                        // pool — the orbit ranking only looks at what
                        // the user might plausibly be looking at.
                        let world = glam::DVec3::new(b.x, b.y, b.z);
                        let Some(sp) =
                            world_to_screen(world, render_origin, view_proj_relative, rect)
                        else {
                            continue;
                        };
                        let dx = sp.x - vp_center.x;
                        let dy = sp.y - vp_center.y;
                        let d_norm = (dx * dx + dy * dy).sqrt() / vp_half_diag.max(1.0);
                        if d_norm > 2.0 {
                            continue;
                        }
                        let viewport_weight = 1.0 / (1.0 + d_norm * d_norm);
                        let mass_factor = (1.0_f32 + b.mass as f32).ln().max(0.0);
                        let influence = mass_factor * viewport_weight;

                        candidates.push((i, primary_idx, influence));
                    }

                    candidates
                        .sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
                    if candidates.len() > self.orbit_top_n {
                        candidates.truncate(self.orbit_top_n);
                    }

                    for (i, primary_idx, _) in &candidates {
                        let Some(el) = self.orbit_smoother.smoothed(
                            bodies,
                            &names_owned,
                            *i,
                            *primary_idx,
                            siblings_for(*primary_idx),
                            g_factor,
                            t_sim,
                        ) else {
                            continue;
                        };
                        let primary = &bodies[*primary_idx];
                        let primary_pos = [primary.x, primary.y, primary.z];
                        let sampled = el.sample_orbit(primary_pos, 64);
                        draw_orbit_polyline(&mut backend, &sampled, project, &bg_style, None, None);
                        draw_orbit_apsides(&mut backend, &el, primary_pos, project, &bg_style);
                    }
                }

                if !self.pinned_orbits.is_empty() {
                    let fg_style = OrbitOverlayStyle::selected_default();
                    let pins: Vec<usize> = self.pinned_orbits.iter().copied().collect();
                    for i in pins {
                        if self.selection.contains(i) {
                            continue;
                        }
                        if is_system_root(bodies, i) {
                            continue;
                        }
                        let primary =
                            self.orbit_hierarchy.primary(i).or_else(|| dominant_primary(bodies, i));
                        let Some(primary_idx) = primary else {
                            continue;
                        };
                        let Some(el) = self.orbit_smoother.smoothed(
                            bodies,
                            &names_owned,
                            i,
                            primary_idx,
                            siblings_for(primary_idx),
                            g_factor,
                            t_sim,
                        ) else {
                            continue;
                        };
                        let primary_b = &bodies[primary_idx];
                        let primary_pos = [primary_b.x, primary_b.y, primary_b.z];
                        let sampled = el.sample_orbit(primary_pos, 96);
                        let body = &bodies[i];
                        let anchor = closest_sample_index(&sampled, [body.x, body.y, body.z]);
                        draw_orbit_polyline(
                            &mut backend,
                            &sampled,
                            project,
                            &fg_style,
                            None,
                            anchor,
                        );
                        draw_orbit_apsides(&mut backend, &el, primary_pos, project, &fg_style);
                    }
                }

                if let Some(idx) = self.selection.single() {
                    if idx < bodies.len() && !is_system_root(bodies, idx) {
                        let primary = self
                            .orbit_hierarchy
                            .primary(idx)
                            .or_else(|| dominant_primary(bodies, idx));
                        if let Some(primary_idx) = primary {
                            if let Some(el) = self.orbit_smoother.smoothed(
                                bodies,
                                &names_owned,
                                idx,
                                primary_idx,
                                siblings_for(primary_idx),
                                g_factor,
                                t_sim,
                            ) {
                                let primary = &bodies[primary_idx];
                                let primary_pos = [primary.x, primary.y, primary.z];
                                let sampled = el.sample_orbit(primary_pos, 128);
                                let style = OrbitOverlayStyle::selected_default();
                                let body = &bodies[idx];
                                let anchor =
                                    closest_sample_index(&sampled, [body.x, body.y, body.z]);
                                draw_orbit_polyline_with_halo(
                                    &mut backend,
                                    &sampled,
                                    project,
                                    &style,
                                    anchor,
                                );
                                draw_orbit_apsides(&mut backend, &el, primary_pos, project, &style);
                            }
                        }
                    }
                }

                // Drop cache entries for bodies removed since the last
                // frame (collision merges, manual deletions). Cheap when
                // the cache is empty / small.
                self.orbit_smoother.prune(&names_owned);
                self.orbit_smoother_last_t = t_sim;
            }

            backend.trail_style = self.trail_style_preset.style(self.trail_width);

            // Apply any COM shift the physics thread accumulated, then drain
            // the per-step trail samples it produced and push them into the
            // ring buffer. Sampling decisions happen on the physics thread;
            // this side is a pure consumer.
            let (shift_x, shift_y) = self.system.take_pending_com_shift();
            self.trail_recorder.apply_com_shift(shift_x, shift_y);
            let samples = self.system.take_trail_samples();
            self.trail_recorder.ingest_with_colors(
                samples,
                self.system.bodies(),
                body_colors_override.as_deref(),
            );

            if self.show_trails {
                let bodies = self.system.bodies();
                let names = self.system.names();
                // Three-tier visibility: per-body override wins; otherwise
                // an authored body (one with a non-empty template name)
                // shows iff its class passes the per-class filter; bodies
                // without a name are off by default. The mass-ratio
                // heuristic is gone — it failed for compact-object
                // systems and multi-system scenes where the dominant
                // mass scale was unrelated to "is this body interesting".
                backend.trail_visibility = Some(
                    bodies
                        .iter()
                        .enumerate()
                        .map(|(i, b)| match self.trail_per_body_override.get(&i) {
                            Some(&explicit) => explicit,
                            None => {
                                let authored = names.get(i).map_or(false, |n| !n.is_empty());
                                authored && self.trail_class_filter.allows(b.class)
                            },
                        })
                        .collect(),
                );
                backend.trail_buffer =
                    Some(std::sync::Arc::new(self.trail_recorder.buffer().clone()));
            } else {
                backend.trail_visibility = None;
                backend.trail_buffer = None;
            };
        }

        // ── Place-mode: click or drag-release to spawn a body ────────────────
        // Editing is locked during a Precision Run (REBOUND-aligned).
        // Place-mode silently ignores input in that state; the UI
        // already greyed out the Add tool in the rail, so user input
        // should not reach here anyway — belt-and-suspenders.
        if self.place_mode && self.is_editing_locked() {
            self.place_drag_start = None;
            ctx.set_cursor_icon(egui::CursorIcon::NotAllowed);
        } else if self.place_mode {
            let pointer = ctx.input(|i| i.pointer.clone());

            // Track drag start (only when pressing on empty space)
            if response.drag_started() {
                if let Some(pos) = pointer.press_origin() {
                    // Only start a place-drag if not clicking an existing body
                    if self.find_body_at(pos, view_proj_relative, rect).is_none() {
                        self.place_drag_start = Some(pos);
                    }
                }
            }

            // Draw velocity arrow while dragging
            if let Some(start) = self.place_drag_start {
                if let Some(cur) = pointer.hover_pos() {
                    let painter = ui.painter();
                    let delta = cur - start;
                    let len = delta.length();
                    if len > 4.0 {
                        // Line
                        painter.line_segment(
                            [start, cur],
                            Stroke::new(1.5, Color32::from_rgba_premultiplied(120, 200, 255, 200)),
                        );
                        // Arrowhead
                        let dir = delta / len;
                        let perp = egui::vec2(-dir.y, dir.x);
                        let tip = cur;
                        let base = tip - dir * 8.0;
                        painter.add(egui::Shape::convex_polygon(
                            vec![tip, base + perp * 4.0, base - perp * 4.0],
                            Color32::from_rgba_premultiplied(120, 200, 255, 200),
                            Stroke::NONE,
                        ));
                        // Velocity preview from the drag delta in
                        // world AU on the ecliptic. Both endpoints
                        // unproject through the same camera so their
                        // difference is dimensional regardless of
                        // zoom or tilt.
                        let speed = match (
                            screen_to_world_on_z_plane(
                                start,
                                view_proj_relative,
                                rect,
                                render_origin,
                            ),
                            screen_to_world_on_z_plane(
                                cur,
                                view_proj_relative,
                                rect,
                                render_origin,
                            ),
                        ) {
                            (Some(s), Some(e)) => {
                                let d = e - s;
                                (d.x * d.x + d.y * d.y).sqrt() * 0.5
                            },
                            _ => 0.0,
                        };
                        painter.text(
                            cur + egui::vec2(8.0, -12.0),
                            egui::Align2::LEFT_CENTER,
                            format!("v={:.3}", speed),
                            egui::FontId::monospace(9.5),
                            Color32::from_rgba_premultiplied(140, 210, 255, 200),
                        );
                    }
                }
            }

            // On release: spawn body
            if response.drag_stopped() || response.clicked() {
                let start = self.place_drag_start.take();
                if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                    // Don't spawn if clicking an existing body
                    if self.find_body_at(cursor, view_proj_relative, rect).is_none() {
                        let spawn_pos = start.unwrap_or(cursor);
                        let Some(spawn_world) = screen_to_world_on_z_plane(
                            spawn_pos,
                            view_proj_relative,
                            rect,
                            render_origin,
                        ) else {
                            self.place_drag_start = None;
                            return;
                        };
                        let wx = spawn_world.x;
                        let wy = spawn_world.y;

                        // Velocity from drag delta — same dimensional
                        // unproject as the preview above.
                        let (vx, vy) = match start
                            .and_then(|s| {
                                screen_to_world_on_z_plane(
                                    s,
                                    view_proj_relative,
                                    rect,
                                    render_origin,
                                )
                            })
                            .zip(screen_to_world_on_z_plane(
                                cursor,
                                view_proj_relative,
                                rect,
                                render_origin,
                            )) {
                            Some((s, e)) => {
                                let d = e - s;
                                (d.x * 0.5, d.y * 0.5)
                            },
                            None => (0.0, 0.0),
                        };

                        let body = apsis::domain::body::Body::from_preset(
                            self.place_preset,
                            self.place_mass,
                        )
                        .at(wx, wy)
                        .with_velocity(vx, vy);

                        self.push_undo(UndoRecord::AddedBodies(1));
                        self.system.add_named_body(apsis::domain::body::NamedBody {
                            body,
                            name: Some(self.place_preset.display_name.to_owned()),
                        });
                    }
                }
            }

            // Override cursor to crosshair while in place-mode (unless dragging a body)
            if self.dragging_body.is_none() {
                ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
            }

            ctx.request_repaint();
        } else {
            // ── Normal / Shift click: select body or toggle multi-select ──────
            self.place_drag_start = None;
            if response.clicked() {
                if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                    let shift = ctx.input(|i| i.modifiers.shift);
                    match self.find_body_at(cursor, view_proj_relative, rect) {
                        Some(idx) => {
                            if shift {
                                // Shift+click: toggle body into/out of the selection.
                                let prev = std::mem::take(&mut self.selection);
                                self.selection = prev.toggle(idx);
                                match &self.selection {
                                    BodySelection::Multi(_) => {
                                        // Multi-select: disable follow and clear the edit form.
                                        self.follow_selected_body = false;
                                        self.selection_form = None;
                                    },
                                    BodySelection::Single(i) => {
                                        // Toggled back to a single body — restore normal state.
                                        let i = *i;
                                        let body = self.system.bodies()[i];
                                        let name = self.system.name(i).to_owned();
                                        self.follow_selected_body = true;
                                        self.selection_form =
                                            Some(SelectionForm::from_body(&body, &name));
                                    },
                                    BodySelection::None => {
                                        self.follow_selected_body = false;
                                        self.selection_form = None;
                                    },
                                }
                            } else {
                                let body = self.system.bodies()[idx];
                                self.selection = BodySelection::select_single(idx);
                                self.follow_selected_body = true;
                                let name = self.system.name(idx).to_owned();
                                self.selection_form = Some(SelectionForm::from_body(&body, &name));

                                let canvas_h = rect.height().max(1.0);
                                let focal_y = canvas_h / (2.0 * (FOV_Y_RAD * 0.5).tan());
                                let r = (body.physical_radius as f32).max(1e-12);
                                let dist = ((r * focal_y) / FRAME_TARGET_PX)
                                    .max(NEAR_PLANE * 5.0)
                                    .max(r * 4.0) as f64;

                                // Capture the click-time pose. Phase-locked
                                // transition lerps every dimension by a
                                // single fraction so distance, pivot, and
                                // orientation reach the same perceptual
                                // progress at the same wall time.
                                self.follow_transition =
                                    Some(crate::app::camera::FollowTransition::capture(
                                        idx,
                                        self.camera.current,
                                        self.camera.current.azimuth,
                                        self.camera.current.elevation,
                                        dist,
                                    ));
                                self.camera.target.distance = dist;
                            }
                        },
                        None => {
                            self.selection = BodySelection::default();
                            self.follow_selected_body = false;
                            self.follow_transition = None;
                            self.selection_form = None;
                        },
                    }
                }
            }
        }

        // ── Template drag-drop ────────────────────────────────────────────────
        // Check if a drag from the template panel was released over this canvas.
        if self.template_drag.is_some() {
            let released = ctx.input(|i| i.pointer.any_released());
            let drop_pos = ctx.input(|i| i.pointer.interact_pos());

            if released {
                if let (Some(build_fn), Some(screen_pos)) = (self.template_drag.take(), drop_pos) {
                    if rect.contains(screen_pos) {
                        // Convert screen pos → world pos by ray-
                        // casting through the camera onto the
                        // ecliptic. Drop is rejected when the camera
                        // can't see that plane at the cursor.
                        let Some(world) = screen_to_world_on_z_plane(
                            screen_pos,
                            view_proj_relative,
                            rect,
                            render_origin,
                        ) else {
                            return;
                        };
                        let wx = world.x;
                        let wy = world.y;
                        let template = build_fn();
                        self.active_units = template.units;
                        let bodies = instantiate_at(&template, wx, wy);
                        self.push_undo(UndoRecord::AddedBodies(bodies.len()));
                        self.system.add_named_bodies(bodies);
                        self.pending_fit = false; // dropped at explicit position — no auto-fit
                    } else {
                        // Released outside canvas — discard
                        self.template_drag = None;
                    }
                } else {
                    self.template_drag = None;
                }
            } else {
                // Still dragging — show a ghost cursor
                if hover_pos.map_or(false, |p| rect.contains(p)) {
                    ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
                }
            }
        }

        // Keep animating while a body is selected (pulsing ring needs repaints)
        if self.selection.is_some() {
            ctx.request_repaint();
        }

        // ── GPU paint callback ────────────────────────────────────────────────
        let device = self.device.as_ref().unwrap().clone();
        let queue = self.queue.as_ref().unwrap().clone();
        let format = self.format.unwrap();

        let view_proj_relative_arr = view_proj_relative.to_cols_array_2d();

        ui.painter().add(egui_wgpu::Callback::new_paint_callback(
            rect,
            CallbackFn {
                backend: self.backend.clone(),
                device,
                queue,
                format,
                screen: [rect.width(), rect.height()],
                viewport_min: [rect.min.x, rect.min.y],
                view_proj_relative: view_proj_relative_arr,
            },
        ));

        // ── Overlay: rings + labels (on top of GPU layer) ─────────────────────
        self.draw_overlay(ui, rect, view_proj_relative, time);

        // ── FPS / frame-time HUD (top-right of canvas, subtle) ────────────────
        if self.show_fps_hud {
            self.draw_fps_overlay(ui, rect);
        }

        // ── 3D camera axis triad (bottom-left, debug) ─────────────────────────
        if self.show_camera_triad {
            self.draw_camera_triad(ui, rect);
        }

        // ── Loading overlay ───────────────────────────────────────────────────
        if self.system.is_loading() {
            self.draw_loading_overlay(ui, rect, time);
            ctx.request_repaint();
        }
    }

    fn draw_fps_overlay(&self, ui: &egui::Ui, rect: egui::Rect) {
        let d = &self.diagnostics;
        if d.is_idle() || d.warming() {
            return;
        }
        let text = format!("{:>4.0} FPS · {:>5.2} ms", d.fps(), d.frame_ms());
        let pos = egui::pos2(rect.right() - 12.0, rect.top() + 10.0);
        let color = crate::app::design::tokens::color::fg::TERTIARY;
        ui.painter().text(pos, egui::Align2::RIGHT_TOP, text, FontId::monospace(10.0), color);
    }

    fn draw_camera_triad(&self, ui: &egui::Ui, rect: egui::Rect) {
        use glam::DVec3;

        const RADIUS_PX: f32 = 28.0;
        const X_COL: Color32 = Color32::from_rgb(212, 102, 102);
        const Y_COL: Color32 = Color32::from_rgb(132, 192, 132);
        const Z_COL: Color32 = Color32::from_rgb(122, 154, 232);

        // The triad only needs the camera orientation, never the eye
        // translation, so the rotation-only view (the same one shaders
        // consume under Floating Origin) is the right choice.
        let view = self.camera.current.view_rotation_only();
        let center = egui::pos2(rect.left() + 44.0, rect.bottom() - 44.0);

        // World axes are direction vectors, so transform_vector3 applies
        // only the rotation part of the view matrix. Camera-space basis
        // is right-handed: +x right, +y up, −z into the screen.
        let mut axes: [(f32, f32, f32, Color32, &str); 3] = [
            project(view.transform_vector3(DVec3::X), X_COL, "x"),
            project(view.transform_vector3(DVec3::Y), Y_COL, "y"),
            project(view.transform_vector3(DVec3::Z), Z_COL, "z"),
        ];
        axes.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        let painter = ui.painter();
        let bg = crate::app::design::tokens::color::fg::TERTIARY;
        painter.circle_stroke(center, RADIUS_PX + 2.0, Stroke::new(0.5, bg.gamma_multiply(0.4)));

        for (dx, dy, depth, color, label) in axes {
            let tip = center + egui::vec2(dx * RADIUS_PX, dy * RADIUS_PX);
            // Axes pointing away from the camera dim — gives the ring a
            // sense of depth without needing a real projection.
            let alpha = if depth > 0.0 { 0.45 } else { 1.0 };
            painter.line_segment([center, tip], Stroke::new(1.5, color.gamma_multiply(alpha)));
            painter.text(
                tip + egui::vec2(dx * 7.0, dy * 7.0),
                egui::Align2::CENTER_CENTER,
                label,
                FontId::monospace(9.5),
                color.gamma_multiply(alpha),
            );
        }
        painter.circle_filled(center, 1.5, bg);
    }

    // ── Overlay ───────────────────────────────────────────────────────────────

    fn draw_overlay(
        &self,
        ui: &egui::Ui,
        rect: egui::Rect,
        view_proj_relative: glam::Mat4,
        time: f32,
    ) {
        let bodies = self.system.bodies();
        let names = self.system.names();

        if bodies.is_empty() {
            return;
        }

        let render_origin = self.camera.current.eye();
        let max_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);
        let canvas_h = rect.height().max(1.0);
        let focal_y = canvas_h / (2.0 * (FOV_Y_RAD * 0.5).tan());
        let mode = self.semantic_scale_mode;
        let min_px: f32 = match mode {
            SemanticScaleMode::Physical => 0.0,
            SemanticScaleMode::Comparative => 3.0,
            SemanticScaleMode::Illustrative => 5.0,
        };

        // Label visibility threshold scales with camera distance so a wide
        // overview shows fewer labels than a close-up framing. Pivot
        // distance comes straight from the orbit camera; the f64
        // subtraction stays precise at AU scale.
        let cam_dist = (render_origin - self.camera.current.pivot).length().max(1e-3) as f32;
        let importance_threshold = (2.0_f64 / cam_dist as f64).clamp(0.001, 1.0);

        let painter = ui.painter();
        let font = FontId::proportional(LABEL_FONT_SIZE);

        // Pulse for selection ring: ±1.5 px at ~3.5 Hz
        let pulse = (time * 3.5).sin() * 1.5_f32;

        for (i, (body, name)) in bodies.iter().zip(names.iter()).enumerate() {
            let world = glam::DVec3::new(body.x, body.y, body.z);
            let Some(body_pos) = world_to_screen(world, render_origin, view_proj_relative, rect)
            else {
                continue;
            };
            let center_rel = RenderRelativeVec3::from_world(world, render_origin);
            let body_dist = center_rel.as_vec3().length().max(1e-6);
            let r_world = radius_world_3d(body.physical_radius, mode, min_px, body_dist, focal_y);
            let visual_r = projected_radius_px(center_rel, r_world, canvas_h);
            let px = body_pos.x;
            let py = body_pos.y;

            // Single body selected → pulsing ring; one of many selected → dim ring.
            let is_primary = self.selection.single() == Some(i);
            let is_in_multi =
                matches!(&self.selection, BodySelection::Multi(_)) && self.selection.contains(i);
            let is_hovered = self.hovered_body == Some(i) && !is_primary && !is_in_multi;

            // ── Hover ring ───────────────────────────────────────────────
            if is_hovered {
                painter.circle_stroke(
                    body_pos,
                    visual_r + RING_GAP - 1.0,
                    Stroke::new(1.0, Color32::from_rgba_premultiplied(160, 160, 200, 90)),
                );
            }

            // ── Multi-select ring (dim, no pulse) ────────────────────────
            if is_in_multi {
                painter.circle_stroke(
                    body_pos,
                    visual_r + RING_GAP,
                    Stroke::new(1.0, Color32::from_rgba_premultiplied(200, 200, 255, 120)),
                );
            }

            // ── Primary selection ring (pulsing) ─────────────────────────
            if is_primary {
                let r1 = visual_r + RING_GAP + pulse;
                let r2 = r1 + 3.5;

                // Outer dim halo
                painter.circle_stroke(
                    body_pos,
                    r2,
                    Stroke::new(0.8, Color32::from_rgba_premultiplied(130, 130, 200, 55)),
                );
                // Main ring
                painter.circle_stroke(
                    body_pos,
                    r1,
                    Stroke::new(1.5, Color32::from_rgba_premultiplied(200, 200, 255, 210)),
                );
                // Inner tick marks at compass points for clarity
                for angle in [0.0_f32, 90.0, 180.0, 270.0] {
                    let rad = angle.to_radians();
                    let dir = egui::vec2(rad.cos(), rad.sin());
                    let inner = body_pos + dir * (r1 - 2.5);
                    let outer = body_pos + dir * (r1 + 2.5);
                    painter.line_segment(
                        [inner, outer],
                        Stroke::new(1.5, Color32::from_rgba_premultiplied(200, 200, 255, 180)),
                    );
                }
            }

            // ── Name label ───────────────────────────────────────────────
            let importance = if max_mass > 0.0 { body.mass / max_mass } else { 0.0 };
            let is_any_selected = is_primary || is_in_multi;
            let show_label =
                visual_r >= 5.0 || importance >= importance_threshold || is_any_selected;

            if show_label {
                // Cap offset so the label never drifts far from the body
                let offset_y = (visual_r + 4.0).min(MAX_LABEL_OFFSET_PX);
                let label_pos = egui::pos2(px, py + offset_y);

                let color = if is_any_selected {
                    Color32::from_rgb(220, 220, 255)
                } else {
                    Color32::from_rgba_premultiplied(175, 175, 195, 195)
                };

                // Shadow
                painter.text(
                    label_pos + egui::vec2(1.0, 1.0),
                    egui::Align2::CENTER_TOP,
                    name.as_str(),
                    font.clone(),
                    Color32::from_black_alpha(110),
                );
                painter.text(
                    label_pos,
                    egui::Align2::CENTER_TOP,
                    name.as_str(),
                    font.clone(),
                    color,
                );
            }
        }
    }

    // ── Loading overlay ───────────────────────────────────────────────────────

    fn draw_loading_overlay(&self, ui: &egui::Ui, rect: egui::Rect, time: f32) {
        let painter = ui.painter();

        // Dim backdrop — subtle, doesn't obliterate the scene.
        painter.rect_filled(rect, 0.0, Color32::from_black_alpha(120));

        let cx = rect.center();

        // ── Spinner: three arcs rotating at different phases ──────────────────
        // Outer ring
        draw_spinner_arc(
            painter,
            cx,
            22.0,
            2.5,
            time,
            1.0,
            Color32::from_rgba_premultiplied(160, 160, 255, 220),
        );
        // Middle ring (counter-rotate, dimmer)
        draw_spinner_arc(
            painter,
            cx,
            14.0,
            2.0,
            -time * 1.4,
            0.75,
            Color32::from_rgba_premultiplied(120, 120, 200, 150),
        );
        // Inner dot
        painter.circle_filled(cx, 3.5, Color32::from_rgba_premultiplied(180, 180, 255, 200));

        // ── "LOADING" label ────────────────────────────────────────────────────
        let label_pos = egui::pos2(cx.x, cx.y + 36.0);
        // Shadow
        painter.text(
            label_pos + egui::vec2(1.0, 1.0),
            egui::Align2::CENTER_TOP,
            "LOADING",
            eframe::egui::FontId::proportional(11.0),
            Color32::from_black_alpha(120),
        );
        painter.text(
            label_pos,
            egui::Align2::CENTER_TOP,
            "LOADING",
            eframe::egui::FontId::proportional(11.0),
            Color32::from_rgba_premultiplied(170, 170, 220, 210),
        );
    }

    // ── Hit-test ──────────────────────────────────────────────────────────────

    fn find_body_at(
        &self,
        cursor: Pos2,
        view_proj_relative: glam::Mat4,
        rect: egui::Rect,
    ) -> Option<usize> {
        let bodies = self.system.bodies();
        let render_origin = self.camera.current.eye();
        let canvas_h = rect.height().max(1.0);
        let focal_y = canvas_h / (2.0 * (FOV_Y_RAD * 0.5).tan());
        let mode = self.semantic_scale_mode;
        let min_px: f32 = match mode {
            SemanticScaleMode::Physical => 0.0,
            SemanticScaleMode::Comparative => 3.0,
            SemanticScaleMode::Illustrative => 5.0,
        };

        // Pick the front-most candidate among bodies whose projected disc
        // covers the cursor. Front-to-back ordering matters in 3D where
        // bodies can occlude each other regardless of insertion order.
        let mut best: Option<(usize, f32)> = None;
        for (i, b) in bodies.iter().enumerate() {
            let world = glam::DVec3::new(b.x, b.y, b.z);
            let Some(screen) = world_to_screen(world, render_origin, view_proj_relative, rect)
            else {
                continue;
            };
            let center_rel = RenderRelativeVec3::from_world(world, render_origin);
            let dist_vec = center_rel.as_vec3();
            let body_dist = dist_vec.length().max(1e-6);
            let r_world = radius_world_3d(b.physical_radius, mode, min_px, body_dist, focal_y);
            let r = projected_radius_px(center_rel, r_world, canvas_h).max(MIN_HIT_PX);
            let dx = cursor.x - screen.x;
            let dy = cursor.y - screen.y;
            if dx * dx + dy * dy > r * r {
                continue;
            }
            let cam_dist_sq = dist_vec.length_squared();
            if best.is_none_or(|(_, prev_d2)| cam_dist_sq < prev_d2) {
                best = Some((i, cam_dist_sq));
            }
        }
        best.map(|(i, _)| i)
    }

    /// Resolve the active [`ColorViewSelection`] against current state.
    ///
    /// Returns `None` when the selection is invalid (unknown ID) or the
    /// body list is empty. Caches the resolved data range in
    /// `self.color_view_range` so the colour-bar UI can render without
    /// re-evaluating the field.
    fn evaluate_color_view(
        &mut self,
        sel: &crate::render::color::ColorViewSelection,
    ) -> Option<Vec<[u8; 3]>> {
        use crate::render::color;

        let bodies = self.system.bodies();
        if bodies.is_empty() {
            self.color_view_range = None;
            return None;
        }

        let field = self.field_registry.get(&sel.field_id)?;
        let normalizer = self.normalizer_registry.get(&sel.normalizer_id)?;
        let colormap = self.colormap_registry.get(&sel.colormap_id)?;

        let ctx = apsis::domain::field::FieldContext {
            bodies,
            accelerations: self.system.accelerations(),
            t: self.system.t(),
            g_factor: self.system.g_factor(),
        };

        let out = color::compute(field, normalizer, colormap, sel.range, &ctx);
        self.color_view_range = Some(out.resolved_range);
        Some(out.colors)
    }
}

// ── Spinner helpers ───────────────────────────────────────────────────────────

/// Draw a rotating arc at `center` with the given `radius` and line `width`.
/// `time` drives rotation speed; `arc_frac` (0..1) is how much of the circle
/// is filled (e.g. 0.75 = 270°). Alpha fades from full at the head to 0 at
/// the tail for a "comet tail" effect.
fn draw_spinner_arc(
    painter: &egui::Painter,
    center: Pos2,
    radius: f32,
    width: f32,
    time: f32,
    arc_frac: f32,
    color: Color32,
) {
    let segments = 32usize;
    let arc_radians = arc_frac * std::f32::consts::TAU;
    let angle_step = arc_radians / segments as f32;
    let base_angle = time * 2.2; // rotation speed

    let [r, g, b, _] = color.to_array();

    for i in 0..segments {
        let t_frac = i as f32 / segments as f32; // 0 = tail, 1 = head
        let alpha = (t_frac * t_frac * 255.0) as u8; // quadratic fade

        let a0 = base_angle + i as f32 * angle_step;
        let a1 = a0 + angle_step * 1.1; // slight overlap to avoid gaps

        let p0 = center + egui::vec2(a0.cos(), a0.sin()) * radius;
        let p1 = center + egui::vec2(a1.cos(), a1.sin()) * radius;

        painter.line_segment(
            [p0, p1],
            Stroke::new(width, Color32::from_rgba_premultiplied(r, g, b, alpha)),
        );
    }
}
