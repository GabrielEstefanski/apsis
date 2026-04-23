use crate::app::render_params::{RenderParams, compute_render_radius};
use crate::app::ui::{SelectionForm, SemanticScaleMode, SimulationApp, UndoRecord};
use gravity_sim_core::physics::orbital::{compute_elements, dominant_primary};
use crate::render::lighting::{LightSpec, SceneLighting};
use crate::render::CallbackFn;
use crate::render::orbit_overlay::{
    OrbitOverlayStyle, draw_orbit_apsides, draw_orbit_polyline,
};
use gravity_sim_core::templates::instantiate_at;
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
/// Follow-mode exponential smoothing rate (1/s). Higher = tighter tracking.
/// At 60 fps, `alpha = 1 - exp(-dt * FOLLOW_RATE) ≈ 0.55` — responsive without
/// feeling jolty on the first-click transition.
const FOLLOW_RATE: f32 = 48.0;
/// If the body has moved farther than this many pixels between frames, snap
/// the camera instead of lerping. Prevents losing fast bodies at high sim
/// speeds (where a single publish can advance the body across the viewport).
const FOLLOW_SNAP_PX: f32 = 120.0;

/// Eccentricity window around 1.0 within which an orbit is considered
/// numerically degenerate (near-parabolic). Chosen deliberately tight so
/// high-eccentricity cometary orbits (e ≈ 0.95 – 0.99) remain visible —
/// those are exactly the orbits worth looking at.
const DEGEN_ECC_WINDOW: f64 = 0.005;

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
    el: &gravity_sim_core::physics::orbital::OrbitalElements,
    primary: &gravity_sim_core::domain::body::Body,
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

        // ── Zoom with inertia ─────────────────────────────────────────────────
        if scroll_y.abs() > 0.0 {
            self.zoom_vel += scroll_y * 0.0004;
            self.zoom_vel = self.zoom_vel.clamp(-0.15, 0.15);
        }

        if self.zoom_vel.abs() > 0.00005 {
            let old_scale = self.scale;
            self.scale = (self.scale * (1.0 + self.zoom_vel)).clamp(0.001, 50_000.0);

            if let Some(mouse) = hover_pos {
                let ratio = self.scale / old_scale;
                self.offset += (mouse - center) * (1.0 - ratio);
                // keep animation target in sync with the scale-adjusted offset
                if let Some(t) = self.camera_anim_target.as_mut() {
                    *t *= ratio;
                }
            }

            self.zoom_vel *= 0.80;
            if self.zoom_vel.abs() < 0.00005 {
                self.zoom_vel = 0.0;
            } else {
                ctx.request_repaint();
            }
        }

        // ── Pan with inertia ──────────────────────────────────────────────────
        if response.dragged() {
            let delta = response.drag_delta();
            self.offset += delta;
            // Cancel smooth-pan if the user takes manual control
            self.camera_anim_target = None;
            self.follow_selected_body = false;
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

        // ── Body-follow tracking ──────────────────────────────────────────────
        // Applied AFTER zoom and pan so neither can drag the followed body off
        // centre. Uses frame-rate-independent exponential smoothing and snaps
        // to the target when the body outran the lerp in a single frame.
        if self.follow_selected_body {
            if let Some(idx) = self.selected_body {
                if let Some(body) = self.system.bodies().get(idx) {
                    let target =
                        egui::vec2(-body.x as f32 * self.scale, -body.y as f32 * self.scale);
                    let delta = target - self.offset;
                    let alpha = 1.0 - (-dt * FOLLOW_RATE).exp();
                    if delta.length_sq() >= FOLLOW_SNAP_PX * FOLLOW_SNAP_PX {
                        // Body moved further than we can smoothly catch — snap.
                        self.offset = target;
                    } else if delta.length_sq() > 0.01 {
                        self.offset += delta * alpha;
                    }
                    ctx.request_repaint();
                } else {
                    self.follow_selected_body = false;
                    self.selected_body = None;
                    self.selection_form = None;
                }
            } else {
                self.follow_selected_body = false;
            }
        }

        // ── Render params ─────────────────────────────────────────────────────
        let render_params = RenderParams {
            world_scale: self.scale,
            mode: self.semantic_scale_mode,
            min_px: match self.semantic_scale_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            },
        };

        // ── Hover detection (before click, drives cursor + ring) ──────────────
        let center_after_pan = rect.center() + self.offset;
        self.hovered_body = hover_pos
            .filter(|p| rect.contains(*p))
            .and_then(|p| self.find_body_at(p, center_after_pan, render_params));

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

            let bodies = self.system.bodies();

            let mut screen_positions = Vec::with_capacity(bodies.len());

            for b in bodies {
                let px = center_after_pan.x + b.x as f32 * self.scale;
                let py = center_after_pan.y + b.y as f32 * self.scale;

                screen_positions.push(([px, py], b));
            }

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
                    .map(|b| LightSpec {
                        world_pos: [b.x as f32, b.y as f32, 0.0],
                        intensity: (b.luminosity / max_lum) as f32,
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
            let r_ref = if let Some(primary) = lights.first() {
                let lx = primary.world_pos[0] as f64;
                let ly = primary.world_pos[1] as f64;
                let (sum_sq, n) = bodies
                    .iter()
                    .filter(|b| !b.is_luminous())
                    .fold((0.0_f64, 0usize), |(acc, k), b| {
                        let dx = b.x - lx;
                        let dy = b.y - ly;
                        (acc + dx * dx + dy * dy, k + 1)
                    });
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

            for (i, (sp, b)) in screen_positions.iter().enumerate() {
                let rgb = match body_colors_override.as_ref() {
                    Some(colors) => colors[i],
                    None => b.color,
                };
                let r = compute_render_radius(b.physical_radius, render_params);

                // Luminous bodies: self-lit disc (emissive carries the colour,
                // albedo stays dark so the unlit side doesn't darken their
                // surface). Non-luminous: pure albedo, no self-emission.
                let base = [
                    rgb[0] as f32 / 255.0,
                    rgb[1] as f32 / 255.0,
                    rgb[2] as f32 / 255.0,
                ];
                let (albedo, emissive) = if b.is_luminous() {
                    ([0.0, 0.0, 0.0, 1.0], [base[0], base[1], base[2], 1.0])
                } else {
                    ([base[0], base[1], base[2], 1.0], [0.0, 0.0, 0.0, 1.0])
                };
                let world_pos = [b.x as f32, b.y as f32, 0.0];

                backend.draw_body(*sp, r, world_pos, albedo, emissive);
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
                let scale = self.scale;
                let cx = center_after_pan.x;
                let cy = center_after_pan.y;
                let project =
                    |p: [f64; 3]| [cx + p[0] as f32 * scale, cy + p[1] as f32 * scale];

                self.orbit_hierarchy.tick(bodies, g_factor);

                // Prune stale pins (e.g. after a collision-merge swap_remove
                // on the sim thread invalidated the index).
                let n_bodies = bodies.len();
                self.pinned_orbits.retain(|&i| i < n_bodies);

                if self.show_orbit_ellipses {
                    let bg_style = OrbitOverlayStyle::background_default();
                    let vp_center = rect.center();
                    let vp_half_diag =
                        (rect.width().powi(2) + rect.height().powi(2)).sqrt() * 0.5;

                    // Collect (idx, primary_idx, elements, influence) for
                    // every body that survives the filter pipeline.
                    let mut candidates: Vec<(usize, usize, gravity_sim_core::physics::orbital::OrbitalElements, f32)> =
                        Vec::with_capacity(bodies.len().min(self.orbit_top_n * 2));

                    for i in 0..bodies.len() {
                        // Selected + pinned bodies draw in a dedicated
                        // foreground pass — skip here to avoid double-draw
                        // and to exempt them from all filter gates.
                        if Some(i) == self.selected_body {
                            continue;
                        }
                        if self.pinned_orbits.contains(&i) {
                            continue;
                        }
                        // Hierarchy must know i's primary; Free bodies
                        // have no elliptical orbit to draw.
                        let Some(level) = self.orbit_hierarchy.level(i) else {
                            continue;
                        };
                        if !level_is_visible(level, &self.orbit_visible_levels) {
                            continue;
                        }
                        let Some(primary_idx) = self.orbit_hierarchy.primary(i) else {
                            continue;
                        };
                        let Some(el) = compute_elements(bodies, i, primary_idx, g_factor)
                        else {
                            continue;
                        };
                        if self.orbit_hide_degenerate
                            && is_degenerate_orbit(&el, &bodies[primary_idx])
                        {
                            continue;
                        }

                        let b = &bodies[i];
                        let sp_x = cx + b.x as f32 * scale;
                        let sp_y = cy + b.y as f32 * scale;
                        let dx = sp_x - vp_center.x;
                        let dy = sp_y - vp_center.y;
                        let d_norm = (dx * dx + dy * dy).sqrt() / vp_half_diag.max(1.0);
                        // Bodies far beyond the viewport contribute nothing.
                        if d_norm > 2.0 {
                            continue;
                        }
                        let viewport_weight = 1.0 / (1.0 + d_norm * d_norm);
                        let mass_factor = (1.0_f32 + b.mass as f32).ln().max(0.0);
                        let influence = mass_factor * viewport_weight;

                        candidates.push((i, primary_idx, el, influence));
                    }

                    // Top-N by influence — partial_sort would be nicer
                    // but N is small and draws are the real cost.
                    candidates.sort_by(|a, b| {
                        b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal)
                    });
                    if candidates.len() > self.orbit_top_n {
                        candidates.truncate(self.orbit_top_n);
                    }

                    for (_i, primary_idx, el, _infl) in &candidates {
                        let primary = &bodies[*primary_idx];
                        let primary_pos = [primary.x, primary.y, 0.0];
                        let sampled = el.sample_orbit(primary_pos, 64);
                        draw_orbit_polyline(&mut backend, &sampled, project, &bg_style);
                        draw_orbit_apsides(&mut backend, el, primary_pos, project, &bg_style);
                    }
                }

                // Pinned orbits — drawn unconditionally in the foreground
                // (selected_default style). Intentionally placed outside the
                // `show_orbit_ellipses` guard so a pin is a user commitment
                // that overrides the global toggle too.
                if !self.pinned_orbits.is_empty() {
                    let fg_style = OrbitOverlayStyle::selected_default();
                    for &i in self.pinned_orbits.iter() {
                        if Some(i) == self.selected_body {
                            continue; // selected pass draws it
                        }
                        let primary = self
                            .orbit_hierarchy
                            .primary(i)
                            .or_else(|| dominant_primary(bodies, i));
                        let Some(primary_idx) = primary else {
                            continue;
                        };
                        let Some(el) = compute_elements(bodies, i, primary_idx, g_factor)
                        else {
                            continue;
                        };
                        let primary_b = &bodies[primary_idx];
                        let primary_pos = [primary_b.x, primary_b.y, 0.0];
                        let sampled = el.sample_orbit(primary_pos, 96);
                        draw_orbit_polyline(&mut backend, &sampled, project, &fg_style);
                        draw_orbit_apsides(&mut backend, &el, primary_pos, project, &fg_style);
                    }
                }

                if let Some(idx) = self.selected_body {
                    if idx < bodies.len() {
                        // Foreground overlay uses the hierarchy's primary
                        // when available; otherwise falls back to the raw
                        // dominant-primary heuristic so the first click on
                        // a body still shows its orbit before the
                        // hierarchy has ticked for this topology.
                        let primary = self
                            .orbit_hierarchy
                            .primary(idx)
                            .or_else(|| dominant_primary(bodies, idx));
                        if let Some(primary_idx) = primary {
                            if let Some(el) =
                                compute_elements(bodies, idx, primary_idx, g_factor)
                            {
                                let primary = &bodies[primary_idx];
                                let primary_pos = [primary.x, primary.y, 0.0];
                                let sampled = el.sample_orbit(primary_pos, 128);
                                let style = OrbitOverlayStyle::selected_default();
                                draw_orbit_polyline(
                                    &mut backend, &sampled, project, &style,
                                );
                                draw_orbit_apsides(
                                    &mut backend, &el, primary_pos, project, &style,
                                );
                            }
                        }
                    }
                }
            }

            backend.center = [center_after_pan.x, center_after_pan.y];
            backend.scale = self.scale;
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
                let dom_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);
                backend.trail_visibility = Some(
                    bodies
                        .iter()
                        .map(|b| dom_mass == 0.0 || b.mass / dom_mass >= self.trail_min_mass_ratio)
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
                    if self.find_body_at(pos, center_after_pan, render_params).is_none() {
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
                        // Velocity label
                        let v_scale = 0.5 / self.scale as f64;
                        let vx = delta.x as f64 * v_scale;
                        let vy = delta.y as f64 * v_scale;
                        let speed = (vx * vx + vy * vy).sqrt();
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
                    if self.find_body_at(cursor, center_after_pan, render_params).is_none() {
                        let spawn_pos = start.unwrap_or(cursor);
                        let wx = (spawn_pos.x - center_after_pan.x) as f64 / self.scale as f64;
                        let wy = (spawn_pos.y - center_after_pan.y) as f64 / self.scale as f64;

                        // Velocity from drag delta
                        let v_scale = 0.5 / self.scale as f64;
                        let (vx, vy) = if let Some(s) = start {
                            let d = cursor - s;
                            (d.x as f64 * v_scale, d.y as f64 * v_scale)
                        } else {
                            (0.0, 0.0)
                        };

                        use gravity_sim_core::domain::materials::density as mat_density;
                        let mut body = gravity_sim_core::domain::body::Body::new(
                            wx,
                            wy,
                            vx,
                            vy,
                            self.place_mass,
                            self.place_material,
                        );
                        body.density = mat_density(self.place_material, self.place_mass);
                        body.sync_physical_properties();

                        self.push_undo(UndoRecord::AddedBodies(1));
                        self.system.add_body(body);
                    }
                }
            }

            // Override cursor to crosshair while in place-mode (unless dragging a body)
            if self.dragging_body.is_none() {
                ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
            }

            ctx.request_repaint();
        } else {
            // ── Normal click: select + pan to center ──────────────────────────
            self.place_drag_start = None;
            if response.clicked() {
                if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                    match self.find_body_at(cursor, center_after_pan, render_params) {
                        Some(idx) => {
                            let body = self.system.bodies()[idx];
                            self.selected_body = Some(idx);
                            self.follow_selected_body = true;
                            let name = self.system.name(idx).to_owned();
                            self.selection_form = Some(SelectionForm::from_body(&body, &name));

                            self.pan_vel = egui::Vec2::ZERO;
                            self.zoom_vel = 0.0;

                            let screen_r =
                                compute_render_radius(body.physical_radius, render_params);
                            if screen_r < 6.0 && body.physical_radius > 1e-30 {
                                let desired_px = 24.0_f32;
                                let new_scale = (desired_px / body.physical_radius as f32)
                                    .clamp(self.scale * 2.0, self.scale * 500.0)
                                    .min(50_000.0);
                                self.scale = new_scale;
                            }

                            self.camera_anim_target = Some(egui::vec2(
                                -body.x as f32 * self.scale,
                                -body.y as f32 * self.scale,
                            ));
                        },
                        None => {
                            self.selected_body = None;
                            self.follow_selected_body = false;
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
                        // Convert screen pos → world pos
                        let wx = (screen_pos.x - center_after_pan.x) as f64 / self.scale as f64;
                        let wy = (screen_pos.y - center_after_pan.y) as f64 / self.scale as f64;
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
        if self.selected_body.is_some() {
            ctx.request_repaint();
        }

        // ── GPU paint callback ────────────────────────────────────────────────
        let device = self.device.as_ref().unwrap().clone();
        let queue = self.queue.as_ref().unwrap().clone();
        let format = self.format.unwrap();

        ui.painter().add(egui_wgpu::Callback::new_paint_callback(
            rect,
            CallbackFn {
                backend: self.backend.clone(),
                device,
                queue,
                format,
                // Canvas dimensions (not full window) so `to_ndc` maps correctly
                // into the canvas-rect viewport that egui_wgpu sets for callbacks.
                screen: [rect.width(), rect.height()],
                viewport_min: [rect.min.x, rect.min.y],
            },
        ));

        // ── Grid labels (AU ticks, on top of GPU grid, below body overlay) ───
        {
            let backend = self.backend.lock().unwrap();
            backend.draw_labels_overlay(ui, rect, &self.physics_cfg.dist_label);
        }

        // ── Overlay: rings + labels (on top of GPU layer) ─────────────────────
        self.draw_overlay(ui, center_after_pan, time);

        // ── Loading overlay ───────────────────────────────────────────────────
        if self.system.is_loading() {
            self.draw_loading_overlay(ui, rect, time);
            ctx.request_repaint();
        }

    }

    // ── Overlay ───────────────────────────────────────────────────────────────

    fn draw_overlay(&self, ui: &egui::Ui, center: Pos2, time: f32) {
        let bodies = self.system.bodies();
        let names = self.system.names();

        if bodies.is_empty() {
            return;
        }

        let max_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);

        // !! Use the SAME min_px as the GPU renderer so rings/labels are anchored
        //    to the visible disc, not to the (possibly sub-pixel) physical radius.
        let render_params = RenderParams {
            world_scale: self.scale,
            mode: self.semantic_scale_mode,
            min_px: match self.semantic_scale_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            },
        };

        // Label visibility: show when body is large enough on screen OR important
        // enough given current zoom. threshold ↑ at low zoom → fewer labels.
        let importance_threshold = (2.0_f64 / self.scale as f64).clamp(0.001, 1.0);

        let painter = ui.painter();
        let font = FontId::proportional(LABEL_FONT_SIZE);

        // Pulse for selection ring: ±1.5 px at ~3.5 Hz
        let pulse = (time * 3.5).sin() * 1.5_f32;

        for (i, (body, name)) in bodies.iter().zip(names.iter()).enumerate() {
            // visual_r = how big the body actually appears on screen (matches GPU)
            let visual_r = compute_render_radius(body.physical_radius, render_params);

            let px = center.x + body.x as f32 * self.scale;
            let py = center.y + body.y as f32 * self.scale;
            let body_pos = egui::pos2(px, py);

            let is_selected = self.selected_body == Some(i);
            let is_hovered = self.hovered_body == Some(i) && !is_selected;

            // ── Hover ring ───────────────────────────────────────────────
            if is_hovered {
                painter.circle_stroke(
                    body_pos,
                    visual_r + RING_GAP - 1.0,
                    Stroke::new(1.0, Color32::from_rgba_premultiplied(160, 160, 200, 90)),
                );
            }

            // ── Selection ring (pulsing) ─────────────────────────────────
            if is_selected {
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
            let show_label = visual_r >= 5.0 || importance >= importance_threshold || is_selected;

            if show_label {
                // Cap offset so the label never drifts far from the body
                let offset_y = (visual_r + 4.0).min(MAX_LABEL_OFFSET_PX);
                let label_pos = egui::pos2(px, py + offset_y);

                let color = if is_selected {
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
        center: Pos2,
        render_params: RenderParams,
    ) -> Option<usize> {
        let bodies = self.system.bodies();

        // Iterate in reverse so top-rendered (last) body wins ties
        for i in (0..bodies.len()).rev() {
            let b = &bodies[i];

            let px = center.x + b.x as f32 * self.scale;
            let py = center.y + b.y as f32 * self.scale;

            // Hit radius: visual size + generous minimum for easy clicking
            let r = compute_render_radius(b.physical_radius, render_params).max(MIN_HIT_PX);

            let dx = cursor.x - px;
            let dy = cursor.y - py;

            if dx * dx + dy * dy <= r * r {
                return Some(i);
            }
        }

        None
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

        let ctx = gravity_sim_core::domain::field::FieldContext {
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
