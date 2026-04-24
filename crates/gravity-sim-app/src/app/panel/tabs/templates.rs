use crate::app::icons;
use crate::app::theme::{BORDER, PANEL_BG, SURFACE_CARD, TEXT_DIM, TEXT_PRI, TEXT_SEC, section};
use crate::app::ui::{SimulationApp, UndoRecord};
use eframe::egui::text::{LayoutJob, TextFormat};
use eframe::egui::{self, Align, Color32, FontId, Frame, Margin, RichText, Stroke};
use gravity_sim_core::templates::{TEMPLATES, TemplateCategory, instantiate_at};

impl SimulationApp {
    // ── Sidebar tab ──────────────────────────────────────────────────────────

    pub(super) fn panel_tab_templates(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        if crate::app::theme::primary_btn(
            ui,
            &format!("{}  Browse all templates", icons::TOOL_TEMPLATES),
        ) {
            self.show_templates_modal = true;
        }

        ui.add_space(3.0);
        ui.label(
            RichText::new("Click to spawn  ·  Drag to place on canvas")
                .size(9.5)
                .color(TEXT_DIM)
                .italics(),
        );

        section(ui, "FEATURED");

        const FEATURED: &[&str] = &["Solar System", "Binary Stars", "TRAPPIST-1"];

        let mut clicked: Option<usize> = None;
        let mut dragged: Option<usize> = None;

        for (global_idx, entry) in
            TEMPLATES.iter().enumerate().filter(|(_, e)| FEATURED.contains(&e.name))
        {
            let (desc, count) = self.templates_meta[global_idx];
            let resp = small_card(ui, entry.name, desc, count);
            if resp.clicked() {
                clicked = Some(global_idx);
            }
            if resp.drag_started() {
                dragged = Some(global_idx);
            }
            ui.add_space(4.0);
        }

        self.apply_template_action(clicked, dragged);

        ui.add_space(2.0);
        let remaining = TEMPLATES.len().saturating_sub(FEATURED.len());
        if remaining > 0 {
            ui.label(
                RichText::new(format!("… {remaining} more in catalog"))
                    .size(9.5)
                    .color(TEXT_DIM)
                    .italics(),
            );
        }
    }

    // ── Full modal catalog ───────────────────────────────────────────────────

    pub(in crate::app) fn draw_templates_modal(&mut self, ctx: &egui::Context) {
        if !self.show_templates_modal {
            return;
        }

        let screen = ctx.screen_rect();
        let modal_w = (screen.width() * 0.55).clamp(460.0, 680.0);
        let modal_h = (screen.height() * 0.78).clamp(380.0, 600.0);

        let query = self.templates_search.trim().to_lowercase();

        // Card index list built each frame — cheap: just usize + &'static refs
        let visible: Vec<(usize, TemplateCategory, &'static str, &'static str, usize)> = TEMPLATES
            .iter()
            .enumerate()
            .filter_map(|(i, e)| {
                let (desc, count) = self.templates_meta[i];
                if query.is_empty()
                    || e.name.to_lowercase().contains(&query)
                    || desc.to_lowercase().contains(&query)
                {
                    Some((i, e.category, e.name, desc, count))
                } else {
                    None
                }
            })
            .collect();

        let mut clicked: Option<usize> = None;
        let mut dragged: Option<usize> = None;
        let mut close = false;

        egui::Window::new("template_catalog")
            .id(egui::Id::new("templates_modal"))
            .title_bar(false)
            .collapsible(false)
            .resizable(false)
            .fixed_size([modal_w, modal_h])
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .frame(
                Frame::NONE
                    .fill(PANEL_BG)
                    .stroke(Stroke::new(1.0, BORDER))
                    .corner_radius(10.0)
                    .inner_margin(Margin::same(16)),
            )
            .show(ctx, |ui| {
                // ── Header ──────────────────────────────────────────────── //
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Templates").size(15.0).strong().color(TEXT_PRI));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new(RichText::new("×").size(18.0).color(TEXT_DIM))
                                    .fill(Color32::TRANSPARENT)
                                    .frame(false),
                            )
                            .on_hover_text("Close (Esc)")
                            .clicked()
                        {
                            close = true;
                        }
                    });
                });

                // Close on Escape
                if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
                    close = true;
                }

                ui.add_space(8.0);

                // ── Search ──────────────────────────────────────────────── //
                ui.horizontal(|ui| {
                    ui.label(RichText::new(icons::SEARCH).size(11.0).color(TEXT_DIM));
                    ui.add_space(2.0);
                    ui.add(
                        egui::TextEdit::singleline(&mut self.templates_search)
                            .desired_width(ui.available_width())
                            .hint_text("Filter by name or description…")
                            .font(FontId::proportional(11.0))
                            .text_color(TEXT_PRI),
                    );
                });

                ui.add_space(4.0);
                ui.label(
                    RichText::new("Click to spawn at origin  ·  Drag to place on canvas")
                        .size(9.5)
                        .color(TEXT_DIM)
                        .italics(),
                );
                ui.add_space(6.0);
                ui.add(egui::Separator::default().horizontal());
                ui.add_space(4.0);

                // ── Scrollable list ──────────────────────────────────────── //
                egui::ScrollArea::vertical().auto_shrink([false; 2]).show(ui, |ui| {
                    if visible.is_empty() {
                        ui.add_space(24.0);
                        ui.label(
                            RichText::new("No templates match your search.")
                                .size(11.0)
                                .color(TEXT_DIM)
                                .italics(),
                        );
                        return;
                    }

                    for cat in [
                        TemplateCategory::Bodies,
                        TemplateCategory::Systems,
                        TemplateCategory::ThreeBodyProblems,
                    ] {
                        let group: Vec<_> = visible.iter().filter(|(_, c, ..)| *c == cat).collect();
                        if group.is_empty() {
                            continue;
                        }

                        section(ui, cat.label());

                        for (global_idx, _, name, desc, count) in &group {
                            let resp = full_card(ui, name, desc, *count);
                            if resp.clicked() {
                                clicked = Some(*global_idx);
                            }
                            if resp.drag_started() {
                                dragged = Some(*global_idx);
                            }
                            ui.add_space(5.0);
                        }

                        ui.add_space(6.0);
                    }
                });
            });

        if close {
            self.show_templates_modal = false;
        }

        self.apply_template_action(clicked, dragged);
        if clicked.is_some() || dragged.is_some() {
            self.show_templates_modal = false;
        }
    }

    // ── Shared spawn logic ───────────────────────────────────────────────────

    fn apply_template_action(&mut self, clicked: Option<usize>, dragged: Option<usize>) {
        if let Some(idx) = clicked {
            let entry = &TEMPLATES[idx];
            let preview = entry.build(self.system.seed());
            let bodies = instantiate_at(&preview, 0.0, 0.0);
            self.push_undo(UndoRecord::AddedBodies(bodies.len()));
            self.system.add_named_bodies(bodies);
            self.pending_fit = true;
            self.reset_drift_peaks();
            if self.sim_name.is_empty() {
                self.sim_name = entry.name.to_owned();
            }
        }

        if let Some(idx) = dragged {
            let seed = self.system.seed();
            let kind = TEMPLATES[idx].kind;
            self.template_drag = Some(Box::new(move || kind.build(seed)));
        }
    }
}

// ── Card widgets ─────────────────────────────────────────────────────────────

/// Full card for the modal: large title + description body + count chip top-right.
fn full_card(ui: &mut egui::Ui, name: &str, description: &str, count: usize) -> egui::Response {
    let mut job = LayoutJob::default();
    job.append(
        name,
        0.0,
        TextFormat {
            font_id: FontId::proportional(12.5),
            color: TEXT_PRI,
            valign: Align::Center,
            ..Default::default()
        },
    );
    job.append(
        &format!("\n{description}"),
        0.0,
        TextFormat { font_id: FontId::proportional(10.0), color: TEXT_SEC, ..Default::default() },
    );

    let w = ui.available_width();
    let resp = ui.add(
        egui::Button::new(job)
            .fill(SURFACE_CARD)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(w, 54.0))
            .corner_radius(7.0)
            .sense(egui::Sense::click_and_drag()),
    );

    ui.painter().text(
        resp.rect.right_top() + egui::vec2(-10.0, 12.0),
        egui::Align2::RIGHT_TOP,
        &format!("{count} bodies"),
        FontId::proportional(9.5),
        TEXT_DIM,
    );

    resp
}

/// Compact card for the sidebar: name + count chip only.
fn small_card(ui: &mut egui::Ui, name: &str, _desc: &str, count: usize) -> egui::Response {
    let w = ui.available_width();
    let resp = ui.add(
        egui::Button::new(RichText::new(name).size(11.0).color(TEXT_PRI))
            .fill(SURFACE_CARD)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(w, 28.0))
            .corner_radius(6.0)
            .sense(egui::Sense::click_and_drag()),
    );

    ui.painter().text(
        resp.rect.right_center() + egui::vec2(-10.0, 0.0),
        egui::Align2::RIGHT_CENTER,
        &format!("{count}"),
        FontId::proportional(9.5),
        TEXT_DIM,
    );

    resp
}
