//! Font registration — embeds IBM Plex Sans (Regular + Medium) and IBM Plex
//! Mono Regular into the egui font registry. Plex Sans becomes the primary
//! Proportional face; Plex Mono becomes the primary Monospace face. The
//! Medium weight is registered as a named family for explicit weight
//! requests by primitives that render section headings.
//!
//! License: SIL Open Font License 1.1. See `assets/fonts/LICENSE.txt`.

use eframe::egui::{FontData, FontDefinitions, FontFamily};
use std::sync::Arc;

const PLEX_SANS_REGULAR: &[u8] = include_bytes!("../../../assets/fonts/IBMPlexSans-Regular.ttf");
const PLEX_SANS_MEDIUM: &[u8] = include_bytes!("../../../assets/fonts/IBMPlexSans-Medium.ttf");
const PLEX_MONO_REGULAR: &[u8] = include_bytes!("../../../assets/fonts/IBMPlexMono-Regular.ttf");

/// Register Plex faces into the given font definitions. Call before
/// `egui_phosphor::add_to_fonts` so the latter inherits the same registry.
pub fn install(fonts: &mut FontDefinitions) {
    use super::tokens::typography::font;

    fonts
        .font_data
        .insert(font::SANS_REGULAR.to_owned(), Arc::new(FontData::from_static(PLEX_SANS_REGULAR)));
    fonts
        .font_data
        .insert(font::SANS_MEDIUM.to_owned(), Arc::new(FontData::from_static(PLEX_SANS_MEDIUM)));
    fonts
        .font_data
        .insert(font::MONO_REGULAR.to_owned(), Arc::new(FontData::from_static(PLEX_MONO_REGULAR)));

    fonts
        .families
        .entry(FontFamily::Proportional)
        .or_default()
        .insert(0, font::SANS_REGULAR.to_owned());

    fonts
        .families
        .entry(FontFamily::Monospace)
        .or_default()
        .insert(0, font::MONO_REGULAR.to_owned());

    fonts.families.insert(
        FontFamily::Name(font::SANS_MEDIUM.into()),
        vec![font::SANS_MEDIUM.to_owned(), font::SANS_REGULAR.to_owned()],
    );
}
