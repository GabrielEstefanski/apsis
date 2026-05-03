//! Design tokens — colours, typography, spacing, border, motion, icon,
//! shape, z-index, and UI scale.
//!
//! Every visual constant in the app sources from this module. Token names
//! mirror the design spec (`color::bg::CANVAS`, `typography::text::BASE`,
//! `space::S5`, etc.) so reviews against the spec are 1:1.
//!
//! Spacing follows a four-pixel grid (multiples of four only). The font
//! size scale skips 12/16/20 by convention. Border widths are device-pixel
//! exact and ignore [`scale`] — they stay crisp at every UI scale.

use egui::Color32;
use std::time::Duration;

/// Colour tokens grouped by semantic role.
pub mod color {
    use super::Color32;

    /// Background levels — five graphite layers, hue ≈ 220, saturation ≤ 5%.
    pub mod bg {
        use super::Color32;
        pub const CANVAS: Color32 = Color32::from_rgb(11, 12, 14);
        pub const SURFACE: Color32 = Color32::from_rgb(18, 19, 23);
        pub const RAISED: Color32 = Color32::from_rgb(25, 27, 32);
        pub const HOVER: Color32 = Color32::from_rgb(31, 34, 40);
        pub const POPOVER: Color32 = Color32::from_rgb(38, 42, 49);
    }

    /// Foreground levels — three contrast tiers.
    pub mod fg {
        use super::Color32;
        pub const PRIMARY: Color32 = Color32::from_rgb(232, 230, 225);
        pub const SECONDARY: Color32 = Color32::from_rgb(171, 171, 165);
        pub const TERTIARY: Color32 = Color32::from_rgb(110, 111, 107);
    }

    /// Single accent (copper-amber) reserved for selection, focus, active state.
    pub mod accent {
        use super::Color32;
        /// Solid accent — selection ring, focus ring, active tab, active tool.
        pub const SOLID: Color32 = Color32::from_rgb(200, 166, 117);
        /// Hairline divider, premultiplied α 0.15.
        pub const HAIRLINE: Color32 = Color32::from_rgba_premultiplied(30, 25, 18, 38);
        /// Bloomberg-flash on changed numeric values, premultiplied α 0.10.
        pub const FLASH: Color32 = Color32::from_rgba_premultiplied(20, 17, 12, 26);
    }

    /// Signal colours — physical state semantics, never decoration.
    pub mod signal {
        use super::Color32;
        pub const LIVE: Color32 = Color32::from_rgb(111, 160, 160);
        pub const WARN: Color32 = Color32::from_rgb(212, 149, 92);
        pub const ERROR: Color32 = Color32::from_rgb(184, 90, 101);
        pub const OK: Color32 = Color32::from_rgb(125, 165, 137);
    }
}

/// Typography tokens.
pub mod typography {
    /// Font family identifiers. Loaded into [`egui::FontDefinitions`] at app entry.
    pub mod font {
        pub const SANS_REGULAR: &str = "PlexSans-Regular";
        pub const SANS_MEDIUM: &str = "PlexSans-Medium";
        pub const MONO_REGULAR: &str = "PlexMono-Regular";
    }

    /// Font sizes (pixels).
    pub mod text {
        pub const XS: f32 = 11.0;
        pub const SM: f32 = 13.0;
        pub const BASE: f32 = 14.0;
        pub const MD: f32 = 15.0;
        pub const LG: f32 = 18.0;
        pub const XL: f32 = 22.0;
    }

    /// Line heights (pixels).
    pub mod lh {
        pub const TIGHT: f32 = 14.0;
        pub const DENSE: f32 = 16.0;
        pub const PROSE: f32 = 20.0;
    }

    /// Letter spacing (pixels).
    pub mod track {
        pub const NORMAL: f32 = 0.0;
        pub const WIDE: f32 = 0.6;
    }
}

/// Spacing scale (pixels) — multiples of four only.
pub mod space {
    pub const S0: f32 = 0.0;
    pub const S1: f32 = 2.0;
    pub const S2: f32 = 4.0;
    pub const S3: f32 = 8.0;
    pub const S4: f32 = 12.0;
    pub const S5: f32 = 16.0;
    pub const S6: f32 = 24.0;
    pub const S7: f32 = 32.0;
    pub const S8: f32 = 48.0;
    pub const S9: f32 = 64.0;
}

/// Border widths and radii. Widths are device-pixel exact; UI [`scale`] does
/// not affect them — they remain 1 px at every zoom level.
pub mod border {
    pub mod width {
        pub const HAIRLINE: f32 = 1.0;
        pub const REGULAR: f32 = 1.0;
    }

    pub mod radius {
        pub const NONE: f32 = 0.0;
        pub const SHARP: f32 = 2.0;
    }
}

/// Motion tokens.
pub mod motion {
    use std::time::Duration;

    pub const FAST: Duration = Duration::from_millis(90);
    pub const DEFAULT: Duration = Duration::from_millis(140);
    pub const SLOW: Duration = Duration::from_millis(220);

    /// Single easing curve in the system: ease-out snap.
    pub const EASING: CubicBezier = CubicBezier { x1: 0.2, y1: 0.0, x2: 0.0, y2: 1.0 };

    pub struct CubicBezier {
        pub x1: f32,
        pub y1: f32,
        pub x2: f32,
        pub y2: f32,
    }
}

/// Icon sizes (pixels) for Phosphor glyphs. Stroke width comes from the
/// Phosphor variant — Regular = 1.5 px.
pub mod icon {
    pub const SM: f32 = 12.0;
    pub const BASE: f32 = 16.0;
    pub const LG: f32 = 20.0;
    pub const XL: f32 = 24.0;
}

/// Primitive shape sizes (pixels) — body colour swatch, status dots, etc.
/// Rendered directly via [`egui::Painter`], not Phosphor.
pub mod shape {
    pub const SWATCH: f32 = 8.0;
    pub const DOT_LIVE: f32 = 6.0;
    pub const DOT_IDLE: f32 = 6.0;
}

/// Z-index for layered surfaces.
pub mod z {
    pub const BASE: i32 = 0;
    pub const SURFACE: i32 = 1;
    pub const RAISED: i32 = 2;
    pub const POPOVER: i32 = 3;
    pub const MODAL: i32 = 4;
    pub const TOAST: i32 = 5;
}

/// UI scale factors. Multiplies [`typography::text`] sizes and [`space`]
/// tokens. Does not multiply [`border::width`], [`icon`] strokes, or any
/// hairline — those remain device-pixel.
pub mod scale {
    pub const COMPACT: f32 = 0.9;
    pub const DEFAULT: f32 = 1.0;
    pub const COMFORTABLE: f32 = 1.1;
    pub const LARGE: f32 = 1.2;
}
