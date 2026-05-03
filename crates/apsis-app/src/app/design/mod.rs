//! Design system — single source of truth for visual and motion values.
//!
//! [`tokens`] holds every constant. Higher layers (primitives, theme bridge)
//! consume tokens; nothing else in the app references raw colours, sizes, or
//! durations directly.

pub mod fonts;
pub mod theme;
pub mod tokens;
