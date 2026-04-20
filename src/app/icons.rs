//! Phosphor icon aliases used throughout the UI.
//!
//! Keeping the full icon set behind named constants means a future swap
//! (e.g. to Lucide) touches only this file. Use `Regular` for inactive
//! states and `Fill` for active/selected states.

pub use egui_phosphor::fill;
pub use egui_phosphor::regular as reg;

// ── Tool rail ────────────────────────────────────────────────────────────────

pub const TOOL_OVERVIEW: &str = reg::GAUGE;
pub const TOOL_ADD: &str = reg::PLUS_CIRCLE;
pub const TOOL_TEMPLATES: &str = reg::STAR;
pub const TOOL_VIEW: &str = reg::EYE;
pub const TOOL_CAMERA: &str = reg::VIDEO_CAMERA;
pub const TOOL_CONFIG: &str = reg::SLIDERS_HORIZONTAL;

pub const TOOL_OVERVIEW_ON: &str = fill::GAUGE;
pub const TOOL_ADD_ON: &str = fill::PLUS_CIRCLE;
pub const TOOL_TEMPLATES_ON: &str = fill::STAR;
pub const TOOL_VIEW_ON: &str = fill::EYE;
pub const TOOL_CAMERA_ON: &str = fill::VIDEO_CAMERA;
pub const TOOL_CONFIG_ON: &str = fill::SLIDERS_HORIZONTAL;

// ── Top bar ──────────────────────────────────────────────────────────────────

pub const MENU: &str = reg::LIST;
pub const SIDEBAR_CLOSE: &str = reg::CARET_DOUBLE_LEFT;
pub const SIDEBAR_OPEN: &str = reg::CARET_DOUBLE_RIGHT;
pub const SETTINGS: &str = reg::GEAR;
pub const HELP: &str = reg::QUESTION;
pub const RECORD: &str = reg::RECORD;
pub const SAVE: &str = reg::FLOPPY_DISK;
pub const LOAD: &str = reg::FOLDER_OPEN;
pub const CLEAR: &str = reg::TRASH;

// ── Playbar ──────────────────────────────────────────────────────────────────

pub const PLAY: &str = reg::PLAY;
pub const PAUSE: &str = reg::PAUSE;
pub const RESET: &str = reg::ARROW_COUNTER_CLOCKWISE;
pub const STEP: &str = reg::SKIP_FORWARD;

// ── Toggles ──────────────────────────────────────────────────────────────────

pub const CHECK_ON: &str = fill::CHECK_SQUARE;
pub const CHECK_OFF: &str = reg::SQUARE;
pub const RADIO_ON: &str = fill::CHECK_CIRCLE;
pub const RADIO_OFF: &str = reg::CIRCLE;

// ── Inspector / misc ────────────────────────────────────────────────────────

pub const FIT_VIEW: &str = reg::FRAME_CORNERS;
pub const TARGET: &str = reg::TARGET;
