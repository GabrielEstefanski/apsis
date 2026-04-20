//! Shared drift-severity classification used by the Overview tab, the
//! playbar stability badge, and the Display / Physics surfaces.
//!
//! The old metrics module also rendered two full diagnostic tables
//! (`panel_metrics_compact` + `panel_diagnostics_detail`) when Advanced
//! still hosted the diagnostics section. After the F3 reorg, every
//! diagnostic value lives in the Overview tab, so only the severity
//! classification remains here.

use crate::app::theme::{ACCENT, DANGER, SUCCESS, TEXT_DIM};
use eframe::egui::Color32;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::app::panel) enum DriftSeverity {
    Excellent,
    Good,
    Acceptable,
    Warning,
    Critical,
}

impl DriftSeverity {
    pub(in crate::app::panel) fn from_peak(peak: f64) -> Self {
        let p = peak.abs();
        if p < 1e-9 { Self::Excellent }
        else if p < 1e-6 { Self::Good }
        else if p < 1e-3 { Self::Acceptable }
        else if p < 1e-1 { Self::Warning }
        else { Self::Critical }
    }

    pub(in crate::app::panel) fn color(self) -> Color32 {
        match self {
            Self::Excellent | Self::Good => SUCCESS,
            Self::Acceptable => TEXT_DIM,
            Self::Warning => ACCENT,
            Self::Critical => DANGER,
        }
    }

    pub(in crate::app::panel) fn dot(self) -> &'static str {
        match self {
            Self::Excellent | Self::Good | Self::Acceptable => "●",
            Self::Warning | Self::Critical => "▲",
        }
    }

    pub(in crate::app::panel) fn label(self) -> &'static str {
        match self {
            Self::Excellent => "excellent",
            Self::Good => "good",
            Self::Acceptable => "acceptable",
            Self::Warning => "warning",
            Self::Critical => "critical",
        }
    }
}
