//! [`ColorView`] — the composition that turns per-body scalars into per-body
//! RGB triples. This is the single value the UI stores as "the active
//! visualization mode" when the user opts in to data-driven colouring.
//!
//! When `ColorView` is `None` (the default), the rest of the render
//! pipeline falls back to each body's material colour — preserving the
//! existing Universe-Sandbox-style categorical palette.

use apsis::domain::field::FieldContext;

/// Selection of which field, normalizer, and colormap are in use. Stored
/// by ID rather than by reference so the app can serialise the choice
/// without tying the UI to the trait-object lifetimes.
#[derive(Clone, Debug)]
pub struct ColorViewSelection {
    pub field_id: String,
    pub normalizer_id: String,
    pub colormap_id: String,
    /// Optional explicit data range. `None` means "auto-detect from the
    /// current body set every frame".
    pub range: Option<(f64, f64)>,
}

impl ColorViewSelection {
    /// Default pick when the user first enables data-driven colouring:
    /// velocity × linear × cool-warm (Universe-Sandbox hot/cold feel).
    pub fn default_velocity() -> Self {
        Self {
            field_id: "velocity".into(),
            normalizer_id: "linear".into(),
            colormap_id: "cool_warm".into(),
            range: None,
        }
    }
}

/// Output of one `compute()` call: an RGB triple per body plus the data
/// range actually used (resolved from `selection.range` or auto-detected).
/// The range is exposed so the UI can render a colour bar.
pub struct ColorViewOutput {
    pub colors: Vec<[u8; 3]>,
    pub resolved_range: (f64, f64),
}

/// Evaluate a [`ColorViewSelection`] against the current frame.
///
/// `field`, `normalizer`, `colormap` are resolved by the caller from the
/// three registries — this function is pure data transformation.
pub fn compute(
    field: &dyn apsis::domain::field::BodyField,
    normalizer: &dyn super::normalizer::Normalizer,
    colormap: &dyn super::colormap::Colormap,
    explicit_range: Option<(f64, f64)>,
    ctx: &FieldContext,
) -> ColorViewOutput {
    let n = ctx.bodies.len();

    // First pass: sample every body. Kept in a scratch vec so the range
    // auto-detection and the normalize pass share one traversal each.
    let mut samples: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        samples.push(field.sample(i, ctx));
    }

    let range = explicit_range.unwrap_or_else(|| auto_range(&samples));

    let mut colors = Vec::with_capacity(n);
    for &v in samples.iter() {
        let t = normalizer.normalize(v, range);
        colors.push(colormap.sample(t));
    }

    ColorViewOutput { colors, resolved_range: range }
}

/// Min/max over finite samples, with safe fallbacks for degenerate cases.
fn auto_range(samples: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in samples {
        if v.is_finite() {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
    }
    if !lo.is_finite() || !hi.is_finite() {
        return (0.0, 1.0);
    }
    if (hi - lo).abs() < 1e-300 {
        return (lo - 0.5, hi + 0.5);
    }
    (lo, hi)
}
