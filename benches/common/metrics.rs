//! Scenario metrics: what a single run produces and how to classify it.
//!
//! Metrics are split into two tiers that drive different baseline
//! tolerance strategies:
//!
//! | Tier      | Examples                          | Default tolerance |
//! |-----------|-----------------------------------|-------------------|
//! | Counter   | substeps, rejections, iters       | `tol_abs = 0`     |
//! | Float     | peak_energy_err, dt stats         | `tol_factor` from 10-run jitter |
//!
//! Counters are expected to be bit-deterministic on a single machine
//! (the integrator is deterministic; rayon is forced single-thread
//! in the harness). Any non-zero delta is a real behavioural change
//! and must be reviewed, not absorbed.
//!
//! Floats are near-deterministic but may flutter by a few ULPs due to
//! fused-multiply-add choices, compiler reassociation, or threshold
//! crossings on the last digits. The recording pass observes the
//! actual jitter across N runs and sizes `tol_factor` to 2× observed
//! — tight by default, but adaptive to real-world noise.

/// Classification that determines the default tolerance strategy when
/// recording a baseline. See the module-level docs for rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricTier {
    /// Discrete, expected bit-deterministic. Recording asserts
    /// `min == max` across runs; otherwise the harness aborts with
    /// a clear error rather than silently widening tolerance.
    Counter,
    /// Floating-point; `tol_factor` derived from observed jitter.
    Float,
}

/// Snapshot of one scenario run. Field names are the stable keys used
/// in the baseline file — renaming requires a baseline update.
#[derive(Debug, Clone)]
pub struct ScenarioMetrics {
    // ── Tier 1: controller counters ──────────────────────────────────
    pub substeps: u64,
    pub rejections_picard: u64,
    pub rejections_truncation: u64,
    pub picard_iters_total: u64,
    pub degraded_total: u64,

    // ── Tier 1: dt profile (controller behaviour) ────────────────────
    pub dt_min: f64,
    pub dt_max: f64,
    pub dt_mean: f64,
    pub dt_p50: f64,
    pub dt_p95: f64,

    // ── Tier 2: numerical quality ────────────────────────────────────
    pub peak_energy_err: f64,
}

impl ScenarioMetrics {
    /// Ordered list of `(name, tier)` for all metrics. The order here
    /// drives the order of entries written to the baseline TOML file.
    /// Keep it stable — reordering would noise up diffs without any
    /// semantic change.
    pub const ALL: &'static [(&'static str, MetricTier)] = &[
        ("substeps", MetricTier::Counter),
        ("rejections_picard", MetricTier::Counter),
        ("rejections_truncation", MetricTier::Counter),
        ("picard_iters_total", MetricTier::Counter),
        ("degraded_total", MetricTier::Counter),
        // dt stats are floats by type but derive from a deterministic
        // sequence of accepted sub-step sizes — treated as Tier 2 so
        // the baseline captures per-run jitter (there should be none,
        // but we measure rather than assume).
        ("dt_min", MetricTier::Float),
        ("dt_max", MetricTier::Float),
        ("dt_mean", MetricTier::Float),
        ("dt_p50", MetricTier::Float),
        ("dt_p95", MetricTier::Float),
        ("peak_energy_err", MetricTier::Float),
    ];

    /// Look up a metric by its stable name. Returns `None` only when
    /// the caller has passed a typo — treat that as a programmer error.
    pub fn get(&self, name: &str) -> Option<f64> {
        match name {
            "substeps" => Some(self.substeps as f64),
            "rejections_picard" => Some(self.rejections_picard as f64),
            "rejections_truncation" => Some(self.rejections_truncation as f64),
            "picard_iters_total" => Some(self.picard_iters_total as f64),
            "degraded_total" => Some(self.degraded_total as f64),
            "dt_min" => Some(self.dt_min),
            "dt_max" => Some(self.dt_max),
            "dt_mean" => Some(self.dt_mean),
            "dt_p50" => Some(self.dt_p50),
            "dt_p95" => Some(self.dt_p95),
            "peak_energy_err" => Some(self.peak_energy_err),
            _ => None,
        }
    }
}

/// Build a [`ScenarioMetrics`] from a vector of `consumed_dt` samples
/// (one per accepted sub-step), the integrator's adaptive stats, and
/// the peak relative energy error observed during the run.
///
/// Kept as a free function rather than a constructor because it's the
/// one place that touches the runtime representation (`AdaptiveStats`),
/// so the struct itself stays as a pure data type.
pub fn assemble(
    dt_samples: &[f64],
    stats: &gravity_sim::physics::integrator::traits::AdaptiveStats,
    peak_energy_err: f64,
) -> ScenarioMetrics {
    let (dt_min, dt_max, dt_mean, dt_p50, dt_p95) = dt_summary(dt_samples);
    ScenarioMetrics {
        substeps: stats.substeps,
        rejections_picard: stats.rejections_picard,
        rejections_truncation: stats.rejections_truncation,
        picard_iters_total: stats.picard_iters,
        degraded_total: stats.degraded,
        dt_min,
        dt_max,
        dt_mean,
        dt_p50,
        dt_p95,
        peak_energy_err,
    }
}

/// Return `(min, max, mean, p50, p95)` for a slice of dt samples.
/// Empty input returns all zeros — a scenario that produced zero
/// substeps is degenerate and will fail other baseline checks first.
fn dt_summary(samples: &[f64]) -> (f64, f64, f64, f64, f64) {
    if samples.is_empty() {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let mut sorted: Vec<f64> = samples.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("dt samples contain NaN"));

    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let mean = sorted.iter().sum::<f64>() / sorted.len() as f64;
    let p50 = percentile(&sorted, 0.50);
    let p95 = percentile(&sorted, 0.95);
    (min, max, mean, p50, p95)
}

/// Nearest-rank percentile on a pre-sorted slice. Coarser than linear
/// interpolation but deterministic and monotonic in sample size — which
/// matters more than smoothness for regression detection.
fn percentile(sorted: &[f64], q: f64) -> f64 {
    debug_assert!((0.0..=1.0).contains(&q));
    let idx = ((sorted.len() as f64 - 1.0) * q).round() as usize;
    sorted[idx]
}
