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

    // ── Tier 2: dt profile (controller behaviour) ────────────────────
    //
    // Includes `dt_p05` in addition to the median (`dt_p50`) because a
    // controller that only subtly relaxes its lower-tail behaviour
    // (think: accepting slightly larger steps through pericenter) can
    // leave `dt_min` unchanged while the aggressive-regime population
    // shifts upward. Reading `p05` and `p50` together makes shape
    // changes visible that a single percentile would hide.
    pub dt_min: f64,
    pub dt_max: f64,
    pub dt_mean: f64,
    pub dt_p05: f64,
    pub dt_p50: f64,
    pub dt_p95: f64,

    // ── Tier 2: numerical quality ────────────────────────────────────
    //
    // `peak_energy_err` catches gross blow-ups but is oscillatory in
    // nature for well-behaved scenes — it can stay flat while energy
    // silently drifts. The extra two metrics close that gap:
    //
    //   * `rel_energy_err_rms`: sqrt(mean(err²)) over all samples.
    //     Penalises *sustained* error rather than isolated spikes.
    //     Useful as a scalar summary of integration quality over the
    //     whole window.
    //   * `energy_drift_slope`: least-squares slope of |err(t)| vs t.
    //     > 0 means the absolute-error envelope is growing (secular
    //     drift, bad). ≈ 0 means the error is bounded / oscillatory
    //     (good — what an adaptive high-order integrator should
    //     produce for non-chaotic scenes). We use |err| rather than
    //     signed err so an oscillation around zero can't mask drift
    //     as "zero slope".
    //
    // Units on the slope are 1/time, so values are scenario-dependent
    // and only meaningful when compared against the same scenario's
    // baseline. That's how regression is checked anyway.
    pub peak_energy_err: f64,
    pub rel_energy_err_rms: f64,
    pub energy_drift_slope: f64,
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
        ("dt_p05", MetricTier::Float),
        ("dt_p50", MetricTier::Float),
        ("dt_p95", MetricTier::Float),
        // Numerical quality metrics. Tier::Float because they involve
        // summations / sqrt / division that can flutter at ULP level
        // across platforms even when the underlying integration is
        // bit-deterministic. On a single machine the recording pass
        // typically derives tol_factor = 1.0 anyway.
        ("peak_energy_err", MetricTier::Float),
        ("rel_energy_err_rms", MetricTier::Float),
        ("energy_drift_slope", MetricTier::Float),
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
            "dt_p05" => Some(self.dt_p05),
            "dt_p50" => Some(self.dt_p50),
            "dt_p95" => Some(self.dt_p95),
            "peak_energy_err" => Some(self.peak_energy_err),
            "rel_energy_err_rms" => Some(self.rel_energy_err_rms),
            "energy_drift_slope" => Some(self.energy_drift_slope),
            _ => None,
        }
    }
}

/// Parallel per-substep samples collected during a scenario run.
///
/// All three vectors are indexed by accepted sub-step: `t[i]` is the
/// simulation time at which sub-step `i` ended, `dt[i]` its size, and
/// `abs_energy_err[i]` the absolute relative energy error at that
/// point. The three lengths must match — [`assemble`] asserts this
/// rather than silently aligning to the shorter one.
#[derive(Debug, Default)]
pub struct RunSamples {
    pub t: Vec<f64>,
    pub dt: Vec<f64>,
    pub abs_energy_err: Vec<f64>,
}

impl RunSamples {
    /// Pre-allocate all three vectors to the same capacity. Tight
    /// upper bounds on substep count are easy to estimate from
    /// scenario duration / dt_budget; getting the capacity right
    /// eliminates the realloc chain that would otherwise land inside
    /// the validation hot path.
    pub fn with_capacity(expected_substeps: usize) -> Self {
        Self {
            t: Vec::with_capacity(expected_substeps),
            dt: Vec::with_capacity(expected_substeps),
            abs_energy_err: Vec::with_capacity(expected_substeps),
        }
    }

    /// Append one sub-step sample. All three vectors grow in lock-step.
    pub fn push(&mut self, t: f64, dt: f64, abs_energy_err: f64) {
        self.t.push(t);
        self.dt.push(dt);
        self.abs_energy_err.push(abs_energy_err);
    }
}

/// Build a [`ScenarioMetrics`] from a batch of per-sub-step samples
/// and the integrator's adaptive stats.
///
/// Kept as a free function rather than a constructor because this is
/// the one place that touches the runtime representation
/// (`AdaptiveStats`), so the struct itself stays a pure data type.
pub fn assemble(
    samples: &RunSamples,
    stats: &gravity_sim_core::physics::integrator::traits::AdaptiveStats,
) -> ScenarioMetrics {
    assert_eq!(
        samples.t.len(),
        samples.dt.len(),
        "RunSamples.t and .dt must have identical lengths"
    );
    assert_eq!(
        samples.t.len(),
        samples.abs_energy_err.len(),
        "RunSamples.t and .abs_energy_err must have identical lengths"
    );

    let dt = dt_summary(&samples.dt);

    let peak_energy_err = samples.abs_energy_err.iter().copied().fold(0.0_f64, f64::max);
    let rel_energy_err_rms = rms(&samples.abs_energy_err);
    let energy_drift_slope = linear_regression_slope(&samples.t, &samples.abs_energy_err);

    ScenarioMetrics {
        substeps: stats.substeps,
        rejections_picard: stats.rejections_picard,
        rejections_truncation: stats.rejections_truncation,
        picard_iters_total: stats.picard_iters,
        degraded_total: stats.degraded,
        dt_min: dt.min,
        dt_max: dt.max,
        dt_mean: dt.mean,
        dt_p05: dt.p05,
        dt_p50: dt.p50,
        dt_p95: dt.p95,
        peak_energy_err,
        rel_energy_err_rms,
        energy_drift_slope,
    }
}

// ── Distribution helpers ─────────────────────────────────────────────────────

struct DtSummary {
    min: f64,
    max: f64,
    mean: f64,
    p05: f64,
    p50: f64,
    p95: f64,
}

/// Summary stats for a slice of dt samples. Empty input returns all
/// zeros — a scenario that produced zero substeps is degenerate and
/// will fail other baseline checks first.
fn dt_summary(samples: &[f64]) -> DtSummary {
    if samples.is_empty() {
        return DtSummary { min: 0.0, max: 0.0, mean: 0.0, p05: 0.0, p50: 0.0, p95: 0.0 };
    }
    let mut sorted: Vec<f64> = samples.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("dt samples contain NaN"));

    DtSummary {
        min: sorted[0],
        max: sorted[sorted.len() - 1],
        mean: sorted.iter().sum::<f64>() / sorted.len() as f64,
        p05: percentile(&sorted, 0.05),
        p50: percentile(&sorted, 0.50),
        p95: percentile(&sorted, 0.95),
    }
}

/// Nearest-rank percentile on a pre-sorted slice. Coarser than linear
/// interpolation but deterministic and monotonic in sample size —
/// which matters more than smoothness for regression detection.
fn percentile(sorted: &[f64], q: f64) -> f64 {
    debug_assert!((0.0..=1.0).contains(&q));
    let idx = ((sorted.len() as f64 - 1.0) * q).round() as usize;
    sorted[idx]
}

// ── Quality metric helpers ───────────────────────────────────────────────────

/// Root mean square of a slice: `sqrt(mean(x²))`. Returns 0 for empty
/// input (as the other metrics do) so callers don't need per-call
/// guards.
fn rms(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = values.iter().map(|v| v * v).sum();
    (sum_sq / values.len() as f64).sqrt()
}

/// Least-squares slope of `y(t)` — the `β` coefficient in `y = α + βt`
/// — computed in the numerically stable `Σ((t - t̄)(y - ȳ))/Σ((t - t̄)²)`
/// form rather than the textbook `Σty − n·t̄·ȳ` which loses precision
/// when the centred covariance is small compared to the raw moments
/// (the common case for drift near the round-off floor).
///
/// Returns 0 for empty input, for inputs with fewer than 2 samples
/// (no line to fit), and for the degenerate case where all `t` values
/// collapse to one point (regression denominator = 0). Treating these
/// as zero is safer than panicking: a scenario that produces such a
/// sample set has bigger problems than a missing slope.
fn linear_regression_slope(t: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(t.len(), y.len(), "t and y must have matching lengths");
    if t.len() < 2 {
        return 0.0;
    }
    let n = t.len() as f64;
    let t_mean: f64 = t.iter().sum::<f64>() / n;
    let y_mean: f64 = y.iter().sum::<f64>() / n;

    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (&ti, &yi) in t.iter().zip(y.iter()) {
        let dt = ti - t_mean;
        num += dt * (yi - y_mean);
        den += dt * dt;
    }
    if den == 0.0 { 0.0 } else { num / den }
}
