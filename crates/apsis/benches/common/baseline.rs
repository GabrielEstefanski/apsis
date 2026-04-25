//! Versioned numerical-regression gate for the IAS15 harness.
//!
//! The file at [`BASELINE_PATH`] is the source of truth: every
//! scenario's expected metric values and per-metric tolerances live
//! there, are committed to git, and evolve only through explicit
//! update runs (`IAS15_BENCH_UPDATE_BASELINE=1 cargo bench`). A bench
//! run without the env var will parse this file, run each scenario,
//! and fail with a detailed diff the moment a metric drifts outside
//! its allowed window.
//!
//! # Tolerance semantics
//!
//! A [`ToleranceSpec`] carries the baseline `value` plus exactly one
//! of:
//!
//! * `tol_abs` — `|measured − value| ≤ tol_abs`. Suited to counters
//!   that are either identical or meaningfully different; `tol_abs=0`
//!   demands bit-exact reproduction.
//! * `tol_factor` — `value/factor ≤ measured ≤ value·factor`. Suited
//!   to floats where legitimate ULP-level jitter may appear; the
//!   recording pass sizes this from observed run-to-run variation.
//!
//! A hard cap [`MAX_ALLOWED_TOL_FACTOR`] prevents silent drift via
//! successive relaxations: anything beyond 1.5× must come from an
//! explicit baseline update whose diff is reviewable in a PR.
//!
//! # Recording pass
//!
//! When invoked with `IAS15_BENCH_UPDATE_BASELINE=1`, the harness
//! runs each scenario [`RECORD_RUNS`] times, then for every metric:
//!
//! * `Counter` — require `min == max` across runs. Any jitter is a
//!   determinism bug (rayon threading, use of `HashMap` ordering,
//!   etc.) and is surfaced as a recording failure rather than
//!   absorbed into a widened tolerance.
//! * `Float` — size `tol_factor = 1 + 2·(max−min)/mean`, clamped to
//!   `[1.0, MAX_ALLOWED_TOL_FACTOR]`. This catches regressions at
//!   ~½× the observed jitter while absorbing genuine noise.

use super::metrics::{MetricTier, ScenarioMetrics};
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Path to the baseline file, resolved from the workspace root.
pub const BASELINE_PATH: &str = "benches/baselines/ias15.toml";

/// Number of runs per scenario during a baseline recording pass.
/// Larger values characterise the jitter distribution more precisely
/// at the cost of proportional wall time; 10 is a practical balance
/// for scenarios in the 100ms–few-seconds range.
pub const RECORD_RUNS: usize = 10;

/// Hard upper bound on any `tol_factor` entry. Values above this are
/// rejected at parse time: if a metric legitimately flutters beyond
/// 1.5× its baseline, something has changed about the underlying
/// computation and deserves an investigation, not a looser tolerance.
pub const MAX_ALLOWED_TOL_FACTOR: f64 = 1.5;

/// Name of the env var that flips the harness from validation mode
/// to recording mode. Using an env var (rather than a CLI flag)
/// keeps the entry point transparent to Criterion's own arg parser.
pub const UPDATE_ENV_VAR: &str = "IAS15_BENCH_UPDATE_BASELINE";

// ── Types ────────────────────────────────────────────────────────────────────

/// Tolerance entry for one metric. Exactly one of `tol_abs` and
/// `tol_factor` must be set; [`Self::validate`] enforces this and
/// also checks the `MAX_ALLOWED_TOL_FACTOR` ceiling.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct ToleranceSpec {
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tol_abs: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tol_factor: Option<f64>,
}

impl ToleranceSpec {
    fn validate(&self) -> Result<(), String> {
        match (self.tol_abs, self.tol_factor) {
            (Some(_), Some(_)) => Err("both tol_abs and tol_factor set; pick one".into()),
            (None, None) => Err("neither tol_abs nor tol_factor set; pick one".into()),
            (Some(a), _) if a < 0.0 || !a.is_finite() => {
                Err(format!("tol_abs {a} must be finite and non-negative"))
            },
            (_, Some(f)) if f < 1.0 || !f.is_finite() => {
                Err(format!("tol_factor {f} must be finite and ≥ 1.0"))
            },
            (_, Some(f)) if f > MAX_ALLOWED_TOL_FACTOR => Err(format!(
                "tol_factor {f} exceeds cap {MAX_ALLOWED_TOL_FACTOR}; \
                 investigate the underlying fluctuation before relaxing further"
            )),
            _ => Ok(()),
        }
    }

    /// Check `measured` against this spec. On violation returns a
    /// human-readable reason that makes the baseline-vs-measured
    /// delta obvious at a glance.
    ///
    /// `tol_factor == 1.0` is interpreted as "bit-exact required" and
    /// uses `to_bits` comparison. The alternative — `value/1 ≤ x ≤
    /// value·1` — is mathematically equivalent to bit-exact *when*
    /// the stored value and the measured value share a bit pattern,
    /// but can fail spuriously if stringification through TOML loses
    /// a round-trip bit that no reviewer would call a regression.
    /// When tol_factor > 1.0, the range check is the well-posed
    /// semantic and we use it unchanged.
    fn check(&self, measured: f64) -> Result<(), String> {
        if let Some(abs) = self.tol_abs {
            let delta = (measured - self.value).abs();
            if delta > abs {
                return Err(format!(
                    "Δ={delta:e} > tol_abs={abs:e} \
                     (baseline={}, measured={})",
                    self.value, measured
                ));
            }
        }
        if let Some(factor) = self.tol_factor {
            if factor == 1.0 {
                if measured.to_bits() != self.value.to_bits() {
                    return Err(format!(
                        "bit-exact mismatch (tol_factor=1.0): \
                         baseline={} (0x{:016x}), measured={} (0x{:016x})",
                        self.value,
                        self.value.to_bits(),
                        measured,
                        measured.to_bits(),
                    ));
                }
            } else {
                let lo = self.value / factor;
                let hi = self.value * factor;
                if measured < lo || measured > hi {
                    return Err(format!(
                        "measured={measured:e} outside [{lo:e}, {hi:e}] \
                         (baseline={}, tol_factor={factor})",
                        self.value
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Per-scenario map of metric name → tolerance spec. `BTreeMap` keeps
/// keys sorted so the serialised TOML has a stable diff-friendly order.
pub type ScenarioBaseline = BTreeMap<String, ToleranceSpec>;

/// Top-level baseline: scenario name → metric-level baseline. Also a
/// `BTreeMap` so scenarios appear in alphabetical order in the file.
#[derive(Debug, Clone, Default)]
pub struct BaselineFile {
    pub scenarios: BTreeMap<String, ScenarioBaseline>,
}

// ── I/O ──────────────────────────────────────────────────────────────────────

/// Load the baseline from [`BASELINE_PATH`]. Errors propagate file
/// I/O, parse, and validation failures with enough context to point
/// at the offending key.
pub fn load() -> Result<BaselineFile, String> {
    let path = baseline_path();
    let content = std::fs::read_to_string(&path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;

    let value: toml::Value =
        toml::from_str(&content).map_err(|e| format!("failed to parse {}: {e}", path.display()))?;

    let root_table =
        value.as_table().ok_or_else(|| format!("{}: root must be a table", path.display()))?;

    let mut scenarios: BTreeMap<String, ScenarioBaseline> = BTreeMap::new();
    for (scenario_name, scenario_val) in root_table {
        let metrics_table = scenario_val.as_table().ok_or_else(|| {
            format!("{}: [{scenario_name}] must be a table of metric entries", path.display())
        })?;

        let mut baseline: ScenarioBaseline = BTreeMap::new();
        for (metric_name, metric_val) in metrics_table {
            let spec: ToleranceSpec = metric_val.clone().try_into().map_err(|e| {
                format!("{}: [{scenario_name}].{metric_name}: parse error — {e}", path.display())
            })?;
            spec.validate()
                .map_err(|e| format!("{}: [{scenario_name}].{metric_name}: {e}", path.display()))?;
            baseline.insert(metric_name.clone(), spec);
        }

        scenarios.insert(scenario_name.clone(), baseline);
    }

    Ok(BaselineFile { scenarios })
}

/// Write the baseline file to [`BASELINE_PATH`], overwriting any
/// existing content. A header comment records when and on which
/// commit the file was regenerated so a reviewer reading a diff
/// immediately knows whether the context of the change matches the
/// PR's claim.
pub fn save(baseline: &BaselineFile, context: &RecordContext) -> Result<(), String> {
    let path = baseline_path();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }

    let mut out = String::new();
    out.push_str("# IAS15 benchmark baseline — versioned regression gate.\n");
    out.push_str("#\n");
    out.push_str("# Regenerate with:  IAS15_BENCH_UPDATE_BASELINE=1 cargo bench\n");
    out.push_str("# Validate with:    cargo bench\n");
    out.push_str("#\n");
    out.push_str(&format!("# Recorded at: {}\n", context.timestamp));
    out.push_str(&format!("# Git commit:  {}\n", context.commit));
    out.push_str(&format!("# Runs per scenario: {}\n", context.runs_per_scenario));
    out.push('\n');

    for (scenario_name, scenario) in &baseline.scenarios {
        out.push_str(&format!("[{scenario_name}]\n"));
        for (metric_name, spec) in scenario {
            out.push_str(&format_tolerance_line(metric_name, spec));
        }
        out.push('\n');
    }

    std::fs::write(&path, out).map_err(|e| format!("failed to write {}: {e}", path.display()))?;
    Ok(())
}

// ── Validation ───────────────────────────────────────────────────────────────

/// Detailed report of a single scenario's failed metrics. Kept as a
/// separate type (rather than `Vec<String>`) so the caller can
/// distinguish "scenario missing from baseline" from "some metrics
/// failed" without string-parsing.
#[derive(Debug)]
pub struct ScenarioDiff {
    pub scenario: String,
    pub failures: Vec<MetricFailure>,
}

#[derive(Debug)]
pub struct MetricFailure {
    pub metric: String,
    pub reason: String,
}

/// Check the measured metrics for `scenario` against the loaded
/// baseline. Missing baseline entries for known metrics and extra
/// metrics (typo guard) both surface as failures.
pub fn check_scenario(
    baseline: &BaselineFile,
    scenario: &str,
    measured: &ScenarioMetrics,
) -> Result<(), ScenarioDiff> {
    let scenario_baseline = match baseline.scenarios.get(scenario) {
        Some(b) => b,
        None => {
            return Err(ScenarioDiff {
                scenario: scenario.into(),
                failures: vec![MetricFailure {
                    metric: "<scenario>".into(),
                    reason: format!(
                        "no baseline entry for scenario '{scenario}'; \
                         run {UPDATE_ENV_VAR}=1 cargo bench to record"
                    ),
                }],
            });
        },
    };

    let mut failures = Vec::new();
    for (metric_name, _tier) in ScenarioMetrics::ALL {
        let spec = match scenario_baseline.get(*metric_name) {
            Some(s) => s,
            None => {
                failures.push(MetricFailure {
                    metric: (*metric_name).into(),
                    reason: format!(
                        "no baseline entry for '{metric_name}' under [{scenario}]; \
                         update baseline to record"
                    ),
                });
                continue;
            },
        };
        let measured_value = measured
            .get(metric_name)
            .expect("ScenarioMetrics::ALL and ScenarioMetrics::get must agree");
        if let Err(reason) = spec.check(measured_value) {
            failures.push(MetricFailure { metric: (*metric_name).into(), reason });
        }
    }

    if failures.is_empty() {
        Ok(())
    } else {
        Err(ScenarioDiff { scenario: scenario.into(), failures })
    }
}

// ── Recording ────────────────────────────────────────────────────────────────

/// Metadata captured at recording time, embedded as comments in the
/// output file to give reviewers context.
pub struct RecordContext {
    pub timestamp: String,
    pub commit: String,
    pub runs_per_scenario: usize,
}

impl RecordContext {
    pub fn capture() -> Self {
        let timestamp = chrono_like_utc();
        let commit = git_head_hash().unwrap_or_else(|| "<unknown>".into());
        Self { timestamp, commit, runs_per_scenario: RECORD_RUNS }
    }
}

/// Build a baseline from a batch of per-scenario metric samples.
/// The samples vector for each scenario must be non-empty; identical
/// scenarios across multiple runs are expected to produce identical
/// counter values (enforced here) and near-identical float values
/// (absorbed into `tol_factor`).
pub fn record(runs: &BTreeMap<String, Vec<ScenarioMetrics>>) -> Result<BaselineFile, String> {
    let mut file = BaselineFile::default();

    for (scenario_name, samples) in runs {
        if samples.is_empty() {
            return Err(format!(
                "scenario '{scenario_name}': no samples recorded \
                 — cannot derive baseline"
            ));
        }

        let mut scenario: ScenarioBaseline = BTreeMap::new();
        for (metric_name, tier) in ScenarioMetrics::ALL {
            let values: Vec<f64> =
                samples.iter().map(|m| m.get(metric_name).expect("metric exists")).collect();
            let spec = derive_tolerance(metric_name, *tier, &values)?;
            scenario.insert((*metric_name).into(), spec);
        }
        file.scenarios.insert(scenario_name.clone(), scenario);
    }

    Ok(file)
}

fn derive_tolerance(
    metric_name: &str,
    tier: MetricTier,
    values: &[f64],
) -> Result<ToleranceSpec, String> {
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = values.iter().sum::<f64>() / values.len() as f64;

    // When all runs are bit-identical we must store that shared value
    // verbatim, not an arithmetic mean of it. `(X * N) / N` drifts by
    // up to 1 ULP via intermediate rounding, which would bite every
    // subsequent validation run (which computes `X` directly and then
    // fails the bit-exact check).
    let value = if min == max { min } else { mean };

    match tier {
        MetricTier::Counter => {
            // Hard determinism requirement. Any run-to-run variation
            // at this tier is a bug in the harness (threading leak,
            // non-deterministic container iteration, etc.), not
            // something to absorb into a tolerance.
            if min != max {
                return Err(format!(
                    "counter metric '{metric_name}' jittered across runs: \
                     min={min}, max={max}, values={values:?} — \
                     harness determinism invariant violated, investigate \
                     before rerunning"
                ));
            }
            Ok(ToleranceSpec { value, tol_abs: Some(0.0), tol_factor: None })
        },
        MetricTier::Float => {
            // Bit-identical across runs → tol_factor = 1.0 (still
            // asserts exact reproduction; ULP drift in future runs
            // will fire). Otherwise scale tol_factor to 2× observed
            // jitter, clamped at the hard cap.
            let tol_factor = if min == max || mean == 0.0 {
                1.0
            } else {
                let rel_range = (max - min).abs() / mean.abs();
                (1.0 + 2.0 * rel_range).min(MAX_ALLOWED_TOL_FACTOR)
            };
            Ok(ToleranceSpec { value, tol_abs: None, tol_factor: Some(tol_factor) })
        },
    }
}

// ── Formatting / utility ─────────────────────────────────────────────────────

fn format_tolerance_line(metric_name: &str, spec: &ToleranceSpec) -> String {
    // Inline-table format so a single metric fits on one line, keeping
    // git diffs tight and readable.
    if let Some(abs) = spec.tol_abs {
        format!(
            "{metric_name} = {{ value = {}, tol_abs = {} }}\n",
            format_f64(spec.value),
            format_f64(abs),
        )
    } else if let Some(factor) = spec.tol_factor {
        format!(
            "{metric_name} = {{ value = {}, tol_factor = {} }}\n",
            format_f64(spec.value),
            format_f64(factor),
        )
    } else {
        unreachable!("ToleranceSpec::validate rejects the both-none case before save")
    }
}

/// Format a float for TOML such that parse-then-format round-trips to
/// the identical f64 bit pattern.
///
/// Rust's default `Display` for f64 emits the shortest decimal that
/// round-trips; `{:e}` does the same in scientific form. We pick
/// between the two based on magnitude (readability only — precision
/// is identical) and add a trailing `.0` for integer-valued floats
/// so TOML parses them as floats, not ints.
///
/// An earlier revision used `{:.6e}` for magnitude extremes, which
/// truncated to 6 decimal digits and broke round-tripping: baselines
/// recorded from a bench run no longer matched the metrics produced
/// by a replay of the same run. Never reintroduce a precision modifier
/// here — the default round-trip guarantee is the whole point.
fn format_f64(x: f64) -> String {
    if x == 0.0 {
        return "0.0".into();
    }
    if !x.is_finite() {
        return format!("{x}");
    }
    let mag = x.abs();
    if !(1e-3..1e6).contains(&mag) {
        format!("{x:e}")
    } else if x == x.trunc() {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}

fn baseline_path() -> PathBuf {
    // CARGO_MANIFEST_DIR is the package root at build time and
    // remains valid at runtime (Criterion benches run with the
    // workspace root as cwd, but relying on an absolute path derived
    // from the manifest is more robust to invocation from sub-dirs).
    let root = std::env::var("CARGO_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("."));
    root.join(BASELINE_PATH)
}

/// UTC timestamp in RFC-3339-ish format without bringing in `chrono`
/// as a dep just for a header comment. Second resolution is enough.
fn chrono_like_utc() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    // Minimal YYYY-MM-DDTHH:MM:SSZ from epoch seconds. Uses standard
    // civil-from-epoch arithmetic (Howard Hinnant).
    let (y, mo, d, h, mi, se) = epoch_to_ymdhms(secs as i64);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{se:02}Z")
}

fn epoch_to_ymdhms(secs: i64) -> (i64, u32, u32, u32, u32, u32) {
    let days = secs.div_euclid(86_400);
    let tod = secs.rem_euclid(86_400);
    let h = (tod / 3600) as u32;
    let mi = ((tod % 3600) / 60) as u32;
    let se = (tod % 60) as u32;
    // Howard Hinnant's civil_from_days:
    let z = days + 719_468;
    let era = if z >= 0 { z / 146_097 } else { (z - 146_096) / 146_097 };
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d, h, mi, se)
}

fn git_head_hash() -> Option<String> {
    let output = std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .current_dir(
            std::env::var("CARGO_MANIFEST_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|_| PathBuf::from(".")),
        )
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let hash = String::from_utf8(output.stdout).ok()?;
    Some(hash.trim().to_string())
}
