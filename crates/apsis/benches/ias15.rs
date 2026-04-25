//! IAS15 full-system benchmark harness.
//!
//! ## Two modes, selected by env var
//!
//! ### Validation (default)
//!
//! ```text
//! cargo bench
//! ```
//!
//! 1. Force single-thread execution so rayon's reduction order cannot
//!    introduce run-to-run noise.
//! 2. Load `benches/baselines/ias15.toml`; abort with a clear message
//!    if missing.
//! 3. Run every scenario once, compare metrics to the baseline, and
//!    fail loud on any out-of-tolerance value *before* Criterion
//!    starts — so a correctness regression surfaces immediately
//!    instead of contaminating timing measurements.
//! 4. Hand control to Criterion for the actual wall-clock benches.
//!
//! ### Recording
//!
//! ```text
//! IAS15_BENCH_UPDATE_BASELINE=1 cargo bench
//! ```
//!
//! Runs each scenario [`common::baseline::RECORD_RUNS`] times,
//! derives per-metric tolerances, and rewrites the TOML baseline.
//! Skips Criterion entirely — recording a baseline is an authoring
//! action, not a measurement one. The rewritten file is expected to
//! be committed as part of whatever change motivated the update.
//!
//! ## Why `harness = false`
//!
//! Declared in `Cargo.toml` so the default libtest harness steps
//! aside and this `main` owns bench startup. Both the single-thread
//! enforcement and the pre-Criterion validation phase require control
//! of the entry point.

mod common;

use common::baseline::{self, BaselineFile, RecordContext, UPDATE_ENV_VAR};
use common::metrics::ScenarioMetrics;
use common::runner;
use common::scenarios::{self, ScenarioSpec};

use criterion::Criterion;
use std::collections::BTreeMap;

const MULTITHREAD_ENV_VAR: &str = "IAS15_BENCH_MULTITHREAD";

fn main() {
    // Multi-thread diagnostic mode: skip single-thread enforcement so
    // rayon uses its default worker count (num_cpus), run only the
    // multi-thread-targeted scenarios, and exit without Criterion.
    // Produces non-deterministic numerical metrics by design — the
    // goal is to answer the allocator-contention question that
    // single-thread runs cannot, not to gate on bit-exactness.
    if std::env::var(MULTITHREAD_ENV_VAR).is_ok() {
        run_multithread_diagnostic();
        return;
    }

    enforce_single_thread();

    if std::env::var(UPDATE_ENV_VAR).is_ok() {
        run_recording_mode();
        return;
    }

    run_validation_mode();
    run_criterion();
}

// ── Modes ────────────────────────────────────────────────────────────────────

/// Multi-thread diagnostic mode. Runs every scenario for which
/// `criterion_bench == false` (today: solar_n641 only — the large-N
/// investigation target), with rayon using its default multi-core
/// worker pool. Prints the phase profile for each.
///
/// Purpose: the standard harness enforces single-thread rayon for
/// bit-exact determinism, which means it cannot observe multi-thread-
/// specific costs (allocator contention, cross-thread cache effects,
/// rayon work-stealing overhead at realistic core counts). This mode
/// deliberately sacrifices determinism to expose those costs.
///
/// No baseline comparison, no Criterion timing. Output is read by a
/// human, compared against the single-thread phase profile, and the
/// delta informs whether multi-thread contention explains the gap
/// between bench numbers and real-app stutter reports.
fn run_multithread_diagnostic() {
    eprintln!("[diagnostic] multi-thread mode — rayon default workers, non-deterministic metrics");
    eprintln!(
        "[diagnostic] running {} scenario(s) flagged criterion_bench=false",
        scenarios::all().iter().filter(|s| !s.criterion_bench).count()
    );
    eprintln!();

    for spec in scenarios::all().into_iter().filter(|s| !s.criterion_bench) {
        // The phase profile is printed by `run_for_validation` itself
        // (when the ias15-profile feature is enabled). Discard the
        // returned metrics — in multi-thread mode their bit pattern
        // is non-deterministic anyway, so comparing against the
        // baseline would be noise.
        let _ = runner::run_for_validation(&spec);
    }
}

fn run_validation_mode() {
    let baseline = match baseline::load() {
        Ok(b) => b,
        Err(err) => {
            eprintln!("\nERROR: failed to load baseline");
            eprintln!("  {err}");
            eprintln!();
            eprintln!("If this is the first run, record an initial baseline with:");
            eprintln!("    {UPDATE_ENV_VAR}=1 cargo bench");
            std::process::exit(1);
        },
    };

    let scenarios_list = scenarios::all();
    let mut any_failed = false;

    for spec in &scenarios_list {
        let metrics = runner::run_for_validation(spec);
        match baseline::check_scenario(&baseline, spec.name, &metrics) {
            Ok(()) => {
                println!("[validation] {}: OK", spec.name);
            },
            Err(diff) => {
                if spec.gate_on_baseline {
                    any_failed = true;
                    eprintln!("\n[validation] {}: REGRESSION", diff.scenario);
                    for failure in &diff.failures {
                        eprintln!("  {}: {}", failure.metric, failure.reason);
                    }
                } else {
                    // Scenario is in diagnostic mode: metrics are
                    // expected to shift across runs as the investigation
                    // proceeds. Report the deltas for awareness but
                    // don't flip any_failed. Flipping `gate_on_baseline`
                    // back to `true` is the explicit action that
                    // re-arms the regression gate.
                    println!("\n[validation] {}: advisory (gate_on_baseline = false)", spec.name);
                    for failure in &diff.failures {
                        println!("  {}: {}", failure.metric, failure.reason);
                    }
                }
            },
        }
    }

    if any_failed {
        eprintln!();
        eprintln!("Baseline validation failed. If the change is intentional, update with:");
        eprintln!("    {UPDATE_ENV_VAR}=1 cargo bench");
        eprintln!(
            "and commit the updated benches/baselines/ias15.toml alongside your code change."
        );
        std::process::exit(2);
    }

    println!();
    println!(
        "[validation] all {} scenarios within tolerance — handing off to Criterion",
        scenarios_list.len()
    );
    println!();
}

fn run_recording_mode() {
    let scenarios_list = scenarios::all();
    let runs_per_scenario = baseline::RECORD_RUNS;

    println!(
        "Recording baseline: {} scenarios × {runs_per_scenario} runs each",
        scenarios_list.len()
    );
    println!();

    let mut runs: BTreeMap<String, Vec<ScenarioMetrics>> = BTreeMap::new();
    for spec in &scenarios_list {
        let samples = collect_samples(spec, runs_per_scenario);
        runs.insert(spec.name.into(), samples);
    }

    let baseline: BaselineFile = baseline::record(&runs).unwrap_or_else(|err| {
        eprintln!("\nERROR: recording failed");
        eprintln!("  {err}");
        std::process::exit(1);
    });

    let context = RecordContext::capture();
    baseline::save(&baseline, &context).unwrap_or_else(|err| {
        eprintln!("\nERROR: failed to write baseline");
        eprintln!("  {err}");
        std::process::exit(1);
    });

    println!();
    println!("Baseline written to {}", baseline::BASELINE_PATH);
    println!("Review the diff and commit alongside the change that motivated the update.");
}

fn collect_samples(spec: &ScenarioSpec, runs: usize) -> Vec<ScenarioMetrics> {
    let mut samples = Vec::with_capacity(runs);
    for i in 0..runs {
        println!("  [{}] run {}/{runs}", spec.name, i + 1);
        samples.push(runner::run_for_validation(spec));
    }
    samples
}

// ── Criterion wiring ─────────────────────────────────────────────────────────

fn run_criterion() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_full_system(&mut criterion);
    criterion.final_summary();
}

fn bench_full_system(c: &mut Criterion) {
    let mut group = c.benchmark_group("ias15_full_system");
    // Throughput is "sub-steps per unit time" — more intuitive than
    // seconds-per-batch when comparing scenarios with different
    // adaptive step sizes.
    group.throughput(criterion::Throughput::Elements(runner::STEPS_PER_ITER as u64));

    for spec in scenarios::all().into_iter().filter(|s| s.criterion_bench) {
        group.bench_function(spec.name, |b| {
            b.iter_batched_ref(
                || runner::bench_setup(&spec),
                runner::step_batch,
                criterion::BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

// ── Infrastructure ───────────────────────────────────────────────────────────

/// Pin rayon to a single thread for the entire bench process. This
/// removes non-determinism from force-evaluation reductions so the
/// metrics compared against the baseline are bit-stable across runs.
///
/// Idempotent in practice: `build_global` can only be called once per
/// process, and this is the very first thing `main` does — subsequent
/// failures would indicate someone smuggled rayon into a lazy-static
/// code path, which we want to hear about loudly.
fn enforce_single_thread() {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global().expect(
        "rayon global thread pool already initialised before bench main; \
             single-thread determinism cannot be guaranteed — aborting",
    );
}
