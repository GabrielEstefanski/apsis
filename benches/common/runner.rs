//! The one place that instantiates `System` and drives it through a
//! scenario. Keeping the System/force-model/integrator wiring in a
//! single module means scenarios stay data, metrics stay data, and
//! there is exactly one code path to audit when runs disagree
//! between validation and Criterion timing.
//!
//! Two entry points, serving two purposes:
//!
//! * [`run_for_validation`] — runs the full scenario, collects every
//!   accepted `consumed_dt` sample + the peak energy error, and
//!   returns a [`ScenarioMetrics`] suitable for baseline comparison
//!   or recording.
//! * [`bench_setup`] + [`step_batch`] — prepare a System at its
//!   initial state and advance it by a fixed number of sub-steps.
//!   Used as the measured closure inside `Criterion::bench_function`
//!   so each timing iteration measures the same unit of work,
//!   independent of the scenario's validation window.

use super::metrics::{self, ScenarioMetrics};
use super::scenarios::ScenarioSpec;
use gravity_sim::core::system::System;
use gravity_sim::physics::integrator::traits::IntegratorKind;

/// Tree-opening parameter for Barnes-Hut. Below the exact O(N²)
/// threshold (all scenarios in this harness are small-N) this value
/// is unused, but the constructor still requires a sensible default.
const THETA: f64 = 0.5;

/// Max Barnes-Hut tree depth; also unused for small-N.
const MAX_DEPTH: usize = 10;

/// Number of accepted sub-steps per Criterion timing iteration. Chosen
/// to amortise setup over a representative window of the controller's
/// behaviour (~few orbits for the Kepler scenarios) without pushing
/// iteration time into the regime where Criterion's warmup becomes
/// ineffective.
pub const STEPS_PER_ITER: usize = 100;

/// Run a scenario from `t=0` through `spec.duration`, collecting the
/// metrics needed for baseline comparison.
///
/// The trail ring buffer is sized to `1` because benches never render
/// trails — sizing it larger would waste work per step without
/// affecting the controller's behaviour we're measuring.
pub fn run_for_validation(spec: &ScenarioSpec) -> ScenarioMetrics {
    let mut sys = build_system(spec);
    let mut dt_samples: Vec<f64> = Vec::with_capacity(4096);
    let mut peak_energy_err = 0.0_f64;

    while sys.t() < spec.duration {
        let t_before = sys.t();
        sys.step();
        let consumed = sys.t() - t_before;
        // Zero consumed_dt would mean the controller stalled at the
        // DT_MIN floor; the IAS15 degraded_total counter catches it
        // separately. Recording a zero here is still correct — it
        // reflects the actual behaviour of the run.
        dt_samples.push(consumed);

        peak_energy_err = peak_energy_err.max(sys.metrics().rel_energy_error.abs());
    }

    let stats = sys
        .metrics()
        .adaptive_stats
        .expect("IAS15 must expose AdaptiveStats; check IntegratorKind::Ias15 was set");

    metrics::assemble(&dt_samples, &stats, peak_energy_err)
}

/// Construct a `System` in the scenario's initial state. Called once
/// per Criterion batch; Criterion excludes this cost from measurement.
pub fn bench_setup(spec: &ScenarioSpec) -> System {
    build_system(spec)
}

/// Advance the given System by [`STEPS_PER_ITER`] sub-steps. This is
/// the closure body Criterion actually times.
///
/// We deliberately do not check `sys.t() < spec.duration` here: the
/// budget-vs-duration relationship is a validation concern, not a
/// timing concern. For the timing bench we only care about the cost
/// of a representative batch of sub-steps.
pub fn step_batch(sys: &mut System) {
    for _ in 0..STEPS_PER_ITER {
        sys.step();
    }
}

// ── Internal ─────────────────────────────────────────────────────────────────

fn build_system(spec: &ScenarioSpec) -> System {
    let mut sys = System::new(
        spec.bodies.clone(),
        THETA,
        spec.dt_budget,
        MAX_DEPTH,
        /* trail_every */ 1,
    );
    sys.set_integrator(IntegratorKind::Ias15);
    sys
}
