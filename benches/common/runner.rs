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

use super::metrics::{self, RunSamples, ScenarioMetrics};
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
/// Sample collection runs outside any Criterion timing iteration —
/// this function is only called from the validation and recording
/// code paths, never from `bench.iter_batched_ref`. The per-substep
/// `RunSamples.push` cost therefore does not contaminate wall-clock
/// measurements.
///
/// The trail ring buffer is sized to `1` because benches never render
/// trails — sizing it larger would waste work per step without
/// affecting the controller's behaviour we're measuring.
pub fn run_for_validation(spec: &ScenarioSpec) -> ScenarioMetrics {
    let mut sys = build_system(spec);

    // Capacity estimate: upper bound on number of accepted substeps.
    // Using `duration / dt_budget` assumes the controller never
    // exceeds the budget (true by construction — it's a cap) and is
    // a loose overestimate when the controller shrinks dt below it.
    // A loose overestimate is exactly what we want: one allocation
    // up front, zero reallocation during the validation loop.
    let capacity = expected_substeps_upper_bound(spec);
    let mut samples = RunSamples::with_capacity(capacity);

    while sys.t() < spec.duration {
        let t_before = sys.t();
        sys.step();
        let t_after = sys.t();
        // Zero consumed_dt would mean the controller stalled at the
        // DT_MIN floor; the IAS15 degraded_total counter catches it
        // separately. Recording a zero here is still correct — it
        // reflects the actual behaviour of the run.
        let consumed = t_after - t_before;
        let abs_err = sys.metrics().rel_energy_error.abs();
        samples.push(t_after, consumed, abs_err);
    }

    let stats = sys
        .metrics()
        .adaptive_stats
        .expect("IAS15 must expose AdaptiveStats; check IntegratorKind::Ias15 was set");

    metrics::assemble(&samples, &stats)
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

/// Loose upper bound on the number of accepted substeps the
/// controller will produce for `spec`. Used only for `Vec::with_capacity`,
/// so overestimating is free and underestimating forces reallocation
/// into the validation hot path — err on the side of generous.
fn expected_substeps_upper_bound(spec: &ScenarioSpec) -> usize {
    // +64 slack to cover the final partial substep, integer rounding,
    // and any transient retry spikes that push us marginally above
    // duration/dt_budget.
    let primary = (spec.duration / spec.dt_budget).ceil() as usize;
    primary.saturating_add(64)
}

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
