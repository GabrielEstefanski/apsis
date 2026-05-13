//! Counter-test — `Continuity::C0` precondition violation.
//!
//! Pairs the existing Mercury / softening counter-test (which exercises the
//! `Exactness::Exact` invariant) with a distinct physical signature for the
//! `Continuity::Smooth` invariant declared by the 1PN correction.
//!
//! ## The experiment
//!
//! Configure an equal-mass two-body orbit in the e = 0.5, a = 1
//! configuration — periapse r = 0.5, apoapse in a truncated Plummer
//! potential close to 1.44 — and register
//! [`PostNewtonian1PN::for_units(UnitSystem::solar_canonical())`] against a
//! [`TruncatedPlummerKernel`] with cutoff `R_c = 1`. The truncated kernel
//! provides `Exactness::Modified + Continuity::C0`; 1PN requires
//! `Exact + Smooth`. Two invariant violations are therefore expected on
//! the single `add_hamiltonian_perturbation` call.
//!
//! With the orbit bound (outside scale α = 0.5 is the maximum-range
//! case that keeps the orbit finite for these parameters) the trajectory
//! crosses the R_c shell twice per period. Symplectic integration assumes
//! a smooth Hamiltonian, so each crossing produces an impulsive
//! energy-error event.
//!
//! The test asserts:
//!
//! 1. **Two diagnostics on registration**, one per violated invariant,
//!    each tagged with `violated_invariant = "Exactness"` or `"Continuity"`.
//! 2. **Causal bijection**: every R_c crossing produces exactly one
//!    energy-error spike and every spike corresponds to an R_c crossing.
//! 3. **Magnitude separation**: spike amplitudes are several orders of
//!    magnitude larger than the smooth-kernel baseline drift.
//!
//! This is the Continuity-invariant analogue of the Mercury / softening
//! counter-test. Together the two tests demonstrate that a single match
//! mechanism catches invariant-specific failures on two formally distinct
//! axes with two distinct observable signatures.

use std::sync::{Arc, Mutex};

use apsis::core::log::{Event, Level, subscribe, unsubscribe};
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::kernel::TruncatedPlummerKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

/// Equal-mass two-body configuration at e = 0.5, a = 1 (COM at origin).
///
/// Body 1 is at `(−a(1−e)/2, 0)` with tangential velocity `+v_peri/2`.
/// Body 2 is the mirror, giving zero COM position and momentum. Both
/// bodies are explicitly unsoftened so that the only precondition
/// violations come from the kernel itself.
fn two_body_eccentric() -> Vec<Body> {
    const A: f64 = 1.0;
    const E: f64 = 0.5;
    const M_TOTAL: f64 = 1.0;
    const M_EACH: f64 = M_TOTAL / 2.0;

    let r_peri = A * (1.0 - E);
    // Relative speed at periapse for Kepler: v² = GM (1+e) / (a(1−e)).
    let v_peri_rel = (M_TOTAL * (1.0 + E) / (A * (1.0 - E))).sqrt();
    let v_each = v_peri_rel / 2.0;

    let body1 = Body::rocky(M_EACH).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_each).unsoftened();
    let body2 = Body::rocky(M_EACH).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_each).unsoftened();
    vec![body1, body2]
}

/// Current pair separation between the two bodies.
fn pair_separation(bodies: &[Body]) -> f64 {
    let dx = bodies[1].pos_x - bodies[0].pos_x;
    let dy = bodies[1].pos_y - bodies[0].pos_y;
    (dx * dx + dy * dy).sqrt()
}

/// Event captured during integration: time, total energy, pair separation.
#[derive(Debug, Clone, Copy)]
struct Sample {
    t: f64,
    e: f64,
    r: f64,
}

/// A zero-crossing where `(r[i-1] − R_c)` has a different sign from
/// `(r[i] − R_c)`. The time is linearly interpolated between the
/// bracketing samples so spikes can be matched against it precisely.
#[derive(Debug, Clone, Copy)]
struct CrossingEvent {
    t: f64,
}

fn find_crossings(samples: &[Sample], r_cut: f64) -> Vec<CrossingEvent> {
    let mut out = Vec::new();
    for i in 1..samples.len() {
        let prev = samples[i - 1].r - r_cut;
        let curr = samples[i].r - r_cut;
        if prev == 0.0 || curr == 0.0 {
            continue;
        }
        if prev.signum() != curr.signum() {
            // Linear interpolation in r to estimate the crossing time.
            let t0 = samples[i - 1].t;
            let t1 = samples[i].t;
            let alpha = prev / (prev - curr);
            let t_cross = t0 + alpha * (t1 - t0);
            out.push(CrossingEvent { t: t_cross });
        }
    }
    out
}

/// A step-level energy-error event where `|ΔE/E|` between successive
/// samples rises above `threshold`. Threshold is set well above the
/// smooth-kernel baseline drift.
#[derive(Debug, Clone, Copy)]
struct SpikeEvent {
    t: f64,
    magnitude: f64,
}

fn find_spikes(samples: &[Sample], threshold: f64) -> Vec<SpikeEvent> {
    let mut out = Vec::new();
    for i in 1..samples.len() {
        let delta = (samples[i].e - samples[i - 1].e).abs();
        let rel = delta / samples[i].e.abs().max(1e-30);
        if rel > threshold {
            out.push(SpikeEvent { t: samples[i].t, magnitude: rel });
        }
    }
    out
}

/// The registration-time invariant-violation check should fire exactly two
/// structured diagnostics: one for Exactness, one for Continuity, each
/// carrying the `violated_invariant` field that identifies which.
#[test]
fn truncated_kernel_plus_1pn_fires_both_exactness_and_continuity_warnings() {
    let captured: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = captured.clone();
    let id = subscribe(move |event: &Event| {
        // Collect every warn-level System event so we can inspect the
        // `violated_invariant` field of each.
        if event.level == Level::Warn {
            sink.lock().unwrap().push(event.clone());
        }
    });

    let kernel = Arc::new(TruncatedPlummerKernel::new(1.0));
    let mut sys = System::new(two_body_eccentric(), UnitSystem::solar_canonical())
        .with_kernel(kernel)
        .with_integrator(IntegratorKind::Yoshida4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(
        UnitSystem::solar_canonical(),
    )));

    let events = captured.lock().unwrap().clone();
    unsubscribe(id);

    // Field values come through the log bus as Debug-formatted strings
    // (warn_diag! uses `format!("{:?}", …)`), so a string literal like
    // "Exactness" is stored as the six characters `"Exactness"` including
    // the enclosing quote marks. Strip them before comparison.
    let invariants: Vec<String> = events
        .iter()
        .filter_map(|ev| {
            ev.fields.iter().find_map(|(k, v)| {
                if *k == "violated_invariant" {
                    Some(v.trim_matches('"').to_string())
                } else {
                    None
                }
            })
        })
        .collect();

    assert_eq!(
        events.len(),
        2,
        "expected exactly two invariant-violation diagnostics, got {}: {:?}",
        events.len(),
        invariants
    );
    assert!(
        invariants.iter().any(|s| s == "Exactness"),
        "missing Exactness diagnostic in {invariants:?}"
    );
    assert!(
        invariants.iter().any(|s| s == "Continuity"),
        "missing Continuity diagnostic in {invariants:?}"
    );
}

/// The central counter-test: every R_c crossing produces exactly one
/// energy-error spike, every spike corresponds to an R_c crossing, and
/// the spike magnitudes are several orders of magnitude above the
/// smooth-kernel baseline drift.
///
/// Running this with the smooth [`PlummerKernel`] on the same bodies
/// would produce no spikes — the pair of `(baseline, violated)` assertions
/// guards against false positives on the detection logic itself.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn truncated_kernel_energy_spikes_are_in_bijection_with_r_cut_crossings() {
    const R_CUT: f64 = 1.0;
    const DT: f64 = 1e-3;
    // With α = 0.8 the truncated orbit (r_peri = 0.5, r_apo ≈ 2.06) has
    // a period close to 2× the Kepler period — ~4π ≈ 12.6 simulation units
    // per revolution. 60_000 steps = 60 simulation units covers ≥ 4 full
    // orbits with 8+ R_c crossings.
    const N_STEPS: usize = 60_000;
    const SPIKE_THRESHOLD: f64 = 1e-6;

    // ── Reference run: smooth PlummerKernel at the same bodies ───────────
    //
    // Establishes the baseline drift amplitude the spike threshold must
    // clearly separate. The smooth kernel produces no impulsive events,
    // so any spike above threshold in the truncated-kernel run is
    // attributable to the discontinuity. Sampling begins AFTER the first
    // step: `System::energy()` returns 0 before any integration has run
    // (the cached kinetic/potential fields are still at their default),
    // and pairing that against the first post-step value would look like
    // a spurious spike.
    let mut sys_smooth = System::new(two_body_eccentric(), UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Yoshida4)
        .with_dt(DT);
    sys_smooth.step();

    let mut baseline_samples: Vec<Sample> = Vec::with_capacity(N_STEPS);
    baseline_samples.push(Sample {
        t: sys_smooth.t(),
        e: sys_smooth.energy(),
        r: pair_separation(sys_smooth.bodies()),
    });
    for _ in 0..N_STEPS {
        sys_smooth.step();
        baseline_samples.push(Sample {
            t: sys_smooth.t(),
            e: sys_smooth.energy(),
            r: pair_separation(sys_smooth.bodies()),
        });
    }

    let baseline_spikes = find_spikes(&baseline_samples, SPIKE_THRESHOLD);
    assert!(
        baseline_spikes.is_empty(),
        "control: smooth kernel produced {} spikes above {SPIKE_THRESHOLD:e} — \
         the threshold needs to be raised or the detection logic is wrong",
        baseline_spikes.len(),
    );

    // ── Truncated run: same bodies, TruncatedPlummerKernel ───────────────
    let kernel = Arc::new(TruncatedPlummerKernel::new(R_CUT));
    let mut sys_trunc = System::new(two_body_eccentric(), UnitSystem::solar_canonical())
        .with_kernel(kernel)
        .with_integrator(IntegratorKind::Yoshida4)
        .with_dt(DT);
    sys_trunc.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(
        UnitSystem::solar_canonical(),
    )));
    sys_trunc.step();

    let mut samples: Vec<Sample> = Vec::with_capacity(N_STEPS);
    samples.push(Sample {
        t: sys_trunc.t(),
        e: sys_trunc.energy(),
        r: pair_separation(sys_trunc.bodies()),
    });
    for _ in 0..N_STEPS {
        sys_trunc.step();
        samples.push(Sample {
            t: sys_trunc.t(),
            e: sys_trunc.energy(),
            r: pair_separation(sys_trunc.bodies()),
        });
    }

    // ── Events ──────────────────────────────────────────────────────────
    let crossings = find_crossings(&samples, R_CUT);
    let spikes = find_spikes(&samples, SPIKE_THRESHOLD);

    // Diagnostics to surface what actually happened if any assertion fails.
    let r_min_trunc = samples.iter().map(|s| s.r).fold(f64::INFINITY, f64::min);
    let r_max_trunc = samples.iter().map(|s| s.r).fold(0.0_f64, f64::max);
    let r_min_smooth = baseline_samples.iter().map(|s| s.r).fold(f64::INFINITY, f64::min);
    let r_max_smooth = baseline_samples.iter().map(|s| s.r).fold(0.0_f64, f64::max);
    let smooth_crossings = find_crossings(&baseline_samples, R_CUT);
    eprintln!(
        "[diag] smooth: r ∈ [{r_min_smooth:.4}, {r_max_smooth:.4}], crossings = {}",
        smooth_crossings.len()
    );
    eprintln!(
        "[diag] trunc:  r ∈ [{r_min_trunc:.4}, {r_max_trunc:.4}], crossings = {}, spikes = {}",
        crossings.len(),
        spikes.len()
    );
    if !spikes.is_empty() {
        let min_mag = spikes.iter().map(|s| s.magnitude).fold(f64::INFINITY, f64::min);
        let max_mag = spikes.iter().map(|s| s.magnitude).fold(0.0_f64, f64::max);
        eprintln!("[diag] spike magnitude range: [{min_mag:.3e}, {max_mag:.3e}]");
    }
    let worst_baseline_delta_dbg = baseline_samples
        .windows(2)
        .map(|w| (w[1].e - w[0].e).abs() / w[1].e.abs().max(1e-30))
        .fold(0.0_f64, f64::max);
    eprintln!("[diag] smooth worst step ΔE/E: {worst_baseline_delta_dbg:.3e}");

    assert!(
        crossings.len() >= 4,
        "expected at least 4 R_c crossings over the integration window, got {}",
        crossings.len(),
    );

    // ── Causal bijection: |spikes| == |crossings|, each paired ──────────
    assert_eq!(
        spikes.len(),
        crossings.len(),
        "spike count ({}) does not match crossing count ({}); \
         spikes without crossings break the causal link",
        spikes.len(),
        crossings.len(),
    );

    // Each spike must be temporally adjacent (within 10·DT) to the
    // matching crossing — the discontinuity produces its numerical
    // signature within a single step of the actual crossing event.
    const MATCHING_TOLERANCE: f64 = 10.0 * DT;
    for (spike, crossing) in spikes.iter().zip(crossings.iter()) {
        let dt = (spike.t - crossing.t).abs();
        assert!(
            dt < MATCHING_TOLERANCE,
            "spike at t = {:.5} lagged crossing at t = {:.5} by {:.5} — \
             expected < {:.5}; non-causal spike indicates the signature is \
             not uniquely tied to the discontinuity",
            spike.t,
            crossing.t,
            dt,
            MATCHING_TOLERANCE,
        );
    }

    // ── Magnitude separation from smooth-kernel baseline ────────────────
    //
    // The spike amplitudes must be several orders of magnitude above the
    // worst step-to-step baseline change in the smooth-kernel run. This
    // rules out generic numerical noise as the source.
    let worst_baseline_delta = baseline_samples
        .windows(2)
        .map(|w| (w[1].e - w[0].e).abs() / w[1].e.abs().max(1e-30))
        .fold(0.0_f64, f64::max);

    let min_spike_magnitude = spikes.iter().map(|s| s.magnitude).fold(f64::INFINITY, f64::min);

    assert!(
        min_spike_magnitude > worst_baseline_delta * 1000.0,
        "smallest spike ({min_spike_magnitude:e}) not separated from smooth-kernel \
         worst step ({worst_baseline_delta:e}) by at least 3 orders of magnitude; \
         the discontinuity signature is indistinguishable from baseline noise",
    );
}
