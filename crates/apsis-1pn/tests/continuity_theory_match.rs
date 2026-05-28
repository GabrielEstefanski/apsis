//! Theory-match harness for the Continuity counter-test — sharpens the
//! `kernel_continuity_counter_test.rs` bijection gate to also pin each
//! spike magnitude against the closed-form bound from paper §3.3.

use std::sync::Arc;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::kernel::TruncatedPlummerKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

const R_CUT: f64 = 1.0;
const ALPHA: f64 = 0.8;
const DT: f64 = 1e-3;
const N_STEPS: usize = 60_000;
const SPIKE_THRESHOLD: f64 = 1e-6;

const A_ORBIT: f64 = 1.0;
const M_TOTAL: f64 = 1.0;
const DELTA_F: f64 = M_TOTAL * (1.0 - ALPHA) / (R_CUT * R_CUT); // 0.2
const EPS_REL_ABS: f64 = M_TOTAL / (2.0 * A_ORBIT); // 0.5

fn two_body_eccentric() -> Vec<Body> {
    const A: f64 = 1.0;
    const E: f64 = 0.5;
    const M_TOTAL: f64 = 1.0;
    const M_EACH: f64 = M_TOTAL / 2.0;

    let r_peri = A * (1.0 - E);
    let v_peri_rel = (M_TOTAL * (1.0 + E) / (A * (1.0 - E))).sqrt();
    let v_each = v_peri_rel / 2.0;

    let body1 = Body::rocky(M_EACH).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_each);
    let body2 = Body::rocky(M_EACH).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_each);
    vec![body1, body2]
}

fn pair_separation(bodies: &[Body]) -> f64 {
    let dx = bodies[1].pos_x - bodies[0].pos_x;
    let dy = bodies[1].pos_y - bodies[0].pos_y;
    (dx * dx + dy * dy).sqrt()
}

fn pair_relative_speed(bodies: &[Body]) -> f64 {
    let dvx = bodies[1].vel_x - bodies[0].vel_x;
    let dvy = bodies[1].vel_y - bodies[0].vel_y;
    (dvx * dvx + dvy * dvy).sqrt()
}

#[derive(Debug, Clone, Copy)]
struct Sample {
    t: f64,
    e: f64,
    r: f64,
    v_rel: f64,
}

#[derive(Debug, Clone, Copy)]
struct CrossingMeasurement {
    t: f64,
    v_cross: f64,
    spike_magnitude: f64,
    e_total: f64,
}

fn measure_crossings(samples: &[Sample]) -> Vec<CrossingMeasurement> {
    const MATCHING_WINDOW_STEPS: usize = 10;
    let mut out = Vec::new();
    for i in 1..samples.len() {
        let prev = samples[i - 1].r - R_CUT;
        let curr = samples[i].r - R_CUT;
        if prev == 0.0 || curr == 0.0 || prev.signum() == curr.signum() {
            continue;
        }
        let alpha_t = prev / (prev - curr);
        let t_cross = samples[i - 1].t + alpha_t * (samples[i].t - samples[i - 1].t);
        let v_cross = samples[i - 1].v_rel + alpha_t * (samples[i].v_rel - samples[i - 1].v_rel);

        let lo = i.saturating_sub(MATCHING_WINDOW_STEPS);
        let hi = (i + MATCHING_WINDOW_STEPS).min(samples.len() - 1);
        let mut max_rel: f64 = 0.0;
        let mut e_at_peak = samples[i].e;
        for j in (lo + 1)..=hi {
            let delta = (samples[j].e - samples[j - 1].e).abs();
            let rel = delta / samples[j].e.abs().max(1e-30);
            if rel > max_rel {
                max_rel = rel;
                e_at_peak = samples[j].e;
            }
        }

        out.push(CrossingMeasurement {
            t: t_cross,
            v_cross,
            spike_magnitude: max_rel,
            e_total: e_at_peak,
        });
    }
    out
}

/// `|ΔE_sys/E_sys| ≤ ΔF·v_cross·dt / |ε_rel|` (paper §3.3).
fn predicted_spike_bound_relative(v_cross: f64, dt: f64, delta_f: f64, eps_rel_abs: f64) -> f64 {
    delta_f * v_cross * dt / eps_rel_abs
}

fn measure_crossing_sequence() -> Vec<CrossingMeasurement> {
    let kernel = Arc::new(TruncatedPlummerKernel::new(R_CUT));
    let mut sys = System::new(two_body_eccentric(), UnitSystem::solar_canonical())
        .with_kernel(kernel)
        .with_integrator(IntegratorKind::Yoshida4)
        .with_dt(DT);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(
        UnitSystem::solar_canonical(),
    )))
    .expect("continuity-theory match: matched UnitSystem");
    sys.step();

    let mut samples: Vec<Sample> = Vec::with_capacity(N_STEPS);
    samples.push(Sample {
        t: sys.t(),
        e: sys.energy(),
        r: pair_separation(sys.bodies()),
        v_rel: pair_relative_speed(sys.bodies()),
    });
    for _ in 0..N_STEPS {
        sys.step();
        samples.push(Sample {
            t: sys.t(),
            e: sys.energy(),
            r: pair_separation(sys.bodies()),
            v_rel: pair_relative_speed(sys.bodies()),
        });
    }

    measure_crossings(&samples)
}

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn continuity_per_crossing_measurements_are_recorded() {
    let crossings = measure_crossing_sequence();
    eprintln!("[continuity-theory] {} crossings recorded", crossings.len());
    eprintln!("[continuity-theory] # | t_cross    | v_cross   | spike |dE/E|  | e_total");
    for (i, c) in crossings.iter().enumerate() {
        eprintln!(
            "[continuity-theory] {:>2} | {:>9.5} | {:>9.5} | {:>13.3e} | {:>10.5}",
            i + 1,
            c.t,
            c.v_cross,
            c.spike_magnitude,
            c.e_total,
        );
    }
    assert!(
        crossings.len() >= 4,
        "expected ≥ 4 crossings over the integration window, got {}",
        crossings.len(),
    );
    let above_threshold = crossings.iter().filter(|c| c.spike_magnitude > SPIKE_THRESHOLD).count();
    assert_eq!(
        above_threshold,
        crossings.len(),
        "{} crossings produced no detectable spike above {SPIKE_THRESHOLD:e} — \
         the bijection guarantee from §3.3 should hold here too",
        crossings.len() - above_threshold,
    );
}

/// `SAFETY_FACTOR = 1.0` — the worst-case envelope is already strict.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn spike_magnitudes_satisfy_jump_bound() {
    let crossings = measure_crossing_sequence();

    const SAFETY_FACTOR: f64 = 1.0;

    let mut worst_ratio: f64 = 0.0;
    for c in &crossings {
        let predicted = predicted_spike_bound_relative(c.v_cross, DT, DELTA_F, EPS_REL_ABS);
        let ratio = c.spike_magnitude / predicted;
        if ratio > worst_ratio {
            worst_ratio = ratio;
        }
        assert!(
            ratio < SAFETY_FACTOR,
            "spike at t = {:.5} (v_cross = {:.5}) violates bound: \
             measured |ΔE/E| = {:.3e} vs predicted {:.3e}; ratio = {:.3} >= {SAFETY_FACTOR}",
            c.t,
            c.v_cross,
            c.spike_magnitude,
            predicted,
            ratio,
        );
    }
    eprintln!(
        "[continuity-theory] {} crossings within bound; worst measured/predicted = {:.3}",
        crossings.len(),
        worst_ratio,
    );
}
