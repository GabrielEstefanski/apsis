//! Continuity counter-test data — the measured side of the §3.3 figure:
//! per-`R_c`-crossing energy-error spike `|ΔE/E|`, the closed-form bound
//! `ΔF·v_cross·δt/|E₀|`, the smooth-kernel floor, and the downsampled
//! separation `r(t)`. Constants mirror
//! `crates/apsis-1pn/tests/truncated_plummer_continuity_validation.rs`.
//!
//! ```text
//! cargo run --release --example continuity_spike_bound -p apsis-1pn -- \
//!     --output paper/figures/data/continuity_spike_bound.csv
//! ```

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
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
const TRAJECTORY_STRIDE: usize = 30;

const A_ORBIT: f64 = 1.0;
const M_TOTAL: f64 = 1.0;
const DELTA_F: f64 = M_TOTAL * (1.0 - ALPHA) / (R_CUT * R_CUT); // 0.2
const E0_ABS: f64 = M_TOTAL / (2.0 * A_ORBIT); // |specific energy| = 0.5

fn two_body_eccentric() -> Vec<Body> {
    const A: f64 = 1.0;
    const E: f64 = 0.5;
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

#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    e: f64,
    r: f64,
    v_rel: f64,
}

struct Crossing {
    t: f64,
    v_cross: f64,
    spike_rel: f64,
}

/// Per-crossing peak `|ΔE/E|` and interpolated crossing speed.
fn measure_crossings(samples: &[Sample]) -> Vec<Crossing> {
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
        for j in (lo + 1)..=hi {
            let delta = (samples[j].e - samples[j - 1].e).abs();
            max_rel = max_rel.max(delta / samples[j].e.abs().max(1e-30));
        }

        out.push(Crossing { t: t_cross, v_cross, spike_rel: max_rel });
    }
    out
}

/// Integrate the pair for `N_STEPS`. `None` kernel = smooth default with no
/// 1PN (the baseline floor run); `Some` = truncated kernel + 1PN.
fn run(kernel: Option<Arc<TruncatedPlummerKernel>>) -> Vec<Sample> {
    let mut sys = System::new(two_body_eccentric(), UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Yoshida4)
        .with_dt(DT);
    if let Some(k) = kernel {
        sys = sys.with_kernel(k);
        sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(
            UnitSystem::solar_canonical(),
        )))
        .expect("continuity figure: matched UnitSystem");
    }
    // energy() is 0 before the first step; sampling it would read as a spike.
    sys.step();

    let mut samples = Vec::with_capacity(N_STEPS + 1);
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
    samples
}

fn worst_step_rel(samples: &[Sample]) -> f64 {
    samples
        .windows(2)
        .map(|w| (w[1].e - w[0].e).abs() / w[1].e.abs().max(1e-30))
        .fold(0.0_f64, f64::max)
}

fn write_crossings(path: &PathBuf, crossings: &[Crossing], smooth_floor: f64) {
    let file = File::create(path).expect("failed to open crossings output");
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "# Continuity counter-test: energy-spike magnitudes vs analytic jump-bound (paper §3.3)"
    )
    .unwrap();
    writeln!(
        w,
        "# protocol: TruncatedPlummerKernel(R_c={R_CUT}, alpha={ALPHA}) + 1PN, Yoshida4, dt={DT:e}, {N_STEPS} steps; e=0.5 a=1 equal-mass"
    )
    .unwrap();
    writeln!(
        w,
        "# bound: |dE/E| <= delta_F * v_cross * dt / |E0|  (delta_F={DELTA_F}, E0={E0_ABS})"
    )
    .unwrap();
    writeln!(
        w,
        "# smooth_floor_rel={smooth_floor:.6e}  (default exact NewtonKernel, Yoshida4, same bodies, no 1PN; worst step-to-step |dE/E|)"
    )
    .unwrap();
    writeln!(w, "crossing,t_cross,v_cross,spike_rel,bound_rel").unwrap();
    for (i, c) in crossings.iter().enumerate() {
        let bound = DELTA_F * c.v_cross * DT / E0_ABS;
        writeln!(w, "{},{:.9e},{:.9e},{:.9e},{:.9e}", i + 1, c.t, c.v_cross, c.spike_rel, bound)
            .unwrap();
    }
    w.flush().unwrap();
}

fn write_trajectory(path: &PathBuf, samples: &[Sample]) {
    let file = File::create(path).expect("failed to open trajectory output");
    let mut w = BufWriter::new(file);
    writeln!(w, "# Pair separation r(t) for the truncated-kernel run (downsampled, stride {TRAJECTORY_STRIDE})").unwrap();
    writeln!(w, "t,r").unwrap();
    for s in samples.iter().step_by(TRAJECTORY_STRIDE) {
        writeln!(w, "{:.6e},{:.6e}", s.t, s.r).unwrap();
    }
    w.flush().unwrap();
}

fn main() {
    let crossings_path = parse_args();
    if let Some(parent) = crossings_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }
    let trajectory_path = crossings_path.with_file_name("continuity_trajectory.csv");

    let smooth_floor = worst_step_rel(&run(None));
    let samples = run(Some(Arc::new(TruncatedPlummerKernel::new(R_CUT))));
    let crossings = measure_crossings(&samples);

    write_crossings(&crossings_path, &crossings, smooth_floor);
    write_trajectory(&trajectory_path, &samples);

    let worst_ratio = crossings
        .iter()
        .map(|c| c.spike_rel / (DELTA_F * c.v_cross * DT / E0_ABS))
        .fold(0.0_f64, f64::max);
    eprintln!(
        "wrote {} crossings to {} (+ r(t) to {}); smooth floor {:.2e}; worst spike {:.0}% of bound",
        crossings.len(),
        crossings_path.display(),
        trajectory_path.display(),
        smooth_floor,
        100.0 * worst_ratio,
    );
}

fn parse_args() -> PathBuf {
    let mut output = PathBuf::from("paper/figures/data/continuity_spike_bound.csv");
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" || arg == "-o" {
            output = PathBuf::from(args.next().expect("--output requires a path argument"));
        }
    }
    output
}
