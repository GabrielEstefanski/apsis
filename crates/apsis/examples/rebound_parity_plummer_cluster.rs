//! REBOUND parity — Plummer cluster (softened kernel), apsis IAS15 side.
//!
//! Reads the committed IC CSV (single source of truth, including `# eps=` in
//! the header), integrates under IAS15 with the softened `NewtonKernel`, and
//! emits a long-format snapshot CSV plus a step-count JSON for the
//! comparator's energy-gate model.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_plummer_cluster -p apsis -- \
//!     --ics validation/rebound-parity/plummer-cluster/ics_n256.csv \
//!     --output validation/rebound-parity/plummer-cluster/out/apsis.csv \
//!     --stats-output validation/rebound-parity/plummer-cluster/out/apsis_stats.json
//! cargo run --release --example rebound_parity_plummer_cluster -p apsis -- --smoke --eps 0.0958
//! ```
//!
//! `--smoke` integrates a single softened pair for t = 1e-6 from rest and
//! prints `a_x = Δvx/t` for the convention check (`smoke_pair.py`).
//!
//! At N > 200 the adaptive-integrator scale advisory fires on stderr —
//! expected and informational for this scenario.
//!
//! ## Protocol
//!
//! `paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md`. Constants
//! here mirror the notebook in lockstep — changes are protocol changes.

use std::env;
use std::fs::{File, create_dir_all, read_to_string};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::NewtonKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

const T_FINAL: f64 = 10.0;
const SAMPLES_PER_TIME_UNIT: u64 = 10;
const DT_INITIAL: f64 = 1.0e-3;

fn main() {
    let cli = parse_cli();

    if cli.smoke {
        run_smoke(cli.eps.expect("--smoke requires --eps"));
        return;
    }

    let ics_path = cli.ics.expect("--ics is required for a cluster run");
    let (bodies, eps) = read_ics(&ics_path);
    let n = bodies.len();

    let mut sys = System::new(bodies, UnitSystem::canonical())
        .with_kernel(Arc::new(NewtonKernel::new(eps)))
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(DT_INITIAL);

    let output_path = cli.output.expect("--output is required for a cluster run");
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity -- Plummer cluster -- apsis IAS15 side").unwrap();
    writeln!(w, "# protocol: paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md")
        .unwrap();
    writeln!(
        w,
        "# n={n}, eps={eps:.18e}, t_final={T_FINAL}, samples_per_tu={SAMPLES_PER_TIME_UNIT}, dt0={DT_INITIAL:.18e}"
    )
    .unwrap();
    writeln!(w, "sample,t,body,x,y,z,vx,vy,vz").unwrap();

    let total_samples = (T_FINAL * SAMPLES_PER_TIME_UNIT as f64).round() as u64 + 1;
    let dt_sample = 1.0 / SAMPLES_PER_TIME_UNIT as f64;

    write_sample(&mut w, 0, &sys);
    for s in 1..total_samples {
        sys.integrate_until(s as f64 * dt_sample);
        write_sample(&mut w, s, &sys);
    }
    w.flush().unwrap();
    eprintln!("wrote {total_samples} samples to {}", output_path.display());

    let stats_path = cli.stats_output.expect("--stats-output is required for a cluster run");
    let mut stats = format!("{{\"substeps_total\": {}", sys.steps());
    if let Some(s) = sys.adaptive_stats() {
        stats.push_str(&format!(
            ", \"rejections\": {}, \"rejections_picard\": {}, \"rejections_truncation\": {}, \
             \"degraded\": {}, \"picard_iters\": {}, \"picard_stagnations\": {}, \
             \"shrink_grow_cycles\": {}",
            s.rejections,
            s.rejections_picard,
            s.rejections_truncation,
            s.degraded,
            s.picard_iters,
            s.picard_stagnations,
            s.shrink_grow_cycles,
        ));
    }
    stats.push('}');
    std::fs::write(&stats_path, stats).expect("failed to write stats json");
    eprintln!("apsis substeps total: {}", sys.steps());
}

fn write_sample(w: &mut BufWriter<File>, sample: u64, sys: &System) {
    let t = sys.t();
    for (i, b) in sys.bodies().iter().enumerate() {
        writeln!(
            w,
            "{sample},{t:.18e},{i},{x:.18e},{y:.18e},{z:.18e},{vx:.18e},{vy:.18e},{vz:.18e}",
            x = b.pos_x,
            y = b.pos_y,
            z = b.pos_z,
            vx = b.vel_x,
            vy = b.vel_y,
            vz = b.vel_z,
        )
        .unwrap();
    }
}

/// Single softened pair from rest, t = 1e-6: prints `a_x = Δvx/t` of the
/// probe for the convention check. Closed form: −m_src/(r²+ε²)^{3/2} at r = 1.
fn run_smoke(eps: f64) {
    let src = Body::rocky(1.0).at_3d(0.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 0.0);
    let probe = Body::rocky(1.0e-12).at_3d(1.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 0.0);
    let mut sys = System::new(vec![src, probe], UnitSystem::canonical())
        .with_kernel(Arc::new(NewtonKernel::new(eps)))
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1.0e-7);
    sys.integrate_until(1.0e-6);
    let a_x = sys.bodies()[1].vel_x / sys.t();
    println!("{a_x:.17e}");
}

fn read_ics(path: &PathBuf) -> (Vec<Body>, f64) {
    let text = read_to_string(path).expect("failed to read IC file");
    let mut eps: Option<f64> = None;
    let mut bodies = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("# eps=") {
            eps = Some(rest.trim().parse().expect("malformed eps header"));
            continue;
        }
        if line.starts_with('#') || line.starts_with("body,") || line.is_empty() {
            continue;
        }
        let f: Vec<f64> =
            line.split(',').skip(1).map(|v| v.parse().expect("malformed IC row")).collect();
        assert!(f.len() >= 7, "malformed IC row: expected 7 fields, got {}", f.len());
        bodies.push(Body::rocky(f[0]).at_3d(f[1], f[2], f[3]).with_velocity_3d(f[4], f[5], f[6]));
    }
    (bodies, eps.expect("IC file missing '# eps=' header"))
}

struct Cli {
    ics: Option<PathBuf>,
    output: Option<PathBuf>,
    stats_output: Option<PathBuf>,
    smoke: bool,
    eps: Option<f64>,
}

fn parse_cli() -> Cli {
    let mut cli = Cli { ics: None, output: None, stats_output: None, smoke: false, eps: None };
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--ics" => cli.ics = Some(PathBuf::from(args.next().expect("--ics requires a path"))),
            "--output" => {
                cli.output = Some(PathBuf::from(args.next().expect("--output requires a path")));
            },
            "--stats-output" => {
                cli.stats_output =
                    Some(PathBuf::from(args.next().expect("--stats-output requires a path")));
            },
            "--smoke" => cli.smoke = true,
            "--eps" => {
                cli.eps = Some(
                    args.next()
                        .expect("--eps requires a float")
                        .parse()
                        .expect("--eps must be a valid float"),
                );
            },
            other => panic!("unknown argument: {other}"),
        }
    }
    if let Some(e) = cli.eps {
        assert!(e >= 0.0, "--eps must be >= 0");
    }
    cli
}
