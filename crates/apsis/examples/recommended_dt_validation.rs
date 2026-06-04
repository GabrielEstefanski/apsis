//! Validation harness — `recommended_dt` heuristic for fixed-step integrators.
//!
//! Iterates 13 templates, computes `recommended_dt` once per scenario via a
//! VV warm-up step, then runs 100 scored substeps under each of VV, Y4, WH
//! at that dt and records per-step energy and angular momentum to CSV.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example recommended_dt_validation -p apsis
//! cargo run --release --example recommended_dt_validation -p apsis -- --output path/to/runs.csv
//! ```
//!
//! Default output path: `validation/recommended-dt/out/runs.csv` relative to
//! the workspace root.
//!
//! ## Heuristic note
//!
//! Formula, gates, and current verdict in
//! [`docs/experiments/2026-05-01-recommended-dt-heuristic.md`](../../../../docs/experiments/2026-05-01-recommended-dt-heuristic.md).
//! Constants here mirror that note; changes flow through it.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::templates::TemplateKind;

// ── Protocol constants ──────────────────────────────────────────────────── //

const N_SUBSTEPS: usize = 100;
const WARMUP_DT_FALLBACK: f64 = 1.0e-3;

/// 13 scenarios from the protocol §Methodology.
const SCENARIOS: &[(TemplateKind, &str)] = &[
    (TemplateKind::BinaryStars, "binary"),
    (TemplateKind::SolarSystem, "solar_system"),
    (TemplateKind::Kepler36, "kepler_36"),
    (TemplateKind::PlutoCharon, "pluto_charon"),
    (TemplateKind::AlphaCentauriAb, "alpha_centauri_ab"),
    (TemplateKind::HotJupiter, "hot_jupiter"),
    (TemplateKind::SunEarthMoon, "sun_earth_moon"),
    (TemplateKind::SunEarthLagrange, "sun_earth_lagrange"),
    (TemplateKind::JupiterTrojans, "jupiter_trojan"),
    (TemplateKind::Hd80606, "hd_80606_b_system"),
    (TemplateKind::Trappist1, "trappist_one"),
    (TemplateKind::ThreeBodyPythagorean, "three_body_pythagorean"),
    (TemplateKind::ThreeBodyFigureEight, "three_body_figure_eight"),
];

const INTEGRATORS: &[(IntegratorKind, &str)] = &[
    (IntegratorKind::VelocityVerlet, "vv"),
    (IntegratorKind::Yoshida4, "y4"),
    (IntegratorKind::WisdomHolman, "wh"),
];

// ── Energy and Lz from current body state ──────────────────────────────── //

/// Total mechanical energy under the exact `1/r` pair potential. Matches
/// the integrator's energy bookkeeping for runs using the default
/// `NewtonKernel::exact()` (ε = 0), which all 13 protocol scenarios do.
fn total_energy(bodies: &[Body], g: f64) -> f64 {
    let ke: f64 = bodies
        .iter()
        .map(|b| 0.5 * b.mass * (b.vel_x * b.vel_x + b.vel_y * b.vel_y + b.vel_z * b.vel_z))
        .sum();
    let mut pe = 0.0;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let dz = bodies[i].pos_z - bodies[j].pos_z;
            let r2 = dx * dx + dy * dy + dz * dz;
            pe -= g * bodies[i].mass * bodies[j].mass / r2.sqrt();
        }
    }
    ke + pe
}

fn lz(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| b.mass * (b.pos_x * b.vel_y - b.pos_y * b.vel_x)).sum()
}

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let cli = parse_cli();
    if let Some(parent) = cli.output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    let file = File::create(&cli.output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# Validation — recommended_dt heuristic for fixed-step integrators").unwrap();
    writeln!(w, "# note: docs/experiments/2026-05-01-recommended-dt-heuristic.md").unwrap();
    writeln!(w, "# scenarios={}, integrators=3, substeps={N_SUBSTEPS}", SCENARIOS.len()).unwrap();
    writeln!(w, "scenario,integrator,sample,t,e_total,lz,dt_recommended").unwrap();

    let mut total_cells = 0usize;
    let mut skipped_scenarios = 0usize;

    for (kind, scenario_name) in SCENARIOS {
        // ── Per-scenario warm-up (VV at template suggested_dt) ──────────── //
        let dt_warmup = warmup_dt(*kind);
        let mut warmup = System::from_template(*kind)
            .with_integrator(IntegratorKind::VelocityVerlet)
            .with_dt(dt_warmup);
        warmup.step();
        let dt_recommended = match warmup.recommended_dt() {
            Some(dt) => dt,
            None => {
                eprintln!(
                    "[skip] {scenario_name}: recommended_dt returned None (single body or degenerate IC)",
                );
                skipped_scenarios += 1;
                continue;
            },
        };
        let g = warmup.metrics().g_factor;

        eprintln!(
            "[scenario] {scenario_name:<26} dt_warmup={dt_warmup:.3e}  dt_recommended={dt_recommended:.3e}",
        );

        // ── Per-cell scored runs (VV, Y4, WH) ─────────────────────────────── //
        for (integrator, integrator_name) in INTEGRATORS {
            let mut sys =
                System::from_template(*kind).with_integrator(*integrator).with_dt(dt_recommended);

            // Sample 0: pre-step state.
            write_sample(&mut w, scenario_name, integrator_name, 0, &sys, g, dt_recommended);

            for n in 1..=N_SUBSTEPS {
                sys.step();
                write_sample(&mut w, scenario_name, integrator_name, n, &sys, g, dt_recommended);
            }
            total_cells += 1;
        }
    }

    w.flush().unwrap();
    eprintln!();
    eprintln!(
        "wrote {} cells ({} scenarios × {} integrators) to {}",
        total_cells,
        SCENARIOS.len() - skipped_scenarios,
        INTEGRATORS.len(),
        cli.output_path.display(),
    );
    if skipped_scenarios > 0 {
        eprintln!("({skipped_scenarios} scenarios skipped — see [skip] log)");
    }
}

// ── Helpers ─────────────────────────────────────────────────────────────── //

fn warmup_dt(kind: TemplateKind) -> f64 {
    let template = kind.build(0);
    template.suggested_dt.unwrap_or(WARMUP_DT_FALLBACK)
}

fn write_sample(
    w: &mut BufWriter<File>,
    scenario: &str,
    integrator: &str,
    sample: usize,
    sys: &System,
    g: f64,
    dt_recommended: f64,
) {
    let bodies = sys.bodies();
    let e = total_energy(bodies, g);
    let l = lz(bodies);
    writeln!(
        w,
        "{scenario},{integrator},{sample},{t:.18e},{e:.18e},{l:.18e},{dt:.18e}",
        t = sys.t(),
        dt = dt_recommended,
    )
    .unwrap();
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

struct Cli {
    output_path: PathBuf,
}

fn parse_cli() -> Cli {
    let mut output_path: Option<PathBuf> = None;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" | "-o" => {
                output_path =
                    Some(PathBuf::from(args.next().expect("--output requires a path argument")));
            },
            other => panic!("unknown argument: {other}"),
        }
    }
    Cli {
        output_path: output_path
            .unwrap_or_else(|| PathBuf::from("validation/recommended-dt/out/runs.csv")),
    }
}
