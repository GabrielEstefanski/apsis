//! Comparator — `recommended_dt` validation.
//!
//! Reads the CSV emitted by `recommended_dt_validation`, groups by
//! (scenario, integrator), computes peak `|ΔE/E_0|` and `|ΔLz/Lz_0|` per
//! cell, and applies the Tier 1 + Tier 2 gates declared in the protocol.
//! Writes a structured JSON report and a stdout summary; exits with code
//! 0 iff every gated cell is within tolerance.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example recommended_dt_compare -p apsis
//! cargo run --release --example recommended_dt_compare -p apsis -- --input runs.csv --output report.json
//! ```
//!
//! ## Tolerances (declared a priori in the protocol notebook)
//!
//! - VV  (gated):  `|ΔE/E_0| ≤ 1e-3`
//! - Y4  (gated):  `|ΔE/E_0| ≤ 1e-6`
//! - VV+Y4 (gated): `|ΔLz/Lz_0| ≤ 1e-10`  (or absolute `1e-10` when |Lz_0| < 1e-12)
//! - WH (informational): no a-priori bound on either metric
//!
//! ## Exit codes
//!
//! - `0` — all gated cells within tolerance.
//! - `1` — input file error.
//! - `2` — at least one gated cell exceeded tolerance.

use std::collections::BTreeMap;
use std::env;
use std::fs::{File, create_dir_all, read_to_string};
use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

// ── Tolerances from the protocol §Hypothesis ────────────────────────────── //

const TOL_REL_E_VV: f64 = 1.0e-3;
const TOL_REL_E_Y4: f64 = 1.0e-6;
const TOL_REL_LZ: f64 = 1.0e-10;
const TOL_ABS_LZ: f64 = 1.0e-10;
const LZ_REL_THRESHOLD: f64 = 1.0e-12;

// ── Records ─────────────────────────────────────────────────────────────── //

#[derive(Debug)]
struct Sample {
    e: f64,
    lz: f64,
    dt_recommended: f64,
}

#[derive(Debug)]
struct CellResult {
    scenario: String,
    integrator: String,
    n_samples: usize,
    e0: f64,
    lz0: f64,
    dt_recommended: f64,
    peak_rel_de: f64,
    peak_lz_drift: f64,
    lz_uses_absolute: bool,
    gated: bool,
    e_gate_passed: Option<bool>,
    e_gate_tolerance: Option<f64>,
    lz_gate_passed: Option<bool>,
    lz_gate_tolerance: Option<f64>,
}

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() -> ExitCode {
    let cli = parse_cli();

    let csv = match read_to_string(&cli.input_path) {
        Ok(s) => s,
        Err(err) => {
            eprintln!("ERROR: failed to read {}: {}", cli.input_path.display(), err);
            return ExitCode::from(1);
        },
    };

    let samples = match load_csv(&csv) {
        Ok(s) => s,
        Err(err) => {
            eprintln!("ERROR: {err}");
            return ExitCode::from(1);
        },
    };

    let cells = analyse(&samples);

    print_report(&cells);

    if let Some(parent) = cli.output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }
    let mut f = File::create(&cli.output_path).expect("failed to open report path");
    write_json(&mut f, &cells).expect("failed to write JSON report");
    eprintln!("\nwrote JSON report to {}", cli.output_path.display());

    let any_gate_fail = cells
        .iter()
        .any(|c| c.gated && (c.e_gate_passed == Some(false) || c.lz_gate_passed == Some(false)));
    if any_gate_fail { ExitCode::from(2) } else { ExitCode::SUCCESS }
}

// ── CSV loading ─────────────────────────────────────────────────────────── //

fn load_csv(s: &str) -> Result<BTreeMap<(String, String), Vec<Sample>>, String> {
    let mut groups: BTreeMap<(String, String), Vec<Sample>> = BTreeMap::new();
    let mut header_seen = false;
    for (lineno, line) in s.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if !header_seen {
            if !line.starts_with("scenario,integrator,sample,") {
                return Err(format!("line {}: unexpected header `{}`", lineno + 1, line));
            }
            header_seen = true;
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() != 7 {
            return Err(format!(
                "line {}: expected 7 columns, got {}: `{}`",
                lineno + 1,
                cols.len(),
                line
            ));
        }
        let scenario = cols[0].to_string();
        let integrator = cols[1].to_string();
        // Validate sample-index and t columns even though we don't store them —
        // mismatched types here indicate a malformed CSV and should fail loudly.
        let _: usize =
            cols[2].parse().map_err(|e| format!("line {}: sample index parse: {e}", lineno + 1))?;
        let _: f64 = cols[3].parse().map_err(|e| format!("line {}: t parse: {e}", lineno + 1))?;
        let e = cols[4].parse::<f64>().map_err(|e| format!("line {}: e parse: {e}", lineno + 1))?;
        let lz =
            cols[5].parse::<f64>().map_err(|e| format!("line {}: lz parse: {e}", lineno + 1))?;
        let dt_rec = cols[6]
            .parse::<f64>()
            .map_err(|e| format!("line {}: dt_recommended parse: {e}", lineno + 1))?;
        groups.entry((scenario, integrator)).or_default().push(Sample {
            e,
            lz,
            dt_recommended: dt_rec,
        });
    }
    if !header_seen {
        return Err("CSV is empty or has no header row".into());
    }
    Ok(groups)
}

// ── Analysis ────────────────────────────────────────────────────────────── //

fn analyse(samples: &BTreeMap<(String, String), Vec<Sample>>) -> Vec<CellResult> {
    let mut out = Vec::with_capacity(samples.len());
    for ((scenario, integrator), rows) in samples {
        let n = rows.len();
        if n == 0 {
            continue;
        }
        let s0 = &rows[0];
        let e0 = s0.e;
        let lz0 = s0.lz;
        let dt_rec = s0.dt_recommended;

        // Peak |ΔE/E_0|.
        let peak_rel_de = if e0.abs() > 0.0 {
            rows.iter().map(|r| ((r.e - e0) / e0).abs()).fold(0.0_f64, f64::max)
        } else {
            // Degenerate; report absolute drift for completeness.
            rows.iter().map(|r| (r.e - e0).abs()).fold(0.0_f64, f64::max)
        };

        // Peak Lz drift — relative when |Lz_0| ≥ threshold, else absolute.
        let lz_uses_absolute = lz0.abs() < LZ_REL_THRESHOLD;
        let peak_lz_drift = if lz_uses_absolute {
            rows.iter().map(|r| (r.lz - lz0).abs()).fold(0.0_f64, f64::max)
        } else {
            rows.iter().map(|r| ((r.lz - lz0) / lz0).abs()).fold(0.0_f64, f64::max)
        };

        // Apply gates per integrator.
        let (gated, e_tol, lz_tol) = match integrator.as_str() {
            "vv" => (
                true,
                Some(TOL_REL_E_VV),
                Some(if lz_uses_absolute { TOL_ABS_LZ } else { TOL_REL_LZ }),
            ),
            "y4" => (
                true,
                Some(TOL_REL_E_Y4),
                Some(if lz_uses_absolute { TOL_ABS_LZ } else { TOL_REL_LZ }),
            ),
            "wh" => (false, None, None),
            other => panic!("unknown integrator label: {other}"),
        };
        let e_gate_passed = e_tol.map(|tol| peak_rel_de <= tol);
        let lz_gate_passed = lz_tol.map(|tol| peak_lz_drift <= tol);

        out.push(CellResult {
            scenario: scenario.clone(),
            integrator: integrator.clone(),
            n_samples: n,
            e0,
            lz0,
            dt_recommended: dt_rec,
            peak_rel_de,
            peak_lz_drift,
            lz_uses_absolute,
            gated,
            e_gate_passed,
            e_gate_tolerance: e_tol,
            lz_gate_passed,
            lz_gate_tolerance: lz_tol,
        });
    }
    out
}

// ── Stdout report ───────────────────────────────────────────────────────── //

fn print_report(cells: &[CellResult]) {
    println!();
    println!("Validation — recommended_dt heuristic — comparison report");
    println!();
    println!(
        "  {:<26} {:<3} {:>13} {:>11} {:>12} {:<7}",
        "scenario", "int", "dt_rec", "|ΔE/E_0|", "|ΔLz|", "verdict"
    );
    println!("  {:-<26} {:-<3} {:->13} {:->11} {:->12} {:-<7}", "", "", "", "", "", "");
    for c in cells {
        let verdict = if !c.gated {
            "info".to_string()
        } else {
            let e_ok = c.e_gate_passed.unwrap_or(true);
            let lz_ok = c.lz_gate_passed.unwrap_or(true);
            if e_ok && lz_ok {
                "pass".to_string()
            } else {
                let mut parts = Vec::new();
                if !e_ok {
                    parts.push("E");
                }
                if !lz_ok {
                    parts.push("Lz");
                }
                format!("FAIL[{}]", parts.join(","))
            }
        };
        let lz_label = if c.lz_uses_absolute { "abs" } else { "rel" };
        println!(
            "  {:<26} {:<3} {:>13.3e} {:>11.3e} {:>9.3e}{}  {}",
            c.scenario,
            c.integrator,
            c.dt_recommended,
            c.peak_rel_de,
            c.peak_lz_drift,
            lz_label,
            verdict,
        );
    }
    println!();

    // Summary counts.
    let total_gated = cells.iter().filter(|c| c.gated).count();
    let passed_gated = cells
        .iter()
        .filter(|c| c.gated && c.e_gate_passed.unwrap_or(true) && c.lz_gate_passed.unwrap_or(true))
        .count();
    let info = cells.iter().filter(|c| !c.gated).count();
    println!("  gated cells: {passed_gated}/{total_gated} pass    informational cells: {info}");
}

// ── JSON emit (manual; no serde dependency) ────────────────────────────── //

fn write_json(f: &mut File, cells: &[CellResult]) -> std::io::Result<()> {
    writeln!(f, "{{")?;
    let any_fail = cells
        .iter()
        .any(|c| c.gated && (c.e_gate_passed == Some(false) || c.lz_gate_passed == Some(false)));
    writeln!(f, "  \"all_passed\": {},", !any_fail)?;
    writeln!(f, "  \"cells\": [")?;
    for (idx, c) in cells.iter().enumerate() {
        writeln!(f, "    {{")?;
        writeln!(f, "      \"scenario\": \"{}\",", c.scenario)?;
        writeln!(f, "      \"integrator\": \"{}\",", c.integrator)?;
        writeln!(f, "      \"n_samples\": {},", c.n_samples)?;
        writeln!(f, "      \"e0\": {:.18e},", c.e0)?;
        writeln!(f, "      \"lz0\": {:.18e},", c.lz0)?;
        writeln!(f, "      \"dt_recommended\": {:.18e},", c.dt_recommended)?;
        writeln!(f, "      \"peak_rel_de\": {:.18e},", c.peak_rel_de)?;
        writeln!(f, "      \"peak_lz_drift\": {:.18e},", c.peak_lz_drift)?;
        writeln!(f, "      \"lz_uses_absolute\": {},", c.lz_uses_absolute)?;
        writeln!(f, "      \"gated\": {},", c.gated)?;
        write_optional_bool(f, "e_gate_passed", c.e_gate_passed)?;
        write_optional_f64(f, "e_gate_tolerance", c.e_gate_tolerance)?;
        write_optional_bool(f, "lz_gate_passed", c.lz_gate_passed)?;
        write_optional_f64_last(f, "lz_gate_tolerance", c.lz_gate_tolerance)?;
        let trailing = if idx + 1 < cells.len() { "," } else { "" };
        writeln!(f, "    }}{}", trailing)?;
    }
    writeln!(f, "  ]")?;
    writeln!(f, "}}")?;
    Ok(())
}

fn write_optional_bool(f: &mut File, key: &str, v: Option<bool>) -> std::io::Result<()> {
    match v {
        Some(b) => writeln!(f, "      \"{key}\": {b},"),
        None => writeln!(f, "      \"{key}\": null,"),
    }
}

fn write_optional_f64(f: &mut File, key: &str, v: Option<f64>) -> std::io::Result<()> {
    match v {
        Some(x) => writeln!(f, "      \"{key}\": {x:.18e},"),
        None => writeln!(f, "      \"{key}\": null,"),
    }
}

fn write_optional_f64_last(f: &mut File, key: &str, v: Option<f64>) -> std::io::Result<()> {
    match v {
        Some(x) => writeln!(f, "      \"{key}\": {x:.18e}"),
        None => writeln!(f, "      \"{key}\": null"),
    }
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

struct Cli {
    input_path: PathBuf,
    output_path: PathBuf,
}

fn parse_cli() -> Cli {
    let mut input_path: Option<PathBuf> = None;
    let mut output_path: Option<PathBuf> = None;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" | "-i" => {
                input_path =
                    Some(PathBuf::from(args.next().expect("--input requires a path argument")));
            },
            "--output" | "-o" => {
                output_path =
                    Some(PathBuf::from(args.next().expect("--output requires a path argument")));
            },
            other => panic!("unknown argument: {other}"),
        }
    }
    Cli {
        input_path: input_path
            .unwrap_or_else(|| PathBuf::from("validation/recommended-dt/out/runs.csv")),
        output_path: output_path
            .unwrap_or_else(|| PathBuf::from("validation/recommended-dt/out/comparison.json")),
    }
}
