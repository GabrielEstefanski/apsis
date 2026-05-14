//! Tier 2 + Tier 4 + Tier 5 driver for the Implicit Midpoint lab notebook
//! (`docs/experiments/2026-05-14-implicit-midpoint-integrator.md`).
//!
//! Outer Solar System (Sun + Jupiter + Saturn) integrated for 10⁶ steps
//! at `dt = 0.05 yr / (2π)` ≈ Gaussian-canonical units. Reports:
//!
//! - peak `|ΔE/E₀|` and time-averaged `|⟨ΔE/E₀⟩|` for IM (Tier 2)
//! - mean iteration count + max-iter-hit fraction (Tier 4)
//! - same scenario under WHFast for cross-integrator parity (Tier 5)
//!
//! Run: `cargo run --release --example implicit_midpoint_outer_solar`.

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

fn outer_solar() -> Vec<Body> {
    // Heliocentric Cartesian, Gaussian-canonical units (G·M_sun = 1).
    // Approximate present-epoch elements; absolute geometry doesn't matter
    // for symplecticity diagnostics, only that the system stays bound.
    vec![
        Body::star(1.0).unsoftened(),
        // Jupiter: a ≈ 5.20 AU, m ≈ 9.55e-4 M_sun, v_circ ≈ √(1/5.2)
        Body::gas_giant(9.547_919e-4)
            .at(5.20, 0.0)
            .with_velocity(0.0, (1.0_f64 / 5.20).sqrt())
            .unsoftened(),
        // Saturn: a ≈ 9.58 AU, m ≈ 2.86e-4 M_sun
        Body::gas_giant(2.858_860e-4)
            .at(0.0, 9.58)
            .with_velocity(-(1.0_f64 / 9.58).sqrt(), 0.0)
            .unsoftened(),
    ]
}

fn total_energy(bodies: &[Body]) -> f64 {
    let mut ke = 0.0;
    let mut pe = 0.0;
    for (i, b) in bodies.iter().enumerate() {
        ke += 0.5 * b.mass * (b.vel_x.powi(2) + b.vel_y.powi(2) + b.vel_z.powi(2));
        for other in bodies.iter().skip(i + 1) {
            let dx = b.pos_x - other.pos_x;
            let dy = b.pos_y - other.pos_y;
            let dz = b.pos_z - other.pos_z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-30);
            pe -= b.mass * other.mass / r;
        }
    }
    ke + pe
}

fn run_one(kind: IntegratorKind, n_steps: u64, dt: f64, sample_every: u64) -> RunStats {
    let mut sys = System::new(outer_solar(), UnitSystem::solar_canonical())
        .with_exact_gravity()
        .with_integrator(kind)
        .with_dt(dt);

    let e0 = total_energy(sys.bodies());
    let mut peak_rel_de = 0.0_f64;
    let mut sum_signed_rel_de = 0.0_f64;
    let mut samples = 0_u64;

    let mut sample_q1 = sys.bodies()[1];
    let mut last_sampled = 0_u64;

    for k in 1..=n_steps {
        sys.step();
        if k.is_multiple_of(sample_every) {
            let e = total_energy(sys.bodies());
            let signed_rel = (e - e0) / e0;
            peak_rel_de = peak_rel_de.max(signed_rel.abs());
            sum_signed_rel_de += signed_rel;
            samples += 1;
            sample_q1 = sys.bodies()[1];
            last_sampled = k;
        }
    }
    let _ = (sample_q1, last_sampled); // silence unused — kept for future Tier 5 trajectory diff

    let stats = sys.metrics().adaptive_stats;
    let e_end = total_energy(sys.bodies());
    RunStats {
        peak_rel_de,
        secular_rel_de: ((e_end - e0) / e0).abs(),
        mean_signed_rel_de: sum_signed_rel_de / samples.max(1) as f64,
        cum_iterations: stats.map(|s| s.picard_iters).unwrap_or(0),
        cum_max_iter_hits: stats.map(|s| s.degraded).unwrap_or(0),
        cum_steps: stats.map(|s| s.substeps).unwrap_or(n_steps),
        final_q_planet1: sys.bodies()[1],
    }
}

#[derive(Debug)]
struct RunStats {
    peak_rel_de: f64,
    secular_rel_de: f64,
    mean_signed_rel_de: f64,
    cum_iterations: u64,
    cum_max_iter_hits: u64,
    cum_steps: u64,
    final_q_planet1: Body,
}

fn main() {
    const N_STEPS: u64 = 1_000_000;
    const DT: f64 = 0.05; // canonical time units
    const SAMPLE_EVERY: u64 = 1_000;

    println!("Implicit Midpoint — Tier 2 + Tier 4 + Tier 5");
    println!(
        "Scenario: Sun + Jupiter + Saturn, N = {N_STEPS} steps, dt = {DT}, sample every {SAMPLE_EVERY}"
    );
    println!();

    println!("Running IM…");
    let im = run_one(IntegratorKind::ImplicitMidpoint, N_STEPS, DT, SAMPLE_EVERY);
    println!("Running WHFast…");
    let wh = run_one(IntegratorKind::WHFast, N_STEPS, DT, SAMPLE_EVERY);

    println!();
    println!("=== Tier 2 — IM symplecticity ===");
    println!("  peak     |ΔE/E₀|  = {:.3e}   (oscillation amplitude)", im.peak_rel_de);
    println!("  endpoint |ΔE/E₀|  = {:.3e}   (e_end − e_start)", im.secular_rel_de);
    println!("  mean     ⟨ΔE/E₀⟩  = {:.3e}   (signed; symplectic ⇒ ≈ 0)", im.mean_signed_rel_de);
    println!();
    println!("=== Tier 4 — IM iteration diagnostic ===");
    println!("  steps                = {}", im.cum_steps);
    println!("  total iterations     = {}", im.cum_iterations);
    println!(
        "  mean iter / step     = {:.2}",
        im.cum_iterations as f64 / im.cum_steps.max(1) as f64
    );
    let n = im.cum_steps.max(1);
    if im.cum_max_iter_hits == 0 {
        println!("  max-iter exhaustions = 0 in {n} steps (≤ {:.0e})", 1.0 / n as f64);
    } else {
        println!(
            "  max-iter exhaustions = {} in {n} ({:.3e})",
            im.cum_max_iter_hits,
            im.cum_max_iter_hits as f64 / n as f64
        );
    }
    println!();
    println!("=== Tier 5 — IM vs WHFast cross-integrator ===");
    let dx = im.final_q_planet1.pos_x - wh.final_q_planet1.pos_x;
    let dy = im.final_q_planet1.pos_y - wh.final_q_planet1.pos_y;
    let dz = im.final_q_planet1.pos_z - wh.final_q_planet1.pos_z;
    let r = (im.final_q_planet1.pos_x.powi(2)
        + im.final_q_planet1.pos_y.powi(2)
        + im.final_q_planet1.pos_z.powi(2))
    .sqrt()
    .max(1.0e-30);
    let rel = (dx * dx + dy * dy + dz * dz).sqrt() / r;
    println!("  |Δr_jupiter| / r       = {rel:.3e}   (reported, no gate)");
    println!("  WHFast peak |ΔE/E₀|    = {:.3e}   (sanity)", wh.peak_rel_de);
}
