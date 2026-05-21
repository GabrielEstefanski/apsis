# Cross-platform bitwise reproducibility — IAS15 adaptive controller via `libm::pow`

**Date:** 2026-05-20
**Subject:** Demonstrate that an apsis record (Cargo.lock + rustc version + source SHA) reproduces f64-bit-identical trajectories across heterogeneous x86_64 hosts (Windows AMD Zen 4 vs Linux Intel Ice Lake) for the four REBOUND-parity scenarios and the Mercury 1PN federation gate, once the only libc transcendental in the IAS15 adaptive controller is routed through the deterministic `libm` crate.

**Baseline commit (pre-fix):** `06bd0a9` (master, post-#152)

**Fix commit:** `ce0f2a9` (`experiment/cross-platform-libm`)

**Tooling:** apsis IAS15 (`crates/apsis/src/physics/integrator/ias15.rs`), `libm = "0.2"` (pure-Rust math), `validation/cross-platform/run_linux_side.sh`, `validation/cross-platform/compare.py`

**Status:** Single bidirectional run executed 2026-05-20. Diagnostic phase identified `f64::powf(1/7)` in the IAS15 step-size controller as the sole bifurcation source. Post-fix run on identical hosts yields byte-identical output files (independent SHA256 verification).

---

## Abstract

The v0.1 paper's central claim — that an apsis record's TOML provenance plus the workspace's `Cargo.lock` plus the rustc version is a complete recipe for reproducing the trajectory it describes — was, at the start of this experiment, defensible only at the physics-equivalent level (conserved-invariant agreement across platforms) and not at the bitwise level. A diagnostic cross-platform run on a c6i.large EC2 instance (Intel Xeon 8375C, Ice Lake) against an identically-configured AMD Ryzen 5 7600X Windows host showed energy preserved to the f64 floor in every parity scenario but trajectories diverging by ~1e+15 ULPs over 100 Kepler orbits — the signature of phase drift with conserved energy.

The signature pointed at the adaptive step-size controller, not the force model. A single line of the IAS15 step-size optimizer (`crates/apsis/src/physics/integrator/ias15.rs:1981`) called `f64::powf(1.0/7.0)` on the dimensionless ratio in the IAS15 7th-root step-size formula (Rein & Spiegel 2015 §2.3 eq. 11). `f64::powf` routes to the platform's libc `pow`; the last-ULP outputs of glibc `pow` and Microsoft UCRT `pow` are not bitwise-equivalent. Replacing this single call with `libm::pow` (the `libm` crate's pure-Rust implementation, deterministic across x86_64 targets) restored full bitwise cross-platform reproducibility for all four parity scenarios (Kepler, figure-8, Pythagorean, retrograde Kepler) and the Mercury 1PN perihelion gate, validated independently by SHA256 of the output files.

The claim that paper artifacts now carry is conditional and scoped: "Bitwise cross-platform reproducibility is achievable when all libc-bound transcendentals in integration-critical paths are routed through deterministic alternatives. This experiment demonstrates the principle for IAS15 + direct summation + Newton kernel + 1PN operator on x86_64; extending it to WHFast (Kepler solver), Mercurius (Hill-radius switching), and the central-force operator requires the analogous audit-and-replace pass, tracked as issues #159–#161."

---

## Motivation

The apsis-record format (PR #92) stores, alongside frame data, a TOML header containing the integrator kind, every registered operator and its declared kernel requirements, the workspace `Cargo.lock` BLAKE3 hash, the `rustc` version, and a per-system reproducibility seed. The stated contract is: `{ source SHA, Cargo.lock, rustc version, header }` is sufficient to regenerate the byte-identical frame stream that the record's trailer hashes.

That contract is trivially provable intra-platform — the same machine running the same commit twice will hash to the same trailer. The CI release-validation job exercises exactly this. The unmet question, until this experiment, is whether the contract holds *across* platforms.

A reviewer familiar with REBOUND, GADGET, AREPO, or PKDGRAV is aware that bitwise cross-platform reproducibility in scientific N-body codes is rare and typically not promised. The principled reasons are well-known: IEEE 754 does not bound transcendental output to the last ULP, libc/libm implementations differ between glibc and MSVC, LLVM does not promise identical instruction selection across targets, FMA changes rounding, and any chaotic regime amplifies microscopic differences exponentially. The expected outcome before this experiment, therefore, was *physics-equivalent reproducibility* (conserved invariants agree at the f64 floor across platforms; trajectories phase-drift) — a strictly weaker, but still publishable, claim.

The experiment reported here found that the weaker outcome was an artifact of a single tractable call site. The stronger claim is recoverable, conditional on a stated methodology and a documented audit of remaining call sites.

### What this experiment is NOT testing

- Not a cross-architecture claim. Targets exercised are both x86_64 with AVX2 + AVX-512 available. ARM, RISC-V, POWER are unaddressed.
- Not a claim about `target-cpu = native` codegen. Both binaries were built with default codegen (`target-cpu = x86-64` baseline). Microarch-specific tuning is unaddressed.
- Not a Mercurius, WHFast, BH, or apsis-central reproducibility claim. The integration-critical paths of those code paths contain additional libc transcendentals (sin/cos/cosh, cbrt, powf respectively) and are tracked separately as issues #159, #160, #161.
- Not a claim that all apsis simulations reproduce bitwise on any pair of x86_64 hosts. The claim is scoped to the parity portfolio above with the libm fix in place.

---

## Setup

### Host A — Windows reference

| Field | Value |
| --- | --- |
| CPU | AMD Ryzen 5 7600X (Zen 4, 6c/12t, 4.7 GHz, AVX2 + AVX-512) |
| OS | Microsoft Windows NT 10.0.26200 |
| rustc | 1.94.1 (e408947bf 2026-03-25) |
| cargo | 1.94.1 |
| target | `x86_64-pc-windows-msvc` (UCRT libc) |
| `Cargo.lock` SHA256 | `F39E4946109916E8C8CCBE9D482502CCFB4166391BC8BCFFCF8302743A6EAFE1` |
| git SHA (pre-fix) | `06bd0a974e6873c70122d613cb22c4e9d3e2e4be` |

### Host B — Linux EC2 instance

| Field | Value |
| --- | --- |
| Instance | c6i.large (Spot, 2 vCPU, 4 GiB RAM) |
| CPU | Intel Xeon Platinum 8375C (Ice Lake, AVX2 + AVX-512) |
| OS | Ubuntu 24.04.4 LTS, kernel 6.17.0-1012-aws |
| glibc | 2.39 |
| rustc | 1.94.1 (pinned via `rustup default 1.94.1` to match Host A) |
| target | `x86_64-unknown-linux-gnu` (glibc) |
| `Cargo.lock` SHA256 | `F39E4946109916E8C8CCBE9D482502CCFB4166391BC8BCFFCF8302743A6EAFE1` (matches Host A) |
| git SHA (post-fix) | `e06ba47` (`experiment/cross-platform-libm`, after `ce0f2a9`) |

Both hosts compile against the same workspace `Cargo.lock` and the same `rustc 1.94.1`, isolating microarch + libc + LLVM target as the experimental axes.

---

## Protocol

### Scenarios

The portfolio is the four published REBOUND-parity scenarios plus the Mercury 1PN gate:

| Scenario | Source | Output |
| --- | --- | --- |
| Kepler e=0.5, 100 orbits | `crates/apsis/examples/rebound_parity_kepler.rs` | 101-row CSV (orbit, t, x/y/v positions, e_total) |
| Figure-8, 50 periods | `crates/apsis/examples/rebound_parity_figure8.rs` | 2001-row CSV (3-body, 12 dynamical + e_total) |
| Pythagorean, 70 t.u. | `crates/apsis/examples/rebound_parity_pythagorean.rs` | 2101-row CSV (chaotic, 3-body) |
| Retrograde Kepler, 10⁴ orbits | `crates/apsis/examples/rebound_parity_retrograde.rs` | 10001-row CSV (long-horizon 2-body, $L_z < 0$) |
| Mercury 1PN perihelion, 500 orbits | `crates/apsis-1pn/examples/mercury_perihelion.rs` | stdout text, `rate = X arcsec/century` |

All five run under `IAS15` with `NewtonKernel::exact()` (ε = 0) and direct O(N²) summation. The Mercury gate additionally registers `PostNewtonian1PN::solar_units()` as a Hamiltonian perturbation.

### Runner

The Linux side is captured by `validation/cross-platform/run_linux_side.sh`. The script captures host metadata (rustc, OS, kernel, CPU flags, glibc, `Cargo.lock` SHA256, git HEAD), builds the workspace in release, runs each scenario writing outputs into `/tmp/apsis-xplat/linux/`, and tarballs the result for transport.

The Windows side runs the same five commands directly via PowerShell, writing into `validation/cross-platform/windows/`.

### Comparator

`validation/cross-platform/compare.py` reads paired outputs from two directories (`--a`, `--b`) and reports per-column ULP-distance (in units of `2⁻⁵² · max(|a|, |b|)`) plus the Mercury rate delta in ppm of the 43 arcsec/century GR prediction. Bit-equality is reported as "bit-equal / total" per column; SHA256 of the file is the independent cross-check.

---

## Phase 1 — pre-fix cross-platform run

Both hosts executed against commit `06bd0a9`, with `rustc 1.94.1`, `Cargo.lock` SHA256 matching. The result was the well-formed signature of *phase drift with conserved energy*:

| Scenario | max ULP drift (position columns) | max ULP drift (e_total) | max \|Δ\| absolute (position) |
| --- | --- | --- | --- |
| kepler | ≈ 4 × 10¹⁵ | 11 | 2.7 × 10⁻⁸ |
| figure8 | ≈ 9 × 10¹⁵ | 5 | 4.1 × 10⁻² |
| pythagorean | ≈ 9 × 10¹⁵ | 1.9 × 10⁵ (chaos amp) | 2.4 × 10⁰ |
| retrograde | ≈ 4 × 10¹⁵ | 120 | 6.8 × 10⁻² |
| mercury_perihelion | n/a (single scalar) | n/a | 42.993 vs 42.991 arcsec/century (Δ = 46 ppm of GR) |

The Kepler scenario alone disambiguates the diagnosis: Kepler at e = 0.5 is integrable, so the trajectory cannot diverge from chaos amplification. Linear last-ULP round-off accumulation would yield ~10⁻¹² over 100 orbits, not ~10⁻⁸. Energy conserved at 11 ULPs (the IAS15 floor) rules out a force-evaluation discrepancy: were the per-step accelerations drifting, energy would track. The signature is consistent only with *step-size scheduling differing across platforms*: same forces, same integrator algebra, different `dt` sequence, integrable trajectory phase-drifts while staying on the same ellipse and conserving its energy.

The adaptive step-size optimizer is the only place IAS15 makes a scheduling decision, and `crates/apsis/src/physics/integrator/ias15.rs:1981` was the only libc transcendental call in that path. The formula is Rein & Spiegel (2015) §2.3 eq. (11), the dimensionless 7th-root rule `dt_required = dt_trial · (ε_b / b̃_6)^(1/7)`:

```rust
fn optimal_dt(&self, dt_current: f64, err: f64) -> f64 {
    if err <= 0.0 { return dt_current * DT_ZERO_ERR_GROWTH; }
    let ratio = (self.epsilon / err).powf(1.0 / 7.0);
    dt_current * DT_SAFETY * ratio
}
```

`f64::powf` resolves to a libc `pow` call. Glibc's `pow` is correctly rounded since 2.18; MSVC's UCRT `pow` is not strictly required to be, and empirically the two implementations differ in the last ULP on the inputs encountered. Each step's `dt` therefore differs by ~1 ULP between platforms, the schedule diverges over thousands of steps, and the trajectory phase-drifts.

The diagnostic phase explicitly ruled out other candidates:

- `f64::mul_add` (FMA): grep across the workspace returned zero call sites.
- `target-cpu = native`: not set anywhere in the workspace (default codegen is `target-cpu = x86_64` baseline, no microarch-specific instruction selection).
- SIMD dispatch (`gravity/simd.rs` AVX2 path): gated on `kernel.is_plummer()`, not exercised by Newton-kernel parity scenarios.
- `sqrt`: hardware instruction (`vsqrtsd`), bit-identical by IEEE 754 spec.
- Other transcendentals in the IAS15 hot path: grep returned none.

---

## Phase 2 — intervention

A single line in `optimal_dt` was changed to call the `libm` crate's pure-Rust `pow` instead of `f64::powf`:

```rust
let ratio = libm::pow(self.epsilon / err, 1.0 / 7.0);
```

`libm` was added to `[workspace.dependencies]` as `libm = "0.2"` and to `crates/apsis/Cargo.toml` as `libm = { workspace = true }`. No other source code changed. The full workspace test suite (610 unit tests across the 4 default-member crates) was re-run; all passed.

### Intra-platform sanity check

Before re-running cross-platform, the Windows side regenerated all five outputs with the fix in place and was compared against the pre-fix Windows baseline (UCRT `pow` vs libm `pow` on the same AMD host). The change in trajectories was ~1e⁻⁸ absolute over 100 Kepler orbits — of the same order as the cross-platform drift seen in Phase 1, confirming that `powf` was the dominant scheduling lever, not a minor contributor. Mercury rate changed from 42.993 to 42.991 arcsec/century on Windows — the same value Linux produced in Phase 1 with glibc `pow`, suggesting that libm crate's `pow` and glibc's `pow` agree at least to the 3-decimal precision of the Mercury gate on these inputs.

---

## Phase 3 — post-fix cross-platform run

Both hosts executed against branch `experiment/cross-platform-libm` (commit `e06ba47`), with `rustc 1.94.1` and the same `Cargo.lock` (the only delta from Phase 1 is the `libm` workspace entry, which does not change any other crate's resolved versions). Output files were compared via `python validation/cross-platform/compare.py` and independently via SHA256.

### Per-column ULP results (`compare.py`)

| Scenario | bit-equal rows / total | max ULP any column |
| --- | --- | --- |
| kepler | 101 / 101 | 0 |
| figure8 | 2001 / 2001 | 0 |
| pythagorean | 2101 / 2101 | 0 |
| retrograde | 10001 / 10001 | 0 |

Every column of every row of every CSV matched bit-for-bit between Windows AMD and Linux Intel.

### Independent SHA256 of full output files

| File | Windows SHA256 (first 16) | Linux SHA256 (first 16) | Match |
| --- | --- | --- | --- |
| kepler.csv | `3C597636FBAF5456…` | `3C597636FBAF5456…` | yes |
| figure8.csv | `A0EFE225475C4092…` | `A0EFE225475C4092…` | yes |
| pythagorean.csv | `87BC699265506400…` | `87BC699265506400…` | yes |
| retrograde.csv | `505C5427EF0C4D4B…` | `505C5427EF0C4D4B…` | yes |

File sizes also match exactly. Three independent signals (per-column ULP, file size, SHA256) converge on the same conclusion.

### Mercury 1PN

| Field | Windows | Linux |
| --- | --- | --- |
| measured Δω (500 orbits) | +2.509906 × 10⁻⁴ rad (+51.7705 arcsec) | identical |
| rate | 42.991 arcsec/century | 42.991 arcsec/century |
| Δ vs Phase 1 Windows (42.993) | −0.002 arcsec/century | n/a |
| Δ between hosts | 0.000 arcsec/century (0 ppm of GR) | — |
| Δ vs GR (43.000) | 9.0 arcsec/century, i.e. 2.1 × 10⁻⁴ relative (the 1PN test-particle truncation error) | same |

The Mercury rate after the fix is bit-identical across platforms. The remaining 2.1 × 10⁻⁴ relative deviation from the 43.0 GR prediction is the test-particle truncation error of the 1PN expansion itself (also present in REBOUND + REBOUNDx 1PN), not a numerical artifact.

---

## Interpretation

The Phase 3 result demonstrates that the apsis-record contract (Cargo.lock + rustc + source SHA + libm crate dependency) is sufficient for byte-identical trajectory reproduction across heterogeneous x86_64 hosts, for the integration regime exercised by the parity portfolio.

The reason this holds, in retrospect, is that every other source of cross-platform non-determinism in this configuration is already controlled by the build setup:

- **Same LLVM IR.** Cargo.lock pins every crate version; rustc 1.94.1 pins the compiler; default codegen pins target-cpu to the x86_64 baseline (no `target-cpu = native`). LLVM emits the same intermediate representation for both targets.
- **Same arithmetic primitives.** IEEE 754 mandates correctly-rounded results for `+ − × ÷ sqrt`. Hardware `sqrtsd` is bit-identical between any IEEE 754 x86_64 CPU.
- **No FMA divergence.** Grep confirmed zero `f64::mul_add` / `fmuladd` call sites. LLVM does not auto-fuse `a*b + c` to FMA under default codegen flags (would require `-Cllvm-args=--fp-contract=fast`).
- **No SIMD reduction order divergence.** The AVX2 leaf-pair kernel in `gravity/simd.rs` is gated on `kernel.is_plummer()`; the Newton-kernel parity scenarios run scalar pairwise. Direct summation over 2–3 bodies has no reduction-order ambiguity.
- **MXCSR defaults are IEEE-correct.** Glibc startup and UCRT startup both leave `MXCSR.DAZ = 0`, `MXCSR.FTZ = 0`.
- **Chaos amplifies differences, not creates them.** Pythagorean is famously sensitive (Lyapunov time ~1 t.u.), but a chaotic system with bit-identical initial conditions and bit-identical update rules produces bit-identical output. With the libm fix removing the only stochastic input, Pythagorean's chaotic trajectory is reproducible.

The `powf` call in the step-size controller was the single point of divergence. Once removed, the entire system is deterministic by construction.

---

## Methodology principle

The result generalizes as a methodological statement, but its application requires a per-integrator audit:

> **Bitwise cross-platform reproducibility is achievable when all libc-bound transcendentals in integration-critical paths are replaced by deterministic alternatives (e.g., the `libm` crate).**

The audit checklist for any new integrator or operator joining the cross-platform validation portfolio:

1. Grep the implementation for `f64::sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `exp`, `ln`, `log`, `powf`, `cbrt`, `sinh`, `cosh`, `tanh`.
2. Any hit in a code path that affects trajectory state (force evaluation, step-size decision, switching boundary, drift solver) must route through `libm` instead of `f64::*`.
3. Hardware-implemented `sqrt`, `recip`, and the basic arithmetic operators do not need replacement.
4. Validate by running the parity scenario both platforms post-change and confirming bit-equal output (SHA256 + per-column ULP).

Pending application of this checklist:

- **WHFast Kepler drift** (`crates/apsis/src/physics/integrator/kepler.rs:36,48`): `sin`, `cos`, `cosh` in the Stumpff-series solver — tracked as **#159**.
- **Mercurius Hill-radius switching** (`crates/apsis/src/physics/integrator/mercurius.rs:312`): `cbrt` in the close-encounter detection threshold — tracked as **#160**.
- **apsis-central force law** (`crates/apsis-central/src/lib.rs:245, 405, 451, 665`): `powf(γ + k)` in the central-force prefactor and potential — tracked as **#161**.

Each requires the analogous one-line replacement plus an extension of the cross-platform parity portfolio with a scenario exercising the integrator or operator. Until all three are completed, the bitwise reproducibility claim scopes only to the configuration tested in this experiment.

---

## Limitations

- **x86_64 only.** ARM and RISC-V hosts have different FP unit microarchitectures and (more importantly) different libc/libm implementations. `libm` crate is target-independent, so the same fix likely extends, but this is untested.
- **Default codegen only.** `cargo build --release` with no `target-cpu` override. Anyone passing `-C target-cpu=native` to one host but not the other voids the guarantee — LLVM will pick microarch-specific instructions that differ across CPUs.
- **Parity portfolio scope.** Five scenarios cover IAS15 + direct summation + Newton kernel + 1PN operator. Other integrators (WHFast, Mercurius, Implicit Midpoint), kernels (Plummer with BH SIMD), and operators (radiation, central force) are not exercised here.
- **Unaudited transcendentals.** The grep was comprehensive for the call sites that affect the parity scenarios, but the workspace contains many `sin`/`cos`/`sqrt`/`powf` calls in display code, IC setup, magnitude conversion, etc. Those do not affect trajectory bits but were not enumerated for this experiment.

---

## References

- Everhart, E. (1985). *An efficient integrator that uses Gauss-Radau spacings.* In Carusi & Valsecchi, eds., *Dynamics of Comets: Their Origin and Evolution*, pp. 185–202. Introduces the 15th-order RADAU integrator on Gauss-Radau spacings that IAS15 builds on; the specific step-size formula used here is Rein & Spiegel's, not Everhart's (R&S 2015 §2.3 explicitly state their controller "is different from and superior to the one proposed by Everhart 1985").
- Rein, H., & Spiegel, D. S. (2015). *IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits.* MNRAS 446, 1424–1437. DOI: [10.1093/mnras/stu2164](https://doi.org/10.1093/mnras/stu2164). §2.3 "Step-size control", eq. (11) is the 7th-root rule `dt_required = dt_trial · (ε_b / b̃_6)^(1/7)` that this experiment exercised.
- IEEE Std 754-2019. *IEEE Standard for Floating-Point Arithmetic.* §5.4.1 specifies correctly-rounded results for basic operations; §9.2 lists recommended-but-not-required correctly-rounded transcendentals.
- `libm` crate: Rust port of MUSL libc's math library. <https://github.com/rust-lang/libm>.
- Apsis-record format: `docs/adr/011-apsis-record.md`.
