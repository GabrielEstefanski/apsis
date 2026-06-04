# ADR-014 — IAU 2015 Nominal GM_sun as the Calibration Primitive

**Status:** Accepted
**Date:** 2026-06-04

**Supersedes (in part):** the `G_SI · MSUN_KG` calibration in
[`crates/apsis/src/units.rs`](../../crates/apsis/src/units.rs)
(`solar_canonical`, `MSUN_KG`).

---

## Context

`UnitSystem::solar_canonical` derived its time scale — and hence the
implied gravitational parameter of the solar-AU system — from the
product `G_SI · MSUN_KG`, with `MSUN_KG = 1.98892 × 10³⁰ kg` (a
textbook solar mass). That makes the gravitational parameter a product
of two independently-measured, ~5-digit constants:

- CODATA 2018 `G = 6.67430(15) × 10⁻¹¹ m³ kg⁻¹ s⁻²` — relative
  uncertainty ~22 ppm, the least-precise constant in physics.
- A textbook `M_sun` known to ~4–5 significant figures.

The product `G · M_sun` inherits both uncertainties and lands **257 ppm
away** from the directly-measured solar gravitational parameter.

This is backwards from how the quantity is actually known. `(GM)_sun`
is determined to ~10 significant figures from planetary ephemerides;
`G` and `M_sun` individually are not. **IAU 2015 Resolution B3**
therefore defines the *nominal solar mass parameter*
`(GM)_sun^N = 1.3271244 × 10²⁰ m³ s⁻²` as an exact conventional
constant, precisely because `GM` is the well-determined quantity and
`G` is the bottleneck. REBOUND (`units.py`: `yr2pi` from
`GM_sun = 1.3271244004193938 × 10²⁰`), Orekit, and JPL all take `GM`
as the primitive.

apsis's `G × M` construction was non-standard, and the audit (finding
**UM**, `validation/audit/ledger.md`) measured the cost:

- The recommended user-facing path `for_units(solar_canonical).c()`
  returned `10064.03`, **128 ppm** off the community value `10065.32`
  (REBOUND) — a reader arriving from REBOUND/Orekit got the *worst* `c`
  on the path the docs point them to.
- The solar-canonical implied `GM` sat 257 ppm off the IAU nominal.
- The apsis-1pn convention gap (Gaussian `for_units` `c` vs the
  `C_SOLAR_UNITS` IAU-julian literal) read ~110 ppm — dominated by the
  wrong `GM`, masking the genuine julian-vs-Gaussian definitional
  difference of ~19 ppm.

## Decision

Take the IAU 2015 nominal `(GM)_sun = 1.3271244 × 10²⁰ m³ s⁻²` as the
calibration primitive.

```rust
pub const GM_SUN_SI: f64 = 1.327_124_4e20;   // IAU 2015 Resolution B3, exact
pub const MSUN_KG:   f64 = GM_SUN_SI / G_SI;  // derived, not an independent literal
```

- `solar_canonical` time scale `= sqrt(AU³ / GM_SUN_SI)`, so `G_code = 1`
  exactly *and* the implied `GM` equals the IAU nominal.
- `MSUN_KG` becomes a derived convenience. `G_code · M_code = GM` holds
  by construction; the kg mass carries `G`'s uncertainty, where it
  belongs.

Effect: `for_units(solar_canonical).c() = 10065.32`, matching REBOUND on
the recommended path. The apsis-1pn convention gap collapses from ~110
to ~19 ppm — the residual being the genuine, irreducible IAU-julian-
vs-Gaussian gap (both conventions valid).

## Validation cross-check

- **Unit suite:** verdicts unchanged. No gate broke — the existing
  tolerances (`solar_g` 1 %, the old `for_units`-gap `1e-3`) absorbed
  the 257 ppm shift. That non-breakage is itself the finding: those
  gates were too loose to see a units-constant move. Closed by
  tightening the `for_units`-gap gate to `5e-5` and adding two pins —
  `solar_canonical_gm_matches_iau_nominal` (`<1e-9`) and
  `msun_kg_is_gm_over_g`.
- **Mercury 1PN gate** (`mercury_precession_gate.rs`): bit-identical
  before/after. It pins `c` to the `C_SOLAR_UNITS` literal (IAU julian,
  `M_sun`-independent) and integrates in code units (`M = 1`), so the
  gated precession does not depend on the kg mass.
- **REBOUND parity (mercurius):** the validation scripts hardcoded the
  old `MSUN_KG`; updated to derive from `GM_SUN_SI` (preserving the
  multiply-then-divide f64 association so the f64 `G` agrees bit-for-bit
  with `solar().g()`). Re-run: Tier-1 conservation parity passes with
  apsis ≡ REBOUND **bit-identical** initial energy; absolute `E₀`
  shifted by the expected 257 ppm (`−4.298790e-3 → −4.297688e-3`).
  Tier-2 (chaotic test-particle elements) fails exactly as the protocol
  notebook already documents — unchanged by this ADR.

## Alternatives rejected

- **Keep `G × M`, just update `M_sun` to a newer textbook value.**
  Any `M_sun` literal still multiplies in `G`'s 22 ppm uncertainty and
  cannot reach the 10-digit `GM`. The problem is the *primitive*, not
  the literal.
- **`GM` primitive for dynamics, textbook `M_sun` for the kg export.**
  Splits the solar mass into two inconsistent values (`G_code · M_code ≠
  GM`) — the "strange, sensitive, unexpected-to-outsiders" model this
  ADR removes.

## Consequences

- `mass_to_si(1 M_sun)` now returns `1.98841 × 10³⁰ kg`, ~257 ppm below
  the textbook `1.989 × 10³⁰`. By design: the GM is the measured
  quantity; the kg inherits CODATA `G`'s uncertainty.
- `paper/notebooks/2026-05-13-rebound-parity-mercurius.md` records an
  old-G `E₀` literal; the cross-impl parity number it feeds the paper
  §Validation table is unaffected (energy/Lz agreement is dimensionless).
  Update tracked with that notebook.
- The density-rendering constant `KG_M3_TO_SOLAR_AU3`
  (`templates/builders.rs`) still quotes the textbook solar mass; it is
  a few-percent radius convenience, independent of the dynamics
  primitive. Reconciled separately.

## References

- IAU 2015 Resolution B3, *Recommended Nominal Conversion Constants for
  Selected Solar and Planetary Properties.* arXiv:1510.07674.
- CODATA 2018: `G = 6.67430(15) × 10⁻¹¹ m³ kg⁻¹ s⁻²`.
- REBOUND `units.py` (GM_sun primitive); Rein & Liu 2012, A&A 537, A128.
- Audit finding UM — `validation/audit/ledger.md`.
