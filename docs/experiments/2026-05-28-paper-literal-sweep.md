# Paper literal sweep — 2026-05-28

**Status:** done (Phase 1+2). Phase 3 (cross-platform §3.5 percentages)
deferred — requires re-running the controller-`pow` ULP analysis on
both Windows and Linux against a fresh Mercury 1PN trajectory.

## Motivation

`paper.md` cites concrete measured numbers in §3.1–§3.4 and Appendix A.
The Tier-1 derivation scaffolding (`feat/tier1-derivations`) surfaced
that one of those literals (the Plummer-violated drift) no longer
reproduces on a fresh run. This sweep audits every other paper-cited
literal against the current code.

## Method

For each literal, identify the harness/script/data file that produces
it, re-run or re-read, capture the fresh value, diff against the
paper text. Volatile fields (`git_sha`, `created_utc`, `crate_hash`,
`rustc_version`) excluded.

## Results

| § | Literal | Paper | Current | Status |
|---|---|---|---|---|
| Abs / 3.2 | Mercury vs GR rel agreement | $2.8\times10^{-5}$ | $2.802\times10^{-5}$ | ✓ |
| 3.2 | measured cumulative $\Delta\omega$ (500 orbits) | $51.7705$ arcsec | $+51.7705$ arcsec | ✓ |
| 3.2 | GR closed-form $\Delta\omega$ | $51.7720$ arcsec | $+51.7720$ arcsec | ✓ |
| 3.2 | per-century rate | $42.991$ arcsec | $42.991$ arcsec | ✓ |
| 3.2 | violated-Plummer drift ($\varepsilon\approx0.02$ AU) | $-83\,128$ arcsec/century | $-2.289\times 10^6$ (post-unwrap) | **superseded** — replaced by closed-form $-2.349\times 10^6$ |
| 3.2 | 4153-orbit residual | $2.2\times10^{-4}$ | $2.177\times10^{-4}$ | ✓ |
| 3.3 | $R_c$ crossings | 11 | 11 | ✓ |
| 3.3 | spike magnitude range | $[4.7\times10^{-6},\,2.0\times10^{-4}]$ | $[4.664\times10^{-6},\,2.013\times10^{-4}]$ | ✓ |
| 3.3 | smooth-kernel baseline floor | $<2.7\times10^{-14}$ | $2.665\times10^{-14}$ | ✓ |
| 3.1 | radiation Burns 1979 rel | $0.7\,\%$ | $1.19\,\%$ | **resolved** — paper updated to 1.2 %; root cause is #133 back-reaction suppression (see Findings) |
| 3.1 | central round-trip rel | $2.7\,\%$ | $2.65\,\%$ | ✓ |
| 3.1 | Kepler baseline drift (no operator) | $<10^{-7}$ | measured ω̇ ~ 2e-17; gate tightened to $<10^{-9}$ | ✓ (paper + gate aligned at 1e-9) |
| 3.4 | Pythagorean $\|\Delta E\|/\|E_0\|$ at $T=70$ | $1.4\times10^{-10}$ | $1.405\times10^{-10}$ | ✓ |
| 3.4 | retrograde $10^4$ orbits $\|\Delta E\|/\|E_0\|$ | $2.6\times10^{-14}$ | $2.583\times10^{-14}$ | ✓ |
| App A | `[unit_system]` `density = "Msun/AU3"` line | present | absent in dump | **stale (schema drift — field removed)** |
| 3.5 | UCRT 96.97 % / libm 95.29 % / $4.4\times10^{-6}$ / $2.8\times10^{-5}$ | — | deferred | unverified (Phase 3) |

## Findings

**Confirmed stale (need paper update):**

1. §3.2 — `-83,128 arcsec/century` → `-136,732`. Root cause is a
   refactor that removed per-body softening: the original measurement
   was with `EPS_BASE = 0.02` per-body and pair-averaged
   ($\varepsilon_\text{eff} \approx 0.0141$ AU); the current API has
   only flat `NewtonKernel::new(0.02)`. Qualitative claim ("three
   orders of magnitude larger than GR, wrong sign") still holds at
   $-136\,732$ vs $+43$ (factor $\sim 3\,180$, wrong sign).
2. Appendix A — `density = "Msun/AU3"` line in `[unit_system]` no
   longer emitted by the record header. The schema dropped this
   field at some point; the embedded example wasn't refreshed.

**Drifted but root cause now known:**

3. §3.1 — radiation Burns 1979 agreement `0.7 %` → `1.19 %`. Root
   cause is PR #133 (test-particle back-reaction suppression at
   extreme mass ratios): the pre-fix Sun absorbed spurious velocity
   from the m=1e-15 dust grain, which partially cancelled the
   constant-r analytic bias. After #133 the Sun is physically fixed
   (correct), and the residual reflects only the constant-r
   approximation — well within the 5 % gate. Paper abstract and
   §3.1 updated to 1.2 %.

**Tighter than paper claims (could tighten or leave conservative):**

4. §3.1 — Kepler baseline drift gate asserts `< 1e-9` (file:
   `crates/apsis-central/tests/round_trip_gate.rs`, function
   `keplerian_baseline_does_not_precess`) but paper says `< 1e-7`.
   Paper is safe but room to tighten.

**Unchanged (9 literals):** §3.2 satisfied case (4 values), §3.2
long-horizon residual, §3.3 (3 values), §3.4 (2 values), §3.1 central
round-trip — all reproduce within rounding.

**Deferred (Phase 3):** §3.5 cross-platform UCRT 96.97 / libm 95.29
percentages depend on re-running the controller-pow oracle analysis
on both Windows and Linux against a fresh Mercury 1PN trajectory.
Inputs (Mercury satisfied case) reproduce exactly per the §3.2 row,
so the underlying trajectory hasn't drifted; the percentages should
hold but are not confirmed.

## Action items

- [x] §3.2 violated case: superseded — the literal-replacement plan was
      dropped in favour of stating the closed-form prediction
      $-2.35\times 10^6$ arcsec/century (the Plummer derivation
      notebook supplies the formula and 5% gate). The aliased
      $-83,128$ / $-136,732$ literals are removed entirely.
- [x] §Future work *Theory-confirmed counter-tests*: closed; both
      §3.2 (Plummer) and §3.3 (continuity) now have closed-form
      backing notebooks.
- [x] Remove `density = "Msun/AU3"` line from Appendix A TOML block —
      done by hand-edit, matches current schema.
- [x] §3.1 radiation Burns: investigated. Root cause is #133
      (test-particle back-reaction suppression at extreme mass ratios):
      the pre-fix Sun accumulated spurious velocity from the dust
      (m_ratio 1e-15) which partially cancelled the constant-r
      analytic bias. After #133 the measurement is physically cleaner;
      1.19 % is the honest residual. Paper abstract and §3.1 updated
      to 1.2 %.
- [x] §3.1 Kepler baseline tightened from `< 10⁻⁷` to `< 10⁻⁹`
      (paper + gate both updated). Measured ω̇ is 2.1e-17, near the
      ULP floor — gate retains ~8 orders of headroom against noise.
      Sweep's "gated at < 10⁻⁹" claim was incorrect at the time
      (gate was 1e-7); it is now true.
- [ ] Phase 3 follow-up: re-run cross-platform ULP analysis once both
      OS images are available; refresh §3.5 percentages.

## Process gap

The general pattern surfaced is that paper-cited literals have no
automated check against current code output. Two of them
(`-83,128`, the Appendix A density unit) silently drifted past
refactors that removed the upstream feature. Future work item: a CI
gate that re-runs each cited harness, captures its literal, and
fails if the value in `paper.md` differs by more than rounding
tolerance.
