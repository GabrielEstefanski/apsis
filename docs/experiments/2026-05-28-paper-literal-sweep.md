# Paper literal sweep — 2026-05-28

**Status:** done (Phases 1–3). Phase 3 (cross-platform §3.5
percentages) closed 2026-07-04: re-measured on the current workload,
paper updated. Full table in the controller-`pow` notebook's
re-measurement section.

## Motivation

`paper.md` cites concrete measured numbers in §3.1–§3.4 and Appendix A.
The Tier-1 derivation scaffolding surfaced
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
| Abs / 3.2 | Mercury vs GR rel agreement | $2.8\times10^{-5}$ | $2.802\times10^{-5}$ | pass |
| 3.2 | measured cumulative $\Delta\omega$ (500 orbits) | $51.7705$ arcsec | $+51.7705$ arcsec | pass |
| 3.2 | GR closed-form $\Delta\omega$ | $51.7720$ arcsec | $+51.7720$ arcsec | pass |
| 3.2 | per-century rate | $42.991$ arcsec | $42.991$ arcsec | pass |
| 3.2 | violated-Plummer drift ($\varepsilon\approx0.02$ AU) | $-83\,128$ arcsec/century | $-2.289\times 10^6$ (post-unwrap) | superseded — replaced by closed-form $-2.349\times 10^6$ |
| 3.2 | 4153-orbit residual | $2.2\times10^{-4}$ | $2.177\times10^{-4}$ | pass |
| 3.3 | $R_c$ crossings | 11 | 11 | pass |
| 3.3 | spike magnitude range | $[4.7\times10^{-6},\,2.0\times10^{-4}]$ | $[4.664\times10^{-6},\,2.013\times10^{-4}]$ | pass |
| 3.3 | smooth-kernel baseline floor | $<2.7\times10^{-14}$ | $2.665\times10^{-14}$ | pass |
| 3.1 | radiation Burns 1979 rel | $0.7\,\%$ | $1.19\,\%$ | resolved — paper updated to 1.2 %; root cause is back-reaction suppression at extreme mass ratios (see Findings) |
| 3.1 | central round-trip rel | $2.7\,\%$ | $2.65\,\%$ | pass |
| 3.1 | Kepler baseline drift (no operator) | $<10^{-7}$ | measured ω̇ ~ 2e-17; gate tightened to $<10^{-9}$ | pass (paper + gate aligned at 1e-9) |
| 3.4 | Pythagorean $\|\Delta E\|/\|E_0\|$ at $T=70$ | $1.4\times10^{-10}$ | $1.405\times10^{-10}$ | pass |
| 3.4 | retrograde $10^4$ orbits $\|\Delta E\|/\|E_0\|$ | $2.6\times10^{-14}$ | $2.583\times10^{-14}$ | pass |
| App A | `[unit_system]` `density = "Msun/AU3"` line | present | absent in dump | **stale (schema drift — field removed)** |
| 3.5 | UCRT / libm oracle-agreement percentages | 96.97 % / 95.29 % (42,662 inputs) | 96.69 % / 95.08 % (43,284 inputs) | resolved — paper updated 2026-07-04; the $4.4\times10^{-6}$ / $2.8\times10^{-5}$ trajectory literals left §3.5 in the 2026-06-16 rewrite |

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
   cause is test-particle back-reaction suppression at
   extreme mass ratios: the pre-fix Sun absorbed spurious velocity
   from the m=1e-15 dust grain, which partially cancelled the
   constant-r analytic bias. After the fix the Sun is physically fixed
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

**Phase 3 (closed 2026-07-04):** re-measured on the current workload
(43,284 unique controller inputs; the count moved from 42,662 with
ADR-014/ADR-015 and the controller-policy alignment): UCRT 96.69 %,
libm 95.08 %, glibc 2.39 96.68 %; no implementation off by more than
1 ULP on any input; libm bit-identical Windows↔Linux on the full set.
Qualitative structure of the original measurement unchanged. Paper
literals updated (count + both percentages, two call-out sites).

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
- [x] §3.1 radiation Burns: investigated. Root cause is
      test-particle back-reaction suppression at extreme mass ratios:
      the pre-fix Sun accumulated spurious velocity from the dust
      (m_ratio 1e-15) which partially cancelled the constant-r
      analytic bias. After the fix the measurement is physically cleaner;
      1.19 % is the honest residual. Paper abstract and §3.1 updated
      to 1.2 %.
- [x] §3.1 Kepler baseline tightened from `< 10⁻⁷` to `< 10⁻⁹`
      (paper + gate both updated). Measured ω̇ is 2.1e-17, near the
      ULP floor — gate retains ~8 orders of headroom against noise.
      Sweep's "gated at < 10⁻⁹" claim was incorrect at the time
      (gate was 1e-7); it is now true.
- [x] Phase 3 follow-up: re-run cross-platform ULP analysis once both
      OS images are available; refresh §3.5 percentages. Done
      2026-07-04 (Windows UCRT + WSL2 glibc 2.39); paper updated.

## Process gap

The general pattern surfaced is that paper-cited literals have no
automated check against current code output. Two of them
(`-83,128`, the Appendix A density unit) silently drifted past
refactors that removed the upstream feature. Future work item: a CI
gate that re-runs each cited harness, captures its literal, and
fails if the value in `paper.md` differs by more than rounding
tolerance.
