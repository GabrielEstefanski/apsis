---
name: Bug report
about: Report a defect — incorrect behavior, regression, numerical drift, build failure
title: ""
labels: bug
---

## What happened

A clear description of the observed behavior.

## What you expected

The behavior the documentation, paper, or code contract leads you to expect.

## Reproduction

Minimal steps to trigger the issue:

```text
commit:      <git rev-parse HEAD>
OS / arch:   <linux x86_64 / macOS arm64 / windows x86_64>
toolchain:   <rustc --version>
integrator:  <ias15 / yoshida4 / mercurius / whfast / implicit_midpoint / verlet>
scenario:    <preset name or short IC description>
```

```rust
// or python; smallest snippet that triggers the issue
```

## Numerical context (if applicable)

For energy drift, orbital element bounds, parity discrepancies, etc.:

- Measurement: `<value with units>`
- Reference / expected: `<value with units>`
- Tolerance / gate: `<value>`

## Impact

- [ ] CI gate failure
- [ ] Regression vs prior tagged version
- [ ] Paper-relevant (validation portfolio, published claim)
- [ ] Cosmetic / docs / DX
