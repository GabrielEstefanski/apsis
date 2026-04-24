# Plummer Softening — Technical Reference

**Module:** `crates/gravity-sim-core/src/domain/body.rs`, `crates/gravity-sim-core/src/physics/gravity/kernel.rs`  

---

## 1. Motivation

In a point-mass N-body simulation, the gravitational potential diverges as two
bodies approach: Φ → −∞ as r → 0. This singularity is unphysical for
finite-size objects and causes two numerical problems:

1. **Force spikes** that blow up the integrator, even with very small Δt.
2. **Artificial scattering** at close encounters, corrupting orbital statistics.

Softening replaces the singular potential with a smooth Plummer kernel that
remains finite at r = 0.

---

## 2. Plummer Kernel

The softened potential for a pair (i, j) is:

```
Φ_ij = −G · m_i · m_j / √(r² + ε²_ij)
```

The corresponding acceleration on body i from body j is:

```
a_i = G · m_j · (x_j − x_i) / (r² + ε²_ij)^(3/2)
```

where `r = |x_j − x_i|` and `ε_ij` is the **effective pairwise softening length**.

### Pairwise softening

Each body carries an individual softening length `ε_i`. For a pair, the pairwise
softening is combined in quadrature (Plummer-equivalent equal-volume averaging):

```
ε²_ij = (ε²_i + ε²_j) / 2
```

Implementation: `crates/gravity-sim-core/src/physics/gravity/kernel.rs::pair_eps2`.

---

## 3. Per-body Default

Each body's softening is derived from its mass via a cubic root scaling, so
bodies with equal density have softening proportional to their physical radius:

```
ε_i = EPS_BASE · m_i^(1/3)     with EPS_BASE = 0.02
```

Implementation: `crates/gravity-sim-core/src/domain/body.rs::default_softening`.

### Physical motivation

For a uniform sphere of density ρ and mass m, the radius is r ∝ m^(1/3). Setting
ε ∝ m^(1/3) ensures the softening volume is proportional to the body's actual
volume, preventing spurious force suppression on small bodies near large ones.

---

## 4. Global Softening Scale

The UI exposes a dimensionless scale factor `s` (default 1.0) that multiplies
all per-body softenings simultaneously:

```
ε_i(effective) = EPS_BASE · m_i^(1/3) · s
```

This is stored in `System::softening_scale` and applied in `set_softening_scale()`.

**Guidelines:**

| Scale | Use case |
|-------|----------|
| 0.1–0.5 | High-accuracy runs, well-separated bodies |
| 1.0 | Default — balanced accuracy/stability |
| 2–5 | Suppressing close-encounter singularities |
| 10 | Nearly "fuzzy" gravity — cosmological toy models |

---

## 5. Validity Criterion

Softening distorts the true gravitational force. The fractional force error for a
pair at separation r with pairwise softening ε is, to leading order:

```
δF/F ≈ (3/2) · (ε/r)²      for ε/r << 1
```

This approximation comes from expanding the Plummer denominator:

```
(r² + ε²)^(−3/2) = r^(−3) · [1 − (3/2)(ε/r)² + …]
```

### Diagnostic: ε/r_min ratio

The simulator tracks the **minimum pairwise separation** `r_min` and the
**maximum effective pairwise softening** `ε_max = max_ij(ε_ij)` at every step
(for N ≤ 512 bodies). Their ratio gives an upper bound on the worst-case
force error in the current configuration:

| ε_max / r_min | Force error (leading order) | Severity |
|:---:|:---:|:---:|
| < 0.1 | < 1.5% | OK (green) |
| 0.1 – 0.3 | 1.5% – 14% | Warning (yellow) |
| > 0.3 | > 14% | Critical (red) |

The warning indicator is shown in the **FORCE ACCURACY** section of the Config
panel in real time.

### When the warning fires

The warning fires whenever two bodies come within ~10× their combined softening
radius. Typical causes:

- A highly eccentric orbit with tight periapsis
- A multi-body close encounter (three-body flyby)
- The softening scale is too large for the system geometry

### Remedies

1. **Reduce `ε scale`** — brings softening closer to the physical radius. Risk:
   the integrator may require a smaller Δt to stay stable.
2. **Reduce Δt** — allows the integrator to resolve close encounters that were
   previously being jumped over.
3. **Switch to Yoshida 4th-order** — its superior energy conservation allows a
   larger Δt margin before instability, giving more headroom near encounters.
4. **Accept the error** — for cosmological or statistical runs where close
   encounters are intentionally suppressed, a large ε is physically justified.

---

## 6. Relationship to Barnes–Hut (θ)

Softening and the Barnes–Hut opening angle θ are independent parameters that
both affect force accuracy:

- **θ** controls the *geometric* multipole error (cell size / distance ratio).
- **ε** controls the *singularity* regularization error at close range.

For publication-quality runs, both should be tuned:
- θ < 0.3 (approaches exact O(N²))
- ε/r_min < 0.1 (< 1.5% force error at closest encounter)

---

## 7. Implementation Files

| File | Role |
|------|------|
| `crates/gravity-sim-core/src/domain/body.rs` | `EPS_BASE`, `default_softening()`, `Body::softening` field |
| `crates/gravity-sim-core/src/physics/gravity/kernel.rs` | `pair_eps2()`, `plummer_acc()`, `plummer_phi()` |
| `crates/gravity-sim-core/src/physics/gravity/engine.rs` | Per-body softening used in BH and exact evaluations |
| `crates/gravity-sim-core/src/core/system/config.rs` | `set_softening_scale()` |
| `crates/gravity-sim-core/src/core/system/helpers.rs` | `compute_closeness()` (r_min tracking) |
| `crates/gravity-sim-core/src/core/metrics.rs` | `Metrics::r_min`, `Metrics::softening_max` |
| `crates/gravity-sim-app/src/app/panel/tabs/config.rs` | ε scale slider + real-time validity indicator |

---

## 8. Known Limitations

1. **r_min tracking is disabled for N > 512.** For large asteroid-belt
   simulations, the O(N²) scan would dominate the frame budget. The diagnostic
   simply shows no warning for these runs.

2. **The force error formula is asymptotic.** For ε/r > 0.5, the leading-order
   approximation underestimates the actual error. The true Plummer factor is
   (1 + (ε/r)²)^(−3/2) times the point-mass force.

3. **Softening does not model real physics.** It is a numerical regularization,
   not a physical collision or tidal model. Merger/fragmentation events require
   separate treatment (planned: Phase 1 — Collisions).
