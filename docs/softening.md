# Plummer Softening — Technical Reference

**Module:** `crates/apsis/src/domain/body.rs`, `crates/apsis/src/physics/gravity/kernel.rs`

---

## 1. Motivation

In a point-mass N-body simulation, the gravitational potential diverges as two
bodies approach: $\Phi \to -\infty$ as $r \to 0$. This singularity is unphysical
for finite-size objects and causes two numerical problems:

1. **Force spikes** that blow up the integrator, even with very small $\Delta t$.
2. **Artificial scattering** at close encounters, corrupting orbital statistics.

Softening replaces the singular potential with a smooth Plummer kernel that
remains finite at $r = 0$.

---

## 2. Plummer kernel

The softened potential for a pair $(i, j)$ is

$$
\Phi_{ij} = -\,\frac{G \, m_i \, m_j}{\sqrt{r^2 + \epsilon_{ij}^2}}.
$$

The corresponding acceleration on body $i$ from body $j$ is

$$
\mathbf{a}_i = \frac{G \, m_j \, (\mathbf{x}_j - \mathbf{x}_i)}{(r^2 + \epsilon_{ij}^2)^{3/2}},
$$

where $r = |\mathbf{x}_j - \mathbf{x}_i|$ and $\epsilon_{ij}$ is the **effective
pairwise softening length**.

### Pairwise softening

Each body carries an individual softening length $\epsilon_i$. For a pair, the
pairwise softening is combined in quadrature (Plummer-equivalent equal-volume
averaging):

$$
\epsilon_{ij}^2 = \frac{\epsilon_i^2 + \epsilon_j^2}{2}.
$$

Implementation: `crates/apsis/src/physics/gravity/kernel.rs::pair_eps2`.

---

## 3. Per-body default

Each body's softening is derived from its mass via a cubic-root scaling, so
bodies with equal density have softening proportional to their physical radius:

$$
\epsilon_i = \mathrm{EPS\_BASE} \cdot m_i^{1/3}, \qquad \mathrm{EPS\_BASE} = 0.02.
$$

Implementation: `crates/apsis/src/domain/body.rs::default_softening`.

### Physical motivation

For a uniform sphere of density $\rho$ and mass $m$, the radius is
$r \propto m^{1/3}$. Setting $\epsilon \propto m^{1/3}$ ensures the softening
volume is proportional to the body's actual volume, preventing spurious force
suppression on small bodies near large ones.

---

## 4. Global softening scale

The UI exposes a dimensionless scale factor $s$ (default $1.0$) that multiplies
all per-body softenings simultaneously:

$$
\epsilon_i^{(\mathrm{eff})} = \mathrm{EPS\_BASE} \cdot m_i^{1/3} \cdot s.
$$

This is stored in `System::softening_scale` and applied in `set_softening_scale()`.

**Guidelines:**

| Scale | Use case |
| --- | --- |
| $0.1$–$0.5$ | High-accuracy runs, well-separated bodies |
| $1.0$ | Default — balanced accuracy/stability |
| $2$–$5$ | Suppressing close-encounter singularities |
| $10$ | Nearly "fuzzy" gravity — cosmological toy models |

---

## 5. Validity criterion

Softening distorts the true gravitational force. The fractional force error for
a pair at separation $r$ with pairwise softening $\epsilon$ is, to leading order,

$$
\frac{\delta F}{F} \approx \frac{3}{2} \, \left(\frac{\epsilon}{r}\right)^2
\qquad \text{for}\ \epsilon / r \ll 1.
$$

This approximation comes from expanding the Plummer denominator,

$$
(r^2 + \epsilon^2)^{-3/2} = r^{-3} \, \left[1 - \tfrac{3}{2}\left(\tfrac{\epsilon}{r}\right)^2 + \ldots\right].
$$

### Diagnostic: $\epsilon / r_\text{min}$ ratio

The simulator tracks the **minimum pairwise separation** $r_\text{min}$ and the
**maximum effective pairwise softening**
$\epsilon_\text{max} = \max_{ij} \epsilon_{ij}$ at every step (for $N \leq 512$
bodies). Their ratio gives an upper bound on the worst-case force error in the
current configuration:

| $\epsilon_\text{max} / r_\text{min}$ | Force error (leading order) | Severity |
| :---: | :---: | :---: |
| $< 0.1$ | $< 1.5\,\%$ | OK (green) |
| $0.1$ – $0.3$ | $1.5\,\%$ – $14\,\%$ | Warning (yellow) |
| $> 0.3$ | $> 14\,\%$ | Critical (red) |

The warning indicator is shown in the **FORCE ACCURACY** section of the Config
panel in real time.

### When the warning fires

The warning fires whenever two bodies come within roughly $10\times$ their
combined softening radius. Typical causes:

- A highly eccentric orbit with tight periapsis;
- A multi-body close encounter (three-body flyby);
- The softening scale set too large for the system geometry.

### Remedies

1. **Reduce $\epsilon$ scale** — brings softening closer to the physical radius.
   Risk: the integrator may require a smaller $\Delta t$ to stay stable.
2. **Reduce $\Delta t$** — allows the integrator to resolve close encounters
   that were previously being jumped over.
3. **Switch to Yoshida 4th-order** — its superior energy conservation allows a
   larger $\Delta t$ margin before instability, giving more headroom near
   encounters.
4. **Accept the error** — for cosmological or statistical runs where close
   encounters are intentionally suppressed, a large $\epsilon$ is physically
   justified.

---

## 6. Relationship to Barnes–Hut ($\theta$)

Softening and the Barnes–Hut opening angle $\theta$ are independent parameters
that both affect force accuracy:

- $\theta$ controls the *geometric* multipole error (cell size / distance ratio).
- $\epsilon$ controls the *singularity* regularisation error at close range.

For publication-quality runs, both should be tuned:

- $\theta < 0.3$ (approaches exact $O(N^2)$);
- $\epsilon / r_\text{min} < 0.1$ ($< 1.5\,\%$ force error at closest encounter).

---

## 7. Implementation files

| File | Role |
| --- | --- |
| `crates/apsis/src/domain/body.rs` | `EPS_BASE`, `default_softening()`, `Body::softening` field |
| `crates/apsis/src/physics/gravity/kernel.rs` | `pair_eps2()`, `plummer_acc()`, `plummer_phi()` |
| `crates/apsis/src/physics/gravity/engine.rs` | Per-body softening used in Barnes–Hut and exact evaluations |
| `crates/apsis/src/core/system/config.rs` | `set_softening_scale()` |
| `crates/apsis/src/core/system/helpers.rs` | `compute_closeness()` ($r_\text{min}$ tracking) |
| `crates/apsis/src/core/metrics.rs` | `Metrics::r_min`, `Metrics::softening_max` |
| `crates/apsis-app/src/app/panel/tabs/config.rs` | $\epsilon$ scale slider + real-time validity indicator |

---

## 8. Known limitations

1. **$r_\text{min}$ tracking is disabled for $N > 512$.** For large
   asteroid-belt simulations, the $O(N^2)$ scan would dominate the frame
   budget. The diagnostic simply shows no warning for these runs.

2. **The force-error formula is asymptotic.** For $\epsilon / r > 0.5$, the
   leading-order approximation underestimates the actual error. The true
   Plummer factor is $(1 + (\epsilon / r)^2)^{-3/2}$ times the point-mass force.

3. **Softening does not model real physics.** It is a numerical regularisation,
   not a physical collision or tidal model. Merger / fragmentation events
   require separate treatment, currently out of scope.
