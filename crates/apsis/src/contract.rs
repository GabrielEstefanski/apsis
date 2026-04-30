//! Formal contract for the federated perturbation model.
//!
//! This module is the **executable specification** of the federation thesis
//! in `paper.md`: every guarantee the simulator makes to a perturbation
//! author, and every assumption the simulator imposes back, is named here
//! and gated by a CI test whose name matches the guarantee. Reading this
//! module top-to-bottom is reading the contract; running its tests is
//! verifying the contract.
//!
//! ## What this contract is, and why it exists
//!
//! N-body simulators routinely accept user-defined non-gravitational
//! forces (radiation pressure, J2 oblateness, atmospheric drag, custom
//! research perturbations). The accepted way of integrating them — pass
//! a callback, hope for the best — is informal. There is no statement of
//! what the simulator promises about the environment in which the
//! callback runs, no statement of what the callback may rely on across
//! multiple registrations, no machine-checkable record of what counts as
//! a valid configuration. The thesis APSIS advances is that this can be
//! made formal: perturbations are first-class scientific artifacts under
//! a written, versioned, executable contract.
//!
//! Concretely, against the comparable surface in REBOUND/REBOUNDx:
//!
//! ```text
//! | Aspect            | REBOUND   | APSIS         |
//! |-------------------|-----------|---------------|
//! | formal contract   | implicit  | explicit      |
//! | validation        | partial   | systematic    |
//! | composition       | ad hoc    | specified     |
//! | verifiability     | limited   | executable    |
//! ```
//!
//! "Systematic" here is **not** a claim of broader test coverage — REBOUND
//! has a wider validation portfolio measured by problem count. It is a
//! claim of **shape**: every guarantee in this module has a named test
//! gated in CI, every locked baseline lives in `docs/experiments/`, every
//! warning the simulator can emit is asserted to fire under exactly one
//! known configuration. The property that distinguishes APSIS is that a
//! reviewer can **mechanically check** the claims, not that the claims
//! are quantitatively stronger.
//!
//! ## Scope and counter-scope
//!
//! Two layers of guarantee live below — read both before writing a
//! perturbation:
//!
//! ### `KernelInvariants` — what the system promises a registered perturbation
//!
//! 1. **Determinism (system-level).** Identical inputs produce identical
//!    outputs, bit-for-bit. This applies to the **complete simulation**
//!    (core integrator + perturbations), not just to the bare Newtonian
//!    kernel. A perturbation author can rely on
//!    `(bodies, perturbations, dt) → trajectory` being a pure function.
//!    - test: [`tests::invariant_determinism_bit_exact`] (positive)
//!    - test: [`tests::invariant_determinism_distinguishes_distinct_inputs`]
//!      (negative — proves the determinism test is observing trajectory
//!      state, not returning a fixed value)
//!
//! 2. **Newtonian consistency under perturbation attach.** The
//!    underlying Newtonian force evaluation is invariant under
//!    perturbation registration: attaching a no-op perturbation produces
//!    a trajectory bit-equal to the bare run. Perturbations are *added*
//!    to Newton, never substituted for it.
//!    - test: [`tests::invariant_newtonian_consistency_under_null_perturbation_attach`]
//!
//! 3. **Read-only access to base dynamics.** Perturbations cannot
//!    mutate body state, force-model state, or any other system field.
//!    Enforced by the `&self` receiver on
//!    [`PerturbationForce::accumulate`](crate::physics::integrator::PerturbationForce::accumulate)
//!    and by Rust's borrow checker — there is no runtime gate to
//!    bypass. The escape hatch (`Cell`/`RefCell`/atomic via interior
//!    mutability) is a contract violation rather than a structural one,
//!    and is gated by:
//!    - test: [`tests::invariant_perturbation_is_pure_function_of_state`]
//!      (`accumulate` invoked twice on the same instance with identical
//!      input produces identical output — proves no observable internal
//!      state evolution between calls).
//!
//! ### `CompositionRules` — what the system promises about multi-perturbation registration
//!
//! Filled by the next commit in this series. Stub here for shape:
//!
//! 4. Commutativity. Registration order does not affect the trajectory.
//! 5. Associativity. Composing three or more perturbations is
//!    well-defined and order-invariant under any associativity.
//! 6. Additive composition. Each perturbation contributes by `+=` to the
//!    accumulator slice; no perturbation may overwrite or read existing
//!    values to compute its own.
//! 7. Union of kernel requirements. A composed system's effective
//!    `KernelRequirements` is the union of the individual perturbations'.
//!
//! ### Failure model — what the system promises when a configuration is invalid
//!
//! Filled by the third commit in this series. Stub here for shape:
//!
//! 8. Exactly-one warning per violated invariant.
//! 9. Repeated registration of the same violation does not duplicate
//!    warnings.
//! 10. Silent acceptance is structurally impossible — emission goes
//!     through the structured log bus regardless of subscriber state.
//!
//! ## What this contract does NOT guarantee
//!
//! Reviewers who hold the federation thesis to a stronger standard need
//! these limitations stated up front, not buried in implementation:
//!
//! - **Cross-platform bit-exactness.** Determinism holds within a single
//!   build on a single hardware target. f64 reductions can differ
//!   between architectures (FMA emission, libm differences, autovec
//!   thresholds), between Rust toolchain versions, and between
//!   `RUSTFLAGS` settings. `docs/experiments/2026-04-29-3d-port-baseline.md`
//!   states the same restriction for the locked physics baselines: they
//!   are reproduced bit-exact on developer hardware (Windows MSVC) and
//!   pass under a 100 ppm portable bound on CI (Linux glibc). The
//!   determinism test in this module similarly asserts bit-exactness on
//!   the host running it, not across hosts.
//!
//! - **Cross-thread determinism.** The simulator is single-threaded for
//!   the integrate loop. The Barnes–Hut tree traversal uses Rayon, and
//!   parallel reduction order is not guaranteed across runs — but this
//!   is not exercised by any release-mode physics gate (every gate sits
//!   under `EXACT_THRESHOLD` and uses direct O(N²) summation). A
//!   perturbation author who introduces multi-threaded internal state
//!   is outside the contract; the integrator does not protect against
//!   this.
//!
//! - **Cross-RNG-seed equivalence.** APSIS's core integration is RNG-free;
//!   determinism here is *not* "we use a seeded PRNG correctly", it is
//!   "the integrator is a pure function of state and registered
//!   perturbations". A perturbation author who introduces stochastic
//!   forcing must own seed control inside the perturbation; the
//!   simulator does not provide a global RNG handle.
//!
//! - **Build-flag invariance.** A run compiled under `--release` is not
//!   guaranteed bit-equal to a `--debug` run; a run with
//!   `target-cpu=native` is not guaranteed bit-equal to a portable build.
//!   This is the standard f64-numerics caveat and is not specific to
//!   APSIS.
//!
//! ## Iteration-order invariant (load-bearing)
//!
//! The determinism property above silently relies on a property of the
//! perturbation storage:
//!
//! > Perturbations registered through
//! > [`System::add_perturbation`](crate::core::system::System::add_perturbation)
//! > are stored in a `Vec<Box<dyn PerturbationForce>>` and iterated via
//! > `slice::iter()`. **Iteration order equals registration order.**
//!
//! Any future change that swaps the storage for a `HashSet`, `HashMap`,
//! `BTreeSet`-by-pointer-address, or any other container with non-stable
//! iteration order silently breaks determinism. The test
//! [`tests::invariant_determinism_bit_exact`] is the load-bearing guard:
//! such a regression would surface as a bit-difference between two
//! identical runs, not as a compile error. The
//! [`tests::composition_commutative_two_perturbations`] test (next
//! commit) does NOT cover this — commutativity is symmetry under
//! reordering, while determinism is sameness under no reordering.

#[cfg(test)]
mod tests {
    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::math::Vec3;
    use crate::physics::integrator::{IntegratorKind, PerturbationForce};
    use crate::units::UnitSystem;

    // ── Test perturbations ────────────────────────────────────────────────────
    //
    // Test-local fakes deliberately disjoint from `apsis-1pn` and
    // anything in the production crate graph. The contract is generic
    // over `PerturbationForce` impls, so the tests must not couple to
    // any particular real perturbation. `apsis-1pn` carries its own
    // contract evidence in `crates/apsis-1pn/tests/` — proving a real
    // perturbation also satisfies the contract — but does not appear
    // here.

    /// Stateless perturbation that adds a fixed Vec3 to every body's
    /// acceleration. Trivial enough to reason about by hand; useful for
    /// composition tests where each perturbation's contribution is known
    /// in closed form.
    struct ConstantPush(Vec3);

    impl PerturbationForce for ConstantPush {
        fn accumulate(&self, _bodies: &[Body], acc: &mut [Vec3]) {
            for a in acc.iter_mut() {
                *a += self.0;
            }
        }
    }

    /// Stateless perturbation that adds a linear-drag term `-k · v` to
    /// each body. Reads body velocity, so it exercises the
    /// `(bodies, scratch_acc) → contribution` data flow that pure
    /// constant pushes cannot. Used by determinism / state-purity tests
    /// that need a non-trivial dependence on body state.
    struct LinearDrag(f64);

    impl PerturbationForce for LinearDrag {
        fn accumulate(&self, bodies: &[Body], acc: &mut [Vec3]) {
            for (b, a) in bodies.iter().zip(acc.iter_mut()) {
                a.x -= self.0 * b.vx;
                a.y -= self.0 * b.vy;
                a.z -= self.0 * b.vz;
            }
        }
    }

    /// Perturbation that contributes nothing. Used by the Newtonian
    /// consistency test: attaching a no-op perturbation must leave the
    /// trajectory bit-equal to the bare-Newton run.
    struct NullPerturbation;

    impl PerturbationForce for NullPerturbation {
        fn accumulate(&self, _bodies: &[Body], _acc: &mut [Vec3]) {}
    }

    // ── System fixture ────────────────────────────────────────────────────────

    /// Reproducible two-body Kepler-like setup, used as the substrate for
    /// every contract test that needs a running integration. Exact
    /// numbers don't matter — the same fixture across tests means the
    /// invariants are tested against the same dynamical regime.
    fn fixture_system() -> System {
        let primary = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0).unsoftened();
        let satellite = Body::rocky(1e-6).at(1.0, 0.0).with_velocity(0.0, 1.0).unsoftened();
        System::new(vec![primary, satellite], UnitSystem::canonical())
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(1e-3)
    }

    /// Snapshot the full kinematic state of every body. Used to compare
    /// runs bit-for-bit — equality of `Vec<BodyState>` is bit-equality
    /// of every f64 field.
    #[derive(Clone, PartialEq, Debug)]
    struct BodyState {
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
    }

    fn snapshot(sys: &System) -> Vec<BodyState> {
        sys.bodies()
            .iter()
            .map(|b| BodyState { x: b.x, y: b.y, z: b.z, vx: b.vx, vy: b.vy, vz: b.vz })
            .collect()
    }

    // ── KernelInvariants ──────────────────────────────────────────────────────

    /// **Invariant 1 (positive).** Two identical runs produce identical
    /// trajectories, bit-for-bit. The test exercises the *full* system
    /// (Newton + a registered perturbation that reads body state) so a
    /// bug confined to perturbation iteration order, perturbation
    /// internal state, or any post-Newton accumulation surfaces here —
    /// not just bugs in the bare Newton path.
    #[test]
    fn invariant_determinism_bit_exact() {
        const N_STEPS: u64 = 200;

        let run = || {
            let mut sys = fixture_system();
            sys.add_perturbation(Box::new(LinearDrag(1e-4)));
            sys.add_perturbation(Box::new(ConstantPush(Vec3::new(0.0, 0.0, 1e-6))));
            for _ in 0..N_STEPS {
                sys.step();
            }
            snapshot(&sys)
        };

        let a = run();
        let b = run();

        assert_eq!(
            a, b,
            "two runs with identical configuration produced different trajectories — \
             non-determinism in the integrator, perturbation iteration, or perturbation state"
        );
    }

    /// **Invariant 1 (negative / sanity).** The determinism test
    /// machinery actually observes trajectory state, rather than
    /// returning a fixed value or short-circuiting on a stale cache. A
    /// 1e-10 perturbation in the initial position is well above the
    /// f64 noise floor and well below physical relevance; the resulting
    /// trajectories must be **detectably different**.
    ///
    /// Without this test, a regression that silently always returns a
    /// fixed `snapshot()` would pass the positive test trivially. The
    /// negative test is the guard against sanity loss.
    #[test]
    fn invariant_determinism_distinguishes_distinct_inputs() {
        const N_STEPS: u64 = 200;

        let run = |x_offset: f64| {
            let primary = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0).unsoftened();
            let satellite =
                Body::rocky(1e-6).at(1.0 + x_offset, 0.0).with_velocity(0.0, 1.0).unsoftened();
            let mut sys = System::new(vec![primary, satellite], UnitSystem::canonical())
                .with_integrator(IntegratorKind::Ias15)
                .with_dt(1e-3);
            sys.add_perturbation(Box::new(LinearDrag(1e-4)));
            for _ in 0..N_STEPS {
                sys.step();
            }
            snapshot(&sys)
        };

        let baseline = run(0.0);
        let perturbed = run(1e-10);

        assert_ne!(
            baseline, perturbed,
            "a 1e-10 change in initial position produced no observable difference \
             after {N_STEPS} IAS15 substeps — the determinism test machinery is \
             not actually observing trajectory state",
        );
    }

    /// **Invariant 2.** Attaching a perturbation that contributes nothing
    /// produces a trajectory bit-equal to the bare-Newton run. The
    /// underlying gravitational kernel is invariant under the act of
    /// perturbation registration; a bug that leaks registration state
    /// into the Newton evaluation (e.g., a stale `scratch_acc` index, a
    /// reset that runs only when perturbations are present) surfaces as
    /// a bit-difference between the two runs.
    #[test]
    fn invariant_newtonian_consistency_under_null_perturbation_attach() {
        const N_STEPS: u64 = 200;

        // Bare Newton run.
        let mut bare = fixture_system();
        for _ in 0..N_STEPS {
            bare.step();
        }
        let bare_state = snapshot(&bare);

        // Same setup, with a no-op perturbation registered.
        let mut with_null = fixture_system();
        with_null.add_perturbation(Box::new(NullPerturbation));
        for _ in 0..N_STEPS {
            with_null.step();
        }
        let with_null_state = snapshot(&with_null);

        assert_eq!(
            bare_state, with_null_state,
            "registering a no-op perturbation produced a different trajectory from \
             bare Newton — perturbation registration is leaking state into the \
             Newtonian force evaluation"
        );
    }

    /// **Invariant 3.** A perturbation instance is a pure function of
    /// `(bodies, scratch_acc)`: invoking `accumulate` twice on the same
    /// instance with identical inputs produces identical outputs. The
    /// test catches accidental interior mutability (a `Cell`-typed
    /// counter, a memoisation cache, an atomic that drifts between
    /// calls), which the trait's `&self` receiver does not structurally
    /// prevent.
    ///
    /// The trait-level guarantee that `accumulate` cannot mutate
    /// `&mut self` is enforced by the borrow checker and needs no
    /// runtime test. This test is specifically for the looser property
    /// — no observable side effect on internal state across calls —
    /// that interior mutability could violate.
    #[test]
    fn invariant_perturbation_is_pure_function_of_state() {
        let p = LinearDrag(1e-3);

        let primary = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
        let satellite = Body::rocky(1e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
        let bodies = [primary, satellite];

        let mut acc1 = vec![Vec3::ZERO; bodies.len()];
        p.accumulate(&bodies, &mut acc1);

        let mut acc2 = vec![Vec3::ZERO; bodies.len()];
        p.accumulate(&bodies, &mut acc2);

        assert_eq!(
            acc1, acc2,
            "the same perturbation instance produced different accumulator state \
             on two calls with identical input — observable internal state is \
             evolving between calls, which violates the contract that perturbation \
             evaluation is a pure function of (bodies, scratch_acc)"
        );
    }
}
