//! Registration-only Exactness assertion for the Plummer cluster protocol
//! (`paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md`
//! §Contract assertion): registering the exactness-requiring 1PN operator
//! on the softened cluster kernel emits exactly one Exactness diagnostic.
//! No integration is performed — integrating 1PN dynamics in cluster
//! units is out of the protocol's scope.

use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

use apsis::core::log::{Event, Level, subscribe, unsubscribe};
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::NewtonKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

const MARKER: &str = "perturbation requires exact 1/r gravity";

#[test]
fn softened_cluster_kernel_under_1pn_emits_one_exactness_diagnostic() {
    // Protocol ε for N = 10³: 0.98 · N^(−0.26) Plummer scale lengths.
    let eps = 0.98 * 1000.0_f64.powf(-0.26) * (3.0 * PI / 16.0);

    let captured: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = captured.clone();
    let id = subscribe(move |event: &Event| {
        if event.message.starts_with(MARKER) {
            sink.lock().unwrap().push(event.clone());
        }
    });

    let bodies = vec![
        Body::rocky(1.0e-3).at_3d(0.5, 0.0, 0.0).with_velocity_3d(0.0, 0.3, 0.0),
        Body::rocky(1.0e-3).at_3d(-0.5, 0.0, 0.0).with_velocity_3d(0.0, -0.3, 0.0),
        Body::rocky(1.0e-3).at_3d(0.0, 0.7, 0.1).with_velocity_3d(0.2, 0.0, 0.0),
    ];
    let mut sys = System::new(bodies, UnitSystem::canonical())
        .with_kernel(Arc::new(NewtonKernel::new(eps)))
        .with_integrator(IntegratorKind::Ias15);
    sys.add_hamiltonian_perturbation(Box::new(
        PostNewtonian1PN::for_units(UnitSystem::canonical()),
    ))
    .expect("registration itself must succeed; the violation is a warning, not an error");

    let events = captured.lock().unwrap().clone();
    unsubscribe(id);

    assert_eq!(events.len(), 1, "exactly one Exactness warning expected");
    assert_eq!(events[0].level, Level::Warn);
    let fields: Vec<&str> = events[0].fields.iter().map(|(k, _)| *k).collect();
    assert!(fields.contains(&"kernel_epsilon"));
    assert!(fields.contains(&"violated_invariant"));
}
